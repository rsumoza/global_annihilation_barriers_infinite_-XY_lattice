#!/usr/bin/env python3
"""Exact arithmetic checks for the restricted infinite initial-state class.

This checks finite block prefixes and the scalar deductions of the manuscript.
It does not execute general critical-state oracles, solve a new global barrier,
or decide whether an arbitrary infinite generator satisfies its input contract.

Input JSON (optional):
  {"base_radius": 4098, "epsilon": "1/1000000", "blocks": [
    {"n": 0, "vertices": [[0, 0, "1/100"]], "norm2": "1/2500"}
  ]}
Blocks must have consecutive indices beginning at zero. All values are exact
integers or rational strings; floating-point JSON numbers are rejected.
The optional norm2 claim is checked against the full-lattice energy norm,
including every bond from the stored support to its zero exterior.
Alternatively give raw_vertices instead of vertices. The raw compact field q_n
is normalized structurally to v_n = w_n q_n / max(1, ||q_n||_X^2), where
w_n = (1/20) 2^(-n-1). Thus its norm budget holds automatically, independently
of the amplitudes of q_n. The support and effective-generation contract remain.
"""

import argparse
import copy
from fractions import Fraction as Q
import json
from pathlib import Path
import sys


BASE_RADIUS = 4098
RADIUS = Q(1, 20)
SIGMA = Q(57, 1000)
TAU = Q(71, 100)
BARRIER_LOWER = Q(567867, 100000000)
BARRIER_UPPER = Q(570116, 100000000)
CRITICAL_NORM_BOUND = Q(6008)
CRITICAL_PROJECTION_CONSTANT = Q(27200000000000)
GLOBAL_CERTIFICATE_SHA256 = (
    "b1cbb7206a0c7492ad6990213a48480c9d25e04a5a21fe1805f0195a871f5999"
)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def integer(value, label):
    require(type(value) is int, label + " must be an integer")
    return value


def rational(value, label):
    require(type(value) is int or isinstance(value, str),
            label + " must be an exact integer or rational string")
    try:
        return Q(value)
    except (ValueError, ZeroDivisionError) as error:
        raise ValueError(label + " is not a valid rational") from error


def full_lattice_norm2(vertices):
    """Each positive-coordinate bond is counted exactly once."""
    anchors = set()
    for x, y in vertices:
        anchors.update(((x, y, 0), (x - 1, y, 0),
                        (x, y, 1), (x, y - 1, 1)))
    total = Q(0)
    for x, y, direction in anchors:
        endpoint = (x + (direction == 0), y + (direction == 1))
        difference = vertices.get(endpoint, Q(0)) - vertices.get((x, y), Q(0))
        total += difference * difference
    return total


def check_prefix(payload):
    require(isinstance(payload, dict), "input must be a JSON object")
    require(set(payload) <= {"base_radius", "epsilon", "blocks"},
            "unknown top-level input field")
    base = integer(payload.get("base_radius", BASE_RADIUS), "base_radius")
    require(base == BASE_RADIUS,
            "the common critical-tail sequence has base radius 4098")
    epsilon = rational(payload.get("epsilon", "1/1000000"), "epsilon")
    require(epsilon > 0, "epsilon must be strictly positive")
    blocks = payload.get("blocks")
    require(isinstance(blocks, list) and len(blocks) > 0,
            "blocks must be a nonempty finite prefix")
    reports = []
    for expected_n, block in enumerate(blocks):
        require(isinstance(block, dict), "each block must be an object")
        require(set(block) <= {"n", "vertices", "raw_vertices", "norm2"},
                "unknown block field")
        n = integer(block.get("n"), "block index")
        require(n == expected_n, "block indices must be consecutive from zero")
        box_radius = base * 2 ** n
        is_raw = "raw_vertices" in block
        require(("vertices" in block) != is_raw,
                "give exactly one of vertices or raw_vertices")
        raw_vertices = block["raw_vertices"] if is_raw else block["vertices"]
        require(isinstance(raw_vertices, list), "vertices must be a list")
        vertices = {}
        for row in raw_vertices:
            require(isinstance(row, list) and len(row) == 3,
                    "each vertex must be [x,y,value]")
            x = integer(row[0], "vertex x")
            y = integer(row[1], "vertex y")
            require(max(abs(x), abs(y)) <= box_radius,
                    "vertex is outside its prescribed support box")
            require((x, y) not in vertices, "duplicate vertex")
            vertices[(x, y)] = rational(row[2], "vertex value")
        norm_bound = RADIUS / 2 ** (n + 1)
        raw_norm2 = full_lattice_norm2(vertices)
        if is_raw:
            divisor = max(Q(1), raw_norm2)
            scale = norm_bound / divisor
            vertices = {vertex: scale * value for vertex, value in vertices.items()}
            norm2 = full_lattice_norm2(vertices)
            require(norm2 == norm_bound ** 2 * raw_norm2 / divisor ** 2,
                    "structural normalization identity failed")
            require(raw_norm2 <= divisor ** 2,
                    "structural normalization does not prove the norm budget")
        else:
            divisor = Q(1)
            norm2 = raw_norm2
        require(norm2 <= norm_bound ** 2, "block exceeds its norm budget")
        if "norm2" in block:
            claimed = rational(block["norm2"], "norm2")
            require(claimed == norm2,
                    "claimed norm2 omits or miscalculates full-lattice bonds")
        reports.append({"n": n, "box_radius": box_radius,
                        "stored_vertices": len(vertices),
                        "norm2": str(norm2),
                        "norm_budget2": str(norm_bound ** 2),
                        "structurally_normalized": is_raw,
                        "raw_norm2": str(raw_norm2),
                        "normalization_divisor": str(divisor)})
    last_n = len(blocks) - 1
    tail = RADIUS / 2 ** (last_n + 1)
    tail_constant = RADIUS * base / 2
    require(tail == tail_constant / (base * 2 ** last_n),
            "geometric tail conversion failed")
    return {
        "status": "FINITE_PREFIX_PASS",
        "base_radius": base,
        "blocks": reports,
        "last_index": last_n,
        "remaining_series_norm_budget": str(tail),
        "perturbation_projection_constant": str(tail_constant),
        "state_projection_constant": str(
            CRITICAL_PROJECTION_CONSTANT + tail_constant),
        "energy_error_budget": check_energy_budget(epsilon),
        "infinite_contract": (
            "Future blocks must also be generated effectively with the stated "
            "support and norm bounds. A checked finite prefix does not prove "
            "that universal property of an arbitrary generator."
        ),
    }


def check_energy_budget(epsilon):
    if not isinstance(epsilon, Q):
        epsilon = rational(epsilon, "epsilon")
    require(epsilon > 0, "epsilon must be strictly positive")
    coefficient = 2 * CRITICAL_NORM_BOUND + RADIUS + 1
    delta = min(Q(1), epsilon / (8 * coefficient))
    compact_approximation_radius = (
        (2 * CRITICAL_NORM_BOUND + RADIUS) * delta + delta ** 2)
    finite_energy_radius_each = epsilon / 16
    barrier_radius = compact_approximation_radius + 2 * finite_energy_radius_each
    require(delta <= 1, "delta must be at most one")
    require(compact_approximation_radius <= epsilon / 8,
            "compact energy errors exceed their allocation")
    require(barrier_radius <= epsilon / 4,
            "barrier radius exceeds its allocation")
    require(2 * barrier_radius <= epsilon / 2,
            "barrier width exceeds its allocation")
    return {
        "epsilon": str(epsilon),
        "compact_state_error_delta": str(delta),
        "compact_energy_error_radius_sum": str(compact_approximation_radius),
        "finite_energy_interval_radius_each": str(finite_energy_radius_each),
        "barrier_interval_radius_upper": str(barrier_radius),
        "barrier_interval_width_upper": str(2 * barrier_radius),
        "scope": "Arithmetic allocation only; critical-state oracles not executed.",
    }


def check_constants():
    segment_ceiling = RADIUS ** 2 / 2
    curvature_floor = SIGMA - TAU * RADIUS
    barrier_floor = BARRIER_LOWER - segment_ceiling
    require(segment_ceiling == Q(1, 800), "wrong segment ceiling")
    require(curvature_floor == Q(43, 2000), "wrong curvature floor")
    require(barrier_floor == Q(442867, 100000000), "wrong barrier floor")
    require(curvature_floor > 0 and barrier_floor > 0,
            "strict class margins must remain positive")
    # tau is a strict rational upper bound on 1/sqrt(2).
    require(2 * TAU ** 2 > 1, "tau must exceed the Hessian Lipschitz constant")
    require(BARRIER_LOWER < BARRIER_UPPER, "invalid seed barrier bracket")
    return {
        "radius": str(RADIUS), "seed_inverse_gap": str(SIGMA),
        "hessian_lipschitz_upper": str(TAU),
        "segment_energy_ceiling_above_minimum": str(segment_ceiling),
        "curvature_floor_on_class_ball": str(curvature_floor),
        "strict_barrier_lower_bound": str(barrier_floor),
        "seed_barrier_lower": str(BARRIER_LOWER),
        "seed_barrier_upper": str(BARRIER_UPPER),
        "critical_norm_bound": str(CRITICAL_NORM_BOUND),
        "critical_projection_constant": str(CRITICAL_PROJECTION_CONSTANT),
        "common_projection_radius_base": BASE_RADIUS,
    }


def demo_payloads():
    single = {"base_radius": BASE_RADIUS, "epsilon": "1/1000000",
              "blocks": [{"n": 0, "vertices": [[0, 0, "1/100"]],
                          "norm2": "1/2500"}]}
    multiscale = {"base_radius": BASE_RADIUS, "epsilon": "1/1000000000000",
                  "blocks": []}
    for n in range(5):
        amplitude = RADIUS / 2 ** (n + 2)
        multiscale["blocks"].append({
            "n": n, "vertices": [[BASE_RADIUS * 2 ** n, 0, str(amplitude)]],
            "norm2": str(4 * amplitude ** 2)})
    translated = {"base_radius": BASE_RADIUS, "epsilon": "1000000",
                  "blocks": [{"n": 0,
                              "vertices": [[-BASE_RADIUS, BASE_RADIUS, "1/100"]],
                              "norm2": "1/2500"}]}
    normalized = {
        "base_radius": BASE_RADIUS, "epsilon": "1/1000000",
        "blocks": [
            {"n": 0, "raw_vertices": [[0, 0, "1000"]]},
            {"n": 1, "raw_vertices": [[BASE_RADIUS * 2, 0, "1/1000"]]},
            {"n": 2, "raw_vertices": [[-BASE_RADIUS * 4, BASE_RADIUS * 4, "1"]]},
            {"n": 3, "raw_vertices": []},
        ]}
    return {"single_compact": single, "multiscale_exact_budget": multiscale,
            "translated_boundary_support": translated,
            "structural_normalization_large_small_and_zero_raw_fields": normalized}


def negative_controls():
    seed = demo_payloads()["single_compact"]
    mutations = []
    data = copy.deepcopy(seed)
    data["blocks"][0]["vertices"][0][2] = "1/10"
    mutations.append(("norm_budget_exceeded", data))
    data = copy.deepcopy(seed)
    data["blocks"][0]["vertices"][0][0] = BASE_RADIUS + 1
    mutations.append(("support_outside_box", data))
    data = copy.deepcopy(seed)
    data["blocks"][0]["norm2"] = "0"
    mutations.append(("omitted_exterior_interface_bonds", data))
    data = copy.deepcopy(seed)
    data["epsilon"] = "0"
    mutations.append(("zero_requested_precision", data))
    data = copy.deepcopy(seed)
    data["blocks"][0]["vertices"].append([0, 0, "1/100"])
    mutations.append(("duplicate_vertex", data))
    data = copy.deepcopy(seed)
    data["blocks"][0]["n"] = 1
    mutations.append(("nonconsecutive_block_index", data))
    data = copy.deepcopy(seed)
    data["blocks"][0]["vertices"][0][2] = 0.01
    mutations.append(("unenclosed_float_input", data))
    data = copy.deepcopy(seed)
    data["base_radius"] = 4097
    mutations.append(("wrong_critical_tail_sequence", data))
    results = []
    for name, data in mutations:
        try:
            check_prefix(data)
        except ValueError as error:
            results.append({"name": name, "status": "REJECTED", "reason": str(error)})
        else:
            raise ValueError("negative control unexpectedly accepted: " + name)
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, help="optional finite-prefix JSON")
    parser.add_argument("--output", type=Path, help="write JSON report")
    parser.add_argument("--self-test", action="store_true",
                        help="reject deliberately corrupted finite inputs")
    arguments = parser.parse_args()
    report = {
        "status": "INFINITE_CLASS_ARITHMETIC_PASS",
        "scope": (
            "Exact scalar consequences and finite rational block-prefix checks. "
            "No new global supervisor execution, general critical-state oracle, "
            "or arbitrary-precision barrier computation is asserted."
        ),
        "seed_certificate_sha256": GLOBAL_CERTIFICATE_SHA256,
        "seed_binding_scope": "Identity reference only; the seed certificate is not re-proved here.",
        "constants": check_constants(),
    }
    if arguments.input:
        report["input_prefix"] = check_prefix(json.loads(arguments.input.read_text()))
    else:
        report["demonstrations"] = {
            name: check_prefix(payload) for name, payload in demo_payloads().items()}
    if arguments.self_test:
        report["negative_controls"] = negative_controls()
    output = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if arguments.output:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(output)
        print("INFINITE_CLASS_ARITHMETIC_PASS: " + str(arguments.output))
    else:
        sys.stdout.write(output)


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError, json.JSONDecodeError) as error:
        print("INFINITE_CLASS_ARITHMETIC_FAIL: " + str(error), file=sys.stderr)
        sys.exit(1)
