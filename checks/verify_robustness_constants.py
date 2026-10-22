#!/usr/bin/env python3
"""Exact rational checks for the analytic coupling-robustness corollary.

This does not solve perturbed equilibria. It verifies inequalities derived from
the original freshly reproduced stationary-state/global-barrier certificate.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path


def exact_json(path):
    return json.loads(path.read_text(), parse_float=str)


def check(condition, message):
    if not condition:
        raise AssertionError(message)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--certificate-root", type=Path,
        default=Path(__file__).resolve().parents[1] / "proof_record")
    parser.add_argument("--output", type=Path,
        default=Path(__file__).with_name("robustness_constants.json"))
    args = parser.parse_args()
    names = ["certificate.json", "fresh_connection/prerequisites/baseline.json",
             "fresh_connection/prerequisites/saddle.json", "caps.json"]
    records = {name: exact_json(args.certificate_root / name) for name in names}
    global_cert, baseline, saddle = (records[name] for name in names[:3])
    caps=records["caps.json"]
    # Bind every subsidiary record to the authenticated global record.
    for name in names[1:]:
        digest = hashlib.sha256((args.certificate_root/name).read_bytes()).hexdigest()
        check(global_cert["output_sha256"].get(name) == digest,
              "Subsidiary record hash mismatch: " + name)
    summary_path = args.certificate_root / "reproduction_summary.json"
    if summary_path.is_file():
        expected = exact_json(summary_path)["certificate_sha256"]
        check(hashlib.sha256((args.certificate_root/"certificate.json").read_bytes()).hexdigest()
              == expected, "Global certificate differs from reproduction record")
    check(global_cert["status"] == "FRESH_OPTIMAL_ANNIHILATION_BARRIER_PASS", "Global certificate status")
    check(global_cert["fresh_results_only"] is True, "Global certificate must use fresh results")
    check(global_cert["antecedent_execution_mode"] == "fresh_subprocesses", "Fresh antecedents required")
    check(baseline["status"] == "INDEPENDENT_CHAIN_PASS", "Minimum certificate status")
    check(saddle["status"] == "INFINITE_SADDLE_INDEX_ONE_PASS", "Saddle certificate status")
    check(saddle["morse_index"] == 1, "Saddle index")
    sigma, tau, delta, radius = F(57,1000), F(71,100), F(1,500), F(1,32)
    eta, rate, hessian_lipschitz = F(1,100), F(13,500), F(711,1000)
    check(tau*tau > F(1,2), "tau must majorize 1/sqrt(2)")
    min_gap = F(baseline["c"]) - F(baseline["rho"])
    saddle_gap = F(saddle["sigma"]) - F(saddle["rho"])
    branch_margin = F(baseline["branch_margin"]) - F(baseline["rho"])
    check(min_gap > sigma, "Minimum exact gap must exceed sigma")
    check(saddle_gap > sigma, "Saddle exact gap must exceed sigma")
    check(branch_margin > F(11,10), "Minimum exact global branch margin")
    gram=caps["G"]
    off=max(abs(F(gram[0][1][0])),abs(F(gram[0][1][1])))
    gram_min=min(F(gram[i][i][0])-off for i in range(2))
    check(gram_min>F(1,4), "Reaction-coordinate functional norm below2")
    plus=next(c for c in caps["caps"] if c["name"]=="plus")
    minus=next(c for c in caps["caps"] if c["name"]=="minus")
    candidate_coordinate_separation=F(plus["candidate_a"][1][0])-F(minus["candidate_a"][1][1])
    check(candidate_coordinate_separation>F(89,100), "Candidate saddles separated in a_2")
    check(F(saddle["rho"])<F(13,1000), "Original saddle correction radius")
    original_saddle_separation=(F(89,100)-4*F(13,1000))/2
    check(original_saddle_separation>F(2,5), "Exact saddles separated by more than0.4")
    continued_saddle_separation=F(2,5)-2*radius
    check(continued_saddle_separation>0, "Distinct continued saddle branches")
    check(1-delta>0, "Positive perturbed couplings")
    mapping_lhs = tau*delta + tau*radius*radius/2
    mapping_rhs = sigma*radius
    mapping_margin = mapping_rhs-mapping_lhs
    check(mapping_margin>0, "Frozen-Newton self-map")
    perturbation = tau*radius+delta/2
    kappa = perturbation/sigma
    gap = sigma-perturbation
    check(kappa<1 and gap>0, "Contraction and index persistence")
    min_charge_margin = F(11,10)-tau*radius
    check(min_charge_margin>0, "Dipole charges survive")
    lower, upper = F(global_cert["barrier_lower_exact"]), F(global_cert["barrier_upper_exact"])
    bridge_excess = radius*radius/2
    check(bridge_excess<lower, "Low-energy bridge cannot cross original separator")
    check(sigma-tau*radius>0, "Original energy convex on continuation ball")
    check(lower-2*delta>0, "Perturbed global barrier stays positive")
    check((1+delta/2)*tau<hessian_lipschitz, "Uniform Hessian Lipschitz bound")
    dynamic_gap=gap-hessian_lipschitz*eta
    check(dynamic_gap>rate, "Uniform exponential X-gradient stability")
    exterior_angle=F(18,1000)+tau*radius
    check(exterior_angle<F(41,1000), "Uniform exterior branch bound throughout the coupling class")
    source_absorption=F(41,1000)**2/6+delta
    check(2*source_absorption<1, "Full ell^1 coupling class Green-tail absorption")
    output = {
        "status":"COUPLING_ROBUSTNESS_RATIONAL_CHECKS_PASS",
        "scope":"Analytic corollary from the freshly reproduced exact XY certificate; no perturbed nonlinear solve",
        "coupling_class":"b in ell^1(edges), norm at most1/500; J_e=1+b_e",
        "parameters":{k:str(v) for k,v in {
            "sigma":sigma,"tau":tau,"delta":delta,"critical_radius":radius,
            "dynamical_radius":eta,"exponential_rate":rate}.items()},
        "derived":{k:str(v) for k,v in {
            "self_map_lhs":mapping_lhs,"self_map_rhs":mapping_rhs,
            "self_map_margin":mapping_margin,"contraction_factor":kappa,
            "uniform_critical_gap":gap,"inverse_upper":1/gap,
            "minimum_charge_margin":min_charge_margin,
            "continued_saddle_separation":continued_saddle_separation,
            "bridge_excess":bridge_excess,
            "barrier_lower_exact":lower-2*delta,
            "barrier_upper_exact":upper+2*delta,
            "dynamical_convexity":dynamic_gap,
            "uniform_exterior_angle":exterior_angle,
            "Green_source_absorption":source_absorption}.items()},
        "displayed_strict_barrier_bracket":["0.00167867","0.00970116"],
        "does_not_claim":["perturbed minimax equals a continued saddle height",
            "dominant mechanism change","new perturbed numerical equilibria",
            "local ell^2 XY flow stability","thermal escape rate",
            "proof of the Green kernel estimate by this script"],
        "fresh_evidence":[{"relative_path":name,
            "sha256":hashlib.sha256((args.certificate_root/name).read_bytes()).hexdigest()}
            for name in names],
        "old_certificate_chain_reexecuted_by_this_script":False,
    }
    check(F(output["displayed_strict_barrier_bracket"][0])<lower-2*delta, "Displayed lower enclosure")
    check(F(output["displayed_strict_barrier_bracket"][1])>upper+2*delta, "Displayed upper enclosure")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(output,indent=2)+"\n")
    print(json.dumps({"status":output["status"],"contraction":str(kappa),"gap":str(gap),
        "displayed_bracket":output["displayed_strict_barrier_bracket"]},indent=2))


if __name__=="__main__":
    main()
