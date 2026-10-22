#!/usr/bin/env python3
"""Exact arithmetic for the new Green-tail localization budget.

This verifies arithmetic consequences of the analytic kernel theorem and
previously certified path inputs. It does not independently validate the
old global certificate, a critical point, or the analytic Green proof.
"""
from fractions import Fraction as F
from pathlib import Path
import json


def fraction(x):
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def check(condition, message):
    """Keep certificate obligations active even under python -O."""
    if not condition:
        raise ValueError(message)


def calculate(epsilon):
    m_critical = F(6007)
    m_path = F(6007)
    r_star = 4097
    delta = F(9, 500)
    a = delta * delta / 6
    q_upper = m_critical**3 / 6
    # 32 sqrt(2) < 46 and sqrt(227) < 16, checked without floating point.
    check(32**2 * 2 < 46**2, "invalid rational kernel majorant")
    check(227 < 16**2, "invalid rational cutoff majorant")
    check(0 <= 2 * a < 1, "tail recurrence is not contractive")
    check(epsilon > 0, "localization tolerance must be positive")
    a_upper = max(r_star * m_critical, 46 * q_upper / (1 - 2 * a))
    projection_constant = max(16 * a_upper, F(480))
    energy_constant = m_path * projection_constant
    n, radius = 0, r_star
    while energy_constant / radius > epsilon / 4:
        n += 1
        radius *= 2
    chi = energy_constant / radius
    check(chi <= epsilon / 4, "localization allocation not met")
    check(n == 0 or energy_constant / (radius // 2) > epsilon / 4,
          "dyadic exponent is not minimal")
    check(radius == r_star * 2**n, "inconsistent dyadic radius")
    return {
        "epsilon": fraction(epsilon),
        "a": fraction(a),
        "q_upper": fraction(q_upper),
        "a_tail_constant_upper": fraction(a_upper),
        "projection_constant_upper": fraction(projection_constant),
        "energy_constant_upper": fraction(energy_constant),
        "dyadic_exponent": n,
        "radius": str(radius),
        "radius_decimal_digits": len(str(radius)),
        "localization_error_upper": fraction(chi),
        "allocated_localization_error": fraction(epsilon / 4),
        "finite_box_materialized": False,
        "new_barrier_interval_claimed": False,
    }


def robust_budget():
    """Rational bounds for the same full ell1 coupling class; no new solve."""
    m, beta, eta, r_star = F(6008), F(1,500), F(41,1000), 4098
    a = eta**2 / 6 + beta
    q = m**3 / 6 + beta
    check(a == F(13681,6000000) and 2*a < 1, "perturbed absorption mismatch")
    check(F(18,1000) + F(71,100)/32 < eta, "continued angle exceeds bound")
    a_upper = max(r_star*m, 46*q/(1-2*a))
    check(a_upper < F(17,10)*10**12, "displayed tail majorant too small")
    check(16*a_upper < F(272,100)*10**13, "displayed projection majorant too small")
    return {
        "coupling_budget": fraction(beta), "critical_norm_upper": fraction(m),
        "exterior_angle": fraction(eta), "common_anchor_radius": r_star,
        "source_absorption": fraction(a), "source_l1_upper": fraction(q),
        "tail_constant_upper": fraction(a_upper),
        "displayed_tail_constant_upper": "1700000000000",
        "displayed_projection_constant_upper": "27200000000000",
        "no_coupling_support_or_moment_requirement": True,
        "new_perturbed_critical_solve": False,
        "conditional_on_analytic_kernel_and_critical_continuation": True
    }


if __name__ == "__main__":
    data = {
        "status": "EXACT_GREEN_TAIL_BUDGET_ARITHMETIC_PASS",
        "analytic_kernel_theorem_proved_by_this_program": False,
        "old_global_certificate_reexecuted_by_this_program": False,
        "scope": "Exact arithmetic only, conditional on the separately proved analytic theorem and certified path inputs.",
        "cases": [calculate(F(1, 100)), calculate(F(1, 10**6))],
        "uniform_perturbed_class": robust_budget(),
    }
    output = Path(__file__).with_name("green_tail_budget.json")
    output.write_text(json.dumps(data, indent=2) + "\n")
    print(json.dumps(data, indent=2))
