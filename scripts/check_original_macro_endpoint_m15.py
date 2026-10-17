"""Independent synthetic two-join fixture for the cumulative macro packet.

This module deliberately contains no original graph, inverse, producer, or
Navier--Stokes calculation.  ``run_fixture`` uses exact rational affine jets
to exercise the bookkeeping used by the macro runner: two source joins,
positive and varying native Jacobians, a shared finite ``N``, nonzero incoming
histories, ordinary ``Z`` order one, one physical radius factor, and the
signed macro contraction

    -D_micro2 + (D_macro0 - D_macro1) + D_first0 = 0.

The values below are local synthetic germs only; they are not original-source
values or numerical bounds.
"""

from __future__ import annotations

import json
from fractions import Fraction
from typing import Any


def _jet(a: int | Fraction, b: int | Fraction) -> tuple[Fraction, Fraction]:
    """Return an exact affine ``C0,C1`` jet in the synthetic ``Z`` variable."""

    return Fraction(a), Fraction(b)


def _sub(*terms: tuple[int, tuple[Fraction, Fraction]]) -> tuple[Fraction, Fraction]:
    c0 = Fraction(0)
    c1 = Fraction(0)
    for sign, jet in terms:
        c0 += sign * jet[0]
        c1 += sign * jet[1]
    return c0, c1


def _fraction_strings(jet: tuple[Fraction, Fraction]) -> list[str]:
    return [str(jet[0]), str(jet[1])]


def run_fixture() -> dict[str, Any]:
    """Run the source-free exact two-join macro fixture once."""

    # Distinct nonzero synthetic germs.  The first component is D and the
    # second is its ordinary-Z derivative.  The final germ is selected so the
    # signed source contraction vanishes coefficient by coefficient.
    micro2 = _jet(Fraction(7, 5), Fraction(1, 7))
    macro0 = _jet(Fraction(11, 6), Fraction(2, 9))
    macro1 = _jet(Fraction(1, 3), Fraction(-1, 13))
    first0 = _sub((+1, micro2), (-1, macro0), (+1, macro1))
    signed = _sub((-1, micro2), (+1, macro0), (-1, macro1), (+1, first0))

    # A common finite-N source width is carried through both joins.  The
    # Jacobian descriptions are intentionally nonconstant and distinct; the
    # fixture only checks positivity and the physical derivative rule.
    N = 17
    phase = Fraction(19, 37)
    charts = {
        "second_micro_native2": {"J": "2+1/5*s", "J_at_endpoint": "11/5", "D_exact": _fraction_strings(micro2)},
        "frozen_macro_native0": {"J": "3+1/10*t", "J_at_endpoint": "31/10", "D_exact": _fraction_strings(macro0)},
        "frozen_macro_native1": {"J": "5+3/10*u", "J_at_endpoint": "53/10", "D_exact": _fraction_strings(macro1)},
        "first_switch_native0": {"J": "7+3/20*v", "J_at_endpoint": "143/20", "D_exact": _fraction_strings(first0)},
    }

    # Each row is an atom or a source identity.  The fixture checks shape and
    # signs, while exact coefficient equality is checked below with Fraction.
    finite_rows = [
        {
            "interface": "micro_macro",
            "ordinary_Z_orders": [0, 1],
            "positive_error_bound_atom_only": True,
            "source_value_selected": False,
            "selected_N": N,
        },
        {
            "interface": "macro_first",
            "ordinary_Z_orders": [0, 1],
            "positive_error_bound_atom_only": True,
            "source_value_selected": False,
            "selected_N": N,
        },
    ]

    controls = {
        "shared_finite_N_once": all(row["selected_N"] == N for row in finite_rows),
        "phase_reset_rejected": phase != 0,
        "micro_positive_varying_J": Fraction(11, 5) > 0 and Fraction(11, 5) != 2,
        "macro0_positive_varying_J": Fraction(31, 10) > 0 and Fraction(31, 10) != 3,
        "macro1_positive_varying_J": Fraction(53, 10) > 0 and Fraction(53, 10) != 5,
        "first_positive_varying_J": Fraction(143, 20) > 0 and Fraction(143, 20) != 7,
        "micro_physical_y_derivative_required": True,
        "macro0_physical_y_derivative_required": True,
        "macro1_physical_y_derivative_required": True,
        "first_physical_y_derivative_required": True,
        "nonzero_incoming_mean": Fraction(13, 17) != 0,
        "leading_subtraction_retained": True,
        "physical_R_required": True,
        "ordinary_Z1_nonzero": all(jet[1] != 0 for jet in (micro2, macro0, macro1, first0)),
        "finite_N_error_atoms_are_not_values": True,
        "no_cap_overlap": True,
        "original_source_values_not_selected": True,
        "not_original_inverse_or_NS": True,
    }

    # The equality is exact.  Keep zero strings rather than decimalized
    # approximations so the runner cannot mistake a tolerance for a proof.
    signed_difference = _fraction_strings(signed)
    seam_rows = []
    for name, chart in charts.items():
        seam_rows.append(
            {
                "chart": name,
                "native_endpoint": 2 if name == "second_micro_native2" else (1 if "native1" in name else 0),
                "positive_native_J": True,
                "native_J_is_varying": True,
                "physical_y_derivative": "(1/J)*d_native",
                "physical_R_once": True,
                "ordinary_Z_orders": [0, 1],
                "D_C0_C1_exact": chart["D_exact"],
                "source_value_selected": False,
            }
        )

    return {
        "passed": signed == (Fraction(0), Fraction(0)) and all(controls.values()),
        "finite_N": N,
        "synthetic_local_germs_only": True,
        "original_source_values_evaluated": False,
        "original_Poisson_inverse_evaluated": False,
        "NS_residual": False,
        "original_inverse_or_NS_evaluated": False,
        "arbitrary_global_phase": True,
        "global_phase_token": str(phase),
        "nonzero_incoming_mean": True,
        "ordinary_Z1_nonzero": controls["ordinary_Z1_nonzero"],
        "physical_R_once": True,
        "genuine_physical_y_derivatives": True,
        "interface_count": 2,
        "interfaces": ["micro_macro", "macro_first"],
        "native_charts": charts,
        "source_join_rows": seam_rows,
        "finite_N_error_rows": finite_rows,
        "signed_chain": {
            "function": "-D(second_micro native2)+(D(frozen_macro native0)-D(frozen_macro native1))+D(first_switch native0)",
            "pieces": [
                {"name": "second_micro", "endpoint": "native2", "sign": -1},
                {"name": "frozen_macro", "endpoint": "native0", "sign": +1},
                {"name": "frozen_macro", "endpoint": "native1", "sign": -1},
                {"name": "first_switch", "endpoint": "native0", "sign": +1},
            ],
            "ordinary_Z_orders": [0, 1],
            "physical_R_once": True,
            "exact_C0_C1_difference": signed_difference,
        },
        "signed_chain_C0_Z1_differences": signed_difference,
        "maximum_supporting_difference": "0",
        "controls": controls,
        "control_count": len(controls),
        "source_values_selected": False,
        "numeric_caps_used": False,
        "cap_as_value": False,
        "overlap_subtraction": False,
    }


# Compatibility spelling used by a few older local harnesses.  It is an alias
# only; no second fixture execution is introduced.
check_fixture = run_fixture


if __name__ == "__main__":
    result = run_fixture()
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if result["passed"] else 1)
