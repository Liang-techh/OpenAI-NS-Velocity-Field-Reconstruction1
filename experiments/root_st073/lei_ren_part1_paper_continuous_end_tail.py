"""Focused diagnostic for the direct continuous end-bump tail.

This reads the latest installed incoming receipt and continuous coefficient
receipt.  It does not build the global field.  The diagnostic checks the
reflected tail primitive near the final bump endpoint and verifies that the
anchored cumulative row carries the materialized terminal balance ``rho``.
The balance is reported as an evaluated residual; it is never replaced by
zero.  Quadrature and incoming-integral uncertainty remain outside the
identity checks.
"""

import json
from pathlib import Path

import mpmath as mp

from lei_ren_part1_paper_axial_correction import from_signed_log, signed_log
from lei_ren_part1_paper_continuous_axial_basis import ContinuousAxialBump
from lei_ren_part1_paper_continuous_axial_solve import ContinuousAxialCorrection


HERE = Path(__file__).parent
INCOMING_RECEIPT = HERE / "lei_ren_part1_paper_continuous_incoming_outer.json"
TERMINAL_RECEIPT = HERE / "lei_ren_part1_paper_continuous_terminal_balance.json"


def _signed(value, digits=100):
    """Serialize an MP value without converting its exponent to a float."""

    value = mp.mpf(value)
    if value == 0:
        return {"sign": 0, "log_abs": None, "arbitrary_exponent_value": "0"}
    return signed_log(value, digits)


def _relative(actual, expected):
    actual = mp.mpf(actual)
    expected = mp.mpf(expected)
    if expected == 0:
        return abs(actual)
    return abs(actual - expected) / abs(expected)


def _finite_difference_tail(basis, lam, s, ell):
    """Differentiate the direct tail while staying inside its support."""

    gap = ell - s
    h = min(mp.mpf("1e-5"), max(gap / 4, mp.mpf("1e-30")))
    if s + h < ell:
        return (basis.tail(lam, s + h) - basis.tail(lam, s - h)) / (2 * h)
    return (basis.tail(lam, s) - basis.tail(lam, s - h)) / h


def run():
    receipt = json.loads(INCOMING_RECEIPT.read_text(encoding="utf-8"))
    incoming = receipt["incoming"]
    solve = receipt["continuous_solve"]
    base_rows = [
        incoming["row_normalization"][key]
        for key in ("scaled_base_m1", "scaled_base_m2")
    ]

    # The saved solve is the latest installed component.  No joined-field
    # construction or incoming regeneration is needed for this diagnostic.
    component = ContinuousAxialCorrection(solve, base_rows)
    precision = component.precision

    with mp.workdps(precision):
        basis = ContinuousAxialBump(precision=precision)
        mu = mp.mpf(solve["input_mu"])
        ell = basis.ell

        # The final bump is centred at -1, so its right support endpoint is
        # offset -0.85.  These gaps are close to ell while keeping the
        # finite-difference derivative check in a resolved endpoint layer.
        gaps = [mp.mpf("5e-2"), mp.mpf("2e-2"), mp.mpf("1e-2"), mp.mpf("5e-3")]
        tail_rows = []
        lam = mp.mpf(".5") - mu
        for gap in gaps:
            s = ell - gap
            partial = basis.primitive(lam, s)
            tail = basis.tail(lam, s)
            full = basis.full(mp.nstr(lam, precision))
            derivative = -mp.exp(lam * s) * basis.values(s)["beta"]
            reflected_derivative = -basis.primitive_jet(-lam, -s)["derivative"]
            finite_difference = _finite_difference_tail(basis, lam, s, ell)
            tail_rows.append(
                {
                    "gap": mp.nstr(gap, 30),
                    "s": mp.nstr(s, 40),
                    "tail": _signed(tail),
                    "tail_positive": bool(tail > 0),
                    "full": _signed(full),
                    "partial_plus_tail_relative_error": mp.nstr(
                        _relative(partial + tail, full), 60
                    ),
                    "direct_derivative": _signed(derivative),
                    "reflected_primitive_derivative": _signed(reflected_derivative),
                    "direct_vs_reflected_derivative_relative_error": mp.nstr(
                        _relative(reflected_derivative, derivative), 60
                    ),
                    "finite_difference_derivative": _signed(finite_difference),
                    "derivative_relative_error": mp.nstr(
                        _relative(finite_difference, derivative), 60
                    ),
                }
            )

        endpoint = ell
        endpoint_tail = basis.tail(lam, endpoint)
        endpoint_partial = basis.primitive(lam, endpoint)
        endpoint_full = basis.full(mp.nstr(lam, precision))

        # The last bump ends at offset -0.85.  Include points immediately
        # below and above that endpoint, where the direct tail must remain
        # positive below and exactly zero at/above the endpoint.
        support_end = -mp.mpf(".85")
        # Decimal .85 and the MP representation of ell can differ by one
        # last bit at 200 dps. Nudge the endpoint outward by an invisible
        # guard so weighted_tail sees the exact post-support branch.
        endpoint_guard = mp.mpf("1e-180")
        offsets = [
            support_end - mp.mpf("1e-6"),
            support_end + endpoint_guard,
            support_end + mp.mpf("1e-6"),
        ]
        rho = component.terminal_balance(1, precision)["value"]
        cumulative_rows = []
        for offset in offsets:
            tail = component.weighted_tail(1, offset)
            cumulative = component.cumulative_row(1, end_offset=offset)
            cumulative_value = from_signed_log(cumulative["nominal"])
            expected = rho - tail
            cumulative_rows.append(
                {
                    "offset": mp.nstr(offset, 40),
                    "weighted_tail": _signed(tail),
                    "cumulative": _signed(cumulative_value),
                    "rho_minus_tail": _signed(expected),
                    "cumulative_identity_relative_error": mp.nstr(
                        _relative(cumulative_value, expected), 60
                    ),
                    "method": cumulative["method"],
                    "quadrature_enclosure_certified": cumulative[
                        "quadrature_enclosure_certified"
                    ],
                }
            )

        endpoint_rows = []
        for row in (1, 2):
            row_rho = component.terminal_balance(row, precision)["value"]
            row_cumulative = from_signed_log(
                component.cumulative_row(row, end_offset=-mp.mpf(".85"))["nominal"]
            )
            endpoint_rows.append(
                {
                    "row": row,
                    "rho": _signed(row_rho),
                    "cumulative_at_last_bump_end": _signed(row_cumulative),
                    "anchor_relative_error": mp.nstr(
                        _relative(row_cumulative, row_rho), 60
                    ),
                    "rho_is_retained": bool(row_rho != 0),
                }
            )

        terminal_receipt = {}
        if TERMINAL_RECEIPT.exists():
            terminal_receipt = json.loads(
                TERMINAL_RECEIPT.read_text(encoding="utf-8")
            )

        report = {
            "diagnostic": "direct continuous end-bump tail and anchored cumulative row",
            "input_receipt": INCOMING_RECEIPT.name,
            "continuous_solve_source": "continuous_solve",
            "precision": precision,
            "input_mu": mp.nstr(mu, precision),
            "bump_ell": mp.nstr(ell, 40),
            "last_bump_center": "-1",
            "last_bump_support_end_offset": "-.85",
            "support_end_endpoint_guard": "1e-180 outward MP guard",
            "tail_definition": "tail(lam,s)=primitive(-lam,-s)",
            "tail_derivative_definition": "-exp(lam*s)*beta(s)",
            "tail_checks": tail_rows,
            "endpoint_identity": {
                "partial_plus_tail_relative_error": mp.nstr(
                    _relative(endpoint_partial + endpoint_tail, endpoint_full), 60
                ),
                "endpoint_tail": _signed(endpoint_tail),
                "endpoint_tail_is_zero": bool(endpoint_tail == 0),
            },
            "anchored_cumulative_rows": cumulative_rows,
            "support_end_anchors": endpoint_rows,
            "terminal_row_one_rho": _signed(rho),
            "terminal_mean_forced_zero": False,
            "materialized_terminal_balance_retained": True,
            "exact_functional_identity": (
                "For complete continuous atoms, M c + b + a p = 0; "
                "the evaluated rho is retained as a materialized residual."
            ),
            "mean_Z_residual_preserved": bool(
                terminal_receipt.get("installed_normalized_mass_Z") is not None
            ),
            "installed_normalized_mass_Z": terminal_receipt.get(
                "installed_normalized_mass_Z"
            ),
            "quadrature_enclosure_certified": False,
            "incoming_uncertainty_enclosed": False,
            "finite_energy_certified": False,
            "scope": (
                "Direct tail and partial-primitive identities only; no global "
                "energy, recursion, or quadrature enclosure claim."
            ),
        }

    with mp.workdps(100):
        small_basis=ContinuousAxialBump(precision=100)
        s=mp.mpf('.1499');lam=mp.mpf('.5')
        direct=small_basis.tail(lam,s)
        subtractive=small_basis.full('.5')-small_basis.primitive(lam,s)
        if direct<=0 or subtractive!=0:
            raise AssertionError('Flat endpoint cancellation demonstration changed')
        report['flat_endpoint_cancellation_example']={
            'precision':100,'s':'.1499','lambda':'.5',
            'direct_tail':_signed(direct,100),
            'full_minus_partial':_signed(subtractive,100),
            'positive_remaining_mass_preserved':True}

    output = HERE / "lei_ren_part1_paper_continuous_end_tail.json"
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "endpoint_tail_is_zero": report["endpoint_identity"][
                    "endpoint_tail_is_zero"
                ],
                "terminal_mean_forced_zero": report["terminal_mean_forced_zero"],
                "mean_Z_residual_preserved": report[
                    "mean_Z_residual_preserved"
                ],
                "support_end_row1_anchor_error": report["support_end_anchors"][
                    0
                ]["anchor_relative_error"],
            }
        ),
        flush=True,
    )
    return report


if __name__ == "__main__":
    run()
