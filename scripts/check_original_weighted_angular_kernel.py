"""Focused independent checks for the pure weighted angular kernel."""

from __future__ import annotations

import json
from pathlib import Path
import sys
from typing import Any, Callable, Mapping, Sequence

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from openai_ns_reconstruction.original_weighted_angular_kernel import (  # noqa: E402
    angular_mode_numerator_epsilon_sectors,
)


def _raises(function: Callable[[], Any]) -> bool:
    try:
        function()
    except (TypeError, ValueError, ArithmeticError):
        return True
    return False


def _ordinary_rows(ctx: Any, function: Callable[[Any], Any], z: Any, order: int) -> list[Any]:
    return [ctx.diff(function, z, k) for k in range(order + 1)]


def _evaluate_sectors(sectors: Mapping[int, Sequence[Any]], epsilon: Any) -> list[Any]:
    if not isinstance(sectors, Mapping) or not sectors:
        raise ValueError("sector polynomial must be a nonempty degree mapping")
    first = next(iter(sectors.values()))
    result = [0 for _ in first]
    for degree, row in sectors.items():
        if len(row) != len(first):
            raise ValueError("sector ordinary orders must match")
        weight = epsilon ** degree
        result = [left + weight * right for left, right in zip(result, row)]
    return result


def _full_numerator(
    ctx: Any,
    z: Any,
    epsilon: Any,
    delta_scale: Any,
    G: Callable[[Any], Any],
    phi: Callable[[Any], Any],
    m: Callable[[Any], Any],
    H: Callable[[Any], Any],
    K: Callable[[Any], Any],
) -> Any:
    """Independent dressed original numerator, retaining F0 explicitly."""
    delta = epsilon * delta_scale
    d = 1 - z * z
    f0 = ctx.exp(-G(z) / epsilon)
    h_dressed = ctx.diff(lambda q: ctx.exp(-G(q) / epsilon) * H(q), z) / f0
    k_dressed = ctx.diff(lambda q: ctx.exp(-G(q) / epsilon) * K(q), z) / f0
    W = 1 - (1 - delta) * z * m(z) - d * ctx.diff(m, z)
    angular = (
        (1 - delta / 2) * H(z)
        - (1 - delta) * z * h_dressed / 2
        - d * k_dressed
        + (2 * delta - 1) * z * K(z)
    )
    return -2 * phi(z) * W + angular


def _angular_difference(
    ctx: Any,
    z: Any,
    epsilon: Any,
    delta_scale: Any,
    G: Callable[[Any], Any],
    phi: Callable[[Any], Any],
    H: Callable[[Any], Any],
    K: Callable[[Any], Any],
) -> Any:
    """Independent angular-only numerator for mode 2."""
    delta = epsilon * delta_scale
    d = 1 - z * z
    f0 = ctx.exp(-G(z) / epsilon)
    h_dressed = ctx.diff(lambda q: ctx.exp(-G(q) / epsilon) * H(q), z) / f0
    k_dressed = ctx.diff(lambda q: ctx.exp(-G(q) / epsilon) * K(q), z) / f0
    return (
        (1 - delta / 2) * H(z)
        - (1 - delta) * z * h_dressed / 2
        - d * k_dressed
        + (2 * delta - 1) * z * K(z)
    )


def _mode_1_numerator(
    ctx: Any,
    z: Any,
    epsilon: Any,
    delta_scale: Any,
    phi: Callable[[Any], Any],
    dM: Callable[[Any], Any],
) -> Any:
    delta = epsilon * delta_scale
    d = 1 - z * z
    return 2 * phi(z) * ((1 - delta) * z * dM(z) + d * ctx.diff(dM, z))


def _box(ctx: Any, interval_ctx: Any, value: Any, radius: str = "1e-90") -> Any:
    width = ctx.mpf(radius)
    return interval_ctx.mpf(
        [ctx.nstr(value - width, ctx.dps), ctx.nstr(value + width, ctx.dps)]
    )


def _inside(point: Any, enclosure: Any) -> bool:
    return point in enclosure


def _nonzero(row: Sequence[Any], tolerance: Any) -> bool:
    return any(abs(value) > tolerance for value in row)


def run_fixture() -> dict[str, Any]:
    ctx = mp.mp.clone()
    ctx.dps = 140
    interval_ctx = MPIntervalContext()
    interval_ctx.dps = 180
    controls: list[dict[str, Any]] = []
    differences: list[Any] = []

    def passed(name: str, condition: bool) -> None:
        if not condition:
            raise AssertionError(name)
        controls.append({"name": name, "passed": True})

    j = ctx.mpf(".37")
    # The source normalization fixes delta_scale=exp(970).  A moderate
    # correlated test value keeps the independent F0=exp(-G/epsilon) oracle
    # numerically well-conditioned; first-switch kappa(t)=1-sigma(t) is not
    # this parameter.
    delta_scale = ctx.exp(ctx.mpf("3.2"))

    def u0(q: Any) -> Any:
        return 4 * q + j

    def phi(q: Any) -> Any:
        return 1 + q / 8 + ctx.exp(-q / 5) / 6 + q**2 / 23

    def H(q: Any) -> Any:
        return phi(q) + ctx.sin(q / 4) / 7 + q**3 / 31

    def Psi(q: Any) -> Any:
        return ctx.exp(q / 6) / 5 + q**2 / 19 + ctx.sin(q / 5) / 11

    def Mpsi(q: Any) -> Any:
        return ctx.sin(q / 3) / 4 + q**2 / 17 + q**3 / 43

    def Kpsi(q: Any) -> Any:
        return ctx.exp(-q / 7) / 3 + q / 13 + q**2 / 29

    def G(q: Any) -> Any:
        return q**2 / 17 + ctx.sin(q / 6) / 5 + q**3 / 47

    def input_rows(z: Any) -> dict[str, list[Any]]:
        return {
            "u0": _ordinary_rows(ctx, u0, z, 3),
            "phi": _ordinary_rows(ctx, phi, z, 3),
            "H": _ordinary_rows(ctx, H, z, 3),
            "Psi": _ordinary_rows(ctx, Psi, z, 3),
            "Mpsi": _ordinary_rows(ctx, Mpsi, z, 3),
            "Kpsi": _ordinary_rows(ctx, Kpsi, z, 3),
            "Gprime": [ctx.diff(G, z, k + 1) for k in range(3)],
        }

    for z_text in ("-.31", ".07", ".29"):
        z = ctx.mpf(z_text)
        rows = input_rows(z)
        decomposition = angular_mode_numerator_epsilon_sectors(
            z,
            rows["u0"],
            rows["phi"],
            rows["H"],
            rows["Psi"],
            rows["Mpsi"],
            rows["Kpsi"],
            rows["Gprime"],
            delta_scale,
        )
        passed(
            f"result_shape_{z_text}",
            set(decomposition) == {"numerator_modes", "denominator_sectors"}
            and len(decomposition["numerator_modes"]) == 3,
        )
        modes = decomposition["numerator_modes"]
        passed(
            f"all_laurent_support_retained_{z_text}",
            set(modes[0]) == {-1, 0, 1, 2}
            and set(modes[1]) == {1, 2}
            and set(modes[2]) == {-1, 0, 1, 2},
        )
        passed(
            f"negative_degree_retained_{z_text}",
            -1 in modes[0]
            and -1 in modes[2]
            and _nonzero(modes[0][-1], ctx.mpf("1e-80"))
            and _nonzero(modes[2][-1], ctx.mpf("1e-80")),
        )
        passed(
            f"common_denominator_support_{z_text}",
            set(decomposition["denominator_sectors"]) == {0, 1},
        )

        for epsilon_text in ("-.23", ".17", ".41"):
            epsilon = ctx.mpf(epsilon_text)
            delta = epsilon * delta_scale

            m0 = lambda q: u0(q) + epsilon * Psi(q)
            h0 = phi
            k0 = lambda q: phi(q) * m0(q)
            dM = lambda q: epsilon * (Mpsi(q) - Psi(q))
            dH = lambda q: H(q) - phi(q)
            dK = lambda q: u0(q) * dH(q) + epsilon * (Kpsi(q) - phi(q) * Psi(q))

            direct_modes = (
                lambda q: _full_numerator(ctx, q, epsilon, delta_scale, G, phi, m0, h0, k0),
                lambda q: _mode_1_numerator(ctx, q, epsilon, delta_scale, phi, dM),
                lambda q: _angular_difference(ctx, q, epsilon, delta_scale, G, phi, dH, dK),
            )
            for mode, direct_function in enumerate(direct_modes):
                sector_row = _evaluate_sectors(modes[mode], epsilon)
                oracle = [ctx.diff(direct_function, z, order) for order in range(3)]
                for order, (value, expected) in enumerate(zip(sector_row, oracle)):
                    difference = abs(value - expected)
                    differences.append(difference)
                    passed(
                        f"mode_{mode}_matches_dressed_mp_diff_{z_text}_{epsilon_text}_Z{order}",
                        difference < ctx.mpf("1e-78"),
                    )

            denominator_row = _evaluate_sectors(
                decomposition["denominator_sectors"], epsilon
            )
            denominator_oracle = [
                ctx.diff(lambda q: 2 * phi(q) * (1 - delta * q * q), z, order)
                for order in range(3)
            ]
            passed(
                f"denominator_matches_direct_{z_text}_{epsilon_text}",
                all(
                    abs(value - expected) < ctx.mpf("1e-78")
                    for value, expected in zip(denominator_row, denominator_oracle)
                ),
            )

            for theta_text in ("-.4", ".23", ".71"):
                theta = ctx.mpf(theta_text)
                m_full = lambda q: m0(q) + theta * dM(q)
                h_full = lambda q: h0(q) + theta**2 * dH(q)
                k_full = lambda q: k0(q) + theta**2 * dK(q)
                reconstructed = [
                    left
                    + theta * middle
                    + theta**2 * right
                    for left, middle, right in zip(
                        _evaluate_sectors(modes[0], epsilon),
                        _evaluate_sectors(modes[1], epsilon),
                        _evaluate_sectors(modes[2], epsilon),
                    )
                ]
                direct_full = [
                    ctx.diff(
                        lambda q: _full_numerator(
                            ctx,
                            q,
                            epsilon,
                            delta_scale,
                            G,
                            phi,
                            m_full,
                            h_full,
                            k_full,
                        ),
                        z,
                        order,
                    )
                    for order in range(3)
                ]
                passed(
                    f"theta_reconstruction_{z_text}_{epsilon_text}_{theta_text}",
                    all(
                        abs(value - expected) < ctx.mpf("1e-76")
                        for value, expected in zip(reconstructed, direct_full)
                    ),
                )

        # Directed containment is checked coefficient by coefficient, including
        # the Laurent gauge sector, with independent narrow intervals.
        iv_rows = {
            name: [_box(ctx, interval_ctx, value) for value in row]
            for name, row in rows.items()
        }
        iv_result = angular_mode_numerator_epsilon_sectors(
            _box(ctx, interval_ctx, z),
            iv_rows["u0"],
            iv_rows["phi"],
            iv_rows["H"],
            iv_rows["Psi"],
            iv_rows["Mpsi"],
            iv_rows["Kpsi"],
            iv_rows["Gprime"],
            _box(ctx, interval_ctx, delta_scale),
        )
        point_modes = decomposition["numerator_modes"]
        passed(
            f"directed_containment_{z_text}",
            all(
                _inside(value, iv_result["numerator_modes"][mode][degree][order])
                for mode in range(3)
                for degree, row in point_modes[mode].items()
                for order, value in enumerate(row)
            )
            and all(
                _inside(value, iv_result["denominator_sectors"][degree][order])
                for degree, row in decomposition["denominator_sectors"].items()
                for order, value in enumerate(row)
            ),
        )

    z = ctx.mpf(".13")
    rows = input_rows(z)
    names = ("u0", "phi", "H", "Psi", "Mpsi", "Kpsi")
    for name in names:
        bad = dict(rows)
        bad[name] = rows[name][:3]
        passed(
            f"missing_{name}_Z3_rejected",
            _raises(
                lambda bad=bad: angular_mode_numerator_epsilon_sectors(
                    z,
                    bad["u0"],
                    bad["phi"],
                    bad["H"],
                    bad["Psi"],
                    bad["Mpsi"],
                    bad["Kpsi"],
                    bad["Gprime"],
                    delta_scale,
                )
            ),
        )
    passed(
        "missing_Gprime_Z2_rejected",
        _raises(
            lambda: angular_mode_numerator_epsilon_sectors(
                z,
                rows["u0"],
                rows["phi"],
                rows["H"],
                rows["Psi"],
                rows["Mpsi"],
                rows["Kpsi"],
                rows["Gprime"][:2],
                delta_scale,
            )
        ),
    )

    return {
        "passed": True,
        "control_count": len(controls),
        "controls": controls,
        "maximum_supporting_difference": ctx.nstr(max(differences), 25),
        "ordinary_formula": (
            "N=-2*phi*W+(1-delta/2)*H-(1-delta)*z*(H_Z+g*H)/2 "
            "-d*(K_Z+g*K)+(2*delta-1)*z*K"
        ),
        "api": {
            "function": "angular_mode_numerator_epsilon_sectors",
            "inputs": "z,u0,phi,H,Psi,Mpsi,Kpsi,Gprime,delta_scale",
            "rows": "u0,phi,H,Psi,Mpsi,Kpsi are raw undressed ordinary Z0..Z3; Gprime is Z0..Z2",
            "returns": "numerator_modes[0..2] and denominator_sectors, each degree -> ordinary Z0..Z2",
            "modes": [
                "m=u0+epsilon*Psi, H=phi, K=phi*(u0+epsilon*Psi)",
                "mDiff=epsilon*(Mpsi-Psi), H=K=0",
                "Hdiff=H-phi, Kdiff=u0*(H-phi)+epsilon*(Kpsi-phi*Psi)",
            ],
            "common_denominator": "normalized 2*phi*L = 2*phi*(1-(epsilon*delta_scale)*z**2)",
            "laurent": "g=-Gprime/epsilon is degree -1; delta=epsilon*delta_scale is degree 1",
            "source_scale": "delta_scale=exp(970) is fixed; first-switch kappa(t)=1-sigma(t) is separate",
        },
        "ancestor_constructor_producer_checker_fixture_run": False,
        "source_loading_or_adapter": False,
        "physical_scales_or_pins": False,
        "quadrature_or_integration": False,
        "temporal_recursion": False,
        "NS_residual": False,
    }


if __name__ == "__main__":
    print(json.dumps(run_fixture(), indent=2))
