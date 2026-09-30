"""Bounded finite-core kinetic-energy scaling receipt.

The computation integrates the actual regular core supplied by
``build_joined_field().inner.core`` over ``0 <= R <= inner.r`` and
``-.5 <= Z <= .5``.  It uses the physical cylindrical Jacobian from the
similarity chart, separates radial from swirl-plus-axial contributions, and
keeps the requested extreme tau values in signed-log form.  This is a finite
core quadrature sample; it is not a global finite-energy or NS-closure claim.
"""

from __future__ import annotations

import json
from pathlib import Path

import mpmath as mp

from lei_ren_part1_paper_exit_field import physical_chart
from lei_ren_part1_paper_joined_outer import build_joined_field


OUT_JSON = Path(__file__).with_suffix(".json")


def _text(value: mp.mpf, digits: int) -> str:
    return mp.nstr(mp.mpf(value), digits)


def signed_log(value: mp.mpf, digits: int) -> dict[str, object]:
    """Serialize a possibly extreme MP number without binary floats."""

    value = mp.mpf(value)
    if value == 0:
        return {
            "sign": 0,
            "log_abs": None,
            "arbitrary_exponent_value": "0",
        }
    return {
        "sign": 1 if value > 0 else -1,
        "log_abs": _text(mp.log(abs(value)), digits),
        "arbitrary_exponent_value": _text(value, digits),
    }


def _relative(left: mp.mpf, right: mp.mpf) -> mp.mpf:
    left = mp.mpf(left)
    right = mp.mpf(right)
    scale = max(abs(left), abs(right))
    return abs(left - right) / scale if scale else abs(left - right)


def _legendre_nodes(order: int) -> tuple[list[mp.mpf], list[mp.mpf]]:
    nodes, weights = mp.gauss_quadrature(order, "legendre")
    return (
        [mp.mpf(nodes[i]) for i in range(order)],
        [mp.mpf(weights[i]) for i in range(order)],
    )


def _samples(core, R_max: mp.mpf, radial_order: int, z_order: int):
    """Return actual core values and mapped quadrature weights."""

    radial_nodes, radial_weights = _legendre_nodes(radial_order)
    z_nodes, z_weights = _legendre_nodes(z_order)
    samples = []
    for x, wx in zip(radial_nodes, radial_weights):
        R = R_max * (x + 1) / 2
        weight_R = R_max * wx / 2
        for xz, wz in zip(z_nodes, z_weights):
            Z = xz / 2
            weight_Z = wz / 2
            state = core.evaluate(R, Z)
            samples.append(
                {
                    "weight": weight_R * weight_Z,
                    "R": R,
                    "Z": Z,
                    "Ur": mp.mpf(state["Ur"]),
                    "F": mp.mpf(state["F"]),
                    "Uz": mp.mpf(state["Uz"]),
                }
            )
    return samples


def _shape_integrals(samples, delta: mp.mpf) -> tuple[mp.mpf, mp.mpf]:
    """Integrate the tau-independent radial and swirl/axial shapes."""

    radial_terms = []
    swirl_axial_terms = []
    for sample in samples:
        Z = sample["Z"]
        R = sample["R"]
        d = 1 - Z * Z
        common = 1 - delta * Z * Z
        radial_terms.append(
            sample["weight"]
            * sample["Ur"]
            * sample["Ur"]
            * common
            * d ** (-(3 - delta) / 2)
        )
        swirl_axial_terms.append(
            sample["weight"]
            * (2 * R * sample["F"] * sample["F"] + sample["Uz"] * sample["Uz"])
            * common
            * d ** (-(3 - 3 * delta) / 2)
        )
    return mp.fsum(radial_terms), mp.fsum(swirl_axial_terms)


def _direct_integral(
    samples,
    *,
    tau: mp.mpf,
    delta: mp.mpf,
    nu: mp.mpf,
) -> tuple[mp.mpf, mp.mpf]:
    """Integrate the physical kinetic energy directly at fixed tau."""

    radial_terms = []
    swirl_axial_terms = []
    for sample in samples:
        Z = sample["Z"]
        R = sample["R"]
        d = 1 - Z * Z
        q = tau / d
        jacobian = (
            2
            * mp.pi
            * nu ** (mp.mpf(3) / 2)
            * q ** ((3 - delta) / 2)
            * (1 - delta * Z * Z)
            / d
        )
        radial_speed_sq = nu * q ** (-1) * sample["Ur"] ** 2
        swirl_axial_speed_sq = nu * q ** (-(1 + delta)) * (
            2 * R * sample["F"] ** 2 + sample["Uz"] ** 2
        )
        radial_terms.append(
            sample["weight"] * mp.mpf(".5") * radial_speed_sq * jacobian
        )
        swirl_axial_terms.append(
            sample["weight"]
            * mp.mpf(".5")
            * swirl_axial_speed_sq
            * jacobian
        )
    return mp.fsum(radial_terms), mp.fsum(swirl_axial_terms)


def _separated_energy(
    shape_radial: mp.mpf,
    shape_swirl_axial: mp.mpf,
    *,
    tau: mp.mpf,
    delta: mp.mpf,
    nu: mp.mpf,
) -> tuple[mp.mpf, mp.mpf]:
    constant = mp.pi * nu ** (mp.mpf(5) / 2)
    radial = constant * tau ** ((1 - delta) / 2) * shape_radial
    swirl_axial = constant * tau ** ((1 - 3 * delta) / 2) * shape_swirl_axial
    return radial, swirl_axial


def _energy_row(
    direct_radial: mp.mpf,
    direct_swirl_axial: mp.mpf,
    separated_radial: mp.mpf,
    separated_swirl_axial: mp.mpf,
    *,
    logtau: mp.mpf,
    delta: mp.mpf,
    baseline_radial: mp.mpf,
    baseline_swirl_axial: mp.mpf,
    baseline_total: mp.mpf,
    digits: int,
) -> dict[str, object]:
    direct_total = direct_radial + direct_swirl_axial
    separated_total = separated_radial + separated_swirl_axial
    radial_gain = direct_radial / baseline_radial
    swirl_gain = direct_swirl_axial / baseline_swirl_axial
    total_gain = direct_total / baseline_total
    expected_radial_log = ((1 - delta) / 2) * logtau
    expected_swirl_log = ((1 - 3 * delta) / 2) * logtau
    row = {
        "logtau": _text(logtau, digits),
        "tau": signed_log(mp.exp(logtau), digits),
        "direct_radial_energy": signed_log(direct_radial, digits),
        "direct_swirl_plus_axial_energy": signed_log(direct_swirl_axial, digits),
        "direct_total_energy": signed_log(direct_total, digits),
        "separated_radial_energy": signed_log(separated_radial, digits),
        "separated_swirl_plus_axial_energy": signed_log(separated_swirl_axial, digits),
        "separated_total_energy": signed_log(separated_total, digits),
        "direct_vs_separated_radial_relative_error": signed_log(
            _relative(direct_radial, separated_radial), digits
        ),
        "direct_vs_separated_swirl_plus_axial_relative_error": signed_log(
            _relative(direct_swirl_axial, separated_swirl_axial), digits
        ),
        "direct_vs_separated_total_relative_error": signed_log(
            _relative(direct_total, separated_total), digits
        ),
        "radial_energy_ratio_to_logtau0": signed_log(radial_gain, digits),
        "swirl_plus_axial_energy_ratio_to_logtau0": signed_log(swirl_gain, digits),
        "total_energy_ratio_to_logtau0": signed_log(total_gain, digits),
        "expected_radial_log_gain": _text(expected_radial_log, digits),
        "expected_swirl_plus_axial_log_gain": _text(expected_swirl_log, digits),
        "measured_minus_expected_radial_log_gain": _text(
            mp.log(abs(radial_gain)) - expected_radial_log, digits
        ),
        "measured_minus_expected_swirl_plus_axial_log_gain": _text(
            mp.log(abs(swirl_gain)) - expected_swirl_log, digits
        ),
    }
    return row


def _physical_integrand_check(inner, core, *, R, Z, tau, nu, delta, precision, digits):
    """Compare one actual physical callable integrand with the separated form."""

    d = 1 - Z * Z
    q = tau / d
    state = core.evaluate(R, Z)
    physical = inner.physical_field(nu=nu, T=0)
    chart = physical_chart(
        {**state, "R": R},
        Z,
        mp.log(q),
        delta=delta,
        nu=nu,
        phi=mp.mpf(".23"),
        precision=precision,
    )
    back = physical.evaluate(*chart["xyz"], -chart["tau"])
    ur, utheta, uz = back["cylindrical_velocity"]
    jacobian = (
        2
        * mp.pi
        * nu ** (mp.mpf(3) / 2)
        * q ** ((3 - delta) / 2)
        * (1 - delta * Z * Z)
        / d
    )
    physical_density = mp.mpf(".5") * (ur * ur + utheta * utheta + uz * uz) * jacobian
    separated_density = (
        mp.pi
        * nu ** (mp.mpf(5) / 2)
        * tau ** ((1 - delta) / 2)
        * state["Ur"] ** 2
        * (1 - delta * Z * Z)
        * d ** (-(3 - delta) / 2)
        + mp.pi
        * nu ** (mp.mpf(5) / 2)
        * tau ** ((1 - 3 * delta) / 2)
        * (2 * R * state["F"] ** 2 + state["Uz"] ** 2)
        * (1 - delta * Z * Z)
        * d ** (-(3 - 3 * delta) / 2)
    )
    formula_velocity = (
        mp.sqrt(nu / q) * state["Ur"],
        mp.sqrt(nu) * q ** (-(1 + delta) / 2) * mp.sqrt(2 * R) * state["F"],
        mp.sqrt(nu) * q ** (-(1 + delta) / 2) * state["Uz"],
    )
    return {
        "R": _text(R, digits),
        "Z": _text(Z, digits),
        "tau": signed_log(tau, digits),
        "region": back.get("region"),
        "R_relative_roundtrip_error": signed_log(
            _relative(back["R"], R), digits
        ),
        "Z_absolute_roundtrip_error": signed_log(abs(back["Z"] - Z), digits),
        "velocity_relative_roundtrip_errors": [
            signed_log(_relative(actual, expected), digits)
            for actual, expected in zip(back["cylindrical_velocity"], formula_velocity)
        ],
        "physical_integrand": signed_log(physical_density, digits),
        "separated_integrand": signed_log(separated_density, digits),
        "physical_vs_separated_relative_error": signed_log(
            _relative(physical_density, separated_density), digits
        ),
    }


def run() -> dict[str, object]:
    print("building joined field", flush=True)
    joined = build_joined_field()
    inner = joined.inner
    core = inner.core
    precision = max(int(inner.precision), 260)
    digits = min(precision - 12, 120)

    with mp.workdps(precision):
        nu = mp.mpf(".01")
        delta = mp.mpf(core.delta)
        R_max = mp.mpf(inner.r)
        z_min = mp.mpf("-.5")
        z_max = mp.mpf(".5")
        logtau_values = [
            mp.mpf(0),
            -mp.mpf(2) / delta,
            -mp.mpf(4) / delta,
            -mp.mpf(6) / delta,
        ]

        orders = ((16, 12), (32, 24))
        quadrature_receipts = []
        finest_data = None
        for radial_order, z_order in orders:
            print(
                f"sampling core radial_order={radial_order} z_order={z_order}",
                flush=True,
            )
            samples = _samples(core, R_max, radial_order, z_order)
            shape_radial, shape_swirl_axial = _shape_integrals(samples, delta)
            raw_rows = []
            for logtau in logtau_values:
                tau = mp.exp(logtau)
                direct_radial, direct_swirl_axial = _direct_integral(
                    samples, tau=tau, delta=delta, nu=nu
                )
                separated_radial, separated_swirl_axial = _separated_energy(
                    shape_radial,
                    shape_swirl_axial,
                    tau=tau,
                    delta=delta,
                    nu=nu,
                )
                raw_rows.append(
                    {
                        "logtau": logtau,
                        "direct_radial": direct_radial,
                        "direct_swirl_axial": direct_swirl_axial,
                        "separated_radial": separated_radial,
                        "separated_swirl_axial": separated_swirl_axial,
                    }
                )
            baseline = raw_rows[0]
            rows = [
                _energy_row(
                    raw["direct_radial"],
                    raw["direct_swirl_axial"],
                    raw["separated_radial"],
                    raw["separated_swirl_axial"],
                    logtau=raw["logtau"],
                    delta=delta,
                    baseline_radial=baseline["direct_radial"],
                    baseline_swirl_axial=baseline["direct_swirl_axial"],
                    baseline_total=baseline["direct_radial"]
                    + baseline["direct_swirl_axial"],
                    digits=digits,
                )
                for raw in raw_rows
            ]
            quadrature_receipts.append(
                {
                    "radial_order": radial_order,
                    "z_order": z_order,
                    "sample_count": len(samples),
                    "shape_radial": signed_log(shape_radial, digits),
                    "shape_swirl_plus_axial": signed_log(shape_swirl_axial, digits),
                    "rows": rows,
                }
            )
            finest_data = {
                "samples": samples,
                "shape_radial": shape_radial,
                "shape_swirl_axial": shape_swirl_axial,
            }

        # One independent physical-callable point check at a moderate time.
        physical_check = _physical_integrand_check(
            inner,
            core,
            R=R_max * mp.mpf(".37"),
            Z=mp.mpf(".2"),
            tau=mp.mpf(1),
            nu=nu,
            delta=delta,
            precision=precision,
            digits=digits,
        )

        fine = quadrature_receipts[-1]
        receipt: dict[str, object] = {
            "kind": "bounded_regular_core_energy_scaling",
            "status": "measured",
            "field": "build_joined_field().inner.core",
            "precision_digits": precision,
            "domain": {
                "R_min": "0",
                "R_max": _text(R_max, digits),
                "R_max_source": "inner.r (regular-core exit R_a)",
                "core_polynomial_limit": _text(mp.mpf("4.1") / core.Lambda, digits),
                "Z_min": _text(z_min, digits),
                "Z_max": _text(z_max, digits),
                "normalized_R": "R / inner.r in [0, 1]",
                "nu": _text(nu, digits),
                "delta": _text(delta, digits),
            },
            "jacobian": (
                "2*pi*nu^(3/2)*q^((3-delta)/2)*(1-delta*Z^2)/(1-Z^2), "
                "q=tau/(1-Z^2)"
            ),
            "energy_decomposition": {
                "radial_speed_squared": "nu*q^(-1)*Ur^2",
                "swirl_plus_axial_speed_squared": (
                    "nu*q^(-(1+delta))*((2*R)*F^2 + Uz^2)"
                ),
                "radial_tau_scale": "tau^((1-delta)/2)",
                "swirl_plus_axial_tau_scale": "tau^((1-3*delta)/2)",
            },
            "logtau_values": [_text(value, digits) for value in logtau_values],
            "quadrature": {
                "orders": "radial 16/32; Z 12/24 Gauss-Legendre",
                "interpretation": (
                    "Bounded sample refinement only; no uniform quadrature or "
                    "global energy proof."
                ),
                "coarse": quadrature_receipts[0],
                "fine": fine,
            },
            "physical_callable_integrand_check": physical_check,
            "scope": {
                "measured": [
                    "actual core.evaluate values Ur, F, and Uz",
                    "finite-core kinetic-energy quadrature",
                    "radial versus swirl-plus-axial tau scaling",
                    "one physical callable velocity/integrand roundtrip",
                ],
                "not_certified": [
                    "global finite energy",
                    "outer or heat-tail contribution",
                    "mean-tail cancellation or modification",
                    "recursive scale transition",
                    "continuous uniform-in-tau bounds",
                    "NS closure",
                ],
            },
        }

    OUT_JSON.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {OUT_JSON}", flush=True)
    print("finite-core rows: 4; global finite energy: not certified", flush=True)
    return receipt


if __name__ == "__main__":
    run()
