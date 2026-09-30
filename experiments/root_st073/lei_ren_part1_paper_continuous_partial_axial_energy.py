"""Partial actual axial-square moments through the source pulse and end bumps.

This provider is deliberately separate from ``ContinuousAxialEnergyMoments``.
It evaluates the actual cumulative axial square from the retained inner join
through ``R_v`` using the same live runtime pulse, bump basis, coefficients and
coefficient tangent.  Partial integrals are positive MP Gauss--Legendre sums;
their quadrature error is therefore measured only, not enclosed.

The pulse and end supports are disjoint.  With

    xi = mu * log(R / R_p),       t = log(R / R_v),

the pulse is supported on ``0 <= xi <= 11`` and the end bumps on
``-3.15 <= t <= -0.85``.  Since ``log(R_v/R_p) = 13/mu``, an end point has
``xi = 13 + mu*t > 11``.  Consequently there is no pulse/end cross term in
the axial square.

This module does not install a second global provider and does not claim a
quadrature enclosure, complete five-moment closure, finite energy, or scale
recursion.  A caller handling ``R >= R_v`` should use the existing complete
provider and its retained terminal residuals.
"""

from __future__ import annotations

import json
from decimal import Decimal
from pathlib import Path
from typing import Any, Callable

import mpmath as mp

from lei_ren_part1_paper_axial_correction import signed_log
from lei_ren_part1_paper_joined_outer import _mp


HERE = Path(__file__).resolve().parent


def _relative(actual: mp.mpf, expected: mp.mpf) -> mp.mpf:
    """Return a scale-safe relative difference."""

    actual = mp.mpf(actual)
    expected = mp.mpf(expected)
    scale = max(abs(actual), abs(expected))
    return abs(actual - expected) / scale if scale else abs(actual - expected)


def _fourth_difference(
    function: Callable[[Any], mp.mpf], center: Any, profile: Any, step: str = "1e-5"
) -> mp.mpf:
    """Differentiate a log-radius primitive on exact Decimal offsets."""

    h = mp.mpf(step)
    values = {
        index: function(profile.log_at(center, mp.nstr(index * h, profile.precision)))
        for index in (-2, -1, 1, 2)
    }
    return (values[-2] - 8 * values[-1] + 8 * values[1] - values[2]) / (12 * h)


class ContinuousPartialAxialEnergy:
    """Actual cumulative axial-square primitive on the incoming/pulse/end path.

    The public ``moments_jet`` method accepts the retained incoming region
    beginning at ``R_h`` and ends at ``R_v``.  The task's ``R_p`` to ``R_v``
    interval is the pulse/end portion; supporting the incoming prefix here
    makes the single inner offset explicit and gives a useful ``R_p``
    boundary value without reconstructing that offset a second time.
    """

    def __init__(self, profile: Any, *, quadrature_order: int = 96) -> None:
        try:
            order = int(quadrature_order)
        except (TypeError, ValueError) as exc:
            raise TypeError("quadrature_order must be an integer") from exc
        if order < 16:
            raise ValueError("quadrature_order must be at least 16")
        self.profile = profile
        self.precision = int(profile.precision)
        self.order = order
        # Keep MP nodes and weights.  Converting a node through binary64 would
        # be unnecessary loss in the positive pulse/end atoms.
        with mp.workdps(self.precision):
            nodes, weights = mp.gauss_quadrature(order, "legendre")
            self._nodes = [mp.mpf(node) for node in nodes]
            self._weights = [mp.mpf(weight) for weight in weights]
        # Runtime objects are shared by ``ContinuousIncomingProfile`` for a
        # fixed Z.  Cache only derived positive atoms; no coefficient or
        # residual is reconstructed here.  The cache is especially useful for
        # the startup pulse, whose point value contains a nested endpoint
        # quadrature.
        self._inner_offset_cache: dict[str, tuple[mp.mpf, mp.mpf]] = {}
        self._pulse_cache: dict[tuple[int, str], tuple[mp.mpf, str]] = {}
        self._pulse_startup_cache: dict[int, mp.mpf] = {}
        self._end_atom_cache: dict[tuple[int, str], tuple[list[mp.mpf], str]] = {}
        self._moment_cache: dict[tuple[str, str], dict[str, Any]] = {}

    def _gauss_integral(
        self, function: Callable[[mp.mpf], mp.mpf], left: mp.mpf, right: mp.mpf
    ) -> mp.mpf:
        """Positive MP Gauss--Legendre quadrature on one finite interval."""

        left = mp.mpf(left)
        right = mp.mpf(right)
        if right <= left:
            return mp.mpf(0)
        midpoint = (left + right) / 2
        half_width = (right - left) / 2
        terms = [
            weight * function(midpoint + half_width * node)
            for node, weight in zip(self._nodes, self._weights)
        ]
        return mp.fsum(terms) * half_width

    def _inner_offsets(self, z: mp.mpf) -> tuple[mp.mpf, mp.mpf]:
        """Return the measured inner axial offsets exactly once."""

        cache_key = mp.nstr(z, self.precision)
        if cache_key in self._inner_offset_cache:
            return self._inner_offset_cache[cache_key]
        source = self.profile.seed_source
        inner = source.inner.evaluate_x(mp.e, z)
        radius = _mp(source.Rh)
        offset = _mp(inner["raw_quadratic_integrals"]["axial"]) - 16 * z * z * radius
        offset_z = (
            _mp(inner["raw_quadratic_integrals"]["axial_Z"])
            - 32 * z * radius
        )
        result = (offset, offset_z)
        self._inner_offset_cache[cache_key] = result
        return result

    def _incoming_prefix(self, log_radius: Any, z: mp.mpf) -> tuple[mp.mpf, mp.mpf]:
        """Return actual inner-seeded axial square at a pre-pulse radius."""

        profile = self.profile
        schedule = profile.schedule
        y = _mp(str(profile.offset(log_radius, schedule.logRref)))
        # The incoming provider is linear in Z.  Evaluating its unit-Z atom
        # and applying ZÂ² keeps the same normalization as the complete live
        # provider and the current coefficient solve.
        unit = profile.incoming_provider.I_uz2(y, mp.mpf(1))
        radius_ref = mp.exp(_mp(str(schedule.logRref)))
        offset, offset_z = self._inner_offsets(z)
        return (
            radius_ref * unit * z * z + offset,
            2 * radius_ref * unit * z + offset_z,
        )

    def _pulse_integral(self, runtime: Any, xi: mp.mpf) -> tuple[mp.mpf, str]:
        """Return the normalized positive pulse square through ``xi``."""

        xi = mp.mpf(xi)
        cache_key = (id(runtime), mp.nstr(xi, self.precision))
        if cache_key in self._pulse_cache:
            return self._pulse_cache[cache_key]
        if xi <= 0:
            result = (mp.mpf(0), "before_pulse_support")
            self._pulse_cache[cache_key] = result
            return result
        if xi >= 11:
            # This is the live complete atom used by the algebraic energy
            # solve.  Reusing it makes the Rv boundary an identity up to MP
            # arithmetic instead of a second nominal atom.
            result = (runtime.Kp, "runtime_complete_pulse_atom")
            self._pulse_cache[cache_key] = result
            return result

        pulse = runtime.pulse

        def weighted_startup(value: mp.mpf) -> mp.mpf:
            source_value = pulse.value_jet(value)["value"]
            return mp.exp(-2 * value) * source_value * source_value

        if xi < mp.mpf(".02"):
            result = (
                self._gauss_integral(weighted_startup, mp.mpf(0), xi),
                "positive_gauss_partial_pulse_startup",
            )
            self._pulse_cache[cache_key] = result
            return result

        # The nested startup primitive is needed only once per live runtime.
        # On the plateau/cutoff, the shared pulse definition is explicit and
        # can be evaluated without invoking that nested primitive again.
        if id(runtime) not in self._pulse_startup_cache:
            self._pulse_startup_cache[id(runtime)] = self._gauss_integral(
                weighted_startup, mp.mpf(0), mp.mpf(".02")
            )

        def weighted_explicit(value: mp.mpf) -> mp.mpf:
            primitive = value - mp.mpf(".01")
            cutoff = pulse.sigma_pair(mp.mpf(11) - value)[0]
            source_value = primitive * cutoff
            return mp.exp(-2 * value) * source_value * source_value

        cuts = [mp.mpf(".02")]
        if mp.mpf(10) < xi:
            cuts.append(mp.mpf(10))
        cuts.append(xi)
        total = self._pulse_startup_cache[id(runtime)]
        for left, right in zip(cuts, cuts[1:]):
            total += self._gauss_integral(weighted_explicit, left, right)
        result = (total, "positive_gauss_partial_pulse")
        self._pulse_cache[cache_key] = result
        return result

    def _end_atoms(
        self, runtime: Any, end: mp.mpf
    ) -> tuple[list[mp.mpf], str]:
        """Return the two positive normalized end atoms through ``end``."""

        end = mp.mpf(end)
        cache_key = (id(runtime), mp.nstr(end, self.precision))
        cached = self._end_atom_cache.get(cache_key)
        if cached is not None:
            return cached
        basis = runtime.basis
        ell = _mp(basis.ell)
        mu = _mp(runtime.mu)
        atoms: list[mp.mpf] = []
        partial = False
        for index, center in enumerate((mp.mpf(-3), mp.mpf(-1))):
            upper = end - center
            if upper <= -ell:
                atom = mp.mpf(0)
            elif upper >= ell:
                # K already contains exp(-26) and the translated
                # exp(-2*mu*center) factor for this same live basis atom.
                atom = _mp(runtime.K[index])
            else:
                partial = True

                def weighted_square(s: mp.mpf) -> mp.mpf:
                    beta = basis.raw(s / ell) / (ell * basis.normalizer)
                    return mp.exp(-2 * mu * s) * beta * beta

                gram = self._gauss_integral(weighted_square, -ell, upper)
                atom = mp.exp(-26 - 2 * mu * center) * gram
            atoms.append(atom)
        result = (
            atoms,
            "positive_gauss_partial_end" if partial else "runtime_complete_end_atom",
        )
        self._end_atom_cache[cache_key] = result
        return result

    def _end_partial(
        self, runtime: Any, coefficients: list[mp.mpf], end: mp.mpf
    ) -> tuple[mp.mpf, str]:
        """Return the normalized end-bump square through ``t = end``."""

        atoms, method = self._end_atoms(runtime, end)
        total = mp.fsum(
            coefficient * coefficient * atom
            for coefficient, atom in zip(coefficients, atoms)
        )
        return total, method

    def _end_bilinear(
        self,
        runtime: Any,
        coefficients: list[mp.mpf],
        tangent: list[mp.mpf],
        end: mp.mpf,
    ) -> mp.mpf:
        """Return the Z derivative of the end square via polarization."""

        atoms, _ = self._end_atoms(runtime, end)
        return 2 * mp.fsum(
            coefficient * derivative * atom
            for coefficient, derivative, atom in zip(coefficients, tangent, atoms)
        )

    def moments_jet(self, log_radius: Any, Z: Any) -> dict[str, Any]:
        """Evaluate the actual cumulative axial square and its Z jet.

        The result is defined for ``R_h <= R <= R_v``.  A point before
        ``R_p`` uses the live incoming prefix; a point after ``R_p`` adds the
        pulse and, when supported, the translated end bumps.  The complete
        provider remains responsible for ``R >= R_v``.
        """

        profile = self.profile
        schedule = profile.schedule
        source = profile.seed_source
        cache_key = (str(log_radius), mp.nstr(_mp(Z), self.precision))
        cached = self._moment_cache.get(cache_key)
        if cached is not None:
            return cached
        with mp.workdps(self.precision):
            z = _mp(Z)
            if abs(z) >= 1:
                raise ValueError("Actual inner seed requires |Z| < 1")
            start = _mp(str(profile.offset(log_radius, schedule.logR_p)))
            end = _mp(str(profile.offset(log_radius, schedule.logR_v)))
            # ``CorrectedProfile.offset`` operates on Decimal checkpoints;
            # the joined source stores ``logRh`` as an MP value.  Preserve
            # the declared digits at this boundary as the other installed
            # moment providers do.
            from_join = _mp(
                str(
                    profile.offset(
                        log_radius,
                        Decimal(mp.nstr(source.logRh, self.precision)),
                    )
                )
            )
            if from_join < 0:
                raise ValueError("Partial axial square starts at the retained inner join Rh")
            if end > 0:
                raise ValueError("Partial axial square stops at Rv; use the complete provider after Rv")

            runtime = profile.runtime(float(z))
            tangent = profile.coefficient_tangent(float(z))
            mu = _mp(runtime.mu)
            if 13 - mp.mpf("3.15") * mu <= 11:
                raise ArithmeticError("Pulse/end support separation failed")

            incoming, incoming_z = self._incoming_prefix(log_radius, z)
            pulse = mp.mpf(0)
            pulse_z = mp.mpf(0)
            end_energy = mp.mpf(0)
            end_energy_z = mp.mpf(0)
            scale = mp.mpf(0)
            lz = -2 * z / (1 + z * z)
            xi = mu * start
            pulse_method = "incoming_prefix"
            end_method = "incoming_prefix"
            region = "incoming"

            if start >= 0:
                ep_log = _mp(
                    schedule.at_log_radius(
                        schedule.logR_p, z
                    )["log_angular_amplitude"]
                )
                scale = mp.exp(_mp(str(schedule.logR_p)) + 2 * ep_log)
                pulse_atom, pulse_method = self._pulse_integral(runtime, xi)
                pulse = runtime.a * runtime.a * pulse_atom / mu
                pulse_z = 2 * runtime.a * tangent["a_Z"] * pulse_atom / mu
                region = "pulse"
                if end >= mp.mpf("-3.15"):
                    end_energy, end_method = self._end_partial(
                        runtime, [_mp(value) for value in runtime.c], end
                    )
                    end_energy_z = self._end_bilinear(
                        runtime,
                        [_mp(value) for value in runtime.c],
                        [_mp(value) for value in tangent["c_Z"]],
                        end,
                    )
                    region = "end_bump" if end < 0 else "terminal_Rv"

            axial = incoming + scale * (pulse + end_energy)
            axial_z = incoming_z + scale * (
                pulse_z + end_energy_z + 2 * lz * (pulse + end_energy)
            )
            angular = profile.angular_moments_jet(log_radius, z)
            z_theta = axial - angular["swirl_energy"] / 2
            z_theta_z = axial_z - angular["swirl_energy_Z"] / 2
            result = {
                "axial_square": axial,
                "axial_square_Z": axial_z,
                "z_theta": z_theta,
                "z_theta_Z": z_theta_z,
                "incoming_axial_square": incoming,
                "incoming_axial_square_Z": incoming_z,
                "pulse_normalized": pulse,
                "pulse_normalized_Z": pulse_z,
                "end_normalized": end_energy,
                "end_normalized_Z": end_energy_z,
                "scale_Rp_Ep_square": scale,
                "xi": xi,
                "end_offset": end,
                "region": region,
                "pulse_atom_method": pulse_method,
                "end_atom_method": end_method,
                "quadrature_order": self.order,
                "positive_partial_quadrature": True,
                "quadrature_error_enclosed": False,
                "actual_inner_offsets_reapplied_once": True,
                "pulse_end_cross_term_zero_by_disjoint_support": True,
                "terminal_residual_retained": True,
                "terminal_residual_forced_zero": False,
                "complete_five_moments": False,
                "finite_energy_certified": False,
                "scale_recursion_certified": False,
            }
            self._moment_cache[cache_key] = result
            return result


def _signed(value: Any, precision: int) -> dict[str, Any]:
    return signed_log(_mp(value), precision)


def _markdown(report: dict[str, Any]) -> str:
    boundary = report["boundary_agreement"]
    rows = report["samples"]
    lines = [
        "# Continuous partial axial-square energy",
        "",
        "This receipt evaluates the actual inner-seeded axial-square primitive through `R_v` with the live continuous pulse, end basis, coefficients, and coefficient tangent.",
        "",
        "## Normalization",
        "",
        "For `xi = mu log(R/R_p)` and `t = log(R/R_v)`, the pulse contribution is",
        "",
        "```text",
        "R_p E_p^2 * a^2/mu * integral_0^xi exp(-2 s) gp(s)^2 ds.",
        "```",
        "",
        "The end contribution is `R_p E_p^2` times the sum of the positive partial bump atoms. Each atom contains `exp(-26 - 2 mu center)` and `integral exp(-2 mu s) beta(s)^2 ds`. The pulse support is `0 <= xi <= 11`; end support is `-3.15 <= t <= -.85`, so `xi = 13 + mu t > 11` on the end band and the cross term is zero by support.",
        "",
        "The incoming prefix is `R_ref I_uz2(y,1) Z^2` plus the measured inner axial offset, with its analytic Z derivative. The offset is applied once and is never rounded to zero.",
        "The normalized Z jet uses `2 a a_Z` for the pulse and `2 sum_i c_i c_i_Z K_i_partial` for the end atoms. The physical scale contributes `scale_Z = 2 LZ scale`, with `LZ = -2 Z/(1+Z^2)`.",
        "",
        "## Numerical checks",
        "",
        f"- quadrature order: `{report['quadrature_order']}` (positive MP Gauss--Legendre; enclosure: `{report['quadrature_error_enclosed']}`)",
        f"- samples: `{', '.join(rows)}`",
        f"- maximum Rv boundary relative difference: `{boundary['maximum_relative_difference']}`",
        f"- maximum checked radial/Z integrand relative difference: `{report['integrand_checks']['maximum_relative_difference']}`",
        "",
        "The boundary compares this partial provider with the existing complete provider at `R_v` for axial square, axial-square Z jet, `z_theta`, and `z_theta` Z jet. Forward checks compare fourth-order log-radius differences with `R U_z^2`, `2 R U_z U_{z,Z}`, `R(U_z^2-U_theta^2/2)`, and `R(2 U_z U_{z,Z}-U_theta U_{theta,Z})`. At the two end centers, axial square and axial-square Z use the separately retained normalized end atoms; direct whole-sum stencils are retained as precision-loss diagnostics because the full accumulated pulse baseline overwhelms the local end atom at the working precision.",
        "",
        "## Limits",
        "",
        "The result is a nominal numerical primitive. Quadrature error, inherited inner/incoming uncertainty, coefficient-input uncertainty, full five-moment closure, heat continuation, finite energy, and scale recursion remain uncertified. The complete provider should be called separately for `R >= R_v`; terminal residuals remain materialized.",
        "",
    ]
    return "\n".join(lines)


def run() -> dict[str, Any]:
    """Build the shared field, run representative checks, and save receipts."""

    from lei_ren_part1_paper_continuous_incoming_outer import ContinuousIncomingOuterField
    from lei_ren_part1_paper_joined_outer import build_joined_field

    folder = Path(__file__).resolve().parent
    source = build_joined_field()
    seeded = json.loads(
        (folder / "lei_ren_part1_paper_seeded_shared_candidate.json").read_text(
            encoding="utf-8"
        )
    )
    atoms = json.loads(
        (folder / "lei_ren_part1_paper_continuous_axial_solve.json").read_text(
            encoding="utf-8"
        )
    )
    field = ContinuousIncomingOuterField(source, prepared={0.3: (seeded, atoms)})
    profile = field.outer
    engine = ContinuousPartialAxialEnergy(profile)

    with mp.workdps(max(field.precision, profile.jet_precision)):
        z = mp.mpf(".3")
        schedule = profile.schedule
        mu = _mp(str(schedule.mu))
        points = {
            "incoming_reference": profile.log_at(schedule.logRref, "0"),
            "pulse_bulk_xi_5": profile.log_at(
                schedule.logR_p, mp.nstr(5 / mu, profile.precision)
            ),
            "end_bump_1_center": profile.log_at(schedule.logR_v, "-3"),
            "end_bump_2_center": profile.log_at(schedule.logR_v, "-1"),
            "terminal_Rv": schedule.logR_v,
        }
        samples = {
            name: engine.moments_jet(radius, z) for name, radius in points.items()
        }

        complete = profile.axial_energy_moment_provider.moments_jet(
            schedule.logR_v, z
        )
        boundary_keys = (
            "axial_square",
            "axial_square_Z",
            "z_theta",
            "z_theta_Z",
        )
        boundary_errors = {
            key: mp.nstr(_relative(samples["terminal_Rv"][key], complete[key]), 50)
            for key in boundary_keys
        }
        boundary_max = max(mp.mpf(value) for value in boundary_errors.values())

        integrand_rows: list[dict[str, Any]] = []
        for name in (
            "incoming_reference",
            "pulse_bulk_xi_5",
            "end_bump_1_center",
            "end_bump_2_center",
        ):
            radius = points[name]
            row = samples[name]
            values = profile.values_with_jets(radius, z)
            physical_radius = mp.exp(_mp(str(radius)))
            targets = {
                "axial_square": physical_radius * values["Uz"] ** 2,
                "axial_square_Z": 2 * physical_radius * values["Uz"] * values["Uz_Z"],
                "z_theta": physical_radius
                * (values["Uz"] ** 2 - values["Utheta"] ** 2 / 2),
                "z_theta_Z": physical_radius
                * (
                    2 * values["Uz"] * values["Uz_Z"]
                    - values["Utheta"] * values["Utheta_Z"]
                ),
            }

            local_end_residual = name.startswith("end_bump_")
            lz = -2 * z / (1 + z * z)
            whole_sum_errors: dict[str, str] = {}
            normalized_end_errors: dict[str, str] = {}
            if local_end_residual:
                # The accumulated axial row is many orders larger than an
                # end atom.  Keep the direct whole-sum stencil as a diagnostic
                # for this precision loss, but validate the end derivative
                # with the owned normalized atoms instead.
                whole_sum_derivatives = {
                    key: _fourth_difference(
                        lambda query, key=key: engine.moments_jet(query, z)[key],
                        radius,
                        profile,
                    )
                    for key in ("axial_square", "axial_square_Z")
                }
                whole_sum_errors = {
                    key: mp.nstr(
                        _relative(whole_sum_derivatives[key], targets[key]), 50
                    )
                    for key in whole_sum_derivatives
                }
                normalized_derivatives = {
                    key: _fourth_difference(
                        lambda query, key=key: engine.moments_jet(query, z)[key],
                        radius,
                        profile,
                    )
                    for key in ("end_normalized", "end_normalized_Z")
                }
                scale = row["scale_Rp_Ep_square"]
                normalized_targets = {
                    "end_normalized": targets["axial_square"] / scale,
                    # end_normalized_Z excludes the Z derivative of the
                    # outer scale, so remove 2 LZ times the square target.
                    "end_normalized_Z": (
                        targets["axial_square_Z"] / scale
                        - 2 * lz * targets["axial_square"] / scale
                    ),
                }
                normalized_end_errors = {
                    key: mp.nstr(
                        _relative(normalized_derivatives[key], normalized_targets[key]),
                        50,
                    )
                    for key in normalized_derivatives
                }
                errors = {
                    "axial_square": normalized_end_errors["end_normalized"],
                    "axial_square_Z": normalized_end_errors["end_normalized_Z"],
                    "z_theta": mp.nstr(
                        _relative(
                            _fourth_difference(
                                lambda query: engine.moments_jet(query, z)["z_theta"],
                                radius,
                                profile,
                            ),
                            targets["z_theta"],
                        ),
                        50,
                    ),
                    "z_theta_Z": mp.nstr(
                        _relative(
                            _fourth_difference(
                                lambda query: engine.moments_jet(query, z)["z_theta_Z"],
                                radius,
                                profile,
                            ),
                            targets["z_theta_Z"],
                        ),
                        50,
                    ),
                }
            else:
                derivatives = {
                    key: _fourth_difference(
                        lambda query, key=key: engine.moments_jet(query, z)[key],
                        radius,
                        profile,
                    )
                    for key in targets
                }
                errors = {
                    key: mp.nstr(_relative(derivatives[key], targets[key]), 50)
                    for key in targets
                }
            integrand_rows.append(
                {
                    "name": name,
                    "log_radius": str(radius),
                    "region": row["region"],
                    "axial_derivative_check_method": (
                        "normalized_end_atoms_with_whole_sum_precision_diagnostic"
                        if local_end_residual
                        else "direct_cumulative_fourth_difference"
                    ),
                    "derivative_relative_errors": errors,
                    "whole_sum_derivative_relative_errors": whole_sum_errors,
                    "normalized_end_derivative_relative_errors": normalized_end_errors,
                }
            )

        max_integrand = max(
            mp.mpf(error)
            for row in integrand_rows
            for error in row["derivative_relative_errors"].values()
        )
        runtime = profile.runtime(float(z))
        terminal = {
            str(index): _signed(
                runtime.component.terminal_balance(index, profile.precision)["value"],
                profile.precision,
            )
            for index in (1, 2)
        }
        offset, offset_z = engine._inner_offsets(z)
        report: dict[str, Any] = {
            "kind": "continuous_partial_axial_square_energy",
            "module": str(Path(__file__).resolve()),
            "Z": ".3",
            "precision": profile.precision,
            "quadrature_order": engine.order,
            "quadrature_error_enclosed": False,
            "schedule": {
                "mu": mp.nstr(mu, profile.precision),
                "logR_p": str(schedule.logR_p),
                "logR_v": str(schedule.logR_v),
                "logRv_minus_logRp": str(schedule.logR_v - schedule.logR_p),
                "end_support": ["-3.15", "-.85"],
                "pulse_support_xi": ["0", "11"],
            },
            "formulas": {
                "incoming": "R_ref*I_uz2(y,1)*Z^2 + (inner_raw_axial - 16*Z^2*R_h)",
                "incoming_Z": "2*R_ref*I_uz2(y,1)*Z + (inner_raw_axial_Z - 32*Z*R_h)",
                "pulse": "R_p*E_p^2*a^2/mu * integral_0^xi exp(-2*s)*gp(s)^2 ds",
                "pulse_Z_normalized": "2*a*a_Z/mu * integral_0^xi exp(-2*s)*gp(s)^2 ds",
                "end": "R_p*E_p^2*sum_i c_i^2*exp(-26-2*mu*center_i)*integral_{-ell}^{t-center_i} exp(-2*mu*s)*beta(s)^2 ds",
                "end_Z_normalized": "2*sum_i c_i*c_i_Z*K_i_partial",
                "physical_scale_Z": "scale_Z = 2*LZ*scale, LZ = -2*Z/(1+Z^2)",
                "z_theta": "axial_square - swirl_energy/2 from the installed angular provider",
            },
            "samples": list(points),
            "sample_rows": {
                name: {
                    "region": row["region"],
                    "xi": mp.nstr(row["xi"], profile.precision),
                    "end_offset": mp.nstr(row["end_offset"], profile.precision),
                    "axial_square": _signed(row["axial_square"], profile.precision),
                    "axial_square_Z": _signed(row["axial_square_Z"], profile.precision),
                    "z_theta": _signed(row["z_theta"], profile.precision),
                    "z_theta_Z": _signed(row["z_theta_Z"], profile.precision),
                    "pulse_normalized": _signed(row["pulse_normalized"], profile.precision),
                    "end_normalized": _signed(row["end_normalized"], profile.precision),
                    "pulse_atom_method": row["pulse_atom_method"],
                    "end_atom_method": row["end_atom_method"],
                }
                for name, row in samples.items()
            },
            "boundary_agreement": {
                "relative_errors": boundary_errors,
                "maximum_relative_difference": mp.nstr(boundary_max, 50),
                "compared_with": "ContinuousAxialEnergyMoments.moments_jet at Rv",
            },
            "integrand_checks": {
                "rows": integrand_rows,
                "maximum_relative_difference": mp.nstr(max_integrand, 50),
                "whole_sum_end_precision_loss_diagnostic": True,
                "radial_targets": "d_logR axial_square = R*Uz^2; d_logR z_theta = R*(Uz^2-Utheta^2/2)",
                "Z_targets": "d_logR axial_square_Z = 2*R*Uz*Uz_Z; d_logR z_theta_Z = R*(2*Uz*Uz_Z-Utheta*Utheta_Z)",
            },
            "inner_offsets": {
                "actual_inner_offsets_reapplied_once": True,
                "axial_offset": _signed(offset, profile.precision),
                "axial_offset_Z": _signed(offset_z, profile.precision),
            },
            "terminal_residuals": {
                "rows": terminal,
                "terminal_residual_retained": True,
                "terminal_residual_forced_zero": False,
            },
            "live_identity": {
                "runtime_pulse_shared": True,
                "runtime_basis_shared": True,
                "runtime_coefficients_shared": True,
                "coefficient_tangent_shared": True,
                "pulse_end_cross_term_zero_by_disjoint_support": True,
            },
            "claims": {
                "partial_axial_energy_implemented": True,
                "quadrature_enclosure_certified": False,
                "complete_five_moments": False,
                "finite_energy_certified": False,
                "scale_recursion_certified": False,
            },
        }

    output_json = Path(__file__).with_suffix(".json")
    output_md = Path(__file__).with_suffix(".md")
    output_json.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    output_md.write_text(_markdown(report), encoding="utf-8")
    print(
        json.dumps(
            {
                "boundary_maximum_relative_difference": report["boundary_agreement"][
                    "maximum_relative_difference"
                ],
                "integrand_maximum_relative_difference": report["integrand_checks"][
                    "maximum_relative_difference"
                ],
                "quadrature_error_enclosed": False,
            }
        ),
        flush=True,
    )
    return report


if __name__ == "__main__":
    run()
