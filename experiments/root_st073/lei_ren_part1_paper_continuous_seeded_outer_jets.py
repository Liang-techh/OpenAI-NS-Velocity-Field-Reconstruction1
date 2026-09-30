"""Bounded Z-jet adapter for the continuous seeded axial exterior.

The adapter differentiates the continuous pulse and end-bump contribution
with respect to Z after the actual Rh-seeded input rows have been transported
through the implicit continuous coefficient solve.  The angular schedule is
retained by identity; its future energy derivative is still the explicit
float-backed finite difference supplied by ``actual_seeded_input_tangents``.

This file is a diagnostic installation.  It does not claim a complete outer
five-moment field, a global Z derivative, finite energy, or scale recursion.
"""

from decimal import Decimal, localcontext
from functools import lru_cache
import json
from pathlib import Path

import mpmath as mp

from lei_ren_part1_paper_axial_correction import from_signed_log, signed_log
from lei_ren_part1_paper_continuous_axial_solve import ContinuousAxialCorrection
from lei_ren_part1_paper_continuous_axial_tangents import ContinuousAxialAlgebra
from lei_ren_part1_paper_continuous_seeded_outer import (
    ContinuousSeededAxialProfile,
    ContinuousSeededOuterField,
)
from lei_ren_part1_paper_joined_outer import (
    _finite_difference,
    _mp,
    build_joined_field,
)
from lei_ren_part1_paper_seeded_input_tangents import (
    actual_seeded_input_tangents,
)


def _signed(value, precision=120):
    """Serialize an MP value without discarding its arbitrary exponent."""

    with mp.workdps(precision):
        return signed_log(_mp(value), precision)


def _signed_logs(values, precision=120):
    return [_signed(value, precision) for value in values]


class ContinuousSeededAxialProfileJets(ContinuousSeededAxialProfile):
    """Continuous seeded profile with analytic axial Z jets.

    ``actual_seeded_input_tangents`` supplies the two normalized incoming-row
    derivatives and the energy-target derivative.  The latter includes an
    explicitly labelled finite-difference angular tail derivative.  The
    continuous axial algebra then supplies ``a_Z`` and ``c_Z``.
    """

    def __init__(
        self,
        source,
        *,
        prepared=None,
        seeded_receipt=None,
        continuous_receipt=None,
        jet_precision=200,
    ):
        super().__init__(source, prepared=prepared)
        self.jet_precision = int(jet_precision)
        self.seeded_receipt = seeded_receipt
        self.continuous_receipt = continuous_receipt

    def _receipts_for(self, Z):
        """Return the seeded and continuous receipts used at one Z value."""

        z = float(Z)
        prepared_key = {float(key) for key in self.prepared}
        if z in prepared_key:
            return self.prepared[z]
        # An explicitly supplied pair is valid only when no prepared Z map is
        # present.  With a map, never reuse the .3 receipt at another Z.
        if (
            not prepared_key
            and self.seeded_receipt is not None
            and self.continuous_receipt is not None
        ):
            return self.seeded_receipt, self.continuous_receipt
        solved = self.seed_solve(z)
        return solved, solved["axial"]

    def _merged_receipt(self, seeded, continuous):
        """Merge top-level actual energy target with continuous atoms."""

        merged = dict(seeded)
        merged["axial"] = continuous
        # ``ContinuousAxialAlgebra`` reads the actual target from this level;
        # preserve it when the continuous receipt itself has no target key.
        if "energy_target" not in merged:
            target = seeded.get("axial", {}).get("energy_target")
            if target is not None:
                merged["energy_target"] = target
        return merged

    @lru_cache(maxsize=64)
    def coefficient_tangent(self, Z):
        """Replay actual input jets and the continuous implicit solve."""

        z = float(Z)
        seeded, continuous = self._receipts_for(z)
        input_jet = actual_seeded_input_tangents(self.seed_source, seeded, z)
        merged = self._merged_receipt(seeded, continuous)
        algebra = ContinuousAxialAlgebra(merged, continuous)
        with mp.workdps(max(self.jet_precision, algebra.precision, self.precision)):
            base_Z = [from_signed_log(row) for row in input_jet["base_Z"]]
            target_Z = from_signed_log(input_jet["energy_target_Z"])
            tangent = algebra.tangent(base_Z, target_Z)
            return {
                "a_Z": tangent["a_Z"],
                "c_Z": list(tangent["c_Z"]),
                "input": input_jet,
                "algebra_precision": algebra.precision,
                "linear_tangent_relative_replay": tangent[
                    "linear_tangent_relative_replay"
                ],
                "energy_tangent_relative_replay": tangent[
                    "energy_tangent_relative_replay"
                ],
            }

    def _base_derivative(self, Z):
        tangent = self.coefficient_tangent(float(Z))
        return [from_signed_log(row) for row in tangent["input"]["base_Z"]]

    def _pulse_integral(self, component, xi):
        """Return the normalized row-one pulse integral through ``xi``."""

        xi = _mp(xi)
        if xi <= 0:
            return mp.mpf(0)
        atom = (
            component.pulse.partial_row(component.mu, 1, mp.nstr(xi, component.precision))
            if xi < 11
            else component.pulse.full_row(component.mu, 1)
        )
        return mp.exp(mp.mpf(atom["log_normalized_pulse_integral"]))

    def _pulse_integral_and_derivative(self, component, xi, tangent):
        xi = _mp(xi)
        if xi <= 0:
            return mp.mpf(0), mp.mpf(0), {"exact_support_zero": True}
        if xi < 11:
            atom = component.pulse.partial_row(
                component.mu, 1, mp.nstr(xi, component.precision)
            )
        else:
            # Keep the pulse atom separate; the translated end-bump
            # primitive is added exactly once below for the current end.
            atom = component.pulse.full_row(component.mu, 1)
        integral = mp.exp(mp.mpf(atom["log_normalized_pulse_integral"]))
        # The pulse integral is an atom of the fixed source shape, so its Z
        # jet is zero; only the solved amplitude contributes here.
        return integral, mp.mpf(0), atom

    def _end_primitive_derivative(self, component, c_Z, end):
        with mp.workdps(max(self.jet_precision, component.precision)):
            lam = mp.mpf(".5") - component.mu
            total = mp.mpf(0)
            for coefficient, center in zip(c_Z, (-3, -1)):
                total += coefficient * mp.exp(lam * center) * component.basis.primitive(
                    lam, _mp(end) - center
                )
            return total

    def _normalization_and_jet(self, logR, Z, *, xi=None, end=None):
        """Return the normalized row-one cumulative atom and its Z jet."""

        component = self.component(float(Z))
        tangent = self.coefficient_tangent(float(Z))
        base_Z = self._base_derivative(float(Z))[0]
        a_Z = tangent["a_Z"]
        if xi is None:
            raise ValueError("The normalized pulse coordinate is required")
        if end is not None and _mp(end)>=mp.mpf('-3.15'):
            terminal=component.terminal_balance(1,max(self.precision,self.jet_precision))
            integral=terminal['pulse_integral']
            N=terminal['value']-component.weighted_tail(1,_mp(end))
            terminal_Z=base_Z+a_Z*integral+self._end_primitive_derivative(component,tangent['c_Z'],0)
            N_Z=terminal_Z-component.weighted_tail(1,_mp(end),coefficients=tangent['c_Z'])
            return N,N_Z,terminal['pulse_atom']
        integral, _, atom = self._pulse_integral_and_derivative(component, xi, tangent)
        N = component.base[0] + component.a * integral
        N_Z = base_Z + a_Z * integral
        primitive = mp.mpf(0)
        primitive_Z = mp.mpf(0)
        if end is not None and _mp(end) >= mp.mpf("-3.15"):
            primitive = component.weighted_primitive(1, _mp(end))
            primitive_Z = self._end_primitive_derivative(component, tangent["c_Z"], end)
            N += primitive
            N_Z += primitive_Z
        return N, N_Z, atom

    def values_with_jets(self, logR, Z):
        """Return values and the continuous axial Z derivative."""

        with mp.workdps(max(self.precision, self.jet_precision)):
            base = self.schedule.at_log_radius(logR, Z)
            inherited = super().values(logR, Z)
            Utheta = _mp(inherited["Utheta"])
            # Through Rv the amplitude has the exact source factor
            # (1+Z^2)^-1. Do not round its derivative through binary64.
            z_mp = _mp(Z)
            LZ = -2*z_mp/(1+z_mp*z_mp) if self.offset(logR,self.schedule.logR_v)<=0 else _mp(base["dlogU_dZ"])
            Uz = _mp(inherited["Uz"])
            start = self.offset(logR, self.schedule.logR_p)
            end = self.offset(logR, self.schedule.logR_v)
            method = "schedule_before_Rp_fallback"
            Uz_Z = _mp(base["Uz_Z"])
            if start >= 0 and end <= 0:
                component = self.component(float(Z))
                tangent = self.coefficient_tangent(float(Z))
                with localcontext() as context:
                    context.prec = max(self.schedule.decimal_precision, self.precision)
                    xi_decimal = self.schedule.mu * start
                xi = _mp(xi_decimal)
                if Decimal("-3.15") <= end <= Decimal("-.85"):
                    jet = component.value_jet(str(end))
                    beta = jet["beta"]
                    beta_Z = mp.mpf(0)
                    for coefficient, center in zip(tangent["c_Z"], (-3, -1)):
                        beta_Z += coefficient * component.basis.values(
                            _mp(end) - center
                        )["beta"]
                    Uz = mp.exp(_mp(base["log_angular_amplitude"])) * beta
                    Uz_Z = mp.exp(_mp(base["log_angular_amplitude"])) * (
                        beta_Z + LZ * beta
                    )
                    method = "continuous_end_bump_analytic_Z_jet"
                else:
                    pulse = component.pulse.value_jet(str(xi))
                    gp = pulse["value"]
                    Uz = mp.exp(_mp(base["log_angular_amplitude"])) * component.a * gp
                    Uz_Z = mp.exp(_mp(base["log_angular_amplitude"])) * (
                        tangent["a_Z"] * gp + LZ * component.a * gp
                    )
                    method = "continuous_pulse_analytic_Z_jet"
            elif start >= 0 and end > 0:
                # The axial source is zero after the end bump.  Its cumulative
                # mean is still retained below, so this branch does not imply
                # a zero terminal mass.
                Uz = mp.mpf(0)
                Uz_Z = mp.mpf(0)
                method = "post_Rv_zero_axial_source"
            return {
                "Utheta": Utheta,
                "Utheta_Z": Utheta * LZ,
                "Uz": Uz,
                "Uz_Z": Uz_Z,
                "base_log_Utheta": str(base["log_angular_amplitude"]),
                "logF_slope": _mp(base["logF_slope"]),
                "logarithmic_slope": _mp(base["logarithmic_slope"]),
                "heat_method": inherited.get("heat_method", base.get("heat_method")),
                "axial_Z_method": method,
                "angular_Z_jet_complete": False,
            }

    values_Z = values_with_jets

    def axial_average_jet(self, logR, Z):
        """Return the retained mean and analytic post-Rp mean Z jet."""

        with mp.workdps(max(self.precision, self.jet_precision)):
            start = self.offset(logR, self.schedule.logR_p)
            if start < 0:
                # The incoming construction is retained exactly through the
                # parent profile. Its Z jet still uses the declared source
                # stencil because the incoming primitive has no jet API.
                parent_receipt = ContinuousSeededAxialProfile.axial_average_receipt(
                    self, logR, Z
                )
                value = from_signed_log(parent_receipt["nominal"])
                z = _mp(Z)
                derivative = _finite_difference(
                    lambda zz: from_signed_log(
                        ContinuousSeededAxialProfile.axial_average_receipt(
                            self, logR, float(zz)
                        )["nominal"]
                    ),
                    z,
                    step=self.seed_source.derivative_step,
                )
                return {
                    "value": _mp(value),
                    "derivative": _mp(derivative),
                    "method": "inherited_incoming_float_Z_fallback_before_Rp",
                    "terminal_mean_forced_zero": False,
                }

            end = self.offset(logR, self.schedule.logR_v)
            component = self.component(float(Z))
            tangent = self.coefficient_tangent(float(Z))
            with localcontext() as context:
                context.prec = max(self.schedule.decimal_precision, self.precision)
                xi = _mp(self.schedule.mu * start)
            if end <= 0:
                N, N_Z, atom = self._normalization_and_jet(
                    logR, Z, xi=xi, end=end
                )
                logE = _mp(
                    self.schedule.at_log_radius(logR, Z)["log_angular_amplitude"]
                )
                lam = mp.mpf(".5") - component.mu
                scale = mp.exp(logE - lam * _mp(end))
                z_mp = _mp(Z)
                LZ = -2*z_mp/(1+z_mp*z_mp)
                value = scale * N
                derivative = scale * (N_Z + LZ * N)
                region = "continuous_seeded_exterior_before_Rv"
            else:
                N, N_Z, atom = self._normalization_and_jet(
                    logR, Z, xi=mp.mpf(11), end=mp.mpf(0)
                )
                at_rv = self.schedule.at_log_radius(self.schedule.logR_v, Z)
                logE = _mp(at_rv["log_angular_amplitude"])
                scale = mp.exp(logE - _mp(end))
                z_mp = _mp(Z)
                LZ = -2*z_mp/(1+z_mp*z_mp)
                value = scale * N
                derivative = scale * (N_Z + LZ * N)
                region = "continuous_seeded_exterior_after_Rv"
            relative_bound = atom.get(
                "log_relative_omitted_absolute_bound",
                atom.get("log_relative_omitted_positive_bound"),
            )
            if relative_bound is not None and "log_normalized_pulse_integral" in atom:
                # This is conditional on fixed nominal coefficients and fixed
                # quadrature atoms; incoming, coefficient, and quadrature
                # uncertainty are not enclosed by this source-kernel bound.
                kernel_bound = mp.exp(
                    mp.mpf(atom["log_normalized_pulse_integral"])
                    + mp.mpf(relative_bound)
                )
                derivative_scale = tangent["a_Z"] + LZ * component.a
                kernel_sign = int(atom.get("omitted_correction_sign", 1))
                derivative_sign = kernel_sign * int(mp.sign(derivative_scale))
                omitted = {
                    "pulse_omitted_correction_sign": kernel_sign,
                    "pulse_omitted_derivative_correction_sign": derivative_sign,
                    "pulse_omitted_absolute_bound": _signed(
                        abs(scale * component.a * kernel_bound), self.precision
                    ),
                    "pulse_omitted_derivative_absolute_bound": _signed(
                        abs(scale * derivative_scale * kernel_bound), self.precision
                    ),
                    "pulse_omitted_bound_condition": (
                        "Conditional on fixed nominal coefficients and fixed "
                        "quadrature atoms; incoming/coefficient/quadrature "
                        "uncertainty is not enclosed."
                    ),
                }
            else:
                omitted = {
                    "pulse_omitted_correction_sign": None,
                    "pulse_omitted_derivative_correction_sign": None,
                    "pulse_omitted_absolute_bound": None,
                    "pulse_omitted_derivative_absolute_bound": None,
                    "pulse_omitted_bound_condition": None,
                }
            return {
                "value": value,
                "derivative": derivative,
                "N": N,
                "N_Z": N_Z,
                "pulse_atom": atom,
                "method": "continuous_mean_analytic_Z_jet",
                "cumulative_representation": (
                    "retained_terminal_balance_minus_direct_end_tail"
                    if _mp(end)>=mp.mpf('-3.15') else "forward_continuous_primitive"
                ),
                "terminal_materialized_residual_retained": True,
                "region": region,
                "terminal_mean_forced_zero": False,
                "angular_Z_jet_complete": False,
                **omitted,
            }

    def axial_average_receipt(self, logR, Z):
        jet = self.axial_average_jet(logR, Z)
        result = dict(
            nominal=_signed(jet["value"], self.precision),
            nominal_Z=_signed(jet["derivative"], self.precision),
            incoming_uncertainty_enclosed=False,
            quadrature_enclosure_certified=False,
            terminal_mean_forced_zero=False,
            continuous_mean_Z_jet_installed=True,
            angular_Z_jet_complete=False,
        )
        for key in (
            "method",
            "cumulative_representation",
            "terminal_materialized_residual_retained",
            "region",
            "N",
            "N_Z",
            "pulse_atom",
            "pulse_omitted_correction_sign",
            "pulse_omitted_derivative_correction_sign",
            "pulse_omitted_absolute_bound",
            "pulse_omitted_derivative_absolute_bound",
            "pulse_omitted_bound_condition",
        ):
            if key in jet:
                result[key] = (
                    _signed(jet[key], self.precision)
                    if key
                    in (
                        "N",
                        "N_Z",
                        "pulse_omitted_absolute_bound",
                        "pulse_omitted_derivative_absolute_bound",
                    )
                    else jet[key]
                )
        return result

    def axial_average(self, logR, Z):
        return self.axial_average_jet(logR, Z)["value"]


class ContinuousSeededOuterFieldJets(ContinuousSeededOuterField):
    """Joined field using the continuous axial value and mean jets."""

    def __init__(
        self,
        source,
        *,
        prepared=None,
        seeded_receipt=None,
        continuous_receipt=None,
        jet_precision=200,
    ):
        super().__init__(source, prepared=prepared)
        self.outer = ContinuousSeededAxialProfileJets(
            source,
            prepared=prepared,
            seeded_receipt=seeded_receipt,
            continuous_receipt=continuous_receipt,
            jet_precision=jet_precision,
        )

    def _outer_values(self, log_radius, z):
        with mp.workdps(self.work_precision):
            row = self.outer.values_with_jets(log_radius, z)
            base = self.schedule.at_log_radius(
                mp.nstr(log_radius, self.precision), z
            )
            Utheta = _mp(row["Utheta"])
            Utheta_Z = _mp(row.get("Utheta_Z", Utheta * _mp(base["dlogU_dZ"])))
            Uz = _mp(row["Uz"])
            Uz_Z = _mp(row["Uz_Z"])
            R = mp.exp(log_radius)
            root = mp.sqrt(2 * R)
            return {
                "Utheta": Utheta,
                "Utheta_Z": Utheta_Z,
                "Uz": Uz,
                "Uz_Z": Uz_Z,
                "R": R,
                "F": Utheta / root,
                "FZ": Utheta_Z / root,
                "logF_slope": _mp(row.get("logF_slope",base["logF_slope"])),
                "logarithmic_slope": _mp(row.get("logarithmic_slope",base["logarithmic_slope"])),
                "heat_method": row.get("heat_method", base.get("heat_method")),
                "axial_Z_method": row.get("axial_Z_method"),
                "angular_Z_jet_complete": False,
            }

    def _outer_average(self, log_radius, z):
        with mp.workdps(self.work_precision):
            start = self.outer.offset(
                mp.nstr(log_radius, self.precision), self.schedule.logR_p
            )
            if start < 0:
                return super()._outer_average(log_radius, z)
            jet = self.outer.axial_average_jet(
                mp.nstr(log_radius, self.precision), z
            )
            return _mp(jet["value"]), _mp(jet["derivative"])

    def _outer_state(self, *args, **kwargs):
        result = super()._outer_state(*args, **kwargs)
        result.update(
            region="continuous_seeded_axial_outer_jets",
            continuous_axial_values_and_means_installed=True,
            continuous_mean_Z_jet_installed=True,
            outer_velocity_jets_complete=False,
            velocity_jet_scope=(
                "Analytic continuous axial pulse/end-bump Uz_Z and mean_Z; "
                "angular correction Z jet remains unavailable."
            ),
            continuous_mean_closure_certified=False,
        )
        return result


def _receipt_row(profile, field, name, logR, z):
    values = profile.values_with_jets(mp.nstr(logR, field.precision), z)
    jet = profile.axial_average_jet(mp.nstr(logR, field.precision), z)
    state = field._outer_state(logR, z)
    return {
        "name": name,
        "logR": mp.nstr(logR, field.precision),
        "axial_Z_method": values["axial_Z_method"],
        "Uz": _signed(values["Uz"], field.precision),
        "Uz_Z": _signed(values["Uz_Z"], field.precision),
        "axial_average": _signed(jet["value"], field.precision),
        "axial_average_Z": _signed(jet["derivative"], field.precision),
        "Ur": _signed(state["Ur"], field.precision),
        "mass": _signed(state["mass"], field.precision),
        "mass_Z": _signed(state["mass_Z"], field.precision),
        "pulse_omitted_absolute_bound": jet.get(
            "pulse_omitted_absolute_bound"
        ),
        "pulse_omitted_derivative_absolute_bound": jet.get(
            "pulse_omitted_derivative_absolute_bound"
        ),
        "pulse_omitted_correction_sign": jet.get(
            "pulse_omitted_correction_sign"
        ),
        "pulse_omitted_derivative_correction_sign": jet.get(
            "pulse_omitted_derivative_correction_sign"
        ),
        "terminal_mean_forced_zero": False,
        "angular_Z_jet_complete": False,
    }


def run():
    print("building retained joined candidate", flush=True)
    source = build_joined_field()
    folder = Path(__file__).parent
    seeded = json.loads(
        (folder / "lei_ren_part1_paper_seeded_shared_candidate.json").read_text(
            encoding="utf-8"
        )
    )
    continuous = json.loads(
        (folder / "lei_ren_part1_paper_continuous_axial_solve.json").read_text(
            encoding="utf-8"
        )
    )
    field = ContinuousSeededOuterFieldJets(
        source,
        prepared={0.3: (seeded, continuous)},
        seeded_receipt=seeded,
        continuous_receipt=continuous,
        jet_precision=int(continuous.get("algebra_precision", 200)),
    )
    profile = field.outer
    with mp.workdps(field.precision):
        z = mp.mpf(".3")
        schedule = field.schedule
        mu = _mp(str(schedule.mu))
        log_rp = _mp(str(schedule.logR_p))
        log_rv = _mp(str(schedule.logR_v))
        points = [
            ("Rp", log_rp),
            ("pulse", log_rp + mp.mpf(5) / mu),
            ("first_end_bump", log_rv - 3),
            ("Rv", log_rv),
            ("after_Rv", log_rv + 1),
        ]
        rows = []
        for name, logR in points:
            print("sampling " + name, flush=True)
            rows.append(_receipt_row(profile, field, name, logR, z))

        pulse_logR = dict(points)["pulse"]
        center = field._outer_state(pulse_logR, z)
        divergence = []
        for h in (mp.mpf("1e-5"), mp.mpf("5e-6")):
            states = {
                i: field._outer_state(pulse_logR + i * h, z)
                for i in (-2, -1, 1, 2)
            }
            radial = (
                states[-2]["Ur"]
                - 8 * states[-1]["Ur"]
                + 8 * states[1]["Ur"]
                - states[2]["Ur"]
            ) / (12 * h)
            axial_y = (
                states[-2]["Uz"]
                - 8 * states[-1]["Uz"]
                + 8 * states[1]["Uz"]
                - states[2]["Uz"]
            ) / (12 * h)
            # This is the analytic axial jet from the pulse coefficient and
            # not a neighboring-Z coefficient solve.
            axial_Z = center["UZ"]
            R = center["R"]
            root = mp.sqrt(2 * R)
            delta = field.delta
            terms = [
                root * radial / R,
                center["Ur"] / root,
                (
                    (1 - z * z) * axial_Z
                    - 2 * z * axial_y
                    - (1 + delta) * z * center["Uz"]
                )
                / (1 - delta * z * z),
            ]
            scale = sum(abs(term) for term in terms)
            divergence.append(
                {
                    "log_radius_step": mp.nstr(h, 20),
                    "mapped_q_divergence": _signed(sum(terms), field.precision),
                    "relative_cancellation": mp.nstr(
                        abs(sum(terms)) / scale if scale else mp.mpf(0), 50
                    ),
                    "axial_Z_source": "analytic continuous pulse jet at fixed Z",
                    "neighboring_Z_coefficient_solves_used": False,
                }
            )

        tangent = profile.coefficient_tangent(z)
        input_jet = tangent["input"]
        tangent_receipt = {
            "a_Z": _signed(tangent["a_Z"], field.outer.jet_precision),
            "c_Z": _signed_logs(tangent["c_Z"], field.outer.jet_precision),
            "algebra_precision": tangent["algebra_precision"],
            "linear_tangent_relative_replay": tangent[
                "linear_tangent_relative_replay"
            ],
            "energy_tangent_relative_replay": tangent[
                "energy_tangent_relative_replay"
            ],
        }
        report = {
            "Z": ".3",
            "source_precision": field.precision,
            "seeded_receipt_incoming_precision": seeded["incoming"].get(
                "precision"
            ),
            "continuous_algebra_precision": continuous.get("algebra_precision"),
            "jet_precision": field.outer.jet_precision,
            "files": {
                "python": str(Path(__file__).resolve()),
                "json": str(Path(__file__).with_suffix(".json").resolve()),
                "markdown": str(Path(__file__).with_suffix(".md").resolve()),
            },
            "schedule_identity_preserved": profile.schedule is source.outer.schedule,
            "angular_identity_preserved": profile.angular is source.outer.angular,
            "tail_identity_preserved": profile.tail is source.outer.tail,
            "pressure_identity_preserved": profile.pressure is source.outer.pressure,
            "actual_input_tangent": input_jet,
            "continuous_coefficient_tangent": tangent_receipt,
            "rows": rows,
            "pulse_radial_divergence": divergence,
            "continuous_axial_values_and_means_installed": True,
            "continuous_mean_Z_jet_installed": True,
            "angular_Z_jet_complete": False,
            "future_angular_derivative_is_finite_difference": True,
            "global_analytic_Z_jet_claim": False,
            "full_five_moments_available": False,
            "moments_were_fabricated": False,
            "terminal_mean_forced_zero": False,
            "finite_energy_certified": False,
            "scale_recursion_certified": False,
            "scope": (
                "Bounded actual seeded continuous axial value/mean Z-jet adapter "
                "at Z=.3. Radial pulse divergence uses the analytic axial jet, "
                "while angular future energy remains float-backed finite difference."
            ),
        }
    path = Path(__file__).with_suffix(".json")
    path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "json": str(path),
                "rows": len(rows),
                "divergence": [
                    row["mapped_q_divergence"] for row in divergence
                ],
                "angular_Z_jet_complete": report["angular_Z_jet_complete"],
            }
        ),
        flush=True,
    )
    return report


if __name__ == "__main__":
    run()
