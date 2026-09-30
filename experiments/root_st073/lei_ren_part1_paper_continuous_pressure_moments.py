"""Bounded cumulative pressure moment on the shared continuous source.

The fifth Lei--Ren moment is evaluated from the actual swirl field,

    Mp(R,Z) = integral_0^R Utheta(R',Z)^2 / (2 R') dR'.

The source inner value at ``Rh`` is the only anchor.  Beyond ``Rh`` the
reference angular schedule, the signed compact angular bump, and the heat
collar/exterior are integrated as separate increments.  The same increments
are added to the actual inner pressure ``P``; a pressure datum is therefore
never silently reset to zero at the join.

Variable angular stages use the shared Gauss nodes from the angular provider
(96 by default).  Constant-slope stages use exact exponential atoms.  Heat
increments use the installed small-x heat polynomial and cap Gauss
quadrature at the finite collar ``0 <= t <= 3``; the exterior polynomial is
integrated analytically.  These are numerical providers, not certified
quadrature or finite-energy enclosures.
"""

from functools import lru_cache
import json
from pathlib import Path

import mpmath as mp
from numpy.polynomial.legendre import leggauss

from lei_ren_part1_paper_axial_correction import signed_log
from lei_ren_part1_paper_continuous_axial_pulse import ContinuousAxialPulse
from lei_ren_part1_paper_joined_outer import _mp


def _signed_value(value, precision):
    """Decode the full-value or older signed-log correction receipts."""

    if isinstance(value, mp.mpf):
        return value
    if not isinstance(value, dict):
        return _mp(value)
    exact = value.get("arbitrary_exponent_value")
    if exact is not None:
        return _mp(exact)
    sign = int(value.get("sign", 0))
    logarithm = value.get("log_abs")
    if sign == 0 or logarithm is None:
        return mp.mpf(0)
    with mp.workdps(max(32, int(precision))):
        return (1 if sign > 0 else -1) * mp.exp(_mp(logarithm))


class ContinuousPressureMoments:
    """Actual cumulative ``Mp`` and pressure on a shared outer profile.

    ``profile`` is the live ``ContinuousIncomingProfile`` used by the
    continuous incoming field.  The class intentionally does not patch that
    profile; the parent integration lane can install this provider on its
    own API when the interface is selected.
    """

    def __init__(self, profile, *, order=96):
        self.profile = profile
        self.schedule = profile.schedule
        self.precision = int(profile.precision)
        self.order = int(order)
        if self.order < 16:
            raise ValueError("At least 16 quadrature nodes required")
        self.angular = getattr(profile, "angular_moment_provider", None)
        if self.angular is None:
            raise TypeError("profile must expose angular_moment_provider")
        self.correction = profile.angular_correction_provider
        self.heat_provider = getattr(self.angular, "heat_provider", None)
        if self.heat_provider is None:
            raise TypeError("angular provider must expose its heat provider")
        with mp.workdps(self.precision):
            # Parse every tiny parameter and local coordinate while the
            # declared precision is active.  Parsing an MP value outside the
            # context would round the ``Rh`` offset before Decimal sees it.
            self.mu = _mp(str(self.schedule.mu))
            self.delta = _mp(str(self.schedule.delta))
            self._rh_y = _mp(str(profile.offset(
                mp.nstr(profile.seed_source.logRh, self.precision),
                self.schedule.logRref)))
            self._tail_y = _mp(str(self.schedule.y_tail))
            nodes, weights = leggauss(self.order)
            self.nodes = [
                ((_mp(str(float(node))) + 1) / 2,
                 _mp(str(float(weight))) / 2)
                for node, weight in zip(nodes, weights)
            ]
        self.stages = list(self.angular.stages)

    @staticmethod
    def atom(rate, left, right):
        """Return ``integral_left^right exp(rate*t) dt`` stably."""

        if right <= left:
            return mp.mpf(0)
        rate = _mp(rate)
        left = _mp(left)
        right = _mp(right)
        if rate:
            return mp.exp(rate * left) * mp.expm1(rate * (right - left)) / rate
        return right - left

    def _target_y(self, log_radius):
        # ``str(mp.mpf('5e152'))`` is intentionally short even inside a
        # high-precision work context.  Convert through ``nstr`` first so a
        # huge absolute source radius does not lose its local offset.
        if isinstance(log_radius, mp.mpf):
            log_radius = mp.nstr(log_radius, self.precision)
        return _mp(str(self.profile.offset(log_radius, self.schedule.logRref)))

    def _base(self, y):
        """Use the shared angular normalization at an exact local offset."""

        return _mp(self.angular._base(mp.nstr(_mp(y), self.precision)))

    def _ratio_jet(self, y, z):
        """Use the shared angular Z ratio and its analytic derivative."""

        return self.angular._ratio_jet(_mp(y), _mp(z))

    def _transfer(self, state, left, end, z, slope):
        """Advance the normalized pressure state over one angular stage.

        The normalized state is ``q = DeltaMp/(U_base(y)^2/2)``.  Thus
        ``q' = G^2 - 2 E' q`` and ``q_Z' = 2 G G_Z - 2 E' q_Z``.
        """

        left = _mp(left)
        end = _mp(end)
        length = end - left
        if length <= 0:
            return list(state)
        with mp.workdps(self.precision):
            if slope is not None:
                if slope == "pulse":
                    rate_s = -mp.mpf(".5") - self.mu
                elif slope == "waiting":
                    rate_s = -(1 + self.delta) / 2
                else:
                    rate_s = _mp(slope)
                rate = 2 * rate_s
                damp = mp.exp(-rate * length)
                source = self._ratio_jet(left, z)
                atom = -mp.expm1(-rate * length) / rate if rate else length
                return [
                    state[0] * damp + source[0] ** 2 * atom,
                    # _ratio_jet returns the logarithmic Z jet LZ, not
                    # dG/dZ.  Since the source is G^2, its derivative is
                    # 2 * LZ * G^2.
                    state[1] * damp + 2 * source[1] * source[0] ** 2 * atom,
                ]

            # Variable stages are deliberately finite.  The very long
            # power/flattening buffers are represented by constant stages in
            # the shared schedule and never reach this Gauss loop.
            Eleft = self._base(left)
            Eend = self._base(end)
            damp = mp.exp(-2 * (Eend - Eleft))
            value = state[0] * damp
            derivative = state[1] * damp
            for unit, weight in self.nodes:
                local = length * unit
                y = left + local
                Et = self._base(y)
                G, LZ = self._ratio_jet(y, z)
                factor = length * weight * mp.exp(-2 * (Eend - Et))
                value += factor * G ** 2
                derivative += factor * 2 * LZ * G ** 2
            return [value, derivative]

    def _preheat_normalized(self, y, z):
        """Integrate the reference pressure increment from ``Rh`` to ``y``."""

        with mp.workdps(self.precision):
            target = min(_mp(y), self._tail_y)
            state = [mp.mpf(0), mp.mpf(0)]
            if target <= self._rh_y:
                return state
            zero = mp.mpf(0)
            if self._rh_y < min(target, zero):
                state = self._transfer(state, self._rh_y, min(target, zero), z, "0.1")
            if target <= zero:
                return state
            bounds = self.schedule._make_stage_bounds()
            for name, slope in self.stages:
                left, right = bounds[name]
                left, right = _mp(str(left)), _mp(str(right))
                if target <= left:
                    break
                end = min(target, right)
                state = self._transfer(state, left, end, z, slope)
                if end == target:
                    break
            return state

    def _reference_increment(self, log_radius, z):
        """Return reference ``(Mp, Mp_Z)`` through the pre-heat schedule."""

        with mp.workdps(self.precision):
            y = self._target_y(log_radius)
            if y > self._tail_y:
                y = self._tail_y
            state = self._preheat_normalized(y, z)
            scale = mp.mpf(".5") * mp.exp(2 * self._base(y))
            return scale * state[0], scale * state[1]

    def _bump_increment(self, log_radius, z):
        """Return physical signed bump pressure increment and Z jet."""

        with mp.workdps(self.precision):
            if isinstance(log_radius, mp.mpf):
                log_radius = mp.nstr(log_radius, self.precision)
            t = _mp(str(self.profile.offset(log_radius, self.schedule.logR_rel)))
            if t <= mp.mpf("-3.15"):
                return {
                    "value": mp.mpf(0),
                    "derivative": mp.mpf(0),
                    "terms": [],
                    "normalized": [mp.mpf(0), mp.mpf(0)],
                }
            lam = -(1 + 2 * self.mu)
            # Use the same exact-Z spelling as the installed continuous
            # coefficient cache.  Passing an MP value that was constructed
            # at a lower ambient precision can otherwise collide with a
            # receipt carrying the full 443-digit key.
            coefficients = self.correction.coefficients(mp.nstr(z, self.precision))
            scale = mp.exp(2 * self._base(_mp(str(self.schedule.y_rel))))
            value = mp.mpf(0)
            derivative = mp.mpf(0)
            terms = []
            for key, center in (("d1", -3), ("d2", -1)):
                center = mp.mpf(center)
                d = _signed_value(coefficients[key], self.precision)
                dz = _signed_value(coefficients[key + "_Z"], self.precision)
                linear = mp.exp(lam * center) * self.correction.weighted_atom(
                    lam, power=1, upper=t - center)
                square = mp.exp(lam * center) * self.correction.weighted_atom(
                    lam, power=2, upper=t - center)
                term = d * linear + d * d * square / 2
                term_z = dz * linear + d * dz * square
                value += term
                derivative += term_z
                terms.append({
                    "center": int(center),
                    "coefficient": signed_log(d, self.precision),
                    "coefficient_Z": signed_log(dz, self.precision),
                    "linear": signed_log(linear, self.precision),
                    "square": signed_log(square, self.precision),
                    "value": signed_log(scale * term, self.precision),
                    "derivative": signed_log(scale * term_z, self.precision),
                })
            return {
                "value": scale * value,
                "derivative": scale * derivative,
                "terms": terms,
                "normalized": [value, derivative],
                "scale": scale,
                "exponent": lam,
                "continuous_bump_atoms_shared": True,
            }

    def _heat_point(self, v, z, heat):
        """Evaluate the installed heat polynomial at one collar node."""

        x = _mp(heat["xi"])
        xz = _mp(heat["xi_Z"])
        h = self.delta / 2
        c1 = h * (1 + h)
        c2 = h * (h + 1) ** 2 * (h + 2)
        sw = ContinuousAxialPulse.sigma_pair(v)[0]
        edge = (3 - v) / 2
        phi = mp.exp(-1 / edge ** 2) if edge > 0 else mp.mpf(0)
        epsilon = _mp(self.schedule.epsilon)
        factor = 1 - epsilon * phi
        K0 = (1 - sw) * (1 - epsilon) + sw * factor
        xx = x * mp.exp(-v)
        xxz = xz * mp.exp(-v)
        dk = sw * factor * (-c1 * xx + c2 * xx * xx / 2)
        dkz = sw * factor * (-c1 + c2 * xx) * xxz
        return {
            "K0": K0,
            "K0_squared": K0 * K0,
            "K0_Z": mp.mpf(0),
            "dk": dk,
            "dk_Z": dkz,
            "reference": K0 * K0 / 2,
            "reference_Z": mp.mpf(0),
            "correction": K0 * dk + dk * dk / 2,
            "correction_Z": (K0 + dk) * dkz,
            "heat_deficit": heat.get("deficit"),
        }

    def heat_increments(self, t, z):
        """Return physical finite heat reference/correction increments."""

        with mp.workdps(self.precision):
            t = _mp(t)
            if t < 0:
                raise ValueError("Heat increments start at Rtail")
            # The shared heat provider now owns the common collar and
            # exterior polynomial atoms.  Use that public helper when it is
            # installed; the local path below remains a compatibility
            # fallback for an older checkout of this experiment.
            shared = getattr(self.heat_provider, "pressure_increments", None)
            if shared is not None:
                row = shared(t, z)
                correction = row.get("heat_correction", row.get("correction", 0))
                correction_z = row.get(
                    "heat_correction_Z", row.get("correction_Z", 0))
                reference = row.get("reference", row.get("heat_reference", 0))
                reference_z = row.get(
                    "reference_Z", row.get("heat_reference_Z", 0))
                normalized = row.get("normalized", {})
                return {
                    "reference": _mp(reference),
                    "reference_Z": _mp(reference_z),
                    "correction": _mp(correction),
                    "correction_Z": _mp(correction_z),
                    "normalized_reference": normalized.get("reference", 0)
                    if isinstance(normalized, dict) else normalized,
                    "normalized_reference_Z": normalized.get("reference_Z", 0)
                    if isinstance(normalized, dict) else 0,
                    "normalized_correction": normalized.get("heat_correction", normalized.get("correction", 0))
                    if isinstance(normalized, dict) else 0,
                    "normalized_correction_Z": normalized.get("heat_correction_Z", normalized.get("correction_Z", 0))
                    if isinstance(normalized, dict) else 0,
                    "scale": row.get("scale"),
                    "heat_deficit": row.get("heat_deficit"),
                    "heat_kernel_shared": True,
                    "collar_order": row.get("quadrature_order", self.order),
                    "collar_end": min(t, mp.mpf(3)),
                    "exterior_polynomial_analytic": bool(t > 3),
                    "shared_heat_provider": True,
                }
            heat = self.profile.angular_schedule_provider.heat_jet(
                self.schedule.logR_tail, z)
            collar_end = min(t, mp.mpf(3))
            lam = 1 + self.delta
            reference = reference_z = correction = correction_z = mp.mpf(0)
            if collar_end > 0:
                for unit, weight in self.nodes:
                    v = collar_end * unit
                    point = self._heat_point(v, z, heat)
                    factor = collar_end * weight * mp.exp(-lam * v)
                    reference += factor * point["reference"]
                    reference_z += factor * point["reference_Z"]
                    correction += factor * point["correction"]
                    correction_z += factor * point["correction_Z"]

            # After the finite collar sigma=1 and phi=0.  The same quadratic
            # heat polynomial is therefore an exact finite polynomial in
            # exp(-t), so every exterior atom is analytic.
            if t > 3:
                h = self.delta / 2
                c1 = h * (1 + h)
                c2 = h * (h + 1) ** 2 * (h + 2)
                x = _mp(heat["xi"])
                xz = _mp(heat["xi_Z"])
                b1 = -c1 * x
                b2 = c2 * x * x / 2
                b1z = -c1 * xz
                b2z = c2 * x * xz
                reference += self.atom(-lam, mp.mpf(3), t)
                correction_coefficients = (b1, b2 + b1 * b1 / 2,
                                            b1 * b2, b2 * b2 / 2)
                correction_coefficients_z = (
                    b1z, b2z + b1 * b1z,
                    b1z * b2 + b1 * b2z, b2 * b2z)
                for j, (coefficient, coefficient_z) in enumerate(
                    zip(correction_coefficients, correction_coefficients_z), 1
                ):
                    atom = self.atom(-(lam + j), mp.mpf(3), t)
                    correction += coefficient * atom
                    correction_z += coefficient_z * atom

            # E_tail is the common source amplitude before K0.  This is the
            # pressure analogue of the heat moment scale, with no extra R.
            tail_log = _mp(str(self.schedule.logR_tail))
            a = (1 + self.delta) / 2
            E_tail = _mp(self.schedule._log_c_inf) - a * tail_log
            scale = mp.exp(2 * E_tail)
            return {
                "reference": scale * reference / 1,
                "reference_Z": scale * reference_z,
                "correction": scale * correction,
                "correction_Z": scale * correction_z,
                "normalized_reference": reference,
                "normalized_reference_Z": reference_z,
                "normalized_correction": correction,
                "normalized_correction_Z": correction_z,
                "scale": scale,
                "heat_deficit": heat.get("deficit"),
                "heat_kernel_shared": True,
                "collar_order": self.order,
                "collar_end": collar_end,
                "exterior_polynomial_analytic": bool(t > 3),
            }

    def _inner_anchor(self, z):
        key = mp.nstr(_mp(z), self.precision)
        return self._inner_anchor_cached(key)

    @lru_cache(maxsize=128)
    def _inner_anchor_cached(self, z_key):
        with mp.workdps(self.precision):
            z = _mp(z_key)
            source = self.profile.seed_source
            row = source.inner.evaluate_x(mp.e, z)
            moments_z = row.get("momentsZ", row.get("moments_Z"))
            return {
                "p": _mp(row["moments"]["p"]),
                "p_Z": _mp(moments_z["p"]),
                "P": _mp(row["P"]),
                "PZ": _mp(row["PZ"]),
            }

    def terminal_pressure_jet(self, Z):
        """Actual P(infinity,Z) and target for the compact bump integral.

        Retain the baseline and heat correction separately; this supplies
        a matching target, not a pressure gauge change or a coefficient solve.
        The heat truncation bound covers the integral value only, not its Z jet.
        """
        with mp.workdps(self.precision):
            z=_mp(Z)
            if abs(z)>=1:
                raise ValueError('Actual inner seed requires |Z|<1')
            anchor=self._inner_anchor(z)
            ref,ref_z=self._reference_increment(self.schedule.logR_tail,z)
            bump=self._bump_increment(self.schedule.logR_tail,z)
            heat=self.heat_provider.complete_pressure_heat_integral(z)
            base=mp.fsum([anchor['P'],ref,heat['reference']])
            base_z=mp.fsum([anchor['PZ'],ref_z])
            correction=heat['heat_correction']
            correction_z=heat['heat_correction_Z']
            required=-mp.fsum([base,correction])
            required_z=-mp.fsum([base_z,correction_z])
            return dict(P_infinity=mp.fsum([base,bump['value'],correction]),
                P_infinity_Z=mp.fsum([base_z,bump['derivative'],correction_z]),
                pressure_baseline=base,pressure_baseline_Z=base_z,
                current_bump=bump['value'],current_bump_Z=bump['derivative'],
                required_bump=required,required_bump_Z=required_z,
                additional_bump_required=required-bump['value'],
                additional_bump_required_Z=required_z-bump['derivative'],
                inner_anchor=anchor,preheat_reference=ref,preheat_reference_Z=ref_z,
                heat=heat,components_retained_separately=True,
                pressure_datum_changed=False,angular_coefficients_changed=False,
                quadrature_error_enclosed=False,arithmetic_error_enclosed=False,
                pressure_Z_truncation_enclosed=False,
                pressure_terminal_compatibility_certified=False)

    def radial_integrand_jet(self, logR, Z):
        """Return the local radial derivatives implied by the actual swirl.

        A finite difference of a cumulative pressure value can lose a tiny
        tail increment when the already accumulated datum is many orders of
        magnitude larger.  This local query keeps the exact pointwise
        integrand available for those source-scale checks:

        ``d_logR Mp = d_logR P = Utheta^2/2`` and
        ``d_logR Mp_Z = d_logR P_Z = Utheta*Utheta_Z``.
        """

        with mp.workdps(self.precision):
            z = _mp(Z)
            values = self.profile.values_with_jets(logR, z)
            return {
                "value": mp.mpf(".5") * values["Utheta"] ** 2,
                "derivative": values["Utheta"] * values["Utheta_Z"],
                "Utheta": values["Utheta"],
                "Utheta_Z": values["Utheta_Z"],
                "finite_difference_precision_limited": True,
            }

    def moments_jet(self, logR, Z, *, include_inner=True):
        """Return ``Mp``, its Z jet, and pressure from ``Rh`` outward."""

        with mp.workdps(self.precision):
            z = _mp(Z)
            if abs(z) >= 1:
                raise ValueError("Actual inner seed requires |Z|<1")
            y = self._target_y(logR)
            if y < self._rh_y:
                raise ValueError("Continuous pressure moments start at Rh")
            anchor = self._inner_anchor(z) if include_inner else {
                "p": mp.mpf(0), "p_Z": mp.mpf(0),
                "P": mp.mpf(0), "PZ": mp.mpf(0)}
            if isinstance(logR, mp.mpf):
                logR = mp.nstr(logR, self.precision)
            if _mp(str(self.profile.offset(logR, self.schedule.logR_tail))) <= 0:
                ref, ref_z = self._reference_increment(logR, z)
                bump = self._bump_increment(logR, z)
                heat = {
                    "reference": mp.mpf(0), "reference_Z": mp.mpf(0),
                    "correction": mp.mpf(0), "correction_Z": mp.mpf(0),
                    "normalized_reference": mp.mpf(0),
                    "normalized_reference_Z": mp.mpf(0),
                    "normalized_correction": mp.mpf(0),
                    "normalized_correction_Z": mp.mpf(0),
                    "heat_kernel_shared": True,
                    "collar_order": self.order,
                    "exterior_polynomial_analytic": False,
                }
                region = "continuous_preheat"
            else:
                ref, ref_z = self._reference_increment(self.schedule.logR_tail, z)
                bump = self._bump_increment(self.schedule.logR_tail, z)
                heat = self.heat_increments(
                    _mp(str(self.profile.offset(logR, self.schedule.logR_tail))), z)
                region = "continuous_heat"
            heat_ref = heat["reference"]
            heat_ref_z = heat["reference_Z"]
            heat_corr = heat["correction"]
            heat_corr_z = heat["correction_Z"]
            total_ref = ref + heat_ref
            total_ref_z = ref_z + heat_ref_z
            total_bump = bump["value"] + heat_corr
            total_bump_z = bump["derivative"] + heat_corr_z
            increment = total_ref + total_bump
            increment_z = total_ref_z + total_bump_z
            local_integrand = self.radial_integrand_jet(logR, z)
            return {
                "p": anchor["p"] + increment,
                "p_Z": anchor["p_Z"] + increment_z,
                "pressure_moment": anchor["p"] + increment,
                "pressure_moment_Z": anchor["p_Z"] + increment_z,
                "P": anchor["P"] + increment,
                "PZ": anchor["PZ"] + increment_z,
                "inner_anchor": anchor,
                "reference_increment": total_ref,
                "reference_increment_Z": total_ref_z,
                "bump_increment": total_bump,
                "bump_increment_Z": total_bump_z,
                "heat_reference_increment": heat_ref,
                "heat_reference_increment_Z": heat_ref_z,
                "heat_correction_increment": heat_corr,
                "heat_correction_increment_Z": heat_corr_z,
                "preheat_reference_increment": ref,
                "preheat_reference_increment_Z": ref_z,
                "preheat_bump_increment": bump["value"],
                "preheat_bump_increment_Z": bump["derivative"],
                "bump_terms": bump.get("terms", []),
                "heat_detail": heat,
                "radial_integrand": local_integrand["value"],
                "radial_integrand_Z": local_integrand["derivative"],
                "region": region,
                "actual_inner_anchor_once": bool(include_inner),
                "pressure_propagated_from_inner": True,
                "full_pressure_asymptote_forced_zero": False,
                "continuous_bump_atoms_shared": True,
                "heat_kernel_shared": True,
                "quadrature_order": self.order,
                "quadrature_error_enclosed": False,
                "heat_polynomial_truncation_enclosed": False,
                "complete_five_moments": False,
                "finite_energy_certified": False,
                "stress_or_PDE_claim": False,
                "scale_recursion_certified": False,
            }

    # Short aliases used by parent integration code and older moment APIs.
    moments = moments_jet
    pressure_jet = moments_jet


ContinuousPressureProvider = ContinuousPressureMoments


def _build_prepared_field():
    from lei_ren_part1_paper_continuous_incoming_outer import ContinuousIncomingOuterField
    from lei_ren_part1_paper_joined_outer import build_joined_field

    source = build_joined_field()
    folder = Path(__file__).parent
    seeded = json.loads((folder / "lei_ren_part1_paper_seeded_shared_candidate.json").read_text())
    atoms = json.loads((folder / "lei_ren_part1_paper_continuous_axial_solve.json").read_text())
    field = ContinuousIncomingOuterField(source, prepared={.3: (seeded, atoms)})
    return source, field


def run():
    """Run radial and separate-correction checks for prepared ``Z=.3``."""

    print("building shared field for cumulative pressure moments", flush=True)
    source, field = _build_prepared_field()
    profile = field.outer
    engine = ContinuousPressureMoments(profile, order=96)
    with mp.workdps(field.precision):
        # Construct the prepared coordinate only after entering the declared
        # source precision.  A default-dps ``mpf('.3')`` carries a rounded
        # binary-looking decimal and cannot match the exact coefficient-cache
        # receipt installed for the prepared seed.
        z = mp.mpf(".3")
        schedule = source.schedule
        rh = source.logRh
        points = [
            ("Rh", rh),
            ("Rref", profile.log_at(schedule.logRref, "0")),
            ("Rp", profile.log_at(schedule.logR_p, "0")),
            ("flatten", profile.log_at(schedule.logR_v, "50")),
            ("bump", profile.log_at(schedule.logR_rel, "-2.96")),
            ("heat_collar", profile.log_at(schedule.logR_tail, ".6")),
            ("heat_exterior", profile.log_at(schedule.logR_tail, "4")),
        ]
        samples = []
        radial_checks = []
        for label, radius in points:
            row = engine.moments_jet(radius, z)
            if label == "Rh":
                anchor = source.inner.evaluate_x(mp.e, z)
                errors = [
                    abs(row["p"] / anchor["moments"]["p"] - 1),
                    abs(row["p_Z"] / anchor["momentsZ"]["p"] - 1),
                    abs(row["P"] / anchor["P"] - 1),
                    abs(row["PZ"] / anchor["PZ"] - 1),
                ]
                if max(errors) > mp.mpf("1e-70"):
                    raise ArithmeticError("Inner pressure/Mp anchor was reapplied incorrectly")
            samples.append({
                "label": label,
                "region": row["region"],
                "p": signed_log(row["p"], field.precision),
                "p_Z": signed_log(row["p_Z"], field.precision),
                "P": signed_log(row["P"], field.precision),
                "PZ": signed_log(row["PZ"], field.precision),
                "reference_increment": signed_log(row["reference_increment"], field.precision),
                "bump_increment": signed_log(row["bump_increment"], field.precision),
                "heat_correction_increment": signed_log(row["heat_correction_increment"], field.precision),
            })

            if label == "Rh":
                continue
            step = mp.mpf("1e-4")
            local = engine.radial_integrand_jet(radius, z)
            targets = [local["value"], local["derivative"]]
            # At Rp and in the long decaying buffers the cumulative value is
            # many orders larger than its local tail increment.  Retain a
            # pointwise analytic check there and run an independent centered
            # finite difference whenever the declared precision can resolve
            # the increment from the accumulated datum.
            errors = [
                abs(row["radial_integrand"] / targets[0] - 1),
                abs(row["radial_integrand_Z"] / targets[1] - 1),
            ]
            finite_difference_errors = None
            if targets[0] and mp.log10(abs(targets[0])) > -mp.mpf("200"):
                neighbors = {
                    shift: engine.moments_jet(
                        profile.log_at(radius, shift * step), z, include_inner=False
                    )
                    for shift in (-2, -1, 1, 2)
                }
                finite_difference_errors = []
                for key, target in (("p", targets[0]), ("p_Z", targets[1]),
                                    ("P", targets[0]), ("PZ", targets[1])):
                    derivative = sum(
                        weight * neighbors[shift][key]
                        for shift, weight in ((-2, 1), (-1, -8), (1, 8), (2, -1))
                    ) / (12 * step)
                    finite_difference_errors.append(abs(derivative / target - 1))
                if max(finite_difference_errors) > mp.mpf("1e-8"):
                    raise ArithmeticError("Pressure/Mp radial derivative differs from actual swirl")
            radial_checks.append({
                "label": label,
                "pointwise_relative_errors": [mp.nstr(error, 40) for error in errors],
                "finite_difference_relative_errors": (
                    None if finite_difference_errors is None
                    else [mp.nstr(error, 40) for error in finite_difference_errors]
                ),
                "finite_difference_precision_limited": finite_difference_errors is None,
            })

        # The signed bump is deliberately checked on its own scale, where a
        # large reference sum cannot hide the tiny compact correction.
        bump_radius = profile.log_at(schedule.logR_rel, "-2.96")
        step = mp.mpf("1e-5")
        bump_neighbors = {
            shift: engine._bump_increment(profile.log_at(bump_radius, shift * step), z)
            for shift in (-2, -1, 1, 2)
        }
        bump_value = engine._bump_increment(bump_radius, z)
        t = mp.mpf("-2.96")
        point = engine.correction.value_jet(t, mp.nstr(z, engine.precision))
        scale = mp.exp(2 * engine._base(_mp(str(schedule.y_rel))))
        lam = -(1 + 2 * engine.mu)
        bump_targets = [
            scale * mp.exp(lam * t) * (point["h"] + point["h"] ** 2 / 2),
            scale * mp.exp(lam * t) * (point["h_Z"] + point["h"] * point["h_Z"]),
        ]
        bump_errors = []
        for key, target in (("value", bump_targets[0]), ("derivative", bump_targets[1])):
            derivative = sum(
                weight * bump_neighbors[shift][key]
                for shift, weight in ((-2, 1), (-1, -8), (1, 8), (2, -1))
            ) / (12 * step)
            bump_errors.append(abs(derivative / target - 1))
        if max(bump_errors) > mp.mpf("1e-8"):
            raise ArithmeticError("Separate signed bump pressure primitive failed")
        if not any(term["value"]["sign"] for term in bump_value["terms"]):
            raise ArithmeticError("Separate signed bump correction was rounded away")

        # Heat corrections are also checked apart from the large reference
        # pressure.  Both collar and analytic exterior are exercised.
        heat_checks = []
        for t in (mp.mpf(".6"), mp.mpf("4")):
            step = mp.mpf("1e-4")
            neighbors = {
                shift: engine.heat_increments(t + shift * step, z)
                for shift in (-2, -1, 1, 2)
            }
            heat = engine.heat_increments(t, z)
            hrow = profile.angular_schedule_provider.heat_jet(schedule.logR_tail, z)
            p0 = engine._heat_point(t, z, hrow)
            E = _mp(schedule._log_c_inf) - (1 + engine.delta) * _mp(str(schedule.logR_tail)) / 2
            target = mp.exp(2 * E) * mp.exp(-(1 + engine.delta) * t)
            targets = [target * p0["correction"], target * p0["correction_Z"]]
            errors = []
            for key, expected in (("correction", targets[0]), ("correction_Z", targets[1])):
                derivative = sum(
                    weight * neighbors[shift][key]
                    for shift, weight in ((-2, 1), (-1, -8), (1, 8), (2, -1))
                ) / (12 * step)
                errors.append(abs(derivative / expected - 1))
            if max(errors) > mp.mpf("1e-8"):
                raise ArithmeticError("Separate heat pressure correction primitive failed")
            heat_checks.append({
                "t": mp.nstr(t, 20),
                "relative_errors": [mp.nstr(error, 40) for error in errors],
                "correction_retained": heat["correction"] != 0,
                "exterior_polynomial_analytic": heat["exterior_polynomial_analytic"],
            })

        report = {
            "samples": samples,
            "radial_integrand_checks": radial_checks,
            "separate_bump_integrand_checks": [mp.nstr(value, 40) for value in bump_errors],
            "separate_heat_integrand_checks": heat_checks,
            "prepared_Z": ".3",
            "inner_pressure_and_Mp_anchor_once": True,
            "actual_Mp_definition": "integral Utheta^2/(2R) dR",
            "pressure_propagation": "inner P/PZ plus the same reference, bump, and heat increments",
            "reference_and_bump_increments_separate": True,
            "heat_reference_and_correction_separate": True,
            "gauss_order_variable_and_heat_collar": 96,
            "constant_stage_exponentials_exact": True,
            "huge_span_variable_gauss_avoided": True,
            "quadrature_error_enclosed": False,
            "heat_polynomial_truncation_enclosed": False,
            "complete_five_moments": False,
            "finite_energy_certified": False,
            "stress_or_PDE_claim": False,
            "scale_recursion_certified": False,
        }
    Path(__file__).with_suffix(".json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "radial_checks": len(radial_checks),
        "bump_errors": report["separate_bump_integrand_checks"],
        "heat_checks": len(heat_checks),
    }), flush=True)
    return report


if __name__ == "__main__":
    run()
