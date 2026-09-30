"""Single-owner live continuous axial atoms and coefficients.

The saved solve receipt remains useful provenance, but this runtime owns one
live continuous bump basis, one live pulse provider, one complete pulse atom
per row, and the numeric coefficient solve derived from those atoms.  The
component adapter shares those objects by identity and keeps the evaluated
terminal residual instead of replacing it by the exact functional identity.
"""

import json
from functools import lru_cache
from pathlib import Path
import types

import mpmath as mp

from lei_ren_part1_paper_axial_correction import from_signed_log, signed_log
from lei_ren_part1_paper_continuous_axial_basis import ContinuousAxialBump
from lei_ren_part1_paper_continuous_axial_pulse import ContinuousAxialPulse
from lei_ren_part1_paper_continuous_axial_solve import ContinuousAxialCorrection
from lei_ren_part1_paper_continuous_axial_tangents import ContinuousAxialAlgebra


HERE = Path(__file__).parent


def _unwrap_receipt(receipt):
    if "continuous_solve" in receipt:
        return receipt["continuous_solve"]
    return receipt


def _decode(value):
    if isinstance(value,mp.mpf):
        return value
    if isinstance(value, dict) and "arbitrary_exponent_value" in value:
        return from_signed_log(value)
    return mp.mpf(value)


def _signed(value, precision):
    value = mp.mpf(value)
    if value == 0:
        return {"sign": 0, "log_abs": None, "arbitrary_exponent_value": "0"}
    return signed_log(value, precision)


def _json_safe(value, precision):
    if isinstance(value, mp.mpf):
        return mp.nstr(value, precision)
    if isinstance(value, dict):
        return {str(key): _json_safe(item, precision) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item, precision) for item in value]
    return value


def _relative(actual, expected):
    actual = mp.mpf(actual)
    expected = mp.mpf(expected)
    if expected == 0:
        return abs(actual)
    return abs(actual - expected) / abs(expected)


class _LiveAxialAlgebra(ContinuousAxialAlgebra):
    """ContinuousAxialAlgebra methods over already evaluated MP atoms."""

    def __init__(self, *, precision, mu, matrix, determinant, pulse, base,
                 target, K, Kp):
        self.precision = precision
        self.mu = mu
        self.M = matrix
        self.det = determinant
        self.p = pulse
        self.base = base
        self.target = target
        self.K = K
        self.Kp = Kp
        self.v = self.affine(self.p)
        self._cached_solution = None

    def solve(self, base=None, target=None):
        if base is None and target is None and self._cached_solution is not None:
            return self._cached_solution
        return super().solve(base=base, target=target)

    def retain_solution(self, amplitude, coefficients):
        self._cached_solution = (amplitude, coefficients)
        self.a = amplitude
        self.c = coefficients


class _RuntimeAxialCorrection(ContinuousAxialCorrection):
    """ContinuousAxialCorrection adapter with shared live atom ownership."""

    def __init__(self, runtime):
        self.runtime = runtime
        self.precision = runtime.precision
        self.basis = runtime.basis
        self.pulse = runtime.pulse
        self.mu = runtime.mu
        self.c = runtime.c
        self.a = runtime.a
        self.base = runtime.base

    def weighted_primitive(self,row,offset):
        """Use the stored full atom whenever a bump is completely passed."""
        with mp.workdps(self.precision):
            if row not in (1,2):raise ValueError('Expected row 1 or 2')
            end=mp.mpf(offset);lam=mp.mpf('.5')-row*self.mu
            atoms=[self.runtime.matrix[row-1][index]
                   if end-center>=self.basis.ell else
                   mp.exp(lam*center)*self.basis.primitive(lam,end-center)
                   for index,center in enumerate((-3,-1))]
            return sum(c*atom for c,atom in zip(self.c,atoms))

    def weighted_tail(self,row,offset,*,coefficients=None):
        """Use the same stored complete atom before each bump's support."""
        with mp.workdps(self.precision):
            if row not in (1,2):raise ValueError('Expected row 1 or 2')
            end=mp.mpf(offset);lam=mp.mpf('.5')-row*self.mu
            c=self.c if coefficients is None else coefficients
            atoms=[self.runtime.matrix[row-1][index]
                   if end-center<=-self.basis.ell else
                   mp.exp(lam*center)*self.basis.tail(lam,end-center)
                   for index,center in enumerate((-3,-1))]
            return sum(value*atom for value,atom in zip(c,atoms))

    @lru_cache(maxsize=8)
    def terminal_balance(self, row, evaluation_precision=None):
        """Evaluate the stored complete atoms without reconstructing JSON."""

        with mp.workdps(max(self.precision, evaluation_precision or self.precision)):
            if row not in (1, 2):
                raise ValueError("Expected row 1 or 2")
            pulse_integral = self.runtime.p[row - 1]
            end = self.runtime.full_end(self.c, row)
            value = self.base[row - 1] + self.a * pulse_integral + end
            atom = self.runtime.pulse_rows[row - 1]
            return {
                "value": value,
                "pulse_integral": pulse_integral,
                "pulse_atom": atom,
                "full_end_primitive": end,
                "runtime_full_end_used": True,
                "materialized_residual_retained": True,
                "exact_functional_identity": (
                    "M c + b + a p = 0 for exact continuous atoms"
                ),
                "functional_identity_certified_for_materialized_coefficients": False,
                "incoming_uncertainty_enclosed": False,
                "quadrature_enclosure_certified": False,
                "terminal_mean_forced_zero": False,
            }


class SharedContinuousAxialRuntime:
    """Own live continuous atoms, the solve, and one compatible component."""

    def __init__(self, seeded, receipt):
        self.seeded = seeded
        source_receipt = receipt
        solve = _unwrap_receipt(receipt)
        self.source_receipt = source_receipt
        self.solve_receipt = solve
        self.precision = int(solve["algebra_precision"])

        with mp.workdps(self.precision):
            self.mu = mp.mpf(solve["input_mu"])
            self.basis = ContinuousAxialBump(precision=self.precision)
            self.pulse = ContinuousAxialPulse(precision=self.precision)

            # Patch only this one pulse object with a complete-row cache.  The
            # original provider remains the implementation for partial rows
            # and value jets, and each cached row keeps its omitted bounds.
            self._pulse_full_cache = {}
            raw_full_row = self.pulse.full_row

            def cached_full_row(pulse_self, mu, row, *, band=48):
                with mp.workdps(self.precision):
                    key = (
                        mp.nstr(_decode(mu), self.precision),
                        int(row),
                        mp.nstr(_decode(band), self.precision),
                    )
                    if key not in self._pulse_full_cache:
                        self._pulse_full_cache[key] = raw_full_row(mu, row, band=band)
                    return self._pulse_full_cache[key]

            self.pulse.full_row = types.MethodType(cached_full_row, self.pulse)

            self.matrix, self.determinant = self.basis.matrix(self.mu)
            self.pulse_rows = [
                self.pulse.full_row(self.mu, row) for row in (1, 2)
            ]
            # Each p_i is evaluated exactly once from its live row atom and is
            # then shared with the algebra and component by list identity.
            self.p = [
                mp.exp(mp.mpf(row["log_normalized_pulse_integral"]))
                for row in self.pulse_rows
            ]

            self.base = self._base_rows(seeded, source_receipt)
            self.target = self._target(seeded, solve)
            self.energy_atom_precision = int(
                solve.get("energy_atom_precision", 100)
            )
            self.Kp = mp.mpf(solve["K_p"])
            gram = self.basis.energy_gram(self.mu)
            self.K = [
                mp.exp(-26 + index * self.mu) * gram for index in (6, 2)
            ]

            self.algebra = _LiveAxialAlgebra(
                precision=self.precision,
                mu=self.mu,
                matrix=self.matrix,
                determinant=self.determinant,
                pulse=self.p,
                base=self.base,
                target=self.target,
                K=self.K,
                Kp=self.Kp,
            )
            # Keep these assignments explicit at the runtime boundary: every
            # downstream tangent and complete-end replay must use these live
            # objects, never atoms reconstructed from the receipt strings.
            self.algebra.M = self.matrix
            self.algebra.p = self.p
            self.algebra.det = self.determinant
            self.algebra.base = self.base
            self.algebra.v = self.algebra.affine(self.algebra.p)
            self.a, self.c = self.algebra.solve()
            self.algebra.retain_solution(self.a, self.c)
            self.component = _RuntimeAxialCorrection(self)
            self.receipt = self._make_receipt(solve)

    @staticmethod
    def _base_rows(seeded, receipt):
        incoming = seeded.get("incoming", {}) if isinstance(seeded, dict) else {}
        if "row_normalization" not in incoming:
            incoming = receipt.get("incoming", {}) if isinstance(receipt, dict) else {}
        rows = incoming.get("row_normalization", {})
        if rows:
            return [
                _decode(rows[key])
                for key in ("scaled_base_m1", "scaled_base_m2")
            ]
        axial = seeded.get("axial", {}) if isinstance(seeded, dict) else {}
        linear = axial.get("linear_rhs_inputs", {})
        if "base" not in linear:
            raise KeyError("No live incoming base rows supplied")
        return [_decode(value) for value in linear["base"]]

    @staticmethod
    def _target(seeded, solve):
        if seeded.get("energy_target") is not None:
            return mp.mpf(seeded["energy_target"])
        if solve.get("energy_target") is not None:
            return mp.mpf(solve["energy_target"])
        return mp.mpf(seeded["axial"]["energy_target"])

    def full_end(self, coefficients=None, row=1):
        """Return the complete end primitive from the stored matrix row."""

        if row not in (1, 2):
            raise ValueError("Expected row 1 or 2")
        values = self.c if coefficients is None else [
            _decode(value) for value in coefficients
        ]
        with mp.workdps(max(self.precision, mp.mp.dps)):
            return sum(
                self.matrix[row - 1][index] * value
                for index, value in enumerate(values)
            )

    @lru_cache(maxsize=4)
    def basis_enclosure(self,*,panels=4096,precision=60):
        """Optional certified basis bounds; pulse/input closure stays open."""
        from lei_ren_part1_paper_continuous_basis_enclosure import ContinuousBasisEnclosure
        encloser=ContinuousBasisEnclosure(precision=precision,
            ell=mp.nstr(self.basis.ell,self.precision))
        return encloser.report(self.solve_receipt['input_mu'],panels=panels,
                              nominal=self.basis)

    @lru_cache(maxsize=8)
    def pulse_enclosure(self,row,*,panels=4096,precision=80,band=48):
        """Complete pulse bounds using this runtime owned pulse provider."""
        from lei_ren_part1_paper_continuous_pulse_enclosure import ContinuousPulseEnclosure
        return ContinuousPulseEnclosure(precision=precision).report(
            self.solve_receipt['input_mu'],row,panels=panels,band=band,nominal=self.pulse)

    @lru_cache(maxsize=4)
    def coefficient_enclosure(self,*,panels=4096,precision=80,Md='.5'):
        """Conditional coefficient bounds, preserving inherited-input gaps."""
        from lei_ren_part1_paper_continuous_solve_enclosure import ContinuousSolveEnclosure
        return ContinuousSolveEnclosure(precision=precision).report(self,panels=panels,Md=Md)

    def _make_receipt(self, solve):
        precision = self.precision
        with mp.workdps(precision):
            row_replay = []
            for index in range(2):
                rhs = self.base[index] + self.a * self.p[index]
                row_replay.append(
                    mp.nstr(
                        abs((self.full_end(self.c, index + 1) + rhs) / rhs),
                        60,
                    )
                )
            energy = self.Kp * self.a * self.a + self.mu * sum(
                key * value * value for key, value in zip(self.K, self.c)
            )
            pulse_rows = _json_safe(self.pulse_rows, precision)
            base_signed = [_signed(value, precision) for value in self.base]
            return {
                "input_mu": mp.nstr(self.mu, precision),
                "a_p": mp.nstr(self.a, precision),
                "c": [_signed(value, precision) for value in self.c],
                "linear_matrix": [
                    [mp.nstr(value, precision) for value in row]
                    for row in self.matrix
                ],
                "linear_rhs_inputs": {
                    "base": base_signed,
                    "pulse": pulse_rows,
                },
                "pulse_rows": pulse_rows,
                "K_p": mp.nstr(self.Kp, self.energy_atom_precision),
                "K_bump": [mp.nstr(value, precision) for value in self.K],
                "u": [_signed(value, precision) for value in self.algebra.affine(self.base)],
                "v": [_signed(value, precision) for value in self.algebra.v],
                "energy_target": mp.nstr(self.target, precision),
                "linear_relative_replay": row_replay,
                "energy_relative_replay": mp.nstr(
                    abs(energy / self.target - 1), 60
                ),
                "energy_atom_precision": self.energy_atom_precision,
                "energy_atom_precision_is_not_certified_accuracy": True,
                "algebra_precision": precision,
                "inherited_base_rows": False,
                "inherited_energy_target": False,
                "incoming_rows_regenerated": True,
                "live_complete_atoms_owned": True,
                "live_basis_and_pulse_shared_by_identity": True,
                "basis_enclosure_available": True,
                "quadrature_enclosure_certified": False,
                "installed_in_global_profile": True,
                "global_mean_closed": False,
                "finite_energy_certified": False,
                "end_component_provider_available": True,
                "terminal_mean_forced_zero": False,
                "exact_functional_identity": (
                    "M c + b + a p = 0 for exact continuous atoms; "
                    "the materialized terminal residual remains observable."
                ),
            }


def _fixture():
    outer = json.loads(
        (HERE / "lei_ren_part1_paper_continuous_incoming_outer.json").read_text(
            encoding="utf-8"
        )
    )
    seeded = json.loads(
        (HERE / "lei_ren_part1_paper_seeded_shared_candidate.json").read_text(
            encoding="utf-8"
        )
    )
    seeded["incoming"] = outer["incoming"]
    solve = outer["continuous_solve"]
    seeded["energy_target"] = solve["energy_target"]
    seeded["axial"]["energy_target"] = solve["energy_target"]
    seeded["axial"]["linear_rhs_inputs"]["base"] = [
        outer["incoming"]["row_normalization"][key]
        for key in ("scaled_base_m1", "scaled_base_m2")
    ]
    return seeded, solve


def run():
    seeded, receipt = _fixture()
    runtime = SharedContinuousAxialRuntime(seeded, receipt)
    precision = runtime.precision
    with mp.workdps(precision):
        component = runtime.component
        identity = {
            "component_runtime": component.runtime is runtime,
            "component_basis": component.basis is runtime.basis,
            "component_pulse": component.pulse is runtime.pulse,
            "component_a": component.a is runtime.a,
            "component_c": component.c is runtime.c,
            "component_base": component.base is runtime.base,
            "algebra_matrix": runtime.algebra.M is runtime.matrix,
            "algebra_p": runtime.algebra.p is runtime.p,
            "algebra_base": runtime.algebra.base is runtime.base,
            "algebra_K": runtime.algebra.K is runtime.K,
            "cached_full_row_1": runtime.pulse.full_row(runtime.mu, 1)
            is runtime.pulse_rows[0],
            "cached_full_row_2": runtime.pulse.full_row(runtime.mu, 2)
            is runtime.pulse_rows[1],
        }

        rows = []
        for row in (1, 2):
            end = runtime.full_end(runtime.c, row)
            direct = runtime.base[row - 1] + runtime.a * runtime.p[row - 1] + end
            terminal = component.terminal_balance(row, precision)
            primitive = component.weighted_primitive(row, 0)
            rows.append(
                {
                    "row": row,
                    "full_end_matrix_dot": _signed(end, precision),
                    "full_end_weighted_primitive": _signed(primitive, precision),
                    "full_end_relative_consistency": mp.nstr(
                        _relative(end, primitive), 60
                    ),
                    "direct_balance": _signed(direct, precision),
                    "component_balance": _signed(terminal["value"], precision),
                    "balance_relative_consistency": mp.nstr(
                        _relative(direct, terminal["value"]), 60
                    ),
                    "nonzero_materialized_residual": bool(direct != 0),
                    "pulse_omitted_bound_present": (
                        terminal["pulse_atom"].get(
                            "log_relative_omitted_positive_bound"
                        )
                        is not None
                    ),
                }
            )

        report = {
            "diagnostic": "live shared continuous axial atom runtime",
            "input_receipt": "lei_ren_part1_paper_continuous_incoming_outer.json",
            "precision": precision,
            "object_identity": identity,
            "rows": rows,
            "runtime_receipt": runtime.receipt,
            "exact_functional_identity": (
                "M c + b + a p = 0 for exact continuous atoms"
            ),
            "materialized_residual_retained": True,
            "terminal_mean_forced_zero": False,
            "quadrature_enclosure_certified": False,
            "incoming_uncertainty_enclosed": False,
            "global_mean_closed": False,
            "finite_energy_certified": False,
            "scope": (
                "One live basis, pulse, complete pulse rows, matrix, and solve; "
                "no fullfield construction or global certificate."
            ),
        }

    basis_bounds=runtime.basis_enclosure(panels=1024,precision=60)
    if not all(basis_bounds['nominal_containment'].values()) or not basis_bounds['determinant_strictly_negative']:
        raise AssertionError('Live basis failed independent midpoint enclosure')
    report['basis_enclosure_summary']=dict(
        panels=basis_bounds['panels'],interval_precision=basis_bounds['interval_precision'],
        nominal_containment=basis_bounds['nominal_containment'],
        determinant_strictly_negative=basis_bounds['determinant_strictly_negative'],
        basis_quadrature_enclosed=True,
        complete_axial_quadrature_enclosed=False)

    pulse_bounds=[runtime.pulse_enclosure(row,panels=1024) for row in (1,2)]
    if not all(b['nominal_center_contained'] for b in pulse_bounds):
        raise AssertionError('Live pulse failed independent range enclosure')
    report['pulse_enclosure_summary']=[dict(row=b['row'],panels=b['panels'],
        nominal_center_contained=b['nominal_center_contained'],
        centered_relative_width=b['centered_relative_width'],
        pulse_quadrature_enclosed=True,complete_axial_closure_enclosed=False)
        for b in pulse_bounds]

    output = HERE / "lei_ren_part1_paper_continuous_axial_runtime.json"
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "all_identity_checks": all(identity.values()),
                "row1_nonzero_residual": rows[0]["nonzero_materialized_residual"],
                "row2_nonzero_residual": rows[1]["nonzero_materialized_residual"],
                "finite_energy_certified": report["finite_energy_certified"],
            }
        ),
        flush=True,
    )
    return report


if __name__ == "__main__":
    run()
