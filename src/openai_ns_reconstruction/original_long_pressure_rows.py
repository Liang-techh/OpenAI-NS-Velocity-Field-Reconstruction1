"""Original long pressure with genuine five-jet inlet and logged mean memory."""
import copy
import hashlib
import json

from .original_core_numeric import ROOT, PREFIX, _bounds
from .original_long_numeric import OriginalLongNumeric
from .original_macro_pressure_rows import OriginalMacroPressureRows
from .original_r110_numeric import OriginalSecondSwitchNumeric
from .original_r110_pressure_rows import OriginalR110PressureRows


class OriginalLongPressureRows(OriginalR110PressureRows):
    def __init__(self, *, long=None, **kwargs):
        if long is not None and kwargs:
            raise ValueError('Supply retained original long source or options')
        self.long = long if long is not None else OriginalLongNumeric(**kwargs)
        if not isinstance(self.long, OriginalLongNumeric):
            raise ValueError('Actual original long source required')
        # Reuse the admitted log/error machinery with the long source cover,
        # without asking the truncated long public velocity for pressure data.
        view = object.__new__(OriginalSecondSwitchNumeric)
        view.__dict__ = dict(self.long.owner.__dict__)
        view._micro_cache = {}; view._upgraded = set(); view.chart = 'long_reshape'
        view._covers = self.long._covers
        OriginalMacroPressureRows.__init__(self, micro=view)
        self.chart, self.derivative_chart = 'long_reshape', 'reshape'
        self._locations = {}; self._long_inputs = {}
        view._coordinate = self._coordinate
        self.r110_rows = OriginalR110PressureRows(micro=self.long.owner)
        self.inlet_rows = self.r110_rows.post_power_field()
        self.hashes.update(self.long.hashes); self.code_hashes.update(self.long.code_hashes)
        self.hashes.update(self.inlet_rows.hashes); self.code_hashes.update(self.inlet_rows.code_hashes)
        path = self.directory/(PREFIX+'compliant_current_generic_shear_loop_jet_bounds.json')
        self.loop = json.loads(path.read_bytes())['current_actual_loop_jet_log_bounds_by_chart']['reshape']
        if self.loop['source_family'] != self.family or self.loop['primitive_units']['B'] != 'original B/Pstar; E=Utheta/Pstar':
            raise ValueError('Actual reshape derivative source and Pstar units required')
        bound = self.read_log(self.loop['slow_phase_held_log_bounds']['A']['y0_Z0'])
        if bound is not None and _bounds(self.audit, self.micro.logN)[0] < _bounds(self.audit, bound)[1]:
            raise ValueError('Original selected N must admit reshape A/N')
        for suffix in ('r110_m15_numeric', 'long_numeric'):
            path = self.directory/(PREFIX+'compliant_current_original_'+suffix+'.json')
            raw = path.read_bytes(); data = json.loads(raw)
            if not data['checks']['passed'] or data['source_family'] != self.family:
                raise ValueError('Checked original long pressure prerequisite required')
            for base, group, target in ((self.directory, data['input_hashes'], self.hashes),
                                         (ROOT, data['implementation_hashes'], self.code_hashes)):
                for name, digest in group.items():
                    if target.get(name, digest) != digest or hashlib.sha256((base/name).read_bytes()).hexdigest() != digest:
                        raise ValueError('Changed original long pressure prerequisite: '+name)
                    target[name] = digest
            self.hashes[path.name] = hashlib.sha256(raw).hexdigest()
        name = 'src/openai_ns_reconstruction/original_long_pressure_rows.py'
        self.code_hashes[name] = hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    def _coordinate(self, coordinate, selected=False):
        c, iv = self.ctx, self.iv
        if selected or type(coordinate) is not str:
            raise ValueError('Original long decimal fraction or y:log-offset required')
        native = coordinate.startswith('y:')
        raw = coordinate[2:] if native else coordinate
        if not c.isfinite(c.mpf(raw)):
            raise ValueError('Finite original long coordinate required')
        if native:
            yp, y = c.mpf(raw), iv.mpf(raw)
            q = y/self.long.T
        else:
            q = iv.mpf(raw); yp, y = self.long.Tp*c.mpf(raw), self.long.T*q
        lo, hi = _bounds(self.audit, q)
        if lo < 0 or hi > 1:
            raise ValueError('Original long coordinate must lie in the finite original interval')
        recipe = dict(kind='exact_decimal_R110_log_offset' if native else 'exact_decimal_long_fraction', value=raw)
        self._locations[q._mpi_] = yp, y, recipe
        return q, False, recipe

    def geometry(self, q):
        c, iv = self.ctx, self.iv
        yp, y, recipe = self._locations.get(q._mpi_,
            (self.long.Tp*sum(_bounds(c, q))/2, self.long.T*q, None))
        proof = []
        theta = self.long.kernels.positive_decay(iv, 1, y, proof)
        if q._mpi_ == iv.mpf(0)._mpi_: theta = iv.mpf(1)
        lo, hi = _bounds(self.audit, 1-theta)
        return dict(logR=iv.log(110)+y, logR_point=c.log(110)+yp,
            length=self.long.T, theta=theta,
            mass=iv.mpf([max(self.audit.mpf(0), lo), min(self.audit.mpf(1), hi)]),
            y=y, y_point=yp, mean_decay_proofs=proof, recipe=recipe,
            defining_length_source='original exact singleton T=400A', length_is_error_cover_only=False)

    def inputs(self, Z, q):
        key = Z, q._mpi_
        if key in self._long_inputs: return self._long_inputs[key]
        c, iv, J, o = self.ctx, self.iv, self.J, self.o
        incoming = self.inlet_rows.inputs(Z, iv.mpf(1)); geo = self.geometry(q)
        point, truth = dict(Psi=incoming['point']['Psi']), dict(Psi=incoming['truth']['Psi'])
        yp, y, proof = geo['y_point'], geo['y'], []
        dp = incoming['point']['Mpsi']-point['Psi']
        dt = incoming['truth']['Mpsi']-truth['Psi']
        upper = _bounds(self.audit, y)[1]
        direct = upper < 10000
        if direct:
            ep, et = c.exp(-yp), iv.exp(-y)
            point['Mpsi'] = incoming['point']['Mpsi'] if not yp else point['Psi']+dp*ep
            truth['Mpsi'] = incoming['truth']['Mpsi'] if upper == 0 else truth['Psi']+dt*et
        else:
            errors = []
            for k in range(6):
                magnitude = o._abs(dt[k])
                log = None if _bounds(self.audit, magnitude)[1] == 0 else iv.log(magnitude)-y
                if log is None:
                    error = iv.mpf(0)
                elif _bounds(self.audit, log)[1] <= _bounds(self.audit, iv.log(iv.mpf('1e-1000')))[0]:
                    error = self.long.owner._cap(log, proof, 'long_pressure_mean_memory_Z'+str(k))
                else:
                    # The combined log is bounded above by log|dt[k]|, so
                    # any positive exponent here is no larger than the
                    # already represented coefficient. A huge negative
                    # exponent is capped by the accepted bound_exp path.
                    error = self.long.kernels.bound_exp(iv, log, proof)
                errors.append(o._symmetric(error))
            point['Mpsi'] = point['Psi']
            truth['Mpsi'] = truth['Psi']+J(iv, errors)
        zero, zeroi = J.constant(c, 0, 5), J.constant(iv, 0, 5)
        point['Psi_y'], truth['Psi_y'] = zero, zeroi
        source = dict(point=zero, truth=zeroi, normalized_value_point=zero,
                      normalized_value_truth=zeroi, parts={}, nonlinear_coefficient_error_norm=iv.mpf(0))
        result = dict(point=point, truth=truth, geometry=geo, source_Uy=source, incoming=incoming,
            mean_memory_proofs=proof, source_mean_memory_exact='(Mpsi_R110-Psi_R110)*exp(-T*q)',
            mean_memory_cap_is_error_only=not direct)
        self._long_inputs[key] = result
        return result

    def at_log_radius(self, Z, y):
        return self.evaluate(Z, 'y:'+y)

    def evaluate(self, Z, coordinate, *, selected_sc_multiple=False):
        out = OriginalMacroPressureRows.evaluate(self, Z, coordinate, selected_sc_multiple=selected_sc_multiple)
        q, _, _ = self._coordinate(coordinate, False); v = self.inputs(Z, q); geo = v['geometry']
        out['leading_U_y_factored'].update(definition='Original long prescribed U_y=0 exactly',
            actual_source='constant R110 axial source', actual_unmodified_U_y_exact_zero=True,
            log_scale_terms=out['leading_U_y_factored']['log_scale_terms'][:1])
        out['geometry'] = dict(log_physical_radius=out['geometry']['log_physical_radius'],
            native_length=self.o._summary(self.long.Tp, self.long.T),
            original_T_singleton_exact_mpf_tuple=self.long.T_tuple,
            native_log_offset=self.o._summary(geo['y_point'], geo['y']),
            R_definition='110*exp(T*q)', native_measure='dy=T dq; no hb factor',
            radius_and_enormous_exponentials_materialized=False,
            direct_decay_only_if_y_below_10000=True, finite_original_length_not_infinite_window=True,
            mean_memory_proofs=self.long._serializable_proof(v['mean_memory_proofs']))
        out['source_cover'].update(generic_derivative_chart='reshape',
            cumulative_history_charts=['first_micro', 'second_micro', 'frozen_macro', 'first_switch', 'second_switch', 'post_power', 'long_reshape'],
            prefix_physical_measure='each chart measure preserved; long T*dq once',
            phase='frac(N*(Y+log(11/10)+T*q-hbB*original_sc/2))')
        for name in ('genuine_macro_modal_Uy_source', 'actual_macro_geometry_and_phase', 'actual_large_Y_mean_relaxation'):
            out['contract'].pop(name, None)
        out['contract'].update(genuine_R110_mean_Psi_fivejets=True, source_mean_ODE_through_Z5=True,
            original_long_T_singleton_retained=True, normalized_decay_coefficients_combined_before_cap=True,
            selected_N_mean_Z3_error_kept_logarithmic=True, no_zero_padding_of_mean_Z3=True,
            reshape_slow_and_fast_derivative_owners=True, complete_upstream_and_long_prefix=True,
            exact_zero_background_Uy_does_not_zero_corrected_errors=True)
        out['limits'].pop('actual_macro_U_Q_Z2_mean_Z3_Uy_Qy_Z1_with_factored_errors', None)
        out['limits'].pop('macro_signed_m15_integrated', None)
        out['limits'].update(actual_long_U_Q_Z2_mean_Z3_Uy_Qy_Z1_with_errors=True,
                             long_signed_m15_integrated=False)
        return out
