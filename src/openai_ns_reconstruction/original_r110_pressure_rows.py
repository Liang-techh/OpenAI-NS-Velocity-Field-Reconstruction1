"""Genuine pressure inlets on constant-U second-switch/post-power charts."""
import copy
import hashlib
import json

from .original_core_numeric import ROOT, PREFIX, _bounds
from .original_first_switch_numeric import OriginalFirstSwitchNumeric
from .original_first_switch_pressure_rows import OriginalFirstSwitchPressureRows
from .original_macro_pressure_rows import OriginalMacroPressureRows
from .original_r110_numeric import OriginalSecondSwitchNumeric, OriginalPostPowerNumeric


class OriginalR110PressureRows(OriginalMacroPressureRows):
    def __init__(self, *, micro=None, first_rows=None, second_rows=None, **kwargs):
        if micro is not None and kwargs:
            raise ValueError('Supply original R110 chart or constructor options')
        micro = micro if micro is not None else OriginalSecondSwitchNumeric(**kwargs)
        if not isinstance(micro, OriginalSecondSwitchNumeric):
            raise ValueError('Original second-switch/post-power source required')
        super().__init__(micro=micro)
        self.chart = micro.chart
        self.derivative_chart = 'switch_power' if self.chart == 'post_power' else 'switch_second'
        path = self.directory/(PREFIX+'compliant_current_generic_shear_loop_jet_bounds.json')
        self.loop = json.loads(path.read_bytes())['current_actual_loop_jet_log_bounds_by_chart'][self.derivative_chart]
        if self.loop['source_family'] != self.family or self.loop['primitive_units']['B'] != 'original B/Pstar; E=Utheta/Pstar':
            raise ValueError('Same original derivative family and Pstar units required')
        bound = self.read_log(self.loop['slow_phase_held_log_bounds']['A']['y0_Z0'])
        if bound is not None and _bounds(self.audit, micro.logN)[0] < _bounds(self.audit, bound)[1]:
            raise ValueError('Original selected N must admit local A/N')
        for suffix in ('first_switch_pressure_rows', 'r110_numeric'):
            path = self.directory/(PREFIX+'compliant_current_original_'+suffix+'.json')
            raw = path.read_bytes(); data = json.loads(raw)
            if not data['checks']['passed'] or data['source_family'] != self.family:
                raise ValueError('Accepted original R110 pressure prerequisites required')
            for base, group, target in ((self.directory, data['input_hashes'], self.hashes),
                                         (ROOT, data['implementation_hashes'], self.code_hashes)):
                for name, digest in group.items():
                    if target.get(name, digest) != digest or hashlib.sha256((base/name).read_bytes()).hexdigest() != digest:
                        raise ValueError('Changed original R110 pressure source: '+name)
                    target[name] = digest
            self.hashes[path.name] = hashlib.sha256(raw).hexdigest()
        if first_rows is None:
            view = object.__new__(OriginalFirstSwitchNumeric)
            view.__dict__ = dict(micro.__dict__)
            view._micro_cache = {}; view._upgraded = set()
            first_rows = OriginalFirstSwitchPressureRows(micro=view)
        if first_rows.ctx is not self.ctx or first_rows.iv is not self.iv or first_rows.o is not self.o:
            raise ValueError('Exact shared source and Taylor contexts required')
        if self.chart == 'post_power' and (second_rows is None or second_rows.chart != 'second_switch'):
            raise ValueError('Explicit genuine second-switch pressure inlet required')
        self.first_rows, self.second_rows = first_rows, second_rows
        self._constant_inputs = {}; self._post_rows = None
        name = 'src/openai_ns_reconstruction/original_r110_pressure_rows.py'
        self.code_hashes[name] = hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    def post_power_field(self):
        if self.chart != 'second_switch':
            raise ValueError('Construct post-power rows from the second-switch endpoint')
        if self._post_rows is None:
            self._post_rows = OriginalR110PressureRows(micro=self.micro.post_power_field(),
                first_rows=self.first_rows, second_rows=self)
        return self._post_rows

    def geometry(self, q):
        c, iv, o = self.ctx, self.iv, self.o
        qp = sum(_bounds(c, q))/2
        post = self.chart == 'post_power'
        if post:
            length = self.micro.post_length
            logR = iv.log(100)+self.micro.post_log*q+iv.mpf([0, _bounds(self.audit, o._width_error(2*(1-q)))[1]])
            logRp = c.log(100)+self.micro.post_log_point*qp
            theta, thetap = iv.exp(-length*q), c.exp(-self.micro.post_log_point*qp)
            if q._mpi_ == iv.mpf(1)._mpi_:
                logR, logRp = iv.log(110), c.log(110)
        else:
            length = iv.mpf([0, _bounds(self.audit, o.hcap)[1]])
            logR = iv.log(100)+iv.mpf([0, _bounds(self.audit, o._width_error(1+q))[1]])
            logRp = c.log(100)
            theta, thetap = iv.exp(-iv.mpf([0, _bounds(self.audit, o._width_error(q))[1]])), c.mpf(1)
        if q._mpi_ == iv.mpf(0)._mpi_:
            theta, thetap = iv.mpf(1), c.mpf(1)
        lo, hi = _bounds(self.audit, 1-theta)
        mass = iv.mpf([max(self.audit.mpf(0), lo), min(self.audit.mpf(1), hi)])
        return dict(logR=logR, logR_point=logRp, length=length, theta=theta,
                    theta_point=thetap, mass=mass, length_is_error_cover_only=not post,
                    defining_width_source='Lp=log(11/10)-2hbS' if post else 'original hbS log atom')

    def inputs(self, Z, q):
        key = Z, q._mpi_
        if key in self._constant_inputs:
            return self._constant_inputs[key]
        c, iv, J = self.ctx, self.iv, self.J
        owner = self.second_rows if self.chart == 'post_power' else self.first_rows
        incoming = owner.inputs(Z, iv.mpf(1)); geo = self.geometry(q)
        point = dict(Psi=incoming['point']['Psi'])
        truth = dict(Psi=incoming['truth']['Psi'])
        point['Mpsi'] = incoming['point']['Mpsi']+(point['Psi']-incoming['point']['Mpsi'])*(1-geo['theta_point'])
        truth['Mpsi'] = incoming['truth']['Mpsi']+(truth['Psi']-incoming['truth']['Mpsi'])*geo['mass']
        for ctx, values in ((c, point), (iv, truth)):
            values['Psi_y'] = J.constant(ctx, 0, 5)
        zero, zeroi = J.constant(c, 0, 5), J.constant(iv, 0, 5)
        source = dict(point=zero, truth=zeroi, normalized_value_point=zero,
                      normalized_value_truth=zeroi, parts={}, nonlinear_coefficient_error_norm=iv.mpf(0))
        out = dict(point=point, truth=truth, geometry=geo, source_Uy=source,
                   incoming=incoming, constant_U_original_source=True)
        self._constant_inputs[key] = out
        return out

    def uniform_background(self, Z, q):
        """U/Q bounds over every native prefix, not just the upper endpoint."""
        iv, J, o = self.iv, self.J, self.o
        v = self.inputs(Z, q); incoming = v['incoming']['truth']; psi = v['truth']['Psi']
        upper = _bounds(self.audit, v['geometry']['mass'])[1]
        mean = incoming['Mpsi']+(psi-incoming['Mpsi'])*iv.mpf([0, upper])
        z = J.variable(iv, iv.mpf(Z), 5); A, D = 1-z*z*o.core.delta_iv, 1-z*z
        u0 = 4*z+o.core.j_iv
        qaxis = (u0*(z*(1+o.core.delta_iv))-4*D)/A
        qpsi = (psi*(2*z)-mean*(z*(1-o.core.delta_iv))-o._derivative(mean)*D)/A
        return dict(U=u0+psi*o.core.epsilon_iv, Q=qaxis+qpsi*o.core.epsilon_iv)

    def evaluate(self, Z, coordinate, *, selected_sc_multiple=False):
        out = super().evaluate(Z, coordinate, selected_sc_multiple=selected_sc_multiple)
        q, _, _ = self.micro._coordinate(coordinate, False)
        post = self.chart == 'post_power'
        out['leading_U_y_factored'].update(definition='Original U_y=0 exactly on this chart',
            actual_source='prescribed post-switch constant U', actual_unmodified_U_y_exact_zero=True,
            log_scale_terms=out['leading_U_y_factored']['log_scale_terms'][:1])
        out['geometry'] = dict(log_physical_radius=out['geometry']['log_physical_radius'],
            original_log_width=self.o._summary(sum(_bounds(self.ctx, self.o.logh))/2, self.o.logh),
            native_length=self.o._summary(self.micro.post_log_point if post else self.ctx.mpf(0), self.geometry(q)['length']),
            native_length_is_error_cover_only=not post, native_hbS_is_parameter=True,
            R_definition=self.micro.radius_definition,
            native_measure='dy=Lp dq; Lp=log(11/10)-2hbS; no extra hbS' if post else 'dy=hbS ds',
            exact_R110_endpoint=post and q._mpi_ == self.iv.mpf(1)._mpi_,
            width_source_identity=copy.deepcopy(self.micro.width_identity))
        charts = ['first_micro', 'second_micro', 'frozen_macro', 'first_switch', 'second_switch']
        if post: charts.append('post_power')
        out['source_cover'].update(generic_derivative_chart=self.derivative_chart,
            cumulative_history_charts=charts, generic_source_domain=['0', '1'] if post else ['1', '2'],
            phase='frac(N*(Y+log(11/10)*q+2hbS*(1-q)-hbB*sc/2))' if post else 'frac(N*(Y+hbS*(1+s)-hbB*sc/2))')
        for name in ('genuine_macro_modal_Uy_source', 'actual_macro_geometry_and_phase', 'actual_large_Y_mean_relaxation'):
            out['contract'].pop(name, None)
        out['contract'].update(genuine_first_switch_fivejet_inlet=True, constant_U_mean_ODE_retained=True,
            source_backgrounds_do_not_override_pressure_inlet=True, exact_zero_Uy_does_not_zero_Qy=True,
            both_finite_N_endpoints_and_fast_phase_retained=True, uniform_prefix_background_available=True,
            actual_local_chart_measure=True)
        out['limits'].pop('actual_macro_U_Q_Z2_mean_Z3_Uy_Qy_Z1_with_factored_errors', None)
        out['limits'].pop('macro_signed_m15_integrated', None)
        out['limits'].update(actual_R110_U_Q_Z2_mean_Z3_Uy_Qy_Z1=True,
                             second_switch_post_power_signed_m15_integrated=False)
        return out
