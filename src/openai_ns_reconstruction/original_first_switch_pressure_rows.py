"""First-switch pressure jets from the original R100 source and controls.

The switch parameter hbS remains distinct from bridge hbB. Their common
numeric width source is guarded upstream. Signed background proxies and
selected-N logarithmic errors remain separate.
"""
import copy
import hashlib
import json

import mpmath as mp

from .original_core_numeric import ROOT, PREFIX, _bounds
from .original_first_switch_numeric import OriginalFirstSwitchNumeric
from .original_macro_numeric import OriginalMacroNumeric
from .original_macro_pressure_rows import OriginalMacroPressureRows
from .original_n1_cutoff_numeric import stable_sigma_point


class OriginalFirstSwitchPressureRows(OriginalMacroPressureRows):
    def __init__(self, *, micro=None, **kwargs):
        if micro is not None and kwargs:
            raise ValueError('Supply retained first-switch source or constructor options')
        micro = micro if micro is not None else OriginalFirstSwitchNumeric(**kwargs)
        if not isinstance(micro, OriginalFirstSwitchNumeric):
            raise ValueError('Actual original first-switch source required')
        super().__init__(micro=micro)
        for suffix in ('macro_pressure_rows', 'first_switch_numeric'):
            path = self.directory/(PREFIX+'compliant_current_original_'+suffix+'.json')
            raw = path.read_bytes(); data = json.loads(raw)
            if not data['checks']['passed'] or data['source_family'] != self.family:
                raise ValueError('Accepted original first-switch prerequisite required: '+suffix)
            for root, group, target in ((self.directory, data['input_hashes'], self.hashes),
                                        (ROOT, data['implementation_hashes'], self.code_hashes)):
                for name, digest in group.items():
                    if (target.get(name, digest) != digest
                            or hashlib.sha256((root/name).read_bytes()).hexdigest() != digest):
                        raise ValueError('Changed original switch prerequisite: '+name)
                    target[name] = digest
            self.hashes[path.name] = hashlib.sha256(raw).hexdigest()
        path = self.directory/(PREFIX+'compliant_current_generic_shear_loop_jet_bounds.json')
        self.loop = json.loads(path.read_bytes())['current_actual_loop_jet_log_bounds_by_chart']['switch_first']
        if (self.loop['source_family'] != self.family or self.loop['chart'] != 'switch_first'
                or self.loop['primitive_units']['B'] != 'original B/Pstar; E=Utheta/Pstar'):
            raise ValueError('Original switch_first derivative source and Pstar units required')
        A = self.read_log(self.loop['slow_phase_held_log_bounds']['A']['y0_Z0'])
        if A is not None and _bounds(self.audit, micro.logN)[0] < _bounds(self.audit, A)[1]:
            raise ValueError('Selected original N must admit first-switch A/N')
        # Independent chart view, shared immutable source/caches and exact
        # Taylor contexts. It avoids virtual switch geometry/coordinates
        # being applied to the genuine upstream macro q=1 evaluation.
        view = object.__new__(OriginalMacroNumeric)
        view.__dict__ = dict(micro.__dict__)
        view._micro_cache = {}; view._upgraded = set()
        self.macro_rows = OriginalMacroPressureRows(micro=view)
        self._switch_inputs = {}
        name = 'src/openai_ns_reconstruction/original_first_switch_pressure_rows.py'
        self.code_hashes[name] = hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    def geometry(self, s):
        c, iv, o = self.ctx, self.iv, self.o
        tiny = o._width_error(s)
        logR = iv.log(100)+iv.mpf([0, _bounds(self.audit, tiny)[1]])
        return dict(logR=logR, logR_point=c.log(100),
                    length=iv.mpf([0, _bounds(self.audit, o.hcap)[1]]),
                    length_is_error_cover_only=True, defining_width_source='original hbS log atom',
                    theta_mass_are_error_covers_only=True,
                    theta=iv.exp(-iv.mpf([0, _bounds(self.audit, tiny)[1]])),
                    mass=iv.mpf([0, _bounds(self.audit, tiny)[1]]), theta_point=c.mpf(1))

    def source_Uy(self, prepared, geometry, order=5, *, phase=None, ell_switch=None):
        if phase is None or ell_switch is None:
            raise ValueError('Actual first-switch phase and angular prefix error required')
        o, iv, J, c = self.o, self.iv, self.J, self.ctx
        raw = super().source_Uy(prepared, geometry, order)
        with mp.workdps(280):
            complement = 1-self.t.sigma_jets(iv, phase)[0]
        sp = sum(_bounds(c, phase))/2
        chip = 1-stable_sigma_point(c, sp)
        # The macro bound already includes exp(ell_macro). The extra
        # first-switch angular exponent is genuinely hbS^2-small.
        extra = o._norm(raw['truth'])*ell_switch*iv.exp(ell_switch)
        truth = (raw['truth']+o._norm_error_jet(extra, order))*complement
        result = dict(raw)
        result.update(point=raw['point']*chip, truth=truth,
            normalized_value_truth=o._norm_error_jet(o._width_error(o._norm(truth)), order),
            normalized_value_point=J.constant(c, 0, order),
            parts={name: row*complement for name, row in raw['parts'].items()},
            nonlinear_coefficient_error_norm=raw['nonlinear_coefficient_error_norm']+extra,
            first_switch_complement=complement, first_switch_complement_point=chip,
            angular_prefix_error=ell_switch)
        return result

    def inputs(self, Z, s):
        key = Z, s._mpi_
        if key in self._switch_inputs:
            return self._switch_inputs[key]
        o, iv = self.o, self.iv
        incoming = self.macro_rows.inputs(Z, iv.mpf(1))
        ell, psi_error, evidence = self.micro._switch_increments(Z, s)
        point = {name: incoming['point'][name] for name in ('Psi', 'Mpsi')}
        truth = dict(Psi=incoming['truth']['Psi']+o._norm_error_jet(psi_error, 5))
        mean_error = o._width_error(s*(o._norm(incoming['truth']['Psi']-incoming['truth']['Mpsi'])+psi_error))
        truth['Mpsi'] = incoming['truth']['Mpsi']+o._norm_error_jet(mean_error, 5)
        geo = self.geometry(s)
        source = self.source_Uy(o._prepared[Z], geo, 5, phase=s, ell_switch=ell)
        point['Psi_y'], truth['Psi_y'] = source['normalized_value_point'], source['normalized_value_truth']
        out = dict(point=point, truth=truth, geometry=geo, source_Uy=source,
                   incoming=incoming, switch_evidence=evidence,
                   switch_Psi_change_error=psi_error, switch_mean_change_error=mean_error)
        self._switch_inputs[key] = out
        return out

    def evaluate(self, Z, coordinate, *, selected_sc_multiple=False):
        out = super().evaluate(Z, coordinate, selected_sc_multiple=selected_sc_multiple)
        s, _, recipe = self.micro._coordinate(coordinate, False)
        v = self.inputs(Z, s)
        width = out['leading_U_y_factored']['log_scale_terms'][1]
        width.update(original_native_parameter_name='hbS', source_log_atom_only=True,
                     bridge_hbB_graph_node_not_aliased=True,
                     numeric_identity_source='accepted shared upstream K1/logh source reboxing')
        out['leading_U_y_factored'].update(
            definition='U_y=epsilon*hbS*(1-sigma(s))*(original normalized current-radius drive with full angular error)',
            actual_source='-hbS*(1-sigma(s))*(Phi/barPhi)*(R*hydro+R*Pstar^2*pressure+R^2*F0^2*swirl)',
            angular_control='1', axial_control='1-sigma(s)',
            angular_phase_measure_power=2, axial_phase_measure_power=2,
            actual_unmodified_U_y_exact_zero_at_s1=s._mpi_ == self.iv.mpf(1)._mpi_)
        out['geometry'] = dict(log_physical_radius=out['geometry']['log_physical_radius'],
            original_log_width=copy.deepcopy(width),
            width_source_identity=copy.deepcopy(self.micro.width_identity),
            R_definition='100*exp(hbS*s)', native_measure='dy=hbS ds; R*dR=R^2*hbS ds',
            original_hbS_is_not_hb_to_a_power=True, defining_width_not_selected_from_cap=True,
            bridge_hbB_switch_hbS_graph_nodes_distinct=True,
            R100_inlet_exact=s._mpi_ == self.iv.mpf(0)._mpi_,
            R1_endpoint_exact=s._mpi_ == self.iv.mpf(1)._mpi_, R1_definition='100*exp(hbS)')
        out['source_cover'].update(generic_derivative_chart='switch_first',
            cumulative_history_charts=['first_micro', 'second_micro', 'frozen_macro', 'first_switch'],
            native_source_domain=['0', '1'],
            macro_prefix_measure_bound_only='Y*dq; actual macro measure=(Y-2hbB)dq; switch measure=hbS ds',
            phase='frac(N*(Y+hbS*s-hbB*original_sc/2))')
        out['first_switch_source_control_evidence'] = copy.deepcopy(v['switch_evidence'])
        for name in ('genuine_macro_modal_Uy_source', 'actual_macro_geometry_and_phase', 'actual_large_Y_mean_relaxation'):
            out['contract'].pop(name, None)
        out['contract'].update(genuine_first_switch_controls_and_modal_source=True,
            complete_macro_mean_and_Psi_inlet_retained=True, shared_exact_Taylor_contexts=True,
            local_angular_exponent_error_retained=True, original_hbS_source_atom_preserved=True,
            hbB_and_hbS_nodes_distinct_with_guarded_numeric_source_identity=True,
            exact_control_zero_at_s1_not_misused_for_finite_N=True)
        out['limits'].pop('actual_macro_U_Q_Z2_mean_Z3_Uy_Qy_Z1_with_factored_errors', None)
        out['limits'].pop('macro_signed_m15_integrated', None)
        out['limits'].update(actual_first_switch_U_Q_Z2_mean_Z3_Uy_Qy_Z1_with_factored_errors=True,
                             first_switch_signed_m15_integrated=False)
        return out
