"""Current finite-width histories on all three original bridge charts.

The unchanged bridge derivative/physical algebra is replayed locally. Only
its three parent-acquisition calls change, to coordinate-labelled current
integral histories. No rounded theta is used to select a microscopic point.
Bounds enclose the nonlinear field; they do not select production values.
"""
import ast
import copy
import json
import math
from pathlib import Path

import mpmath as mp

import lei_ren_part1_paper_compliant_bridge_mixed_C4 as native
from lei_ren_part1_paper_compliant_current_core_physical_assembly import CurrentCorePhysicalAssembly
from lei_ren_part1_paper_compliant_actual_bridge_integrals import (
    named_moments,value_enclosure,bridge_control_bindings)
from lei_ren_part1_paper_compliant_actual_bridge_switch import coefficient_lists,CoreContext
from lei_ren_part1_paper_compliant_actual_switch_mixed_C4 import current_actual_field_bindings
from lei_ren_part1_paper_compliant_actual_reference_restore_mixed_C4 import accepted,sha,HERE,PREFIX
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision
from lei_ren_part1_paper_compliant_inner_bridge_profiles import dress
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

RECEIPT=PREFIX+'current_actual_bridge_mixed_C4_check.json'
THEOREM=PREFIX+'bridge_mixed_C4_check.json'
GATE='current_actual_three_bridge_charts_mixed4_available'
VIEWS={'whole_first':('first',[0,1]),'core_exit':('first',0),
    'first_phase1':('first',1),'whole_second':('second',[1,2]),
    'second_phase1':('second',1),'smoothing_exit':('second',2),
    'whole_macro':('macro',[0,1]),'macro_inlet':('macro',0),'R100_exit':('macro',1)}
SCOPES=('actual_point_moment_history_recovered','full_point_physical_field_evaluation',
    'full_implicit_leading_inputs_recomputed','full_inner_interfaces_certified',
    'full_cartesian_vector_derivatives_certified','global_completed_tensor_admissibility',
    'admissible_stress_lift_constructed','temporal_recursion')


def current_coordinate_replay():
    """Check and replace only actual/inputs/comparison acquisition statements."""
    tree=ast.parse(Path(native.__file__).read_text(encoding='utf8'))
    cls=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=='CompliantBridgeMixedC4')
    original=next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name=='evaluate')
    fn=copy.deepcopy(original)
    expected={'actual':'self.bridge.actual(Z,theta)','inp':'self.bridge.inputs(Z)',
        'comparison':'self.bridge.comparison(Z,theta)'}
    located=[]
    for i,node in enumerate(fn.body):
        if not isinstance(node,ast.Assign) or len(node.targets)!=1:continue
        target=ast.unparse(node.targets[0])
        if target in expected:
            if ast.dump(node.value)!=ast.dump(ast.parse(expected[target],mode='eval').body):
                raise ValueError('Original bridge acquisition changed: '+target)
            located.append(i)
    if len(located)!=3 or located!=list(range(located[0],located[0]+3)):
        raise ValueError('Original three parent calls must remain adjacent and unique')
    replacement=ast.parse('actual,inp,comparison=self.current_parents(Z,value,chart)').body[0]
    fn.body[located[0]:located[-1]+1]=[replacement]
    # The remainder is byte-independent AST equality, including source-factor
    # cancellation and every original derivative formula and coordinate guard.
    restored=copy.deepcopy(fn)
    restored.body[located[0]:located[0]+1]=copy.deepcopy(original.body[located[0]:located[-1]+1])
    if ast.dump(restored)!=ast.dump(original):raise ValueError('Bridge replay changed derivative arithmetic')
    env=dict(native.__dict__)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),
        '<current coordinate-labelled bridge parents>','exec'),env)
    proof=dict(only_three_parent_acquisition_statements_replaced=True,
        unchanged_original_derivative_and_physical_AST=True,
        acquisition='actual,inp,comparison=self.current_parents(Z,value,chart)',
        coordinate_not_recovered_from_rounded_theta=True,passed=True)
    return env['evaluate'],proof


REPLAY,REPLAY_PROOF=current_coordinate_replay()


class CurrentActualBridgeMixedC4(native.CompliantBridgeMixedC4):
    @source_precision
    def __init__(self,owner=None,require_checked=True):
        self.owner=owner if owner is not None else CurrentCorePhysicalAssembly()
        self.anchor=self.owner.dispatch.anchor
        self.switch=self.anchor.patch.reference_mixed.long_mixed.switch_mixed
        self.history=self.switch.history;self.upstream=self.history.upstream
        self.bridge=self.history.bridge;self.core=self.owner.core
        self.ctx=c=self.upstream.ctx;self.family=self.owner.family;self.source=self.owner.source
        self.datum_sha=self.owner.datum_sha;self.hashes=dict(self.owner.hashes)
        self.proofs=[];self.core_bounds=[];self.cache={}
        self.logh=self.upstream.logh;self.h=c.mpf([0,endpoints(self.upstream.cap)[1]])
        self.logF0=self.switch.logF0;self.r=self.bridge.r
        linear=self.core.records['shared_linear_resolvent'];major=self.core.records['core_transfer']
        self.phi_norm=read_interval(c,linear['Phi_model_Xh_norm_upper'])+self.core.correction
        self.psi_norm=read_interval(c,major['fresh_Psi_model_Xh_norm_upper'])+self.core.correction
        self.shared=self.switch.shared_axial_source
        self.comparison_source=self.upstream.comparison.source
        self.core_norm_bindings={name:sha(name) for name in (
            'lei_ren_part1_paper_shared_linear_resolvent.json',PREFIX+'core_transfer.json',
            PREFIX+'core_physical_field.json','lei_ren_part1_paper_shared_analytic_tube.json')}
        self.retained=accepted(THEOREM,self.family,self.source,'bridge_mixed4_available')
        for name,digest in self.retained['input_hashes'].items():
            if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Retained bridge source conflict: '+name)
            self.hashes[name]=digest
        self.hashes[THEOREM]=sha(THEOREM);self.hashes.update(self.core_norm_bindings)
        self.bindings=bridge_control_bindings();self.actual_bindings=current_actual_field_bindings()
        self.graph=self.current_provider_graph()
        if not all(self.graph.values()):raise ValueError('Current finite-width bridge graph differs')
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.acceptance_loaded=False
        if require_checked:
            check=accepted(RECEIPT,self.family,self.source,GATE)
            if check['datum_enclosure_sha256']!=self.datum_sha:raise ValueError('Current bridge datum differs')
            self.hashes.update(check['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.acceptance_loaded=True

    def current_provider_graph(self):
        return dict(same_current_core_context=self.core is self.switch.core is self.bridge.core,
            same_actual_original_core=type(self.core) is CoreContext and self.core.original is self.upstream.core,
            original_profiles_and_full_norm_source_retained=self.core.records is self.upstream.core.records,
            same_native_history=self.history is self.anchor.patch.reference_mixed.long_mixed.history,
            same_native_actual_bridge_integrals=self.upstream is self.history.upstream,
            same_actual_and_comparison_context=self.ctx is self.upstream.comparison.ctx is self.bridge.ctx,
            same_family_source_datum=(self.family,self.source,self.datum_sha)==(
                self.upstream.core.family,self.upstream.core.source,self.upstream.core.datum.datum_sha),
            exact_current_width_source=self.logh is self.upstream.logh,
            same_current_switch_axial_source=self.shared is self.switch.shared_axial_source,
            unchanged_core_rectangular_extension=self.core_rectangular.__func__ is native.CompliantBridgeMixedC4.core_rectangular,
            unchanged_original_coordinate_ledger_operator=native.export_logR_ledger is REPLAY.__globals__['export_logR_ledger'])

    @source_precision
    def current_parents(self,Z,value,chart):
        """One coordinate supplies actual and known-comparison histories."""
        c=self.ctx
        incoming=self.upstream.packet(Z,value,chart)
        inp=self.upstream.prepare(Z)['inputs']
        phi=native.IntervalTaylor(c,incoming['actual_phi_axial5'])
        own=named_moments({name:native.IntervalTaylor(c,row)
            for name,row in incoming['actual_own_six_moments_axial5'].items()})
        actual=dict(F_actual_over_F0_axial5_coefficients=list(phi.coefficients),
            F_actual_true_axial5_divided_by_F0=list(dress(phi,inp['F0_ratios']).coefficients),
            Uz_actual_axial5_coefficients=incoming['actual_raw_V_axial5'],
            actual_moment_shape_axial5_coefficients=coefficient_lists(own),
            actual_Q_axial4_coefficients=incoming['actual_radial_Q_axial4'],
            pressure_axis_axial5_coefficients=incoming['pressure_axis_axial5'],
            pressure_increment_true_axial5_divided_by_R_F0_squared=incoming[
                'actual_pressure_increment_axial5_divided_by_R_F0_squared'],
            current_actual_integral_and_own_feedback_source=incoming,
            exact_coordinate_labelled_source='upstream.packet(Z,value,chart)',
            comparison_moments_substituted=False)
        lo,hi=endpoints(c.mpf(value))
        if chart=='macro':known=self.upstream.comparison.macro(Z,value)
        elif lo==hi:known=self.upstream.comparison.micro(Z,value)
        else:known=self.upstream.comparison_micro_cover(self.upstream.prepare(Z)['data'],lo,hi)
        comparison=dict(phi=value_enclosure(known['phi']),v=value_enclosure(known['V']),
            moments=named_moments({name:value_enclosure(row) for name,row in known['moments'].items()}))
        if (phi.order!=5 or comparison['phi'].order!=6 or comparison['v'].order!=6
                or any(len(row)!=6 for row in incoming['actual_own_six_moments_axial5'].values())):
            raise ValueError('Actual axial5 and independently defined comparison axial6 required')
        return actual,inp,comparison

    @source_precision
    def evaluate(self,Z,value,chart):
        packet=REPLAY(self,Z,value,chart)
        packet.update(datum_enclosure_sha256=self.datum_sha,
            current_coordinate_labelled_finite_width_history_used=True,
            current_actual_bridge_derivative_replay=REPLAY_PROOF,
            current_actual_bridge_mixed4_feedback_proved=True,
            actual_bridge_mixed4_feedback_installed=self.acceptance_loaded,
            **{GATE:self.acceptance_loaded},
            core_bridge_and_R100_local_functional_joins_certified=False,
            R100_functional_join_to_existing_switch_installed=False,
            **dict.fromkeys(SCOPES,False))
        return packet

    def manifest(self):
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            datum_enclosure_sha256=self.datum_sha,current_provider_graph_identity=self.graph,
            actual_parent_acquisition_replay=REPLAY_PROOF,
            actual_source_controls_bindings=self.bindings,actual_finite_width_field_bindings=self.actual_bindings,
            reused_unchanged_operator_fixtures={k:self.retained[k] for k in (
                'core_log_fixture','varying_direction_fixture','controls_fixture','symbolic_checks')},
            exact_current_axial_source=self.shared,exact_current_comparison_source=self.comparison_source,
            current_three_bridge_charts_source_ownership_installed=self.acceptance_loaded,
            **{GATE:self.acceptance_loaded},actual_bridge_mixed4_feedback_installed=self.acceptance_loaded,
            core_bridge_and_R100_local_functional_joins_certified=False,
            R100_functional_join_to_existing_switch_installed=False,
            all_33_current_source_charts_physical_spatial4_time1_mapped=False,
            **dict.fromkeys(SCOPES,False),input_hashes=dict(self.hashes))

    @source_precision
    def report(self):
        result=self.manifest();views={}
        for name,(chart,value) in VIEWS.items():
            views[name]=self.evaluate([-1,1],value,chart)
            print('Current finite-width bridge mixed4: '+name,flush=True)
        result.update(whole_current_bridge_mixed4=views,
            same_core_rectangular_extension_proofs=self.core_bounds,
            factored_positive_width_amplitude_cap_proofs=self.proofs,input_hashes=dict(self.hashes))
        return encode(pack(result))


def run():
    result=CurrentActualBridgeMixedC4(require_checked=False).report()
    Path(__file__).with_suffix('.json').write_bytes((json.dumps(result,indent=2)+'\n').encode('utf8'))
    print('Current three-chart actual bridge mixed4 generated; point/global/interface gates remain open',flush=True)
    return result


if __name__=='__main__':run()
