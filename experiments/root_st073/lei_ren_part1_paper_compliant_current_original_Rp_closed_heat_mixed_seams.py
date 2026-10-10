"""Same-current waiting/collar/exterior joins in original raw mixed units.

Original whole-Z terminal identities identify the functions before directed
enclosures are compared. Native phase rows keep the exact waiting Jacobian.
This adapter does not install a uniform norm or a physical point/time field.
"""
import gzip
import json
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_Rp_pulse_mixed_seams as pulse_seams
from lei_ren_part1_paper_compliant_collar_Gamma_C4_check import functional_source_identities
from lei_ren_part1_paper_compliant_current_heat_source import same_exact_Gamma_future_bindings
from lei_ren_part1_paper_compliant_actual_Rsh_source_join import class_assignment
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

mixed=pulse_seams.mixed
HERE,PREFIX,sha=mixed.HERE,mixed.PREFIX,mixed.sha
NAME=PREFIX+'current_original_Rp_closed_heat_mixed_seams.json.gz'
RECEIPT=PREFIX+'current_original_Rp_closed_heat_mixed_seams_check.json'
GATES=('current_original_Rp_waiting_closed_collar_raw_mixed4_function_seam_installed',
       'current_original_Rp_closed_collar_Gamma_raw_mixed4_function_seam_installed')
OPEN=mixed.OPEN
SEAMS={'waiting_collar':('waiting',1,'heat_collar',0),
       'collar_exterior':('heat_collar',3,'heat_exterior',3)}
OBSERVATIONS={'waiting':('waiting',1),'collar_in':('heat_collar',0),
              'collar_out':('heat_collar',3),'exterior':('heat_exterior',3)}
PAIRS={'waiting_collar':('waiting','collar_in'),'collar_exterior':('collar_out','exterior')}


def closed_output_bindings():
    """Bind the actual closed provider rather than its legacy companion."""
    cls='CurrentOriginalRpSameRepairHeatClosure';stem='current_original_Rp_same_repair_heat_closure'
    specs={'packet':'self.pressure.evaluate(chart,Z,coordinate)',
        'forward':"self.angular.terminal_constants(z)['current_repaired_forward_terminal']",
        'theta':"shape['K_rows'][0]*forward['theta_base']*c.exp(-self.angular.heat.bh*t)",
        'X':"packet['actual_angular_Taylor_after_source_closure']",
        'P0':"packet['original_analytic_P0_Taylor_retained']",
        'P':"-remaining*forward['pressure_scale']",
        'jets':'dict(Mz=zero,Mtheta=theta*X*c.sqrt(2),Mtheta_z=zero,\n            Mztheta=theta*theta*energy,Mp=P-P0,P0=P0,pressure=P,Utheta=theta)'}
    return {key:class_assignment(stem,cls,'evaluate',key,value) for key,value in specs.items()}


class CurrentOriginalRpClosedHeatMixedSeams:
    @source_precision
    def __init__(self,before=None,require_checked=True):
        self.before=before if before is not None else pulse_seams.CurrentOriginalRpPulseMixedSeams()
        if type(self.before) is not pulse_seams.CurrentOriginalRpPulseMixedSeams or not self.before.acceptance_loaded:
            raise ValueError('Accepted typed current pulse mixed seam owner required')
        self.transport=self.before.before;self.owner=self.transport.owner
        self.closed=self.owner.closed;self.post=self.owner.post;self.radius=self.owner.radius
        self.graph=self.transport.graph;self.ctx=self.transport.ctx;self.family_record=self.transport.family_record
        self.canonical=functional_source_identities()
        self.shared_future=same_exact_Gamma_future_bindings()
        self.output_bindings=closed_output_bindings()
        if self.canonical!=self.post.heat_source.functional_proof or not all(self.canonical.values()):
            raise ValueError('Original arbitrary-function heat seam theorem differs')
        if self.shared_future!=self.post.heat_source.bindings['shared_exact_Gamma_future_binding']:
            raise ValueError('Actual current complete Gamma future differs')
        name=PREFIX+'collar_Gamma_C4_check.json';canonical=json.loads((HERE/name).read_bytes())
        flags=('all_passed','waiting_collar_and_collar_Gamma_joins_certified',
               'full_infinite_Gamma_source_and_formal_nonzero_S_retained')
        self.canonical_flags={key:canonical[key] for key in flags}
        if not all(self.canonical_flags.values()) or self.canonical!=canonical['functional_production_source_identities']:
            raise ValueError('Accepted original full-Gamma seam function certificate required')
        self.hashes=dict(self.before.hashes)
        mixed.pulse.radius.post.selected.inlet.add_hashes(self.hashes,canonical['input_hashes'])
        for dep in (name,pulse_seams.NAME,pulse_seams.RECEIPT,Path(__file__).name,
                    PREFIX+'collar_Gamma_C4_check.py',PREFIX+'current_heat_source.py',
                    PREFIX+'actual_Rsh_source_join.py'):
            self.hashes[dep]=sha(dep)
        self.acceptance_loaded=False;self.call_trace=[]
        self.function_transfer=self.source_function_transfer();self.assert_graph()
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not all(receipt[k] for k in GATES) or any(receipt[k] for k in OPEN):
                raise ValueError('Closed heat seam receipt/scope differs')
            if receipt['source_family']!=self.family_record:raise ValueError('Closed heat seam source family differs')
            mixed.pulse.radius.post.selected.inlet.add_hashes(self.hashes,receipt['input_hashes'])
            self.acceptance_loaded=True

    def source_function_transfer(self):
        """Require current exact closure, full future and generic seam rows."""
        a=self.closed.angular.proof;p=self.closed.pressure.proof
        angular=('passed','current_exact_native_angular_branch_identified_by_uniform_uniqueness',
            'same_exact_preheat_collarJ_and_current_JW_function',
            'same_current_exact_inverse_radius_and_full_Gamma_Theta',
            'exact_function_identity_implies_whole_Z_axial5','interval_overlap_not_used_to_zero_constant')
        pressure=('passed','accepted_current_collar_waiting_measure_function_proof',
            'same_reference_down_to_zero_and_native_cumulative_FTC',
            'all_original_raw_stage_integrals_are_one_function_not_box_choices',
            'original_P0_not_redefined_or_pressure_patched',
            'same_exact_inverse_radius_and_Pstar_squared_units_preserved',
            'full_infinite_tail_domination','holomorphic_integral_and_axial5_differentiation_justified')
        required={**{'angular:'+k:a[k] for k in angular},**{'pressure:'+k:p[k] for k in pressure}}
        if not all(required.values()):raise ValueError('Current original forward closure function theorem required')
        if not all(self.shared_future[k] for k in ('passed','both_energy_epsilon_atoms_retained',
                'same_infinite_Gamma_tail_callable','all_three_scaled_Gamma_future_integrands_bound')) \
                or self.shared_future['interval_overlap_used_as_join_proof']:
            raise ValueError('Same original complete Gamma future function required')
        identities={}
        for seam,prefix in (('waiting_collar','waiting_collar'),('collar_exterior','collar_Gamma')):
            for label in ('theta','X','energy','pressure'):
                for k in range(5):
                    for j in range(5-k):
                        key=prefix+'_'+label+'_y%d_Z%d'%(k,j)
                        if not self.canonical.get(key):raise ValueError('Canonical seam derivative missing: '+key)
                        identities[seam+'_'+label+'_y%d_Z%d'%(k,j)]=True
                for j in range(6):
                    if not self.canonical.get(prefix+'_'+label+'_axial'+str(j)):
                        raise ValueError('Canonical whole-Z fifth derivative seam missing')
        return dict(passed=True,current_original_terminal_function_flags=required,
            same_complete_future_and_epsilon_energy_normalization=True,
            same_original_P0_plus_full_native_pressure_integral_identified=True,
            original_forward_Dtheta_and_Cp_identified_as_zero_functions=True,
            canonical_current_primitive_mixed4_identities=identities,
            raw_transfer='Mtheta=sqrt(2)*R^(3/2)*Ev0*theta*X; Mztheta=R*Ev0^2*theta^2*energy; Mp=P-P0; homogeneous meridional terminal histories remain zero. Positive y derivatives use the admitted general cumulative density equations.',
            ordinary_coordinates='D_phase=waiting*D_y; D_t=D_y in collar/exterior. Exact positive Jacobian powers remain separate.',
            mixed_scope='k+j<=4, whole-Z base identities through axial order5',
            independent_P0_not_reset_and_absolute_pressure_not_patched=True,
            interval_overlap_is_not_the_function_proof=True)

    def assert_graph(self):
        e=self.closed.exact;s=self.owner.selected
        result=dict(accepted_current_pulse_and_mixed_owners=self.before.acceptance_loaded
                and all(self.before.assert_graph().values()) and self.transport.acceptance_loaded,
            same_actual_fifteen_chart_provider=self.transport.owner is self.owner and all(self.transport.assert_graph().values()),
            same_closed_current_original_forward_provider=self.closed is self.owner.closed
                and self.closed.acceptance_loaded and all(self.closed.assert_graph().values()),
            same_current_unique_repair=e.repair is self.post.heat.repair is s.seed.exact.repair,
            same_current_full_future=e.future is self.post.steep.future is self.post.heat.future is s.future,
            same_original_independent_P0=s.datum is e.flatten.inlet.datum is self.closed.raw.flat.inlet.datum,
            same_current_context=self.ctx is self.closed.ctx is self.post.ctx,
            same_current_expression_graph=self.graph is self.owner.graph is self.radius.graph,
            current_function_proof=self.function_transfer['passed'],
            original_full_Gamma_certificate=all(self.canonical_flags.values()),
            no_second_repair_or_future_replay=True)
        if not all(result.values()):raise ValueError('Closed heat seam graph differs: '+str(result))
        return result

    def assemble(self,name,Z,left,right):
        if name not in SEAMS:raise ValueError('Current waiting/collar or collar/exterior seam required')
        self.assert_graph();lc,lv,rc,rv=SEAMS[name]
        for view,chart,coordinate in ((left,lc,lv),(right,rc,rv)):
            g=view['geometry'];coord=g['exact_native_coordinate']
            if view['chart']!=chart or coord!={'numerator':coordinate,'denominator':1}:
                raise ValueError('Actual exact seam source route required')
            if mixed.pulse.radius.post.selected.inlet.endpoints(view['Z'])!=mixed.pulse.radius.post.selected.inlet.endpoints(self.ctx.mpf(Z)):
                raise ValueError('Same actual axial source observation required')
            if set(view['log_radius_mixed_rows'])!=set(left['log_radius_mixed_rows']) or len(view['log_radius_mixed_rows'])!=10:
                raise ValueError('Complete ten-output mixed source required')
            for grid in view['log_radius_mixed_rows'].values():
                if len(grid)!=15 or any(type(row) is not mixed.FactorizedMixedSourceRow
                        or row.coefficients.ctx is not self.ctx
                        or any(v.graph is not self.graph for _,v in row.log_scale_parts) for row in grid.values()):
                    raise ValueError('Actual same-graph typed mixed rows required')
            if chart.startswith('heat_'):
                packet=view['original_forward_source_packet']
                if not all(packet[k] for k in ('pressure_is_original_forward_function_after_terminal_identity',
                        'original_axis_pressure_not_replaced_or_tail_patched','whole_Z_pressure_function_identity_consumed')):
                    raise ValueError('Actual closed original heat pressure source required')
                forward=self.closed.angular.terminal_constants(self.ctx.mpf(Z))['current_repaired_forward_terminal']
                if packet['exact_pressure_scale'] is not forward['pressure_scale']:
                    raise ValueError('Same source-bound forward pressure scale required')
        self.call_trace.append(dict(seam=name,actual_typed_current_source_routes=(lc,rc),
            exact_function_identity_before_interval_diagnostics=True))
        return dict(seam=name,Z=self.ctx.mpf(Z),left=left,right=right,
            same_current_original_pressure_and_complete_energy_memory=True,
            native_derivatives_use_exact_waiting_Jacobian_not_an_enclosure_multiplier=True,
            **dict.fromkeys(GATES,self.acceptance_loaded),**dict.fromkeys(OPEN,False))

    @source_precision
    def evaluate(self,name,Z):
        if name not in SEAMS:raise ValueError('Current closed heat seam required')
        lc,lv,rc,rv=SEAMS[name]
        return self.assemble(name,Z,self.transport.evaluate(lc,Z,lv),self.transport.evaluate(rc,Z,rv))


def report(view):
    return {**{k:v for k,v in view.items() if k not in ('left','right')},
        'left':mixed.view_report(view['left']),'right':mixed.view_report(view['right'])}


@source_precision
def run(before=None,observed_transport=None,observed_views=None):
    began=time.monotonic();owner=CurrentOriginalRpClosedHeatMixedSeams(before,require_checked=False)
    if observed_transport is not None:
        if observed_transport is not owner.transport or set(observed_views)!=set(OBSERVATIONS):
            raise ValueError('Complete typed observations of the same actual current transport required')
        views={name:owner.assemble(name,'.521',observed_views[l],observed_views[r]) for name,(l,r) in PAIRS.items()}
    else:views={name:owner.evaluate(name,'.521') for name in SEAMS}
    result=dict(source_family=owner.family_record,candidate_current_closed_heat_mixed_seams_constructed=True,
        actual_source_graph=owner.assert_graph(),original_heat_function_theorem=owner.canonical,
        same_current_complete_Gamma_future_binding=owner.shared_future,
        actual_closed_output_AST_bindings=owner.output_bindings,current_original_function_transfer=owner.function_transfer,
        accepted_full_Gamma_function_certificate_flags=owner.canonical_flags,
        actual_two_closed_heat_mixed_seam_views={name:report(view) for name,view in views.items()},
        actual_source_call_trace=owner.call_trace,fresh_Z='.521',
        scope='Two same-current original raw logR/Z function seams through total mixed order4; no uniform/global/physical-point/time admission',
        **dict.fromkeys(GATES+OPEN,False),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began)
    data=json.dumps(mixed.pulse.raw.packed(result),separators=(',',':'))+'\n'
    (HERE/NAME).write_bytes(gzip.compress(data.encode(),mtime=0))
    return owner,views


if __name__=='__main__':run()
