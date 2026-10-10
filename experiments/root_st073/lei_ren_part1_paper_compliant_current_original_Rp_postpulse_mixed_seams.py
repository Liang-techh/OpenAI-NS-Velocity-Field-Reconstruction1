"""All eight current original raw postpulse mixed function joins.

Six nonheat joins instantiate the original power/steep source theorems.
The two closed heat joins retain their accepted original function proof.
Exact radius/amplitude/Jacobian scales stay separate from directed rows.
Uniform norms, resolved physical points and temporal recursion remain open.
"""
import gzip
import json
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_Rp_closed_heat_mixed_seams as heat_seams
from lei_ren_part1_paper_compliant_current_postpulse_interfaces import flatten_power_current_theorem
from lei_ren_part1_paper_compliant_power_angular_C4_check import functional_source_identities as power_theorem
from lei_ren_part1_paper_compliant_steep_waiting_C4_check import functional_source_identities as steep_theorem
from lei_ren_part1_paper_compliant_current_power_angular_source import exact_prefix_and_decomposition_bindings
from lei_ren_part1_paper_compliant_current_steep_waiting_source import current_source_bindings
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

mixed=heat_seams.mixed
HERE,PREFIX,sha=mixed.HERE,mixed.PREFIX,mixed.sha
NAME=PREFIX+'current_original_Rp_postpulse_mixed_seams.json.gz'
RECEIPT=PREFIX+'current_original_Rp_postpulse_mixed_seams_check.json'
GATES=('current_original_Rp_six_nonheat_postpulse_raw_mixed4_function_seams_installed',
       'current_original_Rp_eight_postpulse_raw_mixed4_function_seams_installed',
       'current_original_Rp_fourteen_interface_source_chain_caller_installed')
OPEN=mixed.OPEN
NONHEAT={'flatten_power':('flatten',100,'outer_power',0),
    'power_angular':('outer_power',1,'outer_angular',-4),
    'angular_steep':('outer_angular',0,'steep_entry',0),
    'entry_power':('steep_entry',1,'steep_power',0),
    'power_exit':('steep_power',1,'steep_exit',0),
    'exit_waiting':('steep_exit',1,'waiting',0)}
SEAMS={**NONHEAT,**heat_seams.SEAMS}
OBSERVATIONS={'flatten_out':('flatten',100),'power_in':('outer_power',0),
    'power_out':('outer_power',1),'angular_in':('outer_angular',-4),
    'angular_out':('outer_angular',0),'entry_in':('steep_entry',0),
    'entry_out':('steep_entry',1),'steep_in':('steep_power',0),
    'steep_out':('steep_power',1),'exit_in':('steep_exit',0),
    'exit_out':('steep_exit',1),'wait_in':('waiting',0)}
PAIRS={'flatten_power':('flatten_out','power_in'),'power_angular':('power_out','angular_in'),
    'angular_steep':('angular_out','entry_in'),'entry_power':('entry_out','steep_in'),
    'power_exit':('steep_out','exit_in'),'exit_waiting':('exit_out','wait_in')}


class CurrentOriginalRpPostpulseMixedSeams:
    @source_precision
    def __init__(self,before=None,require_checked=True):
        self.before=before if before is not None else heat_seams.CurrentOriginalRpClosedHeatMixedSeams()
        if type(self.before) is not heat_seams.CurrentOriginalRpClosedHeatMixedSeams or not self.before.acceptance_loaded:
            raise ValueError('Accepted typed current closed-heat seam owner required')
        self.transport=self.before.transport;self.owner=self.transport.owner
        self.post=self.owner.post;self.radius=self.owner.radius;self.closed=self.owner.closed
        self.graph=self.transport.graph;self.ctx=self.transport.ctx;self.family_record=self.transport.family_record
        self.canonical={'power':power_theorem(),'steep':steep_theorem()}
        self.flatten_theorem=flatten_power_current_theorem()
        self.prefix=exact_prefix_and_decomposition_bindings()
        self.steep_binding=current_source_bindings(self.post.steep_source)
        if self.prefix!=self.post.power.bindings['current_prefix_and_future_decomposition']:
            raise ValueError('Current full C4/C5 future prefix definition differs')
        if self.steep_binding!=self.post.steep_source.bindings:
            raise ValueError('Actual current repaired angular/steep energy and pressure histories differ')
        self.hashes=dict(self.before.hashes);self.canonical_flags={}
        for group,live,gate in (('power',self.post.power,'flatten_power_and_power_angular_joins_certified'),
                ('steep',self.post.steep_source,'angular_steep_and_internal_joins_certified')):
            name=PREFIX+('power_angular_C4' if group=='power' else 'steep_waiting_C4')+'_check.json'
            receipt=json.loads((HERE/name).read_bytes())
            family=(receipt['actual_five_defect_family_sha256'],receipt['implicit_source_sha256'])
            if family!=(self.post.family,self.post.source):raise ValueError('Original postpulse seam theorem family differs')
            self.canonical_flags[group]={key:receipt[key] for key in ('all_passed',gate)}
            if not all(self.canonical_flags[group].values()) or not all(self.canonical[group].values()) \
                    or self.canonical[group]!=receipt['functional_production_source_identities'] \
                    or self.canonical[group]!=live.functional_proof:
                raise ValueError('Accepted actual arbitrary-function '+group+' seam theorem required')
            mixed.pulse.radius.post.selected.inlet.add_hashes(self.hashes,receipt['input_hashes'])
            self.hashes[name]=sha(name)
        density=json.loads((HERE/mixed.RECEIPT).read_bytes())
        self.density_theorem=density['independent_general_cumulative_mixed_equations']
        if not self.density_theorem['passed'] or self.density_theorem['independent_full_cumulative_mixed_equations']!=60:
            raise ValueError('Admitted independent general raw cumulative equations required')
        for name in (heat_seams.NAME,heat_seams.RECEIPT,Path(__file__).name,
                PREFIX+'current_postpulse_interfaces.py',PREFIX+'power_angular_C4_check.py',
                PREFIX+'steep_waiting_C4_check.py',PREFIX+'current_power_angular_source.py',
                PREFIX+'current_steep_waiting_source.py'):
            self.hashes[name]=sha(name)
        self.acceptance_loaded=False;self.call_trace=[]
        self.function_transfer=self.source_function_transfer();self.assert_graph()
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not all(receipt[k] for k in GATES) or any(receipt[k] for k in OPEN):
                raise ValueError('Current postpulse mixed seam receipt/scope differs')
            if receipt['source_family']!=self.family_record:raise ValueError('Current postpulse seam family differs')
            mixed.pulse.radius.post.selected.inlet.add_hashes(self.hashes,receipt['input_hashes'])
            self.acceptance_loaded=True

    def source_function_transfer(self):
        p=self.prefix;b=self.steep_binding
        prefix=('passed','first_five_scaled_angular_and_Gamma_rows_copied_before_appending_fifth',
            'physical_prefix_uses_same_exact_repair_scale','first_five_complete_future_rows_are_actual_C4_prefix',
            'copied_enclosures_do_not_select_point_coefficients')
        steep=('passed','current_waiting_root_and_heat_atoms_not_read_from_saved_packets',
            'current_absolute_pressure_and_nonzero_angular_histories_retained',
            'same_exact_transition_kernel_definitions_with_independent_directed_cells')
        if not all(p[k] for k in prefix) or not all(b[k] for k in steep) \
                or b['source_caps_used_as_defining_field_values'] or b['interval_overlap_used_as_join_proof']:
            raise ValueError('Same complete current future and original forward histories required')
        for key in ('exact_current_kernel_and_parameter_source_bindings',
                'actual_current_angular_s0_future_support_bindings'):
            if not b[key]['passed']:raise ValueError('Current angular/steep exact source binding missing: '+key)
        live=b['actual_current_live_angular_Gamma_history_bindings']
        if set(live)!={'terminal','source','heat','XR','PR'} or not all(live.values()):
            raise ValueError('Actual angular/Gamma forward history AST assignments required')
        native=self.post.power.bindings['current_native_parameter_source_bridge']
        if not native['passed'] or not all(native['actual_native_only_object_graph'].values()):
            raise ValueError('Same original native parameter/datum source required')
        if not all(self.flatten_theorem.values()):raise ValueError('Actual current flatten full future seam theorem missing')
        identities=dict(self.flatten_theorem)
        mappings=(('power_angular','power','new_provider_endpoint'),
            ('angular_steep','steep','angular_entry'),('entry_power','steep','entry_power'),
            ('power_exit','steep','power_exit'),('exit_waiting','steep','exit_waiting'))
        for seam,group,prefix_name in mappings:
            for label in ('theta','X','energy','pressure'):
                for k in range(5):
                    for j in range(5-k):
                        key=prefix_name+'_'+label+'_y%d_Z%d'%(k,j)
                        if not self.canonical[group].get(key):raise ValueError('Actual original seam derivative missing: '+key)
                        identities[seam+'_'+label+'_y%d_Z%d'%(k,j)]=True
        if len(identities)!=360:raise ValueError('Six current postpulse primitive mixed grids incomplete')
        return dict(passed=True,current_six_nonheat_primitive_mixed4_function_identities=identities,
            current_full_future_prefix_flags={key:p[key] for key in prefix},
            current_nonzero_history_and_kernel_flags={key:b[key] for key in steep},
            same_current_native_parameter_source=native,
            complete_flatten_future='EF + q^2 exp(-200mu)/4 * (I(2mu,Lrel) + exp(-2mu Lrel)*(post+signed_angular_change))',
            angular_steep_energy='post = Ein + exp(-1-mu)*(I(2,Ts)+exp(-2Ts)*(Eout+exp(-1-delta/2)*(I(delta,W)+exp(-delta W)*E0/(1-eps)^2)))',
            same_original_nonzero_X_pressure_and_quadratic_memories_retained=True,
            independent_original_P0_not_reset_or_tail_patched=True,
            raw_products_and_positive_radial_density_equations=self.density_theorem,
            accepted_closed_heat_function_joins=self.before.function_transfer,
            coordinate_rule='D_native=J*D_logR; J=(Lrel-4), Ts or true waiting W for phase charts, 1 otherwise',
            functions_identified_before_differentiation_and_interval_comparison=True,
            mixed_scope='ordinary logR/Z total order k+j<=4; no independent order8 grid or uniform norm')

    def assert_graph(self):
        s=self.owner.selected;p=self.post
        result=dict(accepted_same_current_closed_heat_and_mixed_owners=self.before.acceptance_loaded
                and all(self.before.assert_graph().values()) and self.transport.acceptance_loaded,
            same_actual_postpulse_dispatcher=self.post is self.transport.owner.post and all(p.assert_graph().values()),
            same_current_unique_repair=p.heat.repair is s.seed.exact.repair is self.closed.exact.repair,
            same_current_complete_future=p.outer.future is p.steep.future is p.heat.future is s.future,
            same_actual_flatten_and_independent_P0=p.outer.flatten is s.flatten
                and p.outer.flatten.inlet.datum is s.datum,
            same_source_context=self.ctx is p.ctx is self.closed.ctx,
            same_exact_radius_expression_graph=self.graph is self.owner.graph is self.radius.graph,
            same_accepted_pulse_and_Rv_interface_owners=self.before.before.before is self.transport
                and self.transport.before.acceptance_loaded,
            accepted_original_power_steep_function_certificates=all(all(v.values()) for v in self.canonical_flags.values()),
            current_complete_future_and_forward_memory_functions=self.function_transfer['passed'],
            no_repair_future_selection_or_legacy_physical_owner_replay=True)
        if not all(result.values()):raise ValueError('Current postpulse mixed seam graph differs: '+str(result))
        return result

    def assemble(self,name,Z,left,right):
        if name not in SEAMS:raise ValueError('One of the eight current raw postpulse seams required')
        self.assert_graph();lc,lv,rc,rv=SEAMS[name]
        ends=mixed.pulse.radius.post.selected.inlet.endpoints
        expected_labels={'y%d_Z%d'%(k,j) for k in range(5) for j in range(5-k)}
        expected_outputs={'Mz','Mtheta','Mtheta_z','Mztheta','Mp','P0','pressure','Utheta','Uz','Ur'}
        for view,chart,coordinate in ((left,lc,lv),(right,rc,rv)):
            if view['chart']!=chart or view['geometry']['exact_native_coordinate']!={'numerator':coordinate,'denominator':1}:
                raise ValueError('Same exact original postpulse seam route required')
            if ends(view['Z'])!=ends(self.ctx.mpf(Z)):raise ValueError('Same current axial source observation required')
            if set(view['log_radius_mixed_rows'])!=expected_outputs:raise ValueError('All ten original raw outputs required')
            for grid in view['log_radius_mixed_rows'].values():
                if set(grid)!=expected_labels or any(type(row) is not mixed.FactorizedMixedSourceRow
                        or row.coefficients.ctx is not self.ctx
                        or any(v.graph is not self.graph for _,v in row.log_scale_parts) for row in grid.values()):
                    raise ValueError('Complete actual same-graph typed mixed rows required')
            packet=view['original_forward_source_packet']
            if not chart.startswith('heat_'):
                if chart=='flatten':valid=packet['component_units']['velocity']=='U/Ev0; exact formal Ev0/Pstar kept separately'
                else:valid=packet['velocity_units']=='U/Ev0; exact formal Ev0/Pstar retained separately'
                if not valid:raise ValueError('Actual current postpulse velocity source in Ev0 units required')
        if name in heat_seams.SEAMS:
            self.before.assemble(name,Z,left,right)
        self.call_trace.append(dict(seam=name,actual_current_source_routes=(lc,rc),
            exact_whole_function_identity_precedes_interval_diagnostics=True))
        return dict(seam=name,Z=self.ctx.mpf(Z),left=left,right=right,
            same_current_original_pressure_and_complete_energy_memory=True,
            exact_phase_Jacobian_powers_retained_as_scale_functions=True,
            **dict.fromkeys(GATES,self.acceptance_loaded),**dict.fromkeys(OPEN,False))

    @source_precision
    def evaluate(self,name,Z):
        if name not in SEAMS:raise ValueError('Current original postpulse seam required')
        lc,lv,rc,rv=SEAMS[name]
        return self.assemble(name,Z,self.transport.evaluate(lc,Z,lv),self.transport.evaluate(rc,Z,rv))

    def interface_registry(self):
        """Fourteen original adjacent source interfaces, with explicit owners."""
        return {'Rv':dict(owner='accepted CurrentOriginalRpCommonUnitSeam',
                    call='evaluate_interface(Rv,Z)',scope='common Ev0 units and ordinary logR/Z total order4'),
            **{name:dict(owner='accepted CurrentOriginalRpPulseMixedSeams',
                    call='evaluate_interface('+name+',Z)',scope='exact formal pulse boundary and mixed logR/Z total order4')
                for name in heat_seams.pulse_seams.SEAMS},
            **{name:dict(owner='CurrentOriginalRpPostpulseMixedSeams',
                    call='evaluate_interface('+name+',Z)',scope='exact postpulse boundary and mixed logR/Z total order4')
                for name in SEAMS}}

    @source_precision
    def evaluate_interface(self,name,Z):
        """Dispatch to the accepted source owner; retain its original schema."""
        self.assert_graph()
        if name=='Rv':return self.transport.before.evaluate(Z)
        if name in heat_seams.pulse_seams.SEAMS:return self.before.before.evaluate(name,Z)
        if name in SEAMS:return self.evaluate(name,Z)
        raise ValueError('One of the fourteen original adjacent source interfaces required')


report=heat_seams.report


@source_precision
def run(before=None,observed_transport=None,observed_views=None,observed_heat_views=None):
    began=time.monotonic();owner=CurrentOriginalRpPostpulseMixedSeams(before,require_checked=False)
    if observed_transport is not None:
        if observed_transport is not owner.transport or set(observed_views)!=set(OBSERVATIONS) \
                or set(observed_heat_views)!=set(heat_seams.SEAMS):
            raise ValueError('Complete typed observations of the same actual accepted transport required')
        views={name:owner.assemble(name,'.521',observed_views[l],observed_views[r]) for name,(l,r) in PAIRS.items()}
        views.update({name:owner.assemble(name,'.521',view['left'],view['right'])
            for name,view in observed_heat_views.items()})
    else:views={name:owner.evaluate(name,'.521') for name in SEAMS}
    result=dict(source_family=owner.family_record,candidate_current_postpulse_mixed_seams_constructed=True,
        actual_source_graph=owner.assert_graph(),original_power_steep_function_theorems=owner.canonical,
        actual_current_flatten_full_future_function_theorem=owner.flatten_theorem,
        same_actual_C4_C5_complete_future_prefix_binding=owner.prefix,
        same_actual_current_angular_steep_forward_history_binding=owner.steep_binding,
        current_original_function_transfer=owner.function_transfer,
        fourteen_interface_source_chain_registry=owner.interface_registry(),
        accepted_original_function_certificate_flags=owner.canonical_flags,
        actual_eight_postpulse_mixed_seam_views={name:report(view) for name,view in views.items()},
        actual_source_call_trace=owner.call_trace,fresh_Z='.521',
        scope='Eight current original raw logR/Z function joins through total mixed order4; no uniform/global norm or resolved physical point/time field',
        **dict.fromkeys(GATES+OPEN,False),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began)
    data=json.dumps(mixed.pulse.raw.packed(result),separators=(',',':'))+'\n'
    (HERE/NAME).write_bytes(gzip.compress(data.encode(),mtime=0))
    return owner,views


if __name__=='__main__':run()
