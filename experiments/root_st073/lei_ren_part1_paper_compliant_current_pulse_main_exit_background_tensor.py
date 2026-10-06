"""Current full main/exit tensors with all incoming history cross terms.

Source algorithms and current bindings identify functions before enclosures.
The complete producer JSON is stored as deterministic gzip, without pruning.
"""
import ast
import copy
import gzip
import hashlib
import importlib
import json
import math

import mpmath as mp
from lei_ren_part1_paper_compliant_current_pulse_gap_background_tensor import (
    CurrentPulseGapBackgroundTensor,HERE,PREFIX,OPEN,TENSOR_KEYS,SourceAST,
    sha,pack,encode,endpoints,copy_jet,source_precision,accepted,_verify_hashes,
    IntervalTaylor,canonical_tensor_groups)
from lei_ren_part1_paper_compliant_pulse_main_exit_similarity_C4 import (
    main_exit_shapes,partial_linear_kernel,partial_future_energy)
from lei_ren_part1_paper_compliant_pulse_main_exit_physical_C2 import (
    compiled_main_exit_lift,main_exit_to_physical_packet)
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import gp_jets

NAME=PREFIX+'current_pulse_main_exit_background_tensor.json.gz'
RECEIPT=PREFIX+'current_pulse_main_exit_background_tensor_check.json'
GATES=('current_actual_main_exit_full_tensors_available',
    'current_actual_main_exit_full_meridional_decomposition_available',
    'current_actual_main_exit_completed_tensor_join_certified',
    'current_actual_exit_gap_completed_tensor_join_certified')
VIEWS={'whole_main':('pulse_main',(-1,1),('.02','10'),('-3','-1'),None,'1'),
    'whole_exit':('pulse_exit',(-1,1),(10,11),('-3','-1'),None,'1'),
    'main_left':('pulse_main',(-1,1),'.02','-1',None,'1'),
    'main_exit':('pulse_main',(-1,1),10,'-1',None,'1'),
    'exit_gap':('pulse_exit',(-1,1),11,'-1',None,'1'),
    'fresh_main':('pulse_main','.537','4.337','-2.6','.41','.8'),
    'fresh_main_left':('pulse_main','.731','.371','-2.6','.41','.8'),
    'fresh_exit':('pulse_exit','.731','10.537','-2.6','.41','.8')}
SEAMS=('main_exit','exit_gap')
LEGACY_STATEMENT_AST_SHA256={
    ('pulse_main_exit_similarity_C4','source_proof'):(
        'b3d83c33321e054fa156c6168d087380fb842c6a0310baf95861ef346633a6e8',
        '19dcb6fc7aed5ed8d6058d169aeceb7366b2168bf840ae2eff25ffe74f4768fc',
        '427d396a69bc202b1b7df936083fe2c2bbb813a2cab1c3801396bdcba7c08a6a',
        'b3ce9160e22301b085ad3652b1611a8429da3e001d1d3a7a27a2a318262ad7c2',
        'd6abb5ffaac1454fd470438b939ce7aa7fcf378325920144346fb0ee3a83f77f'),
    ('pulse_main_exit_physical_C2','source_and_join_binding'):(
        'b3d83c33321e054fa156c6168d087380fb842c6a0310baf95861ef346633a6e8',
        '755c40a687bfa6c6117eadd499f4af5d9acd609da1b9976022150a61d61d1dfc',
        '2d3ba936fb53725d1ae984fc9a13ec41e4d684cd10f567290741b2d5baa0af08',
        '79734fee33f38c57a6ad881902c5d9716693f8af256fa319b121784ea1386667',
        'd9839671f6a9c19f86a8215b37ce04c3c4c70692cd4cdbd44415a4add43d7ad7',
        '172d6e8b0c595d57c1ed6cefe66187a5178c68f1426917b6b3d95b40f1ed58f5',
        '97e5d7278f315aa438e3b6d21c0d6a7dbb39a32d59d853870f78f9ab04235a30',
        '6b0faaf2377bd931f5f57aad1a9e72d7492eab5f0c0930311956e0949cd23801',
        '56efb4e35ca46b03ee8c6f0c35ea0b2ca16af9a0e4a292bd34d7eedf17d9728c')}


def original_generic_proof(stem,name,operator):
    """Replay mathematical operators; current admission replaces old receipts.

    Only top-level legacy record checks are removed. The original symbolic
    formulas, source AST checks, kernel changes of variable and tensor row
    identities run unchanged. Removed checks are reported explicitly.
    """
    module=importlib.import_module(PREFIX+stem);asts=SourceAST()
    fn=copy.deepcopy(asts.method(stem,name));fn.decorator_list=[]
    legacy={'records','certificate','flat','pressure'};kept=[];removed=[];digests=[]
    for node in fn.body:
        names={v.id for v in ast.walk(node) if isinstance(v,ast.Name)}
        if names.intersection(legacy):
            removed.append(ast.unparse(node).splitlines()[0])
            digests.append(hashlib.sha256(ast.dump(node,include_attributes=False).encode('utf8')).hexdigest())
            continue
        kept.append(node)
    if tuple(digests)!=LEGACY_STATEMENT_AST_SHA256.get((stem,name)):
        raise ValueError('Legacy receipt removal differs from exact reviewed statement allowlist')
    fn.body=kept;env=dict(vars(module));env['operator']=operator
    exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),
        '<original generic main/exit theorem, current admission separate>','exec'),env)
    result=env[name]({})
    if not all(result['identities'].values()):raise ValueError('Original generic main/exit source identity failed')
    result.update(input_hashes={**result['input_hashes'],**asts.hashes},
        omitted_legacy_receipt_checks=removed,omitted_legacy_statement_AST_sha256=digests,
        exact_reviewed_legacy_statement_allowlist_enforced=True,
        scope='Generic original arbitrary-history operator theorem; current source admission is supplied separately.')
    return result


def current_main_exit_source_theorem(field):
    graph=field.history.selected.assert_graph();current=field.atlas.pulse_interfaces.proof
    physical_graph=field.physical.assert_graph()
    native=importlib.import_module(PREFIX+'axial_pulse_field')
    flat=importlib.import_module(PREFIX+'flat_pulse_derivatives')
    high=importlib.import_module(PREFIX+'pulse_high_jets')
    fixed_gp=dict(current_fixed_flat_provider_consumed=physical_graph['original_fixed_flat_gp_beta_derivative_provider'],
        forward_kernel_uses_native_fixed_gp=partial_linear_kernel.__globals__['gp'] is native.gp,
        remaining_energy_uses_native_fixed_gp=partial_future_energy.__globals__['gp'] is native.gp,
        main_shapes_use_current_fixed_gp_jets=main_exit_shapes.__globals__['gp_jets'] is flat.gp_jets,
        current_pulse_main_uses_original_native_dispatch=(field.pulse.main.__func__ is high.CompliantPulseHighJets.main
            and isinstance(field.pulse,native.CompliantAxialPulseField)),
        native_main_uses_same_fixed_gp=native.CompliantAxialPulseField.main.__globals__['gp'] is native.gp,
        current_pulse_mixed_uses_same_fixed_gp_jets=field.pulse._high_packet.__func__.__globals__['gp_jets'] is flat.gp_jets,
        flat_derivative_original_gp_alias_is_same=flat.original_gp is native.gp)
    if not all(fixed_gp.values()):raise ValueError('Current fixed GP source/provider differs from main/exit integrals')
    mixed=current['source_bound_five_current_pulse_mixed4_theorem']
    if not current['passed'] or not mixed['passed']:raise ValueError('Current selected whole pulse function trace theorem required')
    operator=field.gap_tensor.end_tensor.physical_proof
    source=original_generic_proof('pulse_main_exit_similarity_C4','source_proof',operator)
    physical=original_generic_proof('pulse_main_exit_physical_C2','source_and_join_binding',operator)
    asts=SourceAST()
    asts.method('axial_pulse_field','gp');asts.method('flat_pulse_derivatives','gp_jets');asts.method('pulse_high_jets','main')
    asts.expression('pulse_mixed_C4','_high_packet','shape',wanted="gp_jets(c,coord['xi'])")
    asts.expression('current_pulse_main_exit_background_tensor','_data','ap',wanted="copy_jet(c,selected['selected_ap_Taylor'])")
    asts.expression('current_pulse_main_exit_background_tensor','_data','incoming',wanted="[copy_jet(c,value) for value in selected['incoming']['moment_Taylor']]")
    asts.expression('current_pulse_main_exit_background_tensor','_data','future',wanted="copy_jet(c,self.pulse.selection.future.future(Z)['complete_future_energy_Taylor'])/2")
    asts.expression('current_pulse_main_exit_background_tensor','_data','right',wanted='self.gap_tensor.cache[key]')
    asts.expression('current_pulse_main_exit_background_tensor','chart','kernels',wanted="[partial_linear_kernel(c,mu,xi,c.mpf('.5')-i*mu,self.cells,self.window) for i in (1,2)]")
    asts.expression('current_pulse_main_exit_background_tensor','chart','remaining',wanted='partial_future_energy(c,xi,self.cells)')
    if not all(field.gap_tensor.proof['current_exact_velocity_and_remainder_scale_identities'].values()):raise ValueError('Current gap normalization/scale source bridge required')
    return dict(original_generic_full_main_exit_stress_theorem=source,
        original_generic_full_main_exit_physical_and_exit_gap_theorem=physical,
        current_selected_graph=graph,current_fixed_gp_provider_bindings=fixed_gp,
        current_selected_parameter_and_five_history_function_trace_theorem=current,
        current_gap_complete_pressure_and_weight_source_theorem=field.gap_tensor.proof,
        original_full_meridional_paper_stress_AST_theorem=field.gap_tensor.end_tensor.formulas,
        original_full_physical_operator_theorem=operator,unchanged_compiled_main_exit_physical_lift=field.lift_proof,
        current_C5_selected_ap_incoming_controls_future_algorithms_bound=True,
        current_partial_forward_and_remaining_energy_defining_integrals_bound=True,
        current_reduced_absolute_Rv_pressure_and_native_Pin_retained=True,
        current_selected_two_row_inverse_and_energy_equation_instantiates_exit_gap_theorem=True,
        original_main_exit_xi10_same_global_source_function=True,
        original_exit_gap_xi11_uses_current_exact_selection_before_enclosure=True,
        original_kernel_positive_omitted_tail_not_deleted=True,
        all_local_incoming_and_square_meridional_cross_terms_retained=True,
        kernel_cells_only_enclose_defining_functions_not_point_values=True,
        source_function_equality_precedes_common_bounds=True,
        input_hashes={**source['input_hashes'],**physical['input_hashes'],**asts.hashes},passed=True)


def read_producer():
    return json.loads(gzip.decompress((HERE/NAME).read_bytes()))


class CurrentPulseMainExitBackgroundTensor:
    @source_precision
    def __init__(self,gap_tensor=None,require_checked=True,cells=64,window=1000):
        if not isinstance(cells,int) or cells<1 or not isinstance(window,int) or window<1:raise ValueError('Positive directed kernel cells/window required')
        self.gap_tensor=gap_tensor if gap_tensor is not None else CurrentPulseGapBackgroundTensor()
        if not self.gap_tensor.acceptance_loaded:raise ValueError('Checked actual current gap tensor required')
        self.physical=self.gap_tensor.physical;self.history=self.physical.history;self.pulse=self.history.selected.pulse
        self.ctx=self.physical.ctx;self.atlas=self.gap_tensor.atlas
        self.family=self.gap_tensor.family;self.source=self.gap_tensor.source;self.datum_sha=self.gap_tensor.datum_sha
        self.cells=cells;self.window=window;self.assert_graph()
        self.lift,self.lift_proof=compiled_main_exit_lift();self.proof=current_main_exit_source_theorem(self)
        self.cache={};self.hashes=dict(self.gap_tensor.hashes);self.hashes.update(self.proof['input_hashes']);self.hashes.update(self.lift_proof['input_hashes'])
        for stem in ('pulse_main_exit_similarity_C4','pulse_main_exit_physical_C2','current_pulse_main_exit_background_tensor'):
            name=PREFIX+stem+'.py';self.hashes[name]=sha(name)
        self.acceptance_loaded=False
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source,GATES[0]);_verify_hashes(receipt)
            if receipt['datum_enclosure_sha256']!=self.datum_sha or not all(receipt[k] for k in GATES) or any(receipt[k] for k in OPEN):raise ValueError('Current main/exit receipt exceeds scope')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.acceptance_loaded=True

    def assert_graph(self):
        self.gap_tensor.assert_graph()
        if not (self.gap_tensor.acceptance_loaded and self.physical is self.gap_tensor.physical and self.history is self.physical.history
                and self.pulse is self.history.selected.pulse and self.pulse is self.physical.pulse and self.pulse is self.gap_tensor.pulse
                and self.atlas is self.gap_tensor.atlas and self.ctx is self.physical.ctx):
            raise ValueError('Main/exit must share checked current selected pulse, pressure and full histories')

    def _data(self,Z):
        c=self.ctx;key=endpoints(Z)
        if key not in self.cache:
            selected,inlet,u,_,_=self.pulse.data(Z)
            ap=copy_jet(c,selected['selected_ap_Taylor'])
            incoming=[copy_jet(c,value) for value in selected['incoming']['moment_Taylor']]
            energy=copy_jet(c,selected['incoming']['energy_Taylor'])
            future=copy_jet(c,self.pulse.selection.future.future(Z)['complete_future_energy_Taylor'])/2
            if key not in self.gap_tensor.cache:self.gap_tensor.chart('pulse_gap',Z,11)
            right=self.gap_tensor.cache[key]
            if any(value.order!=5 for value in [ap,energy,future]+incoming+[right['J0'],right['P0']]):raise ValueError('Current full selected C5 histories required')
            self.cache[key]=dict(ap=ap,incoming=incoming,incoming_energy=energy,future=future,J0=right['J0'],P0=right['P0'],
                controls=right['controls'],inlet=inlet,u=u,current_selected_result=selected)
        return self.cache[key]

    @source_precision
    def chart(self,chart,Z,coordinate,log_tau='-1',theta=None,viscosity='1'):
        self.assert_graph();c=self.ctx;Z=c.mpf(Z);xi=c.mpf(coordinate);lt=c.mpf(log_tau);nu=c.mpf(viscosity)
        if chart not in ('pulse_main','pulse_exit'):raise ValueError('Current main/exit tensor chart required')
        if not all(mp.isfinite(x) for x in endpoints(Z)+endpoints(xi)+endpoints(lt)+endpoints(nu)) or endpoints(Z)[0]<-1 or endpoints(Z)[1]>1 or endpoints(nu)[0]<=0:raise ValueError('Z[-1,1],finite coordinate/logtau and nu>0 required')
        limits=(c.mpf('.02'),c.mpf(10)) if chart=='pulse_main' else (c.mpf(10),c.mpf(11))
        if endpoints(xi)[0]<endpoints(limits[0])[0] or endpoints(xi)[1]>endpoints(limits[1])[1]:raise ValueError('Original main xi[.02,10] or exit xi[10,11] required')
        if theta is not None and not all(mp.isfinite(x) for x in endpoints(c.mpf(theta))):raise ValueError('Finite theta or all-angle bounds required')
        data=self._data(Z);mu=c.mpf(self.pulse.mu);delta=c.mpf(self.pulse.delta)
        z=IntervalTaylor.variable(c,Z,5);C=IntervalTaylor(c,[1+Z**2,2*Z,1,0,0,0]).reciprocal()
        kernels=[partial_linear_kernel(c,mu,xi,c.mpf('.5')-i*mu,self.cells,self.window) for i in (1,2)]
        remaining=partial_future_energy(c,xi,self.cells)
        rows=main_exit_shapes(c,mu,delta,z,C,c.mpf(self.pulse.Xp),data['ap'],data['incoming'],
            [v['enclosure'] for v in kernels],remaining,data['future'],data['J0'],data['P0'],xi)
        logR=self.physical.logRp+xi/mu
        logB=dict(logPstar=self.physical.logP,actual_log_inlet_U=c.ln(c.mpf(self.gap_tensor.end_tensor.flatten_power.flatten.U)),
            inverse_mu=-xi/(2*mu),finite=-xi)
        logH=-(1-mu)*xi/mu;logD0=c.mpf(self.pulse.logE)
        extras=dict(one=c.mpf(0),incoming1=-(c.mpf('.5')-mu)*xi/mu,incoming2=-(c.mpf('.5')-2*mu)*xi/mu,
            Q=-(1+2*mu)*(13-xi)/mu,end_square=2*logD0)
        def logs(rp,bp,hp=0,extra='one'):
            return dict(source_logR=rp*logR,**{k:bp*value for k,value in logB.items()},signed_original_memory_log=hp*logH,
                exact_source_history_log=extras[extra],normalization=-c.ln(2)/2)
        sectors={}
        for label,parts in rows['stress'].items():
            sectors[label]={}
            for name,part in parts.items():
                rp,bp,dp,hp=part['mode']
                grid={'s%d_Z%d'%(j,n):jet[n]*math.factorial(n) for j,jet in enumerate(part['full_derivative_rows']) for n in range(4-j)}
                sectors[label][name]=dict(mode=part['mode'],extra_source=part['extra_source'],
                    exact_source_log_parts=logs(rp,bp,hp,part['extra_source']),full_stress_mixed3_coefficient_enclosures=grid)
        packet=dict(Z=Z,xi=xi,full_meridional_stress_log_sectors=sectors,exact_source_logs=dict(R=logR,B=logB,H=logH,D0=logD0,extra=extras))
        mapped=main_exit_to_physical_packet(packet,mu)
        velocity=dict(local=rows['velocity_local'],incoming=rows['velocity_incoming'])
        point=self.lift(c,mapped,delta,velocity,lt,theta,nu)
        native_pressure=self.pulse.pressure_moment(Z,xi/self.pulse.mu,data['inlet'],data['u'])
        return dict(point,chart=chart,coordinate=xi,ordinary_radial_derivative='d_logR=d_y=mu*d_xi',
            current_actual_source_stress_packet=mapped,source_full_main_exit_history_rows=rows,
            current_selected_ap=data['ap'],current_selected_controls=data['controls'],current_incoming_moments=data['incoming'],current_incoming_energy=data['incoming_energy'],
            current_complete_C5_future=data['future'],current_full_selected_end_loss=data['J0'],current_signed_absolute_Rv_pressure=data['P0'],
            original_partial_linear_kernel_bounds=kernels,original_partial_future_energy_bound=remaining,
            original_native_forward_absolute_pressure=native_pressure,current_source_three_component_velocity_rows=velocity,
            current_selected_complete_history_graph_retained=True,actual_full_stress_not_local_difference=True,
            all_local_incoming_cross_and_square_radial_terms_retained=True,
            pressure_remaining_reduction_before_enclosure=True,positive_omitted_forward_tails_retained=True,
            source_bounds_not_resolved_physical_point_values=True,**dict.fromkeys(OPEN,False))

    @source_precision
    def interface(self,name,Z=(-1,1),log_tau=('-3','-1'),theta=None,viscosity='1'):
        if name not in SEAMS:raise ValueError('Current main_exit or exit_gap tensor seam required')
        if name=='main_exit':left=self.chart('pulse_main',Z,10,log_tau,theta,viscosity);right=self.chart('pulse_exit',Z,10,log_tau,theta,viscosity)
        else:left=self.chart('pulse_exit',Z,11,log_tau,theta,viscosity);right=self.gap_tensor.chart('pulse_gap',Z,11,log_tau,theta,viscosity)
        a=canonical_tensor_groups(left);b=canonical_tensor_groups(right);common={}
        if set(a)!=set(b):raise ValueError('Full actual current tensor trace layouts differ')
        for key in a:
            values=[endpoints(row['log_absolute_upper'])[1] for row in a[key]+b[key] if not row['exact_zero']]
            common[key]=dict(exact_zero=not values,log_absolute_upper=self.ctx.mpf(max(values))+self.ctx.ln(len(values)) if values else None)
        return dict(seam=name,common_actual_tensor_rows=common,current_common_tensor_contribution_count=len(common),
            current_actual_forward_selected_full_tensor_join_theorem=self.proof,
            source_function_equality_precedes_common_triangle_bounds=True,actual_nonzero_shared_boundary_histories_preserved=True,
            interval_overlap_not_used_as_function_identity=True,**dict.fromkeys(OPEN,False))

    def manifest(self):
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,datum_enclosure_sha256=self.datum_sha,
            current_main_exit_full_history_stress_physical_and_two_join_theorem=self.proof,
            actual_current_tensor_regions_available=['pulse_main','pulse_exit']+self.gap_tensor.manifest()['actual_current_tensor_regions_available'],
            actual_current_completed_tensor_adjacent_interface_count=13,actual_current_completed_tensor_internal_interface_count=4,
            current22_velocity_pressure_interface_inventory=dict(adjacent=14,internal=8),
            source_domain='main xi[.02,10],exit xi[10,11],ordinary logR rows;R>0,|Z|<1,tau>0,constant nu>0,compact logtau',
            directed_partial_kernel_cells=self.cells,directed_partial_kernel_window=self.window,
            producer_storage='Deterministic gzip of complete UTF-8 JSON, no tensor/velocity/history entries removed',
            all_original_full_meridional_cross_terms_and_omitted_forward_tails_retained=True,
            remaining_entrance_core_axis_angular_internal_global_and_temporal_scope_open=True,
            input_hashes=self.hashes,**dict.fromkeys(GATES+OPEN,False))


@source_precision
def run(field=None):
    field=field if field is not None else CurrentPulseMainExitBackgroundTensor(require_checked=False)
    result=field.manifest();result['current_actual_main_exit_tensor_views']={}
    for name,args in VIEWS.items():
        result['current_actual_main_exit_tensor_views'][name]=field.chart(*args)
        print('Current actual main/exit tensor: '+name,flush=True)
    result['current_actual_two_main_exit_tensor_interfaces']={name:field.interface(name) for name in SEAMS}
    payload=(json.dumps(encode(pack(result)),separators=(',',':'))+'\n').encode('utf8')
    (HERE/NAME).write_bytes(gzip.compress(payload,compresslevel=6,mtime=0))
    return result


if __name__=='__main__':run()
