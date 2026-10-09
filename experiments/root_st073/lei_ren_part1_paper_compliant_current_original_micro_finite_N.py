"""Source-bound Section 11 inlet and genuine finite-N prefix through Rm.

The original background before r_minus is unchanged. Changed histories
start with that background's own histories and integrate signed density
increments. Zero initial defect is a source/support definition, not a
background micro increment or a default numerical zero vector.
"""
import ast
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_micro_functions as micro

macro=micro.downstream;switch=macro.downstream
fields,parameters,phase,primitives,bounds=macro.fields,macro.parameters,macro.phase,macro.primitives,macro.bounds
HERE,PREFIX,sha=macro.HERE,macro.PREFIX,macro.sha
NAME=PREFIX+'current_original_micro_finite_N.json.gz'
RECEIPT=PREFIX+'current_original_micro_finite_N_check.json'
GATE='current_original_source_bound_micro_finite_N_prefix_through_Rm_enclosed'
RATES,C0,Z=macro.RATES,macro.C0,macro.Z


def serialized(value):
    if isinstance(value,dict):return {('y%d_Z%d'%key if isinstance(key,tuple) else key):serialized(row) for key,row in value.items()}
    if isinstance(value,(tuple,list)):return [serialized(row) for row in value]
    return micro.serialized(value)


def source_bindings():
    path=HERE/(PREFIX+'current_generic_five_defect_bounds.py')
    tree=ast.parse(path.read_text(encoding='utf8'))
    fn=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='integral_contracts')
    wanted=dict(inlet_defect_exact_zero=True,actual_original_inlet_history_unchanged=True,
        zero_extension_before_y0=True,unchanged_original_axis_pressure_P0=True)
    for key,value in wanted.items():
        got=[n.value.value for n in ast.walk(fn) if isinstance(n,ast.keyword) and n.arg==key and isinstance(n.value,ast.Constant)]
        if got!=[value]:raise ValueError('Original Section11 initial/support contract changed: '+key)
    graph=ast.parse((HERE/(PREFIX+'current_generic_loop_function_sources.py')).read_text(encoding='utf8'))
    flat=next(n for n in ast.walk(graph) if isinstance(n,ast.FunctionDef) and n.name=='flat')
    for key,value in (('flat_when_Delta_ge_eta',True),('active_body_not_evaluated_on_flat_branch',True)):
        got=[n.value.value for n in ast.walk(flat) if isinstance(n,ast.keyword) and n.arg==key and isinstance(n.value,ast.Constant)]
        if got!=[value]:raise ValueError('Original exact flat-loop branch changed: '+key)
    native=ast.parse((HERE/(PREFIX+'current_native_generic_left_inlet.py')).read_text(encoding='utf8'))
    left=next(n for n in ast.walk(native) if isinstance(n,ast.FunctionDef) and n.name=='left_inlet')
    initial=dict(left_loop_q_A_B_exact_zero_by_existing_checked_collar=True,
        generic_inlet_defect_zero_is_declared_Section11_initial_condition=True,
        original_P0_and_all_five_incoming_histories_retained=True)
    for key,value in initial.items():
        got=[n.value.value for n in ast.walk(left) if isinstance(n,ast.keyword) and n.arg==key and isinstance(n.value,ast.Constant)]
        if got!=[value]:raise ValueError('Actual native Section11 initial-condition source changed: '+key)
    assignments=(('current_original_bridge_macro_functions.py','owner',"inlet=frame['second_exit']"),
        ('current_original_bridge_macro_moments.py','owner',"packet=self.fields.records['actual_bridge_integrals']['packets'][label]['second_exit']"),
        ('current_original_bridge_macro_moments.py','evaluate',"source_polys=dict(H=phi.scale(2),M=V,K=(phi*V).scale(2),A=V*V,B=phi*phi,C=phi*phi)"),
        ('current_original_micro_functions.py','histories',"rhs=dict(H=f.scale(phi,2),M=V,K=f.scale(f.multiply(phi,V),2),A=f.multiply(V,V),B=f.multiply(phi,phi),C=f.multiply(phi,phi))"))
    for name,method,text in assignments:
        origin=ast.parse((HERE/(PREFIX+name)).read_text(encoding='utf8'));target=ast.dump(ast.parse(text).body[0])
        matches=[n for fn in ast.walk(origin) if isinstance(fn,ast.FunctionDef) and fn.name==method
            for n in ast.walk(fn) if isinstance(n,ast.Assign) and ast.dump(n)==target]
        if len(matches)!=1:raise ValueError('Actual micro/macro source-definition join changed: '+name)
    own=ast.parse(Path(__file__).read_text(encoding='utf8'))
    select=next(n for n in ast.walk(own) if isinstance(n,ast.FunctionDef) and n.name=='modification_support')
    selector="chart=='first_micro' and hi<=micro.ep(sc/2)[0]"
    if not any(isinstance(n,ast.If) and ast.dump(n.test)==ast.dump(ast.parse(selector,mode='eval').body) for n in ast.walk(select)):
        raise ValueError('Executable current unchanged-left velocity selector changed')
    return dict(original_Section11_initial_and_support_bindings=wanted,
        accepted_native_positive_left_initial_condition_bindings=initial,
        original_micro_macro_inlet_and_six_history_source_assignments=assignments,
        original_exact_q_A_B_flat_branch_bound=True,
        current_coupled_micro_source=micro.source_bindings(),current_macro_source=macro.source_bindings(),
        source_owned_modified_velocity='background for R<=r_minus; E_N=E*exp(A/N), V_N=V+B/N for R>r_minus',
        source_owned_changed_history='background_history_j + integral_0^y exp(-rate_j*(y-t))*delta_density_j(t,Z,frac(N*t)) dt',
        actual_micro_global_phase='frac(N*hb*(s-s_c/2))',
        actual_micro_macro_seam='s=2 and rho=0 are the same R0=Ra*exp(2hb), same phi/V and six original own histories',
        source_equality_by_original_coupled_integral_definition_and_linear_history_ODE_uniqueness=True,
        executable_source_selector='OriginalMicroFiniteN.modification_support; query evaluates unchanged-left/flat-collar or original finite-N right branch',
        zero_inlet_not_inferred_from_Ra_background_increment_or_constructor_default=True)


def weights(series,left,right,end,rate):
    """True microscopic masses with compact source endpoints and offsets."""
    f,c=series.flow,series.c;l,r,e=map(c.mpf,(left,right,end))
    ep=micro.ep
    if any(ep(v)[0]!=ep(v)[1] for v in (l,r,e)) or not 0<=ep(l)[0]<ep(r)[0]<=ep(e)[0]<=2:
        raise ValueError('Ordered exact compact source endpoints required')
    ds=r-l;width=f.h*ds
    if not rate:return width,f.scalar(1),f.scalar(1)
    lam=c.mpf(rate.numerator)/rate.denominator
    # Keep the positive left offset in its own exponential factor rather
    # than losing it by subtracting an astronomical singleton from one.
    decay=series.scalar_series(-lam,r)[0]*series.scalar_series(lam,l)[0]
    tail=series.scalar_series(-lam,e)[0]*series.scalar_series(lam,r)[0]
    average,error=series.scalar_series(-lam,ds,average=True)
    mass=width*average
    return mass,decay,tail


def guard_inlet(family,label,N,P0,width,sc,packet,identity,eta_log,dstar_log):
    switch.guard_saved_source(family,label,N,P0,packet)
    if packet['source_identity_fields']!=identity:raise ValueError('Same actual core family, implicit source and pressure datum required')
    if packet['selected_eta_log_tuple']!=macro.base.encoded(eta_log) or packet['selected_dstar_log_tuple']!=macro.base.encoded(dstar_log):
        raise ValueError('Same exact selected original eta and d_star required')
    if packet['provider_sha256']!=sha(Path(micro.__file__).name):raise ValueError('Same checked current micro source provider required')
    if packet['source_width_log_tuple']!=macro.base.encoded(width):raise ValueError('Same exact source width tuple required')
    if packet['selected_positive_s_c_tuple']!=macro.base.encoded(sc):raise ValueError('Same exact positive phase-origin singleton required')
    if packet['source_owned_zero_defect_definition']!='original background unchanged before r_minus, signed Duhamel lower bound0':
        raise ValueError('Actual current Section11 initial/support definition required')
    if not packet['flat_left_collar_certified_for_same_current_source']:
        raise ValueError('Same-source smooth unchanged/modified velocity join required')
    if set(packet['actual_finite_N_inlet_defect_C0_Z'])!=set(RATES):raise ValueError('All five current initial defects required')
    for pair in packet['actual_finite_N_inlet_defect_C0_Z'].values():
        if len(pair)!=2 or any(not row['exact_zero'] for row in pair):raise ValueError('Source-defined initial value and Z defect must both be exact zero')


class OriginalMicroFiniteN:
    mode='current_original_source_owned_finite_N_rminus_R0_R100_R110_Rm_prefix_enclosures'
    def __init__(self,dps=500,require_checked=True):
        self.micro=micro.OriginalMicroFunctions(dps);self.macro=self.micro.downstream
        self.switch=self.macro.downstream;self.tail=self.switch.downstream
        self.parameters=self.tail.parameters;self.c=self.micro.c;self.family=self.micro.family
        self.hashes=dict(self.micro.hashes);self.cache={};self.inlets={};self.bindings=source_bindings()
        self.source_family=dict(actual_five_defect_family_sha256=self.family,
            implicit_source_sha256=self.micro.origin.source,datum_enclosure_sha256=self.micro.origin.datum)
        self.records={}
        for stem,gate in (
            ('current_generic_shear_loop_domain','current_original_generic_loop_two_sided_collars_and_repair_geometry_certified'),
            ('current_generic_shear_uniform_inputs','current_original_whole_generic_input_margins_and_log_scales_certified'),
            ('current_generic_loop_function_sources','current_original_generic_loop_velocity_and_five_increment_function_graphs_defined'),
            ('current_generic_five_defect_bounds','current_whole_generic_loop_five_defect_Duhamel_log_majorants_certified'),
            ('current_native_generic_left_inlet','current_original_generic_left_inlet_live_source_cover_query_executed')):
            checked_name=PREFIX+stem+'_check.json';checked=json.loads((HERE/checked_name).read_bytes())
            if not checked.get('all_passed') or not checked.get(gate) or checked['source_family']!=self.source_family:
                raise ValueError('Accepted same-family initial/support prerequisite required: '+stem)
            for name,digest in checked['input_hashes'].items():fields.previous.bind(self.hashes,name,digest)
            name=PREFIX+stem+'.json';self.records[stem]=json.loads((HERE/name).read_bytes())
            for path in (name,checked_name,PREFIX+stem+'.py'):fields.previous.bind(self.hashes,path,sha(path))
        self.saved_macro=json.loads(gzip.decompress((HERE/macro.NAME).read_bytes()))
        self.saved_switch=json.loads(gzip.decompress((HERE/switch.NAME).read_bytes()))
        for report,gate in ((self.saved_macro,macro.GATE),(self.saved_switch,switch.GATE)):
            if not report[gate] or report['source_family']!=self.family:raise ValueError('Same accepted current downstream driver required')
        fields.previous.bind(self.hashes,Path(__file__).name,sha(Path(__file__).name))
        if require_checked:
            checked=json.loads((HERE/RECEIPT).read_bytes())
            if not checked.get('all_passed') or checked['source_family']!=self.family:raise ValueError('Accepted current finite-N prefix receipt required')
            for name,digest in checked['input_hashes'].items():fields.previous.bind(self.hashes,name,digest)
            fields.previous.bind(self.hashes,RECEIPT,sha(RECEIPT))

    def inlet(self,label,N=257):
        N=phase.candidate_N(N);key=(label,N)
        if key in self.inlets:return self.inlets[key]
        if N!=self.saved_macro['candidate_N'] or N!=self.saved_switch['candidate_N']:
            raise ValueError('Same actual current finite N throughout prefix and downstream required')
        op=self.micro.owner(label);f,c=op.flow,op.c;sc=self.macro.sc
        ref=self.tail.reference.owner(label);_,_,axis=self.switch.r100.owner(label)
        if switch.canonical_expression(axis['pressure_axis_over_Pstar_squared_axial5'])!=switch.canonical_expression(ref.P0):
            raise ValueError('Same original live analytic P0 required')
        if op.series.flow is not f or ref.flow is not f:raise ValueError('One current live source context, basis and ledger required')
        if f.logs[0]._mpi_!=fields.previous.read_interval(c,self.macro.collar['exact_source_width_log_enclosure'])._mpi_:
            raise ValueError('Exact current micro/collar source width required')
        domain=self.records['current_generic_shear_loop_domain']['current_original_generic_loop_domain']
        scales=self.records['current_generic_shear_uniform_inputs']['current_actual_logarithmic_loop_scales']
        eta=fields.previous.read_interval(c,scales['selected_positive_eta_log'])
        dstar=fields.previous.read_interval(c,scales['logarithmic_selected_positive_lower_constants']['d_star'])
        excess=fields.previous.read_interval(c,domain['left_collar_kappa_minus2_lower'])
        limit=fields.previous.read_interval(c,domain['allowed_eta_log_upper_for_automatic_constant_edges'])
        if eta._mpi_!=self.parameters.eta_log._mpi_ or dstar._mpi_!=self.parameters.dstar_log._mpi_ or micro.ep(excess)[0]<=0:
            raise ValueError('Same selected original eta/d_star and positive strict left excess required')
        if micro.ep(eta)[1]>micro.ep(c.ln(excess)-c.ln(2))[0] or micro.ep(eta)[1]>micro.ep(limit)[0]:
            raise ValueError('Actual eta must satisfy the original strict flat left-edge collar')
        point=sc/2
        if micro.ep(point)[0]!=micro.ep(point)[1] or not 0<micro.ep(point)[0]<1:
            raise ValueError('Exact positive compact r_minus source endpoint required')
        with mp.workdps(c.dps+40):
            left=self.query(label,'first_micro',point,point,N)
            back=left['actual_micro_background_source'];raw=left['original_raw_source']
            recovered=left['original_generic_source']
            if not left['source_owned_modification_support']['increments_exact_zero']:
                raise ValueError('Source-owned unchanged-left selector required at the true inlet')
            if any(not row.zero for part in ('kernels','Z_derivatives') for row in left['original_signed_five_density_C0_Z'][part].values()):
                raise ValueError('Actual inlet velocity and density increments must be exactly zero')
        zero={name:[f.scalar(0),f.scalar(0)] for name in RATES}
        proof=dict(original_velocity_unchanged_for_R_le_rminus_by_source_definition=True,
            original_q_A_B_and_all_derivatives_zero_on_two_sided_left_collar=True,
            source_identity='kappa-2>=left_excess>=2eta implies original q=A=B=0',
            two_sided_collar_phase_fraction=('1/4','3/4'),actual_left_center_phase=point,
            original_background_histories_used_as_changed_incoming_histories=True,
            signed_defect_is_zero_at_its_defining_Duhamel_lower_bound=True,
            executable_piecewise_velocity_source_query_bound_to_actual_current_owner=True,
            background_micro_angular_axial_increment_need_not_be_zero=True,
            zero_initial_Z_rows_follow_from_source_function_identity_not_sampled_smallness=True,
            full_tensor_owner_interface_or_global_N_gate_not_promoted=True)
        packet=dict(source_family=self.family,source_identity_fields=self.source_family,source_frame=label,candidate_N=N,
            exact_common_P0_axial5=ref.P0,provider_sha256=sha(Path(micro.__file__).name),
            source_width_log_tuple=macro.base.encoded(f.logs[0]),selected_positive_s_c_tuple=macro.base.encoded(sc),
            selected_eta_log_tuple=macro.base.encoded(self.parameters.eta_log),selected_dstar_log_tuple=macro.base.encoded(self.parameters.dstar_log),
            actual_left_source_phase=point,exact_left_radius='Ra*exp(hb*s_c/2)',
            retained_positive_left_log_offset=f.h*point,
            source_owned_zero_defect_definition='original background unchanged before r_minus, signed Duhamel lower bound0',
            flat_left_collar_certified_for_same_current_source=True,original_left_flat_join_proof=proof,
            actual_current_left_background=back,original_current_left_raw_source=raw,
            original_current_left_common_source=recovered,actual_finite_N_inlet_defect_C0_Z=zero,
            original_background_inlet_is_not_zeroed=True,**dict.fromkeys(fields.previous.OPEN,False))
        saved=macro.base.encoded(serialized(packet))
        guard_inlet(self.family,label,N,ref.P0,f.logs[0],sc,saved,self.source_family,self.parameters.eta_log,self.parameters.dstar_log)
        self.inlets[key]=packet;return packet

    def modification_support(self,chart,left,right):
        """Executable unchanged-left / original modified-right source selector."""
        c=self.c;sc=self.macro.sc;l,r=c.mpf(left),c.mpf(right)
        if chart not in micro.CHARTS:raise ValueError('Original microscopic source chart required')
        lo,hi=micro.ep(l)[0],micro.ep(r)[1]
        if chart=='first_micro' and hi<=micro.ep(sc/2)[0]:region='unchanged_left'
        elif chart=='first_micro' and lo>=micro.ep(sc/4)[0] and hi<=micro.ep(sc*3/4)[0]:region='original_flat_left_collar'
        elif chart=='first_micro' and lo<micro.ep(sc/2)[1]<hi:region='unchanged_modified_union'
        else:region='original_modified_right'
        return dict(region=region,increments_exact_zero=region in ('unchanged_left','original_flat_left_collar'),
            cut_source='R=Ra*exp(hb*s_c/2)',source_provider_sha256=sha(Path(micro.__file__).name),
            selector_sha256=sha(Path(__file__).name),source_velocity_left='unchanged actual background',
            source_velocity_right='actual E*exp(A/N), actual V+B/N with original Section11 inverse primitives',
            exact_collar_join_by_same_selected_eta_and_strict_kappa_excess=True,
            original_P0_and_background_histories_retained=True)

    def query(self,label,chart,left,right,N=257):
        N=phase.candidate_N(N);op=self.micro.owner(label);f,c=op.flow,op.c
        ref=self.tail.reference.owner(label);_,_,axis=self.switch.r100.owner(label)
        l,r=c.mpf(left),c.mpf(right);key=(label,chart,l._mpi_,r._mpi_,N)
        if key in self.cache:return self.cache[key]
        with mp.workdps(c.dps+40):
            source=op.evaluate(chart,l,r);proxy,raw,recovered,a=switch.recover_source(ref,axis,source)
            recovered.pop('source_frame_conditional_on_same_current_R100_background')
            recovered['source_frame_from_same_actual_Ra_micro_background']=True
            roots,q,proof=switch.general_quotients(proxy,recovered,a,self.parameters.eta_log)
            support=self.modification_support(chart,l,r)
            if support['increments_exact_zero']:
                q={C0:f.scalar(0),Z:f.scalar(0)}
                roots['q']=q[C0]
                got=dict(values={key:f.scalar(0) for key in phase.OUTPUTS},
                    record=dict(original_exact_flat_inverse_A_B_and_Z_zero=True,
                        original_active_inverse_not_evaluated_on_unchanged_or_flat_collar=True))
                proof['original_source_support_q_C0_Z_exact_zero']=True
            else:
                got=primitives.all_u_primitive_bounds(f,roots,q,self.parameters.dstar_log,c.mpf([0,1]))
                got=switch.bounded_exponent_cover(f,got,N)
                if q[C0].zero and q[Z].zero:
                    got['values']={key:f.scalar(0) for key in got['values']}
                    got['record']['original_exact_flat_inverse_A_B_and_Z_zero']=True
                if support['region']=='unchanged_modified_union':
                    for key,row in got['values'].items():
                        lo,hi=micro.ep(row.coefficient)
                        got['values'][key]=type(row)(row.scale,c.mpf([min(mp.mpf(0),lo),max(mp.mpf(0),hi)]),row.ledger)
                    got['record']['source_owned_unchanged_modified_union_includes_exact_zero']=True
            E,V=recovered['common_velocity_E_axial5'],recovered['common_velocity_V_axial5']
            density=phase.densities.density_Z_kernels(E[0],E[1],V[0],V[1],got['values'],N)
            result=dict(source_family=self.family,source_frame=label,candidate_N=N,exact_common_P0_axial5=ref.P0,
                actual_micro_background_source=source,original_raw_source=raw,original_generic_source=recovered,
                original_full_source_quotients=proof,original_roots=roots['roots'],original_q_C0_Z=q,
                original_primitive_values=got['values'],original_primitive_proof=got['record'],
                original_signed_five_density_C0_Z=density,
                actual_source_owned_modified_velocity_C0_Z=density['velocities'],source_owned_modification_support=support,
                actual_source_owned_piecewise_modification_installed=True,
                actual_phase_definition='frac(N*hb*(s-s_c/2))',actual_phase_Z_exact_zero=True,
                actual_selected_positive_s_c=self.macro.sc,actual_positive_left_offset=f.h*(self.macro.sc/2),
                actual_unwrapped_log_radius_phase=source['exact_log_radius_over_Ra']-f.h*(self.macro.sc/2),
                actual_phase_full_period_cover=c.mpf([0,1]),phase_cancellation_not_claimed=True,
                source_intervals_are_enclosures_not_selected_values=True,**dict.fromkeys(fields.previous.OPEN,False))
        self.cache[key]=result;return result

    def seam(self,label):
        op=self.micro.owner(label);f=op.flow;macro_op,_=self.macro.moments.owner(label)
        back=op.evaluate('second_micro',2,2);next_=macro.macro_cell(macro_op,op.series,(0,1),(0,1))
        rows=0
        for name in ('phi','V'):
            for a,b in zip(back['fields'][name],next_['fields'][name]):
                delta=a-b;lo,hi=micro.ep(delta.coefficient)
                if not delta.zero and not lo<=0<=hi:raise ValueError('Same defining micro/macro field source covers must meet')
                rows+=1
        for name in back['histories']:
            for a,b in zip(back['histories'][name],next_['histories'][name]):
                delta=a-b;lo,hi=micro.ep(delta.coefficient)
                if not delta.zero and not lo<=0<=hi:raise ValueError('Same defining micro/macro history source covers must meet')
                rows+=1
        return dict(same_original_source_defining_functions=True,
            original_prescribed_shear_phi_V_integral_definition_shared=True,
            original_six_histories_equal_by_same_initial_atoms_rhs_and_linear_ODE_uniqueness=True,
            independent_source_enclosure_consistency_rows=rows,
            interval_overlap_is_consistency_only_not_the_function_identity_proof=True,
            exact_radius_identity='Ra*exp(hb*2)=Ra*exp(2hb)',
            exact_phase_identity='hb*(2-s_c/2)=(Y-2hb)*0+hb*(2-s_c/2)',
            same_actual_micro_and_macro_source_flow=f is macro_op.flow,
            same_core_atoms_pressure_amplitude_width_and_phase_offset=True)

    def contribution(self,label,N=257):
        N=phase.candidate_N(N);inlet=self.inlet(label,N)
        op=self.micro.owner(label);f,c=op.flow,op.c;ref=self.tail.reference.owner(label)
        local={name:list(rows) for name,rows in inlet['actual_finite_N_inlet_defect_C0_Z'].items()};windows={}
        with mp.workdps(c.dps+40):
            start=self.macro.sc/2
            for chart,partition,end in (('first_micro',(start,c.mpf('.25'),c.mpf('.5'),c.mpf('.75'),c.mpf(1)),c.mpf(1)),
                ('second_micro',(c.mpf(1),c.mpf('1.25'),c.mpf('1.5'),c.mpf('1.75'),c.mpf(2)),c.mpf(2))):
                total={name:[f.scalar(0),f.scalar(0)] for name in RATES};cells=[]
                for left,right in zip(partition,partition[1:]):
                    source=self.query(label,chart,left,right,N);rows={};weight_rows={}
                    for name,rate in RATES.items():
                        mass,decay,tail=weights(op.series,left,right,end,rate)
                        parameters.positive_source(f,mass,'actual_micro_full_cell_Duhamel_mass')
                        density=source['original_signed_five_density_C0_Z']
                        pair=[bounds.symmetric(f,bounds.magnitude(f,density[part][name])*mass*tail) for part in ('kernels','Z_derivatives')]
                        for n,row in enumerate(pair):total[name][n]+=row
                        rows[name]=pair;weight_rows[name]=dict(actual_full_mass=mass,actual_cell_memory=decay,
                            actual_suffix_memory=tail,own_rate=str(rate),true_physical_measure_once=True,
                            actual_compact_left_endpoint=left,actual_compact_right_endpoint=right)
                    cells.append(dict(source=source,actual_signed_cell_driver_C0_Z=rows,actual_own_rate_weights=weight_rows))
                    print('Current original actual micro finite-N cell',label,chart,flush=True)
                memory={name:weights(op.series,partition[0],partition[-1],end,rate)[1] for name,rate in RATES.items()}
                for name in RATES:local[name]=[local[name][n]*memory[name]+total[name][n] for n in range(2)]
                windows[chart]=dict(actual_micro_source_cells=cells,actual_local_signed_driver_C0_Z=total,
                    actual_incoming_memory=memory,actual_propagated_exit_correction_C0_Z={name:list(rows) for name,rows in local.items()})
            joined=self.seam(label);accepted=self.saved_macro['frames'][label];child=self.saved_switch['frames'][label]
            for packet in (accepted,child):switch.guard_saved_source(self.family,label,N,ref.P0,packet)
            restore=lambda row:switch.first.endpoint.restore_row(f,row)
            rows=lambda key,packet:{name:[restore(row) for row in values] for name,values in packet[key].items()}
            memory=lambda key,packet:{name:restore(row) for name,row in packet[key].items()}
            def apply(input_,drive_,memory_):return {name:[input_[name][n]*memory_[name]+drive_[name][n] for n in range(2)] for name in RATES}
            R100=apply(local,rows('actual_micro_exit_R100_local_driver_C0_Z',accepted),memory('retained_micro_exit_R100_incoming_memory',accepted))
            r110_memory={name:f.scalar(c.exp(-c.mpf(rate.numerator)/rate.denominator*c.ln(c.mpf(11)/10))) for name,rate in RATES.items()}
            R110=apply(R100,rows('actual_R100_R110_local_driver_C0_Z',child),r110_memory)
            Rm=apply(local,rows('actual_micro_exit_Rm_local_driver_C0_Z',accepted),memory('retained_micro_exit_Rm_incoming_memory',accepted))
        return dict(source_family=self.family,source_frame=label,candidate_N=N,exact_common_P0_axial5=ref.P0,
            actual_source_bound_Section11_inlet=inlet,actual_micro_source_windows=windows,
            typed_actual_micro_macro_source_join=joined,actual_current_R0_correction_C0_Z=local,
            actual_current_R100_correction_C0_Z=R100,actual_current_R110_correction_C0_Z=R110,
            actual_current_Rm_incoming_correction_C0_Z=Rm,
            actual_current_rminus_R0_R100_R110_Rm_prefix_integrated=True,
            genuine_current_finite_N_R0_R100_R110_Rm_boundary_enclosures_supplied=True,
            zero_before_rminus_is_explicit_source_owned_modification_definition=True,
            local_drivers_composed_with_real_current_prefix_not_arbitrary_zero=True,
            correction_rows_are_signed_function_enclosures_not_selected_point_values=True,
            no_old_N1024_or_ancestor_owner_boundary_transplanted=True,
            conditional_source_frames=('0','.5'),global_N_or_whole_Z_admission=False,
            **dict.fromkeys(fields.previous.OPEN,False))


def run():
    began=time.monotonic();owner=OriginalMicroFiniteN(require_checked=False)
    frames={label:owner.contribution(label) for label in ('0','.5')}
    report=dict(mode=owner.mode,source_family=owner.family,**{GATE:True},candidate_N=257,
        original_source_and_initial_support_bindings=owner.bindings,frames=serialized(frames),
        input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,**dict.fromkeys(fields.previous.OPEN,False))
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(macro.base.encoded(report),separators=(',',':'))+'\n').encode(),compresslevel=6,mtime=0))
    return report


if __name__=='__main__':run()
