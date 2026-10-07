"""Attach the analytic Ra<R<=110 relaxed input to the actual current graph.

Transfer is through current core/bridge/switch callable, integral, six-
history, pressure and exact-width identities, not a family ID or overlap.
The untouched Ra endpoint is the separate stress-free core case. No source
constructor, numerical microscopic width, repair or global N is used.
"""
import ast
import json
import math
from pathlib import Path
from types import SimpleNamespace
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_generic_shear_inputs as inputs

HERE,PREFIX,sha=inputs.HERE,inputs.PREFIX,inputs.sha
NAME=PREFIX+'current_inner_relaxed_inputs.json'
RECEIPT=PREFIX+'current_inner_relaxed_inputs_check.json'
GATE='current_original_inner_Ra_R110_relaxed_generic_input_certified'
OPEN=inputs.OPEN
CHARTS=('bridge_first','bridge_second','bridge_macro','switch_first','switch_second','switch_power')


def exact_switch_source_theorem():
    asts=inputs.packets.recovery.numeric.transport.SourceAST();bindings={}
    def assignment(stem,method,target,wanted):
        fn=asts.method(stem,method);expected=ast.dump(ast.parse(wanted,mode='eval').body)
        if not any(isinstance(n,ast.Assign) and any(ast.unparse(t)==target for t in n.targets)
                   and ast.dump(n.value)==expected for n in ast.walk(fn)):
            raise ValueError('Admitted original inner shear source changed: '+stem+'.'+target)
        bindings[stem+'.'+method+'/'+target+'#'+str(len(bindings))]=wanted
    specs=(
      ('microswitch_mixed_C4','switch_controls','tinyD','[scale_width(row,k+1) for k,row in enumerate(Dbar_y)]'),
      ('microswitch_mixed_C4','switch_controls','a','tinyD'),
      ('microswitch_mixed_C4','switch_controls','a',"[sum((tinyD[j]*(math.comb(k,j)*complement[k-j]) for j in range(k+1)),tinyD[0]*0)+sig[k]*c.mpf('.8') for k in range(4)]"),
      ('microswitch_mixed_C4','switch_controls','Vphase','[quotient*0]*4'),
      ('microswitch_mixed_C4','switch_controls','driver','[drive_y(j,j+2) for j in range(4)]'),
      ('microswitch_mixed_C4','switch_controls','logU','[scale_width(one-a[0],1)/2]+[scale_width(row,1)*(-c.mpf(\'.5\')) for row in a[1:]]'),
      ('microswitch_mixed_C4','postpower','V',"jet(inlet['Uz_actual_axial5_coefficients'])"),
      ('microswitch_mixed_C4','postpower','phi',"phi2*theta**c.mpf('.4')"),
      ('microswitch_mixed_C4','postpower','moments','power_transport(c,theta,as_initial(moments2),phi2,V)'),
      ('microswitch_mixed_C4','postpower','packet',"reshape_mixed(c,Z,self.core.delta,logu,[constant('.1')]+[constant(0)]*3,V,shapes,jet(inlet['pressure_axis_axial5_coefficients']),self.invP2,self.proofs)"))
    for spec in specs:assignment(*spec)
    c=SimpleNamespace(mpf=lambda v:s.Rational(str(v)));env=dict(math=math)
    asts.replay('long_reshape_mixed_C4','exponential_derivatives',env)
    fn=asts.replay('microswitch_mixed_C4','switch_controls',env)
    hb,D,E,R,Fa,Fbar=s.symbols('hb Dbar Ebar R Factual Fbar',positive=True)
    E=s.Symbol('signed_Ebar',real=True);sig=s.Symbol('same_sigma',real=True)
    ds=s.symbols('D0:4');gs=s.symbols('G0:4');cutoff=[sig]+[s.Integer(0)]*4
    drive=lambda j,power:(s.sqrt(R/2)*Fbar*E if j==0 else gs[j])*hb**power
    width=lambda value,power=1:value*hb**power
    controls={branch:fn(c,branch,cutoff,[D,*ds[1:]],drive,Fa/Fbar,width) for branch in ('first','second')}
    checks={}
    def zero(name,left,right):
        if s.cancel(left-right)!=0:raise ArithmeticError('Current inner source identity: '+name)
        checks[name]=True
    for branch,control in controls.items():
        a=1-2*control['logUtheta_phase_derivatives'][0]/hb
        b=2*(control['Uz_positive_phase_derivatives'][0]/hb)/(s.sqrt(2*R)*Fa)
        zero(branch+'_same_angular_a',a,hb*D if branch=='first' else hb*D*(1-sig)+s.Rational(4,5)*sig)
        zero(branch+'_same_signed_axial_b',b,-hb*E*(1-sig) if branch=='first' else 0)
    y=s.Symbol('original_logR',real=True);f2=s.Symbol('same_R2_F',positive=True)
    F=f2*s.exp(-s.Rational(2,5)*y);U=s.sqrt(2*R)*s.exp(y/2)*F
    zero('same_postpower_a_four_fifths',1-2*s.diff(U,y)/U,s.Rational(4,5))
    zero('same_postpower_radial_constant_V',s.diff(s.Symbol('same_R2_V'),y),0)
    Db,Eb,chi=s.symbols('Dbar Ebar chi',positive=True);Eb=s.Symbol('signed_Ebar',real=True)
    et,ez=s.symbols('actual_signed_inertial_error_theta actual_signed_inertial_error_axial',real=True)
    a=chi*Db;b=-chi*Eb;t0=-b/a;p1=Db+et;p2=Eb+ez
    Hbar=(Db*Db+Eb*Eb)/Db
    zero('whole_bridge_H0_independent_of_chi',p1+t0*p2,Hbar+et+(Eb/Db)*ez)
    zero('same_first_switch_kappa',hb*Db+(-hb*Eb*(1-sig))**2/(hb*Db),
         hb*(Db*Db+Eb*Eb*(1-sig)**2)/Db)
    return dict(passed=True,original_switch_and_postpower_AST_bindings=bindings,
        exact_shear_source_identities=checks,
        same_original_micro_controls_and_complete_power_source=True,
        signed_axial_three_scale_drive_and_full_histories_retained=True,input_hashes=asts.hashes)


class CurrentInnerRelaxedInputs:
    def __init__(self,service=None):
        self.service=inputs.packets.CurrentSourcePackets() if service is None else service
        self.ctx=self.service.ctx;self.family=self.service.family;self.conditions=[];self.rows={}
        input_receipt=json.loads((HERE/inputs.RECEIPT).read_bytes())
        if not input_receipt.get('all_passed') or not input_receipt.get(inputs.GATE) or input_receipt.get('source_family')!=self.family:
            raise ValueError('Same checked complete inertial/shear input expressions required')
        self.service.bind_hashes(input_receipt['input_hashes']);self.service.bind_hashes({inputs.RECEIPT:sha(inputs.RECEIPT)})
        def load(stem,checked=True,family=True):
            name=PREFIX+stem+'.json';row=json.loads((HERE/name).read_bytes())
            if checked and not row.get('all_passed'):raise ValueError('Checked current analytic source receipt required: '+stem)
            if family and {k:row[k] for k in inputs.packets.FAMILY_KEYS}!=self.family:
                raise ValueError('Inner source/axis datum differs: '+stem)
            self.service.bind_hashes(row['input_hashes']);self.service.bind_hashes({name:sha(name)})
            self.rows[stem]=row;return row
        def true(row,path,label):
            value=row
            for key in path:value=value[key]
            if value is not True:raise ValueError('Current inner callable/source identity missing: '+label)
            self.conditions.append(label)
        def all_true(row,key,label):
            if not row[key] or any(v is not True for v in row[key].values()):raise ValueError('Exact current source bindings required: '+label)
            self.conditions.append(label)
        collar=load('current_inner_exit_strict_collar_check')
        true(collar,('current_inner_exit_strict_collar_attached_to_current_source_graph_certified',),'current actual core/first/width attached')
        attached=load('current_inner_exit_strict_collar',checked=False)
        source=attached['exact_current_exit_source_attachment']
        true(source,('passed',),'bridge prescribed shear and six integral function identities')
        if len(source['same_function_proof'])<7:raise ValueError('Whole prescribed bridge source-function attachment required')
        for key in ('actual_amplitude_packet_graph','exact_core_atom_and_pressure_functions','exact_scaled_atoms','exact_pressure_primitive','original_direction_equality'):
            value=attached['cross_graph_amplitude_packet_pressure_width_attachment'][key]
            if not value or any(v is False for v in value.values()):raise ValueError('Current bridge common source attachment missing: '+key)
        first=load('current_core_first_interface_check')
        for key in ('current_core_bridge_functional_mixed4_join_certified','all_current_bridge_functional_interfaces_certified'):
            true(first,(key,),'new separate core-first receipt/'+key)
        micro=load('current_microswitch_background_tensor_check')
        mt=micro['current_actual_microswitch_source_and_tensor_theorem']['current_actual_microswitch_source_pressure_units_and_two_endpoint_theorem']
        all_true(mt,'live_original_callable_bindings','micro same original source callables')
        original=mt['current_original_actual_input_bindings']
        for key in ('original_mixed4_callable_identities','original_mixed4_actual_input_AST_bindings'):
            all_true(original,key,'micro actual input/'+key)
        for key in ('actual_switch_object_is_current_history','canonical_pressure_and_ratios_from_current_actual_bridge','source_width_not_a_parameter_value'):
            true(original,(key,),'micro exact inherited source/'+key)
        norm=mt['original_phase_to_ordinary_y_and_full_raw_history_normalization']
        for key in ('passed','original_phase_and_owned_raw_adapter_AST_replayed','five_full_history_ODEs_and_absolute_pressure_datum_retained','inverse_width_never_materialized_or_selected_from_cap'):
            true(norm,(key,),'micro ordinary units/'+key)
        true(norm,('identities','absolute_pressure_y_ODE'),'same full pressure ODE')
        endpoint=mt['original_phase1_and_R2_complete_source_function_theorem']
        for key in ('passed','exact_same_current_R2_phi_V_six_moments_and_absolute_P0','source_functions_compared_before_width_or_amplitude_enclosures'):
            true(endpoint,(key,),'micro same function endpoints/'+key)
        if not endpoint['original_postpower_source_bindings']:raise ValueError('Original postpower function source required')
        switch=load('actual_switch_mixed_C4_check')
        for key in ('actual_R100_functional_mixed4_join_certified','actual_finite_width_bridge_histories_installed_in_switch_mixed4','actual_R100_R110_feedback_mixed4_available'):
            true(switch,(key,),'same current bridge/switch histories/'+key)
        power=load('current_switch_power_background_tensor_check')
        pt=power['current_actual_switch_power_source_and_tensor_theorem']['current_actual_switch_power_source_pressure_and_R110_theorem']
        all_true(pt,'live_original_callable_bindings','power same original source callables')
        adapter=pt['original_output_only_source_locals_adapter']
        for key in ('original_postpower_AST_unchanged_except_one_output_keyword','every_original_nonreturn_statement_retained_exactly','existing_locals_logu_shapes_V_P0_exposed_without_feedback'):
            true(adapter,(key,),'power output-only source exposure/'+key)
        endpoint=pt['actual_original_R110_power_post_and_reshape_boundary_function_theorem']
        for key in ('passed','same_original_R2_moments_V110_F0_and_P0_not_enclosure_overlap','same_reshape_mixed_generator_and_flat_log_jets_imply_R110_full_boundary_mixed4','positive_width_and_G_caps_not_selected_as_functions'):
            true(endpoint,(key,),'power same function/'+key)
        all_true(endpoint,'exact_six_current_power_transport_function_identities','same six current power transport integrals')
        if endpoint['physical_source_function_identities']!=['Utheta','Uz','Ur','P_over_Pstar2','Mtheta','Mtheta_z','Mz','Mztheta','Mp'] or endpoint['total_physical_mixed4_rows_implied']!=135:
            raise ValueError('All original power physical function identities required')
        self.conditions.append('nine current power physical functions and all135 mixed4 source rows')
        pressure=pt['checked_actual_analytic_pressure_function']
        for key in ('same_fourteen_implicit_stage_sources','callable_constructor_and_equations_not_hash_only_proof'):
            true(pressure,(key,),'same analytic pressure/'+key)
        if pressure['exact_analytic_defining_function']!='P0/Pstar^2=-(M2*(1+Z^2)^-2+M0+F_flat(Z))':
            raise ValueError('Original analytic pressure defining function differs')
        global_exit=load('global_exit_certificate',checked=False,family=False)
        if source['admitted_global_exit_prescription']!=global_exit['field_prescription']:
            raise ValueError('Actual prescribed field differs from admitted global exit functions')
        for key in ('actual_whole_axis_Ra_R110_relaxed_cone_analytically_certified','exact_implicit_exit_field_specified'):
            true(global_exit,(key,),'original analytic whole exit/'+key)
        ledger=load('K1_ledger',checked=False,family=False);norms=load('physical_norm_family',checked=False,family=False)
        attachment=attached['cross_graph_amplitude_packet_pressure_width_attachment']
        for key in ('base_analytic_core_family_sha256','uniform_Cstar_family_sha256','admitted_inner_parameter_family_sha256'):
            if attachment[key]!=global_exit[key] or global_exit[key]!=ledger[key]:raise ValueError('Original exit parameter family differs: '+key)
        for key in ('uniform_Cstar_family_sha256','base_analytic_core_family_sha256'):
            if norms[key]!=global_exit[key]:raise ValueError('Current norm/core family differs')
        if not all(ledger['smallness_checks'].values()) or not ledger['h_b_equals_epsilon_b_by_definition']:
            raise ValueError('Exact positive width/norm smallness ledger required')
        self.theorem=exact_switch_source_theorem();self.service.bind_hashes(self.theorem['input_hashes'])
        self.service.bind_hashes({Path(__file__).name:sha(Path(__file__).name)})
        self.proof=self.bounds()

    def bounds(self):
        c=self.ctx;read=lambda v:inputs.packets.interval(c,v)
        ledger=self.rows['K1_ledger'];exit=self.rows['global_exit_certificate'];norms=self.rows['physical_norm_family']
        gamma=read(ledger['gamma']);weak_error=read(ledger['decreasing_K_smallness_bounds']['weak_cone_error'])
        bridge_margin=2*gamma-weak_error
        first_margin=read(exit['first_switch_relaxed_lower_margin']);last_margin=c.mpf(1)
        margin=c.mpf(min(inputs.packets.recovery.endpoints(v)[0] for v in (bridge_margin,first_margin,last_margin)))
        L=read(norms['selected_log_K_upper']);loghb=read(ledger['shared_positive_width_log_enclosure'])
        values=dict(bridge_H0_minus2_lower=bridge_margin,first_switch_H0_minus2_lower=first_margin,
            second_and_power_H0_minus2_lower=last_margin,whole_inner_H0_minus2_lower=margin,
            source_hb_below_half_log_gap=c.ln(c.mpf('.5'))-loghb,
            original_K_exceeds_one_million_log_gap=L-c.ln(1000000),
            original_strong_signed_Q_error_ratio_margin=read(exit['strong_cone_sufficient_ratio_margin']),
            original_first_e2_squared_Dbar_margin=read(exit['first_switch_e2_squared_comparison_D_lower_margin']))
        if any(inputs.packets.recovery.endpoints(v)[0]<=0 for v in values.values()):
            raise ArithmeticError('Current analytic inner input bound unresolved')
        return dict(physical_domain='Ra<R<=110, all Z[-1,1]; original untouched Ra endpoint separate',
            original_source_charts=list(CHARTS),positive_directed_margins=values,
            positive_profile_derivation='positive original f and exact angular exponential/power; no cap representative',
            positive_a_derivation='bridge chi>=hb,Dbar>=1/(2K); first hb*Dbar; second convex interpolation to4/5; power4/5',
            source_log_a_lower_bound=loghb-c.ln(2)-L,
            a_lower_bound_formula='a>=hb/(2Kmax)>0; bound never materializes exact source hb',
            source_log_kappa_upper_bound=10*L,
            source_log_abs_t0_upper_bound=c.ln(4)+2*L,
            t0_bound_formula='|t0|<=|Ebar|/Dbar<=4Kmax^2; first switch adds (1-sigma), second/power t0=0',
            bridge_H0_derivation='H0=Hbar+e_theta+(Ebar/Dbar)*e_z; Hbar>=2+2gamma, |projection error|<=weak_error',
            first_H0_derivation='D_actual>=3.25 and Ebar*E_actual>=-e2^2/4, e2^2<=Dbar; H0-2>=3.25-.25-2',
            second_power_H0_derivation='b=0 and inertial p1=D_actual>3, hence H0-2>1',
            original_relaxed_branches='kappa<=2: H0>2; kappa>2: signed D>0 and signed Q>0 from exact global exit theorem',
            original_core_zero_inlet_case_separate=True,
            full_current_I_and_S_same_units_and_source=True,
            absolute_pressure_and_cumulative_moments_not_reset=True,
            width_or_scalar_caps_not_used_as_source_functions=True,
            p1_p2_whole_inner_norm_bounds_certified=False,
            original_whole_inner_relaxed_input_certified=True,**dict.fromkeys(OPEN,False))

    def query(self,chart,coordinate,Z=(-1,1)):
        c=self.ctx
        if chart not in CHARTS:raise ValueError('Original analytic inner source chart required')
        v=c.mpf(coordinate);z=c.mpf(Z);vl,vh=inputs.packets.recovery.endpoints(v);zl,zh=inputs.packets.recovery.endpoints(z)
        lower,upper=inputs.packets.DOMAINS[chart]
        if not all(mp.isfinite(x) for x in (vl,vh,zl,zh)) or vl<lower or vh>upper or zl<-1 or zh>1:
            raise ValueError('Original finite inner source subbox required')
        includes_Ra=chart=='bridge_first' and vl==0
        return dict(chart=chart,coordinate=v,Z=z,source_family=self.family,
            relaxed_input_certified_for_every_source_point_in_query=not includes_Ra,
            excluded_Ra_stress_free_endpoint_may_be_in_query=includes_Ra,
            exact_zero_Ra_endpoint_only=includes_Ra and vh==0,
            relaxed_or_separate_zero_core_case_certified=True,
            analytic_source_certificate=self.proof,no_source_point_values_evaluated=True,
            **{GATE:not includes_Ra},**dict.fromkeys(OPEN,False))


def run():
    owner=CurrentInnerRelaxedInputs()
    examples={chart:owner.query(chart,('.01','1') if chart=='bridge_first' else inputs.packets.DOMAINS[chart]) for chart in CHARTS}
    examples['closed_first_with_separate_zero']=owner.query('bridge_first',(0,1))
    examples['exact_zero_Ra']=owner.query('bridge_first',0)
    result=dict(source_family=owner.family,current_actual_inner_source_function_attachment_conditions=owner.conditions,
        exact_original_micro_power_and_H0_source_theorem=owner.theorem,
        current_whole_inner_relaxed_input_and_bounds=owner.proof,analytic_source_subbox_examples=examples,
        bound_source_receipts={name:PREFIX+name+'.json' for name in owner.rows},
        **{GATE:True},**dict.fromkeys(OPEN,False),
        strict_completed_tensor_cone_new_regions_admitted=0,
        source_graph_ancestor_constructors_called=False,
        scope='Original current analytic Ra<R<=110 relaxed input and formal/log a,t0 bounds, attached through exact current source identities. Not strict whole cone, p1/p2 norm bounds, full upstream admission, loop, repair or new N.',
        input_hashes=owner.service.hashes)
    (HERE/NAME).write_text(json.dumps(inputs.packets.encode(result),indent=2)+'\n',encoding='utf8')
    print('Actual current analytic inner relaxed input attached: Ra<R<=110',flush=True)
    return result


if __name__=='__main__':run()
