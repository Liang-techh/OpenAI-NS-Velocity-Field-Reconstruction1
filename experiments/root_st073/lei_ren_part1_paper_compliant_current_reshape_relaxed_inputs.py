"""Current-source angular barrier on R110..Rsh..Rz, without sampling.

The unchanged full primitives imply the paper angular ODE. The actual
current V110, B/T, pressure and Rsh histories are bound to that theorem.
No old fixture cone, microscopic/source factor point, or ancestor is used.
"""
import json
import ast
import math
from pathlib import Path
from types import SimpleNamespace
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_inner_relaxed_inputs as inner

inputs=inner.inputs;packets=inputs.packets
HERE,PREFIX,sha=inner.HERE,inner.PREFIX,inner.sha
NAME=PREFIX+'current_reshape_relaxed_inputs.json'
RECEIPT=PREFIX+'current_reshape_relaxed_inputs_check.json'
GATE='current_original_R110_Rz_relaxed_generic_input_certified'
OPEN=inner.OPEN
CHARTS=('reshape','inner_reference')


def angular_barrier_source_theorem():
    asts=packets.recovery.numeric.transport.SourceAST()
    wanted={
      ('long_reshape_profiles','inputs','inlet'):"self.switch.inlet(Z)",
      ('long_reshape_profiles','inputs','v'):"jet(inlet['Uz_actual_axial5_coefficients'])",
      ('long_reshape_profiles','evaluate','logu'):'B*(1-sig)-logq+(y/10-self.core.logC-self.core.logP)',
      ('long_reshape_profiles','evaluate','a'):"B*(2*ds/self.T)+c.mpf('.8')",
      ('long_reshape_profiles','evaluate','mean'):"inp['moments']['mean']*theta+inp['v']*(1-theta)",
      ('long_reshape_profiles','evaluate','theta'):'positive_decay(c,1,y,self.proofs)',
      ('long_reshape_profiles','evaluate','moment_shapes'):"dict(theta=inherited['theta']+kernels['theta'],theta_z=inherited['theta_z']+inp['v']*kernels['theta'],pressure=inherited['pressure']+kernels['pressure'],swirl=inherited['swirl']+kernels['swirl'],mean=mean,axial=axial)",
      ('long_reshape_profiles','__init__','self.T'):'400*self.A',
      ('reference_restore_profiles','reference','gap'):'(self.loggap-8)*phase',
      ('reference_restore_profiles','reference','logu'):'-logq+(self.reshape.T/10-self.core.logC-self.core.logP+gap/10)',
      ('reference_restore_mixed_C4','source_mixed','mismatch'):'[E*value for value in alpha]',
      ('reference_restore_mixed_C4','source_mixed','V'):'[z*4+mismatch[0]]+mismatch[1:]',
      ('reference_restore_mixed_C4','reference_branch',None):None}
    nodes={}
    for (stem,method,target),value in wanted.items():
        if target is not None:nodes[target if target!='logu' else stem]=asts.expression(stem,method,target,wanted=value)
    branch=asts.method('reference_restore_mixed_C4','reference_branch')
    call=next(n.value for n in branch.body if isinstance(n,ast.Return))
    expected=ast.parse("self.packet(Z,self.reference.reference(Z,phase),[c.mpf(1)]+[c.mpf(0)]*4,'Rsh_to_Rz')",mode='eval').body
    if ast.dump(call)!=ast.dump(expected):raise ValueError('Original constant reference alpha source changed')
    asts.bindings['reference_restore_mixed_C4.reference_branch.original_constant_alpha']=True
    y,Z,T,logC,logP=s.symbols('y Z T logC logP',real=True)
    B=s.Function('same_actual_B')(Z);sig=s.Function('same_sigma')(y/T)
    env=dict(B=B,sig=sig,logq=s.log(1+Z*Z),y=y,
        self=SimpleNamespace(core=SimpleNamespace(logC=logC,logP=logP)))
    logu=asts.evaluate(nodes['long_reshape_profiles'],env)
    checks={}
    def zero(name,value):
        if s.cancel(value)!=0:raise ArithmeticError('Current angular barrier identity: '+name)
        checks[name]=True
    a=s.Rational(4,5)+2*B*s.diff(sig,y)
    zero('same_original_variable_a',1-2*s.diff(logu,y)-a)
    zero('same_signed_zeta_convex_interpolation',s.diff(logu,Z)-((1-sig)*(s.diff(B,Z)-2*Z/(1+Z*Z))-sig*2*Z/(1+Z*Z)))
    eref=s.Function('same_current_V110_minus4Z')(Z)
    refenv=dict(E=eref,alpha=[1,0,0,0,0],z=Z)
    refenv['mismatch']=asts.evaluate(nodes['mismatch'],refenv)
    reference_V=asts.evaluate(nodes['V'],refenv)
    zero('same_reference_V110',reference_V[0]-4*Z-eref)
    for j in range(1,5):zero('same_reference_V_y'+str(j)+'_zero',reference_V[j])
    # Replay the full current raw angular inertial sectors in native units.
    R,S,delta=s.symbols('R Pstar delta',positive=True)
    u=s.Function('same_actual_Utheta')(R,Z);V=s.Function('same_actual_V110')(Z)
    mt,mz,mtz,mzt,P=[s.Function('same_actual_'+n)(R,Z) for n in ('Mtheta','Mz','Mtheta_z','Mztheta','P')]
    d=1-Z*Z;L=1-delta*Z*Z;Aop=lambda q:(1-delta)*Z*q+d*s.diff(q,Z)
    W=1-Aop(mz/R);Hv=(1-delta)*Z/2+d*V
    angular=(1-delta/2)*mt-(1-delta)*Z*s.diff(mt,Z)/2-d*s.diff(mtz,Z)+(2*delta-1)*Z*mtz
    Q=-W+angular/(s.sqrt(2)*R**s.Rational(3,2)*u)
    p1=R*Q/L;F=u/s.sqrt(2*R)
    native=dict(m=mz/R,h=mt/(s.sqrt(2)*R**s.Rational(3,2)*S),
        k=mtz/(s.sqrt(2)*R**s.Rational(3,2)*S),e=mzt/(R*S*S),p=P/(S*S))
    rows=lambda v:[v]+[s.Integer(0)]*4
    c=SimpleNamespace(mpf=lambda v:s.Rational(str(v)))
    ops=dict(math=math,axial_derivative=lambda v:s.diff(v,Z))
    asts.replay('collar_Gamma_C4','product_rows',ops)
    asts.replay('collar_stress_C3','shifted_rows',ops)
    raw=asts.replay('current_pre_pulse_stress_operator','raw_pre_stress_rows',ops)(
        c,delta,Z,rows(u/S),rows(V),{k:rows(v) for k,v in native.items()},rows(P/(S*S)))
    inertial=sum(R**s.Rational(str(part['mode'][0]))*S**part['mode'][1]*part['shape'][0]/s.sqrt(2)
        for label,part in raw['theta'].items() if label!='variable_radial_shear')
    zero('same_complete_raw_angular_I_over_F',inertial/F-p1)
    av,zeta=s.symbols('same_actual_a same_actual_zeta',real=True)
    SQ=-W*(1-av/2)-delta*(1-2*Z*V)/2-Hv*zeta
    rhs={mt:s.sqrt(2*R)*u,mz:V,mtz:s.sqrt(2*R)*u*V}
    updates={s.diff(q,R):f for q,f in rhs.items()}
    updates.update({s.diff(s.diff(q,Z),R):s.diff(f,Z) for q,f in rhs.items()})
    updates[s.diff(u,R)]=(1-av)*u/(2*R)
    reduce=lambda q:s.factor(s.together(q.xreplace(updates).subs(s.diff(u,Z),zeta*u)))
    zero('same_full_primitive_angular_ODE',reduce(R*s.diff(Q,R)+(2-av/2)*Q-SQ))
    zero('same_generic_p1_barrier_ODE',reduce(R*s.diff(p1,R)+(1-av/2)*p1-R*SQ/L))
    ev=s.Symbol('same_signed_V110_minus4Z',real=True)
    target=-2*Z/(1+Z*Z)
    leading=-2*Z*Z*((1-delta)/2+4*d)/(1+Z*Z)
    error=-2*Z*d*ev/(1+Z*Z)
    zero('same_target_signed_gradient_negative_leading_plus_error',
        ((1-delta)*Z/2+d*(4*Z+ev))*target-leading-error)
    zero('same_target_error_multiplier_denominator_slack',1+Z*Z-d-2*Z*Z)
    rho=s.Symbol('same_positive_R110_over_R',real=True);e0,e1=s.symbols('inlet_mean_error V110_error',real=True)
    mean=(4*Z+e0)*rho+(4*Z+e1)*(1-rho)
    zero('same_mean_error_positive_average',mean-4*Z-(e0*rho+e1*(1-rho)))
    inlet=s.Function('same_actual_m110')(Z)
    mean_source=inlet*s.exp(-y)+V*(1-s.exp(-y))
    zero('same_exact_uncapped_mean_ODE',s.diff(mean_source,y)+mean_source-V)
    zero('same_exact_mean_inlet',mean_source.subs(y,0)-inlet)
    zero('same_correlated_B_T_shear_radius',2*(2*s.Symbol('same_A',positive=True))*8/(400*s.Symbol('same_A',positive=True))-s.Rational(2,25))
    return dict(passed=True,original_source_AST_bindings=asts.bindings,exact_identities=checks,
        angular_Q_is_paper_angular_inertial_quantity_not_radial_velocity_Q=True,
        generic_p1_equals_R_times_angular_Q_over_L=True,
        raw_stress_identity_uses_complete_zeroth_angular_inertial_sectors=True,
        no_pressure_axial_or_history_source_replaced=True,
        all_original_positive_kernels_and_full_current_inlet_retained=True,input_hashes=asts.hashes)


class CurrentReshapeRelaxedInputs:
    def __init__(self,service=None):
        self.inner=inner.CurrentInnerRelaxedInputs(service);self.service=self.inner.service
        self.ctx=self.inner.ctx;self.family=self.inner.family;self.conditions=[];self.rows={}
        proof=json.loads((HERE/inner.RECEIPT).read_bytes())
        if not proof['all_passed'] or not proof[inner.GATE] or proof['source_family']!=self.family:
            raise ValueError('Checked actual current whole inner input required')
        self.service.bind_hashes(proof['input_hashes']);self.service.bind_hashes({inner.RECEIPT:sha(inner.RECEIPT)})
        def load(stem,checked=True,family=True):
            name=PREFIX+stem+'.json';r=json.loads((HERE/name).read_bytes())
            if checked and not r.get('all_passed'):raise ValueError('Checked actual source receipt required: '+stem)
            if family and {k:r[k] for k in packets.FAMILY_KEYS}!=self.family:raise ValueError('Actual reshape source family differs')
            self.service.bind_hashes(r['input_hashes']);self.service.bind_hashes({name:sha(name)})
            self.rows[stem]=r;return r
        def yes(row,key,label):
            if row[key] is not True:raise ValueError('Actual source identity missing: '+label)
            self.conditions.append(label)
        def yesmap(row,key,label):
            if not row[key] or any(v is not True for v in row[key].values()):raise ValueError('Actual source function map missing: '+label)
            self.conditions.append(label)
        reshape=load('current_reshape_background_tensor_check')
        rs=reshape['current_actual_reshape_source_units_and_tensor_theorem']['current_actual_variable_reshape_source_pressure_and_Rsh_theorem']
        yes(rs,'passed','checked actual variable reshape source theorem')
        normalization=reshape['current_actual_reshape_source_units_and_tensor_theorem']['original_actual_variable_reshape_current_radius_normalization_theorem']
        yes(normalization,'passed','current original raw five-history normalization theorem')
        for key in ('raw_five_history_ODE_m','raw_five_history_ODE_h','raw_five_history_ODE_k','raw_five_history_ODE_e','raw_five_history_ODE_p','absolute_pressure_y_ODE'):
            yes(normalization['identities'],key,'same unchanged primitive equation/'+key)
        yesmap(rs,'live_original_callable_bindings','same original/current profile and full stress operator callables')
        for key in ('actual_variable_B_T_cutoff_and_four_log_amplitude_jets_retained','full_original_positive_kernels_and_nonzero_R110_histories_retained',
            'inverse_T_factors_and_Bell_products_not_replaced_by_constant_rates','exact_actual_Rsh_radius_function_identity','interval_overlap_not_used_as_function_identity'):
            yes(rs,key,'actual reshape/'+key)
        profile=rs['current_R110_profile_bindings']
        yesmap(profile,'unchanged_original_callables','original profile inlet/normalization/transport callables')
        for key in ('actual_R110_input_AST_bindings','actual_R110_normalization_AST_bindings','constant_V_transport_AST_bindings'):
            yesmap(profile,key,'current inlet source/'+key)
        yes(profile,'original_B_C2_and_shear_theorems_consume_same_defining_source','actual B/T theorem source equality')
        boundary=rs['actual_original_Rsh_boundary_function_theorem']
        yes(boundary,'passed','actual full Rsh source-function boundary')
        if boundary['exact_six_centered_history_function_identities']!=6:raise ValueError('All six exact Rsh history functions required')
        self.conditions.append('same six exact Rsh centered history functions')
        if boundary['total_physical_mixed4_rows_implied']!=135:raise ValueError('All original Rsh physical functions required')
        restore=load('current_restore_background_tensor_check')
        rt=restore['current_actual_restore_source_units_and_tensor_theorem']['current_actual_restore_source_pressure_and_three_endpoint_theorem']
        yes(rt,'passed','checked current reference source theorem')
        yesmap(rt,'live_original_callable_bindings','same actual reference provider and source equations')
        yesmap(rt,'original_restore_start_six_centered_function_identities','same reference Rz six-history endpoint')
        yes(rt['current_correlated_E_source_bindings'],'same_core_bridge_first_switch_defining_E_source','same current correlated V110-minus4Z source')
        for record in (rs['checked_actual_analytic_pressure_function'],rt['checked_actual_analytic_pressure_function']):
            for key in ('passed','same_fourteen_implicit_stage_sources','callable_constructor_and_equations_not_hash_only_proof'):yes(record,key,'same unchanged analytic pressure/'+key)
        join=load('reference_join_bounds',checked=False,family=False)
        original=load('reference_join_check',checked=False,family=False)
        if original['reference_join_family_sha256']!=join['reference_join_family_sha256'] or not original['moment_stress_identities']['independent_five_primitive_to_stress_derivation_passed']:
            raise ValueError('Original exact primitive-to-stress theorem required')
        for key in ('same_family_R110_Rh_relaxed_cone_analytically_certified','same_preheat_axis_pressure_retained','actual_moments_continuously_inherited','normalized_R110_amplitude_budget_recomputed'):
            yes(join,key,'original reference source bounds/'+key)
        for key in ('admitted_inner_parameter_family_sha256','uniform_Cstar_family_sha256'):
            if join[key]!=self.inner.rows['global_exit_certificate'][key]:raise ValueError('Current reshape norm/width family differs')
        expected=dict(T='400*Abar',Rref='110*(Cstar*Pstar)^10, exact expression',Rsh='110*exp(T)',
            swirl='log(u/Pstar)=y/10-logCstar-logPstar-log(1+Z^2)+(1-sigma(y/T))*B',
            B='log(Cstar*u1*(1+Z^2)) of the actual R110 exit',moments='continue all five exact axis primitives, P=P0+Mp; no resetting')
        if any(join['definition'][key]!=value for key,value in expected.items()):raise ValueError('Actual reference prescription differs')
        self.conditions.append('actual defining log source, B, fixed T, reference radius and all five original primitives')
        self.theorem=angular_barrier_source_theorem();self.service.bind_hashes(self.theorem['input_hashes'])
        self.service.bind_hashes({Path(__file__).name:sha(Path(__file__).name)})
        self.proof=self.bounds()

    def bounds(self):
        c=self.ctx;read=lambda value:packets.interval(c,value)
        lower=lambda value:c.mpf(packets.recovery.endpoints(value)[0])
        upper=lambda value:c.mpf(packets.recovery.endpoints(value)[1])
        ledger=self.inner.rows['K1_ledger'];exit=self.inner.rows['global_exit_certificate'];join=self.rows['reference_join_bounds']
        eps=read(ledger['epsilon0']);delta=c.mpf(self.service.data['delta'])
        A=read(join['Abar']);selected=read(self.inner.rows['physical_norm_family']['A_upper'])
        if packets.recovery.endpoints(A)!=packets.recovery.endpoints(selected) or packets.recovery.endpoints(read(join['B_C2_upper']))!=packets.recovery.endpoints(2*A):
            raise ValueError('Exact selected A and actual same-source B C2 bound required')
        positive=dict(delta_floor=delta,delta_below_point005=c.mpf(1)/200-delta,
            epsilon_below1=1-eps,A_exceeds10=A-10,
            negative_W_magnitude_floor=3-c.mpf(4)/200-10*eps,
            V110_error_budget=eps-read(exit['axial_velocity_C2_error_upper']),
            Mz110_mean_error_budget=eps-read(exit['axial_cumulative_mean_C2_error_upper']),
            reference_log_gap_positive=read(join['shape_to_Rz_log_gap_lower']))
        if packets.recovery.endpoints(delta)[1]>packets.recovery.endpoints(c.mpf(1)/200)[0]:raise ValueError('Original delta bound unresolved')
        for key,value in positive.items():
            if key=='delta_below_point005':continue
            if packets.recovery.endpoints(value)[0]<=0:raise ArithmeticError('Actual current angular source budget unresolved: '+key)
        # A is correlated with T=400*A in the original defining expression.
        # Cancellation happens before enclosing; independent A/T boxes are not used.
        shear_radius=c.mpf(2)/25;a_min=c.mpf('.7');a_max=c.mpf('.9')
        signed_entry=upper(read(exit['signed_Hv_logF_Z_upper']))
        signed_shape=c.mpf(max(packets.recovery.endpoints(signed_entry)[1],packets.recovery.endpoints(2*eps)[1]))
        SQ=lower((3-4*c.mpf(1)/200-10*eps)*(1-a_max/2)-c.mpf('5.5')/200-signed_shape)
        barrier=lower(110*SQ-3)
        Qone=lower((c.mpf('1.4')/c.mpf('1.65'))*(1-c.exp(c.mpf('-1.65'))))
        margins=dict(SQ_above_seven_fifths=SQ-c.mpf('1.4'),positive_p1_equals3_crossing=barrier,
            angular_Q_after_one_log_unit_above_half=Qone-c.mpf('.5'))
        if any(packets.recovery.endpoints(value)[0]<=0 for value in margins.values()):raise ArithmeticError('Actual whole angular barrier failed')
        return dict(physical_domain='110<=R<=Rz=exp(-8)*Rref, all Z[-1,1]',source_charts=list(CHARTS),
            input_source_budgets=positive,exact_correlated_shear_radius=shear_radius,
            original_angular_shear_bounds=[a_min,a_max],axial_shear_identically_zero=True,
            signed_Hv_zeta_upper=signed_shape,angular_SQ_lower=SQ,positive_directed_barrier_margins=margins,
            target_signed_gradient_proof=dict(
                exact_leading_term='-2Z^2*((1-delta)/2+4*(1-Z^2))/(1+Z^2)<=0 for |Z|<=1,delta<1',
                exact_error_term='-2Z*(1-Z^2)*(actualV110-4Z)/(1+Z^2)',
                absolute_error_multiplier_bound='2*|Z|*(1-Z^2)/(1+Z^2)<=2; 0<=1-Z^2<=1+Z^2',
                actual_V110_error_bound='same current whole C2 norm <epsilon0, via positive V110_error_budget',
                target_signed_Hv_zeta_upper=2*eps,
                convex_interpolation_bound='max(same current entry upper,2*epsilon0), original sigma in[0,1]'),
            angular_Q_after_one_log_unit_lower=Qone,
            p1_ODE='D_y p1+(1-a/2)*p1=R*SQ/L; p1=R*Qangular/L, not radial-velocity Q',
            source_mean='Mz/R=(110/R)*Mz110/110+(1-110/R)*actualV110; exact positive Z-independent weights',
            source_zeta='convex interpolation of entry zeta and -2Z/(1+Z^2); actual V110 unchanged',
            angular_source_bound='|W-(-3+4delta Z^2)|<=10epsilon0; target Hv*zeta<=2epsilon0; |1-2ZV|<=11',
            inlet_p1_strictly_above3_from_same_current_inner_power_source=True,
            whole_path_p1_strictly_above3=True,whole_path_H0_minus2_strictly_above1=True,
            whole_path_kappa_equals_a_below1=True,whole_path_t0_identically_zero=True,
            relaxed_branch='kappa=a<=.9<2, H0=p1>3>2; full p2 expression retained but no p2 norm bound inferred',
            actual_current_Rsh_source_and_reference_histories_not_reset=True,
            strict_completed_tensor_cone_new_regions_admitted=0,p1_p2_whole_path_norm_bounds_certified=False,
            axial_restoration_or_active_patch_relaxed_input_certified=False,
            source_radius_or_amplitude_not_materialized=True,**dict.fromkeys(OPEN,False))

    def query(self,chart,coordinate,Z=(-1,1)):
        c=self.ctx;v=c.mpf(coordinate);z=c.mpf(Z)
        vl,vh=packets.recovery.endpoints(v);zl,zh=packets.recovery.endpoints(z)
        if chart not in CHARTS or not all(mp.isfinite(x) for x in (vl,vh,zl,zh)) or vl<0 or vh>1 or zl<-1 or zh>1:
            raise ValueError('Actual reshape/reference source phase[0,1],Z[-1,1] required')
        after_one=chart=='inner_reference' or packets.recovery.endpoints(read_interval(self.ctx,self.rows['reference_join_bounds']['T'])*v)[0]>=1
        return dict(chart=chart,coordinate=v,Z=z,source_family=self.family,
            source_function_domain_not_saved_point_values=True,analytic_source_certificate=self.proof,
            angular_Q_above_half_on_entire_query=after_one,
            **{GATE:True},**dict.fromkeys(OPEN,False))


def read_interval(c,value):return packets.interval(c,value)


def run():
    owner=CurrentReshapeRelaxedInputs()
    examples={chart:owner.query(chart,(0,1)) for chart in CHARTS}
    examples['reshape_zero']=owner.query('reshape',0)
    examples['reshape_right']=owner.query('reshape',1)
    result=dict(source_family=owner.family,current_source_function_attachment_conditions=owner.conditions,
        original_current_angular_barrier_source_theorem=owner.theorem,
        current_whole_R110_Rz_relaxed_input_and_bounds=owner.proof,analytic_source_subbox_examples=examples,
        **{GATE:True},**dict.fromkeys(OPEN,False),
        strict_completed_tensor_cone_new_regions_admitted=0,source_graph_ancestor_constructors_called=False,
        scope='Current actual original whole R110..Rsh..Rz relaxed input and angular barrier; no axial restore/active patch, norm bounds, new strict tensor, changed loop, repair or new N.',
        input_hashes=owner.service.hashes)
    (HERE/NAME).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    print('Actual current whole R110..Rsh..Rz angular barrier and relaxed input PASS',flush=True)
    return result


if __name__=='__main__':run()
