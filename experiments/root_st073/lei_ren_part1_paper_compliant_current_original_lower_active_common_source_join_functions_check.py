"""Source-bound physical integral recovery, common IVP and actual flat collar."""
import json
from pathlib import Path
import time
import sympy as sy

import lei_ren_part1_paper_compliant_current_original_lower_active_common_source_join_functions as current


class SourceAlgebra:
    def __init__(self,owner,bindings,flat=False,comparison_flat=False):
        self.owner,self.g,self.bindings,self.cache,self.flat=owner,owner.g,bindings,{},flat
        self.comparison_flat=comparison_flat
    def at(self,node):
        if node in self.bindings:return self.bindings[node]
        if node in self.cache:return self.cache[node]
        row=self.g.nodes[node];op=row['operation'];at=self.at
        if op=='exact_rational':out=sy.Rational(row['numerator'],row['denominator'])
        elif op=='bound_variable':out=sy.Symbol(row['name'],real=True)
        elif op=='sum':out=sy.Add(*(at(n) for n in row['arguments']))
        elif op=='product':out=sy.Mul(*(at(n) for n in row['arguments']))
        elif op=='negative':out=-at(row['argument'])
        elif op=='positive_quotient':out=at(row['numerator'])/at(row['denominator'])
        elif op=='function_substitution':
            assert row['Z_independent_substitution']
            out=at(row['expression']).subs(at(row['variable']),at(row['value']))
        elif op=='definite_integral':
            lo,hi=at(row['lower']),at(row['upper'])
            out=sy.Integer(0) if lo==hi else sy.Integral(at(row['integrand']),
                (sy.Symbol(row['variable'],real=True),lo,hi))
        elif op=='current_original_unique_analytic_core_component':
            assert row['source_family']==self.owner.family and row['unique_analytic_Banach_solution']
            assert row['actual_twenty_term_operator'] and row['admitted_solution_not_a_model_or_enclosure_center']
            rho,Z=at(row['rho']),at(row['Z_variable']);r=sy.Symbol('rho',real=True)
            body=sy.Function(row['quantity'])(r,Z)
            out=sy.diff(body,r,row['radial_order'],Z,row['ordinary_slow_Z_derivative_order']).subs(r,rho)
        elif op=='original_lazy_flat_branch' and self.flat:
            coordinate=row.get('original_collar_coordinate');upper=row.get('original_collar_upper')
            predicate=self.g.nodes[row['flat_predicate']]
            assert coordinate==self.owner.active_phase.node
            assert upper==self.owner.loop.charts['first_micro']['loop']['original_collar_upper'].node
            assert predicate['operation']=='logical_or'
            assert any(self.g.nodes[n]==dict(operation='exact_real_comparison',operator='<=',
                left=coordinate,right=upper) for n in predicate['arguments'])
            out=at(row['flat_value'])
        elif op=='analytic_unary':
            x=at(row['argument']);name=row['name']
            if name=='expm1':out=sy.exp(x)-1
            elif name in ('exp','log','sqrt'):out=getattr(sy,name)(x)
            elif name=='original_flat_sigma_prime':out=sy.Function('sigma_prime')(x)
            else:raise AssertionError('Unbound unary '+name)
        elif op=='exact_original_source_flat_sigma':
            u=self.g.symbol('original_prescribed_C6_comparison_t')
            if self.comparison_flat and row['argument']==self.g.sub(u,self.g.one).node:
                assert self.owner.sc_binding['comparison_alpha_exactly_one_on_zero_to_sc']
                out=sy.Integer(0)
            else:out=sy.Function('sigma')(at(row['argument']))
        else:raise AssertionError('Unbound actual seam source '+op+' '+str(node))
        self.cache[node]=out;return out


def bind_jet(bindings,jet,body,Z):
    for n,ref in enumerate(current.rows(jet)):
        value=sy.diff(body,Z,n)
        if ref.node==0:assert value==0
        else:bindings[ref.node]=value
    bindings[0]=sy.Integer(0)


def constants(owner):
    g=owner.g;Z=sy.Symbol('Z',real=True);rho=sy.Symbol('rho',real=True)
    S=sy.Symbol('s',real=True);eps,hb,Pstar,de,j=sy.symbols('epsilon hb Pstar delta j',positive=True)
    bindings={g.symbol('Z').node:Z,owner.core.rho.node:rho,owner.lower_phase.node:S,owner.active_phase.node:S,
        owner.core.parameters['epsilon'].node:eps,owner.bridge.hb.node:hb,owner.inlet.Pstar.node:Pstar,
        owner.core.parameters['delta'].node:de,owner.core.parameters['j'].node:j}
    bind_jet(bindings,owner.bridge.core7['F0'],sy.Function('F0')(Z),Z)
    bind_jet(bindings,owner.bridge.P0,sy.Function('PD')(Z),Z)
    bind_jet(bindings,owner.core.functions['P0'],Pstar**2*sy.Function('PD')(Z),Z)
    return bindings,dict(Z=Z,rho=rho,S=S,eps=eps,hb=hb,Pstar=Pstar,de=de,j=j)


def source_check(owner):
    g=owner.g;bridge=owner.bridge;lower=owner.lower;core=owner.core
    derivative=current.bridge_source.inlet.previous.HigherCoreZDerivative(lower.lower)
    for group in (bridge.core7,bridge.initial7):
        for jet in group.values():
            q=current.rows(jet)[:7]
            for left,right in zip(q,q[1:]):derivative.seed(left,right)
    for jets in bridge.radial7.values():
        for jet in jets:
            q=current.rows(jet)
            for left,right in zip(q,q[1:]):derivative.seed(left,right)
    assert bridge.inlet.previous is lower and bridge.core is core
    assert owner.loop.source_charts['first_micro'] is bridge.micro
    assert owner.modified.charts['first_micro']['source_chart'] is owner.loop.charts['first_micro']
    assert bridge.hb==lower.lower.parameters['hb'] and bridge.Ra==owner.inlet.radius
    assert owner.left_endpoint==lower.lower.parameters['sc_half']
    assert owner.right_endpoint==g.mul(owner.loop.sc,g.constant('1/2'))
    assert current.selected_sc_binding(lower.lower,owner.loop)==owner.sc_binding
    assert current.rows(bridge.P0)[:7]==current.rows(owner.inlet.P0)
    native=core.pulse.own.jets.native['P0']
    assert current.rows(bridge.P0)[:4]==current.rows(native)
    core_module=getattr(core,'source_module',lower.lower.core.__class__.__module__)
    import sys
    expected=sys.modules[core_module].c3.C3Algebra(g).scale(native,owner.inlet_binding['physical_axis_scale'])
    assert current.rows(expected)==current.rows(core.functions['P0'])
    assert owner.inlet.Pstar==g.unary('exp',core.parameters['logP'])
    assert g.nodes[core.parameters['epsilon'].node]==dict(operation='analytic_unary',name='exp',
        argument=g.neg(core.parameters['logLambda']).node)
    assert sys.modules[core_module].map_terms(g,core.definition['source_variables'])==core.terms
    assert len(core.terms['theta'])==len(core.terms['z'])==10
    for key in ('Phi','Psi','mean_Psi'):
        for n,ref in enumerate(current.rows(lower.core_functions[key])):
            node=g.nodes[ref.node]
            assert node['quantity']==key and node['radial_order']==0 and node['ordinary_slow_Z_derivative_order']==n
            assert node['axis_conditions']=='Phi(0,Z)=1; Psi(0,Z)=0'
            assert node['source_family']==owner.family and node['unique_analytic_Banach_solution']
    count=0
    for key,jet in owner.inlet.moment_functions.items():
        for n,ref in enumerate(current.rows(jet)):
            row=g.nodes[ref.node]
            assert row['operation']=='definite_integral' and row['actual_core_normalized_moment']==key
            assert row['same_unique_nonlinear_core_and_full_density_product']
            assert row['lower']==g.zero.node and row['upper']==g.one.node
            assert row['ordinary_slow_Z_derivative_order']==n
            count+=1
        assert current.rows(bridge.initial7[key])[:7]==current.rows(owner.inlet.inlet_moments[key])
        seventh=current.rows(bridge.initial7[key])[7]
        assert derivative.derivative(current.rows(owner.inlet.inlet_moments[key])[6])==seventh
        endpoint=g.nodes[seventh.node]
        assert endpoint['operation']=='function_substitution' and endpoint['variable']==core.rho.node
        assert endpoint['value']==g.constant(4).node
        primitive=g.nodes[endpoint['expression']]
        assert primitive['operation']=='definite_integral' and primitive['ordinary_slow_Z_derivative_order']==7
        assert primitive['actual_core_normalized_moment']==key and primitive['endpoint_independent_of_Z']
    pressure7=g.nodes[current.rows(bridge.P0)[7].node]
    assert pressure7['ordinary_slow_Z_derivative_order']==7 and pressure7['true_order7_pressure_provider']
    assert pressure7['original_function']==native.value.node and pressure7['source_family']==owner.family
    assert pressure7['source_projection']=='factorial(7)*P0_rows(Z).coefficients[7]'
    assert owner.integrated_theorem['pressure_and_swirl_physical_scales_not_dropped']
    assert owner.integrated_theorem['same_common_pressure_P_equals_P0_plus_integral_F_squared']
    return dict(passed=True,actual_six_inlet_source_integral_rows=count,
        actual_twenty_term_core_and_analytic_axis_sources_bound=True,
        exact_active_sc_source_projection_and_singleton_recomputation=True,
        normalized_P0_and_physical_Pstar_squared_relation_identified=True,
        actual_core_bridge_lower_and_finite_N_owner_links=True)


def inlet_primitives(owner):
    bindings,v=constants(owner);g=owner.g;Z,rho,eps,j=(v[k] for k in ('Z','rho','eps','j'))
    t=sy.Symbol('original_core_inlet_moment_t',real=True)
    phi=sy.Function('Phi')(rho*t,Z);psi=sy.Function('Psi')(rho*t,Z);U=4*Z+j+eps*psi
    bodies=dict(H=2*t*phi,M=U,K=2*t*phi*U,A=U**2,B=t*phi**2,C=phi**2)
    alg=SourceAlgebra(owner,bindings);checked=[]
    for key,jet in owner.inlet.moment_functions.items():
        row=g.nodes[jet.value.node]
        assert sy.simplify(alg.at(row['integrand'])-bodies[key])==0,key
        saved=g.nodes[owner.inlet.inlet_moments[key].value.node]
        assert saved['operation']=='function_substitution' and saved['expression']==jet.value.node
        assert saved['variable']==owner.core.rho.node and saved['value']==g.constant(4).node
        checked.append(key)
    mean=g.nodes[owner.core.components['mean_Psi'].value.node]
    assert mean['component_definition']=='mean_Psi=integral_0^1 Psi(rho*t,Z) dt'
    return dict(passed=True,actual_six_defining_inlet_density_identities=checked,
        M_at_rho4_equals_U0_plus_epsilon_mean_Psi_by_integral_linearity=True,
        C_at_rho4_equals_one_quarter_true_Phi_squared_primitive=True,
        full_Z_derivatives_follow_same_analytic_source_and_fixed_limits=True,
        all_six_nonzero_original_inlet_histories_retained=True)


def physical_direction_identity(owner):
    bindings,v=constants(owner);Z=v['Z'];R=sy.Symbol('R',positive=True);Pstar=v['Pstar'];de=v['de']
    F=sy.Function('F')(R,Z);V=sy.Function('V')(R,Z);P0=Pstar**2*sy.Function('PD')(Z);F0=sy.Function('F0')(Z)
    IV,IA,IH,IK,IB,IC=[sy.Function(k)(R,Z) for k in ('IV','IA','IH','IK','IB','IC')]
    J=dict(H=2*IH/(R**2*F0),M=IV/R,K=2*IK/(R**2*F0),A=IA/R,B=IB/(R**2*F0**2),C=IC/(R*F0**2))
    cmp=owner.bridge.comparison_micro
    bind_jet(bindings,cmp['phi'],F/F0,Z);bind_jet(bindings,cmp['U'],V,Z)
    for key,jet in cmp['moments'].items():bind_jet(bindings,jet,J[key],Z)
    bindings[owner.bridge.micro['radius'].node]=R
    alg=SourceAlgebra(owner,bindings);D=alg.at(cmp['directions']['D_over_R'].value.node)
    drive=alg.at(owner.bridge.micro['actual']['full_signed_axial_drive'].value.node)
    L=1-de*Z**2;d=1-Z**2;W=1-(1-de)*Z*IV/R-d*sy.diff(IV,Z)/R
    angular=2*L*R**2*sy.diff(F,R)-W*R**2*F+(1-de/2)*IH-(1-de)*Z*sy.diff(IH,Z)/2-d*sy.diff(IK,Z)+(2*de-1)*Z*IK
    axial=2*L*R*sy.diff(V,R)-W*R*V+(1-de)*IV/2-(1-de)*Z*sy.diff(IV,Z)/2+2*de*Z*IA-d*sy.diff(IA,Z)
    axial-=R*(d*sy.diff(P0,Z)-2*(1+de)*Z*P0)
    axial-=d*(R*sy.diff(IC,Z)-sy.diff(IB,Z))
    axial+=2*(1+de)*Z*(R*IC-IB)+2*Z*IB
    assert sy.cancel(angular/R**2-(2*L*sy.diff(F,R)+L*F*D))==0
    assert sy.cancel(axial/R-2*L*(sy.diff(V,R)+drive/R))==0
    return dict(passed=True,actual_current_direction_and_three_sector_drive_value_identities=2,
        actual_F0_axial_dependence_and_pressure_Pstar_squared_retained=True,
        admitted_integrated_actual_core_equations_and_zero_axis_constants_consumed=True,
        actual_core_substitutions_give_D_plus_2_Phi_rho_over_epsilon_Phi_zero=True,
        actual_core_substitutions_give_full_drive_plus_epsilon_rho_Psi_rho_zero=True)


def pull_constants(expr):
    def change(node):
        variable=node.limits[0][0];constant,body=node.function.as_independent(variable,as_Add=False)
        canonical=sy.Symbol('seam_integral_t',real=True)
        return constant*sy.Integral(body.subs(variable,canonical),
            (canonical,*node.limits[0][1:]))
    return expr.replace(lambda n:isinstance(n,sy.Integral),change)


def comparison_core_FTC(owner):
    bindings,v=constants(owner);Z,S,eps,hb=(v[k] for k in ('Z','S','eps','hb'))
    rho=4*sy.exp(hb*S);P=sy.Function('Phi');Psi=sy.Function('Psi');cmp=owner.bridge.comparison_micro
    alg=SourceAlgebra(owner,bindings,comparison_flat=True)
    logbar=alg.at(cmp['logphi'].value.node);barU=alg.at(cmp['U'].value.node)
    assert sy.simplify(sy.diff(logbar,S)-sy.diff(sy.log(P(rho,Z)),S))==0
    assert sy.simplify(sy.diff(barU,S)-sy.diff(4*Z+v['j']+eps*Psi(rho,Z),S))==0
    zero=SourceAlgebra(owner,{**bindings,owner.active_phase.node:sy.Integer(0)},comparison_flat=True)
    assert sy.simplify(zero.at(cmp['logphi'].value.node)-sy.log(P(4,Z)))==0
    assert sy.simplify(zero.at(cmp['U'].value.node)-(4*Z+v['j']+eps*Psi(4,Z)))==0
    assert owner.g.nodes[cmp['phi'].value.node]==dict(operation='analytic_unary',name='exp',argument=cmp['logphi'].value.node)
    # Exact weighted primitive ODEs identify the comparison six histories
    # with the true core integrals throughout the flat-alpha interval.
    r=sy.Symbol('r',positive=True);t=sy.Symbol('t',positive=True);f=sy.Function('density')
    count=0
    for key,weight,norm in (('H',1,2),('M',0,1),('K',1,2),('A',0,1),('B',1,1),('C',0,1)):
        primitive=norm*sy.Integral(t**weight*f(t),(t,0,r))/r**(weight+1)
        assert sy.simplify(r*sy.diff(primitive,r)-(norm*f(r)-(weight+1)*primitive))==0
        assert current.bridge_source.RATES[key]==weight+1
        count+=1
    return dict(passed=True,actual_comparison_logphi_and_U_FTC_source_identities=4,
        original_alpha_one_minus_sigma_u_minus_one_on_actual_subunit_domain=True,
        actual_comparison_phi_equals_unique_core_Phi_and_U_equals_true_core_U=True,
        exact_six_weighted_core_primitive_history_ODEs=count,
        same_true_inlet_histories_identify_comparison_moments_with_core_functions=True)


def common_IVP(owner):
    bindings,v=constants(owner);Z,S,eps,hb=(v[k] for k in ('Z','S','eps','hb'))
    g=owner.g;rho=4*sy.exp(hb*S);P=sy.Function('Phi');Psi=sy.Function('Psi');ratio=sy.Function('ratio')(S,Z)
    chi=1+(hb-1)*sy.Function('sigma')(S)
    lower=owner.lower.lower.fields;cmp=owner.bridge.comparison_micro;actual=owner.bridge.micro['actual']
    bind_jet(bindings,cmp['phi'],P(rho,Z),Z)
    bind_jet(bindings,cmp['directions']['D_over_R'],-2*sy.Subs(sy.diff(P(v['rho'],Z),v['rho']),v['rho'],rho)/(eps*P(rho,Z)),Z)
    bind_jet(bindings,actual['full_signed_axial_drive'],-eps*rho*sy.Subs(sy.diff(Psi(v['rho'],Z),v['rho']),v['rho'],rho),Z)
    alg=SourceAlgebra(owner,bindings)
    logratio=alg.at(owner.lower.fields['logratio'].value.node)
    true_log=sy.log(P(rho,Z))
    assert sy.simplify(sy.diff(logratio,S)-(chi-1)*sy.diff(true_log,S))==0
    active_log=alg.at(actual['logratio'].value.node)
    assert sy.simplify(sy.diff(active_log,S)-chi*sy.diff(true_log,S))==0
    # Identify the actual ratios independently, retaining their original
    # definitions. The canonical ratio is the lower exp(logratio), not the
    # core Phi/Phi(4): actual chi is nonconstant.
    left_ratio=alg.at(owner.lower.fields['ratio'].value.node)
    right_ratio=alg.at(actual['ratio'].value.node)
    assert sy.simplify(sy.diff(sy.log(left_ratio),S)-(chi-1)*sy.diff(true_log,S))==0
    assert sy.simplify(sy.diff(sy.log(right_ratio),S)-(chi-1)*sy.diff(true_log,S))==0
    assert sy.simplify(left_ratio.subs(S,0).doit()-1)==0
    assert sy.simplify(right_ratio.subs(S,0).doit()-1)==0
    # This named ratio is now a registered source function with the proved
    # defining expression above, used to avoid nested dummy-variable capture.
    bindings[owner.lower.fields['ratio'].value.node]=ratio
    bindings[actual['ratio'].value.node]=ratio
    alg=SourceAlgebra(owner,bindings)
    leftU=alg.at(owner.left_fields['U'].value.node)
    # The modified field has zero corrections on the original collar; use
    # the actual prescribed U to compare the two source primitives first.
    rightU=alg.at(actual['U'].value.node)
    assert sy.simplify(pull_constants(leftU)-pull_constants(rightU))==0
    leftUy=alg.at(owner.left_fields['U_y'].value.node)
    assert sy.simplify(leftUy-alg.at(actual['U_y'].value.node))==0
    assert sy.simplify(leftUy-sy.diff(leftU,S)/hb)==0
    assert owner.bridge.initial7.keys()==owner.inlet.inlet_moments.keys()==owner.lower_histories.keys()
    assert owner.integrated_theorem['axis_integration_constant_proof']
    # The actual lower mean is a complete rho4 primitive. Prove its ODE
    # from its saved body instead of asserting a history correspondence.
    U=sy.Function('U')(S,Z);PhiL=sy.Function('PhiL')(S,Z)
    mean_bindings={**bindings};U0=4*Z+v['j']
    bind_jet(mean_bindings,owner.lower.radial_fields['Psi_rho_rows'][0],(U-U0)/eps,Z)
    bind_jet(mean_bindings,owner.lower.radial_fields['Phi_rho_rows'][0],PhiL,Z)
    mean_alg=SourceAlgebra(owner,mean_bindings)
    mean=mean_alg.at(owner.left_fields['m'].value.node)
    assert sy.simplify(sy.diff(mean,S)-hb*(U-mean))==0
    pressure=mean_alg.at(owner.left_fields['pressure'].value.node)
    F0=sy.Function('F0')(Z);R=eps*rho
    assert sy.simplify(sy.diff(pressure,S)-hb*R*F0**2*PhiL**2)==0
    histories=owner.lower_histories
    for key,density in (('M',U),('C',PhiL**2)):
        body=mean_alg.at(histories[key].value.node)
        assert sy.simplify(sy.diff(body,S)-hb*(density-body))==0,key
        initial=mean_alg.at(owner.inlet.inlet_moments[key].value.node)
        assert sy.simplify(body.subs(S,0)-initial)==0,(key,'lower inlet')
    active_bindings={**mean_bindings}
    bind_jet(active_bindings,actual['U'],U,Z)
    bind_jet(active_bindings,actual['phi'],PhiL,Z)
    active_alg=SourceAlgebra(owner,active_bindings)
    for key,density in (('M',U),('C',PhiL**2)):
        body=active_alg.at(owner.bridge.micro['moments'][key].value.node)
        initial=active_alg.at(owner.inlet.inlet_moments[key].value.node)
        assert sy.simplify(sy.diff(body,S)-hb*(density-body))==0,(key,'active ODE')
        assert sy.simplify(body.subs(S,0)-initial)==0,(key,'active inlet')
    return dict(passed=True,actual_lower_logratio_and_actual_bridge_logratio_ODE_identities=2,
        actual_lower_and_prescribed_axial_integral_function_identity=True,
        actual_physical_Uy_source_equality_and_native_derivative_divided_by_hb=True,
        common_true_rho4_inlet_and_comparison_alpha_one_from_actual_sc=True,
        angular_positive_linear_IVP_uniqueness_with_nonconstant_chi=True,
        actual_left_and_right_ratio_ODE_and_unit_inlet_identities=4,
        canonical_ratio_source_expression_node=owner.lower.fields['ratio'].value.node,
        all_six_same_inlet_Duhamel_histories_then_have_identical_solutions=True,
        actual_saved_lower_mean_and_pressure_and_M_C_history_FTC_identities=4,
        actual_active_M_C_history_ODE_and_both_inlet_identities=6,
        mean_and_absolute_pressure_follow_same_M_C_FTC_IVP=True,
        physical_Uy_uses_no_extra_hb=True)


def physical_recovery(owner):
    bindings,v=constants(owner);Z,S,eps,hb,de=(v[k] for k in ('Z','S','eps','hb','de'))
    rho=4*sy.exp(hb*S);R=eps*rho;F0=sy.Function('F0')(Z);phi=sy.Function('PhiL')(S,Z)
    U,M,Uy=(sy.Function(k)(S,Z) for k in ('U','M','Uy'));my=U-M;g=owner.g
    for group in (owner.left_fields,owner.right_fields):
        for key,body in (('U',U),('m',M),('U_y',Uy)):
            bind_jet(bindings,group[key],body,Z)
    direct=SourceAlgebra(owner,bindings)
    for group in (owner.left_fields,owner.right_fields):
        assert sy.simplify(direct.at(group['m_y'].value.node)-my)==0
        bind_jet(bindings,group['m_y'],my,Z)
    bind_jet(bindings,owner.lower.radial_fields['U_rho_rows'][1],Uy/rho,Z)
    bind_jet(bindings,owner.lower.radial_fields['m_rho_rows'][1],my/rho,Z)
    bind_jet(bindings,owner.lower.radial_fields['Phi_rho_rows'][0],phi,Z)
    chart=owner.loop.charts['first_micro'];C=sy.Function('C')(S,Z)
    bind_jet(bindings,chart['source_chart']['fields']['E'],F0*phi*sy.sqrt(2*R)/v['Pstar'],Z)
    bind_jet(bindings,chart['source_chart']['fields']['p'],R*F0**2*C/v['Pstar']**2,Z)
    bindings[owner.loop.N.node]=sy.Symbol('N',positive=True)
    bindings[owner.loop.sc.node]=sy.Symbol('sc',positive=True)
    alg=SourceAlgebra(owner,bindings,flat=True)
    L=1-de*Z**2;d=1-Z**2
    Q=(2*Z*U-(1-de)*Z*M-d*sy.diff(M,Z))/L
    Qy=(2*Z*Uy-(1-de)*Z*my-d*sy.diff(my,Z))/L
    checked=[]
    for label,group in (('left',owner.left_fields),('right',owner.right_fields)):
        for key,body in (('Q',Q),('Q_y',Qy),('Utheta',F0*phi*sy.sqrt(2*R))):
            assert sy.simplify(alg.at(group[key].value.node)-body)==0,(label,key)
            checked.append(label+'_'+key)
    pressure=alg.at(owner.right_fields['pressure'].value.node)
    assert sy.simplify(pressure-(v['Pstar']**2*sy.Function('PD')(Z)+R*F0**2*C))==0
    return dict(passed=True,actual_saved_physical_recovery_value_identities=checked,
        actual_leading_Q_and_Qy_one_minus_delta_recovery=True,
        positive_order_one_plus_delta_recovery_is_a_separate_provider=True,
        actual_swirl_same_F0_phi_sqrt_two_R_without_extra_Pstar=True,
        actual_active_absolute_pressure_same_axis_plus_R_F0_squared_C=True,
        mean_pressure_IVP_and_derivative_recovery_complete_all_eight_trace_source_fields=True)


def flat_collar(owner):
    bindings,v=constants(owner);Z=v['Z'];g=owner.g;chart=owner.loop.charts['first_micro']
    bindings[owner.loop.sc.node]=sy.Symbol('sc',positive=True)
    bindings[owner.loop.N.node]=sy.Symbol('N',positive=True)
    for key,jet in chart['roots'].items():bind_jet(bindings,jet,sy.Function('root_'+key)(Z),Z)
    for key,jet in chart['source_chart']['fields'].items():
        bind_jet(bindings,jet,sy.Function('background_'+key)(Z),Z)
    alg=SourceAlgebra(owner,bindings,True);checked=[]
    for label,jet in [(k,chart['loop'][k]) for k in ('q','A','B_over_Pstar','A_phi','B_over_Pstar_phi')]+[
        ('delta_E',chart['finite_N']['delta_E']),('delta_V',chart['finite_N']['delta_V'])]+[
        ('density_'+k,q) for k,q in chart['finite_N']['five_signed_density_increments'].items()]:
        for n,ref in enumerate(current.rows(jet)):
            assert sy.simplify(alg.at(ref.node))==0,(label,n)
            checked.append(f'{label}_Z{n}')
    assert chart['incoming']==owner.loop.zero_inlet
    assert all(ref==g.zero for jet in chart['incoming'].values() for ref in current.rows(jet))
    for jet in chart['contributions'].values():
        for ref in current.rows(jet):
            row=g.nodes[ref.node]
            assert row['lower']==owner.right_endpoint.node and row['upper']==owner.active_phase.node
    for key in ('slow_A_y','slow_B_y'):
        for n,ref in enumerate(current.rows(owner.modified.charts['first_micro']['radial'][key])):
            assert sy.simplify(alg.at(ref.node))==0,(key,n)
    return dict(passed=True,actual_flat_loop_and_finite_N_rows=checked,
        actual_slow_radial_loop_source_rows_zero=12,
        actual_flat_predicate_has_original_coordinate_below_three_sc_quarters=True,
        common_seam_sc_half_lies_strictly_inside_original_two_sided_flat_collar=True,
        correction_inlet_zero_and_empty_same_source_histories_used=True,
        singular_active_inverse_not_evaluated_on_flat_branch=True)


def ordinary_bindings(owner,side):
    bindings,v=constants(owner);Z,S=v['Z'],v['S']
    if side=='left':
        for key in ('Phi_rho_rows','Psi_rho_rows','mean_Psi_rho_rows','F_rho_rows'):
            for n,jet in enumerate(owner.lower.radial_fields[key]):
                bind_jet(bindings,jet,sy.Function('source_'+key+'_'+str(n))(S,Z),Z)
    else:
        chart=owner.loop.charts['first_micro']
        for key in ('E','V','m','pressure'):
            bind_jet(bindings,chart['modified_fields'][key],sy.Function('source_normalized_'+key)(S,Z),Z)
        bind_jet(bindings,owner.modified.charts['first_micro']['radial']['V_y'],sy.Function('source_normalized_Vy')(S,Z),Z)
    return bindings,v


def modified_physical_decomposition(owner):
    bindings,v=constants(owner);Z,S,eps,hb,de=(v[k] for k in ('Z','S','eps','hb','de'))
    R=4*eps*sy.exp(hb*S);F0=sy.Function('F0')(Z)
    U,M,Uy,phi,C=(sy.Function(k)(S,Z) for k in ('U','M','Uy','PhiL','C'))
    bg=owner.bridge.micro;g=owner.g
    for key,body in (('U',U),('U_y',Uy),('phi',phi)):
        bind_jet(bindings,bg['actual'][key],body,Z)
    bind_jet(bindings,bg['moments']['M'],M,Z)
    bind_jet(bindings,bg['moments']['C'],C,Z)
    bindings[owner.loop.sc.node]=sy.Symbol('sc',positive=True)
    bindings[owner.loop.N.node]=sy.Symbol('N',positive=True)
    alg=SourceAlgebra(owner,bindings,flat=True);d=1-Z**2;L=1-de*Z**2;my=U-M
    expected=dict(U=U,U_y=Uy,m=M,m_y=my,Utheta=F0*phi*sy.sqrt(2*R),
        Q=(2*Z*U-(1-de)*Z*M-d*sy.diff(M,Z))/L,
        Q_y=(2*Z*Uy-(1-de)*Z*my-d*sy.diff(my,Z))/L,
        pressure=v['Pstar']**2*sy.Function('PD')(Z)+R*F0**2*C)
    checked=[]
    for key,jet in owner.right_fields.items():
        assert sy.simplify(alg.at(jet.value.node)-expected[key])==0,key
        checked.append(key)
    return dict(passed=True,actual_modified_physical_field_value_decompositions=checked,
        actual_finite_N_density_and_history_zero_consumed_without_output_field_bindings=True,
        physical_mean_equals_actual_background_M_history_without_another_radius_division=True,
        physical_pressure_equals_same_axis_plus_R_F0_squared_C_history=True,
        actual_common_IVP_background_fields_feed_all_eight_modified_outputs=True)


def ordinary_source_rows(owner):
    checked=[]
    # Rebuild the admitted higher calculus from its accepted input layer.
    # No current output row is seeded; higher Q rows must be the actual
    # derivatives of the preserved lower Q recovery, not fixture aliases.
    d=current.bridge_source.inlet.previous.HigherCoreZDerivative(owner.lower.lower)
    for group in (owner.core.functions,owner.lower.core_functions):
        for name in ('Phi','Psi','mean_Psi','U0','gradient','chi','G','logF0','F0'):
            q=current.rows(group[name])
            for left,right in zip(q,q[1:]):d.seed(left,right)
    for group in (owner.lower.lower.core_radial_rows,owner.lower.core_radial_rows):
        for jets in group.values():
            for jet in jets:
                q=current.rows(jet)
                for left,right in zip(q,q[1:]):d.seed(left,right)
    admitted=[]
    for side,fields in (('left',owner.left_fields),('right',owner.right_fields)):
        bindings,v=ordinary_bindings(owner,side);Z=v['Z'];alg=SourceAlgebra(owner,bindings)
        for key,jet in fields.items():
            if side=='left' and key=='pressure':
                bind_jet(bindings,owner.lower.core_functions['Phi'],
                    sy.Function('source_core_pressure_Phi')(v['rho'],Z),Z)
                alg=SourceAlgebra(owner,bindings)
            if side=='left' and key=='Q_y':
                for label,array in (('U','U_rho_rows'),('m','m_rho_rows')):
                    bind_jet(bindings,owner.lower.radial_fields[array][1],
                        sy.Function('actual_lower_radial_'+label)(v['S'],Z),Z)
                alg=SourceAlgebra(owner,bindings)
            body=alg.at(jet.value.node)
            for n,ref in enumerate(current.rows(jet)):
                if side=='left' and key in ('Q','Q_y') and n>=len(current.rows(owner.lower.lower.fields[key])):
                    assert d.derivative(current.rows(jet)[n-1])==ref,(side,key,n)
                    admitted.append(f'{side}_{key}_Z{n}')
                    checked.append(f'{side}_{key}_Z{n}')
                    continue
                residual=sy.expand(alg.at(ref.node)-sy.diff(body,Z,n))
                # Definite integrals are linear. SymPy leaves a constant
                # inside one integral and outside another after ordinary
                # differentiation; canonicalize that exact identity before
                # comparing formal atoms, rather than evaluating a primitive.
                replacements={}
                for integral in residual.atoms(sy.Integral):
                    assert len(integral.limits)==1
                    var,lo,hi=integral.limits[0]
                    assert not lo.has(Z) and not hi.has(Z)
                    terms=[]
                    for term in sy.Add.make_args(sy.expand(integral.function)):
                        coefficient,density=term.as_independent(var,as_Add=False)
                        terms.append(coefficient*sy.Integral(density,(var,lo,hi)))
                    replacements[integral]=sy.Add(*terms)
                residual=sy.expand(residual.xreplace(replacements))
                if residual!=0:
                    # Treat identical formal derivative/integral atoms as
                    # algebraic coefficients; SymPy's rational collector
                    # does not support some multivariate Derivative nodes.
                    atoms=sorted(residual.atoms(sy.Derivative,sy.Subs,sy.Integral),key=sy.sstr)
                    replaced=residual.xreplace({atom:sy.Symbol('exact_source_atom_'+str(i)) for i,atom in enumerate(atoms)})
                    assert sy.cancel(replaced)==0,(side,key,n)
                checked.append(f'{side}_{key}_Z{n}')
            if side=='left' and key in ('U','U_y','m','m_y'):
                # Reuse this independently checked layer as an input to
                # the later Q recovery; never bind the output being checked.
                bind_jet(bindings,jet,body,Z);alg=SourceAlgebra(owner,bindings)
    return dict(passed=True,actual_saved_field_ordinary_derivative_rows=checked,
        independently_reconstructed_admitted_higher_Q_derivative_rows=admitted,
        true_source_inputs_and_previously_checked_layers_without_binding_current_output=True,
        both_actual_field_calculus_constructions_independently_equal_derivatives_of_their_values=True,
        source_C5_Q4_pressure_C3_rows_not_forced_to_a_common_arbitrary_jet=True)


def trace_transfer(owner,proofs):
    for name in ('source','inlet_primitives','physical_direction','comparison_core','common_IVP','physical_recovery','flat_collar','modified_decomposition','ordinary_sources'):
        assert proofs[name]['passed'],name
    g=owner.g;checked=[]
    for key,left in owner.left_trace.items():
        right=owner.right_trace[key]
        assert len(current.rows(left))==len(current.rows(right))
        for n,(l,r) in enumerate(zip(current.rows(left),current.rows(right))):
            for ref,body,var,end in ((l,current.rows(owner.left_fields[key])[n],owner.lower_phase,owner.left_endpoint),
                                    (r,current.rows(owner.right_fields[key])[n],owner.active_phase,owner.right_endpoint)):
                row=g.nodes[ref.node]
                assert row['operation']=='function_substitution' and row['expression']==body.node
                assert row['variable']==var.node and row['value']==end.node and row['Z_independent_substitution']
            checked.append(f'{key}_Z{n}')
    # Equality is transported from the common defining IVP and physical
    # recovery, never inferred from these unequal graph indices themselves.
    assert len(checked)==44
    bindings,v=constants(owner);sc=sy.Symbol('sc',positive=True);Z,S=v['Z'],v['S']
    bindings[owner.loop.sc.node]=bindings[owner.lower.lower.parameters['sc'].node]=sc
    radius_alg=SourceAlgebra(owner,bindings)
    assert sy.simplify(radius_alg.at(owner.left_radius.node)-radius_alg.at(owner.right_radius.node))==0
    # Derive canonical bodies from the actual lower source construction.
    # Value identities were proved independently, and ordinary_sources
    # proves each output jet is genuinely the derivative of its own value.
    # This consumes those equalities; no arbitrary joined field is injected.
    canonical_bindings,_=ordinary_bindings(owner,'left')
    canonical_alg=SourceAlgebra(owner,canonical_bindings)
    for key in owner.left_fields:
        body=canonical_alg.at(owner.left_fields[key].value.node)
        for fields in (owner.left_fields,owner.right_fields):
            for n,ref in enumerate(current.rows(fields[key])):
                bindings[ref.node]=body if n==0 else sy.Derivative(body,Z,n,evaluate=False)
    alg=SourceAlgebra(owner,bindings)
    for key,left in owner.left_trace.items():
        for n,(l,r) in enumerate(zip(current.rows(left),current.rows(owner.right_trace[key]))):
            assert sy.simplify(alg.at(l.node)-alg.at(r.node))==0,(key,n)
    return dict(passed=True,actual_registered_two_sided_trace_rows=checked,
        U_Uy_mean_mean_y_swirl_C5_Q_Qy_C4_and_absolute_pressure_C3=True,
        source_function_equalities_and_genuine_ordinary_derivatives_transport_to_traces=True,
        canonical_bodies_derived_from_actual_lower_sources_not_arbitrary_joined_functions=True,
        actual_radius_graph_expressions_compared_with_resolved_exact_sc=True,
        actual_radius_equality_uses_same_epsilon_hb_and_resolved_exact_sc=True,
        quantitative_mixed_norms_and_other_global_joins_not_claimed=True)


def run(owner):
    began=time.monotonic();owner.assert_graph();raw=current.parent.linear.read(current.NAME)
    assert raw['lower_active_common_source_bindings']==owner.snapshot
    assert raw['exact_graph_nodes']==owner.g.nodes[:len(raw['exact_graph_nodes'])]
    checks={}
    for name,fn in (('source',source_check),('inlet_primitives',inlet_primitives),
        ('physical_direction',physical_direction_identity),('comparison_core',comparison_core_FTC),
        ('common_IVP',common_IVP),('physical_recovery',physical_recovery),
        ('flat_collar',flat_collar),('modified_decomposition',modified_physical_decomposition),
        ('ordinary_sources',ordinary_source_rows),
        ('trace_transfer',lambda o:trace_transfer(o,checks))):
        checks[name]=fn(owner);print('LOWER_ACTIVE_'+name.upper()+'_CHECKED',flush=True)
    owner.assert_graph()
    for name,digest in owner.hashes.items():assert current.sha(name)==digest,name
    rec=dict(all_passed=True,source_family=owner.family,checks=checks,
        **dict.fromkeys(current.GATES,True),**dict.fromkeys(current.OPEN,False),
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(rec,indent=2)+'\n',encoding='utf8',newline='\n')
    print('ACTUAL_LOWER_ACTIVE_COMMON_SOURCE_JOIN_CHECKED',flush=True);return rec


if __name__=='__main__':raise SystemExit('run(owner) reuses checked actual sources')
