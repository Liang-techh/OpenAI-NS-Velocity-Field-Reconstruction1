"""Coordinate-changing source seam; construction checks do not admit n1.

The background seam is obtained from the original Duhamel equations with
their actual incoming nodes. No high output derivative is assigned a value.
"""
import json
from pathlib import Path
import time

import sympy as sy

import lei_ren_part1_paper_compliant_current_original_micro_macro_mixed_source_join as current
import lei_ren_part1_paper_compliant_current_original_micro_active_radial4_pressure_C3_functions_check as active
import lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions as bridge

rows,encode=current.outer.rows,current.outer.encode
zero,bind=active.zero,active.bind
CONSTRUCTION_RECEIPT=current.PREFIX+'current_original_micro_macro_mixed_source_join_construction_check.json'


class SeamAlgebra(active.prior.RminusAlgebra):
    """Lexical substitutions and integrals, without an active flat shortcut."""
    def at(self,i):
        i=i.node if hasattr(i,'node') else i
        if i in self.bindings or i in self.cache:return super().at(i)
        row=self.g.nodes[i];op=row['operation'];owner=self.rminus
        if op=='function_substitution':
            assert row['Z_independent_substitution'];scope=dict(self.scope)
            scope[row['variable']]=self.at(row['value'])
            result=SeamAlgebra(owner,self.source_bindings,self.evidence,scope).at(row['expression'])
        elif op=='definite_integral':
            assert row['endpoint_independent_of_Z'];lo,hi=self.at(row['lower']),self.at(row['upper'])
            if lo==hi:self.evidence['empty'].add(i);result=sy.Integer(0)
            else:
                variable=sy.Symbol(row['variable'],real=True);scope=dict(self.scope)
                scope[self.g.symbol(row['variable']).node]=variable
                body=SeamAlgebra(owner,self.source_bindings,self.evidence,scope).at(row['integrand'])
                result=sy.Integral(body,(variable,lo,hi))
        elif op=='original_lazy_flat_branch':
            raise AssertionError(('Active seam cannot use the rminus flat shortcut',i))
        else:return super().at(i)
        self.cache[i]=result;return result


def source_contract(owner):
    owner.assert_graph();g=owner.g
    rec=current.parent.linear.read(current.previous.RECEIPT)
    assert rec['all_passed'] and rec['actual_first_second_micro_contracted_mixed_join_rows']==246
    assert all(rec[key] for key in current.previous.GATES)
    assert owner.base.acceptance_loaded
    assert owner.bridge==owner.rminus.bridge
    assert owner.loop_input==owner.base.loop_input and owner.modified_input==owner.base.modified_input
    assert owner.coordinate.node==owner.bridge['macro']['phase']
    assert owner.jacobian.node==owner.bridge['length']
    assert owner.base.jacobian.node==owner.bridge['hb']!=owner.jacobian.node
    assert owner.loop['second_micro']['coordinate']==owner.bridge['micro']['phase']
    assert owner.loop['frozen_macro']['coordinate']==owner.coordinate.node
    for name,label in (('second_micro','micro'),('frozen_macro','macro')):
        assert owner.loop[name]['radius']==owner.bridge[label]['radius']
        for key in ('E','V'):
            assert owner.loop[name]['roots'][key]==owner.bridge[label]['fields'][key][:6]
    assert owner.modified['frozen_macro']['radial_background']==encode(owner.background_y)
    assert owner.factory.cache[owner.phase.node]==g.mul(owner.N,owner.jacobian)
    for factory in owner.factory.shared['factories'].values():
        assert factory.cache[owner.phase.node]==(g.mul(owner.N,owner.jacobian)
            if factory.variable==owner.coordinate else g.zero)
    assert len(owner.requests)==18
    for request in owner.requests:
        assert owner.factory.cache[request['input']].node==request['native_output']
        expected=g.mul(g.quotient(g.one,owner.jacobian,
            'same positive original frozen-macro physical Jacobian'),owner.ref(request['native_output']))
        assert request['physical_output']==expected.node
    assert len(owner.second_requests)==18
    assert owner.second_factory is not owner.base.factories['second_micro']
    assert owner.second_factory.cache[owner.base.phase.node]==g.mul(owner.N,owner.base.jacobian)
    for request in owner.second_requests:
        assert owner.second_factory.cache[request['input']].node==request['native_output']
        expected=g.mul(g.quotient(g.one,owner.base.jacobian,
            'same positive original micro physical Jacobian'),owner.ref(request['native_output']))
        assert request['physical_output']==expected.node
    assert encode(owner.second_background_y)==owner.modified['second_micro']['radial_background']
    for key,ids in owner.bridge['micro']['fields'].items():
        for node,dy in zip(ids,rows(owner.second_background_y[key])):
            assert owner.second_factory.cache[node]==g.mul(owner.base.jacobian,dy)
    for key in ('U','m','Q'):
        assert encode(owner.charts['frozen_macro']['physical'][key][:2])==[
            owner.modified['frozen_macro']['physical'][key],owner.modified['frozen_macro']['physical'][key+'_y']]
    intermediate_module=current.previous.previous.common_source.previous
    saved=current.parent.linear.read(intermediate_module.NAME)['intermediate_own_C3_chart_bindings']
    support=current.parent.linear.read(intermediate_module.RECEIPT)['checks']
    assert support['source']['actual_preRc_coordinate_radius_domain_and_C3_pressure_source_identity']
    assert support['actual_regional_recovery']['own_V_y_and_mass_y_zero_only_outside_velocity_supports']
    for name in current.CHARTS:
        chart=owner.charts[name]['intermediate']
        assert chart==saved['charts'][name] and name in saved['preRc_names']
        assert chart['region']=='beforeRc' and chart['complete_incoming_mass_not_reset']
        assert chart['fields']['P1']==chart['cumulative']['m13']
        for key in ('Uz1','Utheta1','F1'):assert set(chart['fields'][key])=={0}
        for key in ('m11','m12','m14'):
            assert chart['cumulative'][key]==owner.lower.rin.cut['endpoint'][key][:4]
        assert chart['fixed_pressure']==owner.modified[name]['pressure']
    assert saved['charts']['second_micro']['domain']==[g.one.node,g.constant(2).node]
    assert saved['charts']['frozen_macro']['domain']==[g.zero.node,g.one.node]
    for key,ids in owner.background['fields'].items():
        for node,dy in zip(ids,rows(owner.background_y[key])):
            assert owner.factory.cache[node]==g.mul(owner.jacobian,dy)
    for factory in owner.factory.shared['factories'].values():
        for node in factory.core_extensions:owner.base.guard_core(g.nodes[node])
    assert owner.loop['frozen_macro']['incoming']==owner.loop_input['endpoints']['second_micro']
    assert owner.modified['frozen_macro']['incoming_integrals']==owner.modified_input['endpoints']['second_micro']['integrals']
    return dict(passed=True,actual_new_macro_U_y2_y3_y4_C5_requests=18,
        actual_rebuilt_second_micro_endpoint_U_y2_y3_y4_C5_requests=18,
        second_micro_packet_rederived_from_original_full_chart_not_relabelled_from_near_one=True,
        native_coordinates_and_hb_length_Jacobians_distinct=True,
        actual_macro_phase_lift_is_N_times_length=True,
        actual_saved_y0_y1_and_all_background_radial_seeds_retained=True,
        actual_nonzero_micro_correction_and_pressure_endpoints_are_macro_incoming=True,
        inherited_total_native_calculus_and_periodic_lift_receipt=current.previous.RECEIPT)


def actual_background_recipes(owner):
    """Replay original formulas on exact source nodes, including all inlets."""
    g=owner.g;c=g.constant;a=bridge.BridgeAlgebra(g,7);b=bridge.BridgeAlgebra(g,6)
    saved=owner.bridge;micro=saved['micro'];macro=saved['macro']
    s,m=owner.ref(micro['phase']),owner.ref(macro['phase'])
    hb,length=owner.ref(saved['hb']),owner.ref(saved['length'])
    actual_end={key:b.at(owner.jet(micro['actual'][key]),s,c(2)) for key in ('phi','U')}
    moment_end={key:b.at(owner.jet(ids),s,c(2)) for key,ids in micro['moments'].items()}
    assert encode(actual_end)==saved['actual_micro_end']
    assert encode(moment_end)==saved['actual_moment_end']
    bar=saved['comparison_micro'];bar_end=saved['comparison_micro_end']
    assert bar_end['phi']==encode(a.at(owner.jet(bar['phi']),s,c(2)))
    assert bar_end['U']==encode(a.at(owner.jet(bar['U']),s,c(2)))
    assert bar_end['moments']=={key:encode(a.at(owner.jet(ids),s,c(2))) for key,ids in bar['moments'].items()}
    phi,U=owner.jet(bar_end['phi']),owner.jet(bar_end['U'])
    assert saved['comparison_macro']['phi']==encode(phi)
    assert saved['comparison_macro']['U']==encode(U)
    eq=dict(H=phi,M=U,K=a.mul(phi,U),A=a.mul(U,U),B=a.scale(a.mul(phi,phi),'1/2'),C=a.mul(phi,phi))
    assert bridge.RATES==dict(H=2,M=1,K=2,A=1,B=2,C=1)
    histories={}
    for key,rate in bridge.RATES.items():
        decay=g.unary('exp',g.mul(c(-rate),length,m))
        histories[key]=a.add(a.scale(owner.jet(bar_end['moments'][key]),decay),a.scale(eq[key],g.sub(g.one,decay)))
    assert encode(histories)==saved['comparison_macro']['moments']
    F0=owner.jet(saved['core7']['F0']);P0=owner.jet(saved['P0']);delta=owner.parameters['delta']
    known=bridge.comparison_directions(g,phi,U,histories,F0,P0,delta)
    assert encode(known)==saved['comparison_macro']['directions']
    pstar=owner.rminus.common['inlet_binding']['Pstar'];Pstar=owner.ref(pstar)
    actual=bridge.actual_fields(b,known,actual_end,m,length,hb,owner.radius,b.trim(F0),Pstar,phi,'macro')
    assert encode(actual)==macro['actual']
    ownJ=bridge.histories(b,actual['phi'],actual['U'],moment_end,m,length,'actual_macro')
    assert encode(ownJ)==macro['moments']
    fields,radial,generic=bridge.generic_fields(b,actual,ownJ,b.trim(F0),P0,Pstar,owner.radius,delta)
    assert encode(fields)==macro['fields'] and encode(radial)==macro['radial'] and encode(generic)==macro['generic']
    Ra=owner.ref(saved['Ra'])
    assert macro['radius']==g.mul(Ra,g.unary('exp',g.add(g.mul(c(2),hb),g.mul(length,m)))).node
    assert micro['radius']==g.mul(Ra,g.unary('exp',g.mul(hb,s))).node
    sigma=lambda x:bridge.inlet.previous.lower.own.jets.complete.Gamma.sigma(g,x)
    assert micro['chi']==g.add(g.one,g.mul(g.sub(hb,g.one),sigma(s))).node
    assert macro['chi']==hb.node
    # Consume the real comparison smoothing integral and its actual alpha,
    # rather than declaring an arbitrary comparison profile to be flat.
    log_integrals=[];axial_integrals=[]
    for i,row in enumerate(g.nodes[:len(owner.base.prefix)]):
        if row['operation']=='definite_integral' and row.get('exact_original_comparison_smoothing'):
            assert row['comparison_not_actual_field'] and row['upper']==s.node
            (log_integrals if i<=bar['logphi'][-1] else axial_integrals).append(i)
    assert len(log_integrals)==8 and len(axial_integrals)==8
    # The complete body replay below checks the actual argument s-1 and
    # exactly one hb factor in both comparison primitives.
    core_rho=owner.ref(g.nodes[saved['core7']['Phi'][0]]['rho'])
    t=g.symbol('original_prescribed_C6_comparison_t');rho=g.mul(c(4),g.unary('exp',g.mul(hb,t)))
    atcore=lambda ids:a.at(owner.jet(ids),core_rho,rho)
    alpha=g.sub(g.one,sigma(g.sub(t,g.one)))
    ratio=a.div(atcore(saved['radial7']['Phi'][1]),atcore(saved['core7']['Phi']),
        'same positive unique analytic Phi on rho4..4exp(2hb)')
    logpart=a.integral(a.scale(ratio,g.mul(alpha,hb,rho)),t,g.zero,s,
        exact_original_comparison_smoothing=True,comparison_not_actual_field=True)
    axialpart=a.integral(a.scale(atcore(saved['radial7']['Psi'][1]),
        g.mul(alpha,hb,rho,owner.parameters['epsilon'])),t,g.zero,s,
        exact_original_comparison_smoothing=True,comparison_not_actual_field=True)
    assert encode(logpart)==log_integrals and encode(axialpart)==axial_integrals
    geometry=[]
    H,L,r0=sy.symbols('micro_hb macro_length Ra',positive=True)
    ss,mm=sy.symbols('micro_s macro_m',real=True)
    for n in range(6):
        left=sy.diff(r0*sy.exp(H*ss),ss,n).subs(ss,2)/H**n
        right=sy.diff(r0*sy.exp(2*H+L*mm),mm,n).subs(mm,0)/L**n
        zero(left-right,('actual physical radius seam',n));geometry.append(n)
    return dict(passed=True,actual_source_macro_phi_U_and_all_six_nonzero_inlet_rows_replayed=True,
        actual_frozen_comparison_with_B_half_equilibrium_replayed=True,
        actual_comparison_directions_actual_fields_histories_generic_fields_replayed=True,
        actual_micro_alpha_argument_s_minus_one_and_hb_once_replayed=True,
        actual_chi_micro_sigma_s_and_macro_hb_replayed=True,
        actual_radius_physical_y0_to_y5_equal=geometry,
        original_comparison_primitive_rows=log_integrals+axial_integrals,
        source_background_jet_scope='physical y0..5 as ordinary Z functions; no numerical bounds')


def at_zero(expr,variable):
    expr=expr.subs(variable,0)
    return expr.replace(lambda e:isinstance(e,sy.Integral) and any(
        len(limit)==3 and limit[1]==limit[2] for limit in e.limits),lambda e:sy.Integer(0))


def independent_physical_transport_seam(owner):
    """Original source recipes, arbitrary nonzero inlets, separate native J.

    The fixture is scalar in Z. Its identities hold for arbitrary smooth Z
    functions. The pinned ordinary Z algebra then commutes with this finite
    physical-y induction; no high output row is bound in either graph.
    """
    g=current.source.FunctionTransportGraph();a=bridge.BridgeAlgebra(g,0)
    H,L,R0,P,F,p0,u0,bp,bu,D,A,C,S=sy.symbols(
        'hb length R0 Pstar F0 phi2 U2 barPhi2 barU2 D Hdrive Pdrive Sdrive',positive=True)
    inputs={str(v):g.symbol(str(v)) for v in (H,L,R0,P,F,p0,u0,bp,bu,D,A,C,S)}
    env={ref.node:sy.Symbol(name,positive=True) for name,ref in inputs.items()}
    fixed=lambda ref:bridge.FunctionJet((ref,))
    class Algebra:
        def __init__(self):self.cache=dict(env)
        def at(self,i):
            if i in self.cache:return self.cache[i]
            row=g.nodes[i];op=row['operation'];f=self.at
            if op=='bound_variable':out=sy.Symbol(row['name'],real=True)
            elif op=='exact_rational':out=sy.Rational(row['numerator'],row['denominator'])
            elif op=='sum':out=sy.Add(*(f(v) for v in row['arguments']))
            elif op=='product':out=sy.Mul(*(f(v) for v in row['arguments']))
            elif op=='negative':out=-f(row['argument'])
            elif op=='positive_quotient':out=f(row['numerator'])/f(row['denominator'])
            elif op=='analytic_unary':out=getattr(sy,row['name'])(f(row['argument']))
            elif op=='function_substitution':out=f(row['expression']).subs(f(row['variable']),f(row['value']))
            elif op=='definite_integral':out=sy.Integral(f(row['integrand']),
                (sy.Symbol(row['variable'],real=True),f(row['lower']),f(row['upper'])))
            else:raise AssertionError(op)
            self.cache[i]=out;return out
    y=sy.Symbol('independent_physical_y',real=True);R=R0*sy.exp(y)
    phi=p0*sy.exp(-H*D*R0*(sy.exp(y)-1)/2)
    drive=R*A+R*P**2*C+R**2*F**2*S
    Uy=-H*phi/bp*drive
    phi_rows=[sy.diff(phi,y,n).subs(y,0) for n in range(6)]
    U_rows=[u0]+[sy.diff(Uy,y,n-1).subs(y,0) for n in range(1,6)]
    evidence=[]
    for label,J in (('micro',H),('macro',L)):
        phase=g.symbol('fixture_'+label+'_native');native=sy.Symbol('fixture_'+label+'_native',real=True)
        jac=inputs[str(J)];radius=g.mul(inputs['R0'],g.unary('exp',g.mul(jac,phase)))
        known={key:fixed(inputs[str(body)]) for key,body in
            (('D_over_R',D),('drive_hydro',A),('drive_pressure',C),('drive_swirl',S))}
        got=bridge.actual_fields(a,known,dict(phi=fixed(inputs['phi2']),U=fixed(inputs['U2'])),
            phase,jac,inputs['hb'],radius,fixed(inputs['F0']),inputs['Pstar'],fixed(inputs['barPhi2']),'seam_fixture_'+label)
        incoming={};symbols={}
        for key in bridge.RATES:
            symbol=sy.Symbol('nonzero_incoming_'+key,real=True);ref=g.symbol(str(symbol))
            env[ref.node]=symbol;symbols[key]=symbol;incoming[key]=fixed(ref)
        histories=bridge.histories(a,got['phi'],got['U'],incoming,phase,jac,'seam_fixture_'+label)
        alg=Algebra();actual={key:alg.at(got[key].value.node) for key in ('phi','U')}
        for key,expected in (('phi',phi_rows),('U',U_rows)):
            expression=actual[key]
            for n in range(6):
                zero(at_zero(expression,native)/J**n-expected[n],('original actual field physical jet',label,key,n))
                evidence.append([label,key,n]);expression=sy.diff(expression,native)
        # Histories follow the source ODE. Derive the physical jet sequence
        # from the actual density jets and the genuine incoming value.
        densities=dict(H=[2*v for v in phi_rows],M=U_rows)
        from math import comb
        product=lambda left,right:[sum(comb(n,k)*left[k]*right[n-k] for k in range(n+1)) for n in range(5)]
        densities.update(K=[2*v for v in product(phi_rows,U_rows)],A=product(U_rows,U_rows),
            B=product(phi_rows,phi_rows),C=product(phi_rows,phi_rows))
        for key,rate in bridge.RATES.items():
            expression=alg.at(histories[key].value.node);target=symbols[key]
            for n in range(6):
                zero(at_zero(expression,native)/J**n-target,('original own history physical jet',label,key,n))
                evidence.append([label,key,n]);expression=sy.diff(expression,native)
                if n<5:target=densities[key][n]-rate*target
        # The frozen-comparison equilibrium agrees with the original
        # Duhamel source for arbitrary nonzero incoming values.
        comparison=bridge.histories(a,fixed(inputs['barPhi2']),fixed(inputs['barU2']),incoming,
            phase,jac,'comparison_seam_fixture_'+label)
        forcing=dict(H=2*bp,M=bu,K=2*bp*bu,A=bu**2,B=bp**2,C=bp**2)
        barphi,barU=fixed(inputs['barPhi2']),fixed(inputs['barU2'])
        equilibrium=dict(H=barphi,M=barU,K=a.mul(barphi,barU),A=a.mul(barU,barU),
            B=a.scale(a.mul(barphi,barphi),'1/2'),C=a.mul(barphi,barphi))
        for key,rate in bridge.RATES.items():
            expected=symbols[key]*sy.exp(-rate*J*native)+forcing[key]/rate*(1-sy.exp(-rate*J*native))
            expression=alg.at(comparison[key].value.node)
            for n in range(6):
                zero(at_zero(expression,native)/J**n-sy.diff(expected,native,n).subs(native,0)/J**n,
                    ('original comparison physical jet',label,key,n))
                evidence.append([label,'comparison_'+key,n]);expression=sy.diff(expression,native)
            decay=g.unary('exp',g.mul(g.constant(-rate),jac,phase))
            frozen=a.add(a.scale(incoming[key],decay),a.scale(equilibrium[key],g.sub(g.one,decay)))
            expression=alg.at(frozen.value.node)
            for n in range(6):
                zero(at_zero(expression,native)/J**n-sy.diff(expected,native,n).subs(native,0)/J**n,
                    ('actual frozen-comparison branch physical jet',label,key,n))
                evidence.append([label,'frozen_comparison_'+key,n]);expression=sy.diff(expression,native)
    return dict(passed=True,source_recipe_physical_y0_to_y5_rows=evidence,
        arbitrary_nonzero_phi_U_and_all_six_history_inlets=True,
        separate_symbolic_hb_length_and_common_physical_chi_hb=True,
        full_hydro_pressure_swirl_radius_Pstar_F0_factors_preserved=True,
        B_forcing_phi_squared_and_rate_two_gives_half_equilibrium=True,
        actual_frozen_comparison_specific_half_B_branch_checked_separately=True,
        actual_source_ODE_induction_not_assigned_output_derivatives=True)


def actual_smooth_source_seam(owner):
    """Finite physical-jet induction, followed by the actual loop recipes."""
    g=owner.g;s=owner.base.coordinate;H=owner.base.jacobian
    phase=owner.base.phase
    D=current.previous.PeriodicTotalNativeDerivative(g,s,phase,s,H,owner.N,owner.base.guard_core)
    D.shared['family']=owner.family
    evaluator=current.previous.modified_source.lexical.Algebra
    sigma=bridge.inlet.previous.lower.own.jets.complete.Gamma.sigma
    controls=[]
    for shift in (0,1):
        ref=sigma(g,g.sub(s,g.constant(shift)))
        for n in range(5):
            value=evaluator(g,{s.node:sy.Integer(2)}).at(ref.node)
            zero(value-(1 if n==0 else 0),('actual sigma at micro2',shift,n))
            controls.append([shift,n,ref.node])
            if n<4:ref=D.derivative(ref)
    # alpha=1-sigma(s-1) and every one of its first four derivatives vanish.
    # The comparison primitive derivatives through order five therefore
    # vanish by Leibniz, without assigning any derivative of Phi/Psi.
    from math import comb
    alpha=sy.symbols('alpha0:5');core=sy.symbols('actual_core_factor0:5')
    flat={term:sy.Integer(0) for term in alpha}
    for n in range(5):
        zero(sum(comb(n,k)*alpha[k]*core[n-k] for k in range(n+1)).subs(flat),
            ('actual comparison density derivative from alpha flatness',n))
    # Comparison histories share their inlet and the exact source ODE.
    # Directions are a common rational/Z composition of these histories.
    # Actual phi/U and histories then share one physical ODE jet recurrence.
    # Five background derivatives are required: roots use E_y and V_y.
    loop_source=current.previous.previous.loop_source
    Pstar=owner.ref(owner.rminus.common['inlet_binding']['Pstar'])
    roots=[];finite=[]
    for name,label in (('second_micro','micro'),('frozen_macro','macro')):
        bg=owner.bridge[label];chart=owner.loop[name]
        packet={key:{field:owner.jet(ids) for field,ids in bg[key].items()}
            for key in ('fields','radial','generic')}
        packet['radius']=owner.ref(bg['radius'])
        root=loop_source.generic_roots(g,packet,Pstar)
        assert encode(root)==chart['roots']
        # The pinned loop source already proves ordinary implicit/integral
        # calculus. Consume its actual inverse and source arguments here;
        # do not replay an older module's class aliases in a hydrated graph.
        inverse=g.nodes[chart['loop']['inverse'][0]]
        assert inverse['operation']=='original_monotone_phase_inverse'
        assert inverse['phase_function']==chart['loop']['Phi'][0]
        assert inverse['target_phase']==chart['actual_phase']
        assert inverse['angle_variable']==g.nodes[chart['loop']['angle']]['name']
        assert inverse['true_prescribed_source_C5_implicit_inverse']
        offset=g.sub(g.unary('log',g.quotient(owner.ref(chart['radius']),owner.ref(owner.bridge['Ra']),
            'same positive actual R/Ra')),owner.ref(owner.loop_input['phase_origin']))
        assert offset.node==chart['global_offset']
        assert g.unary('fractional_part',g.mul(owner.N,offset)).node==chart['actual_phase']
        loop={key:owner.jet(chart['loop'][key]) for key in ('A','B_over_Pstar')}
        a=current.outer.Algebra(g,5)
        invN=g.quotient(g.one,owner.N,'same exact positive integer N on every source chart')
        fast=a.scale(loop['A'],invN);ex=a.exp(fast)
        exm1=current.outer.FunctionJet((g.unary('expm1',fast.value),)+rows(ex)[1:])
        E,V=root['E'],root['V'];r=a.mul(E,exm1);v=a.scale(loop['B_over_Pstar'],invN)
        inc=dict(m=v,h=r,k=a.add(a.mul(V,r),a.mul(E,v),a.mul(r,v)),
            e=a.sub(a.add(a.scale(a.mul(V,v),2),a.mul(v,v)),a.add(a.mul(E,r),a.scale(a.mul(r,r),'1/2'))),
            p=a.add(a.mul(E,r),a.scale(a.mul(r,r),'1/2')))
        result=dict(E_N=a.mul(E,ex),V_N=a.add(V,v),delta_E=r,delta_V=v,A_over_N=fast,
            five_signed_density_increments=inc)
        assert encode(result)==chart['finite_N']
        roots.extend([name,key,ids] for key,ids in chart['roots'].items())
        finite.extend([name,key,ids] for key,ids in chart['finite_N'].items() if key!='five_signed_density_increments')
        for key in ('A','B_over_Pstar'):
            assert g.nodes[chart['loop'][key][0]]['operation'] in ('original_lazy_flat_branch','product','positive_quotient')
        assert chart['memory']['p']==g.one.node
    old=current.parent.linear.read(current.previous.RECEIPT)['checks']
    assert old['periodic_lift']['A_and_B_and_all_required_phase_Z_mixed_rows_periodic']
    assert old['independent_total_native_calculus']['angle_binder_holds_phase_while_main_native_factory_has_NJ']
    assert old['actual_leading_recovery_pressure_and_FTC']['actual_mean_ODE_all_five_rows_and_actual_saved_reduced_S_C3']
    return dict(passed=True,actual_sigma_s_and_s_minus_one_value_and_native1_to4_at2=controls,
        actual_comparison_phi_U_physical_y1_to_y5_flat_by_original_alpha_Leibniz=True,
        comparison_six_history_jets_from_common_inlet_and_source_ODE=True,
        actual_comparison_direction_jets_are_common_smooth_Z_compositions=True,
        actual_phi_U_six_histories_and_generic_background_physical_y0_to5_match_by_ODE_induction=True,
        actual_roots_replayed=roots,actual_finite_N_fields_replayed=finite,
        background_y5_supplies_root_y4_no_missing_derivative_padding=True,
        common_global_log_radius_phase_at_micro2_macro0_and_one_N=True,
        original_monotone_inverse_uniqueness_and_periodic_lift_transfer_all_y4_Z5_loop_rows=True,
        actual_micro_complete_correction_is_macro_incoming_so_mean_values_match=True,
        original_mean_ODE_and_one_minus_delta_recovery_give_common_y4_Z4_Q=True,
        theorem_scope='source identities at micro2/macro0; quantitative mixed norms remain open')


def actual_leading_recovery_and_common_N(owner):
    """Discharge leading-source hypotheses before own-output comparisons."""
    Z=sy.Symbol('Z',real=True);y=sy.Symbol('independent_leading_physical_y',real=True)
    U=sy.Function('actual_leading_U')(y,Z);M=sy.Function('actual_leading_mean')(y,Z)
    evidence=[];operators=[]
    previous=current.parent.linear.read(current.previous.RECEIPT)['checks']['actual_leading_recovery_pressure_and_FTC']
    assert previous['actual_N_y0_y1_y2_C3_and_signed_pressure_y1_y2_y3_C3']
    assert previous['actual_mean_ODE_all_five_rows_and_actual_saved_reduced_S_C3']
    for name in current.CHARTS:
        chart=owner.charts[name];env,v=active.prior.active_context(owner.rminus,False);delta=v['delta']
        for n,jet in enumerate(chart['physical']['U']):bind(env,jet,sy.diff(U,y,n),Z)
        bind(env,chart['physical']['m'][0],M,Z)
        alg=SeamAlgebra(owner,env);means=[M]
        for n in range(1,5):means.append(sy.diff(U,y,n-1)-means[-1])
        for n in range(5):
            zero(alg.at(chart['physical']['m'][n].value)-means[n],('actual leading mean recovery',name,n))
            Q=(2*Z*sy.diff(U,y,n)-(1-delta)*Z*means[n]-(1-Z**2)*sy.diff(means[n],Z))/(1-delta*Z**2)
            for k,ref in enumerate(rows(chart['physical']['Q'][n])):
                zero(alg.at(ref)-sy.diff(Q,Z,k),('actual leading Q recovery before seam binding',name,n,k))
                evidence.append([name,n,k,ref.node])
        f={key+suffix:chart['physical'][key][n] for key in ('U','m','Q')
            for n,suffix in enumerate(('', '_y','_yy','_yyy','_yyyy'))}
        R=owner.ref(owner.loop[name]['radius'])
        density=current.rc_source.radial.pressure_density(owner.g,f,owner.parameters['delta'],R)
        higher=current.high.higher_pressure_density(owner.g,f,owner.parameters['delta'],R,density)
        assert encode(density)==encode(chart['density']) and encode(higher)==encode(chart['higher_density'])
        assert encode([density['N1p'],higher['N1p_y'],higher['N1p_yy']])==encode(chart['N_y'])
        # Same fixed parametric operator as in the accepted arbitrary-U/Q
        # theorem. The actual saved reduced-S recipe is separately replayed.
        g=owner.g;a=current.outer.Algebra(g,3)
        z=current.outer.FunctionJet((g.symbol('Z'),g.one,g.zero,g.zero));one=a.fixed(1);de=a.fixed(owner.parameters['delta'])
        d=a.sub(one,a.mul(z,z));L=a.sub(one,a.mul(de,a.mul(z,z)))
        u,q,qy=(a.trim(owner.jet(owner.modified[name]['physical'][key])) for key in ('U','Q','Q_y'))
        qz=a.shift(owner.jet(owner.modified[name]['physical']['Q']))
        div=lambda jet:a.div(jet,L,'same original positive meridional L')
        linear=div(a.add(q,a.scale(a.mul(a.sub(one,de),a.mul(z,qz)),'1/2'),qy))
        quadratic=a.mul(q,a.add(a.scale(q,'1/2'),qy))
        axial=div(a.mul(u,a.add(a.scale(a.mul(z,q),-2),a.mul(d,qz),a.scale(a.mul(z,qy),-2))))
        reduced=a.add(linear,quadratic,axial)
        assert encode(reduced)==owner.modified[name]['S']
        operators.append(dict(chart=name,radius=R.node,N_y=encode(chart['N_y']),saved_S=owner.modified[name]['S']))
    return dict(passed=True,actual_leading_Q_recovery_y0_to_y4_C4_rows=evidence,
        actual_mean_ODE_not_an_assigned_mean_derivative=True,
        actual_N_operators_replayed_on_each_actual_leading_packet=operators,
        arbitrary_U_Q_operator_identity_receipt=current.previous.RECEIPT,
        common_U_and_mean_physical_jets_imply_common_Q_then_common_S_and_N_y0_to2=True,
        common_N_is_derived_before_own_output_input_bindings=True,
        physical_R_y_equals_R_and_both_native_Jacobians_retained=True)


def actual_macro_pressure_and_FTC(owner):
    name='frozen_macro';chart=owner.charts[name];old=owner.modified[name]
    Z=sy.Symbol('Z',real=True);y=sy.Symbol('independent_macro_physical_y',real=True)
    u=sy.Function('independent_actual_macro_U')(y,Z);q=sy.Function('independent_actual_macro_Q')(y,Z)
    env,v=active.prior.active_context(owner.rminus,False);delta=v['delta']
    R=sy.Symbol('actual_macro_R0',positive=True)*sy.exp(y);env[owner.radius.node]=R
    for key,body in (('U',u),('Q',q)):
        for n,jet in enumerate(chart['physical'][key]):bind(env,jet,sy.diff(body,y,n),Z)
    alg=SeamAlgebra(owner,env);d=1-Z**2;L=1-delta*Z**2
    S=(q+(1-delta)*Z*sy.diff(q,Z)/2+sy.diff(q,y))/L+q*(q/2+sy.diff(q,y))+u*(-2*Z*q+d*sy.diff(q,Z)-2*Z*sy.diff(q,y))/L
    N=-S/2+(sy.diff(q,y)+sy.diff(q,y,2))/R;evidence=[]
    for k,node in enumerate(old['S']):zero(alg.at(node)-sy.diff(S,Z,k),('actual macro saved reduced S',k))
    for n,jet in enumerate(chart['N_y']):
        for k,ref in enumerate(rows(jet)):
            zero(alg.at(ref)-sy.diff(N,y,n,Z,k),('actual macro pressure source',n,k))
            evidence.append(['N',n,k,ref.node])
    for n in range(3):
        for key,target in (('pressure_y',R*N),('m15_y',-R**2*N)):
            for k,ref in enumerate(rows(chart[key][n])):
                zero(alg.at(ref)-sy.diff(target,y,n,Z,k),('actual macro pressure derivative',key,n,k))
                zero(alg.at(rows(chart['raw_'+key][n])[k])-alg.at(ref),('actual raw/coupled macro pressure',key,n,k))
                evidence.append([key,n,k,ref.node])
    env,v=active.prior.active_context(owner.rminus,False)
    x=sy.Symbol(owner.g.nodes[owner.coordinate.node]['name'],real=True);J=sy.Symbol('actual_macro_length',positive=True)
    R=sy.Symbol('actual_macro_R0',positive=True)*sy.exp(J*x)
    env[owner.coordinate.node]=x;env[owner.jacobian.node]=J;env[owner.radius.node]=R
    q=sy.Function('actual_macro_saved_Q')(x,Z);S=sy.Function('actual_macro_saved_S')(x,Z)
    bind(env,owner.jet(old['physical']['Q']),q,Z)
    bind(env,owner.jet(old['physical']['Q_y']),sy.diff(q,x)/J,Z)
    bind(env,owner.jet(old['S']),S,Z)
    for key,ids in owner.modified_input['initial'].items():bind(env,owner.jet(ids),sy.Function('actual_global_initial_'+key)(Z),Z)
    for key,ids in old['incoming_integrals'].items():bind(env,owner.jet(ids),sy.Function('actual_nonzero_macro_incoming_'+key)(Z),Z)
    alg=SeamAlgebra(owner,env);Py=-R*S/2+sy.diff(q,x)/J+sy.diff(q,x,2)/J**2;ftc=[]
    for key,target in (('m13',Py),('m15',-R*Py)):
        for k,node in enumerate(old['pressure'][key]):
            zero(sy.diff(alg.at(node),x).doit()/J-sy.diff(target,Z,k),('actual macro original signed FTC',key,k))
            ftc.append([key,k,node])
    for key,power in (('S0',1),('S1',2)):
        for k,node in enumerate(old['local_integrals'][key]):
            row=owner.g.nodes[node]
            assert row['operation']=='definite_integral' and row['lower']==owner.g.zero.node
            assert row['upper']==owner.coordinate.node and row['physical_radius_power']==power
            assert row['physical_Jacobian_applied_once'] and row['ordinary_slow_Z_derivative_order']==k
    return dict(passed=True,actual_macro_N_and_signed_pressure_derivative_C3_rows=evidence,
        actual_original_macro_pressure_primitive_C3_FTC_rows=ftc,
        native_macro_differentiation_divided_by_length_not_hb=True,
        nonzero_incoming_S0_S1_and_global_finite_boundary_retained=True,
        own_pressure_outputs_not_bound=True)


def actual_own_mixed_seam(owner):
    env,Z,slow=active.prior.rin_check.environment(owner.lower.rin,history=True,pressure_theorem=True)
    s=sy.Symbol(owner.g.nodes[owner.base.coordinate.node]['name'],real=True)
    x=sy.Symbol(owner.g.nodes[owner.coordinate.node]['name'],real=True)
    H=env[owner.parameters['hb'].node];J=sy.Symbol('actual_macro_length',positive=True)
    env[owner.jacobian.node]=J
    env[owner.rminus.active_sc.node]=SeamAlgebra(owner,env).at(owner.parameters['sc'])
    R2=sy.Symbol('actual_common_micro2_macro0_radius',positive=True)
    env[owner.base.radius.node]=R2*sy.exp(H*(s-2));env[owner.radius.node]=R2*sy.exp(J*x)
    # Bind defining primitives, including every earlier nonempty chart.
    # Macro local integrals are left in the graph and become empty at m=0.
    for name in ('first_micro','second_micro'):
        for key,ids in owner.base.loop_input['charts'][name]['contributions'].items():
            bind(env,owner.jet(ids),sy.Function('actual_'+name+'_finite_N_primitive_'+key)(s,Z),Z)
        for key,ids in owner.base.modified_input['charts'][name]['local_integrals'].items():
            bind(env,owner.jet(ids),sy.Function('actual_'+name+'_signed_primitive_'+key)(s,Z),Z)
    background_mean=sy.Function('actual_common_background_mean')
    mapped=2+J*x/H
    bind(env,owner.jet(owner.bridge['micro']['fields']['m']),background_mean(s,Z),Z)
    bind(env,owner.jet(owner.bridge['macro']['fields']['m']),background_mean(mapped,Z),Z)
    for key in ('Q','Q_y'):
        bind(env,owner.jet(owner.base.modified_input['charts']['first_micro']['physical'][key]),
            sy.Function('actual_retained_first_micro_leading_'+key)(s,Z),Z)
    # These are proved leading-source inputs to the own pressure operator.
    # They are not the desired own pressure/velocity output rows.
    for name,argument in (('second_micro',s),('frozen_macro',mapped)):
        for key in ('Q','Q_y'):
            bind(env,owner.jet(owner.modified[name]['physical'][key]),sy.Function('common_actual_leading_'+key)(argument,Z),Z)
        for n,jet in enumerate(owner.charts[name]['N_y']):
            bind(env,jet,sy.Function('common_actual_leading_N_y'+str(n))(argument,Z),Z)
    alg=SeamAlgebra(owner,env);flow=[];leading_mean=[]
    for name,value in (('second_micro',owner.g.constant(2)),('frozen_macro',owner.g.zero)):
        coordinate=owner.ref(owner.loop[name]['coordinate'])
        jet=current.outer.Algebra(owner.g,5).at(owner.charts[name]['physical']['m'][0],coordinate,value)
        leading_mean.append([alg.at(ref) for ref in rows(jet)])
    for k,(left,right) in enumerate(zip(*leading_mean)):
        zero(left-right,('actual finite_N leading mean incoming seam',k))
    yy=sy.Symbol('independent_own_log_radius',real=True);delta=alg.at(owner.parameters['delta'])
    for name in current.CHARTS:
        f=owner.charts[name]['own']['fields'];mean=alg.at(f['mean1'][0].value)
        radius=alg.at(owner.ref(owner.loop[name]['radius']))
        model_mean=mean*sy.exp(-yy);model_R=radius*sy.exp(yy)
        model_Q=(-(1+delta)*Z*model_mean-(1-Z**2)*sy.diff(model_mean,Z))/(1-delta*Z**2)
        models=dict(Uz1=sy.Integer(0),Utheta1=sy.Integer(0),F1=sy.Integer(0),mean1=model_mean,
            K1=model_mean,M1=model_R*model_mean,Q1=model_Q,V1=model_R*model_Q)
        for key,body in models.items():
            for n,jet in enumerate(f[key]):
                zero(alg.at(jet.value)-sy.diff(body,yy,n).subs(yy,0),('actual supported own flow',name,key,n))
                flow.append([name,key,n])
        for key in ('M1','mean1','Q1','V1'):
            for new,old in zip(rows(f[key][1]),owner.charts[name]['intermediate']['fields'][key+'_y']):
                zero(alg.at(new)-alg.at(old),('actual original own first derivative',name,key))
    evidence=[]
    for spec in owner.contract:
        for n in range(spec['physical_y_order']+1):
            for k,(left,right) in enumerate(zip(rows(owner.join['left'][spec['group']][spec['field']][n]),
                rows(owner.join['right'][spec['group']][spec['field']][n]))):
                zero(alg.at(left)-alg.at(right),('actual micro2/macro0 own mixed seam',spec['group'],spec['field'],n,k))
                evidence.append(dict(group=spec['group'],field=spec['field'],physical_y_order=n,
                    ordinary_Z_order=k,left_node=left.node,right_node=right.node))
    assert len(evidence)==246
    return dict(passed=True,actual_full_own_micro_macro_mixed_join_rows=evidence,actual_supported_own_flow_rows=flow,
        actual_leading_mean_incoming_value_C5_rows=6,
        actual_micro_nonzero_mean_pressure_and_all_five_histories_transfer_to_macro=True,
        actual_macro_empty_integrals_do_not_zero_predecessor_memory=True,
        actual_own_one_plus_delta_distinguished_from_leading_one_minus_delta=True)


CONSTRUCTION_CHECKS=(source_contract,actual_background_recipes,independent_physical_transport_seam)
CHECKS=CONSTRUCTION_CHECKS+(actual_smooth_source_seam,actual_leading_recovery_and_common_N,
    actual_macro_pressure_and_FTC,actual_own_mixed_seam)


def run_construction(owner):
    began=time.monotonic();checks={}
    for fn in CONSTRUCTION_CHECKS:
        checks[fn.__name__]=fn(owner)
        print('MICRO_MACRO_'+fn.__name__.upper()+'_PASSED',flush=True)
    owner.assert_graph();raw=current.parent.linear.read(current.NAME)
    assert raw['micro_macro_mixed_source_join_bindings']==owner.snapshot
    assert raw['exact_graph_nodes']==owner.g.nodes[:len(raw['exact_graph_nodes'])]
    assert not any(raw[key] for key in current.GATES+current.OPEN)
    rec=dict(construction_checks_passed=True,all_passed=False,source_family=owner.family,checks=checks,
        **dict.fromkeys(current.GATES+current.OPEN,False),
        pending=['actual comparison flatness through physical y5 and derived loop seam',
            'actual macro pressure operator and signed FTC','actual own contracted 246 mixed seam rows'],
        exact_graph_nodes=len(raw['exact_graph_nodes']),
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        numerical_oracle_available=False,global_n1=False,execution_seconds=time.monotonic()-began)
    (current.HERE/CONSTRUCTION_RECEIPT).write_text(json.dumps(rec,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('ACTUAL_MICRO_MACRO_BACKGROUND_CONSTRUCTION_CHECKED_NOT_ACCEPTED',flush=True)
    return rec


def run(owner):
    began=time.monotonic();checks={}
    for fn in CHECKS:
        checks[fn.__name__]=fn(owner)
        print('MICRO_MACRO_'+fn.__name__.upper()+'_PASSED',flush=True)
    owner.assert_graph();raw=current.parent.linear.read(current.NAME)
    assert raw['micro_macro_mixed_source_join_bindings']==owner.snapshot
    assert raw['exact_graph_nodes']==owner.g.nodes[:len(raw['exact_graph_nodes'])]
    rec=dict(all_passed=True,source_family=owner.family,checks=checks,
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        **dict.fromkeys(current.GATES,True),**dict.fromkeys(current.OPEN,False),
        actual_second_micro_frozen_macro_contracted_mixed_join_rows=246,
        exact_graph_nodes=len(raw['exact_graph_nodes']),numerical_oracle_available=False,
        global_n1=False,execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(rec,indent=2)+'\n',encoding='utf-8',newline='\n')
    owner.load_acceptance();print('ACTUAL_MICRO_MACRO_SOURCE_AND_OWN246_ACCEPTED',flush=True)
    return rec


if __name__=='__main__':
    raise SystemExit('run_construction(owner) reuses the live accepted micro ancestor')
