"""Scoped source, variable log-measure, native bump and pressure checks."""
import json
from pathlib import Path
import time

import sympy as sy

import lei_ren_part1_paper_compliant_current_original_true_Rm_patch_finite_N_C5_functions as current
import lei_ren_part1_paper_compliant_current_original_true_Rm_patch_C6_functions_check as leading_check
import lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions_check as radial_check

Algebra=leading_check.Algebra
lexical=radial_check.prior


def jet(g,variables,body,order=5):
    Z=sy.Symbol('Z')
    return current.FunctionJet(tuple(current.source.expression(g,sy.diff(body,Z,k),variables) for k in range(order+1)))


def actual_source(owner):
    g=owner.g;r=current.rows;chart=owner.chart;bg=owner.leading.chart
    free={i:set(v) for i,v in owner.leading.free_source_seeds.items()};name=g.nodes[owner.coordinate.node]['name'];count=0
    assert owner.leading.acceptance_loaded and owner.source_chart is bg and chart['source_chart'] is bg
    assert all(len(r(q))==7 for q in bg['fields'].values()) and all(len(r(q))==6 for q in bg['generic'].values())
    assert owner.incoming is owner.loop.Rm_correction and owner.incoming_integrals is owner.modified.Rm['integrals']
    assert owner.initial is owner.modified.initial and owner.phase_origin is owner.loop.phase_origin
    assert owner.bindings['original_all_N_source']==current.original_transport.source_binding()
    assert owner.bindings['original_log_kernel']==current.parent.active.current.ast_binding(current.original_patch.patch_weights)
    for attr in ('N','eta','dstar','sc'):assert getattr(owner,attr) is getattr(owner.loop,attr)
    offset=g.sub(g.unary('log',g.quotient(bg['radius'],owner.leading.background.bridge.Ra,
        'same original positive actual R/Ra')),owner.phase_origin)
    assert chart['actual_phase']==g.unary('fractional_part',g.mul(owner.N,offset))
    assert lexical.free_variables(g,owner.phase,free)=={name}
    assert chart['loop']['original_collar_upper'] is None and chart['roots']['p2'] is bg['generic']['p2']
    assert all(any(q!=g.zero for q in r(v)) for v in owner.incoming.values())
    assert chart['memory']['p']==g.one
    for key,q in bg['fields'].items():
        for ref,dr in zip(r(q),r(owner.radial_background[key])):
            assert owner.derivative.cache[ref.node]==g.mul(bg['physical_measure'],dr)
    for gamma in owner.leading.gamma:
        assert owner.derivative.cache[gamma[0].node]==gamma[1]
        assert owner.derivative.cache[gamma[1].node]==gamma[2]
    assert owner.derivative.cache[owner.phase.node]==g.zero
    assert owner.derivative.shared['family']==owner.family
    for group,order in ((chart['roots'].values(),5),(chart['history'].values(),5),
        (chart['modified_fields'].values(),5),(chart['finite_N']['five_signed_density_increments'].values(),5),
        (owner.generic.values(),4),(owner.pressure.values(),3),(owner.cumulative_integrals.values(),3)):
        for q in group:
            assert len(r(q))==order+1
            for ref in r(q):assert lexical.free_variables(g,ref,free)<={'Z',name},(ref.node,lexical.free_variables(g,ref,free))
            count+=order+1
    for key in ('A','B_over_Pstar','A_phi','B_over_Pstar_phi','inverse','q'):
        assert len(r(chart['loop'][key]))==6
        for ref in r(chart['loop'][key]):assert lexical.free_variables(g,ref,free)<={'Z',name}
    for key in ('E_y','V_y','m_y','slow_A_y','slow_B_y'):
        assert len(r(owner.radial[key]))==6
        for ref in r(owner.radial[key]):assert lexical.free_variables(g,ref,free)<={'Z',name}
    for key,q in owner.physical.items():assert len(r(q))==(5 if key in ('Q','Q_y') else 6)
    for key,q in chart['contributions'].items():
        for n,ref in enumerate(r(q)):
            row=g.nodes[ref.node]
            assert row['operation']=='definite_integral' and row['ordinary_slow_Z_derivative_order']==n
            assert row['lower']==g.one.node and row['upper']==owner.coordinate.node
            assert row['physical_measure']=='dt/t' and row['original_rate']==str(current.RATES[key])
            assert row['nonzero_original_Rm_correction_inlet_retained'] and row['predecessor_memory_never_reset']
    for key,q in owner.local_integrals.items():
        for n,ref in enumerate(r(q)):
            row=g.nodes[ref.node]
            assert row['operation']=='definite_integral' and row['ordinary_slow_Z_derivative_order']==n
            assert row['lower']==g.one.node and row['upper']==owner.coordinate.node
            assert row['nonzero_preceding_S0_S1_retained'] and row['physical_Jacobian_applied_once']
            assert row['physical_radius_power']==(1 if key=='S0' else 2)
    assert r(chart['modified_fields']['P0'])==r(owner.P0)[:6]
    assert owner.Rh is owner.endpoints['Rh'] and owner.Rh_correction is owner.Rh['correction']
    return dict(passed=True,actual_modified_density_history_generic_pressure_rows=count,
        same_true_background_p2_original_N_phase_and_hbB_origin=True,
        all_nonzero_correction_and_pressure_predecessors_retained=True,
        exact_log_measure_and_rate_zero_pressure_memory=True,
        native_background_and_gamma_first_second_derivative_seeds=True,
        total_fast_phase_radial_C5_and_complete_signed_pressure_C3=True,
        independent_P0_original_rminus_boundary_and_Rh_endpoints=True,
        no_free_angle_old_phase_or_integration_variable=True,
        six_outer_charts_full_preRp_numeric_norms_and_global_claims_not_installed=True)


def log_transport_calculus():
    g=current.source.FunctionTransportGraph();x=g.symbol('fixture_patch_x');Z=sy.Symbol('Z');X=sy.Symbol('fixture_patch_x',positive=True)
    variables={'fixture_patch_x':x,'Z':g.symbol('Z')}
    density={key:(1+Z+Z**5/(11+i))*(1+X**2/(13+i)) for i,key in enumerate(current.RATES)}
    inlet={key:2+Z**2/(17+i)+Z**5/(19+i) for i,key in enumerate(current.RATES)}
    make=lambda body:jet(g,variables,body)
    got,parts,memory=current.transport_log(g,{k:make(v) for k,v in density.items()},
        {k:make(v) for k,v in inlet.items()},x,g.one,'fixture')
    alg=Algebra(g,{x.node:X,g.symbol('Z').node:Z});count=0
    for key,rate in current.RATES.items():
        for n,ref in enumerate(current.rows(got[key])):
            value=alg.at(ref.node)
            assert sy.simplify(X*sy.diff(value,X)+sy.Rational(rate)*value-sy.diff(density[key],Z,n))==0,(key,n,'log ODE')
            assert sy.simplify(value.subs(X,1)-sy.diff(inlet[key],Z,n))==0,(key,n,'nonzero inlet')
            count+=1
    zero={key:make(sy.Integer(0)) for key in current.RATES}
    quiet,_,_=current.transport_log(g,zero,{k:make(v) for k,v in inlet.items()},x,g.one,'quiet_fixture')
    for key,rate in current.RATES.items():
        for n,ref in enumerate(current.rows(quiet[key])):
            assert sy.simplify(alg.at(ref.node)-sy.diff(inlet[key],Z,n)*X**(-sy.Rational(rate)))==0,(key,n,'quiet memory')
            count+=1
    assert memory['p']==g.one
    return dict(passed=True,independent_nonzero_inlet_log_ODE_and_quiet_memory_C5_rows=count,
        all_five_rates_original_dt_over_t_and_pressure_rate_zero_checked=True,
        nonconstant_native_density_exercises_variable_Jacobian=True,mathematical_fixture_only=True)


class InletAlgebra(Algebra):
    def at(self,i):
        if i in self.cache:return self.cache[i]
        row=self.g.nodes[i]
        if row['operation']=='definite_integral' and self.at(row['lower'])==self.at(row['upper']):
            self.cache[i]=sy.Integer(0);return self.cache[i]
        return super().at(i)


def actual_inlet_seam(owner):
    """Compare actual chart root inputs and retained histories at x1.

    Identical root functions/global phase define the same unique original
    loop inverse. Dummy integration-angle names do not change that function.
    """
    g=owner.g;r=current.rows;lead=owner.leading;old=owner.loop.charts['postrestore'];oldbg=old['source_chart']
    Z=sy.Symbol('Z');point=sy.Rational(3,8);P=sy.Rational(3,2);R=sy.Rational(7,3)
    am=sy.exp(-sy.Rational(3,5))/(1+Z*Z)
    bodies={key:2+j+Z/(17+j)+Z**6/(23+j) for j,key in enumerate(('m','h','k','e','p','P0'))}
    bodies.update(E=am,V=4*Z/P,pressure=bodies['P0']+bodies['p'])
    bind={owner.coordinate.node:sy.Integer(1),old['coordinate'].node:sy.Integer(1),
        owner.delta.node:sy.Rational(1,5),lead.background.core.parameters['logP'].node:sy.log(P),
        lead.inlet['radius'].node:R,oldbg['radius'].node:R,lead.bump.normalization.node:sy.Rational(7,5)}
    fill=lambda q,body:bind.update({ref.node:sy.diff(body,Z,n).subs(Z,point) for n,ref in enumerate(r(q))})
    for key,body in bodies.items():
        fill(lead.inlet['fields'][key],body);fill(oldbg['fields'][key],body)
    fill(lead.P0,bodies['P0'])
    for q in lead.control_functions:fill(q,1+Z**6/31)
    fill(oldbg['radial']['E_y'],am/10);fill(oldbg['radial']['V_y'],sy.Integer(0))
    for j,(key,q) in enumerate(owner.incoming.items()):fill(q,5+j+Z**5/(37+j))
    for j,(key,q) in enumerate(owner.incoming_integrals.items()):fill(q,11+j+Z**3/(41+j))
    bind[g.symbol('Z').node]=point;alg=InletAlgebra(g,bind);count=0
    for key,q in lead.chart['fields'].items():
        for n,ref in enumerate(r(q)):
            assert sy.simplify(alg.at(ref.node)-sy.diff(bodies[key],Z,n).subs(Z,point))==0,('leading seam',key,n);count+=1
    for gamma in lead.gamma:
        assert all(alg.at(ref.node)==0 for ref in gamma)
    for key,q in owner.chart['roots'].items():
        for n,ref in enumerate(r(q)):
            assert sy.simplify(alg.at(ref.node)-alg.at(r(old['roots'][key])[n].node))==0,('actual root seam',key,n);count+=1
    for key,q in owner.endpoints['Rm']['correction'].items():
        for n,ref in enumerate(r(q)):
            assert alg.at(ref.node)==alg.at(r(owner.incoming[key])[n].node),('correction seam',key,n);count+=1
    for key,q in owner.endpoints['Rm']['integrals'].items():
        for n,ref in enumerate(r(q)):
            assert alg.at(ref.node)==alg.at(r(owner.incoming_integrals[key])[n].node),('pressure integral seam',key,n);count+=1
    assert owner.bindings['original_loop']==owner.loop.bindings['original_loop']
    return dict(passed=True,actual_leading_root_correction_and_pressure_integral_seam_rows=count,
        true_C6_leading_inputs_and_C5_generic_roots_match_at_x1=True,
        all_three_original_bumps_and_first_second_derivatives_flat_at_x1=True,
        same_root_phase_and_unique_loop_inverse_identify_modified_velocity_seam=True,
        actual_nonzero_C5_corrections_and_C3_integrals_match_by_exact_inlet_identity=True,
        mathematical_source_fixture_not_new_original_scale_norm=True)


def gamma_native_calculus():
    g=current.source.FunctionTransportGraph();x=g.symbol('fixture_gamma_x');phase=g.symbol('fixture_phase')
    W=current.previous.ExactBumpWeights(g);gamma=W.gamma(x,0)
    D=current.radial.NativeDerivative(g,x,phase)
    D.shared['free'].update(W.free_seeds);D.seed(gamma[0],gamma[1]);D.seed(gamma[1],gamma[2])
    body=g.add(g.mul(x,gamma[0]),g.mul(x,x,gamma[1]));derivative=D.derivative(body)
    X=sy.Symbol('X');N=sy.Rational(7,5);beta=sy.exp(-1/(1-(40*(X-sy.Rational(5,4)))**2));G=40*beta/N
    expected=sy.diff(X*G+X*X*sy.diff(G,X),X)
    point=sy.Rational(5,4)+sy.Rational(1,200)
    alg=Algebra(g,{x.node:point,W.normalization.node:N})
    assert sy.simplify(alg.at(derivative.node)-expected.subs(X,point))==0
    partial=W.weight(0,sy.Rational(3,5),upper=x);D.shared['free'].update(W.free_seeds)
    endpoint_derivative=D.derivative(partial);symbolic=Algebra(g,{x.node:X,W.normalization.node:N})
    want=X**sy.Rational(3,5)*40*sy.Function('raw_beta')(40*(X-sy.Rational(5,4)))/N
    assert sy.simplify(symbolic.at(endpoint_derivative.node)-want)==0
    return dict(passed=True,actual_scaled_beta_native_chain_first_second_factors_checked=True,
        exact_moving_partial_weight_upper_endpoint_derivative_checked=True,
        mathematical_fixture_only=True)


def log_pressure_calculus():
    g=current.source.FunctionTransportGraph();x=g.symbol('fixture_pressure_patch_x');X=sy.Symbol('fixture_pressure_patch_x',positive=True);Z=sy.Symbol('Z')
    variables={'fixture_pressure_patch_x':x,'Z':g.symbol('Z')};P=sy.Rational(3,2);Rm=sy.Rational(7,3);de=sy.Rational(1,5)
    m=2+Z+Z**5/11+(3+Z**2+Z**5/13)*X+(1+Z**3+Z**5/17)*X**2
    U=m+X*sy.diff(m,X);Uy=X*sy.diff(U,X)
    physical=current.radial.physical_pressure_fields(g,{'V':jet(g,variables,U/P),'m':jet(g,variables,m/P)},
        {'V_y':jet(g,variables,Uy/P)},g.constant(de),g.constant(P))
    S=current.radial.reduced_pressure_S(g,physical,g.constant(de))
    incoming={key:jet(g,variables,body,3) for key,body in {'S0':7+Z**3/19,'S1':11+Z**3/23}.items()}
    initial={key:jet(g,variables,body,3) for key,body in {'Q_plus_Q_y':13+Z**3/29,'R_Q_y':17+Z**3/31}.items()}
    boundary,local,cumulative,pressure=current.pressure_continuation(g,S,physical,g.constant(Rm),x,incoming,initial)
    L=1-de*Z**2;d=1-Z**2;q=(2*Z*U-(1-de)*Z*m-d*sy.diff(m,Z))/L
    qy=X*sy.diff(q,X);qyy=X*sy.diff(qy,X)
    ss=(q+(1-de)*Z*sy.diff(q,Z)/2+qy)/L+q*(q/2+qy)+U*(-2*Z*q+d*sy.diff(q,Z)-2*Z*qy)/L
    N1p=-ss/2+(qy+qyy)/(Rm*X);point=sy.Rational(3,8);alg=Algebra(g,{x.node:X,g.symbol('Z').node:point});count=0
    for key,multiplier in (('m13',Rm),('m15',-Rm**2*X)):
        for n,ref in enumerate(current.rows(pressure[key])):
            value=alg.at(ref.node)
            assert sy.simplify(sy.diff(value,X)-sy.diff(multiplier*N1p,Z,n).subs(Z,point))==0,(key,n,'direct source')
            ii=alg.at(current.rows(incoming['S0' if key=='m13' else 'S1'])[n].node)
            bb=alg.at(current.rows(boundary['Q_plus_Q_y' if key=='m13' else 'R_Q_y'])[n].node).subs(X,1)
            oo=alg.at(current.rows(initial['Q_plus_Q_y' if key=='m13' else 'R_Q_y'])[n].node)
            want=-ii/2+bb-oo if key=='m13' else ii/2-bb+oo
            assert sy.simplify(value.subs(X,1)-want)==0,(key,n,'inherited pressure');count+=1
    return dict(passed=True,independent_direct_N1p_and_inherited_pressure_C3_rows=count,
        actual_Rm_and_Rm_squared_x_measures_and_m15_endpoint_sign=True,
        nonzero_S0_S1_and_original_global_boundary_retained=True,mathematical_fixture_only=True)


def run(owner):
    began=time.monotonic();owner.assert_graph();raw=current.parent.linear.read(current.NAME)
    assert raw['Rm_patch_finite_N_C5_bindings']==owner.snapshot
    assert raw['exact_graph_nodes']==owner.g.nodes[:len(raw['exact_graph_nodes'])]
    checks={'actual_source':actual_source(owner),'actual_inlet_seam':actual_inlet_seam(owner)};print('TRUE_RM_PATCH_FINITE_N_SOURCE_AND_SEAM_CHECKED',flush=True)
    for name,fn in (('log_transport',log_transport_calculus),('gamma_native',gamma_native_calculus),('log_pressure',log_pressure_calculus)):
        checks[name]=fn();print('TRUE_RM_PATCH_'+name.upper()+'_CHECKED',flush=True)
    owner.assert_graph()
    for name,digest in owner.hashes.items():assert current.sha(name)==digest,name
    receipt=dict(all_passed=True,source_family=owner.family,checks=checks,
        **dict.fromkeys(current.GATES,True),**dict.fromkeys(current.OPEN,False),
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8',newline='\n')
    print('TRUE_RM_PATCH_FINITE_N_C5_AND_PRESSURE_C3_CHECKED',flush=True);return receipt


if __name__=='__main__':raise SystemExit('run(owner) reuses actual checked leading patch')
