"""Scoped compact source, normalization, implicit calculus and seam proofs."""
import json
from pathlib import Path
import time
import sympy as sy

import lei_ren_part1_paper_compliant_current_original_true_compact_C5_functions as current
import lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions_check as prior

Algebra,lexical=prior.Algebra,prior.lexical


def jet(g,body,order=5):return prior.jet(g,{'Z':g.symbol('Z')},sy.sympify(body),order)


def graph(nodes):
    g=current.source.FunctionTransportGraph();g.nodes=[dict(q) for q in nodes]
    g.keys={json.dumps(q,sort_keys=True,separators=(',',':')):i for i,q in enumerate(g.nodes)}
    return g


def actual_source(owner):
    g=owner.g;r=current.rows;free={i:set(v) for i,v in owner.free_source_seeds.items()};native=g.nodes[owner.x.node]['name'];count=0
    for i,row in enumerate(g.nodes):
        if row['operation']=='compact_raw_beta':free[i]=lexical.free_variables(g,row['argument'],free)
    assert owner.outer.acceptance_loaded and owner.A is owner.outer.Rc_background['fields']['E']
    assert owner.incoming_correction is owner.outer.Rc['correction']
    assert owner.incoming_background is owner.outer.Rc_background['history']
    assert owner.incoming_integrals is owner.outer.Rc['integrals'] and owner.initial is owner.outer.initial
    assert all(Algebra(g,{}).at(ref.node)==0 for ref in r(owner.outer.Rc_background['fields']['V']))
    assert all(ref==g.zero for q in owner.outer.patch.loop.zero_inlet.values() for ref in r(q))
    for name in ('N','mu','P0','Pstar','delta'):assert getattr(owner,name) is getattr(owner.outer,name)
    assert tuple(owner.targets)==tuple(owner.residual)==current.compact.ROWS
    for data,order in ((owner.targets.values(),5),(owner.h,5),(owner.background['fields'].values(),6),
        (owner.background['generic'].values(),5),(owner.density.values(),5),(owner.history.values(),5),
        (owner.modified_fields.values(),5),(owner.radial.values(),5),(owner.generic.values(),4),
        (owner.pressure.values(),3),(owner.integrals.values(),3)):
        for q in data:
            assert len(r(q))==order+1
            for ref in r(q):assert lexical.free_variables(g,ref,free)<={'Z',native},ref.node
            count+=order+1
    for q in owner.targets.values():
        for ref in r(q):assert lexical.free_variables(g,ref,free)<={'Z'}
    for n,step in enumerate(owner.implicit['recurrence'],1):
        assert step['order']==n and step['same_full_Jacobian'] and step['only_unknown_highest_row_temporarily_zero']
    assert g.nodes[owner.system.node]['original_branch_transfer']==owner.branch_transfer
    assert g.nodes[owner.system.node]['actual_selected_N']==owner.N.node
    assert owner.branch_transfer['selected_integer']['J']==g.nodes[owner.N.node]['J']
    assert owner.branch_transfer['original_C1_limit_proof']['contraction_at_most']=='1/2'
    for key,q in owner.contributions.items():
        for n,ref in enumerate(r(q)):
            row=g.nodes[ref.node]
            assert row['lower']==g.one.node and row['upper']==owner.x.node
            assert row['finite_N_correction'] and row['original_rate']==str(current.RATES[key])
            assert row['physical_Jacobian_applied_once'] and row['ordinary_slow_Z_derivative_order']==n
    assert owner.memory['p']==g.one
    for key,q in owner.pressure_local.items():
        for ref in r(q):
            row=g.nodes[ref.node]
            assert row['physical_radius_power']==(1 if key=='S0' else 2)
            assert row['physical_Jacobian_applied_once'] and row['nonzero_preceding_S0_S1_retained']
    assert r(owner.modified_fields['P0'])==r(owner.P0)[:6]
    for key,ref in owner.terminal_certificates.items():
        row=g.nodes[ref.node]
        assert row['actual_endpoint_rows']==current.encode(owner.terminal['correction'][key])
        assert row['equivalent_residual_rows']==current.encode(owner.terminal_rhs[key])
        assert row['zero_rows']==[g.zero.node]*6 and row['relative_correction_only']
        assert row['absolute_leading_histories_and_pressure_not_zeroed']
    return dict(passed=True,actual_targets_controls_profile_history_radial_pressure_rows=count,
        actual_Rc_inlets_P0_N_mu_and_nonzero_pressure_history_retained=True,
        actual_control_C5_rows_from_one_full_Jacobian=True,
        original_compact_log_transport_measure_and_relative_only_terminal_theorem=True,
        numerical_high_order_norms_quiet_Rp_global_n1_and_recursion_not_claimed=True)


def original_prefix_and_mu(owner):
    """Compare archived function expressions, not old source range endpoints."""
    old=current.parent.linear.read(current.original_targets.NAME);p=owner.original_prefix
    f=old['actual_C3_density_transport_targets'];og=graph(old['exact_graph_nodes'])
    admission=current.parent.linear.read(current.previous.original_geometry.NAME)
    pg=graph(admission['exact_function_graph_nodes']);pr=admission['exact_current_integral_and_control_graph']
    assert p['actual_C3_target_functions']==f['actual_N_scaled_five_targets']
    assert p['actual_C3_N_scaled_incoming']==f['terminal_N_scaled_histories']
    assert p['actual_C3_amplitude']==f['actual_terminal_amplitude']
    assert p['original_mu_node']==pr['parameters']['mu'] and p['original_N_node']==pr['N']
    # No fixture parameter is used for this original-source recipe identity.
    assert sy.simplify(Algebra(pg,{}).at(p['original_mu_node'])-Algebra(owner.g,{}).at(owner.mu.node))==0
    assert pg.nodes[p['original_N_node']]==owner.g.nodes[owner.N.node]
    Z=sy.Symbol('Z');point=sy.Rational(3,8);mu=sy.Rational(1,13);N=sy.Integer(37)
    A=2+Z+Z**2/17+Z**6/19;H={k:3+i+Z**2/(23+i)+Z**5/(31+i) for i,k in enumerate(current.RATES)}
    nb={owner.g.symbol('Z').node:point,owner.mu.node:mu,owner.N.node:N};ob={}
    for n,ref in enumerate(current.rows(owner.A)):nb[ref.node]=sy.diff(A,Z,n).subs(Z,point)
    for k,q in owner.incoming_correction.items():
        for n,ref in enumerate(current.rows(q)):nb[ref.node]=sy.diff(H[k],Z,n).subs(Z,point)
        for n,i in enumerate(f['terminal_N_scaled_histories'][k]):ob[i]=N*sy.diff(H[k],Z,n).subs(Z,point)
    for n,i in enumerate(f['actual_terminal_amplitude']):ob[i]=sy.diff(A,Z,n).subs(Z,point)
    # The original parameter root is in the shared graph prefix.
    ob[p['original_mu_node']]=mu;ob[p['original_N_node']]=N
    newalg,oldalg=Algebra(owner.g,nb),Algebra(og,ob)
    AA=A*A;expected={'M':N*H['m']/A,'I':N*H['h']/A,'S':N*H['e']/AA,'Cp':N*H['p']/AA,
        current.compact.ROWS[1]:N*(H['k']-A*H['m'])/(mu*AA)}
    for k,q in owner.targets.items():
        for n,ref in enumerate(current.rows(q)):
            want=sy.diff(expected[k],Z,n).subs(Z,point)
            assert sy.simplify(newalg.at(ref.node)-want)==0,(k,n)
            if n<4:assert sy.simplify(oldalg.at(f['actual_N_scaled_five_targets'][k][n])-want)==0,('old prefix',k,n)
    # N-scaled density formula used by the original target transport equals
    # N times the actual physical signed density, with all cross terms.
    E,V,F,B,Ns=sy.symbols('E V F_N B N',nonzero=True)
    dE,dV=F/Ns,B/Ns
    physical=[dV,dE,V*dE+E*dV+dE*dV,2*V*dV+dV**2-E*dE-dE**2/2,E*dE+dE**2/2]
    scaled=[B,F,V*F+E*B+F*B/Ns,2*V*B-E*F+(B**2-F**2/2)/Ns,E*F+F**2/(2*Ns)]
    assert all(sy.expand(Ns*a-b)==0 for a,b in zip(physical,scaled))
    receipt=current.parent.linear.read(current.original_limit.RECEIPT)
    assert receipt['all_passed'] and receipt[current.original_limit.GATE] and receipt['source_family']==owner.family
    return dict(passed=True,original_mu_recipe_and_selected_N_expression_exactly_equal=True,
        actual_new_target_C5_normalization_rows=30,accepted_C3_target_FunctionRef_prefix_rows=20,
        original_amplitude_and_N_scaled_history_expression_identifications_explicit=True,
        original_N_scaled_signed_density_polynomial_identity=True,
        original_profile_and_Duhamel_source_identity_dependency=owner.outer.definition,
        original_unique_branch_and_exact_signed_bump_endpoint_theorem_bound=True,
        no_range_endpoints_or_extra_N_factor_used=True,normalization_calculus_fixture_only=True)


def implicit_calculus():
    g=current.source.FunctionTransportGraph();Z=sy.Symbol('Z');point=sy.Rational(1,7);N=sy.Integer(37);mu=sy.Rational(1,13)
    B=[[1,1,0,0,0],[-sy.Rational(1,4),-sy.Rational(3,4),0,0,0],
       [0,0,1,2,3],[0,0,-2,-5,-9],[0,0,4,9,16]]
    raw=dict(B=B,mu=mu,cross=[sy.Rational(2,3),sy.Rational(3,5)],
        energy=[sy.Rational(5,7),sy.Rational(7,11),sy.Rational(11,13)],
        pressure=[sy.Rational(13,17),sy.Rational(17,19),sy.Rational(19,23)])
    W={k:[[g.constant(v) for v in row] for row in vals] if k=='B' else g.constant(vals) if k=='mu' else [g.constant(v) for v in vals] for k,vals in raw.items()}
    h=[(i+1+Z/(i+5)+Z**2/(i+11)+Z**3/(i+13)+Z**5/(i+17))/1000 for i in range(5)]
    Q=[0,(raw['cross'][0]*h[0]*h[2]+raw['cross'][1]*h[1]*h[4])/mu,0,
        raw['energy'][0]*h[0]**2+raw['energy'][2]*h[1]**2-sum(w*q*q for w,q in zip(raw['energy'],h[2:]))/2,
        sum(w*q*q for w,q in zip(raw['pressure'],h[2:]))/2]
    d=[-sum(b*q for b,q in zip(row,h))-q/N for row,q in zip(B,Q)]
    got,proof=current.implicit_rows(g,W,g.constant(N),[jet(g,q) for q in d],[g.constant(q.subs(Z,point)) for q in h])
    alg=Algebra(g,{g.symbol('Z').node:point})
    for i,q in enumerate(got):
        for n,ref in enumerate(current.rows(q)):assert sy.simplify(alg.at(ref.node)-sy.diff(h[i],Z,n).subs(Z,point))==0,(i,n)
    J=sy.Matrix([[alg.at(q.node) for q in row] for row in proof['Jacobian']])
    inverse=sy.Matrix([[alg.at(q.node) for q in row] for row in proof['inverse']]);assert J*inverse==sy.eye(5)
    for n,pairs in ((4,((4,1,3),(3,2,2))),(5,((5,1,4),(10,2,3)))):
        cross=[g.zero]*5
        for factor,i,j in pairs:
            terms=current.compact.quadratic_action(g,W,[current.rows(q)[i] for q in got],[current.rows(q)[j] for q in got])
            cross=[g.add(a,g.mul(g.constant(factor),b)) for a,b in zip(cross,terms)]
        for k,rhs in enumerate(proof['recurrence'][n-1]['known_rhs']):
            expected=sy.diff(d[k],Z,n).subs(Z,point)+Algebra(g,{g.symbol('Z').node:point}).at(cross[k].node)/N
            assert sy.simplify(Algebra(g,{g.symbol('Z').node:point}).at(rhs.node)-expected)==0,(n,k)
    return dict(passed=True,independent_known_original_quadratic_branch_C5_rows=30,
        exact_full_Jacobian_inverse_and_all_ordinary_cross_terms=True,
        original_fourth_Hessian_factors_4_3_and_fifth_5_10_checked=True,
        fifth_control_rows_nonzero=True,mathematical_fixture_only=True)


def actual_profile_calculus(owner):
    g=owner.g;Z=sy.Symbol('Z');X=sy.Symbol('compact_x',positive=True);point=sy.Rational(3,8)
    A=2+Z+Z**6/19;mu=sy.Rational(1,13);N=sy.Integer(37);h=[1+i+Z**2/(23+i)+Z**5/(31+i) for i in range(5)]
    beta=[1+i+(2+i)*sy.log(X)+(3+i)*sy.log(X)**2/7 for i in range(3)]
    ydiff=lambda q:X*sy.diff(q,X)
    bindings={owner.x.node:X,g.symbol('Z').node:point,owner.mu.node:mu,owner.N.node:N}
    for n,ref in enumerate(current.rows(owner.A)):bindings[ref.node]=sy.diff(A,Z,n).subs(Z,point)
    for i,q in enumerate(owner.h):
        for n,ref in enumerate(current.rows(q)):bindings[ref.node]=sy.diff(h[i],Z,n).subs(Z,point)
    for i,triplet in enumerate(owner.profiles['bumps']):
        for j,ref in enumerate(triplet):bindings[ref.node]=beta[i] if j==0 else ydiff(beta[i]) if j==1 else ydiff(ydiff(beta[i]))
    F=sum(b*q for b,q in zip(beta,h[2:]))/N;G=(beta[0]*h[0]+beta[2]*h[1])/N
    power=X**(-(sy.Rational(1,2)+mu));bodies=dict(F=F,G=G,original_E=A*power,delta_E=A*F,delta_V=A*G,E=A*(power+F),V=A*G)
    alg=Algebra(g,bindings);count=0
    for j,row in enumerate(owner.profiles['radial_y_rows']):
        for k,q in row.items():
            body=bodies[k] if j==0 else ydiff(bodies[k]) if j==1 else ydiff(ydiff(bodies[k]))
            for n,ref in enumerate(current.rows(q)):
                assert sy.simplify(alg.at(ref.node)-sy.diff(body,Z,n).subs(Z,point))==0,(j,k,n);count+=1
    dE,dV=A*F,A*G;E=A*power
    D=dict(m=dV,h=dE,k=(E+dE)*dV,e=dV*dV-E*dE-dE*dE/2,p=E*dE+dE*dE/2)
    for k,q in owner.density.items():
        for n,ref in enumerate(current.rows(q)):assert sy.simplify(alg.at(ref.node)-sy.diff(D[k],Z,n).subs(Z,point))==0,(k,n)
    # Here the actual original raw beta argument is evaluated, without the
    # abstract beta bindings used above. All original supports lie strictly
    # inside (1,2), so these three physical radial rows are flat at both ends.
    for x in (sy.Integer(1),sy.Integer(2)):
        end=Algebra(g,{owner.x.node:x,owner.W['normal'].node:sy.Integer(1)})
        for triplet in owner.profiles['bumps']:
            for ref in triplet:assert sy.simplify(end.at(ref.node))==0,(x,ref.node)
    return dict(passed=True,actual_compact_profile_value_first_second_radial_C5_rows=count,
        actual_five_signed_density_C5_rows=30,actual_original_raw_beta_end_flat_rows=18,
        log_radius_chain_and_1_over_x_and_all_amplitude_cross_terms=True,
        abstract_profile_calculus_fixture_not_original_field_value=True)


def actual_inlet_and_endpoint_coefficient(owner):
    g=owner.g;Z=sy.Symbol('Z');point=sy.Rational(3,8);Rc=sy.Rational(7,3);N=sy.Integer(37);mu=sy.Rational(1,13)
    A=2+Z+Z**6/19;H={k:3+i+Z**2/(23+i)+Z**5/(31+i) for i,k in enumerate(current.RATES)}
    L={k:7+i+Z**2/(41+i)+Z**6/(43+i) for i,k in enumerate(current.RATES)};P0=-11-Z**6/47
    h=[1+i+Z**2/(53+i)+Z**5/(61+i) for i in range(5)]
    bindings={g.symbol('Z').node:point,owner.x.node:sy.Integer(1),owner.Rc.node:Rc,owner.N.node:N,owner.mu.node:mu,owner.W['normal'].node:1}
    bind=lambda q,body:[bindings.__setitem__(ref.node,sy.diff(body,Z,n).subs(Z,point)) for n,ref in enumerate(current.rows(q))]
    bind(owner.A,A);bind(owner.P0,P0)
    for k in current.RATES:bind(owner.incoming_background[k],L[k]);bind(owner.incoming_correction[k],H[k])
    for i,q in enumerate(owner.h):bind(q,h[i])
    for k,q in owner.incoming_integrals.items():bind(q,(13 if k=='S0' else 17)+Z**3/67)
    Q=19+Z**3/71;Qy=23+Z**3/73;B0=29+Z**3/79;B1=31+Z**3/83
    bind(owner.physical['Q'],Q);bind(owner.physical['Q_y'],Qy)
    bind(owner.initial['Q_plus_Q_y'],B0);bind(owner.initial['R_Q_y'],B1)
    alg=Algebra(g,bindings);count=0
    for k,q in owner.inlet['correction'].items():
        for n,ref in enumerate(current.rows(q)):assert sy.simplify(alg.at(ref.node)-sy.diff(H[k],Z,n).subs(Z,point))==0;count+=1
    bodies={k:L[k]+H[k] for k in current.RATES};bodies.update(E=A,V=sy.Integer(0),P0=P0,pressure=P0+L['p']+H['p'])
    for k,q in owner.inlet['modified_fields'].items():
        for n,ref in enumerate(current.rows(q)):assert sy.simplify(alg.at(ref.node)-sy.diff(bodies[k],Z,n).subs(Z,point))==0,(k,n);count+=1
    expected=dict(m13=-(13+Z**3/67)/2+Q+Qy-B0,m15=(17+Z**3/67)/2-Rc*Qy+B1)
    for k,q in owner.inlet['pressure'].items():
        for n,ref in enumerate(current.rows(q)):assert sy.simplify(alg.at(ref.node)-sy.diff(expected[k],Z,n).subs(Z,point))==0,(k,n);count+=1
    # Test the actual endpoint-to-residual expressions against independent
    # signed full-bump moment coefficients, with nonzero arbitrary controls.
    B=[[1,1,0,0,0],[-sy.Rational(1,4),-sy.Rational(3,4),0,0,0],
       [0,0,1,2,3],[0,0,-2,-5,-9],[0,0,4,9,16]]
    C=[sy.Rational(2,3),sy.Rational(3,5)];Ew=[sy.Rational(5,7),sy.Rational(7,11),sy.Rational(11,13)];Pw=[sy.Rational(13,17),sy.Rational(17,19),sy.Rational(19,23)]
    for i,row in enumerate(owner.W['B']):
        for j,ref in enumerate(row):bindings[ref.node]=B[i][j]
    for key,vals in (('cross',C),('energy',Ew),('pressure',Pw)):
        for ref,value in zip(owner.W[key],vals):bindings[ref.node]=value
    QS=Ew[0]*h[0]**2+Ew[2]*h[1]**2-sum(w*q*q for w,q in zip(Ew,h[2:]))/2
    QP=sum(w*q*q for w,q in zip(Pw,h[2:]))/2;cross=C[0]*h[0]*h[2]+C[1]*h[1]*h[4]
    want=dict(m=(H['m']+A*(h[0]+h[1])/N)/2,
        h=(H['h']+A*sum(w*q for w,q in zip(B[2][2:],h[2:]))/N)*2**(-sy.Rational(3,2)),
        e=(H['e']+A*A*sum(w*q for w,q in zip(B[3][2:],h[2:]))/N+A*A*QS/N**2)/2,
        p=H['p']+A*A*sum(w*q for w,q in zip(B[4][2:],h[2:]))/N+A*A*QP/N**2,
        k=(H['k']+A*A*(h[0]+h[1]+mu*(B[1][0]*h[0]+B[1][1]*h[1]))/N+A*A*cross/N**2)*2**(-sy.Rational(3,2)))
    alg=Algebra(g,bindings)
    for k,q in owner.terminal_rhs.items():
        for n,ref in enumerate(current.rows(q)):assert sy.simplify(alg.at(ref.node)-sy.diff(want[k],Z,n).subs(Z,point))==0,('endpoint',k,n)
    return dict(passed=True,actual_nonzero_Rc_correction_modified_pressure_inlet_rows=count,
        actual_terminal_rhs_original_signed_moment_coefficient_C5_rows=30,
        joint_k_minus_A_m_and_mu_divided_row_and_N_factors_retained=True,
        terminal_integral_identity_uses_pinned_original_exact_bump_theorem=True,
        mathematical_fixture_only=True)


def run(owner):
    began=time.monotonic();owner.assert_graph();raw=current.parent.linear.read(current.NAME)
    assert raw['true_compact_C5_bindings']==owner.snapshot and raw['exact_graph_nodes']==owner.g.nodes[:len(raw['exact_graph_nodes'])]
    checks={}
    for name,fn in (('actual_source',actual_source),('original_prefix_and_mu',original_prefix_and_mu),
        ('actual_profiles',actual_profile_calculus),('actual_inlet_and_endpoint',actual_inlet_and_endpoint_coefficient)):
        checks[name]=fn(owner);print('TRUE_COMPACT_'+name.upper()+'_CHECKED',flush=True)
    checks['implicit_calculus']=implicit_calculus();print('TRUE_COMPACT_IMPLICIT_C5_CALCULUS_CHECKED',flush=True)
    owner.assert_graph()
    for name,digest in owner.hashes.items():assert current.sha(name)==digest,name
    receipt=dict(all_passed=True,source_family=owner.family,checks=checks,
        **dict.fromkeys(current.GATES,True),**dict.fromkeys(current.OPEN,False),
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8',newline='\n')
    print('TRUE_COMPACT_C5_AND_2RC_PRESSURE_C3_CHECKED',flush=True);return receipt


if __name__=='__main__':raise SystemExit('run(owner) reuses the actual checked Rc source')
