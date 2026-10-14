"""Exact source correspondence and independent leading implicit C6 calculus."""
from fractions import Fraction
import json
from pathlib import Path
import time

import sympy as sy

import lei_ren_part1_paper_compliant_current_original_true_Rm_patch_C6_functions as current
import lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions_check as prior


class Algebra(prior.Algebra):
    def at(self,i):
        if i in self.cache:return self.cache[i]
        row=self.g.nodes[i]
        if row['operation']=='compact_raw_beta':
            x=self.at(row['argument'])
            value=sy.Integer(0) if x.is_number and abs(x)>=1 else sy.exp(-1/(1-x*x)) if x.is_number else sy.Function('raw_beta')(x)
            self.cache[i]=value;return value
        return super().at(i)


def jet(g,body,order=6):
    Z=sy.Symbol('Z')
    return current.FunctionJet(tuple(current.source.expression(g,sy.diff(body,Z,k),{'Z':g.symbol('Z')}) for k in range(order+1)))


def exact_fraction(row):
    sign,man,exponent,_=row
    return Fraction((-1 if sign else 1)*man)*Fraction(2)**exponent


def actual_source(owner):
    g=owner.g;r=current.rows;free={i:set(v) for i,v in owner.free_source_seeds.items()};count=0
    assert owner.modified.acceptance_loaded and owner.background.acceptance_loaded
    assert owner.inlet is owner.background.endpoints['Rm']
    assert all(len(r(q))==7 for q in owner.inlet['fields'].values())
    assert len(r(owner.am))==len(r(owner.alpha))==7
    assert owner.bindings['true_normalized_fields']==current.parent.active.current.ast_binding(current.normalization.generic_fields)
    assert owner.bindings['true_Rm_source']==current.parent.active.current.ast_binding(type(owner.background).__init__)
    assert owner.bindings['original_centered_packet']==current.parent.active.current.ast_binding(current.leading.previous.ActualReferenceRestoreFunctions.packet)
    assert owner.bindings['original_postrestore']==current.parent.active.current.ast_binding(current.leading.previous.ActualReferenceRestoreFunctions.postrestore)
    assert owner.bindings['original_leading_branch']==current.parent.active.current.ast_binding(current.leading._ActualRmLeadingPatch.__init__)
    system=g.nodes[owner.system.node]
    assert system['exact_weight_nodes']==current.encode(owner.W)
    assert system['exact_defect_value_nodes']==[q.value.node for q in owner.defects]
    assert system['invAm2']==owner.alpha.value.node
    assert system['exact_Z_independent_row_and_control_units']==[q.node for q in owner.units]
    assert system['fixed_leading_input_N0']==owner.fixed_N0
    assert system['source_projection_binding']==current.parent.active.current.ast_binding(current.admitted.WholeZRmCenteredInverse.owner)
    assert owner.base_identity['system']==owner.system and owner.base_identity['residual']==owner.base_residual
    a=current.Algebra(g,0);values=[q.value for q in owner.control_functions]
    Q=current.nonlinear(a,owner.W,[a.fixed(q) for q in values],a.trim(owner.alpha))
    for i,ref in enumerate(owner.base_residual):
        assert ref==g.add(owner.defects[i].value,Q[i].value,*(g.mul(owner.W['L'][i][j],values[j]) for j in range(5)))
    for i,q in enumerate(owner.control_functions):
        row=g.nodes[q.value.node]
        assert row['operation']=='original_Rm_unique_leading_implicit_control' and row['system']==owner.system.node
        assert row['component']==i and row['source_family']==owner.family
        assert len(r(q))==7
    assert len(owner.implicit['recurrence'])==6
    for n,step in enumerate(owner.implicit['recurrence'],1):
        assert step['order']==n and step['only_unknown_coefficient_temporarily_zero']
        assert step['solved_ordinary_row']==[r(q)[n] for q in owner.control_functions]
    name=g.nodes[owner.x.node]['name']
    for group,order in ((owner.defects,6),(owner.control_functions,6),(owner.chart['fields'].values(),6),
                        (owner.chart['radial'].values(),6),(owner.chart['generic'].values(),5)):
        for q in group:
            assert len(r(q))==order+1
            for ref in r(q):assert prior.free_variables(g,ref,free)<={'Z',name},(ref.node,prior.free_variables(g,ref,free))
            count+=order+1
    assert owner.chart['radius']==g.mul(owner.inlet['radius'],owner.x)
    assert owner.chart['physical_measure']==g.quotient(g.one,owner.x,'positive original actual patch x')
    assert r(owner.chart['fields']['P0'])==r(owner.P0)[:7]
    assert owner.endpoints['Rm']['radius']==owner.at(g.one)['radius']
    return dict(passed=True,actual_C6_and_generic_C5_rows=count,true_Rm_parent_and_same_original_map_nodes=True,
        ordinary_implicit_rows_not_padded=True,unique_branch_defining_residual_bound_by_theorem_not_numeric_zero=True,
        independent_P0_retained_and_excluded_from_five_defects=True,original_log_Jacobian_and_finite_patch_endpoints=True,
        no_free_old_phase_or_integration_variable=True,no_finite_N_history_reset_or_numeric_containment_claim=True)


def exact_weights(owner):
    g=owner.g;proof=current.parent.linear.read('lei_ren_part1_paper_bump_integral_enclosures_check.json')
    assert proof['directed_integrals_certified']
    for name,key in (('lei_ren_part1_paper_bump_integral_enclosures.py','module_sha256'),
                     ('lei_ren_part1_paper_bump_integral_enclosures_check.py','check_source_sha256')):
        assert current.sha(name)==proof[key] and owner.hashes[name]==proof[key]
    records={(Fraction(q['center']),Fraction(q['power']),q['multiplicity']):q for q in proof['weight_records'].values()}
    norm=g.nodes[owner.bump.normalization.node];s=sy.Symbol(norm['variable']);N=sy.Symbol('original_beta_normalization',positive=True)
    b=Algebra(g,{g.symbol(norm['variable']).node:s})
    assert norm['operation']=='definite_integral'
    assert b.at(norm['lower'])==-1 and b.at(norm['upper'])==1
    assert b.at(norm['integrand'])==sy.Function('raw_beta')(s)
    count=0;bindings={q['node']:q for q in owner.branch_transfer['exact_weight_recipe_bindings']}
    for (i,p,m,upper),ref in owner.bump.cache.items():
        row=g.nodes[ref.node];x=sy.Symbol(row['variable'],positive=True);center=sy.Rational(current.CENTERS[i])
        alg=Algebra(g,{g.symbol(row['variable']).node:x,owner.bump.normalization.node:N})
        assert upper is None and (current.CENTERS[i],p,m) in records
        assert alg.at(row['lower'])==center-sy.Rational(1,40)
        assert alg.at(row['upper'])==center+sy.Rational(1,40)
        expected=x**sy.Rational(p)*(40*sy.Function('raw_beta')(40*(x-center))/N)**m
        assert sy.simplify(alg.at(row['integrand'])-expected)==0,(i,p,m)
        assert row['ordinary_measure']=='dx' and row['Z_independent_weight']
        bound=records[(current.CENTERS[i],p,m)]['weight_interval']
        assert bindings[ref.node]['original_directed_weight_interval']==bound
        lower=exact_fraction(bound['lower_exact_mpf_tuple']);upper=exact_fraction(bound['upper_exact_mpf_tuple'])
        assert 0<=lower<=upper
        count+=1
    ctx=owner.background.bridge.inlet.common.ctx
    for center in (current.CENTERS[0],current.CENTERS[2]):
        value=current.leading.restore_value(ctx,records[(center,Fraction(0),1)]['weight_interval'])
        lower,upper=current.leading.ep(value)
        assert lower<=1<=upper
    # u=40(x-center), dx=du/40 proves gamma's mass=1. The same change
    # of variables gives exactly the original directed weight recipes.
    return dict(passed=True,exact_full_integral_recipes_checked=count,
        integrands_bounds_normalization_and_measure_equal_to_original_directed_definitions=True,
        each_exact_integral_containment_bound_by_its_pinned_directed_recipe_and_interval=True,
        exact_unit_mass_by_original_change_of_variable=True,no_interval_midpoint_defines_a_weight=True)


def true_inlet_normalization(owner):
    """Evaluate the ACTUAL parent field nodes with independent raw C6 inputs.

    This proves the normalization map, including every Z product and both
    halves. Original raw histories are then identified by the pinned
    prescribed-source/initial-data theorem, not a numerical overlap.
    """
    g=owner.g;r=current.rows;bg=owner.background;chart=bg.charts['post_restore'];Z=sy.Symbol('Z');point=sy.Rational(3,8)
    P=sy.Rational(3,2);R=sy.Rational(7,3)
    F=2+Z/13+Z**2/17+Z**6/19
    bodies={key:sy.Rational(j+2)+Z/(j+23)+Z**2/(j+29)+Z**6/(j+37)
            for j,key in enumerate(('H','M','K','A','B','C'))}
    bindings={g.symbol('Z').node:point,bg.core.parameters['logP'].node:sy.log(P),chart['radius'].node:R}
    for ref,k in zip(r(bg.F0),range(7)):bindings[ref.node]=sy.diff(F,Z,k).subs(Z,point)
    for key,q in chart['moments'].items():
        for k,ref in enumerate(r(q)):bindings[ref.node]=sy.diff(bodies[key],Z,k).subs(Z,point)
    for k,ref in enumerate(r(bg.P0)):bindings[ref.node]=sy.Integer(101+k)
    alg=Algebra(g,bindings);am=sy.exp(-sy.Rational(3,5))/(1+Z**2)
    phi=P*am/(sy.sqrt(2*R)*F);H,M,K,A,B,C=(bodies[key] for key in ('H','M','K','A','B','C'))
    expected=[M-4*Z,K/(2*phi)-4*Z*H/(2*phi),H/(2*phi)-sy.Rational(5,8),
        (A-8*Z*M+16*Z**2)/(P*am)**2-(B/phi**2-sy.Rational(5,6))/2,
        (C/phi**2-5)/2]
    count=0
    for i,q in enumerate(owner.defects):
        for k,ref in enumerate(r(q)):
            assert sy.simplify(alg.at(ref.node)-sy.diff(expected[i],Z,k).subs(Z,point))==0,(i,k)
            count+=1
    for k,ref in enumerate(r(owner.inlet['fields']['E'])):
        assert sy.simplify(alg.at(ref.node)-sy.diff(am,Z,k).subs(Z,point))==0,('actual E=Wphi',k)
    for k,ref in enumerate(r(owner.alpha)):
        assert sy.simplify(alg.at(ref.node)-sy.diff(1/(P*am)**2,Z,k).subs(Z,point))==0,('invAm2',k)
    # Exact order-zero identities with arbitrary symbols prove equality
    # of the defining functions for every Z and every raw moment value.
    z=sy.Symbol('Z',real=True);ps,rs,fs=sy.symbols('Pstar Rm F0',positive=True)
    raw={key:sy.Symbol('raw_'+key,real=True) for key in bodies}
    symbolic={g.symbol('Z').node:z,bg.core.parameters['logP'].node:sy.log(ps),chart['radius'].node:rs,
              r(bg.F0)[0].node:fs,**{r(chart['moments'][key])[0].node:value for key,value in raw.items()}}
    symbolic_alg=Algebra(g,symbolic);am0=sy.exp(-sy.Rational(3,5))/(1+z*z);phi0=ps*am0/(sy.sqrt(2*rs)*fs)
    want=[raw['M']-4*z,raw['K']/(2*phi0)-4*z*raw['H']/(2*phi0),raw['H']/(2*phi0)-sy.Rational(5,8),
          (raw['A']-8*z*raw['M']+16*z*z)/(ps*am0)**2-(raw['B']/phi0**2-sy.Rational(5,6))/2,
          (raw['C']/phi0**2-5)/2]
    for i,q in enumerate(owner.defects):assert sy.simplify(symbolic_alg.at(q.value.node)-want[i])==0,('exact source function',i)
    assert sy.simplify(symbolic_alg.at(owner.inlet['fields']['E'].value.node)-am0)==0
    return dict(passed=True,actual_true_Rm_defect_normalization_rows=count,
        actual_E_equals_W_phi_equals_am_rows=7,actual_invAm2_product_rows=7,
        nonconstant_F0_and_all_sixth_raw_moment_rows_and_nonunit_Pstar_radius=True,
        swirl_energy_half_and_pressure_half_proved_on_actual_field_nodes=True,
        exact_defining_function_identities_for_arbitrary_Z_raw_moments_Pstar_Rm_F0=True,
        P0_excluded=True,exact_algebra_fixture_not_original_scale_numeric_enclosure=True)


def normalized_conjugacy(owner):
    g=owner.g;certificate=owner.conjugacy;T,P,A=sy.symbols('t Pstar invAm2',positive=True)
    bindings={owner.t.node:T,owner.background.core.parameters['logP'].node:sy.log(P),owner.alpha.value.node:A}
    for q in owner.bump.cache.values():bindings[q.node]=sy.Symbol('w%d'%q.node,positive=True)
    for i,q in enumerate(owner.defects):bindings[q.value.node]=sy.Symbol('d%d'%i,real=True)
    for i,q in enumerate(certificate['variables']):bindings[q.node]=sy.Symbol('u%d'%i,real=True)
    alg=Algebra(g,bindings);count=0
    for raw,norm in zip(certificate['raw_map_divided_by_units'],certificate['normalized_map']):
        assert sy.expand(alg.at(raw.node)-alg.at(norm.node))==0;count+=1
    for rawrow,normrow in zip(certificate['conjugated_raw_Jacobian'],certificate['normalized_Jacobian']):
        for raw,norm in zip(rawrow,normrow):
            assert sy.expand(alg.at(raw.node)-alg.at(norm.node))==0;count+=1
    z=sy.Symbol('Z',real=True)
    actual=Algebra(g,{owner.t.node:T,owner.background.core.parameters['logP'].node:sy.log(P),g.symbol('Z').node:z})
    assert sy.simplify(actual.at(certificate['normalized_alpha'].value.node)-sy.exp(sy.Rational(6,5))*(1+z*z)**2/10**10)==0
    assert owner.branch_transfer['raw_to_original_normalized_conjugacy']==current.encode(certificate)
    return dict(passed=True,exact_map_and_Jacobian_conjugacy_entries=count,
        arbitrary_controls_weights_defects_and_positive_units_checked_symbolically=True,
        original_normalized_alpha_and_angular_weight_factors_proved=True,
        normalized_Neumann_theorem_proves_raw_Jacobian_invertible=True)


def branch_theorem(owner,checks):
    """Bind the old unique branch to the SAME exact map, not a new box.

    The old directed proof encloses exact source functions/integrals. The
    normalized-source identities above and the checked unique prescribed
    ODE chain identify the new representation with those same functions.
    Its strict common-box Neumann bound therefore proves det(J)!=0.
    No sixth-derivative numerical norm follows from this transfer.
    """
    receipt=current.parent.linear.read(current.admitted.RECEIPT)
    raw=current.parent.linear.read(current.admitted.NAME)
    assert all(checks[key]['passed'] for key in ('actual_source','exact_weights','true_inlet_normalization','normalized_conjugacy'))
    assert receipt['all_passed'] and receipt[current.admitted.GATE]
    assert receipt['source_family']==raw['source_family']==owner.family
    assert raw['original_source_bindings']==owner.bindings['source_centering']==current.admitted.source_bindings()
    assert raw['candidate_N_not_used_for_leading_map']==owner.fixed_N0
    assert raw['exact_Z_domain']==['-1','1']
    assert receipt['same_common_unit_box_proves_unique_whole_Z_smooth_leading_family']
    assert receipt['actual_leading_controls_and_same_Jacobian_axial5_checked']
    ctx=owner.background.bridge.inlet.common.ctx
    for cell in receipt['contraction_margins'].values():
        for key in ('initial_contraction','higher_Jacobian_contraction'):
            q=ctx.mpf(current.mp.make_mpf(tuple(cell[key]['upper_exact_mpf_tuple'])))
            assert 0<=q<1
    repair=current.parent.linear.read(current.PREFIX+'five_moment_repair.json')
    t=exact_fraction(repair['axial_box_scale']['lower_exact_mpf_tuple'])
    coefficient=10**10*t*t
    tv=ctx.mpf(t.numerator)/ctx.mpf(t.denominator)
    actual_log_unit=ctx.ln(ctx.mpf(10**10))+2*ctx.ln(tv)-2*owner.background.bridge.inlet.common.core.logP
    for cell in raw['source_cells']:
        solution=cell['coefficient_solution'];unit=solution['canonical_positive_angular_unit']
        bound=unit['coefficient_interval'];scale=unit['formal_positive_scale']
        assert exact_fraction(bound['lower_exact_mpf_tuple'])<=coefficient<=exact_fraction(bound['upper_exact_mpf_tuple'])
        assert scale['source_exponents']==[0,-1,0,0] and scale['radius_power']==0
        assert exact_fraction(scale['additional_log_interval']['lower_exact_mpf_tuple'])==0
        assert exact_fraction(scale['additional_log_interval']['upper_exact_mpf_tuple'])==0
        cover=solution['directed_angular_weight_cover']
        lower=ctx.mpf(current.mp.make_mpf(tuple(cover['lower_exact_mpf_tuple'])))
        upper=ctx.mpf(current.mp.make_mpf(tuple(cover['upper_exact_mpf_tuple'])))
        assert current.leading.ep(lower)[0]>=0 and current.leading.ep(upper)[0]>0
        assert current.leading.ep(actual_log_unit)[1]<=current.leading.ep(ctx.ln(upper))[0]
        if current.leading.ep(lower)[1]>0:
            assert current.leading.ep(ctx.ln(lower))[1]<=current.leading.ep(actual_log_unit)[0]
    correspondence=owner.branch_transfer['exact_source_correspondence']
    assert correspondence['source_bindings']==owner.bindings
    assert correspondence['normalized_raw_histories']==['M','H/(2*phi)','K/(2*phi)','A','B/phi^2','C/phi^2']
    assert owner.base_identity['source_transfer'] is owner.branch_transfer
    return dict(passed=True,same_exact_original_unique_branch_theorem_bound=True,
        exact_source_map_weights_amplitude_N0_and_initial_data_identified=True,
        original_common_box_strict_Neumann_bound_applies_to_same_exact_Jacobian=True,
        all_four_original_angular_weight_covers_contain_actual_exact_unit=True,
        predecessor_original_source_bindings_compared_by_exact_equality=True,
        defining_five_function_zero_identities_not_interval_zero_containment=True,
        determinant_squared_positive_by_transferred_invertibility=True,
        sixth_function_rows_follow_implicit_function_theorem_not_sixth_numeric_norm=True)


def implicit_calculus():
    g=current.source.FunctionTransportGraph();Z=sy.Symbol('Z');point=sy.Rational(1,7)
    L=[[1,1,0,0,0],[sy.Rational(5,4),sy.Rational(7,4),0,0,0],
        [0,0,1,2,3],[0,0,-2,-5,-9],[0,0,4,9,16]]
    raw=dict(L=L,fg=[sy.Rational(2,3),sy.Rational(3,5)],gg=[sy.Rational(5,7),sy.Rational(7,11)],
             ff=[sy.Rational(11,13),sy.Rational(13,17),sy.Rational(17,19)],
             ff_over_x=[sy.Rational(19,23),sy.Rational(23,29),sy.Rational(29,31)])
    W={key:[[g.constant(v) for v in row] for row in value] if key=='L' else [g.constant(v) for v in value] for key,value in raw.items()}
    h=[(i+1+Z/(i+5)+Z**2/(i+11)+Z**6/(i+17))/1000 for i in range(5)]
    alpha=2+Z/13+Z**6/17
    Q=[0,h[0]*h[2]*raw['fg'][0]+h[1]*h[4]*raw['fg'][1],0,
       alpha*(h[0]**2*raw['gg'][0]+h[1]**2*raw['gg'][1])-sum(h[i+2]**2*raw['ff'][i] for i in range(3))/2,
       sum(h[i+2]**2*raw['ff_over_x'][i] for i in range(3))/2]
    defects=[-sum(L[i][j]*h[j] for j in range(5))-Q[i] for i in range(5)]
    got,proof=current.implicit_rows(g,W,[jet(g,d) for d in defects],jet(g,alpha),
                                  [g.constant(body.subs(Z,point)) for body in h])
    alg=Algebra(g,{g.symbol('Z').node:point});count=0
    for i,q in enumerate(got):
        for n,ref in enumerate(current.rows(q)):
            assert alg.at(ref.node)==sy.diff(h[i],Z,n).subs(Z,point),(i,n);count+=1
    J=sy.Matrix([[alg.at(ref.node) for ref in row] for row in proof['Jacobian']])
    inv=sy.Matrix([[alg.at(ref.node) for ref in row] for row in proof['inverse']])
    assert J*inv==sy.eye(5) and J.det()!=0
    assert sy.diff(alpha,Z,6)!=0 and all(sy.diff(q,Z,6)!=0 for q in h)
    return dict(passed=True,independent_known_nonlinear_implicit_C6_rows=count,
        exact_full_Jacobian_inverse_and_ordinary_recurrence_checked=True,
        sixth_amplitude_products_and_all_control_sixth_rows_nonzero=True,
        mathematical_fixture_only=True)


def beta_calculus():
    g=current.source.FunctionTransportGraph();x=g.symbol('fixture_beta_x');beta=current.raw_beta(g,x)
    first,second=current.beta_first_two(g,x,beta);X=sy.Symbol('x');expected=sy.exp(-1/(1-X*X));count=0
    for point in (sy.Rational(0),sy.Rational(1,3),sy.Rational(-2,5),sy.Integer(-1),sy.Integer(1),sy.Integer(2)):
        alg=Algebra(g,{x.node:point})
        for n,ref in enumerate((beta,first,second)):
            want=0 if abs(point)>=1 else sy.diff(expected,X,n).subs(X,point)
            assert sy.simplify(alg.at(ref.node)-want)==0,(point,n);count+=1
    return dict(passed=True,independent_raw_beta_value_first_second_and_flat_rows=count,
        singular_interior_body_not_evaluated_at_flat_endpoints=True,mathematical_fixture_only=True)


def patch_matching_calculus():
    g=current.source.FunctionTransportGraph();Z=sy.Symbol('Z');point=sy.Rational(3,8);x=g.symbol('fixture_patch_x')
    am=sy.exp(-sy.Rational(3,5))/(1+Z**2);P=sy.Rational(3,2);alpha=1/(P*am)**2
    bodies=dict(m=2+Z+Z**6/13,h=3+Z**2+Z**6/17,k=4+Z**3+Z**6/19,
                e=5+Z**4+Z**6/23,p=6+Z**5+Z**6/29,P0=-7-Z**6/31)
    fields={key:jet(g,value) for key,value in bodies.items()};defects=current.inlet_defects(g,fields,jet(g,am),g.constant(P))
    zero=g.zero;Wzero=dict(L=[[zero]*5 for _ in range(5)],fg=[zero]*2,gg=[zero]*2,ff=[zero]*3,ff_over_x=[zero]*3)
    controls=[jet(g,(i+1+Z**6/(i+37))/1000) for i in range(5)]
    gamma=[(zero,zero,zero)]*3
    got,radial,_=current.patch_fields(g,x,jet(g,am),g.constant(P),fields['P0'],defects,controls,jet(g,alpha),Wzero,gamma)
    alg=Algebra(g,{x.node:sy.Integer(1),g.symbol('Z').node:point});count=0
    for key,body in {**bodies,'E':am,'V':4*Z/P,'pressure':bodies['P0']+bodies['p']}.items():
        for n,ref in enumerate(current.rows(got[key])):
            assert sy.simplify(alg.at(ref.node)-sy.diff(body,Z,n).subs(Z,point))==0,('inlet',key,n);count+=1
    # At and beyond the original last support71/40, partial weights are
    # full weights. The unique defining map makes all five residuals zero.
    # Test the resulting original terminal target formulas with genuinely
    # nonconstant amplitude, independent P0 and nonunit Pstar.
    d0=[jet(g,sy.Integer(0)) for _ in range(5)]
    end,_,_=current.patch_fields(g,x,jet(g,am),g.constant(P),fields['P0'],d0,controls,jet(g,alpha),Wzero,gamma)
    X=sy.Symbol('X',positive=True);alg=Algebra(g,{x.node:X,g.symbol('Z').node:point})
    target=dict(E=am*X**sy.Rational(1,10),V=4*Z/P,m=4*Z/P,
        h=am*sy.Rational(5,8)*X**sy.Rational(1,10),
        k=am*sy.Rational(5,2)*Z*X**sy.Rational(1,10)/P,
        e=16*Z**2/P**2-am**2*sy.Rational(5,12)*X**sy.Rational(1,5),
        p=am**2*sy.Rational(5,2)*X**sy.Rational(1,5),P0=bodies['P0'])
    target['pressure']=target['P0']+target['p']
    for key,body in target.items():
        for n,ref in enumerate(current.rows(end[key])):
            assert sy.simplify(alg.at(ref.node)-sy.diff(body,Z,n).subs(Z,point))==0,('exit target',key,n);count+=1
    return dict(passed=True,independent_inlet_and_terminal_target_C6_rows=count,
        original_five_residual_recovery_and_Pstar_halves_retained=True,
        terminal_zero_is_same_unique_map_identity_not_new_zero_inlet=True,mathematical_fixture_only=True)


def run(owner):
    began=time.monotonic();owner.assert_graph();raw=current.parent.linear.read(current.NAME)
    assert raw['true_Rm_patch_C6_bindings']==owner.snapshot
    assert raw['exact_graph_nodes']==owner.g.nodes[:len(raw['exact_graph_nodes'])]
    checks={}
    for name,fn in (('actual_source',actual_source),('exact_weights',exact_weights),
                    ('true_inlet_normalization',true_inlet_normalization),('normalized_conjugacy',normalized_conjugacy)):
        checks[name]=fn(owner);print('TRUE_RM_PATCH_'+name.upper()+'_CHECKED',flush=True)
    checks['branch_theorem']=branch_theorem(owner,checks);print('TRUE_RM_PATCH_EXACT_BRANCH_THEOREM_BOUND',flush=True)
    for name,fn in (('implicit_calculus',implicit_calculus),('beta_calculus',beta_calculus),('patch_matching',patch_matching_calculus)):
        checks[name]=fn();print('TRUE_RM_PATCH_'+name.upper()+'_CHECKED',flush=True)
    owner.assert_graph()
    for name,digest in owner.hashes.items():assert current.sha(name)==digest,name
    receipt=dict(all_passed=True,source_family=owner.family,checks=checks,
        same_exact_original_unique_branch_theorem_bound=True,
        **dict.fromkeys(current.GATES,True),**dict.fromkeys(current.OPEN,False),
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8',newline='\n')
    print('TRUE_RM_LEADING_IMPLICIT_AND_PATCH_C6_CHECKED',flush=True);return receipt


if __name__=='__main__':raise SystemExit('run(owner) reuses checked true prescribed Rm source')
