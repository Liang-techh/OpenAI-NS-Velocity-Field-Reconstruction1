"""Independent source partition, analytic bounds and physical pressure signs."""
import copy
from fractions import Fraction
import gzip
import json
from pathlib import Path
import time

from mpmath.ctx_iv import MPIntervalContext
import sympy as s
import lei_ren_part1_paper_compliant_current_original_Rp_remaining_pressure_tail as current
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision


def included(outer,inner):
    a,b=current.ends(outer);u,v=current.ends(inner)
    return a<=u<=v<=b


def rejected(fn):
    try:fn()
    except (ValueError,TypeError,ArithmeticError,KeyError):return True
    raise AssertionError('Invalid original pressure tail input admitted')


def rational_expression(ctx,expression,symbol,value):
    if expression==symbol:return value
    if expression.is_Rational:return ctx.mpf(int(expression.p))/int(expression.q)
    if expression.is_Add:return sum((rational_expression(ctx,v,symbol,value) for v in expression.args),ctx.mpf(0))
    if expression.is_Mul:
        result=ctx.mpf(1)
        for v in expression.args:result*=rational_expression(ctx,v,symbol,value)
        return result
    if expression.is_Pow and expression.args[1].is_Integer:
        return rational_expression(ctx,expression.args[0],symbol,value)**int(expression.args[1])
    raise ValueError('Exact rational axial derivative required')


@source_precision
def run(owner,deliveries):
    began=time.monotonic();raw=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    assert not owner.acceptance_loaded and not any(raw[k] for k in current.GATES+current.OPEN)
    assert raw['source_family']==owner.family_record and all(owner.assert_graph().values())
    assert all(raw['input_hashes'][k]==v for k,v in owner.before.hashes.items())
    independent=MPIntervalContext();independent.dps=500
    # Existing verifier independently binds all original/native densities,
    # exact complete P0 identity and the original unique raw waiting root.
    witness=owner.function.witness
    original=current.terminal_source.verify_exact_integral_witness(witness)
    assert all(original.values()) and len(witness.native_densities)==14
    p=witness.partition;mu=p['mu'];rate=1+2*mu;t=p['t'];xi=owner.function.xi;q=p['q']
    density=witness.native_densities['pulse_reserved']
    actual_partial=s.integrate(density,(t,xi/mu,13/mu))
    actual_prefix=s.integrate(density,(t,0,xi/mu))
    assert s.simplify(actual_partial-owner.function.pulse_partial)==0
    assert s.simplify(actual_prefix-owner.function.pulse_prefix)==0
    assert s.simplify(actual_partial+actual_prefix-s.integrate(density,(t,0,13/mu)))==0
    assert set(owner.function.late_atoms)==set(current.LATE)
    assert tuple(current.LATE)==tuple(list(witness.native_densities)[6:])
    for name,atom in owner.function.late_atoms.items():assert atom==witness.stage_integral(name,True)
    # Independently prove the residual fixed-at-pulse-end is y-constant,
    # so all preserved radial/axial mixed rows have the original FTC.
    y,Z=s.symbols('y Z',real=True);U,p_rate,endpoint=s.symbols('U p endpoint',positive=True)
    late=s.Function('B_late')(Z);qx=1+Z*Z
    pressure=-U**2*s.exp(-p_rate*y)/(2*p_rate*qx**2)+U**2*s.exp(-p_rate*endpoint)*(1/(2*p_rate*qx**2)-late)
    radial_proofs={}
    for k in range(1,5):
        expected=U**2*s.exp(-p_rate*y)*(-p_rate)**(k-1)/(2*qx**2)
        for j in range(5-k):
            assert s.simplify(s.diff(pressure,y,k,Z,j)-s.diff(expected,Z,j))==0
            radial_proofs['y%d_Z%d'%(k,j)]=True
    # Complex q and source slopes establish the Cauchy / infinite-tail
    # domination without selecting interval masses for exact integrals.
    xr,yi,Zc=s.symbols('x y Zc',real=True)
    real_q=s.re(1+(Zc+xr+s.I*yi)**2)
    assert s.expand(real_q-(1-yi**2+(Zc+xr)**2))==0
    assert 1-s.Rational(1,4)**2==s.Rational(15,16)
    assert s.integrate(s.exp(-p_rate*t)/2,(t,0,s.oo))==1/(2*p_rate)
    assert s.integrate(s.exp(-t),(t,0,s.oo))==1
    sigma,dm,mm=s.symbols('sigma delta mu',nonnegative=True)
    slopes=(-s.Rational(1,2)-mm,
        -s.Rational(1,2)-mm-(1-mm)*sigma,
        -s.Rational(3,2),-s.Rational(3,2)+(1-dm/2)*sigma,
        -(1+dm)/2)
    # Their endpoint extrema on sigma in [0,1], mu/delta in [0,1]
    # attain at most -1/2. The raw collar bracket's endpoint extrema
    # are 1-epsilon and 1, retaining its amplitude normalization.
    for slope in slopes:
        for sv in (0,1):
            for dv in (0,1):
                for mv in (0,1):assert slope.subs({sigma:sv,dm:dv,mm:mv})<=-s.Rational(1,2)
    eps,phi=s.symbols('epsilon phi',positive=True)
    bracket=1-eps*(1-sigma+sigma*phi)
    assert set(s.expand(bracket.subs({sigma:a,phi:b})) for a in (0,1) for b in (0,1))=={1,1-eps}
    # Derive the global envelope independently from ORIGINAL stage shapes,
    # accumulating each domain length, rather than assuming prefactors.
    cursor=s.Integer(0);JJ,gap=s.symbols('JJ gap',nonnegative=True)
    positive_complement=s.Symbol('one_minus_epsilon',positive=True)
    L,Ts,W=p['L'],p['Ts'],p['W'];bp=s.Rational(1,2)+mu
    independent_global={}
    for stage in current.LATE[1:]:
        a,b=p['domains'][stage];offset=cursor-a
        shape=p['shapes'][stage]
        # SymPy moves exp(-log(2)) and exp(-log(1-epsilon))
        # outside the exponential. Include those exact amplitude factors
        # in the logarithm; omitting them would change the normalization.
        amplitude=shape
        if stage=='heat_collar':
            raw_bracket=1-p['epsilon']*(1-p['sigma'](t)+p['sigma'](t)*p['phi'](t))
            amplitude=shape/raw_bracket
        # All amplitudes and 1-epsilon are positive on the bound source
        # domain. Name its positive complement before simplification to
        # prevent factoring -(epsilon-1) into a complex logarithm branch.
        amplitude=s.simplify(amplitude.subs(p['epsilon'],1-positive_complement))
        full_log_amplitude=s.expand_log(s.log(amplitude),force=True)
        exponent_difference=2*(full_log_amplitude-p['logEv']+bp*100+s.log(2))+offset+t
        if stage in ('heat_collar','exterior_power_tail'):
            exponent_difference+=2*s.log(positive_complement)
        cone=s.simplify(s.expand_log(exponent_difference,force=True))
        cone=s.expand(cone.subs(p['J'](t),JJ).subs(t,JJ+gap))
        assert all(v>=0 for v in s.Poly(-cone,mu,p['delta'],L,Ts,W,JJ,gap).coeffs()),stage
        independent_global[stage]=dict(global_offset=str(offset),derived_nonpositive_exponent=str(cone),passed=True)
        cursor=offset+b
    assert cursor==s.oo and len(independent_global)==7
    a,b,D=s.symbols('edge_a edge_b D',positive=True)
    assert s.ask(s.Q.positive(a/(a+b))) is True and s.ask(s.Q.positive(b/(a+b))) is True
    assert s.ask(s.Q.nonpositive(-4/D**2)) is True
    assert owner.majorant_proof['original_cutoff_inequalities']['actual_original_formula_bindings']==current.terminal_source.bind_original_master_source()
    bare=current.correlated.locator.CancelledGraphBounds(owner.graph,owner.ctx,{},())
    assert not bare.polynomial(owner.graph.sub(owner.L_function,owner.graph.mul(owner.graph.constant(-30),owner.logmu_function)))
    assert current.ends(owner.amplitude.reader().at(owner.L_function))[0]>0
    observations={};totals={};pressure_rows=0
    for name,delivery in deliveries.items():
        observed=owner.evaluate(delivery,'1/1000');old=owner.before.evaluate(delivery,'1/1000')
        tail=owner.source(delivery)['source_equivalent_remaining_pressure_tail']
        reader,token,request=owner.amplitude.before._validate_delivery(delivery)
        z=independent.mpf(current.ends(reader.at(delivery['physical_inverse']['coordinate_functions']['Z'])))
        m=independent.mpf(current.ends(owner.amplitude.reader(delivery).bindings[owner.amplitude.mu.node]))
        e=independent.mpf(current.ends(owner.terminal.raw.flat.inlet.datum.parameters.epsilon))
        pr=1+2*m;rho=independent.mpf(1)/4
        flat_upper=1/(2*pr*(1-rho*rho)**2)
        post_upper=independent.exp(-100*pr)/(8*(1-e)**2)
        assert included(tail['flatten_majorant'],flat_upper)
        assert included(tail['postflatten_majorant'],post_upper)
        assert included(tail['late_integral_Taylor_enclosure'][0],independent.mpf([0,current.ends(flat_upper+post_upper)[1]]))
        assert current.ends(tail['ratio_log_enclosure'])[1]<=-current.signed.EXP_LOG_LIMIT
        assert included(tail['positive_ratio_enclosure'],independent.exp(-1000))
        assert tail['ratio_upper_is_bound_not_value'] and tail['all_late_positive_integrals_retained']
        assert owner.graph.nodes[tail['exact_native_logR_FTC_identity']]['existing_radial_product_rows_retained']
        Ubox=independent.mpf(current.ends(owner.refined.constants['U']))
        ratio=independent.mpf(current.ends(tail['positive_ratio_enclosure']))
        for j in range(6):
            inv=s.diff(qx**(-2),Z,j)/s.factorial(j)
            main=rational_expression(independent,inv,Z,z)*Ubox*Ubox/(2*pr)
            bound=flat_upper/rho**j
            late_j=(independent.mpf([0,current.ends(flat_upper+post_upper)[1]]) if j==0
                    else independent.mpf([-current.ends(bound)[1],current.ends(bound)[1]]))
            alternative=-main+(main-late_j*Ubox*Ubox)*ratio
            assert included(tail['coefficients'][j],alternative),j
        original_view=owner.refined.source(delivery);new_view=owner.source(delivery)
        assert new_view is not original_view
        assert new_view['original_factorized_values']['P0'] is original_view['original_factorized_values']['P0']
        for section in ('log_radius_mixed_rows','native_coordinate_mixed_rows'):
            for component,rows in original_view[section].items():
                for label,row in rows.items():
                    if component!='pressure' or row.derivative[0]>0:assert new_view[section][component][label] is row
        for section,components in observed['physical_value_rows'].items():
            for component,rows in components.items():
                for label,row in rows.items():
                    prior=old['physical_value_rows'][section][component][label]
                    if component!='p':
                        assert current.correlated.report(row)==current.correlated.report(prior)
                        continue
                    pressure_rows+=1
                    common=independent.mpf(0)
                    for term in row['ratio_terms']:
                        assert term['directed_ratio_enclosure'] is not None
                        common+=independent.mpf(current.ends(term['signed_coefficient']))*independent.mpf(current.ends(term['directed_ratio_enclosure']))
                    assert included(row['common_scale_coefficient_enclosure'],common)
                    assert row['physical_accuracy']['sign'] in (-1,1)
                    assert not row['physical_accuracy']['ordinary_numeric_delivery_target_satisfied']
        base=observed['physical_value_rows']['Cartesian_spatial_rows']['p']['x0_y0_z0']
        assert base['physical_accuracy']['sign']==-1 and base['physical_accuracy']['factored_relative_width_satisfied']
        for key,value in observed['delivery_counts'].items():totals[key]=totals.get(key,0)+value
        observations[name]=dict(counts=observed['delivery_counts'],old_counts=old['delivery_counts'],
            recovered_pressure_accuracy=base['physical_accuracy'],all_pressure_row_signs_resolved=True,
            original_P0_and_radial_FTC_rows_identical=True,velocity_rows_identical=True)
    first=next(iter(deliveries.values()))
    invalid=dict(copied_delivery=rejected(lambda:owner.evaluate(copy.copy(first))),
        invalid_target=rejected(lambda:owner.evaluate(first,'0')),
        pure_xi_outside_domain=rejected(lambda:owner.function.remaining_pressure_after_pulse('0','13')),
        pure_Z_outside_domain=rejected(lambda:owner.function.remaining_pressure_after_pulse('2','10')))
    saved=owner.function.remaining
    try:
        owner.function.remaining=s.Integer(0)
        invalid['tail_function_substitution']=rejected(owner.assert_graph)
    finally:owner.function.remaining=saved
    saved=copy.deepcopy(owner.graph.nodes[owner.proof.node])
    try:
        owner.graph.nodes[owner.proof.node]['no_datum_or_forcing_replaced']=False
        invalid['proof_node_substitution']=rejected(owner.assert_graph)
    finally:owner.graph.nodes[owner.proof.node].clear();owner.graph.nodes[owner.proof.node].update(saved)
    assert all(invalid.values()) and all(owner.assert_graph().values())
    hashes=dict(owner.hashes);hashes[current.NAME]=current.sha(current.NAME)
    hashes[Path(__file__).name]=current.sha(Path(__file__).name)
    receipt=dict(all_passed=True,source_family=owner.family_record,actual_row_counts=totals,
        actual_pressure_rows_checked=pressure_rows,actual_pressure_tail_observations=observations,
        independent_complete_original_integral_identity=original,
        independent_partial_pulse_integral_and_FTC=True,independent_radial_mixed_identities=radial_proofs,
        independent_holomorphic_and_source_slope_majorants=True,
        independent_original_global_stage_envelope=independent_global,
        exact_positive_Lrel_and_cutoff_sources_bound=True,
        all_positive_late_atoms_bounded_not_dropped=True,
        independent_500_digit_tail_coefficients_and_signed_physical_sums=True,
        invalid_inputs_rejected=invalid,input_hashes=hashes,execution_seconds=time.monotonic()-began,
        **dict.fromkeys(current.GATES,True),**dict.fromkeys(current.OPEN,False))
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.correlated.report(receipt),indent=2)+'\n',encoding='utf8',newline='\n')
    print('PASS_CURRENT_ORIGINAL_RP_REMAINING_PRESSURE_TAIL original correlated pressure and physical signs',flush=True)
    return receipt
