"""Independent original cone equivalence and actual signed product bounds."""
import copy
from fractions import Fraction
import gzip
import json
from pathlib import Path
import time

from mpmath.ctx_iv import MPIntervalContext
import sympy as s
import lei_ren_part1_paper_compliant_current_original_Rp_stress_cone as current
import lei_ren_part1_paper_compliant_current_original_Rp_background_stress_check as stress_check
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

bg=current.background


def independent_cone_identities(owner):
    # Derive the paper's cone coordinates before removing positive factors.
    a,U,B,Tt,Tz=s.symbols('a U B Tt Tz',positive=True)
    tilt=-B/(a*U);speed=a*(1+tilt**2)
    first=Tt+tilt*Tz;perpendicular=Tz-tilt*Tt
    quadratic=2*first**2-(speed-2)*perpendicular**2
    H0=(a-2)*a*U**2+B**2;H1=a*U*Tt-B*Tz
    H2=2*a*U**2*H1**2-H0*(a*U*Tz+B*Tt)**2
    checks=dict(v_minus_2_positive_factor=s.cancel(H0-a*U**2*(speed-2))==0,
        first_margin_positive_factor=s.cancel(H1-a*U*first)==0,
        quadratic_margin_positive_factor=s.cancel(H2-a**3*U**4*quadratic)==0)
    Ut=bg.F['Utheta'];derivative=bg.dy(bg.F['Uz'])
    Tt=owner.before.operators['cylindrical_stress']['r_theta'][0][1]
    Tz=owner.before.operators['cylindrical_stress']['r_z'][0][1]
    rate=2+2*bg.MU
    replacements={a:rate,U:Ut,B:2*derivative}
    # Stress symbols must remain separate until the generic identity is formed.
    stress_a,stress_b=s.symbols('stress_a stress_b')
    generic0=(a-2)*a*U**2+B**2
    generic1=a*U*stress_a-B*stress_b
    generic2=2*a*U**2*generic1**2-generic0*(a*U*stress_b+B*stress_a)**2
    expected=dict(a=rate,Utheta_positive=Ut,H0_v_minus_2=generic0.subs(replacements),
        H1_first_margin=generic1.subs(replacements).subs({stress_a:Tt,stress_b:Tz}),
        H2_quadratic_margin=generic2.subs(replacements).subs({stress_a:Tt,stress_b:Tz}))
    checks['actual_original_cone_source_definitions']=all(s.expand(owner.operators[name][0][1]-expr)==0
        and owner.operators[name][0][0]==0 for name,expr in expected.items())
    checks['common_original_physical_stress_prefactor']=owner.before.operators['cylindrical_stress']['r_theta'][0][0]==\
        owner.before.operators['cylindrical_stress']['r_z'][0][0]==-2-bg.DELTA
    # The actual replayed shear has these normalized cone coordinates.
    checks['original_theta_source_law']=owner.before.before.law['identity']=='Utheta_y=-(1/2+mu)*Utheta'
    assert all(checks.values())
    return checks


@source_precision
def run(owner,values,fields):
    began=time.monotonic();raw=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    assert not owner.acceptance_loaded and not raw[current.GATE] and not raw[current.POINT_CONE]
    assert not any(raw[key] for key in current.OPEN)
    assert raw['source_family']==owner.family_record and raw['original_cone_operator_definitions']==owner.definitions
    assert raw['original_cone_equivalence']==owner.equivalence and all(owner.assert_graph().values())
    assert all(raw['input_hashes'][name]==digest for name,digest in owner.before.hashes.items())
    identities=independent_cone_identities(owner)
    c=MPIntervalContext();c.dps=800;observations={};products_checked=bounds_checked=0
    for name,value in values.items():
        report=owner.report(value);record=owner.differential._require(fields[name])
        source=owner.product.before.source(record['delivery'])
        _,_,request=owner.product.amplitude.before._validate_delivery(record['delivery'])
        issued_reader=owner.product.reader(record['delivery'])
        reader=bg.correlated.CorrelatedGraphBounds(bg.correlated.locator.CancelledGraphBounds(
            owner.graph,c,dict(issued_reader.bindings),issued_reader.allowed_exponentials),dict(issued_reader.aliases))
        angle=c.mpf(request.theta.numerator)/request.theta.denominator
        bindings={bg.Z:c.mpf(request.Z.numerator)/request.Z.denominator,
            bg.DELTA:c.mpf(issued_reader.bindings[owner.before.before.delta.node]),
            bg.MU:c.mpf(issued_reader.bindings[owner.before.before.mu.node]),bg.CS:c.cos(angle),bg.SN:c.sin(angle)}
        assert bg.ends(bindings[bg.MU])[0]>0
        assert bg.ends(source['log_radius_mixed_rows']['Utheta']['y0_Z0'].coefficients[0])[0]>0
        signs={};unresolved=0
        for component,result in report['cone_margins'].items():
            groups={}
            for item in result['actual_source_product_ledger']:
                coefficient=stress_check.numeric(s.sympify(item['exact_operator_expression']),c,bindings)
                units={'R':Fraction(item['radial_power']['numerator'],item['radial_power']['denominator'])}
                product=c.mpf(1)
                for atom in item['actual_shared_source_product']:
                    label=atom['source'];k,n=atom['ordinary_derivative'];power=atom['power']
                    row=source['log_radius_mixed_rows'][label]['y%d_Z%d'%(k,n)]
                    assert atom['actual_source_units']==row.source_units
                    assert atom['actual_source_powers']==[bg.signed.rational_record(q) for q in row.powers]
                    product*=c.mpf(row.coefficients[0])**power
                    for unit,q in zip(row.source_units,row.powers):units[unit]=units.get(unit,Fraction(0))+power*q
                assert bg.axial.differential.product.contains(item['signed_product_coefficient'],coefficient*product)
                log=bg.box.pulse.radius.FunctionRef(owner.graph,item['exact_product_log_scale'])
                key=(tuple(sorted(reader.polynomial(log).items())),tuple(sorted(units.items())),item['lambda_exponent'])
                group=groups.setdefault(key,dict(log=log,coefficient=c.mpf(0)))
                group['coefficient']+=coefficient*product;products_checked+=1
            assert not result['exact_zero_enclosure']
            reference=bg.box.pulse.radius.FunctionRef(owner.graph,result['exact_reference_log_scale_function'])
            total=c.mpf(0);failed=False
            for group in groups.values():
                difference=owner.graph.sub(group['log'],reference)
                if not reader.polynomial(difference):ratio=c.mpf(1)
                else:
                    lo,hi=bg.ends(reader.at(difference))
                    if hi>1000:failed=True;continue
                    if hi<=-1000:ratio=c.mpf([0,bg.ends(c.exp(-1000))[1]])
                    elif lo<-1000:ratio=c.mpf([0,bg.ends(c.exp(c.mpf(hi)))[1]])
                    else:ratio=c.exp(reader.at(difference))
                total+=group['coefficient']*ratio
            if failed:
                assert result['common_scale_coefficient_enclosure'] is None
                signs[component]='ratio_unresolved';unresolved+=1
            else:
                assert bg.axial.differential.product.contains(result['common_scale_coefficient_enclosure'],total),component
                signs[component]=result['signed_log_value']['sign'];bounds_checked+=1
            if component=='a':
                assert result['ordinary_numeric_materialized']
                assert bg.axial.differential.product.contains(result['ordinary_numeric_enclosure'],2+2*bindings[bg.MU])
            else:assert not result['ordinary_numeric_materialized']
            assert result['nonlinear_shared_source_products_collected_before_interval_bounds']
        positive=all(sign==1 for sign in signs.values())
        assert positive==report['signed_margins_certify_original_two_component_cone_at_this_point']
        assert report['completed_diagonal_excluded_from_cone'] and not report['regional_cone_certified']
        assert not report[current.POINT_CONE] and all(not report[key] for key in current.OPEN)
        observations[name]=dict(signs=signs,all_strict_cone_margins_positive=positive,
            actual_original_admissible_cone_at_this_point=positive,unresolved_scale_ratio_outputs=unresolved,
            positive_original_mu_swirl_and_common_prefactor=True,regional_cone_remains_open=True)
    value=next(iter(values.values()));entry=owner._fields[id(value)]
    invalid=dict(copied_cone_packet=stress_check.rejected(lambda:owner.report(copy.copy(value))),
        invalid_target=stress_check.rejected(lambda:owner.evaluate(next(iter(fields.values())),'2')))
    node=next(iter(entry[5]));saved=copy.deepcopy(owner.graph.nodes[node])
    try:
        owner.graph.nodes[node]['changed_actual_cone_source']=True
        invalid['changed_cone_source_DAG']=stress_check.rejected(lambda:owner.report(value))
    finally:owner.graph.nodes[node]=saved
    saved=owner.operators['H2_quadratic_margin']
    try:
        owner.operators['H2_quadratic_margin']=[(0,s.Integer(1))]
        invalid['substituted_quadratic_cone_margin']=stress_check.rejected(owner.assert_graph)
    finally:owner.operators['H2_quadratic_margin']=saved
    assert all(invalid.values()) and all(owner.assert_graph().values())
    hashes=dict(owner.hashes);hashes[current.NAME]=current.sha(current.NAME);hashes[Path(__file__).name]=current.sha(Path(__file__).name)
    receipt=dict(all_passed=True,source_family=owner.family_record,original_cone_operator_definitions=owner.definitions,
        original_cone_equivalence=owner.equivalence,independent_original_cone_identities=identities,
        actual_original_cone_point_observations=observations,independently_checked_shared_source_products=products_checked,
        independently_checked_signed_margin_bounds=bounds_checked,independent_800_digit_source_products=True,
        invalid_inputs_rejected=invalid,input_hashes=hashes,execution_seconds=time.monotonic()-began,
        **{current.GATE:True,current.POINT_CONE:all(v['actual_original_admissible_cone_at_this_point'] for v in observations.values())},
        **dict.fromkeys(current.OPEN,False))
    (current.HERE/current.RECEIPT).write_text(json.dumps(bg.correlated.report(receipt),indent=2)+'\n',encoding='utf8',newline='\n')
    print('PASS_CURRENT_ORIGINAL_RP_STRESS_CONE independent point margins',flush=True)
    return receipt
