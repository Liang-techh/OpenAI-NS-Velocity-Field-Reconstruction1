"""Focused source-weighted cone identities, transfer and all-N checks."""
import json
from pathlib import Path
import mpmath as mp
import sympy as s
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_O3_modified_transition_cone as source
from lei_ren_part1_paper_compliant_current_original_cone_operator import cone_margins


def independent_signed_weighted_fixtures():
    """Compare with the original signed D/J/Q implementation at moderate mu."""
    c=MPIntervalContext();c.dps=90;count=positive=negative=0
    with mp.workdps(120):
        for muv in ('.001','.03'):
            mu=c.mpf(muv)
            for gv,hv in ((-1,-4),(-1,4),(0,-1),(0,1),(2,0),(2,-4),(2,4),(3,1)):
                g,h=c.mpf(gv),c.mpf(hv);a=2+mu*g;bs=c.sqrt(mu)*h
                for tv,zv in (('1','-.01'),('1','.01'),('1','10'),('-1','.01')):
                    theta,zeta=c.mpf(tv),c.mpf(zv);tz=zeta/c.sqrt(mu)
                    tested=cone_margins(c,a,bs,theta,tz,mu*(2*g+h*h+mu*g*g)/a)
                    D=theta-bs*tz/a;J=tz+bs*theta/a;Q=2*D*D-(a-2+bs*bs/a)*J*J
                    H=a*theta-h*zeta;M=2*g+h*h+mu*g*g
                    W=(2*a-mu*h*h)*theta*theta-2*a*h*theta*zeta-a*g*zeta*zeta
                    direct=dict(direction=H,quadratic=(a*a+mu*h*h)*W)
                    original=dict(direction=a*D,quadratic=a**3*Q)
                    for key,value in direct.items():
                        lo,hi=source.endpoints(value);u,v=source.endpoints(original[key])
                        if hi<u or v<lo:raise ArithmeticError('Weighted signed cone differs from original '+key)
                    actual=(source.endpoints(a)[0]>0 and source.endpoints(M)[0]>0 and
                            source.endpoints(H)[0]>0 and source.endpoints(W)[0]>0)
                    if actual!=tested['admitted']:raise ArithmeticError('Weighted cone changed original signed admission')
                    positive+=actual;negative+=not actual;count+=1
    if not positive or not negative:raise ValueError('Missing signed success/failure fixtures')
    return count,positive,negative


def independent_source_scale_identities(field):
    R,P,mu,Ad=s.symbols('R Pstar mu Ad',positive=True);count=0
    common=P*s.sqrt(R/2)
    for label,parts in field.errors.model['stress'].items():
        for p in parts.values():
            rate=s.Rational(str(p['mode'][0]));power=p['mode'][1];aa=p['error_Ad_power']
            actual=R**rate*P**power*Ad**aa/s.sqrt(2)/common
            expected=R**(rate-s.Rational(1,2))*P**(power-1)*Ad**aa
            if label=='axial':actual*=s.sqrt(mu);expected*=s.sqrt(mu)
            if s.simplify(actual-expected)!=0:raise ArithmeticError('Relative weighted tensor factor differs')
            count+=1
    return count


def run():
    data=json.loads((source.HERE/source.NAME).read_bytes())
    for name,digest in data['input_hashes'].items():
        if source.sha(name)!=digest:raise ValueError('Changed closed modified O3 source: '+name)
    field=source.CurrentModifiedO3TransitionCone(require_checked=False)
    if source.parameters.encoded(field.theorem)!=data['exact_weighted_source_cone_theorem']:
        raise ValueError('Weighted actual source theorem changed')
    if source.parameters.encoded(field.proof)!=data['whole_modified_closed_O3_cone']:
        raise ValueError('Whole current cone inequality ledger changed')
    if source.parameters.encoded(field.error_rows)!=data['complete_leading_stress_error_log_ledger']:
        raise ValueError('Signed source error sector log ledger changed')
    c=field.ctx
    for name,margin in field.proof['positive_margins'].items():
        if source.endpoints(margin)[0]<=0:raise ArithmeticError('Unresolved whole current margin: '+name)
    if source.endpoints(field.proof['weighted_H_over_theta_lower'])[0]<=source.endpoints(c.mpf('1.9'))[0]:
        raise ArithmeticError('Signed alignment reserve below1.9')
    if source.endpoints(field.proof['original_signed_D_over_theta_lower'])[0]<=source.endpoints(c.mpf('.9'))[0]:
        raise ArithmeticError('Original signed direction reserve below.9')
    if source.endpoints(field.proof['original_signed_Q_over_theta_squared_lower'])[0]<=source.endpoints(c.mpf('1.9'))[0]:
        raise ArithmeticError('Original signed quadratic reserve below1.9')
    if source.endpoints(field.proof['reduced_weighted_quadratic_over_theta_squared_lower'])[0]<=source.endpoints(c.mpf('3.9'))[0]:
        raise ArithmeticError('Weighted quadratic reserve below3.9')
    native_count=0
    for N in (37,10**12):
        rows=field.leading_error_rows(N)
        for label,parts in rows.items():
            for name,row in parts.items():
                ref=field.error_rows[label][name]
                if source.endpoints(row['native_coefficient_absolute_cap'])[1]>source.endpoints(ref['native_coefficient_absolute_cap'])[1]:
                    raise ArithmeticError('Larger N increased source native0 cap: '+name)
                native_count+=1
    specs=(('closed_whole',(0,1),22),('old_zero_shear_seam','0',22),
        ('right_flat_taper',('.25','.5'),22),('right_flat_edge','.5',22),
        ('terminal','1',22),('current_frequency',('0','.5'),10**12))
    for name,t,N in specs:
        if source.parameters.encoded(field.query(t,N))!=data['examples'][name]:
            raise ValueError('Changed source-bound cone query: '+name)
    invalid=((-1,22),(('-.1','.1'),22),(2,22),(0,21),(0,True),(0,1.5),(mp.inf,22))
    for t,N in invalid:
        try:field.query(t,N)
        except (ValueError,TypeError):pass
        else:raise ArithmeticError('Invalid current O3 cone domain admitted')
    # The left flat O2 edge is outside this certificate. Its exact primitive
    # excess remains zero, while the closed O3 seam is now strictly positive.
    left=field.primitive.query('-2',22)
    if left['source_certified_joint_margin_lower_envelope']._mpi_!=c.mpf(0)._mpi_:
        raise ArithmeticError('Left O2 flat edge relabeled as strict')
    if source.endpoints(field.query('0',22)['source_correlated_primitive_lower'])[0]<=0:
        raise ArithmeticError('Old O3 zero-shear seam not repaired')
    for flag in source.OPEN:
        if data.get(flag) is not False:raise ValueError('Scoped O3 proof promoted global or other regional gate')
    fixtures,positive,negative=independent_signed_weighted_fixtures()
    scales=independent_source_scale_identities(field)
    result=dict(all_passed=True,source_family=data['source_family'],
        exact_source_weighted_cone_and_all_N_identities=len(field.theorem['identities']),
        strict_whole_current_source_inequalities=len(field.proof['positive_margins']),
        independent_original_signed_cone_fixtures=fixtures,signed_admitted_fixtures=positive,signed_rejected_fixtures=negative,
        independent_current_Pstar_radius_Ad_weighted_scale_identities=scales,
        larger_N_native_error_comparisons=native_count,current_queries_recomputed=len(specs),invalid_queries_rejected=len(invalid),
        all_N_certificate_uses_exact_monotonicity_not_the_comparison_samples=True,
        all_phases_closed_offset_and_axial_domain_proved_by_source_inequalities=True,
        all_original_pressure_energy_moments_and_cross_sectors_retained=True,
        exact_old_O3_seam_now_strict_left_O2_flat_edge_unadmitted=True,
        current_modified_closed_O3_signed_two_vector_cone_certified=True,
        **{flag:False for flag in source.OPEN},
        input_hashes={**data['input_hashes'],source.NAME:source.sha(source.NAME),Path(__file__).name:source.sha(Path(__file__).name)})
    (source.HERE/source.RECEIPT).write_bytes((json.dumps(result,indent=2)+'\n').encode())
    print('Modified closed O3 cone PASS: all N>=22;',fixtures,'signed fixtures;',len(field.proof['positive_margins']),'strict inequalities',flush=True)
    return result


if __name__=='__main__':run()
