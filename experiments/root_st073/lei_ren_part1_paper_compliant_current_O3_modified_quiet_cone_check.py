"""Focused quiet signed-source, all-N error transfer and join-domain checks."""
import json
from pathlib import Path
import mpmath as mp
import lei_ren_part1_paper_compliant_current_O3_modified_quiet_cone as source
from lei_ren_part1_paper_compliant_current_O3_modified_transition_cone_check import (
    independent_signed_weighted_fixtures,independent_source_scale_identities)


def independent_signed_quiet_shear_fixtures(field):
    """Actual derivative quotient, with positive default axial shear sign.

    Synthetic control vectors test the uniform ball. They are never used
    as the actual current h_N or as another N's implicit solution.
    """
    counts=positive=negative=0
    with mp.workdps(280):
        def point(v):
            lo,hi=source.endpoints(v);return (lo+hi)/2
        R=point(field.repair.R);norm=point(field.repair.normalization)
        ell=mp.mpf(1)/40;ci=[mp.mpf(1)/5,mp.mpf(1)/2,mp.mpf(4)/5]
        G=[source.endpoints(v)[1] for v in field.repair.jets['log_bump_ordinary_derivative_caps']]
        def bump(y,i):
            r=(y-ci[i])/ell
            return mp.exp(-y-1/(1-r*r))/(ell*norm) if abs(r)<1 else mp.mpf(0)
        for mu in (mp.mpf('.001'),mp.mpf('1e-20')):
            alpha=mp.mpf('.5')+mu
            for N in (field.minimum_N,10**12):
                denmin=mp.exp(-alpha)-mu*R*G[0]/N
                if denmin<=0:raise ArithmeticError('Independent fixture denominator nonpositive')
                gerr=2*R*(G[1]+alpha*G[0])/(N*denmin);hcap=2*R*G[1]/(N*denmin)
                for center in ci:
                    for fraction in ('-.75','0','.73'):
                        y=center+ell*mp.mpf(fraction)
                        for values in ((1,-1,1,-1,1),(-1,1,-1,1,-1),('.3','-.7','.8','-.2','.4')):
                            controls=[mp.mpf(v)*R for v in values]
                            def Hs(v):return sum(controls[i+2]*bump(v,i) for i in range(3))
                            def Ha(v):return controls[0]*bump(v,0)+controls[1]*bump(v,2)
                            def U(v):return mp.exp(-alpha*v)+mu*Hs(v)/N
                            def V(v):return mp.sqrt(mu)*Ha(v)/N
                            f=mp.exp(-alpha*y);den=U(y);ey=mp.diff(Hs,y);ay=mp.diff(Ha,y)
                            g=(2*f-(Hs(y)+2*ey)/N)/den;h=2*ay/(N*den)
                            # At these moderate fixtures subtraction is
                            # resolved far beyond both interval caps.
                            actual_a=1-2*mp.diff(U,y)/U(y)
                            actual_b=2*mp.diff(V,y)/U(y)
                            tol=mp.mpf(2)**(-800)
                            if abs((actual_a-2)/mu-g)>tol or abs(actual_b/mp.sqrt(mu)-h)>tol:
                                raise ArithmeticError('Actual positive signed quiet quotient differs')
                            if abs(g-2)>gerr or abs(h)>hcap:raise ArithmeticError('Quiet source quotient exceeds uniform cap')
                            positive+=h>0;negative+=h<0;counts+=1
        if not positive or not negative:raise ValueError('Both axial shear signs must be checked')
    return counts,positive,negative


def run():
    data=json.loads((source.HERE/source.NAME).read_bytes())
    for name,digest in data['input_hashes'].items():
        if source.sha(name)!=digest:raise ValueError('Changed complete quiet cone source: '+name)
    field=source.CurrentModifiedQuietCone(require_checked=False);c=field.ctx
    for key,value in (('exact_signed_quiet_source_and_weighted_cone_theorem',field.theorem),
        ('original_power_local_baseline',field.baseline),('whole_modified_quiet_cone',field.proof),
        ('complete_leading_stress_error_log_ledger',field.error_rows)):
        if source.parameters.encoded(value)!=data[key]:raise ValueError('Quiet source proof/ledger changed: '+key)
    for name,value in field.proof['positive_margins'].items():
        if source.endpoints(value)[0]<=0:raise ArithmeticError('Quiet signed cone inequality unresolved: '+name)
    for key,minimum in (('weighted_H_over_theta_lower','1.9'),('original_signed_D_over_theta_lower','.9'),
        ('original_signed_Q_over_theta_squared_lower','1.9'),('reduced_weighted_quadratic_over_theta_squared_lower','3.9')):
        if source.endpoints(field.proof[key])[0]<=source.endpoints(c.mpf(minimum))[0]:raise ArithmeticError('Quiet signed cone reserve: '+key)
    native_count=0
    for N in (2*field.minimum_N,10**12):
        for label,parts in field.leading_error_rows(N).items():
            for name,row in parts.items():
                if source.endpoints(row['native_coefficient_absolute_cap'])[1]>source.endpoints(field.error_rows[label][name]['native_coefficient_absolute_cap'])[1]:
                    raise ArithmeticError('Quiet N increase enlarged native0 cap')
                native_count+=1
    specs=(('whole_strip',(0,1),field.minimum_N),('entry','0',field.minimum_N),
        ('first_bump',('.175','.225'),field.minimum_N),('first_gap',('.225','.475'),field.minimum_N),
        ('middle_bump',('.475','.525'),field.minimum_N),('last_bump',('.775','.825'),field.minimum_N),
        ('near_last_flat_edge',('.824','.825'),field.minimum_N),('final_exact_strip',('.826','1'),field.minimum_N),
        ('terminal','1',field.minimum_N),('current_frequency',(0,1),10**12))
    for name,y,N in specs:
        q=field.query(y,N)
        if source.parameters.encoded(q)!=data['examples'][name]:raise ValueError('Quiet complete cone query changed: '+name)
        if not q['complete_modified_signed_cone_certified_for_entire_query_box']:raise ArithmeticError('Quiet entry/strip/join missing')
    if not data['examples']['terminal']['exact_full_source_error_zero_after_last_support'] or not data['examples']['final_exact_strip']['exact_full_source_error_zero_after_last_support']:
        raise ArithmeticError('Same implicit terminal source equality missing')
    if data['examples']['entry']['exact_full_source_error_zero_after_last_support'] or data['examples']['first_gap']['exact_full_source_error_zero_after_last_support']:
        raise ArithmeticError('Local flat gap incorrectly erased cumulative history')
    zeros=0
    q,caps=field.errors.inputs('quiet',('.826','1'),field.minimum_N)
    for label,parts in field.errors.model['stress'].items():
        for name,p in parts.items():
            value=source.errors.bound_expression(c,p['expression'],caps,c.mpf([-1,1]),field.errors.delta,0,p['mode'][0])[0,0]
            if source.endpoints(value)!=(0,0):raise ArithmeticError('Post-support original tensor identity lost: '+name)
            zeros+=1
    invalid=((-1,field.minimum_N),(2,field.minimum_N),(('.99','1.01'),field.minimum_N),
        (0,field.minimum_N-1),(0,True),(0,1.5),(mp.inf,field.minimum_N))
    for y,N in invalid:
        try:field.query(y,N)
        except (ValueError,TypeError):pass
        else:raise ArithmeticError('Invalid quiet cone query admitted')
    for key in source.OPEN:
        if data.get(key) is not False:raise ValueError('Quiet scoped source cone promoted global gate')
    signs,pos,neg=independent_signed_quiet_shear_fixtures(field)
    fixtures,admitted,rejected=independent_signed_weighted_fixtures();scales=independent_source_scale_identities(field)
    result=dict(all_passed=True,source_family=data['source_family'],minimum_integer_N=field.minimum_N,
        exact_signed_quiet_source_coordinate_monotonicity_and_weighted_identities=len(field.theorem['identities']),
        strict_whole_quiet_source_inequalities=len(field.proof['positive_margins']),
        independent_actual_positive_signed_quiet_shear_fixtures=signs,positive_axial_shear_fixtures=pos,negative_axial_shear_fixtures=neg,
        independent_original_signed_cone_fixtures=fixtures,signed_admitted_fixtures=admitted,signed_rejected_fixtures=rejected,
        independent_current_Pstar_radius_Ad_weighted_scale_identities=scales,
        exact_post_support_complete_tensor_error_zeros=zeros,larger_N_native_error_comparisons=native_count,
        current_bump_gap_entry_terminal_queries_recomputed=len(specs),invalid_queries_rejected=len(invalid),
        local_power_baseline_restriction_and_all_N_source_monotonicity_not_samples=True,
        same_source_cumulative_gap_history_and_pressure_retained=True,
        current_modified_quiet_signed_cone_certified=True,**{key:False for key in source.OPEN},
        input_hashes={**data['input_hashes'],source.NAME:source.sha(source.NAME),Path(__file__).name:source.sha(Path(__file__).name),
            Path(independent_signed_weighted_fixtures.__code__.co_filename).name:source.sha(Path(independent_signed_weighted_fixtures.__code__.co_filename).name)})
    (source.HERE/source.RECEIPT).write_bytes((json.dumps(source.parameters.encoded(result),indent=2)+'\n').encode())
    print('Modified quiet cone PASS: all N>=',field.minimum_N,';',signs,'positive-sign source fixtures;',zeros,'exact post-support tensor zeros',flush=True)
    return result


if __name__=='__main__':run()
