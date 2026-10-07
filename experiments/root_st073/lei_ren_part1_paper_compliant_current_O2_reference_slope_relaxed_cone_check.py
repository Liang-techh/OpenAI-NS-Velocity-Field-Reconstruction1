"""Focused stronger relaxed admission, full theta units and source seams."""
import json
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_O2_reference_slope_relaxed_cone as source
from lei_ren_part1_paper_compliant_current_original_cone_operator import cone_margins
from lei_ren_part1_paper_compliant_current_pre_pulse_stress_operator import raw_pre_stress_rows
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_compliant_outer_initial import stable_sigma


def independent_paper_relaxed_fixtures():
    c=MPIntervalContext();c.dps=100;count=admitted=rejected=0
    with mp.workdps(130):
        for av in ('.8','1','2','3'):
            for bv in ('0','.7','-2'):
                for tv in ('-.7','.5','3','5'):
                    for zv in ('-3','.2','4'):
                        a,b,theta,axial=map(mp.mpf,(av,bv,tv,zv));k=a+b*b/a
                        # Independent paper dot/cross form with F=1.
                        dot=-a*theta+b*axial;cross=-b*theta-a*axial
                        second=2*dot*dot-(k-2)*cross*cross if k>2 else -dot-a*(2-k)
                        if min(abs(dot),abs(second))<mp.mpf('1e-80'):continue
                        expected=dot<0 and second>0
                        actual=source.relaxed_cone_margins(c,av,bv,tv,zv)
                        if actual['admitted']!=expected:raise ArithmeticError('Signed paper relaxed branch differs')
                        admitted+=expected;rejected+=not expected;count+=1
    if not admitted or not rejected:raise ArithmeticError('Independent relaxed fixtures need both decisions')
    if source.relaxed_cone_margins(c,1,0,'.5',0)['admitted']:
        raise ArithmeticError('Positive theta alone incorrectly admitted kappa<2')
    good=source.relaxed_cone_margins(c,('1.9','2.1'),0,10,1)
    bad=source.relaxed_cone_margins(c,('1.9','2.1'),0,1,100)
    if not good['admitted'] or good['strict_source_shear_for_entire_box'] or bad['admitted']:
        raise ArithmeticError('Crossing kappa2 box branch incorrectly promoted')
    return count,admitted,rejected


def independent_full_reference_slope_fixtures():
    c=MPIntervalContext();c.dps=100;count=0
    with mp.workdps(130):
        for chart,yv in (('Rh_reference','-5'),('Rh_reference','-2.337'),('Rh_reference','0'),
            ('O2_slope','0'),('O2_slope','.537'),('O2_slope','1')):
            y=c.mpf(yv)
            if chart=='Rh_reference':U=c.exp(y/10);X=c.mpf(5)/8;lam=c.mpf(1)/10;mass=None
            else:
                J,mass=source.slope_masses(c,y,256);U=c.exp(y/10-c.mpf(3)/5*J)
                X=(c.mpf(5)/8+mass[0])*c.exp(-c.mpf(3)/2*y)/U
                lam=c.mpf(1)/10-c.mpf(3)/5*stable_sigma(c,y)[0]
            for zv in ('-1','-.373','0','.51','1'):
                constant=lambda v:IntervalTaylor(c,[v]+[c.mpf(0)]*5)
                z=IntervalTaylor(c,[c.mpf(zv),c.mpf(1)]+[c.mpf(0)]*4);zero=z*0
                C=(1+z*z).reciprocal();delta=constant(c.mpf('.001'));L=1-delta*z*z
                A=(1-delta/2)*C+(1-delta)*z*z*C*C
                B=((2*delta*z*z-1)*C+2*(1-z*z)*z*z*C*C)/L
                uu=constant(U)*C;h=uu*X;Pstar=c.mpf(10000);R=c.mpf('1e6')*c.exp(y)
                energy=16*z*z/Pstar**2-(uu*uu*c.mpf(5)/12 if mass is None else C*C*(c.mpf(5)/12+mass[2]/2)*c.exp(-y))
                pressure=(1+z*z)*c.mpf('.017')+(uu*uu*c.mpf('2.5') if mass is None else C*C*(c.mpf('2.5')+mass[1]/2))
                rows=lambda v:[v]+[zero]*4
                u=[uu,uu*lam]+[zero]*3
                raw=raw_pre_stress_rows(c,delta,z,u,rows(4*z),dict(m=rows(4*z),h=rows(h),k=rows(4*z*h),e=rows(energy)),rows(pressure))
                def lift(parts):return sum(R**c.mpf(str(p['mode'][0]))*Pstar**p['mode'][1]*p['shape'][0][0]/c.sqrt(2) for p in parts.values())
                theta,axial=lift(raw['theta']),lift(raw['axial']);F=Pstar*uu[0]/c.sqrt(2*R);a=1-2*lam
                expected=((A*X-C)/L+4*C+4*B*X)[0]*R/C[0]-a
                observed=theta/F
                ol,oh=source.endpoints(observed);el,eh=source.endpoints(expected)
                if oh<el or eh<ol:raise ArithmeticError('Full theta/F source normalization differs')
                relaxed=source.relaxed_cone_margins(c,a,0,observed,axial/F)
                if not relaxed['admitted'] or cone_margins(c,a,0,theta,axial)['admitted']:
                    raise ArithmeticError('Source relaxed/strict scope differs on reference/slope')
                if source.endpoints(observed-(2-a))[0]<=0:
                    raise ArithmeticError('Stronger kappa<=2 direction was omitted')
                count+=1
    return count


def run():
    data=json.loads((source.HERE/source.NAME).read_bytes())
    for name,digest in data['input_hashes'].items():
        if source.sha(name)!=digest:raise ValueError('Changed stronger relaxed source: '+name)
    field=source.CurrentOriginalReferenceSlopeRelaxedCone(require_checked=False);c=field.ctx
    for key,value in (('exact_reference_slope_source_relaxed_cone_theorem',field.theorem),
        ('whole_original_reference_slope_relaxed_cone',field.proof)):
        if source.parameters.encoded(value)!=data[key]:raise ValueError('Original slope proof changed: '+key)
    for name,value in field.proof['positive_margins'].items():
        if source.endpoints(value)[0]<=0:raise ArithmeticError('Stronger relaxed continuous bound missing: '+name)
    specs=(('whole_reference','Rh_reference',(-5,0)),('Rh','Rh_reference',-5),('reference_slope_left','Rh_reference',0),
        ('whole_slope','O2_slope',(0,1)),('reference_slope_right','O2_slope',0),('slope_interior','O2_slope','.537'),('slope_axial_seam','O2_slope',1))
    for name,chart,y in specs:
        q=field.query(chart,y)
        if source.parameters.encoded(q)!=data['examples'][name]:raise ValueError('Original relaxed query changed: '+name)
        if q['complete_modified_signed_cone_certified_for_entire_query_box']:raise ValueError('Original relaxed source admitted strictly')
        if chart=='O2_slope':
            lo,hi=source.endpoints(q['actual_original_h_over_U_enclosure'])
            if lo<source.endpoints(c.mpf(5)/8)[0] or hi>source.endpoints(1-field.proof['source_angular_deficit_lower'])[1]:
                raise ArithmeticError('Source correlated angular deficit enclosure missing')
    left=field.query('Rh_reference',0);right=field.query('O2_slope',0)
    for key in ('native_logR','actual_original_scalar_U_enclosure','actual_original_h_over_U_enclosure'):
        if source.parameters.encoded(left[key])!=source.parameters.encoded(right[key]):raise ArithmeticError('Actual reference/slope source seam differs')
    invalid=(('bad',0),('Rh_reference',-6),('Rh_reference',1),('O2_slope',-1),('O2_slope',2),('O2_slope',('.99','1.01')),('O2_slope',mp.inf))
    for chart,y in invalid:
        try:field.query(chart,y)
        except (ValueError,TypeError):pass
        else:raise ArithmeticError('Invalid original reference/slope box admitted')
    for key in source.OPEN:
        if data.get(key) is not False:raise ValueError('Original relaxed input promoted global gate')
    fixtures,admitted,rejected=independent_paper_relaxed_fixtures();full=independent_full_reference_slope_fixtures()
    result=dict(all_passed=True,source_family=data['source_family'],
        exact_original_source_deficit_geometry_and_relaxed_unit_identities=len(field.theorem['identities']),
        strict_whole_reference_slope_inequalities=len(field.proof['positive_margins']),
        independent_paper_signed_relaxed_branch_fixtures=fixtures,signed_admitted_fixtures=admitted,signed_rejected_fixtures=rejected,
        independent_full_reference_slope_stress_unit_fixtures=full,
        source_queries_recomputed=len(specs),source_seam_scalar_and_radius_equalities=3,invalid_queries_rejected=len(invalid),
        stronger_kappa_le2_lower_not_replaced_by_direction_only=True,
        whole_source_deficit_and_radius_bounds_not_sampled_admission=True,
        current_original_Rh_O2_slope_relaxed_input_cone_certified=True,
        current_original_Rh_O2_slope_strict_cone_certified=False,**{key:False for key in source.OPEN},
        input_hashes={**data['input_hashes'],source.NAME:source.sha(source.NAME),Path(__file__).name:source.sha(Path(__file__).name)})
    (source.HERE/source.RECEIPT).write_bytes((json.dumps(source.parameters.encoded(result),indent=2)+'\n').encode())
    print('Original reference/slope stronger relaxed cone PASS:',full,'full-stress cases;',fixtures,'signed cone decisions',flush=True)
    return result


if __name__=='__main__':run()
