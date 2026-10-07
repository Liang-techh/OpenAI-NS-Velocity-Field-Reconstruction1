"""Check full negative buffer input and prevent strict/global promotion."""
import json
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_O2_relaxed_buffer_cone as source
from lei_ren_part1_paper_compliant_current_original_cone_operator import cone_margins
from lei_ren_part1_paper_compliant_current_pre_pulse_stress_operator import raw_pre_stress_rows
from lei_ren_part1_paper_interval_taylor import IntervalTaylor


def independent_full_stress_degenerate_fixtures():
    c=MPIntervalContext();c.dps=100;count=0
    with mp.workdps(130):
        for tv in ('-11','-7.13','-2'):
            t=c.mpf(tv);us=c.mpf('.4')*c.exp(-t/2)
            ms=c.mpf('.17')*c.exp(-t);ds=c.mpf('.03')*c.exp(-t)
            aq=c.mpf('-.37')-t/2;hf=1+(c.mpf('.87')-1)*c.exp(t)
            for zv in ('-1','-.373','0','.51','1'):
                constant=lambda v:IntervalTaylor(c,[c.mpf(v)]+[c.mpf(0)]*5)
                delta,U,M,D,AZ,AQ,Hf=map(constant,('.001',us,ms,ds,'-.21',aq,hf))
                z=IntervalTaylor(c,[c.mpf(zv),c.mpf(1)]+[c.mpf(0)]*4);zero=z*0
                C=(1+z*z).reciprocal();L=1-delta*z*z
                A=(1-delta/2)*C+(1-delta)*z*z*C*C
                B=((2*delta*z*z-1)*C+2*(1-z*z)*z*z*C*C)/L
                memory=(1+z*z)*c.mpf('.017')
                u=[U*C*(-c.mpf('.5'))**j for j in range(5)]
                hist=dict(m=[M*z*(-1)**j for j in range(5)],
                    h=[U*C*(-D*(-c.mpf('1.5'))**j+(-c.mpf('.5'))**j) for j in range(5)],
                    k=[U*(M-4*D)*z*C*(-c.mpf('1.5'))**j for j in range(5)],
                    e=[U**2*(-1)**j*(z*z*AZ+C*C*(AQ+c.mpf(j)/2)) for j in range(5)])
                pressure=[-U**2*C*C*Hf/2+memory]+[-U**2*C*C*(-1)**j/2 for j in range(1,5)]
                raw=raw_pre_stress_rows(c,delta,z,u,[zero]*5,hist,pressure)
                canonical=sum(p['shape'][0] for name,p in raw['theta'].items() if name!='variable_radial_shear')/U
                expected=(A-C)/L+(-A/L-4*B)*D+(C+B)*M
                if not source.endpoints(canonical[0])[0]<=source.endpoints(expected[0])[1] or not source.endpoints(expected[0])[0]<=source.endpoints(canonical[0])[1]:
                    raise ArithmeticError('Full earlier-buffer theta correlation differs')
                R=c.mpf('1e12')*c.exp(t);Pstar=c.mpf(10000)
                def lift(parts):
                    return sum(R**c.mpf(str(p['mode'][0]))*Pstar**p['mode'][1]*p['shape'][0][0]/c.sqrt(2) for p in parts.values())
                theta,axial=lift(raw['theta']),lift(raw['axial'])
                F=Pstar*u[0][0]/c.sqrt(2*R)
                if source.endpoints(theta)[0]<=0 or source.endpoints(-2*F*theta)[1]>=0:
                    raise ArithmeticError('Independent relaxed source direction missing')
                strict=cone_margins(c,2,0,theta,axial)
                if strict['admitted'] or 'vs_minus2' not in strict['failed']:
                    raise ArithmeticError('Actual nonzero degenerate stress admitted strictly')
                if source.endpoints(strict['margins']['signed_quadratic'])[0]<=0:
                    raise ArithmeticError('Degenerate quadratic direction reserve missing')
                count+=1
    return count


def run():
    data=json.loads((source.HERE/source.NAME).read_bytes())
    for name,digest in data['input_hashes'].items():
        if source.sha(name)!=digest:raise ValueError('Changed relaxed buffer source: '+name)
    field=source.CurrentOriginalO2RelaxedBufferCone(require_checked=False);c=field.ctx
    for key,value in (('exact_original_buffer_relaxed_cone_source_theorem',field.theorem),
        ('whole_original_buffer_baseline',field.baseline),('whole_original_buffer_relaxed_cone',field.proof)):
        if source.parameters.encoded(value)!=data[key]:raise ValueError('Relaxed buffer proof changed: '+key)
    for name,value in field.baseline['positive_margins'].items():
        if source.endpoints(value)[0]<=0:raise ArithmeticError('Relaxed whole-buffer inequality missing: '+name)
    if field.proof['strict_vs_minus2']!=0 or field.proof['current_original_O2_buffer_strict_cone_certified'] is not False:
        raise ValueError('Nonzero kappa2 source relabeled strict')
    zeros=jet_zeros=0
    for N in (1,10**12):
        q,caps=field.errors.inputs('O2',(-11,-2),N)
        for key in ('du','V','dm','dh','dk','de','dp','dr'):
            for jk,value in caps[key].items():
                if source.endpoints(value)!=(0,0):
                    raise ArithmeticError('Pre-support local/history/radial source jet is nonzero: '+key+str(jk))
                jet_zeros+=1
        for label,parts in field.errors.model['stress'].items():
            for name,p in parts.items():
                cap=source.taper.o3.errors.bound_expression(c,p['expression'],caps,c.mpf([-1,1]),field.errors.delta,0,p['mode'][0])[0,0]
                if source.endpoints(cap)!=(0,0):raise ArithmeticError('Pre-support complete source error is not zero: '+name)
                zeros+=1
    specs=(('whole_buffer',(-11,-2),1),('turnoff_to_buffer_seam','-11',1),('interior','-7.13',37),
        ('exact_modulation_flat_edge','-2',22),('current_frequency',(-11,-2),10**12))
    for name,v,N in specs:
        q=field.query(v,N)
        if source.parameters.encoded(q)!=data['examples'][name]:raise ValueError('Relaxed buffer query changed: '+name)
        if q['complete_modified_signed_cone_certified_for_entire_query_box'] or not q['full_original_stress_not_zero']:
            raise ValueError('Relaxed query claimed a strict or zero-stress boundary')
    edge=field.taper.query('-2',22)
    if edge['complete_modified_signed_cone_certified_for_entire_query_box'] or source.endpoints(edge['actual_vs_minus2_source_lower'])!=(0,0):
        raise ValueError('O2 flat seam strict source shear changed')
    invalid=((-12,1),(-1,1),(('-11.001','-2'),1),(-2,0),(-2,True),(-2,1.5),(mp.inf,1))
    for v,N in invalid:
        try:field.query(v,N)
        except (ValueError,TypeError):pass
        else:raise ArithmeticError('Invalid relaxed buffer query admitted')
    for key in source.OPEN:
        if data.get(key) is not False:raise ValueError('Relaxed input promoted global target')
    fixtures=independent_full_stress_degenerate_fixtures()
    result=dict(all_passed=True,source_family=data['source_family'],
        exact_relaxed_buffer_source_and_cone_identities=len(field.theorem['identities']),
        strict_whole_buffer_baseline_inequalities=len(field.baseline['positive_margins']),
        independent_full_stress_relaxed_direction_and_strict_rejection_fixtures=fixtures,
        exact_pre_support_local_cumulative_radial_mixed4_error_zeros=jet_zeros,
        exact_pre_support_complete_tensor_error_zeros=zeros,closed_buffer_and_seam_queries_recomputed=len(specs),
        invalid_queries_rejected=len(invalid),unchanged_modulation_flat_edge_source_shear_verified=True,
        original_full_pressure_energy_moments_and_radial_shear_retained=True,
        current_original_O2_buffer_relaxed_input_cone_certified=True,
        current_original_O2_buffer_strict_cone_certified=False,**{key:False for key in source.OPEN},
        input_hashes={**data['input_hashes'],source.NAME:source.sha(source.NAME),
            Path(__file__).name:source.sha(Path(__file__).name)})
    (source.HERE/source.RECEIPT).write_bytes((json.dumps(source.parameters.encoded(result),indent=2)+'\n').encode())
    print('Original O2 relaxed buffer PASS:',fixtures,'full-stress cases;',zeros,'exact errors; strict admission remains false',flush=True)
    return result


if __name__=='__main__':run()
