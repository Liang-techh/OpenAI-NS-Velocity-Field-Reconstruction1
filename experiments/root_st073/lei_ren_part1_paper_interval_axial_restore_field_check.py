"""Independent partial restoration integrals and inlet continuity check."""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_interval_axial_restore_field import IntervalAxialRestoreField,partial_restore_integrals
from lei_ren_part1_paper_interval_long_reshape_field import RM_LOG_EXACT
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent

def run():
    c=MPIntervalContext();c.dps=80
    with mp.workdps(100):
        t=mp.mpf('.5');actual=partial_restore_integrals(c,c.mpf('.5'))
        def cutoff(s):
            if s<=0:return mp.mpf(1)
            if s>=1:return mp.mpf(0)
            return 1/(1+mp.exp(1/(1-s)**2-1/s**2))
        references=dict(K1=mp.quad(lambda s:mp.exp(-(t-s))*cutoff(s),[0,t/2,t]),
            K16=mp.quad(lambda s:mp.exp(-mp.mpf('1.6')*(t-s))*cutoff(s),[0,t/2,t]),
            K2=mp.quad(lambda s:mp.exp(-(t-s))*cutoff(s)**2,[0,t/2,t]))
        for name,value in references.items():
            lo,hi=endpoints(actual[name])
            if not lo<=value<=hi:raise AssertionError(('independent partial primitive excluded',name))
    field=IntervalAxialRestoreField()
    with mp.workdps(field.calc.precision+60):
        at_zero=field.evaluate_phase(0);start=field.start;count=0
        for key in ('F','Utheta','Utheta_y','Uz','Uz_y','P'):
            for order in (0,1):
                if endpoints(at_zero[key][order])!=endpoints(start[key][order]):raise AssertionError(('inlet field differs',key,order))
                count+=1
        for key in ('z','theta','theta_z','z_theta','p'):
            for order in (0,1):
                if endpoints(at_zero['physical_moments'][key][order])!=endpoints(start['physical_moments'][key][order]):
                    raise AssertionError(('inlet moment differs',key,order))
                count+=1
        one=field.evaluate_phase(1)
        if any(endpoints(one['Uz'][k])!=endpoints((field.z*4)[k]) for k in (0,1)):
            raise AssertionError('axial restoration did not reach 4Z identically')
        if any(endpoints(one['Uz_y'][k])!=(mp.mpf(0),mp.mpf(0)) for k in (0,1)):
            raise AssertionError('terminal axial radial derivative nonzero')
        dispatched=field.evaluate_log_offset(RM_LOG_EXACT)
        terminal=field.evaluate_phase(2)
        if any(endpoints(dispatched['P'][k])!=endpoints(terminal['P'][k]) for k in (0,1)):
            raise AssertionError('unified endpoint dispatcher differs')
    result=dict(independent_partial_restore_primitives_contained=3,
        exact_C1_inlet_field_and_moment_coefficients=count,
        restored_axial_field_exactly_4Z=True,terminal_axial_radial_derivative_zero=True,
        unified_R110_to_Rm_dispatcher_terminal_preserved=True,
        whole_connecting_cone_or_higher_smoothness_certified=False)
    names=(Path(__file__).name,'lei_ren_part1_paper_interval_axial_restore_field.py',
        'lei_ren_part1_paper_interval_long_reshape_field.py','lei_ren_part1_paper_interval_functional_defects.json')
    result['input_hashes']={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names}
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print('Axial restoration checks:',result,flush=True)
    return result

if __name__=='__main__':run()
