"""Representable-scale equivalence and nonzero-swirl checks of scaled step."""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_functional_core_step import initial_rows, advance_one
from lei_ren_part1_paper_logarithmic_core_step import advance_scaled_one
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent

def run():
    c=MPIntervalContext();c.dps=160
    comparisons=0;case_results=[];maximum=c.mpf(0)
    with mp.workdps(200):
        for center in ('0','.5','-.75'):
            z=c.mpf(center);length=10;zero=c.mpf(0)
            ell=[c.mpf('.25')]+[zero]*(length-1)
            # S is a nonzero analytic jet satisfying S'=2 ell S.
            S=[c.mpf('.0625')*c.mpf('.5')**k/mp.factorial(k) for k in range(length)]
            U=[4*z+c.mpf('.125'),c.mpf(4)]+[zero]*(length-2)
            P=[c.mpf(-3),c.mpf('.5'),c.mpf('.25'),c.mpf('-.125')]+[zero]*(length-4)
            fixed=dict(ell_Z_taylor=ell,S_Z_taylor=S,U0_Z_taylor=U,P0_Z_taylor=P)
            for lam_string in ('8','500','1e8'):
                lam=c.mpf(lam_string);eps=1/lam
                scaled=dict(ell_Z_taylor=[eps*x for x in ell],
                    S_Z_taylor=[eps**2*x for x in S],U0_Z_taylor=U,
                    P0_Z_taylor=[eps*x for x in P])
                physical=initial_rows(c,fixed,z,6,required_depth=3)
                direct=initial_rows(c,scaled,z,6,required_depth=3)
                for n in range(6):
                    advance_one(c,fixed,physical,n,z,c.mpf('.125'))
                    advance_scaled_one(c,scaled,direct,n,z,c.mpf('.125'),eps)
                    for name in ('A','Uz','P'):
                        for x,y in zip(physical[name][n+1],direct[name][n+1]):
                            expected=x*eps**(n+1)*(eps if name=='P' else 1)
                            a,b=endpoints(expected);d,e=endpoints(y)
                            if max(a,d)>min(b,e):raise ArithmeticError('nonoverlap: '+name)
                            error=max(abs(a-d),abs(b-e))/max(mp.mpf(1),abs(a),abs(b))
                            if error>mp.mpf('1e-140'):raise ArithmeticError('scaled mismatch')
                            maximum=c.mpf([0,max(endpoints(maximum)[1],error)])
                            comparisons+=1
                if not endpoints(direct['P'][1][0])[0]>0:
                    raise ArithmeticError('nonzero swirl pressure lost')
                case_results.append(dict(center=center,Lambda=lam_string,degree=6,
                    all_coefficients_overlap=True,swirl_pressure_nonzero=True))
        result=dict(cases=case_results,coefficient_comparisons=comparisons,
            maximum_absolute_or_relative_endpoint_difference_upper=str(endpoints(maximum)[1]),
            scaled_step_matches_frozen_physical_step=True,
            scale_definitions=dict(A='physical_A[n]/Lambda**n',Uz='physical_Uz[n]/Lambda**n',
                P='epsilon*physical_P[n]/Lambda**n',ell='epsilon*physical_ell',
                S='epsilon**2*F0**2'),
            Md40_seed_generated=False,degree110_core_generated=False,
            temporal_recursion=False,
            input_hashes={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in (
                Path(__file__).name,'lei_ren_part1_paper_logarithmic_core_step.py',
                'lei_ren_part1_paper_functional_core_step.py',
                'lei_ren_part1_paper_schedule_endpoint_enclosures.py')})
    p=Path(__file__).with_suffix('.json');p.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print('PASS:',comparisons,'coefficient comparisons;',len(case_results),'cases; nonzero swirl retained')
    return result

if __name__=='__main__':run()
