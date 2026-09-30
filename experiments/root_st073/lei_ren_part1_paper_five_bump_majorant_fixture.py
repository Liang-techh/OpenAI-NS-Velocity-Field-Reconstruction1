"""Resolved scalar convergence versus the conditional positive majorant."""
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_five_bump_map import FiveBumpMomentMap
from lei_ren_part1_paper_five_bump_majorant import compute_majorant
from lei_ren_part1_paper_five_bump_response import build_response


def run():
    with mp.workdps(100):
        m=FiveBumpMomentMap(precision=100,order=96)
        zero=compute_majorant(m,0,1,degree=3);threshold=zero['raw']['condition_threshold']
        boundary=compute_majorant(m,threshold,1,degree=3)
        failed=compute_majorant(m,2*threshold,1,degree=3)
        assert zero['raw']['tail_bound']==0 and boundary['raw']['contraction_condition_passed']
        assert abs(boundary['raw']['lipschitz_bound']-mp.mpf('.5'))<mp.mpf('1e-95')
        assert not failed['raw']['contraction_condition_passed']
        e=threshold/1000;d=[e*v/15 for v in (1,-2,3,-4,5)];am=mp.mpf('.5')
        bounds=compute_majorant(m,e,1,degree=3)
        response=build_response(m,am,degree=3);terms=response.evaluate_terms(d)
        approximation=[mp.fsum(v[i] for v in terms.values()) for i in range(5)]
        h=[mp.mpf(0)]*5
        for k in range(32):
            q=m.quadratic(h,am)
            h=[-v for v in m.linear_inverse([a+b for a,b in zip(d,q)])]
        error=mp.fsum(abs(a-b) for a,b in zip(h,approximation))
        assert error<=bounds['raw']['tail_bound']
        report=dict(zero_tail=True,boundary_pass=True,above_threshold_rejected=True,
          resolved_scalar_response_error=mp.nstr(error,40),conditional_tail_bound=mp.nstr(bounds['raw']['tail_bound'],40),
          constants=bounds['metadata'],actual_uniform_defect_norm_verified=False,actual_closure=False)
        Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps({k:v for k,v in report.items() if k!='constants'},indent=2))
        return report

if __name__=='__main__':run()
