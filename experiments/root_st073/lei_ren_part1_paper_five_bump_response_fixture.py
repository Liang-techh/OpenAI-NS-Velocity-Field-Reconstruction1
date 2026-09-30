"""Finite formal response: coefficient replay and independent scalar contraction."""
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_five_bump_map import FiveBumpMomentMap
from lei_ren_part1_paper_five_bump_response import build_response
from lei_ren_part1_paper_axial_dual import AxialDual


def run():
    with mp.workdps(100):
        moment_map=FiveBumpMomentMap(precision=100,order=96)
        am=mp.mpf(5);degree=3
        response=build_response(moment_map,am,degree=degree)
        coeff=response.coefficients
        residual={}
        def add(key,values):
            old=residual.get(key,[mp.mpf(0)]*5)
            residual[key]=[a+b for a,b in zip(old,values)]
        for key,values in coeff.items():
            add(key,[sum(moment_map.matrix[i,j]*values[j] for j in range(5)) for i in range(5)])
        for ka,a in coeff.items():
            for kb,b in coeff.items():
                key=tuple(x+y for x,y in zip(ka,kb))
                if sum(key)<=degree:add(key,moment_map.bilinear(a,b,am))
        for i in range(5):
            key=tuple(int(i==j) for j in range(5));add(key,[mp.mpf(int(i==j)) for j in range(5)])
        coefficient_error=max(abs(v) for rows in residual.values() for v in rows)
        assert coefficient_error<mp.mpf('1e-65')
        d=[mp.mpf(x)*mp.mpf('1e-12') for x in (1,-2,3,-4,5)]
        terms=response.evaluate_terms(d)
        approximated=[mp.fsum(rows[i] for rows in terms.values()) for i in range(5)]
        h=[mp.mpf(0)]*5
        for iteration in range(24):
            q=moment_map.quadratic(h,am)
            h=[-v for v in moment_map.linear_inverse([a+b for a,b in zip(d,q)])]
        error=mp.fsum(abs(a-b) for a,b in zip(approximated,h))/mp.fsum(abs(v) for v in h)
        assert error<mp.mpf('1e-16')
        dz=[mp.mpf(x)*mp.mpf('1e-13') for x in (2,1,-1,3,-2)];amz=mp.mpf('.7')
        dual_response=build_response(moment_map,AxialDual(am,amz,pressure_order=0,width_order=0),degree=degree)
        dual_terms=dual_response.evaluate_terms([AxialDual(a,b,pressure_order=0,width_order=0) for a,b in zip(d,dz)])
        predicted=[mp.fsum(rows[i].tangent.component(0,0) for rows in dual_terms.values()) for i in range(5)]
        step=mp.mpf('1e-10')
        def shifted(k):
            r=build_response(moment_map,am+k*step*amz,degree=degree)
            t=r.evaluate_terms([a+k*step*b for a,b in zip(d,dz)])
            return [mp.fsum(rows[i] for rows in t.values()) for i in range(5)]
        samples={k:shifted(k) for k in (-2,-1,1,2)}
        tangent=[(-samples[2][i]+8*samples[1][i]-8*samples[-1][i]+samples[-2][i])/(12*step) for i in range(5)]
        tangent_error=mp.fsum(abs(a-b) for a,b in zip(predicted,tangent))/mp.fsum(abs(v) for v in tangent)
        assert tangent_error<mp.mpf('1e-25')
        report=dict(formal_defect_degree=degree,coefficient_count=len(coeff),
          coefficientwise_residual_absolute_error=mp.nstr(coefficient_error,30),
          independent_scalar_contraction_relative_error=mp.nstr(error,30),
          independent_first_Z_relative_error=mp.nstr(tangent_error,30),
          quadrature_enclosed=False,degree_remainder_enclosed=False,
          actual_source_closure=False,temporal_recursion=False)
        Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps(report,indent=2),flush=True)
        return report

if __name__=='__main__':run()
