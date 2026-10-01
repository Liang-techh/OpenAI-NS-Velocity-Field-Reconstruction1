"""Locate the finite-core trace defect through its integrated angular equation.

Stored pressure-degree-zero coefficients only; MP diagnostics, not an infinite
remainder or directed source-roundoff certificate. Preserve the original datum.
"""
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_interval_collar_core_inlet import _read_cache, DEFAULT_CACHE
from lei_ren_part1_paper_collar_interval_receipt_reader import read_inlet
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints


def run():
    raw, path, digest = _read_cache(DEFAULT_CACHE)
    with mp.workdps(600):
        rows = raw['coefficients']; degree = len(rows['F'])-1
        z=mp.mpf(raw['Z']); dt=raw['delta']; L=1-dt*z*z; d=1-z*z
        Ra=4/raw['Lambda']; Faxis=rows['F'][0][0].atoms.get(0,mp.mpf(0))
        def scaled(name,k,normalization=1):
            return [row[k].atoms.get(0,mp.mpf(0))*Ra**n/normalization
                    for n,row in enumerate(rows[name])]
        A=scaled('F',0,Faxis); Az=scaled('F',1,Faxis)
        U=scaled('Uz',0); Uz=scaled('Uz',1)
        def add(*polys):
            return [sum(p[n] if n<len(p) else 0 for p in polys)
                    for n in range(max(map(len,polys)))]
        def scale(p,s):return [s*v for v in p]
        def mul(a,b):
            return [sum(a[i]*b[n-i] for i in range(len(a)) if 0<=n-i<len(b))
                    for n in range(len(a)+len(b)-1)]
        W=add([1],[-((1-dt)*z*u+d*uz)/(n+1) for n,(u,uz) in enumerate(zip(U,Uz))])
        H=add([(1-dt)*z/2],scale(U,d))
        rhs=add(mul(W,[(n+1)*v for n,v in enumerate(A)]),
                scale(add(A,scale(mul(U,A),-2*z)),dt/2),mul(H,Az))
        lhs=[2*L*(n+1)*(n+2)*A[n+1] for n in range(degree)]
        defect=add(lhs,scale(rhs,-Ra))
        contributions=[v/(L*(n+2)) for n,v in enumerate(defect)]
        integrated=sum(contributions)
        low=sum(abs(v) for v in contributions[:degree])
        high=sum(contributions[degree:])
        inlet,_=read_inlet();ctx=inlet['ctx']
        ring=inlet['I_theta'].value+2*inlet['R'].value*inlet['F_R'].value
        # The independent existing inlet formula uses the same retained rows,
        # but additionally propagates conditional pressure-integral error.
        direct=ring.component(0,0).value/ctx.mpf(Faxis)
        lo,hi=endpoints(direct)
        if not lo<=integrated<=hi:
            raise AssertionError('Integrated angular defect misses independent inlet interval')
        if low>=abs(high)*mp.mpf('1e-100'):
            raise AssertionError('Trusted recurrence orders fail to cancel at source precision')
        def s(v):return mp.nstr(v,80)
        report=dict(input_path=path.name,input_sha256=digest,precision=600,
            radial_degree=degree,pressure_parameter_degree=0,
            identity='d_R[R T_theta] = R angular_equation_defect / L',
            integrated_normalized_trace_T_over_axis_F=s(integrated),
            independent_trace_interval=[s(lo),s(hi)],
            low_order_contribution_absolute_sum=s(low),
            omitted_recurrence_order_contribution=s(high),
            contributions=[dict(radial_defect_order=n,integrated_trace_over_axis_F=s(v))
                           for n,v in enumerate(contributions)],
            independent_inlet_containment=True,
            original_pressure_datum_preserved=True,
            source_generation_roundoff_enclosed=False,
            infinite_radial_remainder_enclosed=False,
            exact_stress_free_core_certified=False,
            temporal_recursion=False,
            next_dependency='Enclose the analytic continuation beyond retained radial order; use the exact integrated identity, not a fitted trace reset')
        Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
        print('integrated T/Faxis',s(integrated),'low-order absolute sum',s(low),flush=True)
    return report


if __name__=='__main__':run()
