"""Coupled radial continuation with unchanged accepted axis and physical P0.

Pressure-parameter degree zero, at Z=.3. Finite MP diagnostics do not bound
the infinite tail or replace the actual temporal coefficient recursion.
"""
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_interval_collar_core_inlet import _read_cache, DEFAULT_CACHE
from lei_ren_part1_paper_core_recursion import core_coefficients


def trace_diagnostic(core,Ra,Faxis):
    z=core['Z'];dt=core['delta'];L=1-dt*z*z;d=1-z*z
    degree=core['radial_degree']
    def scaled(name,k,normalization=1):
        return [row[k]*Ra**n/normalization for n,row in enumerate(core[name])]
    A=scaled('F',0,Faxis);Az=scaled('F',1,Faxis)
    U=scaled('Uz',0);Uz=scaled('Uz',1)
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
    def moment(p):return sum(v/(n+2) for n,v in enumerate(p))
    inertial=-Ra*sum(W)*sum(A)/L
    inertial+=Ra/L*((1-dt/2)*moment(A)-(1-dt)*z*moment(Az)/2
                   -d*moment(add(mul(Az,U),mul(A,Uz)))
                   +(2*dt-1)*z*moment(mul(A,U)))
    shear=2*sum(n*v for n,v in enumerate(A))
    direct=inertial+shear
    if abs(direct-integrated)>max(abs(integrated),mp.mpf('1e-300'))*mp.mpf('1e-200'):
        raise AssertionError('Moment trace and integrated defect disagree')
    low=sum(abs(v) for v in contributions[:degree])
    if low>=abs(integrated)*mp.mpf('1e-100'):
        raise AssertionError('Computed recurrence orders do not cancel')
    return dict(trace=integrated,low=low,A=A,Az=Az,U=U,Uz=Uz,
                direct_identity_error=abs(direct-integrated),
                first_uncomputed_equation_contribution=contributions[degree])


def run():
    raw,path,digest=_read_cache(DEFAULT_CACHE)
    with mp.workdps(600):
        axis={name:[p.atoms.get(0,mp.mpf(0)) for p in raw['coefficients'][name][0]]
              for name in ('F','Uz','P')}
        maximum_degree=min(map(len,axis.values()))-2
        Ra=4/raw['Lambda'];Faxis=axis['F'][0]
        reports=[];previous=None;baseline=None
        for degree in range(18,maximum_degree+1):
            core=core_coefficients(raw['Z'],raw['delta'],F0_Z_taylor=axis['F'],
                 U0_Z_taylor=axis['Uz'],P0_Z_taylor=axis['P'],
                 radial_degree=degree,precision=600)
            if core['P'][0][0]!=axis['P'][0]:raise AssertionError('P0 changed')
            result=trace_diagnostic(core,Ra,Faxis)
            if baseline is None:
                baseline=result
                earlier=json.loads(Path(__file__).with_name('lei_ren_part1_paper_angular_defect_integral.json').read_text())
                if earlier['input_sha256']!=digest:raise AssertionError('Independent source receipt changed')
                reference=mp.mpf(earlier['integrated_normalized_trace_T_over_axis_F'])
                if abs(result['trace']-reference)>abs(reference)*mp.mpf('1e-75'):
                    raise AssertionError('Regenerated degree18 trace differs from committed source diagnostic')
            prefix_difference=max(abs(result[n][i]-baseline[n][i])
                for n in ('A','Az','U','Uz') for i in range(19))
            row=dict(radial_degree=degree,normalized_trace_over_axis_F=mp.nstr(result['trace'],80),
                     computed_order_absolute_defect_sum=mp.nstr(result['low'],80),
                     independent_moment_identity_error=mp.nstr(result['direct_identity_error'],80),
                     first_uncomputed_equation_contribution=mp.nstr(result['first_uncomputed_equation_contribution'],80),
                     retained_prefix_max_change=mp.nstr(prefix_difference,80),
                     ratio_to_previous=None if previous is None else mp.nstr(abs(result['trace']/previous),80))
            reports.append(row);previous=result['trace']
            print('radial degree',degree,'trace/Faxis',row['normalized_trace_over_axis_F'],flush=True)
        output=dict(input_path=path.name,input_sha256=digest,precision=600,Z='.3',
                    physical_P0_preserved=True,pressure_parameter_degree=0,
                    axis_rows_available={n:len(v) for n,v in axis.items()},results=reports,
                    coupled_swirl_and_pressure_terms_retained=True,
                    committed_degree18_trace_agreement=True,
                    source_generation_roundoff_enclosed=False,infinite_radial_remainder_enclosed=False,
                    extrapolated_tail_bound=False,core_C1_matching_certified=False,temporal_recursion=False,
                    required_next_dependency='Generate higher analytic axis jets and validate an infinite coupled radial remainder; finite-order ratios alone are insufficient')
        Path(__file__).with_suffix('.json').write_text(json.dumps(output,indent=2)+'\n')
        return output


if __name__=='__main__':run()
