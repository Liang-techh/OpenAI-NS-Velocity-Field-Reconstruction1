"""Link the accepted finite collar's first-width atoms to analytic integrals.

This checks stored finite-model coefficients at Z=.3, not true source errors
or the second-width and full nonlinear ODE remainders.
"""
import argparse
import hashlib
import json
import pickle
from pathlib import Path
from types import SimpleNamespace
import mpmath as mp
from lei_ren_part1_paper_pressure_width_second_axial_comparison import SecondAxialPressureWidthComparison
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode


def run(cache_dir):
    cache_dir=Path(cache_dir)
    core_path=cache_dir/'second_axial_core_Z03.pkl'
    source_path=cache_dir/'second_axial_R100_source_Z03.pkl'
    # Only load trusted locally generated project caches.
    with core_path.open('rb') as f: core=pickle.load(f)
    with source_path.open('rb') as f: source=pickle.load(f)
    with mp.workdps(260):
        width=mp.exp(-100-100*mp.mpf('1e152'))
        axis=SimpleNamespace(precision=260,Lambda=core['Lambda'],delta=core['delta'])
        def coefficients(z):
            if mp.mpf(str(z))!=mp.mpf('.3'): raise ValueError('This cache has only Z=.3')
            return core['coefficients']
        component=SimpleNamespace(axis=axis,precision=260,coefficients=coefficients)
        comparison=SecondAxialPressureWidthComparison(
            dict(precision=260,axis=axis,component_pressure_core=component),
            h_b=width,pressure_order=9,width_order=2)
        bar=comparison.evaluate(0,'.3')
    iv=mp.ctx_iv.MPIntervalContext();iv.dps=300
    def exact(value):
        # Direct MP conversion preserves the stored binary number.
        return iv.mpf(value)
    def slots(jet):
        return [[exact(getattr(jet,name).atoms.get((p,0),0)) for p in range(10)]
                for name in ('value','tangent','second')]
    def scale(rows,factor):return [[v*factor for v in row] for row in rows]
    def product(a,b):
        return [sum((a[k]*b[p-k] for k in range(p+1)),iv.mpf(0)) for p in range(10)]
    with mp.workdps(340):
        u=slots(bar['Uz']);u2=[product(u[0],u[0]),scale([product(u[0],u[1])],2)[0],
            [2*(a+b) for a,b in zip(product(u[1],u[1]),product(u[0],u[2]))]]
        constant=lambda c:[[iv.mpf(c) if p==0 and z==0 else iv.mpf(0) for p in range(10)] for z in range(3)]
        expected=[scale(slots(bar['D']),iv.mpf('-0.25')),
            scale(slots(bar['I_z']),-iv.sqrt(exact(comparison.Ra)/2)/2),
            constant(4),scale(u,2),scale(u,4),scale(u2,2),constant(2),constant(2)]
        start=source['functions']['start']
        if start['s']!=2 or start['steps']!=8:raise ValueError('Unexpected collar endpoint or step count')
        rows=[];max_relative=mp.mpf(0);nonzero_zero_targets=0;failed=[]
        for state_name,state,target in zip(('g','u','theta','mz','mixed','axial','swirl','p'),start['normalized_state'],expected):
            for z,name in enumerate(('value','tangent','second')):
                for p in range(10):
                    stored=exact(getattr(state,name).atoms.get((p,1),0))/exact(width)
                    difference=stored-target[z][p];lo,hi=endpoints(difference)
                    error=max(abs(lo),abs(hi));tlo,thi=endpoints(target[z][p])
                    lower=min(abs(tlo),abs(thi)) if tlo*thi>0 else mp.mpf(0)
                    relative=error/lower if lower else None
                    if lower: max_relative=max(max_relative,relative)
                    elif error: nonzero_zero_targets+=1
                    if relative is not None and relative>=mp.mpf('1e-250'):
                        failed.append(dict(state=state_name,Z_derivative=z,pressure_power=p,
                            relative_difference_upper=relative))
                    rows.append(dict(state=state_name,Z_derivative=z,pressure_power=p,
                        stored_width_coefficient_divided_by_h=stored,analytic_expected=target[z][p],
                        difference_interval=difference,absolute_difference_upper=error,
                        relative_difference_upper=relative))
        report=dict(Z='.3',endpoint_s=2,RK_steps=8,pressure_order=9,width_order=2,
            cache_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (core_path,source_path)},
            coefficients_compared=len(rows),maximum_relative_difference_upper=max_relative,
            nonzero_differences_at_zero_targets=nonzero_zero_targets,
            stored_generation_roundoff_tolerance='1e-250',
            within_stored_generation_roundoff=(max_relative<mp.mpf('1e-250') and nonzero_zero_targets==0),
            failed_relative_comparisons=failed,
            analytic_first_width_replacement_coefficients_available=True,
            original_cache_modified=False,
            finite_model_point_Z_only=True,source_driver_error_enclosed=False,
            second_width_ODE_error_enclosed=False,full_ODE_error_enclosed=False,
            temporal_recursion=False,coefficients=rows)
        Path(__file__).with_suffix('.json').write_text(json.dumps(encode(report),indent=2)+'\n')
        print('coefficients',len(rows),'relative difference',mp.nstr(max_relative,18),
              'zero-target differences',nonzero_zero_targets,flush=True)
        # A discrepancy is scientific evidence. Preserve every comparison and
        # the directed analytic replacement instead of blessing the old cache.
        print('coefficients outside relative tolerance',len(failed),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--cache-dir',default='work_paper_cache')
    run(parser.parse_args().cache_dir)
