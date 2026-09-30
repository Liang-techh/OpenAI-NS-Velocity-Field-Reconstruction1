"""Independent resolved root and first-Z checks of incremental inverse."""
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_five_bump_map import FiveBumpMomentMap
from lei_ren_part1_paper_five_bump_inverse import iterate_inverse
from lei_ren_part1_paper_five_bump_majorant import compute_majorant
from lei_ren_part1_paper_axial_dual import AxialDual


def run():
    with mp.workdps(120):
        moment_map=FiveBumpMomentMap(precision=120,order=96)
        z=mp.mpf('.3')
        def inputs(x):
            return tuple(mp.mpf('1e-10')*(i+1)*(1+x*x)/(i+2) for i in range(5))
        def amplitude(x):return 20/(1+x*x)
        def root(x):
            d=inputs(x);a=amplitude(x)
            seed=tuple(-v for v in moment_map.linear_inverse(d))
            return mp.findroot(lambda *h:tuple(v+w for v,w in zip(moment_map.apply(h,a),d)),
                               seed,J=lambda *h:mp.matrix(moment_map.jacobian(h,a)),
                               tol=mp.mpf('1e-220'),maxsteps=20)
        dual_d=[AxialDual(v,mp.diff(lambda x:inputs(x)[i],z)) for i,v in enumerate(inputs(z))]
        dual_a=AxialDual(amplitude(z),mp.diff(amplitude,z))
        # Resolved family has C1 l1 input bound < 1e-8 on |Z|<=1.
        majorant=compute_majorant(moment_map,'1e-8',20)
        result=iterate_inverse(moment_map,dual_d,dual_a,steps=16,majorant=majorant)
        oracle=root(z)
        value_error=max(abs(h.value.evaluate()-oracle[i]) for i,h in enumerate(result['h']))
        step=mp.mpf('1e-15')
        roots={k:root(z+k*step) for k in (-2,-1,1,2)}
        first_Z_error=max(abs(h.tangent.evaluate()-(roots[-2][i]-8*roots[-1][i]
                          +8*roots[1][i]-roots[2][i])/(12*step))
                          for i,h in enumerate(result['h']))
        direct=moment_map.apply(result['h'],dual_a)
        replay_error=max(abs((direct[i]+dual_d[i]-result['terminal_residual'][i]).value.evaluate())
                         for i in range(5))
        print('value / first-Z / replay errors',value_error,first_Z_error,replay_error,flush=True)
        assert value_error<mp.mpf('1e-60') and first_Z_error<mp.mpf('1e-60')
        assert replay_error<mp.mpf('1e-110')
        assert value_error<result['conditional_C1_error_bound']
        report=dict(value_error=mp.nstr(value_error,60),first_Z_error=mp.nstr(first_Z_error,60),
                    residual_replay_error=mp.nstr(replay_error,60),
                    conditional_bound=mp.nstr(result['conditional_C1_error_bound'],60),
                    independent_root_solver=True,first_Z_oracle='four-point stencil of independent root solver, step=1e-15',
                    source_scope='resolved analytic surrogate only',actual_source_certified=False)
        Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps(report,indent=2),flush=True)


if __name__=='__main__':run()
