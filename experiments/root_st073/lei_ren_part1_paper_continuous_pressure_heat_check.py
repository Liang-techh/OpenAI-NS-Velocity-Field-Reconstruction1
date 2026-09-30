"""Check pressure-weighted heat primitives from the common heat owner."""
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_joined_outer import build_joined_field,_mp
from lei_ren_part1_paper_continuous_incoming_outer import ContinuousIncomingOuterField
from lei_ren_part1_paper_continuous_axial_pulse import ContinuousAxialPulse
from lei_ren_part1_paper_axial_correction import signed_log


def run():
    folder=Path(__file__).parent
    print('building shared field for pressure heat primitives',flush=True)
    source=build_joined_field()
    seeded=json.loads((folder/'lei_ren_part1_paper_seeded_shared_candidate.json').read_text())
    atoms=json.loads((folder/'lei_ren_part1_paper_continuous_axial_solve.json').read_text())
    field=ContinuousIncomingOuterField(source,prepared={.3:(seeded,atoms)})
    p=field.outer;s=p.schedule;engine=p.angular_moment_provider.heat_provider
    with mp.workdps(field.precision):
        z=mp.mpf('.3');step=mp.mpf('1e-4');samples=[]
        E=_mp(s._log_c_inf)-(1+_mp(s.delta))*_mp(s.logR_tail)/2
        for t in (mp.mpf('.6'),mp.mpf(4)):
            neighbors={i:engine.pressure_increments(t+i*step,z) for i in (-2,-1,1,2)}
            sw=ContinuousAxialPulse.sigma_pair(t)[0];eps=_mp(s.epsilon)
            edge=(3-t)/2;phi=mp.exp(-1/edge**2) if edge>0 else mp.mpf(0)
            factor=1-eps*phi;K0=(1-sw)*(1-eps)+sw*factor
            heat=p.angular_schedule_provider.heat_jet(p.log_at(s.logR_tail,t),z)
            dk=-sw*factor*heat['deficit'];dkz=sw*factor*heat['H_prime']*heat['xi_Z']
            scale=mp.exp(2*E-(1+_mp(s.delta))*t)/2
            targets=[scale*K0*K0,scale*(2*K0*dk+dk*dk),scale*2*(K0+dk)*dkz]
            errors=[]
            for key,target in zip(('reference','heat_correction','heat_correction_Z'),targets):
                derivative=sum(w*neighbors[i][key] for i,w in ((-2,1),(-1,-8),(1,8),(2,-1)))/(12*step)
                errors.append(abs(derivative/target-1))
            if max(errors)>mp.mpf('1e-8'):raise ArithmeticError('Pressure heat cumulative integrand mismatch')
            samples.append(dict(t=mp.nstr(t,20),integrand_relative_errors=[mp.nstr(v,40) for v in errors]))
        tail=engine.complete_pressure_heat_integral(z)
        if not(tail['reference']>0 and tail['heat_correction']<0):raise ArithmeticError('Pressure tail sign lost')
        neighbors={i:engine.complete_pressure_heat_integral(z+i*step)['heat_correction'] for i in (-2,-1,1,2)}
        dZ=sum(w*neighbors[i] for i,w in ((-2,1),(-1,-8),(1,8),(2,-1)))/(12*step)
        error=abs(dZ/tail['heat_correction_Z']-1)
        if error>mp.mpf('1e-8'):raise ArithmeticError('Complete pressure heat Z jet mismatch')
        origin=engine.pressure_increments(0,z)
        if any(origin[k] for k in ('reference','heat_correction','heat_correction_Z')):
            raise ArithmeticError('Pressure heat increment nonzero at Rtail')
        report=dict(samples=samples,complete_tail_Z_relative_error=mp.nstr(error,40),
            pressure_heat_tail={k:signed_log(tail[k],50) for k in ('reference','heat_correction','heat_correction_Z','analytic_heat_truncation_absolute_bound')},
            Rtail_increment_zero=True,heat_kernel_shared=True,
            infinite_pressure_heat_tail_integrated=True,
            actual_inner_pressure_moment_installed=False,quadrature_error_enclosed=False,
            arithmetic_error_enclosed=False,complete_five_moments=False,
            finite_energy_certified=False,scale_recursion_certified=False)
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'samples':samples,'complete_tail_Z_relative_error':report['complete_tail_Z_relative_error']}),flush=True)
    return report


if __name__=='__main__':run()
