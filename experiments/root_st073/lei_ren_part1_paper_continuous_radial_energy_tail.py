"""Pointwise materialized terminal transport and logarithmic radial energy.

This audits the numerical candidate, not the exact paper solution. Nonzero
rounding residuals are retained; input enclosures are needed to distinguish
exact functional closure from numerical defects. A point is not a Z integral.
"""
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_joined_outer import build_joined_field,_mp
from lei_ren_part1_paper_axial_correction import signed_log


def radial_tail_point(profile,logR,Z,*,nu='1',tau='1'):
    with mp.workdps(profile.precision):
        z=_mp(Z);delta=_mp(profile.schedule.delta);viscosity=_mp(nu);time=_mp(tau)
        if abs(z)>=1 or viscosity<=0 or time<=0:raise ValueError('Require |Z|<1, nu>0, tau>0')
        if profile.offset(logR,profile.schedule.logR_v)<=0:
            raise ValueError('Terminal tail point must be beyond Rv')
        R=mp.exp(_mp(logR));jet=profile.axial_average_jet(logR,z)
        mass=R*jet['value'];massZ=R*jet['derivative']
        d=1-z*z;L=1-delta*z*z
        C=((1-delta)*z*mass+d*massZ)/L
        q=time/d;r=mp.sqrt(2*viscosity*q*R)
        ur=-viscosity*C/r
        dz_dZ=mp.sqrt(viscosity)*time**((1-delta)/2)*L*d**(-(3-delta)/2)
        # 1/2 * 2*pi*r dr dz * ur^2; dr/r = dlogR/2.
        energy_density=mp.pi*viscosity**2*dz_dZ*C*C/2
        return dict(Mz=mass,Mz_Z=massZ,C=C,u_r=ur,r=r,dz_dZ=dz_dZ,
            energy_per_dZ_dlogR=energy_density,terminal_mean_forced_zero=False,
            exact_functional_terminal_closure_certified=False,
            physical_Z_interval_integrated=False,input_error_enclosed=False,
            finite_energy_certified=False,scale_recursion_certified=False)


def run():
    from lei_ren_part1_paper_continuous_incoming_outer import ContinuousIncomingOuterField
    folder=Path(__file__).parent
    print('building current shared field for terminal radial energy audit',flush=True)
    source=build_joined_field()
    seeded=json.loads((folder/'lei_ren_part1_paper_seeded_shared_candidate.json').read_text())
    atoms=json.loads((folder/'lei_ren_part1_paper_continuous_axial_solve.json').read_text())
    field=ContinuousIncomingOuterField(source,prepared={.3:(seeded,atoms)})
    p=field.outer;s=p.schedule
    with mp.workdps(field.precision):
        z=mp.mpf('.3');rows=[];values=[]
        for label,logR in [('post_Rv',p.log_at(s.logR_v,'50')),('heat',p.log_at(s.logR_tail,'4'))]:
            row=radial_tail_point(p,logR,z,nu=1,tau=mp.exp(-4))
            values.append(row)
            rows.append(dict(label=label,**{k:signed_log(row[k],60) for k in ('Mz','Mz_Z','C','u_r','energy_per_dZ_dlogR')}))
        errors={}
        for key in ('Mz','Mz_Z','C','energy_per_dZ_dlogR'):
            scale=max(abs(values[1][key]),abs(values[0][key]))
            errors[key]=mp.nstr(abs(values[1][key]-values[0][key])/scale if scale else 0,40)
        if max(map(mp.mpf,errors.values()))>mp.mpf('1e-60'):
            raise ArithmeticError('Terminal radial transport changed after axial support')
        identity=(abs(values[1]['u_r']*values[1]['r']/(-values[1]['C'])-1)
                  if values[1]['C'] else abs(values[1]['u_r']))
        if identity>mp.mpf('1e-60'):raise ArithmeticError('Physical radial tail scaling mismatch')
        report=dict(samples=rows,Z='.3',nu='1',tau='exp(-4)',
            post_support_constancy_relative_errors=errors,
            physical_ur_r_over_minus_nuC_relative_error=mp.nstr(identity,40),
            density_definition='dE_radial/(dZ dlogR) = pi*nu^2/2 * dz/dZ * C(Z)^2',
            nominal_pointwise_logarithmic_energy_growth_positive=values[1]['energy_per_dZ_dlogR']>0,
            nominal_global_radial_energy_obstruction_retained=values[1]['C']!=0,
            exact_functional_terminal_closure_certified=False,
            physical_Z_interval_integrated=False,input_error_enclosed=False,
            finite_energy_certified=False,scale_recursion_certified=False)
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'constancy':errors,'nonzero_C_retained':report['nominal_global_radial_energy_obstruction_retained'],'finite_energy_certified':False}),flush=True)
    return report


if __name__=='__main__':run()
