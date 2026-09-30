"""Actual complete axial-square primitive after the disjoint end bumps.

The complete pulse and end energy atoms are shared with the live solve.
Measured inner axial energy is transported once. Partial pulse/end energy
and certified quadrature/terminal energy matching remain separate open work.
"""
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_joined_outer import _mp
from lei_ren_part1_paper_axial_correction import signed_log


class ContinuousAxialEnergyMoments:
    def __init__(self,profile):
        self.profile=profile;self.precision=profile.precision

    def moments_jet(self,logR,Z):
        p=self.profile;s=p.schedule
        with mp.workdps(self.precision):
            if p.offset(logR,s.logR_v)<0:
                raise ValueError('Complete axial-square provider starts at Rv; partial energy is not implemented')
            z=_mp(Z)
            if abs(z)>=1:raise ValueError('Actual inner seed requires |Z|<1')
            runtime=p.runtime(float(z));tangent=p.coefficient_tangent(float(z))
            mu=runtime.mu
            if 13-mp.mpf('3.15')*mu<=11:
                raise ArithmeticError('Pulse/end support separation failed')
            source=p.seed_source;inner=source.inner.evaluate_x(mp.e,z)
            Rref=mp.exp(_mp(str(s.logRref)))
            incoming=p.incoming_provider.full_incoming_rows(1)['I_uz2']
            offset=inner['raw_quadratic_integrals']['axial']-16*z*z*source.Rh
            offsetZ=inner['raw_quadratic_integrals']['axial_Z']-32*z*source.Rh
            prior=Rref*incoming*z*z+offset
            priorZ=2*Rref*incoming*z+offsetZ
            Ep=_mp(s.at_log_radius(s.logR_p,z)['log_angular_amplitude'])
            scale=mp.exp(_mp(str(s.logR_p))+2*Ep)
            LZ=-2*z/(1+z*z)
            pulse=runtime.a**2*runtime.Kp/mu
            pulseZ=2*runtime.a*tangent['a_Z']*runtime.Kp/mu
            # K includes exp(-26), the physical Rv Ev^2 / (Rp Ep^2)
            # transport, and the actual owned complete beta-square atoms.
            end=sum(k*c*c for k,c in zip(runtime.K,runtime.c))
            endZ=2*sum(k*c*cz for k,c,cz in zip(runtime.K,runtime.c,tangent['c_Z']))
            value=prior+scale*(pulse+end)
            derivative=priorZ+scale*(pulseZ+endZ+2*LZ*(pulse+end))
            angular=p.angular_moments_jet(logR,z)
            return dict(axial_square=value,axial_square_Z=derivative,
                z_theta=value-angular['swirl_energy']/2,
                z_theta_Z=derivative-angular['swirl_energy_Z']/2,
                incoming_axial_square=prior,incoming_axial_square_Z=priorZ,
                pulse_normalized=pulse,end_normalized=end,
                actual_inner_offsets_reapplied_once=True,owned_complete_atoms=True,
                pulse_end_cross_term_zero_by_disjoint_support=True,
                partial_axial_energy_implemented=False,quadrature_error_enclosed=False,
                complete_five_moments=False,finite_energy_certified=False,
                scale_recursion_certified=False)


def run():
    from lei_ren_part1_paper_continuous_incoming_outer import ContinuousIncomingOuterField
    from lei_ren_part1_paper_joined_outer import build_joined_field
    print('building shared field for complete axial energy',flush=True)
    source=build_joined_field();folder=Path(__file__).parent
    seeded=json.loads((folder/'lei_ren_part1_paper_seeded_shared_candidate.json').read_text())
    atoms=json.loads((folder/'lei_ren_part1_paper_continuous_axial_solve.json').read_text())
    field=ContinuousIncomingOuterField(source,prepared={.3:(seeded,atoms)})
    p=field.outer;engine=p.axial_energy_moment_provider
    with mp.workdps(field.precision):
        z=mp.mpf('.3');center=p.log_at(source.schedule.logR_v,'50')
        row=engine.moments_jet(center,z)
        tail=engine.moments_jet(str(source.schedule.logR_tail),z)
        constancy=[mp.nstr(abs(tail[key]/row[key]-1),40) for key in ('axial_square','axial_square_Z')]
        if max(map(mp.mpf,constancy))>mp.mpf('1e-70'):
            raise ArithmeticError('Axial square primitive changed after its support')
        h=mp.mpf('1e-4')
        neighbors={i:engine.moments_jet(p.log_at(center,i*h),z) for i in (-2,-1,1,2)}
        velocity=p.values_with_jets(center,z);R=mp.exp(_mp(str(center)))
        targets=[R*(velocity['Uz']**2-velocity['Utheta']**2/2),
                 R*(2*velocity['Uz']*velocity['Uz_Z']-velocity['Utheta']*velocity['Utheta_Z'])]
        errors=[]
        for key,target in zip(('z_theta','z_theta_Z'),targets):
            derivative=sum(w*neighbors[i][key] for i,w in ((-2,1),(-1,-8),(1,8),(2,-1)))/(12*h)
            errors.append(mp.nstr(abs(derivative/target-1),40))
        if max(map(mp.mpf,errors))>mp.mpf('1e-8'):
            raise ArithmeticError('Quadratic moment derivative differs from actual velocity integrand')
        report=dict(axial_square=signed_log(row['axial_square'],field.precision),
            axial_square_Z=signed_log(row['axial_square_Z'],field.precision),
            post_support_constancy_relative_errors=constancy,
            quadratic_integrand_relative_errors=errors,
            tail_z_theta=signed_log(tail['z_theta'],field.precision),
            tail_z_theta_Z=signed_log(tail['z_theta_Z'],field.precision),
            owned_complete_atoms=True,actual_inner_offsets_reapplied_once=True,
            partial_axial_energy_implemented=False,quadrature_error_enclosed=False,
            complete_five_moments=False,finite_energy_certified=False,scale_recursion_certified=False)
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'post_support_constancy':constancy,'quadratic_integrand_relative_errors':errors}),flush=True)
    return report


if __name__=='__main__':run()
