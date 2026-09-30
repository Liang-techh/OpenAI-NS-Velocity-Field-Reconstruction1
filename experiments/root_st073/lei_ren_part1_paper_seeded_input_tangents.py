"""Transport measured inner moment jets into the normalized axial equations.

Uses the source factorization Uz proportional to Z and Utheta proportional
to (1+Z^2)^-1, and requires the future angular-energy derivative explicitly.
No missing derivative is silently replaced by zero.
"""
import mpmath as mp
from pathlib import Path
import json
from lei_ren_part1_paper_axial_correction import from_signed_log, signed_log


def seeded_input_tangents(seeded,terminal,Z,*,y_p,mu,future_energy_Z,
                          reference_linear_factors=None,precision=443):
    with mp.workdps(precision):
        z=mp.mpf(Z)
        if z==0 and reference_linear_factors is None:
            raise ValueError('Z=0 requires reference factor coefficients rather than value/Z')
        incoming=seeded['incoming'];transport=incoming['inner_seed_transport']
        R=mp.mpf(terminal['R']);u=mp.sqrt(2*R)*terminal['F']
        uZ=mp.sqrt(2*R)*terminal['FZ']
        theta=mp.mpf(5)/8*R*mp.sqrt(2*R)*u
        thetaZ=mp.mpf(5)/8*R*mp.sqrt(2*R)*uZ
        reference_Z=dict(z=4*R,theta_z=4*theta+4*z*thetaZ,
                         z_theta=32*z*R-mp.mpf(5)/6*R*u*uZ)
        deltaZ={key:terminal['momentsZ'][key]-reference_Z[key] for key in reference_Z}
        delta={key:from_signed_log(value) for key,value in transport['raw_offsets'].items()}
        rp=mp.mpf(transport['log_Rp']);ep=mp.mpf(incoming['log_Ep'])
        LZ=-2*z/(1+z*z);mu=mp.mpf(mu)
        yp=mp.mpf(y_p);dimensions=incoming['dimensionless_integrals']
        # Recover the reference directly, never subtract the large inner
        # seed from an already rounded seeded sum to obtain a tiny input.
        original1=mp.mpf(dimensions['I_z'])*mp.exp(-yp-ep)
        original2=mp.mpf(dimensions['I_theta_z'])*mp.exp(-mp.mpf('1.5')*yp-2*ep)
        original1Z=((1/z-LZ)*original1 if z else
                    mp.mpf(reference_linear_factors['I_z_over_Z'])*mp.exp(-yp-ep))
        original2Z=((1/z-LZ)*original2 if z else
                    mp.mpf(reference_linear_factors['I_theta_z_times_1plusZ2_over_Z'])*mp.exp(-mp.mpf('1.5')*yp-2*ep))
        m1Z=original1Z+mp.exp(-rp-ep)*(deltaZ['z']-LZ*delta['z'])
        m2Z=original2Z+mp.exp(-mp.mpf('1.5')*rp-2*ep)/mp.sqrt(2)*(
            deltaZ['theta_z']-2*LZ*delta['theta_z'])
        # The reference swirl energy normalized by Ep^2 is Z independent.
        # Its derivative cancels analytically, avoiding subtracting its two
        # huge close terms. The axial reference energy has factor Z^2/Ep^2.
        Iuz2=mp.mpf(incoming['dimensionless_integrals']['I_uz2'])
        reference_prior_Z=mu*mp.exp(-yp-2*ep)*Iuz2*(2/z-2*LZ) if z else mp.mpf(0)
        delta_prior_Z=mu*mp.exp(-rp-2*ep)*(deltaZ['z_theta']-2*LZ*delta['z_theta'])
        priorZ=reference_prior_Z+delta_prior_Z
        rows=incoming['row_normalization']
        baseZ=[m1Z*mp.exp(mp.mpf(rows['log_scale1'])),
               m2Z*mp.exp(mp.mpf(rows['log_scale2']))]
        targetZ=-priorZ+mp.mpf(future_energy_Z)
        return dict(base_Z=[signed_log(x,precision) for x in baseZ],
            energy_target_Z=signed_log(targetZ,precision),
            prior_energy_Z=signed_log(priorZ,precision),
            inner_offset_Z={k:signed_log(v,precision) for k,v in deltaZ.items()},
            log_Ep_Z=mp.nstr(LZ,precision),
            future_energy_Z=signed_log(mp.mpf(future_energy_Z),precision),
            input_derivative_enclosure_certified=False,
            unresolved_inner_input_entries=terminal.get('unresolved_input_entries'),
            scope='Actual nominal inner moment jets and source-factorized incoming jets; caller-supplied future energy derivative remains separately sourced.')


def actual_seeded_input_tangents(source,seeded,Z):
    """Actual inner jets with explicitly finite-differenced future energy.

    The remaining angular derivative gap is retained in provenance. No new
    physical coefficient solves at neighboring Z are needed for this path.
    """
    from lei_ren_part1_paper_joined_outer import _finite_difference
    with mp.workdps(source.precision):
        z=mp.mpf(str(Z))
        terminal=source.inner.evaluate_x(mp.e,z)
        tail=source.outer.tail
        def future(zz):
            return mp.mpf(tail.evaluate(float(zz),quadrature_order=source.outer.order)
                          ['energy_target_contribution_nominal'])
        future_jets=[]
        for h in (mp.mpf('1e-4'),mp.mpf('5e-5')):
            future_jets.append(_finite_difference(future,z,step=h))
        reference_factors=None
        if z==0:
            reference_factors=seeded['incoming'].get('continuous_incoming',{}).get('reference_linear_factors')
            if reference_factors is None:
                from lei_ren_part1_paper_axial_incoming import _stage_integrals
                reference=_stage_integrals(source.schedule,.5,order=seeded['incoming']['order'])
                reference_factors=dict(I_z_over_Z=2*reference['I_z'],
                    I_theta_z_times_1plusZ2_over_Z=mp.mpf('2.5')*reference['I_theta_z'])
        result=seeded_input_tangents(seeded,terminal,mp.nstr(z,source.precision),
            y_p=str(source.schedule.y_p),mu=str(source.schedule.mu),
            future_energy_Z=mp.nstr(future_jets[-1],source.precision),
            reference_linear_factors=reference_factors,
            precision=source.precision)
        result['future_energy_derivative_provenance']='Fourth-order float-backed Z difference of actual angular tail nominal target; analytic angular-coefficient derivative not implemented'
        result['future_energy_derivative_refinement']=[signed_log(x,source.precision) for x in future_jets]
        result['future_energy_derivative_absolute_change']=signed_log(abs(future_jets[1]-future_jets[0]),source.precision)
        result['all_input_Z_derivatives_analytic']=False
        return result


def run():
    from lei_ren_part1_paper_joined_outer import build_joined_field
    from lei_ren_part1_paper_continuous_axial_tangents import ContinuousAxialAlgebra
    folder=Path(__file__).parent
    source=build_joined_field()
    seeded=json.loads((folder/'lei_ren_part1_paper_seeded_shared_candidate.json').read_text())
    continuous=json.loads((folder/'lei_ren_part1_paper_continuous_axial_solve.json').read_text())
    print('evaluating actual corrected terminal moment jets',flush=True)
    result=actual_seeded_input_tangents(source,seeded,'.3')
    algebra=ContinuousAxialAlgebra(seeded,continuous)
    with mp.workdps(algebra.precision):
        jet=algebra.tangent([from_signed_log(x) for x in result['base_Z']],
                            from_signed_log(result['energy_target_Z']))
        result['coefficient_Z']=dict(a_Z=signed_log(jet['a_Z'],algebra.precision),
            c_Z=[signed_log(x,algebra.precision) for x in jet['c_Z']],
            linear_tangent_relative_replay=jet['linear_tangent_relative_replay'],
            energy_tangent_relative_replay=jet['energy_tangent_relative_replay'])
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:result['coefficient_Z'][k] for k in
        ('linear_tangent_relative_replay','energy_tangent_relative_replay')}),flush=True)
    return result


if __name__=='__main__':run()
