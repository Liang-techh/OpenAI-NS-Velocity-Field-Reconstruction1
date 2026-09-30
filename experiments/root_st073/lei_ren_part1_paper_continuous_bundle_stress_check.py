"""Exploratory same-candidate five-moment/stress diagnostic after Rv.

Finite output is not admissible stress, PDE residual or closure evidence.
Dependencies must pass their own component checks before acceptance.
"""
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_joined_outer import build_joined_field
from lei_ren_part1_paper_continuous_incoming_outer import ContinuousIncomingOuterField
from lei_ren_part1_paper_continuous_moment_bundle import ContinuousMomentBundle
from lei_ren_part1_paper_axial_correction import signed_log


def run():
    folder=Path(__file__).parent
    print('building shared field for exploratory moment/stress bundle',flush=True)
    source=build_joined_field()
    seeded=json.loads((folder/'lei_ren_part1_paper_seeded_shared_candidate.json').read_text())
    atoms=json.loads((folder/'lei_ren_part1_paper_continuous_axial_solve.json').read_text())
    field=ContinuousIncomingOuterField(source,prepared={.3:(seeded,atoms)})
    p=field.outer;bundle=ContinuousMomentBundle(p);samples=[]
    with mp.workdps(field.precision):
        for label,radius in [('flatten',p.log_at(p.schedule.logR_v,'50')),('heat',p.log_at(p.schedule.logR_tail,'4'))]:
            print('evaluating '+label,flush=True)
            row=bundle.stress(radius,mp.mpf('.3'));stress=row['stress']
            if not all(mp.isfinite(v) for v in stress.values()):
                raise ArithmeticError('Bundle produced nonfinite stress')
            if set(row['moments'])!={'theta','z','theta_z','z_theta','p'}:
                raise ArithmeticError('Missing actual moment')
            shear=abs(stress['S_theta'])+abs(stress['S_z'])
            ratios={k:mp.nstr(abs(stress[k])/shear,40) if shear else None for k in ('T_theta','T_z')}
            samples.append(dict(label=label,moments={k:signed_log(v,60) for k,v in row['moments'].items()},
                moments_Z={k:signed_log(v,60) for k,v in row['moments_Z'].items()},
                P=signed_log(row['P'],60),P_Z=signed_log(row['P_Z'],60),
                stress={k:signed_log(v,60) for k,v in stress.items()},
                total_stress_over_sum_absolute_shears=ratios,
                all_five_moments_evaluated=True))
        report=dict(samples=samples,same_candidate_moment_pressure_velocity_inputs=True,
            analytic_velocity_radial_and_Z_jets=True,exploratory_only=True,
            provider_component_checks_reviewed=True,
            provider_component_scope="nominal anchored primitives and separate correction jets; no enclosures",
            tiny_stress_terms_separately_resolved=False,moment_closure_certified=False,
            stress_cone_certified=False,stress_remainder_decomposed=False,
            full_NS_residual_certified=False,finite_energy_certified=False,scale_recursion_certified=False)
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'stress_ratios':[{r['label']:r['total_stress_over_sum_absolute_shears']} for r in samples],'exploratory_only':True}),flush=True)
    return report


if __name__=='__main__':run()
