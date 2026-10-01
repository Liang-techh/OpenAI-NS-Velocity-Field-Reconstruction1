"""Replace coarse accepted mass enclosures with same-source refinements."""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_candidate_pressure_axis_jets import _interval_from_exact
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode


def run():
    base=Path(__file__).parent
    old_name='lei_ren_part1_paper_candidate_pressure_axis_jets.json'
    refined_name='lei_ren_part1_paper_candidate_pressure_mass_refinement.json'
    data=json.loads((base/old_name).read_text())
    refinement=json.loads((base/refined_name).read_text())
    if data['accepted_schedule_sha256']!=refinement['accepted_schedule_sha256']:
        raise AssertionError('Refined jets change source')
    for name,digest in refinement['input_hashes'].items():
        if hashlib.sha256((base/name).read_bytes()).hexdigest()!=digest:
            raise AssertionError('Refinement dependency changed')
    ctx=MPIntervalContext();ctx.dps=data['precision']
    with mp.workdps(ctx.dps+40):
        mass=_interval_from_exact(ctx,refinement['refined_beta2_mass'])
        mass0=_interval_from_exact(ctx,data['fixed_beta0_mass_upper_normalized'])
        scale=_interval_from_exact(ctx,data['Pstar_squared'])
        for row in data['Taylor_rows']:
            n=row['order']
            qcoef=_interval_from_exact(ctx,data['q_inverse_squared_coefficients'][n])
            fixed2=mass*qcoef
            fixed0=mass0 if n==0 else ctx.mpf(0)
            flatten=_interval_from_exact(ctx,row['flatten_normalized'])
            normalized=-(fixed2+fixed0+flatten)
            row.update(fixed_beta2_normalized=fixed2,
                       fixed_beta0_normalized=fixed0,
                       total_normalized_pressure_coefficient=normalized,
                       physical_pressure_coefficient=scale*normalized)
        data['fixed_beta2_mass_upper_normalized']=mass
        for stage,row in refinement['stages'].items():
            data['stage_mass_intervals'][stage]=row['refined_mass']
        # The inherited old nominal comparison was made against a wider box.
        data.pop('old_axis_pressure_consistency',None)
        data['input_hashes'].update({name:hashlib.sha256((base/name).read_bytes()).hexdigest()
                                    for name in (old_name,refined_name)})
        data['same_source_mass_refinement_applied']=True
        data['refinement_does_not_enclose_original_parameter_errors']=True
        output=base/'lei_ren_part1_paper_candidate_pressure_axis_jets_refined.json'
        output.write_text(json.dumps(encode(data),indent=2)+'\n',encoding='utf-8')
        print('Refined physical pressure jet orders:',len(data['Taylor_rows']),flush=True)
        return data


if __name__=='__main__':run()
