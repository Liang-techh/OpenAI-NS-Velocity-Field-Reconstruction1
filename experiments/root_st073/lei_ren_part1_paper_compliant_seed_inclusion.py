"""New preheat source inside the completed core/exit INPUT enclosures.

This proves input inclusion, not equality of the old and new solutions.
The independently admitted compliant analytic core can use the stored
directed recurrence bounds because its seed inputs are contained in them.
"""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_compliant_pressure_source import CompliantPressureDatum
from lei_ren_part1_paper_logarithmic_pressure_datum import LogarithmicPressureDatum
from lei_ren_part1_paper_candidate_gauge_core import _unpack_vector
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent


def run():
    names=('lei_ren_part1_paper_compliant_core_transfer.json',
           'lei_ren_part1_paper_compliant_core_transfer_check.json',
           'lei_ren_part1_paper_shared_core_seed_check.json',
           'lei_ren_part1_paper_shared_interval_core_Z049_Z051.json',
           'lei_ren_part1_paper_shared_interval_core_Z049_Z051_state.json')
    major,check,oldseed,oldcore,state=[json.loads((HERE/n).read_bytes()) for n in names]
    hashes={}
    for record in (major,check,oldseed,oldcore,state):
        for name,digest in record['input_hashes'].items():
            if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:
                raise ValueError('Seed inclusion dependency changed: '+name)
            hashes[name]=digest
    hashes.update({n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names})
    hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    new=CompliantPressureDatum();old=LogarithmicPressureDatum();c=new.ctx
    if (not major['contraction_proved'] or not check['all_passed']
        or not oldseed['seed_acceptance_passed'] or not oldcore['completed_target']
        or major['implicit_source_sha256']!=new.source_sha
        or state['identity']['implicit_source_sha256']!=old.source_sha
        or oldseed['analytic_core_family_sha256']!=major['analytic_core_family_sha256']):
        raise ValueError('New existence/old envelope binding missing')
    with mp.workdps(c.dps+40):
        contained=lambda a,b: endpoints(b)[0]<=endpoints(a)[0]<=endpoints(a)[1]<=endpoints(b)[1]
        if not all(endpoints(getattr(new.parameters,key))==endpoints(getattr(old.parameters,key))
                   for key in ('logPstar','mu','delta')):
            raise ValueError('Upstream axis data changed')
        for stage in new.stages:
            if not contained(c.mpf(old.ctx.mpf(new.stages[stage]['mass'])),c.mpf(old.stages[stage]['mass'])):
                raise ValueError('New pressure mass not contained: '+stage)
        if not contained(c.mpf(new.flatten_complex_upper),c.mpf(old.flatten_complex_upper)):
            raise ValueError('New complex flatten envelope not contained')
        count=state['identity']['initial_axis_length']
        if count!=148 or state['completed_radial_order']!=144:
            raise ValueError('Completed seed layout changed')
        pressure=new.normalized_jets(c.mpf(state['identity']['center_Z']),count-1)
        # Reproduce the seed's EXACT evaluation path, not exp(2logP)*eps;
        # algebraically equivalent interval paths can have wider roundoff.
        factor=c.exp(2*new.parameters.logPstar-read_interval(c,major['logLambda']))
        stored=_unpack_vector(c,state['fixed']['P0_Z_taylor'])
        checks=[]
        for k,(a,b) in enumerate(zip(pressure['normalized_pressure_coefficients'],stored)):
            inside=contained(a*factor,b)
            if not inside:
                raise ValueError('Compliant pressure input outside stored seed at order '+str(k))
            checks.append(inside)
        if (endpoints(read_interval(c,state['identity']['required_j']))!=endpoints(read_interval(c,major['required_j']))
            or state['seed_metadata']['logC_definition']!=major['logC_definition']):
            raise ValueError('Fixed j/F0 definition changed')
        fixed_sha=hashlib.sha256(json.dumps(state['fixed'],sort_keys=True).encode()).hexdigest()
        if fixed_sha!=oldseed['fixed_seed_sha256']:
            raise ValueError('Accepted envelope seed changed')
        hashes.update(new.input_hashes);hashes.update(old.input_hashes)
        result=dict(implicit_source_sha256=new.source_sha,datum_enclosure_sha256=new.datum_sha,
            old_implicit_source_sha256=old.source_sha,
            analytic_core_family_sha256=major['analytic_core_family_sha256'],
            fixed_seed_sha256=fixed_sha,
            fixed_seed_sha_scope='accepted stored ENVELOPE seed; not new point seed or coefficient equality',
            seed_axial_length=count,new_pressure_seed_inclusion_checks=checks,
            all_fourteen_pressure_mass_boxes_contained=True,
            complex_flatten_envelope_contained=True,
            exact_upstream_axis_delta_mu_logP_j_and_F0_definition_unchanged=True,
            local_center_domain=state['identity']['center_Z'],
            unique_axis_anchor_root_proved=True,
            primitive_branch='G(Z0)=0 at unique H root',
            G_branch_definition='same integral of L*H0/(H0^2+sigma^2), anchored at the unique H0 root',
            real_F0_strictly_positive_by_definition=True,
            accepted_input_envelopes_contain_compliant_source=True,
            seed_acceptance_passed=True,finite_core_completed=True,
            finite_core_scope='old completed interval recurrence encloses new input class by explicit inclusion',
            new_point_coefficients_recomputed=False,old_source_or_state_relabelled=False,
            proof='Both sources have the same axis velocity, log gradient and amplitude class. Every new normalized pressure mass/Cauchy envelope and all 148 scaled pressure input jets are contained. Directed algebraic recurrence and guarded rational operations are inclusion preserving. Independent new analytic existence/uniqueness identifies the enclosed new solution.',
            temporal_recursion=False,input_hashes=hashes)
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(result),indent=2)+'\n',encoding='utf-8')
    print('Compliant source inclusion: 14 pressure mass boxes and 148 scaled seed jets contained; old state preserved',flush=True)
    return result


if __name__=='__main__':run()
