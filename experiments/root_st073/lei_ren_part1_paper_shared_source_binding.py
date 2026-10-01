"""Bind pressure, analytic core admission, and O.2 to one Md40 source."""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent


def run():
    names=('lei_ren_part1_paper_logarithmic_pressure_datum.json',
        'lei_ren_part1_paper_shared_core_majorant.json',
        'lei_ren_part1_paper_shared_core_tail_admission.json',
        'lei_ren_part1_paper_normalized_slow_turnoff.json')
    data={};hashes={}
    for name in names:
        payload=(HERE/name).read_bytes();data[name]=json.loads(payload)
        hashes[name]=hashlib.sha256(payload).hexdigest()
        for dependency,digest in data[name].get('input_hashes',{}).items():
            if hashlib.sha256((HERE/dependency).read_bytes()).hexdigest()!=digest:
                raise ValueError('binding dependency changed: '+dependency)
            hashes[dependency]=digest
    pressure,major,tail,outer=(data[n] for n in names)
    source=pressure['implicit_source_sha256'];datum=pressure['datum_enclosure_sha256']
    if major['analytic_core_family_sha256']!=tail['analytic_core_family_sha256']:
        raise ValueError('Fresh j/core family mismatch')
    if not major['source_j_eta_tol_relation_verified'] or not tail['source_j_eta_tol_relation_verified']:
        raise ValueError('Selected inlet tolerance missing')
    if any(r['implicit_source_sha256']!=source or r['datum_enclosure_sha256']!=datum for r in (major,tail)):
        raise ValueError('pressure source or enclosure fingerprint mismatch')
    definition=pressure['implicit_source_definition'];parameters=pressure['parameter_bounds']
    if definition['Md']!='40' or definition['logPstar']!='exp(Md)+11':
        raise ValueError('source outside the admitted O.2 family')
    con=next(r for r in outer['trials'] if r['Md']=='40')
    if not con['whole_O2_slow_turnoff_conditionally_certified']:
        raise ValueError('complete cutoff plus buffer certificate required')
    if not major['contraction_proved'] or not tail['target_met']:
        raise ValueError('analytic core admission missing')
    c=MPIntervalContext();c.dps=160
    with mp.workdps(200):
        get=lambda r,n:read_interval(c,r[n])
        if pressure['stage_count_total']!=14 or not pressure['complete_preheat_pressure_enclosed']:
            raise ValueError('complete backward pressure required')
        if endpoints(get(parameters['waiting_root'],'root_interval'))[0]<0:
            raise ValueError('nonnegative continuous waiting root required')
        for row in con['full_phase_interval_cells']:
            if endpoints(get(row,'b_squared_upper'))[1]>=4:
                raise ValueError('finite bracket monotonicity requires b_squared<4')
        # The selected matching radius is symbolic, so the enormous logC
        # is never evaluated: logRref=log110+10(logC+logPstar).
        logRref_lower=c.ln(110)+10*(2*get(major,'logLambda')+1000+get(major,'logPstar'))
        if endpoints(logRref_lower)[0]<10 or endpoints(get(major,'Gupper_in_logC_definition'))[0]<0:
            raise ValueError('symbolic matching radius does not meet O.2 threshold')
        result=dict(implicit_source_sha256=source,datum_enclosure_sha256=datum,
            Md='40',analytic_core_family_sha256=major['analytic_core_family_sha256'],
            required_j=major['required_j'],source_j_eta_tol_relation_verified=True,
            full_Section9_parameter_admission=False,source_pressure_condition_discharged_by_fourteen_atom_envelope=True,
            continuous_waiting_root_condition_discharged=True,
            core_analytic_contraction_and_scaled_tail_admitted=True,
            selected_radial_degree=tail['radial_degree'],
            logRref_definition='log(110)+10*(logC+logPstar)',logRref_lower=logRref_lower,
            O2_minimum_finite_strong_bracket=get(con,'finite_strong_bracket_lower'),
            pressure_core_outer_parameter_definitions_consistent=True,
            O2_certificate_bound_to_same_pressure_source=True,
            remaining_O2_inlet_condition='fresh core/transition/functional moment repair must realize exact reference identities',
            fresh_finite_core_generated=False,fresh_reference_inlet_built=False,
            admissible_stress_lift_constructed=False,whole_axis_cone_certified=False,
            actual_heat_pressure_restoration_completed=False,background_assembled=False,
            temporal_recursion=False,input_hashes=hashes)
    result['input_hashes'][Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(result),indent=2)+'\n',encoding='utf-8')
    print('Md40 pressure/core/outer source binding PASS; finite core and repaired inlet remain pending',flush=True)
    return result


if __name__=='__main__':run()
