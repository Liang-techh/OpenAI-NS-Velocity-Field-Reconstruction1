"""Produce new swirl/all-sector source C1; inspect serialized evidence only."""
import argparse
import json
from pathlib import Path
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from openai_ns_reconstruction.original_swirl_source_point import PARENT,PARENT_SHA256,ANCHORED,sha,digest,load
DIRECTORY=ROOT/'experiments/root_st073'
NAME='lei_ren_part1_paper_compliant_current_original_swirl_source_C1_numeric.json'
CODE=('src/openai_ns_reconstruction/original_swirl_source_kernel.py',
      'src/openai_ns_reconstruction/original_swirl_source_point.py',
      'scripts/check_original_swirl_source_kernel.py','scripts/check_original_swirl_transport.py',
      'scripts/run_original_swirl_source_C1.py')


def require(condition,message):
    if not condition:raise AssertionError(message)


def parent():
    require(sha(DIRECTORY/PARENT)==PARENT_SHA256,'Hydro parent bytes changed')
    p=load(DIRECTORY/PARENT);require(p['checks']['passed'] is True,'Passed parent required');return p


def check_packet(r,p):
    require(r['source_family']==p['source_family'] and r['parent_report_sha256']==PARENT_SHA256,'Family/parent changed')
    require(r['source_guards']['passed'] is True and all(v['passed'] for v in r['source_guards']['AST'].values()),'Specific swirl source guard failed')
    require(r['inherited_m15_packet']==p['inherited_m15_packet'] and r['inherited_m15_packet_digest']==p['inherited_m15_packet_digest'],'m15 witness changed')
    require(r['inherited_mass_packet_digest']==p['inherited_mass_packet_digest'],'Saved mass changed')
    require(r['inherited_hydro_source_digest']==digest(p['hydro_source_rows']) and
            r['inherited_hydro_integrals_digest']==digest(p['signed_source_integral_rows']),'Saved hydro changed')
    require(r['inherited_pressure_source_digest']==p['inherited_pressure_source_digest'],'Saved pressure changed')
    require(len(r['swirl_source_rows'])==3 and len(r['all_sector_source_rows'])==3 and len(r['signed_source_integral_rows'])==19,'Missing new queries')
    width_proofs=0;scale_proofs=0
    anchored=load(DIRECTORY/ANCHORED)['anchored_amplitude_packets']
    for row in r['swirl_source_rows']:
        require(row['ordinary_Z_orders']==[0,1,2] and len(row['relative_source_rows'])==3 and len(row['absolute_source_rows'])==3,'Missing B_ZZ')
        require(len(row['relative_amplitude_ordinary_Z0_Z3'])==4 and len(row['mode_delta_pieces'])==6,'Missing relative dressing or mode/delta sector')
        require(row['actual_swirl_set_to_zero'] is False and row['external_amplitude_differentiated_twice'] is False,'Original positive source lost')
        atom=row['absolute_source_atom'];require(atom['selected_value'] is None and atom['actual_amplitude_positive'] is True
            and atom['relative_derivatives_already_dressed'] is True and atom['Gbar_selected_as_G'] is False,'Amplitude selected or differentiated twice')
        original=anchored[row['Z']]
        require(atom['anchored_packet_digest']==digest(original) and atom['anchored_G_source_enclosure']==original['G']
            and atom['selected_logCstar_source_enclosure']==original['selected_logCstar']
            and row['anchored_gradient_source_rows']==original['G_gradient_taylor_coefficients'][:3],
            'Saved actual anchored G/amplitude/gradient changed')
        proof=row['amplitude_log_product_proof'];require(proof['passed'] is True and proof['zero_is_controlled_numerical_approximation_only'] is True
            and proof['actual_positive_amplitude_is_zero'] is False and proof['underflow_used_to_define_field'] is False,'Numerical proxy promoted to field')
        require(all(v['passed'] and v['amplitude_exponentiated_or_selected'] is False for v in proof['orders']),'Amplitude log product failed')
        scale_proofs+=len(proof['orders'])
        tr=row['finite_width_comparison'];require(tr['passed'] is True,'Comparison transport failed')
        for group in tr['ordinary_error_proofs'].values():
            require(all(v['passed'] and v['source_width_selected'] is False for v in group),'Width source selected');width_proofs+=len(group)
        require(row['macro_nonlinear_remainder']['passed'] is True,'Missing expell derivative coupling')
    for row in r['all_sector_source_rows']:
        require(row['components']==['pressure','hydro','swirl'] and row['all_sectors_enclosed_at_query'] is True
            and row['weighted_primitive_complete'] is False,'All-sector source scope changed')
    for row,old in zip(r['signed_source_integral_rows'],p['signed_source_integral_rows']):
        require((row['Z'],row['phase'],row['physical_scale'])==(old['Z'],old['phase'],old['physical_scale']),'Query/physical scale changed')
        require(row['ordinary_Z_orders']==[0,1] and len(row['all_sector_main_component_C0_C1'])==2 and
                len(row['swirl_main_component_relative_C0_C1'])==2,'Missing ordinary C1')
        require(row['mass_Z_independent'] is True and row['raw47_weighted_remainder_retained'] is True and row['complete_weighted_primitive'] is False,'Weighted primitive scope promoted')
    require(r['contract']['accepted_ancestor_replayed'] is False and r['contract']['mass_recomputed'] is False
        and r['contract']['raw213_reduced183_and_all_errors_retained'] is True,'Preservation changed')
    require(all(v is False for v in r['limits'].values()),'Global completion claim introduced')
    return {'passed':True,'swirl_and_all_sector_queries':3,'all_sector_main_C0_C1_queries':19,
            'finite_width_log_error_proofs':width_proofs,'amplitude_log_product_proofs':scale_proofs,
            'positive_absolute_swirl_source_retained':True,'raw213_reduced183_pressure_hydro_and_mass_preserved':True,
            'complete_weighted_primitive_or_global_five_closure':False}


def pins(r):
    count=0
    for base,key in ((DIRECTORY,'input_hashes'),(ROOT,'implementation_hashes')):
        for name,expected in r[key].items():require(sha(base/name)==expected,'Changed pin: '+name);count+=1
    return count


def run():
    started=time.perf_counter()
    from check_original_swirl_source_kernel import run_fixture
    from check_original_swirl_transport import run_fixture as transport_fixture
    from openai_ns_reconstruction.original_swirl_source_point import OriginalSwirlSourcePoint
    fixture=run_fixture();transport=transport_fixture()
    owner=OriginalSwirlSourcePoint();r=owner.build()
    code=dict(owner.implementation_hashes);code.update({name:sha(ROOT/name) for name in CODE})
    r.update(input_hashes=owner.input_hashes,implementation_hashes=code)
    r['checks']={'passed':True,'packet_invariants':check_packet(r,owner.parent),
                 'independent_kernel_fixture':fixture,'independent_transport_fixture':transport,
                 'no_accepted_ancestor_replay':True}
    r['execution_seconds']=round(time.perf_counter()-started,3)
    return r


def inspect(path):
    r=load(path);require(r['checks']['passed'] is True,'Saved report failed')
    for name in ('independent_kernel_fixture','independent_transport_fixture'):
        require(r['checks'][name]['passed'] and all(v['passed'] for v in r['checks'][name]['controls']),'Saved new fixture failed')
    return {'passed':True,'pins':pins(r),'packet_invariants':check_packet(r,parent()),
            'serialized_only':True,'owner_fixture_or_ancestor_run':False}


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);group=ap.add_mutually_exclusive_group(required=True)
    group.add_argument('--write-report',action='store_true');group.add_argument('--inspect-report',action='store_true')
    args=ap.parse_args();path=DIRECTORY/NAME
    if args.write_report:
        r=run();tmp=path.with_suffix('.json.tmp');tmp.write_text(json.dumps(r,indent=2)+'\n',encoding='utf-8',newline='\n');tmp.replace(path)
    else:print(json.dumps(inspect(path),indent=2))
