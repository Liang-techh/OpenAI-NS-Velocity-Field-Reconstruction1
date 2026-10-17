"""Produce new weighted source queries; inspect saved evidence without owners."""
import argparse
import json
from pathlib import Path
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from openai_ns_reconstruction.original_weighted_source_primitive import (
    PARENT,PARENT_SHA256,OriginalWeightedSourcePrimitive,sha,digest,load,
)
DIRECTORY=ROOT/'experiments/root_st073'
NAME='lei_ren_part1_paper_compliant_current_original_weighted_source_primitive_numeric.json'
CODE=(
    'src/openai_ns_reconstruction/original_weighted_pulse_moments.py',
    'src/openai_ns_reconstruction/original_weighted_angular_kernel.py',
    'src/openai_ns_reconstruction/original_weighted_remainder_kernel.py',
    'src/openai_ns_reconstruction/original_weighted_axial_modes.py',
    'src/openai_ns_reconstruction/original_weighted_source_primitive.py',
    'scripts/check_original_weighted_pulse_moments.py',
    'scripts/check_original_weighted_angular_kernel.py',
    'scripts/check_original_weighted_remainder_kernel.py',
    'scripts/run_original_weighted_source_primitive.py',
)


def require(value,message):
    if not value:raise AssertionError(message)


def parent():
    require(sha(DIRECTORY/PARENT)==PARENT_SHA256,'Saved parent bytes changed')
    p=load(DIRECTORY/PARENT);require(p['checks']['passed'] is True,'Passed parent required');return p


def check_packet(r,p):
    require(r['parent_report_sha256']==PARENT_SHA256 and r['source_family']==p['source_family'],'Source family/parent changed')
    require(r['source_guards']['passed'] and all(v['passed'] for v in r['source_guards']['AST'].values()),'Actual source binding failed')
    require(r['source_guards']['source_identity']['angular_control']=='1' and
            r['source_guards']['source_identity']['hbB_and_hbS_distinct_source_atoms'],'Actual angular control changed')
    require(r['inherited_m15_packet']==p['inherited_m15_packet'] and
            r['inherited_m15_packet_digest']==p['inherited_m15_packet_digest'],'m15 ledger witness changed')
    require(r['inherited_mass_packet_digest']==p['inherited_mass_packet_digest'],'Accepted M0 changed')
    require(r['inherited_all_sector_source_digest']==digest(p['all_sector_source_rows']) and
            r['inherited_main_integrals_digest']==digest(p['signed_source_integral_rows']),'Accepted source/main packets changed')
    require(len(r['angular_source_rows'])==len(r['axial_mode_evidence'])==3,'Missing modal queries')
    amplitude_proofs=0
    for row in r['axial_mode_evidence']:
        require(row['passed'] and row['all_original_sectors_retained'] and row['positive_swirl_amplitude_materialized'] is False,'Axial source hydration failed')
        require([(m['sector'],m['mode_j'],m['rate']) for m in row['modes']]==
                [('pressure',0,2),('hydro',0,2),('hydro',1,1),('swirl',0,3),('swirl',1,2),('swirl',2,1)],'Missing axial mode or incorrect radial rate')
        require(len(row['numerical_mode_regressions'])==6 and
                all(v['passed'] and v['numerical_regression_only_not_source_equality_proof'] for v in row['numerical_mode_regressions']),
                'Saved numerical modal reconstruction regression failed')
        for mode in row['swirl_mode_proofs']:
            require(mode['macro_proof']['passed'] and mode['scaled_atom_proof']['passed'],'Per-mode macro/amplitude proof failed')
            require(all(v['passed'] and v['amplitude_exponentiated_or_selected'] is False for v in mode['scaled_atom_proof']['orders']),'Amplitude proxy not proved')
            amplitude_proofs+=len(mode['scaled_atom_proof']['orders'])
    angular_proofs=0
    for row in r['angular_source_rows']:
        require(row['passed'] and row['ordinary_Z_orders']==[0,1,2] and row['no_Laurent_sector_discarded'],'Angular modal source failed')
        require([m['mode_j'] for m in row['modes']]==[0,1,2] and
                [m['radial_rate'] for m in row['modes']]==[1,0,-1],'Angular radial rates changed')
        require(-1 in row['modes'][0]['Laurent_numerator_degrees'] and
                -1 in row['modes'][2]['Laurent_numerator_degrees'],'Gauge Laurent sector lost')
        for group in row['new_H_Kpsi_prefix_width_proofs'].values():
            require(len(group)==4 and all(v['passed'] and v['source_width_selected'] is False for v in group),'New prefix proof failed')
            angular_proofs+=len(group)
    w=r['weighted_pulse_moments'];require(w['passed'] and len(w['rows'])==14 and w['powers']==[1,2],'New weighted moment packet failed')
    require(w['controls']['no_M0_mass_row'] and not w['controls']['old_mass_producer_called'],'M0 replay introduced')
    require(len(r['weighted_remainder_queries'])==19,'Missing weighted queries')
    width_count=0
    for row,old in zip(r['weighted_remainder_queries'],p['signed_source_integral_rows']):
        require((row['Z'],row['phase'])==(old['Z'],old['phase']),'Saved query order changed')
        require(row['ordinary_Z_orders']==[0,1,2] and len(row['raw47_signed_C0_C1'])==2,'Ordinary source Z2/C1 lost')
        require(row['raw46_main_inherited_digest']==digest(old['all_sector_main_component_C0_C1']),'Main component changed')
        for name in ('signed_R_ordinary_rows','factored_R0_rows','factored_R1_rows','relative_swirl_R0_rows','relative_swirl_R1_rows'):
            require(len(row[name])==3 and all('serialized_approximation_absolute_error_upper_exact_mpf_tuple' in v for v in row[name]),'Approximation/truth distinction missing')
        cert=row['finite_width_certificate']
        require(cert['passed'] and cert['h_selected'] is False and cert['hb_materialized'] is False and
                cert['positive_cap_selected_as_value'] is False,'Width selected/materialized')
        for group in ('width_log_proofs','hR1_log_proof_per_order','E3_log_proof_per_order'):
            require(len(cert[group])==3 and all(v['passed'] for v in cert[group]),'Finite-width log proof failed');width_count+=len(cert[group])
        require(row['old_raw47_error_budget_replaced'] is False and row['complete_functional_closure'] is False,'Premature ledger/closure claim')
        atom=row['swirl_absolute_atom'];require(atom['actual_amplitude_positive'] and atom['selected_value'] is None,'Positive swirl source lost')
        if row['phase'] in ('0','0.0'):
            require(all(v['approximation']=='0.0' or v['approximation']=='0' for v in row['signed_R_ordinary_rows']),'Phase-zero exact branch failed')
    require(r['contract']['accepted_ancestor_replayed'] is False and r['contract']['M0_recomputed'] is False and
            r['contract']['raw213_reduced183_and_all_errors_retained'] is True,'Preservation contract changed')
    require(all(v is False for v in r['limits'].values()),'Global completion claim introduced')
    return {'passed':True,'weighted_moment_rows':14,'angular_queries':3,'axial_modal_queries':3,
            'ordinary_R_and_raw47_C0_C1_queries':19,'new_H_Kpsi_width_proofs':angular_proofs,
            'finite_h_log_proofs':width_count,'per_axial_mode_amplitude_log_proofs':amplitude_proofs,
            'numerical_modal_reconstruction_regressions':18,
            'old_raw47_and_all_errors_retained':True,
            'global_functional_primitive_or_five_closure':False}


def pins(r):
    count=0
    for base,key in ((DIRECTORY,'input_hashes'),(ROOT,'implementation_hashes')):
        for name,expected in r[key].items():require(sha(base/name)==expected,'Changed pin: '+name);count+=1
    return count


def run():
    started=time.perf_counter()
    from check_original_weighted_pulse_moments import run_fixture as moment_fixture
    from check_original_weighted_angular_kernel import run_fixture as angular_fixture
    from check_original_weighted_remainder_kernel import run_fixture as remainder_fixture
    fixtures={'weighted_moments':moment_fixture(),'angular_kernel':angular_fixture(),'remainder_kernel':remainder_fixture()}
    owner=OriginalWeightedSourcePrimitive();r=owner.build()
    code=dict(owner.implementation_hashes);code.update({name:sha(ROOT/name) for name in CODE})
    r.update(input_hashes=owner.input_hashes,implementation_hashes=code)
    r['checks']={'passed':True,'packet_invariants':check_packet(r,owner.parent),'new_kernel_fixtures':fixtures,
                 'no_accepted_ancestor_replay':True}
    r['execution_seconds']=round(time.perf_counter()-started,3)
    return r


def inspect(path):
    r=load(path);require(r['checks']['passed'],'Saved weighted report failed')
    for fixture in r['checks']['new_kernel_fixtures'].values():require(fixture['passed'],'Saved new fixture failed')
    return {'passed':True,'pins':pins(r),'packet_invariants':check_packet(r,parent()),
            'serialized_only':True,'owner_fixture_or_ancestor_run':False}


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);g=ap.add_mutually_exclusive_group(required=True)
    g.add_argument('--write-report',action='store_true');g.add_argument('--inspect-report',action='store_true')
    args=ap.parse_args();path=DIRECTORY/NAME
    if args.write_report:
        r=run();tmp=path.with_suffix('.json.tmp');tmp.write_text(json.dumps(r,indent=2)+'\n',encoding='utf-8',newline='\n');tmp.replace(path)
    else:print(json.dumps(inspect(path),indent=2))
