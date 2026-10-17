"""New hydro source/C1 production; serialized inspection never invokes an owner."""
import argparse
import json
from pathlib import Path
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from openai_ns_reconstruction.original_hydro_source_point import PARENT,PARENT_SHA256,sha,digest,load
DIRECTORY=ROOT/'experiments/root_st073'
NAME='lei_ren_part1_paper_compliant_current_original_hydro_source_C1_numeric.json'
CODE=('src/openai_ns_reconstruction/original_hydro_source_kernel.py',
      'src/openai_ns_reconstruction/original_hydro_source_point.py',
      'scripts/check_original_hydro_source_kernel.py','scripts/check_original_hydro_transport.py',
      'scripts/run_original_hydro_source_C1.py')


def require(condition,message):
    if not condition:raise AssertionError(message)


def parent():
    require(sha(DIRECTORY/PARENT)==PARENT_SHA256,'Pressure parent bytes changed')
    p=load(DIRECTORY/PARENT);require(p['checks']['passed'] is True,'Passed parent required');return p


def check_packet(r,p):
    require(r['source_family']==p['source_family'] and r['parent_report_sha256']==PARENT_SHA256,'Source family/parent changed')
    require(r['source_guards']['passed'] is True and all(v['passed'] for v in r['source_guards']['AST'].values()),'Hydro source guards failed')
    require(r['inherited_m15_packet']==p['inherited_m15_packet'] and r['inherited_m15_packet_digest']==p['inherited_m15_packet_digest'],'m15 changed')
    require(r['inherited_mass_packet_digest']==p['inherited_mass_packet_digest'],'Mass changed')
    require(r['inherited_pressure_source_digest']==digest(p['source_jet']) and
            r['inherited_pressure_C1_digest']==digest(p['signed_source_integral']),'Accepted pressure source changed')
    require(len(r['hydro_source_rows'])==3 and len(r['signed_source_integral_rows'])==19,'Missing source/integral queries')
    width_proofs=0
    for row in r['hydro_source_rows']:
        require(row['ordinary_Z_orders']==[0,1,2] and len(row['source_rows'])==3,'Missing ordinary B_ZZ')
        require(row['no_formal_sectors_discarded'] is True and row['moment_difference_baseline_zero_by_shared_source'] is True,'Formal sector loss')
        require(row['full_source_scope'] is False and row['approximation_is_exact_source_value'] is False,'Partial source promoted')
        transport=row['finite_width_comparison'];require(transport['passed'] is True,'Finite-width transport failed')
        proofs=[v for group in transport['ordinary_error_proofs'].values() for v in group]+transport['rho_error_proof']
        require(all(v['passed'] and v['source_width_selected'] is False for v in proofs),'Width/source proof failed')
        width_proofs+=len(proofs)
        require(row['macro_nonlinear_remainder']['passed'] is True,'Macro nonlinear coupling omitted')
    for row,old in zip(r['signed_source_integral_rows'],p['signed_source_integral']['rows']):
        require((row['Z'],row['phase'],row['physical_scale'])==(old['Z'],old['phase'],old['physical_scale']),'Query/physical scale changed')
        require(row['ordinary_Z_orders']==[0,1] and len(row['hydro_main_component_C0_C1'])==2 and
                len(row['pressure_plus_hydro_main_component_C0_C1'])==2,'Missing C1 component')
        require(row['mass_Z_independent'] is True and row['weighted_raw47_remainder_retained'] is True and row['full_source_scope'] is False,'Scope promoted')
    require(all(v is False for v in r['limits'].values()),'Global completion claim introduced')
    require(r['scope']['full_source_scope'] is False and r['scope']['swirl_B_ZZ_recovered'] is False,'Swirl silently completed')
    require(r['contract']['accepted_ancestor_replayed'] is False and r['contract']['mass_recomputed'] is False and
            r['contract']['all_prior_errors_retained'] is True,'Preservation contract changed')
    return {'passed':True,'hydro_queries':3,'hydro_and_pressure_hydro_C0_C1_queries':19,
            'finite_width_log_error_proofs':width_proofs,'raw213_reduced183_and_pressure_preserved':True,
            'full_source_or_global_five_closure':False}


def pins(r):
    count=0
    for base,key in ((DIRECTORY,'input_hashes'),(ROOT,'implementation_hashes')):
        for name,expected in r[key].items():require(sha(base/name)==expected,'Changed pin: '+name);count+=1
    return count


def run():
    started=time.perf_counter()
    from check_original_hydro_source_kernel import run_fixture
    from check_original_hydro_transport import run_fixture as transport_fixture
    from openai_ns_reconstruction.original_hydro_source_point import OriginalHydroSourcePoint
    fixture=run_fixture();transport=transport_fixture()
    owner=OriginalHydroSourcePoint();r=owner.build()
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
