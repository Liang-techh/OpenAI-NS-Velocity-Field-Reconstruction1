"""Run original full first-micro m15 with complete prefix correction errors."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time

import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'src'),str(ROOT/'experiments/root_st073'),str(ROOT/'scripts')]
DIRECTORY=ROOT/'experiments/root_st073'
from openai_ns_reconstruction.original_core_numeric import PREFIX
from openai_ns_reconstruction.original_first_micro_full_m15 import OriginalFirstMicroFullM15,OriginalFirstMicroFullPressureRows
from run_original_first_micro_m15 import independent_S1,independent_original_source,product_errors

NAME=PREFIX+'compliant_current_original_first_micro_full_m15_numeric.json'


def run():
    start=time.perf_counter();source=independent_original_source();S1=independent_S1()
    f=OriginalFirstMicroFullM15(cells=8,precision=80)
    previous=DIRECTORY/(PREFIX+'compliant_current_original_first_micro_m15_numeric.json')
    raw=previous.read_bytes();old=json.loads(raw)
    assert old['source_family']==f.family and old['checks']['passed']
    for n,h in old['implementation_hashes'].items():assert hashlib.sha256((ROOT/n).read_bytes()).hexdigest()==h,n
    f.hashes[previous.name]=hashlib.sha256(raw).hexdigest()
    quarters={p['Z']:p for p in old['samples'] if p['coordinate']==dict(kind='exact_decimal_phase',value='.25') and p['cells_per_piece']==8}
    assert set(quarters)=={'-.3','0','.3'}
    coarse=OriginalFirstMicroFullM15(cells=4,rows=f.rows);products=product_errors(f)
    audit=mp.mp.clone();audit.dps=180
    samples=[];replay=cover=errors=endpoint=zeros=refinement=lower=0
    # Equal boundaries from different axial frames must remain distinct.
    a,b=dict(exact_left=[1,4]),dict(exact_left=[1,4])
    assert len(OriginalFirstMicroFullPressureRows.prefix_cover([a,b],[(a,1)]))==2
    for Z in ('-.3','0','.3'):
        for coordinate in ('.25','.5','.75','1'):
            p=f.evaluate(Z,coordinate);samples.append(p)
            assert all(p['contract'].values()) and len(p['m15_C0_C1_terms'])==13
            assert p['m15_C0_C1_terms'][:6]==p['incoming_q2_m15_terms']
            s,_,_=f.rows.micro._coordinate(coordinate,False)
            local,prefix,frames=f.rows.micro._covers(Z,s)
            cells=OriginalFirstMicroFullPressureRows.prefix_cover(local,prefix)
            metadata=p['finite_N_full_prefix_cover']
            assert metadata['uniform_prefix_source_cells']==len(cells);cover+=1
            for cell,width in prefix:
                assert any(v is cell for v in cells);cover+=1
            assert metadata['correction_errors_cover_full_original_prefix']
            for term in p['m15_C0_C1_terms'][6:]:
                for row in term['coefficient_C0_C1']:
                    ends=[audit.make_mpf(tuple(v)) for v in row['source_enclosure_exact_mpf_tuples']]
                    assert max(abs(v-audit.mpf(row['approximation'])) for v in ends)<=audit.mpf(row['absolute_error_upper']);errors+=1
            if coordinate=='.25':
                assert p['m15_C0_C1_terms']==quarters[Z]['m15_C0_C1_terms'];replay+=1
            if coordinate=='1':
                atom=p['m15_C0_C1_terms'][7]['log_scale_terms'][1]
                assert audit.mpf(atom['approximation'])==0
                assert all(v[1]==0 for v in atom['source_enclosure_exact_mpf_tuples']);endpoint+=1
                sigma,point,truth=f.sigma_atom(f.iv.mpf(1))
                assert point==1 and truth._mpi_==f.iv.mpf(1)._mpi_;endpoint+=1
                assert metadata['uniform_prefix_source_cells']==4*len(frames);cover+=1
            assert p['m15_C0_C1_terms'][8]['actual_function_binding']['endpoint_sign']==-1;lower+=1
            assert not p['limits']['resolved_finite_N_signed_centers']
        p=f.evaluate(Z,'.5',selected_sc_multiple=True);samples.append(p)
        assert p['exact_empty_increment']
        for term in p['m15_C0_C1_terms'][6:]:
            for row in term['coefficient_C0_C1']:
                assert all(v[1]==0 for v in row['source_enclosure_exact_mpf_tuples']);zeros+=1
        fine=f.evaluate(Z,'1');q=coarse.evaluate(Z,'1');samples.append(q)
        for k in (0,1):
            a=fine['m15_C0_C1_terms'][6]['coefficient_C0_C1'][k]
            b=q['m15_C0_C1_terms'][6]['coefficient_C0_C1'][k]
            assert abs(audit.mpf(a['approximation'])-audit.mpf(b['approximation']))<=audit.mpf(a['absolute_error_upper'])+audit.mpf(b['absolute_error_upper']);refinement+=1
        print('Original cumulative m15 reaches complete first_micro endpoint s=1: Z='+Z,flush=True)
    p=f.evaluate('.3','.875');samples.append(p)
    assert p['finite_N_full_prefix_cover']['uniform_prefix_source_cells']>=4
    # New full-prefix envelope preserves all predecessor cell objects even
    # when the endpoint itself is in the final source cell.
    s,_,_=f.rows.micro._coordinate('.875',False)
    local,prefix,_=f.rows.micro._covers('.3',s)
    assert len(f.rows.prefix_cover(local,prefix))>len(local);cover+=1
    rejected=0
    for args in [('1','1'),('.3','1.1'),('.3','0'),('nan','1')]:
        try:f.evaluate(*args)
        except ValueError:rejected+=1
        else:raise AssertionError('Invalid full first-micro query admitted')
    code=dict(f.code_hashes)
    for name in ('scripts/run_original_first_micro_m15.py','scripts/run_original_first_micro_full_m15.py'):
        code[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
    return dict(source_family=f.family,samples=samples,input_hashes=f.hashes,implementation_hashes=code,
        checks=dict(passed=True,independent_actual_Laurent_source=source,independent_actual_S1=S1,
            actual_full_first_micro_packets=len(samples),exact_accepted_quarter_history_replays=replay,
            uniform_full_prefix_and_axial_frame_conditions=cover,serialized_increment_errors=errors,
            exact_upper_sigma_one_conditions=endpoint,exact_empty_increment_rows=zeros,
            retained_lower_sigma_histories=lower,directed_cell_refinement_C0_C1=refinement,
            signed_quadratic_product_error_cases=products,rejected_invalid_queries=rejected),
        execution_seconds=round(time.perf_counter()-start,3),
        limits=dict(full_first_micro_leading_signed_m15_with_complete_factored_finite_N_errors=True,
            resolved_finite_N_signed_centers=False,useful_total_relative_accuracy=False,
            second_micro_m15_integrated=False,remaining_active_compact_signed_pressure=False,
            independent_preRp_signed_seam=False,global_five_moment_numeric_join=False,
            common_global_N=False,exact_heat_exterior=False,stress_cone=False,genuine_temporal_n_ge_2=False,full_corrected_NS=False))


def inspect(path):
    d=json.loads(path.read_text(encoding='utf-8'));assert d['checks']['passed']
    for n,h in d['input_hashes'].items():assert hashlib.sha256((DIRECTORY/n).read_bytes()).hexdigest()==h,n
    for n,h in d['implementation_hashes'].items():assert hashlib.sha256((ROOT/n).read_bytes()).hexdigest()==h,n
    return d


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--write-report',action='store_true');p.add_argument('--inspect-report',action='store_true')
    args=p.parse_args();path=DIRECTORY/NAME;d=inspect(path) if args.inspect_report else run()
    if args.write_report:path.write_text(json.dumps(d,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(checks=d['checks'],execution_seconds=d['execution_seconds']),indent=2),flush=True)
