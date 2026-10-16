"""Assemble original preRp m15 with the signed direct inner/lower inlet."""
import argparse
import copy
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
import time
from types import SimpleNamespace

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'src'),str(ROOT/'experiments/root_st073')]
DIRECTORY=ROOT/'experiments/root_st073'
from openai_ns_reconstruction.original_core_numeric import PREFIX
from openai_ns_reconstruction.original_preRp_signed_inner_m15 import (
    OriginalPreRpSignedInnerM15Numeric,INNER_TERMS,OUTER_TERMS)

NAME=PREFIX+'compliant_current_original_preRp_signed_inner_m15_numeric.json'


def partition_guards(f):
    """Malformed source partitions cannot silently duplicate or relabel history."""
    legacy=f.pre.evaluate('.3','compact','1')
    inner,_=f.inner_lower('.3')
    count=0
    cases=[]
    p=copy.deepcopy(legacy);p['m15_C0_C1_terms'].pop(1);cases.append((p,inner))
    p=copy.deepcopy(legacy);p['m15_C0_C1_terms'].insert(1,copy.deepcopy(p['m15_C0_C1_terms'][0]));cases.append((p,inner))
    for field,value in [('Z','0'),('q','1'),('source_family',{}),('exact_phase_tuple',[0,1,-4,1])]:
        p=copy.deepcopy(inner);p[field]=value;cases.append((legacy,p))
    p=copy.deepcopy(inner);p['m15_C0_C1_terms']=p['m15_C0_C1_terms'][:-1];cases.append((legacy,p))
    for bad_legacy,bad_inner in cases:
        candidate=object.__new__(OriginalPreRpSignedInnerM15Numeric)
        candidate.__dict__.update(f.__dict__)
        candidate._packets={}
        candidate.pre=SimpleNamespace(evaluate=lambda *args,p=bad_legacy:p)
        candidate.inner_lower=lambda Z,p=bad_inner:(p,True)
        try:candidate.evaluate('.3','compact','1')
        except ValueError:count+=1
        else:raise AssertionError('Malformed original source partition admitted')
    return count


def run():
    start=time.perf_counter();f=OriginalPreRpSignedInnerM15Numeric()
    samples=[];source=inner_count=outer_count=units=provenance=seams=zeros=aliases=0
    points=[('compact','1'),('compact','1.25'),('compact','1.4'),('compact','2'),
            ('quiet','0'),('quiet','.5'),('quiet','1')]
    for Z in ('-.3','0','.3'):
        inner,saved=f.inner_lower(Z)
        assert saved
        for chart,x in points:
            old=f.pre.evaluate(Z,chart,x)
            old_before=copy.deepcopy(old)
            p=f.evaluate(Z,chart,x);samples.append(p)
            assert all(p['contract'].values());source+=len(p['contract'])
            assert p['selected_integer']==old['selected_integer'];source+=1
            assert f.pre.evaluate(Z,chart,x)==old_before;source+=1
            terms=p['m15_C0_C1_terms']
            names=tuple(t['name'] for t in terms)
            assert names==INNER_TERMS+OUTER_TERMS+(('quiet_exact_leading',) if chart=='quiet' else ());source+=1
            assert terms[:len(INNER_TERMS)]==inner['m15_C0_C1_terms'];inner_count+=1
            assert terms[len(INNER_TERMS):]==old['m15_C0_C1_terms'][2:];outer_count+=1
            assert not {'inner','lower'}&set(names);inner_count+=1
            for a,b in zip(terms[:6],inner['m15_C0_C1_terms']):
                assert a['log_scale_terms']==b['log_scale_terms'];units+=1
                assert a['coefficient_C0_C1']==b['coefficient_C0_C1'];units+=1
            binding=p['signed_inner_lower_source']
            assert binding['native_q']=='2' and Fraction(binding['actual_Z'])==Fraction(Z);provenance+=1
            assert not binding['archived_FunctionRef_reuse_claimed'];provenance+=1
            assert binding['used_saved_numeric_packet'];provenance+=1
            assert binding['exact_phase_tuple']==inner['exact_phase_tuple'];provenance+=1
            assert not p['limits']['signed_active_and_compact_finite_N_values'];provenance+=1
            assert not p['limits']['independent_signed_numeric_seam_residual'];provenance+=1
            assert not p['limits']['global_five_moment_numeric_join'];provenance+=1
            assert 'source_partition_identity_reference' in p and 'source_partition_identity' not in p;provenance+=1
            if chart=='compact' and x=='1':
                for term in terms[6:]:
                    if term['name']!='active':
                        for row in term['coefficient_C0_C1']:
                            assert all(v[1]==0 for v in row['source_enclosure_exact_mpf_tuples']);zeros+=1
        at=f.evaluate(Z,'compact','2')['m15_C0_C1_terms']
        quiet=f.evaluate(Z,'quiet','0')['m15_C0_C1_terms']
        assert at==quiet[:-1];seams+=1
        assert f.evaluate(Z,'quiet','1')['m15_C0_C1_terms'][:6]==at[:6];seams+=1
        for alias in (Z+'0' if '.' in Z else '0.0',str(Fraction(Z)) if Fraction(Z).denominator==1 else Z):
            a,report=f.inner_lower(alias)
            assert report and a==inner;aliases+=1
        print('Signed inner/lower history consumed through actual compact/quiet toRp: Z='+Z,flush=True)
    guards=partition_guards(f)
    rejected=0
    for args in [('1','compact','1'),('nan','quiet','0'),('.3','compact','.9'),('.3','quiet','1.1'),('.3','other','0')]:
        try:f.evaluate(*args)
        except ValueError:rejected+=1
        else:raise AssertionError('Invalid original external coordinate admitted')
    code=dict(f.code_hashes);code[Path(__file__).relative_to(ROOT).as_posix()]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(source_family=f.family,samples=samples,input_hashes=f.hashes,implementation_hashes=code,
        checks=dict(passed=True,actual_signed_inner_preRp_packets=len(samples),source_partition_conditions=source,
            exact_signed_inner_cumulative_substitutions=inner_count,unchanged_active_compact_own_quiet_terms=outer_count,
            unmodified_inner_physical_units_and_coefficients=units,direct_source_and_open_scope_conditions=provenance,
            compact_quiet_and_Rp_history_seams=seams,exact_empty_compact_increment_zeros=zeros,
            exact_decimal_inner_source_aliases=aliases,rejected_malformed_source_partitions=guards,
            rejected_invalid_queries=rejected),execution_seconds=round(time.perf_counter()-start,3),
        limits=dict(signed_inner_lower_C0_C1_replaces_broad_components=True,
            active_compact_signed_finite_N_pressure=False,independent_preRp_signed_numeric_seam=False,
            resolved_signed_global_centers=False,relative_cancellation_accuracy=False,
            global_five_moment_numeric_join=False,common_global_N=False,exact_heat_exterior=False,
            stress_cone=False,genuine_temporal_n_ge_2=False,full_corrected_NS=False))


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
