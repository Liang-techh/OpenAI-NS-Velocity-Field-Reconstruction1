"""Produce the same-source quiet contraction, or inspect saved bytes only."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from openai_ns_reconstruction.original_quiet_scalar_m15 import (
    PARENT,PARENT_SHA256,AXIAL,QUIET,PRERP,FULL,OriginalQuietScalarM15,
    REMAINDER_OPERATION,contract_quiet_pair,materialize_quiet179_view,sha,digest,load,
    mp,MPIntervalContext,quiet_scalar_difference,
)
DIRECTORY=ROOT/'experiments/root_st073'
NAME='lei_ren_part1_paper_compliant_current_original_quiet_scalar_m15_numeric.json'
CODE=('src/openai_ns_reconstruction/original_quiet_scalar_kernel.py',
      'src/openai_ns_reconstruction/original_quiet_scalar_source.py',
      'src/openai_ns_reconstruction/original_quiet_scalar_m15.py',
      'scripts/run_original_quiet_scalar_m15.py')


def require(value,message):
    if not value:raise AssertionError(message)


def rational_tuple(v):
    s,m,e,b=v
    require(s in (0,1) and m>=0 and m.bit_length()==b and abs(e)<=100000,'Bounded exact binary tuple required')
    return (-1 if s else 1)*F(m)*F(2)**e


def endpoints(row):
    raw=row.get('source_enclosure_exact_mpf_tuples')
    if raw is None:raw=[row[k+'_exact_mpf_tuple'] for k in ('lower','upper')]
    return [rational_tuple(v) for v in raw]


def controls():
    iv=MPIntervalContext();iv.dps=100;rows=[]
    def check(name,value):require(value,name);rows.append({'name':name,'passed':True})
    for i,(a,b,L) in enumerate((((1,2),(3,4),(2,3)),((-2,-1),(3,4),(2,3)),
        ((-2,1),(-1,3),(2,3)),((1,1),(1,1),(2,2)),((3,4),(1,2),(2,3)))):
        result,receipt=quiet_scalar_difference(iv,iv.mpf(a),iv.mpf(b),iv.mpf(L),source_bound=True)
        squares=lambda x:(F(0) if x[0]<=0<=x[1] else min(F(x[0])**2,F(x[1])**2),max(F(x[0])**2,F(x[1])**2))
        aa,bb=squares(a),squares(b);d=(aa[0]-bb[1],aa[1]-bb[0])
        exact=[-F(l)*q/4 for l in L for q in d];lo,hi=map(rational_tuple,result[0]._mpi_)
        check('independent_fraction_range_'+str(i),lo<=min(exact)<=max(exact)<=hi)
        check('ordinary_Z1_exact_zero_'+str(i),result[1]._mpi_==iv.mpf(0)._mpi_)
        check('overlap_not_selected_as_equality_'+str(i),receipt['enclosure_overlap_used_as_equality_or_source_value'] is False)
    for flag in (False,1):
        try:quiet_scalar_difference(iv,iv.mpf(1),iv.mpf(1),iv.mpf(2),source_bound=flag)
        except ValueError:check('reject_false_source_'+str(flag),True)
        else:raise AssertionError('Source proof required')
    try:quiet_scalar_difference(iv,iv.mpf([1,2]),iv.mpf([1,2]),iv.mpf(2),source_bound=True,gap=iv.mpf(0))
    except ValueError:check('unproved_zero_gap_from_overlap_rejected',True)
    else:raise AssertionError('Gap proof required')
    try:quiet_scalar_difference(iv,iv.mpf(1),iv.mpf(2),iv.mpf(0),source_bound=True)
    except ValueError:check('positive_length_required',True)
    else:raise AssertionError('Positive length required')
    correlated,_=quiet_scalar_difference(iv,iv.mpf([1,2]),iv.mpf(['1.5','2.5']),iv.mpf(2),
        source_bound=True,gap=iv.mpf(0),gap_source_bound=True)
    check('source_identity_not_numerical_width_controls_gap',correlated[0]._mpi_==iv.mpf(0)._mpi_)
    # Physical radius and decaying mean must be paired at the same endpoint.
    Rax=F(3);transport=F(17);m=F(2,7);R2=2*transport*Rax;m2=m/(2*transport)
    check('transported_mean_radius_pair',4*R2*m2==4*Rax*m)
    check('wrong_untransported_2Rc_mean_detected',4*R2*m!=4*Rax*m)
    import sympy as s
    y=s.Symbol('y',positive=True);Y=s.exp(40);q=Y-y
    mapped=-s.log(1-q/Y)/40
    check('exact_sigma_argument_substitution',s.simplify(s.expand_log(mapped,force=True)-(1-s.log(y)/40))==0)
    check('exact_kernel_exponent_substitution',s.simplify(s.exp(-q)-s.exp(y-Y))==0)
    check('signed_Jacobian_and_endpoints',s.diff(q,y)==-1 and q.subs(y,Y)==0 and q.subs(y,1)==Y-1)
    logC=s.Symbol('logC');logP=Y+11
    logRax=s.log(110)+10*logC+10*logP+1+(Y-1)
    logRc=s.log(110)+10*logC+11*Y+124
    check('actual_scalar_recipe_Rc_minus_Raxial',s.simplify(logRc-logRax)==14)
    Tw,lp,Fsource=s.symbols('Tw lp F',real=True)
    K_ratio=s.exp(Tw-2)*s.exp(lp)*4*Fsource*s.exp(-Tw-lp-12)
    check('actual_Pstar_kappa_Rp_cancellation',s.simplify(K_ratio-4*s.exp(-14)*Fsource)==0)
    return {'passed':True,'control_count':len(rows),'controls':rows,'accepted_fixture_or_quadrature_replayed':False,
        'original_sigma_or_source_integrals_numerically_executed':False}


def parent():
    require(sha(DIRECTORY/PARENT)==PARENT_SHA256,'Pinned joint parent changed')
    p=load(DIRECTORY/PARENT);require(p['checks']['passed'] is True,'Accepted joint parent required');return p


def packet_check(r,p):
    require(r['source_family']==p['source_family'] and r['parent_report_sha256']==PARENT_SHA256,'Family/parent changed')
    proof=r['quiet_scalar_source'];diag=r['independent_scalar_diagnostic'];rem=r['canonical_nonquiet_remainder']
    require(proof['passed'] is True and all(v['passed'] is True for v in proof['AST'].values()),'Quiet source identity unproved')
    require(proof['gaps']==[] and proof['source_executed'] is False,'Unresolved source pin/identity gap')
    require(all(proof['source_mean_equality'][k] is True for k in ('passed','mstar_equals_F','Ktrue_equals_Kcan','R_axial_pairing')),
        'No cross-producer integral/physical-radius identity')
    require(proof['source_change_of_variables']['passed'] is True and all(proof['physical_scale_guards'].values()),
        'Wrong variable change or physical scale')
    canonical=proof['canonical_quiet']
    require(canonical['quiet_pair']==[31215,0] and canonical['preRp_signed_pair']==[31221,31220], 'Wrong quiet source pair')
    require(canonical['sum_chain']=={'m15':{'node':31221,'children':[31219,31215],'active_plus_compact':31219,
        'quiet':31215,'active_plus_compact_children':[30707,31155]},'C1':{'node':31220,'children':[30708,31156]}},
        'Canonical nonquiet partition changed')
    require(rem['passed'] is True and all(v['passed'] is True for v in rem['AST'].values()),'Canonical nonquiet numerical source unproved')
    require(rem['component_names']==list(rem['components'])==['inner','lower','active','compact'],'Wrong remainder source scope')
    require(rem['canonical_physical_C1_pairs']=={'inner':[53161,53162],'lower':[36592,36593],
        'active':[30707,30708],'compact':[31155,31156]},'Wrong canonical source pairs')
    full=load(DIRECTORY/FULL)['samples']['whole_Z'];pre=load(DIRECTORY/PRERP)['samples']['whole_Z']
    require(rem['enclosure_domain']==full['Z'],'Canonical enclosure domain changed')
    for name in ('inner','lower'):
        require(rem['components'][name]==full['pressure_component_enclosures']['m15'][name]['C1_enclosures'],
            'Saved inner/lower component changed')
    for name in ('active','compact'):
        for k in (0,1):
            source=[endpoints(v['piece_enclosures']['m15'][name][k]) for v in pre['cell_piece_enclosures']]
            wanted=(min(v[0] for v in source),max(v[1] for v in source));lo,hi=endpoints(rem['components'][name][k])
            require(lo<=wanted[0]<=wanted[1]<=hi and hi-lo<=(wanted[1]-wanted[0])*(1+F(1,10**60)),
                'Saved active/compact component hull changed')
    for k,row in enumerate(rem['coefficient_C0_C1']):
        lo,hi=endpoints(row);wanted=(-sum(endpoints(v[k])[1] for v in rem['components'].values()),
            -sum(endpoints(v[k])[0] for v in rem['components'].values()))
        require(lo<=wanted[0]<=wanted[1]<=hi,'Directed remainder bound misses component sum')
        require(hi-lo<=(wanted[1]-wanted[0])*(1+F(1,10**60)),'Unexplained remainder bound inflation')
    require(rem['source_rows_are_enclosures_not_function_values'] and rem['source_errors_retained'],'Remainder source error discarded')
    axial=load(DIRECTORY/AXIAL);means=[v['leading_scalars']['mean'] for v in axial['pressure_row_samples'] if v['native_coordinate']=='1']
    quiet=load(DIRECTORY/QUIET)['actual_quiet_source']
    require(means and all(v==diag['source_mean_summary'] for v in means),'Different axial terminal mean source')
    require(diag['canonical_log_F']==quiet['kernel']['log_F'] and diag['physical_log_length']==quiet['original_quiet_length'],'Different canonical scalar or logarithmic length')
    zero=[[0,0,0,0],[0,0,0,0]]
    require(all(v['source_enclosure_exact_mpf_tuples']==zero for v in diag['same_source_correlated_difference_C0_C1']),'True source difference not exactly zero')
    require(diag['independent_saved_enclosure_difference_C0_C1'][1]['source_enclosure_exact_mpf_tuples']==zero,'Ordinary Z1 not zero')
    require(diag['mean_approximation_or_bound_selected_as_actual_value'] is False and diag['interval_overlap_used_as_source_identity'] is False,'Numerical box used as source theorem')
    require(diag['same_source_integral_theorem_used'] is True and diag['accepted_quadrature_or_owner_replayed'] is False,'Source identity or replay contract changed')
    require(diag['normalized_scalar_exponential_only'] is True and diag['physical_radius_width_or_amplitude_materialized'] is False,'Physical scale was exponentiated')
    require(r['inherited_m15_packet']==p['inherited_m15_packet'] and r['inherited_m15_packet_digest']==p['inherited_m15_packet_digest'],'Original witness changed')
    require(r['inherited_joint180_views_digest']==digest(p['joint_discrepancy_views']) and r['inherited_leading_packets_digest']==p['inherited_leading_packets_digest'],'Accepted source view changed')
    ledger=load(DIRECTORY/r['inherited_m15_packet']['report']);require(len(r['quiet179_discrepancy_views'])==3,'Three saved queries required')
    for view,before in zip(r['quiet179_discrepancy_views'],p['joint_discrepancy_views']):
        require(view['Z']==before['Z'],'Axial query changed')
        sample=next(v for v in ledger['samples'] if v['Z']==view['Z']);original=sample['correlated_reduced_m15_C0_C1_terms']
        require(len(original)==183 and len(sample['raw_m15_C0_C1_terms'])==213,'Original archive scope changed')
        require(digest(original)==view['parent_original_reduced183_digest'] and digest(sample['raw_m15_C0_C1_terms'])==view['raw213_digest'],'Original archive bytes changed')
        require(view['entries'][:-1]==before['entries'][:-2],'A retained forward source changed')
        require(view['view_terms']==179 and view['unchanged_parent_entries']==178 and view['unchanged_original_references']==174,'Wrong quiet view count')
        require(view['parent_joint180_entries_digest']==digest(before['entries']) and view['parent_joint180_expansion_digest']==before['expanded_joint180_digest'],'Parent joint source digest changed')
        require(view['inherited_reserved_target_joint_source_digest']==before['joint_source_digest'] and
            view['nonquiet_remainder_source_digest']==digest(rem),'Current and inherited joint sources confused')
        last=view['entries'][-1]
        require(last['operation']==REMAINDER_OPERATION and last['parent_joint_source_entry_digest']==digest(before['entries'][-1]),'Wrong source group contracted')
        require(last['term']['coefficient_C0_C1']==rem['coefficient_C0_C1'] and last['term']['log_scale_terms']==before['independent_target']['log_scale_terms'],'Physical remainder coefficient/unit changed')
        require(last['source_mean_equality_proof_digest']==digest(proof),'Quiet source proof changed')
        prov=view['exact_zero_quiet_source_provenance']
        require(prov=={'parent_entry':before['entries'][-2],'original_quiet_term':original[181],
            'canonical_quiet_pair':[31215,0],'source_mean_equality_proof_digest':digest(proof),'actual_source_zero':True,
            'old_approximation_enclosure_retained_as_provenance':True,'old_bound_selected_as_value':False},'Quiet source provenance changed')
        expanded=materialize_quiet179_view(original,view['entries'],view['independent_target'])
        require(digest(expanded)==view['expanded_quiet179_digest'] and expanded[:6]==original[:6],'Expansion or six signed inner/lower changed')
        require(view['quiet_scalar_diagnostic_digest']==digest(diag) and view['quiet_pair_consumed_once'] is True,'Quiet pair/diagnostic duplicated')
        for key in ('independent_target','independent_target_digest','selected_N_evidence_digest','largest_unchanged_caps',
            'six_inner_lower_digest','contracted_parent_slots','exact_zero_endpoint_provenance','inherited_contracted_raw44_provenance'):
            require(view[key]==before[key],'Inherited source metadata changed: '+key)
        require(view['external_target_subtraction_required'] is False and view['independent_target_consumed_once'] is True and view['global_terminal_closure'] is False,'Target semantics or global limit changed')
    for key in ('input_hashes','implementation_hashes'):
        require(all(r[key].get(name)==h for name,h in proof[key].items()),'Missing quiet proof source pin')
    true=('same_source_integral_and_physical_K_identity','quiet_pair_consumed_once','target_already_consumed_once',
        'other178_parent_entries_unchanged','original_raw213_reduced183_six_signed_inner_lower_and_N_unchanged',
        'all_finite_N_caps_and_pressure_tail_retained','saved_independent_quadrature_diagnostic_retained')
    false=('numerical_overlap_midpoint_or_cap_used_as_source_value','accepted_ancestor_or_quadrature_replayed',
        'physical_radius_width_amplitude_materialized')
    require(all(r['contract'][k] is True for k in true) and all(r['contract'][k] is False for k in false),'Quiet contraction contract changed')
    require(all(v is False for v in r['limits'].values()),'Global completion introduced')
    return {'passed':True,'quiet179_views':3,'view_terms':179,'unchanged_parent_entries':178,'unchanged_original_references':174,
        'true_quiet_source_pair_zero':True,'independent_numeric_diagnostic_retained':True,
        'original_archives_six_N_and_all_finite_N_caps_retained':True,'global_terminal_closure':False}


def pins(r):
    n=0
    for base,key in ((DIRECTORY,'input_hashes'),(ROOT,'implementation_hashes')):
        for name,h in r[key].items():require(sha(base/name)==h,'Changed source pin: '+name);n+=1
    return n


def run():
    start=time.perf_counter();fixture=controls();o=OriginalQuietScalarM15();r=o.build()
    r.update(input_hashes=o.input_hashes,implementation_hashes=dict(o.implementation_hashes,**{n:sha(ROOT/n) for n in CODE}))
    r['checks']={'passed':True,'independent_new_scalar_controls':fixture,'packet_invariants':packet_check(r,o.parent),
        'no_accepted_ancestor_replay':True}
    r['execution_seconds']=round(time.perf_counter()-start,3);return r


def inspect(path):
    r=load(path);require(r['checks']['passed'] is True and r['checks']['independent_new_scalar_controls']['passed'] is True,'Passed saved report required')
    return {'passed':True,'pins':pins(r),'packet_invariants':packet_check(r,parent()),'serialized_only':True,
        'owner_fixture_or_ancestor_replayed':False}


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);g=ap.add_mutually_exclusive_group(required=True)
    g.add_argument('--write-report',action='store_true');g.add_argument('--inspect-report',action='store_true')
    args=ap.parse_args();path=DIRECTORY/NAME
    if args.write_report:
        r=run();tmp=path.with_suffix('.json.tmp');tmp.write_text(json.dumps(r,indent=2,ensure_ascii=False)+'\n',encoding='utf-8',newline='\n');tmp.replace(path)
    else:print(json.dumps(inspect(path),indent=2))
