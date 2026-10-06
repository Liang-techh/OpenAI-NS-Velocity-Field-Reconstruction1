"""Analytic cover admission, conservative root stops and physical routing."""
import copy
from dataclasses import replace
import itertools
import json
from pathlib import Path
import mpmath as mp
import sympy as s
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_current_global_tensor_cover import (
    CurrentGlobalTensorCover,NAME,RECEIPT,GATES,OPEN,HERE,sha,pack,encode,endpoints,
    source_precision,_verify_hashes,audit_source_hypotheses,SYMBOLS,LOG_HB)
from lei_ren_part1_paper_compliant_current_physical_tensor_locator_operator import implicit_log_coordinate_map
from lei_ren_part1_paper_compliant_current_physical_tensor_locator_check import independent_root,check_refined_view


@source_precision
def nonnarrow_root_enclosures():
    c=MPIntervalContext();c.dps=90;evidence={}
    cases=(('step_limit','.7','-2.6','.8','.37'),('near_delta_one','7','-8','.2','.999'),
        ('zero_delta','-2','-3','1.2','0'),('small_time','1','-1000','.8','.93'))
    for label,z,lt,nu,d in cases:
        out=implicit_log_coordinate_map(c,z,lt,nu,d,max_steps=1)
        lo,hi=endpoints(out['actual_log_lambda']);rlo,rhi=independent_root(z,lt,nu,d)
        if not (mp.isfinite(lo) and mp.isfinite(hi) and lo<=rlo<=rhi<=hi):
            raise ValueError('Limited-step finite bracket lost the original physical root: '+label)
        evidence[label]=dict(independent_direct_equation_root_retained=True,solver_status=out['solver_status'],
            actual_log_lambda=out['actual_log_lambda'],wide_bracket_is_valid_cover=True)
    out=implicit_log_coordinate_map(c,('-2','3'),('-7','-2'),('.3','1.2'),('.1','.999'),max_steps=1)
    lo,hi=endpoints(out['actual_log_lambda'])
    for point in itertools.product(('-2','3'),('-7','-2'),('.3','1.2'),('.1','.999')):
        rlo,rhi=independent_root(*point)
        if not lo<=rlo<=rhi<=hi:raise ValueError('Wide finite parameter bracket lost a corner root')
    evidence['wide_parameter_box']=dict(independent_corner_roots_retained=16,solver_status=out['solver_status'])
    # The root/complement is represented by exact logs; no exp(-10^1000)
    # value or chosen endpoint-field is materialized to claim coverage.
    extreme=implicit_log_coordinate_map(c,'1','-1e1000','.8','.37',max_steps=1)
    if not all(mp.isfinite(v) for v in endpoints(extreme['actual_log_lambda'])):
        raise ValueError('Extreme finite log-time bracket became infinite')
    evidence['extreme_finite_log_time']=dict(finite_bracket_retained=True,
        true_abs_Z_less_than_one=extreme['true_finite_point_has_abs_Z_strictly_below_one'],
        solver_status=extreme['solver_status'],Z=extreme['Z'])
    return evidence


@source_precision
def run(field=None):
    raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
    field=field if field is not None else CurrentGlobalTensorCover(require_checked=False)
    field.assert_graph()
    if encode(pack(field.manifest()))!=raw or any(raw[k] for k in GATES+OPEN):
        raise ValueError('Source-bound global cover manifest differs')
    radial=raw['independent_original_radial_cover_theorem'];axial=raw['independent_global_implicit_coordinate_cover_theorem']
    if not (radial['passed'] and axial['passed'] and field.candidate['passed'] and field.width['passed']
        and len(radial['independent_original_33_radial_sources_and_positive_derivatives'])==33
        and len(radial['independent_all_32_contiguous_source_joins'])==32
        and field.hypotheses['all_original_33_radial_monotonicity_premises_audited']):
        raise ValueError('Independent global analytic cover or actual hypotheses incomplete')
    root_evidence=nonnarrow_root_enclosures();print('Global cover: independent analytic proof and conservative limited-step roots PASS',flush=True)
    views={};families={}
    cases={
      'physical_axis':field.cartesian('0','0','.7','.8',viscosity='.8'),
      'finite_extreme_time_off_axis':field.tensor_log_radius('0','1','-1e1000','.41','.8'),
      'joint_core_boundary':field.correlated_tensor(field.field.radius('Ra'),'.7','-2.6','.41','.8'),
      'unmaterialized_Gamma':field.correlated_tensor(field.field.radius_expression(
          field.field.boundaries['collar_exterior']['expression']+s.exp(SYMBOLS['Md']*10**10)),'.7','-2.6','.41','.8')}
    for label,result in cases.items():
        original=result['original_physical_source_result'];location=original['location']
        if any(result[k] for k in GATES+OPEN) or not result['independent_analytic_cover_applies']:
            raise ValueError('Unchecked fresh cover view scope differs')
        views[label]=[check_refined_view(view,location['actual_log_lambda']) for view in original['candidate_full_tensor_views']]
        families[label]=dict(regions=[row['region'] for row in location['candidates']],
            solver_status=location['solver_status'],actual_log_tau=location['requested_log_tau'],Z=location['Z'],
            original_source_rows_and_actual_lambda_retained=True)
        print('Global cover fresh complete physical view: '+label,flush=True)
    if families['physical_axis']['regions']!=['core_positive_radius'] or families['joint_core_boundary']['regions']!=['core_positive_radius','bridge_first']:
        raise ValueError('Analytic Cartesian axis/joint core boundary cover lost')
    if families['unmaterialized_Gamma']['regions']!=['heat_exterior'] or any(row['nonzero_signed_contributions'] for row in views['unmaterialized_Gamma']):
        raise ValueError('Unbounded exact Gamma cover lost')
    request=field.field.radius('Rref','-7.1');correlated=field.field.locate(request,'.7','-2.6','.41','.8')
    absolute=field.locate_log_radius(correlated['physical_log_radius_enclosure'],'.7','-2.6','.41','.8')
    candidates=absolute['original_physical_source_result']['candidates']
    if len(candidates)<=1 or not absolute['source_dependent_or_absolute_input_uncertainty_retained']:
        raise ValueError('Global cover falsely replaced uncertain absolute region with a single candidate')
    rejected=[];values=field.field._values(field.field.radius('Ra'));c=field.ctx
    for label,symbol,value,delta in (
      ('zero_mu',SYMBOLS['mu'],'0',field.locator.physical.delta),
      ('large_mu',SYMBOLS['mu'],'.3',field.locator.physical.delta),
      ('zero_Md',SYMBOLS['Md'],'0',field.locator.physical.delta),
      ('zero_T',SYMBOLS['T'],'0',field.locator.physical.delta),
      ('short_outer',SYMBOLS['Lrel'],'4',field.locator.physical.delta),
      ('large_hb',LOG_HB,'0',field.locator.physical.delta),
      ('nonfinite_logh',LOG_HB,'-inf',field.locator.physical.delta),
      ('delta_one',None,None,c.mpf(1))):
        changed=dict(values)
        if symbol is not None:changed[symbol]=c.mpf(value)
        try:audit_source_hypotheses(c,changed,delta,field.field._eval)
        except ValueError:rejected.append(label)
        else:raise ValueError('Invalid original global cover premise accepted: '+label)
    location=absolute['original_physical_source_result']
    for label,key,value in (('foreign_family','actual_five_defect_family_sha256','foreign'),
        ('changed_candidate_union','candidates',[location['candidates'][0]]),
        ('changed_physical_z','physical_z',c.mpf(2))):
        forged=copy.deepcopy(location);forged[key]=value
        try:field._attach(forged)
        except ValueError:rejected.append(label)
        else:raise ValueError('Forged global source location accepted: '+label)
    try:field.correlated_tensor(replace(request,source='foreign'),'.7','-2.6','.41','.8')
    except ValueError:rejected.append('foreign_correlated_request')
    else:raise ValueError('Foreign correlated source accepted by global cover')
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,datum_enclosure_sha256=field.datum_sha,
        independent_analytic_cover_proofs_and_actual_source_premises_admitted=True,
        independent_conservative_nonnarrow_root_evidence=root_evidence,
        current_fresh_global_physical_tensor_view_counts=views,current_fresh_global_physical_source_families=families,
        global_cover_retains_absolute_parameter_uncertainty_candidate_count=len(candidates),
        invalid_hypotheses_and_foreign_source_locations_rejected=rejected,
        exact_shared_width_witness_sha256=field.width['witness_sha256'],
        source_union_coverage_is_separate_from_unique_seam_selection=True,
        axis_and_unbounded_Gamma_use_the_original_special_sources=True,
        scope=raw['scope'],input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)},
        all_passed=True,**dict.fromkeys(GATES,True),**dict.fromkeys(OPEN,False))
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Independent source-bound global physical T/E cover PASS; cone/lift/points/NS/energy/flatness/recursion open',flush=True)
    return result


if __name__=='__main__':run()
