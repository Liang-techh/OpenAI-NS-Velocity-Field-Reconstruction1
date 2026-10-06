"""Independent source-radius targets, exact boundaries, crossings and tensors."""
import copy
from dataclasses import replace
import json
from pathlib import Path
import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_current_correlated_tensor_locator import (
    CurrentCorrelatedTensorLocator,SourceRadiusRequest,NAME,RECEIPT,GATES,OPEN,
    HERE,sha,pack,encode,endpoints,source_precision,_verify_hashes)
from lei_ren_part1_paper_compliant_current_correlated_radius_operator import (
    SYMBOLS,OFFSET,LOG_HB,anchors,source_geometry,boundary_sources,canonical_radius)
from lei_ren_part1_paper_compliant_current_physical_tensor_locator_check import check_refined_view


def independent_source_targets():
    # Deliberately specified from shared physical construction anchors, not
    # by calling a native radius routine and immediately inverting its output.
    a=SYMBOLS;h=a['h_bridge'];mu=a['mu'];md=a['Md'];r=a['log_Rref'];p=a['log_Rp'];v=a['log_Rv']
    rel=v+100+a['Lrel'];q=rel+2+a['Ts'];tail=q+a['wait']
    return {
     'core_positive_radius':a['log_epsilon']+s.log(s.Rational(137,1000)),
     'bridge_first':a['log_Ra']+h*s.Rational(417,1000),
     'bridge_second':a['log_Ra']+h*s.Rational(1417,1000),
     'bridge_macro':(a['log_Ra']+2*h+s.log(100))/2,
     'switch_first':s.log(100)+h*s.Rational(417,1000),
     'switch_second':s.log(100)+h*s.Rational(1417,1000),
     'switch_power':(s.log(100)+2*h+s.log(110))/2,
     'reshape':s.log(110)+a['T']*s.Rational(17,40),
     'inner_reference':r-9,'axial_restore':r-s.Rational(15,2),'restore_buffer':r-s.Rational(13,2),
     'actual_patch':r-6+s.log(s.Rational(37,25)),'Rh_reference':r-s.Rational(27,10),
     'O2_slope':r+s.Rational(317,1000),'O2_axial':r+s.exp(md*s.Rational(3,8)),
     'O2_buffer':r+s.exp(md)+s.Rational(37,10),'O3_slope_mu':r+a['log_P']+s.Rational(37,100),
     'O3_power':p-a['Tw']/3,'pulse_entrance':p+1/(100*mu),'pulse_main':p+2/mu,
     'pulse_exit':p+21/(2*mu),'pulse_gap':p+23/(2*mu),
     'pulse_gap_end':v-4-(1/mu-4)*s.Rational(2,5),'pulse_end':v-s.Rational(11,4),
     'flatten':v+37,'outer_power':v+100+(a['Lrel']-4)*s.Rational(3,8),
     'outer_angular':rel-s.Rational(11,4),'steep_entry':rel+s.Rational(37,100),
     'steep_power':rel+1+a['Ts']*s.Rational(3,8),'steep_exit':rel+1+a['Ts']+s.Rational(37,100),
     'waiting':q+a['wait']*s.Rational(3,8),'heat_collar':tail+s.Rational(17,10),
     'heat_exterior':tail+s.Rational(63,10)}


def names(location):return [row['region'] for row in location['candidates']]


@source_precision
def run(field=None):
    raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
    field=field if field is not None else CurrentCorrelatedTensorLocator(require_checked=False)
    field.assert_graph()
    if encode(pack(field.manifest()))!=raw or any(raw[k] for k in GATES+OPEN):raise ValueError('Correlated source geometry manifest differs')
    if canonical_radius(s.log(10)+s.log(11)-s.log(110))!=0 or canonical_radius(SYMBOLS['h_bridge'])!=s.exp(LOG_HB):
        raise ValueError('Exact rational logarithm or positive width source identity lost')
    targets=independent_source_targets()
    if set(targets)!=set(field.geometry):raise ValueError('All 33 independently specified source targets required')
    positions={};target_evidence={}
    for expected,expr in targets.items():
        request=field.radius_expression(expr);location=field.locate(request,'.7','-2.6','.41','.8')
        if names(location)!=[expected] or location['exact_source_boundary_matches']:
            raise ValueError('Independent correlated physical radius branch differs: '+expected+' '+str(names(location)))
        positions[expected]=(request,location)
        target_evidence[expected]=dict(unique_region=True,exact_inverse_source=location['candidates'][0]['exact_native_coordinate_source'])
        print('Check independent correlated radius target: '+expected,flush=True)
    # This is the same source data which previously produced 25 outer
    # candidate regions. The changed input preserves its exact correlation;
    # fixed absolute-coordinate behavior is not claimed to be repaired.
    old=field.locator.locate_log_radius(positions['axial_restore'][1]['physical_log_radius_enclosure'],'.7','-2.6','.41','.8')
    if len(old['candidates'])<=1:raise ValueError('Independent ambiguity comparison fixture no longer exercises common-radius loss')
    boundaries={};geometry=source_geometry()[0];h=SYMBOLS['h_bridge']
    for name,row in field.boundaries.items():
        request=field.radius_expression(row['expression']);location=field.locate(request,'.7','-2.6','.41','.8')
        if location['exact_source_boundary_matches']!=[name]:raise ValueError('Exact original source radius selector differs: '+name)
        expected=([row['left_region'],row['right_region']] if row['kind']=='adjacent' else [row['region']])
        if names(location)!=expected:raise ValueError('Exact boundary neighboring source regions differ: '+name+' '+str(names(location)))
        boundaries[name]=dict(exact_original_radius_identity=True,selected_complete_trace_route=name,candidate_regions=expected)
        # Probe both sides using a strictly positive exact h, despite its
        # numerical cap containing zero. This is never a nominal absolute
        # radius perturbation which could round to the same endpoint.
        if row['kind']=='adjacent':
            for direction,region in ((-1,row['left_region']),(1,row['right_region'])):
                near=field.radius_expression(row['expression']+direction*h*s.Rational(1,10**30))
                located=field.locate(near,'.7','-2.6','.41','.8')
                if names(located)!=[region] or located['exact_source_boundary_matches']:
                    raise ValueError('Exact positive source displacement lost at '+name+' '+str(direction)+' '+str(names(located)))
            boundaries[name]['two_nonzero_hb_source_nudges_unique']=True
        print('Check exact correlated source boundary: '+name,flush=True)
    crossings={}
    for label,expr,box,expected in (
        ('micro',SYMBOLS['log_Ra']+h*(1+OFFSET),('-.1','.1'),['bridge_first','bridge_second']),
        ('restore',SYMBOLS['log_Rref']+OFFSET,('-7.1','-6.9'),['axial_restore','restore_buffer']),
        ('pulse',SYMBOLS['log_Rp']+OFFSET/SYMBOLS['mu'],('9.9','10.1'),['pulse_main','pulse_exit']),
        ('O2_lower',SYMBOLS['log_Rref']+OFFSET,('.9','1.1'),['O2_slope','O2_axial']),
        ('O2_upper',SYMBOLS['log_Rref']+s.exp(SYMBOLS['Md'])+OFFSET,('-.1','.1'),['O2_axial','O2_buffer'])):
        request=field.radius_expression(expr,box);location=field.locate(request,('-1','2'),('-4','-2'),'.41',('.7','1.3'))
        if names(location)!=expected or location['exact_source_boundary_matches']:
            raise ValueError('Crossing source box not retained as union: '+label+' '+str(names(location)))
        crossings[label]=dict(candidate_regions=expected,exact_boundary_not_asserted_for_crossing_box=True)
    # Finite source radius even when its exp cannot be represented by an MPF
    # exponent integer. Only the admitted unbounded exact Gamma tensor is used.
    giant=field.radius_expression(anchors()['Rt']+s.exp(SYMBOLS['Md']*10**10))
    tail=field.tensor(giant,'.7','-2.6','.41','.8')
    if names(tail['location'])!=['heat_exterior']:raise ValueError('Astronomical correlated tail not placed in original Gamma')
    views={}
    for region in ('bridge_second','axial_restore','pulse_gap_end'):
        request,_=positions[region];result=field.tensor(request,'.7','-2.6','.41','.8')
        if names(result['location'])!=[region]:raise ValueError('Fresh correlated tensor changed region')
        views[region]=check_refined_view(result['candidate_full_tensor_views'][0],result['location']['actual_log_lambda'])
    views['unmaterialized_Gamma']=check_refined_view(tail['candidate_full_tensor_views'][0],tail['location']['actual_log_lambda'])
    if views['unmaterialized_Gamma']['nonzero_signed_contributions']:raise ValueError('Exact Gamma zero tensor replaced')
    traces={}
    for name in ('gap_end','patch_support_49','outer_angular_-3_-1'):
        request=field.radius_expression(field.boundaries[name]['expression']);location=field.locate(request,'.7','-2.6','.41','.8')
        value=field.traces(request,location)
        selected=value['exact_original_source_radius_trace_routes']
        if list(selected)!=[name] or not selected[name]['common_actual_tensor_rows']:
            raise ValueError('Real exact source physical trace route failed: '+name)
        traces[name]=len(selected[name]['common_actual_tensor_rows'])
    axis=field.locate(field.radius_expression(-s.oo),'.7','-2.6','.41','.8')
    if names(axis)!=['core_positive_radius'] or endpoints(axis['candidates'][0]['native_coordinate_enclosure'])!=(mp.mpf(0),mp.mpf(0)):
        raise ValueError('Exact axis source radius extension lost')
    rejected=[]
    for label,expression,box in (
        ('unknown_symbol',s.Symbol('foreign_radius'),None),('float_source',s.Float('.1'),None),
        ('unbound_offset',OFFSET,None),('nonfinite_offset',OFFSET,('0','inf')),
        ('invalid_log',s.log(OFFSET),('-1','1')),('zero_touching_log',s.log(OFFSET),('0','1')),
        ('zero_denominator',1/OFFSET,('-1','1'))):
        try:field.radius_expression(expression,box)
        except ValueError:rejected.append(label);continue
        raise ValueError('Invalid correlated source radius request admitted: '+label)
    foreign=replace(positions['bridge_second'][0],source='foreign')
    try:field.locate(foreign,'.7','-2.6','.41','.8')
    except ValueError:rejected.append('foreign_source_request')
    else:raise ValueError('Foreign source-correlated point family admitted')
    request=field.radius_expression(field.boundaries['gap_end']['expression']);location=field.locate(request,'.7','-2.6','.41','.8')
    forged=copy.deepcopy(location);forged['exact_source_boundary_matches']=['core_bridge']
    try:field.traces(request,forged)
    except ValueError:rejected.append('forged_boundary_selector')
    else:raise ValueError('Forged source boundary trace admitted')
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,datum_enclosure_sha256=field.datum_sha,
        current_33_independent_correlated_source_radius_targets=target_evidence,
        absolute_radius_comparison_candidate_count=len(old['candidates']),correlated_radius_comparison_candidate_count=1,
        current_46_exact_source_radius_boundary_selectors=boundaries,current_five_crossing_source_boxes=crossings,
        current_four_fresh_correlated_complete_tensor_views=views,current_fresh_exact_source_trace_counts=traces,
        actual_axis_source_extension_retained=True,
        exact_shared_radius_source_cancellation_precedes_every_candidate_enclosure=True,
        exact_positive_hb_source_not_zero_cap_used_for_strict_signs=True,
        exact_numeric_logarithm_identities_canonicalized_before_enclosure=True,
        requests_are_source_dependent_physical_point_families_not_fixed_absolute_coordinates=True,
        invalid_source_requests_and_foreign_boundary_selection_rejected=rejected,
        scope=raw['scope'],input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)},
        all_passed=True,**dict.fromkeys(GATES,True),**dict.fromkeys(OPEN,False))
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Correlated source radius inverse / 46 exact radius trace selectors PASS; physical global cover/cone/point/recursion open',flush=True)
    return result


if __name__=='__main__':run()
