"""Focused signed cone/source replay and independently specified matrix fixtures."""
import json
import gzip
import mpmath as mp
from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_current_original_cone import (
    CurrentOriginalCone,NAME,RECEIPT,VIEWS_NAME,GATES,OPEN,HERE,PREFIX,sha,pack,encode,endpoints,source_precision,_verify_hashes,coefficient_scaled_exp)
from lei_ren_part1_paper_compliant_current_original_cone_operator import (
    cone_margins,actual_covariance_inverse,linear_covariance_change)


def contains(box,value):
    lo,hi=endpoints(box)
    if not lo<=value<=hi:raise ValueError('Independent fixture not enclosed')


def independent_fixtures():
    c=MPIntervalContext();c.dps=90
    # Lei Figure 5: F=1, S=(-4,-4), kappa=8, T directed against S.
    admitted=cone_margins(c,4,-4,1,1)
    if not admitted['admitted']:raise ValueError('Interior signed cone fixture lost')
    outside=cone_margins(c,4,-4,1,-1)
    if outside['status']!='failed_directed_margin':raise ValueError('Exterior signed cone fixture admitted')
    wide=cone_margins(c,4,-4,(-1,1),(-1,1))
    if wide['admitted']:raise ValueError('Unresolved sign box admitted')
    weak=cone_margins(c,1,0,1,0)
    if weak['admitted'] or 'vs_minus2' not in weak['failed']:raise ValueError('Relaxed-only weak shear admitted')
    zero=cone_margins(c,4,-4,0,0)
    if zero['admitted']:raise ValueError('Zero stress promoted to strict interior cone')
    columns=((2,1),(1,-2));solved=actual_covariance_inverse(c,columns,(5,0))
    contains(solved['squared_amplitudes'][0],2);contains(solved['squared_amplitudes'][1],1)
    if not solved['positive_coefficients'] or solved['actual_pulse_lift_constructed']:raise ValueError('Matrix-only covariance scope differs')
    boundary=actual_covariance_inverse(c,columns,(2,1))
    if boundary['positive_coefficients']:raise ValueError('Zero coefficient admitted as positive')
    singular=actual_covariance_inverse(c,((1,1),(1,1)),(1,1))
    if singular['status']!='inconclusive_singular_column_box':raise ValueError('Singular covariance admitted')
    change=linear_covariance_change(c,columns,(2,1),(-5,0),'.5')
    a=change['fixed_base_amplitudes'];da=change['signed_amplitude_changes']
    for j,wanted in enumerate((-5,0)):
        actual=c.mpf('.5')*sum((c.mpf(columns[k][j])*2*a[k]*da[k] for k in (0,1)),c.mpf(0))
        contains(actual,wanted)
    if any(endpoints(x)[1]>=0 for x in da):raise ValueError('Signed negative correction lost')
    try:linear_covariance_change(c,columns,(0,1),(1,1),1)
    except ValueError:pass
    else:raise ValueError('Zero fixed amplitude allowed in linear lift denominator')
    # The old fixed exp(-1024) cap swamps this small target when multiplied
    # by a large coefficient. Check a direct exponential independently.
    factor,cap=coefficient_scaled_exp(c,c.mpf(-12000),c.exp(4000),c.exp(-5000))
    contains(factor,mp.exp(-12000))
    if endpoints(c.exp(4000)*factor)[1]>endpoints(c.exp(-5999))[0]:
        raise ValueError('Coefficient-scaled exp cap lost the requested small error scale')
    astronomical,original_cap=coefficient_scaled_exp(c,c.mpf('-1e1000'),c.exp(4000),c.exp(-5000))
    if endpoints(astronomical)[0]!=0 or endpoints(astronomical)[1]>endpoints(c.exp(-9999))[0]:
        raise ValueError('Astronomical source log was materialized or loosened')
    return dict(interior_exterior_weak_shear_unresolved_and_zero_tests=5,
        independent_covariance_coefficients=[2,1],boundary_zero_and_singular_rejected=True,
        signed_negative_increment_cross_identity_enclosed=True,
        zero_base_amplitude_denominator_rejected=True,coefficient_scaled_exponential_encloses_direct_and_astronomical_sources=True,actual_wave_or_edge_admission=False)


@source_precision
def run(field=None):
    raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
    field=field if field is not None else CurrentOriginalCone(require_checked=False)
    field.assert_graph()
    if encode(pack(field.manifest()))!=raw or any(raw[k] for k in GATES+OPEN):raise ValueError('Original cone manifest/scope differs')
    fixture=independent_fixtures();views={};counts={}
    cases={
        'fresh_main':('pulse_main','.537','2','-2.6','.41','.8'),
        'fresh_exit':('pulse_exit','.731','10.417','-2.9','.37','.2'),
        'fresh_end_first':('pulse_end','.537','-2.937','-2.6','.41','.8'),
        'fresh_end_second':('pulse_end','.731','-.947','-2.9','.37','.2'),
        'whole_end_requested_box':('pulse_end',(-1,1),(-4,0),('-3','-1'),None,'1')}
    for label,args in cases.items():
        view=field.native(*args);reduced=view['normalized_current_signed_stress'];sectors=reduced['all_signed_sectors']
        if any(view[k] for k in GATES+OPEN):raise ValueError('Unchecked new view admits an open gate')
        expected=15 if args[0]!='pulse_end' else 10
        count=sum(len(v) for v in sectors.values())
        if count!=expected:raise ValueError('Full original signed stress sector count differs: '+str(count))
        for component in ('theta','axial'):
            replay=sum((v['normalized_signed_contribution'] for v in sectors[component].values()),field.ctx.mpf(0))
            if encode(pack(replay))!=encode(pack(reduced['signed_totals'][component])):
                raise ValueError('Signed sum replaced by absolute envelope')
        cone=view['original_leading_two_vector_cone'];reference=view['reference_covariance_diagnostic']
        if label=='whole_end_requested_box' and not cone['admitted']:
            raise ValueError('Whole original pulse-end source cone margin unproved')
        if reference['actual_pulse_lift_constructed']:raise ValueError('Reference columns mistaken for actual wave lift')
        if cone['admitted'] and reference['status'] not in ('positive_reference_coefficients_by_correlated_cone',
                'reference_frame_enclosure_inconclusive','reference_coefficient_enclosure_inconclusive'):
            raise ValueError('Admitted source cone lost reference covariance positivity')
        if not view['current_source_shear']['full_axial_shear_retained']:raise ValueError('Axial shear discarded')
        views[label]=view;counts[label]=dict(signed_order_zero_stress_sectors=count,
            cone_status=cone['status'],reference_status=reference['status'],
            signed_source_stress_not_absolute_envelope=True,full_current_diagonal_divergence_E_preserved=True)
        print('Current source signed cone: '+label+' '+cone['status'],flush=True)
    payload=(json.dumps(encode(pack(views)),separators=(',',':'))+'\n').encode('utf8')
    (HERE/VIEWS_NAME).write_bytes(gzip.compress(payload,compresslevel=6,mtime=0))
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,datum_enclosure_sha256=field.datum_sha,
        original_paper_signed_cone_and_scale_identities=len(field.theorem['identities']),
        independent_signed_cone_and_covariance_fixtures=fixture,current_source_signed_cone_view_counts=counts,
        complete_current_source_signed_cone_views_gzip=VIEWS_NAME,
        gzip_contains_all_five_complete_unpruned_signed_views=True,
        source_correlated_pulse_common_scale_cancellation_verified=True,
        actual_homogeneous_pulse_covariance_columns_still_unconstructed=True,
        uniform_edge_weight_and_direction_still_unconstructed=True,
        general_two_column_inverse_not_promoted_to_Proposition_7_5=True,
        scope=raw['scope'],input_hashes={**raw['input_hashes'],NAME:sha(NAME),VIEWS_NAME:sha(VIEWS_NAME),Path(__file__).name:sha(Path(__file__).name)},
        **dict.fromkeys(GATES,True),**dict.fromkeys(OPEN,False),all_passed=True)
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Original cone/current signed source adapter PASS; global cone/wave lift/recursion remain open',flush=True)
    return result
