"""Independent Decimal source-phase and exact radius-origin checks.

This checks the new lightweight phase path without reconstructing ancestors.
Candidate integers are diagnostic queries, never a selected global frequency.
"""
from decimal import Decimal, localcontext
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_O2_radius_phase_points as current


def read(c,value):
    return c.make_mpf(tuple(value['exact_mpf_tuple'])) if isinstance(value,dict) else c.mpf(value)


def references(owner,manifest):
    c=mp.mp.clone();c.dps=650;maximum=c.mpf(0);wrapped_count=0
    with localcontext() as context:
        context.prec=600
        for query in manifest['original_O2_radius_phase_queries']:
            y=s.Rational(query['original_y_exact']);N=query['explicit_candidate_N']
            dy=Decimal(int(y.p))/Decimal(int(y.q))
            # Independent Decimal exp/ln implementation and arithmetic.
            main=Decimal(N)*(14*(Decimal(40).exp()+11)+(Decimal(110)/4).ln()+dy)
            reference=c.mpf(str(main-main//1))
            point=read(c,query['approximate_original_fractional_phase'])
            arithmetic=read(c,query['phase_arithmetic_absolute_error_upper'])
            discrepancy=abs(reference-point)
            assert discrepancy<=arithmetic+c.mpf('1e-450')
            boxes=[(read(c,box['lower']),read(c,box['upper'])) for box in query['true_original_phase_directed_boxes']]
            assert all(0<=lo<=hi<=1 for lo,hi in boxes)
            assert any(lo<=reference<=hi for lo,hi in boxes)
            assert not query['periodic_projection_full_period']
            offset=query['positive_original_origin_offset']
            log_upper=read(c,offset['error_log_upper']);budget=read(c,offset['representable_absolute_error_budget_upper'])
            assert log_upper<c.ln(budget) and budget>0
            assert offset['strictly_positive'] and offset['sign_in_original_affine_phase']==-1
            assert offset['actual_hb_sc_product_not_materialized_or_set_to_zero']
            exact=query['exact_removed_integer_periods'];dyadic=exact['logCstar']
            assert dyadic['source_exact_mpf_tuple']==list(owner.frame.selected_logCstar_mpf_tuple)
            assert dyadic['exact_multiplier']=={'numerator':10*N,'denominator':1}
            assert dyadic['exact_fraction']=={'numerator':0,'denominator':1}
            assert not dyadic['huge_integer_materialized'] and exact['constant_1000N']==1000*N
            assert query['source_family']==owner.family and query['same_original_radius_phase_bound']
            assert query['source_radius_caps_not_consumed'] and query['candidate_N_not_global_selection']
            assert not query['conditioned_loop_inverse_or_primitives_installed']
            assert not query['numerical_original_source_point_or_integral_oracle_installed']
            origin=query['original_phase_origin']
            assert origin['log_r_minus']=='logRa+hb*s_c/2'
            assert origin['log_positive_log_radius_offset']=='loghb+logsc-log2'
            assert origin['selected_sc_exact_mpf_tuple']==list(owner.sc_tuple)
            maximum=max(maximum,discrepancy);wrapped_count+=len(boxes)
    logP,logC,y,hb,sc=s.symbols('logP logC y hb sc',real=True)
    logR=s.log(110)+10*(logC+logP)+y
    logRa=s.log(4)-4*logP-1000
    logminus=logRa+hb*sc/2
    expected=s.log(s.Rational(110,4))+14*logP+10*logC+1000+y-hb*sc/2
    assert s.expand_log(s.simplify(logR-logminus-expected),force=True)==0
    assert s.diff(logR-logminus,y)==1
    # Independent wrapped-cover check for a retained negative offset across0.
    iv=current.point.MPIntervalContext();iv.dps=80
    wrapped=current.phase.ordinary_mod_one(iv,iv.mpf(['-.001','.0005']))
    assert not wrapped['full_period'] and len(wrapped['boxes'])==2
    lower=[current.point.endpoints(box) for box in wrapped['boxes']]
    assert any(lo<=c.mpf('.9995')<=hi for lo,hi in lower)
    assert any(lo<=c.mpf('.0002')<=hi for lo,hi in lower)
    rejected=0
    for call in (lambda:owner.evaluate(y=-1,N=1),lambda:owner.evaluate(y='nan',N=1),
        lambda:owner.evaluate(y='.5',N=0),lambda:owner.evaluate(y='.5',N=True),
        lambda:owner.evaluate(y='.5',N=1.5),lambda:owner.evaluate(y='.5',N=1<<4096),
        lambda:current.OriginalO2RadiusPhasePoints(dps=20)):
        try:call()
        except ValueError:rejected+=1
    assert rejected==7
    return dict(passed=True,independent_Decimal_original_source_phase_queries=len(manifest['original_O2_radius_phase_queries']),
        Decimal_reference_digits=600,maximum_point_reference_discrepancy=maximum,
        exact_original_logR_logminus_affine_identity=True,original_radial_phase_slope_identity=True,
        retained_positive_origin_error_budget_log_comparisons=True,
        huge_integer_period_removal_exactly_bound_to_original_selected_parameter=True,
        original_phase_directed_cover_boxes=wrapped_count,negative_offset_wrap_union_checked=True,
        invalid_domain_candidate_integer_and_precision_rejections=rejected,
        no_ancestor_constructors_or_source_quadratures_reexecuted=True)


def run():
    began=time.monotonic();owner=current.OriginalO2RadiusPhasePoints()
    manifest=json.loads((current.HERE/current.NAME).read_bytes())
    assert manifest[current.GATE] and manifest['source_family']==owner.family
    for name,digest in manifest['input_hashes'].items():assert current.sha(name)==digest,name
    flags=('conditioned_native_factor_and_loop_evaluation_installed','numerical_original_source_point_or_integral_oracle_installed',
        'actual_five_controls_installed','current_whole_N_selected',*current.point.source.inertial.profiles.loop.OPEN)
    assert all(manifest[key] is False for key in flags)
    result=dict(all_passed=True,**{current.GATE:True},source_family=owner.family,
        independent_phase_origin_and_error_references=references(owner,manifest),
        **dict.fromkeys(flags,False),input_hashes={**manifest['input_hashes'],current.NAME:current.sha(current.NAME),
            Path(__file__).name:current.sha(Path(__file__).name)},execution_seconds=time.monotonic()-began,
        scope='Genuine original O2 radius phase point/error service at candidate N only. Original source-bound positive offset retained. Native factor/loop inverse, primitive evaluation, signed integral, full oracle, controls and global recursion remain open.')
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.point.source.inertial.profiles.loop.encoded(result),indent=2)+'\n',encoding='utf8')
    print('Original O2 true radius phase points: independent Decimal/origin, offset and wrap checks PASS',flush=True)
    return result


if __name__=='__main__':run()
