"""Actual regular n=1 Volterra norm and exact-polynomial radial tail.

This consumes saved, accepted complex source bounds, never their builder.
Canonical equation/axis/uniqueness and forward prefix bindings are required
separately before any completion gate is emitted. No coefficient enclosure
center or physical scale is selected as a function value.
"""
from __future__ import annotations

from fractions import Fraction
import math
from pathlib import Path

from mpmath.ctx_iv import MPIntervalContext
from mpmath.libmp import mpf_pos,mpf_cmp,round_floor

from .original_finite_N_signed_source import load,sha,digest,require
from .original_n1_complex_source_bounds import _interval,_record,_serialization_roundtrip

PARENT='lei_ren_part1_paper_compliant_current_original_C5_terminal_primitive.json'
PARENT_SHA256='588d221fd082d5f84d94f1812eb6342bb6647847b114dac458c0165324a063b1'


def _binding(packet,name):
    require(packet.get('passed') is True and packet.get('parent_report_sha256')==PARENT_SHA256,
        'Actual same-parent '+name+' binding required')
    require(packet.get('checks') and all(v is True for v in packet['checks'].values()),
        'A '+name+' source/equation condition failed')


def build(root,parent,canonical,prefix):
    root=Path(root);path=root/'experiments/root_st073'/PARENT
    require(sha(path)==PARENT_SHA256 and digest(parent)==digest(load(path)),
        'Actual accepted complex-source parent required')
    _binding(canonical,'canonical solution');_binding(prefix,'raw prefix')
    require(all(canonical.get(k) is True for k in ('actual_canonical_owner_bound',
        'actual_zero_axis_condition_bound','actual_raw_equations_bound',
        'actual_raw_B0_B1_g_same_source','actual_axis_weighted_transform_bound',
        'actual_source_normalization_bound','analytic_regular_uniqueness_bound',
        'actual_homogeneous_word_uniqueness_bound','actual_canonical_regular_solution_exists')),
        'Actual canonical owner/equation/axis/uniqueness was not bound')
    prefix_data=prefix['actual_forward_prefix_bound']
    require(prefix_data['passed'] is True and prefix_data['raw_base_radial_degree']==16 and
        prefix_data['weighted_regularized_radial_degree']==15,
        'Actual degree16 value prefix/weighted degree15 map required')
    common_degree=prefix_data['ordinary_Z_domain']['max_rho_degree_with_orders_0_to_3']
    require(type(common_degree) is int and common_degree==14 and
        prefix_data['ordinary_Z_domain']['max_weighted_degree_with_orders_0_to_3']==13,
        'The saved prefix does not have a common ordinary Z0..3 degree14 jet')
    bounds=parent['n1_complex_actual_source_bounds']
    require(bounds['passed'] is True and bounds['actual_source_bounds_proved'] is True and
        bounds['actual_weighted_operator_bounds_proved'] is True and
        bounds['diagnostic_only'] is False and bounds['missing_input_contract']==[],
        'Accepted actual source/operator bounds required')
    require(all(v is True for k,v in bounds['gates'].items()
        if k not in ('solution_norm','radial_Taylor_remainder')),
        'An accepted actual complex domain/source hypothesis is missing')
    green=parent['n1_weighted_complex_Green_kernel']
    require(green['passed'] is True and green['actual_weighted_matrix_entries']['passed'] is True and
        green['exact_controls']['passed'] is True and green['axis_weights']==[2,2,2,2,1,1] and
        green['kinematic_preconditioning_already_in_six_kernels'] is True and
        green['extra_I_plus_J_factor_applied'] is False and
        green['Z_loss_allocation']=='total gap sigma, sigma/k for each of k applications',
        'Accepted exact weighted Green/axis lemma required')
    require(bounds['domain']['rho_bound_radius']=='9/2' and
        bounds['domain']['rho_target_radius']=='4' and
        bounds['domain']['actual_collar_used_as_core_extension'] is False and
        bounds['green_input']['sigma_definition']=='Gg_Z_radius - solution_Z_radius',
        'Actual analytic core radius/Cauchy gap changed')
    ctx=MPIntervalContext();ctx.dps=240
    Msrc=_interval(ctx,bounds['green_input']['Msrc'])
    beta=_interval(ctx,bounds['green_input']['beta'])
    require(Msrc.a>=0 and beta.a>=0 and beta.b<=ctx.mpf(1)/2,
        'Outward source norm and beta<=1/2 required')
    # The exact comparison is with the saved binary endpoints. Decimal
    # displays never enter either the solution norm or its tail bounds.
    Mv=2*Msrc
    radius=ctx.mpf(9)/2;target=ctx.mpf(4);q=target/radius
    zradius=_interval(ctx,bounds['domain']['solution_Z_radius'])
    require(zradius.a>0,'Positive solution Z tube required')
    domain_name='lei_ren_part1_paper_compliant_current_original_n1_inner_analytic_domain.json.gz'
    domain_path=root/'experiments/root_st073'/domain_name
    require(sha(domain_path)==parent['input_hashes'][domain_name],'Actual local complex domain changed')
    domain=load(domain_path)
    require(domain['joint_core_domain']['Z_real_centers']==['-1','1'] and
        domain['joint_core_domain']['rho_complex_modulus_less_than']=='5' and
        domain['source_family']==parent['source_family'],'Different local-center construction family')
    # This is a fresh inner gap after the Green-word Z losses have already
    # been consumed. The disks have each real Z0 in (-1,1) as their center.
    # Saved endpoints can carry more bits than this arithmetic context.
    # Round the positive lower endpoint DOWN before using it as an exact
    # domain parameter; otherwise dividing it can round the target upward.
    safe_endpoint=mpf_pos(zradius._mpi_[0],ctx.prec,round_floor)
    require(mpf_cmp(safe_endpoint,zradius._mpi_[0])<=0,'Working radius exceeds the saved lower endpoint')
    safe_zradius=ctx.make_mpf((safe_endpoint,safe_endpoint))
    zgap=safe_zradius/2
    require(mpf_cmp((zgap+zgap)._mpi_[1],zradius._mpi_[0])<=0,
        'Target disk plus new Cauchy gap leaves the admitted disk')
    degree=common_degree-1
    weighted_tail=Mv*q**common_degree/(1-q)
    base_tail=target*weighted_tail
    odd_tail=ctx.sqrt(target)*weighted_tail
    rows=[]
    for order in range(4):
        zfactor=math.factorial(order)/zgap**order
        rows.append(dict(ordinary_Z_order=order,cauchy_factor=_record(zfactor),
            weighted_v_tail_upper_majorant=_record(weighted_tail*zfactor),
            raw_base_W1_W4_tail_upper_majorant=_record(base_tail*zfactor),
            raw_odd_W5_W6_tail_upper_majorant=_record(odd_tail*zfactor)))
    factor=Fraction(8,9)**common_degree*9
    require(factor==Fraction(8,9)**common_degree/(1-Fraction(8,9)),'Exact radial factor changed')
    value16_tail=Mv*q**16/(1-q)
    # A usable stopping error for future word evaluation, distinct from the
    # rho coefficient truncation. No partial word sum is evaluated here.
    word_tails=[dict(last_retained_word=k,
        omitted_words_upper_majorant=_record(Msrc*beta**(k+1)/(1-beta))) for k in (0,1,4,8)]
    result=dict(passed=True,temporal_order=1,parent_report=PARENT,parent_report_sha256=PARENT_SHA256,
        actual_canonical_binding_sha256=digest(canonical),actual_raw_prefix_binding_sha256=digest(prefix),
        accepted_source_bound_packet_sha256=digest(bounds),accepted_Green_packet_sha256=digest(green),
        actual_canonical_regular_solution_norm_proved=True,
        actual_canonical_raw_radial_Taylor_remainder_certified=True,
        exact_forward_to_raw_prefix_bound=True,
        solution_definition='Normal-convergent sum H^k W0, k>=0; W0 is the preconditioned green_apply(g) output with actual zero-axis data',
        convergence_proof=dict(word_bound='Msrc*(R/2)^k*(C0+k*C1/sigma)^k/k!',
            W0='preconditioned green_apply(g); its six weighted kernels already contain I+J',
            H='preconditioned green_apply(C0*W+B1*dZ(W)); no second I+J factor',
            coefficient_inequality='C0+k*C1/sigma <= k*(C0+C1/sigma), k>=1',
            factorial_inequality='k^k/k! <= e^k',
            geometric_bound='word_k <= Msrc*beta^k',
            normal_convergence_on_compact_subsets=True,
            termwise_derivatives_on_inner_compact_subsets=True,
            solution_norm='||v|| <= Msrc/(1-beta) <= 2*Msrc',
            uniqueness_binding_sha256=digest(canonical)),
        weighted_unknown='(F1/x^2,Uz1/x^2,K1/x^2,P1/x^2,F1x/x,Uz1x/x)',
        axis_transform='W=(rho*v1,rho*v2,rho*v3,rho*v4,x*v5,x*v6)',
        weighted_solution_norm_upper_majorant=_record(Mv),
        Msrc_upper_majorant=_record(Msrc),beta_upper_majorant=_record(beta),
        rho_bound_radius='9/2',rho_target_radius='4',weighted_prefix_degree=degree,
        raw_base_prefix_degree=common_degree,raw_odd_divided_by_x_prefix_degree=degree,
        maximal_value_prefix_degrees=dict(raw_base=16,weighted=15,raw_odd_divided_by_x=15),
        maximal_value_prefix_tail=dict(weighted_v=_record(value16_tail),
            raw_base_W1_W4=_record(target*value16_tail),
            raw_odd_W5_W6=_record(ctx.sqrt(target)*value16_tail)),
        exact_radial_factor=dict(numerator=factor.numerator,denominator=factor.denominator),
        prefix_is_exact_coefficient_polynomial=True,
        coefficient_intervals_are_not_selected_points=True,
        derivative_map='W5/x=2*d_rho F1, W6/x=2*d_rho Uz1; raw degreeN induces weighted degreeN-1',
        axial_domain=dict(real_center_interval_open=['-1','1'],
            local_center='any real Z0 in (-1,1)',
            source_geometry='Union of admitted local disks centered on the real interval; not a single origin disk',
            source_domain_report=domain_name,source_domain_sha256=sha(domain_path),
            solution_radius=_record(zradius),safe_solution_radius_lower=_record(safe_zradius),
            working_radius_endpoint_rounded_down=True,
            target_radius=_record(zgap),
            fresh_Cauchy_gap=_record(zgap),target_in_solution_disk=True),
        ordinary_Z0_Z3_tail_rows=rows,word_truncation_upper_majorants=word_tails,
        certified_values_are_upper_majorants_not_function_values=True,
        computed_word_partial_sum_or_source_point_values=False,
        approximate_coefficient_errors_combined=False,
        full_forward_response_function_identity=False,q_dependent_collar_bridge=False,
        six_core_lower_terms_removed=False,global_five_moment_closure=False,temporal_n_ge_2=False,
        actual_N_or_physical_scale_materialized=False,old_error_caps_modified=False,
        accepted_owner_fixture_recurrence_or_quadrature_replayed=False,
        limitations=['The certified Taylor polynomial uses exact source coefficients; numerical coefficient enclosure errors remain separate',
            'Value prefix reaches radial16, but the saved common ordinary Z0..3 prefix reaches radial14 only',
            'Rescaling this raw bound to forward/physical fields requires the bound scale functions and their ordinary Z jets',
            'The q-dependent collar and six signed source equality are not implied by this core certificate'])
    result['exact_interval_records_roundtripped']=_serialization_roundtrip(ctx,result)
    return result
