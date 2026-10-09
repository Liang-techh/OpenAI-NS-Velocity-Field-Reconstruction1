"""Fresh partial functions, directed masses, kernel/ODE and domain checks."""
from dataclasses import replace
import gzip
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_partial_history_functions as current

ep=current.ep


def interval(c,row):
    return c.mpf((mp.mp.make_mpf(tuple(row['lower_exact_mpf_tuple'])),
                  mp.mp.make_mpf(tuple(row['upper_exact_mpf_tuple']))))


def overlap(a,b):
    lo,hi=ep(a);left,right=ep(b);assert max(lo,left)<=min(hi,right)


def exact_identities():
    x,a,b,d,Z=sy.symbols('x a b d Z',real=True);F=sy.Function('original_density')
    kernel=ODE=0
    for rate in current.transport.RATES.values():
        r=sy.Rational(rate.numerator,rate.denominator)
        assert sy.simplify(sy.exp(-r*(d-b))*sy.exp(-r*(b-x))-sy.exp(-r*(d-x)))==0;kernel+=1
        H=sy.exp(-r*b)*sy.Integral(sy.exp(r*x)*F(x,Z),(x,a,b))
        assert sy.simplify(sy.diff(H,b)-(F(b,Z)-r*H))==0;ODE+=1
    return dict(passed=True,exact_requested_endpoint_kernel_semigroup_identities=kernel,
        exact_original_endpoint_ODE_identities=ODE,
        geometry_endpoints_and_rates_Z_independent=True,
        numerical_interval_overlap_not_used_as_algebraic_identity_proof=True)


def functions(owner,manifest):
    c=owner.c;frames=[];comparisons=masses=zero=micro=0
    p=mp.mp.clone();p.dps=c.dps+150
    for saved in manifest['actual_radial_variable_history_frames']:
        a,b=map(sy.Rational,saved['exact_native_interval']);Z=saved['original_Z_exact'];N=saved['explicit_candidate_N']
        assert sy.Rational(Z)==sy.Rational(37,100) and N in (160,257)
        frame=owner.subinterval(str(a),str(b),Z=Z,N=N);frames.append(frame)
        with mp.workdps(c.dps+40):
            cursor=a
            for piece in frame.record['actual_partial_source_cell_records']:
                left,right=map(sy.Rational,piece['exact_restricted_cell']);pl,pr=map(sy.Rational,piece['exact_parent_cell'])
                assert cursor==left and pl<=left<right<=pr and right<=b;cursor=right
                assert piece['requested_kernel_endpoint']==str(b)
                assert piece['entire_source_and_true_phase_cover_restricted_conservatively']
                source=piece['original_live_source_record'];assert source['source_family']==owner.family
                if piece['chart']=='Rh_reference':
                    assert -5<=pl<pr<=0 and source['candidate_N']==N
                    assert piece['original_exact_N_dependent_coefficients_and_pressure_errors_retained']
                    assert piece['canonical_reference_factors_retained_until_combination']
                else:
                    assert 0<=pl<pr<=1 and tuple(map(sy.Rational,source['exact_Z_range']))==(frame.Z,frame.Z)
                    phase=piece['actual_original_phase_cover']
                    assert phase['explicit_candidate_N']==N and phase['source_graph_sha256']==owner.source_graph_sha256
                    assert piece['old_route_endpoint_masses_and_preintegrated_values_not_reused']
                    assert piece['original_quad_pressure_and_correlated_Z_errors_retained']
                width=right-left;distance=b-right
                wp=p.mpf(int(width.p))/int(width.q);dp=p.mpf(int(distance.p))/int(distance.q)
                for k,rate in current.transport.RATES.items():
                    rr=p.mpf(rate.numerator)/rate.denominator
                    target=wp if not rate else p.exp(-rr*dp)*(-p.expm1(-rr*wp))/rr
                    lo,hi=ep(piece['exact_positive_own_rate_masses'][k]);assert 0<lo<=target<=hi;masses+=1
                    micro+=int(width==sy.Rational(1,10**1000))
            assert cursor==b
            for k in current.transport.RATES:
                assert ep(frame.decay[k])==ep(interval(c,saved['directed_homogeneous_multipliers'][k]))
                for j in ('C0','Z'):
                    assert ep(frame.forcing[k,j])==ep(interval(c,saved['directed_partial_C0_Z_forcing'][k][j]));comparisons+=1
                    if a==b:assert ep(frame.forcing[k,j])==(0,0)
            if a==b:
                assert not frame.record['actual_partial_source_cell_records'];zero+=1
                assert all(ep(v)==(1,1) for v in frame.decay.values())
            assert frame.record['unknown_incoming_functions_retained_at_exact_left_endpoint']
            assert frame.record['separate_original_P0_and_P0_Z_unchanged']
            assert frame.record['reference_exact_N_dependent_order_minus1_and_minus2_functions_retained']
            assert frame.record['O2_full_original_DAG_density_hulls_at_actual_fixed_N_not_all_N_coefficients']
    assert micro==5
    return frames,dict(passed=True,radial_history_query_frames=len(frames),distinct_live_partial_frames=len({id(v) for v in frames}),
        genuine_partial_C0_Z_result_comparisons=comparisons,
        independent_higher_precision_positive_requested_endpoint_mass_checks=masses,
        exact_zero_width_identity_cases=zero,strict_positive_10_to_minus1000_subcell_mass_checks=micro,
        exact_source_partition_has_no_gaps_or_duplicate_mass=True,
        fresh_live_source_owners_not_saved_source_hydration=True)


def whole_and_semigroup(owner):
    c=owner.c;full_checks=prefix_checks=semigroup=0
    with mp.workdps(c.dps+40):
        for N in (160,257):
            ref,slope,live=owner.seed(sy.Rational(37,100),N)
            for index,cell in enumerate(owner.current.route):
                left=current.reference.symbolic_node(owner.current.provider,cell['lower'])[0]
                right=current.reference.symbolic_node(owner.current.provider,cell['upper'])[0]
                local=owner.subinterval(str(left),str(right),Z='.37',N=N);prefix=owner.prefix(str(right),Z='.37',N=N)
                for k,pair in cell['contributions'].items():
                    for j,label in (('C0','value'),('Z','Z')):
                        original=owner.current.integrate(owner.current.built['graph'].nodes[pair[label]],Z='.37',N=N)
                        overlap(local.forcing[k,j],original);full_checks+=1
                        overlap(prefix.forcing[k,j],live.forcing[index,k,j]);prefix_checks+=1
        for a,m,b in (('-5','-2.337','0'),('-5','0','.731'),('0','.137','.731'),('.12','.13','.15'),('.13','.14','1'),('.731','.9','1')):
            left=owner.subinterval(a,m,Z='.37',N=160);right=owner.subinterval(m,b,Z='.37',N=160);whole=owner.subinterval(a,b,Z='.37',N=160)
            for k in current.transport.RATES:
                overlap(left.decay[k]*right.decay[k],whole.decay[k])
                for j in ('C0','Z'):
                    overlap(right.decay[k]*left.forcing[k,j]+right.forcing[k,j],whole.forcing[k,j]);semigroup+=1
    return dict(passed=True,original_full_route_integral_enclosure_consistency_checks=full_checks,
        original_six_route_prefix_operator_consistency_checks=prefix_checks,
        local_semigroup_enclosure_consistency_checks=semigroup,
        overlap_is_supplementary_consistency_not_identity_or_containment=True)


def jets_and_guards(owner,manifest,frames):
    c=owner.c;jets=rejected=0
    for saved in manifest['actual_endpoint_y_Z_jet_records']:
        frame=owner.prefix(saved['exact_native_endpoint'],Z=saved['original_Z_exact'],N=saved['explicit_candidate_N'])
        live=owner.endpoint_jet(frame,chart=saved['endpoint_chart'])
        for k in current.transport.RATES:
            assert ep(live['directed_endpoint_y_derivative_incoming_multipliers'][k])==ep(interval(c,saved['directed_endpoint_y_derivative_incoming_multipliers'][k]))
            for j in ('C0','Z'):
                for field in ('actual_point_C0_Z_density','directed_endpoint_y_derivative_forcing'):
                    assert ep(live[field][k][j])==ep(interval(c,saved[field][k][j]))
                jets+=1
        assert live['conditional_on_unknown_left_history'] and live['higher_velocity_stress_join_not_admitted']
    def reject(call):
        nonlocal rejected
        try:call()
        except (TypeError,ValueError,ArithmeticError,NotImplementedError):rejected+=1
        else:raise AssertionError('Invalid original partial-history request accepted')
    frame=owner.prefix('.731',Z='.37',N=160);seam=owner.prefix(0,Z='.37',N=160)
    for call in (lambda:owner.prefix(2,Z='.37',N=160),lambda:owner.prefix(-6,Z='.37',N=160),
        lambda:owner.subinterval('.5','.4',Z='.37',N=160),lambda:owner.prefix('.7',Z=0,N=160),
        lambda:owner.prefix('.7',Z=2,N=160),lambda:owner.prefix('.7',Z='.37',N=159),
        lambda:owner.prefix('.7',Z='.37',N=True),lambda:owner.endpoint_jet(replace(frame)),
        lambda:owner.endpoint_jet(seam),lambda:owner.endpoint_jet(frame,chart='Rh_reference'),
        lambda:owner.apply_assumed_boundary(frame,{}),lambda:owner.apply_assumed_boundary(replace(frame),{}),
        lambda:owner.mass(sy.Rational(1),sy.Rational(0),sy.Rational(1),current.transport.RATES['m'])):
        reject(call)
    g=owner.current.built['graph'];node=next(iter(owner.current.integrals));row=g.nodes[node];old=row['upper'];row['upper']=row['lower']
    try:reject(lambda:owner.prefix('.731',Z='.37',N=160))
    finally:row['upper']=old
    owner.current.unchanged()
    incoming={role:c.mpf((str((i+1)/1000),str((i+2)/1000))) for i,role in enumerate(frame.forcing)}
    empty=owner.subinterval('.731','.731',Z='.37',N=160);action=owner.apply_assumed_boundary(empty,incoming)
    assert all(ep(action[role])==ep(value) for role,value in incoming.items())
    assert owner.mode!='original_source'
    return dict(passed=True,fresh_endpoint_y_and_mixed_yZ_jet_comparisons=jets,
        invalid_domain_frame_seam_graph_and_boundary_requests_rejected=rejected,
        exact_zero_width_preserves_explicit_nonzero_assumed_boundary=True,
        assumed_boundary_probe_not_original_source_data=True,
        reference_O2_endpoint_jet_seam_requires_explicit_one_sided_chart=True,
        original_graph_unchanged_and_cache_integrity_guard_retained=True)


@current.precise.phase.native.inlet.source_precision
def run():
    began=time.monotonic();raw=gzip.decompress((current.HERE/current.NAME).read_bytes());manifest=json.loads(raw)
    assert manifest[current.GATE]
    for name,digest in manifest['input_hashes'].items():assert current.sha(name)==digest,name
    flags=('actual_original_upstream_history_installed','full_scalar_control_evaluator_compatibility_installed',
        'full_17_chart_source_or_24_cell_integral_oracle_installed','actual_five_controls_installed',
        'actual_terminal_Z_function_closure_installed','current_whole_N_selected',*current.precise.phase.packets.OPEN)
    assert all(manifest[key] is False for key in flags)
    owner=current.OriginalPartialHistoryFunctions();assert owner.family==manifest['source_family']
    identities=exact_identities();frames,partial=functions(owner,manifest)
    checks=dict(exact_kernel_and_endpoint_ODE_identities=identities,
        fresh_radial_functions_source_partitions_and_stable_masses=partial,
        original_endpoints_and_local_semigroup_consistency=whole_and_semigroup(owner),
        original_endpoint_jets_and_partial_scope_guards=jets_and_guards(owner,manifest,frames))
    result=dict(all_passed=True,**{current.GATE:True},source_family=owner.family,source_graph_sha256=owner.source_graph_sha256,
        **checks,**dict.fromkeys(flags,False),compressed_report_raw_json_sha256=hashlib.sha256(raw).hexdigest(),
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began,
        scope='Fresh genuine radial-variable original source/phase range restriction, requested endpoint kernels, stable positive masses, conditional affine histories and true endpoint y/Z jets. Unknown upstream/P0 retained; no full history, whole-Z closure, global N, actual recursion or corrected NS.')
    (current.HERE/current.RECEIPT).write_bytes(json.dumps(current.precise.encode(result),indent=2).encode()+b'\n')
    print('Original arbitrary-radial partial histories and true endpoint jets PASS',flush=True)
    return result


if __name__=='__main__':run()
