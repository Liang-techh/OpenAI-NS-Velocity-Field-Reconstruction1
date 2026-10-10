"""Check sharp bound coordinates and same-source collected Poisson jets."""
import gzip
import json
from pathlib import Path
import time
from unittest.mock import patch
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_original_whole_Z_bridge_sharp_functions as current
import lei_ren_part1_paper_compliant_current_original_whole_Z_bridge_correlated_functions_check as previous_check
from lei_ren_part1_paper_compliant_current_original_bridge_macro_functions import MacroFlow


def encoded(value):return current.encode(current.serialized(value))


def independent_bound_coordinates():
    c = MPIntervalContext(); c.dps = 160
    p = mp.mp.clone(); p.dps = 210
    f = MacroFlow(c,-30,c.ln(4),c.ln(9),c.ln(5),c.ln(100)-c.ln(5),'.01')
    comparisons = rejects = 0
    for offset,powers in ((c.mpf(0),(0,0,0,0,0)),(c.mpf([12,13]),(0,.5,0,0,0))):
        base = f.factor(powers,offset)
        bound = base*2
        for coeff,shift in (([-1,1],0),([1,3],0),([-3,-1],0),([0,3],0),([-3,0],0),([0,0],0),
                (['.75','1.25'],c.ln(2)),(['-1.25','-.75'],c.ln(2)),([-1,1],c.mpf([-2,2]))):
            value = f.factor(powers,offset+shift)*c.mpf(coeff)
            got = current.bound_scale_intersection(value,bound)
            def physical(row):
                l,u = current.ep(row.coefficient)
                a,b = current.ep(row.scale.evaluate())
                return p.mpf(l)*p.exp(p.mpf(a if l>=0 else b)),p.mpf(u)*p.exp(p.mpf(b if u>=0 else a))
            l,u = physical(value); bl,bu = physical(bound); gl,gu = physical(got)
            wantl,wantu = max(l,-bu),min(u,bu)
            assert gl <= wantl+abs(wantl)*p.mpf('1e-150') and gu >= wantu-abs(wantu)*p.mpf('1e-150')
            comparisons += 2
        # Truly disjoint global intervals, rather than an unknown source cap.
        for coeff in ([3,4],[-4,-3]):
            value = f.factor(powers,offset+c.mpf(2))*c.mpf(coeff)
            try:current.bound_scale_intersection(value,bound)
            except ArithmeticError:rejects += 1
            else:raise AssertionError('Disjoint signed source interval accepted')
    huge = c.mpf('1e100')
    for signbox in ([-1,1],[0,1],[-1,0]):
        value = f.factor((0,0,0,0,0),huge)*c.mpf(signbox)
        got = current.bound_scale_intersection(value,f.scalar(2))
        l,u = current.ep(got.finite_interval(max_log=100))
        assert -2 <= l <= u <= 2
        assert l == -2 if signbox[0]<0 else l == 0
        assert u == 2 if signbox[1]>0 else u == 0
        comparisons += 2
    for coeff in ([1,2],[-2,-1]):
        value = f.factor((0,0,0,0,0),-huge)*c.mpf(coeff)
        got = current.bound_scale_intersection(value,f.scalar(2))
        assert not got.zero
        assert current.ep(got.record()['log_absolute_upper'])[1] < -500
        comparisons += 1
    for coeff in ([1,2],[-2,-1]):
        value = f.factor((0,0,0,0,0),huge)*c.mpf(coeff)
        try:current.bound_scale_intersection(value,f.scalar(2))
        except ArithmeticError:rejects += 1
        else:raise AssertionError('Huge disjoint source interval accepted')
    try:current.bound_scale_intersection(f.scalar(1),f.scalar(0))
    except ArithmeticError:rejects += 1
    else:raise AssertionError('Nonzero source accepted against exact zero bound')
    assert current.bound_scale_intersection(f.scalar(0),f.scalar(0)).zero
    return dict(passed=True,independent_signed_endpoint_comparisons=comparisons,
        disjoint_or_zero_bound_rejections=rejects,exp_1e100_clipped_in_bound_units=True,
        nontrivial_formal_basis_and_interval_bound_offset_checked=True,
        tiny_negative_log_tail_retained_not_exact_zero=True)


def independent_inverse():
    def fixture_adapter(source,qrows,dstar,phi,branch):
        # Existing scalar inverse fixtures prescribe arbitrary smooth q,
        # rather than q of the paper cutoff. Differentiate that fixture's
        # same q squared exactly; production queries use cutoff metadata.
        q2 = {current.C0:current.previous.switch.parameters.original.square(qrows[current.C0]),
            current.Z:qrows[current.C0]*qrows[current.Z]*2}
        source = dict(source,original_q_squared_C0_Z=q2,
            original_q_squared_Z_proof=dict(independent_arbitrary_smooth_q_fixture_product_rule=True))
        return current.sharp_inverse_branch(source,qrows,dstar,phi,branch)
    with patch.object(current.previous,'correlated_inverse_branch',fixture_adapter):
        inverse = previous_check.independent_inverse()
        guards = previous_check.flat_symmetry_and_magnitude_guards()
    return dict(passed=True,actual_sharp_inverse_scalar_references=inverse,
        flat_symmetry_magnitude_guards=guards,two_collected_Poisson_factors_in_actual_Z_compiler=True)


def independent_original_q_squared_Z():
    c = MPIntervalContext(); c.dps = 160
    p = mp.mp.clone(); p.dps = 210
    f = MacroFlow(c,-30,c.ln(4),c.ln(9),c.ln(5),c.ln(100)-c.ln(5),'.01')
    eta = p.mpf('.01'); az = p.mpf('.07'); dz = p.mpf('.11')
    compared = 0; families = []
    boxes = (('negative','-.3','-.1'),('negative_transition','-.001','.001'),
        ('transition','.002','.008'),('transition_flat','.008','.02'),('flat','.02','.03'))
    for label,left,right in boxes:
        aa = {current.C0:f.scalar('.8'),current.Z:f.scalar('.07')}
        Delta = f.scalar(c.mpf([left,right])); drow = f.scalar('.11')
        loga = c.ln(c.mpf('.8'))
        original = current.previous.switch.parameters.original.q_enclosure(aa[current.C0],Delta,c.ln(c.mpf('.01')),loga)
        q = original['q']; qr = {current.C0:q,current.Z:f.scalar(0 if q.zero else 1)}
        proof = dict(actual_a_axial5=[aa[current.C0],aa[current.Z]],actual_Delta_axial5=[Delta,drow],
            actual_positive_a=dict(source_log_lower=loga))
        source = dict(q=q,roots=dict(a=aa),original_full_source_quotients=proof)
        rows,record = current.original_q_squared_Z(source,qr)
        assert source['q'] is q and qr[current.C0] is q
        if not q.zero:
            assert record['first_transition_eta_over_eta_cancelled_before_interval_evaluation']
            assert record['no_y_derivative_rows_constructed']
        for point in (p.mpf(left),(p.mpf(left)+p.mpf(right))/2,p.mpf(right)):
            def exact(z):
                delta = point+dz*z; a = p.mpf('.8')+az*z
                x = 1-delta/eta
                sigma = p.mpf(0) if x<=0 else p.mpf(1) if x>=1 else 1/(1+p.exp(1/x**2-1/(1-x)**2))
                return sigma**2*(2*eta-delta)/(2*a)
            for order,key in enumerate((current.C0,current.Z)):
                want = p.diff(exact,p.mpf(0),order)
                lo,hi = current.ep(rows[key].finite_interval(max_log=1000))
                assert lo <= want <= hi,(label,order,want,lo,hi)
                compared += 1
        families.append(label)
    return dict(passed=True,independent_original_cutoff_q_squared_C0_Z_comparisons=compared,
        negative_transition_flat_and_both_mixed_branch_families=families,
        actual_q_C0_object_and_linear_q_Z_not_changed=True,
        independent_scalar_derivative_not_cap_derivative=True)


def forbidden(*args,**kwargs):
    raise AssertionError('Unchanged downstream integrations must not be rerun')


def run():
    began = time.monotonic()
    saved = json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    assert saved[current.GATE]
    for name,digest in saved['input_hashes'].items():
        if current.sha(name) != digest:raise ValueError('Changed source prerequisite: '+name)
    with mp.workdps(540):
        bounds = independent_bound_coordinates()
        inverse = independent_inverse()
        q2 = independent_original_q_squared_Z()
        owner = current.WholeZSharpBridgeFunctions()
        assert owner.identity == saved['source_family'] and owner.N == saved['candidate_N']
        windows_checked = target_pairs = 0
        for row in saved['actual_sharp_bridge_and_Rc_cells']:
            ends = tuple(row['exact_Z_cell'])
            f = owner.owner(ends).flow
            decode = lambda v:current.previous.relative.restore_row(f,v)
            pairs = lambda rows:{key:[decode(v) for v in pair] for key,pair in rows.items()}
            incoming = pairs(row['actual_correlated_bridge_transport']['actual_R100_signed_correction_C0_Z'])
            old = owner.cells[ends]['actual_source_owned_zero_inlet_to_Rc_signed_stages']
            lookup = {}
            for stage,key in (('switch','actual_signed_switch_windows'),('long','actual_signed_long_windows'),('outer','actual_signed_outer_windows')):
                for window in old[stage][key]:lookup[(stage,window['actual_chart'])]=window
            lookup[('patch','actual_patch')] = dict(
                original_incoming_own_rate_memory=old['patch']['original_incoming_own_rate_memory'],
                actual_signed_local_driver_C0_Z=old['patch']['actual_signed_local_driver_at_Rh_C0_Z'])
            for window in row['actual_downstream_retransport_windows']:
                assert encoded(incoming) == window['actual_incoming_correction_C0_Z']
                prior = lookup[(window['stage'],window['actual_chart'])]
                assert prior['actual_signed_local_driver_C0_Z'] == window['unchanged_independent_local_driver_C0_Z']
                assert prior['original_incoming_own_rate_memory'] == window['original_incoming_own_rate_memory']
                memory = {key:decode(v) for key,v in prior['original_incoming_own_rate_memory'].items()}
                local = pairs(prior['actual_signed_local_driver_C0_Z'])
                incoming = {key:[memory[key]*incoming[key][i]+local[key][i] for i in range(2)] for key in current.RATES}
                assert encoded(incoming) == window['actual_exit_signed_correction_C0_Z']
                windows_checked += 1
            assert encoded(incoming) == row['actual_Rc_signed_correction_C0_Z']
            target = current.previous.relative.normalized_targets(f,incoming,
                [decode(v) for v in old['outer']['actual_positive_Rc_normalized_amplitude_C0_Z']],
                f.factor((0,0,0,0,0),owner.source.logmu),owner.source.logmu,owner.N)
            for key,value in target.items():assert encoded(value) == row[key],key
            target_pairs += 5
            assert all(not row[key] for key in current.OPEN)
            assert not row['actual_terminal_controls_installed'] and not row['actual_global_frequency_admitted']
            print('Checked sharp bridge new incoming propagation and Rc targets: '+str(ends),flush=True)
        with patch.object(current.graph.WholeZFullSignedTransport,'full_transport',forbidden), \
                patch.object(current.graph.WholeZFullSignedTransport,'outer_transport',forbidden):
            live = owner.query(current.CELLS[0],'first_micro',(3,4),(1,1))
        expected = saved['actual_sharp_bridge_and_Rc_cells'][0]['actual_correlated_bridge_transport'] \
            ['actual_signed_bridge_windows'][0]['actual_signed_cells'][3]['actual_signed_primitive_C0_Z_phi']
        assert encoded(live['actual_signed_primitive_C0_Z_phi']) == expected
        for alternative in live['actual_signed_inverse_function_alternatives']:
            proof = alternative['original_inverse_and_Z_function_proof']
            if proof['geometry'] != 'flat':
                assert proof['correlated_inverse_identity_binding']['exact_collected_geometry_replacement_counts'] == [1,1]
                assert proof['same_function_active_flat_magnitude_theorem']['tiny_tails_not_moved_to_a_dominating_unrelated_scale']
                assert proof['collected_original_Poisson_geometry_Z_factors']['original_q_p2_Z_chain_retained']
                assert proof['correlated_inverse_identity_binding']['exact_original_q_square_replacement_count'] == 5
                assert proof['original_direct_q_squared_Z_cover']['first_transition_eta_over_eta_cancelled_before_interval_evaluation']
        hashes = dict(saved['input_hashes'])
        for name in (current.NAME,Path(__file__).name):current.bind(hashes,name,current.sha(name))
        receipt = dict(all_passed=True,**{current.GATE:True},source_family=owner.identity,candidate_N=owner.N,
            independent_signed_bound_coordinate_fixture=bounds,independent_collected_inverse_fixture=inverse,
            independent_original_q_squared_Z_fixture=q2,
            exact_new_incoming_transport_windows=windows_checked,exact_Rc_target_pairs=target_pairs,
            actual_first_active_bridge_box_exact_MPI_replayed=True,
            unchanged_downstream_source_integrations_not_repeated=True,
            actual_terminal_controls_installed=False,actual_global_frequency_admitted=False,
            physical_original_exterior_five_targets_closed=False,**dict.fromkeys(current.OPEN,False),
            input_hashes=hashes,execution_seconds=time.monotonic()-began)
        (current.HERE/current.RECEIPT).write_text(json.dumps(current.encode(receipt),indent=2)+'\n',encoding='utf8')
    print('PASS: sharp magnitude and collected Poisson Z bridge',flush=True)
    return receipt


if __name__ == '__main__':run()
