"""Independent scalar identities and narrow correlated-bridge replay.

No unchanged downstream source integration is repeated. Check actual new
inverse arithmetic against the original scalar loop, then exact incoming
propagation, original local driver hydration and target normalization.
"""
import gzip
import json
from pathlib import Path
import time
from types import SimpleNamespace
from unittest.mock import patch
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_original_whole_Z_bridge_correlated_functions as current
import lei_ren_part1_paper_compliant_current_original_Rm_conditioned_phase_density_check as scalar_check
from lei_ren_part1_paper_compliant_current_original_bridge_macro_functions import MacroFlow


def encoded(value):return current.encode(current.serialized(value))


def independent_inverse():
    def tested(source,qr,dstar,phi):
        a,E,t = (source['roots'][key] for key in ('a','E','t0'))
        C = {current.C0:a[current.C0]*E[current.C0],
            current.Z:a[current.Z]*E[current.C0]+a[current.C0]*E[current.Z]}
        B = {current.C0:-C[current.C0]*t[current.C0],
            current.Z:-C[current.Z]*t[current.C0]-C[current.C0]*t[current.Z]}
        source = dict(source,original_C_B_E_C0_Z=dict(C=C,B=B,E=E),original_eta_log=E[current.C0].ctx.mpf(0))
        _,branches,_ = current.backend.signed.signed_u_branches(source,dstar)
        rows = [current.correlated_inverse_branch(source,qr,dstar,phi,b) for b in branches]
        assert all(row['values'] is not None for row in rows)
        return dict(values={name:current.graph.switch_backend.directed_dominating_union(
            [row['values'][name] for row in rows]) for name in current.backend.OUTPUTS},
            record=dict(status='enclosed',geometry=rows[0]['record']['geometry']))
    proxy = SimpleNamespace(**{**vars(scalar_check.current),'Z_FIRST':tested})
    with patch.object(scalar_check,'current',proxy),mp.workdps(180):
        evidence = scalar_check.scalar_fixtures()
    assert len(evidence) == 4 and all(row['passed'] for row in evidence)
    return dict(passed=True,original_scalar_loop_inverse_Z_phi_and_five_density_fixtures=evidence,
        positive_negative_Mobius_small_r_and_zero_crossing_with_nonzero_source_Z=True,
        actual_explicit_phi_terms_and_implicit_inverse_chain_checked=True,
        bounded_magnitude_intersection_in_same_actual_function_evaluation=True)


def independent_quotients():
    c = MPIntervalContext(); c.dps = 180
    p = mp.mp.clone(); p.dps = 230
    f = MacroFlow(c,-30,c.ln(4),c.ln(9),c.ln(5),c.ln(100)-c.ln(5),'.01')
    z = p.mpf('.37')
    a = lambda v:p.mpf('.8')+p.mpf('.07')*v
    b = lambda v:p.mpf('.2')-p.mpf('.05')*v
    E = lambda v:p.mpf('1.3')+p.mpf('.04')*v
    eta = p.mpf('.01')
    pairs = lambda fun:[f.scalar(c.mpf(p.nstr(p.diff(fun,z,n)/p.factorial(n),225))) for n in range(6)]
    aa,ee,cc,bb = [pairs(fun) for fun in (a,E,lambda v:a(v)*E(v),lambda v:b(v)*E(v))]
    sectors = {key:[f.scalar(0)]*4 for key in ('theta_linear','theta_quadratic','axial_linear','axial_quadratic')}
    proxy = SimpleNamespace(flow=f,c=c,Pstar=f.scalar(2),source_radius=f.scalar(3))
    rec = dict(actual_generic_source_numerators=dict(E=ee,C=cc,B=bb),full_signed_inertial_sectors_axial4=sectors)
    roots,qr,proof = current.CORRELATED_QUOTIENTS(proxy,rec,aa,c.ln(c.mpf('.01')))
    delta = lambda v:a(v)+b(v)**2/a(v)-2
    q = lambda v:p.sqrt((2*eta-delta(v))/(2*a(v)))
    comparisons = 0
    for rows,fun in ((proof['actual_t0_axial5'],lambda v:-b(v)/a(v)),
            (proof['actual_kappa_axial5'],lambda v:delta(v)+2),(proof['actual_Delta_axial5'],delta)):
        for n,row in enumerate(rows):
            want = p.diff(fun,z,n)/p.factorial(n)
            lo,hi = current.ep(row.finite_interval(max_log=1000))
            assert lo <= want <= hi,(n,fun)
            comparisons += 1
    for n,key in enumerate((current.C0,current.Z)):
        lo,hi = current.ep(qr[key].finite_interval(max_log=1000))
        assert lo <= p.diff(q,z,n) <= hi
        comparisons += 1
    assert all(v.zero for v in proof['full_signed_p2_axial4'])
    assert roots['roots']['E'][current.C0] is ee[0]
    return dict(passed=True,independent_C_B_E_quotient_and_q_derivative_rows=comparisons,
        nonzero_C_B_E_a_Z_with_original_Taylor_coefficients=True,
        unchanged_signed_inertial_and_one_radius_recipe=True)


def flat_symmetry_and_magnitude_guards():
    c = MPIntervalContext(); c.dps = 160
    f = MacroFlow(c,-30,c.ln(4),c.ln(9),c.ln(5),c.ln(100)-c.ln(5),'.01')
    pair = lambda x,z:{current.C0:f.scalar(x),current.Z:f.scalar(z)}
    roots = dict(a=pair('.8','.07'),t0=pair('-.25','.01'),E=pair('1.3','.04'),p2=pair('.2','.07'))
    C = pair('1.04','.123'); B = pair('.26','.02035')
    qc = c.sqrt(c.mpf('1.17')/c.mpf('1.6'))
    q = pair(qc,qc*(-c.mpf('.070375')/c.mpf('1.17')-c.mpf('.07')/c.mpf('.8'))/2)
    source = dict(roots=roots,q=q[current.C0],original_C_B_E_C0_Z=dict(C=C,B=B,E=roots['E']),original_eta_log=c.ln(c.mpf('.01')))
    _,branches,_ = current.backend.signed.signed_u_branches(source,c.ln(c.mpf('.1')))
    rows = 0
    for phi in ('0','.5','1'):
        for branch in branches:
            got = current.correlated_inverse_branch(source,q,c.ln(c.mpf('.1')),c.mpf(phi),branch)
            assert got['values'] is not None
            for key in ('A','B_over_Pstar','A_Z','B_Z_over_Pstar'):
                assert got['values'][key].zero
                rows += 1
    flatq = pair(0,0)
    flat = dict(source,q=flatq[current.C0])
    _,bb,_ = current.backend.signed.signed_u_branches(flat,c.ln(c.mpf('.1')))
    got = current.correlated_inverse_branch(flat,flatq,c.ln(c.mpf('.1')),c.mpf('.137'),bb[0])
    assert all(v.zero for v in got['values'].values())
    bad = dict(source,original_eta_log=c.mpf(1))
    try:current.correlated_inverse_branch(bad,q,c.ln(c.mpf('.1')),c.mpf('.137'),branches[0])
    except ArithmeticError:pass
    else:raise AssertionError('Magnitude theorem accepted eta>1')
    return dict(passed=True,exact_symmetry_zero_rows=rows,exact_flat_rows=6,
        actual_eta_magnitude_precondition_rejection=True)


def forbidden(*args,**kwargs):
    raise AssertionError('Unchanged downstream source integration is outside this narrow check')


def run():
    began = time.monotonic()
    saved = json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    assert saved[current.GATE]
    for name,digest in saved['input_hashes'].items():
        if current.sha(name) != digest:raise ValueError('Changed source prerequisite: '+name)
    with mp.workdps(540):
        quotient = independent_quotients()
        inverse = independent_inverse()
        guards = flat_symmetry_and_magnitude_guards()
        owner = current.WholeZCorrelatedBridgeFunctions()
        assert owner.identity == saved['source_family'] and owner.N == saved['candidate_N']
        counts = dict(exact_retransport_windows=0,exact_Rc_target_pairs=0,same_function_magnitude_proofs=0)
        for row in saved['actual_correlated_bridge_and_Rc_cells']:
            ends = tuple(row['exact_Z_cell'])
            f = owner.owner(ends).flow
            decode = lambda v:current.relative.restore_row(f,v)
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
                counts['exact_retransport_windows'] += 1
            assert encoded(incoming) == row['actual_Rc_signed_correction_C0_Z']
            target = current.relative.normalized_targets(f,incoming,
                [decode(v) for v in old['outer']['actual_positive_Rc_normalized_amplitude_C0_Z']],
                f.factor((0,0,0,0,0),owner.source.logmu),owner.source.logmu,owner.N)
            for key,value in target.items():assert encoded(value) == row[key],key
            counts['exact_Rc_target_pairs'] += 5
            for chart in row['actual_correlated_bridge_transport']['actual_signed_bridge_windows']:
                for cell in chart['actual_signed_cells']:
                    assert set(cell['actual_signed_primitive_C0_Z_phi']) == set(current.backend.OUTPUTS)
            assert not row['actual_terminal_controls_installed'] and not row['actual_global_frequency_admitted']
            assert all(not row[k] for k in current.OPEN)
            print('Checked actual new incoming propagation and Rc targets: '+str(ends),flush=True)
        # One live first-active source box covers all conditional geometries
        # and the new callback, without repeating all 228 original cells.
        with patch.object(current.graph.WholeZFullSignedTransport,'full_transport',forbidden), \
                patch.object(current.graph.WholeZFullSignedTransport,'outer_transport',forbidden):
            live = owner.query(current.CELLS[0],'first_micro',(3,4),(1,1))
        got = encoded(live['actual_signed_primitive_C0_Z_phi'])
        expected = saved['actual_correlated_bridge_and_Rc_cells'][0]['actual_correlated_bridge_transport'] \
            ['actual_signed_bridge_windows'][0]['actual_signed_cells'][3]['actual_signed_primitive_C0_Z_phi']
        assert got == expected
        for alternative in live['actual_signed_inverse_function_alternatives']:
            proof = alternative['original_inverse_and_Z_function_proof']
            if proof['geometry'] != 'flat':
                bound = proof['same_function_active_flat_magnitude_theorem']
                assert bound['Z_derivative_bounds_unchanged_by_magnitude_theorem'] and bound['derivative_of_cap_not_used']
                before = proof['actual_identity_primitive_before_magnitude_intersection']
                assert before['A_Z'] == proof['original_A_B_first_derivative_enclosures']['A_Z']
                assert before['B_Z_over_Pstar'] == proof['original_A_B_first_derivative_enclosures']['B_Z_over_Pstar']
                counts['same_function_magnitude_proofs'] += 1
        assert counts['same_function_magnitude_proofs'] > 0
        hashes = dict(saved['input_hashes'])
        for name in (current.NAME,Path(__file__).name):current.bind(hashes,name,current.sha(name))
        receipt = dict(all_passed=True,**{current.GATE:True},source_family=owner.identity,candidate_N=owner.N,
            independent_correlated_quotient_fixture=quotient,independent_inverse_primitive_fixture=inverse,
            flat_symmetry_magnitude_guards=guards,
            narrow_actual_source_and_retransport_counts=counts,
            actual_live_first_active_bridge_box_exact_MPI_replayed=True,
            unchanged_downstream_integrations_not_repeated=True,
            actual_terminal_controls_installed=False,actual_global_frequency_admitted=False,
            physical_original_exterior_five_targets_closed=False,**dict.fromkeys(current.OPEN,False),
            input_hashes=hashes,execution_seconds=time.monotonic()-began)
        (current.HERE/current.RECEIPT).write_text(json.dumps(current.encode(receipt),indent=2)+'\n',encoding='utf8')
    print('PASS: correlated bridge functions and actual Rc retransport',flush=True)
    return receipt


if __name__ == '__main__':run()
