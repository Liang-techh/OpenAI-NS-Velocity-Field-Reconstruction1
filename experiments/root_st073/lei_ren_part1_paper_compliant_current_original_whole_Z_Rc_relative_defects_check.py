"""Check new Rc subtraction/normalization with accepted graph reuse."""
import gzip
import json
from pathlib import Path
import time
from types import SimpleNamespace
from unittest.mock import patch
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_original_whole_Z_Rc_relative_defects as current
from lei_ren_part1_paper_compliant_current_original_bridge_macro_functions import MacroFlow


def encoded(value):
    return current.encode(current.serialized(value))


def forbidden(*args, **kwargs):
    raise AssertionError('Unchanged full graph must not be recomputed by a normalization check')


def contains(row, expected):
    c = row.ctx
    ordinary = row.coefficient*c.exp(row.scale.evaluate())
    lo,hi = current.ep(ordinary)
    if not lo <= expected <= hi:
        raise ArithmeticError('Independent analytic normalization outside directed enclosure')


def independent_fixture():
    c = MPIntervalContext(); c.dps = 180
    p = mp.mp.clone(); p.dps = 230
    f = MacroFlow(c,-30,c.ln(4),c.ln(9),c.ln(5),c.ln(100)-c.ln(5),'.01')
    P0 = [f.scalar('.47')]+[f.scalar(0)]*5
    op = SimpleNamespace(flow=f,c=c,P0=P0,Rm_factor=f.scalar(5))
    owner = current.WholeZRcRelativeDefects.__new__(current.WholeZRcRelativeDefects)
    owner.c,owner.N,owner.identity,owner.binding = c,160,{'independent_fixture':True},current.reference_binding()
    owner.weights = current.repair.fresh_weights(c,c.mpf([0,current.ep(c.mpf(1)/6)[1]]))
    owner.matrix = current.repair.fresh_linear_inverse(c,c.mpf([0,current.ep(c.mpf(1)/6)[1]]),owner.weights)
    logmu = c.ln(c.mpf('.07'))
    calls = []
    def packet(ends):
        calls.append(ends)
        z = c.mpf(ends[0]); A = [f.scalar(2+z*z),f.scalar(2*z)]
        D = dict(m=[f.scalar('.13')+f.scalar('.02')*z,f.scalar('.02')],
            h=[f.scalar('.04')-f.scalar('.03')*z+f.scalar('.01')*z*z,f.scalar('-.03')+f.scalar('.02')*z],
            k=[f.scalar('.5')+f.scalar('.1')*z*z,f.scalar('.2')*z],
            e=[f.scalar('.1')+f.scalar('.04')*z,f.scalar('.04')],
            p=[f.scalar('.02')-f.scalar('.01')*z*z,f.scalar('-.02')*z])
        return dict(source_identity=owner.identity,candidate_N=160,exact_Z_cell=list(ends),
            actual_same_source_signed_C0_Z_integral_graph_through_Rc=True,
            exact_common_P0_axial5=P0,exact_Rc_radius=op.Rm_factor*f.factor((0,.5,0,0,0),9),
            actual_Rc_signed_correction_C0_Z=D,
            actual_source_owned_zero_inlet_to_Rc_signed_stages=dict(outer=dict(
                actual_Rc_signed_correction_C0_Z=D,actual_positive_Rc_normalized_amplitude_C0_Z=A)),
            actual_original_Rc_background_C0_Z=object(),
            actual_Rc_signed_complete_own_history_C0_Z=object())
    owner.graph = SimpleNamespace(owner=lambda ends:op,full_transport=packet,source=SimpleNamespace(logmu=logmu))
    owner.cells = {('-.75','-.75'):None,('0','0'):None,('.625','.625'):None}
    def functions(z):
        A = 2+z*z
        d = dict(m=p.mpf('.13')+p.mpf('.02')*z,h=p.mpf('.04')-p.mpf('.03')*z+p.mpf('.01')*z*z,
            k=p.mpf('.5')+p.mpf('.1')*z*z,e=p.mpf('.1')+p.mpf('.04')*z,p=p.mpf('.02')-p.mpf('.01')*z*z)
        return dict(M=d['m']/A,I=d['h']/A,S=d['e']/(A*A),Cp=d['p']/(A*A),
            **{current.repair.ROWS[1]:(d['k']/A**2-d['m']/A)/p.mpf('.07')})
    comparisons = 0
    for ends in owner.cells:
        row = owner.query(ends); z = p.mpf(ends[0]); expect = functions(z)
        for key in current.repair.ROWS:
            derivative = p.diff(lambda q:functions(q)[key],z)
            for i,value in enumerate((expect[key],derivative)):
                contains(row['actual_fixed_N_normalized_target_C0_Z'][key][i],value)
                contains(row['actual_fixed_N_N_scaled_target_C0_Z'][key][i],value*160)
                comparisons += 2
        R = 10*p.exp(9); S = p.mpf(2)
        units = dict(m=R*S,h=p.sqrt(2)*R**p.mpf('1.5')*S,
            k=p.sqrt(2)*R**p.mpf('1.5')*S*S,e=R*S*S,p=S*S)
        D = packet(ends)['actual_Rc_signed_correction_C0_Z']
        for key in current.RATES:
            for i in range(2):
                ordinary = D[key][i].coefficient
                lo,hi = current.ep(ordinary)
                center = (p.mpf(lo)+p.mpf(hi))/2
                contains(row['actual_relative_Rc_physical_five_moment_differences_C0_Z'][key][i],center*units[key])
                comparisons += 1
    # Each query calls the actual source graph once; the additional three
    # packet calls above supply only independent fixture physical units.
    if len(calls) != 6:raise AssertionError('Live callable failed to consume the actual graph once per query')
    try: owner.query(('1','1'))
    except ValueError: pass
    else: raise AssertionError('Unadmitted point-label fallback accepted')
    return dict(passed=True,independent_analytic_rows=comparisons,
        nonzero_amplitude_Z_and_all_five_defect_Z_rows_exercised=True,
        actual_query_calls_graph_once=True,complete_background_rows_not_subtracted_as_independent_hulls=True,
        invalid_domain_rejected=True)


def run():
    began = time.monotonic()
    report = json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    with mp.workdps(540):
        fixture = independent_fixture()
        with patch.object(current.graph.WholeZFullSignedTransport,'full_transport',forbidden):
            owner = current.WholeZRcRelativeDefects()
            if report['source_family'] != owner.identity or report['candidate_N'] != owner.N or not report[current.GATE]:
                raise ValueError('Current relative defect source/report admission changed')
            replay = []
            for ends in current.CELLS:
                replay.append(encoded(owner.accepted_evaluation(ends)))
                print('Checked actual Rc same-input C1 defects/targets: '+str(ends),flush=True)
        if replay != report['actual_Rc_relative_defect_cells']:
            raise ArithmeticError('Actual fixed-N normalized function evaluation exact tuples changed')
        if encoded(owner.weights) != report['original_exact_bump_weight_enclosures'] \
                or encoded(owner.matrix) != report['original_exact_divided_inverse_enclosures']:
            raise ArithmeticError('Fresh integral/inverse constant enclosures changed')
        if current.exact_quotient_theorem() != report['exact_quotient_Z_theorem']:
            raise ArithmeticError('Independent exact quotient theorem changed')
        for row in replay:
            if row['actual_terminal_controls_installed'] or row['actual_global_frequency_admitted'] \
                    or row['physical_original_exterior_five_targets_closed'] or any(row[k] for k in current.OPEN):
                raise AssertionError('Relative defect adapter overclaims remaining reconstruction')
        for name,digest in report['input_hashes'].items():
            if current.sha(name) != digest:raise ValueError('Changed accepted source prerequisite: '+name)
        hashes = dict(report['input_hashes'])
        for name in (current.NAME,Path(__file__).name):current.bind(hashes,name,current.sha(name))
        receipt = dict(all_passed=True,**{current.GATE:True},source_family=owner.identity,candidate_N=owner.N,
            independent_normalization_and_live_callable_fixture=fixture,
            actual_four_Z_cells_exact_fixed_N_normalization_replayed=True,
            accepted_function_evaluation_not_a_support_cap_or_selected_value=True,
            unchanged_full_graph_replay_forbidden_in_this_narrow_check=True,
            fresh_bump_integral_and_divided_inverse_constants_checked=True,
            exact_quotient_Z_identities_checked=True,source_P0_radius_mu_family_and_N_guards_checked=True,
            actual_Rc_relative_defect_functions_installed=True,actual_terminal_controls_installed=False,
            actual_global_frequency_admitted=False,physical_original_exterior_five_targets_closed=False,
            **dict.fromkeys(current.OPEN,False),input_hashes=hashes,execution_seconds=time.monotonic()-began)
        (current.HERE/current.RECEIPT).write_text(json.dumps(current.encode(receipt),indent=2)+'\n',encoding='utf8')
    print('PASS: actual whole-Z fixed-N Rc relative C1 defects/targets',flush=True)
    return receipt


if __name__ == '__main__': run()
