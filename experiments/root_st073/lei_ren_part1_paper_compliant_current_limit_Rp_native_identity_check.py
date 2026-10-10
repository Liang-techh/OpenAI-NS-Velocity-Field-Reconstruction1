"""Check actual recipe projection, pressure provenance and pulse-unit calculus."""
import json
import time
from pathlib import Path
import sympy as sy
import lei_ren_part1_paper_compliant_current_limit_Rp_native_identity as source
from lei_ren_part1_paper_candidate_pressure_function import q_jets


def independent_frame_calculus():
    current=source.current;g=current.source.FunctionTransportGraph();Z=source.Z;S=source.S
    U,M,H,K,EZ,EQ,Pin=sy.symbols('U M H K EZ EQ Pin',real=True);q=1+Z**2
    functions=dict(E=U/q,m=M*Z/S,h=H/q,k=K*Z/(S*q),e=EZ*Z**2+EQ/q**2,p=Pin/q**2,
        P0=sy.Function('P0')(Z))
    symbols={}
    def pair(key):
        symbols[key]=functions[key];symbols[key+'_Z']=sy.diff(functions[key],Z)
        return current.source.C1Function(g.symbol(key),g.symbol(key+'_Z'))
    frame=current.pulse_input(g,pair('E'),{k:pair(k) for k in ('m','h','k','e','p')},pair('P0'))
    def decode(ref):
        row=g.nodes[ref if type(ref) is int else ref.node];op=row['operation']
        if op=='exact_rational':return sy.Rational(row['numerator'],row['denominator'])
        if op=='bound_variable':return symbols[row['name']]
        if op=='sum':return sum(decode(n) for n in row['arguments'])
        if op=='product':return sy.prod(decode(n) for n in row['arguments'])
        if op=='negative':return -decode(row['argument'])
        if op=='positive_quotient':return decode(row['numerator'])/decode(row['denominator'])
        raise AssertionError(row)
    expected=dict(u=U/q,m1=M*Z*q/(S*U),m2=K*Z*q/(S*U**2),X=H/U,
        energy=(EZ*Z**2*q**2+EQ)/U**2,Mp=Pin/q**2,P0=functions['P0'],pressure=functions['P0']+Pin/q**2)
    for key,expr in expected.items():
        assert sy.cancel(decode(frame[key].value)-expr)==0,key
        assert sy.cancel(decode(frame[key].Z)-sy.diff(expr,Z))==0,key+'_Z'
    class ExactContext:
        mpf=staticmethod(sy.Rational)
    rows=q_jets(ExactContext,Z,6)
    for n,row in enumerate(rows):assert sy.cancel(row-sy.diff(q**-2,Z,n)/sy.factorial(n))==0,n
    return dict(actual_current_pulse_input_graph_value_Z_identities=16,
        actual_q_jet_recurrence_order6_true_derivative_identities=7,
        raw_m_k_inverse_Pstar_applied_once=True,independent_pressure_function_not_reset=True)


def run():
    began=time.monotonic();current,old,hashes=source.load_inputs()
    report=json.loads((source.HERE/source.NAME).read_bytes())
    assert report[source.GATE] and report['all_passed'] and report['source_family']==current['source_family']
    for name,digest in report['input_hashes'].items():assert source.sha(name)==digest,name
    projection=source.defining_source_proof()
    assert report['underlying_function_projection']==projection
    assert len(projection['exact_function_identities'])==45
    assert report['current_P0_source_binding']==source.pressure_source_proof(current,old)
    assert report['actual_native_frame_identity']==source.native_frame_proof(projection)
    assert all(row['value_identity'] and row['true_Z_identity'] and row['no_range_value_substitution']
        for row in projection['exact_function_identities'].values())
    calculus=independent_frame_calculus()
    for key in ('actual_native_O4_constructor_consumes_current_limit_frame','actual_numeric_point_values_installed',
        'old_numeric_mu_enclosure_equal_to_exact_mu','patched_Rh_functional_join_proved',
        'physical_original_exterior_five_targets_closed','actual_temporal_scale_recursion_installed'):
        assert report[key] is False,key
    hashes.update(report['input_hashes']);hashes[source.NAME]=source.sha(source.NAME)
    hashes[Path(__file__).name]=source.sha(Path(__file__).name)
    receipt=dict(**{source.GATE:True},all_passed=True,source_family=report['source_family'],
        exact_underlying_outer_identities_checked=45,exact_native_value_Z_identities_checked=12,
        current_to_native_P0_defining_function_identity_bound=True,
        independent_actual_graph_and_pressure_recurrence_calculus=calculus,
        no_numeric_cover_equality_or_inner_Rh_seam_inferred=True,
        actual_O4_constructor_installation_not_claimed=True,input_hashes=hashes,
        execution_seconds=time.monotonic()-began)
    (source.HERE/source.RECEIPT).write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8')
    print('PASS_CURRENT_LIMIT_RP_NATIVE_IDENTITY 45 outer / 12 native / analytic P0 / actual graph calculus',flush=True)
    return receipt


if __name__=='__main__':run()
