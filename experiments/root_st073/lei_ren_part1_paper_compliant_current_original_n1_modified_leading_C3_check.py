"""Independent exact finite-N velocity and own-mean physical identities."""
import json
from pathlib import Path
import time

import sympy as s
import lei_ren_part1_paper_compliant_current_original_n1_modified_leading_C3 as current


class SymbolicGraph:
    zero=s.Integer(0);one=s.Integer(1)
    def constant(self,q):return s.Rational(str(q))
    def add(self,*q):return s.Add(*q)
    def mul(self,*q):return s.Mul(*q)
    def sub(self,a,b):return a-b
    def quotient(self,a,b,proof):return a/b


def independent_fixture():
    z,y,N=s.symbols('z y N',positive=True);t=s.symbols('t',real=True)
    E=1+z+z*z+y;V=2-z+z*z+y*(1+z**3);A=(1+z*z)*y/7;B=(1+z+z*z+z**3)*y*y
    F=N*E*(s.exp(A/N)-1)
    delta_m=s.integrate(s.exp(t-y)*B.subs(y,t)/N,(t,0,y))
    jet=lambda q:current.c3.C3Function(*(s.diff(q,z,k) for k in range(4)))
    functions=current.corrected_functions(SymbolicGraph(),{'E':jet(E),'V':jet(V)},jet(F),jet(B),jet(N*delta_m),jet(0),N)
    truth=dict(E_N=E*s.exp(A/N),V_N=V+B/N,delta_E=E*(s.exp(A/N)-1),delta_V=B/N,
        own_mean_increment_at_right=delta_m,own_mean_increment_at_left=s.Integer(0))
    count=0
    for name,q in truth.items():
        for k,value in enumerate(current.c3.rows(functions[name])):
            assert s.simplify(value-s.diff(q,z,k))==0,(name,k);count+=1
    assert s.simplify(s.diff(delta_m,y)-(B/N-delta_m))==0
    assert s.simplify(delta_m.subs(y,0))==0
    # The rate-one memory is exactly a physical axial radial mean.
    R0,R,x=s.symbols('R0 R x',positive=True)
    direct=s.integrate(B.subs(y,s.log(x/R0))/N,(x,R0,R))/R
    assert s.simplify(direct-delta_m.subs(y,s.log(R/R0)))==0
    delta=s.Rational(1,5);L=1-delta*z*z;d=1-z*z
    Q=(2*z*B/N-(1-delta)*z*delta_m-d*s.diff(delta_m,z))/L
    assert s.simplify(L*(Q+s.diff(Q,y))-(1+delta)*z*B/N+d*s.diff(B/N,z)-2*z*s.diff(B/N,y))==0
    # Incoming history persists in a quiet next chart and remains divided by N once.
    h=s.symbols('h',positive=True)
    quiet=s.exp(-h)*delta_m
    assert s.simplify(s.diff(quiet,h)+quiet)==0 and quiet.subs(h,0)==delta_m
    return count


def rejected(fn):
    try:fn()
    except ValueError:return True
    raise AssertionError('Unsupported modified leading query admitted')


def run(owner,packets):
    began=time.monotonic();fixture=independent_fixture();owner.assert_graph()
    views={k:owner.report(v) for k,v in packets.items()}
    for view in views.values():
        assert all(not view[k] for k in current.OPEN)
        assert view['complete_modified_Q_or_pressure_tail_not_installed']
        assert view['phase_radial_and_mixed_derivatives_not_exported']
        assert view['endpoint_mean_increment_not_whole_cell_mean_oracle']
        assert view['frequency']['N_not_reselected_or_materialized']
        assert set(view['corrected_signed_velocity_and_own_mean_C3_functions'])=={
            'E_N','V_N','delta_E','delta_V','own_mean_increment_at_right','own_mean_increment_at_left'}
        assert all(len(v)==4 for v in view['corrected_signed_velocity_and_own_mean_C3_functions'].values())
        assert all(len(v)==4 for v in view['corrected_source_C3_magnitude_ranges'].values())
    for k in ('E_N','V_N'):
        for row in owner.functions['first_micro'][k].__dict__.values():assert 0<=row.node<len(owner.graph.nodes)
    quiet=owner.functions['O3_power'];primitive=owner.target.primitives['O3_power']
    assert current.c3.encoded(quiet['E_N'])==current.c3.encoded(primitive['source_roots']['E'])
    assert current.c3.encoded(quiet['V_N'])==current.c3.encoded(primitive['source_roots']['V'])
    token=packets['first_micro'];entry=owner._packets[id(token)]
    key='complete_modified_Q_or_pressure_tail_not_installed';saved=entry[1][key];entry[1][key]=False
    rejected(lambda:owner.report(token));entry[1][key]=saved
    saved=owner.functions['first_micro']['V_N'];owner.functions['first_micro']['V_N']=owner.functions['first_micro']['E_N']
    rejected(lambda:owner.report(token));owner.functions['first_micro']['V_N']=saved
    saved=entry[2]['primitive']['proof']['exact_same_source_P0'];entry[2]['primitive']['proof']['exact_same_source_P0']=False
    rejected(lambda:owner.report(token));entry[2]['primitive']['proof']['exact_same_source_P0']=saved
    rejected(lambda:owner.report(current.OriginalModifiedLeadingC3Packet()))
    rejected(lambda:owner.endpoint(('0','.5'),'reference_only_unknown'))
    rejected(lambda:owner.endpoint(('0','.25'),'first_micro'))
    saved=owner.N0;owner.N0+=1;rejected(lambda:owner.report(token));owner.N0=saved
    # A downstream recovery can append new nodes, while existing source nodes stay fixed.
    owner.graph.symbol('downstream_n1_recovery_probe');owner.report(token)
    owner.report(token)
    inputs=dict(owner.hashes);inputs[current.NAME]=current.sha(current.NAME);inputs[Path(__file__).name]=current.sha(Path(__file__).name)
    result=dict(all_passed=True,source_family=owner.family,input_hashes=inputs,
        checks=dict(independent_finite_N_velocity_and_own_mean_Z0_to_Z3_rows=fixture,
            physical_axial_mean_rate_one_and_radial_divergence_identity_verified=True,
            quiet_chart_preserves_nonzero_incoming_mean=True,
            two_same_family_current_micro_endpoint_queries_execute=True,
            same_original_N_C3_prefix_and_common_velocity_unit_retained=True,
            modified_signed_functions_separate_from_directed_magnitude_ranges=True,
            scope_function_and_source_range_mutations_and_foreign_queries_rejected=True,
            source_frequency_and_existing_graph_nodes_protected_while_child_composition_allowed=True,
            mixed_radial_phase_Q_pressure_tail_and_repair_not_admitted=True),
        **{current.GATE:True},**dict.fromkeys(current.OPEN,False),execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_bytes((json.dumps(current.c3.encoded(result),separators=(',',':'))+'\n').encode())
    print('PASS_MODIFIED_LEADING_C3',fixture,flush=True);return result
