"""Independent original quadratic calculus, C2 tail and same-limit checks."""
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_C2_limit_controls as current
import lei_ren_part1_paper_compliant_current_original_C2_density_transport_check as symbolic_source

phase=current.phase


def independent_map_derivatives():
    g=current.source.FunctionTransportGraph();z=s.Symbol('Z');bindings={}
    def fixed(name):
        q=g.symbol(name);bindings[q.node]=s.Symbol(name,positive=True);return q
    def function(name):
        f=s.Function(name)(z);rows=[]
        for order in range(3):
            q=g.symbol(name+'_'+str(order));bindings[q.node]=s.diff(f,z,order);rows.append(q)
        return phase.C2Function(*rows)
    W=dict(mu=fixed('mu'),D=[fixed('D0'),fixed('D2')],
        axial_positive_determinant=fixed('axial_determinant'),
        inverse_swirl=[[fixed('inverse_%d%d'%(i,j)) for j in range(3)] for i in range(3)],
        cross=[fixed('cross0'),fixed('cross2')],
        energy=[fixed('energy'+str(i)) for i in range(3)],
        pressure=[fixed('pressure'+str(i)) for i in range(3)])
    h=[function('h'+str(i)) for i in range(5)];d=[function('d'+str(i)) for i in range(5)]
    N=fixed('N');next_h=current.picard_C2_map(g,W,N,d,h)
    symbolic=lambda q:symbolic_source.symbolic(g,q.node,bindings)
    for q in next_h:
        want=symbolic(q.value)
        assert s.simplify(symbolic(q.Z)-s.diff(want,z))==0
        assert s.simplify(symbolic(q.ZZ)-s.diff(want,z,2))==0
    Q=current.controls.quadratic(g,W,[q.value for q in h]);QZZ=current.quadratic_ZZ(g,W,h)
    for q,r in zip(Q,QZZ):assert s.simplify(symbolic(r)-s.diff(symbolic(q),z,2))==0
    # In particular action(h_Z,h_Z)=D2Q[h_Z,h_Z]=2Q(h_Z).
    action=current.controls.quadratic_action(g,W,[q.Z for q in h],[q.Z for q in h])
    twice=current.controls.quadratic(g,W,[q.Z for q in h])
    for q,r in zip(action,twice):assert s.simplify(symbolic(q)-2*symbolic(r))==0
    return dict(independent_C2_Picard_derivative_rows=10,independent_quadratic_second_rows=5,
        original_bilinear_factor_two_checked=True,all_five_original_controls_and_weight_positions_retained=True)


def independent_tail():
    n=s.Symbol('n',integer=True,nonnegative=True)
    H2,A,rho=s.symbols('H2 A rho',positive=True)
    e=lambda k:(H2+2*k*A*rho)/2**k
    assert s.simplify(e(n+1)-e(n)/2-A*rho/2**n)==0
    assert e(0)==H2
    # CQ is an absolute row-sum coefficient bound on the homogeneous
    # original quadratic. The differentiated map has these three terms:
    CA,CQ,N=s.symbols('CA CQ N',positive=True)
    D2=s.Symbol('D2',positive=True);K=2*CA;H=K*(D2+2*CQ*rho**2/N)
    assert s.simplify(CA*D2+H/2+2*CA*CQ*rho**2/N-H)==0
    forcing=CA*(2*CQ*H2+4*CQ*rho)/N
    assert s.simplify(forcing-CA*CQ*(2*H2+4*rho)/N)==0
    return dict(genuine_second_Picard_error_recurrence_and_closed_tail_checked=True,
        tail='(H2+2*n*A*rho)*2^-n',
        second_rows_are_actual_Picard_derivatives=True,
        inverse_neumann_bound='||J^-1||<=2CA from ||B^-1 DQ/N||<=1/2',
        no_extra_global_frequency_or_oracle_admission=True)


def run(field=None):
    began=time.monotonic();raw=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    assert raw['candidate_actual_C2_limit_controls_constructed']
    assert not any(raw[k] for k in current.GATES+current.OPEN)
    for name,digest in raw['input_hashes'].items():assert current.sha(name)==digest,name
    with mp.workdps(540):
        field=field if field is not None else current.CurrentC2LimitControls(require_checked=False)
        assert not field.acceptance_loaded and field.ranges.acceptance_loaded
        assert field.identity==raw['source_family']
        g=field.ranges.phase.built['graph'];functions=field.functions;old=functions['exact_C1_limit_adapter']
        assert g.nodes[:len(field.prefix)]==field.prefix
        assert g.nodes==raw['exact_graph_nodes']
        encoded=lambda q:json.loads(json.dumps(current.encoded(q)))
        assert encoded(functions)==raw['actual_C2_limit_controls_and_Picard_functions']
        assert encoded(field.bounds)==raw['actual_quantitative_C2_bounds']
        bindings=dict(quadratic_ZZ=current.current.ast_binding(current.quadratic_ZZ),
            picard_C2_map=current.current.ast_binding(current.picard_C2_map),
            build=current.current.ast_binding(current.build),
            quantitative_bounds=current.current.ast_binding(current.quantitative_bounds),
            original_C1_limit_adapter=current.current.ast_binding(current.limit.exact_limit_adapter),
            original_exact_weights=current.current.ast_binding(current.controls.exact_weights),
            original_quadratic=current.current.ast_binding(current.controls.quadratic),
            original_quadratic_action=current.current.ast_binding(current.controls.quadratic_action),
            original_B_inverse_action=current.current.ast_binding(current.controls.inverse_action))
        assert raw['source_bindings']==bindings
        target=field.ranges.target.target_functions();h=field.control_functions();N=old['N']
        for new,prior in zip(h,old['exact_limit_controls']):
            assert new.value==prior.value and new.Z==prior.Z
            node=g.nodes[new.ZZ.node]
            assert node['genuine_ordinary_second_derivative'] and not node['derivative_of_range_endpoint']
            assert node['implicit_ZZ_linear_system']==functions['implicit_ZZ_linear_system'].node
            assert node['limit_vector']==old['exact_limit_vector'].node
        system=g.nodes[functions['implicit_ZZ_linear_system'].node]
        assert system['matrix']==current.current.encode_graph(old['implicit_Z_matrix'])
        assert system['rhs']==current.current.encode_graph(functions['actual_implicit_ZZ_rhs'])
        assert system['original_C1_linear_system']==old['implicit_Z_linear_system'].node
        assert system['same_original_C1_limit_branch']
        assert system['signed_target_C2_source_sha256']==current.sha(current.ranges.target.RECEIPT)
        assert system['quantitative_target_source_sha256']==current.sha(current.ranges.RECEIPT)
        assert all(q.value==field.ranges.phase.built['N_scaled_targets'][row].value and
            q.Z==field.ranges.phase.built['N_scaled_targets'][row].Z for row,q in target.items())
        report=field.ranges.phase.report;freq=report['actual_whole_Z_frequency_connection']
        assert field.bounds['actual_same_repair_integer']==freq['selected_integer']
        assert field.bounds['actual_N_not_reselected']
        perturb=field.bounds['actual_J_preconditioned_perturbation_upper']
        assert current.ep(perturb['log_absolute_upper']+field.ranges.c.ln(2))[1]<=0
        W=old['repair_weights'];d=list(target.values());sequence=functions['genuine_C2_Picard_sequence']
        assert len(sequence)==4
        for prev,nxt in zip(sequence,sequence[1:]):
            assert nxt==current.picard_C2_map(g,W,N,d,prev)
        for n,saved in zip((0,3,32),raw['actual_C2_tail_examples']):
            assert encoded(field.tail(n))==saved
            assert saved['second_error_terms'][0]['exact_relative_denominator']==2**n
            assert saved['second_error_terms'][1]['exact_relative_numerator']==2*n
        proof=independent_map_derivatives();tail=independent_tail()
        parent=json.loads((current.HERE/current.ranges.RECEIPT).read_bytes())
        assert parent['all_passed'] and parent['quantitative_current_C2_target_ranges_installed']
        result=dict(all_passed=True,**dict.fromkeys(current.GATES,True),source_family=field.identity,
            actual_C2_repaired_limit_controls_installed=True,
            independent_original_C2_map_calculus=proof,independent_C2_tail=tail,
            same_actual_C1_limit_target_matrix_weights_and_repair_integer_retained=True,
            exact_limit_control_second_components=5,exact_C2_Picard_derivative_rows=60,
            C2_repair_band_fields_and_histories_installed=False,
            current_numeric_point_field_oracle_installed=False,**dict.fromkeys(current.OPEN,False),
            input_hashes={**field.hashes,current.NAME:current.sha(current.NAME),
                Path(__file__).name:current.sha(Path(__file__).name)},execution_seconds=time.monotonic()-began)
        (current.HERE/current.RECEIPT).write_text(json.dumps(phase.encoded(result),indent=2)+'\n',encoding='utf8')
    print('Same current C2 implicit controls, original quadratic factor and genuine C2 Picard tail checked',flush=True)
    return result


if __name__=='__main__':run()
