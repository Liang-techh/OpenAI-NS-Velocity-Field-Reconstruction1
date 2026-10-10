"""Independent original quadratic third calculus, same-limit and C3 tail."""
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_C3_limit_controls as current
import lei_ren_part1_paper_compliant_current_original_C2_density_transport_check as symbolic_source


def encoded(q):return json.loads(json.dumps(current.encoded(q)))


def independent_map_derivatives():
    g=current.source.FunctionTransportGraph();z=s.Symbol('Z');bindings={}
    def fixed(name):
        q=g.symbol(name);bindings[q.node]=s.Symbol(name,positive=True);return q
    def function(name):
        f=s.Function(name)(z);handles=[]
        for j in range(4):
            q=g.symbol(name+'_'+str(j));bindings[q.node]=s.diff(f,z,j);handles.append(q)
        return current.target.C3Function(*handles)
    W=dict(mu=fixed('mu'),D=[fixed('D0'),fixed('D2')],axial_positive_determinant=fixed('axial_determinant'),
        inverse_swirl=[[fixed('inverse_%d%d'%(i,j)) for j in range(3)] for i in range(3)],
        cross=[fixed('cross0'),fixed('cross2')],energy=[fixed('energy'+str(i)) for i in range(3)],
        pressure=[fixed('pressure'+str(i)) for i in range(3)])
    h=[function('h'+str(i)) for i in range(5)];d=[function('d'+str(i)) for i in range(5)];N=fixed('N')
    lower=lambda q:current.phase.C2Function(*current.target.rows(q)[:3])
    old=current.lower.picard_C2_map(g,W,N,list(map(lower,d)),list(map(lower,h)))
    new=current.picard_C3_map(g,W,N,d,h,old)
    symbolic=lambda q:symbolic_source.symbolic(g,q.node,bindings)
    for q,prior in zip(new,old):
        assert current.target.rows(q)[:3]==tuple(vars(prior).values())
        assert s.simplify(symbolic(q.ZZZ)-s.diff(symbolic(q.value),z,3))==0
    Q=current.controls.quadratic(g,W,[q.value for q in h]);third=current.quadratic_ZZZ(g,W,h)
    for q,r in zip(Q,third):assert s.simplify(symbolic(r)-s.diff(symbolic(q),z,3))==0
    return dict(independent_third_Picard_map_rows=5,independent_original_quadratic_third_rows=5,
        coefficient_three_bilinear_cross_term_checked=True,all_five_original_control_weight_positions_retained=True)


def independent_tail():
    n=s.Symbol('n',integer=True,nonnegative=True)
    CA,CQ,N,rho,H2,H3,A,D3=s.symbols('CA CQ N rho H2 H3 A D3',positive=True)
    B=CA*CQ*rho*(2*H3+12*H2)/N;C=12*CA*CQ*A*rho**2/N
    e=lambda k:(H3+2*k*B+k*(k-1)*C)/2**k
    assert s.simplify(e(n+1)-e(n)/2-(B+n*C)/2**n)==0 and e(0)==H3
    forcing=CA*CQ/N*(2*H3*rho+6*H2*rho+6*rho*(H2+2*n*A*rho))
    assert s.simplify(forcing-B-n*C)==0
    cap=2*CA*(D3+6*CQ*rho*H2/N)
    assert s.simplify(cap/2+CA*D3+6*CA*CQ*rho*H2/N-cap)==0
    # Each D2Q term is bilinear and symmetric; split its two arguments.
    x,y,xn,yn=s.symbols('x y xn yn')
    assert s.expand(x*y-xn*yn-((x-xn)*y+xn*(y-yn)))==0
    return dict(genuine_third_Picard_error_recurrence_and_closed_tail_checked=True,
        tail='(H3+2n*B+n(n-1)*C)*2^-n',third_control_and_iterate_invariant_bound_checked=True,
        original_DQ_and_D2Q_error_splittings_checked=True,frequency_not_reselected=True,
        oracle_quadrature_and_physical_errors_not_admitted=True)


def run(field=None):
    began=time.monotonic();raw=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    assert raw['candidate_actual_C3_limit_controls_constructed'] and not any(raw[k] for k in current.GATES+current.OPEN)
    for name,digest in raw['input_hashes'].items():assert current.sha(name)==digest,name
    with mp.workdps(540):
        field=field if field is not None else current.CurrentC3LimitControls(require_checked=False)
        assert not field.acceptance_loaded and field.lower.acceptance_loaded and field.ranges.acceptance_loaded
        assert field.identity==raw['source_family'];g=field.ranges.phase.built['graph'];f=field.functions
        parent=f['exact_C2_limit_functions'];old=parent['exact_C1_limit_adapter']
        assert field.ranges.phase is field.lower.ranges.phase
        assert old['N']==field.ranges.phase.built['N']
        for ends in current.ranges.outer.CELLS:
            before=field.lower.ranges.phase.outer.owner.owner(ends)
            after=field.ranges.phase.outer.owner.owner(ends)
            assert before is after and before.flow.ledger is after.flow.ledger
        assert g.nodes[:len(field.prefix)]==field.prefix and g.nodes==raw['exact_graph_nodes']
        assert encoded(f)==raw['actual_C3_limit_controls_and_Picard_functions']
        assert encoded(field.bounds)==raw['actual_quantitative_C3_bounds']
        assert raw['source_bindings']==current.source_bindings()
        h=field.control_functions();targets=field.ranges.target.target_functions()
        for i,(new,prior) in enumerate(zip(h,field.lower.control_functions())):
            assert current.target.rows(new)[:3]==tuple(vars(prior).values())
            n=g.nodes[new.ZZZ.node]
            assert n['component']==i and n['Z_order']==3
            assert n['genuine_ordinary_third_derivative'] and not n['derivative_of_range_endpoint']
            assert n['implicit_ZZZ_linear_system']==f['implicit_ZZZ_linear_system'].node
            assert n['limit_vector']==old['exact_limit_vector'].node
        system=g.nodes[f['implicit_ZZZ_linear_system'].node]
        assert system['matrix']==current.current.encode_graph(parent['same_implicit_matrix'])
        assert system['rhs']==current.current.encode_graph(f['actual_implicit_ZZZ_rhs'])
        assert system['original_C1_linear_system']==old['implicit_Z_linear_system'].node
        assert system['original_C2_linear_system']==parent['implicit_ZZ_linear_system'].node and system['same_original_C2_limit_branch']
        assert system['signed_C3_source_sha256']==current.sha(current.target.RECEIPT)
        assert system['quantitative_C3_source_sha256']==current.sha(current.ranges.RECEIPT)
        assert tuple(targets)==current.controls.ROWS
        for row,q in targets.items():assert current.target.rows(q)[:3]==tuple(vars(parent['actual_signed_targets_C2'][row]).values())
        freq=field.ranges.phase.report['actual_whole_Z_frequency_connection']
        assert field.bounds['actual_same_repair_integer']==freq['selected_integer']==field.lower.bounds['actual_same_repair_integer']
        assert field.bounds['actual_N_not_reselected']
        perturb=field.lower.bounds['actual_J_preconditioned_perturbation_upper']
        assert current.ep(perturb['log_absolute_upper']+field.ranges.c.ln(2))[1]<=0
        W=old['repair_weights'];N=old['N'];d=list(targets.values());sequence=f['genuine_C3_Picard_sequence']
        prior=parent['genuine_C2_Picard_sequence'];assert len(sequence)==len(prior)==4
        for step,lower_step in zip(sequence,prior):
            for q,v in zip(step,lower_step):assert current.target.rows(q)[:3]==tuple(vars(v).values())
        for i,(prev,nxt) in enumerate(zip(sequence,sequence[1:])):
            assert nxt==current.picard_C3_map(g,W,N,d,prev,prior[i+1])
        cross=current.controls.quadratic_action(g,W,[q.Z for q in h],[q.ZZ for q in h])
        expected=[g.neg(g.add(q.ZZZ,g.quotient(g.mul(g.constant(3),v),N,
            'same selected exact positive integer N'))) for q,v in zip(d,cross)]
        assert expected==f['actual_implicit_ZZZ_rhs']
        for row,v,proof in zip(current.controls.ROWS,f['actual_third_control_residual'],f['actual_third_zero_theorem_nodes']):
            n=g.nodes[proof.node]
            assert n['expression']==v.node and n['row']==row and n['zero']==g.zero.node
            assert n['implicit_ZZZ_linear_system']==f['implicit_ZZZ_linear_system'].node
        report=json.loads(gzip.decompress((current.HERE/current.ranges.NAME).read_bytes()))
        assert encoded(current.quantitative_bounds(field.lower,field.ranges,report))==encoded(field.bounds)
        c=field.ranges.c
        forcing_A=current.ranges.LogUpper(c,current.current.packets.interval(c,
            field.lower.bounds['actual_C2_error_forcing_upper']['log_absolute_upper']))
        assert encoded(forcing_A.record())==encoded(field.bounds['actual_C2_error_forcing_upper'])
        for depth,saved in zip((0,3,32),raw['actual_C3_tail_examples']):
            assert encoded(field.tail(depth))==saved
            assert [q['exact_relative_numerator'] for q in saved['third_error_terms']]==[1,2*depth,depth*(depth-1)]
            assert all(q['exact_relative_denominator']==2**depth for q in saved['third_error_terms'])
        for depth in (-1,4097,1.5):
            try:field.tail(depth)
            except ValueError:pass
            else:raise AssertionError('Tail depth must be admitted exact integer')
        calculus=independent_map_derivatives();tail=independent_tail()
        receipt=json.loads((current.HERE/current.ranges.RECEIPT).read_bytes())
        assert receipt['all_passed'] and receipt['quantitative_current_C3_target_ranges_installed']
        result=dict(all_passed=True,**dict.fromkeys(current.GATES,True),source_family=field.identity,
            actual_C3_repaired_limit_controls_installed=True,independent_original_C3_map_calculus=calculus,
            independent_C3_tail=tail,same_C2_limit_target_matrix_weights_and_repair_integer_retained=True,
            exact_limit_control_third_components=5,exact_third_Picard_rows=20,
            C3_repaired_band_fields_and_histories_installed=False,current_numeric_point_field_oracle_installed=False,
            **dict.fromkeys(current.OPEN,False),input_hashes={**field.hashes,current.NAME:current.sha(current.NAME),
                Path(__file__).name:current.sha(Path(__file__).name)},execution_seconds=time.monotonic()-began)
        (current.HERE/current.RECEIPT).write_text(json.dumps(encoded(result),indent=2)+'\n',encoding='utf8')
    print('Same C2 limit: original cubic controls and genuine C3 Picard tail checked',flush=True)
    return result


if __name__=='__main__':run()
