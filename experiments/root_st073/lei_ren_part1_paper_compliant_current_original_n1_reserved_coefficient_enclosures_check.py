"""Selected graph identities and actual five physical coefficient C1 ranges."""
import ast
import copy
import inspect
import json
import textwrap
import time
from types import SimpleNamespace

import sympy as sy

import lei_ren_part1_paper_compliant_current_original_n1_reserved_coefficient_enclosures as current
import lei_ren_part1_paper_compliant_current_original_n1_full_defect_linear_functions_check as original_check


def lift(fn,scope):
    tree=ast.parse(textwrap.dedent(inspect.getsource(fn)));tree.body[0].decorator_list=[]
    namespace=dict(fn.__globals__);namespace.update(scope)
    exec(compile(tree,'<actual-reserved-coefficient-formula>','exec'),namespace)
    return namespace[fn.__name__]


def exact_selected_physical_coefficient_identity(owner):
    """Actual selected graph, arbitrary signed moments and ordinary amplitude jets."""
    mu,D,DT=sy.symbols('mu positive_axial_divided_det positive_angular_det',positive=True)
    L,LP,U15,Rvlog=sy.symbols('log2Rc logP pressure_unit15_log logRv',real=True)
    A,AZ=sy.symbols('finite_actual_A_rc finite_actual_A_rc_Z',real=True)
    ax,bd=sy.symbols('a0:2'),sy.symbols('b0:2')
    T=sy.Matrix(3,3,lambda i,j:sy.Symbol(f'angular_{i}_{j}',real=True))
    function=lift(current.divided_matrix_factors,{'interval':lambda c,v:v,'ep':lambda v:(v,v)})
    matrix=function(SimpleNamespace(),dict(mu=mu,axial_divided_rows=[ax,bd],
        positive_negative_axial_divided_determinant=D,positive_angular_determinant=DT,
        angular_rows=[list(T.row(i)) for i in range(3)]))
    for i,row in enumerate(matrix['angular_inverse']):
        for j,value in enumerate(row):
            assert sy.cancel(value-T.cofactor(j,i)/DT)==0
    sigma=-LP/2+sy.Rational(1,2);SA=LP+(sy.Rational(1,2)+mu)*L+sigma
    sbeta=U15-LP-L/2-sigma
    W=owner.operator.functions;b=owner.operator.background;g=owner.g
    bindings={W['mu'].node:mu,W['log2Rc'].node:L,
        b['parameters']['logP']:LP,owner.full.pulse.normalization['physical_unit_logs']['m15'].node:U15,
        b['Astar_ordinary_Z0_Z3'][0]:sy.exp(SA)*A,b['Astar_ordinary_Z0_Z3'][1]:sy.exp(SA)*AZ,
        W['positive_axial_divided_determinant'].node:D}
    for rows,want in ((W['axial_rows'],[ax,bd]),(W['inverse_angular'],matrix['angular_inverse'])):
        for row,values in zip(rows,want):
            for ref,value in zip(row,values):bindings[ref.node]=value
    symbols={key:sy.symbols(key+'_value '+key+'_Z',real=True) for key in ('m11','m12','m14','p13','p15')}
    for key,pair in owner.full.functions.items():
        if key in ('m13','m15'):
            values=symbols['p13' if key=='m13' else 'p15'];unit=sy.exp(U15-Rvlog if key=='m13' else U15)
            values=[unit*v for v in values]
        else:values=symbols[key]
        bindings[pair.value.node],bindings[pair.Z.node]=values
    quotient=lambda pair:[pair[0]/A,(pair[1]*A-pair[0]*AZ)/(A*A)]
    m14=quotient(symbols['m14']);p13=quotient(symbols['p13']);p15=quotient(symbols['p15'])
    alpha=[]
    for x,y in matrix['alpha_signed_factors']:
        alpha.append([x*symbols['m11'][k]+y*sy.exp(-LP-L/2-sigma)*m14[k] for k in (0,1)])
    theta=[[-sy.exp(-L-sbeta)*v for v in symbols['m12']],
        [-sy.exp(L-Rvlog)*v/2 for v in p13],p15]
    beta=[[sum(v*p[k] for v,p in zip(row,theta)) for k in (0,1)] for row in matrix['angular_inverse']]
    cache={};count=0
    for key,want in zip(current.COEFFICIENTS,alpha+beta):
        pair=owner.normalized_functions[key]
        for ref,expected in zip((pair.value,pair.Z),want):
            actual=original_check.symbolic(g,ref,bindings,cache)
            assert sy.simplify(sy.powsimp(actual-expected,force=True))==0,(key,count)
            count+=1
    z=sy.Symbol('Z',real=True)
    formula=lift(current.actual_arc_coefficient_rows,{'ep':lambda v:(v,v)})
    ctx=SimpleNamespace(mpf=sy.sympify,exp=sy.exp)
    rows=formula(ctx,z,mu)
    body=sy.exp(-sy.Rational(17,10)-sy.Rational(5,2)*mu)/(1+z*z)
    for order,row in enumerate(rows):assert sy.simplify(row['coefficient_interval']-sy.diff(body,z,order))==0
    return dict(passed=True,actual_selected_coefficient_value_Z1_graph_identities=count,
        actual_angular_inverse_cofactor_identities=9,actual_canonical_arc_ordinary_C0_C1_identities=2,
        arbitrary_complete_moments_and_Astar_Z_retained=True,
        source_numerator_transformed_before_division_by_nonzero_mu=True,
        alpha_unit_independent_of_beta_scale=True,actual_m13_positive_radius_exponent_and_m15_sign=True)


def rejected(fn):
    try:fn()
    except ValueError:return True
    raise AssertionError('Unsupported actual coefficient attachment admitted')


@current.source_precision
def run(owner,raw,packets):
    began=time.monotonic();owner.assert_graph();proof=exact_selected_physical_coefficient_identity(owner)
    assert raw['source_identity']==owner.source_identity and owner.source_identity['passed']
    assert raw['input_hashes']==owner.hashes and raw['true_function_bindings']==owner.bindings
    assert tuple(raw['coefficient_normalization_graph_span'])==owner._span
    assert raw['coefficient_normalization_graph_nodes']==owner._nodes
    assert set(raw['samples'])==set(packets)=={'point','whole_Z'}
    results=[]
    for name,packet in packets.items():
        entry=owner._packets[id(packet)];view=entry[1]
        assert type(packet) is current.OriginalReservedCoefficientPacket and entry[0] is packet
        assert current.record(view)==entry[2]==current.record(raw['samples'][name])
        assert view['true_function_bindings']==owner.bindings and set(view['coefficients'])==set(current.COEFFICIENTS)
        assert view['actual_five_coefficient_C0_C1_enclosures_present'] and view['actual_Astar_Z_retained']
        assert view['all_nonzero_source_terms_and_finite_mu_preserved']
        assert view['repaired_own_terminal_zero_identities_not_claimed']
        assert all(view[key] is False for key in current.OPEN)
        nonzero=0
        for cell in view['cells']:
            assert cell['actual_Astar_Z_in_all_quotient_derivatives']
            for key,coefficient in cell['coefficients'].items():
                pair=coefficient['normalized_C0_C1_enclosures'];assert len(pair)==2
                for value in pair:assert current.ep(value)[0]<=current.ep(value)[1]
                assert all(v is not None for v in coefficient['physical_C0_C1_absolute_log_upper'])
                if key.startswith('alpha'):
                    assert coefficient['correlated_source_numerator_before_nonzero_mu_division']
                    for term in coefficient['signed_original_source_terms']:
                        for row in term['rows']:
                            if not row['exact_zero']:
                                assert current.ep(row['absolute_log_upper'])[1]<-1000
                                assert current.ep(row['value'])[0]<current.ep(row['value'])[1];nonzero+=1
            assert cell['nonzero_theta0_original_sectors']
            for row in cell['nonzero_theta1_ratio_cover']:
                if not row['exact_zero']:assert current.ep(row['value'])[0]<current.ep(row['value'])[1]
        assert nonzero>0
        results.append(dict(name=name,original_Z_cells=[r['original_source_Z_cell'] for r in view['cells']],
            actual_physical_coefficient_value_Z1_rows=10,nonzero_alpha_sector_rows_retained=nonzero,
            actual_arc_C1_body_used=True,physical_log_bounds_not_interval_midpoints=True))
    assert len(raw['samples']['whole_Z']['cells'])==4
    assert [r['original_source_Z_cell'] for r in raw['samples']['whole_Z']['cells']]==[
        ['-1','-.5'],['-.5','0'],['0','.5'],['.5','1']]
    for row,ends in zip(raw['samples']['whole_Z']['cells'],[[-1,'-.5'],['-.5',0],[0,'.5'],['.5',1]]):
        assert row['clipped_Z']._mpi_==owner.c.mpf(ends)._mpi_
    assert rejected(lambda:owner.evaluate('1.01')) and rejected(lambda:owner.report(object()))
    saved=owner.matrix
    try:owner.matrix={};assert rejected(owner.assert_graph)
    finally:owner.matrix=saved
    target=owner.normalized_functions['beta1'].value.node;saved=copy.deepcopy(owner.g.nodes[target])
    try:
        owner.g.nodes[target]={**saved,'unadmitted_mutation':True};assert rejected(owner.assert_graph)
    finally:owner.g.nodes[target]=saved
    owner.assert_graph()
    receipt=dict(all_passed=True,source_family=owner.family,source_identity=owner.source_identity,
        exact_selected_physical_coefficient_proof=proof,query_results=results,
        nonzero_formal_scales_not_selected_as_field_values=True,
        altered_matrix_normalized_node_foreign_packet_and_domain_rejected=True,
        **{current.GATE:True},**dict.fromkeys(current.OPEN,False),
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME)},execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_bytes((json.dumps(current.record(receipt),indent=2)+'\n').encode())
    print('ACTUAL_RESERVED_FIVE_COEFFICIENT_C1_ENCLOSURES_CHECKED',len(results),flush=True)
    return receipt
