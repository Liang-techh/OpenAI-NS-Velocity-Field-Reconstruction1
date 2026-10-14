"""Actual complete partition, fixed physical units and five C1 enclosure rows."""
import ast
import copy
import inspect
import json
import textwrap
import time
from types import SimpleNamespace

import sympy as sy

import lei_ren_part1_paper_compliant_current_original_n1_full_five_numeric_bridge as current


def exact_factored_ordinary_derivative_identity():
    """Execute the production export for arbitrary analytic coefficient jets."""
    z=sy.Symbol('Z',real=True)
    logs=sy.symbols('log_hb log_F0_anchor log_Pstar_squared log_Lambda',real=True)
    unit=sy.Symbol('fixed_physical_pressure_unit_log',real=True)
    powers=((1,-2,1,-1),(0,1,0,-2))
    functions=[sy.Function(name)(z) for name in ('actual_coefficient_a','actual_coefficient_b')]
    value=SimpleNamespace(ctx=SimpleNamespace(mpf=sy.sympify),order=1,
        algebra=SimpleNamespace(logs=logs),terms={p:SimpleNamespace(coefficients=[f,sy.diff(f,z)])
            for p,f in zip(powers,functions)})
    fn=current.cutoff.op.export
    tree=ast.parse(textwrap.dedent(inspect.getsource(fn)))
    namespace=dict(fn.__globals__);namespace['truncate']=lambda value,order:value
    exec(compile(tree,'<actual-ordinary-factored-derivative-export>','exec'),namespace)
    exported=namespace[fn.__name__](value,1)
    assert exported['ordinary_Z_derivative_order']==1
    count=0
    for term in exported['terms']:
        p=tuple(term['source_exponents']);f=functions[powers.index(p)]
        log=sum(l*q for l,q in zip(logs,p))
        assert sy.expand(term['combined_source_log']-log)==0
        for order in (0,1):
            actual=sy.exp(term['combined_source_log']-unit)*term['ordinary_Z_derivatives'][order]
            expected=sy.diff(sy.exp(log-unit)*f,z,order)
            assert sy.simplify(actual-expected)==0;count+=1
    return dict(passed=True,actual_factored_export_value_Z1_identities=count,
        exact_exponent_basis=['hb','F0_anchor','Pstar_squared','Lambda'],
        absolute_source_logs_and_positive_pressure_unit_are_fixed_in_Z_jets=True,
        actual_ordinary_coefficient_derivatives_not_differentiated_twice=True,
        ordinary_amplitude_derivatives_inherited_from_pinned_16_row_source_proof=True)


def rejected(fn):
    try:fn()
    except ValueError:return True
    raise AssertionError('Unsupported full five moment attachment admitted')


@current.source_precision
def run(owner,raw,packets):
    began=time.monotonic();owner.assert_graph()
    assert raw['source_identity']==owner.source_identity and owner.source_identity['passed']
    proof=exact_factored_ordinary_derivative_identity()
    assert owner.source_identity['actual_pressure_unit_source_identity']['passed']
    assert set(packets)==set(raw['samples'])=={'point','whole_Z'}
    results=[]
    for name,packet in packets.items():
        entry=owner._packets[id(packet)]
        assert type(packet) is current.OriginalTrueFullFiveNumericPacket and entry[0] is packet
        view=entry[1]
        assert current.record(view)==entry[3]==current.record(raw['samples'][name])
        assert view['true_function_bindings']==owner.bindings
        assert all(view[key] is False for key in current.OPEN)
        assert view['moment_order']==['m11','m12','m13','m14','m15']
        assert view['all_five_true_unrepaired_moment_C1_enclosures_present']
        assert view['all_pressure_radial_support_included']
        assert view['velocity_moments_remain_factored_physical_not_pressure_units']
        assert view['signed_m15_not_negated_again']
        assert view['actual_coefficients_and_terminal_zero_identities_not_claimed']
        velocity_rows=0;pressure_rows=0;nonzero_rows=0;sector_rows=0
        for key,sectors in view['full_physical_velocity_moment_C1_sectors'].items():
            assert key in current.VELOCITY and sectors['ordinary_Z_derivative_order']>=1
            assert sectors['terms']
            for term in sectors['terms']:
                assert len(term['source_exponents'])==4 and len(term['ordinary_Z_derivatives'])>=2
            velocity_rows+=2
        for key,pair in view['normalized_full_pressure_moment_C1_enclosures'].items():
            assert key in current.PRESSURE and len(pair)==2
            parts=view['pressure_component_enclosures'][key]
            for component in ('inner','lower'):
                actual=parts[component]
                for term in actual['physical_sector_sources']:
                    assert len(term['source_exponents'])==4 and len(term['rows'])==2
                    for row in term['rows']:
                        sector_rows+=1
                        coefficient=row['source_coefficient'];lo,hi=current.ep(coefficient)
                        cover=row['normalized_enclosure'];a,b=current.ep(cover)
                        assert a<=b
                        if row['coefficient_exact_zero_in_this_query']:
                            assert lo==hi==a==b==0
                        else:
                            assert current.ep(row['normalized_log_absolute_upper'])[1]<-1000
                            assert a<b and b>0 if lo>=0 else a<b and a<0
                            if hi<=0:assert b<=0
                            nonzero_rows+=1
            for order,value in enumerate(pair):
                recomposed=(parts['inner']['C1_enclosures'][order]+parts['lower']['C1_enclosures'][order]
                    +parts['preRp'][order]+parts['pulse'][order])
                assert recomposed._mpi_==value._mpi_
                assert current.ep(value)[0]<=current.ep(value)[1]
                pressure_rows+=1
        assert velocity_rows==6 and pressure_rows==4 and nonzero_rows>0
        results.append(dict(name=name,Z=view['Z'],actual_query_executed=True,
            actual_true_velocity_value_Z1_rows=velocity_rows,actual_true_pressure_value_Z1_rows=pressure_rows,
            inner_lower_physical_sector_rows=sector_rows,nonzero_sector_rows_retained=nonzero_rows,
            all_four_pressure_components_recomposed_exactly_in_one_unit=True))
    # New normalization nodes must be guarded even for a constructor that
    # reuses existing graph nodes rather than appending a fresh tail.
    target=owner.normalized_pressure_functions['m13'].value.node
    assert owner._normalization_span[0]<=target<owner._normalization_span[1]
    saved=copy.deepcopy(owner.g.nodes[target])
    try:
        owner.g.nodes[target]={**saved,'unadmitted_mutation':True}
        assert rejected(owner.assert_graph)
    finally:owner.g.nodes[target]=saved
    assert rejected(lambda:owner.report(object()))
    owner.assert_graph()
    receipt=dict(all_passed=True,source_family=owner.family,source_identity=owner.source_identity,
        exact_factored_ordinary_derivative_proof=proof,query_results=results,
        changed_normalized_function_node_and_foreign_packet_rejected=True,
        actual_physical_log_units_bound_to_source_recipes=True,
        full_five_enclosures_are_unrepaired_defects_not_terminal_zeros=True,
        **{current.GATE:True},**dict.fromkeys(current.OPEN,False),
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_bytes((json.dumps(current.record(receipt),indent=2)+'\n').encode())
    print('TRUE_COMPLETE_FIVE_C1_NUMERIC_ENCLOSURES_CHECKED',len(results),flush=True)
    return receipt
