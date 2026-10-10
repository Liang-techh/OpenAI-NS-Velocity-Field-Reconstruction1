"""Independent Rh Duhamel/FTC, seam primitive and source-range admission."""
import hashlib
import json
import math
from pathlib import Path
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_Rh_C3_continuation as current


def symbolic(g,node,bindings,memo=None):
    memo={} if memo is None else memo
    if node in bindings:return bindings[node]
    if node in memo:return memo[node]
    n=g.nodes[node];op=n['operation'];at=lambda i:symbolic(g,i,bindings,memo)
    if op=='exact_rational':value=s.Rational(n['numerator'],n['denominator'])
    elif op=='bound_variable':value=s.Symbol(n['name'])
    elif op=='sum':value=s.Add(*(at(i) for i in n['arguments']))
    elif op=='product':value=s.Mul(*(at(i) for i in n['arguments']))
    elif op=='negative':value=-at(n['argument'])
    elif op=='positive_quotient':value=at(n['numerator'])/at(n['denominator'])
    elif op=='analytic_unary':
        value={'exp':s.exp,'log':s.log}[n['name']](at(n['argument']))
    elif op=='definite_integral':value=s.Integral(at(n['integrand']),(s.Symbol(n['variable']),at(n['lower']),at(n['upper'])))
    elif op=='function_substitution':value=at(n['expression']).subs(at(n['variable']),at(n['value']))
    else:raise ValueError('Unbound Rh function node '+str(n))
    memo[node]=value;return value


def bind(bindings,q,expr,z):
    for j,row in enumerate(current.target.rows(q)):bindings[row.node]=s.diff(expr,z,j)


def partial_FTC(field):
    g=field.graph;f=field.functions;a=s.Symbol('Rh_reference_endpoint');t=s.Symbol('native_Rh_reference');z=s.Symbol('Z')
    bindings={};count=0
    for key,rate in current.current.RATES.items():
        H=s.Function('incoming_'+key)(z);D=s.Function('density_'+key)(t,z);r=s.Rational(str(rate))
        bind(bindings,f['source_incoming_N_scaled_C3'][key],H,z)
        bind(bindings,f['source_density_N_scaled_C3'][key],D,z)
        expected=s.exp(-r*(a+5))*H+s.Integral(s.exp(-r*(a-t))*D,(t,-5,a))
        for j,row in enumerate(current.target.rows(f['partial_N_scaled_C3'][key])):
            got=symbolic(g,row.node,bindings)
            assert s.expand(got-s.diff(expected,z,j))==0,(key,j)
            assert s.simplify(s.diff(got,a)+r*got-s.diff(D,z,j).subs(t,a))==0,(key,j,'FTC')
            assert s.simplify(got.subs(a,-5)-s.diff(H,z,j))==0,(key,j,'inlet')
            count+=1
        for j,row in enumerate(current.target.rows(f['local_partial_N_scaled_C3'][key])):
            n=g.nodes[row.node]
            assert n['operation']=='definite_integral' and n['variable']=='native_Rh_reference'
            assert n['lower']==g.constant(-5).node and n['upper']==f['endpoint'].node
            assert n['ordinary_slow_Z_derivative_order']==j and n['endpoint_and_kernel_Z_independent']
    assert count==20
    return dict(independent_partial_Duhamel_Z0_Z1_Z2_Z3_rows=count,
        own_rate_radial_FTC_and_nonzero_inlet_checked=True,pressure_rate_zero_memory_one=True)


def signed_density_and_complete(field):
    g=field.graph;f=field.functions;z=s.Symbol('Z');N=s.Symbol('N',positive=True);bindings={field.N.node:N}
    E,V,F,B=(s.Function(key)(z) for key in ('E0','V0','F','B'))
    primitive=field.raw['actual_C3_primitive_functions']['Rh_reference'];w=field.window
    jet=lambda rows:current.target.C3Function(*(current.source.FunctionRef(g,i) for i in rows))
    for rows,expr in ((primitive['source_roots']['E'],E),(primitive['source_roots']['V'],V),(w['F_N'],F),(primitive['B'],B)):
        bind(bindings,jet(rows),expr,z)
    changedE,changedV=E+F/N,V+B/N
    original=dict(m=V,h=E,k=E*V,e=V*V-E*E/2,p=E*E/2)
    changed=dict(m=changedV,h=changedE,k=changedE*changedV,e=changedV**2-changedE**2/2,p=changedE**2/2)
    count=0
    for key,q in f['source_density_N_scaled_C3'].items():
        wanted=s.expand(N*(changed[key]-original[key]))
        for j,row in enumerate(current.target.rows(q)):
            assert s.cancel(symbolic(g,row.node,bindings)-s.diff(wanted,z,j))==0,(key,j)
            count+=1
    bindings={field.N.node:N};complete={}
    for key in current.current.RATES:
        lead=s.Function('lead_'+key)(z);correction=s.Function('correction_'+key)(z)
        bind(bindings,f['actual_leading_history_C3'][key],lead,z)
        bind(bindings,f['partial_N_scaled_C3'][key],correction,z)
        complete[key]=lead+correction/N
        for j,row in enumerate(current.target.rows(f['actual_complete_history_C3'][key])):
            assert s.expand(symbolic(g,row.node,bindings)-s.diff(complete[key],z,j))==0
            count+=1
    X,Y=s.Function('corrected_E')(z),s.Function('corrected_V')(z)
    bind(bindings,f['actual_corrected_E_C3'],X,z);bind(bindings,f['actual_corrected_V_C3'],Y,z)
    expected=dict(m=Y-complete['m'],h=X-3*complete['h']/2,k=X*Y-3*complete['k']/2,
        e=Y*Y-X*X/2-complete['e'],p=X*X/2)
    for key,q in f['actual_complete_history_y_C3'].items():
        for j,row in enumerate(current.target.rows(q)):
            assert s.expand(symbolic(g,row.node,bindings)-s.diff(expected[key],z,j))==0
            count+=1
    P0=s.Function('P0')(z);bind(bindings,f['actual_independent_P0_C3'],P0,z)
    for j,row in enumerate(current.target.rows(f['actual_absolute_pressure_C3'])):
        assert s.expand(symbolic(g,row.node,bindings)-s.diff(P0+complete['p'],z,j))==0;count+=1
    assert count==64
    return dict(independent_full_minus_leading_signed_density_and_complete_history_rows=count,
        correction_added_once_and_original_pressure_not_reset=True)


def canonical_signature(g,node,bindings,chart,memo):
    if node in bindings:return bindings[node]
    if node in memo:return memo[node]
    n=g.nodes[node];op=n['operation'];at=lambda i:canonical_signature(g,i,bindings,chart,memo)
    if op=='exact_rational':value=(op,n['numerator'],n['denominator'])
    elif op=='bound_variable':value=(op,'angle' if n['name']=='loop_angle_'+chart else n['name'])
    elif op in ('sum','product','logical_or'):value=(op,tuple(at(i) for i in n['arguments']))
    elif op=='negative':value=(op,at(n['argument']))
    elif op=='positive_quotient':value=(op,at(n['numerator']),at(n['denominator']))
    elif op=='analytic_unary':value=(op,n['name'],at(n['argument']))
    elif op=='mathematical_pi':value=(op,)
    elif op=='current_original_source_parameter':value=(op,n['name'])
    elif op=='exact_real_comparison':value=(op,n['operator'],at(n['left']),at(n['right']))
    elif op=='original_lazy_flat_branch':value=(op,at(n['active_body']),at(n['flat_predicate']),at(n['flat_value']))
    elif op=='original_monotone_phase_inverse':
        value=(op,at(n['phase_function']),at(n['target_phase']),at(n['angle_lower']),at(n['angle_upper']))
    elif op=='substitute_original_inverse_angle':value=(op,at(n['body']),at(n['inverse_angle']))
    elif op=='definite_integral':value=(op,at(n['integrand']),at(n['lower']),at(n['upper']),
        'angle' if n['variable']=='loop_angle_'+chart else n['variable'])
    else:raise ValueError('Unbound source seam dependency '+str(n))
    value=hashlib.sha256(json.dumps(value,separators=(',',':')).encode()).hexdigest()
    memo[node]=value;return value


def seam_geometry_and_primitives(field):
    g=field.graph;raw=field.connection['exact_current_integral_and_control_graph'];windows={w['chart']:w for w in raw['windows']}
    offsets,domains,x=current.current.exact_geometry();r=s.Symbol('Rm')
    assert s.simplify(offsets['actual_patch'].subs(x,s.E)-offsets['Rh_reference'].subs(x,-5))==0
    assert s.diff(offsets['Rh_reference'],x)==1 and domains['Rh_reference']==(-5,0)
    assert field.join_report['exact_current_function_proof']==current.join.exact_join_proof()
    assert field.join_report['exact_current_function_proof']['implied_ordinary_y_Z_mixed4_rows']==135
    signatures={};counts=0
    for chart in ('actual_patch','Rh_reference'):
        primitive=field.raw['actual_C3_primitive_functions'][chart];w=windows[chart]
        bindings={w['phase']:'same_original_seam_phase',field.N.node:'same_selected_N'}
        for name,rows in primitive['source_roots'].items():
            for j,node in enumerate(rows):bindings[node]='same_leading_'+name+'_Z'+str(j)
        memo={};signature={}
        for key in ('A','B'):
            signature[key]=[canonical_signature(g,node,bindings,chart,memo) for node in primitive[key]]
            counts+=4
        signatures[chart]=signature
    assert signatures['actual_patch']==signatures['Rh_reference']
    incoming=field.correction_functions(-5);outgoing=field.correction_functions(0)
    assert current.target.encoded(incoming)==field.window['incoming']==field.predecessor['outgoing']
    assert current.target.encoded(outgoing)==field.window['outgoing']
    assert any(q.value!=g.zero for q in incoming.values())
    for bad in (-6,1):
        try:field.correction_functions(bad)
        except ValueError:pass
        else:raise AssertionError('Out-of-domain Rh continuation admitted')
    parts=field.window['original_physical_measure_and_partitions']
    assert parts[0]['exact_left']==[-5,1] and parts[-1]['exact_right']==[0,1]
    assert all(a['exact_right']==b['exact_left'] for a,b in zip(parts,parts[1:]))
    a=s.Symbol('Rh_reference_endpoint');t=s.Symbol('native_Rh_reference');z=s.Symbol('Z');partition_rows=0
    for key,rate in current.current.RATES.items():
        H=s.Function('inlet_'+key)(z);D=s.Function('source_'+key)(t,z);rate=s.Rational(str(rate))
        bindings={field.geometry['offset']:t,field.geometry['right_offset']:s.Integer(0),
            field.geometry['left_offset']:s.Integer(-5),field.geometry['width']:s.Integer(5)}
        bind(bindings,field.functions['source_incoming_N_scaled_C3'][key],H,z)
        bind(bindings,field.functions['source_density_N_scaled_C3'][key],D,z)
        for j,node in enumerate(field.window['outgoing'][key]):
            expected=s.exp(-5*rate)*s.diff(H,z,j)+sum((s.Integral(s.exp(rate*t)*s.diff(D,z,j),
                (t,s.Rational(*cell['exact_left']),s.Rational(*cell['exact_right']))) for cell in parts),s.Integer(0))
            assert s.expand(symbolic(g,node,bindings)-expected)==0,(key,j,'partitioned_outgoing')
            partition_rows+=1
    assert partition_rows==20
    # At a=0, finite-integral additivity yields the accepted partitioned outgoing;
    # at a=-5 every local integral vanishes. Neither identity discards the inlet.
    return dict(independent_radius_phase_and_Jacobian_identity=True,
        source_leading_mixed4_identity_and_common_generic_roots_retained=True,
        exact_chart_renamed_C3_primitive_signature_rows=counts,
        exact_nonzero_incoming_and_full_outgoing_alias_rows=40,
        independently_checked_partitioned_outgoing_function_rows=partition_rows,
        whole_reference_integral_partition_additivity_checked=True,domain_guards_checked=True)


def source_projection_and_ranges(field,raw):
    counts=0
    for key,q in field.functions['actual_leading_history_C3'].items():
        for j,ref in enumerate(current.target.rows(q)):
            n=field.graph.nodes[ref.node]
            assert n['quantity']=='history_'+key and n['Z_order']==j and n['Taylor_coefficient_factorial']==math.factorial(j)
            assert n['coordinate']==field.functions['endpoint'].node
            assert n['phase_independent_leading_source']
            assert n['source_projection_binding']==current.current.ast_binding(current.CurrentRhC3Continuation.project_leading_packet)
            counts+=1
    for q in current.target.rows(field.functions['actual_independent_P0_C3']):
        node=field.graph.nodes[q.node]
        assert node['coordinate_independent'] and node['phase_independent_leading_source']
        assert node['exact_source_op_P0_function_not_enclosure_endpoint']
    assert len(field.P0_guards)==4 and field.geometry_binding['source_roots_C3_prefix_equals_geometry_roots']
    c=field.c;rows=[c.mpf(j+1) for j in range(5)];hist={key:list(rows) for key in current.current.RATES}
    packet=dict(exact_common_P0_axial5=rows,original_generic_source=dict(common_original_P0_axial5=rows,common_own_five_histories_axial5=hist))
    projected=field.project_leading_packet(packet)
    for q in projected['histories'].values():
        assert all(current.ep(q[j]-math.factorial(j)*rows[j])==(0,0) for j in range(4))
    for j in range(4):assert current.ep(projected['P0'][j]-math.factorial(j)*rows[j])==(0,0)
    encoded=lambda q:json.loads(json.dumps(current.target.encoded(current.ranges.record(q))))
    cells=0
    for saved in raw['actual_four_Z_Rh_C3_partial_ranges']:
        ends=tuple(saved['exact_Z_cell']);assert encoded(field.quantitative_range(ends))==saved
        prior=field.rows[ends]
        for key in current.current.RATES:
            for bound,old in zip(saved['actual_N_scaled_partial_history_C3_bounds'][key],prior['actual_N_scaled_incoming_C3'][key]):
                if old['exact_zero']:continue
                newlog=current.current.packets.interval(c,bound['log_absolute_upper']);oldlog=current.current.packets.interval(c,old['log_absolute_upper'])
                assert current.ep(newlog)[1]>=current.ep(oldlog)[1],(ends,key)
        assert saved['range_covers_all_partial_endpoints'] and saved['original_nonzero_memory_and_independent_P0_retained']
        cells+=len(saved['original_native_cells'])
    assert counts==20 and cells==20
    return dict(actual_leading_history_source_recipe_rows=counts,ordinary_Taylor_factorials_applied_once=True,
        actual_four_axial_cells=4,actual_reference_native_cells=cells,
        whole_partial_history_majorants_retain_accepted_incoming=True,
        exact_original_reference_masses_and_suffix_decay_at_most_one=True,
        function_definitions_are_not_selected_range_values=True)


def run(field=None):
    began=time.monotonic();raw=current.read(current.NAME)
    for name,digest in raw['input_hashes'].items():assert current.sha(name)==digest,name
    assert raw['candidate_actual_Rh_C3_continuation_constructed'] and not any(raw[k] for k in current.GATES+current.OPEN)
    with mp.workdps(540):
        field=field if field is not None else current.CurrentRhC3Continuation(require_checked=False)
        assert not field.acceptance_loaded and field.graph.nodes==raw['exact_graph_nodes']
        assert field.graph.nodes[:len(field.prefix)]==field.prefix
        encoded=lambda q:json.loads(json.dumps(current.target.encoded(q)))
        assert encoded(field.functions)==raw['actual_Rh_C3_continuation_functions']
        assert raw['source_bindings']==current.source_bindings() and field.identity==raw['source_family']
        partial=partial_FTC(field);signed=signed_density_and_complete(field);seam=seam_geometry_and_primitives(field)
        source_ranges=source_projection_and_ranges(field,raw)
        assert not raw['current_numeric_point_field_oracle_installed']
        result=dict(all_passed=True,**dict.fromkeys(current.GATES,True),source_family=field.identity,
            independent_partial_FTC=partial,independent_signed_density_and_complete_history=signed,
            independent_Rh_seam=seam,actual_source_projection_and_ranges=source_ranges,
            current_numeric_point_field_oracle_installed=False,global_physical_time_Cartesian_heat_cone_and_temporal_recursion_installed=False,
            **dict.fromkeys(current.OPEN,False),input_hashes={**field.hashes,current.NAME:current.sha(current.NAME),
                Path(__file__).name:current.sha(Path(__file__).name)},execution_seconds=time.monotonic()-began)
        (current.HERE/current.RECEIPT).write_text(json.dumps(current.target.encoded(result),indent=2)+'\n',encoding='utf8')
    print('Actual Rh C3 signed partial/complete histories, nonzero seam memory and reference bounds checks passed',flush=True)
    return result


if __name__=='__main__':run()
