"""Relational defining-DAG proof of the actual correction inlet interface.

Fingerprints only align commutative children. Every defining operation and
source parameter is then compared directly; hash equality is not the proof.
"""
import json
from pathlib import Path
import time
import sympy as s
import lei_ren_part1_paper_compliant_current_original_outer_to_repair_C3_bridge as current
import lei_ren_part1_paper_compliant_current_original_Rh_C3_continuation_check as rh_check


class RelationalFunctions:
    def __init__(self,a,b):
        self.a,self.b=a,b;self.memo=set();self.operations={}
        self.align_a=current.ExactFunctionSignatures(a);self.align_b=current.ExactFunctionSignatures(b)

    def equal(self,i,j):
        if (i,j) in self.memo:return
        a,b=self.a[i],self.b[j];op=a['operation']
        assert b['operation']==op,(i,j,'operation')
        one=lambda key:self.equal(a[key],b[key])
        if op in ('sum','product','logical_or'):
            assert len(a['arguments'])==len(b['arguments'])
            left=sorted(a['arguments'],key=self.align_a.at);right=sorted(b['arguments'],key=self.align_b.at)
            for q,r in zip(left,right):self.equal(q,r)
        elif op in ('exact_rational','bound_variable','mathematical_pi',
            'current_original_source_parameter','exact_positive_integer_expression'):
            assert a==b,(i,j,'exact defining leaf')
        elif op=='negative':one('argument')
        elif op=='positive_quotient':
            # Both quotient domains are already admitted by their separately
            # hash-bound source receipts. Different prose certificate labels
            # do not alter their identical positive denominator function.
            assert a['source_positive_certificate'] and b['source_positive_certificate']
            one('numerator');one('denominator')
        elif op=='analytic_unary':assert a['name']==b['name'];one('argument')
        elif op=='exact_real_comparison':assert a['operator']==b['operator'];one('left');one('right')
        elif op=='original_lazy_flat_branch':
            for k in ('active_body','flat_predicate','flat_value'):one(k)
        elif op=='original_monotone_phase_inverse':
            children=('phase_function','target_phase','angle_lower','angle_upper','Delta','eta')
            assert {k:v for k,v in a.items() if k not in children}=={k:v for k,v in b.items() if k not in children}
            for k in children:one(k)
        elif op=='substitute_original_inverse_angle':
            assert a['angle_variable']==b['angle_variable'];one('body');one('inverse_angle')
        elif op=='function_substitution':
            assert a['Z_independent_substitution'] and b['Z_independent_substitution']
            for k in ('expression','variable','value'):one(k)
        elif op=='definite_integral':
            children=('integrand','lower','upper')
            assert {k:v for k,v in a.items() if k not in children}=={k:v for k,v in b.items() if k not in children}
            for k in ('integrand','lower','upper'):one(k)
        elif op=='current_original_leading_function_recipe':
            for k in ('source_family','native_chart','recipe','quantity','Z_order'):
                assert a[k]==b[k],(i,j,'actual source definition '+k)
            for k in ('Taylor_coefficient_factorial','quantity_path','coordinate_independent'):
                default=1 if k=='Taylor_coefficient_factorial' else False if k=='coordinate_independent' else None
                assert a.get(k,default)==b.get(k,default)
            for k in ('phase_independent_leading_source','leading_history_normalization_already_applied'):
                assert a.get(k)==b.get(k),(i,j,'source phase and normalization contract '+k)
            assert a['defining_quantity_not_a_range_value'] and b['defining_quantity_not_a_range_value']
            assert a.get('Taylor_coefficient_factorial',1)==(1,1,2,6)[a['Z_order']]
            one('coordinate');one('Z_variable')
        else:raise AssertionError('Unsupported actual correction dependency '+str(a))
        self.memo.add((i,j));self.operations[op]=self.operations.get(op,0)+1


def exact_inlet_relation(field):
    relation=RelationalFunctions(field.outer_nodes,field.band_nodes)
    relation.equal(field.outer_N,field.band_N)
    rows=[]
    for key,outer_rows in field.outer_terminal.items():
        for j,(left,right) in enumerate(zip(outer_rows,field.band_inlet[key])):
            expression=field.band_nodes[right]
            assert expression['operation']=='product' and len(expression['arguments'])==2
            inverses=[]
            for index in expression['arguments']:
                node=field.band_nodes[index]
                if node['operation']=='positive_quotient' and field.band_nodes[node['numerator']]==field.band_nodes[1]:
                    inverses.append(index)
            assert len(inverses)==1
            inverse=inverses[0];quotient=field.band_nodes[inverse]
            relation.equal(field.outer_N,quotient['denominator'])
            terminal=next(index for index in expression['arguments'] if index!=inverse)
            relation.equal(left,terminal)
            rows.append(dict(history=key,ordinary_Z_order=j,outer_N_scaled_node=left,
                repair_N_scaled_node=terminal,actual_repair_inlet_node=right,
                same_actual_terminal_defining_DAG=True,same_exact_N_divided_once=True))
    assert len(rows)==20
    assert relation.operations['definite_integral']>0
    assert relation.operations['original_monotone_phase_inverse']>0
    assert relation.operations['current_original_leading_function_recipe']>0
    assert relation.operations['analytic_unary']>0
    return dict(independent_original_correction_inlet_rows=rows,
        directly_compared_defining_node_pairs=len(relation.memo),defining_operations=relation.operations,
        fingerprints_only_align_children_and_never_establish_identity=True,
        exact_source_recipe_phase_inverse_Jacobian_integral_limits_and_memory_preserved=True)


def accepted_endpoint_and_units(field):
    receipt=current.outer.rh.read(current.outer.RECEIPT)
    proof=receipt['independent_partial_histories_and_partition']
    assert proof['independently_checked_new_partial_endpoint_equals_old_alias_rows']==100
    assert proof['integral_additivity_on_original_partitions'] and proof['physical_Jacobian_present_once']
    functions=field.outer_report['actual_five_chart_C3_continuation_functions']['O3_power']
    assert functions['exact_native_domain']==[0,2]
    target_report=current.outer.rh.read(current.target.NAME)
    target_rows=target_report['actual_C3_density_transport_targets']['terminal_N_scaled_histories']
    # The current outer Rc aliases belong to this accepted original target
    # graph, not to a reset or newly fitted endpoint packet.
    assert target_rows==field.outer_terminal
    assert tuple(target_rows)==tuple(current.current.RATES)
    inlet=field.correction_inlet_functions()
    assert current.target.encoded(inlet)==field.band_inlet
    assert all(len(rows)==4 for rows in field.band_inlet.values())
    assert not field.quiet_report['current_outer_leading_endpoint_band_seed_function_identity_installed']
    return dict(accepted_new_partial_Rc_endpoint_equals_original_alias=True,
        exact_original_terminal_N_scaled_target_rows_retained=True,
        source_owned_C3_repair_inlet_API=True,ordinary_derivative_factorials_once=True,
        complete_leading_history_interface_still_open=True)


def radius_phase_and_integer_seam(field):
    a=field.outer_nodes;b=field.band_nodes;q=field.quiet_nodes
    relation=RelationalFunctions(a,q)
    relation.equal(field.outer_N,field.quiet_N)
    geo=field.outer_report['actual_geometry_bindings']['O3_power']
    qf=field.quiet_report['actual_C3_repaired_exit_power_and_pulse_functions']
    relation.equal(geo['right_offset'],qf['original_Rc_offset'])
    controls=current.outer.rh.read(current.band.control.NAME)['actual_C3_limit_controls_and_Picard_functions']
    base=controls['exact_C2_limit_functions']['exact_C1_limit_adapter']
    bf=field.band_report['actual_C3_repair_band_functions'];offset=s.Symbol('original_Rc_offset');x=s.Symbol('repair_x')
    bindings={base['original_Rc_offset']:offset,bf['original_band_variable']:x}
    radius=rh_check.symbolic(field.band_graph,bf['original_leading_power']['original_radius_offset'],bindings)
    assert s.simplify(radius.subs(x,1)-offset)==0
    assert s.diff(radius,x).subs(x,1)==1
    assert a[geo['Jacobian']]==dict(operation='exact_rational',numerator=1,denominator=1)
    assert a[geo['domain'][1]]==dict(operation='exact_rational',numerator=2,denominator=1)
    relation2=RelationalFunctions(a,b);relation2.equal(geo['right_offset'],base['original_Rc_offset'])
    phase=a[geo['phase']]
    assert phase['operation']=='analytic_unary' and phase['name']=='fractional_part'
    product=a[phase['argument']]
    assert product['operation']=='product' and sorted(product['arguments'])==sorted((field.outer_N,geo['offset']))
    assert field.source_integer_and_Rc_binding['actual_selected_repair_integer']['J']==a[field.outer_N]['J']
    # R=Rc*x and the verified common logRc give exactly the same absolute
    # phase frac(N*logRc) at outer x=2 and repair x=1. No phase reset occurs.
    return dict(original_outer_x2_equals_repair_x1_log_radius=True,
        original_native_Jacobians_at_Rc_both_one=True,
        exact_same_selected_N_in_outer_repair_and_quiet=True,
        identical_absolute_phase_fractional_part_N_logRc=True,
        absolute_phase_not_reset_at_repair_inlet=True)


def run(field=None):
    began=time.monotonic();raw=current.outer.rh.read(current.NAME)
    for name,value in raw['input_hashes'].items():assert current.sha(name)==value,name
    field=field if field is not None else current.CurrentOuterToRepairC3Bridge(require_checked=False)
    assert not field.acceptance_loaded and field.identity==raw['source_family']
    assert raw['source_bindings']==current.source_bindings()
    assert raw['actual_correction_interface_function_identity']==field.identity_proof
    assert raw['actual_source_integer_and_Rc_binding']==field.source_integer_and_Rc_binding
    assert raw['candidate_actual_outer_to_repair_C3_correction_bridge_constructed']
    assert not any(raw[k] for k in current.GATES+current.OPEN)
    relation=exact_inlet_relation(field);endpoint=accepted_endpoint_and_units(field)
    seam=radius_phase_and_integer_seam(field)
    for key in ('current_outer_leading_endpoint_band_seed_function_identity_installed',
        'current_outer_complete_history_to_repair_inlet_function_identity_installed',
        'selected_native_pulse_constructor_consumes_current_C3_frame'):
        assert not raw[key]
    result=dict(all_passed=True,source_family=field.identity,**dict.fromkeys(current.GATES,True),
        independent_direct_function_relation=relation,accepted_original_endpoint_and_units=endpoint,
        independent_actual_radius_phase_and_integer_seam=seam,
        current_outer_leading_endpoint_band_seed_function_identity_installed=False,
        current_outer_complete_history_to_repair_inlet_function_identity_installed=False,
        selected_native_pulse_constructor_consumes_current_C3_frame=False,
        current_numeric_point_field_oracle_installed=False,
        global_physical_time_Cartesian_heat_cone_and_temporal_recursion_installed=False,
        **dict.fromkeys(current.OPEN,False),
        input_hashes={**field.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('Actual outer-to-repair C3 correction interface defining-DAG checks passed',flush=True)
    return result


if __name__=='__main__':run()
