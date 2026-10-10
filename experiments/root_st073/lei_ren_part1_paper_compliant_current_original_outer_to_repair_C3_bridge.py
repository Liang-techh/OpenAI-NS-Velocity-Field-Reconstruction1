"""Portable exact correction interface from the O3 Rc exit to compact repair.

The outer API is N-scaled. The original repair inlet divides that same
function by the same N. Complete leading-history continuity remains separate.
"""
import hashlib
import gzip
import json
from pathlib import Path
import time
import lei_ren_part1_paper_compliant_current_original_C3_power_to_Rp as quiet

outer,band,target,current,source=quiet.outer,quiet.band,quiet.target,quiet.current,quiet.source
HERE,PREFIX,sha=quiet.HERE,quiet.PREFIX,quiet.sha
NAME=PREFIX+'current_original_outer_to_repair_C3_bridge.json.gz'
RECEIPT=PREFIX+'current_original_outer_to_repair_C3_bridge_check.json'
GATES=('current_original_outer_Rc_to_repair_inlet_C3_correction_function_identity_installed',)
OPEN=quiet.OPEN


def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()


class ExactFunctionSignatures:
    """Fail-closed defining expressions, independent of graph numbering.

Only commutative sum/product order and non-defining proof annotations are
discarded. Source recipes, phase, Jacobian, integration limits, inverse branch,
ordinary derivative order and the exact global integer remain source-defined.
"""
    def __init__(self,nodes):self.nodes,self.memo,self.definitions=nodes,{},{}

    def at(self,i):
        if i in self.memo:return self.memo[i]
        n=self.nodes[i];op=n['operation'];at=self.at
        if op=='exact_rational':body=(op,n['numerator'],n['denominator'])
        elif op=='bound_variable':body=(op,n['name'])
        elif op in ('sum','product','logical_or'):
            body=(op,sorted(at(j) for j in n['arguments']))
        elif op=='negative':body=(op,at(n['argument']))
        elif op=='positive_quotient':body=(op,at(n['numerator']),at(n['denominator']))
        elif op=='analytic_unary':body=(op,n['name'],at(n['argument']))
        elif op=='mathematical_pi':body=(op,)
        elif op in ('current_original_source_parameter','exact_positive_integer_expression'):
            body=n
        elif op=='exact_real_comparison':body=(op,n['operator'],at(n['left']),at(n['right']))
        elif op=='original_lazy_flat_branch':body=(op,at(n['active_body']),at(n['flat_predicate']),at(n['flat_value']))
        elif op=='original_monotone_phase_inverse':
            children=('phase_function','target_phase','angle_lower','angle_upper','Delta','eta')
            body=(op,[at(n[k]) for k in children],{k:v for k,v in n.items() if k not in children})
        elif op=='substitute_original_inverse_angle':body=(op,at(n['body']),at(n['inverse_angle']),n['angle_variable'])
        elif op=='definite_integral':
            body=(op,at(n['integrand']),at(n['lower']),at(n['upper']),n['variable'],
                {k:v for k,v in n.items() if k not in ('integrand','lower','upper','variable')})
        elif op=='function_substitution':
            if not n['Z_independent_substitution']:raise ValueError('Original Z-independent endpoint required')
            body=(op,at(n['expression']),at(n['variable']),at(n['value']),n['Z_independent_substitution'])
        elif op=='current_original_leading_function_recipe':
            j=n['Z_order'];factor=n.get('Taylor_coefficient_factorial',1)
            if j not in range(4) or factor!=(1,1,2,6)[j]:raise ValueError('Ordinary source derivative normalization required')
            if not n['defining_quantity_not_a_range_value']:raise ValueError('Actual original source recipe required')
            body=(op,n['source_family'],n['native_chart'],n['recipe'],n['quantity'],j,factor,
                at(n['coordinate']),at(n['Z_variable']),n.get('quantity_path'),n.get('coordinate_independent',False),
                n.get('phase_independent_leading_source'),n.get('leading_history_normalization_already_applied'))
        else:raise ValueError('Unsupported exact interface dependency '+str(n))
        signature=digest(body);self.memo[i]=signature;self.definitions[signature]=body
        return signature


def accepted_inputs():
    reports,hashes={},{}
    for module in (outer,band,quiet):
        report,receipt=outer.rh.read(module.NAME),outer.rh.read(module.RECEIPT)
        if not receipt['all_passed'] or not all(receipt[k] for k in module.GATES):
            raise ValueError('Accepted actual interface input required '+module.RECEIPT)
        if report['source_family']!=receipt['source_family']:raise ValueError('Report and receipt source identity required')
        for name,value in receipt['input_hashes'].items():
            if name in hashes and hashes[name]!=value:raise ValueError('Conflicting original dependency '+name)
            hashes[name]=value
        hashes[module.NAME]=sha(module.NAME);hashes[module.RECEIPT]=sha(module.RECEIPT)
        reports[module.NAME]=report
    identity=reports[outer.NAME]['source_family']
    if any(row['source_family']!=identity for row in reports.values()):raise ValueError('Same original source required')
    hashes[Path(__file__).name]=sha(Path(__file__).name)
    for name,value in hashes.items():
        if sha(name)!=value:raise ValueError('Changed interface input '+name)
    return reports,hashes


def compare_functions(field):
    a=ExactFunctionSignatures(field.outer_nodes);b=ExactFunctionSignatures(field.band_nodes)
    rows=[]
    if (set(field.outer_terminal)!=set(current.RATES) or set(field.band_inlet)!=set(current.RATES)
        or any(len(field.outer_terminal[key])!=4 or len(field.band_inlet[key])!=4 for key in current.RATES)):
        raise ValueError('Exactly five histories with all four ordinary C3 rows required')
    for key,outer_rows in field.outer_terminal.items():
        for j,(node,inlet) in enumerate(zip(outer_rows,field.band_inlet[key])):
            # Compare the archived repair expression directly to the archived
            # outer terminal / N, without importing an ID from another graph.
            expected=digest(('product',sorted((a.at(node),digest(('positive_quotient',a.at(1),a.at(field.outer_N)))))))
            got=b.at(inlet)
            if expected!=got:raise ValueError('Original correction interface expression differs '+key+' Z'+str(j))
            rows.append(dict(history=key,ordinary_Z_order=j,outer_N_scaled_node=node,
                actual_repair_inlet_node=inlet,exact_outer_N_scaled_divided_by_original_N_signature=expected,
                actual_repair_inlet_signature=got,same_defining_function=True))
    if len(rows)!=20:raise ValueError('All twenty original correction interface rows required')
    return dict(rows=rows,outer_defining_expression_nodes=len(a.memo),repair_defining_expression_nodes=len(b.memo),
        exact_same_global_N_signature=a.at(field.outer_N),
        ordinary_Z_orders=[0,1,2,3],five_history_keys=list(field.outer_terminal),
        outer_N_scaled_equals_N_times_actual_repair_inlet=True,
        preserved_phase_Jacobian_signed_density_and_own_rate_memory=True,
        graph_node_numbering_and_range_overlap_not_used_as_identity=True)


class CurrentOuterToRepairC3Bridge:
    def __init__(self,require_checked=True):
        reports,self.hashes=accepted_inputs();self.identity=reports[outer.NAME]['source_family']
        self.outer_report=reports[outer.NAME];self.band_report=reports[band.NAME]
        self.quiet_report=reports[quiet.NAME];self.outer_nodes=self.outer_report['exact_graph_nodes']
        self.band_nodes=self.band_report['exact_graph_nodes'];self.band_graph=outer.rh.restore_graph(self.band_report)
        functions=self.outer_report['actual_five_chart_C3_continuation_functions']['O3_power']
        self.outer_terminal=functions['original_full_outgoing_N_scaled_C3']
        self.band_inlet=self.band_report['actual_C3_repair_band_functions']['actual_incoming_correction_C3']
        self.outer_N=self.outer_report['actual_geometry_bindings']['O3_power']['selected_N']
        controls=outer.rh.read(band.control.NAME)['actual_C3_limit_controls_and_Picard_functions']
        self.band_N=controls['exact_C2_limit_functions']['exact_C1_limit_adapter']['N']
        a,b=ExactFunctionSignatures(self.outer_nodes),ExactFunctionSignatures(self.band_nodes)
        qf=self.quiet_report['actual_C3_repaired_exit_power_and_pulse_functions']
        self.quiet_N=qf['original_N'];self.quiet_nodes=self.quiet_report['exact_graph_nodes']
        q=ExactFunctionSignatures(self.quiet_nodes)
        if any(nodes[node]['operation']!='exact_positive_integer_expression' for nodes,node in
            ((self.outer_nodes,self.outer_N),(self.band_nodes,self.band_N),(self.quiet_nodes,self.quiet_N))):
            raise ValueError('Original exact positive integer expression required in every consumer')
        if not a.at(self.outer_N)==b.at(self.band_N)==q.at(self.quiet_N):raise ValueError('Same exact source integer required')
        selected=self.band_report['actual_selected_repair_integer']
        integer=self.band_nodes[self.band_N]
        if selected['J']!=integer['J']:raise ValueError('Actual archived integer selection must define this exact N')
        rc=qf['original_Rc_offset'];right=self.outer_report['actual_geometry_bindings']['O3_power']['right_offset']
        if a.at(right)!=q.at(rc):raise ValueError('Same original O3 right radius Rc required')
        self.source_integer_and_Rc_binding=dict(original_N_definition=integer,actual_selected_repair_integer=selected,
            outer_N_node=self.outer_N,band_N_node=self.band_N,quiet_N_node=self.quiet_N,
            outer_Rc_offset_node=right,quiet_Rc_offset_node=rc,same_exact_original_N_and_Rc=True)
        if self.outer_report['source_family']!=self.band_report['source_family']:raise ValueError('Same source family required')
        self.identity_proof=compare_functions(self);self.acceptance_loaded=False
        if require_checked:
            receipt=outer.rh.read(RECEIPT)
            if not receipt['all_passed'] or not all(receipt[k] for k in GATES):raise ValueError('Checked original correction bridge required')
            for name,value in receipt['input_hashes'].items():
                if sha(name)!=value:raise ValueError('Changed correction bridge dependency '+name)
            self.acceptance_loaded=True

    def correction_inlet_functions(self):
        """Actual source-owned band graph handles; no extra normalization."""
        return {key:target.C3Function(*(source.FunctionRef(self.band_graph,i) for i in rows))
            for key,rows in self.band_inlet.items()}


def source_bindings():
    return dict(exact_defining_expressions=current.ast_binding(ExactFunctionSignatures),
        original_cross_graph_identity=current.ast_binding(compare_functions),
        original_outer_chain=current.ast_binding(outer.build_chart),
        original_repair_band=current.ast_binding(band.build),
        original_signed_C3_transport=current.ast_binding(target.build_transport))


def run():
    began=time.monotonic();field=CurrentOuterToRepairC3Bridge(require_checked=False)
    report=dict(candidate_actual_outer_to_repair_C3_correction_bridge_constructed=True,source_family=field.identity,
        actual_correction_interface_function_identity=field.identity_proof,source_bindings=source_bindings(),
        actual_source_integer_and_Rc_binding=field.source_integer_and_Rc_binding,
        correction_units='outer terminal is N-scaled; original repair inlet is the same function divided by N once',
        current_outer_leading_endpoint_band_seed_function_identity_installed=False,
        current_outer_complete_history_to_repair_inlet_function_identity_installed=False,
        selected_native_pulse_constructor_consumes_current_C3_frame=False,
        current_numeric_point_field_oracle_installed=False,
        global_physical_time_Cartesian_heat_cone_and_temporal_recursion_installed=False,
        **dict.fromkeys(GATES+OPEN,False),input_hashes=field.hashes,execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(report,separators=(',',':'))+'\n').encode(),mtime=0))
    print('Original outer-to-repair correction interface constructed',flush=True)
    return field


if __name__=='__main__':run()
