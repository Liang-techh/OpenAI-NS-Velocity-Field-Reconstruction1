"""Signed reference/restore viscous chain with independently saved inlet cap.

Only saved packets are hydrated.  Source functions telescope first; positive
bounds are subsequently attached to the two surviving external functions.
No cap, decimal proxy, or interval midpoint is a signed function evaluation.
"""
import ast
import copy
import hashlib
import json
import sys

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from .original_core_numeric import ROOT, DIRECTORY, PREFIX, _bounds
from .original_active_viscous_m15 import OriginalActiveViscousM15, ActiveEndpointFingerprint
from .original_micro_viscous_m15 import contract_trace_chain


NATIVE = PREFIX+'compliant_current_original_reference_native_m15.json'
NATIVE_SHA = 'e8177ddf61b1a4b139afcdb3fa2b3b95c045a98bc2072ed603ec79ecf8b66293'
ACTIVE = PREFIX+'compliant_current_original_active_viscous_m15_numeric.json'
ACTIVE_SHA = '268f00eea0e118d2afe5b6e9c82aecb3b2c43ceae30887cea4da53965c34bb72'
REPLACEMENT = 'reference_restoration_postrestore_finite_N_viscous_external_endpoints'
OLD_COMBINED = 'restoration_postrestore_finite_N_viscous_external_endpoints'


def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()


class OriginalReferenceEndpointM15(OriginalActiveViscousM15):
    def __init__(self,*,directory=DIRECTORY,precision=150):
        # Deliberately do not call any accepted ancestor constructor.
        self.directory=directory;self.ctx=mp.mp.clone();self.ctx.dps=precision
        self.iv=MPIntervalContext();self.iv.dps=precision+30
        self.hashes={};self.code_hashes={}
        self.native=self.packet(NATIVE,NATIVE_SHA)
        self.family=self.native['source_family']
        self.parent=self.packet(ACTIVE,ACTIVE_SHA)
        self.reference=self.saved('reference_m15_numeric')
        self.restore=self.saved('restore_m15_numeric')
        self.modified=self.archive('modified_Q4_pressure_C3_functions')
        self.loop=self.archive('prescribed_loop_C5_functions')
        self.outer=self.archive('prescribed_R100_Rm_C6_functions')
        self.m=self.modified['modified_Q4_pressure_C3_bindings']
        self.l=self.loop['prescribed_loop_C5_bindings']
        self.background=self.outer['prescribed_R100_Rm_C6_bindings']['charts']
        self.nodes=self.modified['exact_graph_nodes']
        for source in (self.loop,self.outer):
            if self.nodes[:len(source['exact_graph_nodes'])]!=source['exact_graph_nodes']:
                raise ValueError('Common original graph prefixes required')
        self.source_binding=self.bind_sources()
        self.pressure_envelope_binding=self.bind_reference0_envelope()
        self.seams=self.bind_seams()
        name='src/openai_ns_reconstruction/original_reference_endpoint_m15.py'
        self.code_hashes[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    def packet(self,name,expected):
        raw=(self.directory/name).read_bytes()
        if hashlib.sha256(raw).hexdigest()!=expected:raise ValueError('Changed accepted parent '+name)
        data=json.loads(raw)
        if not data['checks']['passed']:raise ValueError('Checked parent required')
        self.merge(data['input_hashes'],self.directory,self.hashes)
        self.merge(data['implementation_hashes'],ROOT,self.code_hashes)
        self.hashes[name]=expected
        return data

    def bind_sources(self):
        requests={
            ('original_macro_m15.py','OriginalMacroM15','correction_densities'):[
                "errors={name:[self.rows.read_log(v) for v in packet['physical_pressure_rows'][name]['finite_N_absolute_error_log_bounds']] for name in ('U','Q','Q_y')}",
                "return dict(endpoint=self.rows.combine(bc,values,1),volume=self.rows.combine(vc,values,1),viscous=errors['Q_y'])"],
            ('original_macro_pressure_rows.py','OriginalMacroPressureRows','evaluate'):[
                "slow_y=[None if v is None else v+lp-self.micro.logN for v in slow_y]",
                "fast_y=[None if v is None else v+lp for v in fast_y]",
                "errors['U_y']=[self.sum_logs(v) for v in zip(slow_y,fast_y)]",
                "errors['Q_y']=self.combine(coeff['Q_y'],errors,1)"],
            ('original_reference_m15.py','OriginalReferenceM15','finite_errors'):[
                "lowpacket=owner.evaluate(Z,'1')",
                "owner=self.rows.long_rows"],
            ('original_restore_m15.py','OriginalRestoreM15','finite_errors'):[
                "lower_packet=owner.evaluate(Z,'1')",
                "upper=self.correction_densities(Z,packet,base)",
                "lower=self.correction_densities(Z,lower_packet,lower_base)",
                "result[name]=[self.rows.sum_logs((None if upper[name][k] is None else upper[name][k]+power*geo['logR'],None if lower[name][k] is None else lower[name][k]+power*lower_geo['logR'])) for k in (0,1)]"],
            ('original_first_micro_pressure_rows.py','OriginalFirstMicroPressureRows','sum_logs'):[
                "values=[v for v in values if v is not None]",
                "high=max(_bounds(self.ctx,v)[1] for v in values)",
                "return self.iv.mpf(high)+self.iv.log(len(values))"],
        }
        rows=[]
        for (file,cls,method),statements in requests.items():
            name='src/openai_ns_reconstruction/'+file;raw=(ROOT/name).read_bytes()
            if hashlib.sha256(raw).hexdigest()!=self.code_hashes.get(name):raise ValueError('Pinned producer required')
            tree=ast.parse(raw);body=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name==cls)
            fn=next(n for n in body.body if isinstance(n,ast.FunctionDef) and n.name==method)
            dumps={ast.dump(n) for n in ast.walk(fn)}
            for statement in statements:
                if ast.dump(ast.parse(statement).body[0]) not in dumps:raise ValueError('Changed producer '+statement)
            rows.append(dict(source=name,source_sha256=self.code_hashes[name],method=cls+'.'+method,
                             exact_source_statements=statements))
        return dict(passed=True,producer_rows=rows,
            function='D=R*(Qy_modified-Qy_original_leading); chart contribution=D(lower)-D(upper)',
            physical_radius_power=1,ordinary_Z_orders=[0,1],
            positive_sum_domination='exp(max(log_bounds)+log(count)) bounds the sum of the nonnegative magnitudes',
            ordinary_Z1_product_rule_retained=True,physical_Pstar_already_in_Qy_error=True,
            slow_B_y_over_N_and_fast_B_phi_both_retained=True,
            separate_inlet_uses_own_reference0_pressure_packet=True,
            old_reference_lower_packet_is_longreshape1_provenance_only=True,
            no_longreshape_reference_alias_authorized=True,
            postrestore_entire_original_155_cap_retained=True,
            inherited_pressure_bound_validity_not_independently_recertified=True)

    def bind_reference0_envelope(self):
        """Identify the error propagator's function, without evaluating it.

        Qy is linear in Uy,U,m,Uz,mz with coefficients depending only on
        delta,Z.  Thus this binding does not require choosing the leading
        numerical pressure proxy or identifying reference0 with long1.
        """
        name='src/openai_ns_reconstruction/original_reference_pressure_rows.py'
        raw=(ROOT/name).read_bytes()
        if hashlib.sha256(raw).hexdigest()!=self.code_hashes[name]:raise ValueError('Pinned reference pressure wrapper required')
        tree=ast.parse(raw);cls=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=='OriginalReferencePressureRows')
        owners=next(n.value for n in cls.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='OWNERS' for t in n.targets))
        if ast.dump(owners)!=ast.dump(ast.parse("dict(reference='inner_reference',restoration='axial_restore',postrestore='restore_buffer')",mode='eval').body):
            raise ValueError('Actual reference derivative owner changed')
        init=next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name=='__init__')
        statements=["self.chart,self.derivative_chart=chart,self.OWNERS[chart]",
            "view._covers=lambda Z,q:self.reference._covers(Z,self.chart,q)",
            "self.loop=json.loads(path.read_bytes())['current_actual_loop_jet_log_bounds_by_chart'][self.derivative_chart]"]
        dumps={ast.dump(n) for n in ast.walk(init)}
        for statement in statements:
            if ast.dump(ast.parse(statement).body[0]) not in dumps:raise ValueError('Changed direct reference envelope owner')
        coeff_name='src/openai_ns_reconstruction/original_first_micro_pressure_rows.py'
        coeff_raw=(ROOT/coeff_name).read_bytes()
        if hashlib.sha256(coeff_raw).hexdigest()!=self.code_hashes[coeff_name]:raise ValueError('Pinned Qy propagator required')
        functions=ast.parse(coeff_raw);fn=next(n for n in functions.body if isinstance(n,ast.FunctionDef) and n.name=='pressure_coefficients')
        ret=next(n for n in ast.walk(fn) if isinstance(n,ast.Return))
        qy=next(k.value for k in ret.value.keywords if k.arg=='Q_y')
        wanted="{'U_y':2*z/L,'U':-z*(1-delta)/L,'m':z*(1-delta)/L,'U_Z':-d/L,'m_Z':d/L}"
        if ast.dump(qy)!=ast.dump(ast.parse(wanted,mode='eval').body):raise ValueError('Actual Qy linear coefficients changed')
        coefficient_statements=["z=J.variable(ctx,ctx.mpf(Z),order)","L=1-z*z*delta","d=1-z*z"]
        for statement in coefficient_statements:
            if ast.dump(ast.parse(statement).body[0]) not in {ast.dump(n) for n in ast.walk(fn)}:
                raise ValueError('Actual Qy coefficient geometry changed')
        combine=next(n for n in ast.walk(functions) if isinstance(n,ast.FunctionDef) and n.name=='combine')
        ordinary="ordinary=coef[j]*(math.factorial(j)*math.comb(k,j))"
        if ast.dump(ast.parse(ordinary).body[0]) not in {ast.dump(n) for n in ast.walk(combine)}:
            raise ValueError('Ordinary Z product rule required')
        guards=self.native['native_derivative_source_guards']['checks']
        formulas=guards['source_ast']['rows']
        if not all(formulas[k]['passed'] for k in ('pressure_formula','radial_formula')):
            raise ValueError('Actual archived physical pressure formula required')
        radial=guards['slow_B_to_total_Vy_to_physical_Uy']['charts']['reference']
        if not radial['passed']:raise ValueError('Actual reference slow/fast derivative binding required')
        sys.setrecursionlimit(max(sys.getrecursionlimit(),12000))
        chart=self.l['charts']['reference']
        f=ActiveEndpointFingerprint(self.nodes,dict(coordinate=chart['coordinate'],loop=self.l['charts']))
        endpoint=chart['domain'][0];scope={chart['coordinate']:f.at(endpoint,{})}
        scoped=[]
        for k in (0,1):
            qy_node=self.m['charts']['reference']['physical']['Q_y'][k]
            rqy_node=self.m['charts']['reference']['finite_boundary']['R_Q_y'][k]
            scoped.append(dict(ordinary_Z_order=k,Qy_function_node=qy_node,RQy_function_node=rqy_node,
                Qy_at_reference0_sha256=f.at(qy_node,scope)[0],RQy_at_reference0_sha256=f.at(rqy_node,scope)[0]))
        return dict(passed=True,method='source_bound_linear_Qy_envelope_at_actual_reference_native0',
            endpoint_scope=dict(coordinate_node=chart['coordinate'],endpoint_node=endpoint),
            graph_function_rows=scoped,reference_wrapper_source=name,reference_wrapper_sha256=self.code_hashes[name],
            exact_reference_owner_statements=statements,numeric_Qy_coefficient_source=coeff_name,
            exact_coefficient_geometry_statements=coefficient_statements,
            exact_numeric_Qy_coefficients=wanted,ordinary_Z_product_rule=ordinary,
            actual_archived_pressure_formula_certificate=copy.deepcopy(formulas['pressure_formula']),
            actual_archived_radial_formula_certificate=copy.deepcopy(formulas['radial_formula']),
            actual_reference_slow_fast_and_Pstar_certificate=copy.deepcopy(radial),
            deduction='Expand my=U-m and my_Z=Uz-mz in the pinned archived Qy formula; coefficients are exactly the numeric linear error propagator. Apply the pinned ordinary-Z product rule to its C1 row.',
            error_inputs=['delta_Uy','delta_U','delta_m','delta_Uz','delta_mz'],
            delta_inputs_from_same_saved_reference_pressure_packet=True,
            own_reference_generic_derivative_owner='inner_reference',
            physical_R_applied_after_Qy_error=True,source_field_values_or_caps_not_selected=True,
            source_function_not_identified_with_longreshape1=True,
            inherited_generic_and_history_enclosures_not_recomputed=True)

    def bind_seams(self):
        recovery=self.native['actual_modified_trace_recovery']
        theorem=self.native['native_derivative_theorem']
        if not recovery['passed'] or not theorem['passed']:raise ValueError('Actual native source theorem required')
        history=self.native['corrected_history_seam']
        if not history['passed'] or not all(history['required_orders'].values()):raise ValueError('History and leading jets required')
        first=[]
        for k in (0,1):
            row=next(r for r in recovery['rows'] if r['quantity']=='R_Q_y' and r['ordinary_Z_order']==k)
            if [row['left_node'],row['right_node']]!=[self.m['charts'][n]['finite_boundary']['R_Q_y'][k] for n in ('reference','restoration')]:
                raise ValueError('Actual RQy function ports required')
            first.append(copy.deepcopy(row))
        second=self.parent['samples'][0]['source_seam_proofs'][0]
        if not second['passed'] or (second['left'],second['right'])!=('restoration','postrestore'):
            raise ValueError('Actual outgoing seam proof required')
        if not second['actual_modified_and_original_leading_traces_equal']:raise ValueError('Leading subtraction required')
        return dict(passed=True,reference_restoration_modified_RQy_rows=first,
            reference_restoration_leading_jet_certificate_sha256=digest(history),
            reference_restoration_native_theorem_sha256=digest(theorem),
            restoration_postrestore_certificate=copy.deepcopy(second),
            delta_definition='R*(Qy_modified-Qy_original_leading)',
            leading_Qy_recovery='Pstar*[2Z*Vy-(1-delta)*Z*(V-m)-(1-Z^2)*(Vz-mz)]/(1-delta*Z^2)',
            leading_inputs_common_through_Z2=True,physical_radius_common_and_Z_independent=True,
            fast_phase_not_reset=True,source_functions_not_error_caps=True)

    def function_ports(self):
        ports=[]
        for name in ('reference','restoration','postrestore'):
            chart=self.l['charts'][name]
            bg=self.background['post_restore' if name=='postrestore' else name]
            visited=set()
            def z_free(node):
                if node in visited:return
                visited.add(node);row=self.nodes[node];op=row['operation']
                if op=='bound_variable':
                    if node!=chart['coordinate']:raise ValueError('Unexpected radius variable')
                elif op in ('original_source_parameter','current_original_source_parameter','exact_positive_integer_expression'):
                    if not row.get('Z_independent'):raise ValueError('Radius must be ordinary-Z independent')
                elif op in ('sum','product'):
                    for value in row['arguments']:z_free(value)
                elif op in ('negative','analytic_unary'):z_free(row['argument'])
                elif op=='positive_quotient':z_free(row['numerator']);z_free(row['denominator'])
                elif op not in ('exact_rational','mathematical_pi'):raise ValueError('Unissued radius operation '+op)
            z_free(chart['radius']);z_free(chart['physical_measure'])
            for node in self.m['charts'][name]['finite_boundary']['R_Q_y'][:2]:
                row=self.nodes[node]
                if row['operation']!='product' or row['arguments'].count(chart['radius'])!=1:
                    raise ValueError('One explicit original R factor required')
            for side,endpoint in zip(('lower','upper'),chart['domain']):
                ports.append(dict(chart=name,side=side,coordinate_node=chart['coordinate'],endpoint_node=endpoint,
                    endpoint_scope=dict(coordinate_node=chart['coordinate'],endpoint_node=endpoint),
                    physical_radius_node=chart['radius'],native_Jacobian_node=chart['physical_measure'],
                    modified_RQy_nodes=self.m['charts'][name]['finite_boundary']['R_Q_y'][:2],
                    leading_normalized_V_nodes=bg['fields']['V'][:3],
                    leading_normalized_mean_nodes=bg['fields']['m'][:3],
                    leading_normalized_Vy_nodes=bg['radial']['V_y'][:2],
                    definition='D=modified_RQy-R*Pstar*leading_normalized_Qy',
                    radius_and_measure_Z_free_dependency_nodes=sorted(visited),
                    physical_R_explicit_once_in_both_Z_rows=True,
                    ordinary_Z_orders=[0,1],signed_coefficient=1 if side=='lower' else -1))
        return ports

    @staticmethod
    def endpoint(data,Z,chart):
        return next(s for s in data['samples'] if s['Z']==Z and s.get('chart')==chart and
                    s['coordinate']==dict(kind='exact_decimal_'+chart,value='1'))

    def inlet_bound(self,Z,selected_integer):
        row=next(s for s in self.reference['pressure_row_samples'] if s['Z']==Z and
                 s['coordinate']==dict(kind='exact_decimal_reference',value='0'))
        if row['selected_integer']!=selected_integer or row['source_family']!=self.family:
            raise ValueError('Own reference0 same-N pressure source required')
        if row['source_cover']['generic_derivative_chart']!='inner_reference':raise ValueError('Own reference derivative owner required')
        if row['source_cover']['cumulative_history_charts'][-1]!='reference':raise ValueError('Reference histories retained')
        if not row['contract']['original_Pstar_once']:raise ValueError('Physical Pstar once required')
        endpoint=row['R_Q_y_endpoint'];errors=row['physical_pressure_rows']['Q_y']['finite_N_absolute_error_log_bounds']
        if errors!=endpoint['finite_N_coefficient_error_log_bounds']:raise ValueError('Same Qy budget required')
        radius=endpoint['finite_N_log_physical_factor']
        if radius!=row['geometry']['log_physical_radius']:raise ValueError('Own physical R factor required')
        weighted=[None if v['exact_zero'] else self.decode(v['log_absolute_upper'])+self.decode(radius) for v in errors]
        if len(weighted)!=2:raise ValueError('Ordinary C0/Z1 bounds required')
        return weighted,dict(source_report=PREFIX+'compliant_current_original_reference_m15_numeric.json',
            source_packet_sha256=digest(row),coordinate=copy.deepcopy(row['coordinate']),
            selected_integer=copy.deepcopy(selected_integer),source_cover=copy.deepcopy(row['source_cover']),
            actual_own_Qy_error_rows=copy.deepcopy(errors),actual_own_log_radius=copy.deepcopy(radius),
            actual_linear_error_input_rows={name:copy.deepcopy(row['physical_pressure_rows'][name]['finite_N_absolute_error_log_bounds'])
                                           for name in ('U','m','U_y')},
            actual_slow_fast_error_components=copy.deepcopy(row['finite_N_axial_radial_error_components']),
            radius_weighted_C0_C1_log_bounds=[self.record(v) if v is not None else None for v in weighted],
            actual_reference0_graph_envelope_binding=self.pressure_envelope_binding,
            physical_R_once=True,physical_Pstar_already_in_error=True,extra_epsilon=False,
            source_function_values_evaluated=False,old_longreshape_lower_packet_not_used=True)

    def evaluate(self,Z):
        parent=next((s for s in self.parent['samples'] if s['Z']==Z),None)
        if parent is None:raise ValueError('Only saved axial queries admitted; ancestors cannot replay')
        raw=parent['raw_m15_C0_C1_terms'];old=parent['correlated_reduced_m15_C0_C1_terms']
        ref=self.endpoint(self.reference,Z,'reference');rest=self.endpoint(self.restore,Z,'restoration');post=self.endpoint(self.restore,Z,'postrestore')
        if len(raw)!=213 or len(old)!=199 or raw[:124]!=ref['m15_C0_C1_terms'] or raw[:140]!=rest['m15_C0_C1_terms'] or raw[:156]!=post['m15_C0_C1_terms']:
            raise ValueError('Exact cumulative source prefixes required')
        if not ref['actual_selected_integer']==rest['actual_selected_integer']==post['actual_selected_integer']:
            raise ValueError('Common actual N required')
        for i,name in ((123,'reference'),(139,'restoration'),(155,'postrestore')):
            if raw[i]['name']!=name+'_finite_N_viscous_error':raise ValueError('Original sector identity required')
        inlet,witness=self.inlet_bound(Z,ref['actual_selected_integer'])
        outgoing=raw[155]
        if len(outgoing['log_scale_terms'])!=1 or not outgoing['log_scale_terms'][0]['positive_error_bound_atom_only']:
            raise ValueError('Entire inherited postrestore aggregate required')
        # No subtraction of either an interval or a positive cap.  Full155
        # bounds its own upper endpoint by the pinned positive producer.
        bound=self.positive_union(inlet+[self.decode(outgoing['log_scale_terms'][0])])
        replacement=copy.deepcopy(outgoing);replacement['name']=REPLACEMENT
        replacement['log_scale_terms']=[self.record(bound)]
        replacement['actual_function_binding']=dict(error_enclosure_only=True,
            function='D(reference native0)-D(postrestore native1=Rm)',delta_definition=self.seams['delta_definition'],
            ordinary_Z1=True,physical_R_once=True,external_endpoints_not_zeroed=True,
            signed_centers_evaluated=False,no_error_cap_subtraction=True,
            independent_own_reference0_budget=True,full_original_postrestore_155_budget_retained=True)
        reduced=[];removed=[]
        for term in old:
            if term['name']==raw[123]['name']:removed.append(term['name']);continue
            if term['name']==OLD_COMBINED:removed.append(term['name']);reduced.append(copy.deepcopy(replacement))
            else:reduced.append(copy.deepcopy(term))
        if set(removed)!={raw[123]['name'],OLD_COMBINED} or len(reduced)!=198:raise ValueError('Exact new198 view required')
        unchanged=[t for t in old if t['name'] not in removed]
        if [t for t in reduced if t['name']!=REPLACEMENT]!=unchanged:raise ValueError('Other197 sectors must stay verbatim')
        chain=contract_trace_chain([
            dict(name='reference',lower='reference0',upper='reference1'),
            dict(name='restoration',lower='restoration0',upper='restoration1'),
            dict(name='postrestore',lower='postrestore0',upper='Rm')],
            dict(restoration0='reference1',postrestore0='restoration1'))
        return dict(Z=Z,source_family=self.family,raw_m15_C0_C1_terms=copy.deepcopy(raw),
            correlated_reduced_m15_C0_C1_terms=reduced,original_raw_terms=213,previous_reduced_terms=199,
            correlated_reduced_terms=198,unchanged_other_sectors=197,signed_trace_chain=chain,
            source_seam_proofs=self.seams,actual_signed_function_ports=self.function_ports(),
            external_endpoint_budget_witness=dict(own_reference0=witness,
                postrestore_full_original_155_aggregate=copy.deepcopy(outgoing),positive_union=self.record(bound),
                source_binding=self.source_binding,postrestore_upper_separate_budget_claimed=False),
            old_budget_terms_retained_as_provenance=[copy.deepcopy(raw[i]) for i in (123,139,155)],
            replacement_outer_endpoints=replacement,independent_target=copy.deepcopy(parent['independent_target']),
            largest_remaining_individual_log_caps=self.largest_caps(reduced,parent['independent_target']),
            contract=dict(raw_213_verbatim=True,other_197_verbatim=True,independent_target_unchanged=True,
                reserved_added_this_stage=0,ordinary_Z1_and_physical_R_once=True,leading_subtraction_retained=True,
                both_external_source_functions_retained=True,no_longreshape_reference_alias=True,
                all_other_volume_width_numerical_errors_retained=True,ancestors_not_reexecuted=True,
                source_bound_signed_trace_contraction=True,error_cap_as_value=False),
            limits=dict(reference_restoration_internal_trace_resolved=True,
                inherited_pressure_bounds_independently_recertified=False,
                external_signed_values_resolved=False,longreshape_reference_seam_proved=False,
                independent_numeric_Rp_seam=False,global_five_moment_common_N_closure=False,
                exact_heat_and_energy=False,admissible_stress_flat_remainder=False,
                temporal_n_recursion_and_summation=False,mean_oscillatory_correction=False,
                independent_corrected_Cartesian_NS=False))

    def build(self):
        return dict(source_family=self.family,parent_report=NATIVE,parent_report_sha256=NATIVE_SHA,
            numerical_parent_report=ACTIVE,numerical_parent_report_sha256=ACTIVE_SHA,
            source_binding=self.source_binding,source_seam_proofs=self.seams,
            reference0_graph_envelope_binding=self.pressure_envelope_binding,
            samples=[self.evaluate(Z) for Z in ('-.3','0','.3')],
            contract=dict(saved_packets_only=True,accepted_ancestor_constructors_not_run=True,
                          cap_is_not_signed_function_value=True),
            limits=dict(numeric_Rp_closure=False,global_all_five_closure=False,full_goal_complete=False))
