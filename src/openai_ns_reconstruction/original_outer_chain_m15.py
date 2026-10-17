"""Source-bound six-interface outer viscous m15 chain.

Only saved archives are hydrated. Accepted ancestor owners never run.
"""
import ast
import hashlib
import json
import sys

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from .original_core_numeric import ROOT,DIRECTORY,PREFIX
from .original_active_viscous_m15 import OriginalActiveViscousM15


PARENT=PREFIX+'compliant_current_original_Rm_inlet_m15_numeric.json'
PARENT_SHA='3f936cbb5f91fd05677f0c292ed783ad72788b0971ed869298f43be9dfc0850d'
SEAMS=('patch_reference','reference_slope','slope_axial','axial_buffer','buffer_transition','transition_power')


def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()).hexdigest()


class OriginalOuterChainM15(OriginalActiveViscousM15):
    def __init__(self,*,directory=DIRECTORY,precision=150):
        self.directory=directory;self.ctx=mp.mp.clone();self.ctx.dps=precision
        self.iv=MPIntervalContext();self.iv.dps=precision+30;self.hashes={};self.code_hashes={}
        raw=(directory/PARENT).read_bytes()
        if hashlib.sha256(raw).hexdigest()!=PARENT_SHA:raise ValueError('Accepted Rm inlet state required')
        self.parent=json.loads(raw);self.family=self.parent['source_family']
        if not self.parent['checks']['passed']:raise ValueError('Checked parent required')
        self.merge(self.parent['input_hashes'],directory,self.hashes)
        self.merge(self.parent['implementation_hashes'],ROOT,self.code_hashes);self.hashes[PARENT]=PARENT_SHA
        self.joined=self.archive('Rm_outer_mixed_source_joins')
        self.j=self.joined['Rm_outer_mixed_source_join_bindings'];self.receipt=self.joined['_receipt']
        self.nodes=self.joined['exact_graph_nodes']
        self.leading=self.archive('true_Rm_patch_C6_functions')
        self.outer=self.archive('true_outer_C6_C5_functions')
        for source in (self.leading,self.outer):
            if self.nodes[:len(source['exact_graph_nodes'])]!=source['exact_graph_nodes']:
                raise ValueError('Actual source archive prefixes must coincide')
        sys.setrecursionlimit(max(sys.getrecursionlimit(),20000))
        from .original_outer_chain_source import build_outer_chain_source
        self.source=build_outer_chain_source(self.nodes,self.joined,self.receipt,root=ROOT,
            hashes=self.hashes,code_hashes=self.code_hashes)
        if not self.source['passed']:raise ValueError('Six actual velocity/pressure source interfaces required')
        from .original_outer_chain_native import build_outer_chain_native
        self.native=build_outer_chain_native(self.nodes,self.joined,root=ROOT)
        if not self.native['passed']:raise ValueError('Actual native Jacobians and phase calculus required')
        from .original_outer_chain_mean import build_outer_chain_mean
        self.mean=build_outer_chain_mean(self.nodes,self.joined,self.receipt,self.leading,self.outer,root=ROOT)
        if not self.mean['passed']:raise ValueError('Actual nonzero mean and original leading germs required')
        self.seam_proofs=self.recover_traces()
        from .original_outer_power_transport import build_outer_power_transport
        power_guard=build_outer_power_transport(self.nodes,self.joined,self.outer,self.leading,root=ROOT)
        self.power_transport=power_guard['power_transport_certificate']
        for row in self.power_transport['source_pins'].values():
            if not row['passed']:raise ValueError('Actual power transport source pin changed')
            path=row.get('path') or row.get('relative')
            if path and path.startswith('experiments/root_st073/'):
                name=path.rsplit('/',1)[-1];sha=row['sha256']
                if name in self.hashes and self.hashes[name]!=sha:raise ValueError('Power source pin disagrees with parent')
                self.hashes[name]=sha
        if self.power_transport.get('passed'):
            self.seam_proofs['power_transport_certificate']=self.power_transport
        for name in ('original_outer_chain_m15.py','original_outer_chain_source.py','original_outer_chain_native.py',
                     'original_outer_chain_mean.py','original_outer_chain_numeric_m15.py','original_outer_power_transport.py'):
            path='src/openai_ns_reconstruction/'+name
            self.code_hashes[path]=hashlib.sha256((ROOT/path).read_bytes()).hexdigest()

    def recover_traces(self):
        """Apply proved function equalities only at their actual endpoint scopes.

        This is theorem application, not a raw graph fingerprint or an
        assignment of numerical values. RQy and D are explicitly local
        recipes whose arguments are archived function/endpoint nodes.
        """
        charts=self.j['charts'];joins=self.j['joins'];nodes=self.nodes
        leading=self.mean['original_leading_Qy_input_rows']
        if not (leading['passed'] and leading['row_count']==96 and leading['radial_row_count']==24):
            raise ValueError('Independent original leading source inputs required')
        source_path='experiments/root_st073/'+PREFIX+'compliant_current_original_modified_Q4_pressure_C3_functions.py'
        source_raw=(ROOT/source_path).read_bytes();source_sha=hashlib.sha256(source_raw).hexdigest()
        if self.hashes.get(source_path.rsplit('/',1)[-1])!=source_sha:
            raise ValueError('Pinned original physical Qy recovery source required')
        module=ast.parse(source_raw.decode('utf-8'))
        function=next(r for r in module.body if isinstance(r,ast.FunctionDef) and r.name=='physical_pressure_fields')
        source_text=ast.get_source_segment(source_raw.decode('utf-8'),function)
        if not all(s in source_text for s in ('my=a.sub(U,m)','Q=div(numerator(U,m));Qy=div(numerator(Uy,my))',
                'b.mul(d,b.shift(mean))','m=a.scale(fields[\'m\'],Pstar)')):
            raise ValueError('Actual physical mean ODE and Qy formula changed')
        qy_formula_guard=dict(passed=True,source=source_path,sha256=source_sha,
            function_AST_sha256=hashlib.sha256(ast.dump(function,include_attributes=False).encode()).hexdigest(),
            modified_m_y_equals_U_minus_m=True,
            Qy='[2Z Uy-(1-delta)Z(U-m)-(1-Z^2)(Uz-mz)]/(1-delta Z^2)',
            leading_substitution='U=Pstar*V; physical_mean=Pstar*m; Uy=Pstar*V_y',
            ordinary_Z1_uses_U_and_m_through_Z2=True,source_executed=False)

        def wrapper(node,expression,chart,endpoint):
            expected=dict(operation='function_substitution',expression=expression,
                variable=chart['coordinate'],value=endpoint,Z_independent_substitution=True)
            if nodes[node]!=expected:raise ValueError('Foreign or unscoped endpoint wrapper '+str(node))
            return dict(node=node,expression=expression,variable=chart['coordinate'],value=endpoint,passed=True)

        # Radius graph has only exact arithmetic, native coordinates, and
        # explicitly Z-independent fixed parameters. This check does not
        # infer fixed parameter values from their enclosing intervals.
        seen=set();stack=[c['radius'] for c in charts.values()];coordinates={c['coordinate'] for c in charts.values()}
        while stack:
            node=stack.pop()
            if node in seen:continue
            seen.add(node);r=nodes[node];op=r['operation']
            if op in ('sum','product'):stack.extend(r['arguments'])
            elif op in ('negative','analytic_unary'):stack.append(r['argument'])
            elif op=='positive_quotient':stack.extend((r['numerator'],r['denominator']))
            elif op=='function_substitution':
                if not r['Z_independent_substitution']:raise ValueError('Moving Z radius substitution')
                stack.extend((r['expression'],r['value']))
            elif op=='exact_rational':pass
            elif op=='bound_variable':
                if node not in coordinates:raise ValueError('Unexpected radius free variable')
            elif op=='current_original_source_parameter':
                if not r.get('Z_independent'):raise ValueError('Z-dependent radius source parameter')
            else:raise ValueError('Unproved radius Z independence: '+op)
        radius_guard=dict(passed=True,checked_node_count=len(seen),ordinary_Z_derivatives_zero=True,
            native_coordinates_not_frozen=True,source_parameter_values_not_selected=True)

        result={}
        for seam in SEAMS:
            proof=self.source['seams'][seam];join=joins[seam]
            if not proof['passed'] or not self.native['seamgeometry'][seam]['passed']:
                raise ValueError('Actual seam source and native guards required')
            wrapper_rows=[];trace_rows=[]
            for quantity,key,y,count,source_key in (
                ('U','U',0,3,'U_rows'),('Uy','U',1,2,'Uy_rows'),
                ('m','m',0,3,None),('Q','Q',0,3,'Q_rows'),('Qy','Q',1,2,'Qy_rows')):
                evidence=self.mean['modified_physical_mean_rows'][seam]['rows'] if source_key is None else proof[source_key]
                if len(evidence)!=count:raise ValueError('Complete ordinary Z source germ required')
                for k,row in enumerate(evidence):
                    if row['ordinary_Z_order']!=k or not row['passed']:raise ValueError('Unproved source derivative')
                    ports={};functions={}
                    for side in ('left','right'):
                        chart=charts[join[side+'_chart']];endpoint=join[side]['native_query']
                        function_node=chart['physical'][key][y][k];endpoint_node=join[side]['physical'][key][y][k]
                        wrapper_rows.append(dict(side=side,quantity=quantity,ordinary_Z_order=k,
                            **wrapper(endpoint_node,function_node,chart,endpoint)))
                        expected=function_node if quantity in ('U','Uy') else endpoint_node
                        if row[side+'_node']!=expected:raise ValueError('Source evidence bound to foreign graph port')
                        ports[side]=endpoint_node;functions[side]=function_node
                    trace_rows.append(dict(quantity=quantity,ordinary_Z_order=k,passed=True,
                        left_node=ports['left'],right_node=ports['right'],function_nodes=functions,
                        common_function_sha256=digest([seam,quantity,k,digest(row)]),
                        source_receipt_sha256=self.source['common_response_receipt_sha256']))
            radius_ports={}
            for side in ('left','right'):
                chart=charts[join[side+'_chart']];radius_node=join[side]['radius'][0]
                wrapper_rows.append(dict(side=side,quantity='R',**wrapper(radius_node,chart['radius'],chart,join[side]['native_query'])))
                radius_ports[side]=radius_node
            for k in (0,1):
                qy=next(row for row in trace_rows if row['quantity']=='Qy' and row['ordinary_Z_order']==k)
                trace_rows.append(dict(quantity='RQy',ordinary_Z_order=k,passed=True,
                    left_recipe=dict(operation='local_product_recipe',arguments=[radius_ports['left'],qy['left_node']]),
                    right_recipe=dict(operation='local_product_recipe',arguments=[radius_ports['right'],qy['right_node']]),
                    recipe_is_archived_graph_node=False,physical_R_once=True,
                    common_radius_premise='accepted original_radius_and_phase; native source guard',
                    common_function_sha256=digest([seam,'R_Qy',k,digest(qy)])))
            background=self.mean['leading_Qy_inputs'][seam]
            if not background['passed']:raise ValueError('Independent original leading germ proof required')
            leading_rows=[]
            for group,key,count in (('input_rows','V',3),('input_rows','m',3),('radial_rows','V_y',2)):
                pair={side:next(row for row in background['side_rows'][side][group] if row['field']==key) for side in ('left','right')}
                for k in range(count):
                    leading_rows.append(dict(quantity=key,ordinary_Z_order=k,
                        left_node=pair['left']['nodes'][k],right_node=pair['right']['nodes'][k],
                        left_scope=proof['left_scope'],right_scope=proof['right_scope'],
                        passed=True,source_node_alias_issued=False))
            common=dict(passed=True,input_rows=leading_rows,input_row_count=8,
                actual_physical_IVP_rows_Z0_to_Z5=background['actual_physical_IVP_rows_Z0_to_Z5'],
                C6_ordinary_Z_differentiation=True,physical_Qy_formula=qy_formula_guard,
                common_Pstar_premise=self.native['scope']['Pstar_nodes'],
                original_leading_controls_not_zeroed=True,output_alias_assigned=False)
            d_rows=[dict(ordinary_Z_order=k,passed=True,physical_R_once=True,
                local_recipe='D=R*(Qy_modified-Qy_leading)',
                common_modified_Qy_ordinary_Z_order=k,common_leading_Qy_ordinary_Z_order=k,
                R_ordinary_Z_derivatives_zero=True,signed_endpoint_values_evaluated=False) for k in (0,1)]
            result[seam]=dict(passed=True,left_chart=proof['left_chart'],right_chart=proof['right_chart'],
                left_endpoint=proof['left_native_query'],right_endpoint=proof['right_native_query'],
                trace_rows=trace_rows,endpoint_wrapper_rows=wrapper_rows,endpoint_wrapper_binding_passed=True,
                leading_Qy_common=common,D_equality_rows=d_rows,
                endpoint_scope_restricted=True,source_values_selected=False,numeric_caps_used=False,
                original_graph_not_mutated=True,incoming_nonzero_mean_history_retained=True)
        return dict(passed=True,joins=result,source_values_selected=False,numeric_caps_used=False,
            physical_Qy_formula_guard=qy_formula_guard,radius_Z_independence_guard=radius_guard,
            trace_row_count=sum(len(row['trace_rows']) for row in result.values()),
            endpoint_wrapper_row_count=sum(len(row['endpoint_wrapper_rows']) for row in result.values()),
            theorem_application_not_raw_graph_fingerprint=True,ancestor_replayed=False)

    def evaluate(self,Z):
        if not hasattr(self,'numeric'):raise ValueError('Build the correlated view first')
        row=next((s for s in self.numeric['samples'] if s['Z']==Z),None)
        if row is None:raise ValueError('Only archived axial queries are admitted')
        return row

    def build(self):
        from .original_outer_chain_numeric_m15 import build_outer_chain_numeric
        self.numeric=build_outer_chain_numeric(self,self.seam_proofs)
        return dict(source_family=self.family,parent_report=PARENT,parent_report_sha256=PARENT_SHA,
            source_derivative_theorem=self.source,native_source_guards=self.native,
            mean_source_certificate=self.mean,source_seam_proofs=self.seam_proofs,
            power_transport_certificate=self.power_transport,
            numerical_source_binding=self.numeric['source_binding'],samples=self.numeric['samples'],
            contract=dict(accepted_ancestor_constructors_not_run=True,source_functions_not_cap_values=True,
                raw_213_and_independent_target_retained=True),
            limits=dict(original_signed_finite_N_endpoint_values=False,numeric_Rp_closure=False,
                global_all_five_common_N_closure=False,exact_heat_energy_stress_flat=False,
                true_temporal_n_recursion_summation=False,mean_oscillatory_correction=False,
                independent_corrected_Cartesian_NS=False,full_goal_complete=False))
