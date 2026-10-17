"""Actual compact/quiet physical seams and the signed tail m15 endpoint.

Only saved source archives are hydrated. No accepted ancestor constructor,
producer or checker runs. Endpoint recipes never become invented graph IDs.
"""
import ast
import hashlib
import json
import sys

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from .original_core_numeric import ROOT,DIRECTORY,PREFIX
from .original_active_viscous_m15 import OriginalActiveViscousM15


PARENT=PREFIX+'compliant_current_original_outer_chain_m15_numeric.json'
PARENT_SHA='cf185c36affc4183984aff5ef9caf119c68ea90d0b8e95535b1f8698e2161e53'
SEAMS=('power_compact','compact_quiet')


def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()).hexdigest()


class OriginalTailM15(OriginalActiveViscousM15):
    def __init__(self,*,directory=DIRECTORY,precision=150):
        self.directory=directory;self.ctx=mp.mp.clone();self.ctx.dps=precision
        self.iv=MPIntervalContext();self.iv.dps=precision+30;self.hashes={};self.code_hashes={}
        raw=(directory/PARENT).read_bytes()
        if hashlib.sha256(raw).hexdigest()!=PARENT_SHA:raise ValueError('Accepted six-outer-seam state required')
        self.parent=json.loads(raw);self.family=self.parent['source_family']
        if not self.parent['checks']['passed']:raise ValueError('Checked parent required')
        self.merge(self.parent['input_hashes'],directory,self.hashes)
        self.merge(self.parent['implementation_hashes'],ROOT,self.code_hashes);self.hashes[PARENT]=PARENT_SHA
        self.outer=self.archive('true_outer_C6_C5_functions')
        self.compact=self.archive('true_compact_C5_functions')
        self.quiet=self.archive('true_quiet_Rp_C5_pressure_C3_functions')
        self.nodes=self.quiet['exact_graph_nodes']
        for archive in (self.outer,self.compact):
            if self.nodes[:len(archive['exact_graph_nodes'])]!=archive['exact_graph_nodes']:
                raise ValueError('Identical actual outer/compact/quiet graph prefixes required')
        self.o=self.outer['true_outer_C6_C5_bindings'];self.c=self.compact['true_compact_C5_bindings']
        self.q=self.quiet['true_quiet_Rp_C5_C3_bindings']
        if self.c['Rc']!=self.o['endpoints']['O3_power']['radius'] or self.q['R2']!=self.c['terminal']['radius']:
            raise ValueError('Actual incoming physical radii required at both seams')
        sys.setrecursionlimit(max(sys.getrecursionlimit(),20000))
        from .original_tail_source import build_tail_source
        self.source=build_tail_source(self.nodes,self.outer,self.compact,self.quiet,root=ROOT,
            hashes=self.hashes,code_hashes=self.code_hashes)
        if not self.source['passed']:raise ValueError('Both actual compact/quiet source identities required')
        from .original_tail_native import build_tail_native
        self.native=build_tail_native(self.nodes,self.outer,self.compact,self.quiet,root=ROOT)
        if not self.native['passed']:raise ValueError('Actual varying-J and physical pressure wiring required')
        from .original_tail_mass_identity import build_tail_mass_identity
        self.mass_identity=build_tail_mass_identity(self.nodes,self.compact,root=ROOT)
        if not self.mass_identity['passed']:raise ValueError('Actual same-system compact mass identities required')
        self.quiet_transport=self.native.get('quiet_transport_certificate',dict(passed=False))
        self.seam_proofs=self.recover_traces()
        self.terminal_zero=self.prove_terminal_D_zero()
        if self.quiet_transport.get('passed'):
            self.seam_proofs['quiet_transport_certificate']=self.quiet_transport
        if self.terminal_zero.get('passed'):
            self.seam_proofs['terminal_D_zero_certificate']=self.terminal_zero
        for name in ('original_tail_m15.py','original_tail_source.py','original_tail_native.py','original_tail_numeric_m15.py',
                     'original_tail_mass_identity.py'):
            path='src/openai_ns_reconstruction/'+name
            self.code_hashes[path]=hashlib.sha256((ROOT/path).read_bytes()).hexdigest()

    def recover_traces(self):
        """Recover pressure by applying its original formula to proved fields."""
        source_path='experiments/root_st073/'+PREFIX+'compliant_current_original_modified_Q4_pressure_C3_functions.py'
        raw=(ROOT/source_path).read_bytes();sha=hashlib.sha256(raw).hexdigest()
        if self.hashes.get(source_path.rsplit('/',1)[-1])!=sha:raise ValueError('Pinned physical pressure source required')
        text=raw.decode('utf-8');tree=ast.parse(text)
        fn=next(r for r in tree.body if isinstance(r,ast.FunctionDef) and r.name=='physical_pressure_fields')
        body=ast.get_source_segment(text,fn)
        if not all(s in body for s in ('my=a.sub(U,m)','Q=div(numerator(U,m));Qy=div(numerator(Uy,my))',
                'b.mul(d,b.shift(mean))','m=a.scale(fields[\'m\'],Pstar)')):
            raise ValueError('Actual physical pressure formula changed')
        formula=dict(passed=True,source=source_path,sha256=sha,
            function_AST_sha256=hashlib.sha256(ast.dump(fn,include_attributes=False).encode()).hexdigest(),
            Q='[2Z U-(1-delta)Z m-(1-Z^2)mz]/(1-delta Z^2)',
            Qy='[2Z Uy-(1-delta)Z(U-m)-(1-Z^2)(Uz-mz)]/(1-delta Z^2)',
            actual_m_y_equals_U_minus_m=True,ordinary_Q_Z2_requires_m_Z3=True,
            ordinary_Qy_Z1_requires_U_Z2_m_Z2_Uy_Z1=True,source_executed=False)
        outer_chart=self.o['charts']['O3_power'];background=self.c['background']
        specs=dict(
            O3_power=dict(chart='O3_power',coordinate=outer_chart['coordinate'],point=48,value=2,
                radius=outer_chart['radius'],physical=outer_chart['physical'],boundary=outer_chart['boundary'],
                fields=outer_chart['modified_fields'],radial=outer_chart['radial'],endpoint=self.o['endpoints']['O3_power']),
            compact_inlet=dict(chart='compact',coordinate=background['phase'],point=1,value=1,
                radius=background['radius'],physical=self.c['physical'],boundary=self.c['boundary'],
                fields=self.c['modified_fields'],radial=self.c['radial'],endpoint=self.c['inlet']),
            compact_terminal=dict(chart='compact',coordinate=background['phase'],point=48,value=2,
                radius=background['radius'],physical=self.c['physical'],boundary=self.c['boundary'],
                fields=self.c['modified_fields'],radial=self.c['radial'],endpoint=self.c['terminal']),
            quiet=dict(chart='quiet',coordinate=self.q['s'],point=0,value=0,
                radius=self.q['R'],physical=self.q['physical'],boundary=self.q['boundary'],
                fields=self.q['fields'],radial=self.q['radial'],endpoint=self.q['inlet']))
        self.tail_specs=specs
        pstar=self.c['Pstar']
        if pstar!=self.q['Pstar'] or self.c['delta']!=self.q['delta']:raise ValueError('Same original physical parameters required')
        factor_rows=[]
        pressure_denominators=set()
        for name,spec in specs.items():
            for key in ('Q','Q_y'):
                row=self.nodes[spec['physical'][key][0]]
                if row.get('operation')!='positive_quotient' or row.get('source_positive_certificate')!='same original L=1-delta*Z^2':
                    raise ValueError('Actual positive original meridional denominator required')
                pressure_denominators.add(row['denominator'])
            for quantity,inputs,count in (('U',spec['fields']['V'],3),('U_y',spec['radial']['V_y'],2),('m',spec['fields']['m'],4)):
                for k in range(count):
                    node=spec['physical'][quantity][k];value=inputs[k];row=self.nodes[node]
                    zero=lambda i:self.nodes[i].get('operation')=='exact_rational' and self.nodes[i].get('numerator')==0
                    valid=(zero(node) and zero(value)) or (row.get('operation')=='product' and sorted(row['arguments'])==sorted([pstar,value]))
                    if not valid:raise ValueError('Physical field factor not bound to actual source '+str((name,quantity,k)))
                    factor_rows.append(dict(chart=name,quantity=quantity,ordinary_Z_order=k,physical_node=node,
                        actual_input_node=value,Pstar_node=pstar,passed=True))
            if not self.native['pressure_ports'][spec['chart']]['modified_mass_ODE']:
                raise ValueError('Actual physical mass ODE required')
        if len(pressure_denominators)!=1:raise ValueError('Same actual pressure denominator across tail charts required')
        formula['actual_common_positive_L_node']=next(iter(pressure_denominators))

        # Prove the physical radius and shared parameters do not depend on Z.
        # Native coordinates remain variable; inner substitutions bind their
        # own actual variables rather than treating them as global constants.
        seen=set();allowed={s['coordinate'] for s in specs.values()}
        def fixed(i,bound):
            key=(i,tuple(sorted(bound)))
            if key in seen:return
            seen.add(key);row=self.nodes[i];op=row['operation']
            if op in ('sum','product'):
                for child in row['arguments']:fixed(child,bound)
            elif op in ('negative','analytic_unary'):fixed(row['argument'],bound)
            elif op=='positive_quotient':fixed(row['numerator'],bound);fixed(row['denominator'],bound)
            elif op=='function_substitution':
                if not row['Z_independent_substitution'] or row['variable']==20:raise ValueError('Non-native radius substitution')
                fixed(row['value'],bound);fixed(row['expression'],bound|{row['variable']})
            elif op=='exact_rational':pass
            elif op=='bound_variable':
                if i not in bound:raise ValueError('Z/free-variable radius dependence')
            elif op in ('current_original_source_parameter','original_source_parameter','exact_positive_integer_expression'):
                if not row.get('Z_independent'):raise ValueError('Z-dependent physical parameter')
            else:raise ValueError('Unproved parameter/radius operation '+op)
        for spec in specs.values():fixed(spec['radius'],allowed)
        fixed(pstar,set());fixed(self.c['delta'],set())
        radius_guard=dict(passed=True,checked_scoped_nodes=len(seen),R_Z_derivatives_zero=True,
            Pstar_and_delta_Z_independent=True,parameter_values_not_selected=True)

        def endpoint_port(spec,quantity,k):
            expression=spec['physical'][quantity][k]
            recipe=dict(operation='local_endpoint_recipe',expression=expression,
                variable=spec['coordinate'],value=spec['point'],Z_independent_substitution=True)
            actual=spec['endpoint'].get('physical',{}).get(quantity)
            if actual is None:return dict(recipe=recipe,recipe_is_archived_node=False,passed=True)
            node=actual[k];expected=dict(operation='function_substitution',expression=expression,
                variable=spec['coordinate'],value=spec['point'],Z_independent_substitution=True)
            if self.nodes[node]!=expected and not (node==expression and self.nodes[node]['operation']=='exact_rational'):
                raise ValueError('Foreign actual physical endpoint wrapper')
            return dict(node=node,recipe=recipe,recipe_is_archived_node=False,actual_wrapper_bound=True,passed=True)

        bg=self.o['background_endpoints']['O3_power']
        leading=dict(
            O3_power=dict(V=bg['fields']['V'],m=bg['fields']['m'],V_y=bg['radial']['V_y']),
            compact_inlet=dict(V=self.c['background']['fields']['V'],m=self.c['background']['fields']['m'],
                V_y=self.c['background']['radial']['V_y']),
            compact_terminal=dict(V=self.c['background']['fields']['V'],m=self.c['background']['fields']['m'],
                V_y=self.c['background']['radial']['V_y']),
            quiet=dict(V=self.q['leading_flow']['profiles_y0_to_y4'][0]['V'],m=self.q['inlet']['leading_fields']['m'],
                V_y=self.q['leading_flow']['profiles_y0_to_y4'][1]['V']))
        joins={}
        for seam,left_name,right_name in (('power_compact','O3_power','compact_inlet'),('compact_quiet','compact_terminal','quiet')):
            source=self.source['seams'][seam];left=specs[left_name];right=specs[right_name]
            if not source['passed']:raise ValueError('Actual endpoint field theorem required')
            for side,spec in (('left',left),('right',right)):
                scope=source[side+'_scope']
                if scope['coordinate_node']!=spec['coordinate'] or scope['native_query_node']!=spec['point']:
                    raise ValueError('Source theorem is scoped to a different endpoint')
            traces=[];bindings=[]
            for quantity,physical_key,count,source_key in (('U','U',3,'modified_U_rows'),('Uy','U_y',2,'modified_Uy_rows'),('m','m',4,'modified_m_rows')):
                evidence=source[source_key]
                if len(evidence)!=count:raise ValueError('Actual source derivative rows missing')
                for k,row in enumerate(evidence):
                    if not row['passed'] or row['ordinary_Z_order']!=k:raise ValueError('Unproved source derivative')
                    ports={}
                    for side,spec in (('left',left),('right',right)):
                        port=endpoint_port(spec,physical_key,k);ports[side]=port
                        expected=port.get('node',spec['physical'][physical_key][k])
                        if row[side+'_node']!=expected:raise ValueError('Source field endpoint port mismatch')
                        bindings.append(dict(side=side,quantity=quantity,ordinary_Z_order=k,**port))
                    # m_Z3 is retained as extra input to Q_Z2, rather than
                    # inflating the three mean trace rows used for Qy_Z1.
                    if quantity!='m' or k<3:
                        traces.append(dict(quantity=quantity,ordinary_Z_order=k,left=ports['left'],right=ports['right'],
                            common_function_sha256=digest([seam,quantity,k,digest(row)]),passed=True))
            lead=source['leading_Qy_inputs']
            if not lead['passed'] or len(lead['rows'])!=8:raise ValueError('Independent original leading inputs required')
            for row in lead['rows']:
                k=row['ordinary_Z_order'];quantity=row['quantity']
                if not row['passed'] or row['left_node']!=leading[left_name][quantity][k] or row['right_node']!=leading[right_name][quantity][k]:
                    raise ValueError('Original leading port namespace mismatch '+str((seam,quantity,k)))
            for quantity,key,count in (('Q','Q',3),('Qy','Q_y',2)):
                for k in range(count):
                    a=endpoint_port(left,key,k);b=endpoint_port(right,key,k)
                    traces.append(dict(quantity=quantity,ordinary_Z_order=k,left=a,right=b,
                        physical_formula=formula,source_field_identity_sha256=digest(source),
                        theorem_application_not_raw_fingerprint=True,passed=True))
                    bindings.extend([dict(side='left',quantity=quantity,ordinary_Z_order=k,**a),dict(side='right',quantity=quantity,ordinary_Z_order=k,**b)])
            for k in (0,1):
                ports={}
                for side,spec in (('left',left),('right',right)):
                    node=spec['endpoint']['boundary']['R_Q_y'][k];expression=spec['boundary']['R_Q_y'][k]
                    expected=dict(operation='function_substitution',expression=expression,variable=spec['coordinate'],
                        value=spec['point'],Z_independent_substitution=True)
                    if self.nodes[node]!=expected:raise ValueError('Foreign actual RQy endpoint product wrapper')
                    ports[side]=dict(node=node,actual_full_chart_product_node=expression,
                        scope=dict(coordinate_node=spec['coordinate'],query_node=spec['point']),physical_R_once=True)
                traces.append(dict(quantity='RQy',ordinary_Z_order=k,left=ports['left'],right=ports['right'],
                    common_radius_and_Qy_source_deduction=True,passed=True))
            common=dict(passed=True,input_rows=lead['rows'],input_row_count=8,original_leading_namespace_bound=True,
                physical_formula=formula,source_values_selected=False,leading_output_alias_assigned=False)
            joins[seam]=dict(passed=True,left_chart=left['chart'],right_chart=right['chart'],
                left_endpoint=left['value'],right_endpoint=right['value'],left_endpoint_node=left['point'],right_endpoint_node=right['point'],
                left_coordinate=left['coordinate'],right_coordinate=right['coordinate'],
                left_scope=source['left_scope'],right_scope=source['right_scope'],
                trace_rows=traces,leading_Qy_common=common,endpoint_binding_rows=bindings,
                endpoint_wrapper_binding_passed=True,compact_missing_physical_wrappers_are_local_recipes=True,
                extra_mean_Z3_input_retained=True,
                D_equality_rows=[dict(ordinary_Z_order=k,passed=True,physical_R_once=True,
                    recipe='R*(Qy_modified-Qy_original_leading)',signed_endpoint_value_evaluated=False) for k in (0,1)],
                incoming_history_retained=True,source_values_selected=False,numeric_caps_used=False)
        return dict(passed=True,joins=joins,physical_pressure_formula=formula,physical_factor_rows=factor_rows,
            radius_parameter_guard=radius_guard,trace_row_count=sum(len(r['trace_rows']) for r in joins.values()),
            source_values_selected=False,numeric_caps_used=False,ancestor_replayed=False,
            original_graph_not_mutated=True,endpoint_recipes_not_registered_as_graph_nodes=True)

    def prove_terminal_D_zero(self):
        source=self.source['terminal_D_zero_inputs'];mass=self.source['compact_terminal_mass_zero']
        if not (source['passed'] and mass['passed'] and self.mass_identity['passed']):
            raise ValueError('Same-system terminal mass source identity required')
        if not (source['relative_m_exact_source_zero'] and source['absolute_modified_m_equals_leading_m']):
            raise ValueError('Relative compact mass zero is not an absolute-history reset')
        if not self.quiet_transport.get('passed'):raise ValueError('Independent actual quiet transport required')
        return dict(passed=True,source_only=True,D_upper_zero=True,ordinary_Z1=True,physical_R_once=True,
            endpoint='compact upper2Rc; identical quiet lower2Rc; propagated through quiet toRp',
            endpoint_node=48,compact_coordinate=self.c['background']['phase'],quiet_coordinate=self.q['s'],
            actual_relative_mass_certificate_node=mass['certificate_node'],
            actual_relative_mass_rows=mass['certificate']['actual_endpoint_rows'][:3],
            compact_modified_mass_endpoint_nodes=self.c['terminal']['modified_fields']['m'][:3],
            compact_leading_mass_endpoint_nodes=self.c['terminal']['background_fields']['m'][:3],
            quiet_modified_physical_mass_endpoint_nodes=self.q['inlet']['physical']['m'][:3],
            quiet_leading_mass_endpoint_nodes=self.q['inlet']['leading_fields']['m'][:3],
            quiet_actual_zero_U_Uy_endpoint_nodes={key:self.q['inlet']['physical'][key][:3 if key=='U' else 2] for key in ('U','U_y')},
            physical_mean_scale_node=self.c['Pstar'],
            compact_quiet_D_identity_sha256=digest(self.seam_proofs['joins']['compact_quiet']),
            independent_mass_row_algebra=self.mass_identity,
            relative_mass_source_identity_sha256=digest(mass),
            velocity_source_identity_sha256=digest(source),
            formula=self.seam_proofs['physical_pressure_formula'],
            deduction='DeltaU=DeltaUy=0 and Delta_m=Delta_mZ=Delta_mZZ=0 imply DeltaQy=DeltaQyZ=0; D=R*DeltaQy; quiet D_y=0',
            absolute_mean_not_zeroed=True,absolute_pressure_history_not_zeroed=True,
            quiet_mean_source='physical m=K*Z/R, actual K source retained',
            original_graph_not_mutated=True,source_values_selected=False,numeric_caps_used=False,ancestor_replayed=False)

    def evaluate(self,Z):
        if not hasattr(self,'numeric'):raise ValueError('Build the correlated view first')
        row=next((s for s in self.numeric['samples'] if s['Z']==Z),None)
        if row is None:raise ValueError('Only archived axial queries admitted')
        return row

    def build(self):
        from .original_tail_numeric_m15 import build_tail_numeric
        self.numeric=build_tail_numeric(self,self.seam_proofs)
        return dict(source_family=self.family,parent_report=PARENT,parent_report_sha256=PARENT_SHA,
            source_derivative_theorem=self.source,native_source_guards=self.native,
            compact_mass_row_identity=self.mass_identity,
            source_seam_proofs=self.seam_proofs,quiet_transport_certificate=self.quiet_transport,
            terminal_D_zero_certificate=self.terminal_zero,
            numerical_source_binding=self.numeric['source_binding'],samples=self.numeric['samples'],
            contract=dict(accepted_ancestor_constructors_not_run=True,source_functions_not_cap_values=True,
                raw_213_and_independent_target_retained=True),
            limits=dict(original_signed_finite_N_inlet_values=False,numeric_Rp_closure=False,
                global_all_five_common_N_closure=False,exact_heat_energy_stress_flat=False,
                true_temporal_n_recursion_summation=False,mean_oscillatory_correction=False,
                independent_corrected_Cartesian_NS=False,full_goal_complete=False))
