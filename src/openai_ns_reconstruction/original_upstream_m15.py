"""Actual upstream source seams and signed m15 inlet continuation.

Hydrate immutable saved graphs only. No accepted ancestor constructor,
producer or checker runs, and endpoint recipes are never archived graph IDs.
"""
import ast
import hashlib
import json
import sys
from fractions import Fraction

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from .original_core_numeric import ROOT,DIRECTORY,PREFIX
from .original_active_viscous_m15 import OriginalActiveViscousM15

PARENT=PREFIX+'compliant_current_original_tail_m15_numeric.json'
PARENT_SHA='a394d3b538e420936e382b36b63dc7405b7a098b57a557c668deba7a2cc60c9a'
SEAMS={'second_power':('second_switch','post_power'), 'power_long':('post_power','long_reshape')}


def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()).hexdigest()


def pressure_congruence(source,formula,*,loop_identity=None,derivative_identity=None):
    """Exact polynomial proof of ordinary pressure-jet congruence.

    The indeterminates are source functions, never numerical cap values.
    L_inverse is differentiated using L=1-delta*Z^2>0.
    """
    def add(*polys):
        out={}
        for poly in polys:
            for key,value in poly.items():out[key]=out.get(key,Fraction(0))+value
        return {key:value for key,value in out.items() if value}
    def scale(poly,c):return {key:value*Fraction(c) for key,value in poly.items() if value*c}
    def mul(*polys):
        out={():Fraction(1)}
        for poly in polys:
            new={}
            for left,a in out.items():
                for right,b in poly.items():
                    key=tuple(sorted(left+right));new[key]=new.get(key,Fraction(0))+a*b
            out={key:value for key,value in new.items() if value}
        return out
    def var(name):return {(name,):Fraction(1)}
    one={():Fraction(1)};z=var('Z');de=var('delta');linv=var('L_inverse')
    def derivative(poly):
        out={}
        for key,value in poly.items():
            for j,name in enumerate(key):
                if name=='Z':d=one
                elif name=='delta':continue
                elif name=='L_inverse':d=scale(mul(de,z,linv,linv),2)
                else:
                    family,order=name.rsplit('_',1);d=var(family+'_'+str(int(order)+1))
                out=add(out,scale(mul({key[:j]+key[j+1:]:Fraction(1)},d),value))
        return out
    def expression(side,kind):
        U=var(side+'U_0');Uy=var(side+'Uy_0');m=var(side+'m_0')
        if kind=='Q':u,mean,mz=U,m,var(side+'m_1')
        else:u,mean,mz=Uy,add(U,scale(m,-1)),add(var(side+'U_1'),scale(var(side+'m_1'),-1))
        numerator=add(scale(mul(z,u),2),scale(mul(add(one,scale(de,-1)),z,mean),-1),
            scale(mul(add(one,scale(mul(z,z),-1)),mz),-1))
        return mul(numerator,linv)
    def common(poly):
        out={}
        for key,value in poly.items():
            normalized=tuple(sorted(name[1:] if name.startswith(('LU_','RU_','LUy_','RUy_','Lm_','Rm_')) else name for name in key))
            out[normalized]=out.get(normalized,Fraction(0))+value
        return {key:value for key,value in out.items() if value}
    def verified_identity(row):
        identity=row.get('actual_function_identity')
        if not isinstance(identity,dict) or identity.get('passed') is not True:
            raise ValueError('Completed actual source-function identity required before pressure identification')
        for side in ('left','right'):
            if identity.get(side+'_node')!=row.get(side+'_node') or identity.get(side+'_scope')!=source.get(side+'_scope'):
                raise ValueError('Pressure input identity has a foreign actual node pair or endpoint scope')
        if identity.get('exact_endpoint_fingerprints') is True:
            if not identity.get('left_endpoint_sha256') or identity.get('left_endpoint_sha256')!=identity.get('right_endpoint_sha256'):
                raise ValueError('Unequal actual endpoint fingerprints cannot be identified')
        else:
            if not isinstance(loop_identity,dict) or not loop_identity.get('passed') or not identity.get('actual_loop_operator_identity_proved') or identity.get('actual_loop_identity_sha256')!=digest(loop_identity):
                raise ValueError('Matching actual issued loop-operator identity required for a non-direct pressure input')
        if identity.get('quantity')=='Uy':
            if not isinstance(derivative_identity,dict) or not derivative_identity.get('passed') or identity.get('issued_derivative_wiring_sha256')!=digest(derivative_identity['charts']):
                raise ValueError('Matching actual native derivative wiring required for a pressure input')
        return digest(identity)
    identities={}
    for family,key,count in (('U','modified_U_rows',3),('Uy','modified_Uy_rows',2),('m','modified_m_rows',4)):
        rows=source[key]
        if len(rows)!=count or not all(row['passed'] for row in rows):raise ValueError('Proved physical source inputs required')
        for k,row in enumerate(rows):identities[family+'_'+str(k)]=verified_identity(row)
    output=[]
    for kind,count in (('Q',3),('Qy',2)):
        left=expression('L',kind);right=expression('R',kind)
        for k in range(count):
            delta=add(left,scale(right,-1));reduced=common(delta)
            if reduced:raise ValueError('Nonzero formal ordinary pressure difference')
            variables={name[1:] for key in delta for name in key if name.startswith(('LU_','RU_','LUy_','RUy_','Lm_','Rm_'))}
            if not variables<=set(identities):raise ValueError('Insufficient ordinary source derivative order')
            output.append(dict(quantity=kind,ordinary_Z_order=k,passed=True,
                actual_input_identity_sha256={name:identities[name] for name in sorted(variables)},
                unreduced_difference_monomials=len(delta),reduced_difference_monomials=0,
                ordinary_derivative_product_rule=True,positive_original_L_shared=True))
            left=derivative(left);right=derivative(right)
    leading_identities={}
    lead=source['leading_Qy_inputs']
    if not lead['passed']:raise ValueError('Independent leading identities required')
    for row in lead['rows']:
        quantity={'V':'U','m':'m','V_y':'Uy','Vy':'Uy'}[row['quantity']]
        leading_identities[quantity+'_'+str(row['ordinary_Z_order'])]=verified_identity(row)
    leading_output=[];left=expression('L','Qy');right=expression('R','Qy')
    for k in (0,1):
        delta=add(left,scale(right,-1));reduced=common(delta)
        variables={name[1:] for key in delta for name in key if name.startswith(('LU_','RU_','LUy_','RUy_','Lm_','Rm_'))}
        if reduced or not variables<=set(leading_identities):raise ValueError('Independent leading pressure congruence failed')
        leading_output.append(dict(quantity='leading_Qy',ordinary_Z_order=k,passed=True,
            actual_input_identity_sha256={name:leading_identities[name] for name in sorted(variables)},
            reduced_difference_monomials=0,common_axially_constant_Pstar=True))
        left=derivative(left);right=derivative(right)
    return dict(passed=True,method='exact rational polynomial ordinary-Z quotient differentiation',rows=output,
        independent_leading_rows=leading_output,
        formula_AST_sha256=formula['function_AST_sha256'],
        inverse_derivative='d_Z L_inverse=2*delta*Z*L_inverse^2',
        independent_leading_application='same linear pressure operator on separately proved V3/m3/Vy2 with common Pstar',
        source_values_selected=False,numeric_caps_used=False)


def pressure_identity_negative_controls():
    controls={}
    for mode in ('missing_identity','foreign_node_pair','unequal_fingerprints'):
        sample={key:[dict(passed=True,ordinary_Z_order=k,left_node=100+k,right_node=200+k) for k in range(n)]
            for key,n in (('modified_U_rows',3),('modified_Uy_rows',2),('modified_m_rows',4))}
        sample.update(left_scope=dict(coordinate_node=1,native_query_node=1),right_scope=dict(coordinate_node=2,native_query_node=0))
        if mode!='missing_identity':
            row=sample['modified_U_rows'][0]
            row['actual_function_identity']=dict(passed=True,left_node=999 if mode=='foreign_node_pair' else row['left_node'],
                right_node=row['right_node'],left_scope=sample['left_scope'],right_scope=sample['right_scope'],
                exact_endpoint_fingerprints=True,left_endpoint_sha256='a'*64,right_endpoint_sha256='b'*64)
        try:pressure_congruence(sample,dict(function_AST_sha256='synthetic_negative_control_only'))
        except ValueError:controls[mode]=True
        else:raise ValueError('Unproved pressure input identification accepted: '+mode)
    return dict(passed=all(controls.values()),control_count=3,controls=controls,source_values_evaluated=False)


class OriginalUpstreamM15(OriginalActiveViscousM15):
    def __init__(self,*,directory=DIRECTORY,precision=150):
        self.pressure_input_guards=pressure_identity_negative_controls()
        self.directory=directory;self.ctx=mp.mp.clone();self.ctx.dps=precision
        self.iv=MPIntervalContext();self.iv.dps=precision+30;self.hashes={};self.code_hashes={}
        raw=(directory/PARENT).read_bytes()
        if hashlib.sha256(raw).hexdigest()!=PARENT_SHA:raise ValueError('Accepted tail state required')
        self.parent=json.loads(raw);self.family=self.parent['source_family']
        if not self.parent['checks']['passed']:raise ValueError('Checked parent required')
        self.merge(self.parent['input_hashes'],directory,self.hashes)
        self.merge(self.parent['implementation_hashes'],ROOT,self.code_hashes);self.hashes[PARENT]=PARENT_SHA
        self.loop=self.archive('prescribed_loop_C5_functions')
        self.modified=self.archive('modified_Q4_pressure_C3_functions')
        self.leading=self.archive('prescribed_R100_Rm_C6_functions')
        self.switch=self.archive('switch_mixed_source_joins')
        self.annular=self.archive('annular_mixed_source_joins')
        self.bridge=self.archive('prescribed_bridge_C6_functions')
        self.nodes=self.annular['exact_graph_nodes']
        for archive in (self.loop,self.modified,self.leading,self.switch,self.bridge):
            if self.nodes[:len(archive['exact_graph_nodes'])]!=archive['exact_graph_nodes']:
                raise ValueError('Identical original upstream graph prefixes required')
        self.l=self.loop['prescribed_loop_C5_bindings'];self.m=self.modified['modified_Q4_pressure_C3_bindings']
        self.b=self.leading['prescribed_R100_Rm_C6_bindings']
        sys.setrecursionlimit(max(sys.getrecursionlimit(),20000))
        from .original_upstream_source import build_upstream_source
        self.source=build_upstream_source(self.nodes,self.loop,self.modified,self.leading,self.switch,self.annular,root=ROOT)
        if not self.source['passed']:raise ValueError('Both actual upstream source identities required')
        from .original_upstream_native import build_upstream_native
        self.native=build_upstream_native(self.nodes,self.loop,self.modified,self.leading,self.switch,root=ROOT)
        self.native_guards=self.native
        if not self.native['passed']:raise ValueError('Actual native and physical pressure wiring required')
        self.complete_source_congruence()
        self.seam_proofs=self.recover_traces()
        for name in ('original_upstream_m15.py','original_upstream_source.py','original_upstream_native.py','original_upstream_numeric_m15.py','original_upstream_derivative.py','original_upstream_loop_identity.py'):
            path='src/openai_ns_reconstruction/'+name
            self.code_hashes[path]=hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
        self.bind_new_pins(self.source)
        self.bind_new_pins(self.native)
        self.bind_new_pins(self.derivative)
        self.bind_new_pins(self.loop_identity)

    def complete_source_congruence(self):
        """Compose actual background germs with their issued loop operators.

        The saved high-source pressure aliases are never used.  Common
        physical background jets, the actual generic-root operator, the
        original unique periodic inverse and exact Q4 derivative wiring
        establish the remaining velocity identities.
        """
        premises={};operator_rows={}
        for name,(left,right) in SEAMS.items():
            source=self.source['seams'][name];join=source['source_join']
            family='annular' if name=='power_long' else 'switch'
            archive=self.annular if family=='annular' else self.switch
            checks=archive['_receipt']['checks'];key='power_reshape' if name=='power_long' else 'second_power'
            lp=checks['actual_loop_phase_transfer' if family=='annular' else 'actual_loop_and_phase_transfer']
            roots_key='actual_original_roots_inverse_phase_and_NJ_bindings' if family=='annular' else 'actual_original_root_inverse_phase_J_bindings'
            phase_seams=lp['actual_radius_global_phase_seams' if family=='annular' else 'actual_radius_and_global_logarithmic_phase_seams']
            smooth_key='source_background_y5_unique_inverse_periodic_lift_then_leading_U_y4_C5' if family=='annular' else 'background_y5_supplies_root_y4_and_original_unique_inverse_periodic_lift_U_y4_C5'
            if not lp['passed'] or not lp[smooth_key] or key not in phase_seams:
                raise ValueError('Actual selected physical-germ/periodic-inverse theorem required')
            recipes=checks['actual_background_recipes']
            recipe_key='actual_five_field_history_generic_recipes_replayed' if family=='annular' else 'actual_field_six_history_and_generic_recipes_replayed'
            if not recipes['passed'] or not all(chart in recipes[recipe_key] for chart in (left,right)):
                raise ValueError('Actual selected background field/history/generic recipes required')
            background=checks['flat_controls_and_background_jets' if family=='annular' else 'physical_controls_and_ODE_jets']
            background_key='actual_background_phi_U_and_six_histories_physical_y0_to5' if family=='annular' else 'background_phi_U_and_six_histories_y0_to5_with_arbitrary_nonzero_inlets'
            if not background['passed'] or not background[background_key]:
                raise ValueError('Actual physical background through fifth y order required')
            rows=[]
            for chart_name in (left,right):
                chart=self.l['charts'][chart_name]
                expected=[chart_name,chart['roots'],chart['actual_phase'],chart['physical_measure']]
                actual=[row for row in lp[roots_key] if row[0]==chart_name]
                if actual!=[expected]:raise ValueError('Saved root/phase theorem bound to foreign chart functions')
                leading=self.b['charts']['post_switch' if chart_name=='post_power' else chart_name]
                for field in ('E','V'):
                    if chart['roots'][field]!=leading['fields'][field][:6]:
                        raise ValueError('Actual generic-root input differs from independent leading source')
                if chart['roots']['p2']!=leading['generic']['p2'][:6]:
                    raise ValueError('Actual root p2 source changed')
                rows.append(dict(chart=chart_name,passed=True,actual_root_phase_J_row=expected,
                    actual_leading_E_V_p2_nodes={field:leading['fields'][field][:6] for field in ('E','V')},
                    leading_p2_nodes=leading['generic']['p2'][:6],
                    coordinate_node=chart['coordinate'],radius_node=chart['radius'],
                    ordinary_Z_rows=list(range(6))))
            contract=join['source_contract']
            successor=contract['successor_incoming_history']
            if not contract['passed'] or not successor['passed']:
                raise ValueError('Actual nonzero predecessor/successor identity required')
            gauge=None
            if family=='annular':
                vg=checks['actual_value_and_gauge_joins']
                if not vg['passed'] or not vg['actual_F0_exp_minus_logC_minus_Lambda_G_read_from_source_not_assigned'] or not vg['actual_B_cancellation_and_exact_restoration_4Z_E_identity']:
                    raise ValueError('Actual R110 anchored F0/B gauge theorem required')
                values=[row for row in vg['actual_whole_Z_phi_U_and_six_history_value_rows'] if row[0]==key]
                if not any(row[:2]==[key,'phi'] for row in values):raise ValueError('Actual R110 phi germ row absent')
                bridge=self.bridge['prescribed_bridge_C6_bindings']['core7']
                if self.b['F0']!=bridge['F0'][:7]:raise ValueError('R110 F0 must be the actual anchored core source')
                endpoint=self.annular['annular_mixed_source_join_bindings']['background_endpoints'][key]
                wrapper_rows=[]
                for side,chart_name,point in (('left',left,1),('right',right,0)):
                    actual=self.b['charts']['post_switch' if chart_name=='post_power' else chart_name]['actual']
                    coordinate=self.l['charts'][chart_name]['coordinate']
                    for field in ('phi','U'):
                        for k in range(3):
                            node=endpoint[side]['actual'][field][k];expression=actual[field][k]
                            expected=dict(operation='function_substitution',expression=expression,variable=coordinate,
                                value=point,Z_independent_substitution=True)
                            if self.nodes[node]!=expected:raise ValueError('R110 actual gauge endpoint wrapper mismatch')
                            wrapper_rows.append(dict(side=side,field=field,ordinary_Z_order=k,node=node,
                                actual_full_function_node=expression,coordinate_node=coordinate,query_node=point,passed=True))
                gauge=dict(passed=True,source_join_key=key,actual_value_rows=values,
                    actual_endpoint_wrapper_rows=wrapper_rows,anchored_core_F0_nodes=bridge['F0'],
                    anchored_core_G_nodes=bridge['G'],anchored_core_logF0_nodes=bridge['logF0'],
                    F0_nodes=self.b['F0'],B_nodes=self.b['B'],
                    left_scope=source['left_scope'],right_scope=source['right_scope'],
                    anchor='post_power native1 / long_reshape native0 at R110',
                    long_native1_premise_used=False,source_values_selected=False)
            premise=dict(passed=True,left_scope=source['left_scope'],right_scope=source['right_scope'],
                source_contract=contract,successor_incoming_equalities=successor,
                background_physical_y0_to5=dict(passed=True,source_receipt_key=background_key,
                    physical_y_orders=list(range(6)),ordinary_Z_commutes_with_source_induction=True,
                    nonzero_original_histories_retained=True),
                unique_inverse_global_phase_source_germ=dict(passed=True,source_receipt_key=smooth_key,
                    actual_root_operator_rows=rows,actual_selected_radius_global_phase_seam=key,
                    periodic_lift_not_fractional_part_derivative=True),
                common_background_operator_inputs=dict(passed=True,
                    actual_recipe_receipt_key=recipe_key,actual_recipe_charts=recipes[recipe_key],
                    physical_background_orders=list(range(6)),source_germ_not_numeric_values=True,
                    actual_R110_gauge=gauge,shared_F0_nodes=self.b['F0'],Pstar_node=self.b['Pstar'],
                    source_recipe_guard='actual field/six-history/generic recipes replayed and equal to these archived ports'),
                source_values_selected=False,numeric_caps_used=False,ancestor_replayed=False)
            premises[name]=premise;operator_rows[name]=dict(passed=True,rows=rows,R110_gauge=gauge)
        from .original_upstream_derivative import build_upstream_derivative
        self.derivative=build_upstream_derivative(self.nodes,self.loop,self.modified,self.leading,{'seams':premises},root=ROOT)
        if not self.derivative['passed']:raise ValueError('Actual issued native derivative source congruence required '+str(self.derivative['checks']['required']))
        from .original_upstream_loop_identity import build_upstream_loop_identity
        derivative_sha=digest(self.derivative)
        for name,premise in premises.items():
            premise['issued_native_derivative']=dict(passed=True,certificate_sha256=derivative_sha,
                certificate_ref='native_derivative_certificate',source_seam=name,source_values_selected=False,numeric_caps_used=False)
        self.loop_identity=build_upstream_loop_identity(self.nodes,self.loop,self.leading,
            dict(seams=premises,native_derivative_certificate=self.derivative,native_derivative_certificate_sha256=derivative_sha),root=ROOT)
        if not self.loop_identity['passed']:raise ValueError('Actual all-root and q/B/Bphi endpoint identities required '+str(self.loop_identity.get('unproved')))
        for name,(left,right) in SEAMS.items():
            source=self.source['seams'][name];direct=source['direct_endpoint_fingerprints']['groups']
            if not all(direct[key]['passed'] for key in ('modified_m','leading_V','leading_m','leading_Vy')):
                raise ValueError('Independent actual mean and leading endpoint identities required')
            if not source['finite_function_source_rows']['passed']:
                raise ValueError('Actual finite-N V_N=V+B/N graph mapping required')
            actual_rows=[]
            for quantity,key,count in (('U','modified_U_rows',3),('Uy','modified_Uy_rows',2),('m','modified_m_rows',4)):
                for k,row in enumerate(source[key]):
                    if k>=count or not row['passed']:raise ValueError('Actual physical field source binding required')
                    field='U_y' if quantity=='Uy' else quantity
                    if row['left_node']!=self.m['charts'][left]['physical'][field][k] or row['right_node']!=self.m['charts'][right]['physical'][field][k]:
                        raise ValueError('Original velocity identity uses foreign actual output node')
                    direct_group=direct['modified_'+quantity]
                    direct_passed=direct_group['rows'][k]['passed']
                    if quantity=='m':
                        if not direct_passed:raise ValueError('Actual modified mean endpoint fingerprint required')
                        method='exact actual endpoint fingerprint; exp(0), empty integral and nonzero incoming memory'
                    elif direct_passed:method='exact actual endpoint fingerprint'
                    else:method='source-operator congruence on actual common physical germ and unique periodic inverse'
                    witness=dict(passed=True,quantity=quantity,ordinary_Z_order=k,
                        left_node=row['left_node'],right_node=row['right_node'],
                        left_scope=source['left_scope'],right_scope=source['right_scope'],method=method,
                        actual_root_operator_sha256=digest(operator_rows[name]),
                        germ_premise_sha256=digest(premises[name]),
                        issued_derivative_wiring_sha256=digest(self.derivative['charts']) if quantity=='Uy' else None,
                        finite_N_operator='Pstar*(V+B_over_Pstar/N)' if quantity=='U' else None,
                        physical_derivative_operator='Pstar*(leading_Vy+slow_B_y/N+B_over_Pstar_phi)' if quantity=='Uy' else None,
                        derivative_holds_fast_phase_fixed_then_adds_Bphi=quantity=='Uy',
                        ordinary_Z_product_chain_rule=True,physical_J_once=True,physical_Pstar_once=True,
                        no_pressure_output_alias=True,source_values_selected=False,numeric_caps_used=False)
                    witness['exact_endpoint_fingerprints']=direct_passed
                    if direct_passed:
                        witness['left_endpoint_sha256']=direct_group['rows'][k]['left_endpoint_sha256']
                        witness['right_endpoint_sha256']=direct_group['rows'][k]['right_endpoint_sha256']
                        if witness['left_endpoint_sha256']!=witness['right_endpoint_sha256']:
                            raise ValueError('Actual direct endpoint identity failed')
                    else:
                        if not self.loop_identity['passed']:raise ValueError('All actual root and loop identities required')
                        witness['actual_loop_operator_identity_proved']=True
                        witness['actual_loop_identity_sha256']=digest(self.loop_identity)
                    row['actual_function_identity']=witness;actual_rows.append(witness)
            for row in source['leading_Qy_inputs']['rows']:
                group={'V':'leading_V','m':'leading_m','V_y':'leading_Vy','Vy':'leading_Vy'}[row['quantity']]
                proof=direct[group]['rows'][row['ordinary_Z_order']]
                if not proof['passed'] or proof['left_node']!=row['left_node'] or proof['right_node']!=row['right_node'] or proof['left_endpoint_sha256']!=proof['right_endpoint_sha256']:
                    raise ValueError('Independent leading input must have its own exact actual endpoint identity')
                row['actual_function_identity']=dict(passed=True,exact_endpoint_fingerprints=True,
                    left_node=row['left_node'],right_node=row['right_node'],
                    left_scope=source['left_scope'],right_scope=source['right_scope'],
                    left_endpoint_sha256=proof['left_endpoint_sha256'],right_endpoint_sha256=proof['right_endpoint_sha256'],
                    namespace='independent original leading',source_values_selected=False,numeric_caps_used=False)
            source['actual_function_identity']=dict(passed=True,rows=actual_rows,
                independent_leading_endpoint_fingerprints={key:direct[key] for key in ('leading_V','leading_m','leading_Vy')},
                source_germ_premise=premises[name],actual_root_operator_binding=operator_rows[name],
                native_derivative_congruence_consumed_by_downstream_adapter=True,
                native_derivative_congruence_unproved_here=False,
                derivative_certificate_sha256=digest(self.derivative),
                source_values_selected=False,numeric_caps_used=False,ancestor_replayed=False)
            source['source_proof_premises']['actual_function_identity']=source['actual_function_identity']
            source['source_proof_premises']['issued_native_derivative_congruence_sha256']=digest(self.derivative)
        self.source['actual_germ_premises']=premises
        self.source['actual_root_operator_bindings']=operator_rows
        self.source['actual_loop_function_identity']=self.loop_identity
        self.source['native_derivative_congruence']=self.derivative
        self.source['actual_function_congruence_completed']=True

    def bind_new_pins(self,packet):
        for row in packet.get('source_pins',{}).values():
            if not isinstance(row,dict):continue
            if row.get('passed') is False:raise ValueError('Changed upstream source pin')
            name=row.get('path') or row.get('relative')
            if not name:continue
            actual=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
            if actual!=row['sha256']:raise ValueError('Upstream source pin mismatch '+name)
            if name.startswith('experiments/root_st073/'):
                key=name.rsplit('/',1)[-1];target=self.hashes
            else:key=name;target=self.code_hashes
            if key in target and target[key]!=actual:raise ValueError('Upstream pin conflicts with parent')
            target[key]=actual

    def recover_traces(self):
        """Recover physical Q/Qy from source-proved U, Uy and mean germs."""
        source_path='experiments/root_st073/'+PREFIX+'compliant_current_original_modified_Q4_pressure_C3_functions.py'
        raw=(ROOT/source_path).read_bytes();sha=hashlib.sha256(raw).hexdigest()
        if self.hashes.get(source_path.rsplit('/',1)[-1])!=sha:raise ValueError('Pinned physical pressure formula required')
        text=raw.decode('utf-8');module=ast.parse(text)
        function=next(r for r in module.body if isinstance(r,ast.FunctionDef) and r.name=='physical_pressure_fields')
        body=ast.get_source_segment(text,function)
        if not all(s in body for s in ('my=a.sub(U,m)','Q=div(numerator(U,m));Qy=div(numerator(Uy,my))',
                'b.mul(d,b.shift(mean))',"m=a.scale(fields['m'],Pstar)")):
            raise ValueError('Original pressure formula changed')
        formula=dict(passed=True,source=source_path,sha256=sha,
            function_AST_sha256=hashlib.sha256(ast.dump(function,include_attributes=False).encode()).hexdigest(),
            Q='[2Z U-(1-delta)Z m-(1-Z^2)mz]/(1-delta Z^2)',
            Qy='[2Z Uy-(1-delta)Z(U-m)-(1-Z^2)(Uz-mz)]/(1-delta Z^2)',
            actual_m_y_equals_U_minus_m=True,ordinary_Q_Z2_requires_m_Z3=True,
            ordinary_Qy_Z1_requires_U_Z2_m_Z2_Uy_Z1=True,source_executed=False)
        pstar=self.b['Pstar'];factor_rows=[];denominators=set();specs={}
        for name in ('second_switch','post_power','long_reshape'):
            loop=self.l['charts'][name];modified=self.m['charts'][name]
            background=self.b['charts']['post_switch' if name=='post_power' else name]
            if loop['coordinate']!=background['phase'] or loop['radius']!=background['radius']:
                raise ValueError('Independent original leading geometry mismatch')
            specs[name]=dict(chart=name,coordinate=loop['coordinate'],radius=loop['radius'],
                physical=modified['physical'],boundary=modified['finite_boundary'],
                fields=loop['modified_fields'],radial=modified['radial'],
                endpoint=self.m['endpoints'][name],
                leading=dict(V=background['fields']['V'],m=background['fields']['m'],Vy=background['radial']['V_y']))
            for quantity,inputs,count in (('U',loop['modified_fields']['V'],3),('U_y',modified['radial']['V_y'],2),('m',loop['modified_fields']['m'],4)):
                for k in range(count):
                    node=modified['physical'][quantity][k];value=inputs[k];r=self.nodes[node]
                    zero=lambda i:self.nodes[i].get('operation')=='exact_rational' and self.nodes[i].get('numerator')==0
                    valid=(zero(node) and zero(value)) or (r.get('operation')=='product' and sorted(r['arguments'])==sorted([pstar,value]))
                    if not valid:raise ValueError('Original physical factor mismatch '+str((name,quantity,k)))
                    factor_rows.append(dict(chart=name,quantity=quantity,ordinary_Z_order=k,physical_node=node,
                        actual_input_node=value,Pstar_node=pstar,passed=True))
            for key in ('Q','Q_y'):
                r=self.nodes[modified['physical'][key][0]]
                if r.get('operation')!='positive_quotient' or r.get('source_positive_certificate')!='same original L=1-delta*Z^2':
                    raise ValueError('Original positive pressure denominator required')
                denominators.add(r['denominator'])
            if not self.native['pressure_ports'][name]['m_y_equals_U_minus_m']:
                raise ValueError('Original physical mass ODE required')
        if len(denominators)!=1:raise ValueError('Same original pressure denominator required')
        formula['actual_common_positive_L_node']=next(iter(denominators))
        # Keep native variables and substitutions scoped while proving R and
        # Pstar independent of the ordinary axial variable Z.
        seen=set();coordinates={spec['coordinate'] for spec in specs.values()}
        def fixed(i,bound):
            key=i,tuple(sorted(bound))
            if key in seen:return
            seen.add(key);row=self.nodes[i];op=row['operation']
            if op in ('sum','product'):
                for child in row['arguments']:fixed(child,bound)
            elif op in ('negative','analytic_unary'):fixed(row['argument'],bound)
            elif op=='positive_quotient':fixed(row['numerator'],bound);fixed(row['denominator'],bound)
            elif op=='function_substitution':
                if not row['Z_independent_substitution'] or row['variable']==20:raise ValueError('Moving axial radius substitution')
                fixed(row['value'],bound);fixed(row['expression'],bound|{row['variable']})
            elif op=='exact_rational':pass
            elif op=='bound_variable':
                if i not in bound:raise ValueError('Free variable in physical radius or parameter')
            elif op in ('current_original_source_parameter','original_source_parameter','exact_positive_integer_expression'):
                if not row.get('Z_independent'):raise ValueError('Axially varying physical parameter')
            else:raise ValueError('Unproved radius operation '+op)
        for spec in specs.values():fixed(spec['radius'],coordinates)
        fixed(pstar,set())
        L=formula['actual_common_positive_L_node'];lr=self.nodes[L]
        negatives=[i for i in lr.get('arguments',[]) if self.nodes[i].get('operation')=='negative']
        if lr.get('operation')!='sum' or len(negatives)!=1 or sorted(lr['arguments'])!=sorted([1,negatives[0]]):
            raise ValueError('Actual L=1-delta*Z^2 source expression required')
        product=self.nodes[self.nodes[negatives[0]]['argument']]
        square=[i for i in product.get('arguments',[]) if self.nodes[i].get('operation')=='product' and self.nodes[i].get('arguments')==[20,20]]
        if product.get('operation')!='product' or len(product['arguments'])!=2 or len(square)!=1 or self.nodes[20]!={'operation':'bound_variable','name':'Z'}:
            raise ValueError('Actual ordinary Z pressure denominator source required')
        delta=next(i for i in product['arguments'] if i!=square[0]);fixed(delta,set())
        formula['actual_delta_node']=delta;formula['actual_Z_node']=20
        formula['delta_Z_independent']=True
        radius_guard=dict(passed=True,checked_scoped_nodes=len(seen),ordinary_R_Z_derivatives_zero=True,
            Pstar_Z_independent=True,delta_node=delta,delta_Z_independent=True,parameter_values_not_selected=True)

        def port(spec,key,k,point,*,boundary=False):
            expression=spec['boundary' if boundary else 'physical'][key][k]
            recipe=dict(operation='local_endpoint_recipe',expression=expression,
                variable=spec['coordinate'],value=point,Z_independent_substitution=True)
            out=dict(recipe=recipe,recipe_is_archived_node=False,full_chart_function_node=expression,passed=True)
            # Existing upper boundary wrappers are bound separately; physical
            # U/Uy/m/Q/Qy wrappers absent from these archives remain recipes.
            if boundary and point==1:
                node=spec['endpoint']['boundary'][key][k]
                expected=dict(operation='function_substitution',expression=expression,variable=spec['coordinate'],
                    value=point,Z_independent_substitution=True)
                if self.nodes[node]!=expected:raise ValueError('Foreign upper endpoint boundary wrapper')
                out.update(node=node,actual_wrapper_bound=True)
            return out

        def source_port(row,side,spec,key,k,point,*,leading=False):
            expression=(spec['leading'] if leading else spec['physical'])[key][k]
            actual=row[side+'_node']
            if actual==expression:return
            expected=dict(operation='function_substitution',expression=expression,variable=spec['coordinate'],
                value=point,Z_independent_substitution=True)
            if self.nodes[actual]!=expected:raise ValueError('Foreign source germ port '+str((side,key,k,actual,expression)))

        result={}
        for name,(left_name,right_name) in SEAMS.items():
            source=self.source['seams'][name];left=specs[left_name];right=specs[right_name]
            if not source['passed']:raise ValueError('Actual upstream field source proof required')
            for side,spec,point in (('left',left,1),('right',right,0)):
                scope=source[side+'_scope']
                if scope['coordinate_node']!=spec['coordinate'] or scope['native_query_node']!=point:
                    raise ValueError('Actual source theorem has foreign endpoint scope')
            traces=[];bindings=[];modified_rows=[]
            for quantity,key,count,source_key in (('U','U',3,'modified_U_rows'),('Uy','U_y',2,'modified_Uy_rows'),('m','m',4,'modified_m_rows')):
                rows=source[source_key]
                if len(rows)!=count:raise ValueError('Complete original modified source germ required')
                for k,row in enumerate(rows):
                    if not row['passed'] or row['ordinary_Z_order']!=k:raise ValueError('Unproved ordinary source derivative')
                    ports={}
                    for side,spec,point in (('left',left,1),('right',right,0)):
                        source_port(row,side,spec,key,k,point)
                        ports[side]=port(spec,key,k,point)
                        bindings.append(dict(side=side,quantity=quantity,ordinary_Z_order=k,**ports[side]))
                    modified_rows.append(dict(quantity=quantity,ordinary_Z_order=k,source_row_sha256=digest(row),
                        left_node=row['left_node'],right_node=row['right_node'],passed=True))
                    if quantity!='m' or k<3:
                        traces.append(dict(quantity=quantity,ordinary_Z_order=k,**ports,
                            common_function_sha256=digest([name,quantity,k,row]),passed=True))
            leading=source['leading_Qy_inputs']
            if not leading['passed'] or len(leading['rows'])!=8:raise ValueError('Independent leading source inputs required')
            leading_rows=[]
            for row in leading['rows']:
                if not row['passed']:raise ValueError('Unproved original leading derivative')
                key={'V_y':'Vy','Vy':'Vy'}.get(row['quantity'],row['quantity']);k=row['ordinary_Z_order']
                for side,spec,point in (('left',left,1),('right',right,0)):
                    source_port(row,side,spec,key,k,point,leading=True)
                leading_rows.append(dict(quantity=key,ordinary_Z_order=k,source_row_sha256=digest(row),
                    left_node=row['left_node'],right_node=row['right_node'],passed=True))
            for quantity,key,count in (('Q','Q',3),('Qy','Q_y',2)):
                for k in range(count):
                    a=port(left,key,k,1);b=port(right,key,k,0)
                    traces.append(dict(quantity=quantity,ordinary_Z_order=k,left=a,right=b,
                        physical_formula=formula,source_field_identity_sha256=digest(source),
                        theorem_application_not_raw_fingerprint=True,passed=True))
            for k in (0,1):
                ports={}
                for side,spec,point in (('left',left,1),('right',right,0)):
                    product=spec['boundary']['R_Q_y'][k];r=self.nodes[product]
                    if r.get('operation')!='product' or sorted(r['arguments'])!=sorted([spec['radius'],spec['physical']['Q_y'][k]]):
                        raise ValueError('Original radius factor must occur once')
                    ports[side]=port(spec,'R_Q_y',k,point,boundary=True)
                traces.append(dict(quantity='RQy',ordinary_Z_order=k,**ports,common_radius_and_Qy_source_deduction=True,passed=True))
            algebra=pressure_congruence(source,formula,loop_identity=self.loop_identity,derivative_identity=self.derivative)
            result[name]=dict(passed=True,left_chart=left_name,right_chart=right_name,left_endpoint=1,right_endpoint=0,
                left_endpoint_node=1,right_endpoint_node=0,left_coordinate=left['coordinate'],right_coordinate=right['coordinate'],
                left_scope=source['left_scope'],right_scope=source['right_scope'],
                trace_rows=traces,modified_source_rows=modified_rows,leading_source_rows=leading_rows,
                pressure_algebra_witness=algebra,
                leading_Qy_common=dict(passed=True,input_rows=leading['rows'],original_leading_namespace_bound=True,
                    physical_formula=formula,source_values_selected=False,leading_output_alias_assigned=False),
                endpoint_binding_rows=bindings,endpoint_wrapper_binding_passed=True,
                missing_physical_wrappers_are_local_recipes=True,extra_mean_Z3_input_retained=True,
                D_equality_rows=[dict(ordinary_Z_order=k,passed=True,physical_R_once=True,
                    recipe='R*(Qy_modified-Qy_original_leading)',signed_endpoint_value_evaluated=False) for k in (0,1)],
                incoming_history_retained=True,source_values_selected=False,numeric_caps_used=False)
        return dict(passed=True,joins=result,physical_pressure_formula=formula,physical_factor_rows=factor_rows,
            radius_parameter_guard=radius_guard,trace_row_count=sum(len(r['trace_rows']) for r in result.values()),
            source_values_selected=False,numeric_caps_used=False,ancestor_replayed=False,
            original_graph_not_mutated=True,endpoint_recipes_not_registered_as_graph_nodes=True)

    def evaluate(self,Z):
        if not hasattr(self,'numeric'):raise ValueError('Build the correlated view first')
        row=next((s for s in self.numeric['samples'] if s['Z']==Z),None)
        if row is None:raise ValueError('Only archived axial queries admitted')
        return row

    def build(self):
        from .original_upstream_numeric_m15 import build_upstream_numeric
        self.numeric=build_upstream_numeric(self,self.seam_proofs)
        return dict(source_family=self.family,parent_report=PARENT,parent_report_sha256=PARENT_SHA,
            source_derivative_theorem=self.source,native_source_guards=self.native,
            source_seam_proofs=self.seam_proofs,pressure_input_identity_guard_controls=self.pressure_input_guards,
            numerical_source_binding=self.numeric['source_binding'],
            samples=self.numeric['samples'],
            contract=dict(accepted_ancestor_constructors_not_run=True,source_functions_not_cap_values=True,
                raw_213_and_independent_target_retained=True),
            limits=dict(original_signed_finite_N_inlet_values=False,numeric_Rp_closure=False,
                global_all_five_common_N_closure=False,exact_heat_energy_stress_flat=False,
                true_temporal_n_recursion_summation=False,mean_oscillatory_correction=False,
                independent_corrected_Cartesian_NS=False,full_goal_complete=False))
