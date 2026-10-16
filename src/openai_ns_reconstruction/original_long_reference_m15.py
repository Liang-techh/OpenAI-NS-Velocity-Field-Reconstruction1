"""Actual long-reshape/reference physical trace and signed m15 continuation.

The source calculus is reused through immutable certificates. New premises
are proved at long native1/reference native0, not borrowed from the distant
reference1/restoration0 seam. Accepted ancestor constructors never run.
"""
import copy
import hashlib
import json
import sys

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from .original_core_numeric import ROOT,DIRECTORY,PREFIX
from .original_active_viscous_m15 import OriginalActiveViscousM15,ActiveEndpointFingerprint
from .original_native_derivative_guards import (
    _phase_binding,_radius_log_jacobian_binding,_slow_b_source_binding,
    _direct_factor,_contains,
    _native_expr,_poly_serial,
)

PARENT=PREFIX+'compliant_current_original_reference_endpoint_m15_numeric.json'
PARENT_SHA='c257265a46adcefba1d00e10c61779a7333d7da411820efd1abc2a2aade95553'
CHARTS=('long_reshape','reference')


def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def long_radius_guard(nodes,chart,N,origin):
    """Admit the actual R=110*exp(T*q) product, without a fictitious sum."""
    coordinate=chart['coordinate'];J=chart['physical_measure'];R=chart['radius']
    radius=nodes[R];exp_rows=[(i,nodes[i]) for i in radius.get('arguments',())
                             if nodes[i].get('operation')=='analytic_unary' and nodes[i].get('name')=='exp']
    if radius['operation']!='product' or len(exp_rows)!=1:raise ValueError('Actual long radius product required')
    exp_id,exponential=exp_rows[0];argument=exponential['argument'];memo={}
    jexpr=_native_expr(nodes,J,coordinate,memo);lexpr=_native_expr(nodes,argument,coordinate,memo)
    other=[_native_expr(nodes,i,coordinate,memo) for i in radius['arguments'] if i!=exp_id]
    offset=nodes[chart['global_offset']];terms=offset.get('arguments',())
    logs=[i for i in terms if nodes[i].get('operation')=='analytic_unary' and nodes[i].get('name')=='log']
    negatives=[i for i in terms if nodes[i].get('operation')=='negative' and nodes[i].get('argument')==origin]
    if offset['operation']!='sum' or len(terms)!=2 or len(logs)!=1 or len(negatives)!=1:
        raise ValueError('Exact original log ratio minus phase origin required')
    ratio=nodes[nodes[logs[0]]['argument']]
    if ratio['operation']!='positive_quotient' or ratio['numerator']!=R:raise ValueError('Actual long radius ratio required')
    ra=ratio['denominator'];raexpr=_native_expr(nodes,ra,coordinate,memo)
    oexpr=_native_expr(nodes,origin,coordinate,memo);nexpr=_native_expr(nodes,N,coordinate,memo)
    phase=_phase_binding(nodes,chart,N)
    flags=dict(log_R_native_derivative_equals_J=bool(jexpr is not None and lexpr is not None and lexpr[1]==jexpr[0]),
        native_J_is_constant=bool(jexpr is not None and not jexpr[1]),
        other_radius_factors_native_independent=bool(other and all(v is not None and not v[1] for v in other)),
        Ra_native_independent=bool(raexpr is not None and not raexpr[1]),
        phase_origin_native_independent=bool(oexpr is not None and not oexpr[1]),
        N_native_independent=bool(nexpr is not None and not nexpr[1]),phase_product=phase['passed'])
    passed=all(flags.values())
    return dict(passed=passed,method='actual_long_radius_product_log_derivative',radius=R,coordinate=coordinate,
        physical_measure=J,exp=exp_id,log_argument=argument,log_argument_operation=nodes[argument]['operation'],
        canonical_J_expression=_poly_serial(jexpr[0]) if jexpr else None,
        canonical_log_R_native_derivative=_poly_serial(lexpr[1]) if lexpr else None,
        global_offset=chart['global_offset'],Ra=ra,phase_origin=origin,checks=flags,
        physical_phase_factor_after_J_division_is_N=passed,
        source_identity='d_native(log R)=J; d_native(N*global_offset)/J=N' if passed else None,
        actual_graph_not_modified=True)


class LongReferenceFingerprint(ActiveEndpointFingerprint):
    """Authorize only the four source-proved slow derivatives in seam scope."""
    def __init__(self,nodes,l,m,theorem,response):
        if not theorem.get('passed') or theorem.get('response_definition_sha256')!=digest(response):
            raise ValueError('Completed actual source theorem and its response required')
        if not theorem.get('actual_root_Z1_rows_consumed')==18 or theorem.get('mixed_yZ_zero_rows_consumed')!=5:
            raise ValueError('Actual root and mixed derivative premises required')
        super().__init__(nodes,dict(coordinate=l['charts']['long_reshape']['coordinate'],loop=l['charts']))
        self.authorizations={};self.alias_hits=set()
        for row in theorem['rows']:
            order=row['ordinary_Z_order']
            if order not in (0,1) or not row['passed']:raise ValueError('Both ordinary source rows required')
            for side,name in zip(('left','right'),CHARTS):
                actual=m['charts'][name]['radial']['slow_B_y'][order]
                if row[side+'_node']!=actual:raise ValueError('Foreign source derivative port')
                chart=l['charts'][name];endpoint=self.one if name=='long_reshape' else self.zero
                self.authorizations[actual]=(chart['coordinate'],endpoint,
                    (self.digest(['proved_long_reference_physical_B_y',digest(response),order]),None))
        if len(self.authorizations)!=4:raise ValueError('Four exact derivative ports required')

    def at(self,node,scope=None):
        if node in self.authorizations:
            coordinate,endpoint,value=self.authorizations[node]
            if scope is not None and scope.get(coordinate)==endpoint:
                self.alias_hits.add(node);return value
        return super().at(node,scope)


class OriginalLongReferenceM15(OriginalActiveViscousM15):
    def __init__(self,*,directory=DIRECTORY,precision=150):
        self.directory=directory;self.ctx=mp.mp.clone();self.ctx.dps=precision
        self.iv=MPIntervalContext();self.iv.dps=precision+30;self.hashes={};self.code_hashes={}
        raw=(directory/PARENT).read_bytes()
        if hashlib.sha256(raw).hexdigest()!=PARENT_SHA:raise ValueError('Accepted endpoint state required')
        self.parent=json.loads(raw);self.family=self.parent['source_family']
        if not self.parent['checks']['passed']:raise ValueError('Checked parent required')
        self.merge(self.parent['input_hashes'],directory,self.hashes)
        self.merge(self.parent['implementation_hashes'],ROOT,self.code_hashes);self.hashes[PARENT]=PARENT_SHA
        self.native=self.saved('reference_native_m15');self.phase=self.saved('reference_phase_lift')
        self.long=self.saved('long_m15_numeric');self.reference=self.saved('reference_m15_numeric');self.restore=self.saved('restore_m15_numeric')
        self.modified=self.archive('modified_Q4_pressure_C3_functions');self.loop=self.archive('prescribed_loop_C5_functions')
        self.outer=self.archive('prescribed_R100_Rm_C6_functions');self.nodes=self.modified['exact_graph_nodes']
        self.m=self.modified['modified_Q4_pressure_C3_bindings'];self.l=self.loop['prescribed_loop_C5_bindings']
        self.background=self.outer['prescribed_R100_Rm_C6_bindings']['charts']
        for source in (self.loop,self.outer):
            if self.nodes[:len(source['exact_graph_nodes'])]!=source['exact_graph_nodes']:
                raise ValueError('Identical original source prefixes required')
        sys.setrecursionlimit(max(sys.getrecursionlimit(),12000))
        self.native_guards=self.build_native_guards()
        self.raw_probe=self.prove_seam('long_reshape','reference',same_native_measure=False,strict=False)
        from .original_long_reference_germ import build_long_reference_germ
        self.germ=build_long_reference_germ(self.nodes,self.l,self.outer,root=ROOT)
        from .original_long_reference_roots import build_long_reference_roots
        self.roots=build_long_reference_roots(self.nodes,self.l,self.germ)
        if not self.germ['passed'] or not self.roots['passed']:raise ValueError('New seam physical germ/root premises incomplete')
        self.response=self.response_definition()
        self.theorem=self.build_theorem()
        self.recovery=self.recover_traces()
        self.seam_proof=dict(passed=self.recovery['passed'],left='long_reshape',right='reference',
            modified_and_leading_traces_both_equal=True,ordinary_Z_orders=[0,1],
            global_fast_phase_not_reset=True,
            physical_R_once=True,nonzero_external_histories_retained=True,
            actual_trace_recovery=self.recovery,source_derivative_theorem=self.theorem,
            leading_jet_certificate_sha256=digest(self.germ),no_caps_used=True)
        for name in ('original_long_reference_m15.py','original_long_reference_germ.py','original_long_reference_roots.py','original_long_reference_numeric_m15.py'):
            path='src/openai_ns_reconstruction/'+name
            self.code_hashes[path]=hashlib.sha256((ROOT/path).read_bytes()).hexdigest()

    def build_native_guards(self):
        inherited=self.native['native_derivative_source_guards']['checks'];ob=self.outer['prescribed_R100_Rm_C6_bindings']
        if not inherited['source_ast']['passed']:raise ValueError('Pinned chart-independent native calculus required')
        rows={};n=self.l['N'];nrow=self.nodes[n]
        if not nrow.get('exact_positive_integer') or not nrow.get('Z_independent'):raise ValueError('Common positive source N required')
        for name in CHARTS:
            chart=self.l['charts'][name];mc=self.m['charts'][name];bg=self.background[name]
            J=chart['physical_measure'];geometry=ob['geometry_certificate']
            positive=(J==ob['T'] and geometry['directed_T_positive']) if name=='long_reshape' else (J==ob['reference_length'] and geometry['directed_reference_gap_positive'])
            phase=_phase_binding(self.nodes,chart,n)
            radius=long_radius_guard(self.nodes,chart,n,self.l['phase_origin']) if name=='long_reshape' else _radius_log_jacobian_binding(self.nodes,chart,n,self.l['phase_origin'])
            slow=_slow_b_source_binding(self.nodes,mc['radial']['slow_B_y'],chart['loop']['B_over_Pstar'],J)
            uy=_direct_factor(self.nodes,mc['physical']['U_y'],ob['Pstar'])
            rqy=_direct_factor(self.nodes,mc['finite_boundary']['R_Q_y'],chart['radius'])
            seeds=all(mc['radial_background'][k]==bg['radial'][v] for k,v in (('E','E_y'),('V','V_y')))
            total=[]
            for k in (0,1):
                vy=mc['radial']['V_y'][k]
                total.append(dict(ordinary_Z_order=k,
                    slow_B_present=_contains(self.nodes,vy,mc['radial']['slow_B_y'][k]),
                    fast_Bphi_present=_contains(self.nodes,vy,chart['loop']['B_over_Pstar_phi'][k])))
            passed=bool(positive and phase['passed'] and radius['passed'] and slow['passed'] and uy['passed'] and rqy['passed'] and seeds and
                        all(t['slow_B_present'] and t['fast_Bphi_present'] for t in total))
            rows[name]=dict(passed=passed,native_Jacobian_node=J,positive_native_J=positive,
                positivity_source='directed_T_positive' if name=='long_reshape' else 'directed_reference_gap_positive',
                phase=phase,radius_log_J=radius,slow_B_source_binding=slow,
                physical_Uy_Pstar_factor=uy,physical_RQy_R_factor=rqy,
                original_radial_seed_binding=seeds,total_Vy_slow_fast_rows=total,
                actual_slow_B_y_nodes=mc['radial']['slow_B_y'][:2],
                physical_formula='Vy_modified=Vy_leading+slow_B_y/N+B_phi')
        if not all(r['passed'] for r in rows.values()):
            failed={name:{key:(value.get('passed') if isinstance(value,dict) else value)
                         for key,value in row.items() if key not in ('total_Vy_slow_fast_rows','actual_slow_B_y_nodes')}
                    for name,row in rows.items() if not row['passed']}
            raise ValueError('Actual long/reference native source guards incomplete '+str(failed))
        return dict(passed=True,charts=rows,common_N_node=n,
            reused_chart_independent_source_calculus=inherited['source_ast'],
            common_N_Z_independent=True,source_only=True,ancestor_constructors_not_run=True,
            no_native_J_equality_claimed=True)

    def response_definition(self):
        result=copy.deepcopy(self.native['common_physical_response'])
        bg=self.background['reference'];chart=self.l['charts']['reference'];p2=self.germ['physical_germ_recipes']['root_p2']
        ports={}
        for name,node,order in (('E',bg['fields']['E'][0],0),('Ez',bg['fields']['E'][1],1),
                                ('p',chart['roots']['p2'][0],0),('pz',chart['roots']['p2'][1],1)):
            ports[name]=dict(kind='original_endpoint_function',original_node=node,ordinary_Z_order=order,
                native_coordinate=chart['coordinate'],native_endpoint=0)
        for name,node,order in (('py',p2['physical_y_Z0_Z1'][0],0),('pyz',p2['physical_y_Z0_Z1'][1],1)):
            ports[name]=dict(kind='local_physical_recipe_function',recipe_node=node,ordinary_Z_order=order,
                recipe_graph='new_long_reference_germ.physical_germ_recipes',not_an_original_graph_node=True)
        ports['q']=dict(kind='original_endpoint_function',original_node=chart['loop']['q'][0],
            ordinary_Z1_node=chart['loop']['q'][1],native_coordinate=chart['coordinate'],native_endpoint=0,
            ordinary_Z1_zero_requires_actual_source_certificate=True)
        ports['dstar']=dict(kind='original_source_parameter',original_node=self.l['dstar'])
        result['ports']=ports
        result['new_germ_certificate_sha256']=digest(self.germ)
        result['endpoint_scope']='long_reshape native1/reference native0; all ordinary Z functions'
        return result

    def build_theorem(self):
        periodic=self.phase['phase_lift_certificate'];tangents=self.native['root_tangent_proof']
        if not periodic['formal_identities']['passed'] or not periodic['source_ast']['passed']:
            raise ValueError('Universal smooth periodic loop proof required')
        if not tangents['passed'] or len(tangents['ordinary_Z1_rows'])!=5:raise ValueError('Universal mixed root tangent identities required')
        passed=bool(self.native_guards['passed'] and self.germ['passed'] and self.roots['passed'])
        rows=[dict(quantity='slow_B_y',ordinary_Z_order=k,left_node=self.m['charts']['long_reshape']['radial']['slow_B_y'][k],
            right_node=self.m['charts']['reference']['radial']['slow_B_y'][k],common_response='G' if k==0 else 'G_Z',passed=passed) for k in (0,1)]
        return dict(passed=passed,rows=rows,response_definition_sha256=digest(self.response),
            actual_root_Z1_rows_consumed=18,mixed_yZ_zero_rows_consumed=5,
            new_germ_premise_sha256=digest(self.germ),new_root_premise_sha256=digest(self.roots),
            universal_root_tangent_identity_sha256=digest(tangents),
            universal_periodic_identity_sha256=digest(periodic['formal_identities']),
            new_active_domain_from_actual_root_functions=True,
            proof_method='pinned_native_chain_rule_induction_with_new_common_physical_two_jets',
            proof=[
                'At long native1 the original flat sigma has zero first and second jets; the actual log swirl hence has physical derivatives1/10 and0. Reference has the same law. This is a seam germ, not a whole long-chart identity.',
                'Common leading fields and Duhamel ODEs give identical physical two-jets and ordinary Z jets. Exact actual root rows give a=4/5,b=t0=0,Delta=-6/5 and their Z1 derivatives zero. Universal tangent identities then give all five y and yZ root tangents zero.',
                'The original physical p2_y and p2_yZ recipes retain genuine Z2 field inputs. The common loop is E*H(p2) at arbitrary held phase; its slow physical derivative is G=E*H/10+E*H_p*p2_y and G_Z retains H_pp*p2_Z*p2_y and H_p*p2_yZ.',
                'The pinned NativeDerivative rules seed J times the physical jets and divide the generated slow derivative by that same J. Structural chain-rule induction gives the common G/G_Z on both actual charts; their native Jacobians are distinct.',
                'Universal Poisson periodic/implicit-inverse identities apply on the new actual positive root domain, including phase integer cuts. Total Vy adds B_phi without another1/N. The retained finite-N mean solves my=Vmodified-mmodified.'
            ],arbitrary_actual_fast_phase=True,no_native_endpoint_as_phase_endpoint=True,
            original_p2_y_not_zeroed=True,original_signed_numeric_values_evaluated=False,
            full_raw_Uy_graph_fingerprint_claimed=False,error_caps_used=False)

    def recover_traces(self):
        f=LongReferenceFingerprint(self.nodes,self.l,self.m,self.theorem,self.response)
        a,b=(self.l['charts'][n] for n in CHARTS);scope={a['coordinate']:f.one,b['coordinate']:f.zero}
        rows=[]
        for group,key,count in (('physical','U',3),('physical','m',3),('physical','m_y',3),
                                ('physical','U_y',2),('physical','Q',2),('physical','Q_y',2),('finite_boundary','R_Q_y',2)):
            for k in range(count):
                left,right=(self.m['charts'][n][group][key][k] for n in CHARTS)
                lh,rh=f.at(left,scope)[0],f.at(right,scope)[0]
                if lh!=rh:raise ArithmeticError('Unproved actual long/reference trace '+str((key,k,left,right)))
                rows.append(dict(quantity=key,ordinary_Z_order=k,left_node=left,right_node=right,common_function_sha256=lh))
        expected={r[side+'_node'] for r in self.theorem['rows'] for side in ('left','right')}
        if f.alias_hits!=expected:raise ValueError('All four actual source derivative ports must be consumed')
        return dict(passed=True,rows=rows,authorized_slow_derivative_ports=sorted(f.alias_hits),
            original_graph_not_mutated=True,endpoint_scope_restricted=True,ordinary_Z2_U_and_mean_retained=True,
            full_raw_Uy_fingerprint=False,source_theorem_not_numeric_accuracy=True)

    def evaluate(self,Z):
        if not hasattr(self,'numeric'):raise ValueError('Build the new correlated view first')
        sample=next((s for s in self.numeric['samples'] if s['Z']==Z),None)
        if sample is None:raise ValueError('Only saved axial queries admitted; no original source evaluation')
        return sample

    def build(self):
        from .original_long_reference_numeric_m15 import build_long_reference_numeric
        self.numeric=build_long_reference_numeric(self,self.seam_proof)
        return dict(source_family=self.family,parent_report=PARENT,parent_report_sha256=PARENT_SHA,
            raw_scope_probe=self.raw_probe,new_physical_germ=self.germ,new_actual_root_certificate=self.roots,
            native_source_guards=self.native_guards,common_physical_response=self.response,
            source_derivative_theorem=self.theorem,actual_trace_recovery=self.recovery,
            source_seam_proof=self.seam_proof,numerical_source_binding=self.numeric['source_binding'],
            samples=self.numeric['samples'],
            contract=dict(accepted_ancestor_constructors_not_run=True,source_functions_not_cap_values=True,
                raw_213_and_independent_target_retained=True),
            limits=dict(original_signed_finite_N_endpoint_values=False,numeric_Rp_closure=False,
                global_all_five_common_N_closure=False,exact_heat_energy_stress_flat=False,
                true_temporal_n_recursion_summation=False,mean_oscillatory_correction=False,
                independent_corrected_Cartesian_NS=False,full_goal_complete=False))
