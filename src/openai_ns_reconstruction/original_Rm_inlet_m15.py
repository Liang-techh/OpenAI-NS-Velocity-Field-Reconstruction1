"""Source-bound postrestore/Rm inlet trace and signed m15 continuation."""
import hashlib
import json
import sys

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from .original_core_numeric import ROOT,DIRECTORY,PREFIX
from .original_active_viscous_m15 import OriginalActiveViscousM15
from .original_Rm_inlet_flat import RmInletFingerprint,build_Rm_inlet_flat_certificate


PARENT=PREFIX+'compliant_current_original_long_reference_m15_numeric.json'
PARENT_SHA='b7f3089fcd0257750dc8e4d082ad4fc59e2badd2360f35618eca7619bf46162b'


def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()


class RmTraceFingerprint(RmInletFingerprint):
    """Actual source-proved velocity/leading-mean ports, in inlet scope."""
    def __init__(self,nodes,l,m,rm,theorem,mean):
        if not theorem.get('passed'):raise ValueError('Completed actual Rm source theorem required')
        if not mean.get('passed'):raise ValueError('Completed actual leading mean source identity required')
        super().__init__(nodes,dict(coordinate=rm['chart']['coordinate'],loop=l['charts'],N=rm['N']))
        self.authorizations={};self.alias_hits=set();self.velocity_factor_rows=[]
        for row in theorem['rows']:
            k=row['ordinary_Z_order']
            if k not in (0,1) or not row['passed']:raise ValueError('Both actual ordinary Uy rows required')
            for side,node,coordinate in (
                ('left',m['charts']['postrestore']['physical']['U_y'][k],l['charts']['postrestore']['coordinate']),
                ('right',rm['physical']['U_y'][k],rm['chart']['coordinate'])):
                if row[side+'_node']!=node:raise ValueError('Foreign physical Uy source port')
                self.authorizations[node]=(coordinate,(self.digest(['proved_Rm_physical_Uy',digest(theorem),k]),None))
        if len(self.authorizations)!=4:raise ValueError('Four actual scoped physical Uy ports required')
        P=l['Pstar'] if 'Pstar' in l else None
        # Verify the accepted U=Pstar*V_N wiring and common Pstar. The actual
        # physical mean derivative is U-m, so only the U ports need aliases.
        for row in theorem['U_rows']:
            k=row['ordinary_Z_order']
            if k not in (0,1,2) or not row['passed']:raise ValueError('Three actual U source rows required')
            for side,chart,physical,coordinate in (
                ('left',l['charts']['postrestore'],m['charts']['postrestore']['physical'],l['charts']['postrestore']['coordinate']),
                ('right',rm['chart'],rm['physical'],rm['chart']['coordinate'])):
                node=physical['U'][k];vn=chart['finite_N']['V_N'][k]
                if row[side+'_node']!=node:raise ValueError('Foreign physical U source port')
                product=nodes[node]
                if product['operation']!='product' or vn not in product['arguments']:
                    raise ValueError('Actual physical U product must contain its V_N function')
                factors=[i for i in product['arguments'] if i!=vn]
                if len(factors)!=1:raise ValueError('Physical U must have one Pstar factor')
                if P is None:P=factors[0]
                if factors!=[P]:raise ValueError('Common physical Pstar required')
                self.authorizations[node]=(coordinate,(self.digest(['proved_Rm_physical_U',digest(theorem),k]),None))
                self.velocity_factor_rows.append(dict(side=side,ordinary_Z_order=k,U_node=node,
                    V_N_node=vn,Pstar_node=P,passed=True))
        for row in mean['rows']:
            k=row['ordinary_Z_order']
            if k not in (0,1,2) or not row['passed']:raise ValueError('Three actual leading mean identities required')
            for side,coordinate in (('left',l['charts']['postrestore']['coordinate']),('right',rm['chart']['coordinate'])):
                self.authorizations[row[side+'_node']]=(coordinate,(self.digest(['proved_Rm_leading_mean',digest(mean),k]),None))

    def at(self,node,scope=None):
        if node in self.authorizations:
            coordinate,value=self.authorizations[node]
            if scope is not None and scope.get(coordinate)==self.one:
                self.alias_hits.add(node);return value
        return super().at(node,scope)


class OriginalRmInletM15(OriginalActiveViscousM15):
    def __init__(self,*,directory=DIRECTORY,precision=150):
        self.directory=directory;self.ctx=mp.mp.clone();self.ctx.dps=precision
        self.iv=MPIntervalContext();self.iv.dps=precision+30;self.hashes={};self.code_hashes={}
        raw=(directory/PARENT).read_bytes()
        if hashlib.sha256(raw).hexdigest()!=PARENT_SHA:raise ValueError('Accepted long/reference state required')
        self.parent=json.loads(raw);self.family=self.parent['source_family']
        if not self.parent['checks']['passed']:raise ValueError('Checked parent required')
        self.merge(self.parent['input_hashes'],directory,self.hashes)
        self.merge(self.parent['implementation_hashes'],ROOT,self.code_hashes);self.hashes[PARENT]=PARENT_SHA
        join_receipt=PREFIX+'compliant_current_original_Rm_outer_mixed_source_joins_check.json'
        receipt_raw=(directory/join_receipt).read_bytes()
        if self.hashes.get(join_receipt)!=hashlib.sha256(receipt_raw).hexdigest():
            raise ValueError('Pinned accepted Rm source-join receipt required')
        receipt=json.loads(receipt_raw)
        if not receipt['all_passed'] or receipt['source_family']!=self.family:
            raise ValueError('Checked same-family Rm source-join receipt required')
        self.merge(receipt['input_hashes'],directory,self.hashes)
        self.modified=self.archive('modified_Q4_pressure_C3_functions')
        self.loop=self.archive('prescribed_loop_C5_functions')
        self.outer=self.archive('prescribed_R100_Rm_C6_functions')
        self.leading=self.archive('true_Rm_patch_C6_functions')
        self.actual=self.archive('true_Rm_patch_finite_N_C5_functions')
        self.nodes=self.actual['exact_graph_nodes'];self.rm=self.actual['Rm_patch_finite_N_C5_bindings']
        self.m=self.modified['modified_Q4_pressure_C3_bindings'];self.l=self.loop['prescribed_loop_C5_bindings']
        self.leading_binding=self.leading['true_Rm_patch_C6_bindings']
        for source in (self.modified,self.loop,self.outer,self.leading):
            if self.nodes[:len(source['exact_graph_nodes'])]!=source['exact_graph_nodes']:
                raise ValueError('Identical actual source prefix required')
        sys.setrecursionlimit(max(sys.getrecursionlimit(),12000))
        self.flat=build_Rm_inlet_flat_certificate(self.nodes,self.rm,self.l,root=ROOT)
        if not self.flat['passed']:raise ValueError('Actual compact beta/incoming history proof incomplete')
        from .original_Rm_native_guard import build_Rm_native_guard
        self.native_guards=build_Rm_native_guard(self.nodes,self.rm,self.leading_binding,self.l,root=ROOT)
        if not self.native_guards['passed']:raise ValueError('Actual variable-J source guard incomplete')
        from .original_Rm_inlet_source_theorem import build_Rm_inlet_source_theorem
        self.theorem=build_Rm_inlet_source_theorem(self.nodes,self.rm,self.leading_binding,self.l,self.m,
            self.flat,self.native_guards,root=ROOT,hashes=self.hashes,code_hashes=self.code_hashes)
        if not self.theorem['passed']:raise ValueError('Actual physical source theorem incomplete')
        from .original_Rm_inlet_mean import build_Rm_inlet_mean_certificate
        self.mean=build_Rm_inlet_mean_certificate(self.nodes,self.leading_binding,self.rm,self.l,self.m,
            self.outer['prescribed_R100_Rm_C6_bindings'],root=ROOT)
        if not self.mean['passed']:raise ValueError('Actual leading mean inlet identity incomplete')
        self.recovery=self.recover_traces()
        self.seam_proof=dict(passed=self.recovery['passed'],left='postrestore',right='Rm_patch',
            left_native_endpoint=1,right_native_endpoint=1,ordinary_Z_orders=[0,1],
            actual_trace_recovery=self.recovery,source_derivative_theorem=self.theorem,
            physical_R_once=True,global_fast_phase_not_reset=True,
            actual_nonconstant_native_radius_germs_not_constant_J_substitutes=True,
            incoming_histories_retained=True,original_leading_Qy_subtraction_common=True,no_caps_used=True)
        coordinate=self.rm['chart']['coordinate'];domain=self.rm['chart']['domain']
        self.seam_proof['inlet_source_seam_scope']=dict(chart='actual_patch',coordinate_node=coordinate,
            endpoint_node=domain[0],native_argument=dict(kind='exact_decimal_Rm_patch',value='1'),
            native_x='1',endpoint='Rm',source_seam_role='postrestore native1 to Rm patch native x=1')
        self.seam_proof['external_outlet_scope']=dict(chart='actual_patch',coordinate_node=coordinate,
            endpoint_node=domain[1],native_argument=dict(kind='exact_exp1',value='exp(1)'),
            native_x='Rh',endpoint='Rh',external_port_role='upper signed endpoint')
        for name in ('original_Rm_inlet_m15.py','original_Rm_inlet_flat.py','original_Rm_native_guard.py',
                     'original_Rm_inlet_source_theorem.py','original_Rm_inlet_mean.py','original_Rm_inlet_numeric_m15.py'):
            path='src/openai_ns_reconstruction/'+name
            self.code_hashes[path]=hashlib.sha256((ROOT/path).read_bytes()).hexdigest()

    def recover_traces(self):
        f=RmTraceFingerprint(self.nodes,self.l,self.m,self.rm,self.theorem,self.mean)
        post=self.l['charts']['postrestore'];patch=self.rm['chart']
        scope={post['coordinate']:f.one,patch['coordinate']:f.one};rows=[]
        for key,count in (('U',3),('m',3),('m_y',3),('U_y',2),('Q',2),('Q_y',2)):
            for k in range(count):
                left=self.m['charts']['postrestore']['physical'][key][k];right=self.rm['physical'][key][k]
                a,b=f.at(left,scope),f.at(right,scope)
                if a[0]!=b[0]:raise ArithmeticError('Unproved actual Rm inlet trace '+str((key,k,left,right)))
                rows.append(dict(quantity=key,ordinary_Z_order=k,left_node=left,right_node=right,
                    common_function_sha256=a[0]))
        for k in (0,1):
            left=self.m['charts']['postrestore']['finite_boundary']['R_Q_y'][k];right=self.rm['boundary']['R_Q_y'][k]
            a,b=f.at(left,scope),f.at(right,scope)
            if a[0]!=b[0]:raise ArithmeticError('Unproved once-R weighted actual Rm inlet trace '+str(k))
            rows.append(dict(quantity='R_Q_y',ordinary_Z_order=k,left_node=left,right_node=right,
                common_function_sha256=a[0]))
        leading_rows=[]
        old=self.outer['prescribed_R100_Rm_C6_bindings']['charts']['post_restore']
        new=self.leading_binding['chart']
        for group,key,count in (('fields','V',3),('fields','m',3),('radial','V_y',2)):
            for k in range(count):
                left=old[group][key][k];right=new[group][key][k]
                a,b=f.at(left,scope),f.at(right,scope)
                if a[0]!=b[0]:raise ArithmeticError('Unproved original leading Qy input '+str((key,k,left,right)))
                leading_rows.append(dict(quantity=key,ordinary_Z_order=k,left_node=left,right_node=right,
                    common_function_sha256=a[0]))
        if f.alias_hits!=set(f.authorizations):raise ValueError('All authorized actual source ports must be consumed')
        uy_ports={r[side+'_node'] for r in self.theorem['rows'] for side in ('left','right')}
        return dict(passed=True,rows=rows,authorized_physical_Uy_ports=sorted(uy_ports),
            all_source_ports_consumed=sorted(f.alias_hits),
            physical_U_V_N_factor_rows=f.velocity_factor_rows,
            leading_Qy_input_rows=leading_rows,
            leading_Qy_source_identity='[2Z Uy-(1-delta)Z(U-m)-(1-Z^2)(Uz-mz)]/(1-delta Z^2); U=Pstar*V, m=Pstar*mean',
            original_leading_Qy_subtraction_common=True,
            original_graph_not_mutated=True,endpoint_scope_restricted=True,
            ordinary_Z2_U_and_mean_retained=True,full_raw_Uy_fingerprint=False,
            source_theorem_not_numeric_accuracy=True)

    def evaluate(self,Z):
        if not hasattr(self,'numeric'):raise ValueError('Build the new correlated view first')
        sample=next((s for s in self.numeric['samples'] if s['Z']==Z),None)
        if sample is None:raise ValueError('Only saved axial queries admitted; original sources not evaluated')
        return sample

    def build(self):
        from .original_Rm_inlet_numeric_m15 import build_Rm_inlet_numeric
        self.numeric=build_Rm_inlet_numeric(self,self.seam_proof)
        return dict(source_family=self.family,parent_report=PARENT,parent_report_sha256=PARENT_SHA,
            compact_flat_and_history_certificate=self.flat,native_source_guards=self.native_guards,
            source_derivative_theorem=self.theorem,leading_mean_certificate=self.mean,actual_trace_recovery=self.recovery,
            source_seam_proof=self.seam_proof,numerical_source_binding=self.numeric['source_binding'],
            samples=self.numeric['samples'],contract=dict(accepted_ancestor_constructors_not_run=True,
                source_functions_not_cap_values=True,raw_213_and_independent_target_retained=True),
            limits=dict(original_signed_finite_N_endpoint_values=False,numeric_Rp_closure=False,
                global_all_five_common_N_closure=False,exact_heat_energy_stress_flat=False,
                true_temporal_n_recursion_summation=False,mean_oscillatory_correction=False,
                independent_corrected_Cartesian_NS=False,full_goal_complete=False))
