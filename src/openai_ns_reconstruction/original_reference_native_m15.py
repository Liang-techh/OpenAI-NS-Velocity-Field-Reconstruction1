"""Physical native-derivative theorem for the original reference inlet.

The proof consumes source functions and pinned calculus rules. It never
selects an original field value from an error cap or sets the fast phase to
an endpoint. Raw m15 provenance and every inherited error budget survive.
"""
from fractions import Fraction
import copy
import hashlib
import json
import sys

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from .original_core_numeric import ROOT,DIRECTORY,PREFIX
from .original_active_viscous_m15 import OriginalActiveViscousM15,ActiveEndpointFingerprint

PHASE=PREFIX+'compliant_current_original_reference_phase_lift.json'
PHASE_SHA='e1640aaa18b99fb376ebdf5e97b9cfa333b6eab567da24020d11d20b0bffdecf'
GERM=PREFIX+'compliant_current_original_reference_physical_germ.json'
GERM_SHA='a29bf10e7fa2fe6f17cc39505790596b6864f788061d8ee36fc68e6dfc18769e'


def canonical(value):return json.dumps(value,ensure_ascii=False,sort_keys=True,separators=(',',':'))
def digest(value):return hashlib.sha256(canonical(value).encode('utf-8')).hexdigest()


def root_tangent_proof():
    """Actual seam two-jets give zero physical tangent of a,b,t0,Delta,q."""
    import sympy as s
    E=s.symbols('E',positive=True);Ey,Eyy,Vy,Vyy=s.symbols('Ey Eyy Vy Vyy',real=True)
    a,b,eta,Delta,q=s.symbols('a b eta Delta q',real=True,nonzero=True)
    ay=-2*(Eyy*E-Ey**2)/E**2
    by=2*(Vyy*E-Vy*Ey)/E**2
    ty=-(by*a-b*ay)/a**2
    dy=ay+(2*b*by*a-b**2*ay)/a**2
    qy=((-dy)*a-(2*eta-Delta)*ay)/(4*a**2*q)
    actual={Ey:E/10,Eyy:E/100,Vy:0,Vyy:0,a:s.Rational(4,5),b:0,Delta:s.Rational(-6,5)}
    expressions=dict(a_y=ay,b_y=by,t0_y=ty,Delta_y=dy,q_y=qy)
    rows=[dict(quantity=k,formula=str(v),after_actual_two_jet_substitution=str(s.simplify(v.subs(actual))),
               exact_zero=s.simplify(v.subs(actual))==0) for k,v in expressions.items()]
    Ez,Eyz,Eyyz,Vyz,Vyyz,az,bz,dz,qz=s.symbols('Ez Eyz Eyyz Vyz Vyyz az bz dz qz',real=True)
    dual=dict(zip((E,Ey,Eyy,Vy,Vyy,a,b,Delta,q),(Ez,Eyz,Eyyz,Vyz,Vyyz,az,bz,dz,qz)))
    actual_Z={Eyz:Ez/10,Eyyz:Ez/100,Vyz:0,Vyyz:0,az:0,bz:0,dz:0,qz:0}
    mixed=[]
    for key,value in expressions.items():
        derivative=sum(s.diff(value,x)*zx for x,zx in dual.items())
        substituted=s.simplify(derivative.subs(actual_Z).subs(actual))
        mixed.append(dict(quantity=key+'Z',ordinary_Z_derivative_taken_before_source_substitutions=True,
            formula=str(derivative),after_actual_source_jets=str(substituted),exact_zero=substituted==0))
    return dict(passed=all(r['exact_zero'] for r in rows+mixed),rows=rows,ordinary_Z1_rows=mixed,
        actual_Ey='E/10',actual_Eyy='E/100',actual_Vy_Vyy='0 by flat original sigma at the seam',
        source_q_active='sigma=1, q=sqrt((2*eta-Delta)/(2*a))>0',
        ordinary_Z_tangents='actual source a_Z,b_Z,t0_Z,Delta_Z,q_Z certificate is required by the owner; mixed yZ rows here use those values explicitly',
        Z_values_selected=False)


def response_definition(germ,l,outer):
    """Explicit common physical derivative, retaining ordinary Z1.

    H=B/E with a=4/5,b=t0=0 and fixed q,eta,d_star. Function recipe IDs
    refer to the accepted germ recipe graph, never the original graph.
    """
    f=outer['prescribed_R100_Rm_C6_bindings']['charts']['reference']['fields']
    roots=l['charts']['reference']['roots'];r=germ['physical_germ_recipes']['root_p2']
    ports={name:dict(kind='original_endpoint_function',original_node=node,ordinary_Z_order=k,
                   native_coordinate=l['charts']['reference']['coordinate'],native_endpoint=1)
           for name,node,k in (('E',f['E'][0],0),('Ez',f['E'][1],1),('p',roots['p2'][0],0),('pz',roots['p2'][1],1))}
    for name,node,k in (('py',r['physical_y_Z0_Z1'][0],0),('pyz',r['physical_y_Z0_Z1'][1],1)):
        ports[name]=dict(kind='accepted_physical_recipe_function',parent_report=GERM,parent_report_sha256=GERM_SHA,
            local_recipe_node=node,ordinary_Z_order=k,not_an_original_graph_node=True)
    ports['q']=dict(kind='original_endpoint_function',original_node=l['charts']['reference']['loop']['q'][0],
                    native_coordinate=l['charts']['reference']['coordinate'],native_endpoint=1,ordinary_Z_order=0,
                    ordinary_Z1_node=l['charts']['reference']['loop']['q'][1],
                    ordinary_Z1_zero_requires_actual_source_certificate=True)
    ports['dstar']=dict(kind='original_source_parameter',original_node=l['dstar'])
    return dict(ports=ports,physical_coordinate='y=log(R/110); d_y=R*d_R',
        held_fast_phase='arbitrary actual fractional_part(N*global_offset), held only in the slow derivative',
        definitions=[
            'a=4/5; b=t0=0; u=p*q/dstar; r=u/sqrt(1+u^2); alpha=2*q/sqrt(1+u^2)',
            't(p,theta)=alpha*(cos(theta)-r)/(1-2*r*cos(theta)+r^2)',
            'K=1/(2*pi*(1+2*q^2)); Phi(p,theta)=K*(theta+int_0^theta t(p,s)^2 ds)',
            'psi=unique inverse Phi(p,psi)=phase; lambda=K*(1+t(p,psi)^2)>0',
            'I_j=int_0^psi partial_p^j t(p,s) ds (j=0,1,2)',
            'Phi_p=K*int_0^psi 2*t*t_p ds; Phi_pp=K*int_0^psi 2*(t_p^2+t*t_pp) ds',
            'psi_p=-Phi_p/lambda; psi_pp=-(Phi_pp+4*K*t*t_p*psi_p+2*K*t*t_theta*psi_p^2)/lambda',
            'H=-a*I_0/(4*pi)',
            'H_p=-a*(I_1+t*psi_p)/(4*pi)',
            'H_pp=-a*(I_2+2*t_p*psi_p+t_theta*psi_p^2+t*psi_pp)/(4*pi)',
            'G=slow_B_y=(E/10)*H+E*H_p*py',
            'G_Z=(Ez/10)*H+(E/10)*H_p*pz+Ez*H_p*py+E*H_pp*pz*py+E*H_p*pyz'],
        ordinary_Z_derivatives_unscaled=True,phase_Z_exact_zero=True,
        true_Z2_field_ports_consumed_by_pyz=True,no_actual_phase_or_source_values_selected=True,
        quantitative_inverse_conditioning_claimed=False)


def evaluate_synthetic_response(*,E,Ez,p,pz,py,pyz,q,dstar,phase,precision=45):
    """Evaluate only caller-supplied moderate parameters for focused checks."""
    c=mp.mp.clone();c.dps=precision
    E,Ez,p,pz,py,pyz,q,dstar,phase=[c.mpf(v) for v in (E,Ez,p,pz,py,pyz,q,dstar,phase)]
    if not(E>0 and q>0 and dstar>0 and 0<phase<1):raise ValueError('Positive synthetic parameters and interior phase required')
    a=c.mpf(4)/5;K=1/(2*c.pi*(1+2*q*q))
    def direction(pp,theta):
        u=pp*q/dstar;root=c.sqrt(1+u*u);r=u/root;alpha=2*q/root
        return alpha*(c.cos(theta)-r)/(1-2*r*c.cos(theta)+r*r)
    t=lambda x:direction(p,x)
    tp=lambda x:c.diff(lambda pp:direction(pp,x),p)
    tpp=lambda x:c.diff(lambda pp:direction(pp,x),p,2)
    Phi=lambda x:K*(x+c.quad(lambda z:t(z)**2,[0,x]))
    psi=c.findroot(lambda x:Phi(x)-phase,(2*c.pi*phase,c.pi),solver='secant')
    if not(0<psi<2*c.pi):raise ArithmeticError('Synthetic inverse left its branch')
    I0=c.quad(t,[0,psi]);I1=c.quad(tp,[0,psi]);I2=c.quad(tpp,[0,psi])
    T,Tp,Ttheta=t(psi),tp(psi),c.diff(t,psi);lam=K*(1+T*T)
    phip=K*c.quad(lambda x:2*t(x)*tp(x),[0,psi])
    phipp=K*c.quad(lambda x:2*(tp(x)**2+t(x)*tpp(x)),[0,psi])
    psip=-phip/lam;psipp=-(phipp+4*K*T*Tp*psip+2*K*T*Ttheta*psip**2)/lam
    H=-a*I0/(4*c.pi);Hp=-a*(I1+T*psip)/(4*c.pi)
    Hpp=-a*(I2+2*Tp*psip+Ttheta*psip**2+T*psipp)/(4*c.pi)
    G=E*H/10+E*Hp*py
    Gz=Ez*H/10+E*Hp*pz/10+Ez*Hp*py+E*Hpp*pz*py+E*Hp*pyz
    return dict(H=H,Hp=Hp,Hpp=Hpp,G=G,Gz=Gz,psi=psi)


class ProvedNativeEndpointFingerprint(ActiveEndpointFingerprint):
    """Use a proved function identity only in its exact native endpoint scope."""
    def __init__(self,nodes,bindings,modified_bindings,theorem,response):
        if not theorem['passed']:raise ValueError('A completed source derivative theorem is required')
        if not theorem.get('actual_root_parameter_Z1_certificate_consumed') or theorem.get('mixed_root_yZ_rows_consumed')!=5:
            raise ValueError('Actual root-Z1 and mixed-yZ premises are required')
        if theorem.get('response_definition_sha256')!=digest(response):raise ValueError('Proved response definition mismatch')
        if len(theorem['rows'])!=2 or {r['ordinary_Z_order'] for r in theorem['rows']}!={0,1}:raise ValueError('Exactly ordinary-Z0/Z1 theorem rows required')
        for row in theorem['rows']:
            k=row['ordinary_Z_order']
            if not row['passed'] or row['quantity']!='slow_B_y':raise ValueError('Unproved function authorization')
            for side,name in (('left','reference'),('right','restoration')):
                if row[side+'_node']!=modified_bindings['charts'][name]['radial']['slow_B_y'][k]:raise ValueError('Wrong original native derivative port')
        super().__init__(nodes,dict(coordinate=bindings['charts']['reference']['coordinate'],loop=bindings['charts']))
        self.authorizations={};self.alias_hits=set();self.outside_endpoint_scope=set()
        response_sha=digest(response)
        for row in theorem['rows']:
            k=row['ordinary_Z_order']
            for side,name in (('left','reference'),('right','restoration')):
                chart=bindings['charts'][name]
                endpoint=self.one if name=='reference' else self.zero
                self.authorizations[row[side+'_node']]=(chart['coordinate'],endpoint,
                    (self.digest(['proved_common_physical_B_y_function',response_sha,k]),None))
    def at(self,node,scope=None):
        if node in self.authorizations:
            coordinate,endpoint,function=self.authorizations[node]
            if scope is not None and scope.get(coordinate)==endpoint:
                self.alias_hits.add(node);return function
            self.outside_endpoint_scope.add(node)
        return super().at(node,scope)


class OriginalReferenceNativeM15(OriginalActiveViscousM15):
    def __init__(self,*,directory=DIRECTORY,precision=100):
        self.directory=directory;self.ctx=mp.mp.clone();self.ctx.dps=precision
        self.iv=MPIntervalContext();self.iv.dps=precision+30;self.hashes={};self.code_hashes={}
        def packet(name,expected):
            raw=(directory/name).read_bytes()
            if hashlib.sha256(raw).hexdigest()!=expected:raise ValueError('Changed accepted packet '+name)
            p=json.loads(raw)
            if not p['checks']['passed']:raise ValueError('Checked packet required')
            self.merge(p['input_hashes'],directory,self.hashes);self.merge(p['implementation_hashes'],ROOT,self.code_hashes)
            self.hashes[name]=expected;return p
        self.phase=packet(PHASE,PHASE_SHA);self.germ=packet(GERM,GERM_SHA)
        self.family=self.phase['source_family']
        self.parent=self.saved('active_viscous_m15_numeric')
        self.modified=self.archive('modified_Q4_pressure_C3_functions');self.loop=self.archive('prescribed_loop_C5_functions')
        self.outer=self.archive('prescribed_R100_Rm_C6_functions');self.nodes=self.modified['exact_graph_nodes']
        self.l=self.loop['prescribed_loop_C5_bindings'];self.m=self.modified['modified_Q4_pressure_C3_bindings']
        for archive in (self.loop,self.outer):
            if self.nodes[:len(archive['exact_graph_nodes'])]!=archive['exact_graph_nodes']:raise ValueError('Common source prefixes required')
        sys.setrecursionlimit(max(sys.getrecursionlimit(),12000))
        from .original_native_derivative_guards import build_native_derivative_guards
        from .original_reference_history_seam import build_reference_history_seam
        self.derivative_guards=build_native_derivative_guards(self.nodes,self.l,self.modified,self.outer,self.germ['source_guards'],ROOT)
        self.history=build_reference_history_seam(self.nodes,self.l,self.modified,self.outer)
        self.tangents=root_tangent_proof();self.response=response_definition(self.germ,self.l,self.outer)
        from .original_reference_root_Z import build_root_Z_certificate
        self.root_Z=build_root_Z_certificate(self.nodes,self.l,self.phase['phase_lift_certificate'],self.germ)
        for name in ('original_reference_native_m15.py','original_native_derivative_guards.py','original_reference_history_seam.py','original_reference_root_Z.py'):
            key='src/openai_ns_reconstruction/'+name;self.code_hashes[key]=hashlib.sha256((ROOT/key).read_bytes()).hexdigest()

    def evaluate(self,*args,**kwargs):
        raise ValueError('Original numerical source values are not evaluated by this function-theorem owner')

    def theorem(self):
        periodic=self.phase['phase_lift_certificate']
        physical=self.germ['physical_two_jet_certificate']
        passed=bool(periodic['passed'] and periodic['global_periodic_lift_proved'] and physical['passed'] and
            self.derivative_guards['passed'] and self.history['passed'] and self.tangents['passed'] and self.root_Z['passed'])
        rows=[]
        for k in (0,1):
            left=self.m['charts']['reference'];right=self.m['charts']['restoration']
            rows.append(dict(quantity='slow_B_y',ordinary_Z_order=k,
                left_node=left['radial']['slow_B_y'][k],right_node=right['radial']['slow_B_y'][k],
                common_response='G' if k==0 else 'G_Z',common_definition_sha256=digest(self.response),
                method='source_chain_rule_naturality_with_common_physical_two_jets',passed=passed))
        return dict(passed=passed,rows=rows,response_definition_sha256=digest(self.response),
            actual_root_parameter_Z1_certificate_consumed=self.root_Z['passed'],
            mixed_root_yZ_rows_consumed=len(self.tangents['ordinary_Z1_rows']),
            proof_method='source-bound calculus deduction, not full raw-graph canonical equality',
            coordinate_maps=dict(reference='y=T+J_ref*q',restoration='y=logref-8+t'),
            native_Jacobians=[self.l['charts']['reference']['physical_measure'],self.l['charts']['restoration']['physical_measure']],
            proof=[
                'The admitted leading input fields agree as slow-Z functions at the seam. Their physical first and second y jets agree by the retained Duhamel ODEs, E_y=E/10 and the flat restoration sigma.',
                'The original root formulas therefore have the same physical first jets. Exact root a,b,t0,Delta,q tangents vanish; p2_y and p2_yZ remain the original nonzero source recipes.',
                'NativeDerivative is the source partial derivative at held actual fast phase, with each source field seeded J times its physical radial derivative. Quotient, moving-limit, substitution and implicit-inverse rules are pinned to that source.',
                'Structural induction through those calculus operations gives D_native(B)/J equal to the explicit common response G for each chart. The smooth periodic lift makes the same identity valid at phase integer cuts; ordinary Z1 commutes and yields G_Z.',
                'The original finite-N m density increment is delta_V and its Duhamel rate is1 with constant J and retained memory; hence corrected m_y=V_modified-m_modified, as used by physical Qy recovery.',
                'No native chart endpoint is substituted for a fast-phase endpoint, no slow derivative is set to zero, and the two Jacobians are never equated.'
            ],direct_full_modified_Uy_fingerprint=False,arbitrary_fast_phase=True,
            original_p2_y_not_replaced_by_zero=True,no_error_caps_used=True)

    def trace_recovery(self,theorem):
        f=ProvedNativeEndpointFingerprint(self.nodes,self.l,self.m,theorem,self.response)
        a,b=(self.l['charts'][n] for n in ('reference','restoration'))
        scope={a['coordinate']:f.one,b['coordinate']:f.zero}
        rows=[]
        for group,key,orders in (('physical','U',3),('physical','m',3),('physical','m_y',3),
                                ('physical','U_y',2),('physical','Q',2),('physical','Q_y',2),
                                ('finite_boundary','R_Q_y',2)):
            for k in range(orders):
                left=self.m['charts']['reference'][group][key][k]
                right=self.m['charts']['restoration'][group][key][k]
                lhash,rhash=f.at(left,scope)[0],f.at(right,scope)[0]
                if lhash!=rhash:raise ArithmeticError('Unproved modified trace after source theorem '+str((key,k,left,right)))
                rows.append(dict(quantity=key,ordinary_Z_order=k,left_node=left,right_node=right,
                    common_endpoint_sha256=lhash,method='endpoint_arithmetic_with_proved_scope_restricted_B_y_identity'))
        expected={r[side+'_node'] for r in theorem['rows'] for side in ('left','right')}
        if f.alias_hits!=expected:raise ArithmeticError('Every proved slow-B derivative port must actually be consumed')
        return dict(passed=True,rows=rows,authorized_derivative_function_ports=sorted(f.alias_hits),
            actual_native_endpoint_scope_required=True,original_graph_not_mutated=True,
            proof_authorized_function_identity_not_unconditional_raw_fingerprint=True,
            original_R_and_Pstar_source_factors_retained=True,ordinary_Z2_U_mean_retained=True,
            empty_integrals=sorted(f.empty),numeric_values_or_caps_used=False)

    def build(self):
        proof=self.theorem()
        if not proof['passed']:raise ValueError('Native derivative premises incomplete')
        recovery=self.trace_recovery(proof)
        return dict(source_family=self.family,parent_report=PHASE,parent_report_sha256=PHASE_SHA,
            native_derivative_source_guards=self.derivative_guards,corrected_history_seam=self.history,
            root_tangent_proof=self.tangents,actual_root_Z1_certificate=self.root_Z,
            common_physical_response=self.response,native_derivative_theorem=proof,
            actual_modified_trace_recovery=recovery,
            unchanged_m15_ledger=copy.deepcopy(self.phase['unchanged_m15_ledger']),
            contract=dict(accepted_ancestor_constructors_not_run=True,raw_213_reduced_199_and_all_budgets_unchanged=True,
                source_function_recipe_not_original_numeric_value=True,no_error_cap_as_value=True,
                actual_phase_not_selected_or_assumed_zero=True),
            limits=dict(actual_native_B_y_naturality_proved=proof['passed'],
                actual_modified_Uy_Qy_RQy_seam_proved=recovery['passed'],new_budget_contraction=False,
                independent_numeric_Rp_seam=False,global_five_moment_common_N_closure=False,
                exact_heat_and_energy=False,admissible_stress_flat_remainder=False,
                temporal_n_recursion_and_summation=False,mean_oscillatory_correction=False,
                independent_corrected_Cartesian_NS=False))
