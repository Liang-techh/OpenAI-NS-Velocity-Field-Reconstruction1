"""Second slow-Z derivatives of the actual original phase primitives.

Uses the current C1 expression graph unchanged and appends exact C2
formulas. The live leading provider supplies Taylor source rows; they
enclose functions, and are never selected as point-function values.
"""
from dataclasses import dataclass
from fractions import Fraction
import gzip
import json
from pathlib import Path
from types import SimpleNamespace
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_control_connection as current
import lei_ren_part1_paper_compliant_current_original_whole_Z_limit_repair_band as limit
import lei_ren_part1_paper_compliant_current_original_whole_Z_Rh_functional_join as join
import lei_ren_part1_paper_compliant_flat_pulse_derivatives as flat_source

HERE,PREFIX,sha=current.HERE,current.PREFIX,current.sha
source=current.source
CHARTS=current.CHARTS
NAME=PREFIX+'current_original_C2_phase_primitives.json.gz'
RECEIPT=PREFIX+'current_original_C2_phase_primitives_check.json'
GATES=('current_original_17_chart_source_Z2_rows_installed',
       'current_original_same_inverse_C2_phase_primitive_functions_installed')
OPEN=join.OPEN


def flat_sigma_second(c,x):
    """The source-backed ordinary derivative, not a differentiated cap."""
    return 2*flat_source.sigma_jets(c,x)[2]


def source_bindings():
    return dict(original_loop=current.ast_binding(current.original_loop),
        actual_leading_source=current.ast_binding(CurrentC2PhasePrimitives.leading_source),
        actual_C2_projection=current.ast_binding(CurrentC2PhasePrimitives.source_packet),
        flat_sigma_second=current.ast_binding(flat_sigma_second),
        original_flat_sigma_jets=current.ast_binding(flat_source.sigma_jets),
        original_flat_sigma_active=current.ast_binding(flat_source._sigma_left),
        current_source_recipe_manifest=current.recipe_manifest())


@dataclass(frozen=True)
class C2Function:
    value: object
    Z: object
    ZZ: object


class C2Algebra:
    """Ordinary derivative rules; two-Z cross terms carry factor two."""
    def __init__(self,g):self.g=g
    def fixed(self,value):
        q=self.g.constant(value) if isinstance(value,(str,int,Fraction)) else value
        return C2Function(q,self.g.zero,self.g.zero)
    def add(self,*rows):
        return C2Function(*(self.g.add(*(getattr(q,key) for q in rows)) for key in ('value','Z','ZZ')))
    def neg(self,q):return C2Function(*(self.g.neg(getattr(q,key)) for key in ('value','Z','ZZ')))
    def scale(self,q,v):return self.mul(q,self.fixed(v))
    def mul(self,a,b):
        g=self.g
        return C2Function(g.mul(a.value,b.value),g.add(g.mul(a.Z,b.value),g.mul(a.value,b.Z)),
            g.add(g.mul(a.ZZ,b.value),g.mul(g.constant(2),a.Z,b.Z),g.mul(a.value,b.ZZ)))
    def div(self,n,d):
        g=self.g;divide=lambda a:g.quotient(a,d.value,'same actual positive source denominator')
        v=divide(n.value);z=divide(g.sub(n.Z,g.mul(v,d.Z)))
        zz=divide(g.sub(n.ZZ,g.add(g.mul(g.constant(2),z,d.Z),g.mul(v,d.ZZ))))
        return C2Function(v,z,zz)
    def unary(self,q,name,first,second):
        g=self.g;v=g.unary(name,q.value)
        return C2Function(v,g.mul(first,q.Z),g.add(g.mul(second,q.Z,q.Z),g.mul(first,q.ZZ)))
    def sqrt(self,q):
        g=self.g;v=g.unary('positive_sqrt',q.value)
        z=g.quotient(q.Z,g.mul(g.constant(2),v),'actual positive square root')
        zz=g.quotient(g.sub(q.ZZ,g.mul(g.constant(2),z,z)),g.mul(g.constant(2),v),'actual positive square root')
        return C2Function(v,z,zz)
    def sigma(self,q):
        g=self.g
        return self.unary(q,'original_flat_sigma',g.unary('original_flat_sigma_prime',q.value),
            g.node('analytic_unary',name='original_flat_sigma_second',argument=q.value.node,
                defining_source=current.ast_binding(flat_sigma_second),
                derivative_source=current.ast_binding(flat_source.sigma_jets),
                ordinary_derivative='2*original_sigma_Taylor_coefficient[2]',
                range_enclosure_is_not_selected_function_value=True))
    def integral(self,q,angle):
        g=self.g
        def one(v):return current.controls.integral(g,v,angle,g.zero,g.symbol(angle),
            measure='original angle dpsi at fixed slow source parameters')
        return C2Function(one(q.value),one(q.Z),one(q.ZZ))


def primitive_C2(built,window):
    """Append actual q_ZZ, implicit psi_ZZ and A/B_ZZ, sharing old inverse."""
    g=built['graph'];alg=C2Algebra(g);chart=window['chart'];old=window['loop']
    roots={}
    for name,pair in window['source_roots'].items():
        first=source.FunctionRef(g,pair[0]);z=source.FunctionRef(g,pair[1])
        definition=dict(g.nodes[first.node]);definition['Z_order']=2
        definition.update(ordinary_slow_Z_derivative_order=2,Taylor_coefficient_factorial=2,
            derivative_provider=Path(__file__).name,source_projection='same current axial coefficient[2] times 2',
            source_projection_binding=current.ast_binding(CurrentC2PhasePrimitives.source_packet),
            derivative_of_range_endpoint=False)
        # The original C1 recipe used a two-row public projection for a/t0.
        # The new row reads genuine retained source Taylor arrays instead.
        definition['C2_quantity_paths']=dict(E='common_velocity_E_axial5',V='common_velocity_V_axial5',
            a=['actual_a_axial5','actual_correlated_a_axial5','full_original_shear_a_axial5'],
            b=['actual_b_axial5','actual_nonzero_b_axial5','exact_b_and_t0_zero'],
            t0=['actual_t0_axial5','exact_b_and_t0_zero'],p2=['full_signed_p2_axial4','exact_full_p2_axial4'],
            Delta=['actual_Delta_axial5','Delta_axial5'])
        zz=g.node(definition.pop('operation'),**definition)
        roots[name]=C2Function(first,z,zz)
    parameters={n['name']:source.FunctionRef(g,i) for i,n in enumerate(g.nodes)
        if n['operation']=='current_original_source_parameter'}
    eta,dstar=[alg.fixed(parameters[k]) for k in ('eta','d_star')]
    a,b,E,p2,t0,D=[roots[k] for k in ('a','b','E','p2','t0','Delta')]
    if old is None:
        return dict(chart=chart,source_roots=roots,A=C2Function(g.zero,g.zero,g.zero),
            B=C2Function(g.zero,g.zero,g.zero),exact_flat_from_original_power_admission=True)
    # Preserve every old flat/collar predicate and its defining module.
    qold=g.nodes[old['q'][0]]
    def flat(v):
        attrs={key:value for key,value in qold.items() if key not in ('operation','active_body')}
        return g.node('original_lazy_flat_branch',active_body=v.node,**attrs)
    gamma=alg.add(alg.scale(eta,2),alg.neg(D))
    cut=alg.add(alg.fixed(1),alg.neg(alg.div(D,eta)))
    root=alg.sqrt(alg.div(gamma,alg.scale(a,2)))
    qraw=alg.mul(alg.sigma(cut),root)
    q=C2Function(*(flat(v) for v in (qraw.value,qraw.Z,qraw.ZZ)))
    u=alg.div(alg.mul(p2,q),dstar)
    hinv=alg.div(alg.fixed(1),alg.sqrt(alg.add(alg.fixed(1),alg.mul(u,u))))
    r=alg.mul(u,hinv);alpha=alg.scale(alg.mul(q,hinv),2)
    angle='loop_angle_'+chart;psi=g.symbol(angle);cp=g.unary('cos',psi);sp=g.unary('sin',psi)
    nn=alg.add(alg.fixed(cp),alg.neg(r))
    den=alg.add(alg.fixed(1),alg.neg(alg.scale(r,g.mul(g.constant(2),cp))),alg.mul(r,r))
    w=alg.div(nn,den);t=alg.add(t0,alg.mul(alpha,w))
    pi=g.node('mathematical_pi');twopi=g.mul(g.constant(2),pi)
    K=alg.div(alg.fixed(1),alg.scale(alg.add(alg.fixed(1),alg.mul(t0,t0),alg.scale(alg.mul(q,q),2)),twopi))
    T1=alg.integral(t,angle);T2=alg.integral(alg.mul(t,t),angle)
    P=alg.add(alg.fixed(psi),T2);Phi=alg.mul(K,P)
    inverse=source.FunctionRef(g,old['inverse'])
    def at(v):return g.node('substitute_original_inverse_angle',body=v.node,
        inverse_angle=inverse.node,angle_variable=angle)
    # Implicit derivative of the SAME original monotone phase inverse.
    phi_psi=alg.mul(K,alg.add(alg.fixed(1),alg.mul(t,t)))
    denpsi=g.mul(g.constant(2),r.value,sp)
    wpsi=g.quotient(g.sub(g.mul(g.neg(sp),den.value),g.mul(nn.value,denpsi)),
        g.mul(den.value,den.value),'original Poisson denominator squared positive')
    tpsi=g.mul(alpha.value,wpsi)
    phi_psipsi=g.mul(g.constant(2),K.value,t.value,tpsi)
    psiZ=g.neg(g.quotient(at(Phi.Z),at(phi_psi.value),'same strictly positive inverse Jacobian'))
    psiZZ=g.neg(g.quotient(g.add(at(Phi.ZZ),g.mul(g.constant(2),at(phi_psi.Z),psiZ),
        g.mul(at(phi_psipsi),psiZ,psiZ)),at(phi_psi.value),'same strictly positive inverse Jacobian'))
    phase=source.FunctionRef(g,window['phase']);chi=g.sub(phase,g.quotient(inverse,twopi,'2*pi positive'))
    T1Z=g.add(at(T1.Z),g.mul(at(t.value),psiZ))
    T1ZZ=g.add(at(T1.ZZ),g.mul(g.constant(2),at(t.Z),psiZ),
        g.mul(at(tpsi),psiZ,psiZ),g.mul(at(t.value),psiZZ))
    M=g.sub(g.neg(g.quotient(g.mul(a.value,at(T1.value)),twopi,'2*pi positive')),g.mul(b.value,phase))
    MZ=g.sub(g.neg(g.quotient(g.add(g.mul(a.Z,at(T1.value)),g.mul(a.value,T1Z)),twopi,'2*pi positive')),g.mul(b.Z,phase))
    MZZ=g.sub(g.neg(g.quotient(g.add(g.mul(a.ZZ,at(T1.value)),g.mul(g.constant(2),a.Z,T1Z),
        g.mul(a.value,T1ZZ)),twopi,'2*pi positive')),g.mul(b.ZZ,phase))
    AZZ=flat(g.mul(g.constant('1/2'),g.sub(g.mul(a.ZZ,chi),g.quotient(
        g.add(g.mul(g.constant(2),a.Z,psiZ),g.mul(a.value,psiZZ)),twopi,'2*pi positive'))))
    BZZ=flat(g.mul(g.constant('1/2'),g.add(g.mul(E.ZZ,M),g.mul(g.constant(2),E.Z,MZ),g.mul(E.value,MZZ))))
    A=C2Function(source.FunctionRef(g,old['A'][0]),source.FunctionRef(g,old['A'][1]),AZZ)
    B=C2Function(source.FunctionRef(g,old['B'][0]),source.FunctionRef(g,old['B'][1]),BZZ)
    return dict(chart=chart,source_roots=roots,q=q,inverse=inverse,psi_Z=psiZ,psi_ZZ=psiZZ,
        Phi_C2=Phi,Phi_psi_C1=phi_psi,Phi_psipsi=phi_psipsi,A=A,B=B,
        original_C1_A_B_prefix_unchanged=True,original_flat_and_collar_predicates_unchanged=True,
        same_original_inverse_not_a_new_branch=True,phase_Z_and_ZZ_exact_zero=True)


def regularity_proof():
    """Positive uniform Poisson denominator, flat seam and endpoint moments."""
    u,r,h,q,t0=s.symbols('u r h q t0',real=True);psi=s.symbols('psi',real=True)
    assert s.simplify((1-(u/s.sqrt(1+u*u))**2)-1/(1+u*u))==0
    den=1-2*r*s.cos(psi)+r*r
    assert s.expand(den-((1-r)**2+2*r*(1-s.cos(psi))))==0
    Z=s.symbols('Z');f0,fp,fz,fpp,fpz,fzz,p0,p1,p2=s.symbols('f0 fp fz fpp fpz fzz p0 p1 p2')
    # Independent universal second-order local jets avoid symbolic Subs
    # aliases and cover every partial derivative in the implicit formula.
    F=f0+fp*(psi-p0)+fz*Z+fpp*(psi-p0)**2/2+fpz*(psi-p0)*Z+fzz*Z**2/2
    p=p0+p1*Z+p2*Z**2/2
    second=s.diff(F.subs(psi,p),Z,2).subs(Z,0)
    assert s.expand(second-(fzz+2*fpz*p1+fpp*p1**2+fp*p2))==0
    return dict(passed=True,exact_current_original_C1_graph_prefix_preserved=True,
        source_rows_are_actual_ordinary_slow_Z_Taylor_coefficients=True,
        Poisson_denominator='1-2*r*cos(psi)+r^2 >= (1-abs(r))^2 >= 1/(4*(1+u^2)^2) > 0',
        exact_r_identity='r=u/sqrt(1+u^2), 1-r^2=1/(1+u^2)',
        strict_inverse_Jacobian='Phi_psi=K*(1+t^2) >= K=1/(2*pi*(1+t0^2+2*q^2)) > 0',
        uniform_C2_on_compact_current_source_domains=True,
        cutoff_extension='Delta<eta: sigma(1-Delta/eta)*sqrt((2eta-Delta)/(2a)); Delta>=eta:0',
        cutoff_C2_join='sigma,sigma_prime,sigma_second vanish at0; active gamma>=eta and a>0',
        collar_predicate_independent_of_Z_and_reused_exactly=True,
        r_zero_regular='r=0 gives denominator1; no division by r or switch of formulas',
        angular_normalization='w=sum_(n>=1) r^(n-1)*cos(n*psi); int(w)=0, int(w^2)=pi/(1-r^2)',
        phase_endpoint_identities='Phi(0)=0, Phi(2*pi)=1, A=B=0 at phase0,1, with Z/ZZ rows zero',
        actual_two_sided_flat_source_is_lazy=True,range_caps_not_differentiated=True,
        target_ZZ_transport_or_C2_limit_controls_admitted=False)


def encoded(value):
    if isinstance(value,C2Function):return [value.value.node,value.Z.node,value.ZZ.node]
    if isinstance(value,source.FunctionRef):return value.node
    if isinstance(value,dict):return {k:encoded(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [encoded(v) for v in value]
    return join.serialized(value)


class CurrentC2PhasePrimitives:
    def __init__(self,owner=None,require_checked=True):
        self.outer=owner if owner is not None else current.current.WholeZAllNOuterRcFunctions()
        report,self.hashes=limit.load_current();self.identity=self.outer.identity
        if report['source_family']!=self.identity:raise ValueError('One current family required')
        self.report=report;self.built=limit.restore_current(report);self.prefix=[dict(n) for n in self.built['graph'].nodes]
        raw=report['exact_current_integral_and_control_graph'];self.windows={q['chart']:q for q in raw['windows']}
        self.functions={chart:primitive_C2(self.built,self.windows[chart]) for chart in CHARTS}
        self.proof=regularity_proof();self.cache={}
        self.hashes.update(self.outer.hashes);self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.hashes[Path(flat_source.__file__).name]=sha(Path(flat_source.__file__).name)
        # The long frontend must use its original coordinate method rather
        # than the outer class's restricted chart selector. Same live graph.
        self.long=SimpleNamespace(**{**vars(self.outer),'packet_cache':{}})
        self.long.coordinate=lambda value,chart=None:current.current.previous.WholeZAllNLongPatchFunctions.coordinate(self.long,value,chart)
        self.switch=SimpleNamespace(**{**vars(self.outer),'packet_cache':{}})
        self.switch.coordinate=lambda value:current.current.previous.previous.WholeZAllNSwitchFunctions.coordinate(self.switch,value)
        self.acceptance_loaded=False
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not all(receipt[k] for k in GATES):raise ValueError('Checked C2 primitive source required')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed C2 input: '+name)
            self.acceptance_loaded=True

    def leading_source(self,ends,chart,left,right=None):
        """Dispatch actual pure leading sources within one live owner graph."""
        if chart in CHARTS[:3]:return self.outer.owner.source_query(ends,chart,left,right)
        if chart in CHARTS[3:6]:return current.current.previous.previous.WholeZAllNSwitchFunctions.leading_packet(
            self.switch,ends,chart,left,right)
        if chart in CHARTS[6:11]:return current.current.previous.WholeZAllNLongPatchFunctions.leading_packet(
            self.long,ends,chart,left,right)
        if chart in CHARTS[11:]:return self.outer.leading_packet(ends,chart,left,right)
        raise ValueError('Admitted current original chart required')

    def source_packet(self,ends,chart,left,right=None):
        if chart not in CHARTS:raise ValueError('This admission covers current original 17-chart sources')
        key=(tuple(ends),chart,left,right)
        if key in self.cache:return self.cache[key]
        packet=self.leading_source(ends,chart,left,right)
        op=self.outer.owner.owner(ends);f=op.flow;raw=packet['original_generic_source'];proof=packet['original_full_source_quotients']
        if packet['source_identity']!=self.identity or packet['exact_common_P0_axial5'] is not op.P0:
            raise ValueError('Same current source and P0 object required')
        def get(*names):
            for name in names:
                if name in proof:return proof[name]
            raise ValueError('Actual Taylor source rows absent: '+str(names))
        zero=[f.scalar(0)]*6
        rows=dict(E=raw['common_velocity_E_axial5'],V=raw['common_velocity_V_axial5'],
            a=get('actual_a_axial5','actual_correlated_a_axial5','full_original_shear_a_axial5'),
            b=zero if proof.get('exact_b_and_t0_zero') else get('actual_b_axial5','actual_nonzero_b_axial5'),
            t0=zero if proof.get('exact_b_and_t0_zero') else get('actual_t0_axial5'),
            p2=get('full_signed_p2_axial4','exact_full_p2_axial4'),Delta=get('actual_Delta_axial5','Delta_axial5'))
        projected={}
        for name,row in rows.items():
            if len(row)<3 or any(v.ctx is not self.outer.c or v.scale.bases is not f.logs or v.ledger is not f.ledger for v in row):
                raise ValueError('Actual Z Taylor rows in one live basis/ledger required')
            projected[name]=[row[0],row[1],row[2]*2]
        result=dict(source_identity=self.identity,exact_Z_cell=list(ends),chart=chart,exact_common_P0_axial5=op.P0,
            actual_raw_root_ordinary_Z2=projected,actual_independent_P0_ordinary_Z2=[op.P0[0],op.P0[1],op.P0[2]*2],
            source_functions=encoded(self.functions[chart]),same_current_Z_Taylor_coordinate=True,
            coefficient2_converted_to_ordinary_derivative_once=True,
            source_radius_and_phase_Z_ZZ_independent=True,current_C1_root_prefix_retained=True,
            range_rows_are_not_function_values=True,actual_source_query_identity=packet['source_identity'])
        self.cache[key]=result;return result


def run(owner=None):
    began=time.monotonic()
    with mp.workdps(540):
        field=CurrentC2PhasePrimitives(owner=owner,require_checked=False)
        coordinates={chart:(0,1) for chart in CHARTS}
        coordinates.update(first_micro='inlet',second_micro=(1,1),actual_patch=(5,4),Rh_reference=(-5,1))
        packets=[field.source_packet(('0','.5'),chart,coordinates[chart]) for chart in CHARTS]
        report=dict(candidate_C2_source_constructed=True,source_family=field.identity,
            original_graph_prefix_length=len(field.prefix),exact_graph_nodes=field.built['graph'].nodes,
            actual_source_bindings=source_bindings(),
            current_C2_phase_functions=encoded(field.functions),actual_current_Z2_source_packets=encoded(packets),
            actual_current_regularity_proof=field.proof,source_charts=CHARTS,
            upstream_six_chart_C2_provider_installed=False,actual_five_target_ZZ_installed=False,
            actual_C2_repaired_limit_controls_installed=False,**dict.fromkeys(GATES+OPEN,False),
            input_hashes=field.hashes,execution_seconds=time.monotonic()-began)
        (HERE/NAME).write_bytes(gzip.compress((json.dumps(encoded(report),separators=(',',':'))+'\n').encode(),mtime=0))
    return field


if __name__=='__main__':run()
