"""Directed actual C2 source/inverse/density/five-target magnitude ranges.

All large quantities remain logarithmic upper bounds, not field values.
The original signed C2 graphs define the functions. Bounds use their exact
derivative identities and the current live source cells, including mixed
active/flat boxes and nonzero predecessor memories.
"""
from dataclasses import dataclass
from fractions import Fraction
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_C2_density_transport as target
import lei_ren_part1_paper_compliant_current_generic_shear_source_bounds as log_source

phase=target.phase;current=phase.current;outer=current.current
HERE,PREFIX,sha,ep=phase.HERE,phase.PREFIX,phase.sha,current.ep
LogUpper=log_source.LogUpper
NAME=PREFIX+'current_original_C2_target_ranges.json.gz'
RECEIPT=PREFIX+'current_original_C2_target_ranges_check.json'
GATES=('current_original_actual_17_chart_C2_primitive_density_ranges_installed',
    'current_original_actual_four_Z_C2_five_target_ranges_installed')
OPEN=phase.OPEN


@dataclass(frozen=True)
class JetBound:
    value: object
    Z: object
    ZZ: object


class Bounds:
    def __init__(self,c):self.c=c;self.zero=self.constant(0);self.one=self.constant(1)
    def constant(self,v):return LogUpper.constant(self.c,v)
    def row(self,v):
        if v.ctx is not self.c:raise ValueError('Same live source context required')
        return LogUpper(self.c,v.record()['log_absolute_upper'])
    def sum(self,*v):return LogUpper.add(self.c,v)
    def mul(self,*v):
        result=self.one
        for q in v:result=result*q
        return result
    def scale(self,v,k):return self.constant(k)*v
    def div(self,v,lower):return v.divide_positive(lower)
    def power(self,v,n):return self.mul(*([v]*n))
    def minimum(self,*values):
        if any(v.log is None for v in values):return self.zero
        return LogUpper(self.c,self.c.mpf(min(ep(v.log)[1] for v in values)))
    def fixed(self,v):return JetBound(v,self.zero,self.zero)
    def add(self,*rows):return JetBound(*(self.sum(*(getattr(q,k) for q in rows)) for k in ('value','Z','ZZ')))
    def product(self,a,b):
        return JetBound(a.value*b.value,self.sum(a.Z*b.value,a.value*b.Z),
            self.sum(a.ZZ*b.value,self.scale(a.Z*b.Z,2),a.value*b.ZZ))
    def scaled(self,a,v):return JetBound(*(v*getattr(a,k) for k in ('value','Z','ZZ')))
    def quotient(self,n,d,lower):
        value=self.div(n.value,lower)
        z=self.div(self.sum(n.Z,value*d.Z),lower)
        zz=self.div(self.sum(n.ZZ,self.scale(z*d.Z,2),value*d.ZZ),lower)
        return JetBound(value,z,zz)


def positive_lower(row):
    """Lower bound from a strictly positive enclosure of this same source."""
    lo,hi=ep(row.coefficient)
    if lo<=0:raise ValueError('Strict current positive source enclosure required')
    return row.ctx.mpf(ep(row.scale.evaluate()+row.ctx.ln(row.ctx.mpf(lo)))[0])


def regularity_bound_theorem():
    u=s.Symbol('u',real=True);r=u/s.sqrt(1+u*u);h=1/s.sqrt(1+u*u)
    assert s.simplify(s.diff(r,u)-(1+u*u)**s.Rational(-3,2))==0
    assert s.simplify(s.diff(r,u,2)+3*u*(1+u*u)**s.Rational(-5,2))==0
    assert s.simplify(s.diff(h,u,2)-(2*u*u-1)*(1+u*u)**s.Rational(-5,2))==0
    # For any real u, h<=1, |u|h^3<=1 and |u|h^5<=1;
    # |h_uu|<=(2u^2+1)/(1+u^2)^(5/2)<=3.
    return dict(passed=True,
        derivative_convention='ordinary Z derivatives; raw Taylor row2 converted exactly once upstream',
        original_cutoff='q=sigma(1-Delta/eta)*sqrt((2eta-Delta)/(2a)) on Delta<eta; 0 otherwise',
        active_constraints='a>0; kappa=a+b^2/a>=0; Delta=kappa-2<eta<=1/2; eta<=gamma=2eta-Delta<=3',
        sigma_global_first_bound=32,sigma_global_second_bound=1792,
        sigma_bound_proof='Original flat-tail theorem C_n*exp(4-m^-2)*m^(-3n), m=min(x,1-x)<=1/2. '
            'For n=1,2 the expression increases on (0,1/2]; C1=4,C2=28 give 32,1792. Outside [0,1] derivatives are zero.',
        flat_endpoint_second_derivatives_zero=True,body_branch_avoids_unnecessary_eta_reciprocal=True,
        original_signed_u_uniform_bounds='|r|<=1,h<=1; |r_Z|<=|u_Z|; '
            '|r_ZZ|<=3|u_Z|^2+|u_ZZ|; |h_Z|<=|u_Z|; |h_ZZ|<=3|u_Z|^2+|u_ZZ|',
        Poisson_denominator_lower='D>=1/(4*(1+u^2)^2)',
        Poisson_angle_derivative_bound='|w_psi|<=8/D_min^2',
        original_phase_normalization='int_0^(2pi)(1+t^2)=2pi*nu; nu=1+t0^2+2q^2>=1',
        actual_inverse_Jacobian_lower='Phi_psi>=1/(2pi*nu_upper)>0',
        original_A_uniform_bound='Flat A=0; active a<2+eta<=5/2 and |phi-psi/(2pi)|<=1 imply |A|<=5/4',
        exponential_bound='N>=N0>=160; |A|<=5/4 => |exp(A/N)|,|exprel(A/N)|<=exp(5/(4N0))',
        source_owned_collar_and_quiet_power_are_exact_flat=True,
        bounds_are_on_original_functions_and_not_point_values=True,
        current_C2_function_regularity_proof=phase.regularity_proof(),
        original_sigma_source_bindings=phase.source_bindings()['original_flat_sigma_jets'])


def primitive_bounds(field,ends,chart,left,right):
    c=field.c;op=field.phase.outer.owner.owner(ends);f=op.flow;bd=Bounds(c)
    packet=field.phase.source_packet(ends,chart,left,right)
    original=field.phase.leading_source(ends,chart,left,right)
    raw=packet['actual_raw_root_ordinary_Z2']
    rows={k:JetBound(*(bd.row(v) for v in values)) for k,values in raw.items()}
    zero=bd.fixed(bd.zero)
    eta_log=field.phase.outer.owner.source.eta_log;dstar_log=field.phase.outer.owner.source.dstar_log
    if ep(eta_log)[1]>ep(c.ln(c.mpf('.5')))[0]:raise ValueError('Original eta<=1/2 required')
    a_lower=positive_lower(raw['a'][0]);eta=f.factor((0,0,0,0,0),eta_log)
    difference=raw['Delta'][0]-eta
    flat=chart=='O3_power' or original.get('source_owned_collar_q_C0_Z_exact_zero',False) or ep(difference.coefficient)[0]>=0
    proof=dict(exact_same_source_P0=original['exact_common_P0_axial5'] is op.P0,
        actual_positive_a_log_lower=a_lower,actual_eta_log=eta_log,actual_dstar_log=dstar_log,
        actual_lazy_branch_difference=difference.record(),source_owned_collar_flat=original.get('source_owned_collar_q_C0_Z_exact_zero',False),
        original_quiet_power_admission=field.phase.windows[chart].get('quiet_source_proof'),
        actual_C2_function_handles=phase.encoded(field.phase.functions[chart]),
        actual_source_ordinary_Z2_rows={k:[v.record() for v in values] for k,values in raw.items()},
        actual_independent_P0_rows=[v.record() for v in packet['actual_independent_P0_ordinary_Z2']],
        source_context_basis_and_ledger_retained=True,full_original_phase_box=[0,1])
    if not proof['exact_same_source_P0']:raise ValueError('Same current P0 required')
    if flat:
        proof.update(branch='exact_original_flat',all_original_A_B_second_rows_exact_zero=True)
        return dict(A=zero,B=zero,E=rows['E'],V=rows['V'],q=zero,proof=proof)
    a,b,E,V,t0,D,p2=[rows[k] for k in ('a','b','E','V','t0','Delta','p2')]
    # Body has sigma=1. Mixed boxes use an active branch with gamma>=eta,
    # then include exact flat zero; no positive gamma theorem on the flat part.
    body=ep(raw['Delta'][0].coefficient)[1]<=0
    gamma=eta*2-raw['Delta'][0]
    gamma_lower=positive_lower(gamma) if body else c.mpf(ep(eta_log)[0])
    root0=LogUpper(c,(c.ln(3)-c.ln(2)-a_lower)/2)
    L1=bd.scale(bd.sum(bd.div(D.Z,gamma_lower),bd.div(a.Z,a_lower)),c.mpf('.5'))
    L2=bd.scale(bd.sum(bd.div(D.ZZ,gamma_lower),bd.div(bd.power(D.Z,2),2*gamma_lower),
        bd.div(a.ZZ,a_lower),bd.div(bd.power(a.Z,2),2*a_lower)),c.mpf('.5'))
    sig1=bd.zero if body else bd.scale(bd.div(D.Z,eta_log),32)
    sig2=bd.zero if body else bd.sum(bd.scale(bd.div(bd.power(D.Z,2),2*eta_log),1792),bd.scale(bd.div(D.ZZ,eta_log),32))
    q=JetBound(root0,root0*bd.sum(sig1,L1),
        root0*bd.sum(sig2,bd.scale(sig1*L1,2),L2,bd.power(L1,2)))
    u=bd.quotient(bd.product(p2,q),bd.fixed(LogUpper(c,dstar_log)),dstar_log)
    curvature=bd.sum(bd.scale(bd.power(u.Z,2),3),u.ZZ)
    h=JetBound(bd.one,u.Z,curvature);r=JetBound(bd.one,u.Z,curvature)
    alpha=bd.scaled(bd.product(q,h),bd.constant(2))
    invD=bd.scale(bd.power(bd.sum(bd.one,bd.power(u.value,2)),2),4)
    den1=bd.scale(r.Z,4);den2=bd.sum(bd.scale(bd.power(r.Z,2),2),bd.scale(r.ZZ,4))
    w0=bd.scale(invD,2)
    w1=bd.sum(r.Z,w0*den1)*invD
    w2=bd.sum(r.ZZ,bd.scale(w1*den1,2),w0*den2)*invD
    w=JetBound(w0,w1,w2);t=bd.add(t0,bd.product(alpha,w));tpsi=alpha.value*bd.scale(bd.power(invD,2),8)
    nu=bd.add(bd.fixed(bd.one),bd.product(t0,t0),bd.scaled(bd.product(q,q),bd.constant(2)))
    two_pi=bd.constant(2*c.pi);inv_two_pi=bd.constant(1/(2*c.pi))
    K=JetBound(inv_two_pi,inv_two_pi*nu.Z,inv_two_pi*bd.sum(nu.ZZ,bd.scale(bd.power(nu.Z,2),2)))
    P=JetBound(two_pi*nu.value,bd.scale(two_pi*t.value*t.Z,2),
        bd.scale(two_pi*bd.sum(bd.power(t.Z,2),t.value*t.ZZ),2))
    Phi=bd.product(K,P)
    inverse_J=two_pi*nu.value
    PhiPsiZ=bd.sum(K.Z*bd.sum(bd.one,bd.power(t.value,2)),bd.scale(K.value*t.value*t.Z,2))
    PhiPsiPsi=bd.scale(K.value*t.value*tpsi,2)
    psi1=inverse_J*Phi.Z
    psi2=inverse_J*bd.sum(Phi.ZZ,bd.scale(PhiPsiZ*psi1,2),PhiPsiPsi*bd.power(psi1,2))
    T1=JetBound(two_pi*t.value,bd.sum(two_pi*t.Z,t.value*psi1),
        bd.sum(two_pi*t.ZZ,bd.scale(t.Z*psi1,2),tpsi*bd.power(psi1,2),t.value*psi2))
    active_a=JetBound(bd.minimum(a.value,bd.constant(c.mpf('2.5'))),a.Z,a.ZZ)
    M=bd.add(bd.scaled(bd.product(active_a,T1),inv_two_pi),b)
    A=JetBound(bd.constant(c.mpf('1.25')),
        bd.scale(bd.sum(a.Z,active_a.value*psi1*inv_two_pi),c.mpf('.5')),
        bd.scale(bd.sum(a.ZZ,bd.scale(a.Z*psi1*inv_two_pi,2),active_a.value*psi2*inv_two_pi),c.mpf('.5')))
    B=bd.scaled(bd.product(E,M),bd.constant(c.mpf('.5')))
    proof.update(branch='body_sigma_one' if body else 'active_flat_union',
        active_gamma_positive_log_lower=gamma_lower,active_gamma_upper=3,
        q_root_upper=root0.record(),q_bounds=record(q),same_original_inverse_psi_Z_upper=psi1.record(),
        same_original_inverse_psi_ZZ_upper=psi2.record(),actual_inverse_Jacobian_reciprocal_upper=inverse_J.record(),
        original_Poisson_denominator_reciprocal_upper=invD.record(),
        flat_part_included_as_exact_zero=True,caps_not_differentiated=True)
    return dict(A=A,B=B,E=E,V=V,q=q,proof=proof)


def density_bounds(c,primitive,N0):
    if type(N0) is not int or N0<160 or N0.bit_length()>4096:raise ValueError('Admitted exact leading-input N0 required')
    bd=Bounds(c);A,B,E,V=[primitive[k] for k in ('A','B','E','V')]
    factor=bd.constant(c.exp(c.mpf('1.25')/N0));epsilon=bd.constant(c.mpf(1)/N0)
    F=JetBound(E.value*A.value*factor,
        bd.sum(E.Z*A.value,E.value*A.Z)*factor,
        bd.sum(E.ZZ*A.value,bd.scale(E.Z*A.Z,2),E.value*bd.sum(A.ZZ,bd.power(A.Z,2)*epsilon))*factor)
    EF,VF,EB,VB=[bd.product(a,b) for a,b in ((E,F),(V,F),(E,B),(V,B))]
    F2,B2=bd.product(F,F),bd.product(B,B);zero=bd.fixed(bd.zero)
    first=dict(m=B,h=F,k=bd.add(VF,EB),e=bd.add(bd.scaled(VB,bd.constant(2)),EF),p=EF)
    second=dict(m=zero,h=zero,k=bd.product(F,B),e=bd.add(B2,bd.scaled(F2,bd.constant(c.mpf('.5')))),p=bd.scaled(F2,bd.constant(c.mpf('.5'))))
    total={key:bd.add(first[key],bd.scaled(second[key],epsilon)) for key in current.RATES}
    return dict(first=first,second=second,N_scaled=total,F=F,
        full_exp_exprel_factor_upper=factor.record(),epsilon_upper=epsilon.record(),
        signed_coefficient_functions_enclosed_by_magnitude_bounds=True)


def record(value):
    if isinstance(value,JetBound):return [value.value.record(),value.Z.record(),value.ZZ.record()]
    if isinstance(value,LogUpper):return value.record()
    if isinstance(value,dict):return {k:record(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [record(v) for v in value]
    return value


class CurrentC2TargetRanges:
    def __init__(self,target_field=None,owner=None,require_checked=True):
        self.target=target_field if target_field is not None else target.CurrentC2DensityTransport(owner=owner)
        if not self.target.acceptance_loaded:raise ValueError('Accepted actual C2 five-target function source required')
        self.phase=self.target.phase;self.c=self.phase.outer.c;self.identity=self.target.identity
        self.hashes=dict(self.target.hashes)
        for name in (target.RECEIPT,target.NAME,Path(__file__).name,Path(log_source.__file__).name):self.hashes[name]=sha(name)
        self.cache={};self.proof=regularity_bound_theorem();self.acceptance_loaded=False
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not all(receipt[k] for k in GATES):raise ValueError('Checked actual C2 target ranges required')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed C2 range source '+name)
            self.acceptance_loaded=True

    def weights(self,ends,chart,left,right,rate):
        owner=self.phase.outer
        if chart in phase.CHARTS[:3]:
            mass,decay,suffix=owner.owner.weights(ends,chart,left,right,rate)
            return dict(original_mass=mass,cell_decay=decay,suffix_decay=suffix)
        if chart in phase.CHARTS[3:6]:
            native=owner.owner.switch_owner(ends)
            mass,decay,suffix=outer.previous.previous.original.own_weights(native.first,chart,left,right,rate)
            return dict(original_mass=mass,cell_decay=decay,suffix_decay=suffix)
        if chart in phase.CHARTS[6:11]:return outer.previous.WholeZAllNLongPatchFunctions.weights(owner,ends,chart,left,right,rate)
        return owner.weights(ends,chart,left,right,rate)

    def cell(self,ends,chart,left,right):
        key=(tuple(ends),chart,left,right)
        if key not in self.cache:
            primitive=primitive_bounds(self,tuple(ends),chart,left,right)
            density=density_bounds(self.c,primitive,self.phase.outer.N0)
            self.cache[key]=dict(primitive=primitive,density=density)
        return self.cache[key]

    def transport(self,ends):
        ends=tuple(ends)
        if ends not in outer.CELLS:raise ValueError('Original admitted Z cell required')
        c=self.c;bd=Bounds(c);history={key:bd.fixed(bd.zero) for key in current.RATES};windows=[]
        partitions=current.native_partitions()
        for chart in phase.CHARTS:
            partition=partitions[chart];local={key:bd.fixed(bd.zero) for key in current.RATES};cells=[]
            for left,right in zip(partition,partition[1:]):
                query=self.cell(ends,chart,left,right);weights={};contributions={}
                for key,rate in current.RATES.items():
                    weight=self.weights(ends,chart,left,right,rate)
                    mass=bd.row(weight['original_mass']);suffix=bd.row(weight['suffix_decay'])
                    contribution=bd.scaled(query['density']['N_scaled'][key],mass*suffix)
                    local[key]=bd.add(local[key],contribution);contributions[key]=contribution
                    weights[key]={k:v.record() for k,v in weight.items() if k in ('original_mass','cell_decay','suffix_decay')}
                cells.append(dict(exact_left=left,exact_right=right,actual_source_and_primitive_range=record(query['primitive']),
                    actual_density_ranges=record(query['density']),original_positive_own_rate_weights=weights,
                    actual_N_scaled_contribution_C2=record(contributions),physical_Jacobian_applied_once=True))
            before=history;memories={key:self.weights(ends,chart,partition[0],partition[-1],rate)['cell_decay'] for key,rate in current.RATES.items()}
            history={key:bd.add(bd.scaled(before[key],bd.row(memories[key])),local[key]) for key in current.RATES}
            windows.append(dict(chart=chart,cells=cells,actual_N_scaled_incoming_C2=record(before),
                original_own_rate_memory={k:v.record() for k,v in memories.items()},
                actual_N_scaled_local_C2=record(local),actual_N_scaled_outgoing_C2=record(history),
                second_row_predecessor_never_reset=True,global_phase_cover=[0,1]))
            print('Actual C2 magnitude transport',ends,chart,flush=True)
        endpoint=self.phase.source_packet(ends,'O3_power',(2,1));amplitude=endpoint['actual_raw_root_ordinary_Z2']['E']
        A=JetBound(*(bd.row(v) for v in amplitude));logA=positive_lower(amplitude[0]);A2=bd.product(A,A)
        mu=LogUpper(c,self.phase.outer.logmu);targets={}
        for row,key,degree in (('M','m',1),('I','h',1),('S','e',2),('Cp','p',2)):
            targets[row]=bd.quotient(history[key],A if degree==1 else A2,degree*logA)
        joint=bd.add(history['k'],bd.product(A,history['m']))
        targets[current.controls.ROWS[1]]=bd.quotient(joint,bd.scaled(A2,mu),2*logA+self.phase.outer.logmu)
        targets={key:targets[key] for key in current.controls.ROWS}
        maximum=LogUpper(c,c.mpf(max(ep(q.ZZ.log)[1] for q in targets.values() if q.ZZ.log is not None)))
        return dict(source_identity=self.identity,exact_Z_cell=ends,actual_17_chart_C2_transports=windows,
            actual_normalized_N_scaled_histories_C2=record(history),actual_terminal_amplitude_C2=[v.record() for v in amplitude],
            actual_positive_terminal_amplitude_log_lower=logA,actual_positive_mu_log=self.phase.outer.logmu,
            actual_N_scaled_five_target_C2_ranges=record(targets),actual_joint_angular_numerator_C2_range=record(joint),
            actual_five_target_ZZ_maximum=maximum.record(),actual_five_target_C2_function_handles=phase.encoded(self.target.target_functions()),
            original_independent_P0_in_source_and_not_added_twice=True,
            uniform_for_all_integer_N_ge_N0=True,epsilon_domain='[0,1/N0]',
            range_caps_not_point_function_values=True,actual_C2_repaired_limit_controls_installed=False)


def run(target_field=None,owner=None,field=None):
    began=time.monotonic()
    with mp.workdps(540):
        field=field if field is not None else CurrentC2TargetRanges(target_field=target_field,owner=owner,require_checked=False)
        field.hashes[Path(__file__).name]=sha(Path(__file__).name)
        rows=[field.transport(ends) for ends in outer.CELLS]
        report=dict(candidate_actual_C2_target_ranges_constructed=True,source_family=field.identity,
            actual_current_regular_majorant_theorem=field.proof,actual_four_Z_C2_target_transports=rows,
            actual_C2_function_source_bindings=dict(primitive_bounds=current.ast_binding(primitive_bounds),
                density_bounds=current.ast_binding(density_bounds),actual_source_projection=current.ast_binding(phase.CurrentC2PhasePrimitives.source_packet),
                original_log_bound_arithmetic=current.ast_binding(LogUpper),
                original_phase_density_transport_source=target.RECEIPT,original_phase_density_transport_sha256=sha(target.RECEIPT)),
            actual_C2_repaired_limit_controls_installed=False,**dict.fromkeys(GATES+OPEN,False),
            input_hashes=field.hashes,execution_seconds=time.monotonic()-began)
        (HERE/NAME).write_bytes(gzip.compress((json.dumps(phase.encoded(record(report)),separators=(',',':'))+'\n').encode(),mtime=0))
    return field


if __name__=='__main__':run()
