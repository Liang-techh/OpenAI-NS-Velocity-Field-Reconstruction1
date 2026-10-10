"""Directed same-source ordinary C3 primitive/density/five-target bounds.

Magnitude caps stay logarithmic. Signed functions are the separately accepted
C3 source/inverse/transport graph, never interval endpoints or cap derivatives.
"""
from dataclasses import dataclass
import gzip
import json
import math
from pathlib import Path
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_C3_source_targets as target
import lei_ren_part1_paper_compliant_current_original_C2_target_ranges as lower

phase,current,outer=lower.phase,lower.current,lower.outer
HERE,PREFIX,sha,ep=lower.HERE,lower.PREFIX,lower.sha,lower.ep
LogUpper,positive_lower=lower.LogUpper,lower.positive_lower
NAME=PREFIX+'current_original_C3_target_ranges.json.gz'
RECEIPT=PREFIX+'current_original_C3_target_ranges_check.json'
GATES=('current_original_actual_17_chart_C3_primitive_density_ranges_installed',
    'current_original_actual_four_Z_C3_five_target_ranges_installed')
OPEN=lower.OPEN
KEYS=('value','Z','ZZ','ZZZ')


@dataclass(frozen=True)
class JetBound:
    value: object
    Z: object
    ZZ: object
    ZZZ: object


def rows(q):return tuple(getattr(q,k) for k in KEYS)


class Bounds(lower.Bounds):
    def fixed(self,v):return JetBound(v,self.zero,self.zero,self.zero)
    def add(self,*values):return JetBound(*(self.sum(*(getattr(q,k) for q in values)) for k in KEYS))
    def product(self,a,b):
        aa,bb=rows(a),rows(b)
        return JetBound(*(self.sum(*(self.scale(aa[k]*bb[j-k],math.comb(j,k)) for k in range(j+1))) for j in range(4)))
    def scaled(self,a,v):return JetBound(*(v*q for q in rows(a)))
    def quotient(self,n,d,log_lower):
        return self.reciprocal_quotient(n,d,LogUpper(self.c,-log_lower))
    def reciprocal_quotient(self,n,d,inverse_upper):
        nn,dd=rows(n),rows(d);out=[]
        for j in range(4):
            out.append(self.sum(nn[j],*(self.scale(dd[k]*out[j-k],math.comb(j,k)) for k in range(1,j+1)))*inverse_upper)
        return JetBound(*out)


def sigma_third_cap(c):return phase.flat_source.sigma_tail_bound(c,3,c.mpf('.5'))


def regularity_bound_theorem():
    u=s.Symbol('u',real=True);r=u/s.sqrt(1+u*u);h=1/s.sqrt(1+u*u)
    assert s.simplify(s.diff(r,u,3)-(12*u*u-3)/(1+u*u)**s.Rational(7,2))==0
    assert s.simplify(s.diff(h,u,3)-(9*u-6*u**3)/(1+u*u)**s.Rational(7,2))==0
    return dict(passed=True,ordinary_derivative_convention='Z3=6*actual Taylor coefficient[3] once in the source provider',
        inherited_positive_source_theorem=lower.regularity_bound_theorem(),
        signed_u_third_uniform_bound=15,
        signed_u_third_chain_bound='15*|u_Z|^3+9*|u_Z*u_ZZ|+|u_ZZZ|',
        signed_u_bound_proof='|12u^2-3|/(1+u^2)^(7/2)<=15; '
            '|9u-6u^3|/(1+u^2)^(7/2)<=15 using |u|/sqrt(1+u^2)<=1.',
        sigma_third_source_envelope='256*exp(4-m^-2)*m^-9, 0<m<=1/2',
        sigma_third_global_cap='256*exp(-1/2)*(3/sqrt(2))^9; peak m=sqrt(2)/3',
        sigma_third_provider=current.ast_binding(sigma_third_cap),
        original_sigma_tail_provider=current.ast_binding(phase.flat_source.sigma_tail_bound),
        inverse_third_formula='J^-1*(Phi_ZZZ+3Phi_psiZZ*psi_Z+3Phi_psipsiZ*psi_Z^2+'
            'Phi_psipsipsi*psi_Z^3+3*(Phi_psiZ+Phi_psipsi*psi_Z)*psi_ZZ)',
        fixed_angle_cross_partials_and_composed_T1_terms_retained=True,
        original_flat_third_rows_exact_zero=True,
        genuine_signed_C3_function_source=target.RECEIPT,
        genuine_signed_C3_function_receipt_sha256=sha(target.RECEIPT),
        bounds_are_not_function_values_or_derivatives_of_caps=True)


def primitive_bounds(field,ends,chart,left,right):
    c=field.c;bd=Bounds(c);op=field.phase.outer.owner.owner(ends);f=op.flow
    packet=field.target.source_packet(ends,chart,left,right)
    original=field.phase.leading_source(ends,chart,left,right)
    raw=packet['actual_raw_root_ordinary_Z3'];source_rows={k:JetBound(*(bd.row(q) for q in values)) for k,values in raw.items()}
    zero=bd.fixed(bd.zero);eta_log=field.phase.outer.owner.source.eta_log;dstar_log=field.phase.outer.owner.source.dstar_log
    if ep(eta_log)[1]>ep(c.ln(c.mpf('.5')))[0]:raise ValueError('Original eta<=1/2 required')
    a_lower=positive_lower(raw['a'][0]);eta=f.factor((0,0,0,0,0),eta_log);difference=raw['Delta'][0]-eta
    flat=chart=='O3_power' or original.get('source_owned_collar_q_C0_Z_exact_zero',False) or ep(difference.coefficient)[0]>=0
    proof=dict(exact_same_source_P0=original['exact_common_P0_axial5'] is op.P0,
        actual_positive_a_log_lower=a_lower,actual_eta_log=eta_log,actual_dstar_log=dstar_log,
        actual_lazy_branch_difference=difference.record(),source_owned_collar_flat=original.get('source_owned_collar_q_C0_Z_exact_zero',False),
        original_quiet_power_admission=field.phase.windows[chart].get('quiet_source_proof'),
        actual_C3_function_handles=target.encoded(field.target.primitives[chart]),
        actual_source_ordinary_Z3_rows={k:[v.record() for v in values] for k,values in raw.items()},
        actual_independent_P0_rows=[v.record() for v in packet['actual_independent_P0_ordinary_Z3']],
        source_context_basis_and_ledger_retained=True,full_original_phase_box=[0,1],
        live_source_phase_Z_exact_zero=original['actual_phase_Z_exact_zero'])
    if not proof['exact_same_source_P0'] or not proof['live_source_phase_Z_exact_zero']:raise ValueError('Same source P0/phase required')
    if flat:
        proof.update(branch='exact_original_flat',all_original_A_B_third_rows_exact_zero=True)
        return dict(A=zero,B=zero,E=source_rows['E'],V=source_rows['V'],q=zero,proof=proof)
    a,b,E,V,t0,D,p2=[source_rows[k] for k in ('a','b','E','V','t0','Delta','p2')]
    body=ep(raw['Delta'][0].coefficient)[1]<=0
    gamma=eta*2-raw['Delta'][0];gamma_lower=positive_lower(gamma) if body else c.mpf(ep(eta_log)[0])
    root0=LogUpper(c,(c.ln(3)-c.ln(2)-a_lower)/2)
    L1=bd.scale(bd.sum(bd.div(D.Z,gamma_lower),bd.div(a.Z,a_lower)),c.mpf('.5'))
    L2=bd.scale(bd.sum(bd.div(D.ZZ,gamma_lower),bd.div(bd.power(D.Z,2),2*gamma_lower),
        bd.div(a.ZZ,a_lower),bd.div(bd.power(a.Z,2),2*a_lower)),c.mpf('.5'))
    L3=bd.scale(bd.sum(bd.div(D.ZZZ,gamma_lower),bd.scale(bd.div(D.Z*D.ZZ,2*gamma_lower),3),
        bd.scale(bd.div(bd.power(D.Z,3),3*gamma_lower),2),bd.div(a.ZZZ,a_lower),
        bd.scale(bd.div(a.Z*a.ZZ,2*a_lower),3),bd.scale(bd.div(bd.power(a.Z,3),3*a_lower),2)),c.mpf('.5'))
    root=JetBound(root0,root0*L1,root0*bd.sum(L2,bd.power(L1,2)),
        root0*bd.sum(L3,bd.scale(L1*L2,3),bd.power(L1,3)))
    cut1,cut2,cut3=[bd.div(q,eta_log) for q in (D.Z,D.ZZ,D.ZZZ)]
    sig1=bd.zero if body else bd.scale(cut1,32)
    sig2=bd.zero if body else bd.sum(bd.scale(bd.power(cut1,2),1792),bd.scale(cut2,32))
    sig3=bd.zero if body else bd.sum(bd.constant(sigma_third_cap(c))*bd.power(cut1,3),
        bd.scale(cut1*cut2,3*1792),bd.scale(cut3,32))
    q=bd.product(JetBound(bd.one,sig1,sig2,sig3),root)
    u=bd.quotient(bd.product(p2,q),bd.fixed(LogUpper(c,dstar_log)),dstar_log)
    curvature=bd.sum(bd.scale(bd.power(u.Z,2),3),u.ZZ)
    cubic=bd.sum(bd.scale(bd.power(u.Z,3),15),bd.scale(u.Z*u.ZZ,9),u.ZZZ)
    h=JetBound(bd.one,u.Z,curvature,cubic);r=JetBound(bd.one,u.Z,curvature,cubic)
    alpha=bd.scaled(bd.product(q,h),bd.constant(2))
    invD=bd.scale(bd.power(bd.sum(bd.one,bd.power(u.value,2)),2),4)
    den1=bd.scale(r.Z,4);den2=bd.sum(bd.scale(bd.power(r.Z,2),2),bd.scale(r.ZZ,4))
    den3=bd.sum(bd.scale(r.Z*r.ZZ,6),bd.scale(r.ZZZ,4))
    den=JetBound(bd.constant(4),den1,den2,den3)
    w0=bd.scale(invD,2);w1=bd.sum(r.Z,w0*den1)*invD
    w2=bd.sum(r.ZZ,bd.scale(w1*den1,2),w0*den2)*invD
    w3=bd.sum(r.ZZZ,bd.scale(w2*den1,3),bd.scale(w1*den2,3),w0*den3)*invD
    w=JetBound(w0,w1,w2,w3)
    # |D_psi| and |D_psipsi| have the same Z-row bounds 2|r_j|.
    dpsi=bd.scaled(r,bd.constant(2))
    wpsi=bd.reciprocal_quotient(bd.add(bd.fixed(bd.one),bd.product(w,dpsi)),den,invD)
    wpsipsi=bd.reciprocal_quotient(bd.add(bd.fixed(bd.one),
        bd.scaled(bd.product(wpsi,dpsi),bd.constant(2)),bd.product(w,dpsi)),den,invD)
    t=bd.add(t0,bd.product(alpha,w));tpsi=bd.product(alpha,wpsi);tpsipsi=bd.product(alpha,wpsipsi)
    nu=bd.add(bd.fixed(bd.one),bd.product(t0,t0),bd.scaled(bd.product(q,q),bd.constant(2)))
    two_pi=bd.constant(2*c.pi);inv_two_pi=bd.constant(1/(2*c.pi))
    K=bd.quotient(bd.fixed(inv_two_pi),nu,c.mpf(0))
    tt=bd.product(t,t);P=JetBound(two_pi*nu.value,*(two_pi*v for v in rows(tt)[1:]))
    Phi=bd.product(K,P);PhiPsi=bd.product(K,bd.add(bd.fixed(bd.one),tt))
    PhiPsiPsi=bd.scaled(bd.product(K,bd.product(t,tpsi)),bd.constant(2))
    PhiPsiPsiPsi=bd.scale(K.value*bd.sum(bd.power(tpsi.value,2),t.value*tpsipsi.value),2)
    inverse_J=two_pi*nu.value;psi1=inverse_J*Phi.Z
    psi2=inverse_J*bd.sum(Phi.ZZ,bd.scale(PhiPsi.Z*psi1,2),PhiPsiPsi.value*bd.power(psi1,2))
    psi3=inverse_J*bd.sum(Phi.ZZZ,bd.scale(PhiPsi.ZZ*psi1,3),
        bd.scale(PhiPsiPsi.Z*bd.power(psi1,2),3),PhiPsiPsiPsi*bd.power(psi1,3),
        bd.scale(bd.sum(PhiPsi.Z,PhiPsiPsi.value*psi1)*psi2,3))
    T1=JetBound(two_pi*t.value,bd.sum(two_pi*t.Z,t.value*psi1),
        bd.sum(two_pi*t.ZZ,bd.scale(t.Z*psi1,2),tpsi.value*bd.power(psi1,2),t.value*psi2),
        bd.sum(two_pi*t.ZZZ,bd.scale(t.ZZ*psi1,3),bd.scale(tpsi.Z*bd.power(psi1,2),3),
            tpsipsi.value*bd.power(psi1,3),bd.scale(t.Z*psi2,3),bd.scale(tpsi.value*psi1*psi2,3),t.value*psi3))
    active_a=JetBound(bd.minimum(a.value,bd.constant(c.mpf('2.5'))),a.Z,a.ZZ,a.ZZZ)
    M=bd.add(bd.scaled(bd.product(active_a,T1),inv_two_pi),b)
    chi=JetBound(bd.one,psi1*inv_two_pi,psi2*inv_two_pi,psi3*inv_two_pi)
    AA=bd.scaled(bd.product(active_a,chi),bd.constant(c.mpf('.5')))
    A=JetBound(bd.constant(c.mpf('1.25')),AA.Z,AA.ZZ,AA.ZZZ)
    B=bd.scaled(bd.product(E,M),bd.constant(c.mpf('.5')))
    proof.update(branch='body_sigma_one' if body else 'active_flat_union',
        active_gamma_positive_log_lower=gamma_lower,active_gamma_upper=3,q_root_upper=root0.record(),
        root_log_derivative_bounds=[v.record() for v in (L1,L2,L3)],sigma_third_global_upper=bd.constant(sigma_third_cap(c)).record(),
        q_bounds=record(q),same_original_inverse_psi_C3_upper=[v.record() for v in (psi1,psi2,psi3)],
        actual_inverse_Jacobian_reciprocal_upper=inverse_J.record(),original_Poisson_denominator_reciprocal_upper=invD.record(),
        fixed_angle_t_C3_upper=record(t),fixed_angle_tpsi_C3_upper=record(tpsi),fixed_angle_tpsipsi_C3_upper=record(tpsipsi),
        fixed_angle_Phi_C3_upper=record(Phi),fixed_angle_PhiPsi_C3_upper=record(PhiPsi),
        fixed_angle_PhiPsiPsi_C3_upper=record(PhiPsiPsi),fixed_angle_PhiPsiPsiPsi_upper=PhiPsiPsiPsi.record(),
        composed_T1_C3_upper=record(T1),flat_part_included_as_exact_zero=True,caps_not_differentiated=True)
    return dict(A=A,B=B,E=E,V=V,q=q,proof=proof)


def density_bounds(c,primitive,N0):
    if type(N0) is not int or N0<160 or N0.bit_length()>4096:raise ValueError('Exact leading-input N0 required')
    bd=Bounds(c);A,B,E,V=[primitive[k] for k in ('A','B','E','V')]
    factor=bd.constant(c.exp(c.mpf('1.25')/N0));epsilon=bd.constant(c.mpf(1)/N0)
    F=JetBound(E.value*A.value*factor,bd.sum(E.Z*A.value,E.value*A.Z)*factor,
        bd.sum(E.ZZ*A.value,bd.scale(E.Z*A.Z,2),E.value*bd.sum(A.ZZ,bd.power(A.Z,2)*epsilon))*factor,
        bd.sum(E.ZZZ*A.value,bd.scale(E.ZZ*A.Z,3),bd.scale(E.Z*bd.sum(A.ZZ,bd.power(A.Z,2)*epsilon),3),
            E.value*bd.sum(A.ZZZ,bd.scale(A.Z*A.ZZ*epsilon,3),bd.power(A.Z,3)*bd.power(epsilon,2)))*factor)
    EF,VF,EB,VB=[bd.product(a,b) for a,b in ((E,F),(V,F),(E,B),(V,B))]
    F2,B2=bd.product(F,F),bd.product(B,B);zero=bd.fixed(bd.zero)
    first=dict(m=B,h=F,k=bd.add(VF,EB),e=bd.add(bd.scaled(VB,bd.constant(2)),EF),p=EF)
    second=dict(m=zero,h=zero,k=bd.product(F,B),e=bd.add(B2,bd.scaled(F2,bd.constant(c.mpf('.5')))),p=bd.scaled(F2,bd.constant(c.mpf('.5'))))
    total={key:bd.add(first[key],bd.scaled(second[key],epsilon)) for key in current.RATES}
    return dict(first=first,second=second,N_scaled=total,F=F,full_exp_exprel_factor_upper=factor.record(),
        epsilon_upper=epsilon.record(),signed_coefficient_functions_enclosed_by_magnitude_bounds=True)


def record(value):
    if isinstance(value,JetBound):return [q.record() for q in rows(value)]
    if isinstance(value,LogUpper):return value.record()
    if isinstance(value,dict):return {k:record(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [record(v) for v in value]
    return value


class CurrentC3TargetRanges:
    def __init__(self,target_field=None,owner=None,require_checked=True):
        self.target=target_field if target_field is not None else target.CurrentC3SourceTargets(owner=owner)
        if not self.target.acceptance_loaded:raise ValueError('Accepted signed C3 five-target functions required')
        self.phase=self.target.phase;self.c=self.phase.outer.c;self.identity=self.target.identity
        self.hashes=dict(self.target.hashes)
        for name in (target.NAME,target.RECEIPT,Path(__file__).name,Path(lower.__file__).name,lower.RECEIPT):self.hashes[name]=sha(name)
        self.cache={};self.proof=regularity_bound_theorem();self.acceptance_loaded=False
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not all(receipt[k] for k in GATES):raise ValueError('Checked C3 ranges required')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed C3 range source '+name)
            self.acceptance_loaded=True

    def weights(self,ends,chart,left,right,rate):return lower.CurrentC2TargetRanges.weights(self,ends,chart,left,right,rate)

    def cell(self,ends,chart,left,right):
        key=(tuple(ends),chart,left,right)
        if key not in self.cache:
            primitive=primitive_bounds(self,tuple(ends),chart,left,right)
            self.cache[key]=dict(primitive=primitive,density=density_bounds(self.c,primitive,self.phase.outer.N0))
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
                    weight=self.weights(ends,chart,left,right,rate);mass=bd.row(weight['original_mass']);suffix=bd.row(weight['suffix_decay'])
                    contribution=bd.scaled(query['density']['N_scaled'][key],mass*suffix)
                    local[key]=bd.add(local[key],contribution);contributions[key]=contribution
                    weights[key]={k:v.record() for k,v in weight.items() if k in ('original_mass','cell_decay','suffix_decay')}
                cells.append(dict(exact_left=left,exact_right=right,actual_source_and_primitive_range=record(query['primitive']),
                    actual_density_ranges=record(query['density']),original_positive_own_rate_weights=weights,
                    actual_N_scaled_contribution_C3=record(contributions),physical_Jacobian_applied_once=True))
            before=history;memories={key:self.weights(ends,chart,partition[0],partition[-1],rate)['cell_decay'] for key,rate in current.RATES.items()}
            history={key:bd.add(bd.scaled(before[key],bd.row(memories[key])),local[key]) for key in current.RATES}
            windows.append(dict(chart=chart,cells=cells,actual_N_scaled_incoming_C3=record(before),
                original_own_rate_memory={k:v.record() for k,v in memories.items()},
                actual_N_scaled_local_C3=record(local),actual_N_scaled_outgoing_C3=record(history),
                third_row_predecessor_never_reset=True,global_phase_cover=[0,1]))
            print('Actual C3 magnitude transport',ends,chart,flush=True)
        endpoint=self.target.source_packet(ends,'O3_power',(2,1));amplitude=endpoint['actual_raw_root_ordinary_Z3']['E']
        A=JetBound(*(bd.row(v) for v in amplitude));logA=positive_lower(amplitude[0]);A2=bd.product(A,A)
        mu=LogUpper(c,self.phase.outer.logmu);targets={}
        for row,key,degree in (('M','m',1),('I','h',1),('S','e',2),('Cp','p',2)):
            targets[row]=bd.quotient(history[key],A if degree==1 else A2,degree*logA)
        joint=bd.add(history['k'],bd.product(A,history['m']))
        targets[current.controls.ROWS[1]]=bd.quotient(joint,bd.scaled(A2,mu),2*logA+self.phase.outer.logmu)
        targets={key:targets[key] for key in current.controls.ROWS}
        maximum=LogUpper(c,c.mpf(max(ep(q.ZZZ.log)[1] for q in targets.values() if q.ZZZ.log is not None)))
        return dict(source_identity=self.identity,exact_Z_cell=ends,actual_17_chart_C3_transports=windows,
            actual_normalized_N_scaled_histories_C3=record(history),actual_terminal_amplitude_C3=[v.record() for v in amplitude],
            actual_positive_terminal_amplitude_log_lower=logA,actual_positive_mu_log=self.phase.outer.logmu,
            actual_N_scaled_five_target_C3_ranges=record(targets),actual_joint_angular_numerator_C3_range=record(joint),
            actual_five_target_ZZZ_maximum=maximum.record(),actual_five_target_C3_function_handles=target.encoded(self.target.target_functions()),
            original_independent_P0_in_source_and_not_added_twice=True,uniform_for_all_integer_N_ge_N0=True,
            epsilon_domain='[0,1/N0]',range_caps_not_point_function_values=True,actual_C3_repaired_limit_controls_installed=False)


def source_bindings():
    return dict(primitive_bounds=current.ast_binding(primitive_bounds),density_bounds=current.ast_binding(density_bounds),
        source_projection=current.ast_binding(target.CurrentC3SourceTargets.source_packet),
        original_log_bound_arithmetic=current.ast_binding(LogUpper),
        original_own_rate_weights=current.ast_binding(lower.CurrentC2TargetRanges.weights),
        third_sigma_cap=current.ast_binding(sigma_third_cap),
        original_signed_C3_receipt=target.RECEIPT,original_signed_C3_receipt_sha256=sha(target.RECEIPT))


def run(target_field=None,owner=None,field=None):
    began=time.monotonic()
    with mp.workdps(540):
        field=field if field is not None else CurrentC3TargetRanges(target_field=target_field,owner=owner,require_checked=False)
        field.hashes[Path(__file__).name]=sha(Path(__file__).name)
        transports=[field.transport(ends) for ends in outer.CELLS]
        report=dict(candidate_actual_C3_target_ranges_constructed=True,source_family=field.identity,
            actual_current_regular_majorant_theorem=field.proof,actual_four_Z_C3_target_transports=transports,
            actual_C3_function_source_bindings=source_bindings(),actual_C3_repaired_limit_controls_installed=False,
            **dict.fromkeys(GATES+OPEN,False),input_hashes=field.hashes,execution_seconds=time.monotonic()-began)
        (HERE/NAME).write_bytes(gzip.compress((json.dumps(phase.encoded(record(report)),separators=(',',':'))+'\n').encode(),mtime=0))
    return field


if __name__=='__main__':run()
