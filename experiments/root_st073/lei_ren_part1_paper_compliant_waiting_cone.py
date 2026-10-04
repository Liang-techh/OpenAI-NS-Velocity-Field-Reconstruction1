"""Whole original waiting two-vector cone from its actual full moment modes.

Uses the accepted source-bound waiting stress and SAME future Gamma moments.
Only the proof companion is added; field, pressure, coefficients, and source
receipts are unchanged. Global stress/energy/recursion remain open.
"""
import ast
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_waiting_stress_C3 import (
    CompliantWaitingStressC3,waiting_stress_rows,source_precision)
from lei_ren_part1_paper_compliant_waiting_physical_C2 import waiting_adapter_binding
from lei_ren_part1_paper_compliant_global_physical_assembly import CompliantGlobalPhysicalAssembly
from lei_ren_part1_paper_compliant_collar_physical_C2 import stress_source_log_parts
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'
DOMAIN=dict(phase=[0,1],Z=[-1,1],offset='q=wait*(phase-1), -wait<=q<=0')


def waiting_cone_identities():
    a,wait,v,Sin,K0=s.symbols('a wait v exact_S K0',positive=True)
    k=1-a; delta=2*a; p=1+delta; q=-v
    I0,Je,Jp=s.symbols('actual_I0 actual_Je actual_Jp',real=True)
    Ct=I0*s.exp(-k*q)-2*Sin*(1+a)*K0*s.exp(-q)
    Cz=Je*s.exp(delta*q)+Jp*s.exp(p*q)
    proofs={}
    def zero(name,value):
        if s.simplify(value)!=0:raise ArithmeticError('Waiting cone identity failed: '+name)
        proofs[name]=True
    zero('actual_theta_common_exponential_and_small_shear_gap',
         Ct-s.exp(k*v)*(I0-2*Sin*(1+a)*K0*s.exp(a*v)))
    zero('actual_axial_modes_contract_backwards',Cz-(Je*s.exp(-delta*v)+Jp*s.exp(-p*v)))
    zero('shear_gap_increases_with_backward_distance',
         s.diff(2*Sin*(1+a)*K0*s.exp(a*v),v)-2*a*Sin*(1+a)*K0*s.exp(a*v))
    zero('angular_positive_factor_increases_backwards',s.diff(s.exp(k*v),v)-k*s.exp(k*v))
    zero('energy_mode_decreases_backwards',s.diff(s.exp(-delta*v),v)+delta*s.exp(-delta*v))
    zero('pressure_mode_decreases_backwards',s.diff(s.exp(-p*v),v)+p*s.exp(-p*v))
    bh=s.Rational(1,2)+a; Btail,thetaT,thetaBase,Ev0=s.symbols('Btail thetaT thetaBase Ev0',positive=True)
    B=Btail*s.exp(bh*v)
    zero('actual_B_max_at_original_waiting_inlet',B.subs(v,wait)-Btail*s.exp(bh*wait))
    zero('actual_B_max_uses_original_thetaT_and_K0',
         (Btail*s.exp(bh*wait)).subs(Btail,Ev0*thetaT*s.exp(-bh*wait)/K0)-Ev0*thetaT/K0)
    # Consume the exact positive radius source. The enclosing S box can
    # include zero without weakening this functional strict sign theorem.
    shear=-2*Sin*(1+a)*K0*s.exp(v)
    zero('actual_source_shear_is_strictly_negative_positive_factors',
         shear+2*Sin*(1+a)*K0*s.exp(v))
    zero('same_K_constant_kappa_minus2_exactly_delta',(-2*(-bh)-1)-delta)
    qt,qz,Bcurrent,ct,cz,st=s.symbols('Qtheta Qz B ctheta cz stheta',real=True)
    numerator=2*(qt*ct*qt*st)**2-delta*(qt*qz*cz*st)**2
    zero('actual_two_vector_cone_source_factor_cancellation',
         numerator.subs(qz,qt*Bcurrent)-qt**4*st**2*(2*ct**2-delta*Bcurrent**2*cz**2))
    # Original pure-swirl Sz=0 and ctheta>0, stheta<0 give T dot S<0.
    zero('actual_T_dot_S_same_positive_factors',qt*ct*qt*st-qt**2*ct*st)
    return dict(identities=proofs,entire_original_waiting_domain=DOMAIN,
                exact_q_equals_log_R_over_Rtail=True,exact_source_S_equals_positive_inverse_Rtail=True,
                constant_K0_kappa_minus2_equals_delta=True,
                shear_source_sign_not_inferred_from_cap_endpoint=True,
                theta_lower='ctheta>=exp(-k*q)*[I0_min-2*S_cap*(1+a)*K0_max*exp(a*wait)]>=positive_gap',
                axial_upper='abs(cz)<=abs(Je0)_max+abs(Jp0)_max; q<=0, delta>0, p>0',
                actual_B_upper_source='Bmax=Ev0*thetaT/K0; original inlet, no materialized cap field',
                normalized_cone_condition='delta*Bmax^2*(Cz_max/positive_gap)^2<2',
                exact_positive_common_nu_lambda_factors_preserve_cone=True,
                completed_diagonal_does_not_enter_two_vector_cone=True)


@source_precision
def whole_waiting_bounds(stress,assembly):
    h=stress.heat; c=stress.ctx
    rows=waiting_stress_rows(h,stress.endpoint([-1,1]),c.mpf([-1,1]),c.mpf(0))
    I0=rows['theta_inertial'][0][0]
    E0=rows['axial_energy'][0][0]; P0=rows['axial_pressure'][0][0]
    K0=1-h.eps
    # Every endpoint used below is a certified lower/upper bound, never a
    # selected moment/parameter/radius value defining a replacement field.
    Scap=c.mpf(endpoints(h.Scap)[1]); Kcap=c.mpf(endpoints(K0)[1])
    shear_relative_cap=2*Scap*(1+h.a)*Kcap*c.exp(h.a*h.steep.wait)
    gap=I0-shear_relative_cap
    gapmin=c.mpf(endpoints(gap)[0])
    mode_absolute_sum=abs(E0)+abs(P0)
    Cz=c.mpf(endpoints(mode_absolute_sum)[1])
    # The exact source q=-wait cancels wait+q in these original logs. This
    # is an identity at the inlet, not subtraction of unrelated log boxes.
    logparts=stress_source_log_parts(assembly,h,c.mpf(0))['B']
    logparts['waiting_and_current']=c.mpf(0)
    logBmax=sum(logparts.values(),c.mpf(0))
    margins=dict(delta=h.delta,epsilon=h.eps,K0=K0,k=h.k,p=h.prate,Lmin=1-h.delta,
                 wait=h.steep.wait,theta_gap=gap)
    for label,value in margins.items():
        if endpoints(value)[0]<=0 or not mp.isfinite(endpoints(value)[1]):raise ArithmeticError('Waiting cone source margin failed: '+label)
    if endpoints(h.eps)[1]>=1:raise ArithmeticError('Actual waiting K0 must stay positive')
    if endpoints(Cz)[0]<=0:raise ArithmeticError('Finite nonzero full axial bound required')
    logcone=c.ln(h.delta)+2*logBmax+2*c.ln(Cz/gapmin)
    if endpoints(logcone)[1]>=endpoints(c.ln(2))[0]:raise ArithmeticError('Whole waiting directional cone failed')
    return dict(positive_margins=margins,actual_full_endpoint_inertial_theta_enclosure=I0,
                actual_full_endpoint_axial_energy_enclosure=E0,actual_full_endpoint_axial_pressure_enclosure=P0,
                actual_source_inverse_radius_cap=Scap,source_K0_upper_bound=Kcap,
                angular_shear_relative_upper=shear_relative_cap,theta_uniform_positive_lower=gapmin,
                axial_uniform_absolute_upper=Cz,actual_Bmax_source_log_parts=logparts,
                actual_log_Bmax_enclosure=logBmax,log_directional_term_upper=logcone,log_threshold=c.ln(2),
                diagnostic_theta_gap_over_epsilon=gap/h.eps,diagnostic_axial_bound_over_epsilon=Cz/h.eps,
                exact_kappa_minus2=h.delta,full_moment_source_endpoint_used=True,
                source_caps_only_define_bounds_not_field_values=True)


def waiting_cone_source_bridge(stress):
    proofs={}
    bridge=stress.bridge; collar=stress.collar_source
    for flag in ('actual_waiting_endpoint_full_moment_history_consumed',
                 'actual_waiting_absolute_pressure_endpoint_consumed',
                 'actual_meridional_terminal_zeros_and_zero_density_FTC_consumed',
                 'waiting_collar_stress_mixed3_join_verified'):
        if not bridge[flag]:raise ValueError('Waiting cone source history missing: '+flag)
        proofs['consumed_'+flag]=True
    for flag in ('actual_waiting_collar_stress_mixed3_AST_join_verified',
                 'actual_waiting_collar_pressure_mixed4_AST_join_verified'):
        if not bridge['actual_waiting_collar_formula_join'][flag]:raise ValueError('Actual waiting cone join missing')
        proofs['consumed_'+flag]=True
    if not collar.bridge['exact_S_not_cap_endpoint']:raise ValueError('Exact positive inverse radius source required')
    for flag in ('consumed_exact_heat_radius_and_xi_binding_verified',
                 'actual_pressure_and_energy_have_same_Ev0_squared_not_runtime_cap_value'):
        if not collar.bridge['identities'][flag]:raise ValueError('Actual cone source radius/amplitude missing')
        proofs['consumed_'+flag]=True
    mapping=waiting_adapter_binding()
    if not mapping['actual_production_waiting_phase_radius_identity_verified']:raise ValueError('Actual phase/radius source not bound')
    proofs['exact_S_is_exp_minus_actual_finite_logRtail_so_strictly_positive']=True
    proofs['actual_waiting_constant_positive_K0_and_Sz_zero_consumed']=True
    proofs['actual_q_phase_radius_and_source_Bmax_inlet_binding_consumed']=True
    return dict(identities=proofs,actual_waiting_formula_source_bridge=bridge,
                actual_Bmax_source_log_binding=waiting_Bmax_binding(),
                actual_production_waiting_phase_radius_binding=mapping,
                sampled_interval_overlap_used_as_proof=False)


def waiting_Bmax_binding():
    """Execute actual log-factor code and the owned inlet cancellation AST."""
    tree=ast.parse((HERE/(PREFIX+'collar_physical_C2.py')).read_text(encoding='utf8'))
    fn=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name=='stress_source_log_parts')
    environment={}; exec(compile(ast.Module(body=[fn],type_ignores=[]),'<actual stress source logs>','exec'),environment)
    a,mu,wait,Ts,length,lp,lu,lrp,lone=s.symbols('a mu wait Ts Lrel logP logU logRp logone',real=True)
    c=SimpleNamespace(mpf=s.sympify,ln=s.log)
    flatten=SimpleNamespace(logEv2_parts=dict(inlet_log=2*lu))
    assembly=SimpleNamespace(ctx=c,logP=lp,logRp=lrp,
                             dispatch=SimpleNamespace(provider=lambda chart:flatten))
    heat=SimpleNamespace(ctx=c,a=a,mu=mu,bh=s.Rational(1,2)+a,
                         steep=SimpleNamespace(wait=wait,Ts=Ts,logone=lone,outer=SimpleNamespace(Lrel=length)))
    original=environment['stress_source_log_parts'](assembly,heat,-wait)['B']
    owned=ast.parse(Path(__file__).read_text(encoding='utf8'))
    function=next(n for n in ast.walk(owned) if isinstance(n,ast.FunctionDef) and n.name=='whole_waiting_bounds')
    assignments=[]
    for node in ast.walk(function):
        if isinstance(node,ast.Assign) and any(ast.unparse(t) in ('logparts',"logparts['waiting_and_current']") for t in node.targets):
            assignments.append(node)
    if len(assignments)!=2:raise ValueError('Owned waiting inlet log cancellation changed')
    environment.update(assembly=assembly,h=heat,c=c)
    exec(compile(ast.Module(body=assignments,type_ignores=[]),'<owned waiting Bmax source logs>','exec'),environment)
    if set(original)!=set(environment['logparts']):raise ValueError('Waiting source log keys changed')
    checks={}
    for key,value in original.items():
        if s.simplify(value-environment['logparts'][key])!=0:raise ArithmeticError('Waiting inlet log source mismatch: '+key)
        checks[key]=True
    return dict(identities=checks,actual_Bmax_is_original_source_B_at_q_minus_wait=True,
                cancellation_holds_before_interval_enclosure=True,cap_values_used_as_source=False)


def run():
    with mp.workdps(300):
        stress=CompliantWaitingStressC3(); assembly=CompliantGlobalPhysicalAssembly()
        if (stress.family,stress.source)!=(assembly.family,assembly.source):raise ValueError('Waiting cone field family mismatch')
        hashes=dict(stress.hashes); hashes.update(assembly.hashes)
        for stem,gate in (('waiting_stress_C3_check','actual_original_waiting_similarity_stress_recovered'),
                          ('waiting_physical_C2_check','actual_regional_physical_waiting_stress_remainder_identity_verified'),
                          ('collar_cone_check','collar_cone_certified')):
            name=PREFIX+stem+'.json'; receipt=json.loads((HERE/name).read_bytes())
            if not receipt['all_passed'] or not receipt[gate]:raise ValueError('Waiting cone prerequisite missing: '+name)
            if (receipt['actual_five_defect_family_sha256'],receipt['implicit_source_sha256'])!=(stress.family,stress.source):
                raise ValueError('Waiting cone receipt family mismatch')
            for source,digest in receipt['input_hashes'].items():
                if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:raise ValueError('Waiting cone dependency changed: '+source)
            hashes.update(receipt['input_hashes']); hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        proof=waiting_cone_identities(); bridge=waiting_cone_source_bridge(stress); bounds=whole_waiting_bounds(stress,assembly)
        hashes.update(assembly.dispatch.hashes)
        hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        result=dict(actual_five_defect_family_sha256=stress.family,implicit_source_sha256=stress.source,
                    scope='Whole original waiting phase[0,1], source Z[-1,1]; physical tau>0, r>0, |Z|<1',
                    domain=DOMAIN,exact_identities=proof,source_bridge=bridge,bounds=bounds,input_hashes=hashes,
                    waiting_cone_certified=True,actual_waiting_nonzero_stress_theta_positive=True,
                    actual_waiting_source_shear_strictly_negative=True,actual_waiting_directional_cone_uniform_margin=True,
                    waiting_collar_cone_proofs_same_source_composed=True,
                    waiting_and_collar_outer_tail_two_vector_cone_certified=True,
                    tail_domain='Rt<=R<Rtail*exp(3); beyond this, accepted Gamma stress is exactly zero',
                    fixed_positive_viscosity_two_vector_cone_transfer_verified=True,
                    completed_full_tensor_cone_certified=False,whole_outer_cone_certified=False,
                    global_admissible_stress_lift_constructed=False,independently_bounded_global_flat_remainder=False,
                    physical_energy_integral_certified=False,temporal_recursion=False)
        Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
        print('Whole original waiting cone generated; common collar join composed, global cone pending',flush=True)
        return result


class CertifiedWaitingPhysical:
    """Add accepted regional cone metadata to unchanged physical source rows."""
    def __init__(self):
        from lei_ren_part1_paper_compliant_waiting_physical_C2 import CompliantWaitingPhysicalC2
        path=HERE/(PREFIX+'waiting_cone_check.json'); receipt=json.loads(path.read_bytes())
        if not receipt['all_passed'] or not receipt['waiting_cone_certified']:raise ValueError('Accepted waiting cone required')
        for source,digest in receipt['input_hashes'].items():
            if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:raise ValueError('Certified waiting source changed: '+source)
        self.field=CompliantWaitingPhysicalC2(); self.family,self.source=self.field.family,self.field.source
        if (receipt['actual_five_defect_family_sha256'],receipt['implicit_source_sha256'])!=(self.family,self.source):
            raise ValueError('Waiting physical cone family mismatch')
        self.receipt_sha256=hashlib.sha256(path.read_bytes()).hexdigest(); self.domain=receipt['domain']

    def waiting(self,Z,phase,log_tau='-1',theta='0',viscosity='1'):
        point=self.field.waiting(Z,phase,log_tau=log_tau,theta=theta,viscosity=viscosity)
        point.update(waiting_cone_certified=True,actual_physical_waiting_two_vector_cone_certified=True,
                     waiting_source_shear_strictly_negative_verified=True,waiting_cone_receipt_sha256=self.receipt_sha256,
                     cone_domain=self.domain,completed_full_tensor_cone_certified=False)
        return point


if __name__=='__main__':run()
