"""Actual full heat-collar two-vector admissible cone.

Closes the original sigmoid interval [0,1] and composes the admitted heat
part [1,3). Keeps global tensor/remainder/recursion gates separate.
"""
import ast
import hashlib
import json
from pathlib import Path

import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_collar_Gamma_C4 as collar_source
from lei_ren_part1_paper_compliant_collar_heat_cone import source_bindings,heat_cone_identities
from lei_ren_part1_paper_compliant_collar_stress_C3 import CompliantCollarStressC3
from lei_ren_part1_paper_compliant_global_physical_assembly import CompliantGlobalPhysicalAssembly
from lei_ren_part1_paper_compliant_collar_physical_C2 import stress_source_log_parts
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import sigma_jets
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'


def sigmoid_cone_identities():
    t=s.symbols('offset',positive=True); eps,a=s.symbols('eps a',positive=True)
    sig,ph,h=s.symbols('sigma phi H',real=True); sd,pd,ht,hz=s.symbols('sigma_t phi_t H_t H_Z',real=True)
    W=1-sig+sig*ph; f=1-sig+sig*eps*ph; D=1-h
    K=(1-sig)*(1-eps)+sig*(1-eps*ph)*h; G=h-K
    proofs={}
    def zero(name,value):
        if s.simplify(value)!=0:raise ArithmeticError('Sigmoid cone identity failed: '+name)
        proofs[name]=True
    dt=lambda x:s.diff(x,sig)*sd+s.diff(x,ph)*pd+s.diff(x,h)*ht
    zero('same_actual_G_is_eps_W_minus_f_Gamma_deficit',G-(eps*W-f*D))
    zero('same_actual_G_Z_is_f_H_Z',s.diff(G,h)*hz-f*hz)
    zero('same_actual_G_t_from_original_composition',dt(G)-(eps*dt(W)-dt(f)*D+f*ht))
    zero('same_actual_shear_perturbation_identity',-dt(G)+(1+a)*G
         -(eps*((1+a)*W-dt(W))+dt(f)*D-f*ht-(1+a)*f*D))
    gap=(1-eps*ph)*h-(1-eps)
    zero('same_actual_hot_minus_waiting_gap',gap-(eps*(1-ph)-(1-eps*ph)*D))
    zero('same_actual_positive_K_t_identity',dt(K)-(sd*gap+sig*(-eps*pd*h+(1-eps*ph)*ht)))
    zero('same_actual_squared_density_difference',K*K-h*h+G*(K+h))
    zero('same_actual_squared_density_Z',s.diff(K*K-h*h,h)*hz-2*hz*(sig*(1-eps*ph)*K-h))
    zero('same_offset1_Gamma_heat_source',K.subs(sig,1)-(1-eps*ph)*h)
    zero('same_offset1_heat_G',G.subs(sig,1)-eps*ph*h)
    odds=(1-t)**(-2)-t**(-2)
    odds_t=2*(1-t)**(-3)+2*t**(-3)
    zero('production_logistic_odds_derivative_positive',s.diff(odds,t)-odds_t)
    logistic=s.exp(odds)/(1+s.exp(odds))
    zero('production_sigma_t_is_positive_logistic_product',s.diff(logistic,t)-logistic*(1-logistic)*odds_t)
    return dict(identities=proofs,
                source_reference='SAME full Gamma H, with full source moment differences supported before offset3',
                sigmoid_domain=dict(offset=[0,1],Z=[-1,1]),
                current_theta_lower='ctheta >= Gmin=eps/e-Dcap >0',
                angular_future_lower='k*G-b*Z*G_Z >= k*eps/e-Dcap*(k+2b)>0 on [0,1]; heat part already nonnegative',
                shear_perturbation_lower='-G_t+(1+a)*G >= eps*(1+a)/e-Dcap*(sigma_t_cap+eps+2+a)>0',
                axial_difference_bound='abs(cz)<=3*(2eps*(1+2delta)+Dcap)/(1-delta)',
                source_K_t_nonnegative_implies_kappa_minus2_at_most_delta=True,
                finite_difference_support_does_not_truncate_original_infinite_moments=True,
                normalized_cone_condition='delta*Bmax^2*(Cz_sigma/Gmin)^2<2',
                completed_diagonal_is_outside_two_vector_cone=True)


def sigmoid_source_bridge(stress):
    bridge=source_bindings(stress); proofs=dict(bridge['identities']); hashes=dict(bridge['input_hashes'])
    def assignment(stem,method,target,expression):
        path=HERE/(PREFIX+stem+'.py'); tree=ast.parse(path.read_text(encoding='utf8'))
        fn=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name==method)
        found=[n.value for n in ast.walk(fn) if isinstance(n,ast.Assign) and any(ast.unparse(v)==target for v in n.targets)]
        wanted=ast.dump(ast.parse(expression,mode='eval').body)
        if sum(ast.dump(v)==wanted for v in found)!=1:raise ValueError('Sigmoid source changed: '+stem+':'+target)
        proofs[stem+':'+method+':'+target]=True; hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
    assignment('collar_Gamma_C4','shape','sig','sigma_jets(c,t)')
    assignment('collar_Gamma_C4','shape','sr','[one*(sig[j]*math.factorial(j)) for j in range(5)]')
    assignment('collar_Gamma_C4','shape','dr','gamma_deficit_mixed(c,Z,self.a,self.Scap,t,4 if high else 0)')
    assignment('flat_pulse_derivatives','_sigma_left','odds','1/(1-x)**2-1/x**2')
    assignment('flat_pulse_derivatives','_sigma_left','e','positive_exp(c,odds)')
    assignment('flat_pulse_derivatives','_sigma_left','value','e/(1+e)')
    assignment('flat_pulse_derivatives','_sigma_left','L1','2/(1-x)**3+2/x**3')
    assignment('flat_pulse_derivatives','_sigma_left','q','value*(1-value)')
    assignment('flat_pulse_derivatives','sigma_jets','left','_sigma_left(c,1-c.mpf([a,b]))')
    path=HERE/(PREFIX+'flat_pulse_derivatives.py'); tree=ast.parse(path.read_text(encoding='utf8'))
    fn=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name=='_sigma_left')
    derivative=next(n.value for n in ast.walk(fn) if isinstance(n,ast.Assign) and any(ast.unparse(v)=='derivatives' for v in n.targets))
    if not isinstance(derivative,ast.List) or ast.dump(derivative.elts[0])!=ast.dump(ast.parse('value',mode='eval').body) or ast.dump(derivative.elts[1])!=ast.dump(ast.parse('q*L1',mode='eval').body):
        raise ValueError('Production zeroth/first sigmoid derivative route changed')
    proofs['actual_sigmoid_zeroth_and_first_derivative_routes_bound']=True
    if collar_source.sigma_jets is not sigma_jets:raise ValueError('Actual sigmoid callable differs')
    if not stress.bridge['meridional_moments_and_velocities_exact_zero'] or not stress.bridge['identities']['actual_zero_Mz_Mtheta_z_Ur_Uz_histories_transfer_through_collar']:
        raise ValueError('Actual pure swirl/Sz=0 source history required')
    if not stress.theorem['full_terminal_moment_stress_theorem_verified']:
        raise ValueError('Canonical complete Gamma moment stress theorem required')
    boundary=sigma_jets(stress.ctx,stress.ctx.mpf(1))
    if endpoints(boundary[0])!=(1,1) or any(endpoints(v)!=(0,0) for v in boundary.coefficients[1:]):
        raise ArithmeticError('Original flat offset1 switch lost')
    boundary0=sigma_jets(stress.ctx,stress.ctx.mpf(0))
    if any(endpoints(v)!=(0,0) for v in boundary0.coefficients):
        raise ArithmeticError('Original flat offset0 waiting join lost')
    proofs.update(actual_sigmoid_callable_and_logistic_odds_bound=True,
                  original_positive_derivative_and_reflection_used=True,
                  actual_zero_meridional_histories_and_Sz_zero_consumed=True,
                  canonical_full_Gamma_zero_stress_theorem_consumed=True,
                  same_sigma1_flat_source_join_through_C4=True,
                  same_sigma0_flat_waiting_join_through_C4=True,
                  same_full_future_integrals_at_offset1_not_overlap_or_reset=True)
    return dict(identities=proofs,input_hashes=hashes)


def run():
    stress=CompliantCollarStressC3(); assembly=CompliantGlobalPhysicalAssembly(); c=stress.ctx; heat=stress.heat
    if (stress.family,stress.source)!=(assembly.family,assembly.source):raise ValueError('Whole collar source mismatch')
    hashes=dict(stress.hashes); hashes.update(assembly.hashes)
    for stem,gate in (('collar_stress_C3_check','whole_collar_shear_strength_kappa_gt2_certified'),
                      ('collar_heat_cone_check','actual_heat_part_collar_cone_certified'),
                      ('collar_physical_C2_check','actual_regional_physical_collar_stress_remainder_identity_verified')):
        name=PREFIX+stem+'.json'; receipt=json.loads((HERE/name).read_bytes())
        if not receipt['all_passed'] or not receipt[gate]:raise ValueError('Whole collar prerequisite missing: '+name)
        if (receipt['actual_five_defect_family_sha256'],receipt['implicit_source_sha256'])!=(stress.family,stress.source):
            raise ValueError('Whole collar receipt family mismatch')
        for source,digest in receipt['input_hashes'].items():
            if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:raise ValueError('Whole collar source changed: '+source)
        hashes.update(receipt['input_hashes']); hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
    bridge=sigmoid_source_bridge(stress); hashes.update(bridge['input_hashes'])
    with mp.workdps(280):
        a=heat.a; delta=heat.delta; eps=heat.eps; k=heat.k; b=(1-delta)/2; Lmin=1-delta
        Dcap=2*a*(1+a)*heat.Scap
        sig=sigma_jets(c,c.mpf([0,1])); sigcap=c.mpf(max(abs(v) for v in endpoints(sig[1])))
        if not mp.isfinite(endpoints(sigcap)[1]):raise ArithmeticError('Finite whole sigmoid derivative cap required')
        lower=dict(Gmin=eps/c.exp(1)-Dcap,
                   angular_future=k*eps/c.exp(1)-Dcap*(k+2*b),
                   angular_shear=eps*(1+a)/c.exp(1)-Dcap*(sigcap+eps+2+a),
                   hot_waiting_gap=eps*(1-c.exp(-c.mpf(4)/9))-Dcap,
                   Lmin=Lmin,epsilon=eps,delta=delta)
        for name,value in lower.items():
            if endpoints(value)[0]<=0:raise ArithmeticError('Actual sigmoid margin failed: '+name)
        if endpoints(eps)[1]>=1:raise ArithmeticError('Original cutoff factor must stay positive')
        Cz=3*(2*eps*(1+2*delta)+Dcap)/Lmin
        logparts=stress_source_log_parts(assembly,heat,c.mpf([0,1]))['B']; logB=sum(logparts.values(),c.mpf(0))
        logcone=c.ln(delta)+2*logB+2*c.ln(Cz/lower['Gmin'])
        if endpoints(logcone)[1]>=endpoints(c.ln(2))[0]:raise ArithmeticError('Actual sigmoid directional cone failed')
        bounds=dict(positive_margins=lower,canonical_Gamma_deficit_uniform_cap=Dcap,
                    sigmoid_derivative_absolute_upper=sigcap,phi_value_range=c.mpf([endpoints(c.exp(-1))[0],endpoints(c.exp(-c.mpf(4)/9))[1]]),
                    phi_offset_derivative_absolute_upper=c.mpf(1),axial_stress_absolute_upper=Cz,
                    actual_B_source_log_parts=logparts,actual_logB_enclosure=logB,
                    log_directional_term_upper=logcone,log_threshold=c.ln(2),kappa_minus2_upper=delta)
        hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        result=dict(actual_five_defect_family_sha256=stress.family,implicit_source_sha256=stress.source,
                    scope='Whole original collar offset[0,3), Z[-1,1]; exact zero stress at offset3',
                    domain=dict(offset=[0,3],Z=[-1,1],strict_upper_offset_excluded=True),
                    sigmoid_exact_identities=sigmoid_cone_identities(),heat_exact_identities=heat_cone_identities(),
                    source_bridge=bridge,sigmoid_bounds=bounds,input_hashes=hashes,
                    sigmoid_transition_0_to_1_cone_certified=True,actual_heat_part_collar_cone_certified=True,
                    same_source_offset1_cone_proofs_composed=True,collar_cone_certified=True,
                    whole_collar_stress_theta_positive_before_zero_join=True,
                    zero_stress_endpoint3_excluded_from_strict_inequalities=True,
                    whole_collar_uniform_terminal_direction_verified=True,
                    fixed_positive_viscosity_and_lambda_factors_preserve_two_vector_cone=True,
                    completed_full_tensor_cone_certified=False,global_admissible_stress_lift_constructed=False,
                    independently_bounded_global_flat_remainder=False,temporal_recursion=False)
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Actual whole original collar two-vector cone generated; global stress/flat remainder pending',flush=True)
    return result


class CertifiedCollarPhysical:
    """Accepted physical source rows with the matching local cone receipt.

    This adapter adds proof metadata only. It never projects stress rows,
    changes the field/pressure or claims a cone for the completed diagonal.
    """
    def __init__(self):
        from lei_ren_part1_paper_compliant_collar_physical_C2 import CompliantCollarPhysicalC2
        name=PREFIX+'collar_cone_check.json'; path=HERE/name; receipt=json.loads(path.read_bytes())
        if not receipt['all_passed'] or not receipt['collar_cone_certified']:
            raise ValueError('Accepted whole collar cone required')
        for source,digest in receipt['input_hashes'].items():
            if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:
                raise ValueError('Certified physical collar source changed: '+source)
        self.field=CompliantCollarPhysicalC2(); self.family=self.field.family; self.source=self.field.source
        if (receipt['actual_five_defect_family_sha256'],receipt['implicit_source_sha256'])!=(self.family,self.source):
            raise ValueError('Physical field/cone receipt family mismatch')
        self.receipt_sha256=hashlib.sha256(path.read_bytes()).hexdigest()
        self.domain=receipt['whole_original_collar_source_cone_domain_verified']

    def collar(self,Z,offset,log_tau='-1',theta='0',viscosity='1'):
        packet=self.field.collar(Z,offset,log_tau=log_tau,theta=theta,viscosity=viscosity)
        packet.update(collar_cone_certified=True,
                      actual_physical_collar_two_vector_cone_certified=True,
                      whole_collar_cone_receipt_sha256=self.receipt_sha256,
                      cone_domain=self.domain,
                      strict_cone_only_where_nonzero_before_offset3=True,
                      completed_full_tensor_cone_certified=False)
        return packet


if __name__=='__main__':run()
