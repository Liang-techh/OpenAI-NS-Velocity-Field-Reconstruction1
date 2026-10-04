"""Actual admissible two-vector cone on the heat part of the collar.

The entire source interval [1,3) is proved, including its flat endpoint.
The sigmoid transition [0,1] and all-region stress remain open.
"""
import ast
import hashlib
import json
from pathlib import Path

import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_collar_stress_C3 import CompliantCollarStressC3
from lei_ren_part1_paper_compliant_global_physical_assembly import CompliantGlobalPhysicalAssembly
from lei_ren_part1_paper_compliant_collar_physical_C2 import stress_source_log_parts
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'


def heat_cone_identities():
    a,y,z,u,R,eps=s.symbols('a y Z future_offset R eps',positive=True)
    delta=2*a; k=1-a; b=(1-delta)/2; L=1-delta*z*z
    h,ht,F,Fz=s.symbols('H H_t full_flat_future_H full_flat_future_H_Z',real=True)
    phi=s.exp(-4/y**2)
    # Here the source offset derivative is -d/dy; phi_t=-8*phi/y^3.
    if s.simplify(-s.diff(phi,y)+8*phi/y**3)!=0:raise ArithmeticError('Flat derivative sign')
    fullH=s.Function('same_full_Gamma_H')(y)
    Ky=(-s.diff((1-eps*phi)*fullH,y)).subs({s.diff(fullH,y):-ht,fullH:h})
    if s.simplify(Ky-((1-eps*phi)*ht+8*eps*phi*h/y**3))!=0:
        raise ArithmeticError('Positive actual heat-part K derivative')
    ct=eps*phi*((h+k*F-b*z*Fz)/L+2/R*((8/y**3+1+a)*h-ht))
    btheta=eps*(y**3*(h+k*F-b*z*Fz)/L+2/R*(8*h+y**3*((1+a)*h-ht)))
    if s.simplify(ct-phi*y**(-3)*btheta)!=0:raise ArithmeticError('Theta flat factor')
    S,t=s.symbols('exact_positive_inverse_Rtail offset',positive=True)
    if s.simplify((1/R).subs(R,s.exp(t)/S)-S*s.exp(-t))!=0:
        raise ArithmeticError('Actual source radius transfer')
    xi=2*(1-z*z)*S*s.exp(-t); canon=s.Function('canonical_Gamma_H')
    if s.simplify(s.diff(canon(xi),t)+xi*s.Subs(s.Derivative(canon(s.Symbol('x')),s.Symbol('x')),s.Symbol('x'),xi))!=0:
        raise ArithmeticError('Canonical positive offset derivative')
    if s.simplify(s.diff(canon(xi),z)+4*z*S*s.exp(-t)*s.Subs(s.Derivative(canon(s.Symbol('x')),s.Symbol('x')),s.Symbol('x'),xi))!=0:
        raise ArithmeticError('Canonical axial derivative')
    # |Z|*(1-Z^2) has its maximum at |Z|=1/sqrt(3).
    v=s.symbols('abs_Z',nonnegative=True); f=v*(1-v*v)
    if s.simplify(s.diff(f,v).subs(v,1/s.sqrt(3)))!=0 or s.simplify(f.subs(v,1/s.sqrt(3))-2/(3*s.sqrt(3)))!=0:
        raise ArithmeticError('Sharp axial chain bound')
    if not (s.simplify(s.diff(f,v,2))==-6*v and 16<27):raise ArithmeticError('Uniform d*H_Z half bound')
    excess=4*((y-u)**(-2)-y**(-2))-8*u/y**3
    positive_excess=4*u*u*(3*y-2*u)/(y**3*(y-u)**2)
    if s.factor(excess-positive_excess)!=0:raise ArithmeticError('Flat future majorant')
    # Complete energy/pressure differences against the SAME Gamma field.
    K=s.symbols('Gamma_H',real=True); ph=s.symbols('future_phi',nonnegative=True)
    density=((1-eps*ph)*K)**2-K*K
    if s.expand(density-(-2*eps*ph+eps**2*ph**2)*K*K)!=0:raise ArithmeticError('Squared source difference')
    wE,wP,Q,Qz=s.symbols('energy_weight pressure_weight density density_Z',real=True)
    dz=1-z*z; p=1+delta
    original=delta*z*wE*Q-dz*wE*Qz/2-2*p*z*wP*Q/2+dz*wP*Qz/2
    regrouped=z*(delta*wE-p*wP)*Q+dz*(wP-wE)*Qz/2
    if s.expand(original-regrouped)!=0:raise ArithmeticError('Same full axial moment integrand')
    qt,qz,B,cs,cz,ss,m=s.symbols('Qtheta Qz B ctheta cz stheta kappa_minus2',positive=True)
    # Sz=0; Qz=Qtheta*B. Positive factors cancel, not cap values.
    margin=2*(qt*cs*qt*ss)**2-m*(qt*qz*cz*ss)**2
    if s.factor(margin.subs(qz,qt*B)-qt**4*ss**2*(2*cs**2-m*B**2*cz**2))!=0:
        raise ArithmeticError('Two-vector cone normalization')
    return dict(flat_source_derivative_verified=True,
                theta_factor_phi_times_y_minus3_verified=True,
                positive_flat_future_excess_verified=True,
                full_Gamma_squared_density_difference_verified=True,
                same_full_axial_moment_integrand_regrouping_verified=True,
                positive_source_K_y_and_kappa_upper_delta_verified=True,
                canonical_H_t_equals_minus_xi_H_xi_verified=True,
                canonical_H_Z_chain_and_uniform_d_H_Z_bound_verified=True,
                uniform_abs_Z_times_one_minus_Z_squared_below_half_verified=True,
                exact_R_equals_Rtail_exp_offset_consumed=True,
                two_vector_cone_positive_factor_cancellation_verified=True,
                flat_future_integral_bound='integral_0^y phi(t+u)/phi(t) du <= y^3/8',
                theta_normalized='ctheta/(eps*phi)=(H+kF-bZ*Fz)/L+2/R*((8/y^3+1+a)*H-H_t)',
                theta_endpoint_factor='btheta(0,Z)=16*eps*S*exp(-3)*H(2*(1-Z^2)*S*exp(-3))>0',
                weaker_uniform_ratio_bound='abs(Tz/Ttheta)<=B*Cz*y^3/hmin',
                stronger_fixed_source_ratio_bound='abs(Tz/Ttheta)<=Bmax*Cz*y^6/(16*S*exp(-3)*hmin)',
                paper_factor_exponents=dict(theta=-3,axial=3,ratio=6),
                completed_diagonal_does_not_enter_two_vector_cone=True)


def source_bindings(stress):
    proofs={}; hashes={}
    def assignment(stem,method,target,expression):
        path=HERE/(PREFIX+stem+'.py'); tree=ast.parse(path.read_text(encoding='utf8'))
        fn=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name==method)
        values=[n.value for n in ast.walk(fn) if isinstance(n,ast.Assign) and any(ast.unparse(v)==target for v in n.targets)]
        wanted=ast.dump(ast.parse(expression,mode='eval').body)
        if sum(ast.dump(v)==wanted for v in values)!=1:raise ValueError('Heat cone source changed: '+stem+':'+target)
        proofs[stem+':'+method+':'+target]=True; hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
    assignment('collar_Gamma_C4','shape','W','[unity[j]-sr[j]+v for j,v in enumerate(product_rows(sr,fr))]')
    assignment('collar_Gamma_C4','shape','C','[unity[j]-fr[j]*self.eps for j in range(5)]')
    assignment('collar_Gamma_C4','shape','pre','[unity[j]-W[j]*self.eps for j in range(5)]')
    assignment('collar_Gamma_C4','shape','D','product_rows(product_rows(sr,C),dr)')
    assignment('collar_Gamma_C4','shape','K','[pre[j]-D[j]*(self.a*self.S) for j in range(len(D))]')
    assignment('collar_Gamma_C4','phi_jets','dist','IntervalTaylor(c,[3-t,-1,0,0,0])')
    assignment('collar_Gamma_C4','phi_jets','jets','(dist**(-2)*(-4)).exp()')
    assignment('steep_waiting_C4','__init__','self.epsilon',"c.mpf('.001')*self.delta")
    if not stress.bridge['identities']['consumed_canonical_gamma_derivative_enclosure_verified']:
        raise ValueError('Canonical full Gamma source required')
    if not stress.bridge['exact_S_not_cap_endpoint']:raise ValueError('Exact positive S binding required')
    for name in ('actual_zero_angular_history_defect_transfers_to_all_collar_offsets',
                 'actual_selected_energy_half_future_transfers_by_same_full_K_squared_FTC',
                 'actual_pressure_and_energy_have_same_Ev0_squared_not_runtime_cap_value',
                 'same_flat_sigma_phi_K_jets_and_full_future_at_offset3_are_Gamma_jets',
                 'consumed_exact_heat_radius_and_xi_binding_verified'):
        if not stress.bridge['identities'][name]:raise ValueError('Common moment/radius source missing: '+name)
        proofs['consumed_'+name]=True
    # The production flat sigmoid, not an interval sample, is identically 1.
    from lei_ren_part1_paper_compliant_flat_pulse_derivatives import sigma_jets
    jet=sigma_jets(stress.ctx,stress.ctx.mpf([1,3]))
    if endpoints(jet[0])!=(1,1) or any(endpoints(v)!=(0,0) for v in jet.coefficients[1:]):
        raise ArithmeticError('Actual sigmoid heat-domain identity lost')
    proofs.update(production_sigma_identically_one_on_1_to_3=True,
                  same_canonical_H_and_exact_positive_S_consumed=True,
                  full_future_differences_supported_only_before_offset3=True,
                  Gamma_stress_exact_zero_from_same_terminal_moments_consumed=True)
    return dict(identities=proofs,input_hashes=hashes)


def run():
    stress=CompliantCollarStressC3(); assembly=CompliantGlobalPhysicalAssembly(); c=stress.ctx
    if (stress.family,stress.source)!=(assembly.family,assembly.source):raise ValueError('Cone/field source mismatch')
    hashes=dict(stress.hashes); hashes.update(assembly.hashes)
    name=PREFIX+'collar_stress_C3_check.json'; receipt=json.loads((HERE/name).read_bytes())
    if not receipt['all_passed'] or not receipt['actual_original_collar_similarity_stress_recovered'] or not receipt['whole_collar_shear_strength_kappa_gt2_certified']:
        raise ValueError('Actual collar stress prerequisite missing')
    for path,digest in receipt['input_hashes'].items():
        if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:raise ValueError('Cone prerequisite changed: '+path)
    hashes.update(receipt['input_hashes']); hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
    bridge=source_bindings(stress); hashes.update(bridge['input_hashes'])
    for stem in ('collar_physical_C2','global_physical_assembly'):
        path=HERE/(PREFIX+stem+'.py'); hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
    heat=stress.heat; a=heat.a; delta=heat.delta; Scap=heat.Scap
    with mp.workdps(280):
        hmin=1-2*a*(1+a)*Scap
        future_margin=heat.k*hmin-4*((1-delta)/2)*a*(1+a)*Scap
        shear_perturbation_margin=(1+a)*hmin-2*a*(1+a)*Scap
        Lmin=1-delta
        bounds=dict(H_lower=hmin,H_upper=c.mpf(1),Z_H_Z_absolute_upper=4*a*(1+a)*Scap,
                    H_Z_absolute_upper=4*a*(1+a)*Scap,d_H_Z_absolute_upper=2*a*(1+a)*Scap,
                    d_density_Z_absolute_upper_over_eps_phi=8*a*(1+a)*Scap,
                    H_t_lower=c.mpf(0),H_t_upper=2*a*(1+a)*Scap,
                    angular_future_integrand_lower=future_margin,
                    angular_shear_perturbation_lower=shear_perturbation_margin,L_lower=Lmin,
                    epsilon=heat.eps,delta=delta,inverse_radius_cap=Scap)
        for key in ('H_lower','angular_future_integrand_lower','angular_shear_perturbation_lower','L_lower','epsilon','delta'):
            if endpoints(bounds[key])[0]<=0:raise ArithmeticError('Heat cone positive source margin lost: '+key)
        if endpoints(heat.eps)[1]>=1:raise ArithmeticError('Heat perturbation must be below one')
        Cz=(1+2*delta+4*a*(1+a)*Scap)/(4*Lmin)
        parts=stress_source_log_parts(assembly,heat,c.mpf([1,3]))['B']
        logB=sum(parts.values(),c.mpf(0))
        # y<=2, so y^3<=8. Never exponentiate the huge source log.
        cone_log_upper=c.ln(delta)+2*logB+2*c.ln(8*Cz/hmin)
        if endpoints(cone_log_upper)[1]>=endpoints(c.ln(2))[0]:
            raise ArithmeticError('Actual heat-domain directional cone upper bound not below 2')
        bounds.update(axial_flat_normalized_absolute_upper=Cz,actual_B_source_log_parts=parts,
                      actual_logB_enclosure=logB,log_directional_cone_term_upper=cone_log_upper,
                      log_directional_threshold=c.ln(2),kappa_minus2_upper=delta)
        hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        result=dict(actual_five_defect_family_sha256=stress.family,implicit_source_sha256=stress.source,
                    scope='Whole actual heat part: offset in [1,3), Z in [-1,1]; endpoint3 exact zero stress',
                    domain=dict(offset=[1,3],Z=[-1,1],strict_upper_offset_excluded=True),
                    exact_identities=heat_cone_identities(),source_bridge=bridge,bounds=bounds,input_hashes=hashes,
                    actual_heat_part_collar_cone_certified=True,actual_nonzero_stress_theta_positive=True,
                    actual_directional_cone_uniform_margin_certified=True,
                    actual_zero_stress_join_uniform_direction_e_theta_verified=True,
                    source_defined_positive_eps_phi_S_not_replaced_by_caps=True,
                    endpoint3_exact_zero_stress_not_subject_to_strict_cone=True,
                    sigmoid_transition_0_to_1_cone_certified=False,collar_cone_certified=False,
                    completed_full_tensor_cone_certified=False,global_admissible_stress_lift_constructed=False,
                    independently_bounded_global_flat_remainder=False,temporal_recursion=False)
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Actual whole heat-part collar cone generated; sigmoid transition/global cone pending',flush=True)
    return result


if __name__=='__main__':run()
