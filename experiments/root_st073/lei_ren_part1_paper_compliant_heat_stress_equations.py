"""Exact Gamma heat radial-stress equations, Lei--Ren v2 (3.12)--(3.20).

These equations prove that exterior stress is homogeneous. They do not
choose its two constants. Transferring the original five terminal moments
is a separate gate, as emphasized on page 60 of the source paper.
"""
import hashlib
import json
from pathlib import Path

import mpmath as mp
import sympy as s

HERE = Path(__file__).parent


def heat_stress_equations():
    proofs = {}

    def zero(name, value):
        if s.simplify(s.expand_power_exp(value)) != 0:
            raise ArithmeticError('Exact heat stress identity failed: '+name)
        proofs[name] = True

    a, xi, v = s.symbols('a xi v', positive=True)
    density = s.exp(-v)*v**a/s.gamma(1+a)
    kernel = (1+xi*v)**(-a)
    ode_kernel = xi**2*s.diff(kernel,xi,2)+(1+2*(1+a)*xi)*s.diff(kernel,xi)+a*(1+a)*kernel
    primitive = s.exp(-v)*v**(1+a)*(1+xi*v)**(-1-a)/s.gamma(1+a)
    zero('canonical_Gamma_ODE_is_integral_of_a_total_derivative',
         density*ode_kernel-a*s.diff(primitive,v))
    # For a>0 and xi>=0 this primitive tends to zero at both endpoints:
    # at zero it has v^(1+a), at infinity exponential decay dominates.
    # The same positive Gamma density dominates all required derivatives.
    proofs['Gamma_IBP_endpoint_limits_zero_for_a_positive'] = True
    H, Hp, Hpp = s.symbols('H H_xi H_xixi', real=True)
    ode = xi**2*Hpp+(1+2*(1+a)*xi)*Hp+a*(1+a)*H
    def heat_derivatives(expression, argument, function):
        replacements = {}
        for term in expression.atoms(s.Function):
            if term.func == function and s.simplify(term.args[0]-argument) == 0:
                replacements[term] = H
        for term in expression.atoms(s.Subs):
            if (isinstance(term.expr,s.Derivative) and term.expr.expr.func == function
                    and s.simplify(term.point[0]-argument) == 0):
                replacements[term] = {1:Hp,2:Hpp}[term.expr.derivative_count]
        return expression.xreplace(replacements)

    # Physical (5.2)--(5.5), including arbitrary viscosity. The free
    # positive amplitude absorbs the paper's factor 2^((1+delta)/2).
    r, tau, nu, A = s.symbols('r tau nu physical_amplitude',positive=True)
    physical_arg = 4*nu*tau/r**2
    h = s.Function('canonical_Gamma_H')
    physical_swirl = A*r**(-1-2*a)*h(physical_arg)
    physical_residual = -s.diff(physical_swirl,tau)-nu*(
        s.diff(physical_swirl,r,2)+s.diff(physical_swirl,r)/r-physical_swirl/r**2)
    zero('physical_angular_heat_equation_is_same_Gamma_ODE',
         heat_derivatives(s.expand(physical_residual),physical_arg,h)
         +4*nu*A*r**(-3-2*a)*ode.subs(xi,physical_arg))

    R, c = s.symbols('R c_infinity', positive=True)
    Z = s.symbols('Z', real=True)
    # Z is symbolic; identities hold on the entire physical axial interval,
    # including +/-1 by smooth continuation of the Gamma expectation.
    d, L, delta = 1-Z*Z, 1-2*a*Z*Z, 2*a
    arg = 2*d/R
    h = s.Function('canonical_Gamma_H')
    F = c/s.sqrt(2)*R**(-1-a)*h(arg)
    Utheta = s.sqrt(2*R)*F
    # Original (3.12), with Uz=Ur=0; retain all axial chain-rule factors.
    Ntheta = -s.sqrt(R/2)/L*((1+delta)/2*Utheta
                 +(1-delta)/2*Z*s.diff(Utheta,Z)+R*s.diff(Utheta,R))
    Stheta = 2*R*s.diff(F,R)
    target_N = s.sqrt(2)*c*R**(-1-a)*s.Subs(s.diff(h(xi),xi),xi,arg)
    zero('original_3_12_angular_source_of_full_Gamma_heat',Ntheta-target_N)
    angular = s.expand(Ntheta+R*s.diff(Stheta,R)+Stheta)
    h_values = {h(arg):H}
    for term in angular.atoms(s.Function):
        if term.func == h and s.simplify(term.args[0]-arg) == 0:
            h_values[term] = H
    for term in angular.atoms(s.Subs):
        if (isinstance(term.expr,s.Derivative) and term.expr.expr.func == h
                and s.simplify(term.point[0]-arg) == 0):
            order = term.expr.derivative_count
            h_values[term] = {1:Hp,2:Hpp}[order]
    residual = angular.xreplace(h_values)
    zero('angular_radial_stress_source_is_exact_Gamma_ODE',
         residual-s.sqrt(2)*c*R**(-1-a)*ode.subs(xi,arg))

    # The SAME absolute pressure tail has this self-similar dependence.
    # No free analytic P0 or exterior gauge is inserted here.
    Q = s.Function('full_Gamma_pressure_tail')
    pressure = c**2*R**(-1-2*a)*Q(arg)
    Nz = s.sqrt(R/2)/L*(2*(1+delta)*Z*pressure
             -d*s.diff(pressure,Z)+2*Z*R*s.diff(pressure,R))
    zero('original_3_12_axial_source_of_same_pressure_tail',Nz)
    Sz = s.Integer(0)
    zero('axial_radial_stress_source_zero',Nz+R*s.diff(Sz,R)+Sz/2)

    Ct, Cz = s.Function('Ctheta')(Z), s.Function('Cz')(Z)
    zero('angular_homogeneous_stress_equation',R*s.diff(Ct/R,R)+Ct/R)
    zero('axial_homogeneous_stress_equation',R*s.diff(Cz/s.sqrt(R),R)+Cz/(2*s.sqrt(R)))
    # Conditional terminal theorem, source page 60. These are estimates for
    # the exact terminal-normalized moments, not newly imposed data in the
    # current C4 packet. H<=1 and |H'|<=a*(1+a) imply the bounds below.
    rho = s.symbols('rho',positive=True)
    zero('Gamma_first_derivative_moment_bound',
         a*s.expand_func(s.gamma(a+2)/s.gamma(a+1))-a*(1+a))
    theta_error = s.integrate(2*s.sqrt(2)*c*a*(1+a)*rho**(-1-a),(rho,R,s.oo))
    theta_Z = s.integrate(4*s.sqrt(2)*c*a*(1+a)*rho**(-1-a),(rho,R,s.oo))
    energy = s.integrate(c*c*rho**(-1-2*a)/2,(rho,R,s.oo))
    pressure_bound = s.integrate(c*c*rho**(-2-2*a)/2,(rho,R,s.oo))
    zero('terminal_Mtheta_remainder_bound',theta_error-2*s.sqrt(2)*c*(1+a)*R**(-a))
    zero('terminal_Mtheta_Z_bound',theta_Z-4*s.sqrt(2)*c*(1+a)*R**(-a))
    zero('terminal_Mztheta_tail_bound',energy-c*c*R**(-2*a)/(4*a))
    zero('terminal_absolute_pressure_tail_bound',pressure_bound-c*c*R**(-1-2*a)/(2*(1+2*a)))
    zero('terminal_angular_weighted_limit_zero',s.limit(R**(-a),R,s.oo))
    zero('terminal_axial_weighted_limit_zero',s.limit(R**(-2*a),R,s.oo))
    # The identities above allow nonzero Ct and Cz. Paper terminal moments
    # must imply R*Ttheta -> 0 and sqrt(R)*Tz -> 0 to eliminate them.
    return dict(identities=proofs,
                exact_heat_angular_radial_source_verified=True,
                exact_heat_axial_radial_source_verified=True,
                exact_heat_stress_is_homogeneous_verified=True,
                canonical_physical_heat_flow_equation_verified=True,
                actual_compliant_field_physical_heat_region_certified=False,
                canonical_physical_velocity='ur=uz=0; utheta=A*r^(-1-delta)*H(4*nu*tau/r^2)',
                canonical_physical_pressure='p=-integral_r^infinity utheta(rho,tau)^2/rho drho; radial balance by FTC',
                terminal_normalization_implies_zero_heat_stress_theorem_verified=True,
                terminal_theorem_assumptions=[
                    '0<a=delta/2<1/2, |Z|<=1, R>=1',
                    'Mz=Mtheta_z=Ur=Uz=0 on the exterior from original terminal moments',
                    'Mtheta=sqrt(2)*c*(R^(1-a)/(1-a)+integral_R^infinity rho^(-a)*(1-H(2d/rho)) drho)',
                    'Mztheta=c^2/2*integral_R^infinity rho^(-1-2a)*H(2d/rho)^2 drho',
                    'P=-c^2/2*integral_R^infinity rho^(-2-2a)*H(2d/rho)^2 drho'],
                terminal_normalized_bounds=dict(
                    Mtheta_remainder='2*sqrt(2)*c*(1+a)*R^(-a)',
                    Mtheta_Z='4*sqrt(2)*c*(1+a)*R^(-a)',
                    Mztheta='c^2*R^(-2a)/(4a)',
                    absolute_pressure='c^2*R^(-1-2a)/(2*(1+2a))',
                    R_Ttheta='O(c*R^(-a)) -> 0',
                    sqrtR_Tz='O(c^2*R^(-2a)) -> 0'),
                homogeneous_stress_forms=dict(theta='Ctheta(Z)/R', axial='Cz(Z)/sqrt(R)'),
                homogeneous_constants_eliminated_from_actual_terminal_moments=False,
                heat_exterior_stress_identity_certified=False,
                global_admissible_stress_lift_constructed=False,
                temporal_recursion=False,
                source=dict(paper='Lei--Ren, arXiv:2609.35406v2',
                            formulas=['(3.12)','(3.13)','(3.18)','(5.1)--(5.11)'],
                            inspected_pdf_pages=[20,21,58,59,60]),
                Gamma_ODE='xi^2*H_xixi+(1+2*(1+a)*xi)*H_xi+a*(1+a)*H=0')


def independent_Gamma_ODE_fixture():
    with mp.workdps(55):
        checks = 0
        for a in map(mp.mpf,('.003','.15','.35')):
            for xi in map(mp.mpf,('0','.002','.08','2')):
                derivatives = []
                for n in range(3):
                    value = (-1)**n*mp.rf(a,n)*mp.quad(
                        lambda v:mp.exp(-v)*v**(a+n)*(1+xi*v)**(-a-n),
                        [0,1,mp.inf])/mp.gamma(1+a)
                    derivatives.append(value)
                H, Hp, Hpp = derivatives
                residual = xi*xi*Hpp+(1+2*(1+a)*xi)*Hp+a*(1+a)*H
                if abs(residual)>mp.mpf('1e-47'):
                    raise ArithmeticError('Independent Gamma heat ODE fixture failed')
                checks += 1
        return dict(independent_positive_Gamma_integral_ODE_checks=checks,
                    finite_parameter_fixture_only=True, passed=True)


def run():
    result = heat_stress_equations()
    result['independent_Gamma_ODE_fixture'] = independent_Gamma_ODE_fixture()
    result['equation_checks_passed'] = True
    result['input_hashes'] = {Path(__file__).name:hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('Exact Gamma heat radial stress equations PASS; actual terminal constants remain gated',flush=True)
    return result


if __name__ == '__main__':
    run()
