"""Exact symbolic substitution from the core equations to the scaled map.

Checks all nonlinear terms, including the pressure derivative and swirl.
The algebra is parameter-independent and does not assert a matching result.
"""
import hashlib
import json
from pathlib import Path
import sympy as sp


def run():
    s,z,dt,L,d,eps,g,S=sp.symbols('s z delta L d epsilon g S')
    U,U_z,P,P_z,H0,W0=sp.symbols('U0 U0_z P0 P0_z H0 W0')
    phi,phi_s,phi_z=sp.symbols('Phi Phi_s Phi_z')
    psi,psi_s,psi_z,M,M_z,C,C_z=sp.symbols('Psi Psi_s Psi_z M M_z Pcal Pcal_z')
    W=W0-eps*((1-dt)*z*M+d*M_z)
    H=H0+eps*d*psi
    u=U+eps*psi
    # Original angular RHS divided by Lambda*L*F0.
    original_theta=eps*(W*(s*phi_s+phi)+H*phi_z
        +dt*phi/2-dt*z*phi*u)/L-g*H*phi/L
    beta=-(W0+dt/2-dt*z*U)/L
    chi=g*H0/L
    theta_terms=[-beta*phi,W0/L*s*phi_s,H0/L*phi_z,
        -eps*(1-dt)*z/L*M*phi,
        -eps*(1-dt)*z/L*M*s*phi_s,
        -eps*d/L*M_z*phi,-eps*d/L*M_z*s*phi_s,
        -eps*dt*z/L*psi*phi,eps*d/L*psi*phi_z,
        -g*d/L*psi*phi]
    theta_difference=sp.expand(original_theta+chi*phi-eps*sum(theta_terms))
    # p=P0+epsilon*Pcal; Pcal=F0²*integral_0^s Phi².
    original_z=(W*eps*s*psi_s+H*(U_z+eps*psi_z)
        +(1+dt)*(u-2*z*u*u)/2+d*(P_z+eps*C_z)
        -2*(1+dt)*z*(P+eps*C)-2*z*eps*s*S*phi*phi)/L
    B=(H0*U_z+(1+dt)*(U-2*z*U*U)/2+d*P_z-2*(1+dt)*z*P)/L
    z_terms=[W0/L*s*psi_s,
        (d*U_z+(1+dt)*(1-4*z*U)/2)/L*psi,
        H0/L*psi_z,-eps*(1-dt)*z/L*M*s*psi_s,
        -eps*d/L*M_z*s*psi_s,-eps*(1+dt)*z/L*psi*psi,
        eps*d/L*psi*psi_z,d/L*C_z,-2*(1+dt)*z/L*C,
        -2*z/L*s*S*phi*phi]
    z_difference=sp.expand(original_z-B-eps*sum(z_terms))
    if theta_difference!=0 or z_difference!=0:
        raise AssertionError('Scaled map does not reproduce the core equations')
    base=Path(__file__).parent
    names=['lei_ren_part1_paper_core_recursion.py',
           'lei_ren_part1_paper_amplitude_factored_core.py',
           'lei_ren_part1_paper_nonlinear_map_majorant.py']
    report=dict(input_hashes={n:hashlib.sha256((base/n).read_bytes()).hexdigest() for n in names},
        angular_symbolic_difference=str(theta_difference),axial_symbolic_difference=str(z_difference),
        angular_term_count=len(theta_terms),axial_term_count=len(z_terms),
        scaled_equations=['2(s Phi_ss+2 Phi_s)+chi Phi=epsilon Etheta',
                          '2(s Psi_ss+Psi_s)=B+epsilon Ez'],
        substitutions=['s=Lambda R','epsilon=1/Lambda','F=F0 Phi',
                       'Uz=U0+epsilon Psi','P=P0+epsilon Pcal',
                       'F0_Z/F0=-Lambda g','g=L H0/(H0²+sigma²)',
                       'M=(integral_0^s Psi)/s','Pcal=F0² integral_0^s Phi²'],
        beta_definition=str(beta),chi_definition=str(chi),
        normalized_equations_exactly_reproduce_core=True,
        coefficient_uniqueness='At radial order n, known coefficients through n determine the next angular/axial rows by division by 2L(n+1)(n+2) and 2L(n+1)²; pressure follows by integration. L is zero-free on the admitted analytic tube.',
        scope='Exact formal equation identity; analytic existence and directed input inclusion are separate dependencies.',
        original_parameter_errors_enclosed=False,matching_certified=False,temporal_recursion=False)
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print('Exact angular and axial scaled identities:0; retained terms:10+10',flush=True)
    return report


if __name__=='__main__':run()
