"""Independent moderate-parameter checks for the common physical B_y response."""
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import mpmath as mp
from openai_ns_reconstruction.original_reference_native_m15 import evaluate_synthetic_response,root_tangent_proof


def run():
    c=mp.mp.clone();c.dps=55
    E,Ez,p,pz,py,pyz,q,d,phase=map(c.mpf,('1.3','.12','.7','-.23','.21','.08','.6','.9','.31'))
    a=c.mpf(4)/5;K=1/(2*c.pi*(1+2*q*q))
    def H(pp,ph=phase):
        u=pp*q/d;r=u/c.sqrt(1+u*u);alpha=2*q/c.sqrt(1+u*u)
        t=lambda theta:alpha*(c.cos(theta)-r)/(1-2*r*c.cos(theta)+r*r)
        Phi=lambda theta:K*(theta+c.quad(lambda x:t(x)**2,[0,theta]))
        psi=c.findroot(lambda theta:Phi(theta)-ph,2*c.pi*ph,
            df=lambda theta:K*(1+t(theta)**2),solver='newton')
        return -a*c.quad(t,[0,psi])/(4*c.pi)
    got=evaluate_synthetic_response(E=E,Ez=Ez,p=p,pz=pz,py=py,pyz=pyz,q=q,dstar=d,phase=phase,precision=50)
    expected=dict(H=H(p),Hp=c.diff(H,p),Hpp=c.diff(H,p,2))
    tolerance=c.mpf('1e-38')
    for key,value in expected.items():assert abs(value-got[key])<tolerance,(key,str(value-got[key]))
    # Independent differentiation of the actual two-variable family, not the
    # stored G/Gz expression. Phase is held here, then varied explicitly.
    def B(y,z=0,ph=phase):
        return (E+Ez*z)*c.exp(y/10)*H(p+pz*z+(py+pyz*z)*y,ph)
    G=c.diff(lambda y:B(y),0)
    Gz=c.diff(lambda z:c.diff(lambda y:B(y,z),0),0)
    assert abs(G-got['G'])<tolerance
    assert abs(Gz-got['Gz'])<tolerance
    J=c.mpf('3.7')
    native_ref=c.diff(lambda s:B(J*(s-1)),1)
    native_restore=c.diff(lambda t:B(t),0)
    assert abs(native_ref/J-native_restore)<tolerance
    assert abs(native_ref-native_restore)>c.mpf('1e-4')
    assert abs(got['G'])>c.mpf('1e-4')
    wrong_Gz=Ez*got['H']/10+E*got['Hp']*pz/10+Ez*got['Hp']*py+E*got['Hpp']*pz*py
    assert abs(wrong_Gz-got['Gz'])>c.mpf('1e-5')
    N=c.mpf(17)
    total=c.diff(lambda y:B(y,ph=phase+N*y)/N,0)
    Bphi=c.diff(lambda ph:B(0,ph=ph),phase)
    assert abs(total-(G/N+Bphi))<tolerance
    assert abs(total-(G/N+Bphi/N))>c.mpf('1e-4')
    tangent=root_tangent_proof();assert tangent['passed']
    return dict(passed=True,synthetic_parameters_only=True,original_source_values_evaluated=False,
        arbitrary_fast_phase=str(phase),H_Hp_Hpp_independent_inverse_differentiation=True,
        physical_y_and_ordinary_Z1_independent_differentiation=True,
        unequal_J_naturality=True,slow_derivative_nonzero=True,
        omitted_p2_yZ_control_rejected=True,wrong_J_control_rejected=True,
        total_fast_phase_factor_N_control_rejected=True,root_tangents=tangent,
        maximum_absolute_difference=str(max(abs(expected[k]-got[k]) for k in expected)))


if __name__=='__main__':print(json.dumps(run(),indent=2))
