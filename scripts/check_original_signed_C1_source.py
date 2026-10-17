"""Independent analytic checks for the new fixed-phase source integral C1."""
from pathlib import Path
import sys

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from openai_ns_reconstruction.original_signed_C1_source import (
    correlated_phi_ratio, macro_nonlinear_enclosure, main_coefficient,
    pressure_drive_ordinary, pressure_source_ordinary,
)


def _raises(function):
    try:
        function()
    except ValueError:
        return True
    return False


def run_fixture():
    c=mp.mp.clone(); c.dps=100
    iv=MPIntervalContext();iv.dps=120
    controls=[];differences=[]
    def passed(name,condition):
        if not condition:raise AssertionError(name)
        controls.append({"name":name,"passed":True})
    def b(z):return c.exp(z/4)+z**3/3-c.sin(z)
    delta=c.mpf('.02')
    for text in ('-.3','0','.23'):
        z=c.mpf(text);jet=[c.diff(b,z,k) for k in range(3)]
        for mass_text in ('0','.1249','.4793','.5'):
            mass=c.mpf(mass_text)
            def oracle(q):return mass*((1-delta)*q*b(q)+(1-q*q)*c.diff(b,q))/(1-delta*q*q)
            result=main_coefficient(z,delta,jet,mass,mass_Z_independent=True)
            for k,value in enumerate(result):
                difference=abs(value-c.diff(oracle,z,k));differences.append(difference)
                passed(text+'_'+mass_text+'_ordinary_Z'+str(k),difference<c.mpf('1e-90'))
            boxes=[iv.mpf([c.nstr(v-c.mpf('1e-95'),100),c.nstr(v+c.mpf('1e-95'),100)]) for v in jet]
            truth=main_coefficient(iv.mpf(text),iv.mpf('.02'),boxes,iv.mpf(mass_text),mass_Z_independent=True)
            passed('directed_'+text+'_'+mass_text,all(iv.mpf(c.nstr(p,100)) in t for p,t in zip(result,truth)))
    passed('missing_Bzz_rejected',_raises(lambda:main_coefficient(z,delta,jet[:2],mass,mass_Z_independent=True)))
    passed('unproved_mass_independence_rejected',_raises(lambda:main_coefficient(z,delta,jet,mass,mass_Z_independent=False)))
    passed('None_Bzz_rejected',_raises(lambda:main_coefficient(z,delta,[jet[0],jet[1],None],mass,mass_Z_independent=True)))
    def p0(q):return c.exp(q/3)+q**5/7-c.sin(q)/9
    def ratio(q):return 1+q/7+q*q/11
    for text in ('-.3','0','.23'):
        z=c.mpf(text);pjet=[c.diff(p0,z,k) for k in range(4)]
        qjet=[c.diff(ratio,z,k) for k in range(3)]
        drive=pressure_drive_ordinary(z,delta,pjet)
        source=pressure_source_ordinary(z,delta,pjet,qjet,c.mpf('3.7'))
        def pd(q):return (2*(1+delta)*q*p0(q)-(1-q*q)*c.diff(p0,q))/(2*(1-delta*q*q))
        for label,result,oracle in (('pressure_drive',drive,pd),('pressure_source',source,lambda q:-c.mpf('3.7')*ratio(q)*pd(q))):
            for k,value in enumerate(result):
                difference=abs(value-c.diff(oracle,z,k));differences.append(difference)
                passed(label+'_'+text+'_ordinary_Z'+str(k),difference<c.mpf('1e-90'))
    passed('p0_Z3_missing_rejected',_raises(lambda:pressure_drive_ordinary(z,delta,pjet[:3])))
    z=c.mpf('.23');h=c.mpf('.000001');hcap=iv.mpf('.001')
    def phi(q):return 2+q+q*q/2
    def variation(q):return c.mpf('.2')+q/10+q*q/20
    def exact_ratio(q):return phi(q)/(phi(q)+h*variation(q))
    pjet=[c.diff(phi,z,k) for k in range(3)]
    djet=[c.diff(variation,z,k) for k in range(3)]
    # Deliberately wider Phi uncertainty tests the common-source identity.
    # The same Phi must cancel algebraically before interval evaluation.
    boxes=[iv.mpf([c.nstr(v-c.mpf('.01'),100),c.nstr(v+c.mpf('.01'),100)]) for v in pjet]
    caps=[iv.mpf(c.nstr(abs(v),100))*hcap*iv.mpf('1.1') for v in djet]
    point,truth,proof=correlated_phi_ratio(c,iv,boxes,caps,hcap,iv.log(iv.mpf('.000001')),
        shared_phi_identity_proved=True,cap_linear_width_proved=True,error_threshold='1e-5')
    for k,box in enumerate(truth):
        actual=c.diff(exact_ratio,z,k)
        passed('nonzero_h_shared_Phi_ratio_ordinary_Z'+str(k),iv.mpf(c.nstr(actual,100)) in box)
    passed('ratio_proxy_not_actual_function',abs(exact_ratio(z)-point[0])>c.mpf('1e-10')
           and proof['ratio_approximation_is_exact_source_value'] is False)
    passed('shared_Phi_identity_required',_raises(lambda:correlated_phi_ratio(c,iv,boxes,caps,hcap,iv.log(iv.mpf('.000001')),
        shared_phi_identity_proved=False,cap_linear_width_proved=True,error_threshold='1e-5')))
    passed('linear_width_bound_required',_raises(lambda:correlated_phi_ratio(c,iv,boxes,caps,hcap,iv.log(iv.mpf('.000001')),
        shared_phi_identity_proved=True,cap_linear_width_proved=False,error_threshold='1e-5')))
    passed('unadmitted_ratio_error_cap_rejected',_raises(lambda:correlated_phi_ratio(c,iv,boxes,caps,hcap,iv.log(iv.mpf('.000001')),
        shared_phi_identity_proved=True,cap_linear_width_proved=True,error_threshold='1e-50')))
    # A nonconstant ell creates derivative error even when B_Z=B_ZZ=0.
    # Componentwise relative widening would incorrectly leave those rows zero.
    simple=[iv.mpf(1),iv.mpf(0),iv.mpf(0)]
    truth,proof=macro_nonlinear_enclosure(c,iv,simple,iv.mpf('.1'),iv.mpf('.001'),macro_norm_bound_proved=True)
    for k,box in enumerate(truth):
        actual=c.diff(lambda q:c.exp(c.mpf('.003')*(q-z)),z,k)
        passed('macro_expell_derivative_coupling_Z'+str(k),iv.mpf(c.nstr(actual,100)) in box)
    passed('macro_ell_proof_required',_raises(lambda:macro_nonlinear_enclosure(c,iv,simple,
        iv.mpf('.1'),iv.mpf('.001'),macro_norm_bound_proved=False)))
    return {"passed":True,"control_count":len(controls),"controls":controls,
            "maximum_supporting_difference":c.nstr(max(differences),25),
            "independent_oracle":"mpmath differentiation of explicit exp/sin/polynomial pressure and quotient, C(B)*M, and a shared-Phi ratio with nonzero finite width",
            "accepted_ancestor_or_fixture_run":False,"original_NS_residual":False}


if __name__=='__main__':
    import json
    print(json.dumps(run_fixture(),indent=2))
