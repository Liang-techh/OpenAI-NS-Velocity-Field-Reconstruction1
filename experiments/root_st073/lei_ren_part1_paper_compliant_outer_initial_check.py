"""Independent normalized five-primitive and interface check for F31."""
# Recomputed for the distinct compliant pressure/moment family.
# Formula origin: lei_ren_part1_paper_shared_outer_initial_check.py; old .01 receipts remain unchanged.

import hashlib
import json
from pathlib import Path
import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_outer_initial import SharedOuterInitial,turnoff_kernels
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_'


def symbolic():
    y,Z,pinv=s.symbols('y Z pinv',real=True);q=1+Z**2
    target=lambda u,V,m,h,k,e,p:[V-m,u-s.Rational(3,2)*h,u*V-s.Rational(3,2)*k,
                                V**2*pinv-u**2/2-e,u**2/2]
    def test(u,V,m,h,k,e,p,rules):
        values=(m,h,k,e,p)
        residuals=[s.simplify((s.diff(v,y)-rhs).subs(rules)) for v,rhs in zip(values,target(u,V,m,h,k,e,p))]
        if any(v!=0 for v in residuals):raise ArithmeticError('Normalized moment ODE failed: '+str(residuals))
        return [str(v) for v in residuals]
    u=s.exp(y/10)/q;V=4*Z;m=V;h=s.Rational(5,8)*u;k=4*Z*h
    e=16*Z**2*pinv-s.Rational(5,12)*u**2;p=s.Rational(5,2)*u**2
    reference=test(u,V,m,h,k,e,p,{})
    J,I1,I2,I3=[s.Function(n)(y) for n in ('J','I1','I2','I3')];sig=s.Function('sig')(y)
    u=s.exp(y/10-s.Rational(3,5)*J)/q;V=4*Z;m=V
    h=(s.Rational(5,8)+I1)*s.exp(-3*y/2)/q;k=4*Z*h
    e=16*Z**2*pinv-(s.Rational(5,12)+I3/2)*s.exp(-y)/q**2
    p=(s.Rational(5,2)+I2/2)/q**2
    rules={s.diff(J,y):sig,s.diff(I1,y):s.exp(8*y/5-s.Rational(3,5)*J),
           s.diff(I2,y):s.exp(y/5-s.Rational(6,5)*J),s.diff(I3,y):s.exp(6*y/5-s.Rational(6,5)*J)}
    slope=test(u,V,m,h,k,e,p,rules)
    u1,m1,h1,k1,e1,p1=s.symbols('u1 m1 h1 k1 e1 p1',real=True)
    K,K2,B=[s.Function(n)(y) for n in ('K','K2','B')];t=y-1
    d=s.exp(-t);d3=s.exp(-3*t/2);dr=s.exp(-t/2)
    u=u1*dr;V=4*Z*B;m=m1*d+4*Z*K
    h=h1*d3+u1*(dr-d3);k=k1*d3+4*Z*u1*dr*K
    e=e1*d+16*Z**2*pinv*K2-u1**2*t*d/2;p=p1+u1**2*(1-d)/2
    axial=test(u,V,m,h,k,e,p,{s.diff(K,y):B-K,s.diff(K2,y):B**2-K2})
    R,Pstar=s.symbols('R Pstar',positive=True)
    # Re-derive the normalized ODE directly from physical cumulative RHS.
    funcs=[s.Function(n)(y) for n in ('m','h','k','e','p')]
    UU,VV=s.symbols('UU VV')
    scales=[R,s.sqrt(2)*R**s.Rational(3,2)*Pstar,s.sqrt(2)*R**s.Rational(3,2)*Pstar,R*Pstar**2,Pstar**2]
    rates=[1,s.Rational(3,2),s.Rational(3,2),1,0]
    physical_rhs=[R*VV,R*s.sqrt(2*R)*Pstar*UU,R*s.sqrt(2*R)*Pstar*UU*VV,
                  R*(VV**2-Pstar**2*UU**2/2),Pstar**2*UU**2/2]
    wanted=target(UU,VV,*funcs)
    physical=[s.simplify(rhs/scale-rate*f-w).subs(pinv,1/Pstar**2)
              for rhs,scale,rate,f,w in zip(physical_rhs,scales,rates,funcs,wanted)]
    if any(s.simplify(v)!=0 for v in physical):raise ArithmeticError('Physical normalization mismatch')
    return dict(reference_five_ODE_residuals=reference,slope_five_ODE_residuals=slope,
                axial_five_ODE_residuals=axial,physical_normalization_residuals=[str(s.simplify(v)) for v in physical])


def overlaps(a,b):
    a0,a1=endpoints(a);b0,b1=endpoints(b);return max(a0,b0)<=min(a1,b1)


def run():
    result=symbolic();field=SharedOuterInitial();c=field.ctx
    raw=json.loads((HERE/(PREFIX+'compliant_outer_initial.json')).read_bytes())
    for name,digest in raw['input_hashes'].items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:raise ValueError('Outer producer dependency changed: '+name)
    keys=('Utheta_over_Pstar','Utheta_y_over_Pstar','Uz','Uz_y','Mz_over_R',
          'Mtheta_over_sqrt2_R_3half_Pstar','Mtheta_z_over_sqrt2_R_3half_Pstar',
          'Mztheta_over_R_Pstar_squared','Mp_over_Pstar_squared','P_over_Pstar_squared')
    interfaces=[]
    with mp.workdps(210):
        for Z in ('-1','0','.5','1'):
            pairs=[('reference/slope',field.reference(Z,0),field.slope(Z,0)),
                   ('slope/turnoff',field.slope(Z,1),field.axial(Z,phase=0)),
                   ('turnoff/buffer',field.axial(Z,phase=1),field.axial(Z,buffer_offset=0))]
            for name,left,right in pairs:
                checks={k:all(overlaps(a,b) for a,b in zip(left[k],right[k])) for k in keys}
                checks['Ur']=overlaps(left['Ur_over_sqrt_R_over_2'],right['Ur_over_sqrt_R_over_2'])
                if not all(checks.values()):raise ArithmeticError('Interface mismatch '+name+' '+Z+str(checks))
                interfaces.append(dict(Z=Z,interface=name,checks=checks))
        whole=field.axial([-1,1],phase='.5',cells=64)
        if not all(len(whole[k])==2 for k in keys):raise ArithmeticError('Whole-axis C1 API failed')
        buffer=field.axial('.5',buffer_offset=11)
        if endpoints(buffer['Uz'][0])!=(mp.mpf(0),mp.mpf(0)):
            raise ArithmeticError('Final buffer is not zero axial velocity')
        if endpoints(buffer['Mz_over_R'][0])[1]<=0:
            raise ArithmeticError('Actual axial moment history was lost')
        K=turnoff_kernels(c,c.exp(40),40,cells=256)
        if endpoints(K['retained_far_tail'])[1]<=0:
            raise ArithmeticError('Positive far kernel tail discarded')
        # Closed form check of the small-y kernel mass against independently
        # integrated exact exponential weight for B<=1, B^2<=B.
        for y in ('1','2','10'):
            K=turnoff_kernels(c,c.mpf(y),40,cells=128)
            for name in ('B_mass','B_squared_mass'):
                if endpoints(K[name])[0]<0 or endpoints(K[name])[1]>endpoints(1-c.exp(-(c.mpf(y)-1)))[1]+mp.mpf('1e-150'):
                    raise ArithmeticError('Normalized positive kernel range failed')
    result.update(interface_checks=interfaces,total_independent_symbolic_identities=20,
        exact_interface_values_and_first_derivative_enclosures_checked=True,
        smooth_all_order_interfaces_follow_from_original_flat_cutoff=True,
        all_order_interface_jets_numerically_evaluated=False,
        whole_axis_C1_outer_initial_API_checked=True,nonzero_axial_history_after_zero_Uz_checked=True,
        positive_far_kernel_tail_retained=True,actual_five_defect_family_sha256=field.family,
        complete_corrected_outer_built=False,heat_exterior_matched=False,temporal_recursion=False,
        input_hashes=dict(raw['input_hashes']))
    for name in (Path(__file__).name,PREFIX+'compliant_outer_initial.json'):
        result['input_hashes'][name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print('Same-family initial outer: 20 physical/moment identities and 12 interfaces PASS; heat and later repairs pending',flush=True)
    return result


if __name__=='__main__':run()
