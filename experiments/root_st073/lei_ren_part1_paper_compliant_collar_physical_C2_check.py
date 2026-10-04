"""Focused physical-collar momentum/tensor/remainder validation."""
import hashlib
import json
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_collar_physical_C2 import physical_collar_identities,physical_source_row
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'


def independent_cartesian_fixture():
    """Full moments, implicit physical coordinates and native differentiation.

    Tensor first jets differentiate the defining full moments under their
    integrals and use the FTC for the moving lower limit. Cartesian
    differentiation checks the moving basis; higher tensor jets do not
    enter divergence. No production coefficient helper supplies these jets.
    """
    with mp.workdps(40):
        a,S,eps,B0=map(mp.mpf,('.15','.004','.002','1.3'))
        delta=2*a; k=1-a; bh=mp.mpf('.5')+a; p=1+delta; b=(1-delta)/2
        cache={}; hcache={}
        def H(x,n=0):
            if not x:return (-1)**n*mp.rf(a,n)*mp.rf(1+a,n)
            key=(x._mpf_,n,mp.mp.prec)
            if key not in hcache:
                hcache[key]=(-1)**n*mp.rf(a,n)*mp.rf(1+a,n)*x**(-1-a-n)*mp.hyperu(1+a+n,2,1/x)
            return hcache[key]
        def sigma(v):
            if v<=0:return mp.mpf(0)
            if v>=1:return mp.mpf(1)
            return 1/(1+mp.exp(1/v**2-1/(1-v)**2))
        def K(v,z,n=0):
            sig=sigma(v); phi=mp.exp(-4/(3-v)**2) if v<3 else mp.mpf(0)
            x=2*(1-z*z)*S*mp.exp(-v)
            if n==0:return (1-sig)*(1-eps)+sig*(1-eps*phi)*H(x)
            xz=-4*z*S*mp.exp(-v)
            if n==1:return sig*(1-eps*phi)*H(x,1)*xz
            if n==2:return sig*(1-eps*phi)*(H(x,2)*xz*xz-4*S*mp.exp(-v)*H(x,1))
            raise ValueError('Independent K derivative order unsupported')
        def square(v,z,n):
            if n==0:return K(v,z,0)**2
            if n==1:return 2*K(v,z,0)*K(v,z,1)
            if n==2:return 2*K(v,z,1)**2+2*K(v,z,0)*K(v,z,2)
            raise ValueError('Independent squared K derivative order unsupported')
        def angular_gamma(z):
            x=2*(1-z*z)*S*mp.exp(-3)
            return (H(x)+x*((1+a)*H(x)+x*H(x,1)))/k
        def moment(v,z,label,n=0):
            key=(v._mpf_,z._mpf_,label,n,mp.mp.prec)
            if key in cache:return cache[key]
            breaks=[v]+([mp.mpf(1)] if v<1 else [])+[mp.mpf(3)]
            if label=='A':
                value=mp.quad(lambda q:mp.exp(k*(q-v))*((1-K(q,z)) if n==0 else -K(q,z,1)),breaks)
                value+=mp.exp(k*(3-v))*(mp.diff(angular_gamma,z,n)-(1/k if n==0 else 0))
                value+=1/k if n==0 else 0
            else:
                rate=delta if label=='E' else p; factor=1 if label=='E' else mp.mpf('.5')
                value=mp.quad(lambda q:mp.exp(-rate*q)*square(q,z,n)*factor,breaks)
                def tail(q):
                    if not q:return factor/rate if n==0 else mp.mpf(0)
                    return square(3-mp.log(q)/rate,z,n)*factor/rate
                value+=mp.exp(-3*rate)*mp.quad(tail,[0,1]); value*=mp.exp(rate*v)
            cache[key]=value; return value
        def coordinates(r,z,tau,nu):
            zs=z/mp.sqrt(nu)
            lam=mp.findroot(lambda q:q*q-zs*zs*q**(2*delta)-tau,mp.sqrt(tau)+abs(zs),verify=True)
            Z=zs/lam**(1-delta); R=r*r/(2*nu*lam*lam)
            return lam,mp.log(R*S),Z,R
        def swirl(r,z,tau,nu):
            lam,v,Z,R=coordinates(r,z,tau,nu)
            return mp.sqrt(nu)*lam**(-1-delta)*B0*mp.exp(-bh*v)*K(v,Z)
        def pressure(r,z,tau,nu):
            lam,v,Z,R=coordinates(r,z,tau,nu)
            return -nu*lam**(-2-2*delta)*B0*B0*mp.exp(-p*v)*moment(v,Z,'P')
        def pressure_gradient(r,z,tau,nu):
            # The full pressure FTC was independently accepted by the C4
            # pressure checker. Reuse it, and differentiate the implicit
            # coordinate map natively instead of re-differentiating its
            # quadrature at twice the working precision for every axis.
            lam,v,Z,R=coordinates(r,z,tau,nu); P=moment(v,Z,'P'); Pz=moment(v,Z,'P',1)
            Py=p*P-K(v,Z)**2/2; eta=-2-2*delta; B2=B0*B0*mp.exp(-p*v)
            result=[]
            for axis,base in ((0,r),(1,z)):
                def mapped(q,index):
                    return coordinates(q,z,tau,nu)[index] if axis==0 else coordinates(r,q,tau,nu)[index]
                dl,dv,dZ=[mp.diff(lambda q:mapped(q,index),base) for index in range(3)]
                result.append(-nu*lam**eta*B2*(eta*P*dl/lam+(Py-p*P)*dv+Pz*dZ))
            return result
        def stress_first_jets(r,z,tau,nu):
            lam,v,Z,R=coordinates(r,z,tau,nu); B=B0*mp.exp(-bh*v); L=1-delta*Z*Z; d=1-Z*Z
            A=moment(v,Z,'A'); Az=moment(v,Z,'A',1)
            E=moment(v,Z,'E'); Ez=moment(v,Z,'E',1); P=moment(v,Z,'P'); Pz=moment(v,Z,'P',1)
            Ezz=moment(v,Z,'E',2); Pzz=moment(v,Z,'P',2)
            K0=K(v,Z); Kz=K(v,Z,1); Ky=mp.diff(lambda q:K(q,Z),v); Kyy=mp.diff(lambda q:K(q,Z),v,2)
            Qt=mp.sqrt(R/2)*B; Qz=mp.sqrt(R/2)*B*B
            It=(k*A-b*Z*Az-K0)/L; C=2/R*(Ky-(1+a)*K0)
            Cy=2/R*(Kyy-(2+a)*Ky+(1+a)*K0)
            theta_y=-It-(Ky+b*Z*Kz)/L+Cy-a*C
            N=delta*Z*E-d*Ez/2-2*p*Z*P+d*Pz; J=N/L
            axial_y=-J/2+(d*Pz-2*Z*(p*P-K0*K0/2))/L
            Nz=delta*E+(delta+1)*Z*Ez-d*Ezz/2-2*p*P-(2*p+2)*Z*Pz+d*Pzz
            Jz=Nz/L+N*2*delta*Z/(L*L); beta=-2-delta
            scale=nu*lam**beta
            Ttheta=scale*Qt*(It+C); Tz=scale*Qz*J
            Ttheta_r=scale*Qt*2/r*theta_y; Tz_r=scale*Qz*2/r*axial_y
            Tz_z=mp.sqrt(nu)*lam**(beta+delta-1)*Qz*(beta*Z*J+d*Jz-2*Z*axial_y)/L
            return Ttheta,Tz,Ttheta_r,Tz_r,Tz_z
        residuals=[]; divergences=[]; radial=[]; nonzero=0
        for offset,nu in ((mp.mpf('.5'),mp.mpf('.01')),(mp.mpf('2'),mp.mpf('.7'))):
            print('Independent Cartesian collar fixture offset '+str(offset)+' nu '+str(nu),flush=True)
            Z=mp.mpf('.4'); tau=mp.mpf('.7'); lam=mp.sqrt(tau/(1-Z*Z)); angle=mp.mpf('.6')
            r=mp.sqrt(nu)*lam*mp.sqrt(2*mp.exp(offset)/S); z=mp.sqrt(nu)*lam**(1-delta)*Z
            point=[r*mp.cos(angle),r*mp.sin(angle),z]
            def vel(x,y,z,ta=tau):
                rr=mp.sqrt(x*x+y*y); g=swirl(rr,z,ta,nu)
                return [-y/rr*g,x/rr*g,mp.mpf(0)]
            def pres(x,y,z):return pressure(mp.sqrt(x*x+y*y),z,tau,nu)
            Tt,Tz,Ttr,Tzr,Tzz=stress_first_jets(r,z,tau,nu)
            print('Complete future stress moment first jets recovered',flush=True)
            pr,pz=pressure_gradient(r,z,tau,nu)
            pressure_cart=[mp.cos(angle)*pr,mp.sin(angle)*pr,pz]
            # Complete the independently obtained first stress jets. The
            # unneeded derivatives of r*Tz_z cancel structurally in divergence.
            def tensor(x,y,zz):
                rr=mp.sqrt(x*x+y*y); cs=x/rr; sn=y/rr
                # Ttheta_z never enters this tensor's Cartesian divergence.
                tt=Tt+Ttr*(rr-r); tz=Tz+Tzr*(rr-r)+Tzz*(zz-z); dd=rr*Tzz
                return [[-2*cs*sn*tt+sn*sn*dd,(cs*cs-sn*sn)*tt-cs*sn*dd,cs*tz],
                        [(cs*cs-sn*sn)*tt-cs*sn*dd,2*cs*sn*tt+cs*cs*dd,sn*tz],[cs*tz,sn*tz,mp.mpf(0)]]
            def partial(fn,j,n=1):
                def evaluate(q):
                    args=list(point); args[j]=q; return fn(*args)
                return mp.diff(evaluate,point[j],n)
            u=vel(*point); gzz=mp.diff(lambda q:swirl(r,q,tau,nu),z,2)
            Etheta=-nu*gzz; E=[-mp.sin(angle)*Etheta,mp.cos(angle)*Etheta,mp.mpf(0)]
            nonzero+=abs(Etheta)>mp.mpf('1e-20')
            cart_res=[]
            for i in range(3):
                ut=-mp.diff(lambda ta:vel(*point,ta)[i],tau)
                conv=sum(u[j]*partial(lambda x,y,z:vel(x,y,z)[i],j) for j in range(3))
                lap=sum(partial(lambda x,y,z:vel(x,y,z)[i],j,2) for j in range(3))
                momentum=ut+conv+pressure_cart[i]-nu*lap
                div=sum(partial(lambda x,y,z:tensor(x,y,z)[i][j],j) for j in range(3))
                error=momentum+div-E[i]
                if abs(error)>mp.mpf('1e-29'):raise ArithmeticError('Independent Cartesian stress/remainder decomposition failed: '+str((offset,nu,i,error)))
                residuals.append(mp.nstr(error,12)); cart_res.append(momentum)
            divu=sum(partial(lambda x,y,z:vel(x,y,z)[j],j) for j in range(3))
            radial_res=mp.cos(angle)*cart_res[0]+mp.sin(angle)*cart_res[1]
            if abs(divu)>mp.mpf('1e-29') or abs(radial_res)>mp.mpf('1e-29'):raise ArithmeticError('Independent divergence/radial momentum failed')
            divergences.append(mp.nstr(divu,12)); radial.append(mp.nstr(radial_res,12))
        if nonzero!=2:raise ArithmeticError('Independent fixture did not exercise nonzero axial viscosity')
        return dict(Cartesian_stress_remainder_identity_components_checked=6,
                    Cartesian_decomposition_errors=residuals,divergence_errors=divergences,radial_momentum_errors=radial,
                    nonzero_axial_viscosity_remainders_exercised=int(nonzero),
                    stress_jets_from_full_integral_derivatives_and_FTC=True,
                    pressure_gradient_from_accepted_full_integral_FTC_and_native_implicit_map=True,
                    full_future_moments_not_radially_truncated=True, moderate_fixture_only=True,passed=True)


def run():
    name=PREFIX+'collar_physical_C2.json'; record=json.loads((HERE/name).read_bytes()); hashes=dict(record['input_hashes'])
    for source,digest in hashes.items():
        if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:raise ValueError('Physical collar dependency changed: '+source)
    proof=physical_collar_identities()
    if proof!=record['physical_collar_identities']:raise ValueError('Physical collar exact identity changed')
    c=MPIntervalContext(); c.dps=280; finite=zeros=0
    def check_row(row):
        nonlocal finite,zeros
        lo,hi=endpoints(read_interval(c,row['signed_coefficient']))
        if not mp.isfinite(lo) or not mp.isfinite(hi):raise ArithmeticError('Nonfinite physical signed coefficient')
        if row['exact_zero']:
            if lo or hi or row['log_absolute_upper'] is not None:raise ArithmeticError('Invalid exact physical zero row')
            zeros+=1
        else:
            if not all(mp.isfinite(v) for v in endpoints(read_interval(c,row['log_absolute_upper']))):raise ArithmeticError('Nonfinite physical logarithmic bound')
            finite+=1
    for point in record['samples']+[record['whole_source_collar']]:
        if not point['actual_regional_physical_collar_stress_remainder_identity_verified']:raise ValueError('Regional physical identity missing')
        for flag in ('global_admissible_stress_lift_constructed','collar_cone_certified','whole_outer_cone_certified',
                     'independently_bounded_global_flat_remainder','physical_energy_integral_certified','full_background_NS_validation','temporal_recursion'):
            if point[flag]:raise ValueError('Physical collar scope overclaimed: '+flag)
        for flag in ('exact_physical_radial_momentum_residual','exact_completed_tensor_radial_divergence',
                     'exact_physical_radial_remainder','exact_physical_axial_remainder','exact_physical_divergence'):
            if endpoints(read_interval(c,point[flag]))!=(0,0):raise ArithmeticError('Nonzero structural pure-swirl component: '+flag)
        for key,count in (('physical_cylindrical_stress_mixed3',10),('physical_cylindrical_stress_divergence_mixed2',6)):
            for grid in point[key].values():
                if len(grid)!=count:raise ValueError('Physical derivative grid incomplete')
                for row in grid.values():check_row(row)
        for key in ('completed_theta_theta_stress_mixed2','physical_angular_axial_viscosity_remainder_mixed2'):
            if len(point[key])!=6:raise ValueError('Physical mixed2 completion/remainder missing')
            for row in point[key].values():check_row(row)
        for key in ('completed_background_stress_tensor_cartesian_components','physical_completed_stress_divergence_cartesian',
                    'physical_remainder_cartesian','physical_momentum_residual_decomposition_cartesian'):
            for rows in point[key].values():
                for row in rows:check_row(row)
        if point['Gamma_endpoint_physical_remainder_exact_zero']:
            if any(not row['exact_zero'] for row in point['physical_angular_axial_viscosity_remainder_mixed2'].values()):
                raise ValueError('Same-source Gamma endpoint remainder not zero through mixed2')
    # Check the stress/remainder viscosity derivative factors without
    # materializing a huge amplitude or a source cap as an exact value.
    units=0
    for base in (1,mp.mpf('.5')):
        for order in range(4):
            row=physical_source_row(c,c.mpf(2),dict(example=c.mpf(3)),c.mpf(-4),c.mpf(-1),c.mpf('.01'),order,nu_base=c.mpf(base))
            if endpoints(row['physical_viscosity_exponent'])!=endpoints(c.mpf(base)-c.mpf(order)/2):raise ArithmeticError('Physical derivative viscosity power wrong')
            units+=1
    fixture=independent_cartesian_fixture()
    hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest(); hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result=dict(actual_five_defect_family_sha256=record['actual_five_defect_family_sha256'],implicit_source_sha256=record['implicit_source_sha256'],
                all_passed=True,actual_regional_physical_collar_stress_remainder_identity_verified=True,
                actual_finite_physical_signed_rows_checked=finite,actual_exact_zero_rows_checked=zeros,
                viscosity_derivative_unit_rows_checked=units,physical_collar_identities=proof,
                independent_cartesian_fixture=fixture,global_admissible_stress_lift_constructed=False,
                independently_bounded_global_flat_remainder=False,collar_cone_certified=False,temporal_recursion=False,input_hashes=hashes)
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('Actual physical collar completed stress/remainder and independent Cartesian fixture PASS',flush=True)
    return result


if __name__=='__main__':run()
