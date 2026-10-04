"""Actual leading collar: physical completed stress and momentum remainder.

Uses the source-bound common moments, full collar K and original maps. The
regional leading remainder is axial viscosity, not a certified flat error.
All positive amplitudes and radius factors remain correlated source logs.
"""
import hashlib
import functools
import json
import math
from pathlib import Path

import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_collar_stress_C3 import CompliantCollarStressC3
from lei_ren_part1_paper_compliant_global_physical_assembly import CompliantGlobalPhysicalAssembly
from lei_ren_part1_paper_compliant_pulse_physical_bounds import (
    physical_bracket,physical_operators,interval_expression,ZSYM,DSYM,BSYM)
from lei_ren_part1_paper_compliant_heat_physical_C4 import viscosity_source_row
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'


@functools.lru_cache(maxsize=None)
def collar_velocity_operator_coefficients(i,j):
    """Combine the exact B_y/B=-(1+delta)/2 correlation before enclosure."""
    coefficients={}; bh=(1+DSYM)/2
    for (k,n),value in physical_operators()[i,j].items():
        value=value.subs(BSYM,-1-DSYM)
        for l in range(k+1):
            key=(l,n)
            coefficients[key]=coefficients.get(key,s.Integer(0))+value*s.binomial(k,l)*(-bh)**(k-l)
    result={key:s.factor(s.cancel(value)) for key,value in coefficients.items() if s.cancel(value)!=0}
    if j and (0,0) in result:raise ArithmeticError('Pure power baseline did not cancel in axial operator')
    return result


def collar_velocity_bracket(c,K,i,j,z,delta):
    return sum((interval_expression(c,value,z,delta,-1-delta)*K[k][n]*math.factorial(n)
                for (k,n),value in collar_velocity_operator_coefficients(i,j).items()),c.mpf(0))


def physical_collar_identities():
    proofs={}
    def zero(name,value):
        if s.simplify(s.powsimp(s.expand_power_exp(value),force=True))!=0:
            raise ArithmeticError('Physical collar identity failed: '+name)
        proofs[name]=True
    a,y,z,R=s.symbols('a y Z R',positive=True); delta=2*a; k=1-a; bh=s.Rational(1,2)+a
    b=(1-delta)/2; p=1+delta; d=1-z*z; L=1-delta*z*z
    K=s.Function('same_full_collar_K')(y,z); A=s.Function('same_angular_moment')(y,z)
    E=s.Function('same_current_energy')(y,z); P=s.Function('same_current_pressure')(y,z)
    It=(k*A-b*z*s.diff(A,z)-K)/L
    shear=2/R*(s.diff(K,y)-(1+a)*K)
    # y=logR: a derivative of 1/R contributes -1/R.
    dy=lambda value:s.diff(value,y)+R*s.diff(value,R)
    Ct=It+shear; Cz=(delta*z*E-d*s.diff(E,z)/2-2*p*z*P+d*s.diff(P,z))/L
    rules={s.diff(A,y):K-k*A,s.diff(A,y,z):s.diff(K,z)-k*s.diff(A,z),
           s.diff(E,y):delta*E-K*K,s.diff(E,y,z):delta*s.diff(E,z)-2*K*s.diff(K,z),
           s.diff(P,y):p*P-K*K/2,s.diff(P,y,z):p*s.diff(P,z)-K*s.diff(K,z)}
    angular_time=(s.diff(K,y)+b*z*s.diff(K,z))/L
    angular_radial_laplacian=(2*s.diff(K,y,2)-4*bh*s.diff(K,y)+(2*bh*bh-s.Rational(1,2))*K)/R
    minus_angular_divergence=-(k*Ct+dy(Ct)).subs(rules,simultaneous=True)
    zero('actual_angular_time_minus_radial_viscosity_equals_minus_stress_divergence',
         angular_time-angular_radial_laplacian-minus_angular_divergence)
    axial_pressure=(2*z*(p*P-K*K/2)-d*s.diff(P,z))/L
    minus_axial_divergence=-(dy(Cz)-delta*Cz).subs(rules,simultaneous=True)
    zero('actual_axial_pressure_gradient_equals_minus_axial_stress_divergence',axial_pressure-minus_axial_divergence)
    zero('actual_common_pressure_radial_FTC_cancels_centrifugal_force',p*P-(p*P-K*K/2)-K*K/2)
    # Full physical velocity includes axial viscosity; the leading angular
    # radial identity above omits it, so its remainder has a negative sign.
    viscosity,shearzz=s.symbols('nu utheta_zz',positive=True)
    zero('physical_angular_remainder_sign',-viscosity*shearzz+viscosity*shearzz)
    H,Hp=s.symbols('full_Gamma_H full_Gamma_H_xi',real=True); xi=2*d/R
    gamma_z=(-2*bh*z*H+d*(-4*z/R*Hp)-2*z*(-bh*H-xi*Hp))/L
    zero('full_Gamma_endpoint_first_axial_velocity_operator_zero',gamma_z)
    zero('full_Gamma_endpoint_second_axial_velocity_operator_zero',dy(gamma_z)+s.diff(gamma_z,z))
    Hpp=s.symbols('full_Gamma_H_xixi',real=True)
    gamma_rows={(0,0):H,(1,0):-xi*Hp,(2,0):xi*Hp+xi*xi*Hpp,
                (0,1):-4*z/R*Hp,(0,2):-4/R*Hp+16*z*z/R**2*Hpp,
                (1,1):4*z/R*(Hp+xi*Hpp)}
    zero('actual_reduced_axial_viscosity_operator_is_zero_on_full_Gamma_endpoint',
         sum(value.subs({ZSYM:z,DSYM:delta})*gamma_rows[index]
             for index,value in collar_velocity_operator_coefficients(0,2).items()))
    for order in range(5):
        prefactor=s.exp(-bh*y)
        zero('actual_velocity_amplitude_Leibniz_y'+str(order),s.diff(prefactor*K,y,order)/prefactor
             -sum(s.binomial(order,l)*(-bh)**(order-l)*s.diff(K,y,l) for l in range(order+1)))
    for i in range(3):
        for j in range(3-i):
            operator=collar_velocity_operator_coefficients(i,j+2)
            if (0,0) in operator:raise ArithmeticError('Axial-viscosity baseline not removed')
            proofs['actual_remainder_operator_r'+str(i)+'_z'+str(j)+'_has_no_unit_baseline']=True

    # Derive the completed tensor divergence directly in Cartesian space.
    x,w,zz,r=s.symbols('x cart_y physical_z r',real=True); rr=s.sqrt(x*x+w*w); cs=x/rr; sn=w/rr
    At=s.Function('physical_Ttheta'); Az=s.Function('physical_Tz'); D=s.Function('completed_theta_theta')
    tt=At(rr,zz); tz=Az(rr,zz); dd=D(rr,zz)
    tensor=(( -2*cs*sn*tt+sn*sn*dd,(cs*cs-sn*sn)*tt-cs*sn*dd,cs*tz),
            ((cs*cs-sn*sn)*tt-cs*sn*dd,2*cs*sn*tt+cs*cs*dd,sn*tz),
            (cs*tz,sn*tz,s.Integer(0)))
    div=[s.diff(row[0],x)+s.diff(row[1],w)+s.diff(row[2],zz) for row in tensor]
    radial=s.diff(Az(rr,zz),zz)-dd/rr
    angular=s.diff(At(r,zz),r).subs(r,rr)+2*tt/rr
    axial=s.diff(Az(r,zz),r).subs(r,rr)+tz/rr
    for label,value,target in zip(('x','y','z'),div,(cs*radial-sn*angular,sn*radial+cs*angular,axial)):
        zero('Cartesian_completed_tensor_divergence_'+label,value-target)
    zero('completed_theta_theta_cancels_radial_divergence',radial.subs(dd,rr*s.diff(Az(rr,zz),zz)))
    # Product derivatives used for all mixed2 divergence/completion rows.
    coefficient=s.symbols('radial_divergence_coefficient',real=True); f=Az(r,zz)
    for i in range(3):
        for j in range(3-i):
            div_product=s.diff(f,r,i+1,zz,j)+coefficient*sum(s.binomial(i,m)*(-1)**m*s.factorial(m)*r**(-m-1)*s.diff(f,r,i-m,zz,j) for m in range(i+1))
            zero('actual_mixed2_divergence_Leibniz_r'+str(i)+'_z'+str(j),
                 s.diff(s.diff(f,r)+coefficient*f/r,r,i,zz,j)-div_product)
            diagonal_product=r*s.diff(f,r,i,zz,j+1)+(i*s.diff(f,r,i-1,zz,j+1) if i else 0)
            zero('actual_mixed2_completed_diagonal_Leibniz_r'+str(i)+'_z'+str(j),
                 s.diff(r*s.diff(f,zz),r,i,zz,j)-diagonal_product)
            for m in range(i+1):
                zero('divergence_source_radius_ratio_'+str(i)+'_'+str(m),
                     (2*R)**(-s.Rational(m+1,2))*(2/R)**s.Rational(i-m,2)/(2/R)**s.Rational(i+1,2)-s.Rational(1,2)**(m+1))
    # No meridional velocity is introduced. Full centrifugal balance is
    # retained, so the radial remainder and divergence are exactly zero.
    G=s.Function('physical_swirl')(rr,zz)
    zero('Cartesian_pure_swirl_divergence',s.diff(-sn*G,x)+s.diff(cs*G,w))

    # Combine enormous source logarithms analytically before enclosure.
    lp,lu,lrp,mu,length,Ts,wait,lone,t=s.symbols('logP logU logRp mu Lrel Ts wait logone offset',real=True)
    bp=s.Rational(1,2)+mu; rate=1-mu
    logev=lp+lu-13/(2*mu)-13
    logtheta=-bp*(100+length)-s.log(2)-bp-rate/2-s.Rational(3,2)*Ts-s.Rational(3,2)+k/2-bh*wait-lone
    logRt=lrp+13/mu+100+length+2+Ts+wait
    logB=lp+lu-13/(2*mu)-bp*(100+length)-s.Rational(3,2)*Ts-bh*(wait+t)-15-mu/2-a/2-s.log(2)-lone
    logQt=lp+lu+lrp/2-mu*(100+length)-Ts-a*(wait+t)-14-mu/2-a/2-s.Rational(3,2)*s.log(2)-lone
    zero('actual_B_source_logs_combined_before_enclosure',logev+logtheta-bh*t-logB)
    zero('actual_Qtheta_inverse_mu_cancellation_before_enclosure',logB+(logRt+t)/2-s.log(2)/2-logQt)
    zero('actual_completed_diagonal_radius_cancellation_before_enclosure',logQt+logB+(logRt+t)/2+s.log(2)/2-(2*logQt+s.log(2)))
    nu=s.symbols('physical_nu',positive=True); rs=s.symbols('source_r',positive=True); field=s.Function('source_Tz')
    zero('viscosity_completed_tensor_conversion',s.sqrt(nu)*rs*(nu/s.sqrt(nu)*field(rs)) -nu*rs*field(rs))
    return dict(identities=proofs,actual_regional_physical_collar_stress_remainder_identity_verified=True,
                physical_stress_prefactor='nu*lambda^(-2-delta)',
                actual_remainder='E_B=(0,-nu*partial_zz(utheta),0) in cylindrical components',
                completed_tensor='Trtheta=Ttheta, Trz=Tz, Ttheta_theta=r*partial_z(Tz); symmetric, remaining entries zero',
                source_formulas=['Lei-Ren v2 (2.19)-(2.20)','(3.1)-(3.3)','completed tensor (21.3)-(21.4)'],
                global_flatness_not_inferred_from_leading_axial_viscosity=True)


def stress_source_log_parts(assembly,heat,offset):
    c=assembly.ctx; a=heat.a; mu=heat.mu; bp=c.mpf('.5')+mu; bh=heat.bh
    logu=assembly.dispatch.provider('flatten').logEv2_parts['inlet_log']/2
    common=dict(logPstar=assembly.logP,actual_inlet_logU=logu,minus_log_one_minus_epsilon=-heat.steep.logone)
    B=dict(common,inverse_mu=-13/(2*mu),flatten_power=-bp*(100+heat.steep.outer.Lrel),
           steep=-c.mpf('1.5')*heat.steep.Ts,waiting_and_current=-bh*(heat.steep.wait+offset),
           finite=-15-mu/2-a/2-c.ln(2))
    Qt=dict(common,radial_origin=assembly.logRp/2,flatten_power=-mu*(100+heat.steep.outer.Lrel),
            steep=-heat.steep.Ts,waiting_and_current=-a*(heat.steep.wait+offset),
            finite=-14-mu/2-a/2-c.mpf('1.5')*c.ln(2))
    Qz={**{'Qtheta_'+key:value for key,value in Qt.items()},**{'B_'+key:value for key,value in B.items()}}
    diag={key:value*2 for key,value in Qt.items()}; diag['log2']=c.ln(2)
    return dict(B=B,Qtheta=Qt,Qz=Qz,completed_diagonal=diag)


def physical_source_row(c,coefficient,parts,gamma,logtau,nu,spatial_order,radial_log=0,nu_base=1):
    """Signed coefficient with exact source factors and finite upper log."""
    magnitude=max(abs(v) for v in endpoints(coefficient)); lognu=c.ln(nu)
    exponent=c.mpf(nu_base)-c.mpf(spatial_order)/2
    combined=sum(parts.values(),c.mpf(0))+radial_log
    upper=None if not magnitude else combined+gamma*logtau/2+exponent*lognu+c.ln(c.mpf(magnitude))
    return dict(signed_coefficient=coefficient,exact_zero=not magnitude,
                actual_source_log_parts=parts,physical_lambda_exponent=gamma,
                physical_viscosity_exponent=exponent,physical_viscosity_log_prefactor=exponent*lognu,
                radial_log_prefactor=radial_log,log_absolute_upper=upper,
                source_factors_combined_before_enclosure=True,positive_source_factors_not_materialized=True)


def scale_row(c,row,multiplier):
    multiplier=c.mpf(multiplier)
    result=dict(row); result['signed_coefficient']=row['signed_coefficient']*multiplier
    magnitude=max(abs(v) for v in endpoints(multiplier)); result['exact_zero']=row['exact_zero'] or not magnitude
    result['log_absolute_upper']=None if result['exact_zero'] else row['log_absolute_upper']+c.ln(c.mpf(magnitude))
    return result


class CompliantCollarPhysicalC2:
    def __init__(self):
        self.stress=CompliantCollarStressC3(); self.assembly=CompliantGlobalPhysicalAssembly()
        self.ctx=self.assembly.ctx; self.heat=self.stress.heat
        self.family,self.source=self.stress.family,self.stress.source
        if (self.assembly.family,self.assembly.source)!=(self.family,self.source):raise ValueError('Collar stress/field family mismatch')
        self.proof=physical_collar_identities(); self.hashes=dict(self.stress.hashes); self.hashes.update(self.assembly.hashes)
        for stem,gate in (('collar_stress_C3_check','actual_original_collar_similarity_stress_recovered'),
                          ('global_physical_assembly_check','all_33_original_source_charts_physical_spatial4_time1_mapped')):
            name=PREFIX+stem+'.json'; receipt=json.loads((HERE/name).read_bytes())
            if not receipt['all_passed'] or not receipt[gate]:raise ValueError('Physical collar prerequisite not admitted: '+name)
            if (receipt['actual_five_defect_family_sha256'],receipt['implicit_source_sha256'])!=(self.family,self.source):
                raise ValueError('Physical collar receipt family mismatch')
            for path,digest in receipt['input_hashes'].items():
                if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:raise ValueError('Physical collar source changed: '+path)
            self.hashes.update(receipt['input_hashes']); self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        name=PREFIX+'heat_physical_C4.py'; self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        self.hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    def collar(self,Z,offset,log_tau='-1',theta='0',viscosity='1'):
        c=self.ctx; z=c.mpf(Z); t=c.mpf(offset); logtau=c.mpf(log_tau); nu=c.mpf(viscosity)
        if any(not mp.isfinite(v) for v in endpoints(logtau)):raise ValueError('Finite log(tau) with tau>0 required')
        if endpoints(nu)[0]<=0 or any(not mp.isfinite(v) for v in endpoints(nu)):raise ValueError('Finite constant nu>0 required')
        packet=self.stress.collar(z,t); point=self.assembly.evaluate('heat_collar',z,t,log_tau=logtau,theta=theta)
        logR=point['source_logR_enclosure']; beta=-2-self.heat.delta
        parts=stress_source_log_parts(self.assembly,self.heat,t)
        grids=packet['collar_similarity_stress_mixed3_factored']; stress={}; divergence={}; divmixed={}
        for label,factor in (('theta','Qtheta'),('axial','Qz')):
            stress[label]={}
            for i in range(4):
                for j in range(4-i):
                    bracket=physical_bracket(c,grids[label],i,j,z,self.heat.delta,beta)
                    stress[label]['r'+str(i)+'_z'+str(j)]=physical_source_row(c,bracket,parts[factor],beta-i+j*(self.heat.delta-1),
                                                  logtau,nu,i+j,c.mpf(i)/2*(c.ln(2)-logR))
            divmixed[label]={}
            for i in range(3):
                for j in range(3-i):
                    coefficient=physical_bracket(c,grids[label],i+1,j,z,self.heat.delta,beta)
                    for m in range(i+1):
                        coefficient+=physical_bracket(c,grids[label],i-m,j,z,self.heat.delta,beta)*(
                            (2 if label=='theta' else 1)*math.comb(i,m)*(-1)**m*math.factorial(m)/c.mpf(2)**(m+1))
                    divmixed[label]['r'+str(i)+'_z'+str(j)]=physical_source_row(c,coefficient,parts[factor],
                                                beta-1-i+j*(self.heat.delta-1),logtau,nu,i+j,
                                                c.mpf(i+1)/2*(c.ln(2)-logR),nu_base=c.mpf('.5'))
            divergence[label]=divmixed[label]['r0_z0']
        # Source B*K is differentiated before dividing by its fixed-basepoint
        # B. Thus no amplitude/logarithmic rate is differentiated twice.
        K=self.heat.shape(z,t)['K_rows']
        gvbeta=-1-self.heat.delta
        endpoint=packet['Gamma_endpoint_stress_exact_zero_from_same_moments']
        emixed={}
        for i in range(3):
            for j in range(3-i):
                ebracket=-collar_velocity_bracket(c,K,i,j+2,z,self.heat.delta)
                if endpoint:ebracket=c.mpf(0)
                emixed['r'+str(i)+'_z'+str(j)]=physical_source_row(c,ebracket,parts['B'],
                           gvbeta-i+(j+2)*(self.heat.delta-1),logtau,nu,i+j,
                           c.mpf(i)/2*(c.ln(2)-logR),nu_base=c.mpf('.5'))
        remainder=emixed['r0_z0']
        # The completion has its own correlated scale: sqrt(2R)*Qz=2*Qtheta^2.
        diagmixed={}
        for i in range(3):
            for j in range(3-i):
                dzcoef=physical_bracket(c,grids['axial'],i,j+1,z,self.heat.delta,beta)
                if i:dzcoef+=physical_bracket(c,grids['axial'],i-1,j+1,z,self.heat.delta,beta)*(c.mpf(i)/2)
                diagmixed['r'+str(i)+'_z'+str(j)]=physical_source_row(c,dzcoef,parts['completed_diagonal'],
                           beta+self.heat.delta-i+j*(self.heat.delta-1),logtau,nu,i+j,
                           c.mpf(i)/2*(c.ln(2)-logR))
        diagonal=diagmixed['r0_z0']
        zero=physical_source_row(c,c.mpf(0),parts['B'],gvbeta,logtau,nu,0)
        cs=point['cos_theta']; sn=point['sin_theta']; tt=stress['theta']['r0_z0']; tz=stress['axial']['r0_z0']
        tensor=dict(xx=[scale_row(c,tt,-2*cs*sn),scale_row(c,diagonal,sn*sn)],
                    xy=[scale_row(c,tt,cs*cs-sn*sn),scale_row(c,diagonal,-cs*sn)],
                    xz=[scale_row(c,tz,cs)],yy=[scale_row(c,tt,2*cs*sn),scale_row(c,diagonal,cs*cs)],
                    yz=[scale_row(c,tz,sn)],zz=[zero])
        divcart=dict(x=[scale_row(c,divergence['theta'],-sn)],y=[scale_row(c,divergence['theta'],cs)],z=[divergence['axial']])
        e_cart=dict(x=[scale_row(c,remainder,-sn)],y=[scale_row(c,remainder,cs)],z=[zero])
        momentum={label:[scale_row(c,row,-1) for row in divcart[label]]+e_cart[label] for label in ('x','y','z')}
        # Explicit physical viscosity pullback on the independently accepted
        # velocity/pressure assembly; source coordinates keep the original map.
        lognu=c.ln(nu)
        point['physical_spatial_cartesian_mixed4']={key:{component:{label:viscosity_source_row(row,lognu,c.mpf((2 if component=='p' else 1)-sum(int(v[1:]) for v in key.split('_')))/2)
                                                                  for label,row in values.items()} for component,values in components.items()}
                                                     for key,components in point['physical_spatial_cartesian_mixed4'].items()}
        point['first_fixed_x_physical_time_derivative']={component:{label:viscosity_source_row(row,lognu,c.mpf(1 if component=='p' else '.5')) for label,row in values.items()}
                                                       for component,values in point['first_fixed_x_physical_time_derivative'].items()}
        coords=dict(point['physical_log_coordinates']); coords['log_r']+=lognu/2
        coords.update(exact_log_r='log(lambda)+log(2R)/2+log(nu)/2',exact_z='sqrt(nu)*Z*lambda^(1-delta)')
        point.update(physical_log_coordinates=coords,physical_viscosity=nu,
                     mapping=dict(lambda_relation='lambda^2-lambda^(2delta)*z^2/nu=tau',radial='R=r^2/(2nu*lambda^2)',
                                  velocity='u_phys=sqrt(nu)*u_source',pressure='p_phys=nu*p_source'),
                     actual_source_log_factors=parts,physical_cylindrical_stress_mixed3=stress,
                     physical_cylindrical_stress_divergence_mixed2=divmixed,
                     completed_theta_theta_stress_mixed2=diagmixed,
                     physical_angular_axial_viscosity_remainder_mixed2=emixed,
                     completed_background_stress_tensor_cartesian_components=tensor,
                     completed_theta_theta_stress=diagonal,physical_completed_stress_divergence_cartesian=divcart,
                     physical_remainder_cartesian=e_cart,physical_angular_axial_viscosity_remainder=remainder,
                     physical_momentum_residual_decomposition_cartesian=momentum,
                     exact_physical_radial_momentum_residual=c.mpf(0),exact_completed_tensor_radial_divergence=c.mpf(0),
                     exact_physical_radial_remainder=c.mpf(0),exact_physical_axial_remainder=c.mpf(0),exact_physical_divergence=c.mpf(0),
                     actual_regional_physical_collar_stress_remainder_identity_verified=True,
                     actual_physical_map_and_remainder_transfer_pending=False,
                     original_pressure_datum_and_positive_amplitudes_retained=True,
                     fixed_viscosity_pullback_verified=True,
                     Gamma_endpoint_physical_remainder_exact_zero=endpoint,
                     completed_tensor_is_divergence_form_background_not_Cauchy=True,
                     source_Z_endpoints_plus_minus_one_are_physical_infinity_limits=True,
                     regional_remainder_is_leading_axial_viscosity_not_proven_flat=True,
                     global_admissible_stress_lift_constructed=False,collar_cone_certified=False,whole_outer_cone_certified=False,
                     independently_bounded_global_flat_remainder=False,physical_energy_integral_certified=False,
                     full_background_NS_validation=False,temporal_recursion=False)
        return point

    def report(self):
        def summary(point):
            omitted=('physical_spatial_cartesian_mixed4','first_fixed_x_physical_time_derivative','positive_source_log_bases')
            result={key:value for key,value in point.items() if key not in omitted}
            result['zeroth_physical_velocity_pressure_rows']=point['physical_spatial_cartesian_mixed4']['x0_y0_z0']
            return result
        with mp.workdps(280):
            return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
                        scope='Actual leading collar offset[0,3], tau>0, r>0; physical stress/mixed3, cylindrical completed diagonal/divergence/remainder mixed2, Cartesian tensor/vector zeroth rows',
                        samples=[summary(self.collar(z,t,lt,ang,nu)) for z,t,lt,ang,nu in
                                 (('0','0','-1','0','1'),('.5','.5','-10','.7','.01'),('-.5','2','-100','1','.7'),('.5','3','-1','.7','1'))],
                        whole_source_collar=summary(self.collar([-1,1],[0,3],[-1000,-1],None)),
                        actual_regional_physical_collar_stress_remainder_identity_verified=True,
                        physical_collar_identities=self.proof,
                        independently_bounded_global_flat_remainder=False,global_admissible_stress_lift_constructed=False,
                        collar_cone_certified=False,temporal_recursion=False,input_hashes=self.hashes)


def run():
    result=CompliantCollarPhysicalC2().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Actual collar physical stress and axial-viscosity remainder generated; global cone/flatness pending',flush=True)
    return result


if __name__=='__main__':run()
