"""Same-source rooted velocity, pressure and nonsingular physical vorticity.

H(a)=0 is retained in fresh coupled radial seeds. Infinite model tails vanish
at the shared root to the requested axial orders; admitted nonlinear tails
remain. Transport uses full analytic norms, not a fitted or midpoint field.
Physical results are signed sums of separately retained logarithmic factors.
"""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
from lei_ren_part1_paper_compliant_root_centered_peak import CompliantRootCenteredPeak
from lei_ren_part1_paper_compliant_core_physical_field import symmetric,intersection
from lei_ren_part1_paper_compliant_core_uniform_bounds import embedding
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_functional_core_step import initial_rows
from lei_ren_part1_paper_logarithmic_core_step import advance_scaled_one
from lei_ren_part1_paper_analytic_radial_tail import tail_factor
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_candidate_exact_amplitude import _encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'


class CompliantRootedCoreField:
    def __init__(self):
        self.peak=CompliantRootCenteredPeak();self.core=self.peak.core;self.ctx=c=self.core.ctx
        self.rebuild=self.peak.amplitude.rebuild;self.epsilon=self.peak.epsilon
        self.hashes=dict(self.peak.hashes);name=PREFIX+'root_centered_peak_check.json'
        receipt=json.loads((HERE/name).read_bytes())
        if not receipt['all_passed'] or receipt['implicit_source_sha256']!=self.core.source:
            raise ValueError('Accepted current-source root peak required')
        for path,digest in receipt['input_hashes'].items():
            if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:raise ValueError('Root peak changed: '+path)
        self.hashes.update(receipt['input_hashes'])
        for path in (name,Path(__file__).name,'lei_ren_part1_paper_functional_core_step.py',
                     'lei_ren_part1_paper_logarithmic_core_step.py','lei_ren_part1_paper_analytic_radial_tail.py'):
            self.hashes[path]=hashlib.sha256((HERE/path).read_bytes()).hexdigest()
        self.root_rows=self.build_root_rows();self.cache={}

    def build_root_rows(self,degree=24,depth=3):
        c=self.ctx;order=degree+depth;a=self.peak.anchor
        t=IntervalTaylor.variable(c,0,order);one=IntervalTaylor.constant(c,1,order);z=one*a+t
        # Exact H(a)=0 BEFORE any enclosure of the correlated anchor.
        H=t*self.peak.hp+t*t*self.peak.q1-4*t**3
        L=one-z*z*self.core.delta;u=4*z+self.core.j
        gradient=L*H/(H*H+self.core.sigma**2)
        S=[]
        for k in range(order+1):
            bound=self.rebuild.S_bound/(self.rebuild.eta/2)**k
            S.append(c.mpf([0,endpoints(bound)[1]]) if k==0 else symmetric(c,bound))
        pressure=self.core.datum.normalized_jets(a,order)
        scale=c.exp(2*self.core.logP-self.core.logLambda)
        fixed=dict(ell_Z_taylor=[-value for value in gradient.coefficients],S_Z_taylor=S,
            U0_Z_taylor=list(u.coefficients),P0_Z_taylor=[scale*v for v in pressure['normalized_pressure_coefficients']])
        rows=initial_rows(c,fixed,a,degree,required_depth=depth)
        for n in range(degree):advance_scaled_one(c,fixed,rows,n,a,self.core.delta,self.epsilon)
        if endpoints(fixed['ell_Z_taylor'][0])!=(mp.mpf(0),mp.mpf(0)):
            raise ArithmeticError('Exact root gradient identity lost')
        return dict(Z=a,radial_degree=degree,axial_depth=depth,rows=rows,fixed=fixed,
            H_Z_taylor=list(H.coefficients),same_pressure_Cstar_and_positive_swirl_source=True,
            exact_shared_H_root_used_before_enclosure=True,old_finite_rows_read=False,
            local_radial_generation_not_temporal_recursion=True)

    def rooted_profile(self,xi,rho,i=0,k=0):
        c=self.ctx;r=c.mpf(rho);s=c.mpf(xi);N=self.root_rows['radial_degree']
        if min(i,k)<0 or i>2 or k>2 or i+k>3:raise ValueError('Local mixed derivatives total<=3 and i,k<=2 required')
        if endpoints(r)[0]<0 or endpoints(r)[1]>endpoints(c.mpf('4.1'))[1] or endpoints(abs(s))[1]>4:
            raise ValueError('Original rho in[0,4.1], rooted xi in[-4,4] required')
        key=(s._mpi_,r._mpi_,i,k)
        if key in self.cache:return self.cache[key]
        rows=self.root_rows['rows'];distance=abs(self.peak.b*s);rmax=c.mpf(endpoints(r)[1])
        factor=tail_factor(c,degree=N,radial_order=i,axial_order=k,radius=rmax,h=self.core.h)['tail_per_Xh_norm'] if endpoints(rmax)[1]>0 else c.mpf(0)
        # chi(a+t) has order>=2 in t. Thus chi^n cannot contribute to
        # axial order k when 2n>k: the Bessel radial tail here is EXACTLY 0.
        if 2*(N+1)<=k:raise ArithmeticError('Root Bessel tail valuation not applicable')
        tail=self.core.correction*factor
        values={};details={}
        for label,norm in (('Phi',self.rebuild.phi_norm),('Psi',self.rebuild.psi_norm),('MeanPsi',self.rebuild.psi_norm)):
            polynomial=c.mpf(0)
            for n in range(i,N+1):
                if label!='Phi' and n==0:continue
                coefficient=rows['A' if label=='Phi' else 'Uz'][n][k]
                if label!='Phi':coefficient=coefficient/self.epsilon
                if label=='MeanPsi':coefficient=coefficient/(n+1)
                polynomial+=coefficient*(math.factorial(n)//math.factorial(n-i))*math.factorial(k)*r**(n-i)
            transport=norm*embedding(c,self.core.h,rmax,i,k+1)*distance
            value=polynomial+symmetric(c,tail+transport)
            if endpoints(r)==(mp.mpf(0),mp.mpf(0)) and i==0:
                value=c.mpf(1 if label=='Phi' and k==0 else 0)
            values[label]=value;details[label]=dict(root_finite_polynomial=polynomial,
                root_nonlinear_radial_tail_upper=tail,full_analytic_axial_transport_upper=transport)
        affine=c.mpf(0)
        if i==0:
            if k==0:affine=4*self.peak.anchor+self.core.j+4*self.peak.b*s
            elif k==1:affine=c.mpf(4)
        values['Uz']=affine+self.epsilon*values['Psi'];values['Mz_over_R']=affine+self.epsilon*values['MeanPsi']
        result=dict(xi=s,rho=r,radial_order=i,axial_order=k,source_jets=values,error_details=details,
            root_model_radial_tail_exact_zero=True,normalized_Psi_built_before_affine_recombination=True,
            Mz_over_R_is_same_Uz_radial_average=True,full_analytic_source_transport=True)
        self.cache[key]=result;return result

    def pressure_integral(self,xi,rho):
        c=self.ctx;r=c.mpf(rho);s=c.mpf(xi)
        if endpoints(r)[1]==0:return dict(V=c.mpf(0),finite_convolution=c.mpf(0),error_upper=c.mpf(0))
        N=self.root_rows['radial_degree'];A=[self.root_rows['rows']['A'][n][0] for n in range(N+1)]
        finite=sum((A[i]*A[j]*r**(i+j+1)/(i+j+1) for i in range(N+1) for j in range(N+1)),c.mpf(0))
        cover=self.rooted_profile(s,c.mpf([0,endpoints(r)[1]]))
        bounds=cover['error_details']['Phi'];E=bounds['root_nonlinear_radial_tail_upper']+bounds['full_analytic_axial_transport_upper']
        Pmax=sum((abs(coefficient)*c.mpf(endpoints(r)[1])**n for n,coefficient in enumerate(A)),c.mpf(0))
        error=r*(2*Pmax+E)*E
        lower=c.mpf(endpoints(r)[0])*c.mpf(endpoints(self.core.phi_floor)[0])**2
        upper=c.mpf(endpoints(r)[1])*c.mpf(endpoints(self.core.phi_ceiling)[1])**2
        V=intersection(c,finite+symmetric(c,error),c.mpf([endpoints(lower)[0],endpoints(upper)[1]]))
        return dict(V=V,finite_convolution=finite,error_upper=error,
            original_primitive='V=int_0^rho Phi(s,Z)^2 ds; V_rho=Phi^2, V(0,Z)=0',
            same_nonlinear_Phi_and_positive_primitive=True)

    @staticmethod
    def term(coefficient,*logs):
        return dict(coefficient_enclosure=coefficient,positive_scale_log_terms=list(logs),
            definition='coefficient * product(exp(log) for log in separate scale terms)',source_enclosure_not_point_selection=True)

    def field(self,xi,rho,log_tau='-1',theta='0'):
        c=self.ctx;s=c.mpf(xi);r=c.mpf(rho);t=c.mpf(log_tau);angle=c.mpf(theta)
        if not all(mp.isfinite(v) for value in (s,r,t,angle) for v in endpoints(value)):
            raise ValueError('Finite source coordinates/log-time/angle required')
        peak=self.peak.evaluate(s);axis=endpoints(r)==(mp.mpf(0),mp.mpf(0))
        packets={};jets={}
        for i,k in ((0,0),(1,0),(0,1),(0,2),(1,1)):
            packet=self.rooted_profile(s,r,i,k);packets[str((i,k))]=packet;jets[(i,k)]=packet['source_jets']
        z=self.peak.anchor+s*self.peak.b;D=1-z*z;L=1-self.core.delta*z*z
        ell=(t-c.ln(D))/2;delta=self.core.delta;A=peak['normalized_F0_relative_to_shared_anchor']
        P,U,M=jets[(0,0)]['Phi'],jets[(0,0)]['Uz'],jets[(0,0)]['Mz_over_R']
        Pr,Ur,Mr=jets[(1,0)]['Phi'],jets[(1,0)]['Uz'],jets[(1,0)]['Mz_over_R']
        Pz,Uz,Mz=jets[(0,1)]['Phi'],jets[(0,1)]['Uz'],jets[(0,1)]['Mz_over_R']
        Mzz=jets[(0,2)]['Mz_over_R'];Mrz=jets[(1,1)]['Mz_over_R']
        numerator=2*z*U-(1-delta)*z*M-D*Mz
        Q=numerator/L;Qr=(2*z*Ur-(1-delta)*z*Mr-D*Mrz)/L
        Qz=(2*U+2*z*Uz-(1-delta)*M+(1+delta)*z*Mz-D*Mzz)/L+2*delta*z*numerator/L**2
        zero=self.term(c.mpf(0));sqrt_half=c.sqrt(r/2);sqrt_two=c.sqrt(2*r)
        radial=[zero] if axis else [self.term(sqrt_half*Q,-self.core.logLambda/2,-ell)]
        swirl=[zero] if axis else [self.term(sqrt_two*A*P,-self.core.logC,-self.core.logLambda/2,-ell,-delta*ell)]
        axial=[self.term(U,-ell,-delta*ell)]
        primitive=self.pressure_integral(s,r);datum=self.core.datum.normalized_jets(z,0)['normalized_pressure_coefficients'][0]
        pressure=[self.term(datum,2*self.core.logP,-2*ell,-2*delta*ell),
                  self.term(A*A*primitive['V'],-self.core.logLambda,-2*self.core.logC,-2*ell,-2*delta*ell)]
        if axis:wr=wt=[zero]
        else:
            # sqrt(epsilon)/b=1/sigma is canceled analytically. Never form K'/b.
            regular=(2+delta)*z*P-D*Pz+2*z*r*Pr
            wr=[self.term(sqrt_two*A*regular/L,-self.core.logC,-self.core.logLambda/2,-2*ell),
                self.term(sqrt_two*A*D*peak['K_derivatives'][0]*P/(self.core.sigma*L),-self.core.logC,-2*ell)]
            wt=[self.term(sqrt_half*(D*Qz-2*z*(Q+r*Qr))/L,-self.core.logLambda/2,-2*ell,delta*ell),
                self.term(-2*sqrt_half*jets[(1,0)]['Psi'],-self.core.logLambda/2,-2*ell,-delta*ell)]
        wz=[self.term(2*A*(P+r*Pr),-self.core.logC,-2*ell,-delta*ell)]
        co=c.cos(angle);si=c.sin(angle)
        def rotate(radial,angular):
            scale=lambda terms,m:[self.term(m*term['coefficient_enclosure'],*term['positive_scale_log_terms']) for term in terms]
            return dict(x=scale(radial,co)+scale(angular,-si),y=scale(radial,si)+scale(angular,co))
        velocity=rotate(radial,swirl);velocity['z']=axial
        vorticity=rotate(wr,wt);vorticity['z']=wz
        div=Q+r*Qr+(D*Uz-(1+delta)*z*U-2*z*r*Ur)/L
        if not endpoints(div)[0]<=0<=endpoints(div)[1]:raise ArithmeticError('Moment-bound divergence diagnostic excludes zero')
        return dict(xi=s,rho=r,log_tau=t,theta=angle,axis=axis,log_lambda=ell,
            root_and_offset=dict(shared_root=self.peak.anchor,separate_offset_width=self.peak.b,xi=s),
            physical_coordinates=dict(radial=self.term(sqrt_two,-self.core.logLambda/2,ell),
                x=self.term(co*sqrt_two,-self.core.logLambda/2,ell),y=self.term(si*sqrt_two,-self.core.logLambda/2,ell),
                z_terms=[self.term(self.peak.anchor,ell,-delta*ell),self.term(self.peak.b*s,ell,-delta*ell)]),
            cylindrical_velocity=dict(ur=radial,utheta=swirl,uz=axial),pressure_terms=pressure,
            cartesian_velocity=velocity,cylindrical_vorticity=dict(omega_r=wr,omega_theta=wt,omega_z=wz),
            cartesian_vorticity=vorticity,radial_recovery=dict(Q=Q,Q_rho=Qr,Q_Z=Qz),
            normalized_divergence_interval=div,divergence_scale_log=-2*ell,
            structural_divergence_exact_by_same_moment_primitive=True,
            pressure_radial_primitive=primitive,source_mixed_jet_packets=packets,
            vorticity_amplitude_inverse_width_canceled_before_enclosure=True,
            full_local_rooted_velocity_pressure_vorticity_enclosures=True,
            exact_source_axis_regular=True,absolute_swirl_not_materialized=True,
            full_point_physical_field_evaluation=False,measured_blowup_dynamics=False,temporal_recursion=False)

    def report(self):
        points={}
        for xi in ('-1','0','1'):
            for time in ('-1','-10'):
                points[xi+':'+time]=[self.field(xi,rho,time,'.7') for rho in ('0','2','4.1')]
        return dict(actual_five_defect_family_sha256=self.core.family,implicit_source_sha256=self.core.source,
            datum_enclosure_sha256=self.core.datum.datum_sha,rooted_fresh_core_packet=self.root_rows,
            physical_field_packets=points,same_source_rooted_vector_and_pressure_available=True,
            same_source_rooted_vorticity_available=True,structural_divergence_exact_by_same_moment_primitive=True,
            positive_factored_original_pressure_increment=True,uncapped_original_scales_preserved=True,
            full_point_physical_field_evaluation=False,all_annular_source_values_resolved=False,
            measured_blowup_dynamics=False,whole_vortex_aspect_ratio_measured=False,
            physical_energy_integral_certified=False,admissible_stress_lift_constructed=False,temporal_recursion=False,
            input_hashes=self.hashes)


def run():
    with mp.workdps(400):result=CompliantRootedCoreField().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(_encode(result),indent=2)+'\n',encoding='utf8')
    print('Same-source rooted velocity, pressure, Cartesian vector and physical vorticity generated',flush=True)
    return result


if __name__=='__main__':run()
