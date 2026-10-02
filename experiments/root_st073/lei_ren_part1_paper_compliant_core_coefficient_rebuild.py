"""Fresh radial coefficient enclosures of the admitted compliant core.

Seeds use the SAME selected Cstar and analytic preheat datum. No old finite
coefficient rows are read. The coupled radial equations retain all source
terms. Analytic Xh norms bound infinite radial tails. Outputs are directed
source enclosures, not midpoint solutions or temporal coefficient recursion.
"""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
from lei_ren_part1_paper_compliant_core_physical_field import (
    CompliantCorePhysicalField,symmetric,model_phi_jets,gridkey)
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_functional_core_step import initial_rows
from lei_ren_part1_paper_logarithmic_core_step import advance_scaled_one
from lei_ren_part1_paper_analytic_radial_tail import tail_factor
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_compliant_five_moment_repair import pack

HERE=Path(__file__).parent


class CompliantCoreCoefficientRebuild:
    def __init__(self):
        self.core=CompliantCorePhysicalField();self.ctx=c=self.core.ctx
        self.hashes=dict(self.core.hashes)
        for name in (Path(__file__).name,'lei_ren_part1_paper_logarithmic_core_step.py',
                     'lei_ren_part1_paper_functional_core_step.py',
                     'lei_ren_part1_paper_interval_taylor.py',
                     'lei_ren_part1_paper_analytic_radial_tail.py'):
            self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        major=self.core.records['core_transfer'];tube=self.core.records['shared_analytic_tube']
        self.eta=read_interval(c,tube['complex_tube_radius'])
        self.phi_norm=read_interval(c,major['Phi_ball_norm_upper'])
        self.psi_norm=read_interval(c,major['Psi_ball_norm_upper'])
        required=self.core.Lambda*self.core.Gbar+2*self.core.logLambda+1000
        if endpoints(self.core.logC-required)[0]<=0:
            raise ValueError('Selected Cstar must satisfy the admitted complex amplitude guard')
        # For this selected Cstar, |epsilon^2 F0^2| on the SAME complex
        # tube is <= exp(-6 logLambda-2000). This representable positive
        # bound is not a selected value of the implicit F0 source.
        self.S_bound=c.exp(-6*self.core.logLambda-2000)
        if endpoints(self.S_bound)[1]<=0:raise ArithmeticError('Positive swirl source bound lost')
        self.model_tail_cache={}

    def seed(self,Z,degree,depth):
        c=self.ctx;z0=c.mpf(Z);order=degree+depth
        if endpoints(z0)[0]<-1 or endpoints(z0)[1]>1:raise ValueError('Core axial center outside [-1,1]')
        z=IntervalTaylor.variable(c,z0,order);one=IntervalTaylor.constant(c,1,order)
        u=4*z+self.core.j;L=one-z*z*self.core.delta
        H=z*((1-self.core.delta)/2)+(one-z*z)*u
        denominator=H*H+self.core.sigma**2
        coefficients=list(denominator.coefficients)
        coefficients[0]=H[0]**2+self.core.sigma**2
        gradient=L*H/IntervalTaylor(c,coefficients)
        S=[]
        for k in range(order+1):
            bound=self.S_bound/(self.eta/2)**k
            S.append(c.mpf([0,endpoints(bound)[1]]) if k==0 else symmetric(c,bound))
        pressure=self.core.datum.normalized_jets(z0,order)
        scale=c.exp(2*self.core.logP-self.core.logLambda)
        fixed=dict(ell_Z_taylor=[-v for v in gradient.coefficients],S_Z_taylor=S,
            U0_Z_taylor=list(u.coefficients),
            P0_Z_taylor=[scale*v for v in pressure['normalized_pressure_coefficients']])
        return fixed,initial_rows(c,fixed,z0,degree,required_depth=depth)

    def rebuild(self,Z,degree=24,depth=5):
        if not isinstance(degree,int) or degree<4 or not isinstance(depth,int) or depth<1:
            raise ValueError('Integer radial degree>=4 and axial depth>=1 required')
        c=self.ctx;fixed,rows=self.seed(Z,degree,depth)
        for n in range(degree):
            advance_scaled_one(c,fixed,rows,n,c.mpf(Z),self.core.delta,self.core.epsilon)
        return dict(Z=c.mpf(Z),radial_degree=degree,axial_depth=depth,rows=rows,fixed=fixed,
            row_units=dict(A='Phi radial coefficient in rho=Lambda*R',
                           Uz='physical Uz radial coefficient in rho',
                           P='epsilon*physical P radial coefficient in rho'),
            implicit_F0_source='F0=exp(-selected_logCstar-Lambda*G), G anchored at original H root',
            scaled_swirl_source='S=epsilon^2*F0^2; positive implicit source, nonzero Cauchy enclosure',
            S_complex_upper=self.S_bound,S_cauchy_radius=self.eta/2,
            actual_selected_logCstar=self.core.logC,
            fresh_compliant_axis_pressure_seed=True,old_finite_coefficient_rows_read=False,
            fresh_coupled_radial_coefficient_enclosures_recomputed=True,
            point_parameter_representatives_selected=False,temporal_recursion=False)

    def profile(self,packet,rho,radial_order=0,axial_order=0):
        c=self.ctx;r=c.mpf(rho);i=radial_order;k=axial_order;N=packet['radial_degree']
        if not isinstance(i,int) or not isinstance(k,int) or min(i,k)<0 or i>N or k>packet['axial_depth']:
            raise ValueError('Requested derivative outside finite coefficient depth')
        if endpoints(r)[0]<0 or endpoints(r)[1]>endpoints(c.mpf('4.1'))[1]:
            raise ValueError('Original scaled core radius rho in [0,4.1] required')
        fields={};tails={};polynomials={}
        for label,norm in (('Phi',self.phi_norm),('Uz',self.core.epsilon*self.psi_norm),
                           ('Mz_over_R',self.core.epsilon*self.psi_norm)):
            rows=packet['rows']['A' if label=='Phi' else 'Uz']
            value=c.mpf(0)
            for n in range(i,N+1):
                coefficient=rows[n][k]/(n+1) if label=='Mz_over_R' else rows[n][k]
                value+=coefficient*(math.factorial(n)//math.factorial(n-i))*math.factorial(k)*r**(n-i)
            # The affine U0 is entirely present in row zero. Only epsilon*Psi
            # contributes to the Uz/Mz radial tail. Averaging divides each
            # coefficient by n+1, so the same tail bound remains conservative.
            if endpoints(r)[1]==0:tail=c.mpf(0)
            else:
                factor=tail_factor(c,degree=N,radial_order=i,axial_order=k,
                    radius=c.mpf(endpoints(r)[1]),h=self.core.h)['tail_per_Xh_norm']
                if i+k<=5:
                    # The explicit Phi Bessel model has factorial decay.
                    # Its nonlinear correction has the admitted Xh tail.
                    # Psi's model is linear in rho and thus has no tail here.
                    if label=='Phi':
                        key=(packet['Z']._mpi_,r._mpi_,N)
                        if key not in self.model_tail_cache:
                            chi=self.core.axis_inputs(packet['Z'])['chi']
                            self.model_tail_cache[key]=model_phi_jets(c,chi,r,degree=N)[1]
                        tail=self.model_tail_cache[key][gridkey(i,k)]+self.core.correction*factor
                    else:tail=self.core.epsilon*self.core.correction*factor
                else:tail=norm*factor
            polynomials[label]=value;tails[label]=tail;fields[label]=value+symmetric(c,tail)
        return dict(rho=r,radial_order=i,axial_order=k,finite_polynomials=polynomials,
            infinite_radial_tail_bounds=tails,source_profile_enclosures=fields,
            tail_decomposition='Explicit Bessel model plus admitted nonlinear correction for total order<=5; full analytic norms at higher orders',
            analytic_tail_norms_admitted=True,finite_polynomial_is_not_complete_solution=True)

    def values(self,packet,rho):
        c=self.ctx;base=self.profile(packet,rho);axial=self.profile(packet,rho,axial_order=1)
        z=packet['Z'];L=1-self.core.delta*z*z;d=1-z*z
        values=base['source_profile_enclosures'];Q=(2*z*values['Uz']-(1-self.core.delta)*z*values['Mz_over_R']
            -d*axial['source_profile_enclosures']['Mz_over_R'])/L
        r=c.mpf(rho);N=packet['radial_degree']
        pfinite=sum((packet['rows']['P'][n][0]*r**n for n in range(N+1)),c.mpf(0))
        # From |Phi_n|<=B/(20^n(n+1)^2), the convolution coefficient is
        # <=B^2*(n+1)/20^n. P_scaled[n+1]=S*(Phi^2)_n/(n+1).
        ptail=self.S_bound*self.phi_norm**2*r**(N+1)/(20**N*(1-r/20))
        return dict(rho=r,Z=z,Phi=values['Phi'],Uz=values['Uz'],Mz_over_R=values['Mz_over_R'],
            radial_recovery_Q=Q,P_scaled=pfinite+symmetric(c,ptail),
            P_scaled_finite_polynomial=pfinite,P_scaled_infinite_tail_bound=ptail,
            profile_derivatives=[base,axial],
            physical_fields='Ur=sqrt(R/2)*Q; F=F0*Phi; Utheta=sqrt(2R)*F; P=P_scaled/epsilon',
            actual_F0_not_materialized=True,full_point_physical_field_evaluation=False,
            temporal_recursion=False)

    def report(self):
        packets={}
        for Z in ('.3','.5','-.5','0'):
            packet=self.rebuild(Z)
            packet['value_packets']=[self.values(packet,r) for r in ('0','2','4','4.1')]
            packets[Z]=packet
            print('Fresh compliant core radial rows: Z='+Z+', degree=24',flush=True)
        return dict(actual_five_defect_family_sha256=self.core.family,implicit_source_sha256=self.core.source,
            datum_enclosure_sha256=self.core.datum.datum_sha,recomputed_center_packets=packets,
            fresh_coupled_radial_coefficients_recomputed=True,infinite_radial_tails_bound=True,
            same_selected_Cstar_and_pressure_source=True,old_finite_coefficient_rows_read=False,
            full_point_physical_field_evaluation=False,measured_blowup_dynamics=False,
            admissible_stress_lift_constructed=False,temporal_recursion=False,input_hashes=self.hashes)


def run():
    with mp.workdps(280):result=CompliantCoreCoefficientRebuild().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()
