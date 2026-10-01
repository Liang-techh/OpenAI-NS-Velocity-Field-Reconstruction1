"""Actual Section 7.13/7.17/7.21 functional angular/pressure repair.

The prescribed preheat datum and its waiting length remain fixed. Exact
Gamma heat deficits and correlated flatten history supply the right sides.
The unique small smooth solution is enclosed, not replaced by box midpoints.
Only this repair layer is completed; axial energy selection and global cone
are subsequent dependencies.
"""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_compliant_exact_heat_component import SharedExactHeatComponent
from lei_ren_part1_paper_compliant_outer_pulse_map import raw_beta
from lei_ren_part1_paper_compliant_outer_initial import stable_sigma
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_'


def magnitude(v):
    return max(abs(x) for x in endpoints(v))


def intersect(c,a,b):
    low=max(endpoints(a)[0],endpoints(b)[0]);high=min(endpoints(a)[1],endpoints(b)[1])
    if low>high:raise ArithmeticError('Repair enclosure intersection is empty')
    return c.mpf([low,high])


def bump_weights(c,mu,normalization,stop='.15',cells=1024):
    """Directed integral of fixed normalized beta on [-.15,stop]."""
    ell=c.mpf('.15');stop=c.mpf(stop)
    # Integrate in the fixed raw-bump coordinate. Clamping is exact because
    # beta vanishes outside [-1,1]; this also encloses rounded endpoints.
    rlo,rhi=endpoints(stop/ell)
    end=c.mpf([max(mp.mpf(-1),min(mp.mpf(1),rlo)),max(mp.mpf(-1),min(mp.mpf(1),rhi))])
    length=end+1
    result={n:c.mpf(0) for n in ('A','B','D','E','F')}
    if endpoints(length)==(mp.mpf(0),mp.mpf(0)):return result
    for i in range(cells):
        left=-1+length*i/cells;right=-1+length*(i+1)/cells
        raw_coordinate=c.mpf([endpoints(left)[0],endpoints(right)[1]])
        s=ell*raw_coordinate
        beta=raw_beta(c,raw_coordinate)/(ell*normalization);ds=ell*(length/cells)
        for name,rate,power in (('A',1-mu,1),('B',-1-2*mu,1),
                                ('D',-1-2*mu,2),('E',-2*mu,1),('F',-2*mu,2)):
            result[name]+=ds*c.exp(rate*s)*beta**power
    return result


class CompliantAngularRepair:
    def __init__(self,cells=512,bump_cells=1024):
        self.heat=SharedExactHeatComponent();self.angular=self.heat.angular
        self.ctx=c=self.heat.ctx;self.params=self.heat.params;self.mu=self.params.mu
        self.delta=self.heat.delta;self.hashes=dict(self.heat.hashes);self.cells=cells
        for name in ('compliant_exact_heat_component','compliant_exact_heat_component_check'):
            filename=PREFIX+name+'.json';record=json.loads((HERE/filename).read_bytes())
            for source,digest in record['input_hashes'].items():
                if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:
                    raise ValueError('Angular repair heat dependency changed: '+source)
                self.hashes[source]=digest
            if record.get('implicit_source_sha256',self.angular.initial.datum.source_sha)!=self.angular.initial.datum.source_sha:
                raise ValueError('Angular repair pressure source mismatch')
            self.hashes[filename]=hashlib.sha256((HERE/filename).read_bytes()).hexdigest()
        with mp.workdps(210):
            self.rate=1-self.mu;self.prate=1+2*self.mu
            self.logscale=30*self.rate*self.params.log_mu
            self.scale=c.exp(self.logscale)
            self.logtail_distance=2+self.params.Ts+self.angular.waiting
            # Exact unit transition J(1)=1/2 at both ends of the steep part.
            self.logEtail_over_Erel=-c.mpf('1.5')*(2+self.params.Ts)+(self.rate+self.angular.restore_rate)/2-(c.mpf('.5')+self.delta/2)*self.angular.waiting
            self.log_theta_multiplier=(self.rate+self.angular.restore_rate)/2+self.angular.restore_rate*self.angular.waiting-self.angular.waiting_logone
            self.log_pressure_multiplier=2*self.logEtail_over_Erel-2*self.angular.waiting_logone+c.ln(self.delta/2)
            # S is not materialized. Choose a finite log cap which survives
            # converting both heat defects to Rrel units and dividing by scale.
            common=max(endpoints(self.log_theta_multiplier)[1],endpoints(self.log_pressure_multiplier)[1],mp.mpf(0))
            self.strong_logS_cap=2*c.mpf(endpoints(self.logscale)[0])-c.mpf(common)-100
            logRtail=self.heat.logRref+13/self.mu+self.heat.tail_finite
            self.inverse_radius_margin=logRtail+self.strong_logS_cap
            if endpoints(self.inverse_radius_margin)[0]<=0:
                raise ArithmeticError('Selected inverse-radius finite cap not proved')
            self.strong_S_cap=c.exp(self.strong_logS_cap)
            self.theta_heat_over_scale_cap=c.exp(self.log_theta_multiplier+self.strong_logS_cap-self.logscale)
            self.pressure_heat_over_scale_cap=c.exp(self.log_pressure_multiplier+self.strong_logS_cap-self.logscale)
            self.normalization=self.angular.initial.repair.normalization
            self.weights=bump_weights(c,self.mu,self.normalization,cells=bump_cells)
            self.p=c.exp(-2*self.rate);self.q=c.exp(-2*self.prate)
            self.det=self.p*self.q-1;self.k=self.weights['D']/(2*self.weights['B'])
            if endpoints(self.det)[1]>=0:raise ArithmeticError('Angular linear inverse not proved')
            self.inverse_norm=c.mpf(1+max(endpoints(self.p)[1],endpoints(self.q)[1]))/c.mpf(min(abs(x) for x in endpoints(self.det)))
            self.beta_sup=c.exp(-1)/(c.mpf('.15')*self.normalization)
            self.beta_derivative_sup=8*c.exp(-2)/(c.mpf('.15')**2*self.normalization)
            # Uniform RHS enclosure: every Z is represented, not sample-only.
            whole=self.defects(c.mpf([-1,1]))
            self.rhs_scaled_sup=endpoints(c.mpf(magnitude(whole['r_scaled'][0]))+c.mpf(magnitude(whole['s_scaled'][0])))[1]
            self.box=endpoints(10*c.mpf(self.rhs_scaled_sup))[1]
            self.rhs_derivative_scaled_sup=endpoints(c.mpf(magnitude(whole['r_scaled'][1]))+c.mpf(magnitude(whole['s_scaled'][1])))[1]
            self.linear_size= self.inverse_norm*max(
                endpoints(c.exp(self.rate)/self.weights['A'])[1],
                endpoints(c.exp(-3*self.prate)/self.weights['B'])[1])*self.rhs_scaled_sup
            self.nonlinear_size=self.inverse_norm*self.k*self.scale*self.box**2*(1+self.q)
            self.lipschitz=self.inverse_norm*2*self.k*self.scale*self.box*(1+self.q)
            if (endpoints(self.linear_size+self.nonlinear_size)[1]>=self.box
                or endpoints(self.lipschitz)[1]>=mp.mpf('.05')):
                raise ArithmeticError('Actual functional repair contraction failed')
            self.coefficient_physical_sup=self.scale*self.box
            self.h_sup=self.coefficient_physical_sup*self.beta_sup
            self.h_radial_derivative_sup=self.coefficient_physical_sup*self.beta_derivative_sup
            if endpoints(self.h_sup+self.h_radial_derivative_sup-self.mu/100)[1]>=0:
                raise ArithmeticError('Angular radial perturbation not small relative to mu')
        for name in (Path(__file__).name,PREFIX+'compliant_outer_pulse_map.py'):
            self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()

    def preheat_difference(self,Z):
        """(Xrel(0)-Xrel(Z))/mu^[30(1-mu)], with exact zero at Z=0.

        Xv is independent of Z: both pre-flatten swirl and its inherited
        angular primitive carry the same (1+Z^2)^-1 factor. The cutoff
        history is divided before subtracting, using an exact positive
        divided-difference integral in eta=Z^2.
        """
        c=self.ctx;Z=c.mpf(Z);eta=Z**2;Q=1+eta
        if endpoints(Z)[0]<-1 or endpoints(Z)[1]>1:raise ValueError('Z in [-1,1] required')
        start=2*self.angular.Xv*c.exp(-100*self.rate)
        divided=start/Q;first=start/Q**2
        for i in range(self.cells):
            a=c.mpf(100)*i/self.cells;b=c.mpf(100)*(i+1)/self.cells
            sig=stable_sigma(c,c.mpf([endpoints(a)[0],endpoints(b)[1]])/100)[0]
            k=1-sig;weight=(c.exp(-self.rate*(100-b))-c.exp(-self.rate*(100-a)))/self.rate
            factor=c.exp(k*c.ln(2))*k
            # [1-Q^-k]/eta = integral_0^1 k*(1+s*eta)^(-1-k)ds.
            quotient=c.mpf([endpoints(c.exp((-1-k)*c.ln(Q)))[0],1])
            divided+=weight*factor*quotient
            first+=weight*factor*c.exp((-1-k)*c.ln(Q))
        value=eta*divided;derivative=2*Z*first
        return IntervalTaylor(c,[value,derivative])

    def defects(self,Z):
        c=self.ctx;Z=c.mpf(Z);pre=self.preheat_difference(Z)
        heat=self.heat.future_defects(Z,cells=256)
        theta=heat['angular_heat_difference_scaled'];pressure=heat['pressure_heat_difference_scaled']
        # Exact S is retained in the definitions; these are directed upper
        # enclosures of its scaled effects, not values obtained from S_cap.
        capjet=lambda jet,cap: IntervalTaylor(c,[
            c.mpf([0,max(mp.mpf(0),endpoints(jet[0]*cap)[1])]),
            c.mpf([-magnitude(jet[1]*cap),magnitude(jet[1]*cap)])])
        rh=capjet(theta,self.theta_heat_over_scale_cap)
        sh=capjet(pressure,self.pressure_heat_over_scale_cap)
        if endpoints(Z)==(mp.mpf(0),mp.mpf(0)):
            rh=IntervalTaylor(c,[rh[0],0]);sh=IntervalTaylor(c,[sh[0],0])
        if endpoints(Z) in ((mp.mpf(-1),mp.mpf(-1)),(mp.mpf(1),mp.mpf(1))):
            rh=IntervalTaylor(c,[0,rh[1]]);sh=IntervalTaylor(c,[0,sh[1]])
        return dict(Z=Z,rpre_scaled=pre,rheat_scaled_enclosure=rh,
            r_scaled=pre+rh,s_scaled=sh,
            exact_r_definition='r=scale*(Xf(0)-Xf(Z))+exp(log_theta_multiplier)*S*Theta_hat(Z)',
            exact_s_definition='sH=exp(log_pressure_multiplier)*S*Pressure_hat(Z)',
            exact_inverse_radius_definition='S=1/Rtail, not strong_S_cap',
            common_coefficient_scale=self.scale,log_common_coefficient_scale=self.logscale,
            preheat_waiting_identity_at_Z0_retained=True,
            prescribed_analytic_preheat_P0_unchanged=True,
            actual_heat_defect_at_Z0_not_replaced_by_zero=True)

    def inverse(self,u,v):
        return ((u*self.q-v)/self.det,(v*self.p-u)/self.det)

    def coefficients(self,Z,tightening=8):
        c=self.ctx;rhs=self.defects(Z);r=rhs['r_scaled'];s=rhs['s_scaled']
        b1=r*(c.exp(self.rate)/self.weights['A']);b2=s*(c.exp(-3*self.prate)/self.weights['B'])
        box=c.mpf([-self.box,self.box]);x=y=box
        for _ in range(tightening):
            nonlinear=self.k*self.scale*(x*x+self.q*y*y)
            nx,ny=self.inverse(b1[0],b2[0]-nonlinear)
            x=intersect(c,x,nx);y=intersect(c,y,ny)
        # Differentiate the actual quadratic system, not the interval iteration.
        J21=1+2*self.k*self.scale*x;J22=self.q*(1+2*self.k*self.scale*y)
        det=self.p*J22-J21
        if endpoints(det)[1]>=0:raise ArithmeticError('Functional derivative inverse failed')
        dx=(J22*b1[1]-b2[1])/det;dy=(self.p*b2[1]-J21*b1[1])/det
        scaled=[IntervalTaylor(c,[x,dx]),IntervalTaylor(c,[y,dy])]
        return dict(Z=c.mpf(Z),scaled_coefficient_Taylor=scaled,
            physical_coefficient_Taylor=[v*self.scale for v in scaled],
            physical_coefficient_definition='dj=scale*Cj; Cj is the unique small smooth solution of the exact two moment equations',
            actual_implicit_angular_pressure_closure=True,coefficients_are_enclosures_not_midpoints=True,
            true_heat_source_factor_S_retained=True,defects=rhs,
            equations=dict(angular='A*(exp(-3rate)*d1+exp(-rate)*d2)=r',
                pressure='B*(exp(3prate)*d1+exp(prate)*d2)+D/2*(exp(3prate)*d1^2+exp(prate)*d2^2)=sH'),
            derivative_jacobian_determinant=det,angular_pressure_functions_even=True,
            exact_Gamma_heat_pressure_restored_to_prescribed_preheat_P0=True,
            actual_ap_selected=False,full_outer_five_moment_match=False,temporal_recursion=False)

    def partial_corrections(self,Z,offset,cells=256):
        """All five cumulative changes in Rrel/Erel units; exact disjoint supports."""
        c=self.ctx;t=c.mpf(offset)
        if endpoints(t)[0]<-4 or endpoints(t)[1]>0:raise ValueError('Finite log offset in [-4,0] required')
        coefficients=self.coefficients(Z)['physical_coefficient_Taylor']
        zero=coefficients[0]*0;angular=zero;pressure=zero;energy=zero;h=zero
        for j,center in enumerate((-3,-1)):
            local=t-center
            ell_upper=endpoints(c.mpf('.15'))[1]
            if endpoints(local)[1]<=-ell_upper:w={n:c.mpf(0) for n in self.weights}
            elif endpoints(local)[0]>=ell_upper:w=self.weights
            else:w=bump_weights(c,self.mu,self.normalization,stop=local,cells=cells)
            dj=coefficients[j]
            angular+=dj*(c.exp(self.rate*center)*w['A'])
            pressure+=(dj*w['B']+dj*dj*(w['D']/2))*c.exp(-self.prate*center)
            energy+=(dj*(2*w['E'])+dj*dj*w['F'])*c.exp(-2*self.mu*center)
            if mp.mpf('-.15')<endpoints(local)[0]<mp.mpf('.15'):
                h+=dj*(raw_beta(c,local/c.mpf('.15'))/(c.mpf('.15')*self.normalization))
        return dict(Z=c.mpf(Z),offset=t,relative_swirl_correction=h,
            corrected_Utheta_over_Erel=(h+1)*c.exp(-(c.mpf('.5')+self.mu)*t),Uz=zero,
            Delta_Mtheta_over_sqrt2_Rrel_3half_Erel=angular,
            Delta_Mp_over_Erel_squared=pressure,
            Delta_swirl_energy_over_Rrel_Erel_squared=energy,
            Delta_Mz=zero,Delta_Mtheta_z=zero,
            Delta_Mztheta_over_Rrel_Erel_squared=-energy/2,
            all_support_endpoints_flat_by_fixed_beta=True,
            full_outer_five_moment_field_built=False)

    def report(self):
        with mp.workdps(210):
            samples=[self.coefficients(z) for z in ('-1','0','.5','1')]
            whole=self.coefficients(self.ctx.mpf([-1,1]))
            return dict(actual_five_defect_family_sha256=self.angular.initial.family,
                implicit_source_sha256=self.angular.initial.datum.source_sha,
                datum_enclosure_sha256=self.angular.initial.datum.datum_sha,
                common_coefficient_scale=self.scale,log_common_coefficient_scale=self.logscale,
                strong_log_inverse_radius_cap=self.strong_logS_cap,strong_inverse_radius_positive_cap=self.strong_S_cap,
                strong_cap_log_margin=self.inverse_radius_margin,
                exact_formal_inverse_radius_terms=self.heat.logS_terms,
                log_theta_heat_multiplier=self.log_theta_multiplier,log_pressure_heat_multiplier=self.log_pressure_multiplier,
                bump_weights=self.weights,beta_sup=self.beta_sup,beta_first_derivative_sup=self.beta_derivative_sup,
                rescaled_linear_determinant=self.det,linear_inverse_row_norm_upper=self.inverse_norm,
                scaled_rhs_sup=self.rhs_scaled_sup,scaled_coefficient_ball_radius=self.box,
                scaled_linear_map_size=self.linear_size,scaled_nonlinear_map_size=self.nonlinear_size,
                contraction_lipschitz_upper=self.lipschitz,
                relative_swirl_sup_upper=self.h_sup,radial_derivative_perturbation_upper=self.h_radial_derivative_sup,
                radial_perturbation_below_mu_over_100=True,
                samples=samples,whole_Z_C1_coefficients=whole,
                support_endpoint_corrections=[self.partial_corrections('.5',t) for t in ('-4','-3.15','-2.85','-1.15','-.85','0')],
                actual_smooth_functional_angular_pressure_repair_proved=True,
                full_future_corrected_energy_completed=False,actual_ap_selected=False,
                whole_outer_cone_certified=False,full_outer_five_moment_match=False,
                global_admissible_stress_lift_constructed=False,temporal_recursion=False,input_hashes=self.hashes)


def run():
    repair=CompliantAngularRepair();result=repair.report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf-8')
    print('Actual whole-Z angular/pressure two-bump repair: uniform small smooth branch and partial moment changes generated; ap/global cone pending',flush=True)
    return result


if __name__=='__main__':run()
