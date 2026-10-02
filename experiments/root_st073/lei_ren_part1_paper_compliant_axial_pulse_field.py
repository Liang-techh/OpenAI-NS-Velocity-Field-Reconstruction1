"""Selected axial pulse and five partial primitives in segmented coordinates.

Main coordinate xi=mu*log(R/Rp); end coordinate s=log(R/Rv) retains the
two narrow supports. Super-small factors are formal exact sources with
directed positive caps. Terminal zeros follow from selected moment equations
and empty backward supports, not deletion of incoming moment histories.
"""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_compliant_axial_amplitude_selection import CompliantAxialAmplitude
from lei_ren_part1_paper_compliant_outer_pulse_map import raw_beta
from lei_ren_part1_paper_compliant_outer_angular_repair import magnitude
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_'


def sigma_enclosure(c,x):
    """Original flat sigma, without exponentiating exp(-1/mu^2)."""
    def point(v):
        if v<=0:return c.mpf(0)
        if v>=1:return c.mpf(1)
        y=c.mpf(v);odds=-1/y**2+1/(1-y)**2
        if endpoints(odds)[1]<-1000:return c.mpf([0,endpoints(c.exp(-1000))[1]])
        if endpoints(odds)[0]>1000:return c.mpf([endpoints(1-c.exp(-1000))[0],1])
        e=c.exp(odds);return e/(1+e)
    lo,hi=endpoints(c.mpf(x));a=point(lo);b=point(hi)
    return c.mpf([max(mp.mpf(0),endpoints(a)[0]),min(mp.mpf(1),endpoints(b)[1])])


def gp(c,xi,cells=128):
    """True gp and its derivative: primitive(sigma(50v))*sigma(11-xi)."""
    xi=c.mpf(xi);lo,hi=endpoints(xi)
    if hi<=0 or lo>=11:return dict(value=c.mpf(0),derivative=c.mpf(0))
    def primitive(x):
        if x<=0:return c.mpf(0)
        y=c.mpf(x)
        if x>=mp.mpf('.02'):return y-c.mpf('.01')
        total=c.mpf(0)
        for i in range(cells):
            a=y*i/cells;b=y*(i+1)/cells
            total+=(b-a)*sigma_enclosure(c,50*c.mpf([endpoints(a)[0],endpoints(b)[1]]))
        return total
    P=c.mpf([max(mp.mpf(0),endpoints(primitive(lo))[0]),max(mp.mpf(0),endpoints(primitive(hi))[1])])
    cut=sigma_enclosure(c,11-xi)
    # 0<=sigma'<=8; it is exactly zero off its open transition.
    slope=c.mpf(0) if hi<=10 or lo>=11 else c.mpf([0,8])
    derivative=sigma_enclosure(c,50*xi)*cut-P*slope
    return dict(value=P*cut,derivative=derivative,exact_shape_retained=True)


def gp_energy(c,xi,K,cells=1024):
    xi=c.mpf(xi);lo,hi=endpoints(xi)
    if lo==hi==0:return c.mpf(0)
    if lo>=11:return K
    result=c.mpf(0)
    for i in range(cells):
        a=xi*i/cells;b=xi*(i+1)/cells
        shape=gp(c,c.mpf([endpoints(a)[0],endpoints(b)[1]]))['value']
        result+=c.exp(-2*a)*decay_integral(c,2,xi/cells)*shape**2
    return result


def decay_integral(c,rate,length):
    """Integral_0^length exp(-rate*v)dv, including sub-precision lengths."""
    rate=c.mpf(rate);length=c.mpf(length)
    if endpoints(length)[0]<0:raise ValueError('Positive integration length required')
    if endpoints(rate*length)[1]<mp.mpf('1e-20'):
        return c.mpf([endpoints(length*c.exp(-rate*length))[0],endpoints(length)[1]])
    return (1-c.exp(-rate*length))/rate


def gp_future_energy(c,xi,cells=512):
    """Positive remaining main energy; direct future integral near xi=11."""
    xi=c.mpf(xi);length=11-xi
    if endpoints(length)==(mp.mpf(0),mp.mpf(0)):return c.mpf(0)
    total=c.mpf(0)
    for i in range(cells):
        a=xi+length*i/cells;b=xi+length*(i+1)/cells
        shape=gp(c,c.mpf([endpoints(a)[0],endpoints(b)[1]]))['value']
        total+=c.exp(-2*a)*decay_integral(c,2,length/cells)*shape**2
    return total


def backward_bump_weights(c,mu,normalization,start,cells=256):
    """Positive future beta integrals in fixed raw coordinates; no 1-minus-CDF."""
    ell=c.mpf('.15');start=c.mpf(start)
    rlo,rhi=endpoints(start/ell)
    begin=c.mpf([max(mp.mpf(-1),min(mp.mpf(1),rlo)),max(mp.mpf(-1),min(mp.mpf(1),rhi))])
    length=1-begin;out=[c.mpf(0),c.mpf(0),c.mpf(0)]
    if endpoints(length)==(mp.mpf(0),mp.mpf(0)):return out
    for i in range(cells):
        a=begin+length*i/cells;b=begin+length*(i+1)/cells
        r=c.mpf([endpoints(a)[0],endpoints(b)[1]]);s=ell*r
        beta=raw_beta(c,r)/(ell*normalization);ds=ell*length/cells
        for row in (1,2):out[row-1]+=ds*c.exp((c.mpf('.5')-row*mu)*s)*beta
        out[2]+=ds*c.exp(-2*mu*s)*beta**2
    return out


class CompliantAxialPulseField:
    def __init__(self):
        self.selection=CompliantAxialAmplitude();self.ctx=c=self.selection.ctx
        self.pulse=self.selection.pulse;self.mu=self.selection.mu;self.delta=self.selection.future.delta
        self.hashes=dict(self.selection.hashes)
        for stem in ('compliant_axial_amplitude_selection','compliant_axial_amplitude_selection_check'):
            name=PREFIX+stem+'.json';r=json.loads((HERE/name).read_bytes())
            for source,digest in r['input_hashes'].items():
                if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:
                    raise ValueError('Pulse selected-source dependency changed: '+source)
                self.hashes[source]=digest
            if not r['actual_ap_selected']:raise ValueError('Actual selected amplitude required')
            self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        with mp.workdps(210):
            self.rate=1-self.mu;self.prate=1+2*self.mu
            self.logE=self.selection.log_end_scale
            self.logEcap=self.selection.future.params.log_mu-1000
            self.logE2cap=2*self.selection.future.params.log_mu-1000
            if endpoints(self.logEcap-self.logE)[0]<=0 or endpoints(self.logE2cap-2*self.logE)[0]<=0:
                raise ArithmeticError('Pulse end-factor cap not proved')
            self.Ecap=c.mpf([0,endpoints(c.exp(self.logEcap))[1]])
            self.E2cap=c.mpf([0,endpoints(c.exp(self.logE2cap))[1]])
            inlet=self.pulse.buffer.power('0',1)
            self.Xp=inlet['Mtheta_over_sqrt2_R_3half_Pstar'][0]/inlet['Utheta_over_Pstar'][0]
            self.beta_derivative_sup=self.selection.future.repair.beta_derivative_sup
            # The selector weights contain exp(-26) in Rp units. Partial
            # backward energy is in Rv units and uses the full local Gram.
            self.full_end_energy_weights=[c.exp(-2*self.mu*center)*self.pulse.basis['energy_gram']
                                          for center in (-3,-1)]
        self.hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    def factor(self,logvalue,cut=-1000):
        """Outward exponential enclosure; the exact logarithm remains source data."""
        c=self.ctx;cut=c.mpf(cut);low,high=endpoints(logvalue)
        if high<=endpoints(cut)[0]:return c.mpf([0,endpoints(c.exp(cut))[1]])
        # Very broad coordinate boxes can include t=0; only materialize
        # finite endpoints. The unrepresentable lower tail remains positive.
        lower=c.mpf(0) if low<-mp.mpf('1e25') else c.exp(c.mpf(low))
        if high>mp.mpf('1e25'):raise ValueError('Use a segmented coordinate with a finite upper log factor')
        upper=c.exp(c.mpf(high))
        return c.mpf([max(mp.mpf(0),endpoints(lower)[0]),endpoints(upper)[1]])

    def data(self,Z):
        c=self.ctx;selected=self.selection.select(Z);inlet=self.pulse.buffer.power(Z,1)
        jet=lambda key:IntervalTaylor(c,inlet[key])
        u=jet('Utheta_over_Pstar');e=jet('Mztheta_over_R_Pstar_squared')/(u*u)
        rows=[IntervalTaylor(c,j.coefficients) for j in selected['actual_incoming_moment_Taylor_enclosures']]
        return selected,inlet,u,e,rows

    def pressure_moment(self,Z,t,inlet,u,log_parts=None):
        c=self.ctx
        if log_parts is None:
            log_parts=dict(time_term=-self.prate*t)
        log_decay=sum(log_parts.values(),c.mpf(0));decay=self.factor(log_decay)
        kernel=(decay_integral(c,self.prate,t) if endpoints(self.prate*t)[1]<mp.mpf('1e-20')
                else (1-decay)/self.prate)
        p=IntervalTaylor(c,inlet['Mp_over_Pstar_squared'])+u*u*(kernel/2)
        p0rows=self.selection.future.angular.initial.datum.normalized_jets(c.mpf(Z),1)['normalized_pressure_coefficients']
        p0=IntervalTaylor(c,p0rows)
        return dict(Mp_over_Pstar_squared=p,P0_over_Pstar_squared=p0,P_over_Pstar_squared=p+p0,
            exact_pressure_source_unchanged=True,pressure_enclosure_is_absolute_not_a_relative_cone_bound=True,
            retained_swirl_pressure_decay_log_parts=log_parts,
            retained_swirl_pressure_decay_log_enclosure=log_decay)

    def radial(self,Z,B,m1):
        # Paper (3.8)-(3.9), where R is a similarity coordinate. Utheta_Z/
        # Utheta=-2Z/(1+Z^2) throughout O.4, so Mz_Z/(R Utheta) is
        # m1_Z-2Z*m1/(1+Z^2). This is not the cylindrical-radius formula.
        c=self.ctx;Z=c.mpf(Z);q=1+Z**2;d=1-Z**2;L=1-self.delta*Z**2
        return (2*Z*B[0]-(1-self.delta)*Z*m1[0]-d*(m1[1]-2*Z*m1[0]/q))/L

    def main(self,Z,xi,entrance_t=None,cells=256,window=1000):
        """Forward normalized partial primitives on the entire main support."""
        c=self.ctx;xi=c.mpf(xi)
        if endpoints(xi)[0]<0 or endpoints(xi)[1]>11:raise ValueError('Main xi in [0,11] required')
        selected,inlet,u,e0,rows=self.data(Z);ap=selected['selected_ap_Taylor']
        shape=gp(c,xi);B=ap*shape['value'];By=ap*(self.mu*shape['derivative'])
        t=xi/self.mu if entrance_t is None else c.mpf(entrance_t)
        length=min(endpoints(t)[0],mp.mpf(window));m=[]
        for row in (1,2):
            lam=c.mpf('.5')-row*self.mu
            kernel=c.mpf(0)
            for i in range(cells):
                a=c.mpf(length)*i/cells;b=c.mpf(length)*(i+1)/cells
                v=c.mpf([endpoints(a)[0],endpoints(b)[1]])
                kernel+=(c.exp(-lam*a)-c.exp(-lam*b))/lam*gp(c,xi-self.mu*v)['value']
            if endpoints(t)[1]>length:
                # All omitted main history is bounded by its positive
                # exponential mass; it is not asserted to be zero.
                kernel+=c.mpf([0,endpoints(11*c.exp(-lam*c.mpf(length))/lam)[1]])
            norm=sum((c.mpf(magnitude(v)) for v in rows[row-1].coefficients),c.mpf(0))
            cut=-1100-c.mpf(max(mp.mpf(0),endpoints(c.ln(norm))[1]))
            incoming=rows[row-1]*self.factor(-lam*t,cut)
            m.append(incoming+ap*kernel)
        X=1/self.rate+(self.Xp-1/self.rate)*self.factor(-self.rate*t)
        K=gp_energy(c,xi,self.selection.K)
        emu=(e0*self.mu+ap*ap*K-decay_integral(c,2,xi)/2)*c.exp(2*xi)
        e=emu/self.mu;energy_chart='forward from actual incoming energy'
        if endpoints(xi)[0]>=10:
            # Use the selected energy identity rather than subtracting two
            # independent enclosures of ap^2*K and the total swirl energy.
            D=13-xi;future=self.selection.future.future(Z)['complete_future_energy_Taylor']/2
            C=selected['selected_scaled_end_coefficient_Taylor'];end_energy=C[0]*0
            for Cj,Kj in zip(C,self.full_end_energy_weights):end_energy+=Cj*Cj*Kj
            e=(future*c.exp(-2*D)+decay_integral(c,2,D)/(2*self.mu)
               -ap*ap*(c.exp(2*xi)*gp_future_energy(c,xi)/self.mu)
               -end_energy*(c.exp(-2*D)*self.E2cap))
            energy_chart='backward from selected positive terminal target'
        return dict(stage='selected O.4 main pulse',Z=c.mpf(Z),coordinate=dict(kind='main',xi=xi,entrance_t=entrance_t),
            Uz_over_Utheta=B,Uz_y_over_Utheta=By-B*(c.mpf('.5')+self.mu),
            Mz_over_R_Utheta=m[0],Mtheta_z_over_sqrt2_R_3half_Utheta_squared=m[1],
            Mtheta_over_sqrt2_R_3half_Utheta=X,Mztheta_over_R_Utheta_squared=e,
            Ur_over_sqrt_R_over_2_Utheta=self.radial(Z,B,m[0]),Ur_Z_available=False,
            pressure=self.pressure_moment(Z,t,inlet,u),partial_gp_energy=K,partial_energy_chart=energy_chart,
            formal_log_Utheta_over_Pstar=dict(inlet_log=c.ln(u[0]),inverse_mu_term=-xi/(2*self.mu),finite_offset=-xi),
            exact_velocity_definition='Utheta=Ep*exp(-(.5+mu)t); Uz=Utheta*ap(Z)*gp(mu*t)',
            exact_forward_linear_moment_definition='mi(t)=exp(-lambda_i*t)*mi(0)+ap*integral_0^t exp(-lambda_i*(t-v))*gp(mu*v)dv',
            allfive_partial_primitives_callable=True,incoming_and_omitted_tails_retained=True,
            divergence_preserved_by_Mz_recovery=True,retained_axial_order=1,
            whole_outer_cone_certified=False,temporal_recursion=False)

    def entrance(self,Z,t):
        c=self.ctx;t=c.mpf(t)
        if endpoints(t)[0]<0:raise ValueError('Entrance t>=0 required')
        return self.main(Z,self.mu*t,entrance_t=t)

    def _gap(self,Z,D,leading_log,coordinate,pressure_log_parts):
        c=self.ctx;selected,inlet,u,e0,rows=self.data(Z)
        C=selected['selected_scaled_end_coefficient_Taylor'];zero=C[0]*0
        m=[];moment_logs=[];normal=self.pulse.initial.repair.normalization
        full_weights=backward_bump_weights(c,self.mu,normal,-4)
        weights=[full_weights,full_weights]
        for row in (1,2):
            lam=c.mpf('.5')-row*self.mu;total=zero
            logfactor=leading_log-row*D
            # leading_log is already reduced: no cancellation of two
            # rounded inverse-mu logs is used to recover this factor.
            factor=self.factor(logfactor,self.logEcap)
            for Cj,center,w in zip(C,(-3,-1),weights):total-=Cj*(c.exp(lam*center)*w[row-1])
            m.append(total*factor);moment_logs.append(logfactor)
        end_energy=zero
        for Cj,Kj in zip(C,self.full_end_energy_weights):end_energy+=Cj*Cj*Kj
        future=self.selection.future.future(Z)['complete_future_energy_Taylor']/2
        e=(future*c.exp(-2*D)+decay_integral(c,2,D)/(2*self.mu)
           -end_energy*(c.exp(-2*D)*self.E2cap))
        time_log=13/self.mu-D/self.mu
        angular_log_parts={key:value*(self.rate/self.prate) for key,value in pressure_log_parts.items()}
        X=1/self.rate+(self.Xp-1/self.rate)*self.factor(sum(angular_log_parts.values(),c.mpf(0)))
        if coordinate['kind']=='gap_main':
            xi=coordinate['xi'];velocity_log=dict(inlet_log=c.ln(u[0]),inverse_mu_term=-xi/(2*self.mu),finite_offset=-xi)
        else:
            s=coordinate['offset_from_Rv'];velocity_log=dict(inlet_log=c.ln(u[0]),inverse_mu_term=-13/(2*self.mu),finite_offset=-13-(c.mpf('.5')+self.mu)*s)
        return dict(stage='selected O.4 inactive gap',Z=c.mpf(Z),coordinate=coordinate,
            Uz_over_Utheta=zero,Uz_y_over_Utheta=zero,
            Mz_over_R_Utheta=m[0],Mtheta_z_over_sqrt2_R_3half_Utheta_squared=m[1],
            Mtheta_over_sqrt2_R_3half_Utheta=X,Mztheta_over_R_Utheta_squared=e,
            Ur_over_sqrt_R_over_2_Utheta=self.radial(Z,zero,m[0]),Ur_Z_available=False,
            pressure=self.pressure_moment(Z,time_log,inlet,u,pressure_log_parts),
            formal_log_Utheta_over_Pstar=velocity_log,retained_angular_memory_log_parts=angular_log_parts,
            exact_combined_moment_log_enclosures=moment_logs,
            exact_gap_moment_definition='mi=-sum_j Cj*Wji*exp(logE+lambda_i*(13-xi)/mu+lambda_i*center_j)',
            exact_gap_energy_definition='e=exp(-2D)*e_v+(1-exp(-2D))/(4mu)-exp(-2D)*exp(2logE)*sum_j Cj^2*Kj',
            inactive_gap_is_zero_velocity_not_zero_history=True,allfive_partial_primitives_callable=True,
            divergence_preserved_by_Mz_recovery=True,retained_axial_order=1,
            whole_outer_cone_certified=False,temporal_recursion=False)

    def gap(self,Z,xi):
        """Main-scale gap chart; xi=mu*log(R/Rp), 11<=xi<13."""
        c=self.ctx;xi=c.mpf(xi);D=13-xi
        if endpoints(xi)[0]<11 or endpoints(xi)[1]>=13 or endpoints(D/self.mu)[0]<4:
            raise ValueError('Gap xi>=11 and distance from Rv >=4 required; use gap_from_end near Rv')
        L=self.pulse.rows['saddle_L'];u0=self.pulse.rows['saddle_u0']
        finite=-3*L+2*c.ln(u0)-c.ln(6)/2-2*c.ln(self.mu)
        leading=(D/2-1)/self.mu+finite
        return self._gap(Z,D,leading,dict(kind='gap_main',xi=xi),
                         dict(inverse_mu_term=-self.prate*xi/self.mu))

    def gap_from_end(self,Z,s):
        """End-scale gap chart, -1/mu<=s<=-4; keeps finite offsets exact."""
        c=self.ctx;s=c.mpf(s)
        if endpoints(s)[1]>-4 or endpoints(s+1/self.mu)[0]<0:
            raise ValueError('Gap end chart requires -1/mu<=s<=-4')
        D=-self.mu*s;L=self.pulse.rows['saddle_L'];u0=self.pulse.rows['saddle_u0']
        finite=-3*L+2*c.ln(u0)-c.ln(6)/2-2*c.ln(self.mu)
        leading=-1/self.mu-s/2+finite
        return self._gap(Z,D,leading,dict(kind='gap_end',offset_from_Rv=s),
                         dict(inverse_mu_term=-13*self.prate/self.mu,finite_offset=-self.prate*s))

    def end(self,Z,s,cells=256):
        """Actual backward end field, s=log(R/Rv) in [-4,0]."""
        c=self.ctx;s=c.mpf(s)
        if endpoints(s)[0]<-4 or endpoints(s)[1]>0:raise ValueError('End s in [-4,0] required')
        selected,inlet,u,e0,rows=self.data(Z);C=selected['selected_scaled_end_coefficient_Taylor']
        zero=C[0]*0;Bhat=zero;Byhat=zero;mhat=[zero,zero];end_energy=zero
        ell=c.mpf('.15');normal=self.pulse.initial.repair.normalization
        for Cj,center in zip(C,(-3,-1)):
            local=s-center
            beta=raw_beta(c,local/ell)/(ell*normal)
            if endpoints(local)[1]<=-endpoints(ell)[1] or endpoints(local)[0]>=endpoints(ell)[1]:
                beta_y=c.mpf(0)
            elif endpoints(local)==(mp.mpf(0),mp.mpf(0)):beta_y=c.mpf(0)
            else:beta_y=c.mpf([-endpoints(self.beta_derivative_sup)[1],endpoints(self.beta_derivative_sup)[1]])
            Bhat+=Cj*beta;Byhat+=Cj*beta_y
            w=backward_bump_weights(c,self.mu,normal,local,cells)
            for row in (1,2):mhat[row-1]-=Cj*(c.exp((c.mpf('.5')-row*self.mu)*(center-s))*w[row-1])
            end_energy+=Cj*Cj*(c.exp(-2*self.mu*center)*w[2])
        B=Bhat*self.Ecap;m=[v*self.Ecap for v in mhat]
        future=self.selection.future.future(Z)['complete_future_energy_Taylor']/2
        # The selected energy equation fixes this positive terminal target.
        # Cancellation-free positive integral for the very small mu*s term.
        distance=-s
        baseline=c.mpf([endpoints(distance*c.exp(-2*self.mu*distance)/2)[0],endpoints(distance/2)[1]])
        e=future*c.exp(2*self.mu*s)+baseline-end_energy*(c.exp(2*self.mu*s)*self.E2cap)
        t=13/self.mu+s
        X=1/self.rate+(self.Xp-1/self.rate)*self.factor(-13*self.rate/self.mu-self.rate*s)
        terminal=endpoints(s)==(mp.mpf(0),mp.mpf(0))
        return dict(stage='selected O.4 end bumps',Z=c.mpf(Z),coordinate=dict(kind='end',offset_from_Rv=s),
            Uz_over_Utheta=B,Uz_y_over_Utheta=(Byhat-Bhat*(c.mpf('.5')+self.mu))*self.Ecap,
            formal_scaled_Uz_over_Utheta=Bhat,formal_scaled_Mz_mixed_moments=mhat,common_log_end_coefficient_scale=self.logE,
            Mz_over_R_Utheta=m[0],Mtheta_z_over_sqrt2_R_3half_Utheta_squared=m[1],
            Mtheta_over_sqrt2_R_3half_Utheta=X,Mztheta_over_R_Utheta_squared=e,
            Ur_over_sqrt_R_over_2_Utheta=self.radial(Z,B,m[0]),Ur_Z_available=False,
            pressure=self.pressure_moment(Z,t,inlet,u,
                dict(inverse_mu_term=-13*self.prate/self.mu,finite_offset=-self.prate*s)),
            positive_terminal_future_energy_Taylor=future,
            formal_log_Utheta_over_Pstar=dict(inlet_log=c.ln(u[0]),inverse_mu_term=-13/(2*self.mu),finite_offset=-13-(c.mpf('.5')+self.mu)*s),
            exact_backward_linear_moment_definition='mi(s)=-exp(logE)*sum_j Cj*integral_s^0 exp(lambda_i*(v-s))*beta(v-center_j)dv',
            exact_backward_energy_definition='e(s)=exp(2mu*s)*e(0)+integral_s^0 exp(2mu*(s-v))/2 dv-exp(2logE)*sum_j Cj^2*integral_s^0 exp(2mu*(s-v))*beta_j(v)^2 dv',
            terminal_linear_zero_by_empty_future_supports=terminal,positive_energy_target_at_Rv=terminal,
            allfive_partial_primitives_callable=True,terminal_identity_uses_selected_actual_moment_equations=True,
            divergence_preserved_by_Mz_recovery=True,retained_axial_order=1,
            whole_outer_cone_certified=False,temporal_recursion=False)

    def report(self):
        with mp.workdps(210):
            return dict(actual_five_defect_family_sha256=self.selection.future.angular.initial.family,
                implicit_source_sha256=self.selection.future.angular.initial.datum.source_sha,
                selected_mu=self.mu,selected_delta=self.delta,
                bump_normalization=self.pulse.initial.repair.normalization,
                actual_ap_selected=True,actual_main_and_end_pulse_partial_fields_callable=True,
                main_samples=[self.main('.5',x) for x in ('0','1','10.5','11')],
                entrance_samples=[self.entrance('.5',x) for x in ('0','1')],
                gap_samples=[self.gap('.5',x) for x in ('11','12','12.1','12.9')],
                gap_end_samples=[self.gap_from_end('.5',s) for s in (-self.ctx.mpf('.9')/self.mu,'-1000','-4')],
                end_samples=[self.end('.5',s) for s in ('-4','-3','-2','-1','0')],
                whole_Z_terminal=self.end([-1,1],0),
                positive_end_factor_cap=self.Ecap,positive_end_square_factor_cap=self.E2cap,
                exact_formal_end_factor_log=self.logE,
                full_end_energy_weights_in_Rv_units=self.full_end_energy_weights,
                selector_end_energy_weights_in_Rp_units=self.pulse.basis['end_energy_weights'],
                inactive_gap_partial_fields_callable=True,full_O4_interval_composed=True,
                gap_chart_coverage='gap xi in[11,12] plus gap_from_end s in[-1/mu,-4], joined at xi=12; end s in[-4,0]',
                radial_recovery_paper_equations='(3.8)-(3.9); R is a similarity coordinate',
                full_outer_five_moment_match=False,whole_outer_cone_certified=False,
                global_admissible_stress_lift_constructed=False,temporal_recursion=False,input_hashes=self.hashes)


def run():
    field=CompliantAxialPulseField();result=field.report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf-8')
    print('Actual selected O.4 pulse: entrance/main/gap/end, five partial primitives and paper (3.9) radial recovery generated; full outer/cone pending',flush=True)
    return result


if __name__=='__main__':run()
