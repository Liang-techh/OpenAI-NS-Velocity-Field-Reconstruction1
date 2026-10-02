"""Corrected post-Rv swirl, five primitives and direct infinite tail integrals.

The pulse's terminal axial/mixed zeros are inherited, never imposed afresh.
Angular and pressure histories remain forward cumulative quantities. A
backward pressure integral is reported independently until the absolute
preheat-datum identity has been composed over the entire matched field.
"""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_compliant_axial_pulse_field import CompliantAxialPulseField, sigma_enclosure
from lei_ren_part1_paper_compliant_outer_angular_candidate import unit_kernels
from lei_ren_part1_paper_compliant_outer_pulse_map import raw_beta
from lei_ren_part1_paper_compliant_outer_buffer import decay_integral as _decay_integral
from lei_ren_part1_paper_compliant_outer_angular_repair import magnitude, bump_weights
from lei_ren_part1_paper_interval_functional_defects import logjet
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_'


def decay_integral(c,rate,length):
    return _decay_integral(c,c.mpf(rate),c.mpf(length))


def future_bump_weights(c,mu,normalization,start,cells=256):
    """Direct remaining A/B/D/E/F integrals; no full-minus-past cancellation."""
    ell=c.mpf('.15');lo,hi=endpoints(c.mpf(start)/ell)
    begin=c.mpf([max(mp.mpf(-1),min(mp.mpf(1),lo)),max(mp.mpf(-1),min(mp.mpf(1),hi))])
    length=1-begin;out={key:c.mpf(0) for key in ('A','B','D','E','F')}
    for i in range(cells):
        a=begin+length*i/cells;b=begin+length*(i+1)/cells
        r=c.mpf([endpoints(a)[0],endpoints(b)[1]]);v=ell*r
        beta=raw_beta(c,r)/(ell*normalization);ds=ell*length/cells
        for key,rate,power in (('A',1-mu,1),('B',-1-2*mu,1),('D',-1-2*mu,2),('E',-2*mu,1),('F',-2*mu,2)):
            out[key]+=ds*c.exp(rate*v)*beta**power
    return out


def transition_tails(c,t,mu,delta,kind,cells=256):
    """Direct remaining energy/pressure kernels on one steep transition."""
    t=c.mpf(t);rate=(1-mu if kind=='in' else 1-delta/2)
    J,_=unit_kernels(c,t,rate,kind,cells)
    length=1-t;energy=c.mpf(0);pressure=c.mpf(0)
    for i in range(cells):
        a=t+length*i/cells;b=t+length*(i+1)/cells
        v=c.mpf([endpoints(a)[0],endpoints(b)[1]]);ds=length/cells
        nextJ=J+ds*sigma_enclosure(c,v);jc=c.mpf([endpoints(J)[0],endpoints(nextJ)[1]])
        if kind=='in':
            energy+=ds*c.exp(-2*mu*v-2*rate*jc)
            pressure+=ds*c.exp(-(1+2*mu)*v-2*rate*jc)/2
        else:
            energy+=ds*c.exp(-2*v+2*rate*jc)
            pressure+=ds*c.exp(-3*v+2*rate*jc)/2
        J=nextJ
    return dict(energy=energy,pressure=pressure)


class CompliantCorrectedOuterField:
    def __init__(self,cells=256):
        self.pulse=CompliantAxialPulseField();self.future=self.pulse.selection.future
        self.repair=self.future.repair;self.heat=self.future.heat;self.angular=self.future.angular
        self.ctx=c=self.future.ctx;self.params=self.future.params
        self.mu=self.future.mu;self.delta=self.future.delta;self.a=self.delta/2
        self.eps=self.heat.epsilon;self.rate=1-self.mu;self.prate=1+2*self.mu
        self.k=1-self.a;self.phrate=1+self.delta;self.bp=c.mpf('.5')+self.mu;self.bh=c.mpf('.5')+self.a
        self.Lrel=self.future.Lrel;self.cells=cells;self.cache={};self.hashes=dict(self.pulse.hashes)
        for stem in ('compliant_axial_pulse_field','compliant_axial_pulse_field_check'):
            name=PREFIX+stem+'.json';r=json.loads((HERE/name).read_bytes())
            for source,digest in r['input_hashes'].items():
                if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:raise ValueError('Outer pulse source changed: '+source)
                self.hashes[source]=digest
            if not r.get('full_O4_interval_composed',r.get('actual_O4_partial_pulse_field_installed',False)):
                raise ValueError('Actual complete O.4 pulse required')
            self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        with mp.workdps(210):
            self.Xv=self.pulse.end('0',0)['Mtheta_over_sqrt2_R_3half_Utheta']
            self.Nf=c.exp(-200*self.mu)/4;self.Pf=c.exp(-100*self.prate)/4
            self.Nrel=self.Nf*c.exp(-2*self.mu*self.Lrel);self.Prel=self.Pf*c.exp(-self.prate*self.Lrel)
            self.Ns=self.future.Ns;self.Ps=c.exp(-2-self.mu)
            self.Nq=self.future.Nq;self.Pq=self.Ps*c.exp(-3*self.params.Ts)
            self.Nt=self.future.Nt;self.Pt=self.Pq*c.exp(-2-self.a)
            self.Ntail=self.future.Ntail;self.Ptail=self.Pt*c.exp(-self.phrate*self.angular.waiting)
            self.tailmult=c.exp(-2*self.angular.waiting_logone)
            self.Lf=-100*self.bp-c.ln(2);self.LR=self.Lf-self.bp*self.Lrel
            self.Ls=self.LR-self.bp-self.rate/2;self.Lq=self.Ls-c.mpf('1.5')*self.params.Ts
            self.Lt=self.Lq-c.mpf('1.5')+self.k/2;self.Ltail=self.Lt-self.bh*self.angular.waiting
            self.fullin=transition_tails(c,0,self.mu,self.delta,'in',cells)
            self.fullout=transition_tails(c,0,self.mu,self.delta,'out',cells)
            _,self.Iin=unit_kernels(c,1,self.rate,'in',cells)
            _,self.Iout=unit_kernels(c,1,self.k,'out',cells)
            exactlog=2*(self.angular.loguRp0-13/(2*self.mu)-13)
            self.logEv0cap=2*self.params.log_mu-1000
            if endpoints(self.logEv0cap-exactlog)[0]<=0:raise ArithmeticError('Ev0 physical scale cap failed')
            self.Ev0Pstar2=c.mpf([0,endpoints(c.exp(self.logEv0cap))[1]])
            self.logEv0Pstar2_parts=dict(inlet_log=2*self.angular.loguRp0,inverse_mu_term=-13/self.mu,finite_offset=c.mpf(-26))
        self.hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    def jet(self,value):return IntervalTaylor(self.ctx,[self.ctx.mpf(value),0])

    def Zkey(self,Z):return tuple(endpoints(self.ctx.mpf(Z)))

    def phase(self,value):
        t=self.ctx.mpf(value)
        if endpoints(t)[0]<0 or endpoints(t)[1]>1:raise ValueError('Phase in [0,1] required')
        return t

    def flatten_raw(self,Z,phase):
        c=self.ctx;phase=self.phase(phase);Z=c.mpf(Z);t=100*phase
        q=IntervalTaylor(c,[1+Z**2,2*Z]);lq=logjet(q);ratio=lq-c.ln(2)
        sig=sigma_enclosure(c,phase);F=(ratio*sig).exp();integral=q*0
        for i in range(self.cells):
            a=t*i/self.cells;b=t*(i+1)/self.cells;v=c.mpf([endpoints(a)[0],endpoints(b)[1]])
            integral+=(ratio*sigma_enclosure(c,v/100)).exp()*((c.exp(-self.rate*(t-b))-c.exp(-self.rate*(t-a)))/self.rate)
        X=(integral+self.Xv*c.exp(-self.rate*t))/F
        L=ratio*sig-lq-self.bp*t
        if endpoints(phase)==(mp.mpf(1),mp.mpf(1)):L=self.jet(self.Lf)
        return X,L,t

    def flatten_tails(self,Z,t):
        c=self.ctx;Z=c.mpf(Z);q=IntervalTaylor(c,[1+Z**2,2*Z]);ratio=logjet(q)-c.ln(2)
        length=100-t;energy=q*0;pressure=q*0
        for i in range(self.cells):
            a=t+length*i/self.cells;b=t+length*(i+1)/self.cells
            v=c.mpf([endpoints(a)[0],endpoints(b)[1]]);ds=length/self.cells
            F=(ratio*(2*sigma_enclosure(c,v/100))).exp()/(q*q)
            energy+=F*(ds*c.exp(-2*self.mu*v));pressure+=F*(ds*c.exp(-self.prate*v)/2)
        return energy,pressure

    def heat_local(self,Z,t):
        """Pure Gamma exterior in current-radius units, even at huge offsets."""
        c=self.ctx;Z=c.mpf(Z);t=c.mpf(t);d=1-Z**2
        Sbox=self.repair.strong_S_cap*self.pulse.factor(-t)
        Sbox=c.mpf([0,endpoints(Sbox)[1]]);B=(1+self.a)**2*(2+self.a)
        F=c.mpf([endpoints(1+self.a-B*Sbox)[0],endpoints(1+self.a)[1]])
        G=c.mpf([endpoints(1+self.a-2*B*Sbox)[0],endpoints(1+self.a)[1]])
        D=IntervalTaylor(c,[2*d*F,-4*Z*G]);K=self.jet(1)-D*(self.a*Sbox)
        absZ=max(abs(v) for v in endpoints(Z))
        def enclosed(lead,err,derivative,errZ):
            value=c.mpf([max(mp.mpf(0),endpoints(lead-err)[0]),max(mp.mpf(0),endpoints(lead)[1])])
            return IntervalTaylor(c,[value,derivative+c.mpf([-endpoints(errZ)[1],endpoints(errZ)[1]])])
        th=enclosed(2*d*(1+self.a),2*self.a*(1+self.a)*(2+self.a)*d**2*Sbox,
                    -4*Z*(1+self.a),8*self.a*(1+self.a)*(2+self.a)*absZ*Sbox)
        def square(rate):
            return enclosed(4*d*(1+self.a)/(rate+1),4*(B+self.a*(1+self.a)**2)*d**2*Sbox/(rate+2),
                            -8*Z*(1+self.a)/(rate+1),16*(B+self.a*(1+self.a)**2)*absZ*Sbox/(rate+2))
        ph=square(self.phrate)/2;eh=square(self.delta)
        X=(th*Sbox+1/self.k)/K
        Pbar=(ph*(-self.a*Sbox)+1/(2*self.phrate))/(K*K)
        Ebar=(eh*(-self.a*Sbox)+1/self.delta)/(K*K)
        return dict(K=K,X_target=X,pressure_tail_over_Utheta_squared=Pbar,
                    energy_tail_over_R_Utheta_squared=Ebar,theta_hat=th,pressure_hat=ph,energy_hat=eh,
                    current_inverse_radius_cap=Sbox,
                    exact_current_inverse_radius='S_current=(1/Rtail)*exp(-offset), never its cap')

    def collar_shape(self,Z,t):
        c=self.ctx;t=c.mpf(t);lo,hi=endpoints(t)
        if lo<0 or hi>3:raise ValueError('Collar offset in[0,3] required')
        sig=sigma_enclosure(c,t)
        def phi_point(v):
            if v>=3:return c.mpf(0)
            logv=-4/(3-c.mpf(v))**2
            return self.pulse.factor(logv)
        phi=c.mpf([endpoints(phi_point(hi))[0],endpoints(phi_point(lo))[1]])
        W=1-sig+sig*phi;C=1-self.eps*phi;pre=1-self.eps*W
        D=IntervalTaylor(c,self.heat.deficit(Z,t)['deficit_scaled_Taylor'].coefficients)*(sig*C)
        S=c.mpf([0,endpoints(self.repair.strong_S_cap)[1]])
        return dict(K=self.jet(pre)-D*(self.a*S),D=D,W=W,pre=pre,sigma=sig,C=C)

    def heat_tails(self,Z,t):
        """Arbitrary collar/pure-heat backward targets, in Rtail units."""
        c=self.ctx;t=c.mpf(t)
        if endpoints(t)[0]<0:raise ValueError('Heat offset>=0 required')
        if endpoints(t)[0]>=3:
            loc=self.heat_local(Z,t);K=loc['K'];S=loc['current_inverse_radius_cap']
            E=(loc['energy_hat']*(-self.a*S)+1/self.delta)*self.pulse.factor(-self.delta*t)
            P=(loc['pressure_hat']*(-self.a*S)+1/(2*self.phrate))*self.pulse.factor(-self.phrate*t)
            return dict(K=K,X_target=loc['X_target'],energy=E,pressure=P,pure_heat=loc,
                        theta_target_local=loc['X_target']*K)
        loc3=self.heat_local(Z,3);S=c.mpf([0,endpoints(self.repair.strong_S_cap)[1]])
        # Exterior terms are normalized at Rtail; S_current(3)=S*exp(-3).
        theta_hat=loc3['theta_hat']*c.exp(-3*self.a)
        pressure_hat=loc3['pressure_hat']*c.exp(-3*(self.phrate+1))
        energy_hat=loc3['energy_hat']*c.exp(-3*(self.delta+1))
        atoms={key:c.mpf(0) for key in ('JW','PW','PW2','EW','EW2')};length=3-t
        zero=theta_hat*0
        for i in range(self.cells):
            left=t+length*i/self.cells;right=t+length*(i+1)/self.cells
            # Clamp outward rounding at the flat endpoint; every true cell
            # stays in [t,3], and the integrands are defined there.
            v=c.mpf([max(mp.mpf(0),endpoints(left)[0]),min(mp.mpf(3),endpoints(right)[1])]);ds=length/self.cells
            sh=self.collar_shape(Z,v);W=sh['W'];D=sh['D']
            sq=zero+2*sh['pre']-D*(self.a*S)
            theta_hat+=D*(ds*self.a*c.exp(self.k*v))
            pressure_hat+=D*sq*(ds*c.exp(-self.phrate*v)/2)
            energy_hat+=D*sq*(ds*c.exp(-self.delta*v))
            atoms['JW']+=ds*c.exp(self.k*v)*W
            for prefix,rate in (('P',self.phrate),('E',self.delta)):
                atoms[prefix+'W']+=ds*c.exp(-rate*v)*W;atoms[prefix+'W2']+=ds*c.exp(-rate*v)*W**2
        K=self.collar_shape(Z,t)['K']
        E=(energy_hat*(-self.a*S)+c.exp(-self.delta*t)/self.delta
           -2*self.eps*atoms['EW']+self.eps**2*atoms['EW2'])
        P=(pressure_hat*(-self.a*S)+c.exp(-self.phrate*t)/(2*self.phrate)
           -self.eps*atoms['PW']+self.eps**2*atoms['PW2']/2)
        theta=theta_hat*(S*c.exp(-self.k*t))+1/self.k+self.eps*atoms['JW']*c.exp(-self.k*t)
        return dict(K=K,X_target=theta/K,theta_target_local=theta,
                    energy=E,pressure=P,scaled_heat_deficit_tails=dict(theta=theta_hat,pressure=pressure_hat,energy=energy_hat),
                    separate_epsilon_future_atoms=atoms)

    def data(self,Z):
        key=self.Zkey(Z)
        if key in self.cache:return self.cache[key]
        c=self.ctx;Z=c.mpf(Z);coeff=self.repair.coefficients(Z)['physical_coefficient_Taylor']
        coeff=[IntervalTaylor(c,v.coefficients) for v in coeff]
        zero=coeff[0]*0;fullA=zero
        for dj,center in zip(coeff,(-3,-1)):fullA+=dj*(c.exp(self.rate*center)*self.repair.weights['A'])
        Xf=self.flatten_raw(Z,1)[0];eq=1/self.rate
        XR=(Xf-eq)*c.exp(-self.rate*self.Lrel)+fullA+eq
        Xs=(XR+self.Iin)*c.exp(-self.rate/2);Xq=Xs+self.params.Ts
        Xt=(Xq+self.Iout)*c.exp(-self.k/2)
        Xtail=(Xt-1/self.k)*c.exp(-self.k*self.angular.waiting)+1/self.k
        heat0=self.heat_tails(Z,0)
        Epost=(heat0['energy']*(self.Ntail*self.tailmult)
               +self.fullin['energy']+self.Ns*decay_integral(c,2,self.params.Ts)
               +self.Nq*self.fullout['energy']+self.Nt*decay_integral(c,self.delta,self.angular.waiting))
        Ppost=(heat0['pressure']*(self.Ptail*self.tailmult)
               +self.fullin['pressure']+self.Ps*decay_integral(c,3,self.params.Ts)/2
               +self.Pq*self.fullout['pressure']+self.Pt*decay_integral(c,self.phrate,self.angular.waiting)/2)
        terminal=self.pulse.end(Z,0)
        item=dict(coeff=coeff,Xf=Xf,XR=XR,Xs=Xs,Xq=Xq,Xt=Xt,Xtail=Xtail,heat0=heat0,Epost=Epost,Ppost=Ppost,
                  pulse_terminal=terminal,angular_tail_constant_defect=(Xtail*(1-self.eps)-heat0['theta_target_local']))
        self.cache[key]=item
        future=self.correction(Z,-4,backward=True)
        Ef=Epost*self.Nrel+future['E']*self.Nrel+self.Nf*decay_integral(c,2*self.mu,self.Lrel)
        Pf=Ppost*self.Prel+future['P']*self.Prel+self.Pf*decay_integral(c,self.prate,self.Lrel)/2
        fe,fp=self.flatten_tails(Z,c.mpf(0));item['Erv']=Ef+fe;item['Prv']=Pf+fp
        return item

    def correction(self,Z,s,backward=False):
        c=self.ctx;s=c.mpf(s);coeff=self.data(Z)['coeff'];zero=coeff[0]*0
        h=zero;A=zero;P=zero;E=zero
        for dj,center in zip(coeff,(-3,-1)):
            local=s-center;ell=c.mpf('.15');lo,hi=endpoints(local)
            # raw_beta handles interval intersections with the entire
            # support, including boxes whose lower endpoint is outside it.
            h+=dj*(raw_beta(c,local/ell)/(ell*self.repair.normalization))
            if backward:w=future_bump_weights(c,self.mu,self.repair.normalization,local,self.cells)
            elif hi<=-endpoints(ell)[1]:w={n:c.mpf(0) for n in self.repair.weights}
            elif lo>=endpoints(ell)[1]:w=self.repair.weights
            else:w=bump_weights(c,self.mu,self.repair.normalization,local,self.cells)
            A+=dj*(c.exp(self.rate*center)*w['A'])
            P+=(dj*w['B']+dj*dj*w['D']/2)*c.exp(-self.prate*center)
            E+=(dj*(2*w['E'])+dj*dj*w['F'])*c.exp(-2*self.mu*center)
        return dict(h=h,A=A,P=P,E=E)

    def packet(self,Z,stage,coordinate,y,logtheta,X,Eglobal,Pglobal,local_heat=None):
        c=self.ctx;data=self.data(Z);zero=self.jet(0)
        theta=logtheta.exp() if local_heat is None else local_heat['theta']
        if local_heat is None:
            energy=Eglobal/(theta*theta)*(c.exp(-y)/2)
            pressureback=-Pglobal/(theta*theta)
        else:
            energy=local_heat['energy'];pressureback=local_heat['pressure']
        if endpoints(Eglobal[0])[0]<0 or endpoints(Pglobal[0])[0]<0:
            raise ArithmeticError('Positive direct outer tail lost')
        if (endpoints(Eglobal[0])[0]==0 or endpoints(Pglobal[0])[0]==0) and local_heat is None:
            raise ArithmeticError('Direct tail positivity unresolved outside the formal pure-heat chart')
        if local_heat is not None and (endpoints(energy[0])[0]<=0 or endpoints(pressureback[0])[1]>=0):
            raise ArithmeticError('Exact positive pure-heat factors lost')
        # Same actual history, not a replacement Mp obtained from P=-tail.
        prefix=data['Prv']-Pglobal
        if endpoints(prefix[0])[1]<0:raise ArithmeticError('Positive accumulated pressure history has a disjoint negative enclosure')
        if coordinate['origin']=='Rv' and endpoints(coordinate['offset'])==(mp.mpf(0),mp.mpf(0)):
            prefix=zero
        else:
            prefix=IntervalTaylor(c,[c.mpf([max(mp.mpf(0),endpoints(prefix[0])[0]),endpoints(prefix[0])[1]]),prefix[1]])
        terminal=data['pulse_terminal']['pressure']
        Mp=IntervalTaylor(c,terminal['Mp_over_Pstar_squared'].coefficients)+prefix*self.Ev0Pstar2
        P0=IntervalTaylor(c,terminal['P0_over_Pstar_squared'].coefficients);Pforward=Mp+P0
        return dict(stage=stage,Z=c.mpf(Z),coordinate=coordinate,
            formal_log_Utheta_over_Ev0=logtheta,Utheta_over_Ev0_enclosure=theta,
            Uz=zero,Ur=zero,Mz_over_R_Utheta=zero,Mtheta_z_over_sqrt2_R_3half_Utheta_squared=zero,
            Mtheta_over_sqrt2_R_3half_Utheta=X,Mztheta_over_R_Utheta_squared=energy,
            Mp_over_Pstar_squared=Mp,P0_over_Pstar_squared=P0,P_over_Pstar_squared=Pforward,
            independent_backward_pressure_over_Utheta_squared=pressureback,
            remaining_swirl_energy_in_Rv_Ev0_squared=Eglobal,remaining_pressure_integral_in_Ev0_squared=Pglobal,
            forward_Mp_increment_in_Ev0_squared=prefix,
            physical_Ev0_over_Pstar_squared_scale_log_parts=self.logEv0Pstar2_parts,
            exact_pressure_definition='P=P0+Mp; backward integral is independent until global constant offset is proven zero',
            exact_energy_definition='Mztheta=.5*integral_R^infinity Utheta_corrected^2 drho, inherited from selected pulse',
            allfive_partial_primitives_callable=True,post_pulse_zero_histories_inherited=True,
            divergence_preserved_by_paper_3_9=True,Ur_Z_available=True,
            retained_axial_order=1,whole_outer_cone_certified=False,
            exact_tail_integrals_positive_even_if_the_cap_includes_zero=True,
            full_global_pressure_terminal_identity_verified=False,full_outer_five_moment_match=False,temporal_recursion=False)

    def flatten(self,Z,phase):
        c=self.ctx;X,L,t=self.flatten_raw(Z,phase);data=self.data(Z)
        fe,fp=self.flatten_tails(Z,t);fc=self.correction(Z,-4,True)
        E=fe+self.Nf*decay_integral(c,2*self.mu,self.Lrel)+(data['Epost']+fc['E'])*self.Nrel
        P=fp+self.Pf*decay_integral(c,self.prate,self.Lrel)/2+(data['Ppost']+fc['P'])*self.Prel
        return self.packet(Z,'corrected O.5 flatten',dict(origin='Rv',phase=c.mpf(phase),offset=t),t,L,X,E,P)

    def power_buffer(self,Z,phase):
        c=self.ctx;phase=self.phase(phase);t=(self.Lrel-4)*phase;data=self.data(Z);eq=1/self.rate
        X=(data['Xf']-eq)*c.exp(-self.rate*t)+eq;L=self.jet(self.Lf-self.bp*t)
        fc=self.correction(Z,-4,True)
        E=(data['Epost']+fc['E'])*self.Nrel+self.Nf*c.exp(-2*self.mu*t)*decay_integral(c,2*self.mu,self.Lrel-t)
        P=(data['Ppost']+fc['P'])*self.Prel+self.Pf*c.exp(-self.prate*t)*decay_integral(c,self.prate,self.Lrel-t)/2
        return self.packet(Z,'corrected O.6 power buffer',dict(origin='Rf',phase=phase,offset=t),100+t,L,X,E,P)

    def angular_end(self,Z,s):
        c=self.ctx;s=c.mpf(s)
        if endpoints(s)[0]<-4 or endpoints(s)[1]>0:raise ValueError('Angular end offset in[-4,0] required')
        data=self.data(Z);past=self.correction(Z,s);future=self.correction(Z,s,True);eq=1/self.rate
        Xbase=(data['Xf']-eq)*c.exp(-self.rate*(self.Lrel+s))+eq
        X=(Xbase+past['A']*c.exp(-self.rate*s))/(past['h']+1)
        L=logjet(past['h']+1)+(self.LR-self.bp*s)
        E=(data['Epost']+future['E']+c.exp(-2*self.mu*s)*decay_integral(c,2*self.mu,-s))*self.Nrel
        P=(data['Ppost']+future['P']+c.exp(-self.prate*s)*decay_integral(c,self.prate,-s)/2)*self.Prel
        return self.packet(Z,'corrected O.6 angular bumps',dict(origin='Rrel',offset=s),100+self.Lrel+s,L,X,E,P)

    def steep_in(self,Z,t):
        c=self.ctx;t=self.phase(t);data=self.data(Z);J,I=unit_kernels(c,t,self.rate,'in',self.cells)
        X=(data['XR']+I)*c.exp(-self.rate*(t-J));L=self.jet(self.LR-self.bp*t-self.rate*J)
        tails=transition_tails(c,t,self.mu,self.delta,'in',self.cells)
        E=(data['heat0']['energy']*(self.Ntail*self.tailmult)+tails['energy']
           +self.Ns*decay_integral(c,2,self.params.Ts)+self.Nq*self.fullout['energy']
           +self.Nt*decay_integral(c,self.delta,self.angular.waiting))*self.Nrel
        P=(data['heat0']['pressure']*(self.Ptail*self.tailmult)+tails['pressure']
           +self.Ps*decay_integral(c,3,self.params.Ts)/2+self.Pq*self.fullout['pressure']
           +self.Pt*decay_integral(c,self.phrate,self.angular.waiting)/2)*self.Prel
        return self.packet(Z,'corrected O.7 steep entry',dict(origin='Rrel',offset=t),100+self.Lrel+t,L,X,E,P)

    def steep_power(self,Z,phase):
        c=self.ctx;phase=self.phase(phase);t=self.params.Ts*phase;left=self.params.Ts*(1-phase);data=self.data(Z)
        X=data['Xs']+t;L=self.jet(self.Ls-c.mpf('1.5')*t)
        E=(data['heat0']['energy']*(self.Ntail*self.tailmult)+self.Nq*self.fullout['energy']
           +self.Nt*decay_integral(c,self.delta,self.angular.waiting)+self.Ns*c.exp(-2*t)*decay_integral(c,2,left))*self.Nrel
        P=(data['heat0']['pressure']*(self.Ptail*self.tailmult)+self.Pq*self.fullout['pressure']
           +self.Pt*decay_integral(c,self.phrate,self.angular.waiting)/2+self.Ps*c.exp(-3*t)*decay_integral(c,3,left)/2)*self.Prel
        return self.packet(Z,'corrected O.7 steep power',dict(origin='Rs',phase=phase,offset=t),101+self.Lrel+t,L,X,E,P)

    def steep_out(self,Z,t):
        c=self.ctx;t=self.phase(t);data=self.data(Z);J,I=unit_kernels(c,t,self.k,'out',self.cells)
        X=(data['Xq']+I)*c.exp(-self.k*J);L=self.jet(self.Lq-c.mpf('1.5')*t+self.k*J)
        tails=transition_tails(c,t,self.mu,self.delta,'out',self.cells)
        E=(data['heat0']['energy']*(self.Ntail*self.tailmult)+self.Nt*decay_integral(c,self.delta,self.angular.waiting)+self.Nq*tails['energy'])*self.Nrel
        P=(data['heat0']['pressure']*(self.Ptail*self.tailmult)+self.Pt*decay_integral(c,self.phrate,self.angular.waiting)/2+self.Pq*tails['pressure'])*self.Prel
        return self.packet(Z,'corrected O.7 steep exit',dict(origin='Rq',offset=t),101+self.Lrel+self.params.Ts+t,L,X,E,P)

    def waiting(self,Z,phase):
        c=self.ctx;phase=self.phase(phase);t=self.angular.waiting*phase;left=self.angular.waiting*(1-phase);data=self.data(Z)
        X=(data['Xt']-1/self.k)*c.exp(-self.k*t)+1/self.k;L=self.jet(self.Lt-self.bh*t)
        E=(data['heat0']['energy']*(self.Ntail*self.tailmult)+self.Nt*c.exp(-self.delta*t)*decay_integral(c,self.delta,left))*self.Nrel
        P=(data['heat0']['pressure']*(self.Ptail*self.tailmult)+self.Pt*c.exp(-self.phrate*t)*decay_integral(c,self.phrate,left)/2)*self.Prel
        return self.packet(Z,'corrected O.7 waiting',dict(origin='Rt',phase=phase,offset=t),102+self.Lrel+self.params.Ts+t,L,X,E,P)

    def collar_heat(self,Z,t):
        c=self.ctx;t=c.mpf(t);data=self.data(Z);heat=self.heat_tails(Z,t);K=heat['K']
        X=heat['X_target']+data['angular_tail_constant_defect']*(self.pulse.factor(-self.k*t))/K
        L=logjet(K)+(self.Ltail-self.angular.waiting_logone-self.bh*t)
        E=heat['energy']*(self.Nrel*self.Ntail*self.tailmult)
        P=heat['pressure']*(self.Prel*self.Ptail*self.tailmult)
        y=102+self.Lrel+self.params.Ts+self.angular.waiting+t
        if endpoints(t)[0]>=3:
            local=heat['pure_heat'];theta=K*self.pulse.factor(self.Ltail-self.angular.waiting_logone-self.bh*t)
            overrides=dict(theta=theta,energy=local['energy_tail_over_R_Utheta_squared']/2,
                           pressure=-local['pressure_tail_over_Utheta_squared'])
        else:overrides=None
        packet=self.packet(Z,'corrected collar / exact Gamma exterior',dict(origin='Rtail',offset=t),y,L,X,E,P,overrides)
        packet.update(heat_profile_definition=self.heat.heat_definition,
            heat_bracket_Taylor=K,
            segmented_log_Utheta_over_Ev0=dict(waiting_endpoint_log=self.Ltail,
                minus_log_one_minus_epsilon=-self.angular.waiting_logone,
                current_power_log=-self.bh*t,log_heat_bracket_Taylor=logjet(K)),
            exact_current_inverse_radius_log_parts=dict(**self.heat.logS_terms,current_offset=-t),
            exact_Gamma_profile_retained=True,angular_constant_defect_not_reset=data['angular_tail_constant_defect'],
            forward_angular_target_equality_verified=False,exact_heat_velocity_from_offset_3=True,
            entire_infinite_tail_integrals_included=True)
        return packet

    def report(self):
        with mp.workdps(210):
            samples=[self.flatten('.5',p) for p in (0,'.5',1)]
            samples +=[self.power_buffer('.5',p) for p in (0,1)]
            samples +=[self.angular_end('.5',s) for s in (-4,-3,-2,-1,0)]
            samples +=[self.steep_in('.5',p) for p in (0,1)]
            samples +=[self.steep_power('.5',p) for p in (0,'.5',1)]
            samples +=[self.steep_out('.5',p) for p in (0,1)]
            samples +=[self.waiting('.5',p) for p in (0,1)]
            samples +=[self.collar_heat('.5',t) for t in (0,1,3,4,1000)]
            whole=[self.flatten([-1,1],'.5'),self.angular_end([-1,1],-3),self.collar_heat([-1,1],3)]
            whole_Rv=self.flatten([-1,1],0)
            wide_support=self.angular_end('.5',[-3.2,-3])
            data=self.data('.5');pulse=data['pulse_terminal'];Z=self.ctx.mpf('.5')
            q=IntervalTaylor(self.ctx,[1+Z**2,2*Z])
            return dict(actual_five_defect_family_sha256=self.angular.initial.family,
                implicit_source_sha256=self.angular.initial.datum.source_sha,samples=samples,whole_Z_examples=whole,
                angular_interval_crossing_support=wide_support,
                whole_Z_Rv_interface=whole_Rv,
                whole_Z_pulse_terminal=self.data([-1,1])['pulse_terminal'],
                source_pulse_terminal=pulse,post_Rv_energy_tail_in_Ev0_units=data['Erv'],
                pulse_future_energy_in_Ev0_units=pulse['positive_terminal_future_energy_Taylor']*2/(q*q),
                future_pressure_at_Rv_in_Ev0_units=data['Prv'],
                formal_pressure_infinity_offset_atoms=dict(
                    P0_over_Pstar_squared=pulse['pressure']['P0_over_Pstar_squared'],
                    Mp_Rv_over_Pstar_squared=pulse['pressure']['Mp_over_Pstar_squared'],
                    remaining_pressure_over_Ev0_squared=data['Prv'],
                    physical_Ev0_over_Pstar_squared_scale_log_parts=self.logEv0Pstar2_parts),
                angular_tail_constant_defect=data['angular_tail_constant_defect'],
                formal_pressure_infinity_offset='P0+Mp(Rv)+Ev0^2*Prv; not yet proven zero by an absolute source identity',
                post_Rv_corrected_velocity_and_five_primitives_callable=True,
                direct_arbitrary_radius_energy_pressure_heat_tails_callable=True,
                post_Rv_axial_and_radial_zero_inherited_from_actual_pulse=True,
                post_Rv_energy_moment_terminal_zero_from_delta_positive=True,
                full_global_pressure_terminal_identity_verified=False,forward_angular_heat_target_equality_verified=False,
                full_outer_five_moment_match=False,C4_outer_bounds_certified=False,whole_outer_cone_certified=False,
                global_admissible_stress_lift_constructed=False,temporal_recursion=False,input_hashes=self.hashes)


def run():
    field=CompliantCorrectedOuterField();result=field.report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf-8')
    print('Corrected post-Rv outer fields and five primitives through infinite Gamma heat tails generated; absolute pressure/angular terminal identities and C4/cone pending',flush=True)
    return result


if __name__=='__main__':run()
