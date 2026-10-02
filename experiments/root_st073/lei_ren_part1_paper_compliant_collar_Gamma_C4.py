"""Original flat heat collar and full Gamma exterior, leading mixed C4.

Gamma derivatives are enclosed directly by positive expectation derivatives.
Finite Taylor composition computes derivative coefficients only; it does
not define the heat function by a finite inverse-radius expansion.
"""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
from lei_ren_part1_paper_compliant_steep_waiting_C4 import CompliantSteepWaitingC4
from lei_ren_part1_paper_compliant_power_angular_C4 import quotient_log_rates
from lei_ren_part1_paper_compliant_angular_high_jets import integrated_gamma_tails,pochhammer,log_taylor
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import sigma_jets,flat_tail_exp,positive_exp,polynomial_add,polynomial_mul
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_'


def phi_polynomials():
    rows=[[1]]
    for n in range(4):
        p=rows[-1]; derivative=[j*p[j] for j in range(1,len(p))] or [0]
        rows.append(polynomial_add(polynomial_add([v*(-8) for v in p],polynomial_mul([0,0,3*n],p)),
            polynomial_mul([0,0,0,-1],derivative)))
    return rows


PHI_POLYNOMIALS=phi_polynomials()


def phi_jets(c,t):
    """Original exp(-4/(3-t)^2), flat zero for t>=3; y Taylor0..4."""
    t=c.mpf(t); lo,hi=endpoints(t)
    if lo<0:raise ValueError('Nonnegative collar offset required')
    if lo>=3:return IntervalTaylor.constant(c,0,4)
    # Perform endpoint subtraction inside the directed context even when
    # this public provider is called outside report()'s mp.workdps scope.
    D=c.mpf(endpoints(c.mpf(3)-c.mpf(lo))[1])
    log_upper=-4/D**2; flat=endpoints(log_upper)[1]<-1000
    upper=flat_tail_exp(c,log_upper,D) if flat else c.exp(log_upper)
    lower=c.mpf(0) if hi>=3 or flat else positive_exp(c,-4/(c.mpf(3)-c.mpf(hi))**2)
    value=c.mpf([endpoints(lower)[0],endpoints(upper)[1]])
    if hi>=3 or flat:
        rows=[value]
        for n,p in enumerate(PHI_POLYNOMIALS[1:],1):
            coefficient=sum(abs(v)*3**j for j,v in enumerate(p))
            peak=c.sqrt(c.mpf(8)/(3*n)); rho=D if endpoints(D)[1]<endpoints(peak)[0] else peak
            bound=flat_tail_exp(c,c.ln(coefficient)-4/rho**2-3*n*c.ln(rho),D)
            radius=endpoints(bound)[1]; rows.append(c.mpf([-radius,radius])/math.factorial(n))
        return IntervalTaylor(c,rows)
    dist=IntervalTaylor(c,[3-t,-1,0,0,0])
    jets=(dist**(-2)*(-4)).exp()
    return IntervalTaylor(c,[value]+list(jets.coefficients[1:]))


def positive_moment_derivative(c,a,xi,n,divided=False):
    """Signed derivative bounds of exact positive Gamma expectations.

    divided=True: F=(1-H)/(a*xi)=E[v int_0^1(1+q*xi*v)^(-a-1)dq].
    Otherwise n>=1: H^(n)/a. Bounds follow from 0<=1-(1+x)^(-b)<=b*x,
    integrated against the full Gamma density. No truncated H definition.
    """
    if divided:
        M=pochhammer(c,1+a,n)*pochhammer(c,1+a,n+1)/(n+1)
        nextM=pochhammer(c,1+a,n+1)*pochhammer(c,1+a,n+2)/(n+2)
        sign=(-1)**n
    else:
        if n<1:raise ValueError('Positive H derivative order required')
        M=pochhammer(c,1+a,n-1)*pochhammer(c,1+a,n)
        nextM=pochhammer(c,1+a,n)*pochhammer(c,1+a,n+1)
        sign=(-1)**n
    low=max(mp.mpf(0),endpoints(M-nextM*xi)[0]); high=endpoints(M)[1]
    return c.mpf([low,high])*sign


def stirling_second(n,k):
    if n==k==0:return 1
    if k==0 or k>n:return 0
    return stirling_second(n-1,k-1)+k*stirling_second(n-1,k)


def gamma_deficit_mixed(c,Z,a,S_cap,t,yorder=4):
    """Ordinary y derivatives0..4 of Dhat, each true axial Taylor0..5.

    Dhat=(1-H(2*d*S*exp(-t)))/(a*S), d=1-Z^2. For k>=1, powers
    xi^l/S are cancelled symbolically as (2*d*exp(-t))^l*S^(l-1),
    so no division by a numerical cap containing zero occurs.
    """
    Z=c.mpf(Z); t=c.mpf(t)
    if endpoints(Z)[0]<-1 or endpoints(Z)[1]>1 or endpoints(t)[0]<0:
        raise ValueError('Actual Gamma Z/offset domain required')
    decay=c.exp(-t); d=IntervalTaylor(c,[1-Z**2,-2*Z,-1,0,0,0])
    S=c.mpf([0,endpoints(S_cap)[1]]); base=d*(2*decay); xi=base*S
    dx=IntervalTaylor(c,[0]+list(xi.coefficients[1:])); powers=[IntervalTaylor.constant(c,1,5)]
    for j in range(1,6):powers.append(powers[-1]*dx)
    F=sum((powers[n]*(positive_moment_derivative(c,a,xi[0],n,True)/math.factorial(n)) for n in range(6)),base*0)
    rows=[base*F]
    G={l:sum((powers[n]*(positive_moment_derivative(c,a,xi[0],l+n)/math.factorial(n)) for n in range(6)),base*0)
        for l in range(1,yorder+1)}
    for k in range(1,yorder+1):
        row=base*0
        for l in range(1,k+1):row+=base**l*G[l]*(stirling_second(k,l)*(-1)**(k+1)*S**(l-1))
        rows.append(row)
    return rows


def product_rows(left,right):
    return [sum((left[j]*right[k-j]*math.comb(k,j) for j in range(k+1)),left[0]*0) for k in range(len(left))]


class CompliantCollarGammaC4:
    def __init__(self,cells=64):
        if not isinstance(cells,int) or cells<1:raise ValueError('Positive cell count required')
        self.steep=CompliantSteepWaitingC4(128); self.outer=self.steep.outer; self.ctx=c=self.steep.ctx
        self.mu=self.steep.mu; self.delta=self.steep.delta; self.a=self.delta/2
        self.eps=self.steep.epsilon; self.k=self.steep.k; self.bh=self.steep.bh; self.prate=1+self.delta
        self.family=self.steep.family; self.source=self.steep.source; self.Scap=self.steep.S
        self.S=c.mpf([0,endpoints(self.Scap)[1]]); self.cells=cells; self.hashes=dict(self.steep.hashes)
        self.cache={}; self.shape_cache={}; self.tail_cache={}; self.gamma_cache={}
        name=PREFIX+'compliant_steep_waiting_C4_check.json'; receipt=json.loads((HERE/name).read_bytes())
        if not receipt['all_passed'] or not receipt['angular_steep_and_internal_joins_certified']:
            raise ValueError('Accepted same-source O7 and waiting terminal required')
        for source,digest in receipt['input_hashes'].items():
            if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:raise ValueError('Collar prerequisite changed: '+source)
        self.hashes.update(receipt['input_hashes']); self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        self.theta_base=self.steep.thetaT*c.exp(-self.bh*self.steep.wait-self.steep.logone)
        self.Ev2=self.outer.flatten.Ev2; self.pressure_scale=self.Ev2*self.theta_base**2
        self.hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    def shape(self,Z,t,high=True):
        c=self.ctx; Z=c.mpf(Z); t=c.mpf(t); key=(tuple(endpoints(Z)),tuple(endpoints(t)),high)
        if key in self.shape_cache:return self.shape_cache[key]
        sig=sigma_jets(c,t); phi=phi_jets(c,t); zero=IntervalTaylor.constant(c,0,5); one=zero+1
        sr=[one*(sig[j]*math.factorial(j)) for j in range(5)]
        fr=[one*(phi[j]*math.factorial(j)) for j in range(5)]
        unity=[one]+[zero]*4
        W=[unity[j]-sr[j]+v for j,v in enumerate(product_rows(sr,fr))]
        # The original W is in[0,1]; use its exact range for C0 rather
        # than lose sigmoid/phi correlation in a full collar box.
        W[0]=IntervalTaylor.constant(c,c.mpf([max(mp.mpf(0),endpoints(W[0][0])[0]),min(mp.mpf(1),endpoints(W[0][0])[1])]),5)
        C=[unity[j]-fr[j]*self.eps for j in range(5)]
        pre=[unity[j]-W[j]*self.eps for j in range(5)]
        dr=gamma_deficit_mixed(c,Z,self.a,self.Scap,t,4 if high else 0)
        if high:D=product_rows(product_rows(sr,C),dr)
        else:D=[sr[0]*C[0]*dr[0]]
        K=[pre[j]-D[j]*(self.a*self.S) for j in range(len(D))]
        if endpoints(K[0][0])[0]<=0:raise ArithmeticError('Actual Gamma collar bracket positivity lost')
        out=dict(K_rows=K,D_rows=D,W_rows=W,pre_rows=pre,sigma=sig,phi=phi)
        self.shape_cache[key]=out; return out

    def local_Gamma(self,Z,t):
        c=self.ctx; Z=c.mpf(Z); t=c.mpf(t); key=(tuple(endpoints(Z)),tuple(endpoints(t)))
        if key in self.gamma_cache:return self.gamma_cache[key]
        Sc=self.S*c.exp(-t)
        tails=integrated_gamma_tails(c,Z,self.a,Sc,5,0)
        D=gamma_deficit_mixed(c,Z,self.a,self.Scap,t)
        one=IntervalTaylor.constant(c,1,5); K=[one-D[0]*(self.a*self.S)]+[-v*(self.a*self.S) for v in D[1:]]
        angular=one/self.k+tails['theta']*Sc
        energy=one/self.delta-tails['energy']*(self.a*Sc)
        pressure=one/(2*self.prate)-tails['pressure']*(self.a*Sc)
        out=dict(K_rows=K,angular_numerator=angular,energy_numerator=energy,pressure_numerator=pressure,
            current_inverse_radius_cap=Sc,full_infinite_Gamma_defect_jets=tails)
        self.gamma_cache[key]=out; return out

    def collar_tails(self,Z,t):
        c=self.ctx; Z=c.mpf(Z); t=c.mpf(t); key=(tuple(endpoints(Z)),tuple(endpoints(t)))
        if key in self.tail_cache:return self.tail_cache[key]
        if endpoints(t)[0]<0 or endpoints(t)[1]>3:raise ValueError('Original collar t in[0,3] required')
        # The complete pure Gamma exterior is integrated analytically to
        # infinity in the derivative enclosures, not stopped at a cutoff.
        tails=integrated_gamma_tails(c,Z,self.a,self.Scap,5,3)
        angular=tails['theta']; pressure=tails['pressure']; energy=tails['energy']
        atoms={name:c.mpf(0) for name in ('JW','PW','PW2','EW','EW2')}
        length=3-t
        if endpoints(length)[1]>0:
            ds=length/self.cells
            for i in range(self.cells):
                a=t+length*i/self.cells; b=t+length*(i+1)/self.cells
                v=c.mpf([max(mp.mpf(0),endpoints(a)[0]),min(mp.mpf(3),endpoints(b)[1])])
                shape=self.shape(Z,v,False); D=shape['D_rows'][0]; pre=shape['pre_rows'][0]
                W=shape['W_rows'][0][0]; square=D*(pre*2-D*(self.a*self.S))
                angular+=D*(self.a*ds*c.exp(self.k*v))
                pressure+=square*(ds*c.exp(-self.prate*v)/2)
                energy+=square*(ds*c.exp(-self.delta*v))
                atoms['JW']+=ds*c.exp(self.k*v)*W
                for label,rate in (('P',self.prate),('E',self.delta)):
                    atoms[label+'W']+=ds*c.exp(-rate*v)*W
                    atoms[label+'W2']+=ds*c.exp(-rate*v)*W**2
        one=IntervalTaylor.constant(c,1,5)
        A=one/self.k+(angular*self.S+atoms['JW']*self.eps)*c.exp(-self.k*t)
        E=one*(c.exp(-self.delta*t)/self.delta-2*self.eps*atoms['EW']+self.eps**2*atoms['EW2'])-energy*(self.a*self.S)
        P=one*(c.exp(-self.prate*t)/(2*self.prate)-self.eps*atoms['PW']+self.eps**2*atoms['PW2']/2)-pressure*(self.a*self.S)
        out=dict(angular_numerator=A,remaining_energy_in_Rtail_units=E,
            remaining_pressure_in_Rtail_units=P,scaled_full_future_Gamma_defects=dict(theta=angular,pressure=pressure,energy=energy),
            separate_epsilon_atoms=atoms)
        self.tail_cache[key]=out; return out

    def data(self,Z):
        c=self.ctx; Z=c.mpf(Z); key=tuple(endpoints(Z))
        if key in self.cache:return self.cache[key]
        terminal=self.steep.waiting(Z,1); heat0=self.collar_tails(Z,0)
        Xtail=terminal['angular_Taylor']; Ptail=terminal['pressure_over_Pstar_squared_Taylor']
        defect=Xtail*(1-self.eps)-heat0['angular_numerator']
        pressure3=self.forward_pressure(Z,3,Ptail)
        out=dict(Xtail=Xtail,Ptail=Ptail,angular_tail_constant_defect=defect,pressure3=pressure3,waiting_terminal=terminal)
        self.cache[key]=out; return out

    def forward_pressure(self,Z,t,Ptail):
        c=self.ctx; t=c.mpf(t); integral=Ptail*0
        if endpoints(t)[1]>0:
            ds=t/self.cells
            for i in range(self.cells):
                a=t*i/self.cells; b=t*(i+1)/self.cells
                v=c.mpf([max(mp.mpf(0),endpoints(a)[0]),min(mp.mpf(3),endpoints(b)[1])])
                K=self.shape(Z,v,False)['K_rows'][0]
                integral+=K*K*(ds*c.exp(-self.prate*v)/2)
        return Ptail+integral*self.pressure_scale

    def packet(self,Z,t,shape,X,energy,pressure,kind,extra):
        c=self.ctx; K=shape['K_rows']; rates=quotient_log_rates(K)
        rates[0]+=self.mu-self.a
        theta=K[0]*(self.theta_base*c.exp(-self.bh*t))
        out=self.outer._packet(Z,theta,X,energy,pressure,rates,dict(origin='Rtail',offset=t),
            dict(stage=kind,original_flat_collar_and_exact_Gamma_source_retained=True,
                same_waiting_terminal_angular_pressure_histories_retained=True,
                full_infinite_Gamma_tail_included=True,actual_Gamma_not_defined_by_finite_S_series=True,
                exact_relative_velocity_log_parts=dict(flatten_exit=-100*self.steep.bp-c.ln(2),
                    power_length=-self.steep.bp*self.outer.Lrel,entry_offset=-self.steep.bp-self.steep.rate/2,
                    steep_power=-c.mpf('1.5')*self.steep.Ts,exit_offset=-c.mpf('1.5')+self.k/2,
                    waiting_offset=-self.bh*self.steep.wait,minus_log_one_minus_epsilon=-self.steep.logone,
                    current_power=-self.bh*t,log_actual_Gamma_collar_bracket_Taylor=log_taylor(K[0])),
                **extra))
        out.update(heat_bracket_y_derivatives=K,unbounded_exterior_covered=kind=='exact Gamma exterior')
        return out

    def collar(self,Z,t):
        c=self.ctx; t=c.mpf(t); data=self.data(Z); tails=self.collar_tails(Z,t); shape=self.shape(Z,t)
        K=shape['K_rows'][0]
        X=(tails['angular_numerator']+data['angular_tail_constant_defect']*c.exp(-self.k*t))/K
        energy=tails['remaining_energy_in_Rtail_units']*c.exp(self.delta*t)/(K*K*2)
        pressure=self.forward_pressure(Z,t,data['Ptail'])
        return self.packet(Z,t,shape,X,energy,pressure,'original heat collar',
            dict(original_collar_phi_y_Taylor=shape['phi'],original_sigma_y_Taylor=shape['sigma'],
                scaled_Gamma_deficit_y_derivatives=shape['D_rows'],complete_future_tail=tails,
                retained_angular_tail_constant_defect=data['angular_tail_constant_defect']))

    def exterior(self,Z,t):
        c=self.ctx; t=c.mpf(t)
        if endpoints(t)[0]<3:raise ValueError('Exact Gamma exterior t>=3 required')
        data=self.data(Z); local=self.local_Gamma(Z,t); K=local['K_rows'][0]
        X=(local['angular_numerator']+data['angular_tail_constant_defect']*c.exp(-self.k*t))/K
        energy=local['energy_numerator']/(K*K*2)
        loc3=self.local_Gamma(Z,3)
        pressure=data['pressure3']+(loc3['pressure_numerator']*c.exp(-3*self.prate)
            -local['pressure_numerator']*c.exp(-self.prate*t))*self.pressure_scale
        return self.packet(Z,t,local,X,energy,pressure,'exact Gamma exterior',
            dict(full_local_Gamma_future_source=local,
                retained_angular_tail_constant_defect=data['angular_tail_constant_defect'],
                infinite_offset_handled_by_nonnegative_decay_intervals=True))

    def report(self):
        with mp.workdps(270):
            return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
                selected_mu=self.mu,selected_delta=self.delta,actual_epsilon=self.eps,
                epsilon_relation='epsilon=.001*delta',original_heat_definition='H(xi)=E_Gamma(1+a)[(1+xi*v)^(-a)]',
                exact_Gamma_source_derivative_enclosure='positive expectation derivative and Lipschitz bounds; finite derivative composition, no finite S-series definition',
                collar_samples=[self.collar(z,t) for z in ('-1','0','.5','1') for t in ('0','1','2','3')],
                exterior_samples=[self.exterior(z,t) for z in ('-1','0','.5','1') for t in ('3','4','10','1000')],
                whole_Z_collar_box=self.collar([-1,1],[0,3]),whole_Z_collar_inlet=self.collar([-1,1],0),
                whole_Z_collar_terminal=self.collar([-1,1],3),
                whole_Z_exterior_inlet=self.exterior([-1,1],3),whole_Z_unbounded_exterior=self.exterior([-1,1],[3,mp.inf]),
                phi_endpoint_crossing=self.collar([-1,1],['2.99','3']),
                entire_collar_and_unbounded_Gamma_high_mixed_derivatives_available=True,
                waiting_collar_and_collar_Gamma_joins_certified=False,full_pulse_C4_installed=True,
                full_outer_C4_certified=False,physical_energy_integral_certified=False,
                whole_outer_cone_certified=False,temporal_recursion=False,input_hashes=self.hashes)


def run():
    result=CompliantCollarGammaC4().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Original collar and full unbounded Gamma exterior: source mixed derivatives and retained histories generated',flush=True)
    return result


if __name__=='__main__':run()
