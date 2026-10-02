"""Independent physical/source checks for original microswitch phase4."""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
import sympy as s
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_microswitch_mixed_C4 import (
    CompliantMicroswitchMixedC4,comparison_radial_directions,switch_controls,phase_physical,
    IntervalTaylor,MTH,MTHZ,MZ,MZT,MP,FactoredAlgebra)
from lei_ren_part1_paper_compliant_long_reshape_mixed_C4 import scaled_positive_source,CompliantLongReshapeMixedC4
from lei_ren_part1_paper_compliant_long_reshape_mixed_C4_check import canonical_source
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import sigma_jets
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_compliant_five_moment_repair import pack

HERE=Path(__file__).parent
NAME='lei_ren_part1_paper_compliant_microswitch_mixed_C4.json'


def sigma(t):
    if t<=0:return mp.mpf(0)
    if t>=1:return mp.mpf(1)
    a=mp.exp(-1/t**2);b=mp.exp(-1/(1-t)**2)
    return a/(a+b)


def jet(c,fn,z,tol,order=5):
    return IntervalTaylor(c,[c.mpf([q-tol,q+tol])/math.factorial(n)
        for n in range(order+1) for q in (mp.diff(fn,z,n),)])


def comparison_fixture():
    """Closed physical comparison moments, independent symbolic y/Z diff."""
    with mp.workdps(70):
        c=MPIntervalContext();c.dps=105;tol=mp.mpf('1e-50');y0=mp.mpf('.2');z0=mp.mpf('.3')
        y,z=s.symbols('y z',real=True);delta=s.Rational(1,2000)
        phi=s.exp(z/30+z*z/100);v=4*z+s.Rational(1,50)-z**3/200
        F0=s.exp(z/20+z*z/60);pressure=s.Rational(1,40)+z*z/80+z**6/50000
        initial=[s.Rational(i+1,5)+z/100+z*z/50+(i+1)*z**5/1000 for i in range(6)]
        H=initial[0]*s.exp(-2*y)+phi*(1-s.exp(-2*y))
        m=initial[1]*s.exp(-y)+v*(1-s.exp(-y))
        K=initial[2]*s.exp(-2*y)+phi*v*(1-s.exp(-2*y))
        A=initial[3]*s.exp(-y)+v*v*(1-s.exp(-y))
        B=initial[4]*s.exp(-2*y)+phi*phi*(1-s.exp(-2*y))/2
        P=initial[5]*s.exp(-y)+phi*phi*(1-s.exp(-y))
        point=lambda expr:s.lambdify((y,z),expr,'mpmath')
        zjet=lambda expr:jet(c,lambda zz:point(expr)(y0,zz),z0,tol,6)
        Fbase=point(F0)(y0,z0)
        ratios=lambda power:[c.mpf(mp.diff(lambda zz:point(F0**power)(y0,zz),z0,n)/Fbase**power) for n in range(7)]
        moments={MTH:zjet(H),MZ:zjet(m),MTHZ:zjet(K),MZT:dict(axial=zjet(A),swirl=zjet(B)),MP:zjet(P)}
        bounded=comparison_radial_directions(c,c.mpf(z0),c.mpf(delta.p)/int(delta.q),zjet(phi),zjet(v),moments,zjet(pressure),ratios(1),ratios(2))
        d=1-z*z;L=1-delta*z*z;W=1-(1-delta)*z*m-d*s.diff(m,z)
        f=F0*phi;h=F0*H;k=F0*K;b=F0*F0*B;p=F0*F0*P
        angular=(1-delta/2)*h-(1-delta)*z*s.diff(h,z)/2-d*s.diff(k,z)+(2*delta-1)*z*k
        formulas=dict(D_over_R=(-W+angular/(2*f))/L,
            hydro=(-W*v+(1-delta)*(m-z*s.diff(m,z))/2+2*delta*z*A-d*s.diff(A,z))/(2*L),
            pressure=(2*(1+delta)*z*pressure-d*s.diff(pressure,z))/(2*L),
            swirl=(-2*delta*z*b+d*s.diff(b,z)+2*(1+delta)*z*p-d*s.diff(p,z))/(2*L))
        count=0
        for name,rows in bounded.items():
            for k,row in enumerate(rows):
                if row.order!=5:raise ValueError('Comparison direction did not consume exactly one axial order')
                for n in range(6):
                    expected=point(s.diff(formulas[name],y,k,z,n))(y0,z0)/(Fbase**2 if name=='swirl' else 1)
                    lo,hi=endpoints(row[n]*math.factorial(n))
                    if not lo-tol*10000<=expected<=hi+tol*10000:raise ArithmeticError('Independent comparison direction failed: '+name)
                    count+=1
        return dict(independent_comparison_ordinary_y_Z_derivatives=count,true_F0_and_F0_squared_derivatives=True,
            axial6_inputs_direction_axial5_outputs=True,finite_fixture_only=True,actual_source_admission=False,passed=True)


def controls_fixture():
    with mp.workdps(70):
        c=MPIntervalContext();c.dps=105;tol=mp.mpf('1e-50');h=mp.mpf('.025');z=mp.mpf('.3');Pstar=mp.mpf('2.4')
        D=lambda zz:mp.mpf('.003')+zz/1000+zz**3/2000
        H=lambda zz:mp.mpf('.02')+zz/1000
        P=lambda zz:mp.mpf('.0002')+zz**2/50000
        S=lambda zz:mp.mpf('.00003')+zz**3/100000
        F0=lambda zz:mp.exp(zz/20+zz*zz/60)
        phi100=lambda zz:1+zz/50+zz**3/500
        barphi=lambda zz:1+zz/100+zz*zz/200
        radius=lambda t:100*mp.exp(h*t)
        def a(t,zz,branch):
            base=h*radius(t)*D(zz)
            return base if branch=='first' else base*(1-sigma(t-1))+mp.mpf('.8')*sigma(t-1)
        first_phi=lambda t,zz:phi100(zz)*mp.exp(-h*D(zz)*(radius(t)-100)/2)
        count=0
        for branch,phase in (('first',mp.mpf('.35')),('second',mp.mpf('1.65'))):
            R=radius(phase);Fbase=F0(z)
            if branch=='first':current_phi=lambda zz:first_phi(phase,zz)
            else:
                J=mp.quad(lambda t:radius(t)*(1-sigma(t-1)),[1,phase]);mass=mp.quad(lambda t:sigma(t),[0,phase-1])
                current_phi=lambda zz:first_phi(1,zz)*mp.exp(-h*h*D(zz)*J/2-mp.mpf('.4')*h*mass)
            quotient=jet(c,lambda zz:current_phi(zz)/barphi(zz),z,tol)
            Drows=[jet(c,lambda zz:R*D(zz),z,tol)]*4
            algebra=FactoredAlgebra(c,[c.mpf(mp.log(h)),2*c.mpf(mp.log(Pstar)),2*c.mpf(mp.log(Fbase)),c.mpf(0)],[])
            scale=algebra.width
            def driver(k,power):
                return scale(algebra.lift(jet(c,lambda zz:R*H(zz),z,tol))
                    +algebra.shift(jet(c,lambda zz:R*P(zz),z,tol),(0,1,0,0))
                    +algebra.shift(jet(c,lambda zz:2**k*R*R*F0(zz)**2*S(zz)/Fbase**2,z,tol),(0,0,1,0)),power)
            control=algebra.tree(switch_controls(c,branch,sigma_jets(c,c.mpf(phase-(0 if branch=='first' else 1))),Drows,driver,algebra.lift(quotient),scale))
            rhs=lambda t,zz:-h*h*(1-sigma(t))*(first_phi(t,zz)/barphi(zz))*(radius(t)*H(zz)+radius(t)*Pstar**2*P(zz)+radius(t)**2*F0(zz)**2*S(zz))
            formulas={'a_phase_derivatives':lambda t,zz:a(t,zz,branch),
                'logF_phase_derivatives':lambda t,zz:-h*a(t,zz,branch)/2,
                'logUtheta_phase_derivatives':lambda t,zz:h*(1-a(t,zz,branch))/2,
                'Uz_positive_phase_derivatives':rhs if branch=='first' else lambda t,zz:mp.mpf(0)}
            for name,fn in formulas.items():
                for k,row in enumerate(control[name]):
                    for n in range(6):
                        expected=mp.diff(fn,(phase,z),(k,n));lo,hi=endpoints(row[n]*math.factorial(n))
                        if not lo-tol*10000<=expected<=hi+tol*10000:raise ArithmeticError('Original switch control Leibniz derivative failed: '+name)
                        count+=1
        return dict(independent_original_switch_control_derivatives=count,both_branches_and_all_drive_scales=True,
            finite_fixture_only=True,actual_source_admission=False,passed=True)


def physical_fixture():
    """Full original physical primitive RHS integrals for varying phase V/u."""
    with mp.workdps(65):
        c=MPIntervalContext();c.dps=95;tol=mp.mpf('1e-45');h=mp.mpf('.025');t=mp.mpf('.3');z=mp.mpf('.2');Pstar=mp.mpf('2.4');delta=mp.mpf('.0005')
        R=lambda q:100*mp.exp(h*q)
        ell=lambda q,zz:mp.mpf('-.8')-mp.log(1+zz*zz)+(mp.mpf('.03')+zz/100+zz**3/1000)*(q+q*q/2-q**3/3+q**4/24)
        u=lambda q,zz:Pstar*mp.exp(ell(q,zz))
        pol=[mp.mpf(v) for v in (1,-1,1)]+[-mp.mpf(1)/3,mp.mpf(1)/6]
        V=lambda q,zz:4*zz+(mp.mpf('.01')+zz*zz/500)*sum(a*q**k for k,a in enumerate(pol))
        initial=lambda i,zz:mp.mpf(i+1)/20+zz/1000+zz**2/500+zz**5/10000
        p0=lambda zz:mp.mpf('.03')+zz/100+zz**4/1000
        rhs=dict(Mtheta=lambda q,zz:h*mp.sqrt(2)*R(q)**mp.mpf('1.5')*u(q,zz),
            Mtheta_z=lambda q,zz:h*mp.sqrt(2)*R(q)**mp.mpf('1.5')*u(q,zz)*V(q,zz),
            Mz=lambda q,zz:h*R(q)*V(q,zz),axial=lambda q,zz:h*R(q)*V(q,zz)**2,
            swirl=lambda q,zz:h*R(q)*u(q,zz)**2,Mp=lambda q,zz:h*u(q,zz)**2/2)
        values={name:[mp.diff(lambda zz:initial(i,zz),z,n)/math.factorial(n)+mp.quad(
            lambda q:mp.diff(lambda zz:fn(q,zz),z,n),[0,t])/math.factorial(n) for n in range(6)]
            for i,(name,fn) in enumerate(rhs.items())}
        cover=lambda row:IntervalTaylor(c,[c.mpf([q-tol,q+tol]) for q in row])
        data={name:cover(row) for name,row in values.items()};uu=u(t,z);rr=R(t)
        uj=jet(c,lambda zz:u(t,zz),z,tol);amp=uj/c.mpf(uu)
        shapes=dict(theta=data['Mtheta']/(uj*(c.sqrt(2)*c.mpf(rr)**c.mpf('1.5'))),
            theta_z=data['Mtheta_z']/(uj*(c.sqrt(2)*c.mpf(rr)**c.mpf('1.5'))),mean=data['Mz']/c.mpf(rr),
            axial=data['axial']/c.mpf(rr),swirl=data['swirl']/(uj*uj*c.mpf(rr)),pressure=(data['Mp']*2)/(uj*uj))
        algebra=FactoredAlgebra(c,[c.mpf(mp.log(h)),2*c.mpf(mp.log(Pstar)),c.mpf(0),2*c.mpf(mp.log(uu/Pstar))],[])
        width=algebra.width
        logU=[jet(c,lambda zz,k=k:mp.diff(lambda q:ell(q,zz),t,k),z,tol) for k in range(1,5)]
        Vs=[jet(c,lambda zz,k=k:mp.diff(lambda q:V(q,zz),t,k),z,tol) for k in range(5)]
        packet=phase_physical(c,c.mpf(z),c.mpf(delta),algebra.width(1),amp,c.mpf(mp.log(uu/Pstar)),[algebra.lift(row) for row in logU],Vs,shapes,jet(c,p0,z,tol),algebra.shift(1,(0,-1,0,0)),width,[])
        def integral_exp_poly(q):
            return sum(a*(mp.exp(h*q)*sum((-1)**j*mp.factorial(i)/mp.factorial(i-j)*q**(i-j)/h**(j+1)
                for j in range(i+1))-(-1)**i*mp.factorial(i)/h**(i+1)) for i,a in enumerate(pol))
        def mean(q,zz):
            mass=initial(2,zz)+100*(4*zz*(mp.exp(h*q)-1)+h*(mp.mpf('.01')+zz*zz/500)*integral_exp_poly(q))
            return mass/R(q)
        def ur(q,zz):
            mm=mean(q,zz);mz=mp.diff(lambda a:mean(q,a),zz)
            return mp.sqrt(R(q)/2)*(2*zz*V(q,zz)-(1-delta)*zz*mm-(1-zz*zz)*mz)/(1-delta*zz*zz)
        scales={'Utheta_over_current_Utheta':uu,'Uz':1,'Ur_over_current_sqrt_R_over_2':mp.sqrt(rr/2),'P_over_Pstar2':1,
            'Mtheta_over_current_sqrt2_R_1p5_Utheta':mp.sqrt(2)*rr**mp.mpf('1.5')*uu,
            'Mtheta_z_over_current_sqrt2_R_1p5_Utheta':mp.sqrt(2)*rr**mp.mpf('1.5')*uu,'Mz_over_current_R':rr,
            'Mztheta_over_current_R_Pstar2':rr*Pstar**2,'Mp_over_Pstar2':Pstar**2}
        names={'Mtheta_over_current_sqrt2_R_1p5_Utheta':'Mtheta','Mtheta_z_over_current_sqrt2_R_1p5_Utheta':'Mtheta_z',
            'Mz_over_current_R':'Mz','Mp_over_Pstar2':'Mp'}
        count=0
        for group in ('physical_velocity_pressure_phase_Z_mixed4','physical_five_primitive_phase_Z_mixed4'):
            for name,grid in packet[group].items():
                for key,value in grid.items():
                    k,n=[int(q[1:]) for q in key.split('_')]
                    if name in ('Utheta_over_current_Utheta','Uz','Ur_over_current_sqrt_R_over_2'):
                        fn={'Utheta_over_current_Utheta':u,'Uz':V,'Ur_over_current_sqrt_R_over_2':ur}[name];expected=mp.diff(fn,(t,z),(k,n))
                    elif k==0:
                        if name=='P_over_Pstar2':expected=mp.diff(p0,z,n)+values['Mp'][n]*math.factorial(n)/Pstar**2
                        elif name.startswith('Mztheta'):expected=(values['axial'][n]-values['swirl'][n]/2)*math.factorial(n)
                        else:expected=values[names[name]][n]*math.factorial(n)
                    else:
                        if name=='P_over_Pstar2':fn=lambda q,zz:rhs['Mp'](q,zz)/Pstar**2
                        elif name.startswith('Mztheta'):fn=lambda q,zz:rhs['axial'](q,zz)-rhs['swirl'](q,zz)/2
                        else:fn=rhs[names[name]]
                        expected=mp.diff(fn,(t,z),(k-1,n))
                    expected/=scales[name];lo,hi=endpoints(value)
                    if not lo-tol*10000<=expected<=hi+tol*10000:raise ArithmeticError('Independent physical phase derivative failed: '+name+' '+key)
                    count+=1
        return dict(independent_full_physical_phase_Z_derivatives=count,varying_axial_velocity_and_nonzero_histories=True,
            true_radial_prefactors_and_original_pressure=True,finite_fixture_only=True,actual_source_admission=False,passed=True)


def factored_scale_fixture():
    """Extreme factors whose combined source is moderate: premature cap fails."""
    with mp.workdps(70):
        c=MPIntervalContext();c.dps=105
        algebra=FactoredAlgebra(c,[c.mpf(-3000),c.mpf(6001),c.mpf(-4000),c.mpf(5999)],[])
        z=IntervalTaylor.variable(c,c.mpf('.3'),5);q=z*z+2
        sources=[algebra.width(algebra.shift(q,(0,1,0,0)),2),
            algebra.shift(algebra.width(q,2),(0,0,0,1)),
            algebra.shift(algebra.shift(q,(0,1,0,0)),(0,0,1,0))]
        # Verify against exp(1), exp(-1), exp(2001), without tiny exponentials.
        for source,logbase in zip(sources,(1,-1,2001)):
            row=algebra.resolve(source)
            for n in range(6):
                expected=c.exp(c.mpf(logbase))*q[n]
                if endpoints(row[n])!=endpoints(expected):raise ArithmeticError('Source factors capped before the final combined log')
        return dict(extreme_source_products_checked=3,ordinary_axial_coefficients_checked=18,
            catches_premature_width_or_amplitude_caps=True,finite_fixture_only=True,actual_source_admission=False,passed=True)


def functional_identities():
    h=s.symbols('h',positive=True);D=s.symbols('D0:4');q=s.symbols('q0:5');v,m=s.symbols('v m')
    sig=s.symbols('sig0:4');first=[h**(k+1)*D[k] for k in range(4)]
    second=[sum(s.binomial(k,j)*first[j]*((1-sig[0]) if k-j==0 else -sig[k-j]) for j in range(k+1))
        +s.Rational(4,5)*sig[k] for k in range(4)]
    left={sig[k]:0 for k in range(4)};right={sig[k]:(1 if k==0 else 0) for k in range(4)}
    checks=0
    def zero(expr):
        nonlocal checks
        if s.simplify(expr)!=0:raise ArithmeticError('Original switch functional source identity failed')
        checks+=1
    for k in range(4):zero(second[k].subs(left)-first[k]);zero(second[k].subs(right)-(s.Rational(4,5) if k==0 else 0))
    # V derivatives vanish at first-switch terminal since every complement
    # derivative vanishes. At phase1, second angular control equals first.
    for k in range(4):zero(sum(s.binomial(k,j)*((1-sig[0]) if j==0 else -sig[j])*q[k-j] for j in range(k+1)).subs(right))
    # At R2, logU_s=h/10 and its higher derivatives vanish. Physical
    # primitives retain their fixed basepoint normalization on both sides.
    rates={'theta':s.Rational(8,5),'theta_z':s.Rational(8,5),'swirl':s.Rational(6,5),'pressure':s.Rational(1,5)}
    for rate in rates.values():
        for k in range(1,5):zero(h*(h*rate)**(k-1)-h**k*rate**(k-1))
    for k in range(5):
        zero((h/s.Integer(10))**k-h**k/s.Integer(10)**k)
        zero(sum(s.binomial(k,j)*h**j*q[j]*(h/2)**(k-j) for j in range(k+1))
            -h**k*sum(s.binomial(k,j)*q[j]/2**(k-j) for j in range(k+1)))
    for k in range(1,5):
        zero(h*(h**(k-1))*v-h**k*v)
        zero(h*(h**(k-1))*v*v-h**k*v*v)
        zero(h**k*(-1)**(k-1)*(v-m)-h**k*((-1)**(k-1)*(v-m)))
    y,J=s.symbols('y JD');logF2=-h/5-h*h*J/2
    zero(logF2-s.Rational(2,5)*(y-2*h)-(-s.Rational(2,5)*y+s.Rational(3,5)*h-h*h*J/2))
    f=s.symbols('fraction');zero((2*h+f*(y-2*h)).subs(f,1)-y)
    # IBP weight sum and positive radial kernel masses, exact for all Z.
    t=s.symbols('t',real=True);theta=s.symbols('theta',positive=True)
    for rate,mass in ((2,1-theta**2),(1,1-theta),(2,(1-theta**2)/2)):
        multiplier=2 if rate==2 and mass==1-theta**2 else 1
        zero(s.integrate(multiplier*s.exp(-rate*(y-t)),(t,0,y)).subs(s.exp(-y),theta)-mass)
    alpha=s.Function('alpha');g=s.Function('g')
    averaged=alpha(y)*g(y)+s.Integral(-s.diff(alpha(t),t)*g(t),(t,0,y))
    zero(s.diff(averaged,y)-alpha(y)*s.diff(g(y),y))
    zero(s.diff(alpha(y)+s.Integral(-s.diff(alpha(t),t),(t,0,y)),y))
    a,b,c=(s.Function(name) for name in ('a','b','c'))
    integrand=a(y)*b(y)*c(y)
    zero(s.diff(s.Integral(integrand,(y,0,t)),t)-integrand.subs(y,t))
    for k in range(4):
        rule=sum(s.factorial(k)/s.factorial(i)/s.factorial(j)/s.factorial(k-i-j)
            *s.diff(a(y),y,i)*s.diff(b(y),y,j)*s.diff(c(y),y,k-i-j)
            for i in range(k+1) for j in range(k-i+1))
        zero(s.diff(integrand,y,k)-rule)
    return dict(symbolic_functional_source_identities=checks,
        phase1_join='same actual histories/P0; original sigma endpoint jets; both V positive-order derivatives zero',
        R2_join='same actual R2 histories/P0; all phase rows equal hb^k times original postpower rows',
        R110_join='same actual R2 transport; exact theta=R2/R; angular identity with +.6hb; accepted original post(110) source',
        includes_IBP_and_positive_kernel_masses=True,interval_overlap_is_not_the_join_proof=True,passed=True)


def run():
    with mp.workdps(280):
        raw=json.loads((HERE/NAME).read_bytes());provider=CompliantMicroswitchMixedC4();c=provider.ctx
        for path,digest in raw['input_hashes'].items():
            if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:raise ValueError('Microswitch source changed: '+path)
        packets=[raw[name] for name in ('whole_first','whole_second','actual_R100_inlet','first_side_phase1','second_side_phase1','actual_R2_exit')]+raw['interior_packets']
        counts={name:0 for name in ('physical_velocity_pressure_phase_Z_mixed4','physical_five_primitive_phase_Z_mixed4')};ledgers=0
        for packet in packets:
            if not packet['source_width_and_amplitude_products_capped_only_after_final_derivatives']:raise ValueError('Premature phase source capping')
            bases=[read_interval(c,row) for row in packet['factored_source_log_bases']]
            width=read_interval(c,packet['width_enclosure_is_not_source'])
            if endpoints(width)[0]!=0 or endpoints(bases[0])[1]>endpoints(c.ln(c.mpf(endpoints(width)[1])))[0]:raise ValueError('Formal positive hb not enclosed by the numerical width interval')
            rows={row['physical_row']:row for row in packet['final_factored_physical_row_ledgers']}
            for group in counts:
                for name,grid in packet[group].items():
                    if set(grid)!={'s'+str(k)+'_Z'+str(n) for k in range(5) for n in range(5-k)}:raise ValueError('Incomplete phase mixed grid')
                    for key,value in grid.items():
                        lo,hi=endpoints(read_interval(c,value))
                        if not lo<=hi or not mp.isfinite(lo) or not mp.isfinite(hi):raise ArithmeticError('Nonfinite phase derivative')
                        countrow=packet['physical_logR_Z_mixed4_log_bound_ledger'][group][name][key.replace('s','y',1)]
                        k=int(key.split('_')[0][1:]);size=max(abs(lo),abs(hi))
                        if countrow['source_phase_order']!=k:raise ValueError('Incorrect formal logR derivative ledger')
                        if size:
                            expected=c.ln(c.mpf(size))-k*provider.logh;stored=read_interval(c,countrow['log_absolute_upper_from_capped_phase'])
                            if endpoints(stored)!=endpoints(expected):raise ArithmeticError('Formal hb^-k log derivative conversion changed')
                        counts[group]+=1;ledgers+=1
                        final=rows[name+'/'+key];reconstructed=c.mpf(0);y_logs=[]
                        for term in final['terms']:
                            exponents=term['source_exponents']
                            combined=sum((bases[i]*p for i,p in enumerate(exponents) if p),c.mpf(0))
                            if endpoints(combined)!=endpoints(read_interval(c,term['combined_positive_source_log'])):raise ArithmeticError('Source factors lost in final derivative log')
                            if not term['all_Bell_Leibniz_and_axial_factors_already_combined']:raise ValueError('Intermediate cap claimed as final derivative')
                            coefficient=read_interval(c,term['final_ordinary_coefficient'])
                            reconstructed+=scaled_positive_source(c,combined,IntervalTaylor.constant(c,coefficient),[])[0]
                            magnitude=max(abs(v) for v in endpoints(coefficient))
                            if magnitude:
                                powers=list(exponents);powers[0]-=k
                                y_logs.append(sum((bases[i]*p for i,p in enumerate(powers) if p),c.mpf(0))+c.ln(c.mpf(magnitude)))
                        if endpoints(reconstructed)!=(lo,hi):raise ArithmeticError('Final physical row differs from its factored source sum')
                        if countrow['exact_zero']!=(not y_logs) or countrow['nonzero_factored_terms']!=len(y_logs) or not countrow['formal_width_powers_cancelled_before_log_bound']:raise ValueError('Premature phase cap polluted logR derivative bound')
                        if y_logs:
                            upper=c.mpf(max(endpoints(row)[1] for row in y_logs))+c.ln(len(y_logs))
                            if endpoints(read_interval(c,countrow['log_absolute_upper']))!=endpoints(upper):raise ArithmeticError('Factored logR triangle bound changed')
            if not packet['inverse_hb_not_materialized'] or not packet['phase_derivatives_do_not_differentiate_width_caps']:raise ValueError('Exact width replaced by a cap')
            if any(packet[key] for key in ('bridge_switch_inlet_join_certified','full_inner_interfaces_certified','full_cartesian_vector_derivatives_certified','admissible_stress_lift_constructed','temporal_recursion')):raise ValueError('Unbuilt spatial/temporal scope promoted')
        for group in counts:
            if canonical_source(raw['first_side_phase1'][group])!=canonical_source(raw['second_side_phase1'][group]):raise ArithmeticError('Exact phase1 two-sided physical derivatives differ')
        capcount=0
        for proof in raw['factored_width_and_amplitude_cap_proofs']:
            if endpoints(read_interval(c,proof['log_magnitude_upper']))[1]>endpoints(read_interval(c,proof['log_cap']))[0] or not proof['exact_source_not_replaced']:raise ArithmeticError('Factored actual width/source cap invalid')
            capcount+=1
        distance=read_interval(c,raw['comparison_smoothing_log_distance_to_R100'])
        if endpoints(distance)[0]<=endpoints(2*provider.switch.cap)[1]:raise ValueError('Original comparison smoothing reaches switching region')
        if endpoints(read_interval(c,raw['comparison_smoothing_continuation_log_distance']))[0]<=endpoints(2*provider.switch.cap)[1]:raise ValueError('Comparison smoothing leaves its covered core continuation')
        enclosure=raw['comparison_moment_enclosure_source']
        if enclosure['axial_Taylor_orders']!=list(range(7)) or enclosure['exact_moment_point_values_reconstructed']:raise ValueError('Comparison moment bounds promoted to exact point values')
        shared=json.loads((HERE/'lei_ren_part1_paper_compliant_long_reshape_mixed_C4.json').read_bytes())['shared_exact_axial_source']
        if canonical_source(shared)!=canonical_source(raw['shared_exact_axial_source']):raise ValueError('Microswitch axial source graph differs from accepted reshape source')
        for packet in packets:
            expected=(dict(op='sum',args=shared['V100']['args']+[
                dict(formal_integral='I_first_switch',shared_source=shared['shared_source_namespace'],bounds=['0','s'],
                    source_integrand=shared['formal_signed_integrals']['I_first_switch'])]) if packet['branch']=='first' else shared['V110'])
            if canonical_source(packet['actual_axial_function_source'])!=canonical_source(expected):raise ValueError('Actual partial first-switch integral source changed')
            binding=packet['first_switch_derivative_source']
            if binding['integrand']!='-hb^2*(1-sigma(s))*(phi_actual/barphi)*(R*hydro+R*Pstar^2*pressure+R^2*F0^2*swirl)':raise ValueError('Actual first-switch derivative sign/scale changed')
        R2=raw['actual_R2_exit']['actual_parent_axial5_packet'];power0=raw['actual_R2_power_inlet']
        if canonical_source(R2['actual_moment_shape_axial5_coefficients'])!=canonical_source(power0['actual_postswitch_moment_shapes']):raise ArithmeticError('Actual R2 moments reset at power inlet')
        for name,right in (('F_actual_over_F0_axial5_coefficients','actual_postswitch_phi_axial5'),('Uz_actual_axial5_coefficients','actual_postswitch_V_axial5'),('pressure_axis_axial5_coefficients','original_P0_axial5')):
            if canonical_source(R2[name])!=canonical_source(power0[right]):raise ArithmeticError('Actual R2 source/P0 changed')
        postcounts={group:0 for group in ('physical_velocity_pressure_y_Z_mixed4','physical_five_primitive_y_Z_mixed4')}
        for name in ('whole_postswitch_power','actual_R2_power_inlet','actual_R110_power_exit'):
            packet=raw[name]
            if packet['source_positive_R2']!='100*exp(2*hb)' or packet['source_theta']!='exp(-zeta)=R2/R; no rounded R2':raise ValueError('Original formal R2/transport source changed')
            if canonical_source(packet['actual_axial_function_source'])!=canonical_source(shared['V110']):raise ValueError('Postpower axial source changed')
            width=read_interval(c,packet['width_enclosure_is_not_source']);logh=read_interval(c,packet['exact_positive_width_log'])
            if endpoints(width)[0]!=0 or endpoints(logh)[1]>endpoints(c.ln(c.mpf(endpoints(width)[1])))[0]:raise ValueError('Formal R2 width outside its admitted interval')
            fraction=read_interval(c,packet['source_fraction']);length=c.ln(c.mpf(110)/100)-2*width
            zeta=length*fraction;theta=c.exp(-zeta)
            if endpoints(read_interval(c,packet['zeta_enclosure_only']))!=endpoints(zeta) or endpoints(read_interval(c,packet['theta_enclosure_only']))!=endpoints(theta):raise ValueError('Formal transport log/ratio enclosure changed')
            # Source radius is the convex log combination 2hb*(1-f)+f*log1.1.
            # hb>0 and 2hb<log1.1 establish the [100,110] clipping analytically.
            if endpoints(length)[0]<=0 or not 0<=endpoints(fraction)[0]<=endpoints(fraction)[1]<=1:raise ValueError('Formal postpower radius ordering not established')
            originalR=c.mpf(100)*c.exp(2*width+zeta)
            if endpoints(fraction)==(1,1):expectedR=c.mpf(110)
            else:expectedR=c.mpf([max(mp.mpf(100),endpoints(originalR)[0]),min(mp.mpf(110),endpoints(originalR)[1])])
            if endpoints(read_interval(c,packet['R_enclosure_only']))!=endpoints(expectedR):raise ValueError('Formal source radius not contained in numerical enclosure')
            if packet['formal_log_radius_tree']!=dict(op='sum',args=['log(100)','2hb',dict(op='mul',args=['fraction','log(110/100)-2hb'])]):raise ValueError('Formal positive R2 source tree changed')
            for group in postcounts:
                for grid in packet[group].values():
                    if set(grid)!={'y'+str(k)+'_Z'+str(n) for k in range(5) for n in range(5-k)}:raise ValueError('Incomplete full R2..110 power mixed4')
                    for value in grid.values():
                        lo,hi=endpoints(read_interval(c,value))
                        if not lo<=hi or not mp.isfinite(lo) or not mp.isfinite(hi):raise ArithmeticError('Nonfinite original full postpower derivative')
                        postcounts[group]+=1
        original=encode(pack(provider.switch.post([-1,1],110)))
        end=raw['actual_R110_power_exit']
        for originalname,newname in (('Uz_actual_axial5_coefficients','actual_postswitch_V_axial5'),('pressure_axis_axial5_coefficients','original_P0_axial5')):
            if canonical_source(original[originalname])!=canonical_source(end[newname]):raise ValueError('Original R110 V/P0 changed')
        if original['exact_post_log_identity']!=end['source_angular_identity'].replace('log(F/F100)=-.4','log(F/F100)=-(2/5)'):raise ValueError('Original R110 angular identity changed')
        if endpoints(read_interval(c,end['R_enclosure_only']))!=(110,110):raise ValueError('Exact R110 endpoint offset not cancelled')
        for endpoint in (0,1):
            cutoff=sigma_jets(c,c.mpf(endpoint))
            if endpoints(cutoff[0])!=(endpoint,endpoint) or any(endpoints(q)!=(0,0) for q in cutoff[1:]):raise ValueError('Original flat cutoff endpoint jets changed')
        result=dict(actual_five_defect_family_sha256=provider.family,implicit_source_sha256=provider.source,
            comparison_fixture=comparison_fixture(),controls_fixture=controls_fixture(),physical_fixture=physical_fixture(),
            factored_scale_fixture=factored_scale_fixture(),functional_identities=functional_identities(),
            actual_phase_mixed_bounds_checked=counts,exact_formal_logR_derivative_ledgers_checked=ledgers,
            actual_complete_postpower_mixed_bounds_checked=postcounts,
            final_factored_physical_row_source_sums_checked=ledgers,comparison_own_history_axial6_enclosure_proof=True,
            actual_partial_first_switch_integrand_bound_to_shared_source=True,formal_hb_R2_zeta_theta_and_R110_containments_checked=True,
            source_width_and_amplitude_caps_deferred_until_final_physical_rows=True,
            actual_factored_positive_width_amplitude_caps_checked=capcount,exact_phase1_two_sided_physical_grids=True,
            actual_R2_moments_velocity_pressure_retained=True,original_comparison_smoothing_finished_before_R100=True,
            switch_phase_mixed4_available=True,switch_logR_mixed4_log_ledger_available=True,
            complete_postswitch_power_installed=True,bridge_switch_inlet_join_certified=False,
            phase1_R2_and_R110_functional_mixed4_joins_certified=True,
            full_inner_interfaces_certified=False,full_cartesian_vector_derivatives_certified=False,
            admissible_stress_lift_constructed=False,temporal_recursion=False,all_passed=True,
            input_hashes={**raw['input_hashes'],NAME:hashlib.sha256((HERE/NAME).read_bytes()).hexdigest(),Path(__file__).name:hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('Original microswitch phase4, exact hb logR scales and actual source histories PASS',flush=True)
    return result


if __name__=='__main__':run()
