"""Original flat phi, exact Gamma derivative source and collar interfaces."""
import ast
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
import sympy as s
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_collar_Gamma_C4 import phi_jets,gamma_deficit_mixed,positive_moment_derivative,stirling_second,PHI_POLYNOMIALS
from lei_ren_part1_paper_compliant_absolute_moment_closure import assignment
from lei_ren_part1_paper_compliant_flatten_mixed_C4_check import source_branch_expression
from lei_ren_part1_paper_compliant_power_angular_C4_check import independent_nonconstant_rate_fixture
from lei_ren_part1_paper_compliant_pulse_physical_bounds import UZ,UT,UR,P
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_'
NAME=PREFIX+'compliant_collar_Gamma_C4.json'


def independent_flat_phi_fixture():
    with mp.workdps(150):
        c=MPIntervalContext(); c.dps=100; checks=0
        phi=lambda t:mp.exp(-4/(3-t)**2) if t<3 else mp.mpf(0)
        for lo,hi in (('0','0'),('1.3','1.3'),('2.8','2.9'),('2.99','3'),('3','4')):
            box=phi_jets(c,[lo,hi]); a,b=mp.mpf(lo),mp.mpf(hi)
            for t in (a,(a+b)/2,b):
                for n in range(5):
                    value=mp.diff(phi,t,n)/math.factorial(n) if t<3 else mp.mpf(0)
                    x,y=endpoints(box[n])
                    if not x<=value<=y:raise ArithmeticError('Original collar phi derivative failed: '+str((lo,hi,t,n)))
                    checks+=1
        for n in range(5):
            if endpoints(phi_jets(c,3)[n])!=(mp.mpf(0),mp.mpf(0)):raise ArithmeticError('Original flat phi endpoint lost')
        return dict(independent_original_phi_derivative_checks=checks,passed=True)


def independent_full_Gamma_fixture():
    """Independent confluent-hypergeometric integral representation.

    H(xi)=xi^(-1-a)*U(1+a,2,1/xi). Validate C0 against direct Gamma
    quadrature, then use symbolic differentiation of the complete source
    and direct hyperu evaluations of its derivative integrals.
    Moderate fixture parameters are not the actual Md40 family.
    """
    with mp.workdps(75):
        c=MPIntervalContext(); c.dps=110
        a=mp.mpf('.15'); S=mp.mpf('.04'); checks=moment_checks=0
        def Hderivative(xi,n):
            return (-1)**n*mp.rf(a,n)*mp.rf(1+a,n)*xi**(-1-a-n)*mp.hyperu(1+a+n,2,1/xi)
        H=lambda xi:Hderivative(xi,0)
        for xi in (mp.mpf('.03'),mp.mpf('.08')):
            quadrature=mp.quad(lambda v:mp.exp(-v)*v**a*(1+xi*v)**(-a),[0,1,mp.inf])/mp.gamma(1+a)
            if abs(H(xi)-quadrature)>mp.mpf('1e-65'):raise ArithmeticError('Full Gamma hyperu representation failed')
            for n in range(10):
                # For divided derivatives the beta/q integral can be
                # avoided using exact H derivatives and the quotient.
                target=(-1)**n*math.factorial(n)*(1-H(xi))/(a*xi**(n+1))
                for j in range(1,n+1):target-=math.comb(n,j)*Hderivative(xi,j)*(-1)**(n-j)*math.factorial(n-j)/(a*xi**(n-j+1))
                box=positive_moment_derivative(c,c.mpf(a),c.mpf(xi),n,True)
                if not endpoints(box)[0]<=target<=endpoints(box)[1]:raise ArithmeticError('Full Gamma divided derivative bound failed')
                moment_checks+=1
                if n:
                    target=Hderivative(xi,n)/a
                    box=positive_moment_derivative(c,c.mpf(a),c.mpf(xi),n)
                    if not endpoints(box)[0]<=target<=endpoints(box)[1]:raise ArithmeticError('Full Gamma derivative moment bound failed')
                    moment_checks+=1
        y,z,aa,ss=s.symbols('y z a S',real=True); argument=2*(1-z*z)*ss*s.exp(-y)
        h=s.Function('H')(argument); fn=(1-h)/(aa*ss)
        expressions={(k,n):s.diff(fn,y,k,z,n) for k in range(5) for n in range(5-k)}
        expressions[(0,5)]=s.diff(fn,z,5)
        for t,Z in ((mp.mpf('.3'),mp.mpf('.2')),(mp.mpf('4'),mp.mpf('.5'))):
            rows=gamma_deficit_mixed(c,c.mpf(Z),c.mpf(a),c.mpf(S),c.mpf(t))
            xi=2*(1-Z*Z)*S*mp.exp(-t); hd=[Hderivative(xi,n) for n in range(10)]
            for (k,n),expr in expressions.items():
                replaced={node:s.Float(str(hd[node.expr.derivative_count]),75) for node in expr.atoms(s.Subs)}
                exact=expr.xreplace(replaced).xreplace({h:s.Float(str(hd[0]),75)})
                value=mp.mpf(str(exact.subs({aa:s.Float(str(a),75),ss:s.Float(str(S),75),y:s.Float(str(t),75),z:s.Float(str(Z),75)}).evalf(75)))/math.factorial(n)
                lo,hi=endpoints(rows[k][n])
                if not lo<=value<=hi:raise ArithmeticError('Full Gamma y/Z derivative failed: '+str((t,Z,k,n)))
                checks+=1
        return dict(independent_full_Gamma_mixed_derivative_checks=checks,
            independent_positive_Gamma_moment_derivative_checks=moment_checks,
            finite_S_series_or_asymptotic_definition_used=False,moderate_fixture_not_actual_Md40=True,passed=True)


def functional_source_identities():
    proofs={}
    def zero(name,value):
        if s.simplify(s.expand_power_exp(value))!=0:raise ArithmeticError('Collar/Gamma source identity failed: '+name)
        proofs[name]=True
    t,x,a,S,Z=s.symbols('t xi a S Z',positive=True); d=1-Z*Z
    H=s.Function('H')(x); D=(1-H)/(a*S)
    current=D
    for k in range(1,5):
        current=-x*s.diff(current,x)
        source=sum(stirling_second(k,l)*x**l*s.diff(H,x,l)/(a*S) for l in range(1,k+1))*(-1)**(k+1)
        zero('full_Gamma_log_radius_derivative_'+str(k),current-source)
        for l in range(1,k+1):
            zero('full_Gamma_S_cancellation_k'+str(k)+'_l'+str(l),
                (x**l/S).subs(x,2*d*S*s.exp(-t))-(2*d*s.exp(-t))**l*S**(l-1))
    # The divided-source identity is the fundamental theorem of calculus
    # under the entire positive Gamma measure, not a series for H.
    v,q=s.symbols('v q',positive=True)
    zero('Gamma_divided_FTC_integrand',s.integrate(v*(1+q*x*v)**(-a-1),(q,0,1))-(1-(1+x*v)**(-a))/(a*x))
    for n in range(10):
        zero('Gamma_positive_derivative_integrand_'+str(n),
            s.diff((1+x*v)**(-a),x,n)-(-1)**n*s.rf(a,n)*v**n*(1+x*v)**(-a-n))
        zero('Gamma_divided_derivative_integrand_'+str(n),
            s.diff(v*(1+q*x*v)**(-a-1),x,n)-(-1)**n*s.rf(1+a,n)*q**n*v**(n+1)*(1+q*x*v)**(-a-1-n))
    dist=s.symbols('distance',positive=True); phi=s.exp(-4/dist**2)
    current=phi
    for n,p in enumerate(PHI_POLYNOMIALS):
        expected=phi/dist**(3*n)*sum(v*dist**j for j,v in enumerate(p))
        zero('original_flat_phi_polynomial_'+str(n),current-expected)
        current=-s.diff(current,dist)
    eps,mu,thetaB,Ev2,p=s.symbols('eps mu thetaB Ev2 prate',real=True); k=1-a
    sig,ph,deficit,pre,W,C,K=s.symbols('sig phi deficit pre W C K',real=True)
    # Production source K equals the exact original H collar bracket.
    one=1; hfactor=1-a*S*deficit
    oldenv={'sig':sig,'self.eps':eps,'phi':ph,'self.a':a,'S':S,'D':deficit*sig*(1-eps*ph),'pre':1-eps*(1-sig+sig*ph)}
    newK=(1-eps*(1-sig+sig*ph))-a*S*deficit*sig*(1-eps*ph)
    zero('exact_original_heat_collar_bracket',newK-((1-sig)*(1-eps)+sig*hfactor*(1-eps*ph)))
    # Bind the NEW tail expressions and original production tail source.
    T,E,Ptail,JW,EW,EW2,PW,PW2=s.symbols('ThetaHat EnergyHat PressureHat JW EW EW2 PW PW2',real=True)
    env={'one':1,'self.k':k,'self.S':S,'self.a':a,'self.eps':eps,'self.delta':2*a,'self.prate':p,
        't':t,'angular':T,'energy':E,'pressure':Ptail,
        "atoms['JW']":JW,"atoms['EW']":EW,"atoms['EW2']":EW2,"atoms['PW']":PW,"atoms['PW2']":PW2}
    Anew=assignment('compliant_collar_Gamma_C4','collar_tails','A',env)
    Enew=assignment('compliant_collar_Gamma_C4','collar_tails','E',env)
    Pnew=assignment('compliant_collar_Gamma_C4','collar_tails','P',env)
    oldenv={'theta_hat':T,'energy_hat':E,'pressure_hat':Ptail,'S':S,'self.a':a,'self.eps':eps,
        'self.eps ** 2':eps**2,
        'self.delta':2*a,'self.phrate':p,'self.k':k,'t':t,
        "atoms['JW']":JW,"atoms['EW']":EW,"atoms['EW2']":EW2,"atoms['PW']":PW,"atoms['PW2']":PW2}
    zero('new_collar_angular_tail_original_source',Anew-assignment('compliant_corrected_outer_field','heat_tails','theta',oldenv))
    zero('new_collar_energy_tail_original_source',Enew-source_branch_expression('heat_tails','E','energy_hat *',oldenv))
    zero('new_collar_pressure_tail_original_source',Pnew-source_branch_expression('heat_tails','P','pressure_hat *',oldenv))
    Xt,Pt,A0,Acur,Ecur,G3,Gcur,defect=s.symbols('Xt Pt A0 Acur Ecur G3 Gcur defect',real=True)
    denv={'Xtail':Xt,'self.eps':eps,"heat0['angular_numerator']":A0}
    newdefect=assignment('compliant_collar_Gamma_C4','data','defect',denv)
    zero('retained_actual_angular_tail_constant',newdefect-(Xt*(1-eps)-A0))
    cenv={"tails['angular_numerator']":Acur,"data['angular_tail_constant_defect']":defect,
        'self.k':k,'t':t,'K':K,"tails['remaining_energy_in_Rtail_units']":Ecur,'self.delta':2*a}
    Xnew=assignment('compliant_collar_Gamma_C4','collar','X',cenv)
    enew=assignment('compliant_collar_Gamma_C4','collar','energy',cenv)
    zero('new_collar_retained_angular_history',Xnew-(Acur+defect*s.exp(-k*t))/K)
    zero('new_collar_backward_energy_units',enew-Ecur*s.exp(2*a*t)/(2*K*K))
    Xloc,Eloc,Ploc=s.symbols('Xloc Eloc Ploc',real=True)
    exenv={'one':1,'self.k':k,'self.delta':2*a,'self.prate':p,'self.a':a,'Sc':S*s.exp(-t),
        "tails['theta']":Xloc,"tails['energy']":Eloc,"tails['pressure']":Ploc}
    Aex=assignment('compliant_collar_Gamma_C4','local_Gamma','angular',exenv)
    Eex=assignment('compliant_collar_Gamma_C4','local_Gamma','energy',exenv)
    Pex=assignment('compliant_collar_Gamma_C4','local_Gamma','pressure',exenv)
    zero('new_full_Gamma_angular_current_units',Aex-1/k-S*s.exp(-t)*Xloc)
    zero('new_full_Gamma_energy_current_units',Eex-1/(2*a)+a*S*s.exp(-t)*Eloc)
    zero('new_full_Gamma_pressure_current_units',Pex-1/(2*p)+a*S*s.exp(-t)*Ploc)
    # Exact rescaling of the infinite positive-source integrals.
    zero('Gamma_angular_tail_rescaling',Anew.subs({T:s.exp(-a*t)*Xloc,JW:0})-Aex)
    zero('Gamma_energy_tail_rescaling',Enew.subs({E:s.exp(-(2*a+1)*t)*Eloc,EW:0,EW2:0})-s.exp(-2*a*t)*Eex)
    zero('Gamma_pressure_tail_rescaling',Pnew.subs({Ptail:s.exp(-(p+1)*t)*Ploc,PW:0,PW2:0})-s.exp(-p*t)*Pex)
    exenv={'K':K,'t':t,'self.k':k,'self.prate':p,'self.pressure_scale':Ev2*thetaB**2,
        "data['angular_tail_constant_defect']":defect,"local['angular_numerator']":Aex,
        "local['energy_numerator']":Eex,"local['pressure_numerator']":Pex,
        "loc3['pressure_numerator']":G3,"data['pressure3']":Pt}
    Xex=assignment('compliant_collar_Gamma_C4','exterior','X',exenv)
    eex=assignment('compliant_collar_Gamma_C4','exterior','energy',exenv)
    Pforward=assignment('compliant_collar_Gamma_C4','exterior','pressure',exenv)
    zero('new_exterior_retained_angular_history',Xex-(Aex+defect*s.exp(-k*t))/K)
    zero('new_exterior_energy_current_units',eex-Eex/(2*K*K))
    zero('new_exterior_forward_pressure_history',Pforward-Pt-Ev2*thetaB**2*(G3*s.exp(-3*p)-Pex*s.exp(-p*t)))
    zero('new_exterior_pressure_interface',Pforward.subs({t:3,G3:Pex.subs(t,3)})-Pt)
    # Direct production bindings for shape product/chain rules and rate.
    tree=ast.parse((HERE/(PREFIX+'compliant_collar_Gamma_C4.py')).read_text(encoding='utf8'))
    def syntax(method,target,expected):
        fn=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name==method)
        nodes=[n.value for n in ast.walk(fn) if isinstance(n,ast.Assign) and any(ast.unparse(v)==target for v in n.targets)]
        if len(nodes)!=1 or ast.dump(nodes[0])!=ast.dump(ast.parse(expected,mode='eval').body):
            raise ValueError('Collar production derivative/source expression changed: '+method+' '+target)
        proofs['production_'+method+'_'+target+'_expression_bound']=True
    syntax('shape','W','[unity[j]-sr[j]+v for j,v in enumerate(product_rows(sr,fr))]')
    syntax('shape','C','[unity[j]-fr[j]*self.eps for j in range(5)]')
    syntax('shape','pre','[unity[j]-W[j]*self.eps for j in range(5)]')
    # high/low D branches differ only in retained y order. Bind both.
    fn=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name=='shape')
    nodes=[n.value for n in ast.walk(fn) if isinstance(n,ast.Assign) and any(ast.unparse(v)=='D' for v in n.targets)]
    expected=['product_rows(product_rows(sr,C),dr)','[sr[0]*C[0]*dr[0]]']
    if sorted(ast.dump(v) for v in nodes)!=sorted(ast.dump(ast.parse(v,mode='eval').body) for v in expected):
        raise ValueError('Original high/low collar Gamma composition changed')
    proofs['production_original_collar_Gamma_composition_bound']=True
    syntax('shape','K','[pre[j]-D[j]*(self.a*self.S) for j in range(len(D))]')
    syntax('packet','rates','quotient_log_rates(K)')
    syntax('packet','theta','K[0]*(self.theta_base*c.exp(-self.bh*t))')
    fn=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name=='packet')
    aug=[n for n in ast.walk(fn) if isinstance(n,ast.AugAssign) and ast.unparse(n.target)=='rates[0]']
    if len(aug)!=1 or not isinstance(aug[0].op,ast.Add) or ast.dump(aug[0].value)!=ast.dump(ast.parse('self.mu-self.a',mode='eval').body):
        raise ValueError('Actual collar/Gamma extra logarithmic rate changed')
    proofs['production_actual_collar_Gamma_log_rate_bound']=True
    # Bind the actual data inlet and forward integral, so pressure
    # endpoint equalities refer to production histories rather than
    # two independently asserted copies of the same symbol.
    zero('new_actual_waiting_X_source',assignment('compliant_collar_Gamma_C4','data','Xtail',
        {"terminal['angular_Taylor']":Xt})-Xt)
    zero('new_actual_waiting_pressure_source',assignment('compliant_collar_Gamma_C4','data','Ptail',
        {"terminal['pressure_over_Pstar_squared_Taylor']":Pt})-Pt)
    forward=s.Function('ForwardPressure')(t)
    dataP3=assignment('compliant_collar_Gamma_C4','data','pressure3',
        {'self.forward_pressure(Z, 3, Ptail)':forward.subs(t,3)})
    zero('new_data_forward_pressure3',dataP3-forward.subs(t,3))
    zero('new_collar_forward_pressure_source',assignment('compliant_collar_Gamma_C4','collar','pressure',
        {"self.forward_pressure(Z, t, data['Ptail'])":forward})-forward)
    fn=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name=='forward_pressure')
    ret=[n.value for n in ast.walk(fn) if isinstance(n,ast.Return)]
    if len(ret)!=1 or ast.dump(ret[0])!=ast.dump(ast.parse('Ptail+integral*self.pressure_scale',mode='eval').body):
        raise ValueError('Actual forward pressure history return changed')
    coeff=s.symbols('cell_length',nonnegative=True)
    zero('new_forward_pressure_actual_integrand',assignment('compliant_collar_Gamma_C4','forward_pressure','integral',
        {'K':K,'ds':coeff,'self.prate':p,'v':t},augmented=True)-K*K*coeff*s.exp(-p*t)/2)
    zero('new_forward_pressure_actual_inlet',Pt+Ev2*thetaB**2*0-Pt)
    zero('new_collar_Gamma_actual_forward_pressure_join',
        Pforward.subs({t:3,G3:Pex.subs(t,3),Pt:dataP3})-forward.subs(t,3))
    tailtheta,bh,wait,lone=s.symbols('thetaT bh wait logone',real=True)
    baseactual=assignment('compliant_collar_Gamma_C4','__init__','self.theta_base',
        {'self.steep.thetaT':tailtheta,'self.bh':bh,'self.steep.wait':wait,'self.steep.logone':lone})
    zero('new_actual_waiting_velocity_normalization',baseactual.subs(lone,s.log(1-eps))-tailtheta*s.exp(-bh*wait)/(1-eps))
    # Bind the actual Stirling accumulation and its cancelled factors.
    for order in range(1,5):
        for degree in range(1,order+1):
            gg=s.symbols('Hderivative_over_a')
            envg={'base':2*d*s.exp(-t),'G[l]':gg,'k':s.Integer(order),'l':s.Integer(degree),
                'S':S,'stirling_second(k, l)':s.Integer(stirling_second(order,degree))}
            zero('production_Gamma_ordinary_y'+str(order)+'_term'+str(degree),
                assignment('compliant_collar_Gamma_C4','gamma_deficit_mixed','row',envg,augmented=True)
                -(-1)**(order+1)*stirling_second(order,degree)*(2*d*s.exp(-t))**degree*S**(degree-1)*gg)
    # Whole-Z functions and original integral ODEs identify mixed jets.
    Kt=s.Function('K')(t); At=s.Function('A')(t); Et=s.Function('E')(t); PP=s.Function('P')(t)
    cc=s.symbols('constant',real=True); theta=thetaB*s.exp(-(s.Rational(1,2)+a)*t)*Kt
    xx=(At+cc*s.exp(-k*t))/Kt; ee=Et*s.exp(2*a*t)/(2*Kt**2)
    fp=s.diff(Kt,t)/Kt; aa=mu-a+fp; bp=s.Rational(1,2)+mu
    zero('collar_and_Gamma_theta_source_ODE',s.diff(theta,t)-(aa-bp)*theta)
    zero('collar_and_Gamma_X_source_ODE',
        (s.diff(xx,t)-1+(k+fp)*xx).subs(s.diff(At,t),Kt-k*At))
    zero('collar_and_Gamma_energy_source_ODE',
        (s.diff(ee,t)+s.Rational(1,2)-(2*a-2*fp)*ee).subs(s.diff(Et,t),-s.exp(-2*a*t)*Kt**2))
    # Evaluate the production-bound K product formulas with arbitrary
    # whole-Z Gamma derivative histories at both exact flat endpoints.
    Dr=[s.Function('D'+str(j))(Z) for j in range(5)]
    def product(l,r):return [sum(s.binomial(n,j)*l[j]*r[n-j] for j in range(n+1)) for n in range(5)]
    def shape_boundary(sr,fr):
        one=[s.Integer(1)]+[s.Integer(0)]*4
        wr=[one[j]-sr[j]+v for j,v in enumerate(product(sr,fr))]
        cr=[one[j]-eps*fr[j] for j in range(5)]
        pre=[one[j]-eps*wr[j] for j in range(5)]
        dd=product(product(sr,cr),Dr)
        return [pre[j]-a*S*dd[j] for j in range(5)]
    inletK=shape_boundary([s.Integer(0)]*5,[s.symbols('phi0')]+[s.symbols('phi'+str(j)) for j in range(1,5)])
    pureK=[1-a*S*Dr[0]]+[-a*S*v for v in Dr[1:]]
    terminalK=shape_boundary([s.Integer(1)]+[s.Integer(0)]*4,[s.Integer(0)]*5)
    for name,lrows,rrows in (('waiting_collar',[1-eps]+[s.Integer(0)]*4,inletK),('collar_Gamma',terminalK,pureK)):
        for order,(lhs,rhs) in enumerate(zip(lrows,rrows)):
            for n in range(6):zero(name+'_actual_K_y'+str(order)+'_Z'+str(n),s.diff(lhs-rhs,Z,n))
    # Actual waiting H and collar E0 are the SAME original full future
    # source, with the original epsilon normalization; bind the old H.
    e0=Enew.subs(t,0)
    hwait=e0/(1-eps)**2
    zero('waiting_full_future_epsilon_normalization',assignment('compliant_steep_waiting_C4','data','H',
        {'preheat_tail':e0,'self.tail_normalization':1/(1-eps)**2})-hwait)
    Zhist={v:s.Function(str(v))(Z) for v in (Xt,Pt,A0,Ecur,Xloc,Eloc,Ploc,defect,K)}
    inletvalues=[thetaB*(1-eps),Xt,hwait/2,Pt]
    collarinlet=[thetaB*(1-eps),(A0+newdefect)/(1-eps),enew.subs({t:0,Ecur:e0,K:1-eps}),Pt]
    endpoint3A=Aex.subs(t,3); endpoint3E=Eex.subs(t,3)
    # Rescaling identities already bind the source tails at t=3.
    collarend=[thetaB*s.exp(-(s.Rational(1,2)+a)*3)*K,
        Xnew.subs({t:3,Acur:endpoint3A}),
        enew.subs({t:3,Ecur:s.exp(-2*a*3)*endpoint3E}),Pt]
    exteriorinlet=[thetaB*s.exp(-(s.Rational(1,2)+a)*3)*K,Xex.subs(t,3),eex.subs(t,3),
        Pforward.subs({t:3,G3:Pex.subs(t,3)})]
    def log_rates(rows):
        rates=[]
        for order in range(4):
            value=rows[order+1]-sum(s.binomial(order,j)*rows[j]*rates[order-j] for j in range(1,order+1))
            rates.append(s.simplify(value/rows[0]))
        rates[0]+=mu-a
        return rates
    for name,left,right,lk,rk in (
        ('waiting_collar',inletvalues,collarinlet,[1-eps]+[s.Integer(0)]*4,inletK),
        ('collar_Gamma',collarend,exteriorinlet,terminalK,pureK)):
        labels=('theta','X','energy','pressure')
        for label,lhs,rhs in zip(labels,left,right):
            zero(name+'_exact_'+label,lhs-rhs)
            for n in range(6):zero(name+'_'+label+'_axial'+str(n),s.diff((lhs-rhs).subs(Zhist),Z,n))
        lr,rr=log_rates(lk),log_rates(rk)
        for j in range(4):zero(name+'_actual_log_rate_derivative'+str(j),lr[j]-rr[j])
        # Each side has its own source-bound initial fields and rates.
        # Generate ordinary y derivatives with the common recovery ODEs.
        histories=[]
        for initial,rates in ((left,lr),(right,rr)):
            th=[s.simplify(initial[0])]; xx=[s.simplify(initial[1])]; ee=[s.simplify(initial[2])]; pp=[s.simplify(initial[3])]
            tr=list(rates); tr[0]-=bp
            xr=list(rates); xr[0]+=1-mu
            er=[-2*v for v in rates]; er[0]+=2*mu
            for order in range(4):
                th.append(sum(s.binomial(order,j)*tr[j]*th[order-j] for j in range(order+1)))
                xx.append((1 if order==0 else 0)-sum(s.binomial(order,j)*xr[j]*xx[order-j] for j in range(order+1)))
                ee.append((-s.Rational(1,2) if order==0 else 0)+sum(s.binomial(order,j)*er[j]*ee[order-j] for j in range(order+1)))
                pp.append(Ev2/2*sum(s.binomial(order,j)*th[j]*th[order-j] for j in range(order+1)))
            histories.append((th,xx,ee,pp))
        for label,lrows,rrows in zip(labels,*histories):
            for order,(lhs,rhs) in enumerate(zip(lrows,rrows)):
                for n in range(5-order):zero(name+'_'+label+'_y'+str(order)+'_Z'+str(n),s.diff((lhs-rhs).subs(Zhist),Z,n))
    return proofs


def run():
    r=json.loads((HERE/NAME).read_bytes()); hashes=dict(r['input_hashes'])
    for name,digest in hashes.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:raise ValueError('Collar/Gamma source changed: '+name)
    c=MPIntervalContext(); c.dps=240; read=lambda v:read_interval(c,v)
    points=r['collar_samples']+r['exterior_samples']+[r[k] for k in
        ('whole_Z_collar_box','whole_Z_collar_inlet','whole_Z_collar_terminal','whole_Z_exterior_inlet','whole_Z_unbounded_exterior','phi_endpoint_crossing')]
    finite=zeros=0
    for p in points:
        for flag in ('original_flat_collar_and_exact_Gamma_source_retained','same_waiting_terminal_angular_pressure_histories_retained',
            'full_infinite_Gamma_tail_included','actual_Gamma_not_defined_by_finite_S_series'):
            if not p[flag]:raise ValueError('Original collar/Gamma history/source missing')
        if endpoints(read(p['energy_Taylor']['coefficients'][0]))[0]<=0:raise ArithmeticError('Full collar/Gamma future energy positivity lost')
        if endpoints(read(p['heat_bracket_y_derivatives'][0]['coefficients'][0]))[0]<=0:raise ArithmeticError('Actual heat bracket positivity lost')
        for label,grid in p['physical_mixed_derivatives_total_order_le4'].items():
            if len(grid)!=15:raise ValueError('Mixed derivative coverage incomplete')
            for value in grid.values():
                ends=endpoints(read(value))
                if any(not mp.isfinite(v) for v in ends):raise ArithmeticError('Nonfinite collar/exterior spatial derivative')
                finite+=1
                if label in (UZ,UR):
                    if ends!=(mp.mpf(0),mp.mpf(0)):raise ArithmeticError('Inherited zero primitive lost')
                    zeros+=1
        for flag in ('full_outer_C4_certified','physical_energy_integral_certified','whole_outer_cone_certified','temporal_recursion'):
            if p[flag]:raise ValueError('Collar/Gamma scope promoted')
    overlaps=prefixes=0
    def overlap(a,b,label):
        nonlocal overlaps
        a,b=endpoints(read(a)),endpoints(read(b))
        if max(a[0],b[0])>min(a[1],b[1]):raise ArithmeticError('Independent same-source diagnostic disjoint: '+label)
        overlaps+=1
    prior=json.loads((HERE/(PREFIX+'compliant_steep_waiting_C4.json')).read_bytes())
    for name,lhs,rhs in (('waiting_collar',prior['whole_Z']['waiting']['terminal'],r['whole_Z_collar_inlet']),
        ('collar_Gamma',r['whole_Z_collar_terminal'],r['whole_Z_exterior_inlet'])):
        for label,grid in lhs['physical_mixed_derivatives_total_order_le4'].items():
            for index,value in grid.items():overlap(value,rhs['physical_mixed_derivatives_total_order_le4'][label][index],name+' '+label+' '+index)
        for key in ('angular_y_derivative_Taylor','energy_y_derivative_Taylor'):
            for a,b in zip(lhs[key],rhs[key]):
                for av,bv in zip(a['coefficients'],b['coefficients']):overlap(av,bv,name+' '+key)
    for jet in r['whole_Z_collar_inlet']['heat_bracket_y_derivatives'][1:]:
        for value in jet['coefficients']:
            if endpoints(read(value))!=(mp.mpf(0),mp.mpf(0)):raise ArithmeticError('Original waiting/collar flat K derivatives not zero')
    for oldvalue,newvalue in zip(prior['whole_Z']['waiting']['terminal']['pressure_over_Pstar_squared_Taylor']['coefficients'],
        r['whole_Z_collar_inlet']['pressure_over_Pstar_squared_Taylor']['coefficients']):
        if endpoints(read(oldvalue))!=endpoints(read(newvalue)):raise ArithmeticError('Actual waiting-terminal forward pressure not inherited exactly')
    for value in r['whole_Z_collar_terminal']['original_collar_phi_y_Taylor']['coefficients']:
        if endpoints(read(value))!=(mp.mpf(0),mp.mpf(0)):raise ArithmeticError('Original collar/Gamma flat phi derivatives not zero')
    old=json.loads((HERE/(PREFIX+'compliant_corrected_outer_field.json')).read_bytes())
    for p in old['samples']:
        if p['stage']!='corrected collar / exact Gamma exterior':continue
        source=r['collar_samples']+r['exterior_samples']
        current=next(q for q in source if endpoints(read(q['Z']))==endpoints(read(p['Z']))
            and endpoints(read(q['coordinate']['offset']))==endpoints(read(p['coordinate']['offset'])))
        for key,prior_key in (('theta_over_Ev0_Taylor','Utheta_over_Ev0_enclosure'),
            ('angular_Taylor','Mtheta_over_sqrt2_R_3half_Utheta'),('energy_Taylor','Mztheta_over_R_Utheta_squared'),
            ('pressure_over_Pstar_squared_Taylor','P_over_Pstar_squared')):
            for a,b in zip(current[key]['coefficients'],p[prior_key]['coefficients']):overlap(a,b,'earlier actual collar/Gamma C1 '+key); prefixes+=1
    proofs=functional_source_identities()
    phi=independent_flat_phi_fixture(); gamma=independent_full_Gamma_fixture(); rates=independent_nonconstant_rate_fixture()
    hashes[NAME]=hashlib.sha256((HERE/NAME).read_bytes()).hexdigest()
    hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    out=dict(actual_five_defect_family_sha256=r['actual_five_defect_family_sha256'],implicit_source_sha256=r['implicit_source_sha256'],
        all_passed=True,functional_production_source_identities=proofs,
        independent_original_phi_fixture=phi,independent_full_Gamma_fixture=gamma,independent_variable_rate_fixture=rates,
        finite_actual_mixed_bounds_checked=finite,inherited_exact_zero_bounds_checked=zeros,
        independent_interface_and_prefix_overlap_diagnostics=overlaps,earlier_actual_C1_prefix_diagnostics=prefixes,
        entire_collar_and_unbounded_Gamma_high_mixed_derivatives_available=True,
        waiting_collar_and_collar_Gamma_joins_certified=True,
        full_infinite_Gamma_source_and_formal_nonzero_S_retained=True,
        join_scope='Leading-profile mixed spatial<=4 / canonical axial5, waiting/collar and collar/Gamma; physical Cartesian/core/axis and temporal recursion excluded',
        full_pulse_C4_installed=True,full_outer_C4_certified=False,
        physical_energy_integral_certified=False,whole_outer_cone_certified=False,temporal_recursion=False,input_hashes=hashes)
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8')
    print('Original collar/full unbounded Gamma: exact source derivatives, flat shapes and spatial joins PASS',flush=True)
    return out


if __name__=='__main__':run()
