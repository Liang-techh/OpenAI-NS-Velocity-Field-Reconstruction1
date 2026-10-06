"""Actual quiet-power bump jets, partial primitives and functional exit bridge."""
import math
import ast
from types import SimpleNamespace
import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_current_O3_independent_repair_operator import (
    SourceAST,exp_average,endpoints)
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import beta_jets,BETA_POLYNOMIALS
from lei_ren_part1_paper_compliant_outer_pulse_map import raw_beta

def bump_rows(c,y,center,ell,normalization):
    raw=beta_jets(c,(y-center)/ell)
    beta=[raw[k]*math.factorial(k)/(ell**(k+1)*normalization) for k in range(5)]
    return [c.exp(-y)*sum((math.comb(n,k)*(-1)**(n-k)*beta[k] for k in range(n+1)),c.mpf(0))
      for n in range(5)]

def partial_weight(c,y,center,ell,normalization,power,multiplicity,full,cells,mu=None):
    """Enclose the same continuous partial integral, exact full/zero branches.

    A non-None mu selects the signed divided axial weight instead of a
    positive x^power weight. No inverse of an interval meeting an edge.
    """
    if mu is not None and (endpoints(mu)[0]<=0 or endpoints(center-ell)[0]<=0 or endpoints(full)[1]>=0):
        raise ValueError('Signed D clamp requires positive mu, positive whole logx support and negative full weight')
    def endpoint(value):
        r=(c.mpf(value)-center)/ell;lo,hi=endpoints(r)
        if hi<=-1:return c.mpf(0)
        if lo>=1:return full
        stop=c.mpf([max(mp.mpf(-1),lo),min(mp.mpf(1),hi)])
        length=stop+1
        length=c.mpf([max(mp.mpf(0),endpoints(length)[0]),endpoints(length)[1]])
        total=c.mpf(0)
        for i in range(cells):
            a=-1+length*i/cells;b=-1+length*(i+1)/cells
            box=c.mpf([endpoints(a)[0],endpoints(b)[1]]);raw=raw_beta(c,box)
            logx=center+ell*box
            mass=raw*(length/cells)/normalization
            if mu is not None:
                value=mass*(-logx)*exp_average(c,-mu*logx)
            elif multiplicity==1:
                value=mass*c.exp(power*logx)
            else:
                value=raw**2*(length/cells)/(ell*normalization**2)*c.exp((power-1)*logx)
            total+=value
        if mu is None:
            return c.mpf([max(mp.mpf(0),endpoints(total)[0]),min(endpoints(full)[1],endpoints(total)[1])])
        return c.mpf([max(endpoints(full)[0],endpoints(total)[0]),min(mp.mpf(0),endpoints(total)[1])])
    lo,hi=endpoints(y);left=endpoint(lo);right=endpoint(hi)
    return c.mpf([min(endpoints(left)[0],endpoints(right)[0]),max(endpoints(left)[1],endpoints(right)[1])])

def partial_weights(repair,y,cells):
    c=repair.ctx;mu=repair.mu;W=repair.weights;ci=W['centers'];ell=W['radius'];normalization=repair.normalization
    alpha=c.mpf('.5')+mu;H=W['H'];matrix=repair.matrix
    out={key:[] for key in ('mass','D','I','S','Cp','cross','energy','pressure')}
    for i,center in enumerate(ci):
        for name,p,full in (('mass',c.mpf(0),c.mpf(1)),
          ('I',c.mpf('.5'),H['I']*c.exp(center/2)),
          ('S',-alpha,H['S']*c.exp(-alpha*center)),
          ('Cp',-alpha-1,H['Cp']*c.exp((-alpha-1)*center))):
            out[name].append(partial_weight(c,y,center,ell,normalization,p,1,full,cells))
        for name,p,full in (('energy',c.mpf(0),matrix['energy_weights'][i]),
          ('pressure',c.mpf(-1),matrix['pressure_weights'][i])):
            out[name].append(partial_weight(c,y,center,ell,normalization,p,2,full,cells))
        if i in (0,2):
            j=0 if i==0 else 1
            out['D'].append(partial_weight(c,y,center,ell,normalization,c.mpf(0),1,W['divided_axial_rows'][j],cells,mu=mu))
            out['cross'].append(partial_weight(c,y,center,ell,normalization,c.mpf('.5'),2,matrix['cross_weights'][j],cells))
    return out

def partial_repair_primitives(c,mu,N,h,W):
    a=(h[0],h[1]);e=h[2:];sa=c.sqrt(mu)/N;se=mu/N
    mass=W['mass'];D=W['D'];I=W['I'];S=W['S'];Cp=W['Cp'];cross=W['cross'];energy=W['energy'];pressure=W['pressure']
    M=sa*(a[0]*mass[0]+a[1]*mass[2])
    J=sa*(a[0]*(mass[0]+mu*D[0])+a[1]*(mass[2]+mu*D[1]))
    J+=sa*se*(a[0]*e[0]*cross[0]+a[1]*e[2]*cross[1])
    Itotal=se*sum((e[i]*I[i] for i in range(3)),sa*0)
    Stotal=sa**2*(a[0]**2*energy[0]+a[1]**2*energy[2])-se*sum((e[i]*S[i] for i in range(3)),sa*0)
    Stotal-=se**2/2*sum((e[i]**2*energy[i] for i in range(3)),sa*0)
    Ptotal=se*sum((e[i]*Cp[i] for i in range(3)),sa*0)+se**2/2*sum((e[i]**2*pressure[i] for i in range(3)),sa*0)
    return dict(M=M,J=J,I=Itotal,S=Stotal,Cp=Ptotal)

def actual_repaired_history_theorem():
    asts=SourceAST();checks={}
    def zero(name,a,b):
        if s.simplify(s.expand(a-b))!=0:raise ArithmeticError('Repaired source identity: '+name)
        checks[name]=True
    r,y,ell,B0,center,p=s.symbols('raw_coordinate logx ell same_raw_mass center p',real=True)
    raw=s.exp(-1/(1-r*r))
    for n,polynomial in enumerate(BETA_POLYNOMIALS):
        poly=sum(v*r**i for i,v in enumerate(polynomial))
        zero('actual_reused_raw_beta_derivative_'+str(n),s.diff(raw,r,n),raw*poly/(1-r*r)**(2*n))
    asts.method('flat_pulse_derivatives','beta_jets');asts.method('flat_pulse_derivatives','beta_polynomials')
    beta=s.Function('same_normalized_beta')(y-center)
    for n in range(5):
        zero('actual_log_translated_bump_ordinary_row_'+str(n),s.diff(s.exp(-y)*beta,y,n),
          s.exp(-y)*sum(s.binomial(n,k)*(-1)**(n-k)*s.diff(beta,y,k) for k in range(n+1)))
    zero('actual_single_partial_log_FTC_density',s.exp((p+1)*y)*(s.exp(-y)*beta),s.exp(p*y)*beta)
    zero('actual_product_partial_log_FTC_density',s.exp((p+1)*y)*(s.exp(-y)*beta)**2,s.exp((p-1)*y)*beta**2)
    t=s.Symbol('edge_distance',positive=True)
    for n in range(9):
        if s.limit(t**(-n)*s.exp(-1/t),t,0,dir='+')!=0:raise ArithmeticError('Raw flat jet limit failed')
        checks['actual_flat_edge_polynomial_exp_limit_'+str(n)]=True
    mu,N=s.symbols('mu N',positive=True);ctx=SimpleNamespace(sqrt=s.sqrt)
    h=s.symbols('a0 a2 e0 e1 e2');D=s.symbols('D0 D2');li=s.symbols('I0 I1 I2');ls=s.symbols('S0 S1 S2');lp=s.symbols('P0 P1 P2')
    cross=s.symbols('C0 C2');energy=s.symbols('E0 E1 E2');pressure=s.symbols('Q0 Q1 Q2')
    W=dict(mass=[s.Integer(1)]*3,D=D,I=li,S=ls,Cp=lp,cross=cross,energy=energy,pressure=pressure)
    fn=asts.replay('current_O3_repaired_histories_operator','partial_repair_primitives',dict())
    primitives=fn(ctx,mu,N,h,W);sa=s.sqrt(mu)/N;se=mu/N
    B=[[1,1,0,0,0],[D[0],D[1],0,0,0],[0,0,*li],[0,0,*[-v for v in ls]],[0,0,*lp]]
    qfn=asts.replay('current_O3_independent_repair_operator','scaled_quadratic',dict())
    Q=qfn(ctx,mu,N,dict(cross_weights=cross,energy_weights=energy,pressure_weights=pressure),h)
    scaled=(primitives['M']/sa,(primitives['J']-primitives['M'])/(mu*sa),
      primitives['I']/se,primitives['S']/se,primitives['Cp']/se)
    mapped=[sum(B[i][j]*h[j] for j in range(5))+Q[i] for i in range(5)]
    for i,name in enumerate(('M','D','I','S','Cp')):zero('actual_full_partial_integrals_reproduce_checked_repair_map_'+name,scaled[i],mapped[i])
    # Actual parent source defects, not defects chosen as minus the map.
    # The accepted previous source bridge identifies D before enclosure.
    defects=s.symbols('same_actual_scaled_dM same_actual_scaled_dD same_actual_scaled_dI same_actual_scaled_dS same_actual_scaled_dCp')
    incoming=dict(M=sa*defects[0],J=sa*(defects[0]+mu*defects[1]),
      I=se*defects[2],S=se*defects[3],Cp=se*defects[4])
    residual=[mapped[i]+defects[i] for i in range(5)]
    factored=dict(M=sa*residual[0],J=sa*(residual[0]+mu*residual[1]),
      I=se*residual[2],S=se*residual[3],Cp=se*residual[4])
    for name in ('M','J','I','S','Cp'):
        zero('actual_terminal_delta_is_checked_implicit_equation_residual_'+name,
          primitives[name]+incoming[name],factored[name])
    checks['actual_zero_exit_uses_checked_same_source_residuals_not_chosen_defects']=True
    # Bind the actual field/control lineage and actual repair inlet scalar
    # definition, including the parent transport from t=2 to t=2+y.
    asts.expression('current_O3_repaired_histories','__init__','self.histories',wanted='self.repair.histories')
    asts.expression('current_O3_repaired_histories','history','parent',wanted='self.histories.history(region,Z,coordinate)')
    asts.expression('current_O3_repaired_histories','history','h',wanted='self.repair.controls')
    asts.expression('current_O3_independent_repair','__init__','self.defects',wanted='correlated_scaled_defects(self.histories)')
    asts.expression('current_O3_independent_repair','__init__','self.controls',
      wanted='coefficient_enclosure(c,self.mu,self.N,self.matrix,self.defects,self.certificate)')
    asts.expression('current_O3_independent_repair_operator','correlated_scaled_defects','source',wanted='field.scalars(2)')
    checks['actual_same_field_inlet_defects_and_unique_control_vector_program_bound']=True
    inlet=s.symbols('same_Bm2 same_Bh2 same_Bk2 same_Be2 same_Bp2')
    f2=s.Symbol('actual_f2',positive=True)
    source_env=dict(N=N,c=SimpleNamespace(sqrt=s.sqrt),mu=mu,source=dict(zip(('m','h','k','e','p'),inlet)),f2=f2)
    dm_program=asts.evaluate(asts.expression('current_O3_independent_repair_operator','correlated_scaled_defects','dM'),source_env)
    zero('actual_inlet_M_scaled_defect_from_parent_source',dm_program,inlet[0]/(f2*sa))
    source_node=asts.method('current_O3_independent_repair_operator','correlated_scaled_defects')
    source_return=next(n.value for n in ast.walk(source_node) if isinstance(n,ast.Return))
    scaled_defect=asts.evaluate(source_return,{**source_env,'dM':dm_program,'dD':defects[1]})
    for i,name in ((2,'I'),(3,'S'),(4,'Cp')):
        expected_inlet={'I':inlet[1]/(f2*se),'S':inlet[3]/(f2*f2*se),'Cp':inlet[4]/(f2*f2*se)}[name]
        zero('actual_inlet_'+name+'_scaled_defect_from_parent_source',scaled_defect[i],expected_inlet)
    # The prior accepted exact source correlation is consumed at runtime;
    # never enclose (Bk/f2^2-Bm/f2)/mu by independent boxes.
    checks['actual_inlet_J_uses_prior_source_correlated_D_identity']=True
    asts.expression('current_O3_modulated_histories','history','t',wanted='1+v')
    asts.expression('current_O3_modulated_histories','history','scalars',wanted='self.scalars(t)')
    cumulative=asts.method('current_O3_modulated_histories_operator','cumulative_scalar_enclosures')
    asts.expression('current_O3_modulated_histories_operator','cumulative_scalar_enclosures','upper_t',wanted="min(hi,mp.mpf('.5'))")
    asts.expression('current_O3_modulated_histories_operator','cumulative_scalar_enclosures','lower_t',wanted="min(lo,mp.mpf('.5'))")
    checks['actual_inlet_and_exit_share_same_fixed_full_modulation_integral']=True
    transport=next(n.value for n in ast.walk(cumulative) if isinstance(n,ast.AugAssign)
      and ast.unparse(n.target)=='result[name]' and ast.dump(n.value)==ast.dump(ast.parse('c.exp(-c.mpf(str(rate))*t)',mode='eval').body))
    expctx=SimpleNamespace(exp=s.exp,mpf=s.Rational)
    for name,rate in (('m',1),('h',s.Rational(3,2)),('k',s.Rational(3,2)),('e',1),('p',0)):
        if rate:
            exit_transport=asts.evaluate(transport,dict(c=expctx,rate=rate,t=2+y))
            inlet_transport=asts.evaluate(transport,dict(c=expctx,rate=rate,t=s.Integer(2)))
            zero('actual_parent_post_support_transport_from_inlet_'+name,exit_transport/inlet_transport,s.exp(-rate*y))
        else:checks['actual_parent_post_support_Cp_is_same_constant']=True
    # Bind the full partial weights to the very same checked map parameters.
    partial_fn=asts.method('current_O3_repaired_histories_operator','partial_weight')
    guard=ast.parse('lo>=1',mode='eval').body
    if not any(isinstance(n,ast.If) and ast.dump(n.test)==ast.dump(guard)
      and any(isinstance(v,ast.Return) and isinstance(v.value,ast.Name) and v.value.id=='full' for v in n.body)
      for n in ast.walk(partial_fn)):raise ValueError('Same exact full-support partial branch changed')
    asts.expression('current_O3_independent_repair_operator','bump_weights','centers',
      wanted='[c.mpf(1)/5,c.mpf(1)/2,c.mpf(4)/5]')
    asts.expression('current_O3_independent_repair_operator','bump_weights','ell',wanted='c.mpf(1)/40')
    fixed_ci=(s.Rational(1,5),s.Rational(1,2),s.Rational(4,5))
    for i,ci in enumerate(fixed_ci):
        if (1-ci)/s.Rational(1,40)<=1:raise ArithmeticError('Exit precedes a full bump support')
        checks['actual_q2_full_raw_support_guard_'+str(i)]=True
    H=s.symbols('same_HI same_HS same_HCp same_Haxial')
    repair_stub=SimpleNamespace(ctx=SimpleNamespace(mpf=s.Rational,exp=s.exp),mu=mu,
      normalization=s.Symbol('same_positive_raw_normalization'),weights=dict(centers=fixed_ci,radius=s.Rational(1,40),
        H=dict(zip(('I','S','Cp','axial'),H)),divided_axial_rows=D),
      matrix=dict(cross_weights=cross,energy_weights=energy,pressure_weights=pressure))
    partial_fn=asts.replay('current_O3_repaired_histories_operator','partial_weights',
      dict(partial_weight=lambda c,y,center,ell,norm,power,multiplicity,full,cells,mu=None:full))
    actual_full=partial_fn(repair_stub,s.Integer(1),128)
    for i,ci in enumerate(fixed_ci):
        zero('actual_full_mass_same_normalized_weight_'+str(i),actual_full['mass'][i],1)
        for name,power,hfactor in (('I',s.Rational(1,2),H[0]),('S',-s.Rational(1,2)-mu,H[1]),
          ('Cp',-s.Rational(3,2)-mu,H[2])):
            zero('actual_full_'+name+'_same_map_weight_'+str(i),actual_full[name][i],hfactor*s.exp(power*ci))
    for name,expected_weights in (('D',D),('cross',cross),('energy',energy),('pressure',pressure)):
        for i,value in enumerate(expected_weights):zero('actual_full_'+name+'_same_map_weight_'+str(i),actual_full[name][i],value)
    f2,x=s.symbols('actual_f2 actual_x',positive=True)
    M,I,J,S,P=s.symbols('actual_partial_M actual_partial_I actual_partial_J actual_partial_S actual_partial_Cp')
    old_scalars=dict(zip(('m','h','k','e','p'),s.symbols('same_old_Bm same_old_Bh same_old_Bk same_old_Be same_old_Bp')))
    scalar_program=asts.expression('current_O3_repaired_histories','history','scalar',wanted=
      "dict(m=previous['m']+self.f2*primitives['M']/x,h=previous['h']+self.f2*primitives['I']/x**c.mpf('1.5'),k=previous['k']+self.f2**2*primitives['J']/x**c.mpf('1.5'),e=previous['e']+self.f2**2*primitives['S']/x,p=previous['p']+self.f2**2*primitives['Cp'])")
    scalar_values=asts.evaluate(scalar_program,dict(c=SimpleNamespace(mpf=s.Rational),
      self=SimpleNamespace(f2=f2),previous=old_scalars,primitives=dict(M=M,I=I,J=J,S=S,Cp=P),x=x))
    expected=dict(m=old_scalars['m']+f2*M/x,h=old_scalars['h']+f2*I/x**s.Rational(3,2),
      k=old_scalars['k']+f2*f2*J/x**s.Rational(3,2),e=old_scalars['e']+f2*f2*S/x,
      p=old_scalars['p']+f2*f2*P)
    for name in expected:zero('actual_callable_source_partial_physical_units_'+name,scalar_values[name],expected[name])
    inlet_from_actual_defects=dict(m=f2*sa*defects[0],h=f2*se*defects[2],
      k=f2*f2*sa*(defects[0]+mu*defects[1]),e=f2*f2*se*defects[3],p=f2*f2*se*defects[4])
    transported_inlet={name:value*s.exp(-rate*y) for name,value,rate in (
      ('m',inlet_from_actual_defects['m'],1),('h',inlet_from_actual_defects['h'],s.Rational(3,2)),
      ('k',inlet_from_actual_defects['k'],s.Rational(3,2)),('e',inlet_from_actual_defects['e'],1),
      ('p',inlet_from_actual_defects['p'],0))}
    terminal_values=asts.evaluate(scalar_program,dict(c=SimpleNamespace(mpf=s.Rational),
      self=SimpleNamespace(f2=f2),previous=transported_inlet,primitives=primitives,x=s.exp(y)))
    terminal_factored=dict(m=f2*s.exp(-y)*sa*residual[0],
      h=f2*s.exp(-s.Rational(3,2)*y)*se*residual[2],
      k=f2*f2*s.exp(-s.Rational(3,2)*y)*sa*(residual[0]+mu*residual[1]),
      e=f2*f2*s.exp(-y)*se*residual[3],p=f2*f2*se*residual[4])
    for name in terminal_values:
        zero('actual_callable_exit_history_factors_same_checked_source_root_'+name,
          terminal_values[name],terminal_factored[name])
    Ahat=s.Symbol('same_A_over_Pstar');bump=[[s.Symbol('g'+str(i)+'_row'+str(j)) for j in range(5)] for i in range(3)]
    env=dict(c=SimpleNamespace(mpf=s.Rational),Ahat=Ahat,se=mu/N,sa=s.sqrt(mu)/N,h=h,rows=bump)
    du=asts.evaluate(asts.expression('current_O3_repaired_histories','history','du'),env)
    uz=asts.evaluate(asts.expression('current_O3_repaired_histories','history','Vnew'),env)
    for j in range(5):
        zero('actual_same_constant_controls_theta_source_row_'+str(j),du[j],Ahat*mu/N*sum(h[i+2]*bump[i][j] for i in range(3)))
        zero('actual_same_constant_controls_axial_source_row_'+str(j),uz[j],Ahat*s.sqrt(mu)/N*(h[0]*bump[0][j]+h[1]*bump[2][j]))
    F,G,Ebase=s.symbols('same_F same_G actual_base_power',real=True)
    # Differentiate the actual normalized primitive changes using their
    # own source FTC densities, not differences of numerical endpoint boxes.
    for name,factor,rate,density in (('m',f2,1,f2*G),
      ('h',f2,s.Rational(3,2),f2*F),
      ('k',f2*f2,s.Rational(3,2),f2*f2*(Ebase+F)*G),
      ('e',f2*f2,1,f2*f2*(G*G-Ebase*F-F*F/2))):
        primitive=s.Function('same_FTC_primitive_'+name)(y)
        value=factor*s.exp(-rate*y)*primitive
        derived=s.diff(value,y).subs(s.diff(primitive,y),s.exp(rate*y)*density/factor)
        zero('actual_repair_normalized_transport_ODE_'+name,derived,density-rate*value)
    pressure_primitive=s.Function('same_Cp_partial')(y)
    zero('actual_repair_pressure_FTC',s.diff(f2*f2*pressure_primitive,y).subs(
      s.diff(pressure_primitive,y),Ebase*F+F*F/2),f2*f2*(Ebase*F+F*F/2))
    return dict(identities=checks,passed=True,input_hashes=asts.hashes,
      same_implicit_controls_used_in_profiles_partial_primitives_and_exit=True,
      checked_parent_D_source_identity_required='actual_D_is_same_signed_history_difference_after_transport_and_scaling',
      exact_exit_reduction_requires_all_full_support_weights_and_checked_defining_equation=True,
      original_incoming_histories_not_zeroed=True,
      primitive_partial_definitions='single int beta(w)*exp(p*(c+w))dw; product int beta(w)^2*exp((p-1)*(c+w))dw; D int beta(w)*(exp(-mu*(c+w))-1)/mu dw',
      fixed_ordinary_logR_and_axial_orders=(4,5),
      physical_factors='R=R0*x; same A=Pstar*Ua*C*f2, f2=exp(-1-3mu/2)',
      outside_repair_five_delta_histories_follow_same_homogeneous_ODE_after_exact_exit=True,
      complete_modified_tensor_cone_physical_NS_remain_open=True)
