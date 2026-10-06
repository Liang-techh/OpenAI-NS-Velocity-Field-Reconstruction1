"""New quiet-power five-bump map with correlated axial divided difference.

No legacy repair matrix is reused. Fixed log-translated compact bumps and
the actual source defect correlation give a uniformly invertible scaled
map, even though the unscaled first two rows nearly coincide.
"""
from types import SimpleNamespace
import ast
import hashlib
from pathlib import Path
import math
import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_current_O3_modulated_histories_operator import MESH
from lei_ren_part1_paper_compliant_current_O3_finite_frequency_profiles import (
    SourceAST,cutoff_rows,endpoints,source_precision)
from lei_ren_part1_paper_compliant_outer_pulse_map import raw_beta

CENTERS=('1/5','1/2','4/5')
RADIUS='1/40'
ROWS=('M','(J-M)/mu','I','S','Cp')
CONTROLS=('axial0','axial2','swirl0','swirl1','swirl2')

def exact_repair_theorem():
    checks={};asts=SourceAST()
    def zero(name,a,b):
        if s.simplify(s.expand(a-b))!=0:raise ArithmeticError('Independent repair identity: '+name)
        checks[name]=True
    mu,N=s.symbols('mu N',positive=True)
    t,q,w=s.symbols('same_logR_offset quiet_q local_logx',real=True)
    ctx=SimpleNamespace(exp=s.exp,mpf=s.Rational)
    f=asts.evaluate(asts.expression('pre_pulse_mixed_C4','power','f',
      wanted="c.exp((-c.mpf('.5')-mu)*t)"),dict(c=ctx,mu=mu,t=q))
    zero('actual_quiet_power_exponent',f,s.exp(-(s.Rational(1,2)+mu)*q))
    same_endpoint_J=asts.evaluate(asts.expression('outer_buffer','transition_kernels','J',
      wanted="c.mpf('.5')"),dict(c=ctx))
    zero('actual_reflected_transition_unit_J_source',same_endpoint_J,s.Rational(1,2))
    ft=asts.evaluate(asts.expression('pre_pulse_mixed_C4','slope_mu','factor',
      wanted="c.exp(-t/2-mu*K['J'])"),dict(c=ctx,t=s.Integer(1),mu=mu,K=dict(J=same_endpoint_J)))
    f2=s.exp(-1-s.Rational(3,2)*mu)
    zero('actual_repair_inlet_f2',ft*f.subs(q,1),f2)
    parent=s.Symbol('same_power_inlet_Utheta_over_Pstar')
    power_u=asts.evaluate(asts.expression('pre_pulse_mixed_C4','power','u',wanted='u1*f'),dict(u1=parent,f=f))
    zero('actual_power_full_amplitude_function',power_u,parent*s.exp(-(s.Rational(1,2)+mu)*q))
    bump_lower=asts.expression('outer_pulse_map','raw_beta','lower')
    expected=ast.parse('c.exp(-1/(1-c.mpf(far)**2)) if far<1 else c.mpf(0)',mode='eval').body
    if ast.dump(bump_lower)!=ast.dump(expected):raise ValueError('Actual raw flat bump source changed')
    checks['actual_same_raw_flat_bump_program_bound']=True
    # The old normalization is reusable, but none of its repair matrix is.
    # Bind its exact unweighted raw-integral definition and normalized ratio.
    normalizer_path=Path(__file__).parent/'lei_ren_part1_paper_bump_integral_enclosures.py'
    normalizer_tree=ast.parse(normalizer_path.read_text(encoding='utf8'))
    normalizer_class=next(n for n in normalizer_tree.body if isinstance(n,ast.ClassDef) and n.name=='PaperBumpIntegralEnclosures')
    normalizer_init=next(n for n in normalizer_class.body if isinstance(n,ast.FunctionDef) and n.name=='__init__')
    wanted_normalizer=ast.parse('self._normalization, normalization_record = self._raw_integral(Fraction(5,4), Fraction(0), 1)').body[0]
    if not any(ast.dump(n)==ast.dump(wanted_normalizer) for n in ast.walk(normalizer_init)):
        raise ValueError('Same exact unweighted raw bump normalization definition changed')
    normalized_weight=next(n for n in normalizer_class.body if isinstance(n,ast.FunctionDef) and n.name=='weight')
    if not any(isinstance(n,ast.If) and ast.unparse(n.test)=='multiplicity == 1'
      and any(isinstance(v,ast.Return) and ast.dump(v.value)==ast.dump(ast.parse('raw/self._normalization',mode='eval').body)
        for v in n.body) for n in ast.walk(normalized_weight)):
        raise ValueError('Same normalized raw-weight source changed')
    asts.hashes[normalizer_path.name]=hashlib.sha256(normalizer_path.read_bytes()).hexdigest()
    checks['actual_exact_unweighted_raw_bump_normalization_program_bound']=True
    raw_mass,ell,raw=s.symbols('same_exact_positive_raw_mass ell same_raw_beta',positive=True)
    zero('actual_log_bump_Jacobian_cancels_radius',raw/(ell*raw_mass)*ell,raw/raw_mass)
    zero('actual_normalized_log_bump_H0_equals_one',raw_mass/raw_mass,s.Integer(1))
    J,Abar,beta=s.symbols('same_cutoff_J actual_Abar actual_beta',real=True)
    im=s.exp(-1+s.Rational(3,2)*mu)*s.exp(t/2-mu*J)*beta/N
    ij=s.exp(-1+3*mu)*s.exp(t/2-2*mu*J)*beta*s.exp(mu*Abar/N)/N
    exponent=mu*(s.Rational(3,2)-J+Abar/N)
    zero('actual_axial_defects_common_integrand_ratio',ij,im*s.exp(exponent))
    zero('actual_axial_defect_divided_difference_before_enclosure',
         (ij-im)/mu,im*(s.exp(exponent)-1)/mu)
    # Replay the defining existing scalar m/k density and common transport
    # programs before taking a divided difference of their exact integrals.
    scalar_source=asts.method('current_O3_modulated_histories_operator','scalar_integrand_enclosure')
    scalar_return=next(n for n in ast.walk(scalar_source) if isinstance(n,ast.Return))
    returned={node.arg:node.value for node in scalar_return.value.keywords}
    env=dict(c=ctx,t=t,mu=mu,J=J,beta=beta,N=N,G=s.exp(mu*Abar/N))
    source_m=asts.evaluate(returned['m'],env);source_k=asts.evaluate(returned['k'],env)
    cumulative=asts.method('current_O3_modulated_histories_operator','cumulative_scalar_enclosures')
    transport=next(n for n in ast.walk(cumulative) if isinstance(n,ast.AugAssign)
      and ast.unparse(n.target)=='result[name]'
      and ast.dump(n.value)==ast.dump(ast.parse('c.exp(-c.mpf(str(rate))*t)',mode='eval').body))
    # SHAPES is a module-level source definition, not a fitted table.
    history_tree=ast.parse((Path(__file__).parent/'lei_ren_part1_paper_compliant_current_O3_modulated_histories_operator.py').read_text(encoding='utf8'))
    shapes_node=next(n.value for n in history_tree.body if isinstance(n,ast.Assign)
      and any(isinstance(k,ast.Name) and k.id=='SHAPES' for k in n.targets))
    expected_shapes=ast.parse("{'m':('C',1,1),'h':('C',1,'3/2'),'k':('C^2',2,'3/2'),'e':('C^2',2,1),'p':('C^2',2,0)}",mode='eval').body
    if ast.dump(shapes_node)!=ast.dump(expected_shapes):raise ValueError('Actual history source shapes/transport rates changed')
    tm=asts.evaluate(transport.value,dict(c=ctx,rate=1,t=s.Integer(2)))
    tk=asts.evaluate(transport.value,dict(c=ctx,rate=s.Rational(3,2),t=s.Integer(2)))
    zero('actual_history_m_repair_inlet_transport',tm,s.exp(-2))
    zero('actual_history_k_repair_inlet_transport',tk,s.exp(-3))
    zero('actual_scaled_M_density_bound_to_existing_history_source',tm*source_m/f2,im)
    zero('actual_scaled_J_density_bound_to_existing_history_source',tk*source_k/f2**2,ij)
    for stem in ('current_O3_modulated_histories_operator','current_O3_independent_repair_operator'):
        fn='scalar_integrand_enclosure' if stem.endswith('modulated_histories_operator') else 'correlated_scaled_defects'
        asts.expression(stem,fn,'chi',wanted='cutoff_rows(c,t)[0]')
        asts.expression(stem,fn,'phase',wanted='N*t')
    symbolic_trig=SimpleNamespace(pi=s.pi,sin=s.sin,cos=s.cos,sqrt=s.sqrt)
    phase_env=dict(c=symbolic_trig,chi=s.Symbol('chi'),phase=s.Symbol('phase'),mu=mu)
    actual_A=asts.evaluate(asts.expression('current_O3_modulated_histories_operator','scalar_integrand_enclosure','A',
      wanted='-mu*chi**2*c.sin(4*c.pi*phase)/(8*c.pi)'),phase_env)
    actual_Abar=asts.evaluate(asts.expression('current_O3_independent_repair_operator','correlated_scaled_defects','Abar',
      wanted='-chi**2*c.sin(4*c.pi*phase)/(8*c.pi)'),phase_env)
    zero('actual_Abar_is_existing_history_A_divided_by_mu',actual_Abar,actual_A/mu)
    actual_beta=asts.evaluate(asts.expression('current_O3_modulated_histories_operator','scalar_integrand_enclosure','beta',
      wanted='-c.sqrt(mu)*chi*c.cos(2*c.pi*phase)/(2*c.pi)'),phase_env)
    scaled_beta=asts.evaluate(asts.expression('current_O3_independent_repair_operator','correlated_scaled_defects','beta_scaled',
      wanted='-chi*c.cos(2*c.pi*phase)/(2*c.pi)'),phase_env)
    zero('actual_scaled_beta_is_existing_history_beta_divided_by_sqrt_mu',scaled_beta,actual_beta/s.sqrt(mu))
    new_q=asts.evaluate(asts.expression('current_O3_independent_repair_operator','correlated_scaled_defects','q',
      wanted="c.mpf('1.5')-J+Abar/N"),dict(c=ctx,J=J,Abar=Abar,N=N))
    zero('actual_correlated_defect_q_from_same_source_quotient',new_q,s.Rational(3,2)-J+Abar/N)
    zero('actual_D_is_same_signed_history_difference_after_transport_and_scaling',
      N/s.sqrt(mu)*(tk*source_k/f2**2-tm*source_m/f2)/mu,
      s.exp(-1+s.Rational(3,2)*mu)*s.exp(t/2-mu*J)*(beta/s.sqrt(mu))*(s.exp(exponent)-1)/mu)
    checks['actual_correlation_uses_same_cutoff_phase_and_scalar_history_programs']=True
    Ua,C,P,R0,Bm,Bh,Bk,Be,Bp=s.symbols('Ua C Pstar R0 Bm Bh Bk Be Bp',positive=True)
    Amp=P*Ua*C*f2
    physical=(R0*P*Ua*C*Bm,s.sqrt(2)*R0**s.Rational(3,2)*P**2*Ua**2*C**2*Bk,
      s.sqrt(2)*R0**s.Rational(3,2)*P*Ua*C*Bh,R0*P**2*Ua**2*C**2*Be,
      P**2*Ua**2*C**2*Bp)
    units=(Amp*R0,s.sqrt(2)*Amp**2*R0**s.Rational(3,2),
      s.sqrt(2)*Amp*R0**s.Rational(3,2),Amp**2*R0,Amp**2)
    for name,value,unit,expected in zip(('M','J','I','S','Cp'),physical,units,
      (Bm/f2,Bk/f2**2,Bh/f2,Be/f2**2,Bp/f2**2)):
        zero('actual_full_Z_and_Pstar_cancel_'+name,value/unit,expected)
    x,F,G=s.symbols('x F G',positive=True);alpha=s.Rational(1,2)+mu
    Ebase=Amp*x**(-alpha);dE=Amp*F;dU=Amp*G;R=R0*x
    direct=(dU*R0,s.sqrt(2*R)*(Ebase+dE)*dU*R0,s.sqrt(2*R)*dE*R0,
      (dU**2-((Ebase+dE)**2-Ebase**2)/2)*R0,
      ((Ebase+dE)**2-Ebase**2)/(2*R)*R0)
    densities=(G,s.sqrt(x)*(x**(-alpha)+F)*G,s.sqrt(x)*F,
      G*G-x**(-alpha)*F-F*F/2,(x**(-alpha)*F+F*F/2)/x)
    for name,value,unit,density in zip(('M','J','I','S','Cp'),direct,units,densities):
        zero('actual_complete_physical_repair_density_'+name,value/unit,density)
    zero('actual_J_linear_weight',s.sqrt(x)*x**(-alpha),x**(-mu))
    zero('actual_S_linear_weight',x**(-alpha),x**(-s.Rational(1,2)-mu))
    zero('actual_Cp_linear_weight',x**(-alpha)/x,x**(-s.Rational(3,2)-mu))
    a,e=s.symbols('scaled_axial scaled_swirl',real=True)
    Fscaled=mu*e/N;Gscaled=s.sqrt(mu)*a/N
    zero('scaled_cross_quadratic_coefficient',Fscaled*Gscaled/(mu*s.sqrt(mu)/N),a*e/N)
    zero('scaled_kinetic_quadratic_coefficient',Gscaled**2/(mu/N),a*a/N)
    zero('scaled_swirl_square_quadratic_coefficient',Fscaled**2/(2*mu/N),mu*e*e/(2*N))
    p,center,ell=s.symbols('weight_power center ell',real=True)
    zero('log_translated_single_bump_weight',s.exp(p*(center+w)),s.exp(p*center)*s.exp(p*w))
    zero('log_translated_product_bump_weight',s.exp((p-1)*(center+w)),
      s.exp((p-1)*center)*s.exp((p-1)*w))
    H,r0,r2,dc=s.symbols('same_Hneg same_row0 same_row2 center_gap',real=True)
    zero('correlated_axial_divided_det',
      H*s.exp(-mu*s.Rational(1,5))*(s.exp(-mu*dc)-1)/mu,
      (H*s.exp(-mu*(s.Rational(1,5)+dc))-H*s.exp(-mu*s.Rational(1,5)))/mu)
    nodes=s.symbols('r0 r1 r2');V=s.Matrix([[1,r,r*r] for r in nodes])
    zero('equal_log_centers_swirl_Vandermonde_det',V.det(),
      (nodes[1]-nodes[0])*(nodes[2]-nodes[0])*(nodes[2]-nodes[1]))
    return dict(identities=checks,passed=True,input_hashes=asts.hashes,
      fixed_logx_centers=CENTERS,fixed_logx_radius=RADIUS,
      actual_bump='g_i(x)=beta(log(x)-c_i)/x; beta(w)=raw_beta(w/ell)/(ell*same_exact_normalization)',
      repair_band='q in[1,2]; x=exp(q-1); R0=Rw*exp(1)',
      source_swirl='E0=A(Z)*x^(-1/2-mu); A=Pstar*Ua*C*exp(-1-3mu/2)',
      source_axial='original Uz=0; original incoming moments are not zeroed',
      scaled_profiles='deltaE=A*mu/N*sum(e_i*g_i); deltaUz=A*sqrt(mu)/N*sum(a_i*g_i)',
      actual_unscaled_five_map=dict(M='integral G dx',J='integral sqrt(x)*(x^(-alpha)+F)*G dx',
        I='integral sqrt(x)*F dx',S='integral (G^2-x^(-alpha)*F-F^2/2) dx',
        Cp='integral (x^(-alpha)*F+F^2/2)/x dx'),
      scaled_row_order=ROWS,scaled_control_order=CONTROLS,
      correlated_defect_D='exp(-1+3mu/2)*integral exp(s/2-mu*J)*(beta/sqrt(mu))*q(s)*integral_0^1 exp(mu*r*q(s))dr ds; q=3/2-J+Abar/N',
      constants_are_Z_independent_after_actual_common_amplitude_normalization=True,
      complete_repaired_field_source_installed=False)

def upper(c,x):return c.mpf(endpoints(c.mpf(x))[1])

def exp_average(c,x):
    """Enclose the same exact integral_0^1 exp(r*x)dr without cancellation."""
    lo,hi=endpoints(x);return c.exp(c.mpf([min(mp.mpf(0),lo),max(mp.mpf(0),hi)]))

def bump_weights(c,mu,normalization,cells=256):
    ell=c.mpf(1)/40;alpha=c.mpf('.5')+mu
    powers=(c.mpf('.5'),-alpha,-alpha-1,-mu)
    H=[c.mpf(0) for _ in powers];K=[c.mpf(0) for _ in range(3)]
    D=[c.mpf(0),c.mpf(0)];centers=[c.mpf(1)/5,c.mpf(1)/2,c.mpf(4)/5]
    for i in range(cells):
        y=c.mpf([-1+mp.mpf(2)*i/cells,-1+mp.mpf(2)*(i+1)/cells]);w=ell*y
        raw=raw_beta(c,y);mass=raw*(c.mpf(2)/cells)/normalization
        square_mass=raw**2*(c.mpf(2)/cells)/(ell*normalization**2)
        for j,p in enumerate(powers):H[j]+=mass*c.exp(p*w)
        for j,p in enumerate((c.mpf('-.5'),c.mpf(-1),c.mpf(-2))):K[j]+=square_mass*c.exp(p*w)
        for j,center in enumerate((centers[0],centers[2])):
            v=center+w;D[j]+=mass*(-v)*exp_average(c,-mu*v)
    if any(endpoints(v)[0]<=0 for v in H+K):raise ArithmeticError('Positive common bump weights required')
    return dict(H=dict(zip(('I','S','Cp','axial'),H)),gram=dict(zip(('cross','energy','pressure'),K)),
      divided_axial_rows=D,centers=centers,radius=ell,cells=cells)

def inverse_and_bounds(c,mu,W):
    H=W['H'];ci=W['centers'];gap=ci[2]-ci[0];step=ci[1]-ci[0]
    detA=-gap*H['axial']*c.exp(-mu*ci[0])*exp_average(c,-mu*gap)
    D=W['divided_axial_rows']
    inverseA=[[D[1]/detA,-1/detA],[-D[0]/detA,1/detA]]
    powers=(c.mpf('.5'),-c.mpf('.5')-mu,-c.mpf('1.5')-mu)
    signs=(1,-1,1)
    E=[[signs[i]*H[k]*c.exp(powers[i]*x) for x in ci] for i,k in enumerate(('I','S','Cp'))]
    factors=[signs[i]*H[k]*c.exp(powers[i]*ci[0]) for i,k in enumerate(('I','S','Cp'))]
    nodes=[c.exp(p*step) for p in powers];detE=factors[0]*factors[1]*factors[2]
    for i,j in ((0,1),(0,2),(1,2)):
        difference=c.exp(powers[i]*step)*(powers[j]-powers[i])*step*exp_average(c,(powers[j]-powers[i])*step)
        detE*=difference
    if endpoints(detA)[1]>=0 or endpoints(detE)[0]<=0:raise ArithmeticError('New fixed axial/swirl determinants unresolved')
    inverseE=[]
    for i in range(3):
        row=[]
        for j in range(3):
            ii=[k for k in range(3) if k!=j];jj=[k for k in range(3) if k!=i]
            cofactor=(E[ii[0]][jj[0]]*E[ii[1]][jj[1]]-E[ii[0]][jj[1]]*E[ii[1]][jj[0]])*((-1)**(i+j))
            row.append(cofactor/detE)
        inverseE.append(row)
    B=[[c.mpf(1),c.mpf(1)]+[c.mpf(0)]*3,D+[c.mpf(0)]*3]
    B += [[c.mpf(0)]*2+row for row in E]
    inv=[[c.mpf(0)]*5 for _ in range(5)]
    for i in range(2):inv[i][:2]=inverseA[i]
    for i in range(3):inv[i+2][2:]=inverseE[i]
    norm=max(endpoints(sum((upper(c,abs(v)) for v in row),c.mpf(0)))[1] for row in inv)
    gram=W['gram']
    cross=[gram['cross']*c.exp(-ci[k]/2) for k in (0,2)]
    energy=[gram['energy']*c.exp(-x) for x in ci]
    pressure=[gram['pressure']*c.exp(-2*x) for x in ci]
    CQ=c.mpf(max(endpoints(value)[1] for value in
      (sum(cross),energy[0]+energy[2]+mu*sum(energy)/2,mu*sum(pressure)/2)))
    return dict(linear=B,inverse=inv,axial_divided_determinant=detA,swirl_determinant=detE,
      inverse_infinity_norm_upper=upper(c,norm),quadratic_norm_times_N_upper=CQ,
      cross_weights=cross,energy_weights=energy,pressure_weights=pressure,
      original_unscaled_axial_determinant_is_mu_times_divided_determinant=True)

def scaled_quadratic(c,mu,N,matrix,h):
    a0,a2,e0,e1,e2=h;zero=a0*0;energy=matrix['energy_weights'];pressure=matrix['pressure_weights']
    return [zero,(matrix['cross_weights'][0]*a0*e0+matrix['cross_weights'][1]*a2*e2)/N,
      zero,(energy[0]*a0*a0+energy[2]*a2*a2-mu*(energy[0]*e0*e0+energy[1]*e1*e1+energy[2]*e2*e2)/2)/N,
      mu*(pressure[0]*e0*e0+pressure[1]*e1*e1+pressure[2]*e2*e2)/(2*N)]

def source_defect_bounds(c,mu):
    L=c.mpf('2.5');a=1/(8*c.pi);b=1/(2*c.pi);q=c.mpf('1.5')+a
    M=L*c.exp(-c.mpf('.75')+c.mpf('1.5')*mu)*b
    return [upper(c,M),upper(c,M*q*c.exp(mu*q)),
      upper(c,L*c.exp(-c.mpf('1.5')+c.mpf('1.5')*mu)*a*c.exp(mu*a)),
      upper(c,L*c.exp(3*mu)*(b*b+a*c.exp(2*mu*a))),
      upper(c,L*c.exp(4+3*mu)*a*c.exp(2*mu*a))]

@source_precision
def correlated_scaled_defects(field):
    c=field.ctx;mu=field.mu;N=c.mpf(field.N);f2=c.exp(-1-c.mpf('1.5')*mu)
    source=field.scalars(2);dM=N/c.sqrt(mu)*source['m']/f2
    result=c.mpf(0)
    for left,right in zip(MESH[:-1],MESH[1:]):
        a=c.mpf(str(left));b=c.mpf(str(right))
        for i in range(field.cells):
            x=a+(b-a)*i/field.cells;y=a+(b-a)*(i+1)/field.cells
            t=c.mpf([endpoints(x)[0],endpoints(y)[1]])
            chi=cutoff_rows(c,t)[0];J=c.mpf(0) if endpoints(t)[1]<=0 else c.mpf([0,endpoints(t)[1]])
            phase=N*t;beta_scaled=-chi*c.cos(2*c.pi*phase)/(2*c.pi)
            Abar=-chi**2*c.sin(4*c.pi*phase)/(8*c.pi)
            q=c.mpf('1.5')-J+Abar/N
            result+=(y-x)*c.exp(t/2-mu*J)*beta_scaled*q*exp_average(c,mu*q)
    dD=c.exp(-1+c.mpf('1.5')*mu)*result
    return [dM,dD,N/mu*source['h']/f2,N/mu*source['e']/f2**2,N/mu*source['p']/f2**2]

def contraction_certificate(c,mu,N,matrix):
    bounds=source_defect_bounds(c,mu);D=c.mpf(max(endpoints(v)[1] for v in bounds));CA=matrix['inverse_infinity_norm_upper'];CQ=matrix['quadratic_norm_times_N_upper']
    radius=upper(c,2*CA*D);threshold=int(mp.ceil(endpoints(upper(c,4*CA*CQ*radius))[1]))+1
    if N<threshold:raise ValueError('Finite frequency below the independent repair-only bound')
    lipschitz=upper(c,2*CA*CQ*radius/N)
    image=upper(c,CA*D+CA*CQ*radius**2/N)
    if endpoints(lipschitz)[1]>=1 or endpoints(image)[1]>=endpoints(radius)[0]:
        raise ArithmeticError('Strict scaled contraction/inclusion failed')
    return dict(scaled_defect_abs_bounds_all_integer_N_at_least1=bounds,
      scaled_controls_ball_radius=radius,repair_only_sufficient_integer_N=threshold,
      actual_finite_N=N,strict_contraction_upper=lipschitz,strict_image_radius_upper=image,
      unique_exact_implicit_controls_exist=True,exact_controls_are_Z_independent=True,
      no_source_defect_midpoint_or_zero_substitution=True,common_cone_frequency_admitted=False)

def coefficient_enclosure(c,mu,N,matrix,defects,certificate,iterations=12):
    radius=certificate['scaled_controls_ball_radius'];r=endpoints(radius)[1]
    box=[c.mpf([-r,r]) for _ in range(5)];inverse=matrix['inverse']
    for _ in range(iterations):
        Q=scaled_quadratic(c,mu,N,matrix,box)
        mapped=[-sum((inverse[i][j]*(defects[j]+Q[j]) for j in range(5)),c.mpf(0)) for i in range(5)]
        next_box=[]
        for old,new in zip(box,mapped):
            lo=max(endpoints(old)[0],endpoints(new)[0]);hi=min(endpoints(old)[1],endpoints(new)[1])
            if lo>hi:raise ArithmeticError('Exact implicit repair enclosure is empty')
            next_box.append(c.mpf([lo,hi]))
        box=next_box
    return box
