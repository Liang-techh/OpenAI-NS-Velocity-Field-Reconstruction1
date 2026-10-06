"""Original leading two-vector cone and covariance algebra, not a wave solve.

OpenAI (4.11), (4.20)-(4.23), (7.1), (7.24), (7.28), (7.31);
Lei--Ren (3.18), (3.21). Strict signed tests precede any square root.
"""
import ast
import functools
import sympy as s
from lei_ren_part1_paper_compliant_current_O3_transition_background_tensor import SourceAST
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

PAPERS={
 'openai':dict(title='Finite Time Blowup for Navier-Stokes',pages=166,
  url='https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf',
  sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
  references={'profile_shear':'pp.26-27, (4.8)-(4.11)',
   'cone':'pp.30-32, (4.20)-(4.23), Lemma 4.5',
   'uniform_edges':'p.33, Theorem 4.6(iii)-(iv), (4.26)-(4.27)',
   'shear_repair':'pp.39-40, (4.38)-(4.40); Appendix C',
   'frame':'p.74, (7.1)', 'actual_covariance':'pp.82-84, Proposition 7.5, (7.24)-(7.30)',
   'signed_linear_lift':'p.84, Proposition 7.6, (7.31)-(7.34)'}),
 'lei_ren':dict(title='Finite-Time Blowup for Navier-Stokes with Smooth Forcing, Part I',
  supplied_filename='2609.35406v2.pdf',pages=245,
  sha256='8396b998dcf737cd6a7c12ff1019c6b2fd308a7e6f0907c0c6647a9b2ec64e8d',
  references={'scales':'pp.11,18, (2.2)-(2.5), (3.1)-(3.3)',
   'shear':'p.21, (3.18)', 'cone':'p.22, (3.21)-(3.23)'}),
 'duraiswami':dict(title='Self-similar swirl between contracting porous walls: the GD1998 exact Navier-Stokes solution revisited',
  supplied_filename='2609.17642v1.pdf',pages=31,
  sha256='a38af23a5873e7671ecf0bcd7f3d4c45851c67689a7551ae0e51203ba166fe95',
  role='Comparison/numerical-method reference; not the original OpenAI paper or cone admission')}


@functools.lru_cache(maxsize=1)
def original_cone_theorem():
    F,a=s.symbols('F a',positive=True);b,tt,tz=s.symbols('bs Ttheta Tz',real=True)
    St,Sz=-F*a,F*b;ts=-b/a;vs=a+b*b/a
    dot=tt*St+tz*Sz;cross=-tt*Sz+tz*St
    D=tt+ts*tz;J=tz-ts*tt
    proofs={}
    def zero(name,expr):
        if s.cancel(s.expand(expr))!=0:raise ValueError('Original cone identity differs: '+name)
        proofs[name]=True
    zero('Lei_kappa_equals_OpenAI_vs',-(St**2+Sz**2)/(F*St)-vs)
    zero('negative_Lei_dot_equals_Fa_OpenAI_direction',-dot-F*a*D)
    zero('Lei_cross_equals_minus_Fa_OpenAI_cross',cross+F*a*J)
    zero('same_signed_quadratic_cone',2*dot**2-(vs-2)*cross**2-(F*a)**2*(2*D**2-(vs-2)*J**2))
    I1,I2=s.symbols('Itheta Iz',real=True);ps=(I1/F,I2/F)
    zero('Ttheta_equals_F_ps_minus_s',I1+St-F*(ps[0]-a))
    zero('Tz_equals_F_ps_minus_s',I2+Sz-F*(ps[1]+b))
    zero('stress_coordinate_Pc_minus_vs',D.subs({tt:I1+St,tz:I2+Sz})-F*(ps[0]+ts*ps[1]-vs))
    zero('stress_coordinate_Jc',J.subs({tt:I1+St,tz:I2+Sz})-F*(ps[1]-ts*ps[0]))
    nu,fac=s.symbols('nu lambda_factor',positive=True)
    zero('physical_common_positive_stress_scale',
        2*(nu*fac*D)**2-(vs-2)*(nu*fac*J)**2-(nu*fac)**2*(2*D**2-(vs-2)*J**2))
    h,delta,lam=s.symbols('h delta lambda',positive=True)
    zero('q_power_equals_lambda_power_exponent',2*(-(s.Rational(1,2)+delta/2)-s.Rational(1,2))-(-2-delta))
    # Bind the shear signs to the actual original repository stress source.
    asts=SourceAST();R,Ut=s.symbols('R Utheta',positive=True);Uty,Uzy=s.symbols('Utheta_y Uz_y',real=True)
    env=dict(a=Ut,ay=Uty,by=Uzy,root=s.sqrt(2*R),R=R)
    from pathlib import Path
    import hashlib
    path=Path(__file__).with_name('lei_ren_part1_paper_mp_stress.py')
    tree=ast.parse(path.read_text(encoding='utf8'));fn=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name=='evaluate_mp_stress')
    shear={}
    for target in ('Stheta','Sz'):
        nodes=[n.value for n in ast.walk(fn) if isinstance(n,ast.Assign) and any(ast.unparse(t)==target for t in n.targets)]
        if len(nodes)!=1 or not isinstance(nodes[0],ast.IfExp):raise ValueError('Original default shear branch changed')
        shear[target]=asts.evaluate(nodes[0].body,env)
    FF=Ut/s.sqrt(2*R)
    zero('actual_default_theta_shear_normalization',-shear['Stheta']/FF-(1-2*Uty/Ut))
    zero('actual_default_axial_shear_normalization',shear['Sz']/FF-2*Uzy/Ut)
    asts.hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
    # Positive reference coefficients and the actual two-column inverse.
    Ac,u,hp,hm=s.symbols('Ac u hplus hminus',positive=True)
    TN,TK=s.symbols('TN TK',real=True)
    Hp=s.Matrix([-hp*Ac,-hp*u]);Hm=s.Matrix([-hm*Ac,hm*u]);H=Hp.row_join(Hm)
    yp=(-TN/Ac-TK/u)/(2*hp);ym=(-TN/Ac+TK/u)/(2*hm)
    for k,x in enumerate(H*s.Matrix([yp,ym])-s.Matrix([TN,TK])):zero('reference_covariance_solve_'+str(k),x)
    H11,H12,H21,H22=s.symbols('H11 H12 H21 H22',real=True)
    HH=s.Matrix([[H11,H12],[H21,H22]]);det=HH.det()
    yy=s.Matrix([(H22*tt-H12*tz)/det,(-H21*tt+H11*tz)/det])
    for k,x in enumerate(HH*yy-s.Matrix([tt,tz])):zero('actual_covariance_inverse_identity_'+str(k),x)
    ap,am,eps=s.symbols('amplitude_plus amplitude_minus epsilon',positive=True)
    sp,sz=s.symbols('Sigma_theta Sigma_z',real=True)
    dd=HH.inv()*s.Matrix([sp,sz])/eps;da=s.Matrix([dd[0]/(2*ap),dd[1]/(2*am)])
    for k,x in enumerate(eps*HH*s.diag(2*ap,2*am)*da-s.Matrix([sp,sz])):zero('linear_cross_covariance_identity_'+str(k),x)
    return dict(paper_versions=PAPERS,identities=proofs,input_hashes=asts.hashes,
      variable_map=dict(X='R',eta='Z',q='lambda^2',h='delta/2',E='Utheta',U='Uz',
        F='Utheta/sqrt(2R)',s='(a,-bs)=-S/F',a='1-2*Utheta_y/Utheta',
        bs='2*Uz_y/Utheta',ps='(Itheta,Iz)/F',ts='-bs/a=Sz/Stheta',vs='a+bs^2/a=Lei kappa'),
      physical_normalization='u=sqrt(nu)*lambda^beta*U; p=nu*lambda^(-2-2delta)*P; T=nu*lambda^(-2-delta)*Tprofile; cone is homogeneous in T',
      exact_signed_cone='a>0; vs-2>0; D=Ttheta-(bs/a)*Tz>0; 2*D^2-(vs-2)*(Tz+(bs/a)*Ttheta)^2>0',
      covariance_scope='Reference columns are not actual pulses. Actual H requires (7.27) integration, finite error control, common slow boxes, edge weight, and smooth amplitudes.',
      higher_order_scope='The original positive target is order zero. Proposition 7.6 lifts later signed increments linearly about fixed positive amplitudes; it does not require each increment to lie in the leading cone.',
      zero_stress_scope='Strict interior cone excludes T=0. Edge direction/weight and exact stress-free core/exterior must be handled separately.',passed=True)


def cone_margins(c,a,bs,theta,axial,vs_minus2=None):
    """Signed interval test. A box containing zero is inconclusive, not failed."""
    a,bs,theta,axial=(c.mpf(v) for v in (a,bs,theta,axial))
    if endpoints(a)[0]<=0:return dict(status='shear_unresolved_or_nonpositive',a=a,admitted=False)
    km=c.mpf(vs_minus2) if vs_minus2 is not None else a-2+bs**2/a
    ts=-bs/a;D=theta+ts*axial;J=axial-ts*theta
    Q=2*D**2-km*J**2
    margins=dict(a=a,vs_minus2=km,signed_direction=D,signed_cross=J,signed_quadratic=Q)
    failed=[key for key in ('vs_minus2','signed_direction','signed_quadratic') if endpoints(margins[key])[1]<=0]
    uncertain=[key for key in ('vs_minus2','signed_direction','signed_quadratic') if endpoints(margins[key])[0]<=0<endpoints(margins[key])[1]]
    return dict(status='failed_directed_margin' if failed else 'inconclusive_enclosure' if uncertain else 'strict_two_vector_cone',
        margins=margins,failed=failed,inconclusive=uncertain,admitted=not failed and not uncertain,
        completed_tensor_or_global_wave_admission=False)


def actual_covariance_inverse(c,columns,target):
    """Algebraic interval inverse of supplied columns; no pulse provenance claim."""
    (x1,z1),(x2,z2)=((c.mpf(x),c.mpf(z)) for x,z in columns);tt,tz=(c.mpf(x) for x in target)
    det=x1*z2-x2*z1
    lo,hi=endpoints(det)
    if lo<=0<=hi:return dict(status='inconclusive_singular_column_box',actual_pulse_lift_constructed=False)
    y=((z2*tt-x2*tz)/det,(-z1*tt+x1*tz)/det)
    positive=all(endpoints(x)[0]>0 for x in y)
    return dict(status='positive_algebraic_coefficients' if positive else 'nonpositive_or_unresolved_coefficients',
        determinant=det,squared_amplitudes=y,amplitudes=tuple(c.sqrt(x) for x in y) if positive else None,
        positive_coefficients=positive,actual_pulse_lift_constructed=False,
        source_scope='Supplied matrix algebra only; no (7.27) pulse integration or smooth edge extension inferred')


def reference_covariance(c,a,bs,km,theta,axial):
    """Diagnostic reference decomposition from (7.28), e_sigma omitted explicitly."""
    tested=cone_margins(c,a,bs,theta,axial,km)
    if not tested['admitted']:return dict(status='reference_not_admitted',actual_pulse_lift_constructed=False)
    norm=c.sqrt(a*a+bs*bs);N=(-a/norm,bs/norm);K=(-N[1],N[0])
    TN=theta*N[0]+axial*N[1];TK=theta*K[0]+axial*K[1];c0=-c.sqrt(km/2)
    ratio=abs(c0*TK/TN);upper=endpoints(ratio)[1]
    if upper>=1:return dict(status='reference_frame_enclosure_inconclusive',actual_pulse_lift_constructed=False)
    # A diagnostic local u, not the uniform global phase parameter of (7.1).
    rho=(c.mpf(upper)+1)/2;u=rho/c.sqrt(1-rho*rho);Ac=-c0*c.sqrt(1+u*u)
    columns=tuple(tuple(-Ac*N[j]-sign*u*K[j] for j in (0,1)) for sign in (1,-1))
    solved=actual_covariance_inverse(c,columns,(theta,axial))
    # H's interval columns lose their shared frame/Ac correlation. Retain
    # that inverse diagnostic, but bound reference coefficients using the
    # exact ratio identity |TK/u|/(-TN/Ac) = ratio/rho.
    R=-TN/Ac;relative=ratio/rho
    lower=endpoints(R*(1-relative)/2)[0];upper_y=endpoints(R*(1+relative)/2)[1]
    positive=lower>0
    y_bounds=(c.mpf([lower,upper_y]),)*2
    return dict(status='positive_reference_coefficients_by_correlated_cone' if positive else 'reference_coefficient_enclosure_inconclusive',
        N=N,K=K,c0=c0,local_diagnostic_u=u,
        reference_columns=columns,reference_inverse=solved,cone_ratio=ratio,
        correlated_reference_squared_amplitude_bounds=y_bounds,
        exact_reference_coefficient_identity='y_plus=(-TN/Ac-TK/u)/2; y_minus=(-TN/Ac+TK/u)/2; h_plus=h_minus=1',
        strict_reference_positive_coefficients=positive,
        actual_pulse_lift_constructed=False,omitted_covariance_column_errors=True,
        diagnostic_local_phase_not_uniform_global_choice=True)


def linear_covariance_change(c,columns,positive_squared_amplitudes,Sigma,epsilon):
    """(7.31) algebra for a signed increment; no constructed pulse/edge claim."""
    epsilon=c.mpf(epsilon);y=tuple(c.mpf(x) for x in positive_squared_amplitudes)
    if endpoints(epsilon)[0]<=0 or any(endpoints(x)[0]<=0 for x in y):
        raise ValueError('Fixed positive epsilon and base squared amplitudes required')
    inverse=actual_covariance_inverse(c,columns,tuple(c.mpf(x)/epsilon for x in Sigma))
    if 'squared_amplitudes' not in inverse:raise ValueError('Invertible covariance column box required')
    d=inverse['squared_amplitudes'];amplitudes=tuple(c.sqrt(x) for x in y)
    return dict(signed_amplitude_changes=tuple(dj/(2*aj) for dj,aj in zip(d,amplitudes)),
        fixed_base_amplitudes=amplitudes,epsilon=epsilon,
        negative_increment_permitted=True,actual_linear_pulse_lift_constructed=False)
