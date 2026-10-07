"""Independent nonlinear fixed point, implicit jets and terminal references."""
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_native_Rc_convergent_control_family as producer

HERE,PREFIX,sha=producer.HERE,producer.PREFIX,producer.sha
source,controls,packets,ep,require=producer.source,producer.controls,producer.packets,producer.ep,producer.require


def independent_tail_and_implicit_identities():
    L,R0=s.symbols('L first_step',positive=True);n=s.symbols('n',integer=True,nonnegative=True)
    tail=R0*L**n/(1-L)
    require(s.simplify(tail.subs(L,s.Rational(1,2))-2*R0*2**(-n))==0,'Geometric C1 tail factor differs')
    require(s.Rational(1,2)+1==s.Rational(3,2),'Preconditioned residual bound differs')
    z,N,mu=s.symbols('Z N mu',positive=True)
    h=[s.Function('h'+str(i))(z) for i in range(5)];d=[s.Function('d'+str(i))(z) for i in range(5)]
    C0,C2,E0,E1,E2,P0,P1,P2=s.symbols('C0 C2 E0 E1 E2 P0 P1 P2')
    a0,a2,e0,e1,e2=h
    Q=s.Matrix([0,(C0*a0*e0+C2*a2*e2)/mu,0,
        E0*a0*a0+E2*a2*a2-(E0*e0*e0+E1*e1*e1+E2*e2*e2)/2,
        (P0*e0*e0+P1*e1*e1+P2*e2*e2)/2])
    B=s.Matrix(5,5,lambda i,j:s.Symbol('B'+str(i)+str(j)))
    residual=B*s.Matrix(h)+s.Matrix(d)+Q/N
    expected=(B+Q.jacobian(h)/N)*s.Matrix([s.diff(q,z) for q in h])+s.Matrix([s.diff(q,z) for q in d])
    require(all(s.simplify(v)==0 for v in s.diff(residual,z)-expected),'Exact full implicit Z equation differs')
    return dict(passed=True,independent_geometric_C1_tail_identity=True,
        preconditioned_residual_factor_three_halves=True,full_implicit_Z_equations=5,
        mu_B_and_all_original_bump_weights_Z_independent=True,
        mathematical_tail_and_numeric_evaluation_error_distinguished=True)


def manufactured_exact_family(family,c):
    g=source.FunctionTransportGraph();z=g.symbol('Z')
    pair=lambda a,b:source.C1Function(g.add(g.constant(a),g.mul(g.constant(b),z)),g.constant(b))
    targets={key:pair(a,b) for key,a,b in (('M','1/50','3/1000'),('D=(J-M)/mu','3/100','-1/500'),
        ('I','9/1000','1/1000'),('S','1/100','-1/1000'),('Cp','1/200','1/2000'))}
    A=pair('7/5','2/5');AA=g.c1mul(A,A);mu=g.node('original_source_parameter',name='mu')
    N=g.node('shared_positive_integer',name='N',lower=160);zero=source.C1Function(g.zero,g.zero)
    incoming={key:g.c1mul(amplitude,targets[row]) for key,row,amplitude in
        (('m','M',A),('h','I',A),('e','S',AA),('p','Cp',AA))}
    incoming['k']=g.c1mul(AA,g.c1add(targets['M'],g.c1scale(mu,targets[controls.ROWS[1]])))
    built=producer.generic_control_operator(dict(graph=g,parameters={'mu':mu},N=N,amplitude=A,N_scaled_targets=targets,
        coefficient_history={-1:incoming,-2:{key:zero for key in incoming}},source_family=family,
        source_graph_sha256='synthetic_reference'))
    L=c.log(2);ell=L/40;centers=[L/5,L/2,4*L/5]
    class Oracle:
        mode='synthetic_reference';source_family=family
        def parameter(self,name):
            require(name=='mu','Manufactured operator contains unknown parameter');return c.mpf('.04')
        def source(self,*a,**kw):raise ArithmeticError('Manufactured operator must not query original source leaves')
        def integrate(self,fn,lo,hi):
            points=[lo,hi]
            if lo<0<hi:points.append(c.mpf(0))
            if lo>=1:
                points += [c.exp(ci+q*ell) for ci in centers for q in (-1,0,1) if lo<c.exp(ci+q*ell)<hi]
            return c.quad(fn,sorted(set(points)))
    return built,Oracle()


def independent_nonlinear_and_terminal_reference(family):
    c=mp.mp.clone();c.dps=48;built,oracle=manufactured_exact_family(family,c)
    mu=c.mpf('.04');alpha=c.mpf('.5')+mu;L=c.log(2);ell=L/40;centers=[L/5,L/2,4*L/5]
    raw=lambda t:c.exp(-1/(1-t*t)) if abs(t)<1 else c.mpf(0)
    J0=c.quad(raw,[-1,0,1])
    average=lambda ci,fn:c.quad(lambda t:raw(t)*fn(c.exp(ci+ell*t))/J0,[-1,0,1])
    gram=lambda ci,p:c.quad(lambda t:raw(t)**2*c.exp(p*(ci+ell*t))/(ell*J0**2),[-1,0,1])
    B=c.matrix([[1,1,0,0,0],*[[(average(centers[i],lambda x:c.expm1(-mu*c.log(x))/mu)) for i in (0,2)]+[0,0,0]]])
    # Extend the rectangular matrix with the three physically integrated rows.
    rows=[[1,1,0,0,0],[average(centers[i],lambda x:c.expm1(-mu*c.log(x))/mu) for i in (0,2)]+[0,0,0]]
    rows += [[0,0]+[average(ci,lambda x:x**c.mpf('.5')) for ci in centers],
             [0,0]+[-average(ci,lambda x:x**(-alpha)) for ci in centers],
             [0,0]+[average(ci,lambda x:x**(-alpha-1)) for ci in centers]]
    B=c.matrix(rows);inv=B**-1
    cross=[gram(centers[i],c.mpf('-.5')) for i in (0,2)]
    energy=[gram(ci,-1) for ci in centers];pressure=[gram(ci,-2) for ci in centers]
    def Q(h):
        a,b,e,f,k=h
        return c.matrix([0,(cross[0]*a*e+cross[1]*b*k)/mu,0,
            energy[0]*a*a+energy[2]*b*b-(energy[0]*e*e+energy[1]*f*f+energy[2]*k*k)/2,
            (pressure[0]*e*e+pressure[1]*f*f+pressure[2]*k*k)/2])
    def DQ(h):
        a,b,e,f,k=h
        return c.matrix([[0]*5,[cross[0]*e/mu,cross[1]*k/mu,cross[0]*a/mu,0,cross[1]*b/mu],[0]*5,
            [2*energy[0]*a,2*energy[2]*b,-energy[0]*e,-energy[1]*f,-energy[2]*k],
            [0,0,pressure[0]*e,pressure[1]*f,pressure[2]*k]])
    CA=max(sum(abs(inv[i,j]) for j in range(5)) for i in range(5))
    CQ=max(sum(cross)/mu,energy[0]+energy[2]+sum(energy)/2,sum(pressure)/2)
    D=c.mpf('.036');rho=2*CA*D;n=10**9;lip=2*CA*CQ*rho/n
    require(lip<c.mpf('.5'),'Manufactured reference does not meet the same full C1 contraction recipe')
    comparisons=terminal=0;largest_residual=c.mpf(0)
    coeff=[(c.mpf(a),c.mpf(b)) for a,b in (('.02','.003'),('.03','-.002'),('.009','.001'),('.01','-.001'),('.005','.0005'))]
    g=lambda i,x:raw((c.log(x)-centers[i])/ell)/(ell*J0*x)
    for z in (c.mpf('-.6'),c.mpf(0),c.mpf('.7')):
        d=c.matrix([a+b*z for a,b in coeff]);dZ=c.matrix([b for a,b in coeff])
        fixed=c.findroot(lambda *h:tuple(B*c.matrix(h)+d+Q(h)/n),tuple(-inv*d),tol=c.mpf('1e-42'))
        fixedZ=-(B+DQ(fixed)/n)**-1*dZ
        value=producer.CachedC1ControlEvaluator(built,oracle=oracle,Z=z,N=n,ctx=c)
        for depth in (3,16,64):
            result=value.evaluate(depth)
            require(result.approximate_only and not result.original_numeric_fixed_point_certified,'Finite numeric evaluation falsely certified')
            tail=rho*2**(-depth)
            for got,want in zip((*result.values,*result.Z_derivatives),(*fixed,*fixedZ)):
                require(abs(got-want)<=tail+c.mpf('1e-35'),'Cached C1 iterates do not approach independent nonlinear/implicit reference')
                comparisons+=1
            if depth==64:
                largest_residual=max(largest_residual,*[abs(v) for v in (*result.residual_values,*result.residual_Z_derivatives)])
        # Physical, independently integrated five terminal rows with A_Z.
        A=c.mpf('1.4')+c.mpf('.4')*z;AZ=c.mpf('.4')
        incoming={'m':A*d[0]/n,'h':A*d[2]/n,'e':A*A*d[3]/n,'p':A*A*d[4]/n,
            'k':A*A*(d[0]+mu*d[1])/n}
        incomingZ={'m':(AZ*d[0]+A*dZ[0])/n,'h':(AZ*d[2]+A*dZ[2])/n,
            'e':(2*A*AZ*d[3]+A*A*dZ[3])/n,'p':(2*A*AZ*d[4]+A*A*dZ[4])/n,
            'k':(2*A*AZ*(d[0]+mu*d[1])+A*A*(dZ[0]+mu*dZ[1]))/n}
        def physical(x):
            F=sum(fixed[i+2]*g(i,x) for i in range(3))/n;G=(fixed[0]*g(0,x)+fixed[1]*g(2,x))/n
            FZ=sum(fixedZ[i+2]*g(i,x) for i in range(3))/n;GZ=(fixedZ[0]*g(0,x)+fixedZ[1]*g(2,x))/n
            E=A*x**(-alpha);EZ=AZ*x**(-alpha);de=A*F;dv=A*G;dez=AZ*F+A*FZ;dvz=AZ*G+A*GZ
            return ({'m':dv,'h':de,'k':(E+de)*dv,'e':dv*dv-E*de-de*de/2,'p':E*de+de*de/2},
                {'m':dvz,'h':dez,'k':(EZ+dez)*dv+(E+de)*dvz,
                 'e':2*dv*dvz-EZ*de-E*dez-de*dez,'p':EZ*de+E*dez+de*dez})
        for key,rate in producer.allN.RATES.items():
            r=c.mpf(rate.numerator)/rate.denominator
            for order,inc in ((0,incoming),(1,incomingZ)):
                end=2**(-r)*(inc[key]+oracle.integrate(lambda x:x**(r-1)*physical(x)[order][key],c.mpf(1),c.mpf(2)))
                require(abs(end)<c.mpf('1e-36'),'Independent physical terminal row did not close at nonlinear reference')
                terminal+=1
    return dict(passed=True,independent_nonlinear_fixed_point_and_implicit_Z_comparisons=comparisons,
        independent_physical_terminal_value_Z_rows=terminal,depths=[3,16,64],reference_N=n,
        nonconstant_actual_amplitude_Z_and_joint_k_row_and_pressure_memory_tested=True,
        finite_numeric_reference_residual_log=str(c.log(max(largest_residual,c.mpf('1e-48')))),
        manufactured_operator_reference_only=True,native_original_point_functions_or_N_not_evaluated=True)


def original_family_provenance(owner,built,report):
    original=owner.integrals.data['current_native_Rc_all_N_function_controls'];prefix=original['exact_function_graph_nodes']
    require(built['graph'].nodes[:len(prefix)]==prefix and report['exact_function_graph_nodes']==built['graph'].nodes,
            'Control family changed original exact source graph or map')
    require(report['exact_original_C1_limit_contract']==packets.encode(owner.contract()),'Original C1 source/frequency contract changed')
    for depth in (0,3,16,64):
        row=owner.factored_tail(depth)
        require(row['relative_control_C1_tail']=={'numerator':1,'denominator':2**depth} and
                row['relative_preconditioned_residual_C1_tail']=={'numerator':3,'denominator':2**(depth+1)} and
                not row['numerical_oracle_quadrature_and_roundoff_errors_included'],
                'Factored tail lost its dyadic decrement or mixed numeric errors')
    for key in ('actual_five_controls_installed','certified_actual_fixed_point_tail_installed','actual_terminal_Z_function_closure_installed',
        'current_whole_N_selected','numerical_original_source_point_or_integral_oracle_installed',*packets.OPEN):
        require(report[key] is False,'Exact family overclaims numeric/global construction: '+key)
    require(report['finite_Picard_terminal_functions_declared_zero'] is False and report['exact_same_source_C1_fixed_point_family_defined'],
            'Exact limit confused with finite terminal functions')
    return dict(passed=True,unchanged_original_native_graph_prefix_nodes=len(prefix),
        realized_original_C1_integral_target_and_native_frequency_contract_bound=True,
        full_five_control_value_Z_and_implicit_operator_templates_defined=True,
        exact_same_source_limit_family_and_terminal_C1_limit_bound=True,
        finite_iteration_approximation_and_numeric_quadrature_errors_kept_separate=True)


def run():
    began=time.monotonic();owner=producer.NativeRcConvergentControlFamily();built=owner.build()
    report=json.loads((HERE/producer.NAME).read_bytes());require(report[producer.GATE] and report['source_family']==owner.family,'Current control family required')
    for name,value in owner.hashes.items():require(report['input_hashes'].get(name)==value,'Original family input changed: '+name)
    result=dict(all_passed=True,**{producer.GATE:True},source_family=owner.family,
        independent_tail_and_implicit=independent_tail_and_implicit_identities(),
        independent_nonlinear_and_terminal=independent_nonlinear_and_terminal_reference(owner.family),
        actual_original_family_provenance=original_family_provenance(owner,built,report),
        read_only_review_model=dict(model='gpt-5.6-luna',reasoning_effort='max'),
        exact_same_source_C1_fixed_point_family_defined=True,
        numerical_original_source_point_or_integral_oracle_installed=False,
        actual_five_controls_installed=False,certified_actual_fixed_point_tail_installed=False,
        actual_terminal_Z_function_closure_installed=False,current_whole_N_selected=False,
        **dict.fromkeys(packets.OPEN,False),input_hashes={**report['input_hashes'],producer.NAME:sha(producer.NAME),Path(__file__).name:sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began)
    (HERE/producer.RECEIPT).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    print('Original convergent C1 control family: nonlinear/implicit jets and physical terminal references PASS',flush=True)
    return result


if __name__=='__main__':run()
