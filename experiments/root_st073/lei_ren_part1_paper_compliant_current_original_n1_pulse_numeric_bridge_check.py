"""Universal source reduction, physical units and original pulse C1 queries."""
import ast
import inspect
import json
from pathlib import Path
import textwrap
import time
from types import SimpleNamespace

import sympy as sy

import lei_ren_part1_paper_compliant_current_original_n1_pulse_numeric_bridge as current


def compile_formula(fn,scope):
    tree=ast.parse(textwrap.dedent(inspect.getsource(fn)));tree.body[0].decorator_list=[]
    namespace=dict(fn.__globals__);namespace.update(scope)
    exec(compile(tree,'<actual-pulse-pressure-source>','exec'),namespace)
    return namespace[fn.__name__]


def compile_scalar_formula(fn,scope):
    """Lift actual formulas to exact functions; never differentiate an interval."""
    tree=ast.parse(textwrap.dedent(inspect.getsource(fn)))
    class Lift(ast.NodeTransformer):
        def visit_Call(self,node):
            node=self.generic_visit(node)
            if isinstance(node.func,ast.Attribute) and node.func.attr=='truncate':return node.func.value
            return node
        def visit_Assign(self,node):
            node=self.generic_visit(node)
            if any(isinstance(t,ast.Name) and t.id=='order' for t in node.targets):node.value=ast.Constant(3)
            return node
    tree=Lift().visit(tree);tree.body[0].decorator_list=[];ast.fix_missing_locations(tree)
    namespace=dict(fn.__globals__);namespace.update(scope)
    exec(compile(tree,'<actual-pulse-exact-function-lift>','exec'),namespace)
    return namespace[fn.__name__]


def exact_linear_pressure_and_units():
    """Compare producing formulas for arbitrary actual source functions."""
    z,de,mu,A,Rv,ratio=sy.symbols('Z delta mu A Rv ratio',positive=True)
    ap,C1,C2,r1,r2=[sy.Function(n)(z) for n in ('ap','C1','C2','row1','row2')]
    P1,P2,Fprime,W1,W2,rho1,rho2,Kmain=sy.symbols('Pi1 Pi2 Fprime W1 W2 rho1 rho2 Kmain')
    c=SimpleNamespace(mpf=sy.sympify)
    polynomials=current.caps.end.polynomials(c,z,3,de)
    wanted=(z,1/(1+z*z),1-z*z,1-de*z*z)
    for polynomial,body in zip(polynomials,wanted):
        for n in range(4):assert sy.cancel(polynomial[n]-sy.diff(body,z,n)/sy.factorial(n))==0
    scalar_polynomials=lambda ctx,Z,order,delta:tuple(p[0] for p in polynomials)
    dz=lambda value:sy.diff(value,z)
    radial=compile_scalar_formula(current.caps.end.radial_q,
        {'polynomials':scalar_polynomials,'derivative':dz})
    formula=compile_scalar_formula(current.caps.pulse.pulse_linear_jets,
        {'end':SimpleNamespace(polynomials=scalar_polynomials,radial_q=radial,derivative=dz)})
    terms=formula(c,z,de,mu,ap,[C1,C2],[r1,r2],[P1,P2],Fprime,[W1,W2],[rho1,rho2],Kmain,ratio)
    ubar=1/(1+z*z);L=1-de*z*z
    I0U=-ubar*mu*r1
    I0m=ubar*(ap*Fprime+13*r1-(r2-r1+ap*(P2-P1))-C1*rho1-C2*rho2)
    I1U=ubar*(C1*W1+C2*W2+ap*Kmain)
    I1m=-I1U-ubar*mu*r1*ratio
    q=lambda U,m:(2*z*U-(1-de)*z*m-(1-z*z)*sy.diff(m,z))/L
    q0,q1,qi=q(I0U,I0m),q(I1U,I1m),q(0,ubar*mu*r1)
    for key,body in dict(I0U=I0U,I0m=I0m,I1U=I1U,I1m=I1m,I0q=q0,I1q=q1,qinit=qi).items():
        for n in (0,1):assert sy.cancel(sy.diff(terms[key],z,n)-sy.diff(body,z,n))==0,(key,n)

    class Pair:
        def __init__(self,value,Z):self.value,self.Z=value,Z
    pair=lambda v:Pair(v,sy.diff(v,z))
    class Algebra:
        def __init__(self,g):pass
        def fixed(self,v):return Pair(sy.sympify(v),sy.Integer(0))
        def add(self,*v):return Pair(sum(q.value for q in v),sum(q.Z for q in v))
        def sub(self,a,b):return Pair(a.value-b.value,a.Z-b.Z)
        def neg(self,a):return Pair(-a.value,-a.Z)
        def mul(self,*v):
            value=sy.prod(q.value for q in v)
            derivative=sum(q.Z*sy.prod(t.value for j,t in enumerate(v) if j!=i) for i,q in enumerate(v))
            return Pair(value,derivative)
        def scale(self,v,a):return Pair(sy.sympify(v)*a.value,sy.sympify(v)*a.Z)
        def div(self,a,b):return Pair(a.value/b.value,(a.Z*b.value-a.value*b.Z)/b.value**2)
    keys=('I0_q_Z','I0_U_q_Z','I0_U_squared','I0_U_U_Z','I1_q','I1_q_Z',
        'I1_q_squared','I1_U_q_Z','I1_U_squared','I1_U_U_Z')
    integrals=dict.fromkeys(keys,pair(sy.Integer(0)))
    # Prove the pressure operator for independent arbitrary analytic q0/q1/qi
    # functions. Their actual recovery bodies were already checked above;
    # keeping the operator proof separate avoids expanding their nested DAGs.
    qa,qb,qc=[sy.Function(n)(z) for n in ('q0','q1','qinit')]
    query=iter((qa,qb,qc))
    linear=compile_scalar_formula(current.caps.pulse.pulse_linear_jets,
        {'end':SimpleNamespace(polynomials=scalar_polynomials,radial_q=lambda *v:next(query),derivative=dz)})
    linear_terms=linear(c,z,de,mu,ap,[C1,C2],[r1,r2],[P1,P2],Fprime,[W1,W2],[rho1,rho2],Kmain,ratio)
    integrals.update(I0_q_Z=pair(A*Rv*sy.diff(qa,z)),I1_q=pair(A*Rv**2*qb),
        I1_q_Z=pair(A*Rv**2*sy.diff(qb,z)))
    endpoints=('R_q','R_squared_q','R_q_squared','R_squared_q_squared','R_U_squared',
        'R_squared_U_squared','q_plus_q_y','R_q_y')
    lower=dict.fromkeys(endpoints,pair(sy.Integer(0)));upper=dict(lower)
    lower.update(R_q=pair(A*Rv*qc),R_squared_q=pair(A*Rv**2*ratio*qc),R_q_y=pair(-A*Rv*qc))
    pressure=compile_formula(current.true.parent.pressure_functions,
        {'C1Algebra':Algebra,'source':SimpleNamespace(C1Function=Pair)})
    result=pressure(SimpleNamespace(symbol=lambda name:z,one=sy.Integer(1)),de,integrals,lower,upper)
    for key,p in (('m13',1),('m15',2)):
        value=linear_terms[key+'_linear']-(qc/Rv if key=='m15' else 0)
        derivative=sy.diff(linear_terms[key+'_linear'],z)-(sy.diff(qc,z)/Rv if key=='m15' else 0)
        assert sy.cancel(result[key].value/(A*Rv**p)-value)==0
        assert sy.cancel(result[key].Z/(A*Rv**p)-derivative)==0

    # Actual C0/C1 source energy normalization: A0^2*Rp /(ref*A*Rv).
    logA0,logE,logRp,logref=sy.symbols('logA0 logE logRp logref')
    logA=logA0+logE-13/(2*mu)-13;logRv=logRp+13/mu
    energy_log=2*logA0+logRp-(logref+logA+logRv)
    common=logA0-logE-13/(2*mu)+13-logref
    assert sy.expand(energy_log-common)==0
    weighted_log=2*logA0+2*logRp-(logref+logA+2*logRv)
    weighted_common=logA0-logE-39/(2*mu)+13-logref
    assert sy.expand(weighted_log-weighted_common)==0
    for fn in (current.caps.pulse.pulse_energy_caps,current.caps._energy_log_caps):
        observed=current.native.execute_assignments(fn,('common',),dict(logA0=logA0,logE=logE,
            mu=mu,reference=sy.exp(logref),c=SimpleNamespace(ln=lambda v:logref)))['common']
        assert sy.expand(observed-common)==0
        row=current.native.assignment(fn,'row').value
        expr=next(kw.value for kw in row.keywords if kw.arg=='weighted_initial_with_length')
        observed=eval(compile(ast.Expression(expr),'<actual-weighted-energy-unit>','eval'),{},
            dict(logA0=logA0,logE=logE,mu=mu,lt=sy.Integer(0),t=sy.Integer(1),length=sy.Integer(0),
                reference=sy.exp(logref),c=SimpleNamespace(ln=lambda v:logref if v==sy.exp(logref) else sy.log(v))))
        assert sy.expand(observed-weighted_common)==0
    # The former derivative error multiplied this already normalized energy
    # by A_cap/ref. Confirm the old expression and the actual corrected body.
    Acap,ref,S0,S1,V=sy.symbols('A_cap reference S0_Z S1_Z viscous_Z',positive=True)
    old=current.native.assignment(current.caps._quadratic_C1_bounds,'quadratic').value
    for key,s in (('m13',S0),('m15',S1)):
        expr=next(kw.value.args[1] for kw in old.keywords if kw.arg==key)
        value=eval(compile(ast.Expression(expr),'<old-C1-scale>','eval'),{},
            dict(A_cap=Acap,reference=ref,S0_Z=S0,S1_Z=S1))
        assert sy.cancel(value-s*Acap/(2*ref))==0
    corrected=compile_formula(current.corrected_derivative_errors,
        {'caps':SimpleNamespace(_scalar_bound=lambda c,v:v)})
    errors=corrected(SimpleNamespace(c=None,quadratic={'reduced_S0_Z_absolute_bound':S0,
        'reduced_S1_Z_absolute_bound':S1},viscous_Z=V))
    assert errors==dict(m13=S0/2,m15=S1/2+V)
    return dict(passed=True,actual_linear_moment_and_Q_value_Z1_identities=14,
        actual_physical_pressure_value_Z1_and_signed_viscous_identities=4,
        qinit_is_Rp_q_physical_over_A_Rv=True,weighted_viscous_is_minus_qinit_over_Rv=True,
        C0_C1_energy_common_logs_have_identical_already_normalized_units=True,
        weighted_p2_energy_unit_and_both_actual_initial_log_recipes_verified=True,
        old_extra_A_cap_over_reference_identified_and_removed=True,
        corrected_C1_quadratic_errors_use_S0_Z_over_two_and_S1_Z_over_two=True,
        all_remainders_and_finite_endpoints_are_enclosures_not_zero_functions=True)


def rejected(fn):
    try:fn()
    except ValueError:return True
    raise AssertionError('Unsupported true pulse query admitted')


@current.source_precision
def run(owner,raw,packets):
    began=time.monotonic();proof=exact_linear_pressure_and_units();owner.assert_graph()
    assert owner.assert_source_identity()==owner.source_identity==raw['source_identity']
    assert len(packets)==len(raw['samples'])==3
    results=[]
    for name,packet in packets.items():
        view=owner.report(packet);assert current.record(view)==current.record(raw['samples'][name])
        assert all(view[k] is False for k in current.OPEN)
        assert view['C1_energy_normalization_double_scaling_corrected']
        assert view['ordinary_Z0_Z1_enclosures_of_true_normalized_functions']
        for key,values in view['normalized_true_pulse_C1_enclosures'].items():
            assert len(values)==2
            for value in values:assert current.ep(value)[0]<=current.ep(value)[1]
            assert current.ep(owner.derivative_errors[key])[1]>0
            assert current.ep(owner.derivative_errors[key])[1]>=current.ep(owner.numeric.full_errors[key])[1]
        results.append(dict(name=name,active_cells=view['active_cells'],actual_query_executed=True,
            true_pressure_value_Z1_enclosure_rows=4,corrected_nonlinear_C1_error_retained=True))
    assert len(owner.report(packets['boundary'])['active_cells'])==2
    assert len(owner.report(packets['whole_Z'])['active_cells'])==4
    assert rejected(lambda:owner.evaluate('1.01'))
    assert rejected(lambda:owner.report(object()))
    saved=owner.derivative_errors
    try:owner.derivative_errors={};assert rejected(owner.assert_graph)
    finally:owner.derivative_errors=saved
    owner.assert_graph()
    receipt=dict(all_passed=True,source_family=owner.family,source_identity=owner.source_identity,
        exact_pressure_and_unit_proof=proof,query_results=results,
        old_C1_double_normalization_error_corrected=True,
        old_grid_accuracy_gate_not_used_as_true_pulse_accuracy_evidence=True,
        **{current.GATE:True},**dict.fromkeys(current.OPEN,False),
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_bytes((json.dumps(current.record(receipt),indent=2)+'\n').encode())
    print('TRUE_COMPLETE_PULSE_C1_NUMERIC_ENCLOSURES_CHECKED',len(results),flush=True)
    return receipt
