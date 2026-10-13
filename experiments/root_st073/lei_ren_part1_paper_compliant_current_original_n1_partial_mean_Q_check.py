"""Independent moving-endpoint mean, radial recovery and pressure identities."""
import json
from pathlib import Path
import time

import sympy as s
import lei_ren_part1_paper_compliant_current_original_n1_partial_mean_Q as current
import lei_ren_part1_paper_compliant_current_original_n1_radial_free_pressure_integrals as pressure


class SymbolicGraph:
    def __init__(self):self.zero=self.ref(s.Integer(0));self.one=self.ref(s.Integer(1))
    def ref(self,value):return current.source.FunctionRef(self,value)
    def constant(self,value):return self.ref(s.Rational(str(value)))
    def symbol(self,name):return self.ref(s.Symbol(name))
    def ids(self,values):assert all(v.graph is self for v in values)
    def add(self,*values):return self.ref(s.Add(*(v.node for v in values)))
    def neg(self,v):return self.ref(-v.node)
    def sub(self,a,b):return self.ref(a.node-b.node)
    def mul(self,*values):return self.ref(s.Mul(*(v.node for v in values)))
    def quotient(self,a,b,proof):return self.ref(a.node/b.node)
    def unary(self,name,v):return self.ref({'exp':s.exp}[name](v.node))
    def node(self,operation,**attrs):
        if operation=='substitute_original_native_coordinate':
            return self.ref(attrs['body'].subs(s.Symbol(attrs['native_variable']),attrs['query_coordinate']))
        if operation=='definite_integral':
            return self.ref(s.integrate(attrs['integrand'],(s.Symbol(attrs['variable']),attrs['lower'],attrs['upper'])))
        raise AssertionError(operation)


def independent_mean_fixture():
    g=SymbolicGraph();z=s.Symbol('Z');t=s.Symbol('native_first_micro');q=s.Symbol('partial_mean_query_first_micro')
    N=s.Integer(7);h=s.Rational(2,5);left=s.Rational(1,3)
    B=(1+z+z*z+z**3)*(t*t+s.sin(N*h*(t-left)))
    incoming=1-z+z**3
    jet=lambda f:current.c3.C3Function(*(g.ref(s.diff(f,z,k)) for k in range(4)))
    base=dict(chart='first_micro',offset=h*(t-left),left_offset=s.Integer(0),domain=[left,s.Integer(2)],
        Jacobian=h,actual_current_source_cells=[])
    got=current.partial_mean_functions(g,base,{'density':{'m':jet(B)},'incoming':{'m':jet(incoming)}},
        {'delta_V':jet(B/N)},g.ref(N))
    # Physical R/Rleft=exp(h*(native-left)); exact normalized radial mean.
    truth=s.exp(-h*(q-left))*incoming/N+s.integrate(s.exp(h*(t-q))*B*h/N,(t,left,q))
    count=0
    for k,value in enumerate(current.c3.rows(got['mean'])):
        assert s.simplify(value.node-s.diff(truth,z,k))==0;count+=1
    for k,value in enumerate(current.c3.rows(got['mean_y'])):
        assert s.simplify(value.node-s.diff(s.diff(truth,q)/h,z,k))==0;count+=1
    assert s.simplify(got['mean'].value.node.subs(q,left)-incoming/N)==0
    assert s.simplify(got['mean'].value.node-got['mean'].value.node.subs(q,2))!=0
    delta=s.Rational(1,5);Q=current.recover_Q_increment(g,got['delta_V'],got['mean'],g.ref(delta))
    d=1-z*z;L=1-delta*z*z
    expected=(2*z*B.subs(t,q)/N-(1-delta)*z*truth-d*s.diff(truth,z))/L
    for k,value in enumerate((Q.value,Q.Z,Q.ZZ)):
        assert s.simplify(value.node-s.diff(expected,z,k))==0;count+=1
    return count


def independent_radial_free_identity():
    z,y,delta=s.symbols('z y delta');R=s.exp(y);m=s.Function('m')(y,z)
    a=1-delta;d=1-z*z;L=1-delta*z*z;U=m+s.diff(m,y)
    q=(2*z*U-a*z*m-d*s.diff(m,z))/L;qz=s.diff(q,z);qy=s.diff(q,y);uz=s.diff(U,z)
    direct=(q+a*z*qz/2+qy)/L+q*(q/2+qy)+U*(-2*z*q+d*qz-2*z*qy)/L
    # Differentiate the proposed finite-interval primitives; arbitrary m,
    # rather than a fitted profile, verifies the whole functional identity.
    C0=R*q/L+R*q*q/2-2*z*z*R*U*U/(L*L)
    C1=R*R*q/L+R*R*q*q/2-2*z*z*R*R*U*U/(L*L)
    local0=a*z*qz/(2*L)+d*U*qz/L-2*delta*z*z*U*U/(L*L)+2*z*d*U*uz/(L*L)
    local1=-q/L+a*z*qz/(2*L)-q*q/2+d*U*qz/L+2*a*z*z*U*U/(L*L)+2*z*d*U*uz/(L*L)
    assert s.simplify(s.expand(R*direct-s.diff(C0,y)-R*local0))==0
    assert s.simplify(s.expand(R*R*direct-s.diff(C1,y)-R*R*local1))==0
    # Verify the implementation coefficients independently with free moments.
    names0=('q_Z','U_q_Z','U_squared','U_U_Z');names1=('q','q_Z','q_squared','U_q_Z','U_squared','U_U_Z')
    I0={k:s.Symbol('I0_'+k) for k in names0};I1={k:s.Symbol('I1_'+k) for k in names1}
    B={k:s.Symbol('B_'+k) for k in ('R_q','R_q_squared','R_U_squared','R_squared_q',
        'R_squared_q_squared','R_squared_U_squared','q_plus_R_q_R','R_squared_q_R')}
    result=pressure.reduced_S_integrals(z,delta,dict(unweighted=I0,weighted=I1),B)
    substitutions0={qz:I0['q_Z'],U*qz:I0['U_q_Z'],U*U:I0['U_squared'],U*uz:I0['U_U_Z']}
    # Pointwise arbitrary-mean proof above establishes the original formulas;
    # free coefficient maps check no accidental weighting or boundary rename.
    expected0=B['R_q']/L+B['R_q_squared']/2-2*z*z*B['R_U_squared']/L**2
    expected0+=a*z*I0['q_Z']/(2*L)+d*I0['U_q_Z']/L-2*delta*z*z*I0['U_squared']/L**2+2*z*d*I0['U_U_Z']/L**2
    expected1=B['R_squared_q']/L+B['R_squared_q_squared']/2-2*z*z*B['R_squared_U_squared']/L**2
    expected1+=-I1['q']/L+a*z*I1['q_Z']/(2*L)-I1['q_squared']/2+d*I1['U_q_Z']/L+2*a*z*z*I1['U_squared']/L**2+2*z*d*I1['U_U_Z']/L**2
    assert s.simplify(result['S']-expected0)==0 and s.simplify(result['R_S']-expected1)==0
    fixed=pressure.fixed_pressure_integrals(z,delta,dict(unweighted=I0,weighted=I1),B)
    assert s.simplify(fixed['N1p']+result['S']/2-B['q_plus_R_q_R'])==0
    assert s.simplify(fixed['R_N1p']+result['R_S']/2-B['R_squared_q_R'])==0
    assert s.simplify(s.diff(q+s.diff(q,y),y)/R-(2*s.diff(q,y)/R+(s.diff(q,y,2)-s.diff(q,y))/R))==0
    assert s.simplify(s.diff(R*s.diff(q,y),y)/R-(s.diff(q,y)+s.diff(q,y,2)))==0
    return True


def rejected(fn):
    try:fn()
    except ValueError:return True
    raise AssertionError('Unsupported partial mean/Q query admitted')


@current.source_precision
def run(owner,packets):
    began=time.monotonic();count=independent_mean_fixture();identity=independent_radial_free_identity()
    for packet in packets.values():
        view=owner.report(packet)
        assert all(not view[k] for k in current.OPEN)
        assert view['whole_cell_range_not_outgoing_endpoint_value']
        assert view['modified_Q_increment_only_not_complete_Q']
        assert len(view['whole_cell_partial_own_mean_increment_Z3_magnitude'])==4
        assert len(view['whole_cell_Q_increment_Z2_magnitude'])==3
        assert view['complete_pressure_source_tail_support_and_repair_not_installed']
    token=packets['first_micro'];entry=owner._packets[id(token)]
    key='modified_Q_increment_only_not_complete_Q';old=entry[1][key];entry[1][key]=False
    rejected(lambda:owner.report(token));entry[1][key]=old
    old=owner.delta;owner.delta=owner.c.mpf('.1');rejected(lambda:owner.report(token));owner.delta=old
    rejected(lambda:owner.range_cell(('0','.5'),'first_micro',(0,1),(1,2)))
    rejected(lambda:owner.range_cell(('0','.5'),'second_micro',(3,2),(1,1)))
    rejected(lambda:owner.report(current.OriginalPartialMeanQPacket()))
    owner.g.symbol('next_pressure_integral_child');owner.report(token)
    hashes=dict(owner.hashes)
    for name in (current.NAME,Path(__file__).name,Path(pressure.__file__).name):hashes[name]=current.sha(name)
    result=dict(all_passed=True,source_family=owner.family,input_hashes=hashes,
        checks=dict(independent_free_endpoint_mean_Z3_mean_y_Z3_and_Q_Z2_rows=count,
            nonzero_incoming_physical_radial_mean_and_free_endpoint_preserved=True,
            arbitrary_smooth_leading_mean_radial_free_pressure_integral_identities=identity,
            finite_viscous_boundary_signs_and_weights_verified=True,
            two_actual_micro_whole_cell_partial_mean_and_Q_increment_queries_execute=True,
            same_N_delta_source_scope_and_parent_composition_guards=True,
            full_pressure_tails_support_and_n1_repair_not_admitted=True),
        **{current.GATE:True},**dict.fromkeys(current.OPEN,False),execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_bytes((json.dumps(current.c3.encoded(result),separators=(',',':'))+'\n').encode())
    print('PASS_PARTIAL_MEAN_Q_RADIAL_FREE_PRESSURE',count,flush=True);return result
