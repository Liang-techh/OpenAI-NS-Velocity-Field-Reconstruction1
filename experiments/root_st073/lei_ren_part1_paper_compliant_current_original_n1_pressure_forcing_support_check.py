"""Actual zero field factory, empty end support and exact future coverage."""
import ast
import inspect
import json
from pathlib import Path
import textwrap

import mpmath as mp
import sympy as s
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
import lei_ren_part1_paper_compliant_current_original_n1_pressure_forcing_support as current


def actual_zero_factory():
    c=mp.iv;c.dps=80
    z=IntervalTaylor.variable(c,c.mpf(['-.5','.5']),5)
    theta=2+z+z*z;one=IntervalTaylor.constant(c,1,5)
    got=current.mixed.flatten_mixed(c,c.mpf('.1'),one,[c.mpf('.2')]*5,
                                  theta,one*3,one*4,one*5,c.mpf('.7'))
    checked=0
    for label in (current.mixed.UZ,current.mixed.UR):
        for value in got['physical_mixed_derivatives_total_order_le4'][label].values():
            assert current.physical.ep(value)==(mp.mpf(0),mp.mpf(0));checked+=1
    assert current.physical.ep(got['physical_mixed_derivatives_total_order_le4'][current.mixed.UT]['y0_Z0'])[0]>0
    assert current.physical.ep(got['physical_mixed_derivatives_total_order_le4'][current.mixed.P]['y0_Z0'])[0]>0
    # Both end bumps are strictly to the left of this entire open end collar.
    for center in (-3,-1):
        local=c.mpf(['-.8','0'])-center
        assert current.physical.ep(current.axial.raw_beta(c,local/c.mpf('.15')))==(mp.mpf(0),mp.mpf(0))
        assert all(current.physical.ep(value)==(mp.mpf(0),mp.mpf(0)) for value in
                   current.axial.backward_bump_weights(c,c.mpf(['.01','.1']),c.mpf([1,2]),local,4))
    tree=ast.parse(textwrap.dedent(inspect.getsource(current.axial.CompliantAxialPulseField.end)))
    pairs=[node.iter for node in ast.walk(tree) if isinstance(node,ast.For)]
    assert any(ast.unparse(value)=='zip(C, (-3, -1))' for value in pairs)
    return checked


def future_coverage(owner):
    data=owner.reports[current.radius.NAME];f=data['exact_waiting_functions'];nodes=owner.nodes
    rp,mu,L,T,W,v=s.symbols('logRp mu Lrel Ts waiting v',positive=True)
    binding={data['original_logRp_node']:rp,f['mu']:mu,f['Lrel']:L,f['Ts']:T,f['waiting']:W}
    cache={}
    def at(i):
        if i in binding:return binding[i]
        if i in cache:return cache[i]
        row=nodes[i];op=row['operation']
        if op=='exact_rational':value=s.Rational(row['numerator'],row['denominator'])
        elif op=='bound_variable':
            assert row['name']=='current_Rp_radius_coordinate';value=v
        elif op=='sum':value=s.Add(*(at(j) for j in row['arguments']))
        elif op=='product':value=s.Mul(*(at(j) for j in row['arguments']))
        elif op=='negative':value=-at(row['argument'])
        elif op=='positive_quotient':value=at(row['numerator'])/at(row['denominator'])
        else:raise AssertionError('Unexpected scalar radius dependency '+op)
        cache[i]=value;return value
    span={'flatten':(0,100),'outer_power':(0,1),'outer_angular':(-4,0),
          'steep_entry':(0,1),'steep_power':(0,1),'steep_exit':(0,1),
          'waiting':(0,1),'heat_collar':(0,3),'heat_exterior':(3,None)}
    previous=rp+13/mu;rows=[]
    for chart in current.FUTURE:
        expression=at(owner.maps[chart]['logR']);a,b=span[chart]
        lo=s.expand(expression.subs(v,a));assert s.expand(lo-previous)==0,(chart,'radius join')
        hi=None if b is None else s.expand(expression.subs(v,b))
        rows.append(dict(chart=chart,exact_lower=s.sstr(lo),exact_upper='infinity' if hi is None else s.sstr(hi)))
        if hi is not None:previous=hi
    assert at(owner.maps['flatten']['pulse_term'])==13/mu
    return rows


def independent_forcing_zero():
    R,z,delta=s.symbols('R Z delta',positive=True)
    # For physical meridional velocity Ur, q=Ur/sqrt(R/2).
    Ur,U,m=s.Integer(0),s.Integer(0),s.Integer(0)
    q=(2*z*U-(1-delta)*z*m-(1-z*z)*s.diff(m,z))/(1-delta*z*z)
    qr=(s.diff(Ur,s.Symbol('y'))-Ur/2)/(R*s.sqrt(R/2))
    assert q==qr==0
    S=(q+(1-delta)*z*s.diff(q,z)/2+R*s.diff(q,R))/(1-delta*z*z)
    S+=q*(q/2+R*s.diff(q,R))+U*(-2*z*q+(1-z*z)*s.diff(q,z)-2*z*R*s.diff(q,R))/(1-delta*z*z)
    force=-S/2+2*qr+R*s.diff(qr,R)
    assert S==force==R*force==0
    assert all(value==0 for value in (R*q,R*q*q,R*U*U,R*R*q,R*R*q*q,R*R*U*U,q+R*qr,R*R*qr))
    return True


def actual_known_pressure_zero():
    c=mp.iv;op=current.operator;algebra=op.N1LogAlgebra(c,(c.mpf(0),)*4,[])
    zero=algebra.lift(op.Jet.constant(c,0,4))
    angular=algebra.lift(op.Jet.constant(c,1,4))
    z=algebra.lift(op.Jet.variable(c,c.mpf(['-.5','.5']),4))
    rows=dict(F=[angular,zero,zero],Uz=[zero]*3,Q=[zero]*3)
    out=op.known_forcing(algebra,rows,z,c.mpf(2),c.mpf('.2'))
    assert all(current.physical.ep(value)==(mp.mpf(0),mp.mpf(0))
               for jet in out['N1p'].terms.values() for value in jet.coefficients)
    return True


def run(owner):
    count=actual_zero_factory();coverage=future_coverage(owner);independent_forcing_zero()
    actual_known_pressure_zero()
    for chart in current.FUTURE:
        view=owner.future_zero(chart);assert all(not view[k] for k in current.OPEN)
        assert not view['full_velocity_support_or_angular_heat_pressure_closure_implied']
    try:owner.future_zero('O3_power')
    except ValueError:pass
    else:raise AssertionError('Loop flatness incorrectly admitted as physical meridional support')
    hashes=dict(owner.hashes)
    for name in (current.NAME,Path(__file__).name):hashes[name]=current.sha(name)
    receipt=dict(all_passed=True,source_family=owner.family,input_hashes=hashes,
        checks=dict(actual_Uz_Ur_mixed_y_Z_zero_rows=count,
                    actual_native_end_beta_and_backward_histories_zero_on_open_Rv_collar=True,
                    exact_same_source_nine_chart_log_radius_coverage=coverage,
                    known_physical_S_N1p_and_all_upper_boundary_terms_exact_zero=True,
                    actual_original_n1_known_pressure_operator_zero_with_nonzero_angular_field=True,
                    original_nonzero_angular_and_pressure_rows_retained=True,
                    local_loop_flatness_not_used_as_support=True),
        **{current.GATE:True},**dict.fromkeys(current.OPEN,False))
    (current.HERE/current.RECEIPT).write_bytes((json.dumps(receipt,separators=(',',':'))+'\n').encode())
    print('PASS_N1_PRESSURE_FORCING_SUPPORT',count,len(coverage),flush=True)
    return receipt
