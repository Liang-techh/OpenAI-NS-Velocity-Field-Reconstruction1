"""Independent regular kernels, interval integration, and source guards."""
import json
import math
from pathlib import Path
import time

import sympy as s
import lei_ren_part1_paper_compliant_current_original_n1_initial_integral as current
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

op=current.operator


def contains(got,want):
    a,b=current.ep(got);c,d=current.ep(want)
    return a<=c and b>=d


def rejected(fn):
    try:fn()
    except ValueError:return True
    raise AssertionError('Unsupported integral scope admitted')


def exact_kernels():
    """Derive the integrated kinematic feedback kernels independently."""
    x,t,u=s.symbols('x t u',positive=True)
    f=s.integrate((u/t)**3,(t,u,x))
    uz=s.integrate(u/t,(t,u,x))
    k=-s.integrate(t*t*(u/t),(t,u,x))/(x*x)
    assert s.simplify(f-u*(1-(u/x)**2)/2)==0
    assert s.simplify(uz-u*s.log(x/u))==0
    assert s.simplify(k+u*(1-(u/x)**2)/2)==0
    return dict(exact_double_regular_inverse_kinematic_kernels=3)


def polynomial_integration(c):
    algebra=op.N1LogAlgebra(c,(c.mpf(0),)*4,[])
    J=op.Jet;z=algebra.lift(J.variable(c,c.mpf('.25'),1))
    a=2+z;b=-3+z*2;v=1+z*z;endpoint=c.mpf(2);N=32
    total={key:algebra.lift(0,1) for key in current.FIELDS}
    for j in range(N):
        cell=c.mpf([c.mpf(j).a*2/N,c.mpf(j+1).b*2/N])
        values=current.integrate_cell(algebra,dict(g={'4':a*cell,'5':b,'6':v}),cell,endpoint,c.mpf(2)/N)
        total={key:total[key]+values[key] for key in current.FIELDS}
    expected=dict(F1=b/2,Uz1=v,K1=-v/2,P1=a*2,d_x_F1=b/2,d_x_Uz1=v)
    count=0
    for field,value in total.items():
        for k in range(2):
            assert contains(value.terms[(0,0,0,0)][k],expected[field].terms[(0,0,0,0)][k]),field
            count+=1
    # A separate finite-width fixture checks measure cancellation against
    # independently integrated exact collar kernels. Its rational hb is
    # ONLY a test parameter, never a production source representative.
    hb=c.mpf(1)/7;phase=c.mpf('.5');endpoint=2*c.exp(hb*phase/2)
    total={key:algebra.lift(0,1) for key in current.FIELDS}
    for j in range(N):
        interval=c.mpf([c.mpf(j).a/(2*N),c.mpf(j+1).b/(2*N)])
        cell=2*c.exp(hb*interval/2)
        g={key:algebra.width(value,-1) for key,value in (('4',a*cell),('5',b),('6',v))}
        values=current.integrate_cell(algebra,dict(g=g),cell,endpoint,cell/(4*N),True)
        total={key:total[key]+values[key] for key in current.FIELDS}
    quadratic=endpoint*endpoint/8-1+2/(endpoint*endpoint)
    expected=dict(F1=b*(quadratic/hb),Uz1=v*((endpoint*endpoint/4-2*c.ln(endpoint/2)-1)/hb),
        K1=-v*(quadratic/hb),P1=a*((endpoint*endpoint-4)/(2*hb)),
        d_x_F1=b*((endpoint**4-16)/(4*hb*endpoint**3)),d_x_Uz1=v*((endpoint*endpoint-4)/(2*hb*endpoint)))
    for field,value in total.items():
        assert set(value.terms)=={(0,0,0,0)}
        for k in range(2):
            assert contains(value.terms[(0,0,0,0)][k],expected[field].terms[(0,0,0,0)][k]),field
            count+=1
    return dict(exact_core_and_finite_width_collar_polynomial_integral_rows=count,
        signed_source_and_ordinary_Z_derivatives_checked=True,collar_inverse_width_cancelled_in_test=True)


def kernel_bounds(c):
    count=0;precise=type(c)();precise.dps=c.dps*2
    for a,b in ((0,1),(0,2),('.7','.8'),('1.9',2)):
        cell=c.mpf([a,b]);bound=current.log_kernel(c,cell,c.mpf(2))
        for j in range(9):
            value=precise.mpf(a) if j==0 else precise.mpf(b) if j==8 else precise.mpf(a)+(precise.mpf(b)-precise.mpf(a))*j/8
            want=precise.mpf(0) if j==0 and a==0 else value*precise.ln(2/value)
            assert contains(bound,want),(a,b,j)
            count+=1
    return dict(axis_and_stationary_point_log_kernel_containment_rows=count)


@source_precision
def run(owner,packets):
    began=time.monotonic();owner.assert_graph();c=owner.c
    checks=dict(**exact_kernels(),**polynomial_integration(c),**kernel_bounds(c));terms=0
    for name,packet in packets.items():
        view=owner.report(packet);assert view['computed_iterate']=='Wbase=(I+J)Gg'
        assert view['rigorous_entire_cell_interval_integration'] and view['coupled_B0_B1_feedback_still_open']
        assert all(not view[key] for key in current.OPEN)
        for value in owner._packets[id(packet)][2].values():
            assert value.order==1 and all(key[0]>=0 for key in value.terms)
            terms+=len(value.terms)
    checks['production_computed_field_source_sectors']=terms
    axis=owner.evaluate('.371','0');assert all(not value.terms for value in owner._packets[id(axis)][2].values())
    checks['all_six_exact_zero_axis_traces']=True
    inlet=owner.evaluate('.371','0','first')
    assert op.snapshot(owner._packets[id(inlet)][2])==op.snapshot(owner._packets[id(packets['core_exit'])][2])
    checks['source_integral_core_exit_equals_exact_first_inlet']=True
    checks['guards_rejected']={name:rejected(fn) for name,fn in (
        ('core_extension',lambda:owner.evaluate('.371','4.01')),
        ('collar_outside_Rin',lambda:owner.evaluate('.371','.51','first')),
        ('wrong_chart',lambda:owner.evaluate('.371','0','macro')),
        ('invalid_partition',lambda:owner.evaluate('.371','4',core_cells=0)),
        ('interval_endpoint',lambda:owner.evaluate('.371',c.mpf([1,2]))),
        ('foreign_packet',lambda:owner.report(current.OriginalN1InitialIntegralPacket())))}
    packet=packets['actual_first_keep'];entry=owner._packets[id(packet)];original=entry[2]['P1']
    entry[2]['P1']=original+1
    try:checks['live_integral_mutation_rejected']=rejected(lambda:owner.report(packet))
    finally:entry[2]['P1']=original
    hashes=dict(owner.hashes)
    for name in (current.NAME,Path(__file__).name):hashes[name]=current.sha(name)
    result=dict(all_passed=True,source_family=owner.family,input_hashes=hashes,checks=checks,
        **{current.GATE:True},**dict.fromkeys(current.OPEN,False),execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(op.snapshot(result),separators=(',',':'))+'\n',encoding='utf-8',newline='\n')
    print('PASS_N1_INITIAL_REGULAR_INTEGRALS',checks['exact_core_and_finite_width_collar_polynomial_integral_rows'],terms,flush=True)
    return result
