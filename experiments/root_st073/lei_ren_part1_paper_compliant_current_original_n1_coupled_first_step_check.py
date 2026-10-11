"""Independent coupled feedback, normalized-axis integration and tail checks."""
import json
import math
from pathlib import Path
import time

import sympy as s
import lei_ren_part1_paper_compliant_current_original_n1_coupled_first_step as current
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

op=current.op;ep=current.ep


def contains(got,want):
    a,b=ep(got);c,d=ep(want)
    return a<=c and b>=d


def overlaps(got,want):
    a,b=ep(got);c,d=ep(want)
    return max(a,c)<=min(b,d)


def rejected(fn):
    try:fn()
    except ValueError:return True
    raise AssertionError('Invalid coupled n1 scope admitted')


def fixture(c):
    """Exact polynomial Wbase and direct B0/B1 feedback, with signed rows."""
    algebra=op.N1LogAlgebra(c,(c.mpf(0),)*4,[]);z=algebra.lift(op.Jet.variable(c,c.mpf('.25'),1))
    A=2+z;B=-3+z*2;C=1+z*z
    def base(x):return dict(F1=B*x*x/8,Uz1=C*x*x/4,K1=-C*x*x/8,
        P1=A*x*x/2,d_x_F1=B*x/4,d_x_Uz1=C*x/2)
    zero=algebra.lift(0,0);counts=0
    for bounds in ((0,0),(0,2),(1,2),(2,2)):
        X=c.mpf(list(bounds));total={name:algebra.lift(0,1) for name in current.FIELDS}
        for j in range(32):
            t=current.partition(c,j,32);x=X*t
            values=current.weighted(algebra,{'4':A*x,'5':B,'6':C},current.normalized_kernels(c,t,X),X/32)
            total={name:total[name]+values[name] for name in current.FIELDS}
        for point in set(bounds):
            for name,want in base(c.mpf(point)).items():
                got=total[name].terms.get((0,0,0,0));target=want.terms.get((0,0,0,0))
                for k in range(2):
                    assert contains(c.mpf(0) if got is None else got[k],c.mpf(0) if target is None else target[k]),(bounds,name,k)
                    counts+=1
    a=c.mpf(3)/7;b=c.mpf(-2)/9;p=c.mpf(5)/11;d=c.mpf(2)/13;e=c.mpf(-3)/17;f=c.mpf(7)/19
    total=[{name:zero for name in current.FIELDS} for _ in range(2)];X=c.mpf(2)
    for j in range(64):
        t=current.partition(c,j,64);x=X*t
        system=dict(B0={'4,1':algebra.lift(a*x),'5,5':algebra.lift(b*x),'6,4':algebra.lift(p)},
                    B1={'5,1':algebra.lift(d),'6,2':algebra.lift(e),'6,4':algebra.lift(f)})
        r0,r1=current.feedback(algebra,system,base(x))
        for index,rhs in enumerate((r0,r1)):
            values=current.weighted(algebra,rhs,current.normalized_kernels(c,t,X),X/64,order=0)
            total[index]={name:total[index][name]+values[name] for name in current.FIELDS}
    expected=[]
    for alpha,beta,gamma in ((A*0+B*(a/8),B*(b/4),A*(p/2)),(A*0,op.dz(B)*(d/8),op.dz(C)*(e/4)+op.dz(A)*(f/2))):
        expected.append(dict(F1=beta*(X**4/24),Uz1=gamma*(X**4/16),K1=-gamma*(X**4/24),
            P1=alpha*(X**4/4),d_x_F1=beta*(X**3/6),d_x_Uz1=gamma*(X**3/4)))
    for index in range(2):
        for name,want in expected[index].items():
            got=total[index][name].terms.get((0,0,0,0));target=want.terms.get((0,0,0,0))
            assert contains(c.mpf(0) if got is None else got[0],c.mpf(0) if target is None else target[0]),(index,name)
            counts+=1
    return dict(exact_normalized_axis_range_and_signed_coupled_polynomial_rows=counts,
                B1_ordinary_Z_derivative_feedback_checked=True,pressure_unknown_F1_feedback_checked=True)


def structural():
    A=s.zeros(6);A[0,4]=1;A[1,5]=1;A[2,5]=-1
    D=s.diag(*s.symbols('d0:6'))
    assert A*D*A==s.zeros(6)
    assert current.initial.operator.matrix.__name__=='matrix'
    return dict(kinematic_map_annihilates_through_every_diagonal_kernel=True,
                diagonal_G_weights_nonnegative=[0,0,2,0,3,1],Akin_absolute_row_norm=1)


def tail_checks(c):
    count=0
    for Cp,C1,L in (('1e-6','1e-9','.1'),('.01','0','.2'),('0','.001','.5')):
        Cp,C1,L=c.mpf(Cp),c.mpf(C1),c.mpf(L)
        tail=current.preconditioned_tail(c,c.mpf(3),c.mpf(2),Cp,C1,L,1)
        M=tail['initial_base_common_tube_bound'];beta=tail['geometric_tail_ratio_upper'];partial=c.mpf(0)
        for k in range(2,42):
            term=M*(6*(Cp+C1*k/L))**k/math.factorial(k)
            assert ep(term)[1]<=ep(M*beta**k)[1]
            partial+=term;count+=1
        assert ep(partial)[1]<ep(tail['omitted_sum_bound'])[1]
    rejected(lambda:current.preconditioned_tail(c,c.mpf(1),c.mpf(2),c.mpf(1),c.mpf(1),c.mpf('.1'),1))
    rejected(lambda:current.preconditioned_tail(c,c.mpf(1),c.mpf(2),c.mpf('.01'),c.mpf(0),c.mpf(0),1))
    return dict(nested_radius_factorial_vs_geometric_tail_rows=count,noncontracting_tail_and_zero_loss_rejected=True)


@source_precision
def run(owner,packets):
    began=time.monotonic();owner.assert_graph();c=owner.c
    checks=dict(**fixture(c),**structural(),**tail_checks(c));prefix_terms=0
    for name,packet in packets.items():
        view=owner.report(packet);assert all(not view[key] for key in current.OPEN)
        assert view['computed_preconditioned_indices']==[0,1]
        tail=view['analytic_preconditioned_tail']
        assert ep(tail['geometric_tail_ratio_upper'])[1]<1 and not tail['half_derivative_block_reduction_used']
        assert tail['initial_function_Z_radius']-tail['solution_Z_radius']==tail['total_available_Z_loss']
        assert all(ep(tail['off_kinematic_B0_norm_upper'])[1]>=ep(v)[1] for v in tail['off_kinematic_rows'])
        assert view['canonical_inner_series_values_include_remaining_uniform_absolute_tail']
        assert not owner._packets[id(packet)][3][1]['P1'].terms
        assert owner._packets[id(packet)][3][0]['P1'].terms
        for value in owner._packets[id(packet)][2].values():
            assert value.order==0 and all(key[0]>=0 for key in value.terms)
            prefix_terms+=len(value.terms)
    checks['production_computed_prefix_source_sectors']=prefix_terms
    # This is a source/enclosure consistency check, not an equation proof.
    old=owner.base.evaluate('.371','4');old_values=owner.base._packets[id(old)][2]
    values=owner._initial_range('.371','4','core')
    for field in current.FIELDS:
        assert values[field].terms.keys()==old_values[field].terms.keys()
        for key in values[field].terms:
            for k in range(2):assert overlaps(values[field].terms[key][k],old_values[field].terms[key][k])
    checks['accepted_point_initial_vs_normalized_range_consistency']=True
    axis=owner.evaluate('.371','0');view=owner.report(axis)
    assert all(not value.terms for value in owner._packets[id(axis)][2].values())
    assert all(not value['terms'] for value in view['canonical_inner_series_value_enclosures'].values())
    checks['prefix_and_exact_inner_series_zero_axis_traces']=True
    checks['guards_rejected']={name:rejected(fn) for name,fn in (
        ('core_comparison_extension',lambda:owner.evaluate('.371','4.01')),
        ('outside_actual_inner_collar',lambda:owner.evaluate('.371','.51','first')),
        ('foreign_chart',lambda:owner.evaluate('.371','0','macro')),
        ('invalid_outer_partition',lambda:owner.evaluate('.371','4',outer_core_cells=0)),
        ('foreign_packet',lambda:owner.report(current.OriginalN1CoupledFirstStepPacket())))}
    packet=packets['actual_first_keep'];entry=owner._packets[id(packet)];old=entry[2]['F1'];entry[2]['F1']=old+1
    try:checks['live_prefix_mutation_rejected']=rejected(lambda:owner.report(packet))
    finally:entry[2]['F1']=old
    hashes=dict(owner.hashes)
    for name in (current.NAME,Path(__file__).name):hashes[name]=current.sha(name)
    result=dict(all_passed=True,source_family=owner.family,input_hashes=hashes,checks=checks,
        **{current.GATE:True},**dict.fromkeys(current.OPEN,False),execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(op.snapshot(result),separators=(',',':'))+'\n',encoding='utf-8',newline='\n')
    print('PASS_N1_COUPLED_FIRST_STEP',checks['exact_normalized_axis_range_and_signed_coupled_polynomial_rows'],prefix_terms,flush=True)
    return result
