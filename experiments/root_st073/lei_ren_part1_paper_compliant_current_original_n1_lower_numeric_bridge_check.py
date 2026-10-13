"""Exact collar-difference/recovery/density identities and genuine C1 queries."""
import json
from pathlib import Path
import time
from types import SimpleNamespace

import sympy as sy

import lei_ren_part1_paper_compliant_current_original_n1_lower_numeric_bridge as current
import lei_ren_part1_paper_compliant_current_original_n1_inner_numeric_bridge_check as mathcheck

prior,op,ep=current.prior,current.op,current.ep


def exact_lower_pressure_identities():
    """Execute producing ASTs with arbitrary analytic source functions."""
    z,r=sy.symbols('Z rho',real=True,positive=True)
    eps,hb,de,lo,hi=sy.symbols('epsilon hb delta phase_lower phase_upper',positive=True)
    psi=sy.Function('Psi_core')(r,z);phi=sy.Function('Phi_core')(r,z)
    mean=sy.Function('mean_Psi_core')(r,z)
    ratio,chi,chir,dpsi,dmean=sy.symbols('ratio chi chi_rho dPsi dmeanPsi')
    corepsi=[sy.diff(psi,r,i) for i in range(3)]
    corephi=[sy.diff(phi,r,i) for i in range(3)]
    difference=prior.collar_derivative_corrections(corepsi,corephi,chi,chir,ratio,dpsi,dmean,r)
    wanted=dict(dpsi1=(ratio*chi-1)*corepsi[1],
        dpsi2=ratio*chir*corepsi[1]+ratio*chi*(chi-1)*(corephi[1]/corephi[0])*corepsi[1]
            +(ratio*chi-1)*corepsi[2],
        dmean1=(dpsi-dmean)/r,
        dmean2=((ratio*chi-1)*corepsi[1]-2*(dpsi-dmean)/r)/r)
    for key,want in wanted.items():assert sy.cancel(difference[key]-want)==0,key
    # Derive the second radial transport independently from the ratio ODE.
    rr,cc=sy.Function('ratio')(r,z),sy.Function('chi')(r)
    transported=rr*cc*sy.diff(psi,r)
    differentiated=sy.diff(transported,r).subs(sy.diff(rr,r),(cc-1)*sy.diff(phi,r)*rr/phi)
    true_second=rr*(sy.diff(cc,r)*sy.diff(psi,r)+cc*(cc-1)*sy.diff(phi,r)*sy.diff(psi,r)/phi
        +cc*sy.diff(psi,r,2))
    assert sy.cancel(differentiated-true_second)==0

    class ScalarAlgebra:
        def __init__(self):self.ctx=SimpleNamespace(mpf=lambda v:sy.Rational(str(v)))
        def lift(self,v,order=None):return sy.sympify(v)
        def shift(self,v,p):return v*eps**(-p[3])
        def width(self,v,power=1):return v*hb**power
    algebra=ScalarAlgebra()
    fakeop=SimpleNamespace(dz=lambda v:sy.diff(v,z),truncate=lambda v,n:v)
    leading=mathcheck.compile_mathematical(prior.mean.factored_leading_rows,{'op':fakeop})
    uaxis=sy.Function('Uaxis')(z);U=uaxis+eps*psi;M=uaxis+eps*mean
    rows=dict(Psi=corepsi,meanPsi=[sy.diff(mean,r,i) for i in range(3)],
        axis=dict(z=z,u=uaxis,d=1-z*z,L=1-de*z*z))
    fields=leading(algebra,rows,SimpleNamespace(delta=de))
    Q=(2*z*U-(1-de)*z*M-(1-z*z)*sy.diff(M,z))/(1-de*z*z)
    for i in range(3):
        for name,want in (('Uz',U),('mean',M),('Q',Q)):
            assert sy.cancel(fields[name][i]-sy.diff(want,r,i))==0,(name,i)

    # Compare the actual true C1 differential-pressure AST with the legacy
    # forcing AST for arbitrary Q,U. Every Z derivative is symbolic and
    # ordinary. No numerical fixture or derivative of a cap is involved.
    class Pair:
        def __init__(self,value,Z):self.value,self.Z=value,Z
    class C1:
        def __init__(self,g):pass
        def fixed(self,v):return Pair(sy.sympify(v),sy.Integer(0))
        def add(self,*v):return Pair(sum(q.value for q in v),sum(q.Z for q in v))
        def sub(self,a,b):return Pair(a.value-b.value,a.Z-b.Z)
        def scale(self,c,v):
            c=sy.sympify(c);return Pair(c*v.value,c*v.Z)
        def mul(self,*values):
            result=self.fixed(1)
            for b in values:
                result=Pair(result.value*b.value,result.Z*b.value+result.value*b.Z)
            return result
        def div(self,a,b):return Pair(a.value/b.value,(a.Z*b.value-a.value*b.Z)/b.value**2)
    q=sy.Function('Q')(r,z);u=sy.Function('U')(r,z)
    qy=r*sy.diff(q,r);qyy=r*sy.diff(qy,r)
    jet=lambda v:SimpleNamespace(value=v,Z=sy.diff(v,z),ZZ=sy.diff(v,z,2))
    fakeg=SimpleNamespace(constant=lambda v:sy.Rational(str(v)),symbol=lambda v:z,one=sy.Integer(1),
        zero=sy.Integer(0),nodes=[{'name':'phase'}],mul=lambda *v:sy.Mul(*v),quotient=lambda a,b,*v:a/b)
    lowerfields=dict(U=jet(u),Q=jet(q),Q_y=jet(qy),Q_yy=jet(qyy),radius=eps*r,
        physical_measure=eps*hb*r,weighted_physical_measure=eps**2*hb*r**2,
        query=SimpleNamespace(node=0),domain=dict(lower=lo,upper=hi))
    controls=SimpleNamespace(integral=lambda g,body,var,low,high,**kw:body*(high-low))
    actual=mathcheck.compile_mathematical(current.true.direct_pressure_functions,
        {'parent':SimpleNamespace(C1Algebra=C1),'source':SimpleNamespace(C1Function=Pair),
         'controls':controls})
    _,bindings=actual(fakeg,lowerfields,dict(delta=de))
    axial=mathcheck.compile_mathematical(op.axial_second,{'dz':lambda v:sy.diff(v,z)})
    legacy=mathcheck.compile_mathematical(op.known_forcing,
        {'dz':lambda v:sy.diff(v,z),'axial_second':axial,'truncate':lambda v,n:v})
    result=legacy(algebra,dict(F=[sy.Integer(0)]*3,Uz=[sy.diff(u,r,i) for i in range(3)],
        Q=[sy.diff(q,r,i) for i in range(3)]),z,r,de)['N1p']
    assert sy.cancel(bindings['N1p'].value-result)==0
    assert sy.cancel(bindings['N1p'].Z-sy.diff(result,z))==0
    for name,measure in (('m13',eps*hb*r),('m15',-eps**2*hb*r*r)):
        pair=bindings['signed_physical_integrands'][name]
        assert sy.cancel(pair.value-measure*result)==0
        assert sy.cancel(pair.Z-measure*sy.diff(result,z))==0

    # The interval integrator is a positive-measure RANGE enclosure. Only
    # for a constant radius would this product be the exact integral.
    covers=prior.physical_pressure_range_integrals(algebra,result,r,lo,hi)
    assert sy.cancel(covers['m13']-eps*hb*r*(hi-lo)*result)==0
    assert sy.cancel(covers['m15']+eps**2*hb*r*r*(hi-lo)*result)==0
    return dict(passed=True,universal_source_difference_equations=4,
        ratio_ODE_implies_true_second_radial_transport=True,
        leading_physical_U_mean_Q_rho0_rho1_rho2_identities=9,
        actual_true_and_numeric_N1p_value_Z1_identities=2,
        both_signed_physical_density_value_Z1_identities=4,
        positive_measure_range_cover_factor_identities=2,
        range_product_equality_to_actual_integral_not_asserted=True,
        leading_one_minus_delta_and_both_viscous_derivatives_verified=True,
        accepted_production_mathematical_AST_executed=True)


def rejected(fn):
    try:fn()
    except ValueError:return True
    raise AssertionError('Unsupported true lower numerical query admitted')


@current.source_precision
def run(owner,raw,packets):
    began=time.monotonic();proof=exact_lower_pressure_identities();owner.assert_graph()
    identity=owner.assert_source_identity()
    assert identity==owner.source_identity and identity['passed']
    assert raw['true_function_bindings']==owner.bindings
    assert len(packets)==len(raw['samples'])==3
    results=[]
    for name,packet in packets.items():
        view=owner.report(packet);assert op.snapshot(view)==op.snapshot(raw['samples'][name])
        assert all(view[k] is False for k in current.OPEN)
        assert view['ordinary_Z0_Z1_enclosures_of_true_functions']
        assert view['range_product_encloses_integral_not_a_selected_value']
        source=view['original_collar_source_enclosures'];error=source['actual_source_reference_and_error_proof']
        assert error['reference_core_is_not_selected_as_actual_field'] and error['actual_comparison_errors_retained']
        assert ep(error['Cauchy_Taylor_gap'])[0]>0
        for row in error['original_sigma_bounds']:
            assert row['source_hb_power']==2 and ep(row['nonzero_coefficient_upper'])[0]>0
        assert view['true_function_graph_and_numerical_enclosure_context_kept_separate']
        for key,row in view['true_lower_pressure_C1_enclosures'].items():
            assert row['ordinary_Z_derivative_order']==1 and row['terms']
            assert all(len(term['ordinary_Z_derivatives'])==2 for term in row['terms'])
            assert all(term['source_exponents'][0]>=1 for term in row['terms'])
            assert all(term['source_exponents'][1]==0 for term in row['terms'])
        combined=view['inner_plus_lower_prefix']
        if name=='point':
            assert combined is not None and set(combined['physical_five_prefix_moments_Z1'])==set(prior.moments.MOMENTS)
            assert combined['pressure_m13_m15_stop_at_rminus'] and combined['preRp_and_pulse_pressure_not_added']
            ce=owner.cutoff._packets[id(owner._packets[id(packet)][3])][4][3]
            lv=owner._packets[id(packet)][4][0]
            for key in prior.moments.MOMENTS:
                expected=ce[key]+lv[key] if key in lv else ce[key]
                assert op.snapshot(op.export(expected,1))==op.snapshot(combined['physical_five_prefix_moments_Z1'][key])
        else:assert combined is None
        results.append(dict(name=name,Z=view['Z'],actual_original_query_executed=True,
            true_signed_pressure_value_Z1_rows=4,
            whole_Z_range_enclosure=name=='whole_Z',Rin_to_rminus_integrals_only=True,
            composed_with_actual_Rin_history=combined is not None))
    assert ep(owner.report(packets['whole_Z'])['Z'])==(-1,1)
    assert rejected(lambda:owner.evaluate('1.01'))
    assert rejected(lambda:owner.report(prior.OriginalUnchangedPressurePrefixPacket()))
    npacket=owner._packets[id(packets['point'])][2]
    assert rejected(lambda:owner.evaluate('.211',reuse_packet=npacket))
    wrong_Rin=next(entry[0] for entry in owner.cutoff._packets.values() if entry[1]['chart']=='first'
        and entry[1]['coordinate']._mpi_==owner.c.mpf(1)._mpi_)
    assert rejected(lambda:owner.evaluate('.211',reuse_packet=owner._packets[id(packets['independent_point'])][2],
        Rin_packet=wrong_Rin))
    original=owner._bindings
    try:
        owner._bindings={};assert rejected(owner.assert_graph)
    finally:owner._bindings=original
    owner.assert_graph();pins={**owner.hashes,current.NAME:current.sha(current.NAME)}
    receipt=dict(all_passed=True,source_family=owner.family,source_identity=identity,
        exact_mathematical_proof=proof,query_results=results,input_hashes=pins,
        checks=dict(actual_lower_C1_functions_enclosed_by_same_source_comparison_errors=True,
            two_original_points_and_whole_Z_range_execute=True,
            physical_source_sectors_and_all_nonzero_errors_retained=True,
            actual_Rin_plus_lower_prefix_composition_verified=True,
            foreign_wrong_Z_and_changed_binding_packets_rejected=True,
            full_external_pressure_and_repaired_family_not_promoted=True),
        **{current.GATE:True},**dict.fromkeys(current.OPEN,False),
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_bytes((json.dumps(op.snapshot(receipt),indent=2)+'\n').encode('utf8'))
    print('TRUE_LOWER_PRESSURE_C1_NUMERICAL_BRIDGE_CHECKED',len(results),flush=True)
    return receipt


if __name__=='__main__':
    proof=exact_lower_pressure_identities()
    print(json.dumps(proof,indent=2))
