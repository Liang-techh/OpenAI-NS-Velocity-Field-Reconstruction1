"""Closed quotient/q derivatives, true whole radial cells and branch guards."""
import gzip
import json
import math
from pathlib import Path
from types import SimpleNamespace
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_Rm_positive_quotients_q as current
from lei_ren_part1_paper_compliant_current_original_bridge_macro_functions import MacroFlow
from lei_ren_part1_paper_interval_taylor import IntervalTaylor

fields,ep=current.fields,current.ep


def closed_fixture(zvalue):
    c=MPIntervalContext();c.dps=200;z=sy.Symbol('Z',real=True);R=sy.Rational;zz=mp.mpf(zvalue)
    E=1+z/4+z*z/10;a=R(4,5)+R(3,100)*z*z;b=R(1,10)-z/25+z*z/50
    p1bar=R(4,5)+z/5-z**3/10;p2bar=-R(1,5)+z/7+R(3,100)*z*z
    eta=R(7,1000);dstar=R(1,200);kap=a+b*b/a;Delta=kap-2
    q2=(2*eta-Delta)/(2*a);q=sy.sqrt(q2);P0=1+z/5
    f=MacroFlow(c,c.mpf(-30),c.ln(9),c.ln(4),c.ln(5),c.ln(100)-c.ln(5),c.mpf('.01'))
    tol=mp.mpf('1e-140')
    def jet(expr):
        rows=[]
        for n in range(6):
            value=sy.lambdify(z,sy.diff(expr,z,n)/math.factorial(n),'mpmath')(zz)
            rows.append(c.mpf([value-tol,value+tol]))
        return f.jet(IntervalTaylor(c,rows))
    op=SimpleNamespace(flow=f,c=c,P0=jet(P0),Pstar=f.factor((0,.5,0,0,0)),Rm_factor=f.scalar(7))
    source=dict(common_original_P0_axial5=op.P0,source_geometry=dict(point=True,exact_x=[5,4]),
        original_fixed_source_log_bases=f.logs,source_ledger_is_same_object=True,actual_patch_source_owner=type(op).__name__,
        actual_generic_source_numerators=dict(E=jet(E),C=jet(E*a),B=jet(E*b),positive_denominator=jet(E*E*a)),
        full_signed_inertial_sectors_axial4=dict(theta_linear=jet(E*p1bar/3)[:5],
            theta_quadratic=current.scale(jet(E*p1bar*R(2,3))[:5],f.factor((0,-.5,0,0,0))),
            axial_linear=jet(E*p2bar*R(3,5))[:5],
            axial_quadratic=current.scale(jet(E*p2bar*R(2,5))[:5],f.factor((0,-.5,0,0,0)))),
        original_physical_radius=f.scalar(c.mpf(35)/4))
    result=current.recover_quotients(op,source,c.ln(c.mpf(7)/1000),c.ln(c.mpf(1)/200));targets=[]
    for key,expr in dict(a=a,b=b,t0=-b/a,kappa=kap,kappa_minus2=Delta).items():
        targets.append((result['shear_quotients_axial5'][key],expr))
    for key,expr in dict(p1=p1bar,p2=p2bar).items():
        targets.append((result['full_signed_inertial_quotients_before_shared_R_axial4'][key],expr))
        targets.append((result['actual_factored_p1_p2_axial4'][key],expr*R(35,4)))
    targets.extend(((result['H0_over_shared_R_axial4'],p1bar-p2bar*b/a),
        (result['J_over_shared_R_axial4'],p2bar+p1bar*b/a),
        (result['original_shear_q']['q_axial_coefficients'],q),
        (result['original_shear_q']['q_squared_axial_coefficients'],q2),
        (result['original_signed_u_before_shared_R_axial4'],p2bar*q/dstar)))
    comparisons=0
    for rows,expr in targets:
        for n,row in enumerate(rows):
            value=sy.lambdify(z,sy.diff(expr,z,n)/math.factorial(n),'mpmath')(zz)
            cover=row.coefficient*c.exp(row.scale.evaluate());lo,hi=ep(cover)
            assert lo<=value<=hi,(zvalue,n,str(expr),mp.nstr(value,15));comparisons+=1
    assert result['original_shear_q']['original_C0_lazy_cutoff']['branch']=='active'
    assert result['common_original_P0_axial5'] is op.P0
    return dict(passed=True,fixture_Z=zvalue,independent_closed_quotient_q_and_signed_u_derivative_comparisons=comparisons,
        source_tolerance_only_for_finite_reference_fixture=True,no_native_factor_materialization=True)


def guard_checks():
    c=MPIntervalContext();c.dps=120
    f=MacroFlow(c,c.mpf(-30),c.ln(9),c.ln(4),c.ln(5),c.ln(100)-c.ln(5),c.mpf('.01'))
    rejects=0
    for row in (f.scalar(0),f.scalar(-1),f.scalar(c.mpf([-1,1]))):
        try:current.positive_source(f,row,'invalid source')
        except ValueError:rejects+=1
        else:raise AssertionError('Unproved denominator accepted')
    a=[f.scalar(c.mpf('.8'))]+[f.scalar(0)]*5;proof=current.positive_source(f,a[0],'a')
    try:current.quotient(f,a,[f.scalar(c.mpf('.8'))]+a[1:],proof)
    except ValueError:rejects+=1
    else:raise AssertionError('Different denominator source object accepted')
    flat=current.shear_q_jets(f,a,[f.scalar(c.mpf('.2'))]+a[1:],c.ln(c.mpf('.1')),proof)
    assert flat['original_C0_lazy_cutoff']['branch']=='flat'
    assert all(row.zero for row in flat['q_axial_coefficients']+flat['q_squared_axial_coefficients'])
    cutoff=current.shear_q_jets(f,a,[f.scalar(c.mpf('.05'))]+a[1:],c.ln(c.mpf('.1')),proof)
    assert cutoff['original_C0_lazy_cutoff']['branch']=='active'
    assert not cutoff['analytic_axial_jets_available'] and cutoff['q_axial_coefficients'] is None
    try:current.zero_mode_sqrt(f.factor((0,.5,0,0,0)))
    except ValueError:rejects+=1
    else:raise AssertionError('Unverified nonzero source powers collapsed by square root')
    zero=[f.scalar(0)]*6;E=[f.scalar(1)]+zero[1:]
    op=SimpleNamespace(flow=f,c=c,P0=E,Pstar=f.scalar(3),Rm_factor=f.scalar(7),controls=[zero]*5)
    geometry=dict(point=True,exact_x=[5,4])
    source=dict(common_original_P0_axial5=op.P0,source_geometry=geometry,
        original_fixed_source_log_bases=f.logs,source_ledger_is_same_object=True,actual_patch_source_owner=type(op).__name__,
        actual_generic_source_numerators=dict(E=E,C=a,B=zero,positive_denominator=a),
        full_signed_inertial_sectors_axial4=dict(theta_linear=zero[:5],theta_quadratic=zero[:5],axial_linear=zero[:5],axial_quadratic=zero[:5]),
        original_physical_radius=f.scalar(c.mpf(35)/4))
    patch=dict(geometry=geometry,original_P0_normalized_axial5=op.P0,actual_H_x_derivative_axial5=[E],actual_gamma_ordinary_x_derivatives=[[c.mpf(0)]*5]*3)
    foreign=MacroFlow(c,c.mpf(-30),c.ln(9),c.ln(4),c.ln(5),c.ln(100)-c.ln(5),c.mpf('.01'))
    bad_sector={**source['full_signed_inertial_sectors_axial4'],'theta_linear':[foreign.scalar(0)]+zero[1:5]}
    binding_rejects=0
    for changed,changed_patch in (
        ({**source,'actual_patch_source_owner':'different owner'},patch),
        ({**source,'original_fixed_source_log_bases':foreign.logs},patch),
        ({**source,'full_signed_inertial_sectors_axial4':bad_sector},patch),
        ({**source,'original_physical_radius':f.scalar(-1)},patch),
        ({**source,'original_physical_radius':foreign.scalar(c.mpf(35)/4)},patch),
        ({**source,'original_physical_radius':f.scalar(9)},patch),
        (source,{**patch,'geometry':dict(point=True,exact_x=[3,2])})):
        try:current.recover_quotients(op,changed,c.ln(c.mpf('.1')),c.ln(c.mpf('.005')),patch=changed_patch)
        except ValueError:binding_rejects+=1
        else:raise AssertionError('Foreign owner, geometry, source algebra or radius accepted')
    return dict(passed=True,unproved_or_mismatched_source_denominators_rejected=rejects,
        owner_geometry_source_algebra_and_positive_radius_binding_rejections=binding_rejects,
        unverified_nonzero_factor_square_root_rejected=True,
        exact_flat_q_derivative_zero_checks=12,cutoff_derivatives_not_substituted_for_source_jets=True)


def correlated_shear_fixture(zvalue):
    x,z=sy.symbols('x Z',positive=True);R=sy.Rational
    change=sy.Function('f')(x,z);am=sy.Function('am')(z);H=x**R(1,10)+change;E=am*H
    identity=(E-2*x*sy.diff(E,x))/E-(R(4,5)-2*(x*sy.diff(change,x)-change/10)/H)
    assert sy.simplify(identity)==0
    c=MPIntervalContext();c.dps=200;xx=mp.mpf(5)/4;zz=mp.mpf(zvalue);tol=mp.mpf('1e-140')
    f=MacroFlow(c,c.mpf(-30),c.ln(9),c.ln(4),c.ln(5),c.ln(100)-c.ln(5),c.mpf('.01'))
    epsilon=R(1,1000)*(1+z+z**3);H=x**R(1,10)+epsilon*x**3
    expected=R(4,5)-2*(3*epsilon*x**3-epsilon*x**3/10)/H
    def jet(expr):
        values=[]
        for n in range(6):
            value=sy.lambdify((x,z),sy.diff(expr,z,n)/math.factorial(n),'mpmath')(xx,zz)
            values.append(c.mpf([value-tol,value+tol]))
        return f.jet(IntervalTaylor(c,values))
    zero=jet(sy.Integer(0));op=SimpleNamespace(flow=f,c=c,controls=[zero,zero,jet(epsilon),zero,zero])
    packet=dict(geometry=dict(point=True,exact_x=[5,4]),actual_H_x_derivative_axial5=[jet(H)],
        actual_gamma_ordinary_x_derivatives=[[c.mpf(xx)**3,3*c.mpf(xx)**2,6*c.mpf(xx),c.mpf(6),c.mpf(0)],
            [c.mpf(0)]*5,[c.mpf(0)]*5])
    actual,evidence=current.correlated_shear_a(op,packet)
    for n,row in enumerate(actual):
        value=sy.lambdify((x,z),sy.diff(expected,z,n)/math.factorial(n),'mpmath')(xx,zz)
        lo,hi=ep(row.coefficient*c.exp(row.scale.evaluate()));assert lo<=value<=hi
    assert evidence['common_axial_amplitude_canceled_analytically']
    return dict(passed=True,fixture_Z=zvalue,exact_amplitude_and_reference_shear_identity=True,
        independent_correlated_shear_axial_derivative_comparisons=6)


def native(owner,report):
    proofs=jets=P0=branches=units=small=terminal_zeros=0
    for label in ('0','.5'):
        for name,call in (('whole_patch',lambda:owner.cell(label,(1,1),'Rh')),
            ('first_support',lambda:owner.cell(label,(6,5),(13,10))),
            ('terminal',lambda:owner.cell(label,(71,40),'Rh')),
            ('active_point',lambda:owner.evaluate(label,(5,4)))):
            data=call();assert current.base.encoded(fields.serialized(data))==report['frames'][label][name]
            op=owner.upstream.upstream.owner(label).op;f=op.flow
            assert data['common_original_P0_axial5'] is op.P0;P0+=1
            for proof in data['source_positive_theorems'].values():
                assert ep(proof['source_row'].coefficient)[0]>0
                assert proof['strict_positive_complete_source_enclosure'] and proof['source_value_not_selected'];proofs+=1
            q=data['original_shear_q'];assert q['original_C0_lazy_cutoff']['branch']=='active'
            assert data['actual_correlated_shear_source']['common_axial_amplitude_canceled_analytically']
            assert q['analytic_axial_jets_available'] and ep(data['shear_quotients_axial5']['kappa_minus2'][0].coefficient)[1]<0;branches+=1
            assert len(q['q_axial_coefficients'])==len(q['q_squared_axial_coefficients'])==6
            assert q['q0_square_root_has_no_nontrivial_source_powers_to_collapse']
            assert q['q_squared_axial_coefficients'][0].scale.powers==(0,0,0,0,0)
            assert ep(q['q_axial_coefficients'][0].coefficient)[0]>0
            for rows in (data['shear_quotients_axial5']['a'],q['q_axial_coefficients'],q['q_squared_axial_coefficients']):
                for row in rows[1:]:
                    if not row.zero:assert row.scale.powers[1]<=-1;small+=1
                    if name=='terminal':assert row.zero;terminal_zeros+=1
            for rows in (*data['shear_quotients_axial5'].values(),q['q_axial_coefficients'],q['q_squared_axial_coefficients'],
                         *data['full_signed_inertial_quotients_before_shared_R_axial4'].values(),data['original_signed_u_before_shared_R_axial4']):
                assert all(row.scale.bases is f.logs and row.ledger is f.ledger for row in rows);jets+=len(rows)
            R=data['exact_same_shared_positive_radius_factor']
            assert data['exact_same_shared_radius_positive_theorem']['source_row'] is R
            assert ep(R.coefficient)[0]>0 and R.scale.bases is f.logs and R.ledger is f.ledger
            assert data['source_context_basis_ledger_geometry_and_P0_bound_to_same_actual_owner']
            for key,rows in data['full_signed_inertial_quotients_before_shared_R_axial4'].items():
                for n,row in enumerate(rows):
                    assert current.base.encoded(fields.serialized(row*R))==current.base.encoded(fields.serialized(data['actual_factored_p1_p2_axial4'][key][n]));units+=1
            assert data['inertial_quotients_differentiated_before_attaching_one_radius_factor']
            assert data['original_selected_parameters_do_not_transfer_old_owner_margin_theorems']
            assert not data['whole_axis_provider_or_new_dstar_cone_margin_certified']
            assert not data['actual_phase_inverse_or_finite_N_density_integrals_installed']
            assert all(data[key] is False for key in fields.previous.OPEN)
    return dict(passed=True,actual_positive_complete_source_denominator_theorems=proofs,
        actual_whole_cell_and_point_quotient_q_u_axial_coefficients=jets,
        exact_same_P0_object_checks=P0,certified_original_sigma1_active_branch_checks=branches,
        one_shared_radius_algebraic_prefactor_checks=units,
        actual_small_Pstar_inverse_squared_high_Z_source_modes=small,
        exact_terminal_a_q_q_squared_positive_Z_derivative_zero_checks=terminal_zeros,
        whole_radial_patch_at_two_conditional_axial_frames_only=True)


def run():
    began=time.monotonic()
    with mp.workdps(180):
        fixtures=[closed_fixture(z) for z in ('.3','-.4')]
        correlated=[correlated_shear_fixture(z) for z in ('.3','-.4')];guards=guard_checks()
    print('Closed quotient/q derivatives and source/cutoff guards PASS',flush=True)
    report=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    owner=current.OriginalRmPositiveQuotientsQ(require_checked=False);evidence=native(owner,report)
    assert report[current.GATE] and report['source_family']==owner.family
    for name,digest in report['input_hashes'].items():fields.previous.bind(owner.hashes,name,digest)
    for name in (current.NAME,Path(__file__).name):fields.previous.bind(owner.hashes,name,current.sha(name))
    result=dict(all_passed=True,**{current.GATE:True},source_family=owner.family,
        closed_source_quotient_q_signed_u_fixtures=fixtures,denominator_and_lazy_branch_guards=guards,
        exact_correlated_reference_shear_fixtures=correlated,
        actual_live_source=evidence,**dict.fromkeys(fields.previous.OPEN,False),
        input_hashes=owner.hashes,execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.base.encoded(result),indent=2)+'\n',encoding='utf8')
    print('Actual Rm whole radial positive quotients and original q PASS',flush=True);return result


if __name__=='__main__':run()
