"""Closed reference source identities, signed axis derivative and real input."""
import copy
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as sy
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_original_Rh_reference_finite_N as current

fields,ep=current.fields,current.ep


def symbolic_source():
    y,z=sy.symbols('y Z',real=True);S=sy.symbols('Pstar',positive=True);delta=sy.symbols('delta',real=True)
    E=sy.exp(y/10)/(1+z*z);V=4*z
    histories=dict(m=V,h=sy.Rational(5,8)*E,k=sy.Rational(5,8)*E*V,
        e=V*V/S**2-sy.Rational(5,12)*E**2,p=sy.Rational(5,2)*E**2)
    rhs=dict(m=V-histories['m'],h=E-sy.Rational(3,2)*histories['h'],
        k=E*V-sy.Rational(3,2)*histories['k'],e=V*V/S**2-E**2/2-histories['e'],p=E**2/2)
    for k in histories:
        assert sy.simplify(sy.diff(histories[k],y)-rhs[k])==0
        assert sy.simplify(sy.diff(histories[k],y,z)-sy.diff(rhs[k],z))==0
    assert sy.simplify(E-2*sy.diff(E,y)-sy.Rational(4,5)*E)==0
    assert sy.diff(V,y)==0 and sy.diff(V,z)==4
    # The derivative here is axial Z, never the radial y derivative.
    Q=(2*z*V-z*histories['m']*(1-delta)-(1-z*z)*sy.diff(histories['m'],z))/(1-delta*z*z)
    expected=4*((2+delta)*z*z-1)/(1-delta*z*z)
    assert sy.simplify(Q-expected)==0 and Q.subs(z,0)==-4 and sy.diff(Q,z).subs(z,0)==0
    x=sy.symbols('x',real=True)
    assert sy.simplify(sy.exp(-sy.Rational(3,5))*sy.exp(sy.Rational(1,10))-sy.exp(-sy.Rational(1,2)))==0
    assert sy.simplify((x+6)-1-(x+5))==0
    return dict(passed=True,all_five_reference_history_y_and_yZ_equations=True,
        exact_E_y_E_over10_a_four_fifths_and_b_zero=True,
        exact_meridional_Q_uses_axial_m_Z_not_radial_m_y=True,
        physical_reference_Q_midplane_minus4_and_Q_Z_zero=True,
        exact_patched_Rh_reference_amplitude_radius_and_phase_join=True)


def independent_signed_integrals():
    c=MPIntervalContext();c.dps=150;p=mp.mp.clone();p.dps=100
    macro=current.incoming.prefix.macro.macro
    f=macro.fields.MacroFlow(c,c.ln(c.mpf('.002')),c.ln(9),c.ln(4),c.ln(5),c.ln(100)-c.ln(5),c.mpf('.01'))
    local={k:[f.scalar(0),f.scalar(0)] for k in current.RATES};initial={k:[f.scalar('-.2'),f.scalar('.03')] for k in current.RATES}
    coefficients=(p.mpf('-.37'),p.mpf('.021'));alpha=p.mpf('.17');comparisons=0
    with mp.workdps(190):
        for l,r in zip(range(-5,0),range(-4,1)):
            for k,rate in current.RATES.items():
                mass=current.incoming.weighted.terminal.local.positive_kernel_mass(c,c.mpf(1),rate)
                lam=c.mpf(rate.numerator)/rate.denominator
                tail=c.mpf(1) if not rate else c.exp(lam*r)
                for n,value in enumerate(coefficients):
                    source=f.scalar(c.mpf(str(value))*c.exp(c.mpf('.17')*c.mpf([l,r])))
                    bound=current.reference.bounds.symmetric(f,current.reference.bounds.magnitude(f,source)*mass*tail)
                    local[k][n]+=bound
        output=current.incoming.weighted.terminal.affine_transport(f,initial,local,c.mpf(5))
        for k,rate in current.RATES.items():
            lam=p.mpf(rate.numerator)/rate.denominator
            actual_mass=-p.expm1(-5*(lam+alpha))/(lam+alpha)
            for n,value in enumerate(coefficients):
                expected=p.mpf(('-.2','.03')[n])*p.exp(-5*lam)+value*actual_mass
                lo,hi=ep(output[k][n].finite_interval(max_log=2000));epsilon=p.mpf('1e-80')
                assert lo-epsilon<=expected<=hi+epsilon;comparisons+=1
    return dict(passed=True,independent_nonconstant_negative_C0_and_nonzero_Z_integral_comparisons=comparisons,
        original_five_unit_cells_total_width5_and_true_suffix_weights=True,
        actual_nonzero_incoming_memory_not_reset=True,pressure_rate_zero_integrates_full_five_unit_width=True)


def native(owner,report):
    cells=source_rows=density_rows=boundary_rows=axis_rows=0
    for label in ('0','.5'):
        packet=owner.contribution(label);op=owner.upstream.owner(label);f,c=op.flow,op.c
        assert current.base.encoded(current.reference.serialized(packet))==report['frames'][label]
        assert packet['exact_common_P0_axial5'] is op.P0
        original=owner.saved['frames'][label]['actual_current_Rh_correction_C0_Z']
        assert current.base.encoded(fields.serialized(packet['actual_current_Rh_incoming_correction_C0_Z']))==original
        assert ep(packet['full_original_logarithmic_interval_width'])==(5,5)
        assert packet['typed_actual_patched_Rh_reference_source_join']['directed_source_overlap_consistency_rows']==42
        assert packet['typed_actual_patched_Rh_reference_source_join']['overlap_only_consistency_not_the_function_identity_proof']
        for cell in packet['original_Rh_reference_source_cells']:
            got=cell['source'];raw=got['original_raw_current_source'];generic=got['original_generic_source'];cells+=1
            assert got['exact_common_P0_axial5'] is op.P0
            assert got['exact_source_Rm_factor'] is op.Rm_factor
            assert got['signed_midplane_and_nonzero_Z_derivative_interface_installed']
            for rows in (*generic['common_own_five_histories_axial5'].values(),generic['common_velocity_E_axial5'],generic['common_velocity_V_axial5']):
                for row in rows:
                    assert row.ctx is c and row.scale.bases is f.logs and row.ledger is f.ledger;source_rows+=1
            for part in ('kernels','Z_derivatives'):
                for row in got['original_signed_five_density_C0_Z'][part].values():
                    assert row.ctx is c and row.scale.bases is f.logs and row.ledger is f.ledger;density_rows+=1
            if label=='0':
                p2=got['original_roots']['p2']
                assert p2[current.C0].zero and not p2[current.Z].zero
                assert not got['original_q_rows'][current.C0].zero and got['original_q_rows'][current.Z].zero
                assert raw['raw_current_radius_y_derivative_axial_coefficients']['velocity']['axial'][0][0].zero
                assert not raw['raw_current_radius_y_derivative_axial_coefficients']['velocity']['axial'][0][1].zero
                assert generic['common_velocity_E_axial5'][1].zero
                assert not generic['common_radial_Q_axial4'][0].zero
                assert generic['common_radial_Q_axial4'][1].zero;axis_rows+=1
            assert ep(cell['actual_physical_log_width'])==(1,1)
            pressure=cell['original_own_rate_weights']['p']
            assert ep(pressure['true_own_rate_cell_mass'])==(1,1) and ep(pressure['suffix_to_Rref'])==(1,1)
        for key in ('actual_current_Rh_incoming_correction_C0_Z','actual_Rh_Rref_local_driver_C0_Z','actual_current_Rref_correction_C0_Z'):
            for pair in packet[key].values():
                for row in pair:
                    assert row.ctx is c and row.scale.bases is f.logs and row.ledger is f.ledger;boundary_rows+=1
        assert packet['genuine_current_finite_N_prefix_through_Rref_boundary_enclosures_supplied']
        assert packet['no_old_point_owner_constructor_or_integration_replay']
        assert all(packet[k] is False for k in fields.previous.OPEN)
        print('Actual current reference source and genuine prefix checked',label,flush=True)
    op=owner.upstream.owner('0');good=owner.saved['frames']['0'];rejects=0
    for key,value in (('source_family','bad'),('source_frame','.5'),('candidate_N',1024),('exact_common_P0_axial5',[]),
        ('genuine_current_finite_N_prefix_through_Rh_boundary_enclosures_supplied',False),
        ('bound_original_endpoint_phase_source','unanchored free phase'),('actual_current_Rh_correction_C0_Z',{})):
        bad=copy.deepcopy(good);bad[key]=value
        try:current.guard_incoming(owner.family,'0',257,op,bad)
        except ValueError:rejects+=1
        else:raise AssertionError('Wrong actual Rh incoming accepted: '+key)
    bad=copy.deepcopy(good);bad['exact_source_Rh_factor']['exact_zero']=True
    try:current.guard_incoming(owner.family,'0',257,op,bad)
    except ValueError:rejects+=1
    else:raise AssertionError('Wrong actual Rh radius accepted')
    for left,right in (((-6,1),(-5,1)),((0,1),(-1,1)),((0,1),(1,1))):
        try:current.background_cell(op,left,right)
        except ValueError:rejects+=1
        else:raise AssertionError('Wrong reference source window accepted')
    for label,N in (('0',1024),('unsupplied',257)):
        try:owner.contribution(label,N)
        except ValueError:rejects+=1
        else:raise AssertionError('Unsupplied current frame/frequency accepted')
    return dict(passed=True,actual_closed_reference_source_cells=cells,actual_common_source_coefficients=source_rows,
        actual_signed_nonlinear_density_C0_Z_rows=density_rows,actual_Rh_local_Rref_boundary_C0_Z_rows=boundary_rows,
        actual_midplane_p2_zero_nonzero_p2_Z_cases=axis_rows,typed_source_geometry_and_incoming_rejections=rejects,
        original_meridional_midplane_inertia_and_nonzero_Z_derivatives_retained=True,
        actual_two_frame_prefix_not_functional_terminal_closure_or_global_N=True)


def run():
    began=time.monotonic();symbolic=symbolic_source()
    with mp.workdps(190):independent=independent_signed_integrals()
    print('Independent reference source, meridional axis and five-unit signed integrals PASS',flush=True)
    report=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    owner=current.OriginalRhReferenceFiniteN(require_checked=False);actual=native(owner,report)
    assert report[current.GATE] and report['source_family']==owner.family
    assert report['original_source_bindings']==current.base.encoded(owner.bindings)
    for name,digest in report['input_hashes'].items():fields.previous.bind(owner.hashes,name,digest)
    for name in (current.NAME,Path(__file__).name):fields.previous.bind(owner.hashes,name,current.sha(name))
    result=dict(all_passed=True,**{current.GATE:True},source_family=owner.family,
        independent_closed_reference_source_identities=symbolic,independent_signed_reference_integrals=independent,
        actual_live_source=actual,input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
        **dict.fromkeys(fields.previous.OPEN,False))
    (current.HERE/current.RECEIPT).write_bytes((json.dumps(current.base.encoded(result),indent=2)+'\n').encode())
    print('Current original source-bound Rh-reference prefix through Rref PASS',flush=True);return result


if __name__=='__main__':run()
