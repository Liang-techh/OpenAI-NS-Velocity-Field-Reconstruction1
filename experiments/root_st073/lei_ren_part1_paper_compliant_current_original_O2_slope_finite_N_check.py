"""Original scalar mass order, O2 history identities and true current input."""
import copy
import gzip
import json
import math
from pathlib import Path
import time
import mpmath as mp
import sympy as sy
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_original_O2_slope_finite_N as current

fields,ep=current.fields,current.ep


def symbolic_source():
    y,z=sy.symbols('y Z',real=True);S=sy.symbols('Pstar',positive=True)
    J=sy.Function('J')(y);Mh=sy.Function('Mh')(y);Mp=sy.Function('Mp')(y);Me=sy.Function('Me')(y);sig=sy.Function('sigma')(y)
    derivative={sy.diff(J,y):sig,sy.diff(Mh,y):sy.exp(sy.Rational(8,5)*y-sy.Rational(3,5)*J),
        sy.diff(Mp,y):sy.exp(y/5-sy.Rational(6,5)*J),sy.diff(Me,y):sy.exp(sy.Rational(6,5)*y-sy.Rational(6,5)*J)}
    qi=1/(1+z*z);E=qi*sy.exp(y/10-sy.Rational(3,5)*J);V=4*z
    h=qi*(sy.Rational(5,8)+Mh)*sy.exp(-sy.Rational(3,2)*y)
    hist=dict(m=V,h=h,k=h*V,e=V*V/S**2-qi**2*(sy.Rational(5,12)+Me/2)*sy.exp(-y),p=qi**2*(sy.Rational(5,2)+Mp/2))
    rhs=dict(m=V-hist['m'],h=E-sy.Rational(3,2)*h,k=E*V-sy.Rational(3,2)*hist['k'],e=V*V/S**2-E**2/2-hist['e'],p=E**2/2)
    for key in hist:
        difference=(sy.diff(hist[key],y)-rhs[key]).subs(derivative)
        assert sy.simplify(difference)==0 and sy.simplify(sy.diff(difference,z))==0
    assert sy.simplify((1-2*sy.diff(E,y)/E).subs(derivative)-(sy.Rational(4,5)+sy.Rational(6,5)*sig))==0
    eta=sy.symbols('eta',positive=True);a=sy.symbols('a',positive=True)
    q_squared=(2*eta-(a-2))/(2*a)
    assert sy.simplify(q_squared.subs(a,2)-eta/2)==0
    return dict(passed=True,original_all_five_slope_history_y_and_yZ_equations=True,
        original_J_theta_pressure_energy_mass_order_verified=True,variable_shear_a_four_fifths_plus_six_fifths_sigma=True,
        exact_Rref_join_by_zero_J_and_masses=True,positive_eta_endpoint_q_squared_eta_over2_not_flat_zero=True)


def independent_scalar_integrals():
    c=MPIntervalContext();c.dps=100
    # Independent explicit sigma and mesh quadrature are diagnostics only.
    # Certification remains the original monotone directed rectangles.
    def sigma(x):
        if x<=0:return 0.
        if x>=1:return 1.
        exponent=1/(x*x)-1/((1-x)*(1-x))
        if exponent>700:return 0.
        if exponent<-700:return 1.
        return 1/(1+math.exp(exponent))
    count=8192;step=1/count;J=0.;mass=[0.,0.,0.];old_sigma=0.;old_integrands=[1.,1.,1.];snapshots={}
    for i in range(1,count+1):
        y=i*step;now_sigma=sigma(y);J+=step*(old_sigma+now_sigma)/2
        integrands=[math.exp(rate*y-.6*power*J) for rate,power in ((1.6,1),(.2,2),(1.2,2))]
        mass=[v+step*(l+r)/2 for v,l,r in zip(mass,old_integrands,integrands)]
        old_sigma=now_sigma;old_integrands=integrands
        if i in (count//4,count//2,count):snapshots[i]=[J,*mass]
    comparisons=0
    for i,values in snapshots.items():
        y=c.mpf(i)/count
        Jbox,massbox=current.original.slope_masses(c,y,current.SCALAR_CELLS)
        for box,value in zip([Jbox,*massbox],values):
            lo,hi=ep(box);assert lo-mp.mpf('1e-7')<=value<=hi+mp.mpf('1e-7');comparisons+=1
        sig=current.original.sigma_jets(c,y)[0]
        lo,hi=ep(sig);assert lo-mp.mpf('1e-12')<=sigma(i/count)<=hi+mp.mpf('1e-12');comparisons+=1
    J0,M0=current.original.slope_masses(c,c.mpf(0),current.SCALAR_CELLS)
    J1,M1=current.original.slope_masses(c,c.mpf(1),current.SCALAR_CELLS)
    assert ep(J0)==(0,0) and all(ep(v)==(0,0) for v in M0)
    assert ep(J1)==(mp.mpf('.5'),mp.mpf('.5'))
    return dict(passed=True,independent_scalar_sigma_J_and_three_mass_comparisons=comparisons,
        numeric_mesh8192_error_allowance1e_minus7_diagnostic_only=True,
        original_directed_rectangles_are_certification_not_numeric_mesh=True,
        exact_zero_inlet_and_symmetric_J1_half_verified=True)


def native(owner,report):
    cells=coefficients=density_rows=boundary_rows=axis_cases=0
    for label in ('0','.5'):
        packet=owner.contribution(label);op=owner.owner(label);f,c=op.flow,op.c
        assert current.base.encoded(current.reference.serialized(packet))==report['frames'][label]
        assert packet['exact_common_P0_axial5'] is op.P0
        saved=owner.saved['frames'][label]
        assert current.base.encoded(fields.serialized(packet['actual_current_Rref_incoming_correction_C0_Z']))==saved['actual_current_Rref_correction_C0_Z']
        assert ep(packet['full_original_logarithmic_interval_width'])==(1,1)
        assert packet['typed_actual_Rref_O2_source_join']['source_overlap_consistency_rows']==42
        for cell in packet['actual_O2_slope_source_cells']:
            source=cell['source'];generic=source['original_generic_source'];back=source['original_closed_O2_slope_background'];cells+=1
            assert source['exact_common_P0_axial5'] is op.P0 and back['exact_source_Rm_factor'] is op.Rm_factor
            assert back['original_scalar_partition']==current.SCALAR_CELLS
            for rows in (*generic['common_own_five_histories_axial5'].values(),generic['common_velocity_E_axial5'],generic['common_velocity_V_axial5']):
                for row in rows:
                    assert row.ctx is c and row.scale.bases is f.logs and row.ledger is f.ledger;coefficients+=1
            for part in ('kernels','Z_derivatives'):
                for row in source['original_signed_five_density_C0_Z'][part].values():
                    assert row.ctx is c and row.scale.bases is f.logs and row.ledger is f.ledger;density_rows+=1
            a=back['actual_a_axial5']
            with mp.workdps(c.dps+40):
                tolerance=mp.mpf('1e-480');lo,hi=ep(a[0].coefficient)
                assert mp.mpf(4)/5-tolerance<=lo<=hi<=2+tolerance
            assert all(row.zero for row in a[1:])
            assert ep(cell['actual_original_log_width'])==(mp.mpf('.25'),mp.mpf('.25'))
            pressure=cell['original_own_rate_weights']['p']
            assert ep(pressure['actual_own_rate_mass'])==(mp.mpf('.25'),mp.mpf('.25')) and ep(pressure['actual_suffix_to_slope_exit'])==(1,1)
            if label=='0':
                p2=source['original_roots']['p2'];assert p2[current.C0].zero and not p2[current.Z].zero
                assert not generic['common_radial_Q_axial4'][0].zero and generic['common_radial_Q_axial4'][1].zero
                assert source['original_q_rows'][current.Z].zero;axis_cases+=1
        endpoint=packet['actual_O2_slope_exit_background_source'];back=endpoint['original_closed_O2_slope_background']
        assert ep(back['original_J'])==(mp.mpf('.5'),mp.mpf('.5')) and ep(back['original_sigma'])==(1,1)
        assert ep(back['actual_a_axial5'][0].coefficient)==(2,2)
        # Original q(a=2)^2=eta/2 is positive; never force the endpoint flat.
        assert not endpoint['original_q_rows'][current.C0].zero and endpoint['original_q_rows'][current.Z].zero
        for key in ('actual_current_Rref_incoming_correction_C0_Z','actual_O2_slope_local_driver_C0_Z','actual_current_O2_slope_exit_correction_C0_Z'):
            for pair in packet[key].values():
                for row in pair:
                    assert row.ctx is c and row.scale.bases is f.logs and row.ledger is f.ledger;boundary_rows+=1
        assert packet['genuine_current_finite_N_prefix_through_O2_slope_exit_supplied']
        assert packet['no_old_O2_owner_or_source_replay'] and all(packet[k] is False for k in fields.previous.OPEN)
        print('Actual current O2 slope source and genuine prefix checked',label,flush=True)
    label='0';op=owner.owner(label);good=owner.saved['frames'][label];rejects=0
    for key,value in (('source_family','bad'),('source_frame','.5'),('candidate_N',1024),('exact_common_P0_axial5',[]),
        ('genuine_current_finite_N_prefix_through_Rref_boundary_enclosures_supplied',False),('actual_current_Rref_correction_C0_Z',{})):
        bad=copy.deepcopy(good);bad[key]=value
        try:current.guard_incoming(owner.family,label,257,op,bad)
        except ValueError:rejects+=1
        else:raise AssertionError('Wrong actual Rref input accepted')
    bad=copy.deepcopy(good);bad['exact_source_Rref_factor']['exact_zero']=True
    try:current.guard_incoming(owner.family,label,257,op,bad)
    except ValueError:rejects+=1
    else:raise AssertionError('Wrong Rref radius accepted')
    for left,right,count in (((-1,1),(0,1),128),((1,1),(0,1),128),((0,1),(2,1),128),((0,1),(1,1),0)):
        try:current.background_cell(op,left,right,count)
        except ValueError:rejects+=1
        else:raise AssertionError('Wrong original slope cell or scalar count accepted')
    for label,N in (('0',1024),('unsupplied',257)):
        try:owner.contribution(label,N)
        except ValueError:rejects+=1
        else:raise AssertionError('Unsupplied current frame/frequency accepted')
    return dict(passed=True,actual_closed_O2_slope_source_cells=cells,actual_common_source_coefficients=coefficients,
        actual_signed_nonlinear_density_C0_Z_rows=density_rows,actual_Rref_local_slope_exit_boundary_C0_Z_rows=boundary_rows,
        midplane_zero_p2_nonzero_p2_Z_cases=axis_cases,typed_actual_source_and_incoming_rejections=rejects,
        original_endpoint_positive_eta_q_and_Z_rows_retained=True,
        true_current_two_frame_prefix_not_functional_closure_or_global_N=True)


def run():
    began=time.monotonic();symbolic=symbolic_source()
    with mp.workdps(150):independent=independent_scalar_integrals()
    print('Independent original O2 history identities and directed scalar reference comparisons PASS',flush=True)
    report=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    owner=current.OriginalO2SlopeFiniteN(require_checked=False);actual=native(owner,report)
    assert report[current.GATE] and report['source_family']==owner.family
    assert report['original_source_bindings']==current.base.encoded(owner.bindings)
    for name,digest in report['input_hashes'].items():fields.previous.bind(owner.hashes,name,digest)
    for name in (current.NAME,Path(__file__).name):fields.previous.bind(owner.hashes,name,current.sha(name))
    result=dict(all_passed=True,**{current.GATE:True},source_family=owner.family,
        independent_original_O2_slope_source_identities=symbolic,independent_scalar_integrals=independent,
        actual_live_source=actual,input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
        **dict.fromkeys(fields.previous.OPEN,False))
    (current.HERE/current.RECEIPT).write_bytes((json.dumps(current.base.encoded(result),indent=2)+'\n').encode())
    print('Current original genuine finite-N prefix through O2 slope exit PASS',flush=True);return result


if __name__=='__main__':run()
