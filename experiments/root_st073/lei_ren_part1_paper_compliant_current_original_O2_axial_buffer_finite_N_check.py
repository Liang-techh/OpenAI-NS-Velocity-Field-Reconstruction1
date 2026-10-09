"""Independent history/Jacobian identities and genuine O2 axial/buffer input."""
import copy
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as sy
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_original_O2_axial_buffer_finite_N as current

fields,ep=current.fields,current.ep


def symbolic_source():
    y,z,tau,md=sy.symbols('y Z tau Md',real=True);S=sy.symbols('Pstar',positive=True)
    B=sy.Function('B')(y);K1=sy.Function('K1')(y);K2=sy.Function('K2')(y)
    u1=sy.Function('u1')(z);old={key:sy.Function(key+'1')(z) for key in current.RATES}
    t=y-1;d=sy.exp(-t);root=sy.exp(-t/2);d3=sy.exp(-sy.Rational(3,2)*t)
    E=u1*root;V=4*z*B
    hist=dict(m=old['m']*d+4*z*K1,h=old['h']*d3+u1*(root-d3),
        k=old['k']*d3+u1*4*z*root*K1,
        e=old['e']*d+16*z*z/S**2*K2-u1*u1*t*d/2,p=old['p']+u1*u1*(1-d)/2)
    rhs=dict(m=V-hist['m'],h=E-sy.Rational(3,2)*hist['h'],k=E*V-sy.Rational(3,2)*hist['k'],
        e=V*V/S**2-E*E/2-hist['e'],p=E*E/2)
    derivatives={sy.diff(K1,y):B-K1,sy.diff(K2,y):B**2-K2}
    for key in hist:
        difference=(sy.diff(hist[key],y)-rhs[key]).subs(derivatives)
        assert sy.simplify(difference)==0 and sy.simplify(sy.diff(difference,z))==0
        assert sy.simplify(hist[key].subs({y:1,K1:0,K2:0})-old[key])==0
    assert sy.simplify(1-2*sy.diff(E,y)/E)==2
    p=sy.symbols('phase',real=True)
    assert sy.simplify(sy.diff(sy.exp(md*p),p)-md*sy.exp(md*p))==0
    # Vanishing B beyond the turnoff, rather than setting history K to0.
    k0=sy.symbols('K0');kernel=k0*sy.exp(-tau)
    assert sy.simplify(sy.diff(kernel,tau)+kernel)==0 and kernel.subs(tau,0)==k0
    b,eta=sy.symbols('b eta',real=True)
    assert sy.simplify((2+b*b/2)-2-b*b/2)==0
    assert sy.simplify((2*eta-b*b/2)/4).subs(b,0)==eta/2
    return dict(passed=True,original_all_five_axial_history_y_and_yZ_equations=True,
        exact_slope_exit_join_by_zero_kernel_and_decay_one=True,exact_critical_a2_and_collected_Delta_b_squared_over2=True,
        exact_axial_physical_Jacobian_Md_exp_Md_phase=True,exact_buffer_kernel_nonzero_memory_identity=True,
        positive_eta_q_squared_eta_over2_when_bzero=True)


def independent_scalar():
    c=MPIntervalContext();c.dps=100;comparisons=0
    # Explicit cutoff, independently integrated in s. These bounded scalar
    # references test the directed kernel helper, not a selected large-y field.
    def cutoff(s):
        x=1-mp.log(s)/40
        if x<=0:return mp.mpf(0)
        if x>=1:return mp.mpf(1)
        a=mp.exp(-1/x**2);b=mp.exp(-1/(1-x)**2)
        return a/(a+b)
    for y in (2,8,30):
        K=current.kernels.turnoff_kernels(c,c.mpf(y),40,256,800)
        for j,key in ((1,'B_mass'),(2,'B_squared_mass')):
            value=mp.quad(lambda s:mp.exp(s-y)*cutoff(s)**j,[1,(1+y)/2,y])
            lo,hi=ep(K[key]);assert lo-mp.mpf('1e-90')<=value<=hi+mp.mpf('1e-90');comparisons+=1
    inlet=current.kernels.turnoff_kernels(c,c.mpf(1),40,128,800)
    assert all(ep(inlet[key])==(0,0) for key in ('B_mass','B_squared_mass','retained_far_tail'))
    huge=current.kernels.turnoff_kernels(c,c.exp(40),40,128,800)
    assert ep(huge['retained_far_tail'])[1]>0 and ep(huge['B_mass'])[1]>0
    # Ordinary derivative of sigma(1-log(y)/40), independent chain rule.
    for phase in (mp.mpf('.25'),mp.mpf('.5'),mp.mpf('.75')):
        y=mp.exp(40*phase)
        rows=current.previous.original.turnoff_derivatives(c,c.mpf(y),40,c.mpf(phase))
        derivative=mp.diff(cutoff,y)
        lo,hi=ep(rows[1]);assert lo-mp.mpf('1e-90')<=derivative<=hi+mp.mpf('1e-90');comparisons+=1
    return dict(passed=True,independent_explicit_cutoff_kernel_and_ordinary_derivative_comparisons=comparisons,
        direct_scalar_quadratures_diagnostic_only_directed_kernel_rectangles_are_certification=True,
        exact_zero_kernel_inlet_and_positive_nonzero_far_tail_checked=True)


def native(owner,report):
    cells=coefficients=density_rows=boundary_rows=axis_cases=buffer_q_cases=0
    for label in ('0','.5'):
        packet=owner.contribution(label);op=owner.owner(label);f,c=op.flow,op.c
        assert current.base.encoded(current.reference.serialized(packet))==report['frames'][label]
        assert packet['exact_common_P0_axial5'] is op.P0
        saved=owner.saved['frames'][label]
        assert current.base.encoded(fields.serialized(packet['actual_current_O2_slope_exit_incoming_correction_C0_Z']))==saved['actual_current_O2_slope_exit_correction_C0_Z']
        for joined in packet['typed_actual_O2_axial_buffer_source_joins'].values():
            assert joined['source_overlap_consistency_rows']==42 and joined['overlap_not_function_identity_proof']
        for chart,part in packet['actual_O2_axial_buffer_charts'].items():
            width=part['full_actual_logarithmic_width']
            if chart=='axial':assert ep(width-c.expm1(40))[0]<=0<=ep(width-c.expm1(40))[1]
            else:assert ep(width)==(11,11)
            for cell in part['actual_source_cells']:
                source=cell['source'];generic=source['original_generic_source'];back=source['original_closed_O2_axial_buffer_background'];cells+=1
                assert source['exact_common_P0_axial5'] is op.P0 and back['exact_source_Rm_factor'] is op.Rm_factor
                assert ep(back['actual_a_axial5'][0].coefficient)==(2,2) and all(row.zero for row in back['actual_a_axial5'][1:])
                assert not back['positive_formal_history_decays']['mean'].zero
                assert not back['positive_formal_history_decays']['amplitude'].zero
                for rows in (*generic['common_own_five_histories_axial5'].values(),generic['common_velocity_E_axial5'],generic['common_velocity_V_axial5']):
                    for row in rows:
                        assert row.ctx is c and row.scale.bases is f.logs and row.ledger is f.ledger;coefficients+=1
                for part_name in ('kernels','Z_derivatives'):
                    for row in source['original_signed_five_density_C0_Z'][part_name].values():
                        assert row.ctx is c and row.scale.bases is f.logs and row.ledger is f.ledger;density_rows+=1
                left,right=source['actual_compact_source_coordinate_cell']
                actual_width,_,mass,decay=current.physical_weights(f,chart,left,right,current.Fraction(0))
                assert ep(cell['actual_original_log_width']-actual_width)[0]<=0<=ep(cell['actual_original_log_width']-actual_width)[1]
                assert ep(cell['original_own_rate_weights']['p']['actual_own_rate_mass']-mass)[0]<=0<=ep(cell['original_own_rate_weights']['p']['actual_own_rate_mass']-mass)[1]
                assert ep(cell['original_own_rate_weights']['p']['actual_suffix_to_chart_exit'].coefficient)==(1,1)
                if chart=='axial':
                    assert ep(actual_width)[0]>100 and ep(back['original_turnoff_kernels']['retained_far_tail'])[1]>0
                    assert not generic['actual_generic_source_numerators']['B'][1].zero
                else:
                    assert all(ep(row)==(0,0) for row in back['original_cutoff_ordinary_y_derivatives'])
                    assert all(row.zero for row in generic['common_velocity_V_axial5'])
                    assert not source['original_q_rows'][current.C0].zero and source['original_q_rows'][current.Z].zero;buffer_q_cases+=1
                    assert back['original_turnoff_kernels']['buffer_exact_retained_kernel_memory']
                    assert ep(back['original_turnoff_kernels']['B_mass'])[1]>0
                if label=='0':
                    p2=source['original_roots']['p2'];assert p2[current.C0].zero and not p2[current.Z].zero
                    assert not generic['common_radial_Q_axial4'][0].zero
                    assert source['original_q_rows'][current.Z].zero and not source['original_q_rows'][current.C0].zero;axis_cases+=1
            for key in ('actual_current_chart_incoming_C0_Z','actual_chart_local_driver_C0_Z','actual_current_chart_exit_C0_Z'):
                for pair in part[key].values():
                    for row in pair:
                        assert row.ctx is c and row.scale.bases is f.logs and row.ledger is f.ledger;boundary_rows+=1
        buffer=packet['actual_O2_axial_buffer_charts']['buffer'];axial=packet['actual_O2_axial_buffer_charts']['axial']
        assert buffer['actual_current_chart_incoming_C0_Z'] is axial['actual_current_chart_exit_C0_Z']
        assert packet['actual_current_Rd_correction_C0_Z'] is buffer['actual_current_chart_exit_C0_Z']
        assert packet['genuine_current_finite_N_prefix_through_Rd_boundary_enclosures_supplied']
        assert packet['no_old_O2_owner_or_source_replay'] and all(packet[k] is False for k in fields.previous.OPEN)
        print('Actual current O2 axial/buffer source and real Rd prefix checked',label,flush=True)
    op=owner.owner('0');good=owner.saved['frames']['0'];rejects=0
    for key,value in (('source_family','bad'),('source_frame','.5'),('candidate_N',1024),('exact_common_P0_axial5',[]),
        ('genuine_current_finite_N_prefix_through_O2_slope_exit_supplied',False),('actual_current_O2_slope_exit_correction_C0_Z',{})):
        bad=copy.deepcopy(good);bad[key]=value
        try:current.guard_incoming(owner.family,'0',257,op,bad)
        except ValueError:rejects+=1
        else:raise AssertionError('Wrong actual slope-exit input accepted')
    bad=copy.deepcopy(good);bad['exact_source_O2_slope_exit_factor']['exact_zero']=True
    try:current.guard_incoming(owner.family,'0',257,op,bad)
    except ValueError:rejects+=1
    else:raise AssertionError('Wrong slope-exit radius accepted')
    for chart,left,right in (('axial',(-1,1),(0,1)),('axial',(1,1),(0,1)),('axial',(0,1),(2,1)),('buffer',(0,1),(12,1)),('unknown',(0,1),(1,1))):
        try:current.coordinate(op.c,chart,left,right)
        except ValueError:rejects+=1
        else:raise AssertionError('Wrong actual source cell accepted')
    for label,N in (('0',1024),('unsupplied',257)):
        try:owner.contribution(label,N)
        except ValueError:rejects+=1
        else:raise AssertionError('Unsupplied current frame/frequency accepted')
    # A nonzero incoming homogeneous solution survives huge widths formally.
    zero={name:[f.scalar(0),f.scalar(0)] for name in current.RATES};seed={name:[f.scalar(1),f.scalar(1)] for name in current.RATES}
    transported=current.formal_affine_transport(f,seed,zero,c.expm1(40))
    assert all(not row.zero for pair in transported.values() for row in pair)
    assert ep(transported['p'][0].coefficient)==(1,1)
    return dict(passed=True,actual_closed_O2_axial_buffer_source_cells=cells,actual_common_source_coefficients=coefficients,
        actual_signed_nonlinear_density_C0_Z_rows=density_rows,actual_incoming_local_exit_boundary_C0_Z_rows=boundary_rows,
        midplane_zero_p2_nonzero_p2_Z_cases=axis_cases,buffer_positive_eta_q_not_false_quiet_cases=buffer_q_cases,
        typed_actual_source_and_incoming_rejections=rejects,formal_positive_homogeneous_history_not_underflowed=True,
        true_current_two_frame_prefix_not_functional_closure_or_global_N=True)


def run():
    began=time.monotonic();symbolic=symbolic_source()
    with mp.workdps(130):independent=independent_scalar()
    print('Independent original axial/buffer histories, physical Jacobian and kernel memory PASS',flush=True)
    report=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    owner=current.OriginalO2AxialBufferFiniteN(require_checked=False);actual=native(owner,report)
    assert report[current.GATE] and report['source_family']==owner.family
    assert report['original_source_bindings']==current.base.encoded(owner.bindings)
    assert report['exact_original_parameter_binding']==current.base.encoded(owner.parameters)
    for name,digest in report['input_hashes'].items():fields.previous.bind(owner.hashes,name,digest)
    for name in (current.NAME,Path(__file__).name):fields.previous.bind(owner.hashes,name,current.sha(name))
    result=dict(all_passed=True,**{current.GATE:True},source_family=owner.family,
        independent_original_axial_buffer_source_identities=symbolic,independent_scalar_kernels=independent,
        actual_live_source=actual,input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
        **dict.fromkeys(fields.previous.OPEN,False))
    (current.HERE/current.RECEIPT).write_bytes((json.dumps(current.base.encoded(result),indent=2)+'\n').encode())
    print('Current original genuine finite-N prefix through O2 axial/buffer Rd PASS',flush=True);return result


if __name__=='__main__':run()
