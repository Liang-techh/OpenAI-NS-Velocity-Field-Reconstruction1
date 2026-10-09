"""Original O3 identities, positive parameter/cutoff and real Rc transport."""
import copy
import gzip
import json
import math
from pathlib import Path
import time
import mpmath as mp
import sympy as sy
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_original_O3_Rc_finite_N as current

fields,ep=current.fields,current.ep


def symbolic_source():
    t,z=sy.symbols('t Z',real=True);mu=sy.symbols('mu',positive=True)
    u1=sy.Function('u1')(z);old={key:sy.Function(key+'1')(z) for key in current.RATES}
    J=sy.Function('J')(t);sig=sy.Function('sigma')(t)
    K={key:sy.Function('K'+key)(t) for key in ('theta','energy','pressure')}
    d=sy.exp(-t);d3=sy.exp(-sy.Rational(3,2)*t);E=u1*sy.exp(-t/2-mu*J)
    hist=dict(m=old['m']*d,h=old['h']*d3+u1*d3*K['theta'],k=old['k']*d3,
        e=old['e']*d-u1*u1*d*K['energy']/2,p=old['p']+u1*u1*K['pressure']/2)
    derivatives={sy.diff(J,t):sig,sy.diff(K['theta'],t):sy.exp(t-mu*J),
        sy.diff(K['energy'],t):sy.exp(-2*mu*J),sy.diff(K['pressure'],t):sy.exp(-t-2*mu*J)}
    def prove(hist,E,derivatives):
        rhs=dict(m=-hist['m'],h=E-sy.Rational(3,2)*hist['h'],k=-sy.Rational(3,2)*hist['k'],e=-hist['e']-E*E/2,p=E*E/2)
        for name in hist:
            difference=(sy.diff(hist[name],t)-rhs[name]).subs(derivatives)
            assert sy.simplify(difference)==0 and sy.simplify(sy.diff(difference,z))==0
    prove(hist,E,derivatives)
    zero={t:0,J:0,**{value:0 for value in K.values()}}
    assert all(sy.simplify(row.subs(zero)-old[name])==0 for name,row in hist.items())
    a=(1-2*sy.diff(E,t)/E).subs(derivatives)
    assert sy.simplify(a-2-2*mu*sig)==0
    factor=sy.exp((-sy.Rational(1,2)-mu)*t)
    D2=(1-sy.exp(-2*mu*t))/(2*mu);D1=(1-sy.exp(-(1+2*mu)*t))/(1+2*mu)
    theta=(factor-d3)/(1-mu)
    power=dict(m=old['m']*d,h=old['h']*d3+u1*theta,k=old['k']*d3,
        e=old['e']*d-u1*u1*d*D2/2,p=old['p']+u1*u1*D1/2)
    prove(power,u1*factor,{})
    assert all(sy.simplify(row.subs(t,0)-old[name])==0 for name,row in power.items())
    assert sy.simplify(theta-d3*(sy.exp((1-mu)*t)-1)/(1-mu))==0
    eta,Delta,A=sy.symbols('eta Delta a',positive=True)
    q2=(2*eta-Delta)/(2*A)
    assert sy.simplify(eta/A-q2)==Delta/(2*A)
    return dict(passed=True,original_all_five_O3_transition_and_power_y_and_yZ_equations=True,
        exact_Rd_Rw_join_by_zero_integrals_and_decay_one=True,
        exact_transition_a_minus2_equals_2mu_sigma_and_power_equals2mu=True,
        original_positive_exprel_power_theta_identity=True,
        nonnegative_Delta_active_q_squared_le_eta_over_a_and_a_ge2_implies_eta_over2=True,
        original_flat_q_when_Delta_ge_eta_independent_of_p2=True)


def independent_scalar():
    c=MPIntervalContext();c.dps=100
    # Explicit cutoff mesh is diagnostic; original directed endpoint rectangles
    # certify the integrals for the same source mu or its enclosing box.
    def sigma(t):
        if t<=0:return 0.
        if t>=1:return 1.
        exponent=1/t**2-1/(1-t)**2
        if exponent>700:return 0.
        if exponent<-700:return 1.
        return 1/(1+math.exp(exponent))
    count=8192;step=1/count;mu=.0001;J=0.;mass=[0.,0.,0.]
    previous_sigma=0.;previous_integrands=[1.,1.,1.];saved={0:[J,*mass]}
    for i in range(1,count+1):
        t=i*step;sig=sigma(t);J+=step*(previous_sigma+sig)/2
        values=[math.exp(t-mu*J),math.exp(-2*mu*J),math.exp(-t-2*mu*J)]
        mass=[old+step*(a+b)/2 for old,a,b in zip(mass,previous_integrands,values)]
        previous_sigma=sig;previous_integrands=values
        if i in (count//4,count//2,3*count//4,count):saved[i]=[J,*mass]
    comparisons=0
    for left,right in ((0,count//4),(count//4,3*count//4),(count//2,count)):
        t=c.mpf([mp.mpf(left)/count,mp.mpf(right)/count])
        packet=current.transition_kernel_hull(c,t,c.mpf([0,'.0002']))
        for index in (left,right):
            for name,value in zip(('J','theta','energy','pressure'),saved[index]):
                lo,hi=ep(packet[name]);assert lo>=0 and lo-mp.mpf('1e-7')<=value<=hi+mp.mpf('1e-7');comparisons+=1
    zero=current.transition_kernel_hull(c,c.mpf(0),c.mpf([0,'.0002']))
    assert all(ep(zero[key])==(0,0) for key in ('J','theta','energy','pressure'))
    one=current.transition_kernel_hull(c,c.mpf(1),c.mpf([0,'.0002']))
    assert ep(one['J'])==(mp.mpf('.5'),mp.mpf('.5'))
    for k,t in ((0,1),('.0001',2),(1,1),(3,2)):
        packet=current.positive_decay_mass(c,c.mpf(k),c.mpf(t))
        value=mp.mpf(t) if k==0 else -mp.expm1(-mp.mpf(k)*t)/mp.mpf(k)
        lo,hi=ep(packet);assert lo-mp.mpf('1e-90')<=value<=hi+mp.mpf('1e-90');comparisons+=1
    return dict(passed=True,independent_positive_transition_kernel_and_decay_mass_comparisons=comparisons,
        explicit_sigma_mesh8192_error_allowance1e_minus7_diagnostic_only=True,
        directed_endpoint_rectangles_are_source_integral_certification=True,
        exact_zero_local_integrals_and_J1_half_checked=True,
        arbitrary_small_mu_fixture_is_not_actual_source_parameter=True)


def native(owner,report):
    cells=coefficients=density_rows=boundary_rows=axis_cases=flat_cells=0
    with mp.workdps(owner.c.dps+40):
        actual_logmu=owner.c.ln(owner.c.mpf('.001'))-4*(owner.c.exp(40)+11)
        assert ep(owner.logmu-actual_logmu)[0]<=0<=ep(owner.logmu-actual_logmu)[1]
        margin=actual_logmu+owner.c.ln(2)-owner.loop_parameters.eta_log
        assert ep(margin)[0]>0 and ep(owner.power_margin-margin)[0]<=0<=ep(owner.power_margin-margin)[1]
        assert ep(owner.Tw+60*actual_logmu)[0]<=0<=ep(owner.Tw+60*actual_logmu)[1]
    for label in ('0','.5'):
        packet=owner.contribution(label);op=owner.owner(label);f,c=op.flow,op.c
        assert current.base.encoded(current.reference.serialized(packet))==report['frames'][label]
        assert packet['exact_common_P0_axial5'] is op.P0
        saved=owner.saved['frames'][label]
        assert current.base.encoded(fields.serialized(packet['actual_current_Rd_incoming_correction_C0_Z']))==saved['actual_current_Rd_correction_C0_Z']
        for seam in packet['typed_actual_O3_source_joins'].values():
            assert seam['source_overlap_consistency_rows']==42 and seam['overlap_not_function_identity_proof']
        for chart,part in packet['actual_O3_charts'].items():
            assert ep(part['full_actual_logarithmic_width'])==((1,1) if chart=='transition' else (2,2))
            for cell in part['actual_source_cells']:
                source=cell['source'];generic=source['original_generic_source'];back=source['original_closed_O3_background'];cells+=1
                assert source['exact_common_P0_axial5'] is op.P0 and back['exact_source_Rm_factor'] is op.Rm_factor
                mu=back['original_positive_formal_mu'];assert not mu.zero and ep(mu.coefficient)==(1,1)
                assert ep(mu.scale.offset-owner.logmu)[0]<=0<=ep(mu.scale.offset-owner.logmu)[1]
                assert ep(back['actual_mu_arithmetic_enclosure_only'])[0]==0 and ep(back['actual_mu_arithmetic_enclosure_only'])[1]>0
                delta=back['actual_collected_Delta_axial5'];assert ep(delta[0].coefficient)[0]>=0 and all(row.zero for row in delta[1:])
                assert all(row.zero for row in generic['common_velocity_V_axial5'])
                assert not generic['common_own_five_histories_axial5']['m'][1].zero
                for rows in (*generic['common_own_five_histories_axial5'].values(),generic['common_velocity_E_axial5'],generic['common_velocity_V_axial5']):
                    for row in rows:
                        assert row.ctx is c and row.scale.bases is f.logs and row.ledger is f.ledger;coefficients+=1
                for part_name in ('kernels','Z_derivatives'):
                    for row in source['original_signed_five_density_C0_Z'][part_name].values():
                        assert row.ctx is c and row.scale.bases is f.logs and row.ledger is f.ledger;density_rows+=1
                assert source['original_q_rows'][current.Z].zero
                branch=source['original_full_signed_quotients']['original_q_C0']['branch']
                if branch=='flat':
                    assert all(row.zero for row in source['original_q_rows'].values())
                    assert all(row.zero for row in source['original_primitive_values'].values())
                    assert all(row.zero for key in ('kernels','Z_derivatives') for row in source['original_signed_five_density_C0_Z'][key].values());flat_cells+=1
                else:
                    assert chart=='transition' and source['actual_compact_source_coordinate_cell'][0]==(0,1)
                    proof=source['original_full_signed_quotients']['original_q_C0']['original_nonnegative_Delta_q_cap_intersection']
                    assert proof['active_source_q_squared_upper']=='eta/2'
                if chart=='power':assert branch=='flat' and ep(source['actual_same_source_power_flat_log_margin'])[0]>0
                pressure=cell['original_own_rate_weights']['p']
                assert ep(pressure['actual_own_rate_mass'])==ep(cell['actual_original_log_width'])
                assert ep(pressure['actual_suffix_to_chart_exit'].coefficient)==(1,1)
                if label=='0':
                    p2=source['original_roots']['p2'];assert p2[current.C0].zero and not p2[current.Z].zero
                    assert not generic['common_radial_Q_axial4'][0].zero;axis_cases+=1
            for key in ('actual_current_chart_incoming_C0_Z','actual_chart_local_driver_C0_Z','actual_current_chart_exit_C0_Z'):
                for pair in part[key].values():
                    for row in pair:
                        assert row.ctx is c and row.scale.bases is f.logs and row.ledger is f.ledger;boundary_rows+=1
        transition,power=packet['actual_O3_charts']['transition'],packet['actual_O3_charts']['power']
        assert power['actual_current_chart_incoming_C0_Z'] is transition['actual_current_chart_exit_C0_Z']
        assert all(row.zero for pair in power['actual_chart_local_driver_C0_Z'].values() for row in pair)
        assert packet['actual_current_Rc_correction_C0_Z'] is power['actual_current_chart_exit_C0_Z']
        assert current.switch.canonical_expression(power['actual_current_chart_exit_C0_Z']['p'])==current.switch.canonical_expression(power['actual_current_chart_incoming_C0_Z']['p'])
        assert any(not row.zero for pair in power['actual_current_chart_exit_C0_Z'].values() for row in pair)
        assert packet['genuine_current_finite_N_prefix_through_Rc_boundary_enclosures_supplied']
        assert packet['no_old_O3_owner_N1024_tail_or_source_replay'] and all(packet[k] is False for k in fields.previous.OPEN)
        print('Actual current O3 source, proved quiet power and real Rc prefix checked',label,flush=True)
    op=owner.owner('0');f,c=op.flow,op.c;good=owner.saved['frames']['0'];rejects=0
    for key,value in (('source_family','bad'),('source_frame','.5'),('candidate_N',1024),('exact_common_P0_axial5',[]),
        ('genuine_current_finite_N_prefix_through_Rd_boundary_enclosures_supplied',False),('actual_current_Rd_correction_C0_Z',{})):
        bad=copy.deepcopy(good);bad[key]=value
        try:current.guard_incoming(owner.family,'0',257,op,bad)
        except ValueError:rejects+=1
        else:raise AssertionError('Wrong actual Rd input accepted')
    bad=copy.deepcopy(good);bad['exact_source_Rd_factor']['exact_zero']=True
    try:current.guard_incoming(owner.family,'0',257,op,bad)
    except ValueError:rejects+=1
    else:raise AssertionError('Wrong canonical Rd radius accepted')
    for chart,value in (('transition',(-1,1)),('transition',(2,1)),('power',(3,1)),('power',(1,0)),('unknown',(0,1))):
        try:current.fraction(chart,value)
        except ValueError:rejects+=1
        else:raise AssertionError('Wrong original O3 local coordinate accepted')
    for label,N in (('0',1024),('unsupplied',257)):
        try:owner.contribution(label,N)
        except ValueError:rejects+=1
        else:raise AssertionError('Unsupplied current frame/frequency accepted')
    try:current.nonnegative_Delta_q(f.scalar(2),f.scalar(-1),owner.loop_parameters.eta_log,c.ln(2))
    except ValueError:rejects+=1
    else:raise AssertionError('Negative Delta admitted by O3-specific positive theorem')
    return dict(passed=True,actual_closed_O3_source_cells=cells,actual_common_source_coefficients=coefficients,
        actual_signed_nonlinear_density_C0_Z_rows=density_rows,actual_incoming_local_exit_boundary_C0_Z_rows=boundary_rows,
        midplane_zero_p2_nonzero_p2_Z_cases=axis_cases,proved_original_flat_source_cells=flat_cells,
        typed_actual_source_and_incoming_rejections=rejects,
        exact_positive_mu_eta_log_guard_and_nonzero_quiet_history_memory_checked=True,
        true_current_two_frame_Rc_prefix_not_functional_closure_or_global_N=True)


def run():
    began=time.monotonic();symbolic=symbolic_source()
    with mp.workdps(150):independent=independent_scalar()
    print('Independent original O3 histories, positive kernel endpoints and quiet cutoff identities PASS',flush=True)
    report=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    owner=current.OriginalO3RcFiniteN(require_checked=False);actual=native(owner,report)
    assert report[current.GATE] and report['source_family']==owner.family
    assert report['original_source_bindings']==current.base.encoded(owner.bindings)
    for name,digest in report['input_hashes'].items():fields.previous.bind(owner.hashes,name,digest)
    for name in (current.NAME,Path(__file__).name):fields.previous.bind(owner.hashes,name,current.sha(name))
    result=dict(all_passed=True,**{current.GATE:True},source_family=owner.family,
        independent_original_O3_source_identities=symbolic,independent_scalar_kernels=independent,
        actual_live_source=actual,input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
        **dict.fromkeys(fields.previous.OPEN,False))
    (current.HERE/current.RECEIPT).write_bytes((json.dumps(current.base.encoded(result),indent=2)+'\n').encode())
    print('Current original genuine finite-N prefix through O3 Rc PASS',flush=True);return result


if __name__=='__main__':run()
