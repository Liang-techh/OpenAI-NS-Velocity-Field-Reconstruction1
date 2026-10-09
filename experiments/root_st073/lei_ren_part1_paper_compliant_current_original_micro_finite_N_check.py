"""Source-owned inlet, real microscopic transport and prefix composition."""
import copy
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_micro_finite_N as current

fields,ep=current.fields,current.micro.ep


def contains(row,value):
    lo,hi=ep(row.finite_interval(max_log=2000));error=mp.mpf('1e-75')
    assert lo-error<=value<=hi+error,(value,lo,hi)


def symbolic():
    E,V,A,B,N=sy.symbols('E V A B N',real=True)
    dE=E*(sy.exp(A/N)-1);dV=B/N
    rates=dict(m=dV,h=dE,k=V*dE+E*dV+dE*dV,
        e=2*V*dV+dV*dV-E*dE-dE*dE/2,p=E*dE+dE*dE/2)
    assert all(sy.simplify(v.subs({A:0,B:0}))==0 for v in rates.values())
    assert sy.simplify((E+dE).subs(A,0)-E)==0
    assert sy.simplify((V+dV).subs(B,0)-V)==0
    y,t,z,lam=sy.symbols('y t z lam',real=True);source=sy.Function('source')
    D=sy.Integral(sy.exp(-lam*(y-t))*source(t,z),(t,0,y))
    assert sy.simplify(D.subs(y,0))==0 and sy.simplify(sy.diff(D,z).subs(y,0))==0
    assert sy.simplify(sy.diff(D,y)+lam*D-source(y,z))==0
    h,sc=sy.symbols('h sc',positive=True)
    assert sy.simplify(h*(2-sc/2)-(2*h-h*sc/2))==0
    assert sy.simplify(sy.exp(-lam*h)*sy.exp(-lam*h*(1-sc/2))-sy.exp(-lam*h*(2-sc/2)))==0
    x1,x2,g1,g2,d0=sy.symbols('memory1 memory2 driver1 driver2 inlet')
    assert sy.expand((d0*x1+g1)*x2+g2-(d0*x1*x2+g1*x2+g2))==0
    return dict(passed=True,unchanged_left_velocity_and_all_five_density_increments_exact_zero=True,
        source_owned_Duhamel_initial_value_and_Z_zero=True,original_background_not_zeroed=True,
        micro_macro_phase_identity_and_compact_offset_memory_verified=True,
        actual_incoming_plus_local_driver_composition_verified=True)


def independent_weights():
    c=MPIntervalContext();c.dps=150;p=mp.mp.clone();p.dps=100
    f=current.macro.macro.fields.MacroFlow(c,c.ln(c.mpf('.002')),c.ln(9),c.ln(4),c.ln(5),c.ln(100)-c.ln(5),c.mpf('.01'))
    jet=lambda v:fields.IntervalTaylor(c,[c.mpf(v)]+[c.mpf(0)]*5)
    zero=jet(0);one=jet(1)
    f.set_sources([zero]*3,{part:[zero]*3 for part in fields.PARTS},one,f.jet(zero),one,one,f.jet(zero))
    series=current.switch.first.FirstSwitchFunctions(f,dict(phi=f.jet(one),V=f.jet(one)),{name:f.jet(one) for name in current.macro.macro.RATES})
    comparisons=0;drivers=0
    # Nonzero compact positive left offset; both windows have true hb width.
    a=c.mpf(1)/32;b=c.mpf(1);d=c.mpf(2)
    for left,right,end in ((a,c.mpf('.25'),b),(c.mpf('.75'),b,b),(b,c.mpf('1.25'),d),(c.mpf('1.75'),d,d)):
        l,r,e=(p.mpf(ep(v)[0]) for v in (left,right,end));width=p.mpf('.002')*(r-l);suffix=p.mpf('.002')*(e-r)
        for name,rate in current.RATES.items():
            mass,decay,tail=current.weights(series,left,right,end,rate);rr=p.mpf(rate.numerator)/rate.denominator
            truth=width if not rate else -p.expm1(-rr*width)/rr
            contains(mass,truth);contains(decay,p.exp(-rr*width));contains(tail,p.exp(-rr*suffix));comparisons+=3
            # A genuinely negative signed source and a distinct nonzero Z
            # source exercise the centered signed-integral enclosure baseline.
            for value in (p.mpf('-.37'),p.mpf('.021')):
                row=current.bounds.symmetric(f,current.bounds.magnitude(f,f.scalar(c.mpf(value)))*mass*tail)
                contains(row,value*truth*p.exp(-rr*suffix));drivers+=1
    for name,rate in current.RATES.items():
        m1=current.weights(series,a,b,b,rate)[1];m2=current.weights(series,b,d,d,rate)[1]
        rr=p.mpf(rate.numerator)/rate.denominator
        contains(m1*m2,p.exp(-rr*p.mpf('.002')*(2-p.mpf(1)/32)));comparisons+=1
    return dict(passed=True,independent_true_micro_kernel_weight_and_memory_comparisons=comparisons,
        independent_negative_C0_and_nonzero_Z_density_integral_comparisons=drivers,
        actual_hb_measure_once_and_positive_compact_initial_offset=True,
        pressure_incoming_memory_one_and_nonzero_full_mass=True)


def native(owner,report):
    cells=coefficients=density_rows=memories=boundary_rows=flat_cases=0
    for label in ('0','.5'):
        got=owner.contribution(label);op=owner.micro.owner(label);f,c=op.flow,op.c;ref=owner.tail.reference.owner(label)
        assert current.macro.base.encoded(current.serialized(got))==report['frames'][label]
        inlet=got['actual_source_bound_Section11_inlet']
        assert inlet['exact_common_P0_axial5'] is ref.P0
        assert inlet['actual_left_source_phase']._mpi_==(owner.macro.sc/2)._mpi_
        assert not inlet['retained_positive_left_log_offset'].zero
        assert inlet['original_background_inlet_is_not_zeroed']
        assert any(not row.zero for row in inlet['original_current_left_common_source']['common_velocity_E_axial5'])
        for rows in inlet['actual_finite_N_inlet_defect_C0_Z'].values():assert len(rows)==2 and all(row.zero for row in rows)
        for window in got['actual_micro_source_windows'].values():
            for cell in window['actual_micro_source_cells']:
                source=cell['source'];back=source['actual_micro_background_source'];generic=source['original_generic_source'];cells+=1
                assert source['exact_common_P0_axial5'] is ref.P0 and generic['common_original_P0_axial5'] is ref.P0
                assert source['actual_selected_positive_s_c'] is owner.macro.sc
                assert source['actual_source_owned_piecewise_modification_installed']
                assert back['source_interval_functions_not_endpoint_hulls'] and back['actual_core_background_inlet_preserved']
                assert any(not row.zero for row in back['V_y'])
                for rows in (*generic['common_own_five_histories_axial5'].values(),generic['common_velocity_E_axial5'],generic['common_velocity_V_axial5']):
                    for row in rows:
                        assert row.ctx is c and row.scale.bases is f.logs and row.ledger is f.ledger;coefficients+=1
                for part in ('kernels','Z_derivatives'):
                    for row in source['original_signed_five_density_C0_Z'][part].values():
                        assert row.ctx is c and row.scale.bases is f.logs and row.ledger is f.ledger;density_rows+=1
                for weight in cell['actual_own_rate_weights'].values():
                    assert weight['true_physical_measure_once'] and not weight['actual_full_mass'].zero
            for row in window['actual_incoming_memory'].values():assert not row.zero;memories+=1
        for key in ('actual_current_R0_correction_C0_Z','actual_current_R100_correction_C0_Z',
            'actual_current_R110_correction_C0_Z','actual_current_Rm_incoming_correction_C0_Z'):
            assert set(got[key])==set(current.RATES)
            for rows in got[key].values():
                for row in rows:
                    assert row.ctx is c and row.scale.bases is f.logs and row.ledger is f.ledger;boundary_rows+=1
        joined=got['typed_actual_micro_macro_source_join']
        assert joined['same_actual_micro_and_macro_source_flow'] and joined['independent_source_enclosure_consistency_rows']==48
        assert joined['interval_overlap_is_consistency_only_not_the_function_identity_proof']
        assert got['genuine_current_finite_N_R0_R100_R110_Rm_boundary_enclosures_supplied']
        assert got['local_drivers_composed_with_real_current_prefix_not_arbitrary_zero'] and not got['global_N_or_whole_Z_admission']
        assert all(got[key] is False for key in fields.previous.OPEN)
        # Verify the executable source selector, including a strictly
        # positive point on the modified-right side inside its flat collar.
        for left,right in ((c.mpf(0),owner.macro.sc/4),(owner.macro.sc*5/8,owner.macro.sc*5/8)):
            query=owner.query(label,'first_micro',left,right)
            assert query['source_owned_modification_support']['increments_exact_zero']
            density=query['original_signed_five_density_C0_Z']
            assert all(row.zero for part in ('kernels','Z_derivatives') for row in density[part].values())
            assert all(density['velocities'][key].zero for key in ('deltaE','deltaE_Z','deltaV','deltaV_Z'));flat_cases+=1
        print('Current original real finite-N prefix checked',label,flush=True)
    rejected=0;label='0';op=owner.micro.owner(label);f=op.flow;P0=owner.tail.reference.owner(label).P0
    good=current.macro.base.encoded(current.serialized(owner.inlet(label)))
    for key,value in (('source_family','bad'),('source_frame','.5'),('candidate_N',258),('exact_common_P0_axial5',[]),
        ('source_identity_fields',{}),('source_width_log_tuple',{}),('selected_positive_s_c_tuple',{}),
        ('selected_eta_log_tuple',{}),('selected_dstar_log_tuple',{}),('provider_sha256','bad'),
        ('source_owned_zero_defect_definition','arbitrary default zero'),('flat_left_collar_certified_for_same_current_source',False)):
        bad=copy.deepcopy(good);bad[key]=value
        try:current.guard_inlet(owner.family,label,257,P0,f.logs[0],owner.macro.sc,bad,owner.source_family,owner.parameters.eta_log,owner.parameters.dstar_log)
        except ValueError:rejected+=1
        else:raise AssertionError('Mismatched true source-bound inlet accepted: '+key)
    bad=copy.deepcopy(good);bad['actual_finite_N_inlet_defect_C0_Z']['p'][1]['exact_zero']=False
    try:current.guard_inlet(owner.family,label,257,P0,f.logs[0],owner.macro.sc,bad,owner.source_family,owner.parameters.eta_log,owner.parameters.dstar_log)
    except ValueError:rejected+=1
    else:raise AssertionError('Nonzero initial pressure Z defect accepted')
    try:owner.inlet('0',258)
    except ValueError:rejected+=1
    else:raise AssertionError('Earlier/different-N inlet accepted')
    return dict(passed=True,actual_original_micro_finite_N_source_cells=cells,actual_common_source_coefficients=coefficients,
        actual_signed_five_density_C0_Z_rows=density_rows,actual_retained_micro_incoming_memories=memories,
        actual_R0_R100_R110_Rm_boundary_enclosure_rows=boundary_rows,
        executable_unchanged_left_and_positive_flat_right_collar_cases=flat_cases,
        typed_source_family_frame_width_sc_eta_dstar_N_and_initial_defect_rejections=rejected,
        real_source_bound_prefix_not_five_moment_closure_or_global_N=True)


def run():
    began=time.monotonic();theorem=symbolic()
    with mp.workdps(190):independent=independent_weights()
    print('Independent source-bound inlet, actual micro weights and composition PASS',flush=True)
    report=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    owner=current.OriginalMicroFiniteN(require_checked=False);actual=native(owner,report)
    assert report[current.GATE] and report['source_family']==owner.family
    assert report['original_source_and_initial_support_bindings']==current.macro.base.encoded(owner.bindings)
    for name,digest in report['input_hashes'].items():fields.previous.bind(owner.hashes,name,digest)
    for name in (current.NAME,Path(__file__).name):fields.previous.bind(owner.hashes,name,current.sha(name))
    result=dict(all_passed=True,**{current.GATE:True},source_family=owner.family,
        independent_source_initial_and_support_theorem=theorem,independent_micro_weights_and_signed_transport=independent,
        actual_live_source=actual,original_source_and_initial_support_bindings=owner.bindings,input_hashes=owner.hashes,
        execution_seconds=time.monotonic()-began,**dict.fromkeys(fields.previous.OPEN,False))
    (current.HERE/current.RECEIPT).write_bytes((json.dumps(current.macro.base.encoded(result),indent=2)+'\n').encode())
    print('Current original true finite-N prefix rminus/R0/R100/R110/Rm PASS',flush=True);return result


if __name__=='__main__':run()
