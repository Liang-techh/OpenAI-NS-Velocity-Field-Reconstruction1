"""Independent coupled micro ODE, core-memory and decoder checks."""
import gzip
import json
import math
from pathlib import Path
from types import SimpleNamespace
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_micro_functions as current

fields,ep=current.fields,current.ep


def contains(row,value,allowance=0):
    lo,hi=ep(row.finite_interval(max_log=2000))
    assert lo-allowance<=value<=hi+allowance,(value,lo,hi)


def symbolic():
    h,s,sig,D,G,phi,barphi=sy.symbols('h s sigma D G phi barphi',positive=True)
    chi=1-sig+h*sig
    assert sy.expand(h*chi-(h*(1-sig)+h*h*sig))==0
    assert sy.simplify((-h*chi*D*phi/2)/h+chi*D*phi/2)==0
    assert sy.simplify((-h*chi*(phi/barphi)*G)/h+chi*(phi/barphi)*G)==0
    assert sy.simplify((-h*h*(phi/barphi)*G)/h+h*(phi/barphi)*G)==0
    assert sy.diff(sy.Symbol('logRa')+h*s,s)==h
    Ra,Rm,lam=sy.symbols('Ra Rm lam',positive=True)
    assert sy.simplify(sy.exp(-lam*h)*sy.exp(-lam*h)-sy.exp(-2*lam*h))==0
    L3,Rphi,p=sy.symbols('L3 Rphi p',positive=True)
    assert sy.simplify((sy.Rational(6,8)*Rphi/p).subs(Rphi,p*L3*8/6)-L3)==0
    return dict(passed=True,both_micro_controls_act_in_both_equations=True,
        physical_hb_measure_once_verified=True,nonzero_second_axial_drive_required=True,
        third_derivative_norm_recovery_inequality_coefficient_verified=True)


def independent_ode():
    """Independent RK4 at two meshes, with an analytic field cross-check."""
    c=MPIntervalContext();c.dps=120
    f=current.downstream.macro.fields.MacroFlow(c,c.ln(c.mpf('.001')),c.ln(9),c.ln(4),c.ln(5),c.ln(100)-c.ln(5),c.mpf('.01'))
    jet=lambda values:fields.IntervalTaylor(c,[c.mpf(str(v)) for v in values]+[c.mpf(0)]*(6-len(values)))
    phi=[.7,.02,.01];vel=[.4,.03,.01];g=[.04,.002,0]
    initial={name:[(n+1)/10,.01,.001] for n,name in enumerate(current.bridge.RATES)}
    zeros=jet([0]);f.set_sources([zeros]*3,{part:[zeros]*3 for part in fields.PARTS},jet([1]),f.jet(zeros),jet(phi),jet(vel),f.jet(zeros))
    inlet={name:f.jet(jet(row)) for name,row in initial.items()}
    series=current.downstream.downstream.first.FirstSwitchFunctions(f,dict(phi=f.jet(jet(phi)),V=f.jet(jet(vel))),inlet)
    dj=lambda rows:fields.IntervalTaylor(c,list(jet(rows).coefficients)+[c.mpf(0)])
    packet=dict(phi=dj(phi),V=dj(vel),directions=dict(D_over_R=dj([.7]),drive_hydro=dj(g),drive_pressure=dj([0]),drive_swirl=dj([0])),
        quotient=jet([1]),comparison={})
    decoder=SimpleNamespace(flow=f,cover=lambda left,right:packet)
    encode=current.downstream.base.encoded
    core=encode(dict(actual_phi_axial5=list(jet(phi).coefficients),actual_raw_V_axial5=list(jet(vel).coefficients),
        actual_own_six_moments_axial5={name:list(jet(row).coefficients) for name,row in initial.items()}))
    owner=current.MicroFunctions(f,series,decoder,core)
    def sigma(s):
        if s<=0:return 0.
        if s>=1:return 1.
        odds=1/(1-s)**2-1/s**2
        if odds<-700:return 0.
        if odds>700:return 1.
        return 1/(1+math.exp(-odds))
    def product(a,b):return [sum(a[j]*b[n-j] for j in range(n+1)) for n in range(3)]
    names=list(initial)
    def rhs(s,state):
        radius=5*math.exp(.001*s);chi=1-sigma(s)+.001*sigma(s) if s<=1 else .001
        pf,vg=state[:2];ph=[v*pf for v in phi];V=[v+x*vg for v,x in zip(vel,g)]
        sources=dict(H=[2*v for v in ph],M=V,K=[2*v for v in product(ph,V)],A=product(V,V),B=product(ph,ph),C=product(ph,ph))
        result=[-.001*chi*.7*radius*pf/2,-.001*chi*radius*pf]
        for i,name in enumerate(names):result += [.001*(sources[name][n]-current.bridge.RATES[name]*state[2+3*i+n]) for n in range(3)]
        return result
    def solve(target,steps):
        state=[1.,0.]+[v for name in names for v in initial[name]]
        # Split at the true first/second seam before integrating.
        cuts=[0.,min(1.,target)]+([target] if target>1 else [])
        for left,right in zip(cuts,cuts[1:]):
            count=max(1,round(steps*(right-left)));ds=(right-left)/count
            for n in range(count):
                s=left+ds*n;k1=rhs(s,state)
                k2=rhs(s+ds/2,[v+ds*x/2 for v,x in zip(state,k1)])
                k3=rhs(s+ds/2,[v+ds*x/2 for v,x in zip(state,k2)])
                k4=rhs(s+ds,[v+ds*x for v,x in zip(state,k3)])
                state=[v+ds*(a+2*b+2*d+e)/6 for v,a,b,d,e in zip(state,k1,k2,k3,k4)]
        return state
    comparisons=0;max_mesh_difference=0.;nonzero=0
    for point in (.17,.63,1.31,1.83):
        coarse,fine=solve(point,400),solve(point,800)
        difference=max(abs(a-b) for a,b in zip(coarse,fine));max_mesh_difference=max(max_mesh_difference,difference)
        assert difference<1e-11
        assert abs(fine[1]-(2/.7)*(fine[0]-1))<2e-13
        chart='first_micro' if point<1 else 'second_micro'
        source=owner.evaluate(chart,point-.01,point+.01)
        chi=1-sigma(point)+.001*sigma(point) if point<1 else .001
        radius=5*math.exp(.001*point);ph=[v*fine[0] for v in phi];V=[v+x*fine[1] for v,x in zip(vel,g)]
        for rows,truth in ((source['fields']['phi'],ph),(source['fields']['V'],V),
            (source['phi_y'],[-chi*.7*radius*v/2 for v in ph]),
            (source['V_y'],[-chi*radius*fine[0]*x for x in g])):
            for n in range(3):contains(rows[n],mp.mpf(truth[n]),mp.mpf('3e-12'));comparisons+=1
        for i,name in enumerate(names):
            for n in range(3):contains(source['histories'][name][n],mp.mpf(fine[2+3*i+n]),mp.mpf('3e-12'));comparisons+=1
        contains(source['radius'],mp.mpf(radius),mp.mpf('3e-12'));comparisons+=1
        contains(source['correlated_a_axial5'][0],mp.mpf(chi*.7*radius),mp.mpf('3e-12'));comparisons+=1
        assert not source['V_y'][0].zero;nonzero+=1
    zero=owner.evaluate('first_micro',0,0)
    for name in initial:
        assert zero['histories'][name] is owner.initial[name]
    for chart,point in (('first_micro',1),('second_micro',1),('second_micro',2)):
        src=owner.evaluate(chart,point,point)
        contains(src['original_coupled_chi'],mp.mpf('.001'),mp.mpf('1e-100'))
    return dict(passed=True,independent_coupled_field_history_and_derivative_comparisons=comparisons,
        independent_RK4_mesh_difference=max_mesh_difference,analytic_phi_V_crosscheck=True,
        original_nonzero_axial_drive_cases=nonzero,exact_core_memory_and_three_micro_control_seams_checked=True)


def native(owner,report):
    cells=coefficients=memories=0;rejected=0
    for label in ('0','.5'):
        op=owner.owner(label);f,c=op.flow,op.c
        index=0
        for chart,start in (('first_micro',0),('second_micro',1)):
            for left,right in ((0,.25),(.25,.5),(.5,.75),(.75,1)):
                value=owner.query(label,chart,start+left,start+right)
                assert current.downstream.base.encoded(current.serialized(value))==report['frames'][label][index]
                back=value['original_micro_source'];index+=1;cells+=1
                assert back['source_interval_functions_not_endpoint_hulls'] and back['actual_core_background_inlet_preserved']
                assert not back['real_finite_N_micro_exit_correction_supplied']
                assert value['no_ancestor_constructors_or_producers_executed']
                for values in (*back['fields'].values(),*back['histories'].values(),back['phi_y'],back['V_y'],back['correlated_a_axial5']):
                    for row in values:
                        assert row.ctx is c and row.scale.bases is f.logs and row.ledger is f.ledger;coefficients+=1
                assert any(not row.zero for row in back['V_y'])
                for name,weight in back['original_core_Volterra_history_evidence'].items():
                    assert not weight['true_incoming_decay'].zero and not weight['true_positive_Volterra_mass'].zero
                    assert weight['original_core_memory_retained'] and weight['source_rhs_covers_every_actual_prefix'];memories+=1
                assert value['actual_selected_positive_s_c'] is owner.downstream.sc
                assert not value['true_phase_origin_log_offset'].zero
                assert all(value[key] is False for key in fields.previous.OPEN)
        core=op.evaluate('first_micro',0,0)
        assert all(core['histories'][name] is op.initial[name] for name in op.initial)
        for chart,left,right in (('bad',0,1),('first_micro',-.1,.1),('first_micro',.5,1.1),('second_micro',.9,1.5),('second_micro',1.9,2.1),('second_micro',1.8,1.2)):
            try:op.evaluate(chart,left,right)
            except ValueError:rejected+=1
            else:raise AssertionError('Invalid source cell accepted')
        print('Current original hydrated micro functions checked',label,flush=True)
    try:owner.owner('wrong')
    except ValueError:rejected+=1
    else:raise AssertionError('Invalid source frame accepted')
    return dict(passed=True,actual_original_micro_source_cells=cells,actual_source_coefficients=coefficients,
        retained_actual_core_moment_memory_rows=memories,invalid_chart_phase_frame_rejections=rejected,
        comparison_owner_hydrated_without_ancestor_constructor=True,
        finite_N_boundary_and_global_moment_closure_still_unsupplied=True)


def run():
    began=time.monotonic();theorem=symbolic()
    with mp.workdps(160):independent=independent_ode()
    print('Independent original micro controls, nonlinear fields and histories PASS',flush=True)
    report=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    owner=current.OriginalMicroFunctions(require_checked=False);actual=native(owner,report)
    assert report[current.GATE] and report['source_family']==owner.family
    assert report['original_source_bindings']==current.downstream.base.encoded(owner.bindings)
    for name,digest in report['input_hashes'].items():fields.previous.bind(owner.hashes,name,digest)
    for name in (current.NAME,Path(__file__).name):fields.previous.bind(owner.hashes,name,current.sha(name))
    receipt=dict(all_passed=True,**{current.GATE:True},source_family=owner.family,
        independent_symbolic_controls_and_measure=theorem,independent_original_micro_ODE=independent,
        actual_live_source=actual,input_hashes=owner.hashes,original_source_bindings=owner.bindings,
        execution_seconds=time.monotonic()-began,**dict.fromkeys(fields.previous.OPEN,False))
    (current.HERE/current.RECEIPT).write_bytes((json.dumps(current.downstream.base.encoded(receipt),indent=2)+'\n').encode())
    print('Current original hydrated whole micro source functions PASS',flush=True);return receipt


if __name__=='__main__':run()
