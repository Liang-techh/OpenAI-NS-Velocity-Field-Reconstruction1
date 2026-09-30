"""Recompute one shared-parameter candidate; no full source certificate."""
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_core_adapter import build_source_core
from lei_ren_part1_paper_axis_norm_bounds import axis_norm_bounds, _json_value
from lei_ren_part1_paper_connection_scale_gate import run as input_gates
from lei_ren_part1_paper_exit_comparison import Section923Comparison
from lei_ren_part1_paper_exit_tangents import ExitTangents
from lei_ren_part1_paper_exit_continuation import ExitContinuation
from lei_ren_part1_paper_exit_switches import ExitSwitches
from lei_ren_part1_paper_exit_field import cone_receipt


def run(precision=260):
    with mp.workdps(precision):
        params=dict(j='1e-14', Lambda='1e36', logC='5e151',
                    logPstar='14', delta='1e-200')
        bundle=build_source_core(precision=precision, degree=18, **params)
        print('shared pressure and core rebuilt', flush=True)
        # These are provisional budgets, not certified mixed C3/K bounds.
        logK=mp.mpf('1e152'); loghb=-100-100*logK
        hb=mp.exp(loghb)
        receipt=dict(candidate_name='shared-axis-bound-candidate-1',
            shared_parameters=bundle['shared_parameters'], precision=precision,
            radial_degree=18, schedule_decimal_precision=bundle['profile'].schedule.decimal_precision,
            A_provisional='1e150', logK_provisional_upper=str(logK),
            log_hb=mp.nstr(loghb, precision), hb_equals_epsilon=True,
            epsilon_positive=bool(hb>0), source_parameter_regime_certified=False,
            mixed_core_norm_certified=False, K_upper_certified=False,
            outer_connection_complete=False, scale_recursion_established=False,
            global_finite_energy_certified=False,
            axis_pressure_anchor=bundle['pressure'].anchor_receipt(),
            axis_norm_bounds=_json_value(axis_norm_bounds(bundle['axis'])),
            necessary_input_gates=input_gates(precision, h_b=hb,
                                              write_receipt=False, **params))
        receipt['samples']=[]
        for Z in (bundle['axis'].Z0, mp.mpf('.3')):
            for s in ('1','4','4.1'):
                v=bundle['core'].evaluate(mp.mpf(s)/bundle['Lambda'],Z)
                receipt['samples'].append(dict(Z=mp.nstr(Z,80),s=s,
                    F_positive=bool(v['F']>0), F_R_negative=bool(v['F_R']<0),
                    Uz=mp.nstr(v['Uz'],60),log_abs_F=mp.nstr(mp.log(abs(v['F'])),precision)))
        path=Path(__file__).with_suffix('.json')
        def save():
            path.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
        save()
        comparison=Section923Comparison(bundle,h_b=hb,transition_steps=16)
        provider=ExitContinuation(ExitTangents(comparison,epsilon=hb,steps=16,
                                             derivative_step='1e-45'))
        switches=ExitSwitches(provider,steps=16)
        receipt['connection_samples']=[]
        Z=mp.mpf('.3')
        probes=[('stage1_mid',lambda:switches.evaluate_switch_phase(1,'.5',Z)),
                ('stage2_mid',lambda:switches.evaluate_switch_phase(2,'.5',Z)),
                ('stage2_end',lambda:switches.evaluate_switch_phase(2,'1',Z)),
                ('R110',lambda:switches.evaluate_R('110',Z))]
        for name,evaluate in probes:
            print('computing '+name,flush=True)
            v=evaluate()
            row=dict(name=name,Z='.3',R=mp.nstr(v['R'],60),
                a=mp.nstr(v['a'],60), b=mp.nstr(v['b'],60),
                F_positive=bool(v['F']>0),Uz=mp.nstr(v['Uz'],60),
                Ur=mp.nstr(v['Ur'],60),cone=cone_receipt(v),
                explicit_phase_coordinate=v.get('explicit_phase_coordinate',False),
                physical_radius_offset_resolved=v.get('physical_radius_offset_resolved'))
            receipt['connection_samples'].append(row); save()
            print(json.dumps(row),flush=True)
        receipt['connection_run_complete']=True
        save()
        return receipt


if __name__=='__main__':
    run()
