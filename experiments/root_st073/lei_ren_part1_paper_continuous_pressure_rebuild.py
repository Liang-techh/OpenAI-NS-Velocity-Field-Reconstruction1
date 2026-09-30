"""Rebuild the nonlinear inner field with a continuous pre-Rv pressure datum.

The dominant separated axis jet omits the post-Rv tail as explicitly marked.
No completed field's pressure is shifted. This pilot measures the actual
terminal pressure of the newly rebuilt shared field, not full stress closure.
"""
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_joined_outer import build_joined_field
from lei_ren_part1_paper_continuous_incoming_outer import ContinuousIncomingOuterField
from lei_ren_part1_paper_continuous_moment_bundle import ContinuousMomentBundle
from lei_ren_part1_paper_axial_correction import signed_log


def run():
    path=Path(__file__);folder=path.parent
    print('rebuilding inner core from continuous pressure anchor',flush=True)
    source=build_joined_field(continuous_pressure=True)
    seeded=json.loads((folder/'lei_ren_part1_paper_seeded_shared_candidate.json').read_text())
    atoms=json.loads((folder/'lei_ren_part1_paper_continuous_axial_solve.json').read_text())
    field=ContinuousIncomingOuterField(source,prepared={.3:(seeded,atoms)})
    p=field.outer
    with mp.workdps(field.precision):
        z=mp.mpf('.3');p.pressure_moments_jet(p.schedule.logR_tail,z)
        engine=p.pressure_moment_provider
        if not engine.MP_nodes or engine.order!=192:
            raise ArithmeticError('Public rebuilt field did not select coherent MP pressure quadrature')
        result=engine.terminal_pressure_jet(z)
        bundle=source.provider.reshape.switches.comparison.bundle
        inner=source.inner.evaluate_x(mp.e,z)
        baseline=json.loads((folder/'lei_ren_part1_paper_continuous_pressure_target.json').read_text())
        previous=mp.mpf(baseline['values']['P_infinity']['arbitrary_exponent_value'])
        report=dict(Z='.3',continuous_axis_pressure_installed_before_core=True,
            nonlinear_core_and_inner_reconstructed=True,pressure_gauge_shift_applied=False,
            pressure_stage_MP_order=192,axis_pressure_MP_order=192,
            axis_K=signed_log(bundle['K'],field.precision),
            P_infinity=signed_log(result['P_infinity'],field.precision),
            P_infinity_Z=signed_log(result['P_infinity_Z'],field.precision),
            absolute_defect_relative_to_previous=mp.nstr(abs(result['P_infinity']/previous),50),
            pressure_transport_coefficient=signed_log(2*(1+mp.mpf('1e-200'))*z*result['P_infinity']-(1-z*z)*result['P_infinity_Z'],field.precision),
            inner_join_stress={k:signed_log(inner['stress'][k],80) for k in ('T_theta','T_z')},
            unresolved_inner_input_entries=inner['unresolved_input_entries'],
            post_Rv_axis_pressure_tail_installed=False,
            quadrature_error_enclosed=False,pressure_terminal_compatibility_certified=False,
            finite_energy_certified=False,stress_cone_certified=False,scale_recursion_certified=False)
        stress_rows=[]
        coupled=ContinuousMomentBundle(p)
        for label,radius in [('flatten',p.log_at(p.schedule.logR_v,'50')),('heat',p.log_at(p.schedule.logR_tail,'4'))]:
            print('rebuilt coupled stress '+label,flush=True)
            actual=coupled.stress(radius,z)
            stress_rows.append(dict(label=label,
                stress={k:signed_log(v,80) for k,v in actual['stress'].items()},
                stress_components={k:signed_log(v,80) for k,v in actual['stress_components'].items()},
                stress_cone_certified=False))
        report['coupled_stress_samples']=stress_rows
        legacy=json.loads((folder/'lei_ren_part1_paper_continuous_bundle_stress_check.json').read_text())
        comparisons=[]
        for current,previous in zip(stress_rows,legacy['samples']):
            if current['label']!=previous['label']:
                raise ArithmeticError('Legacy comparison points differ')
            ratios={k:mp.nstr(abs(mp.mpf(current['stress'][k]['arbitrary_exponent_value'])/
                mp.mpf(previous['stress'][k]['arbitrary_exponent_value'])),40) for k in ('T_theta','T_z')}
            comparisons.append(dict(label=current['label'],absolute_total_stress_relative_to_legacy=ratios,
                legacy_pressure_stage_nodes='binary64 96',rebuilt_pressure_stage_nodes='MP 192',
                stress_cone_certified=False))
        report['legacy_comparison']=comparisons
        if not mp.isfinite(result['P_infinity']) or not mp.isfinite(result['P_infinity_Z']):
            raise ArithmeticError('Rebuilt terminal pressure is nonfinite')
    path.with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:report[k] for k in ('P_infinity','P_infinity_Z','absolute_defect_relative_to_previous')},indent=2),flush=True)
    return report


if __name__=='__main__':run()
