"""Extract the actual pressure-restoration target including the infinite heat tail."""
import json
import sys
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_joined_outer import build_joined_field
from lei_ren_part1_paper_continuous_incoming_outer import ContinuousIncomingOuterField
from lei_ren_part1_paper_axial_correction import signed_log


def run(*,mp_node_audit=False,mp_orders=()):
    path=Path(__file__);folder=path.parent
    print('building shared candidate for actual terminal pressure target',flush=True)
    seeded=json.loads((folder/'lei_ren_part1_paper_seeded_shared_candidate.json').read_text())
    atoms=json.loads((folder/'lei_ren_part1_paper_continuous_axial_solve.json').read_text())
    field=ContinuousIncomingOuterField(build_joined_field(),prepared={.3:(seeded,atoms)})
    p=field.outer
    with mp.workdps(field.precision):
        z=mp.mpf('.3');p.pressure_moments_jet(p.schedule.logR_tail,z)
        row=p.pressure_moment_provider.terminal_pressure_jet(z)
        keys=('P_infinity','P_infinity_Z','pressure_baseline','pressure_baseline_Z',
              'current_bump','current_bump_Z','required_bump','required_bump_Z',
              'additional_bump_required','additional_bump_required_Z')
        report=dict(Z='.3',values={k:signed_log(row[k],field.precision) for k in keys},
            complete_heat_tail_included=True,components_retained_separately=True,
            pressure_datum_changed=False,angular_coefficients_changed=False,
            coefficient_feasibility_checked=False,quadrature_error_enclosed=False,
            pressure_terminal_compatibility_certified=False,stress_cone_certified=False,
            finite_energy_certified=False,scale_recursion_certified=False)
        if mp_node_audit:
            provider=p.pressure_moment_provider
            old_nodes=provider.nodes
            nodes,weights=mp.gauss_quadrature(provider.order,'legendre')
            try:
                provider.nodes=[((x+1)/2,w/2) for x,w in zip(nodes,weights)]
                audited=provider.terminal_pressure_jet(z)
            finally:
                provider.nodes=old_nodes
            report['MP_pressure_stage_node_audit']=dict(order=provider.order,
                scope='pressure variable-stage nodes only; inner anchor, angular schedule, bump and heat inputs unchanged',
                values={k:signed_log(audited[k],field.precision) for k in keys},
                pressure_shift=signed_log(audited['P_infinity']-row['P_infinity'],field.precision),
                pressure_Z_shift=signed_log(audited['P_infinity_Z']-row['P_infinity_Z'],field.precision),
                quadrature_error_enclosed=False,inner_datum_precision_audited=False)
        if mp_orders:
            provider=p.pressure_moment_provider
            old_nodes=provider.nodes
            comparisons=[]
            try:
                for order in mp_orders:
                    print('pressure-only MP stage quadrature order '+str(order),flush=True)
                    nodes,weights=mp.gauss_quadrature(order,'legendre')
                    provider.nodes=[((x+1)/2,w/2) for x,w in zip(nodes,weights)]
                    evaluated=provider.terminal_pressure_jet(z)
                    comparisons.append(dict(order=order,
                        P_infinity=signed_log(evaluated['P_infinity'],field.precision),
                        P_infinity_Z=signed_log(evaluated['P_infinity_Z'],field.precision)))
            finally:
                provider.nodes=old_nodes
            changes=[]
            for previous,current in zip(comparisons,comparisons[1:]):
                changes.append(dict(orders=[previous['order'],current['order']],
                    P_absolute_change=mp.nstr(abs(mp.mpf(previous['P_infinity']['arbitrary_exponent_value'])-mp.mpf(current['P_infinity']['arbitrary_exponent_value'])),50),
                    PZ_absolute_change=mp.nstr(abs(mp.mpf(previous['P_infinity_Z']['arbitrary_exponent_value'])-mp.mpf(current['P_infinity_Z']['arbitrary_exponent_value'])),50)))
            report['MP_pressure_stage_order_comparison']=dict(samples=comparisons,
                successive_absolute_changes=changes,
                scope='only pressure variable-stage quadrature order/nodes change; one shared field',
                quadrature_error_enclosed=False,inner_datum_precision_audited=False)
        # A correction target must be the negative actual terminal mismatch.
        for total,target in [('P_infinity','additional_bump_required'),('P_infinity_Z','additional_bump_required_Z')]:
            scale=max(abs(row[total]),abs(row[target]))
            if scale and abs(row[total]+row[target])/scale>mp.mpf('1e-100'):
                raise ArithmeticError('Terminal matching target has inconsistent sign or anchoring')
    path.with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:report['values'][k] for k in ('P_infinity','P_infinity_Z')},indent=2),flush=True)
    return report


if __name__=='__main__':run(mp_node_audit='--mp-node-audit' in sys.argv,
    mp_orders=(96,128,192) if '--mp-order-audit' in sys.argv else ())
