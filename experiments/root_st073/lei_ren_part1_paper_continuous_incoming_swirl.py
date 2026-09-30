"""Regenerate incoming reference swirl energy from the installed schedule.

Uses the same normalized angular propagation as actual cumulative moments,
before applying actual inner offsets. No inherited float swirl row is reused.
Quadrature/input uncertainty remains unenclosed.
"""
from types import SimpleNamespace
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_joined_outer import _mp
from lei_ren_part1_paper_axial_correction import signed_log


class ContinuousIncomingSwirl:
    def __init__(self,source,*,order=96):
        from lei_ren_part1_paper_continuous_angular_moments import ContinuousAngularMoments
        correction=source.outer.angular._continuous_provider
        self.source=source;self.precision=source.precision
        proxy=SimpleNamespace(schedule=source.schedule,precision=source.precision,
            angular_correction_provider=correction,log_at=source.outer.log_at,
            offset=source.outer.offset)
        self.engine=ContinuousAngularMoments(proxy,order=order,
            mp_nodes=getattr(source,'continuous_pressure_anchor',False))

    def reference(self,Z):
        with mp.workdps(self.precision):
            z=_mp(Z);s=self.source.schedule
            y=mp.nstr(_mp(s.y_p),self.precision)
            state=self.engine._baseline(y,mp.nstr(z,self.precision))
            E=self.engine._base(y)
            scale=mp.exp(_mp(s.y_p)+2*E)/2
            return dict(I_swirl=scale*state[1],I_swirl_Z=scale*state[3],
                normalized_reference_swirl=state[1],
                inner_offset_included=False,shared_angular_propagation=True,
                quadrature_order=self.engine.order,quadrature_error_enclosed=False)


def run():
    from lei_ren_part1_paper_joined_outer import build_joined_field
    from lei_ren_part1_paper_continuous_incoming_outer import ContinuousIncomingOuterField
    print('building shared field with regenerated incoming swirl',flush=True)
    source=build_joined_field();folder=Path(__file__).parent
    seeded=json.loads((folder/'lei_ren_part1_paper_seeded_shared_candidate.json').read_text())
    atoms=json.loads((folder/'lei_ren_part1_paper_continuous_axial_solve.json').read_text())
    field=ContinuousIncomingOuterField(source,prepared={.3:(seeded,atoms)})
    p=field.outer;s=source.schedule
    with mp.workdps(field.precision):
        z=mp.mpf('.3');receipt=p.seed_solve(float(z))
        dimensions=receipt['incoming']['dimensionless_integrals']
        reference=source._continuous_incoming_swirl_reference.reference(z)
        row=p.angular_moments_jet(str(s.logR_p),z)
        Rref=mp.exp(_mp(s.logRref))
        offset=row['inner_offsets'][1];offsetZ=row['inner_offsets'][3]
        expected=2*Rref*mp.mpf(dimensions['I_swirl'])+offset
        expectedZ=2*Rref*reference['I_swirl_Z']+offsetZ
        errors=[abs(expected/row['swirl_energy']-1),abs(expectedZ/row['swirl_energy_Z']-1)]
        if max(errors)>mp.mpf('1e-70'):raise ArithmeticError('Regenerated incoming swirl disagrees with actual angular primitive')
        change=abs(mp.mpf(dimensions['I_swirl'])/mp.mpf(seeded['incoming']['dimensionless_integrals']['I_swirl'])-1)
        factor_error=abs(reference['I_swirl_Z']/(-4*z/(1+z*z)*reference['I_swirl'])-1)
        if factor_error>mp.mpf('1e-70'):raise ArithmeticError('Incoming swirl Z factorization lost')
        runtime=p.runtime(float(z));tangent=p.coefficient_tangent(float(z))
        report=dict(I_swirl=signed_log(reference['I_swirl'],80),
            inherited_swirl_relative_change=mp.nstr(change,40),
            actual_angular_Rp_relative_matching=[mp.nstr(v,40) for v in errors],
            reference_Z_factorization_relative_error=mp.nstr(factor_error,40),
            scaled_linear_rhs_inputs=receipt['axial']['linear_rhs_inputs'].get('base'),
            solve_receipt_linear_relative_replay=runtime.solve_receipt.get('linear_relative_replay'),
            live_coefficient_tangent_evaluated=bool(tangent),
            inherited_swirl_replaced=True,inner_offset_applied_once=True,
            input_error_enclosed=False,complete_five_moments=False,
            finite_energy_certified=False,scale_recursion_certified=False)
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'swirl_change':report['inherited_swirl_relative_change'],'matching':report['actual_angular_Rp_relative_matching']}),flush=True)
    return report


if __name__=='__main__':run()
