"""Pressure compatibility audit of the saved post-axial-support probe.

This audit uses rounded receipts, not a new quadrature or an error enclosure.
With Uz=Uz_y=Uz_Z=0, N_z is exactly the axial pressure contribution
sqrt(R/2) * [2(1+delta) Z P - (1-Z^2) P_Z] / (1-delta Z^2).
The separate linear/quadratic contributions are not recoverable when their
difference from this term falls below the receipt's serialization precision.
"""
import json
from pathlib import Path
import mpmath as mp


def run():
    path=Path(__file__)
    probe=json.loads(path.with_name('lei_ren_part1_paper_continuous_bundle_stress_check.json').read_text())
    rows=[]
    with mp.workdps(100):
        z=mp.mpf('.3');delta=mp.mpf('1e-200');d=1-z*z
        for sample in probe['samples']:
            if sample['label'] not in ('flatten','heat'):
                raise ValueError('Audit requires known post-axial-support samples')
            value=lambda row:mp.mpf(row['arbitrary_exponent_value'])
            P=value(sample['P']);PZ=value(sample['P_Z'])
            H=2*(1+delta)*z*P-d*PZ
            Iz=value(sample['stress']['I_z']);Nz=value(sample['stress']['N_z'])
            scale=max(abs(Iz),abs(Nz))
            rows.append(dict(label=sample['label'],P=mp.nstr(P,50),P_Z=mp.nstr(PZ,50),
                pressure_transport_coefficient=mp.nstr(H,50),
                required_P_Z_for_zero_pressure_transport=mp.nstr(2*(1+delta)*z*P/d,50),
                axial_inertia_minus_pressure_relative_at_receipt_precision=mp.nstr(abs(Iz-Nz)/scale,30),
                subdominant_moment_terms_resolved=False,
                nonzero_pressure_transport_detected=bool(H)))
        if not all(row['nonzero_pressure_transport_detected'] for row in rows):
            raise ArithmeticError('Expected materialized pressure mismatch was not reproduced')
    report=dict(samples=rows,Z='.3',delta='1e-200',
        source='continuous_bundle_stress_check.json; post-Rv Uz and jets vanish by support',
        audit_of_rounded_receipts=True,independent_full_profile_replay=False,
        pressure_gauge_reset_applied=False,pressure_terminal_compatibility_certified=False,
        finite_energy_certified=False,stress_cone_certified=False,scale_recursion_certified=False,
        next_action='Restore actual inner-seeded pressure total with angular correction; do not subtract a Z-dependent pressure datum without updating the inner construction.')
    path.with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(rows,indent=2))
    return report


if __name__=='__main__':run()
