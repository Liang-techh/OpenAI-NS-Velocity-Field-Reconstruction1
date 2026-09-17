"""Axisymmetric volume quadrature; diagnostic only, never candidate promotion."""
import json
from dataclasses import replace
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
from .eq45_supported_delivery import Eq45SupportedDeliveryField


def run(output='artifacts/bipolar_energy/report.json'):
    cfg = json.loads(Path('configs/constraints.json').read_text())
    reports = {}
    for name in ('eq45_supported', 'eq45_bipolar'):
        field = Eq45SupportedDeliveryField.load_candidate(
            f'artifacts/delivery/{name}/candidate.json')
        rows = []
        previous = {}
        for order in (24, 48, 96, 192):
            nodes, weights = leggauss(order)
            r, z = np.meshgrid(nodes + 1, 2 * nodes, indexing='ij')
            volume = 2 * np.pi * r * weights[:, None] * (2 * weights[None, :])
            for t in cfg['validation']['times']:
                velocity = field.velocity(r, 0*r, z, t)
                energy = float(.5 * np.sum(volume * np.sum(velocity**2, axis=-1)))
                rows.append(dict(order=order, time=t, energy=energy,
                    relative_change=None if t not in previous else
                    abs(energy-previous[t])/max(abs(energy), 1e-300)))
                previous[t] = energy
        reports[name] = dict(candidate_sha256=field.sha256, rows=rows,
            reference_energy_error=abs(previous[.25]-cfg['nontriviality']['reference_energy']),
            energy_range_satisfied=all(.1 <= e <= 10 for e in previous.values()),
            reference_normalization_satisfied=abs(previous[.25]-1) <= .001,
            finest_quadrature_change_satisfied=all(
                row['relative_change'] <= .001 for row in rows if row['order']==192))
    # Uniform scaling preserves signs and parity. It does NOT preserve momentum
    # balance: convection is quadratic whereas temporal/viscous terms are linear.
    source = field.candidate
    source_energy = previous[.25]
    scale = float(1 / np.sqrt(source_energy))
    basis = source.parent.profile_basis
    normalized = replace(source, parent=replace(source.parent, profile_basis=replace(
        basis, phi_coefficients=tuple(scale*c for c in basis.phi_coefficients),
        swirl_coefficients=tuple(scale*c for c in basis.swirl_coefficients))))
    candidate_path = Path(output).parent / 'normalized_candidate.json'
    candidate_path.parent.mkdir(parents=True, exist_ok=True)
    normalized.save_json(candidate_path)
    normalized_rows = []
    for t in cfg['validation']['times']:
        values = normalized.velocity_xyz(r, 0*r, z, t)
        energy = float(.5*np.sum(volume*np.sum(values**2, axis=-1)))
        normalized_rows.append(dict(time=t, energy=energy))
    times = np.linspace(.25, .75, 21)
    tau = 1-times
    probe_r, probe_z = .1*np.sqrt(tau), .1*tau**.495
    core = normalized.velocity_xyz(probe_r, 0*times, probe_z, times)
    scaled = core*np.column_stack((np.sqrt(tau), tau**.505, tau**.505))
    core_drift = float(np.max(np.linalg.norm(scaled-scaled[0], axis=1))/np.linalg.norm(scaled[0]))
    report = dict(normalized_bipolar=dict(candidate_sha256=normalized.sha256,
        source_sha256=source.sha256, scale=scale, rows=normalized_rows,
        core_probe_count=21, scaled_core_profile_relative_drift=core_drift,
        core_drift_gate_satisfied=core_drift <= .05,
        core_signs_satisfied=bool(np.all(core[:,0]<0) and np.all(core[:,1:]>0)),
        candidate_path=candidate_path.as_posix(), pde_validated=False,
        scope='Explicit new experimental candidate; momentum diagnostics for its unnormalized parent do not transfer.'),
        method='Gauss-Legendre on r in [0,2], z in [-2,2]; volume weight 2*pi*r; exact axisymmetry assumed from candidate construction',
        scope='Finite-window energy diagnostic with separately serialized energy-normalized child. No threshold changes, PDE acceptance or source identification.',
        candidates=reports)
    path=Path(output); path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: {a:b for a,b in v.items() if a!='rows'} for k,v in reports.items()},indent=2))


if __name__ == '__main__':
    run()
