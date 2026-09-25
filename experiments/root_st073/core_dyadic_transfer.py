"""Measure the ST073-V inner kernel's actual dyadic profile transfer.

This tests a necessary ingredient for scale recursion: whether velocity
profiles at fixed similarity coordinates remain controlled as tau halves.
It is not an inner/outer construction or a Navier--Stokes acceptance test.
"""

import json

import numpy as np

from radial_continuation import FullRadialField, ROOT


def sample(field, k, X, eta):
    tau = 0.5 * 2.0 ** (-k)
    xx, ee = np.meshgrid(X, eta, indexing="ij")
    x, e = xx.ravel(), ee.ravel()
    points = field.from_similarity(x, e, tau)
    data = field.evaluate(points, tau)
    q = tau / (1.0 - e * e)
    normalized = data["velocity"].copy() / np.sqrt(field.nu)
    normalized[:, 0] *= q ** 0.5
    normalized[:, 1:] *= (q ** field.A)[:, None]
    normalized_residual = data["residual"].copy() / np.sqrt(field.nu)
    normalized_residual[:, 0] *= q ** 1.5
    normalized_residual[:, 1:] *= (q ** (field.A + 1.0))[:, None]
    return dict(
        k=k,
        tau=tau,
        normalized_velocity=normalized,
        normalized_residual=normalized_residual,
        physical_residual_max=float(np.max(np.linalg.norm(data["residual"], axis=1))),
        normalized_residual_max=float(np.max(np.linalg.norm(normalized_residual, axis=1))),
        divergence_max=float(np.max(np.abs(data["divergence"]))),
        peak_azimuthal=float(np.max(np.abs(data["velocity"][:, 1]))),
        min_radial_speed=float(np.min(data["velocity"][:, 0])),
    )


def run():
    path = ROOT / "NS_ST073_Full_Local_Recurrence/data/ST073-V.json"
    field = FullRadialField.load(path)
    X = np.linspace(0.001, 1.0 / 64.0, 12)
    eta = np.array([-0.3, -0.15, 0.0, 0.15, 0.3])
    samples = [sample(field, k, X, eta) for k in range(7)]
    transfers = []
    for before, after in zip(samples[:-1], samples[1:]):
        old = before["normalized_velocity"]
        new = after["normalized_velocity"]
        difference = new - old
        transfers.append(dict(
            from_k=before["k"], to_k=after["k"],
            relative_profile_rms=float(np.linalg.norm(difference) / np.linalg.norm(old)),
            max_abs_profile_change=float(np.max(np.abs(difference))),
            peak_azimuthal_growth=float(after["peak_azimuthal"] / before["peak_azimuthal"]),
            normalized_residual_max_before=before["normalized_residual_max"],
            normalized_residual_max_after=after["normalized_residual_max"],
        ))
    report = dict(
        source=str(path.relative_to(ROOT)),
        source_model="ST073-V local full-momentum radial recurrence",
        X_nodes=X.tolist(), eta_nodes=eta.tolist(),
        k_definition="tau=0.5*2**(-k)",
        dyadic_scaling=dict(
            h=float(field.h), radial_length_ratio=float(2.0 ** -0.5),
            axial_length_ratio=float(2.0 ** -field.D),
            axial_to_radial_aspect_growth=float(2.0 ** field.h),
            angular_axial_speed_growth=float(2.0 ** field.A),
            six_step_aspect_growth=float(2.0 ** (6.0 * field.h)),
        ),
        velocity_normalization="radial q**(1/2), angular/axial q**A, divided by sqrt(nu); q=tau/(1-eta**2)",
        residual_normalization="radial q**(3/2), angular/axial q**(A+1), divided by sqrt(nu)",
        scales=[{key: value for key, value in row.items()
                 if key not in ("normalized_velocity", "normalized_residual")}
                for row in samples],
        transfers=transfers,
        scope="Finite inner-core similarity grid X<=1/64, |eta|<=.3, registered k=0..6; no outer matching, recursive fixed point, global energy, or complete-field PDE acceptance.",
        scale_recursion_established=False,
        pde_validated=False,
    )
    output = ROOT / "core_dyadic_transfer.json"
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(dict(output=str(output), transfers=transfers), indent=2))
    return report


if __name__ == "__main__":
    run()
