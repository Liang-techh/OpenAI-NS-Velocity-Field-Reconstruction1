import numpy as np

from openai_ns_reconstruction.constrained_eq45_bipolar_phi04_capacity import (
    audit_bipolar_phi04_capacity,
)


def test_bipolar_phi04_capacity_rejects_even_mode_as_axial_reach_fix():
    report = audit_bipolar_phi04_capacity()

    assert report["truth_boundary"]["velocity_changed"] is False
    assert report["truth_boundary"]["candidate_artifact_changed"] is False
    assert report["truth_boundary"]["pde_validated"] is False
    assert report["held_out_pde_residual"]["evaluated"] is False

    assert report["base_profile_parameter_count"] == 12
    assert report["eta4_tensor_parameter_count"] == 20
    assert report["newly_active_parameters_in_proposed_minimal_extension"] == 1
    assert report["zero_padding_embedding_max_abs_velocity_error"] < 1e-12

    assert report["numerical_rank"] == 2
    assert 1.0 < report["condition_number"] < 5.0
    assert report["phi04_novelty_fraction_vs_phi12"] > 0.85
    assert report["finite_difference_refinement_relative_change"] < 1e-10

    phi12 = report["local_response"]["Phi(1,2)"]
    phi04 = report["local_response"]["Phi(0,4)"]
    assert phi04["midplane_response_rms_per_unit_coefficient"] < 1e-12
    assert phi04["axial_tip_response_rms_per_unit_coefficient"] > 15.0 * phi12[
        "axial_tip_response_rms_per_unit_coefficient"
    ]
    assert phi04["axial_tip_to_shoulder_response_ratio"] > 0.15
    assert phi12["axial_tip_to_shoulder_response_ratio"] < 0.03

    # The source-aligned bipolar child gets its desired opposite-z axial outflow
    # from an odd-eta Phi(0,1) seed.  The proposed eta-even Phi(0,4) direction
    # injects a purely z-even axial-velocity response instead of preserving it.
    assert phi04["parity"]["even_fraction_of_axial_parity_response"] > 0.999999999
    assert phi04["parity"]["axial_velocity_odd_z_response_rms"] < 1e-12

    morphology = report["whole_domain_vorticity_morphology"]
    for key in morphology["base"]:
        assert np.isfinite(morphology["base"][key])
        assert abs(morphology["plus"][key] - morphology["minus"][key]) < 1e-12

    # A bounded +/-0.25 Phi(0,4) perturbation does not move the main axial
    # 25%-core reach on the current bipolar field and barely moves smooth RMS
    # morphology, despite the local Jacobian direction being independent.
    assert abs(morphology["plus"]["core_axial_q99"] - morphology["base"]["core_axial_q99"]) < 1e-12
    axial_relative = abs(
        morphology["plus"]["weighted_axial_rms"] - morphology["base"]["weighted_axial_rms"]
    ) / morphology["base"]["weighted_axial_rms"]
    radial_relative = abs(
        morphology["plus"]["weighted_radial_rms"] - morphology["base"]["weighted_radial_rms"]
    ) / morphology["base"]["weighted_radial_rms"]
    assert axial_relative < 5e-4
    assert radial_relative < 1e-4
