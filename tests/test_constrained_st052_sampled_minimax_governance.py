from copy import deepcopy

import pytest

from openai_ns_reconstruction.constrained_st052_sampled_minimax_governance import (
    audit_st052_sampled_minimax_governance,
    load_constraints,
    load_contract,
    load_project_status,
)


def test_st052_sampled_minimax_governance_passes_frozen_contract():
    receipt = audit_st052_sampled_minimax_governance()
    assert receipt["audit"] == "passed"
    assert receipt["sampled_momentum_max_improved_on_two_fresh_seeds"] is True
    assert receipt["sampled_momentum_L2_improved_on_two_fresh_seeds"] is False
    assert receipt["registered_cr001_seed_replayed"] is False
    assert receipt["pde_validated"] is False
    assert receipt["st052_is_canonical_candidate"] is False
    assert receipt["canonical_velocity_export_ready"] is True


def test_rejects_max_only_improvement_being_promoted_to_pde_or_unqualified_improvement():
    contract = load_contract()
    contract["observed_fresh_validation"]["pde_validated"] = True
    with pytest.raises(ValueError, match="ST052 PDE promotion"):
        audit_st052_sampled_minimax_governance(contract=contract)

    contract = load_contract()
    contract["observed_fresh_validation"]["unqualified_candidate_improved_claim_allowed"] = True
    with pytest.raises(ValueError, match="unqualified improvement claim"):
        audit_st052_sampled_minimax_governance(contract=contract)


def test_rejects_fresh_seed_validation_replacing_preregistered_cr001_sample():
    contract = load_contract()
    contract["observed_fresh_validation"]["fresh_seed_results_are_the_preregistered_CR001_acceptance_sample"] = True
    with pytest.raises(ValueError, match="fresh-seed/registered-sample identity"):
        audit_st052_sampled_minimax_governance(contract=contract)

    contract = load_contract()
    contract["observed_fresh_validation"]["registered_CR001_seed_replayed_by_pr508"] = True
    with pytest.raises(ValueError, match="registered seed replay claim"):
        audit_st052_sampled_minimax_governance(contract=contract)


def test_rejects_feasible_checkpoint_being_laundered_as_optimizer_convergence():
    contract = load_contract()
    contract["optimization_and_selection"]["optimizer_stages_converged"] = True
    with pytest.raises(ValueError, match="optimizer convergence claim"):
        audit_st052_sampled_minimax_governance(contract=contract)

    contract = load_contract()
    contract["optimization_and_selection"]["restored_checkpoint_is_local_or_global_optimum"] = True
    with pytest.raises(ValueError, match="optimizer optimality laundering"):
        audit_st052_sampled_minimax_governance(contract=contract)


def test_rejects_post_holdout_retuning_and_validation_recycling():
    contract = load_contract()
    contract["optimization_and_selection"]["post_holdout_retuning_under_same_experiment_identity_allowed"] = True
    with pytest.raises(ValueError, match="post-holdout retuning"):
        audit_st052_sampled_minimax_governance(contract=contract)

    contract = load_contract()
    contract["optimization_and_selection"][
        "fresh_validation_samples_may_be_recycled_for_model_selection_and_still_called_held_out"
    ] = True
    with pytest.raises(ValueError, match="held-out sample recycling"):
        audit_st052_sampled_minimax_governance(contract=contract)


def test_rejects_cr001_registered_time_or_threshold_drift():
    constraints = load_constraints()
    constraints["validation"]["times"] = [0.25, 0.35, 0.45, 0.55, 0.65, 0.75]
    with pytest.raises(ValueError, match="validation times"):
        audit_st052_sampled_minimax_governance(constraints=constraints)

    constraints = load_constraints()
    constraints["validation"]["thresholds"]["pde_residual_L2"] = 0.01
    with pytest.raises(ValueError, match="momentum L2 threshold"):
        audit_st052_sampled_minimax_governance(constraints=constraints)


def test_rejects_free_force_or_collapsed_velocity_shortcuts():
    constraints = load_constraints()
    constraints["forcing"]["restriction"] = "fit residual-defined pointwise force"
    with pytest.raises(ValueError, match="free-force guard"):
        audit_st052_sampled_minimax_governance(constraints=constraints)

    constraints = load_constraints()
    constraints["nontriviality"]["enforcement"] = "allow amplitude collapse"
    with pytest.raises(ValueError, match="nontriviality guard"):
        audit_st052_sampled_minimax_governance(constraints=constraints)


def test_rejects_st052_silent_canonical_replacement_or_truth_state_promotion():
    contract = load_contract()
    contract["canonical_delivery"]["st052_replaces_canonical_candidate"] = True
    with pytest.raises(ValueError, match="ST052 canonical replacement"):
        audit_st052_sampled_minimax_governance(contract=contract)

    status = deepcopy(load_project_status())
    status["states"]["pde_validated"] = True
    with pytest.raises(ValueError, match="project-status pde_validated"):
        audit_st052_sampled_minimax_governance(project_status=status)


def test_rejects_sampled_training_or_edge_max_as_continuum_certificate():
    contract = load_contract()
    contract["optimization_and_selection"]["sampled_epigraph_is_continuum_bound"] = True
    with pytest.raises(ValueError, match="sampled epigraph continuum promotion"):
        audit_st052_sampled_minimax_governance(contract=contract)

    contract = load_contract()
    contract["sampled_evidence_scope"]["independent_random_or_edge_grid_max_is_continuum_upper_bound"] = True
    with pytest.raises(ValueError, match="sampled evidence promotion"):
        audit_st052_sampled_minimax_governance(contract=contract)
