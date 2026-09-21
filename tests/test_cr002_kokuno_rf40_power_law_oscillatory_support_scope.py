from __future__ import annotations

import copy
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
AUDITOR_PATH = ROOT / "src/openai_ns_reconstruction/audit_kokuno_rf40_power_law_oscillatory_support_scope.py"
CONTRACT_PATH = ROOT / "configs/kokuno_rf40_power_law_oscillatory_support_scope.json"

spec = importlib.util.spec_from_file_location("cr002_osc_support_auditor", AUDITOR_PATH)
assert spec and spec.loader
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def _load_contract():
    return json.loads(CONTRACT_PATH.read_text(encoding="utf-8"))


def _write_temp_root(tmp_path: Path, contract: dict) -> Path:
    (tmp_path / "configs").mkdir()
    (tmp_path / "configs/kokuno_rf40_power_law_oscillatory_support_scope.json").write_text(
        json.dumps(contract, indent=2) + "\n", encoding="utf-8"
    )
    source = ROOT / "configs/constraints.json"
    (tmp_path / "configs/constraints.json").write_bytes(source.read_bytes())
    return tmp_path


def test_repository_contract_audits():
    mod.audit(ROOT)


def test_contract_has_exact_zero_power_law_and_nontrivial_inner_split():
    c = _load_contract()
    a2 = c["upstream_scope"]["agent2_1010"]
    assert a2["inner_oscillatory_nontrivial_gate"] is True
    assert a2["registered_power_law_oscillatory_abs_max_gate"] == 0.0
    locked = c["machine_locked_distinctions"]
    assert locked["inner_oscillatory_signal_required_nontrivial"] is True
    assert locked["registered_power_law_oscillatory_velocity_required_zero"] is True
    assert locked[
        "power_law_composite_label_implies_nontrivial_oscillatory_contribution_on_power_law_interval"
    ] is False


def test_zero_summand_still_allows_composite_toy_witness():
    c = _load_contract()
    toy = c["scope_logic_witness"]["toy_values"]
    total = [a + b for a, b in zip(toy["u_lead_at_probe"], toy["u_osc_at_probe"])]
    assert toy["u_osc_at_probe"] == [0.0, 0.0, 0.0]
    assert total == toy["u_total_at_probe"]
    assert any(abs(v) > 0 for v in total)


@pytest.mark.parametrize("key", mod.REQUIRED_FALSE)
def test_required_false_promotions_are_fail_closed(tmp_path, key):
    c = _load_contract()
    c["machine_locked_distinctions"][key] = True
    root = _write_temp_root(tmp_path, c)
    with pytest.raises(AssertionError):
        mod.audit(root)


@pytest.mark.parametrize("key", mod.REQUIRED_TRUE)
def test_required_true_scope_facts_cannot_be_erased(tmp_path, key):
    c = _load_contract()
    c["machine_locked_distinctions"][key] = False
    root = _write_temp_root(tmp_path, c)
    with pytest.raises(AssertionError):
        mod.audit(root)


def test_inner_nontriviality_cannot_be_relabelled_as_power_law_nontriviality(tmp_path):
    c = _load_contract()
    c["machine_locked_distinctions"][
        "power_law_composite_label_implies_nontrivial_oscillatory_contribution_on_power_law_interval"
    ] = True
    root = _write_temp_root(tmp_path, c)
    with pytest.raises(AssertionError, match="power_law_composite_label"):
        mod.audit(root)


def test_repository_zero_support_cannot_be_promoted_to_public_source_fact(tmp_path):
    c = _load_contract()
    c["machine_locked_distinctions"]["zero_power_law_oscillatory_increment_is_public_source_support_fact"] = True
    root = _write_temp_root(tmp_path, c)
    with pytest.raises(AssertionError, match="public_source"):
        mod.audit(root)


def test_a4_scoped_divergence_cannot_be_promoted_to_oscillatory_structure_validation(tmp_path):
    c = _load_contract()
    c["machine_locked_distinctions"][
        "power_law_independent_divergence_audit_implies_nontrivial_oscillatory_structure_validated_there"
    ] = True
    root = _write_temp_root(tmp_path, c)
    with pytest.raises(AssertionError, match="oscillatory_structure"):
        mod.audit(root)


def test_cr001_threshold_mutation_is_rejected(tmp_path):
    c = _load_contract()
    c["cr001_snapshot"]["momentum_max"] = 0.01
    root = _write_temp_root(tmp_path, c)
    with pytest.raises(AssertionError, match="CR001 snapshot"):
        mod.audit(root)


def test_free_force_laundering_is_rejected(tmp_path):
    c = _load_contract()
    c["cr001_snapshot"]["residual_defined_free_forcing_forbidden"] = False
    root = _write_temp_root(tmp_path, c)
    with pytest.raises(AssertionError, match="CR001 snapshot"):
        mod.audit(root)


def test_eq45_delivery_cannot_be_downgraded_by_kokuno_scope(tmp_path):
    c = _load_contract()
    c["canonical_delivery_independence"]["velocity_export_ready"] = False
    root = _write_temp_root(tmp_path, c)
    with pytest.raises(AssertionError, match="canonical Eq45"):
        mod.audit(root)


def test_scientific_change_laundering_is_rejected(tmp_path):
    c = _load_contract()
    c["scientific_changes"]["velocity_or_profile_bytes_changed"] = True
    root = _write_temp_root(tmp_path, c)
    with pytest.raises(AssertionError, match="scientific"):
        mod.audit(root)
