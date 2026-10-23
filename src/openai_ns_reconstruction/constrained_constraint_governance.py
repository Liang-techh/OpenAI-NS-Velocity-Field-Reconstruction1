"""Machine-readable drift audit for the constrained reconstruction contract.

This module checks that basis/optimizer experiments do not silently change the
preregistered CR001 acceptance contract.  It does not evaluate Navier--Stokes
residuals and must not be used as PDE acceptance evidence.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any, Mapping


@dataclass(frozen=True)
class ConstraintDrift:
    config: str
    path: str
    baseline: Any
    observed: Any


@dataclass(frozen=True)
class ConstraintGovernanceAudit:
    acceptance_contract_pass: bool
    metadata_identity_pass: bool
    known_findings_match: bool
    seed_separation_pass: bool
    protected_drift: tuple[ConstraintDrift, ...]
    duplicate_experiment_ids: Mapping[str, tuple[str, ...]]
    config_count: int


def _load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"{path} must contain a JSON object")
    return value


def _lookup(mapping: Mapping[str, Any], dotted_path: str) -> Any:
    value: Any = mapping
    for key in dotted_path.split("."):
        if not isinstance(value, Mapping) or key not in value:
            raise KeyError(f"missing protected path {dotted_path!r}")
        value = value[key]
    return value


def _known_duplicate_map(governance: Mapping[str, Any]) -> dict[str, tuple[str, ...]]:
    expected: dict[str, tuple[str, ...]] = {}
    findings = governance.get("known_findings", [])
    if not isinstance(findings, list):
        raise ValueError("known_findings must be a list")
    for finding in findings:
        if not isinstance(finding, Mapping):
            raise ValueError("each known finding must be an object")
        if finding.get("code") != "duplicate_experiment_id":
            continue
        value = finding.get("value")
        paths = finding.get("paths")
        if not isinstance(value, str) or not isinstance(paths, list) or not all(
            isinstance(path, str) for path in paths
        ):
            raise ValueError("malformed duplicate_experiment_id finding")
        expected[value] = tuple(sorted(paths))
    return expected


def audit_constraint_directory(config_dir: str | Path) -> ConstraintGovernanceAudit:
    """Audit all ``constraints*.json`` experiment configs in ``config_dir``.

    ``constraint_governance.json`` defines protected dotted paths.  The base
    experiment is compared against every variant.  Candidate parameter bounds
    and experiment-specific nested blocks are deliberately not protected so a
    representation lane can evolve without silently changing the predeclared
    PDE/domain/forcing/nontriviality/validation contract.
    """

    config_dir = Path(config_dir)
    governance = _load_json(config_dir / "constraint_governance.json")
    baseline_name = governance.get("baseline")
    protected_paths = governance.get("protected_paths")
    seed_rule = governance.get("seed_separation")
    if not isinstance(baseline_name, str):
        raise ValueError("baseline must be a string")
    if not isinstance(protected_paths, list) or not all(
        isinstance(path, str) and path for path in protected_paths
    ):
        raise ValueError("protected_paths must be a nonempty string list")
    if not isinstance(seed_rule, Mapping):
        raise ValueError("seed_separation must be an object")

    config_paths = sorted(
        path
        for path in config_dir.glob("constraints*.json")
        if path.name != "constraint_governance.json"
    )
    if not config_paths:
        raise ValueError("no constraint configs found")
    baseline_path = config_dir / baseline_name
    if baseline_path not in config_paths:
        raise ValueError(f"baseline {baseline_name!r} not found")

    configs = {path.name: _load_json(path) for path in config_paths}
    baseline = configs[baseline_name]

    drift: list[ConstraintDrift] = []
    for name, config in configs.items():
        if name == baseline_name:
            continue
        for dotted_path in protected_paths:
            baseline_value = _lookup(baseline, dotted_path)
            observed_value = _lookup(config, dotted_path)
            if observed_value != baseline_value:
                drift.append(
                    ConstraintDrift(
                        config=name,
                        path=dotted_path,
                        baseline=baseline_value,
                        observed=observed_value,
                    )
                )

    training_path = seed_rule.get("training_path")
    validation_path = seed_rule.get("validation_path")
    must_differ = seed_rule.get("must_differ")
    if not isinstance(training_path, str) or not isinstance(validation_path, str):
        raise ValueError("seed paths must be strings")
    if must_differ is not True:
        raise ValueError("this audit only supports must_differ=true")
    seed_separation_pass = all(
        _lookup(config, training_path) != _lookup(config, validation_path)
        for config in configs.values()
    )

    id_to_paths: dict[str, list[str]] = {}
    for name, config in configs.items():
        experiment_id = config.get("experiment_id")
        if not isinstance(experiment_id, str) or not experiment_id:
            raise ValueError(f"{name} has no nonempty experiment_id")
        id_to_paths.setdefault(experiment_id, []).append(name)
    duplicates = {
        experiment_id: tuple(sorted(paths))
        for experiment_id, paths in id_to_paths.items()
        if len(paths) > 1
    }
    expected_duplicates = _known_duplicate_map(governance)

    return ConstraintGovernanceAudit(
        acceptance_contract_pass=not drift and seed_separation_pass,
        metadata_identity_pass=not duplicates,
        known_findings_match=duplicates == expected_duplicates,
        seed_separation_pass=seed_separation_pass,
        protected_drift=tuple(drift),
        duplicate_experiment_ids=duplicates,
        config_count=len(configs),
    )
