"""Constraint-governance audit for constrained reconstruction experiment configs.

This module does not judge PDE validity. It checks that later experiment
configs preserve the preregistered CR001 contract unless they explicitly
fork it, and it surfaces metadata ambiguity such as duplicate experiment IDs.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import argparse
import json
from pathlib import Path
from typing import Any, Iterable


PROTECTED_PATHS: tuple[tuple[str, ...], ...] = (
    ("schema_version",),
    ("source_map",),
    ("units",),
    ("equation",),
    ("nu",),
    ("domain", "physical"),
    ("domain", "evaluation_box"),
    ("domain", "support"),
    ("domain", "boundary"),
    ("domain", "time_interval"),
    ("domain", "initial_condition"),
    ("domain", "pressure_gauge"),
    ("structure", "h"),
    ("structure", "tau"),
    ("structure", "radial_scale"),
    ("structure", "axial_scale"),
    ("structure", "leading_swirl_axial_scale"),
    ("structure", "symmetry"),
    ("structure", "core_probe"),
    ("structure", "core_sign_requirements"),
    ("forcing",),
    ("nontriviality",),
    ("optimization", "seed"),
    ("optimization", "sampling"),
    ("optimization", "interior_points"),
    ("optimization", "maximum_function_evaluations"),
    ("optimization", "loss_weights"),
    ("optimization", "hard_constraints"),
    ("validation", "seed"),
    ("validation", "held_out_points"),
    ("validation", "times"),
    ("validation", "operator"),
    ("validation", "derivative_steps"),
    ("validation", "quadrature_orders_per_axis"),
    ("validation", "vary_one_at_a_time"),
    ("validation", "norms"),
    ("validation", "thresholds"),
    ("validation", "failure_policy"),
    ("validation", "time_derivatives"),
    ("excluded_claims",),
)


@dataclass(frozen=True)
class ConstraintGovernanceReport:
    canonical: str
    files: tuple[str, ...]
    protected_contract_pass: bool
    protected_mismatches: tuple[str, ...]
    semantic_errors: tuple[str, ...]
    identity_collisions: dict[str, tuple[str, ...]]
    repeated_lineage_metadata: dict[str, tuple[str, ...]]
    training_validation_seeds_separate: bool
    governance_pass: bool

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _read_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"{path}: top-level JSON must be an object")
    return value


def _get(data: dict[str, Any], path: tuple[str, ...]) -> Any:
    current: Any = data
    for key in path:
        if not isinstance(current, dict) or key not in current:
            raise KeyError(".".join(path))
        current = current[key]
    return current


def _lineage_key(data: dict[str, Any]) -> str | None:
    parent = data.get("parent_experiment")
    reason = data.get("change_reason")
    if parent is None and reason is None:
        return None
    return json.dumps([parent, reason], ensure_ascii=True, sort_keys=True)


def audit_constraint_family(
    paths: Iterable[Path | str],
    *,
    canonical_path: Path | str | None = None,
) -> ConstraintGovernanceReport:
    normalized = tuple(sorted((Path(p) for p in paths), key=lambda p: p.name))
    if not normalized:
        raise ValueError("at least one constraint config is required")

    if canonical_path is None:
        canonical = next((p for p in normalized if p.name == "constraints.json"), normalized[0])
    else:
        canonical = Path(canonical_path)
    if canonical not in normalized:
        raise ValueError("canonical_path must be included in paths")

    configs = {path: _read_json(path) for path in normalized}
    base = configs[canonical]

    mismatches: list[str] = []
    errors: list[str] = []
    seed_separation = True

    for path, data in configs.items():
        for protected in PROTECTED_PATHS:
            try:
                expected = _get(base, protected)
                actual = _get(data, protected)
            except KeyError as exc:
                errors.append(f"{path.name}: missing protected field {exc.args[0]}")
                continue
            if actual != expected:
                mismatches.append(f"{path.name}: {'.'.join(protected)} differs from canonical")

        status = data.get("status")
        if not isinstance(status, str) or "validated" not in status:
            errors.append(f"{path.name}: status must explicitly encode validation state")

        nu = data.get("nu")
        if isinstance(nu, bool) or not isinstance(nu, (int, float)) or nu <= 0:
            errors.append(f"{path.name}: nu must be positive")

        try:
            t0, t1 = _get(data, ("domain", "time_interval"))
            times = _get(data, ("validation", "times"))
            training_seed = _get(data, ("optimization", "seed"))
            validation_seed = _get(data, ("validation", "seed"))
            thresholds = _get(data, ("validation", "thresholds"))
            reference_energy = _get(data, ("nontriviality", "reference_energy"))
            reference_time = _get(data, ("nontriviality", "reference_time"))
        except (KeyError, TypeError, ValueError):
            errors.append(f"{path.name}: malformed time/validation/nontriviality contract")
            continue

        if not (isinstance(t0, (int, float)) and isinstance(t1, (int, float)) and t0 < t1):
            errors.append(f"{path.name}: invalid time interval")
        elif (
            not isinstance(times, list)
            or not times
            or any(not isinstance(t, (int, float)) or t < t0 or t > t1 for t in times)
        ):
            errors.append(f"{path.name}: validation times must stay inside the preregistered window")
        elif times[0] != t0 or times[-1] != t1:
            errors.append(f"{path.name}: validation times must include both window endpoints")

        if training_seed == validation_seed:
            seed_separation = False
            errors.append(f"{path.name}: training and validation seeds must differ")

        if isinstance(reference_energy, bool) or not isinstance(reference_energy, (int, float)) or reference_energy <= 0:
            errors.append(f"{path.name}: reference energy must be positive")
        if reference_time != t0:
            errors.append(f"{path.name}: energy reference time must equal the window start")

        if not isinstance(thresholds, dict) or not thresholds:
            errors.append(f"{path.name}: validation thresholds must be a nonempty object")
        else:
            for key, value in thresholds.items():
                if isinstance(value, bool) or not isinstance(value, (int, float)) or value <= 0:
                    errors.append(f"{path.name}: threshold {key} must be positive")

        force = data.get("forcing")
        if not isinstance(force, dict):
            errors.append(f"{path.name}: forcing contract missing")
        else:
            restriction = str(force.get("restriction", "")).lower()
            if "no residual-dependent" not in restriction or "pointwise free force" not in restriction:
                errors.append(f"{path.name}: forcing restriction no longer forbids residual-defined/free forcing")

        bounds = data.get("optimization", {}).get("candidate_parameter_bounds", {})
        amplitude_rule = bounds.get("amplitude") if isinstance(bounds, dict) else None
        if not isinstance(amplitude_rule, str) or "not optimized" not in amplitude_rule or "reject outside" not in amplitude_rule:
            errors.append(f"{path.name}: amplitude-collapse guard is missing")

    ids: dict[str, list[str]] = {}
    for path, data in configs.items():
        experiment_id = data.get("experiment_id")
        if not isinstance(experiment_id, str) or not experiment_id:
            errors.append(f"{path.name}: experiment_id must be a nonempty string")
            continue
        ids.setdefault(experiment_id, []).append(path.name)
    collisions = {
        key: tuple(sorted(names))
        for key, names in ids.items()
        if len(names) > 1
    }

    lineage_groups: dict[str, list[str]] = {}
    for path, data in configs.items():
        key = _lineage_key(data)
        if key is not None:
            lineage_groups.setdefault(key, []).append(path.name)
    repeated_lineage = {
        key: tuple(sorted(names))
        for key, names in lineage_groups.items()
        if len(names) > 1
    }

    protected_pass = not mismatches
    governance_pass = protected_pass and not errors and not collisions

    return ConstraintGovernanceReport(
        canonical=canonical.name,
        files=tuple(path.name for path in normalized),
        protected_contract_pass=protected_pass,
        protected_mismatches=tuple(sorted(mismatches)),
        semantic_errors=tuple(sorted(errors)),
        identity_collisions=collisions,
        repeated_lineage_metadata=repeated_lineage,
        training_validation_seeds_separate=seed_separation,
        governance_pass=governance_pass,
    )


def _default_paths(root: Path) -> list[Path]:
    return sorted((root / "configs").glob("constraints*.json"))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Audit constrained CR001 config governance.")
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--strict-identity", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)

    report = audit_constraint_family(_default_paths(args.root))
    payload = json.dumps(report.to_dict(), indent=2, sort_keys=True)
    if args.output is None:
        print(payload)
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload + "\n", encoding="utf-8")

    if report.protected_mismatches or report.semantic_errors:
        return 2
    if args.strict_identity and report.identity_collisions:
        return 3
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
