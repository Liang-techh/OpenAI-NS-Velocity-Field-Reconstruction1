"""Fail-closed CR002 audit for strict-inner vorticity versus velocity identity.

This audit does not alter any candidate.  It binds the exact A5 #884 registration
surface and canonical CR001 constraints, then machine-locks the narrower
representation rule that a local curl/vorticity artifact is not a substitute
for the velocity-field identity required by the project deliverable.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Callable, Iterable, Sequence


TASK = "CR002-KOKUNO-VORTICITY-VELOCITY-IDENTITY-081"
CONFIG_REL = Path("configs/kokuno_vorticity_velocity_identity_scope.json")
UPSTREAM_REL = Path(
    "src/openai_ns_reconstruction/"
    "kokuno_a5_strict_inner_vorticity_artifact_ingest_contract.py"
)
CONSTRAINTS_REL = Path("configs/constraints.json")


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _git_blob_sha1(data: bytes) -> str:
    header = f"blob {len(data)}\0".encode("ascii")
    return hashlib.sha1(header + data).hexdigest()


def _load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _field_a(x: float, y: float, z: float) -> tuple[float, float, float]:
    del z
    return (-0.5 * y, 0.5 * x, 0.0)


def _field_b(x: float, y: float, z: float) -> tuple[float, float, float]:
    ux, uy, _ = _field_a(x, y, z)
    return (ux, uy, 0.3)


def _centered_jacobian(
    field: Callable[[float, float, float], Sequence[float]],
    point: Sequence[float],
    h: float = 1.0e-6,
) -> tuple[tuple[float, float, float], ...]:
    base = [float(v) for v in point]
    rows = [[0.0, 0.0, 0.0] for _ in range(3)]
    for axis in range(3):
        plus = list(base)
        minus = list(base)
        plus[axis] += h
        minus[axis] -= h
        fp = field(*plus)
        fm = field(*minus)
        for component in range(3):
            rows[component][axis] = (float(fp[component]) - float(fm[component])) / (2.0 * h)
    return tuple(tuple(row) for row in rows)


def _divergence(jacobian: Sequence[Sequence[float]]) -> float:
    return float(jacobian[0][0] + jacobian[1][1] + jacobian[2][2])


def _curl(jacobian: Sequence[Sequence[float]]) -> tuple[float, float, float]:
    return (
        float(jacobian[2][1] - jacobian[1][2]),
        float(jacobian[0][2] - jacobian[2][0]),
        float(jacobian[1][0] - jacobian[0][1]),
    )


def _max_abs(values: Iterable[float]) -> float:
    return max(abs(float(value)) for value in values)


def mechanics_witness() -> dict:
    points = [
        (0.0, 0.0, 0.0),
        (0.2, -0.3, 0.1),
        (-0.4, 0.1, -0.2),
        (0.35, 0.25, 0.4),
    ]
    rows = []
    max_divergence_difference = 0.0
    max_curl_difference = 0.0
    min_velocity_difference = float("inf")
    max_expected_curl_error = 0.0
    for point in points:
        ja = _centered_jacobian(_field_a, point)
        jb = _centered_jacobian(_field_b, point)
        da = _divergence(ja)
        db = _divergence(jb)
        ca = _curl(ja)
        cb = _curl(jb)
        ua = _field_a(*point)
        ub = _field_b(*point)
        dv = tuple(float(b - a) for a, b in zip(ua, ub))
        max_divergence_difference = max(max_divergence_difference, abs(da - db))
        max_curl_difference = max(
            max_curl_difference, _max_abs(a - b for a, b in zip(ca, cb))
        )
        min_velocity_difference = min(
            min_velocity_difference, sum(value * value for value in dv) ** 0.5
        )
        max_expected_curl_error = max(
            max_expected_curl_error,
            _max_abs((ca[0], ca[1], ca[2] - 1.0, cb[0], cb[1], cb[2] - 1.0)),
        )
        rows.append(
            {
                "point": list(point),
                "divergence_A": da,
                "divergence_B": db,
                "curl_A": list(ca),
                "curl_B": list(cb),
                "velocity_difference": list(dv),
            }
        )

    tolerance = 1.0e-8
    if max_divergence_difference > tolerance:
        raise AssertionError("same-curl witness divergence equality drifted")
    if max_curl_difference > tolerance:
        raise AssertionError("same-curl witness curl equality drifted")
    if max_expected_curl_error > tolerance:
        raise AssertionError("same-curl witness no longer has curl=(0,0,1)")
    if min_velocity_difference < 0.299999999:
        raise AssertionError("same-curl witness lost its nonzero velocity ambiguity")

    return {
        "sample_count": len(points),
        "max_divergence_difference": max_divergence_difference,
        "max_curl_difference": max_curl_difference,
        "max_expected_curl_error": max_expected_curl_error,
        "min_velocity_difference": min_velocity_difference,
        "interpretation": (
            "two smooth axis-regular axisymmetric-with-swirl local fields have the same "
            "nonzero vorticity and the same zero divergence but different velocity"
        ),
    }


def audit(*, repo_root: Path | None = None) -> dict:
    root = _repo_root() if repo_root is None else Path(repo_root)
    config = _load_json(root / CONFIG_REL)
    constraints = _load_json(root / CONSTRAINTS_REL)
    upstream_bytes = (root / UPSTREAM_REL).read_bytes()
    constraints_bytes = (root / CONSTRAINTS_REL).read_bytes()

    if config.get("task") != TASK:
        raise AssertionError("task identity drift")
    if config.get("status") != "governance_only_not_scientific_admission":
        raise AssertionError("governance scope promoted")

    audited = config["audited_upstream"]
    if audited["head"] != "1d947d4c71812b59b6729befb2be9c46914852cb":
        raise AssertionError("audited #884 exact head drift")
    if _git_blob_sha1(upstream_bytes) != audited["source_git_blob_sha1"]:
        raise AssertionError("audited #884 source blob drift")
    upstream_text = upstream_bytes.decode("utf-8")
    for token in (
        '"vorticity_surface_registered": True',
        '"vorticity_sha256_registered": True',
        '"eligible_for_velocity_export_ready": False',
        '"usable_as_whole_domain_vorticity_morphology_evidence": False',
        '"complete_candidate_api_ready": False',
    ):
        if token not in upstream_text:
            raise AssertionError(f"upstream truth-boundary token missing: {token}")

    cr001 = config["canonical_cr001_binding"]
    if _git_blob_sha1(constraints_bytes) != cr001["git_blob_sha1"]:
        raise AssertionError("canonical CR001 blob drift")
    if constraints["nu"] != cr001["nu"]:
        raise AssertionError("nu drift")
    if constraints["domain"]["physical"] != cr001["physical_domain"]:
        raise AssertionError("physical-domain drift")
    if constraints["domain"]["evaluation_box"] != cr001["evaluation_box"]:
        raise AssertionError("evaluation-box drift")
    if constraints["domain"]["support"] != cr001["support"]:
        raise AssertionError("support drift")
    if constraints["domain"]["time_interval"] != cr001["time_interval"]:
        raise AssertionError("time-window drift")
    if constraints["forcing"]["mode"] != cr001["forcing_mode"]:
        raise AssertionError("forcing-family drift")
    thresholds = constraints["validation"]["thresholds"]
    for key in ("divergence_max", "divergence_L2", "pde_residual_max", "pde_residual_L2"):
        if thresholds[key] != cr001[key]:
            raise AssertionError(f"CR001 threshold drift: {key}")
    if "No residual-dependent basis or pointwise free force" not in constraints["forcing"]["restriction"]:
        raise AssertionError("free-force firewall drift")
    if "reject collapsed candidates" not in constraints["nontriviality"]["enforcement"]:
        raise AssertionError("amplitude-collapse firewall drift")
    if "changing thresholds requires a new experiment version" not in constraints["validation"]["failure_policy"]:
        raise AssertionError("post-hoc threshold firewall drift")

    provenance = config["provenance_classes"]
    if set(provenance) != {"user_requirement", "public_source_fact", "autonomous_design", "pending_unknown"}:
        raise AssertionError("four-way provenance classification drift")
    if any(not provenance[key] for key in provenance):
        raise AssertionError("empty provenance class")

    boundary = config["identity_boundary"]
    required_false = (
        "vorticity_sha256_is_velocity_candidate_identity",
        "matching_vorticity_alone_identifies_velocity",
        "matching_vorticity_and_divergence_locally_identifies_velocity_without_global_boundary_data",
        "strict_inner_vorticity_consistency_can_promote_whole_domain_velocity_identity",
        "strict_inner_vorticity_consistency_can_promote_visual_correspondence",
        "strict_inner_vorticity_consistency_can_promote_pde_validation",
    )
    for key in required_false:
        if boundary[key] is not False:
            raise AssertionError(f"forbidden identity promotion enabled: {key}")
    if boundary["velocity_identity_must_remain_bound_to_the_velocity_artifact_identity"] is not True:
        raise AssertionError("velocity artifact identity no longer authoritative")

    route = config["route_state"]
    for key in (
        "kokuno_velocity_export_ready",
        "kokuno_visual_correspondence_verified",
        "kokuno_pde_validated",
        "kokuno_paper_exact",
        "kokuno_openai_field_identified",
    ):
        if route[key] is not False:
            raise AssertionError(f"route truth state promoted: {key}")
    if route["canonical_eq45_velocity_export_ready_unaffected"] is not True:
        raise AssertionError("scoped Kokuno audit incorrectly blocks canonical Eq45 delivery")

    witness = mechanics_witness()
    return {
        "task": TASK,
        "status": "pass",
        "upstream_source_git_blob_sha1": audited["source_git_blob_sha1"],
        "canonical_constraints_git_blob_sha1": cr001["git_blob_sha1"],
        "mechanics_witness": witness,
        "scientific_admission_changed": False,
        "candidate_bytes_changed": False,
        "thresholds_changed": False,
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=None)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args(argv)
    receipt = audit(repo_root=args.repo_root)
    body = json.dumps(receipt, sort_keys=True, indent=2, allow_nan=False) + "\n"
    if args.output is None:
        print(body, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(body, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
