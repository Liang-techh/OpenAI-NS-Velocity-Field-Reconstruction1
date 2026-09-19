"""CR002 audit for ST052 sampled-grid tensor semantics.

The integrated exporter already binds the numerical payload to a whole-candidate
identity.  This audit guards a narrower representation boundary: because x, y,
and z use equal 33-point coordinate vectors, array shape alone cannot establish
which tensor dimension is the physical x/y/z axis (or how a downstream consumer
interprets u/v/w as Cartesian components).

This module changes no candidate and evaluates no Navier--Stokes acceptance gate.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from . import st052_grid_export as grid_export


CONTRACT_FILENAME = "st052m_grid_tensor_semantics_contract.json"


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _load_json(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise AssertionError(f"expected JSON object: {path}")
    return data


def _assert_equal(actual: Any, expected: Any, label: str) -> None:
    if actual != expected:
        raise AssertionError(f"{label} drifted: {actual!r} != {expected!r}")


def axis_permutation_witness() -> dict[str, Any]:
    """Return a deterministic witness that shape equality does not bind axis meaning."""
    x = np.linspace(-2.0, 2.0, 33, dtype=float)
    y = x.copy()
    z = x.copy()
    xx, yy, zz = np.meshgrid(x, y, z, indexing="ij")
    field = xx + 10.0 * yy + 100.0 * zz
    swapped = np.swapaxes(field, 0, 1)

    if field.shape != swapped.shape or field.shape != (33, 33, 33):
        raise AssertionError("permutation witness must preserve the frozen cubic shape")
    if not (np.array_equal(x, y) and np.array_equal(y, z)):
        raise AssertionError("witness requires equal frozen coordinate vectors")

    probe = (3, 9, 17)
    original_value = float(field[probe])
    swapped_value = float(swapped[probe])
    difference = abs(original_value - swapped_value)
    if not difference > 0.0:
        raise AssertionError("axis permutation witness unexpectedly became symmetric")

    return {
        "shape": list(field.shape),
        "coordinate_vectors_equal": True,
        "permutation": "swap x and y tensor axes",
        "probe_index": list(probe),
        "original_value": original_value,
        "permuted_value": swapped_value,
        "absolute_difference": difference,
        "shape_only_semantics_sufficient": False,
    }


def audit_grid_tensor_semantics(repo_root: str | Path | None = None) -> dict[str, Any]:
    """Fail closed on live exporter, truth-state, or CR001 semantic drift."""
    root = Path(repo_root) if repo_root is not None else _repo_root()
    contract = _load_json(root / "configs" / CONTRACT_FILENAME)
    constraints = _load_json(root / "configs" / "constraints.json")

    live = contract["live_export_contract"]
    _assert_equal(grid_export.SCHEMA, live["schema"], "grid export schema")
    _assert_equal(grid_export.DEFAULT_GRID_POINTS, live["grid_points_per_axis"], "grid points")
    _assert_equal(list(grid_export.DEFAULT_TIMES), live["times"], "grid times")
    _assert_equal(list(grid_export.DEFAULT_DOMAIN), live["domain"], "grid domain")

    # These strings are part of the currently integrated exporter manifest and
    # are the machine-readable semantic statements a downstream consumer must
    # preserve rather than infer from shape alone.
    source_text = (root / "src" / "openai_ns_reconstruction" / "st052_grid_export.py").read_text(
        encoding="utf-8"
    )
    if '"array_layout": "u,v,w: [time,x,y,z]"' not in source_text:
        raise AssertionError("live export no longer declares [time,x,y,z] layout")
    if '"meshgrid_indexing": "ij"' not in source_text:
        raise AssertionError("live export no longer declares meshgrid indexing=ij")
    if '"grid_export_payload_checksum_bound_to_whole_candidate": True' not in source_text:
        raise AssertionError("live export no longer declares payload/candidate checksum binding")

    state = contract["verified_state"]
    for key in (
        "export_manifest_axis_layout_declared",
        "export_manifest_meshgrid_indexing_declared",
        "npz_mat_payload_checksum_bound",
    ):
        if state.get(key) is not True:
            raise AssertionError(f"verified positive fact must remain true: {key}")
    for key in (
        "shape_only_consumer_axis_semantics_verified",
        "matlab_octave_consumer_layout_metadata_enforced",
        "component_basis_machine_bound_in_consumer",
        "axis_permutation_equivalent_to_same_physical_field",
        "derived_grid_vorticity_semantics_verified",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        if state.get(key) is not False:
            raise AssertionError(f"unverified semantic/scientific claim promoted: {key}")

    frozen = contract["cr001_immutable"]
    _assert_equal(constraints["nu"], frozen["nu"], "nu")
    _assert_equal(constraints["domain"]["physical"], frozen["physical_domain"], "physical domain")
    _assert_equal(constraints["domain"]["evaluation_box"], frozen["evaluation_box"], "evaluation box")
    _assert_equal(constraints["domain"]["support"], frozen["support"], "support")
    _assert_equal(constraints["domain"]["time_interval"], frozen["time_interval"], "time interval")
    _assert_equal(constraints["forcing"]["mode"], frozen["forcing_mode"], "forcing mode")
    _assert_equal(constraints["forcing"]["parameters"], frozen["forcing_parameter_bounds"], "forcing bounds")
    if "No residual-dependent basis or pointwise free force" not in constraints["forcing"]["restriction"]:
        raise AssertionError("residual-defined free forcing prohibition drifted")
    _assert_equal(constraints["nontriviality"]["reference_energy"], frozen["reference_energy"], "reference energy")
    _assert_equal(
        constraints["nontriviality"]["reference_energy_abs_tolerance"],
        frozen["reference_energy_abs_tolerance"],
        "reference energy tolerance",
    )
    validation = constraints["validation"]
    _assert_equal(validation["seed"], frozen["validation_seed"], "validation seed")
    _assert_equal(validation["held_out_points"], frozen["held_out_points"], "held-out points")
    _assert_equal(validation["times"], frozen["validation_times"], "validation times")
    _assert_equal(validation["derivative_steps"], frozen["derivative_steps"], "derivative ladder")
    _assert_equal(
        validation["quadrature_orders_per_axis"],
        frozen["quadrature_orders_per_axis"],
        "quadrature ladder",
    )
    thresholds = validation["thresholds"]
    for key in ("divergence_max", "divergence_L2", "pde_residual_max", "pde_residual_L2"):
        _assert_equal(thresholds[key], frozen[key], key)
    if "reject collapsed candidates" not in constraints["nontriviality"]["enforcement"]:
        raise AssertionError("amplitude-collapse prohibition drifted")

    witness = axis_permutation_witness()
    return {
        "schema": contract["schema"],
        "task_id": contract["task_id"],
        "live_export_contract_checked": True,
        "cr001_checked_unchanged": True,
        "axis_permutation_witness": witness,
        "truth_boundary": state,
    }
