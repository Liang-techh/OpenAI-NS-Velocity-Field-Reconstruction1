"""Frozen reference adapter with an initial compact curl velocity correction.

This module is deliberately a small bridge between the assembled two-patch
field and :class:`SolenoidalScaleTransport`.  ``localized_drift_1800_fit``
stores its patch controls as potential *time derivatives*; those controls are
loaded by :func:`global_two_patch_candidate.load` and therefore vanish at its
frozen ``tau0``.  The separate ``scale_reference_velocity_step.json`` record
contains an explicitly labelled initial velocity correction.  Its coefficients
are inserted directly into ``basis_data`` at ``tau0`` and are never multiplied
by a time increment.

The canonical step record has the following small schema::

    {
      "status": "completed",
      "sources": {
        "two_patch_candidate": {"path": "localized_drift_1800_fit.json",
                                 "sha256": "..."},
        "geometry": {"path": "full_wave_frozen_cache.json",
                      "sha256": "..."}
      },
      "reference": {
        "tau0": 0.00024408986266228504,
        "degree": 2, "modes": [0, 1, 2],
        "patch_count": 2, "velocity_columns_per_patch": 135,
        "pressure_columns_per_patch": 45,
        "control_semantics": "initial_velocity",
        "velocity_coefficients_patch1": [[...135 values...]],
        "velocity_coefficients_patch2": [[...135 values...]],
        "pressure_coefficients_patch1": [[...45 values...]],
        "pressure_coefficients_patch2": [[...45 values...]]
      }
    }

For coordination with the grouped worker, top-level ``reference`` fields,
``selected`` fields, and the aliases ``initial_velocity_coefficients`` or
``velocity_correction_coefficients`` are accepted.  The provenance, geometry,
semantics, and coefficient dimensions remain mandatory.  A descriptor that
mentions a derivative or time slope is rejected rather than silently being
interpreted as an initial field.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path
from typing import Any, Iterable, Mapping

import numpy as np


ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from global_two_patch_candidate import load as load_two_patch  # noqa: E402
from supported_fourier_basis import basis_data  # noqa: E402


SOURCE_PATH = ROOT / "localized_drift_1800_fit.json"
GEOMETRY_PATH = ROOT / "full_wave_frozen_cache.json"
PRESSURE_PATH = ROOT / "scale_generator_pressure_projection.json"
STEP_PATH = ROOT / "scale_reference_velocity_step.json"
DEGREE = 2
MODES = (0, 1, 2)
Q = (DEGREE + 1) ** 2
CONTROL_COUNT_PER_PATCH = 180
VELOCITY_COUNT_PER_PATCH = 135
PRESSURE_COUNT_PER_PATCH = 45


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _mapping(value: Any, name: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be an object")
    return value


def _read_json(path: Path, name: str) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"Missing {name}: {path}")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise ValueError(f"Could not read {name}: {path}") from exc
    return dict(_mapping(value, name))


def _close(left: Any, right: Any, name: str, rtol: float = 2.0e-12,
           atol: float = 2.0e-14) -> None:
    a = np.asarray(left, dtype=float)
    b = np.asarray(right, dtype=float)
    if a.shape != b.shape or not np.allclose(a, b, rtol=rtol, atol=atol):
        raise ValueError(f"Step {name} does not match frozen geometry")


def _path_name(value: Any, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"Missing {name} path")
    return Path(value).name


def _descriptor(container: Mapping[str, Any], keys: Iterable[str], name: str) -> Mapping[str, Any]:
    """Find one provenance descriptor while keeping a strict required shape."""

    sources = container.get("sources")
    if isinstance(sources, Mapping):
        for key in keys:
            value = sources.get(key)
            if isinstance(value, Mapping):
                return value
        # Permit an exact basename key used by compact worker reports.
        for value in sources.values():
            if isinstance(value, Mapping) and Path(str(value.get("path", ""))).name == name:
                return value
    source_hashes = container.get("source_hashes")
    if isinstance(source_hashes, Mapping):
        for key in keys:
            value = source_hashes.get(key)
            if isinstance(value, Mapping):
                return value
            if isinstance(value, str) and Path(key).name == name:
                return {"path": key, "sha256": value}
        for value in source_hashes.values():
            if isinstance(value, Mapping) and Path(str(value.get("path", ""))).name == name:
                return value
    return {}


def _check_descriptor(step: Mapping[str, Any], keys: Iterable[str], expected: Path,
                      label: str) -> dict[str, str]:
    descriptor = _descriptor(step, keys, expected.name)
    # Also accept the un-nested source descriptor, but never accept a path or
    # hash by itself: both are needed to make this loader reproducible.
    if not descriptor:
        direct_path = step.get(f"{label}_path")
        direct_hash = step.get(f"{label}_sha256")
        if direct_path is not None or direct_hash is not None:
            descriptor = {"path": direct_path, "sha256": direct_hash}
    path_name = _path_name(descriptor.get("path"), f"{label} provenance")
    digest = descriptor.get("sha256")
    if not isinstance(digest, str) or len(digest) != 64:
        raise ValueError(f"Missing {label} provenance sha256")
    if path_name != expected.name:
        raise ValueError(f"Step {label} path {path_name!r} is not {expected.name!r}")
    actual = _sha256(expected)
    if digest.lower() != actual:
        raise ValueError(f"Step {label} hash mismatch for {expected.name}")
    return {"path": expected.name, "sha256": actual}


def _sections(step: Mapping[str, Any]) -> list[Mapping[str, Any]]:
    result: list[Mapping[str, Any]] = [step]
    for key in ("reference", "selected", "correction", "initial", "inputs"):
        value = step.get(key)
        if isinstance(value, Mapping):
            result.append(value)
    return result


def _first(sections: Iterable[Mapping[str, Any]], keys: Iterable[str], default: Any = None) -> Any:
    for section in sections:
        for key in keys:
            if key in section:
                return section[key]
    return default


def _normalise_semantics(value: Any) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError("Step is missing control_semantics")
    text = " ".join(value.lower().replace("_", " ").replace("-", " ").split())
    # The current worker report says that time-tangent controls are *not*
    # reused.  That wording is safe; an affirmative derivative/slope meaning
    # is not.
    negated = any(
        phrase in text
        for phrase in ("not reused", "not a derivative", "not derivative", "not a slope")
    )
    forbidden = ("derivative", "time slope", "time rate", "tangent", "d u d t", "du dt")
    if any(word in text for word in forbidden) and not negated:
        raise ValueError("Velocity correction controls must be initial amplitudes, not time derivatives")
    if "initial" not in text or "velocity" not in text:
        raise ValueError("control_semantics must identify an initial velocity correction")
    return text


def _as_control_block(value: Any, name: str) -> np.ndarray:
    array = np.asarray(value, dtype=float)
    if array.shape != (VELOCITY_COUNT_PER_PATCH,):
        raise ValueError(f"{name} must contain exactly {VELOCITY_COUNT_PER_PATCH} real controls")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} contains non-finite controls")
    return array.copy()


def _extract_velocity_blocks(sections: list[Mapping[str, Any]]) -> tuple[np.ndarray, np.ndarray]:
    aliases = (
        "velocity_coefficients",
        "initial_velocity_coefficients",
        "velocity_correction_coefficients",
        "initial_velocity_correction_coefficients",
        "reference_velocity_coefficients",
    )
    raw = _first(sections, aliases)
    if raw is not None:
        array = np.asarray(raw, dtype=float)
        if array.shape == (2, VELOCITY_COUNT_PER_PATCH):
            blocks = (array[0].copy(), array[1].copy())
        elif array.shape == (2 * VELOCITY_COUNT_PER_PATCH,):
            blocks = (array[:VELOCITY_COUNT_PER_PATCH].copy(),
                      array[VELOCITY_COUNT_PER_PATCH:].copy())
        else:
            raise ValueError("velocity_coefficients must have shape (2,135) or (270,)")
    else:
        blocks_list: list[np.ndarray] = []
        for index in (1, 2):
            value = _first(
                sections,
                (
                    f"patch{index}_velocity_coefficients",
                    f"patch{index}_initial_velocity_coefficients",
                    f"patch{index}_velocity_correction_coefficients",
                    f"velocity_coefficients_patch{index}",
                    f"initial_velocity_coefficients_patch{index}",
                    f"velocity_correction_coefficients_patch{index}",
                ),
            )
            if value is None:
                raise ValueError("Step is missing velocity correction coefficients")
            blocks_list.append(_as_control_block(value, f"patch{index}_velocity_coefficients"))
        blocks = (blocks_list[0], blocks_list[1])
    for index, block in enumerate(blocks, 1):
        if not np.all(np.isfinite(block)):
            raise ValueError(f"patch{index} velocity correction contains non-finite values")
    return blocks


def _unpack_pressure_only(value: Any, name: str) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
    """Decode optional 45-real pressure-only controls for one patch."""

    if value is None:
        return None
    array = np.asarray(value, dtype=float)
    if array.shape != (PRESSURE_COUNT_PER_PATCH,):
        raise ValueError(
            f"{name} must contain {PRESSURE_COUNT_PER_PATCH} real pressure controls when supplied"
        )
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} contains non-finite controls")
    cursor = 0
    result: list[np.ndarray] = []
    for mode in MODES:
        if mode == 0:
            result.append(array[cursor:cursor + Q].astype(complex))
            cursor += Q
        else:
            result.append(array[cursor:cursor + 2 * Q:2] + 1j * array[cursor + 1:cursor + 2 * Q:2])
            cursor += 2 * Q
    if cursor != len(array):
        raise ValueError(f"{name} unpack failed")
    return tuple(result)  # type: ignore[return-value]


def _extract_pressure_blocks(sections: list[Mapping[str, Any]]) -> tuple[tuple[np.ndarray, ...], ...] | None:
    aliases = (
        "pressure_coefficients",
        "initial_pressure_coefficients",
        "pressure_correction_coefficients",
    )
    raw = _first(sections, aliases)
    if raw is None:
        per_patch = []
        found = False
        for index in (1, 2):
            value = _first(
                sections,
                (
                    f"patch{index}_pressure_coefficients",
                    f"patch{index}_initial_pressure_coefficients",
                    f"patch{index}_pressure_correction_coefficients",
                    f"pressure_coefficients_patch{index}",
                    f"initial_pressure_coefficients_patch{index}",
                    f"pressure_correction_coefficients_patch{index}",
                ),
            )
            if value is not None:
                found = True
            per_patch.append(value)
        if not found:
            return None
        if any(value is None for value in per_patch):
            raise ValueError("Both pressure correction patch blocks are required when pressure is supplied")
        raw = per_patch
    array = np.asarray(raw, dtype=float)
    if array.shape == (2, PRESSURE_COUNT_PER_PATCH):
        values = [array[0], array[1]]
    elif array.shape == (2 * PRESSURE_COUNT_PER_PATCH,):
        values = [array[:PRESSURE_COUNT_PER_PATCH], array[PRESSURE_COUNT_PER_PATCH:]]
    else:
        raise ValueError(
            f"pressure_coefficients must have shape (2,{PRESSURE_COUNT_PER_PATCH}) or "
            f"({2 * PRESSURE_COUNT_PER_PATCH},)"
        )
    unpacked = [_unpack_pressure_only(value, f"pressure patch {i}") for i, value in enumerate(values, 1)]
    return (tuple(unpacked[0]), tuple(unpacked[1]))  # type: ignore[arg-type]


def _pressure_baseline(path: Path) -> tuple[tuple[np.ndarray, ...], ...]:
    report = _read_json(path, "pressure projection")
    if report.get("status") != "completed":
        raise ValueError("Pressure projection is not completed")
    fit = _mapping(report.get("fit"), "pressure projection fit")
    coefficients = fit.get("coefficients")
    if coefficients is None:
        first = fit.get("first_patch_coefficients")
        second = fit.get("second_patch_coefficients")
        if first is None or second is None:
            raise ValueError("Pressure projection has no coefficients")
        coefficients = list(first) + list(second)
    array = np.asarray(coefficients, dtype=float)
    if array.shape != (2 * PRESSURE_COUNT_PER_PATCH,):
        raise ValueError("Pressure projection coefficients must contain 90 values")
    first = _unpack_pressure_only(array[:PRESSURE_COUNT_PER_PATCH], "pressure baseline patch 1")
    second = _unpack_pressure_only(array[PRESSURE_COUNT_PER_PATCH:], "pressure baseline patch 2")
    return (tuple(first), tuple(second))  # type: ignore[arg-type]


def _unpack_velocity_only(value: Any, name: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Decode the 135-real velocity layout used by the step worker.

    Mode 0 has 27 real potential-component columns.  Modes 1 and 2 each have
    27 complex columns packed as real/imag pairs, yielding
    ``27 + 54 + 54 = 135`` values per patch.  Pressure columns are separate
    and are intentionally not included in this velocity decoder.
    """

    array = np.asarray(value, dtype=float)
    if array.shape != (VELOCITY_COUNT_PER_PATCH,):
        raise ValueError(
            f"{name} must contain exactly {VELOCITY_COUNT_PER_PATCH} real velocity controls"
        )
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} contains non-finite controls")
    cursor = 0
    result: list[np.ndarray] = []
    for mode in MODES:
        if mode == 0:
            result.append(array[cursor:cursor + 3 * Q].astype(complex))
            cursor += 3 * Q
        else:
            count = 3 * Q
            result.append(array[cursor:cursor + 2 * count:2]
                          + 1j * array[cursor + 1:cursor + 2 * count:2])
            cursor += 2 * count
    if cursor != len(array):
        raise ValueError(f"{name} unpack failed")
    return tuple(result)  # type: ignore[return-value]


def _check_step_geometry(step: Mapping[str, Any], sections: list[Mapping[str, Any]],
                         source_report: Mapping[str, Any], geometry: Mapping[str, Any]) -> dict[str, Any]:
    reference = geometry.get("inputs")
    reference = _mapping(reference, "geometry inputs")
    mean = _mapping(reference.get("mean"), "geometry mean")
    wave = _mapping(reference.get("wave"), "geometry wave")
    tau0 = float(mean["tau"])
    degree = _first(sections, ("degree", "patch_degree"))
    modes = _first(sections, ("modes", "patch_modes"))
    patch_count = _first(sections, ("patch_count",))
    per_patch = _first(sections, ("control_count_per_patch", "controls_per_patch"))
    if per_patch is None:
        velocity_count = _first(sections, ("velocity_columns_per_patch",))
        pressure_count = _first(sections, ("pressure_columns_per_patch",))
        if velocity_count is not None and pressure_count is not None:
            per_patch = int(velocity_count) + int(pressure_count)
    if degree is None or int(degree) != DEGREE:
        raise ValueError("Step degree must be 2")
    if modes is None or tuple(int(v) for v in modes) != MODES:
        raise ValueError("Step modes must be [0, 1, 2]")
    if patch_count is None or int(patch_count) != 2:
        raise ValueError("Step patch_count must be 2")
    if per_patch is None or int(per_patch) != CONTROL_COUNT_PER_PATCH:
        raise ValueError("Step control_count_per_patch must be 180")
    step_tau = _first(sections, ("tau0", "reference_tau0", "initial_tau"))
    if step_tau is None:
        raise ValueError("Step is missing frozen tau0")
    if not np.isclose(float(step_tau), tau0, rtol=2.0e-12, atol=2.0e-14):
        raise ValueError("Step tau0 does not match frozen geometry")
    # If the worker records patch geometry, require exact agreement with the
    # immutable source report.  The actual values remain sourced from that
    # report so a step cannot silently move the support.
    expected_patches = (
        (source_report["inputs"]["first_patch_center"], source_report["inputs"]["first_patch_widths"]),
        (source_report["inputs"]["second_patch_center"], source_report["inputs"]["second_patch_widths"]),
    )
    patches = _first(sections, ("patches",))
    if patches is not None:
        if isinstance(patches, Mapping):
            patches = [patches.get("first"), patches.get("second")]
        if not isinstance(patches, (list, tuple)) or len(patches) != 2:
            raise ValueError("Step patches must contain two patch descriptors")
        for index, (patch, (center, widths)) in enumerate(zip(patches, expected_patches), 1):
            patch = _mapping(patch, f"patch {index}")
            if "center" not in patch or "widths" not in patch:
                raise ValueError(f"Step patch {index} is incomplete")
            _close(patch["center"], center, f"patch {index} center")
            _close(patch["widths"], widths, f"patch {index} widths")
    geometry_carrier = wave.get("carrier")
    step_carrier = _first(sections, ("carrier", "wave_carrier"))
    if step_carrier is not None:
        _close(step_carrier, geometry_carrier, "carrier")
    return {"tau0": tau0, "degree": DEGREE, "modes": list(MODES),
            "patch_count": 2, "control_count_per_patch": CONTROL_COUNT_PER_PATCH}


def _check_step(step_path: Path, source_path: Path, geometry_path: Path,
                pressure_path: Path, source_report: Mapping[str, Any],
                geometry: Mapping[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    step = _read_json(step_path, "scale reference step")
    if step.get("status") != "completed":
        raise ValueError("Scale reference step is not completed")
    if step.get("accepted") is True:
        raise ValueError("An accepted recursive candidate cannot be imported by this diagnostic adapter")
    source_hash = _check_descriptor(
        step,
        (
            "two_patch_candidate", "candidate", "localized_drift_1800_fit",
            "drift_report", "reference",
        ),
        source_path,
        "two_patch_candidate",
    )
    geometry_hash = _check_descriptor(
        step,
        ("geometry", "full_wave_frozen_cache", "frozen_geometry"),
        geometry_path,
        "geometry",
    )
    sections = _sections(step)
    semantics = _normalise_semantics(
        _first(
            sections,
            ("control_semantics", "coefficient_semantics", "correction_kind", "semantics"),
        )
    )
    geometry_meta = _check_step_geometry(step, sections, source_report, geometry)
    velocity_blocks = _extract_velocity_blocks(sections)
    pressure_blocks = _extract_pressure_blocks(sections)
    pressure_descriptor = _descriptor(
        step, ("pressure_report", "pressure_projection", "pressure_source"),
        pressure_path.name,
    )
    pressure_baseline = None
    pressure_hash = None
    if pressure_descriptor:
        pressure_hash = _check_descriptor(
            step, ("pressure_report", "pressure_projection", "pressure_source"),
            pressure_path, "pressure_report",
        )
        pressure_baseline = _pressure_baseline(pressure_path)
    elif pressure_blocks is not None:
        raise ValueError("Selected pressure correction requires pressure baseline provenance")
    parsed = {
        "step_path": str(step_path),
        "source_hashes": {
            source_path.name: source_hash["sha256"],
            geometry_path.name: geometry_hash["sha256"],
        },
        "control_semantics": semantics,
        "geometry": geometry_meta,
        "velocity_coefficients": [block.tolist() for block in velocity_blocks],
        "pressure_coefficients_present": pressure_blocks is not None,
        "pressure_baseline_present": pressure_baseline is not None,
    }
    if pressure_hash is not None:
        parsed["source_hashes"][pressure_path.name] = pressure_hash["sha256"]
    return step, {
        "meta": parsed,
        "velocity_blocks": velocity_blocks,
        "pressure_blocks": pressure_blocks,
        "pressure_baseline": pressure_baseline,
    }


def _points(points: Any) -> np.ndarray:
    array = np.asarray(points, dtype=float)
    if array.ndim == 1 and array.shape == (3,):
        array = array[None, :]
    if array.ndim != 2 or array.shape[1] != 3:
        raise ValueError("points must have shape (N, 3)")
    if not np.all(np.isfinite(array)):
        raise ValueError("points must be finite")
    return array


class _BaseReferenceView:
    """Expose the assembled base with the optional pressure projection only."""

    def __init__(self, owner: "ScaleReferenceCandidate"):
        self._owner = owner
        self.tau0 = owner.tau0
        self.h = owner.h
        self.nu = owner.nu

    def velocity(self, points: Any) -> np.ndarray:
        points_array = _points(points)
        return np.asarray(
            self._owner.base_field.fields(points_array, self.tau0)[0], dtype=float
        )

    def pressure(self, points: Any) -> np.ndarray:
        points_array = _points(points)
        base = np.asarray(
            self._owner.base_field.fields(points_array, self.tau0)[1], dtype=float
        )
        return base + self._owner._pressure_from_modes(
            points_array, self._owner._pressure_baseline_modes)

    def fields(self, points: Any, tau: float | None = None) -> tuple[np.ndarray, np.ndarray]:
        if tau is not None and not np.isclose(float(tau), self.tau0, rtol=0.0, atol=2.0e-15):
            raise ValueError("Base reference view is frozen at tau0")
        points_array = _points(points)
        return self.velocity(points_array), self.pressure(points_array)


class ScaleReferenceCandidate:
    """Velocity-only frozen reference with optional pressure diagnostics."""

    def __init__(self, base_field: Any, localized: Any, snapshot: Mapping[str, Any],
                 candidate: Mapping[str, Any], step: Mapping[str, Any], parsed: Mapping[str, Any]):
        self.base_field = base_field
        self.localized = localized
        self.snapshot = snapshot
        self.candidate = candidate
        self.step = step
        self.tau0 = float(parsed["meta"]["geometry"]["tau0"])
        self.nu = float(snapshot["inputs"]["mean"]["nu"])
        self.h = float(localized.inner.h)
        self.source_hashes = dict(parsed["meta"]["source_hashes"])
        self.control_semantics = str(parsed["meta"]["control_semantics"])
        self.degree = DEGREE
        self.modes = MODES
        self.center = tuple(float(v) for v in candidate["inputs"]["first_patch_center"]), tuple(float(v) for v in candidate["inputs"]["second_patch_center"])
        self.widths = tuple(float(v) for v in candidate["inputs"]["first_patch_widths"]), tuple(float(v) for v in candidate["inputs"]["second_patch_widths"])
        carrier = np.asarray(snapshot["inputs"]["wave"]["carrier"], dtype=float)
        self.carriers = ({0: np.zeros(2), 1: carrier, 2: 2.0 * carrier},
                         {0: np.zeros(2), 1: carrier, 2: 2.0 * carrier})
        blocks = parsed["velocity_blocks"]
        self._velocity_modes = tuple(_unpack_velocity_only(block, f"patch {index} velocity")
                                     for index, block in enumerate(blocks, 1))
        pressure_blocks = parsed["pressure_blocks"]
        pressure_baseline = parsed.get("pressure_baseline")
        self._pressure_baseline_modes = pressure_baseline
        if pressure_blocks is None:
            self._pressure_modes = pressure_baseline
        elif pressure_baseline is None:
            self._pressure_modes = pressure_blocks
        else:
            self._pressure_modes = tuple(
                tuple(
                    pressure_baseline[patch][mode] + pressure_blocks[patch][mode]
                    for mode in range(len(MODES))
                )
                for patch in range(2)
            )
        self.velocity_coefficients = tuple(np.asarray(block, dtype=float) for block in blocks)
        # This view is useful for a replay that wants to separate the new
        # initial velocity correction from the already fitted pressure defect.
        self.base_reference = _BaseReferenceView(self)

    def _correction_velocity(self, points: np.ndarray) -> np.ndarray:
        result = np.zeros((len(points), 3), dtype=float)
        for patch_index in range(2):
            modes = self._velocity_modes[patch_index]
            for mode in MODES:
                values = basis_data(
                    points, self.center[patch_index], self.widths[patch_index],
                    mode, DEGREE, self.carriers[patch_index][mode],
                )[0]
                result += np.einsum("ncq,q->nc", values, modes[mode]).real
        return result

    def _correction_pressure(self, points: np.ndarray) -> np.ndarray:
        return self._pressure_from_modes(points, self._pressure_modes)

    def _pressure_from_modes(self, points: np.ndarray, pressure_modes) -> np.ndarray:
        result = np.zeros(len(points), dtype=float)
        if pressure_modes is None:
            return result
        for patch_index in range(2):
            modes = pressure_modes[patch_index]
            for mode in MODES:
                values = basis_data(
                    points, self.center[patch_index], self.widths[patch_index],
                    mode, DEGREE, self.carriers[patch_index][mode],
                )[1]
                result += np.einsum("nq,q->n", values, modes[mode]).real
        return result

    def velocity(self, points: Any) -> np.ndarray:
        points_array = _points(points)
        base_velocity = np.asarray(self.base_field.fields(points_array, self.tau0)[0], dtype=float)
        if base_velocity.shape != (len(points_array), 3):
            raise ValueError("Base field returned an invalid velocity shape")
        return base_velocity + self._correction_velocity(points_array)

    def pressure(self, points: Any) -> np.ndarray:
        points_array = _points(points)
        base_pressure = np.asarray(self.base_field.fields(points_array, self.tau0)[1], dtype=float)
        if base_pressure.shape != (len(points_array),):
            raise ValueError("Base field returned an invalid pressure shape")
        return base_pressure + self._correction_pressure(points_array)

    def fields(self, points: Any, tau: float | None = None) -> tuple[np.ndarray, np.ndarray]:
        if tau is not None and not np.isclose(float(tau), self.tau0, rtol=0.0, atol=2.0e-15):
            raise ValueError("ScaleReferenceCandidate is frozen at tau0; no time derivative is assigned")
        points_array = _points(points)
        return self.velocity(points_array), self.pressure(points_array)


def load_reference(step: Path | str = STEP_PATH, source: Path | str = SOURCE_PATH,
                   geometry: Path | str = GEOMETRY_PATH,
                   pressure: Path | str = PRESSURE_PATH) -> ScaleReferenceCandidate:
    """Load and validate the frozen corrected reference field.

    The step JSON is mandatory.  All coefficients are initial amplitudes at
    ``tau0``; this function intentionally does not expose a time derivative for
    them.
    """

    step_path = Path(step)
    source_path = Path(source)
    geometry_path = Path(geometry)
    pressure_path = Path(pressure)
    source_report = _read_json(source_path, "two-patch source")
    geometry_report = _read_json(geometry_path, "frozen geometry")
    step_report, parsed = _check_step(
        step_path, source_path, geometry_path, pressure_path,
        source_report, geometry_report,
    )
    base_field, localized, snapshot, candidate, _ = load_two_patch(source_path)
    _close(snapshot["inputs"]["mean"]["tau"], parsed["meta"]["geometry"]["tau0"], "tau0")
    return ScaleReferenceCandidate(
        base_field, localized, snapshot, candidate, step_report, parsed
    )


load = load_reference


if __name__ == "__main__":
    # Do not perform a grid replay here.  This command only exercises strict
    # provenance and reports the adapter metadata for the next worker.
    reference = load_reference()
    print(json.dumps({
        "status": "loaded",
        "tau0": reference.tau0,
        "control_semantics": reference.control_semantics,
        "source_hashes": reference.source_hashes,
        "pressure_correction": reference._pressure_modes is not None,
        "accepted": False,
        "scale_recursion_established": False,
    }, sort_keys=True))
