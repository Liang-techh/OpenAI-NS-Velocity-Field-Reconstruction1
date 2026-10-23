"""Truth-bounded ParaView PVD time-series indexing for velocity snapshots.

This module does not change or validate a velocity field. It only groups
already-written VTK-family snapshot files behind one deterministic PVD time
index so visualization tools can recover the declared physical time values.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
import xml.etree.ElementTree as ET


CLAIM_SCOPE = "visualization_time_index_only"


@dataclass(frozen=True)
class PVDCollectionArtifact:
    pvd_path: Path
    metadata_path: Path
    pvd_sha256: str
    times: tuple[float, ...]
    relative_files: tuple[str, ...]
    candidate_id: str
    claim_scope: str = CLAIM_SCOPE
    visualization_ready: bool = False
    visual_correspondence_verified: bool = False
    pde_validated: bool = False
    paper_exact: bool = False
    openai_field_identified: bool = False
    blowup_proved: bool = False


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _validated_times(times: object) -> tuple[float, ...]:
    try:
        values = tuple(float(value) for value in times)  # type: ignore[arg-type]
    except (TypeError, ValueError) as exc:
        raise ValueError("times must be a finite one-dimensional numeric sequence") from exc
    if not values:
        raise ValueError("at least one time value is required")
    if not all(math.isfinite(value) for value in values):
        raise ValueError("times must be finite")
    if any(b <= a for a, b in zip(values, values[1:])):
        raise ValueError("times must be strictly increasing")
    return values


def _portable_relative_files(snapshot_paths: object, output_parent: Path) -> tuple[tuple[Path, ...], tuple[str, ...]]:
    try:
        paths = tuple(Path(value) for value in snapshot_paths)  # type: ignore[arg-type]
    except TypeError as exc:
        raise ValueError("snapshot_paths must be a sequence of paths") from exc
    if not paths:
        raise ValueError("at least one snapshot path is required")

    output_root = output_parent.resolve()
    resolved: list[Path] = []
    relative: list[str] = []
    for path in paths:
        full = path.resolve()
        if not full.is_file():
            raise ValueError(f"snapshot does not exist as a file: {path}")
        try:
            rel = full.relative_to(output_root)
        except ValueError as exc:
            raise ValueError("all snapshots must live at or below the PVD output directory") from exc
        if full in resolved:
            raise ValueError("snapshot paths must be unique")
        resolved.append(full)
        relative.append(rel.as_posix())
    return tuple(resolved), tuple(relative)


def _pvd_bytes(times: tuple[float, ...], relative_files: tuple[str, ...]) -> bytes:
    root = ET.Element(
        "VTKFile",
        {
            "type": "Collection",
            "version": "0.1",
            "byte_order": "LittleEndian",
        },
    )
    collection = ET.SubElement(root, "Collection")
    for time, filename in zip(times, relative_files, strict=True):
        ET.SubElement(
            collection,
            "DataSet",
            {
                "timestep": format(time, ".17g"),
                "group": "",
                "part": "0",
                "file": filename,
            },
        )
    ET.indent(root, space="  ")
    body = ET.tostring(root, encoding="utf-8", xml_declaration=True)
    return body + b"\n"


def write_pvd_time_collection(
    snapshot_paths: object,
    times: object,
    output_path: str | Path,
    *,
    candidate_id: str,
) -> PVDCollectionArtifact:
    """Write a deterministic single-part PVD time index plus a truth-boundary sidecar.

    The snapshot files are treated as immutable visualization inputs. They are
    not parsed as velocity data and no numerical or PDE property is inferred
    from a successful index write.
    """

    if not isinstance(candidate_id, str) or not candidate_id.strip():
        raise ValueError("candidate_id must be a non-empty string")

    pvd_path = Path(output_path)
    if pvd_path.suffix.lower() != ".pvd":
        raise ValueError("output_path must end in .pvd")
    pvd_path.parent.mkdir(parents=True, exist_ok=True)
    metadata_path = pvd_path.with_suffix(pvd_path.suffix + ".meta.json")
    if pvd_path.exists() or metadata_path.exists():
        raise FileExistsError("refusing to overwrite an existing PVD collection or sidecar")

    time_values = _validated_times(times)
    snapshots, relative_files = _portable_relative_files(snapshot_paths, pvd_path.parent)
    if len(snapshots) != len(time_values):
        raise ValueError("snapshot_paths and times must have the same length")

    pvd_path.write_bytes(_pvd_bytes(time_values, relative_files))
    pvd_sha256 = _sha256(pvd_path)

    metadata = {
        "schema_version": 1,
        "claim_scope": CLAIM_SCOPE,
        "candidate_id": candidate_id.strip(),
        "pvd_file": pvd_path.name,
        "pvd_sha256": pvd_sha256,
        "series_index_ready": True,
        "snapshots": [
            {
                "time": time,
                "file": filename,
                "sha256": _sha256(path),
            }
            for time, filename, path in zip(time_values, relative_files, snapshots, strict=True)
        ],
        "velocity_changed": False,
        "visualization_ready": False,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
    }
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    return PVDCollectionArtifact(
        pvd_path=pvd_path,
        metadata_path=metadata_path,
        pvd_sha256=pvd_sha256,
        times=time_values,
        relative_files=relative_files,
        candidate_id=candidate_id.strip(),
    )
