import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET

import pytest

from openai_ns_reconstruction.constrained_pvd_collection import write_pvd_time_collection


def _snapshot(path: Path, payload: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(payload, encoding="utf-8")
    return path


def test_writes_time_aware_single_part_collection_and_truth_sidecar(tmp_path: Path):
    frames = tmp_path / "frames"
    snapshots = [
        _snapshot(frames / "velocity_000.vtk", "frame 0\n"),
        _snapshot(frames / "velocity_001.vtk", "frame 1\n"),
        _snapshot(frames / "velocity_002.vtk", "frame 2\n"),
    ]
    times = [0.25, 0.5, 0.75]

    artifact = write_pvd_time_collection(
        snapshots,
        times,
        tmp_path / "velocity_series.pvd",
        candidate_id="eq45-seed-sha256:abc123",
    )

    root = ET.parse(artifact.pvd_path).getroot()
    assert root.attrib == {"type": "Collection", "version": "0.1", "byte_order": "LittleEndian"}
    datasets = root.findall("./Collection/DataSet")
    assert [float(node.attrib["timestep"]) for node in datasets] == times
    assert [node.attrib["part"] for node in datasets] == ["0", "0", "0"]
    assert [node.attrib["file"] for node in datasets] == [
        "frames/velocity_000.vtk",
        "frames/velocity_001.vtk",
        "frames/velocity_002.vtk",
    ]

    metadata = json.loads(artifact.metadata_path.read_text(encoding="utf-8"))
    assert metadata["claim_scope"] == "visualization_time_index_only"
    assert metadata["candidate_id"] == "eq45-seed-sha256:abc123"
    assert metadata["series_index_ready"] is True
    assert metadata["velocity_changed"] is False
    for claim in (
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        assert metadata[claim] is False
    assert metadata["pvd_sha256"] == hashlib.sha256(artifact.pvd_path.read_bytes()).hexdigest()
    for record, snapshot in zip(metadata["snapshots"], snapshots, strict=True):
        assert record["sha256"] == hashlib.sha256(snapshot.read_bytes()).hexdigest()


def test_output_is_deterministic_for_same_series(tmp_path: Path):
    frames = tmp_path / "frames"
    snapshots = [
        _snapshot(frames / "a.vtk", "a\n"),
        _snapshot(frames / "b.vtk", "b\n"),
    ]
    left = write_pvd_time_collection(
        snapshots,
        [0.25, 0.75],
        tmp_path / "left.pvd",
        candidate_id="candidate",
    )
    right = write_pvd_time_collection(
        snapshots,
        [0.25, 0.75],
        tmp_path / "right.pvd",
        candidate_id="candidate",
    )
    assert left.pvd_path.read_bytes() == right.pvd_path.read_bytes()


def test_fail_closed_on_bad_series_inputs(tmp_path: Path):
    frame = _snapshot(tmp_path / "frame.vtk", "frame\n")
    outside = _snapshot(tmp_path.parent / f"{tmp_path.name}-outside.vtk", "outside\n")

    with pytest.raises(ValueError, match="same length"):
        write_pvd_time_collection([frame], [0.25, 0.5], tmp_path / "a.pvd", candidate_id="c")
    with pytest.raises(ValueError, match="strictly increasing"):
        write_pvd_time_collection([frame, frame], [0.5, 0.25], tmp_path / "b.pvd", candidate_id="c")
    with pytest.raises(ValueError, match="finite"):
        write_pvd_time_collection([frame], [float("nan")], tmp_path / "c.pvd", candidate_id="c")
    with pytest.raises(ValueError, match="below the PVD output directory"):
        write_pvd_time_collection([outside], [0.25], tmp_path / "d.pvd", candidate_id="c")
    with pytest.raises(ValueError, match="non-empty"):
        write_pvd_time_collection([frame], [0.25], tmp_path / "e.pvd", candidate_id="")
    with pytest.raises(ValueError, match=".pvd"):
        write_pvd_time_collection([frame], [0.25], tmp_path / "e.xml", candidate_id="c")

    first = write_pvd_time_collection([frame], [0.25], tmp_path / "existing.pvd", candidate_id="c")
    assert first.pvd_path.exists()
    with pytest.raises(FileExistsError, match="overwrite"):
        write_pvd_time_collection([frame], [0.25], tmp_path / "existing.pvd", candidate_id="c")
