"""Extract bounded admitted source ports from saved evidence, without math."""
from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "experiments/root_st073"
PREFIX = "lei_ren_part1_paper_compliant_current_original_"
ARCHIVE = PREFIX + "Rm_loop_density_identity.json.gz"
ARCHIVE_HASH = "29c4fec78d3c622b36bb10867792b3d4df1b62f332cf3bc93cd5ed684706fe29"
DECODED_HASH = "ab2286f3b3142494b88d903d72c965d3d9f15f34f1d52478fb7cfd20919f67c1"
ANCESTOR = PREFIX + "Rm_generic_native_identity.json"
ANCESTOR_HASH = "71db4b9121e591878b4047025ed6d98e2a9a78d3659d82430ecfc294c0cafbf2"
OUTPUT = BASE / "original_Rm_loop_ports_admission.json"


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                    ensure_ascii=True, allow_nan=False).encode()).hexdigest()


def extract():
    compressed = (BASE / ARCHIVE).read_bytes()
    assert hashlib.sha256(compressed).hexdigest() == ARCHIVE_HASH
    decoded = gzip.decompress(compressed)
    assert hashlib.sha256(decoded).hexdigest() == DECODED_HASH
    data = json.loads(decoded)
    ancestor_bytes = (BASE / ANCESTOR).read_bytes()
    assert hashlib.sha256(ancestor_bytes).hexdigest() == ANCESTOR_HASH
    ancestor = json.loads(ancestor_bytes)
    result, phase = data["actual_Rm_loop_density_identity"], data["actual_Rm_common_phase_source"]
    assert result["passed"] and phase["passed"] and not result["unresolved_gates"]
    assert len(result["rows"]) == 22 and all(row["status"] == "proved" for row in result["rows"])
    inherited = ("actual_Rm_generic_native_identity", "actual_Rm_history_native_source", "actual_Rm_absolute_m15",
                 "corrected_own_swirl_provider", "corrected_reserved_point_numeric", "current_Rm_pressure_inlet_identity")
    assert all(data[key] == ancestor[key] for key in inherited)
    rows = []
    for row in result["rows"]:
        rows.append(dict(field=row["field"], order=row["order"], status=row["status"],
            postrestore={key:row["postrestore"][key] for key in ("node", "value")},
            patch={key:row["patch"][key] for key in ("node", "value")},
            complete_row_digest=digest(row)))
    packet = dict(schema="admitted_Rm_loop_source_ports_v1", passed=True, source_only=True,
        source_archive=ARCHIVE, source_archive_sha256=ARCHIVE_HASH, decoded_archive_sha256=DECODED_HASH,
        ancestor_payload_report=ANCESTOR, ancestor_payload_sha256=ANCESTOR_HASH,
        actual_loop_result_digest=digest(result), actual_loop_checkpoint_sha256="08ee51c4e8657017dad693a8bc1ae8c08a503c820957b7db4335b52e9d74fc1d",
        nodes_digest=result["nodes_digest"], graph_node_count=205200, ordinary_orders=[0, 1], rows=rows,
        common_phase_source=phase, loop_input_packet=data["actual_Rm_loop_density_inputs"],
        inherited_packet_digests={key:digest(data[key]) for key in inherited},
        input_hashes=data["input_hashes"], implementation_hashes=data["implementation_hashes"],
        extraction_implementation_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        exact_saved_data_projection=True, mathematical_operator_executed=False,
        numerical_values_computed=False, native_loop_seam_proved=False, full_goal_complete=False)
    raw = (json.dumps(packet, indent=2, ensure_ascii=True, allow_nan=False) + "\n").encode()
    with OUTPUT.open("xb") as stream:
        stream.write(raw)
    return dict(passed=True, output=str(OUTPUT.relative_to(ROOT)), bytes=len(raw),
                sha256=hashlib.sha256(raw).hexdigest(), admitted_rows=len(rows), mathematics_executed=False)


if __name__ == "__main__":
    print(json.dumps(extract()))
