"""Prepare and inspect the original outer-pressure enclosure handoff.

The source preflight and the later ``--write-report``/``--inspect-report``
stages are data-only: they authenticate saved reports, extract the two
ordinary-Z endpoint rows and accepted signed D rows, and validate immutable
JSON checkpoints.  Only a root invocation of ``--compute-new`` imports the
three reviewed new kernels after that preflight and runs the 1600-bit-context
production chain.  Accepted producers, checkers, fixtures, and helper
mathematics are never replayed.
"""
from __future__ import annotations

import argparse
import ast
import copy
import hashlib
import json
from pathlib import Path
import time


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "experiments/root_st073"

PARENT_NAME = "lei_ren_part1_paper_compliant_current_original_pressure_prefix_contraction.json"
PARENT_PATH = "experiments/root_st073/" + PARENT_NAME
PARENT_BYTES = 1_018_242
PARENT_SHA256 = "a10279cd63a6268da3269d33adc1e6de9f8171fd8fc13b8fb7f142d1845077b3"

NATIVE_NAME = "lei_ren_part1_paper_compliant_current_original_native_directed_controls_v2.json"
NATIVE_PATH = "experiments/root_st073/" + NATIVE_NAME
NATIVE_BYTES = 42_177_499
NATIVE_SHA256 = "10abe68b98fdf87709a7a7c9f27103523701e746c357ae62386038700052490a"

REPORT_NAME = "lei_ren_part1_paper_compliant_current_original_outer_pressure_enclosure.json"
REPORT_PATH = BASE / REPORT_NAME
RUNNER_PATH = "scripts/run_original_outer_pressure_enclosure.py"

CHECKPOINT_DIR = ROOT / "tmp/original_outer_pressure_enclosure"
SOURCE_CHECKPOINT = CHECKPOINT_DIR / "new_source.json"
OUTER_CHECKPOINT = CHECKPOINT_DIR / "new_outer.json"
SIGNED_MACRO_CHECKPOINT = CHECKPOINT_DIR / "new_signed_macro.json"
TAIL_CHECKPOINT = CHECKPOINT_DIR / "new_tail.json"
H_CHECKPOINT = CHECKPOINT_DIR / "new_H.json"
ABSOLUTE_CHECKPOINT = CHECKPOINT_DIR / "new_absolute.json"

KERNEL_PATH = "src/openai_ns_reconstruction/original_outer_pressure_enclosure.py"
KERNEL_SHA256 = "f332c4615f8135bdfdd821f8abdfe98194da1869165648f20b65bce67210ccbe"
KERNEL_BYTES = 30_073
FIXTURE_PATH = "experiments/root_st073/original_outer_pressure_enclosure_fixture.json"
FIXTURE_SHA256 = "35ef9725fafce721b9f8d2bdc57e51b6b22dec38684be45ac567a40c84a78006"
FIXTURE_BYTES = 27_815
SIGNED_MACRO_KERNEL_PATH = "src/openai_ns_reconstruction/original_outer_signed_macro_scale.py"
SIGNED_MACRO_KERNEL_SHA256 = "cef0ec74b65926c28350e5550fd704545df4889a1acabbe363b64fbf6b3b95e2"
SIGNED_MACRO_KERNEL_BYTES = 21_664
SIGNED_MACRO_FIXTURE_PATH = "experiments/root_st073/original_outer_signed_macro_scale_fixture.json"
SIGNED_MACRO_FIXTURE_SHA256 = "ac9e636bbb4f58d0281453d4bb52199946ab96910e4c08cf9cd243baf7e7f8e6"
H_KERNEL_PATH = "src/openai_ns_reconstruction/original_H_m15_Rm_transport.py"
H_KERNEL_SHA256 = "2d6120b2f06926f4df44b78484ff815a32385d3667930af85c48a94ca57ae9ba"
H_KERNEL_BYTES = 15_013
H_FIXTURE_PATH = "experiments/root_st073/original_H_m15_Rm_transport_fixture.json"
H_FIXTURE_SHA256 = "eeb629e63f534e528306e293898a378e170985682c5f4ed5e8c7b785eb638a06"
H_FIXTURE_BYTES = 4_678

# The strict logarithmic threshold and the wider working interval precision
# are separate contracts.  The latter is the only 1600-bit setting here.
RADIUS_BITS = 1024
CONTEXT_PRECISION_BITS = 1600

INSPECTION_REPAIR_INPUT_PINS = {
    "experiments/root_st073/original_outer_pressure_runner_execution_source.txt":
        "a50fce081283981a13b3be487a5b16162d918d8cbac6a8a1e8ae5083c482507e",
}
INSPECTION_FORMAT_REPAIR = {
    "executed_runner_preimage_sha256": "a50fce081283981a13b3be487a5b16162d918d8cbac6a8a1e8ae5083c482507e",
    "initial_unaccepted_report_sha256": "bb6dbfa90925bfc12677358c4eca10d49217b8b1abab57e558decf27343cade0",
    "initial_saved_inspection_passed": False,
    "reason": "Authenticated source selected_Z_interval is a coordinate descriptor without directed flag, not an emitted numeric interval.",
    "repair_scope": "Inspector descriptor classification and saved-only report provenance; mathematical function ASTs unchanged.",
    "mathematical_components_replayed": False,
}

ORDINARY_Z_ORDERS = [0, 1]
Z_INTERVAL = ["37/100", "93/250"]
OPAQUE_ATOMS = [[183935, 184677], [183936, 184678]]
H_ANCHOR_DIGEST = "ec934fbfce5aa21708541d5d65a87a920c3239e599f333370aef1dc650985e17"

SOURCE_NODE_HASHES = {
    "prefix": "2d6efeb13efe91dd95d4f24fa03eb4626c7685d04e5f29f2bdc6532e101994e3",
    "pressure": "dfbba64e1f106fa42c233af28a096f8dfad1e35580163c54cad9eb523d21271e",
    "correlated": "83b8594b5b775520cd38c9da1e6633e800ebf170dfb48b6b12eb34f903514fc0",
}
ADOPTED_ROOTS = [186096, 186099]


class DataOnlyBlocked(RuntimeError):
    """Raised when a root-only production artifact is not available."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def sha256_bytes(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_digest(value: object) -> str:
    raw = json.dumps(value, sort_keys=True, separators=(",", ":"),
                     ensure_ascii=True, allow_nan=False).encode("utf-8")
    return sha256_bytes(raw)


def check_lf_single_eof(raw: bytes, label: str) -> None:
    require(b"\r" not in raw, label + " contains CR bytes")
    require(raw.endswith(b"\n"), label + " must end in one LF")
    require(not raw.endswith(b"\n\n"), label + " has more than one terminal LF")


def authenticated_json(path: Path, *, expected_bytes: int, expected_sha256: str,
                       label: str) -> dict:
    """Hash and authenticate bytes before decoding JSON."""
    raw = path.read_bytes()
    check_lf_single_eof(raw, label)
    require(len(raw) == expected_bytes,
            f"{label} byte count changed: {len(raw)} != {expected_bytes}")
    actual = sha256_bytes(raw)
    require(actual == expected_sha256,
            f"{label} SHA-256 changed: {actual} != {expected_sha256}")
    value = json.loads(raw)
    require(isinstance(value, dict), label + " must decode to an object")
    return value


def authenticated_file(path: Path, *, expected_bytes: int, expected_sha256: str,
                       label: str) -> None:
    """Authenticate a pinned implementation/fixture before any decode."""
    raw = path.read_bytes()
    check_lf_single_eof(raw, label)
    require(len(raw) == expected_bytes,
            f"{label} byte count changed: {len(raw)} != {expected_bytes}")
    actual = sha256_bytes(raw)
    require(actual == expected_sha256,
            f"{label} SHA-256 changed: {actual} != {expected_sha256}")


def source_ast_preflight(path: Path, label: str) -> None:
    raw = path.read_bytes()
    try:
        tree = ast.parse(raw.decode("utf-8"), filename=label)
    except (UnicodeDecodeError, SyntaxError) as exc:
        raise ValueError(label + " AST/source preflight failed") from exc
    require(isinstance(tree, ast.Module), label + " AST is not a module")


def authenticated_inputs() -> tuple[dict, dict]:
    parent = authenticated_json(ROOT / PARENT_PATH, expected_bytes=PARENT_BYTES,
                                expected_sha256=PARENT_SHA256, label=PARENT_PATH)
    native = authenticated_json(ROOT / NATIVE_PATH, expected_bytes=NATIVE_BYTES,
                                expected_sha256=NATIVE_SHA256, label=NATIVE_PATH)
    authenticate_inherited_pins((parent, native))
    return parent, native


def _safe_inherited_relative(name: str) -> Path:
    path = Path(name)
    require(not path.is_absolute(), "absolute inherited path: " + name)
    if ".." in path.parts:
        require(name == "../../src/openai_ns_reconstruction/schedule_pressure.py",
                "unexpected inherited traversal: " + name)
        return (BASE / path).resolve()
    if path.parts and path.parts[0] in {"experiments", "src", "scripts", "tmp"}:
        return (ROOT / path).resolve()
    return (BASE / path.name).resolve()


def authenticate_inherited_pins(reports: tuple[dict, ...]) -> None:
    """Hash every inherited pin without decoding any pinned artifact."""
    for report in reports:
        for field in ("input_hashes", "implementation_hashes"):
            pins = report.get(field)
            require(isinstance(pins, dict), "missing inherited " + field)
            for name, expected in pins.items():
                require(isinstance(name, str) and isinstance(expected, str),
                        "malformed inherited " + field + " pin")
                path = _safe_inherited_relative(name)
                require(path.is_file(), field + " path missing: " + name)
                actual = sha256_file(path)
                require(actual == expected, field + " changed: " + name)


def checkpoint_json(path: Path) -> dict:
    raw = path.read_bytes()
    check_lf_single_eof(raw, str(path))
    value = json.loads(raw)
    require(isinstance(value, dict), str(path) + " must decode to an object")
    return value


def write_immutable_json(path: Path, value: dict) -> str:
    """Create one checkpoint without ever replacing an existing file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    encoded = (json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False)
               + "\n").encode("utf-8")
    check_lf_single_eof(encoded, str(path))
    try:
        with path.open("xb") as stream:
            stream.write(encoded)
    except FileExistsError as exc:
        raise DataOnlyBlocked("immutable checkpoint already exists: " + str(path)) from exc
    return sha256_bytes(encoded)


def _exact_mpf_tuple(value: object, label: str) -> list[int]:
    require(isinstance(value, list) and len(value) == 4,
            label + " must be a four-field exact mpf tuple")
    require(all(type(item) is int for item in value),
            label + " has non-integer exact mpf fields")
    sign, mantissa, exponent, precision = value
    require(sign in (0, 1) and mantissa >= 0 and exponent == int(exponent),
            label + " has malformed sign/mantissa/exponent")
    if mantissa == 0:
        require(value == [0, 0, 0, 0], label + " has a noncanonical zero")
    else:
        require(mantissa & 1, label + " has a non-normalized mantissa")
    require(precision == mantissa.bit_length(),
            label + " precision does not match the exact mantissa")
    return value


def raw_iv(value: object, label: str) -> dict:
    require(isinstance(value, dict), label + " must be a directed interval object")
    require("lower_exact_mpf_tuple" in value and "upper_exact_mpf_tuple" in value,
            label + " is missing exact endpoints")
    _exact_mpf_tuple(value["lower_exact_mpf_tuple"], label + ".lower")
    _exact_mpf_tuple(value["upper_exact_mpf_tuple"], label + ".upper")
    require(value.get("directed") is True, label + " is not directed")
    return copy.deepcopy(value)


def canonical_exact_tuple_records(value: object, path: str = "") -> int:
    """Validate exact tuples, distinguishing the pinned source Z descriptor."""
    if isinstance(value, dict):
        count = 0
        if "lower_exact_mpf_tuple" in value or "upper_exact_mpf_tuple" in value:
            if path == ".source_consumer_checkpoint.H_source_basis_snapshot.selected_Z_interval":
                require(set(value) == {"lower", "upper", "width", "lower_exact_mpf_tuple",
                                       "upper_exact_mpf_tuple"},
                        "source Z descriptor shape changed")
                _exact_mpf_tuple(value["lower_exact_mpf_tuple"], "source Z descriptor lower")
                _exact_mpf_tuple(value["upper_exact_mpf_tuple"], "source Z descriptor upper")
            else:
                raw_iv(value, "canonical exact tuple record" + path)
            count = 1
        return count + sum(canonical_exact_tuple_records(item, path + "." + str(key))
                           for key, item in value.items())
    if isinstance(value, list):
        return sum(canonical_exact_tuple_records(item, path + "[" + str(i) + "]")
                   for i, item in enumerate(value))
    return 0


def _raw_to_interval(c: object, value: dict):
    """Convert one authenticated raw interval only inside root production."""
    import mpmath as mp

    lower = tuple(value["lower_exact_mpf_tuple"])
    upper = tuple(value["upper_exact_mpf_tuple"])
    return c.mpf([mp.make_mpf(lower), mp.make_mpf(upper)])


def _raw_rows_to_intervals(c: object, rows: list[dict]) -> list[object]:
    return [_raw_to_interval(c, row) for row in rows]


def _serialize_math(value: object) -> object:
    """Serialize mpmath values without rounding or importing math at module load."""
    rows = getattr(value, "_mpi_", None)
    if rows is not None:
        if len(rows) == 4 and all(type(part) is int for part in rows):
            rows = (rows, rows)
        elif len(rows) == 1:
            rows = (rows[0], rows[0])
        require(len(rows) == 2, "malformed produced interval")
        return {
            "lower_exact_mpf_tuple": list(tuple(rows[0])),
            "upper_exact_mpf_tuple": list(tuple(rows[1])),
            "directed": True,
        }
    point = getattr(value, "_mpf_", None)
    if point is not None:
        return list(tuple(point))
    if isinstance(value, dict):
        return {key: _serialize_math(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_serialize_math(item) for item in value]
    return value


def _production_checkpoint(value: dict, stage: str, **links: str) -> dict:
    """Attach a self-authenticating root-only checkpoint before validation."""
    result = copy.deepcopy(value)
    result["checkpoint"] = {
        "schema": "original_outer_pressure_enclosure_production_checkpoint_v1",
        "stage": stage,
        "written_once": True,
    }
    result.update(links)
    result["checkpoint_sha256"] = canonical_digest(result)
    return result


def _checkpoint_with_seconds(value: dict, seconds: float) -> dict:
    result = copy.deepcopy(value)
    result["execution_seconds"] = max(0.0, round(float(seconds), 6))
    result["checkpoint_sha256"] = canonical_digest(
        {key: item for key, item in result.items()
         if key != "checkpoint_sha256"})
    return result


def _checkpoint_sha(value: dict) -> str:
    expected = value.get("checkpoint_sha256")
    require(isinstance(expected, str), "production checkpoint hash missing")
    body = {key: item for key, item in value.items()
            if key != "checkpoint_sha256"}
    require(expected == canonical_digest(body),
            "production checkpoint hash changed")
    return expected


def explicit_rows(rows: object, label: str, field: str = "value") -> list[dict]:
    """Select exactly ordinary Z rows 0 and 1 without a q/operator loop."""
    require(isinstance(rows, list) and len(rows) >= 2,
            label + " must contain explicit ordinary-Z rows 0 and 1")
    first, second = rows[0], rows[1]
    require(isinstance(first, dict) and isinstance(second, dict),
            label + " rows must be objects")
    require(first.get("ordinary_Z_order") == 0 and second.get("ordinary_Z_order") == 1,
            label + " first two rows are not ordinary Z 0 and 1")
    require(field in first and field in second, label + " rows lack " + field)
    return [raw_iv(first[field], label + "[0]." + field),
            raw_iv(second[field], label + "[1]." + field)]


def explicit_log_rows(rows: object, label: str) -> list[dict]:
    require(isinstance(rows, list) and len(rows) >= 2,
            label + " must contain explicit ordinary-Z log rows 0 and 1")
    first, second = rows[0], rows[1]
    require(isinstance(first, dict) and isinstance(second, dict),
            label + " rows must be objects")
    require(first.get("ordinary_Z_order") == 0 and second.get("ordinary_Z_order") == 1,
            label + " first two rows are not ordinary Z 0 and 1")
    require("log_absolute_upper_bound" in first and "log_absolute_upper_bound" in second,
            label + " rows lack log_absolute_upper_bound")
    return [raw_iv(first["log_absolute_upper_bound"], label + "[0].log_absolute_upper_bound"),
            raw_iv(second["log_absolute_upper_bound"], label + "[1].log_absolute_upper_bound")]


def path_value(value: dict, path: str) -> object:
    current: object = value
    for key in path.split("."):
        require(isinstance(current, dict) and key in current,
                "missing saved input path: " + path)
        current = current[key]
    return current


def source_preconditions(parent: dict, native: dict) -> None:
    authenticated_file(ROOT / KERNEL_PATH, expected_bytes=KERNEL_BYTES,
                       expected_sha256=KERNEL_SHA256,
                       label=KERNEL_PATH)
    source_ast_preflight(ROOT / KERNEL_PATH, KERNEL_PATH)
    authenticated_file(ROOT / FIXTURE_PATH, expected_bytes=FIXTURE_BYTES,
                       expected_sha256=FIXTURE_SHA256,
                       label=FIXTURE_PATH)
    authenticated_file(ROOT / SIGNED_MACRO_KERNEL_PATH,
                       expected_bytes=SIGNED_MACRO_KERNEL_BYTES,
                       expected_sha256=SIGNED_MACRO_KERNEL_SHA256,
                       label=SIGNED_MACRO_KERNEL_PATH)
    source_ast_preflight(ROOT / SIGNED_MACRO_KERNEL_PATH, SIGNED_MACRO_KERNEL_PATH)
    authenticated_file(ROOT / SIGNED_MACRO_FIXTURE_PATH,
                       expected_bytes=14_908,
                       expected_sha256=SIGNED_MACRO_FIXTURE_SHA256,
                       label=SIGNED_MACRO_FIXTURE_PATH)
    require(isinstance(H_KERNEL_BYTES, int) and H_KERNEL_BYTES > 0,
            "reviewed H kernel byte count is missing")
    authenticated_file(ROOT / H_KERNEL_PATH,
                       expected_bytes=H_KERNEL_BYTES,
                       expected_sha256=H_KERNEL_SHA256,
                       label=H_KERNEL_PATH)
    source_ast_preflight(ROOT / H_KERNEL_PATH, H_KERNEL_PATH)
    if H_FIXTURE_SHA256 != "REVIEW_PENDING":
        require(isinstance(H_FIXTURE_BYTES, int) and H_FIXTURE_BYTES > 0,
                "reviewed H fixture byte count is missing")
        authenticated_file(ROOT / H_FIXTURE_PATH,
                           expected_bytes=H_FIXTURE_BYTES,
                           expected_sha256=H_FIXTURE_SHA256,
                           label=H_FIXTURE_PATH)
    prefix = path_value(parent, "actual_pressure_prefix_contraction")
    source = path_value(parent, "prefix_source_component")
    require(prefix.get("passed") is True and source.get("passed") is True,
            "accepted pressure prefix/source component is not passed")
    require(prefix.get("ordinary_Z_orders") == ORDINARY_Z_ORDERS,
            "pressure prefix ordinary-Z scope changed")
    require(prefix.get("normalization") == "own_Rm^2",
            "pressure prefix normalization changed")
    require(prefix.get("opaque_atoms") == [
        {"ordinary_Z_order": 0, "handles": OPAQUE_ATOMS[0],
         "reason": "unresolved first-micro lower RQy and outer-micro/macro half-S1"},
        {"ordinary_Z_order": 1, "handles": OPAQUE_ATOMS[1],
         "reason": "unresolved first-micro lower RQy and outer-micro/macro half-S1"},
    ], "pressure opaque atoms changed")
    anchor = source.get("accepted_H_anchor_receipt")
    require(isinstance(anchor, dict), "accepted H anchor receipt missing")
    require(anchor.get("shifted_H_anchor_digest") == H_ANCHOR_DIGEST,
            "accepted H anchor digest changed")
    require(anchor.get("do_not_add_raw_increment_to_H") is True,
            "raw pressure increment may not be added to H")
    require(anchor.get("first_lower_RQy_remains_opaque") is True,
            "first lower RQy opacity was removed")
    patch = source.get("graph_patch_receipt")
    require(isinstance(patch, dict), "source graph patch receipt missing")
    require(patch.get("restoration_base_sha256") == SOURCE_NODE_HASHES["prefix"]
            and patch.get("restoration_full_sha256") == SOURCE_NODE_HASHES["pressure"]
            and patch.get("correlated_full_sha256") == SOURCE_NODE_HASHES["correlated"],
            "source graph namespace hash changed")
    require(source.get("source_contract", {}).get("cumulative_terms_preserved") == 212,
            "cumulative term count changed")
    require(source.get("source_contract", {}).get("signed_core_lower_terms_preserved") == 6,
            "signed core lower term count changed")
    require(source.get("source_contract", {}).get("reserved_n1_own_swirl_registration_count_preserved") == 1,
            "reserved own-swirl count changed")
    require(source.get("source_contract", {}).get("new_reserved_term") is False,
            "new reserved term was introduced")
    require(source.get("parent_source_partition", {}).get("adopted_roots") == ADOPTED_ROOTS,
            "adopted pressure roots changed")

    q_packet = path_value(native, "actual_micro_Q_increment_bounds")
    micro_packet = path_value(native, "actual_micro_background_I1_numeric")
    macro_packet = path_value(native, "actual_macro_background_I1_numeric")
    macro_finite = path_value(native, "actual_macro_finite_N_I1_numeric")
    first_micro_seam = path_value(native, "actual_first_micro_seam")
    reserved_m15 = path_value(native, "actual_reserved_point_m15")
    correlated = path_value(native, "actual_micro_macro_correlated_m15")
    for packet, label in ((q_packet, "native micro Q increment"),
                          (micro_packet, "native micro background"),
                          (macro_packet, "native macro background"),
                          (macro_finite, "native macro finite-N")):
        require(packet.get("passed") is True, label + " is not passed")
        require(packet.get("accepted_ancestors_replayed") is False,
                label + " replayed accepted ancestors")
    q_contract = q_packet.get("source_contract", {})
    require(q_contract.get("same_source_Q_replacement_not_extra_ledger_contribution") is True,
            "micro Q replacement ledger guard missing")
    require(q_contract.get("full_recovered_Q_equals_original_by_same_U_and_mean_function_identity") is True,
            "micro Q same-source identity missing")
    require(q_contract.get("phase_Z_independent_and_total_y_phase_derivative_N") is True,
            "micro Q phase/Z independence guard missing")
    require(q_contract.get("accepted_Q_recovery_controls_consumed") is True,
            "accepted Q recovery control receipt missing")
    graph = q_contract.get("source_graph_handles", {}).get("first_micro", {})
    require(graph.get("original_Q_source_nodes") == [129312, 129316, 129322],
            "first-micro original Q source ports changed")
    require(graph.get("physical_Q_ordinary_Z0_Z2") == [184686, 184690, 184696],
            "first-micro physical Q handles changed")
    require(graph.get("physical_Q_ordinary_Z0_Z2") != graph.get("original_Q_source_nodes"),
            "namespace guard collapsed physical and original Q handles")
    require(q_packet.get("domain", {}).get("Z") == Z_INTERVAL,
            "native physical Z domain changed")
    require(q_packet.get("new_linear_Q_controls", {}).get("physical_Pstar_applied_once") is True,
            "physical Pstar-once guard missing")
    require(q_packet.get("new_linear_Q_controls", {}).get("operator_linearity_exact") is True,
            "Q operator-linearity guard missing")
    require(micro_packet.get("normalized_by_W") is True
            and micro_packet.get("finite_N_error_added_once_as_normalized_radius") is True,
            "micro W normalization/error-once guard changed")
    require(macro_packet.get("normalized_by_Wmacro") is True
            and macro_packet.get("Wmacro_definition") == "Wmacro=10000",
            "macro W normalization changed")
    require(macro_finite.get("physical_error_log_bounds", {}).get("error_center_not_selected") is True,
            "macro error center became a selected value")
    require(reserved_m15.get("passed") is True
            and reserved_m15.get("accepted_ancestors_replayed") is False
            and reserved_m15.get("source_executed") is False,
            "reserved M15 source receipt is not accepted source data")
    require(reserved_m15.get("inner_common_unit_cover_is_not_a_function_value") is True
            and reserved_m15.get("ordinary_Z1_is_physical_derivative_divided_by_fixed_unit") is True,
            "reserved M15 common-unit semantics changed")
    raw_iv(reserved_m15.get("common_m15_unit_log"),
           "native.actual_reserved_point_m15.common_m15_unit_log")
    for order in ORDINARY_Z_ORDERS:
        h_row = first_micro_seam.get("ordinary_Z", {}).get(str(order), {})
        require(h_row.get("common_unit") == "reference*A*Rv^2",
                "first-micro H common unit changed")
        require(h_row.get("caps_not_selected_as_function_values") is True,
                "first-micro H log cap became a selected value")
        raw_iv(h_row.get("absolute_common_unit_log_upper"),
               f"native.actual_first_micro_seam.ordinary_Z[{order}].absolute_common_unit_log_upper")
    geometry = macro_packet.get("geometric_measure", {})
    require(geometry.get("radius_definition") == "R=Ra*exp(2*h_b+ell*m)"
            and geometry.get("length_definition") == "ell=Y-2*h_b; Y=log(100)-log(Ra)",
            "macro radius/length identity changed")
    require(geometry.get("R_upper") == 100 and geometry.get("Wmacro") == 10000,
            "macro endpoint identity changed")
    require(geometry.get("physical_R_squared_length_applied_once") is True
            and geometry.get("R_positive_for_0_lt_m_le_1") is True,
            "macro physical radius/length receipt changed")
    target = macro_finite.get("source_contract", {}).get(
        "source_graph_handles", {}).get("target_from_parent", {})
    require(target.get("passed") is True
            and target.get("physical_fields", {}).get("Q") == [140780, 140784, 140790],
            "macro Q target source namespace changed")
    require(target.get("physical_measure", {}).get("R_squared_length_applied_once") is True,
            "macro target physical measure receipt changed")
    joined = correlated.get("actual_joined_pressure_function", {})
    require(joined.get("reduced_boundary_B", {}).get("outer") == [184671, 184672],
            "outer B source handles changed")
    require(joined.get("reduced_volume_D", {}).get("outer") == [184673, 184674],
            "outer D source handles changed")
    require(joined.get("shared_U_Q_B_seam", {}).get("outer_B_plus_D") == [184675, 184676],
            "outer S1 source handles changed")
    require(joined.get("ordinary_Z_orders") == ORDINARY_Z_ORDERS
            and joined.get("pressure_increment") == [184685, 184686],
            "outer pressure join receipt changed")
    require(first_micro_seam.get("passed") is True
            and first_micro_seam.get("accepted_ancestors_replayed") is False,
            "actual first-micro seam receipt is not accepted")
    require(first_micro_seam.get("ordinary_Z_orders") == ORDINARY_Z_ORDERS,
            "actual first-micro seam ordinary-Z scope changed")
    require(first_micro_seam.get("archived_preRp_endpoint_has_extra_Pstar") is True,
            "archived pre-Rp extra-Pstar receipt changed")
    require(first_micro_seam.get("archived_preRp_m15_is_not_corrected_physical_source") is True,
            "archived pre-Rp source distinction changed")
    require(first_micro_seam.get("no_second_F0_derivative_applied") is True,
            "second F0 derivative guard changed")
    require(first_micro_seam.get("old_six_numeric_terms_and_all_error_caps_unchanged") is True,
            "old six-term/error-cap preservation changed")
    graph = first_micro_seam.get("corrected_physical_endpoint_graph", {})
    require(first_micro_seam.get("corrected_physical_endpoint_graph_digest")
            == "19322a78adef610dd8298a81c7157bce2ecd2743c71ef0311340a8f22f65b235",
            "corrected physical pre-Rp graph digest changed")
    require(graph.get("original_pairs", {}).get("preRp_m15") == [31221, 31220],
            "archived pre-Rp pair changed")
    require(graph.get("corrected_pairs", {}).get("preRp_m15") == [31250, 31251],
            "corrected physical pre-Rp pair changed")
    require(graph.get("old_graph_unchanged") is True,
            "corrected physical graph was treated as current pressure graph")
    for order in ORDINARY_Z_ORDERS:
        rows = first_micro_seam.get("ordinary_Z", {}).get(str(order), {})
        require(len(rows.get("signed_approximation", [])) == 12,
                "first-micro seam signed approximation row count changed")
        require(len(rows.get("error_radius", [])) == 61,
                "first-micro seam error-radius row count changed")
        require(len(rows.get("collar_axis_inertial_model", [])) == 1,
                "collar axis inertial model was duplicated or omitted")
        require(len(rows.get("collar_inertial_error_radius", [])) == 21,
                "collar inertial error-radius row count changed")


def build_source_consumer_packet(parent: dict, native: dict) -> dict:
    """Extract raw saved inputs; this function performs no interval arithmetic."""
    source_preconditions(parent, native)
    micro = path_value(native, "actual_micro_background_I1_numeric")
    q_packet = path_value(native, "actual_micro_Q_increment_bounds")
    macro = path_value(native, "actual_macro_background_I1_numeric")
    macro_finite = path_value(native, "actual_macro_finite_N_I1_numeric")
    first_micro_seam = path_value(native, "actual_first_micro_seam")
    reserved_m15 = path_value(native, "actual_reserved_point_m15")

    micro_rows = micro["signed_source_rows"]
    macro_fields = macro["actual_macro_fields"]
    micro_errors = q_packet["charts"]["first"]["log_caps"]
    macro_errors = macro_finite["physical_error_log_bounds"]

    micro_d = micro["reduced_micro_volume_D"]
    macro_d = macro_finite["D_full"]
    micro_logw = micro["physical_scale"]["log_W"]
    macro_w = macro["geometric_measure"]["Wmacro_record"]
    macro_r_upper = macro["geometric_measure"]["R_bound"]
    require(macro["geometric_measure"].get("R_upper") == 100,
            "macro exact R upper is not 100")
    require(macro["geometric_measure"].get("Wmacro") == 10000,
            "macro exact W is not 10000")

    delta = raw_iv(q_packet["actual_delta"], "delta")
    log_rm = raw_iv(parent["actual_pressure_prefix_log_transport"]["logRm"],
                    "logRm")

    packet = {
        "schema": "original_outer_pressure_enclosure_source_consumer_v1",
        "passed": True,
        "source_only": True,
        "source_executed": False,
        "accepted_ancestors_replayed": False,
        "parent_report": PARENT_PATH,
        "parent_report_sha256": PARENT_SHA256,
        "native_report": NATIVE_PATH,
        "native_report_sha256": NATIVE_SHA256,
        "source_family": {
            "namespace_prefix": 185704,
            "namespace_prefix_sha256": SOURCE_NODE_HASHES["prefix"],
            "pressure_full_node": 186086,
            "pressure_full_sha256": SOURCE_NODE_HASHES["pressure"],
            "correlated_full_node": 186100,
            "correlated_full_sha256": SOURCE_NODE_HASHES["correlated"],
            "adopted_roots": copy.deepcopy(ADOPTED_ROOTS),
            "namespace_guard": "original_Q_source_nodes_are_not_physical_Q_graph_handles",
        },
        "ordinary_Z_orders": copy.deepcopy(ORDINARY_Z_ORDERS),
        "z": copy.deepcopy(Z_INTERVAL),
        "z_source": "native.actual_micro_Q_increment_bounds.domain.Z",
        "domain": {
            "Z": copy.deepcopy(Z_INTERVAL),
            "Z_source": "native.actual_micro_Q_increment_bounds.domain.Z",
            "physical_U_Q_units": True,
            "Pstar_applied_once": True,
            "nonreset_memory": True,
            "phase_Z_independent": True,
        },
        "opaque_atoms": copy.deepcopy(OPAQUE_ATOMS),
        "accepted_H_anchor_receipt": copy.deepcopy(
            parent["prefix_source_component"]["accepted_H_anchor_receipt"]),
        "H_anchor_scope": {
            "digest": H_ANCHOR_DIGEST,
            "role": "moment_inlet_M15_anchor",
            "not_axis_preheat_H_equals_1": True,
            "numeric_H_installed_here": False,
            "raw_pressure_increment_added": False,
        },
        "outer_graph_source_receipt": {
            "source_path": "native.actual_micro_macro_correlated_m15.actual_joined_pressure_function",
            "outer_B": [184671, 184672],
            "outer_D": [184673, 184674],
            "outer_S1": [184675, 184676],
            "outer_half_S1": [184677, 184678],
            "pressure_increment": [184685, 184686],
            "outer_graph_namespace": "accepted_correlated_joined_pressure_function",
            "current_pressure_graph_namespace_not_used": True,
        },
        "macro_Q_target_source_receipt": {
            "source_path": "native.actual_macro_finite_N_I1_numeric.source_contract.source_graph_handles.target_from_parent",
            "original_macro_Q_source_nodes": [140780, 140784, 140790],
            "original_micro_Q_source_nodes": [129312, 129316, 129322],
            "same_source_Q_replacement_not_extra_ledger_contribution": True,
            "target_namespace_separate_from_outer_graph": True,
        },
        "first_micro_seam_source_receipt": {
            "source_path": "native.actual_first_micro_seam",
            "ordinary_Z_orders": copy.deepcopy(ORDINARY_Z_ORDERS),
            "signed_approximation_rows_per_Z": 12,
            "error_radius_rows_per_Z": 61,
            "collar_axis_inertial_model_rows_per_Z": 1,
            "collar_inertial_error_radius_rows_per_Z": 21,
            "archived_preRp_endpoint_has_extra_Pstar": True,
            "archived_preRp_m15_is_not_corrected_physical_source": True,
            "no_second_F0_derivative_applied": True,
            "old_six_numeric_terms_and_all_error_caps_unchanged": True,
            "structural_channels_total": 6,
            "signed_contributors": 5,
            "core_swirl_coupling_unsigned_error_in_error_radius": True,
            "no_sixth_signed_value": True,
            "independent_target_links_closed": False,
            "archived_preRp_pair": [31221, 31220],
            "corrected_physical_preRp_pair": [31250, 31251],
            "corrected_physical_graph_sha256":
                "19322a78adef610dd8298a81c7157bce2ecd2743c71ef0311340a8f22f65b235",
            "current_pressure_graph_namespace_not_used": True,
        },
        "endpoints": {
            "micro_lower": {
                "U_bg": explicit_rows(micro_rows["U"], "micro_lower.U_bg"),
                "Q_bg": explicit_rows(micro_rows["Q"], "micro_lower.Q_bg"),
                "delta_U_logs": explicit_log_rows(
                    micro_errors["physical_delta_U"], "micro_lower.delta_U_logs"),
                "delta_Q_logs": explicit_log_rows(
                    micro_errors["physical_delta_Q"], "micro_lower.delta_Q_logs"),
                "logR_upper": raw_iv(q_packet["physical_R_upper_log"],
                                      "micro_lower.logR_upper"),
                "source_paths": {
                    "U_bg": "native.actual_micro_background_I1_numeric.signed_source_rows.U",
                    "Q_bg": "native.actual_micro_background_I1_numeric.signed_source_rows.Q",
                    "delta_U_logs": "native.actual_micro_Q_increment_bounds.charts.first.log_caps.physical_delta_U",
                    "delta_Q_logs": "native.actual_micro_Q_increment_bounds.charts.first.log_caps.physical_delta_Q",
                    "logR_upper": "native.actual_micro_Q_increment_bounds.physical_R_upper_log",
                },
            },
            "macro_upper": {
                "U_bg": explicit_rows(macro_fields["U_b"], "macro_upper.U_bg"),
                "Q_bg": explicit_rows(macro_fields["Q_b"], "macro_upper.Q_bg"),
                "delta_U_logs": explicit_log_rows(
                    macro_errors["delta_U"], "macro_upper.delta_U_logs"),
                "delta_Q_logs": explicit_log_rows(
                    macro_errors["delta_Q"], "macro_upper.delta_Q_logs"),
                "logR_upper": {
                    "raw_interval": None,
                    "radius_upper_bound": raw_iv(macro_r_upper, "macro_upper.R_upper"),
                    "radius_upper_integer": 100,
                    "radius_upper_integer_is_exact_geometry": True,
                    "formula": "ln(100)",
                    "materialize_in_root_kernel": True,
                },
                "source_paths": {
                    "U_bg": "native.actual_macro_background_I1_numeric.actual_macro_fields.U_b",
                    "Q_bg": "native.actual_macro_background_I1_numeric.actual_macro_fields.Q_b",
                    "delta_U_logs": "native.actual_macro_finite_N_I1_numeric.physical_error_log_bounds.delta_U",
                    "delta_Q_logs": "native.actual_macro_finite_N_I1_numeric.physical_error_log_bounds.delta_Q",
                    "R_upper": "native.actual_macro_background_I1_numeric.geometric_measure.R_bound",
                },
            },
        },
        "volume_terms": {
            "micro": {
                "D_rows": explicit_rows(micro_d["normalized_values"],
                                          "micro.volume_terms.D_rows"),
                "normalization": "own_W",
                "logW": raw_iv(micro_logw, "micro.volume_terms.logW"),
                "source_path": "native.actual_micro_background_I1_numeric.reduced_micro_volume_D.normalized_values",
                "logW_source_path": "native.actual_micro_background_I1_numeric.physical_scale.log_W",
            },
            "macro": {
                "D_rows": explicit_rows(macro_d["D_rows"],
                                          "macro.volume_terms.D_rows"),
                "normalization": "own_Wmacro",
                "logW": {
                    "raw_interval": None,
                    "weight_exact": raw_iv(macro_w, "macro.volume_terms.Wmacro"),
                    "weight_integer": 10000,
                    "formula": "ln(10000)",
                    "materialize_in_root_kernel": True,
                },
                "source_path": "native.actual_macro_finite_N_I1_numeric.D_full.D_rows",
                "weight_source_path": "native.actual_macro_background_I1_numeric.geometric_measure.Wmacro_record",
            },
        },
        "signed_macro_scale_contract": {
            "module_path": SIGNED_MACRO_KERNEL_PATH,
            "function": "assemble",
            "fixture_path": SIGNED_MACRO_FIXTURE_PATH,
            "reviewed_module_sha256": SIGNED_MACRO_KERNEL_SHA256,
            "reviewed_fixture_sha256": SIGNED_MACRO_FIXTURE_SHA256,
            "arguments": [
                "c", "G_micro", "G_macro", "D_micro", "D_macro",
                "logR_micro_upper", "logW_micro", "logRm",
            ],
            "macro_endpoint_identity": True,
            "R_macro_m1_integer": 100,
            "R_macro_squared_integer": 10000,
            "Wmacro_integer": 10000,
            "radius_formula": "R(m=1)=R_a*exp(2*h_b+ell)=100",
            "length_formula": "ell=Y-2*h_b; Y=log(100)-log(R_a)",
            "radius_formula_is_not_cap": True,
            "B_outer_over_Wmacro":
                "G_macro - bounded_micro_weighted_G",
            "S1_outer_over_Wmacro":
                "B_outer_over_Wmacro + D_macro + bounded_D_micro_over_Wmacro",
            "G_jets_source": "new_outer.metadata.G_jets",
            "recompute_G": False,
            "coarse_outer_bound_preserved": True,
            "signed_B_source_handle": 184671,
            "signed_S1_source_handle": 184675,
            "signed_macro_checkpoint_stage": "after_additional_component_before_validation",
            "macro_radius_source_receipt": {
                "source_path": "native.actual_macro_background_I1_numeric.geometric_measure",
                "radius_definition": macro["geometric_measure"]["radius_definition"],
                "length_definition": macro["geometric_measure"]["length_definition"],
                "R_macro_m1_integer": 100,
                "R_macro_squared_integer": 10000,
                "Wmacro_integer": 10000,
                "identity_is_not_cap": True,
            },
        },
        "H_transport_contract": {
            "module_path": H_KERNEL_PATH,
            "function": "transport",
            "reviewed_module_sha256": H_KERNEL_SHA256,
            "fixture_path": H_FIXTURE_PATH,
            "reviewed_fixture_sha256": H_FIXTURE_SHA256,
            "arguments": ["c", "Hlogbounds", "commonunitlog", "logRm"],
            "Hlogbounds_source": [
                "native.actual_first_micro_seam.ordinary_Z['0'].absolute_common_unit_log_upper",
                "native.actual_first_micro_seam.ordinary_Z['1'].absolute_common_unit_log_upper",
            ],
            "commonunitlog_source": "native.actual_reserved_point_m15.common_m15_unit_log",
            "logRm_source": "parent.actual_pressure_prefix_log_transport.logRm",
            "ordinary_Z_saved_tube_only": True,
            "common_unit": "reference*A*Rv^2",
            "common_m15_unit_is_U_scale": True,
            "scale_Z_independent": True,
            "radius_bits": RADIUS_BITS,
            "context_precision_bits": CONTEXT_PRECISION_BITS,
            "numeric_H_installed_here": False,
            "threshold_failure_is_terminal": True,
            "no_fallback_or_replay": True,
        },
        "H_source_basis_snapshot": {
            "source_path": "native.actual_first_micro_seam",
            "source_log_basis_names": copy.deepcopy(first_micro_seam["source_log_basis_names"]),
            "source_log_bases": copy.deepcopy(first_micro_seam["source_log_bases"]),
            "ordinary_Z": copy.deepcopy(first_micro_seam["ordinary_Z"]),
            "selected_Z_interval": copy.deepcopy(first_micro_seam["selected_Z_interval"]),
            "decompositions_not_added_again": True,
            "source_rows_evaluated_again": False,
        },
        "H_transport_inputs": {
            "Hlogbounds": [
                raw_iv(first_micro_seam["ordinary_Z"]["0"][
                    "absolute_common_unit_log_upper"], "Hlogbounds[0]"),
                raw_iv(first_micro_seam["ordinary_Z"]["1"][
                    "absolute_common_unit_log_upper"], "Hlogbounds[1]"),
            ],
            "commonunitlog": raw_iv(reserved_m15["common_m15_unit_log"],
                                     "commonunitlog"),
            "source_only": True,
        },
        "delta": delta,
        "delta_source_path": "native.actual_micro_Q_increment_bounds.actual_delta",
        "logRm": log_rm,
        "logRm_source_path": "parent.actual_pressure_prefix_log_transport.logRm",
        "formula": {
            "outer_B": "[R^2*G]_macro - [R^2*G]_micro",
            "G": "Q/L + Q^2/2 - 2*Z^2*U^2/L^2",
            "endpoint_pair_scope": True,
            "transport_scope": "D only",
            "half_factor": "outer S1 is divided by 2 exactly once in the linear tail",
        },
        "kernel_contract": {
            "module_path": KERNEL_PATH,
            "function": "enclose",
            "arguments": ["c", "endpoints", "volume_terms", "z", "delta", "logRm", "radius_bits1024"],
            "radius_bits": RADIUS_BITS,
            "context_precision_bits": CONTEXT_PRECISION_BITS,
            "returns": ["outerS1_over_Rm2_raw_iv_pair", "rawB_optional", "rawD_optional", "metadata"],
            "reviewed_sha256": KERNEL_SHA256,
            "fixture_path": FIXTURE_PATH,
            "fixture_sha256": FIXTURE_SHA256,
        },
        "contract": {
            "source_preflight_only": True,
            "new_kernel_not_called": True,
            "accepted_producers_checkers_fixtures_not_called": True,
            "ordinary_Z_rows_explicit": True,
            "raw_intervals_preserved": True,
            "opaque_lower_RQy_retained": True,
            "no_raw_pressure_increment_added_to_H": True,
            "absolute_M15_is_H_anchor_plus_shifted_tail": True,
            "full_pressure_numeric": False,
            "absolute_M15_numeric": False,
            "terminal_match": False,
            "global_closure": False,
            "NS_claim": False,
        },
    }
    return packet


def checkpoint_payload_digest(value: dict) -> str:
    body = {key: item for key, item in value.items() if key != "checkpoint"}
    return canonical_digest(body)


def attach_checkpoint(value: dict, stage: str) -> dict:
    value = copy.deepcopy(value)
    value["checkpoint"] = {
        "schema": "original_outer_pressure_enclosure_checkpoint_v1",
        "stage": stage,
        "written_once": True,
        "payload_sha256": checkpoint_payload_digest(value),
    }
    return value


def validate_source_checkpoint(source: dict, parent: dict, native: dict) -> None:
    require(source.get("schema") == "original_outer_pressure_enclosure_source_consumer_v1",
            "source checkpoint schema changed")
    require(source.get("checkpoint", {}).get("stage") == "source_consumer",
            "source checkpoint stage changed")
    require(source["checkpoint"].get("payload_sha256") == checkpoint_payload_digest(source),
            "source checkpoint payload digest changed")
    require(source.get("parent_report_sha256") == PARENT_SHA256
            and source.get("native_report_sha256") == NATIVE_SHA256,
            "source checkpoint input hashes changed")
    source_preconditions(parent, native)
    require(source.get("opaque_atoms") == OPAQUE_ATOMS, "source opaque atoms changed")
    require(source.get("accepted_H_anchor_receipt", {}).get("shifted_H_anchor_digest") == H_ANCHOR_DIGEST,
            "source H anchor changed")
    require(source.get("source_only") is True and source.get("source_executed") is False,
            "source checkpoint execution status changed")
    expected = build_source_consumer_packet(parent, native)
    actual_payload = {key: item for key, item in source.items()
                      if key != "checkpoint"}
    require(actual_payload == expected,
            "source checkpoint inputs or source receipts differ from authenticated inputs")


def _numeric_inputs(c: object, source: dict) -> tuple[dict, dict, object, object, object, object]:
    """Materialize only authenticated source records for the new kernels."""
    endpoints = {}
    for endpoint_name in ("micro_lower", "macro_upper"):
        packet = source["endpoints"][endpoint_name]
        endpoint = {
            "U_bg": _raw_rows_to_intervals(c, packet["U_bg"]),
            "Q_bg": _raw_rows_to_intervals(c, packet["Q_bg"]),
            "delta_U_logs": _raw_rows_to_intervals(c, packet["delta_U_logs"]),
            "delta_Q_logs": _raw_rows_to_intervals(c, packet["delta_Q_logs"]),
        }
        if endpoint_name == "micro_lower":
            endpoint["logR_upper"] = _raw_to_interval(c, packet["logR_upper"])
        else:
            endpoint["logR_upper"] = c.ln(c.mpf(100))
        endpoints[endpoint_name] = endpoint
    volume_terms = {
        "micro": {
            "D_rows": _raw_rows_to_intervals(
                c, source["volume_terms"]["micro"]["D_rows"]),
            "logW": _raw_to_interval(c, source["volume_terms"]["micro"]["logW"]),
        },
        "macro": {
            "D_rows": _raw_rows_to_intervals(
                c, source["volume_terms"]["macro"]["D_rows"]),
            "logW": c.ln(c.mpf(10_000)),
        },
    }
    z = c.mpf([c.mpf(37) / 100, c.mpf(93) / 250])
    delta = _raw_to_interval(c, source["delta"])
    log_rm = _raw_to_interval(c, source["logRm"])
    return endpoints, volume_terms, z, delta, log_rm, c


def _load_new_kernel(path: Path, module_name: str):
    """Load one reviewed new kernel after source preflight has completed."""
    import importlib.util

    spec = importlib.util.spec_from_file_location(module_name, path)
    require(spec is not None and spec.loader is not None,
            "cannot load reviewed new kernel: " + str(path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _outer_kernel_checkpoint(c: object, source: dict, source_sha: str) -> dict:
    module = _load_new_kernel(ROOT / KERNEL_PATH, "st073_new_outer_pressure_kernel")
    endpoints, volume_terms, z, delta, log_rm, _ = _numeric_inputs(c, source)
    result = module.enclose(c, endpoints, volume_terms, z, delta, log_rm,
                            radius_bits=RADIUS_BITS)
    payload = _serialize_math(result)
    require(isinstance(payload, dict), "outer kernel result is not an object")
    metadata = payload.get("metadata")
    require(isinstance(metadata, dict), "outer kernel metadata is missing")
    endpoints_meta = metadata.get("endpoints")
    require(isinstance(endpoints_meta, dict), "outer kernel endpoint metadata is missing")
    metadata["G_jets"] = {
        "micro_lower": copy.deepcopy(
            endpoints_meta["micro_lower"]["G_signed_rows"]),
        "macro_upper": copy.deepcopy(
            endpoints_meta["macro_upper"]["G_signed_rows"]),
        "source": "new_outer.metadata.endpoints[*].G_signed_rows",
    }
    payload.update({
        "schema": "original_outer_pressure_enclosure_kernel_v1",
        "source_executed": True,
        "accepted_ancestors_replayed": False,
        "source_checkpoint_sha256": source_sha,
        "context_precision_bits": CONTEXT_PRECISION_BITS,
        "radius_bits": RADIUS_BITS,
    })
    return _production_checkpoint(payload, "outer_kernel")


def _signed_macro_checkpoint(c: object, source: dict, outer: dict,
                             source_sha: str, outer_sha: str) -> dict:
    module = _load_new_kernel(ROOT / SIGNED_MACRO_KERNEL_PATH,
                              "st073_new_signed_macro_scale")
    metadata = outer.get("metadata", {})
    jets = metadata.get("G_jets")
    require(isinstance(jets, dict), "saved outer G jets are missing")
    g_micro = _raw_rows_to_intervals(c, jets["micro_lower"])
    g_macro = _raw_rows_to_intervals(c, jets["macro_upper"])
    d_micro = _raw_rows_to_intervals(
        c, source["volume_terms"]["micro"]["D_rows"])
    d_macro = _raw_rows_to_intervals(
        c, source["volume_terms"]["macro"]["D_rows"])
    log_r_micro = _raw_to_interval(
        c, source["endpoints"]["micro_lower"]["logR_upper"])
    log_w_micro = _raw_to_interval(
        c, source["volume_terms"]["micro"]["logW"])
    log_rm = _raw_to_interval(c, source["logRm"])
    result = module.assemble(
        c, g_micro, g_macro, d_micro, d_macro, log_r_micro, log_w_micro,
        log_rm, radius_bits=RADIUS_BITS, macro_endpoint_identity=True,
        geometry_Z_independent=True)
    payload = _serialize_math(result)
    require(isinstance(payload, dict), "signed macro result is not an object")
    payload["schema"] = "original_outer_signed_macro_scale_checkpoint_v1"
    payload["source_executed"] = True
    payload["accepted_ancestors_replayed"] = False
    payload["source_checkpoint_sha256"] = source_sha
    payload["outer_checkpoint_sha256"] = outer_sha
    payload["context_precision_bits"] = CONTEXT_PRECISION_BITS
    payload["radius_bits"] = RADIUS_BITS
    payload.setdefault("metadata", {})["context_precision_bits"] = CONTEXT_PRECISION_BITS
    return _production_checkpoint(payload, "signed_macro_scale")


def _tail_checkpoint(c: object, parent: dict, signed: dict,
                     outer_sha: str, signed_sha: str) -> dict:
    raw_result = _result_field(signed, ("rawIVresult",),
                               "signed macro raw S1/Rm^2")
    signed_rows = _raw_rows_to_intervals(c, raw_result)
    old_rows = [raw_iv(row, "parent.actual_pressure_prefix_contraction.new_partial_rows[" + str(order) + "]")
                for order, row in enumerate(
                    parent["actual_pressure_prefix_contraction"]["new_partial_rows"])]
    require(len(old_rows) == 2,
            "parent.actual_pressure_prefix_contraction.new_partial_rows must have two rows")
    old_intervals = _raw_rows_to_intervals(c, old_rows)
    shifted = [old + row / 2 for old, row in zip(old_intervals, signed_rows)]
    payload = {
        "schema": "original_outer_pressure_enclosure_linear_tail_v2",
        "passed": True,
        "source_executed": True,
        "accepted_ancestors_replayed": False,
        "old_parent_new_partial_rows": copy.deepcopy(
            parent["actual_pressure_prefix_contraction"]["new_partial_rows"]),
        "shifted_tail_Rm2": [_serialize_math(row) for row in shifted],
        "outerS1_factor": {"numerator": 1, "denominator": 2},
        "outer_checkpoint_sha256": outer_sha,
        "signed_macro_checkpoint_sha256": signed_sha,
        "signed_macro_scale_used": True,
        "coarse_outer_bound_only": True,
        "context_precision_bits": CONTEXT_PRECISION_BITS,
        "metadata": {
            "schema": "original_outer_pressure_enclosure_linear_tail_v2",
            "half_factor_applied_once": True,
            "H_half_factor_applied": False,
        },
    }
    return _production_checkpoint(payload, "linear_tail")


def _h_checkpoint(c: object, source: dict, source_sha: str) -> dict:
    module = _load_new_kernel(ROOT / H_KERNEL_PATH, "st073_new_H_m15_Rm_transport")
    inputs = source["H_transport_inputs"]
    h_logs = _raw_rows_to_intervals(c, inputs["Hlogbounds"])
    common_log = _raw_to_interval(c, inputs["commonunitlog"])
    log_rm = _raw_to_interval(c, source["logRm"])
    result = module.transport(c, h_logs, common_log, log_rm,
                              radius_bits=RADIUS_BITS,
                              scale_Z_independent=True)
    payload = _serialize_math(result)
    require(isinstance(payload, dict), "H transport result is not an object")
    payload.update({
        "schema": "original_H_m15_Rm_transport_checkpoint_v1",
        "source_executed": True,
        "accepted_ancestors_replayed": False,
        "source_checkpoint_sha256": source_sha,
        "context_precision_bits": CONTEXT_PRECISION_BITS,
        "radius_bits": RADIUS_BITS,
    })
    return _production_checkpoint(payload, "H_transport")


def _absolute_checkpoint(c: object, h_packet: dict, tail: dict,
                         h_sha: str, tail_sha: str) -> dict:
    h_rows = _raw_rows_to_intervals(
        c, _result_field(h_packet, ("rawIVresult",), "H raw H/Rm^2"))
    tail_rows = _raw_rows_to_intervals(c, tail["shifted_tail_Rm2"])
    absolute = [h + row for h, row in zip(h_rows, tail_rows)]
    payload = {
        "schema": "original_Rm_absolute_m15_checkpoint_v1",
        "passed": True,
        "source_executed": True,
        "accepted_ancestors_replayed": False,
        "new_absolute_rows": [_serialize_math(row) for row in absolute],
        "H_checkpoint_sha256": h_sha,
        "tail_checkpoint_sha256": tail_sha,
        "H_half_factor_applied": False,
        "outer_half_factor_already_in_tail": True,
        "context_precision_bits": CONTEXT_PRECISION_BITS,
        "metadata": {
            "schema": "original_Rm_absolute_m15_checkpoint_v1",
            "identity": "H_m15_inlet/Rm^2 + shifted_tail_Rm^2",
            "H_m15_inlet": "M15(sc/2)+firstlowerRQy",
            "numeric_H_source_bound": True,
            "terminal_target": "Rp (not closed); current_radius=Rm (intermediate)",
            "global_closure": False,
            "NS_claim": False,
        },
    }
    return _production_checkpoint(payload, "absolute_M15")


def validate_raw_outer_packet(outer: dict) -> None:
    require(outer.get("schema") == "original_outer_pressure_enclosure_kernel_v1",
            "new outer checkpoint schema changed")
    require(outer.get("source_executed") is True
            and outer.get("accepted_ancestors_replayed") is False,
            "new outer execution provenance changed")
    result = outer.get("result", outer)
    key = next((candidate for candidate in
                ("outerS1_over_Rm2", "outerS1_norm", "outerS1_Rm2",
                 "rawIVresult") if candidate in result), None)
    require(key is not None, "new outer result lacks outer S1/Rm^2")
    value = result[key]
    require(isinstance(value, list) and len(value) == 2,
            "new outer S1/Rm^2 must have explicit ordinary-Z rows 0 and 1")
    raw_iv(value[0], "new outer S1/Rm^2[0]")
    raw_iv(value[1], "new outer S1/Rm^2[1]")
    metadata = outer.get("metadata", result.get("metadata"))
    if not isinstance(metadata, dict) and isinstance(result, dict):
        metadata = result.get("auditablemetadata")
    require(isinstance(metadata, dict), "new outer metadata missing")
    require(metadata.get("radius_bits") == RADIUS_BITS,
            "new outer radius threshold precision changed")
    require(metadata.get("normalized_by") == "Rm^2"
            or metadata.get("normalization") == "own_Rm^2"
            or metadata.get("rawIVresult_is_outerS1norm") is True,
            "new outer normalization changed")
    jets = metadata.get("G_jets")
    require(isinstance(jets, dict)
            and jets.get("source") ==
            "new_outer.metadata.endpoints[*].G_signed_rows",
            "new outer saved G-jet handoff missing")
    _validate_two_row_result(jets.get("micro_lower"),
                             "new outer G_jets.micro_lower")
    _validate_two_row_result(jets.get("macro_upper"),
                             "new outer G_jets.macro_upper")
    require(_checkpoint_sha(outer) == outer.get("checkpoint_sha256"),
            "new outer checkpoint hash changed")


def _result_field(packet: dict, names: tuple[str, ...], label: str) -> object:
    result = packet.get("result") if isinstance(packet.get("result"), dict) else packet
    for name in names:
        if name in packet:
            return packet[name]
        if isinstance(result, dict) and name in result:
            return result[name]
    raise ValueError(label + " is missing")


def _validate_two_row_result(value: object, label: str) -> None:
    require(isinstance(value, list) and len(value) == 2,
            label + " must have explicit ordinary-Z rows 0 and 1")
    raw_iv(value[0], label + "[0]")
    raw_iv(value[1], label + "[1]")


def validate_signed_macro_packet(signed: dict, outer: dict, source: dict) -> None:
    """Validate the separately checkpointed signed macro-scale component.

    The component is deliberately checked as a consumer of the first outer
    kernel's saved G jets.  This validator never recomputes G or any scale.
    The worker's final schema may wrap the result, so both top-level and
    ``result`` fields are accepted while the raw interval names stay exact.
    """
    schema = signed.get("schema")
    require(schema in {
        "original_outer_signed_macro_scale_v1",
        "original_outer_signed_macro_scale_checkpoint_v1",
    }, "new signed macro-scale checkpoint schema changed")
    require(signed.get("source_executed") is True
            and signed.get("accepted_ancestors_replayed") is False,
            "new signed macro-scale execution provenance changed")
    _validate_two_row_result(
        _result_field(signed, ("rawS1_macro_scale",),
                      "signed macro-scale raw S1"),
        "rawS1_macro_scale")
    _validate_two_row_result(
        _result_field(signed, ("rawB_macro_scale",),
                      "signed macro-scale raw B"),
        "rawB_macro_scale")
    _validate_two_row_result(
        _result_field(signed, ("rawIVresultS1_Rm2_BOUND", "rawIVresult"),
                      "signed macro-scale S1/Rm^2 bound"),
        "rawIVresultS1_Rm2_BOUND")
    metadata = signed.get("metadata")
    if not isinstance(metadata, dict):
        result = signed.get("result")
        metadata = result.get("metadata") if isinstance(result, dict) else None
    require(isinstance(metadata, dict), "new signed macro-scale metadata missing")
    require(metadata.get("macro_endpoint_identity") is True,
            "signed macro-scale macro endpoint identity missing")
    require(metadata.get("macro_endpoint_source_identity") ==
            "R(m=1)=100 exactly; Wmacro=R^2=10000 exactly",
            "signed macro-scale exact endpoint source identity changed")
    require(metadata.get("R_macro_m1_integer", 100) == 100,
            "signed macro-scale R macro endpoint changed")
    require(metadata.get("Wmacro_integer", 10000) == 10000,
            "signed macro-scale Wmacro changed")
    require(metadata.get("radius_bits") == RADIUS_BITS,
            "signed macro-scale radius threshold precision changed")
    require(metadata.get("macro_endpoint_numeric_cap_not_selected", True) is True,
            "macro endpoint cap became a selected value")
    recomputed = metadata.get("recomputed_G",
                              metadata.get("G_recomputed"))
    if recomputed is None:
        recomputed = not metadata.get("source_contract", {}).get(
            "no_source_rows_recomputed", False)
    require(recomputed is False,
            "signed macro-scale recomputed coarse G jets")
    require(metadata.get("coarse_outer_bound_preserved", True) is True,
            "coarse outer bound was not preserved")
    require(metadata.get("outer_B_macro_scale", {}).get("signed") is True
            and metadata.get("outer_S1_macro_scale", {}).get("signed") is True,
            "signed macro-scale B/S1 rows are not marked signed")
    require(metadata.get("outer_S1_macro_scale", {}).get(
        "pressure_half_factor_applied") is False,
            "signed macro-scale applied the outer half factor")
    outer_hash = (signed.get("outer_checkpoint_sha256")
                  or metadata.get("outer_checkpoint_sha256"))
    if outer_hash is not None:
        require(outer_hash == outer.get("checkpoint_sha256"),
                "signed macro-scale outer checkpoint link changed")
    require(source.get("signed_macro_scale_contract", {}).get(
        "G_jets_source") == "new_outer.metadata.G_jets",
            "source G-jet handoff contract changed")
    require(_checkpoint_sha(signed) == signed.get("checkpoint_sha256"),
            "signed macro checkpoint hash changed")


def validate_h_packet(h_packet: dict, source: dict) -> None:
    require(h_packet.get("schema") ==
            "original_H_m15_Rm_transport_checkpoint_v1",
            "H transport checkpoint schema changed")
    require(h_packet.get("source_executed") is True
            and h_packet.get("accepted_ancestors_replayed") is False,
            "H transport execution provenance changed")
    _validate_two_row_result(
        _result_field(h_packet, ("rawIVresult",), "H raw H/Rm^2"),
        "H raw H/Rm^2")
    metadata = h_packet.get("metadata")
    require(isinstance(metadata, dict), "H transport metadata missing")
    require(metadata.get("schema") == "original_H_m15_Rm_transport_v1"
            and metadata.get("normalization") ==
            "absolute_H_over_common_unit_to_Rm2",
            "H transport normalization changed")
    require(metadata.get("ordinary_Z_orders") == ORDINARY_Z_ORDERS
            and metadata.get("scale_Z_independent") is True
            and metadata.get("common_unit_positive") is True
            and metadata.get("common_unit_definition") ==
            "U=reference*A*Rv^2",
            "H transport scale receipt changed")
    require(metadata.get("radius_bits") == RADIUS_BITS,
            "H transport radius threshold changed")
    rows = metadata.get("rows")
    require(isinstance(rows, list) and len(rows) == 2,
            "H transport row metadata changed")
    for order, row in enumerate(rows):
        require(row.get("ordinary_Z_order") == order
                and row.get("original_signed_H_unknown") is True
                and row.get("strict_upper_below_threshold") is True
                and row.get("emitted_kind") == "BOUND"
                and row.get("emitted_not_value_or_zero") is True,
                "H transport strict bound receipt changed")
        require(row.get("source_atom") ==
                "native.actual_first_micro_seam.ordinary_Z[" + str(order)
                + "].absolute_common_unit_log_upper",
                "H transport source atom changed")
    contract = metadata.get("source_contract", {})
    require(contract.get("source_actual_first_micro_seam") is True
            and contract.get("input_is_absolute_common_unit_log_upper") is True
            and contract.get("input_is_H_over_U") is True
            and contract.get("signed_H_unknown_covered") is True
            and contract.get("source_rows_recomputed") is False
            and contract.get("additional_local_to_common_conversion") is False
            and contract.get("F0_derivative_added") is False
            and contract.get("basis_evaluation") is False
            and contract.get("requadrature") is False
            and contract.get("ratio_materialized") is False
            and contract.get("huge_exponentiation") is False,
            "H transport source contract changed")
    require(source.get("H_transport_contract", {}).get(
        "scale_Z_independent") is True,
        "source H scale contract changed")
    require(_checkpoint_sha(h_packet) == h_packet.get("checkpoint_sha256"),
            "H transport checkpoint hash changed")


def validate_absolute_packet(absolute: dict, h_packet: dict,
                             tail: dict) -> None:
    require(absolute.get("schema") == "original_Rm_absolute_m15_checkpoint_v1",
            "absolute M15 checkpoint schema changed")
    require(absolute.get("source_executed") is True
            and absolute.get("accepted_ancestors_replayed") is False,
            "absolute M15 execution provenance changed")
    _validate_two_row_result(absolute.get("new_absolute_rows"),
                             "new absolute M15 rows")
    require(absolute.get("H_half_factor_applied") is False
            and absolute.get("outer_half_factor_already_in_tail") is True,
            "absolute M15 half-factor accounting changed")
    require(absolute.get("H_checkpoint_sha256") ==
            h_packet.get("checkpoint_sha256")
            and absolute.get("tail_checkpoint_sha256") ==
            tail.get("checkpoint_sha256"),
            "absolute M15 checkpoint links changed")
    metadata = absolute.get("metadata", {})
    require(metadata.get("identity") ==
            "H_m15_inlet/Rm^2 + shifted_tail_Rm^2"
            and metadata.get("H_m15_inlet") ==
            "M15(sc/2)+firstlowerRQy"
            and metadata.get("numeric_H_source_bound") is True
            and metadata.get("terminal_target") ==
            "Rp (not closed); current_radius=Rm (intermediate)"
            and metadata.get("global_closure") is False
            and metadata.get("NS_claim") is False,
            "absolute M15 scope changed")
    require(_checkpoint_sha(absolute) == absolute.get("checkpoint_sha256"),
            "absolute M15 checkpoint hash changed")


def validate_tail_packet(tail: dict, source: dict, parent: dict,
                         outer: dict, signed: dict) -> None:
    require(tail.get("schema") in {
        "original_outer_pressure_enclosure_linear_tail_v1",
        "original_outer_pressure_enclosure_linear_tail_v2",
    },
            "new tail checkpoint schema changed")
    require(tail.get("source_executed") is True
            and tail.get("accepted_ancestors_replayed") is False,
            "new tail execution provenance changed")
    old = parent["actual_pressure_prefix_contraction"]["new_partial_rows"]
    require(tail.get("old_parent_new_partial_rows") == old,
            "tail did not retain parent new_partial_rows verbatim")
    require(tail.get("outerS1_factor") == {"numerator": 1, "denominator": 2},
            "tail outer S1 half factor changed")
    shifted = tail.get("shifted_tail_Rm2")
    require(isinstance(shifted, list) and len(shifted) == 2,
            "shifted tail lacks explicit ordinary-Z rows")
    raw_iv(shifted[0], "shifted_tail_Rm2[0]")
    raw_iv(shifted[1], "shifted_tail_Rm2[1]")
    outer_hash = (tail.get("outer_checkpoint_sha256")
                  or tail.get("coarse_outer_checkpoint_sha256"))
    require(outer_hash == outer.get("checkpoint_sha256"),
            "tail does not point to the saved outer checkpoint")
    signed_hash = (tail.get("signed_macro_checkpoint_sha256")
                   or tail.get("signed_macro_scale_checkpoint_sha256"))
    require(signed_hash == signed.get("checkpoint_sha256"),
            "tail does not use the saved signed macro-scale result")
    require((tail.get("signed_macro_scale_used")
             or tail.get("signed_macro_scale_consumed")) is True,
            "tail signed macro-scale consumption flag missing")
    require((tail.get("coarse_outer_bound_only")
             or tail.get("coarse_outer_result_preserved_as_bound_only")) is True,
            "tail did not preserve coarse outer result as bound-only")
    require(tail.get("outerS1_factor") == {"numerator": 1, "denominator": 2},
            "tail consumed outer S1 more than once or with the wrong factor")
    require(_checkpoint_sha(tail) == tail.get("checkpoint_sha256"),
            "tail checkpoint hash changed")


def production_artifact(path: Path) -> dict:
    require(path.exists(), "root-only production checkpoint is missing: " + str(path))
    value = checkpoint_json(path)
    _checkpoint_sha(value)
    return value


def implementation_hashes() -> dict:
    return {
        RUNNER_PATH: sha256_file(ROOT / RUNNER_PATH),
        KERNEL_PATH: (sha256_file(ROOT / KERNEL_PATH)
                      if (ROOT / KERNEL_PATH).exists() else KERNEL_SHA256),
        SIGNED_MACRO_KERNEL_PATH: (sha256_file(ROOT / SIGNED_MACRO_KERNEL_PATH)
                                   if (ROOT / SIGNED_MACRO_KERNEL_PATH).exists()
                                   else SIGNED_MACRO_KERNEL_SHA256),
        H_KERNEL_PATH: (sha256_file(ROOT / H_KERNEL_PATH)
                        if (ROOT / H_KERNEL_PATH).exists()
                        else H_KERNEL_SHA256),
    }


def merge_pins(*maps: dict) -> dict:
    result = {}
    for mapping in maps:
        require(isinstance(mapping, dict), "pin map is not an object")
        for name, value in mapping.items():
            require(name not in result or result[name] == value,
                    "conflicting source pin: " + str(name))
            result[name] = value
    return result


def build_report(parent: dict, native: dict, source: dict, outer: dict, signed: dict,
                 tail: dict, h_packet: dict, absolute: dict) -> dict:
    validate_source_checkpoint(source, parent, native)
    validate_raw_outer_packet(outer)
    validate_signed_macro_packet(signed, outer, source)
    validate_tail_packet(tail, source, parent, outer, signed)
    validate_h_packet(h_packet, source)
    validate_absolute_packet(absolute, h_packet, tail)
    report = copy.deepcopy(parent)
    report.update({
        "stage": "outer_pressure_enclosure_and_H_anchored_absolute_m15",
        "parent_report": PARENT_NAME,
        "parent_report_sha256": PARENT_SHA256,
        "inspection_format_repair": copy.deepcopy(INSPECTION_FORMAT_REPAIR),
        "source_consumer_checkpoint": copy.deepcopy(source),
        "source_checkpoint_sha256": sha256_bytes(SOURCE_CHECKPOINT.read_bytes()),
        "actual_outer_pressure_enclosure": outer,
        "actual_outer_signed_macro_scale": signed,
        "actual_H_m15_Rm_transport": h_packet,
        "actual_Rm_shifted_tail": tail,
        "actual_Rm_absolute_m15": {
            "schema": "original_Rm_absolute_m15_v1",
            "passed": True,
            "source_only": False,
            "source_executed": True,
            "accepted_H_anchor_receipt": copy.deepcopy(source["accepted_H_anchor_receipt"]),
            "H_anchor_digest": H_ANCHOR_DIGEST,
            "first_lower_RQy_retained_opaque": True,
            "pressure_increment_adopted_over_Rm2": {
                "identity": "first_micro_lower_RQy/Rm^2 + shifted_tail_Rm^2",
                "numeric": False,
            },
            "absolute_M15_over_Rm2": {
                "identity": "H_m15_inlet/Rm^2 + shifted_tail_Rm^2",
                "H_m15_inlet": "M15(sc/2)+firstlowerRQy",
                "numeric": False,
                "numeric_bound": True,
                "numeric_enclosure": True,
            },
            "new_partial_rows": copy.deepcopy(tail["shifted_tail_Rm2"]),
            "new_absolute_rows": copy.deepcopy(absolute["new_absolute_rows"]),
            "H_transport_checkpoint_sha256": h_packet["checkpoint_sha256"],
            "absolute_checkpoint_sha256": absolute["checkpoint_sha256"],
            "raw_pressure_increment_added_to_H": False,
            "terminal_target": "Rp (not closed); current_radius=Rm (intermediate)",
            "global_closure": False,
            "NS_claim": False,
        },
        "new_partial_rows": copy.deepcopy(tail["shifted_tail_Rm2"]),
        "new_absolute_rows": copy.deepcopy(absolute["new_absolute_rows"]),
        "checkpoint_hashes": {
            "source": sha256_bytes(SOURCE_CHECKPOINT.read_bytes()),
            "outer": sha256_bytes(OUTER_CHECKPOINT.read_bytes()),
            "signed_macro": sha256_bytes(SIGNED_MACRO_CHECKPOINT.read_bytes()),
            "linear_tail": sha256_bytes(TAIL_CHECKPOINT.read_bytes()),
            "H_transport": sha256_bytes(H_CHECKPOINT.read_bytes()),
            "absolute_M15": sha256_bytes(ABSOLUTE_CHECKPOINT.read_bytes()),
        },
        "execution_seconds": {
            "outer_kernel": float(outer.get("execution_seconds", 0.0)),
            "signed_macro_scale": float(signed.get("execution_seconds", 0.0)),
            "linear_tail": float(tail.get("execution_seconds", 0.0)),
            "H_transport": float(h_packet.get("execution_seconds", 0.0)),
            "absolute_M15": float(absolute.get("execution_seconds", 0.0)),
        },
        "input_hashes": merge_pins(
            parent.get("input_hashes", {}),
            native.get("input_hashes", {}), {
            PARENT_NAME: PARENT_SHA256,
            NATIVE_NAME: NATIVE_SHA256,
            FIXTURE_PATH: FIXTURE_SHA256,
            SIGNED_MACRO_FIXTURE_PATH: SIGNED_MACRO_FIXTURE_SHA256,
            H_FIXTURE_PATH: H_FIXTURE_SHA256,
        }, INSPECTION_REPAIR_INPUT_PINS),
        "implementation_hashes": merge_pins(
            parent.get("implementation_hashes", {}),
            native.get("implementation_hashes", {}),
            implementation_hashes()),
        "source_only": False,
        "accepted_ancestors_replayed": False,
        "source_bindings_snapshot": {
            "outer_graph_source_receipt": copy.deepcopy(
                source["outer_graph_source_receipt"]),
            "macro_Q_target_source_receipt": copy.deepcopy(
                source["macro_Q_target_source_receipt"]),
            "first_micro_seam_source_receipt": copy.deepcopy(
                source["first_micro_seam_source_receipt"]),
            "H_transport_contract": copy.deepcopy(
                source["H_transport_contract"]),
            "H_basis_rows_per_Z": {
                "signed_approximation": 12,
                "error_radius": 61,
            },
        },
        "execution_audit": {
            "accepted_ancestors_replayed": False,
            "accepted_producers_checkers_fixtures_replayed": False,
            "new_kernel_called_once": True,
            "new_kernel_checkpoint_saved_before_validation": True,
            "new_signed_macro_component_called_once": True,
            "new_signed_macro_checkpoint_saved_before_validation": True,
            "new_tail_checkpoint_saved_before_validation": True,
            "new_H_transport_called_once": True,
            "new_H_checkpoint_saved_before_validation": True,
            "new_absolute_checkpoint_saved_before_validation": True,
            "execution_seconds": {
                "outer_kernel": float(outer.get("execution_seconds", 0.0)),
                "signed_macro_scale": float(signed.get("execution_seconds", 0.0)),
                "linear_tail": float(tail.get("execution_seconds", 0.0)),
                "H_transport": float(h_packet.get("execution_seconds", 0.0)),
                "absolute_M15": float(absolute.get("execution_seconds", 0.0)),
            },
            "record_clock_seconds": 0,
        },
    })
    contract = copy.deepcopy(parent.get("scope", {}))
    contract.update({
        "outer_pressure_enclosure": True,
        "outer_signed_macro_scale": True,
        "H_m15_Rm_transport": True,
        "Rm_shifted_tail": True,
        "absolute_M15_bound": True,
        "full_pressure_numeric": False,
        "absolute_M15_numeric": False,
        "terminal_match": False,
        "global_five_moment_closure": False,
        "NS_claim": False,
    })
    report["scope"] = contract
    report["checks"] = {
        "passed": True,
        "source_preflight": True,
        "kernel_result_shape": True,
        "signed_macro_scale_result_shape": True,
        "tail_shape": True,
        "checkpoint_hashes_match": True,
        "opaque_lower_retained": True,
        "half_factor_consumed_once": True,
        "old_212_six_reserved_counts_preserved": True,
        "reviewed_kernel_sha256": KERNEL_SHA256,
        "reviewed_fixture_sha256": FIXTURE_SHA256,
        "reviewed_signed_macro_kernel_sha256": SIGNED_MACRO_KERNEL_SHA256,
        "reviewed_signed_macro_fixture_sha256": SIGNED_MACRO_FIXTURE_SHA256,
        "reviewed_H_kernel_sha256": H_KERNEL_SHA256,
        "reviewed_H_fixture_sha256": H_FIXTURE_SHA256,
    }
    return report


def compute_new() -> dict:
    """Run the root-only new-kernel chain after the immutable source handoff.

    This function is intentionally never called by the data-only preflight
    tests.  Every producer result is checkpointed before its validator runs;
    a threshold failure therefore stops the chain without replay or fallback.
    """
    parent, native = authenticated_inputs()
    source = attach_checkpoint(build_source_consumer_packet(parent, native), "source_consumer")
    digest = write_immutable_json(SOURCE_CHECKPOINT, source)
    source = checkpoint_json(SOURCE_CHECKPOINT)
    validate_source_checkpoint(source, parent, native)

    # Numerical imports occur only after authenticated source setup and the
    # source checkpoint.  The inspection/report stages never import them.
    import mpmath as mp
    from mpmath.ctx_iv import MPIntervalContext

    old_precision = mp.mp.prec
    mp.mp.prec = 53
    c = MPIntervalContext()
    c.prec = CONTEXT_PRECISION_BITS
    try:
        stage_started = time.perf_counter()
        outer = _outer_kernel_checkpoint(c, source, digest)
        outer = _checkpoint_with_seconds(
            outer, time.perf_counter() - stage_started)
        write_immutable_json(OUTER_CHECKPOINT, outer)
        outer = production_artifact(OUTER_CHECKPOINT)
        validate_raw_outer_packet(outer)

        stage_started = time.perf_counter()
        signed = _signed_macro_checkpoint(
            c, source, outer, digest, outer["checkpoint_sha256"])
        signed = _checkpoint_with_seconds(
            signed, time.perf_counter() - stage_started)
        write_immutable_json(SIGNED_MACRO_CHECKPOINT, signed)
        signed = production_artifact(SIGNED_MACRO_CHECKPOINT)
        validate_signed_macro_packet(signed, outer, source)

        stage_started = time.perf_counter()
        tail = _tail_checkpoint(
            c, parent, signed, outer["checkpoint_sha256"],
            signed["checkpoint_sha256"])
        tail = _checkpoint_with_seconds(
            tail, time.perf_counter() - stage_started)
        write_immutable_json(TAIL_CHECKPOINT, tail)
        tail = production_artifact(TAIL_CHECKPOINT)
        validate_tail_packet(tail, source, parent, outer, signed)

        stage_started = time.perf_counter()
        try:
            h_packet = _h_checkpoint(c, source, digest)
        except Exception as exc:
            failure = _production_checkpoint({
                "schema": "original_H_m15_Rm_transport_checkpoint_v1",
                "passed": False,
                "source_executed": True,
                "accepted_ancestors_replayed": False,
                "failure": {
                    "stage": "H_transport",
                    "threshold_failure_is_terminal": True,
                    "fallback_or_replay": False,
                    "error_type": type(exc).__name__,
                    "error": str(exc),
                },
            }, "H_transport_failure", source_checkpoint_sha256=digest)
            failure = _checkpoint_with_seconds(
                failure, time.perf_counter() - stage_started)
            write_immutable_json(H_CHECKPOINT, failure)
            raise DataOnlyBlocked(
                "H transport failed; immutable failure checkpoint written") from exc
        h_packet = _checkpoint_with_seconds(
            h_packet, time.perf_counter() - stage_started)
        write_immutable_json(H_CHECKPOINT, h_packet)
        h_packet = production_artifact(H_CHECKPOINT)
        validate_h_packet(h_packet, source)

        stage_started = time.perf_counter()
        absolute = _absolute_checkpoint(
            c, h_packet, tail, h_packet["checkpoint_sha256"],
            tail["checkpoint_sha256"])
        absolute = _checkpoint_with_seconds(
            absolute, time.perf_counter() - stage_started)
        write_immutable_json(ABSOLUTE_CHECKPOINT, absolute)
        absolute = production_artifact(ABSOLUTE_CHECKPOINT)
        validate_absolute_packet(absolute, h_packet, tail)
    finally:
        mp.mp.prec = old_precision
    return {
        "passed": True,
        "source_checkpoint_sha256": sha256_bytes(SOURCE_CHECKPOINT.read_bytes()),
        "outer_checkpoint_sha256": sha256_bytes(OUTER_CHECKPOINT.read_bytes()),
        "signed_macro_checkpoint_sha256": sha256_bytes(
            SIGNED_MACRO_CHECKPOINT.read_bytes()),
        "tail_checkpoint_sha256": sha256_bytes(TAIL_CHECKPOINT.read_bytes()),
        "H_checkpoint_sha256": sha256_bytes(H_CHECKPOINT.read_bytes()),
        "absolute_checkpoint_sha256": sha256_bytes(ABSOLUTE_CHECKPOINT.read_bytes()),
        "record_clock_seconds": 0,
    }


def write_report() -> dict:
    require(REPORT_PATH.exists() is False, "report exists; use --inspect-report")
    parent, native = authenticated_inputs()
    source = checkpoint_json(SOURCE_CHECKPOINT)
    validate_source_checkpoint(source, parent, native)
    outer = production_artifact(OUTER_CHECKPOINT)
    signed = production_artifact(SIGNED_MACRO_CHECKPOINT)
    tail = production_artifact(TAIL_CHECKPOINT)
    h_packet = production_artifact(H_CHECKPOINT)
    absolute = production_artifact(ABSOLUTE_CHECKPOINT)
    report = build_report(parent, native, source, outer, signed, tail,
                          h_packet, absolute)
    encoded = (json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False)
               + "\n").encode("utf-8")
    check_lf_single_eof(encoded, REPORT_NAME)
    try:
        with REPORT_PATH.open("xb") as stream:
            stream.write(encoded)
    except FileExistsError as exc:
        raise DataOnlyBlocked("report appeared during write: " + str(REPORT_PATH)) from exc
    return {"passed": True, "bytes": len(encoded), "sha256": sha256_bytes(encoded)}


def inspect_report() -> dict:
    parent, native = authenticated_inputs()
    require(REPORT_PATH.exists(), "report is missing: " + str(REPORT_PATH))
    report_raw = REPORT_PATH.read_bytes()
    check_lf_single_eof(report_raw, REPORT_NAME)
    report = json.loads(report_raw)
    require(isinstance(report, dict), "report must decode to an object")
    source = checkpoint_json(SOURCE_CHECKPOINT)
    outer = production_artifact(OUTER_CHECKPOINT)
    signed = production_artifact(SIGNED_MACRO_CHECKPOINT)
    tail = production_artifact(TAIL_CHECKPOINT)
    h_packet = production_artifact(H_CHECKPOINT)
    absolute = production_artifact(ABSOLUTE_CHECKPOINT)
    validate_source_checkpoint(source, parent, native)
    validate_raw_outer_packet(outer)
    validate_signed_macro_packet(signed, outer, source)
    validate_tail_packet(tail, source, parent, outer, signed)
    validate_h_packet(h_packet, source)
    validate_absolute_packet(absolute, h_packet, tail)
    require(report.get("actual_outer_pressure_enclosure") == outer,
            "report and outer checkpoint differ")
    require(report.get("source_consumer_checkpoint") == source,
            "report and source checkpoint differ")
    require(report.get("actual_outer_signed_macro_scale") == signed,
            "report and signed macro-scale checkpoint differ")
    require(report.get("actual_H_m15_Rm_transport") == h_packet,
            "report and H transport checkpoint differ")
    require(report.get("actual_Rm_shifted_tail") == tail,
            "report and tail checkpoint differ")
    absolute_report = report.get("actual_Rm_absolute_m15", {})
    require(absolute_report.get("new_absolute_rows") ==
            absolute.get("new_absolute_rows")
            and absolute_report.get("new_partial_rows") ==
            tail.get("shifted_tail_Rm2"),
            "report and absolute M15 checkpoint rows differ")
    require(report.get("input_hashes", {}).get(PARENT_NAME) == PARENT_SHA256
            and report.get("input_hashes", {}).get(NATIVE_NAME) == NATIVE_SHA256,
            "report input hashes changed")
    expected_inputs = merge_pins(
        parent.get("input_hashes", {}), native.get("input_hashes", {}), {
            PARENT_NAME: PARENT_SHA256,
            NATIVE_NAME: NATIVE_SHA256,
            FIXTURE_PATH: FIXTURE_SHA256,
            SIGNED_MACRO_FIXTURE_PATH: SIGNED_MACRO_FIXTURE_SHA256,
            H_FIXTURE_PATH: H_FIXTURE_SHA256,
        }, INSPECTION_REPAIR_INPUT_PINS)
    require(report.get("input_hashes") == expected_inputs,
            "report inherited/new input pins changed")
    require(report.get("implementation_hashes") == merge_pins(
        parent.get("implementation_hashes", {}),
        native.get("implementation_hashes", {}),
        implementation_hashes()),
        "report inherited/new implementation pins changed")
    require(report.get("implementation_hashes", {}).get(RUNNER_PATH) == sha256_file(ROOT / RUNNER_PATH),
            "report runner implementation hash changed")
    expected_checkpoint_hashes = {
        "source": sha256_bytes(SOURCE_CHECKPOINT.read_bytes()),
        "outer": sha256_bytes(OUTER_CHECKPOINT.read_bytes()),
        "signed_macro": sha256_bytes(SIGNED_MACRO_CHECKPOINT.read_bytes()),
        "linear_tail": sha256_bytes(TAIL_CHECKPOINT.read_bytes()),
        "H_transport": sha256_bytes(H_CHECKPOINT.read_bytes()),
        "absolute_M15": sha256_bytes(ABSOLUTE_CHECKPOINT.read_bytes()),
    }
    require(report.get("checkpoint_hashes") == expected_checkpoint_hashes,
            "report checkpoint hashes changed")
    require(report.get("execution_seconds") == {
        "outer_kernel": float(outer.get("execution_seconds", 0.0)),
        "signed_macro_scale": float(signed.get("execution_seconds", 0.0)),
        "linear_tail": float(tail.get("execution_seconds", 0.0)),
        "H_transport": float(h_packet.get("execution_seconds", 0.0)),
        "absolute_M15": float(absolute.get("execution_seconds", 0.0)),
    }, "report execution seconds changed")
    tree = ast.parse((ROOT / RUNNER_PATH).read_text(encoding="utf-8"), filename=RUNNER_PATH)
    require(isinstance(tree, ast.Module), "runner AST unavailable")
    return {
        "passed": True,
        "serialized_only": True,
        "bytes": len(report_raw),
        "sha256": sha256_bytes(report_raw),
        "report_sha256": sha256_bytes(report_raw),
        "source_checkpoint_sha256": sha256_bytes(SOURCE_CHECKPOINT.read_bytes()),
        "outer_checkpoint_sha256": sha256_bytes(OUTER_CHECKPOINT.read_bytes()),
        "signed_macro_checkpoint_sha256": sha256_bytes(SIGNED_MACRO_CHECKPOINT.read_bytes()),
        "tail_checkpoint_sha256": sha256_bytes(TAIL_CHECKPOINT.read_bytes()),
        "H_checkpoint_sha256": sha256_bytes(H_CHECKPOINT.read_bytes()),
        "absolute_checkpoint_sha256": sha256_bytes(ABSOLUTE_CHECKPOINT.read_bytes()),
        "accepted_producers_checkers_fixtures_replayed": False,
        "ast_checked": True,
        "canonical_exact_tuple_records": canonical_exact_tuple_records(report),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--compute-new", action="store_true",
                       help="run the root-only new-kernel chain after source preflight")
    group.add_argument("--write-report", action="store_true")
    group.add_argument("--inspect-report", action="store_true")
    args = parser.parse_args()
    started = time.perf_counter()
    try:
        if args.compute_new:
            result = compute_new()
        elif args.write_report:
            result = write_report()
        else:
            result = inspect_report()
    except DataOnlyBlocked as exc:
        print(json.dumps({"passed": False, "blocked": True,
                          "reason": str(exc),
                          "record_clock_seconds": 0}, ensure_ascii=False))
        return 2
    result.setdefault("record_clock_seconds", 0)
    result["elapsed_seconds_not_a_recorded_math_clock"] = round(time.perf_counter() - started, 6)
    print(json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
