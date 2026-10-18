"""Transport the saved four chart pressure window into the ``R_m^2`` scale.

This runner is deliberately a saved-data adapter.  It authenticates the four
saved parent reports and their inherited source pins before decoding them,
passes only the saved directed ``S1`` rows to the reviewed transport kernel,
and records the resulting partial pressure-root group.  ``--write-report`` and
``--inspect-report`` never import a numerical provider; inspection uses only
JSON, hashes, ASTs, and exact integer tuple checks.
"""
from __future__ import annotations

import argparse
import ast
import copy
import hashlib
import json
from collections.abc import Mapping, Sequence
from pathlib import Path
import sys
import time


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "experiments/root_st073"
RUNNER = "scripts/run_original_pressure_window_transport.py"
KERNEL = "src/openai_ns_reconstruction/original_pressure_window_transport.py"
FIXTURE = "original_pressure_window_transport_fixture.json"
FIXTURE_BYTES = 34984
FIXTURE_SHA = "6b70690ad1d5801529331a6d71c40741f3dc9f44e2a6aea201e29205646a38a1"
ADOPTION = "lei_ren_part1_paper_compliant_current_original_restoration_pressure_adoption.json"
FINITE_N = "lei_ren_part1_paper_compliant_current_original_directed_finite_N_controls.json"
NATIVE = "lei_ren_part1_paper_compliant_current_original_native_directed_controls_v2.json"
BACKGROUND = "lei_ren_part1_paper_compliant_current_original_restoration_background_I1.json"
REPORT_NAME = "lei_ren_part1_paper_compliant_current_original_pressure_window_transport.json"
REPORT = BASE / REPORT_NAME
CHECKPOINT = ROOT / "tmp/original_pressure_window_transport/new_math.json"
ADOPTION_BYTES = 674832
ADOPTION_SHA = "2f10516f66168e513cca3ea16ce25d82bc78a4850c37c30ed867ffda0f558234"
FINITE_N_BYTES = 1368720
FINITE_N_SHA = "152d211cd137e7e0498fd157a0ac1467541dc6214ae59c13af62844112003d9a"
NATIVE_BYTES = 42177499
NATIVE_SHA = "10abe68b98fdf87709a7a7c9f27103523701e746c357ae62386038700052490a"
BACKGROUND_BYTES = 31700329
BACKGROUND_SHA = "1c962bf1b2deb48971edf8bf8cc19df9678c2f7b09b483096b035fe8cb29f3fd"
REVIEWED_KERNEL_SHA = "a9efc7f6808455a6af7dad26c43a7fdf94062d90371b14cb3a8cbc18c98d077e"
CHARTS = ("long_reshape", "reference", "restoration", "post_restore")
ORDERS = (0, 1)
ZERO = (0, 0, 0, 0)
RADIUS = (0, 1, -1024, 1)
NEG_RADIUS = (1, 1, -1024, 1)
FALSE_SCOPE = (
    "full_pressure_numeric",
    "complete_pressure_numeric",
    "global_closure",
    "global_five_moment_closure",
    "independent_reserved_full_pulse_m15_matched",
    "full_goal_complete",
)
HALF_S1_HANDLES = {
    0: (185674, 185676, 186070, 186072),
    1: (185675, 185677, 186071, 186073),
}
EXTERNAL_RM_HANDLES = {0: 183696, 1: 183697}
LOCAL_S1_HANDLES = {
    "long_reshape": (166093, 166094),
    "reference": (171813, 171814),
    "restoration": (178037, 178038),
    "post_restore": (183640, 183641),
}
PARENTS = {
    ADOPTION: (ADOPTION_BYTES, ADOPTION_SHA),
    FINITE_N: (FINITE_N_BYTES, FINITE_N_SHA),
    NATIVE: (NATIVE_BYTES, NATIVE_SHA),
    BACKGROUND: (BACKGROUND_BYTES, BACKGROUND_SHA),
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def digest(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def sha(path: Path) -> str:
    return digest(path.read_bytes())


def _safe_relative(name: str) -> Path:
    path = Path(name)
    require(not path.is_absolute(), "Absolute path is not allowed: " + name)
    require(".." not in path.parts, "Path traversal is not allowed: " + name)
    return path


def _input_path(name: str) -> Path:
    path = Path(name)
    require(not path.is_absolute(), "Absolute inherited input path: " + name)
    if ".." in path.parts:
        require(name == "../../src/openai_ns_reconstruction/schedule_pressure.py",
                "Unexpected inherited input traversal: " + name)
        return (BASE / path).resolve()
    if path.parts and path.parts[0] == "experiments":
        return (ROOT / path).resolve()
    return (BASE / path.name).resolve()


def _implementation_path(name: str) -> Path:
    path = _safe_relative(name)
    return (ROOT / path).resolve()


def _authenticated_json(name: str, expected_bytes: int, expected_sha: str) -> tuple[dict, bytes]:
    path = BASE / name
    raw = path.read_bytes()
    require(len(raw) == expected_bytes and digest(raw) == expected_sha,
            "Authenticated parent changed before decode: " + name)
    value = json.loads(raw)
    require(isinstance(value, dict), "Authenticated report is not an object: " + name)
    return value, raw


def _authenticate_parent_files() -> tuple[dict[str, dict], dict[str, bytes]]:
    # Read and authenticate all bytes before decoding any parent JSON.
    raws: dict[str, bytes] = {}
    for name, (expected_bytes, expected_sha) in PARENTS.items():
        raw = (BASE / name).read_bytes()
        require(len(raw) == expected_bytes and digest(raw) == expected_sha,
                "Authenticated parent changed before decode: " + name)
        raws[name] = raw
    values = {name: json.loads(raw) for name, raw in raws.items()}
    require(all(isinstance(value, dict) for value in values.values()),
            "Authenticated parent is not an object")
    adoption = values[ADOPTION]
    _authenticate_inherited_pins(adoption)
    return values, raws


def _authenticate_inherited_pins(adoption: Mapping[str, object]) -> None:
    inputs = adoption.get("input_hashes")
    implementations = adoption.get("implementation_hashes")
    require(isinstance(inputs, dict) and isinstance(implementations, dict),
            "Adoption inherited source pins are missing")
    for name, expected in inputs.items():
        require(isinstance(name, str) and isinstance(expected, str),
                "Malformed inherited input pin")
        path = _input_path(name)
        require(path.is_file() and sha(path) == expected,
                "Inherited input pin changed: " + name)
    for name, expected in implementations.items():
        require(isinstance(name, str) and isinstance(expected, str),
                "Malformed inherited implementation pin")
        path = _implementation_path(name)
        require(path.is_file() and sha(path) == expected,
                "Inherited implementation pin changed: " + name)


def _preflight_file(relative: str, expected_sha: str | None = None) -> dict[str, object]:
    path = ROOT / _safe_relative(relative)
    raw = path.read_bytes()
    require(b"\r" not in raw and raw.endswith(b"\n") and not raw.endswith(b"\n\n"),
            "LF/EOF preflight failed: " + relative)
    tree = ast.parse(raw.decode("utf-8"), filename=relative)
    actual_sha = digest(raw)
    if expected_sha is not None:
        require(actual_sha == expected_sha, "Reviewed source changed: " + relative)
    return {
        "path": relative,
        "bytes": len(raw),
        "sha256": actual_sha,
        "ast_parse": True,
        "source_executed": False,
        "ast_nodes": len(list(ast.walk(tree))),
    }


def _fixture() -> dict:
    path = BASE / FIXTURE
    raw = path.read_bytes()
    require(len(raw) == FIXTURE_BYTES and digest(raw) == FIXTURE_SHA,
            "Saved transport fixture changed")
    value = json.loads(raw)
    require(isinstance(value, dict) and value.get("passed") is True,
            "Saved transport fixture did not pass")
    require(value.get("module_sha256") == REVIEWED_KERNEL_SHA,
            "Saved transport fixture/kernel pin mismatch")
    checks = value.get("checks")
    require(isinstance(checks, dict) and checks and all(item is True for item in checks.values()),
            "Saved transport fixture checks are incomplete")
    return value


def preflight() -> dict[str, object]:
    runner = _preflight_file(RUNNER)
    kernel = _preflight_file(KERNEL, REVIEWED_KERNEL_SHA)
    fixture = _fixture()
    return {"runner": runner, "kernel": kernel, "fixture_sha256": FIXTURE_SHA,
            "fixture_bytes": FIXTURE_BYTES, "fixture_passed": fixture["passed"]}


def _endpoint(value: object, label: str = "endpoint") -> tuple[int, int, int, int]:
    require(isinstance(value, (list, tuple)) and len(value) == 4,
            label + " must be a four component exact MPF tuple")
    require(all(type(part) is int for part in value), label + " tuple is not integer")
    result = tuple(value)
    sign, mantissa, exponent, bit_count = result
    require(sign in (0, 1) and mantissa >= 0 and bit_count == mantissa.bit_length(),
            label + " tuple is noncanonical")
    require((mantissa == 0 and result == ZERO) or (mantissa > 0 and mantissa & 1),
            label + " tuple is noncanonical")
    return result


def _cmp_abs(left: tuple[int, int, int, int], right: tuple[int, int, int, int]) -> int:
    _, lm, le, _ = left
    _, rm, re, _ = right
    if lm == 0 or rm == 0:
        return (lm > 0) - (rm > 0)
    ltop, rtop = lm.bit_length() + le, rm.bit_length() + re
    if ltop != rtop:
        return (ltop > rtop) - (ltop < rtop)
    if le >= re:
        lv, rv = lm << (le - re), rm
    else:
        lv, rv = lm, rm << (re - le)
    return (lv > rv) - (lv < rv)


def _cmp_exact(left: tuple[int, int, int, int], right: tuple[int, int, int, int]) -> int:
    if left[1] == 0:
        if right[1] == 0:
            return 0
        return 1 if right[0] else -1
    if right[1] == 0:
        return -1 if left[0] else 1
    if left[0] != right[0]:
        return -1 if left[0] else 1
    result = _cmp_abs(left, right)
    return -result if left[0] else result


def _record(value: object, label: str = "record") -> dict[str, object]:
    require(isinstance(value, dict), label + " must be an object")
    lower = _endpoint(value.get("lower_exact_mpf_tuple"), label + ".lower")
    upper = _endpoint(value.get("upper_exact_mpf_tuple"), label + ".upper")
    require(_cmp_exact(lower, upper) <= 0, label + " is reversed")
    require(value.get("directed") is True, label + " must be directed")
    return {"lower_exact_mpf_tuple": list(lower), "upper_exact_mpf_tuple": list(upper),
            "directed": True}


def _record_equal(left: object, right: object) -> bool:
    if not isinstance(left, dict) or not isinstance(right, dict):
        return False
    try:
        a = _record(left)
        b = _record(right)
    except (TypeError, ValueError):
        return False
    return a == b


def _records(value: object) -> int:
    if isinstance(value, dict):
        count = 0
        if "lower_exact_mpf_tuple" in value or "upper_exact_mpf_tuple" in value:
            _record(value)
            count = 1
        return count + sum(_records(item) for item in value.values()
                           if isinstance(item, (dict, list)))
    if isinstance(value, list):
        return sum(_records(item) for item in value if isinstance(item, (dict, list)))
    return 0


def _row_record(row: object, label: str, order: int | None = None) -> dict[str, object]:
    require(isinstance(row, dict), label + " row is not an object")
    if order is not None:
        require(row.get("ordinary_Z_order") == order, label + " row order changed")
    value = row.get("value", row)
    return _record(value, label + ".value")


def _row_pair(value: object, label: str) -> list[list[tuple[int, int, int, int]]]:
    require(isinstance(value, list) and len(value) == 2, label + " must be a Z0/Z1 pair")
    result = []
    for order, item in enumerate(value):
        if isinstance(item, dict) and "value" in item:
            item = item["value"]
        rec = _record(item, label + "[Z" + str(order) + "]")
        result.append([tuple(rec["lower_exact_mpf_tuple"]), tuple(rec["upper_exact_mpf_tuple"])])
    return result


def _saved_s1_rows(finite: Mapping[str, object]) -> dict[str, list[list[tuple[int, int, int, int]]]]:
    packet = finite.get("actual_directed_finite_N_controls")
    require(isinstance(packet, dict) and packet.get("passed") is True,
            "Finite-N packet is not a passed saved packet")
    charts = packet.get("charts")
    require(isinstance(charts, dict) and set(charts) == set(CHARTS),
            "Finite-N chart set changed")
    result = {}
    for chart in CHARTS:
        data = charts[chart]
        require(isinstance(data, dict), "Finite-N chart is not an object: " + chart)
        full = data.get("S1_full")
        require(isinstance(full, dict) and full.get("passed") is True,
                "Missing passed S1_full packet: " + chart)
        require(full.get("computed_from_B_and_D_once") is True
                and full.get("separate_S1_radius_not_added") is True
                and full.get("no_extra_radius") is True,
                "S1 radius provenance changed: " + chart)
        rows = full.get("rows")
        require(isinstance(rows, list) and len(rows) == 2, "S1 row count changed: " + chart)
        result[chart] = []
        for order, row in enumerate(rows):
            record = _row_record(row, chart + ".S1_full", order)
            result[chart].append([
                tuple(record["lower_exact_mpf_tuple"]),
                tuple(record["upper_exact_mpf_tuple"]),
            ])
        require(len(result[chart]) == 2, "S1 pair shape changed: " + chart)
        bfull, dfull = data.get("B_full"), data.get("D_full")
        require(isinstance(bfull, dict) and isinstance(dfull, dict),
                "B_full/D_full missing: " + chart)
        require(bfull.get("passed") is True and dfull.get("passed") is True,
                "B_full/D_full did not pass: " + chart)
        require(bfull.get("error_radius_is_not_a_selected_value") is True
                and dfull.get("computed_from_full_I1_once") is True
                and dfull.get("separate_D_error_not_added") is True
                and dfull.get("no_extra_D_radius") is True,
                "B/D radius provenance changed: " + chart)
        radius = bfull.get("symmetric_radius_added_once")
        require(_record_equal(radius, {
            "lower_exact_mpf_tuple": list(RADIUS),
            "upper_exact_mpf_tuple": list(RADIUS),
            "directed": True,
        }), "B_full radius changed: " + chart)
        for key in ("B_full", "D_full"):
            rows2 = data[key].get("rows")
            require(isinstance(rows2, list) and len(rows2) == 2,
                    key + " row count changed: " + chart)
            for order, row in enumerate(rows2):
                _row_record(row, chart + "." + key, order)
    return result


def _native_arithmetic(native: Mapping[str, object]) -> dict[str, dict[str, object]]:
    long_packet = native.get("actual_long_reference_background_I1_numeric")
    rest_packet = native.get("actual_restoration_post_jet_repair")
    require(isinstance(long_packet, dict) and isinstance(rest_packet, dict),
            "Native arithmetic packets are missing")
    long_charts = long_packet.get("charts")
    corrected = rest_packet.get("corrected_background")
    require(isinstance(long_charts, dict) and isinstance(corrected, dict),
            "Native chart packets are missing")
    rest_charts = corrected.get("charts")
    require(isinstance(rest_charts, dict), "Corrected native chart map is missing")
    result = {}
    for chart in CHARTS:
        source = long_charts if chart in ("long_reshape", "reference") else rest_charts
        item = source.get(chart)
        require(isinstance(item, dict) and isinstance(item.get("arithmetic"), dict),
                "Native arithmetic missing: " + chart)
        result[chart] = item["arithmetic"]
    return result


def _validate_baseline_bindings(finite: Mapping[str, object], native: Mapping[str, object]) -> None:
    packet = finite["actual_directed_finite_N_controls"]
    arithmetic = _native_arithmetic(native)
    for chart in CHARTS:
        baseline = packet["charts"][chart].get("baseline_centres")
        require(isinstance(baseline, dict), "Baseline centres missing: " + chart)
        native_rows = arithmetic[chart]
        for key in ("six_integrals", "B_over_W", "D", "S1_over_W"):
            require(baseline.get(key) == native_rows.get(key),
                    "Saved baseline/native arithmetic mismatch: " + chart + "." + key)
        require(baseline.get("signed_centres_retained_exactly") is True,
                "Signed baseline centres were not retained: " + chart)


def _reference_length_binding(native: Mapping[str, object], background: Mapping[str, object]) -> dict[str, object]:
    raw_tree = native.get("actual_native_directed_controls_v2", {}).get("raw_interval_tree", {})
    reference = raw_tree.get("long_reference", {}).get("charts", {}).get("reference", {})
    scale = reference.get("scale", {})
    length = scale.get("reference_length")
    accepted = background.get("actual_long_reference_background_I1_numeric", {}).get("charts", {}).get("reference", {}).get("length")
    require(isinstance(length, dict) and isinstance(accepted, dict),
            "Saved reference length record is missing")
    _record(length, "native reference_length")
    _record(accepted, "accepted reference background length")
    require(length == accepted, "Native reference_length differs from accepted reference length")
    return {
        "native_path": "actual_native_directed_controls_v2.raw_interval_tree.long_reference.charts.reference.scale.reference_length",
        "accepted_path": "actual_long_reference_background_I1_numeric.charts.reference.length",
        "record": copy.deepcopy(length),
        "exact_record_equal": True,
    }


def _axial_domain_binding(native: Mapping[str, object]) -> dict[str, object]:
    raw_tree = native.get("actual_native_directed_controls_v2", {}).get("raw_interval_tree", {})
    shared = raw_tree.get("shared", {})
    z = shared.get("z")
    require(isinstance(z, dict), "Native shared axial domain record is missing")
    canonical = _record(z, "native shared axial domain")
    require(z.get("lower_display") == "0.37" and z.get("upper_display") == "0.372",
            "Native shared axial domain changed")
    return {
        "source_path": "actual_native_directed_controls_v2.raw_interval_tree.shared.z",
        "record": copy.deepcopy(z),
        "canonical_record": canonical,
        "whole_Z": True,
    }


def _external_binding(adoption: Mapping[str, object]) -> dict[str, object]:
    """Bind the coefficient -1 source before any new arithmetic is imported."""
    packet = adoption.get("actual_restoration_pressure_adoption")
    require(isinstance(packet, dict), "Adoption nested packet is missing")
    external = packet.get("external_Rm_normalized_C1")
    require(isinstance(external, dict), "Adoption external Rm packet is missing")
    require(external.get("status") == "BOUND"
            and external.get("coefficient") == -1
            and external.get("source_handles") == [183696, 183697],
            "Adoption external Rm source changed")
    enclosure = external.get("enclosure")
    canonical = _record(enclosure, "adoption external Rm enclosure")
    require(enclosure.get("ordinary_Z_orders") == [0, 1]
            and enclosure.get("enclosure_not_actual_zero_or_value") is True,
            "Adoption external Rm enclosure metadata changed")
    # The new operation is authorized only by this saved bound.  Fail before
    # importing the provider if the inherited packet differs from the required
    # directed +/- 2^-1024 enclosure.
    require(canonical["lower_exact_mpf_tuple"] == list(NEG_RADIUS)
            and canonical["upper_exact_mpf_tuple"] == list(RADIUS),
            "Adoption external Rm enclosure is not +/- 2^-1024")
    return {
        "source_packet": copy.deepcopy(external),
        "enclosure": canonical,
        "source_report": ADOPTION,
    }


def _source_partition(background: Mapping[str, object]) -> dict[str, object]:
    certificate = background.get("actual_restoration_correlated_pressure", {}).get("exact_linearization")
    require(isinstance(certificate, dict) and certificate.get("passed") is True,
            "Restoration exact linearization certificate is missing")
    rows = certificate.get("rows")
    require(isinstance(rows, list) and len(rows) == 2,
            "Restoration exact linearization row count changed")
    partition_rows = []
    for order, row in enumerate(rows):
        require(isinstance(row, dict) and row.get("ordinary_Z_order") == order,
                "Linearization ordinary-Z order changed")
        require(row.get("reduced_handle") == (186096 if order == 0 else 186099),
                "Adopted reduced root changed")
        atoms = row.get("signed_atoms")
        require(isinstance(atoms, dict), "Signed atom certificate missing")
        atom_values = {int(key): value for key, value in atoms.items()}
        expected_group = HALF_S1_HANDLES[order]
        expected_external = EXTERNAL_RM_HANDLES[order]
        require(all(atom_values.get(handle) == 1 for handle in expected_group),
                "Half-S1 signed atom group changed")
        require(atom_values.get(expected_external) == -1,
                "External Rm coefficient changed")
        remaining = {handle: coefficient for handle, coefficient in atom_values.items()
                     if handle not in expected_group and handle != expected_external}
        require(len(remaining) == 7, "Opaque signed atom count changed")
        partition_rows.append({
            "ordinary_Z_order": order,
            "adopted_reduced_root": row["reduced_handle"],
            "half_S1_atoms": {str(handle): 1 for handle in expected_group},
            "half_S1_group_coefficient": 1,
            "external_Rm_atom": {str(expected_external): -1},
            "remaining_signed_atoms": {str(handle): coefficient for handle, coefficient in remaining.items()},
            "remaining_atoms_opaque": True,
        })
    return {
        "certificate_path": "background.actual_restoration_correlated_pressure.exact_linearization",
        "full_saved_certificate": copy.deepcopy(certificate),
        "adopted_roots": [186096, 186099],
        "ordinary_Z_orders": [0, 1],
        "rows": partition_rows,
        "local_S1_handles": {chart: list(handles) for chart, handles in LOCAL_S1_HANDLES.items()},
        "partial_factor": {"numerator": 1, "denominator": 2, "formula": "half_S1"},
        "remaining_atoms_are_opaque": True,
        "no_opaque_parent_root_window_added": True,
    }


def _source_scale_proof(native: Mapping[str, object], reference_binding: Mapping[str, object]) -> dict[str, object]:
    return {
        "source": "native_directed_controls_v2 C6 endpoint identities",
        "source_paths": [
            "actual_native_directed_controls_v2.physical_geometry.formula_source.reference",
            "actual_restoration_native_endpoints.source_guards",
            "actual_native_directed_controls_v2.raw_interval_tree.restoration.charts.*.scales",
        ],
        "C6_definitions": {
            "Rsh": "110*exp(T)",
            "Rz": "110*exp(logref-8)",
            "Rrestore": "110*exp(logref-7)",
            "Rm": "110*exp(logref-6)",
            "Lref": "logref-T-8",
        },
        "exact_algebra": {
            "Rsh_squared_over_Rm_squared": "exp(-2*(Lref+2))",
            "Rz_squared_over_Rm_squared": "exp(-4)",
            "Rrestore_squared_over_Rm_squared": "exp(-2)",
        },
        "exact_source_strings": {
            "Rsh": "Rsh=110exp(T)",
            "Rz": "Rz=110exp(logref-8)",
            "Rrestore": "Rrestore=110exp(logref-7)",
            "Rm": "Rm=110exp(logref-6)",
            "Lref": "Lref=logref-T-8",
            "Rsh_squared_over_Rm_squared": "Rsh^2/Rm^2=e^-2(Lref+2)",
            "Rz_squared_over_Rm_squared": "Rz^2/Rm^2=e^-4",
            "Rrestore_squared_over_Rm_squared": "Rrestore^2/Rm^2=e^-2",
        },
        "geometry_Z_independent": True,
        "reference_length_record_reused": True,
        "reference_length_binding": copy.deepcopy(reference_binding),
        "rounded_log_R_subtraction": False,
        "huge_exponential_materialized": False,
    }


def _setup() -> dict[str, object]:
    values, raws = _authenticate_parent_files()
    adoption, finite, native, background = (values[ADOPTION], values[FINITE_N], values[NATIVE], values[BACKGROUND])
    nested_adoption = adoption.get("actual_restoration_pressure_adoption")
    require(isinstance(nested_adoption, dict)
            and nested_adoption.get("schema") == "original_restoration_pressure_adoption_v1",
            "Adoption nested packet schema changed")
    if adoption.get("schema") is not None:
        require(adoption.get("schema") == nested_adoption.get("schema"),
                "Adoption parent schema changed")
    require(adoption.get("source_only") is True and adoption.get("accepted_ancestors_replayed") is False,
            "Adoption parent scope changed")
    _validate_baseline_bindings(finite, native)
    s1_rows = _saved_s1_rows(finite)
    reference_binding = _reference_length_binding(native, background)
    axial_domain = _axial_domain_binding(native)
    external_binding = _external_binding(adoption)
    partition = _source_partition(background)
    scale_proof = _source_scale_proof(native, reference_binding)
    namespace = nested_adoption.get("namespace")
    require(isinstance(namespace, dict), "Adoption source namespace is missing")
    return {
        "values": values,
        "raws": raws,
        "adoption": adoption,
        "finite": finite,
        "native": native,
        "background": background,
        "s1_rows": s1_rows,
        "reference_binding": reference_binding,
        "axial_domain": axial_domain,
        "external_binding": external_binding,
        "source_namespace": copy.deepcopy(namespace),
        "source_partition": partition,
        "source_scale_proof": scale_proof,
    }


def _saved_lref(setup: Mapping[str, object]) -> list[tuple[int, int, int, int]]:
    binding = setup["reference_binding"]
    record = binding["record"]
    return [tuple(record["lower_exact_mpf_tuple"]), tuple(record["upper_exact_mpf_tuple"])]


def _serialize_interval(value: object, label: str = "interval") -> dict[str, object]:
    rows = getattr(value, "_mpi_", None)
    require(rows is not None, label + " is not an interval")
    if len(rows) == 4 and all(type(part) is int for part in rows):
        rows = (rows, rows)
    require(len(rows) == 2, label + " has malformed endpoints")
    lower, upper = _endpoint(tuple(rows[0]), label + ".lower"), _endpoint(tuple(rows[1]), label + ".upper")
    require(_cmp_exact(lower, upper) <= 0, label + " is reversed")
    return {"lower_exact_mpf_tuple": list(lower), "upper_exact_mpf_tuple": list(upper), "directed": True}


def _compute_new(setup: Mapping[str, object], audit: Mapping[str, object]) -> dict[str, object]:
    require(not CHECKPOINT.exists() and not REPORT.exists(),
            "Saved pressure window calculation exists; replay/overwrite refused")
    sys.path.insert(0, str(ROOT / "src"))
    from mpmath import make_mpf
    from mpmath.ctx_iv import MPIntervalContext
    from openai_ns_reconstruction import original_pressure_window_transport as provider

    context = MPIntervalContext()
    context.prec = 1600
    saved_rows = setup["s1_rows"]
    rows_by_chart = {
        chart: tuple(context.mpf([make_mpf(pair[0]), make_mpf(pair[1])])
                     for pair in saved_rows[chart])
        for chart in CHARTS
    }
    external_binding = setup["external_binding"]
    external_source = external_binding["source_packet"]
    external_enclosure = external_binding["enclosure"]
    started = time.monotonic()
    transported = provider.transport(context, rows_by_chart, _saved_lref(setup),
                                      geometry_Z_independent=True, radius_bits=1024)
    # The required partial group is half of the transported S1 sum plus the
    # retained coefficient -1 external Rm bound [-2^-1024,+2^-1024].
    external = context.mpf([make_mpf(tuple(external_enclosure["lower_exact_mpf_tuple"])),
                             make_mpf(tuple(external_enclosure["upper_exact_mpf_tuple"]))])
    half = context.mpf(1) / context.mpf(2)
    raw_rows = transported["rawIVresult"]
    require(isinstance(raw_rows, Sequence) and len(raw_rows) == 2,
            "Transport kernel did not return two ordinary rows")
    partial_rows = tuple(row * half + external for row in raw_rows)
    packet = {
        "schema": "original_pressure_window_transport_v1",
        "passed": True,
        "ordinary_Z_orders": [0, 1],
        "charts": list(CHARTS),
        "normalization": "own_Rm^2",
        "input_S1_rows": {chart: [_serialize_interval(item, chart + ".S1")
                                   for item in rows_by_chart[chart]] for chart in CHARTS},
        "transport_raw_rows": [_serialize_interval(item, "transport_raw[Z" + str(order) + "]")
                               for order, item in enumerate(raw_rows)],
        "transport_metadata": copy.deepcopy(transported["metadata"]),
        "partial_window_rows": [_serialize_interval(item, "partial[Z" + str(order) + "]")
                                for order, item in enumerate(partial_rows)],
        "external_Rm": {
            **copy.deepcopy(external_source),
            "source_report": external_binding["source_report"],
            "normalized_by": "Wpost=Rm^2",
            "enclosure": {
                **copy.deepcopy(external_enclosure),
                "ordinary_Z_orders": [0, 1],
                "enclosure_not_actual_zero_or_value": True,
            },
            "rows": [{
                "ordinary_Z_order": order,
                "value": {
                    "lower_exact_mpf_tuple": list(external_enclosure["lower_exact_mpf_tuple"]),
                    "upper_exact_mpf_tuple": list(external_enclosure["upper_exact_mpf_tuple"]),
                    "directed": True,
                },
            } for order in ORDERS],
            "coefficient_minus_one_preserved": True,
        },
        "source_namespace": copy.deepcopy(setup["source_namespace"]),
        "axial_domain": copy.deepcopy(setup["axial_domain"]),
        "source_partition": copy.deepcopy(setup["source_partition"]),
        "source_scale_proof": copy.deepcopy(setup["source_scale_proof"]),
        "functional_identity": {
            "formula": "adopted_root/Rm^2 = opaque_prefix/Rm^2 + partial_window_Rm_normalized_C1",
            "ordinary_Z_orders": [0, 1],
            "geometry": "frozen",
            "opaque_prefix_retained": True,
            "partial_window_is_function_enclosure": True,
            "point_centers_selected": False,
            "full_pressure_numeric": False,
        },
        "source_contract": {
            "saved_S1_rows_consumed_exactly": True,
            "finite_N_radius_added_to_B_once": True,
            "D_from_full_I1_once": True,
            "S1_from_B_plus_D_once": True,
            "source_rows_recomputed": False,
            "accepted_producer_or_helper_replayed": False,
            "opaque_parent_root_window_added": False,
            "geometry_Z_independent": True,
            "new_values_are_directed_function_enclosures": True,
            "point_centers_selected": False,
        },
        "claims": {
            "partial_pressure_root_group": True,
            "full_pressure_closure": False,
            "global_five_moment_closure": False,
            "global_closure": False,
            "full_goal_complete": False,
        },
    }
    # Persist the only new arithmetic before validating or inspecting its output.
    checkpoint = {
        "payload": packet,
        "interval_precision_bits": context.prec,
        "execution_seconds": time.monotonic() - started,
        "kernel_sha256": REVIEWED_KERNEL_SHA,
        "computation_runner_sha256": sha(ROOT / RUNNER),
        "parent_report": ADOPTION,
        "parent_report_sha256": ADOPTION_SHA,
        "source_report_sha256": {
            ADOPTION: ADOPTION_SHA,
            FINITE_N: FINITE_N_SHA,
            NATIVE: NATIVE_SHA,
            BACKGROUND: BACKGROUND_SHA,
        },
        "accepted_ancestors_replayed": False,
        "source_execution": False,
    }
    _write_json(CHECKPOINT, checkpoint)
    return checkpoint


def _write_json(path: Path, value: Mapping[str, object]) -> None:
    require(not path.exists(), "Saved output exists; replay/overwrite refused: " + str(path))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
                    encoding="utf-8", newline="\n")


def _expected_inputs(adoption: Mapping[str, object]) -> dict[str, str]:
    result = dict(adoption["input_hashes"])
    for name, expected in ((ADOPTION, ADOPTION_SHA), (FINITE_N, FINITE_N_SHA), (NATIVE, NATIVE_SHA),
                           (BACKGROUND, BACKGROUND_SHA), (FIXTURE, FIXTURE_SHA)):
        if name in result:
            require(result[name] == expected, "Inherited input pin conflicts: " + name)
        result[name] = expected
    return result


def _expected_implementations(adoption: Mapping[str, object]) -> dict[str, str]:
    result = dict(adoption["implementation_hashes"])
    for name, expected in ((KERNEL, REVIEWED_KERNEL_SHA), (RUNNER, sha(ROOT / RUNNER))):
        if name in result:
            require(result[name] == expected, "Inherited implementation pin conflicts: " + name)
        result[name] = expected
    return result


def _validate_packet(packet: object, setup: Mapping[str, object]) -> int:
    require(isinstance(packet, dict) and packet.get("schema") == "original_pressure_window_transport_v1",
            "Pressure window packet schema changed")
    require(packet.get("passed") is True and packet.get("ordinary_Z_orders") == [0, 1],
            "Pressure window packet did not pass")
    require(packet.get("charts") == list(CHARTS), "Pressure window chart order changed")
    input_rows = packet.get("input_S1_rows")
    require(isinstance(input_rows, dict) and set(input_rows) == set(CHARTS),
            "Saved S1 input rows missing")
    for chart in CHARTS:
        expected = [{
            "lower_exact_mpf_tuple": list(item[0]),
            "upper_exact_mpf_tuple": list(item[1]),
            "directed": True,
        } for item in setup["s1_rows"][chart]]
        actual = input_rows[chart]
        require(isinstance(actual, list) and len(actual) == 2,
                "Serialized S1 row count changed: " + chart)
        require(all(_record_equal(row, exp) for row, exp in zip(actual, expected)),
                "Serialized S1 rows changed: " + chart)
    for key in ("transport_raw_rows", "partial_window_rows"):
        rows = packet.get(key)
        require(isinstance(rows, list) and len(rows) == 2, key + " shape changed")
        for order, row in enumerate(rows):
            _record(row, key + "[Z" + str(order) + "]")
    metadata = packet.get("transport_metadata")
    require(isinstance(metadata, dict)
            and metadata.get("geometry_Z_independent") is True
            and metadata.get("radius_bits") == 1024,
            "Transport metadata changed")
    require(metadata.get("source_contract", {}).get("long_coefficient_exponentiated") is False,
            "Long coefficient was exponentiated")
    external = packet.get("external_Rm")
    expected_external = setup["external_binding"]["enclosure"]
    require(isinstance(external, dict) and external.get("coefficient") == -1
            and external.get("coefficient_minus_one_preserved") is True
            and external.get("status") == "BOUND"
            and external.get("source_handles") == [183696, 183697]
            and external.get("source_report") == ADOPTION,
            "External Rm source changed")
    enclosure = external.get("enclosure")
    require(isinstance(enclosure, dict)
            and enclosure.get("lower_exact_mpf_tuple") == expected_external["lower_exact_mpf_tuple"]
            and enclosure.get("upper_exact_mpf_tuple") == expected_external["upper_exact_mpf_tuple"]
            and enclosure.get("enclosure_not_actual_zero_or_value") is True,
            "External Rm bound changed")
    ext_rows = external.get("rows")
    require(isinstance(ext_rows, list) and len(ext_rows) == 2,
            "External Rm row count changed")
    for order, row in enumerate(ext_rows):
        _row_record(row, "external_Rm", order)
        require(row["value"]["lower_exact_mpf_tuple"] == expected_external["lower_exact_mpf_tuple"]
                and row["value"]["upper_exact_mpf_tuple"] == expected_external["upper_exact_mpf_tuple"],
                "External Rm per-row bound changed")
    require(packet.get("source_namespace") == setup["source_namespace"],
            "Pressure source namespace changed")
    require(packet.get("axial_domain") == setup["axial_domain"],
            "Axial domain binding changed")
    partition = packet.get("source_partition")
    require(partition == setup["source_partition"], "Source partition changed")
    scale_proof = packet.get("source_scale_proof")
    require(scale_proof == setup["source_scale_proof"], "Source scale proof changed")
    identity = packet.get("functional_identity")
    require(isinstance(identity, dict)
            and identity.get("ordinary_Z_orders") == [0, 1]
            and identity.get("geometry") == "frozen"
            and identity.get("opaque_prefix_retained") is True
            and identity.get("partial_window_is_function_enclosure") is True
            and identity.get("point_centers_selected") is False,
            "Partial functional identity changed")
    contract = packet.get("source_contract")
    require(isinstance(contract, dict)
            and contract.get("saved_S1_rows_consumed_exactly") is True
            and contract.get("source_rows_recomputed") is False
            and contract.get("accepted_producer_or_helper_replayed") is False
            and contract.get("opaque_parent_root_window_added") is False,
            "Source contract changed")
    claims = packet.get("claims")
    require(isinstance(claims, dict) and claims.get("full_pressure_closure") is False
            and claims.get("global_five_moment_closure") is False
            and claims.get("global_closure") is False
            and claims.get("full_goal_complete") is False,
            "Pressure window claims overreached")
    return _records(packet)


def _validate_checkpoint(checkpoint: Mapping[str, object], setup: Mapping[str, object], audit: Mapping[str, object]) -> int:
    require(checkpoint.get("kernel_sha256") == REVIEWED_KERNEL_SHA,
            "Checkpoint kernel pin changed")
    require(checkpoint.get("parent_report_sha256") == ADOPTION_SHA
            and checkpoint.get("parent_report") == ADOPTION,
            "Checkpoint parent pin changed")
    require(checkpoint.get("source_report_sha256") == {
        ADOPTION: ADOPTION_SHA, FINITE_N: FINITE_N_SHA,
        NATIVE: NATIVE_SHA, BACKGROUND: BACKGROUND_SHA,
    }, "Checkpoint source pins changed")
    require(checkpoint.get("interval_precision_bits") == 1600
            and isinstance(checkpoint.get("execution_seconds"), (int, float))
            and checkpoint.get("execution_seconds") >= 0,
            "Checkpoint execution metadata changed")
    require(checkpoint.get("accepted_ancestors_replayed") is False
            and checkpoint.get("source_execution") is False,
            "Checkpoint execution audit changed")
    require(checkpoint.get("computation_runner_sha256") == sha(ROOT / RUNNER),
            "Checkpoint runner pin changed")
    packet = checkpoint.get("payload")
    return _validate_packet(packet, setup)


def _package(setup: Mapping[str, object], audit: Mapping[str, object]) -> dict[str, object]:
    require(not REPORT.exists(), "Report overwrite refused")
    require(CHECKPOINT.exists(), "New pressure window checkpoint is missing")
    checkpoint = json.loads(CHECKPOINT.read_bytes())
    require(isinstance(checkpoint, dict), "Checkpoint is not an object")
    _validate_checkpoint(checkpoint, setup, audit)
    adoption = setup["adoption"]
    report = copy.deepcopy(adoption)
    packet = checkpoint["payload"]
    report.update({
        "schema": "original_pressure_window_transport_v1",
        "stage": "current_original_pressure_window_transport",
        "parent_report": ADOPTION,
        "parent_report_sha256": ADOPTION_SHA,
        "input_hashes": _expected_inputs(adoption),
        "implementation_hashes": _expected_implementations(adoption),
        "actual_pressure_window_transport": packet,
        "partial_window_rows": copy.deepcopy(packet["partial_window_rows"]),
        "source_namespace": copy.deepcopy(packet["source_namespace"]),
        "axial_domain": copy.deepcopy(packet["axial_domain"]),
        "source_partition": copy.deepcopy(packet["source_partition"]),
        "source_scale_proof": copy.deepcopy(packet["source_scale_proof"]),
        "reference_length_binding": copy.deepcopy(setup["reference_binding"]),
        "authenticated_source_reports": {
            ADOPTION: {"bytes": ADOPTION_BYTES, "sha256": ADOPTION_SHA},
            FINITE_N: {"bytes": FINITE_N_BYTES, "sha256": FINITE_N_SHA},
            NATIVE: {"bytes": NATIVE_BYTES, "sha256": NATIVE_SHA},
            BACKGROUND: {"bytes": BACKGROUND_BYTES, "sha256": BACKGROUND_SHA},
            FIXTURE: {"bytes": FIXTURE_BYTES, "sha256": FIXTURE_SHA},
        },
        "execution_audit": {
            "accepted_ancestors_replayed": False,
            "source_execution": False,
            "new_saved_data_production_attempts": 1,
            "checkpoint_sha256": sha(CHECKPOINT),
            "computation_runner_sha256": checkpoint["computation_runner_sha256"],
            "kernel_sha256": checkpoint["kernel_sha256"],
        },
        "source_only": False,
        "accepted_ancestors_replayed": False,
        "execution_seconds": checkpoint["execution_seconds"],
        "interval_precision_bits": 1600,
    })
    scope = dict(adoption.get("scope", {}))
    scope.update({
        "pressure_window_transport_partial_group_computed": True,
        "partial_pressure_root_group_numeric": True,
        "new_values_are_directed_function_enclosures": True,
        "point_centers_selected": False,
    })
    for key in FALSE_SCOPE:
        scope[key] = False
    report["scope"] = scope
    report["limits"] = list(adoption.get("limits", [])) + [
        "Partial pressure-root group only; the seven remaining signed atoms in each row remain opaque.",
        "The four chart S1 rows are consumed as saved directed enclosures; no source or accepted parent math is replayed.",
        "The coefficient -1 external Rm source remains a directed [-2^-1024,+2^-1024] bound, never a selected center or zero.",
        "The partial identity is adopted_root/Rm^2 = opaque_prefix/Rm^2 + partial_window_Rm_normalized_C1 with frozen geometry.",
        "Macro current-N I1_full and D_full sources remain passed; mixed B_outer/S1 and three RQy boundary terms remain outside this partial transport.",
    ]
    return report


def package() -> None:
    audit = preflight()
    setup = _setup()
    report = _package(setup, audit)
    _write_json(REPORT, report)
    raw = REPORT.read_bytes()
    print(json.dumps({"report": str(REPORT), "bytes": len(raw), "sha256": digest(raw)}))


def inspect() -> None:
    audit = preflight()
    setup = _setup()
    require(REPORT.exists(), "Pressure window report is missing")
    require(CHECKPOINT.exists(), "Pressure window checkpoint is missing")
    checkpoint_raw = CHECKPOINT.read_bytes()
    checkpoint = json.loads(checkpoint_raw)
    require(isinstance(checkpoint, dict), "Checkpoint is not an object")
    count = _validate_checkpoint(checkpoint, setup, audit)
    raw = REPORT.read_bytes()
    report = json.loads(raw)
    require(isinstance(report, dict), "Report is not an object")
    expected_inputs = _expected_inputs(setup["adoption"])
    expected_implementations = _expected_implementations(setup["adoption"])
    require(report.get("input_hashes") == expected_inputs
            and report.get("implementation_hashes") == expected_implementations,
            "Report source pins changed")
    require(report.get("actual_post_restore_logu_C1") == setup["adoption"].get("actual_post_restore_logu_C1"),
            "Inherited post-restore packet changed")
    require(report.get("actual_restoration_pressure_adoption") == setup["adoption"].get("actual_restoration_pressure_adoption"),
            "Inherited adoption packet changed")
    require(report.get("preserved") == setup["adoption"].get("preserved"),
            "Inherited preservation ledger changed")
    require(report.get("actual_pressure_window_transport") == checkpoint.get("payload"),
            "Checkpoint packet does not match report")
    require(report.get("partial_window_rows") == checkpoint["payload"].get("partial_window_rows")
            and report.get("source_namespace") == checkpoint["payload"].get("source_namespace")
            and report.get("axial_domain") == checkpoint["payload"].get("axial_domain"),
            "Report partial rows/source namespace/domain changed")
    require(report.get("source_partition") == checkpoint["payload"]["source_partition"]
            and report.get("source_scale_proof") == checkpoint["payload"]["source_scale_proof"],
            "Report source partition/scale proof changed")
    execution = report.get("execution_audit")
    require(isinstance(execution, dict)
            and execution.get("checkpoint_sha256") == digest(checkpoint_raw)
            and execution.get("accepted_ancestors_replayed") is False
            and execution.get("source_execution") is False,
            "Report execution audit changed")
    require(report.get("source_only") is False
            and report.get("accepted_ancestors_replayed") is False
            and report.get("interval_precision_bits") == 1600
            and report.get("execution_seconds") == checkpoint.get("execution_seconds"),
            "Report execution metadata changed")
    scope = report.get("scope")
    require(isinstance(scope, dict) and all(scope.get(key) is False for key in FALSE_SCOPE),
            "Report scope flags overclaimed")
    require(scope.get("pressure_window_transport_partial_group_computed") is True
            and scope.get("partial_pressure_root_group_numeric") is True
            and scope.get("point_centers_selected") is False,
            "Report partial scope flags changed")
    print(json.dumps({
        "passed": True,
        "report": REPORT_NAME,
        "sha256": digest(raw),
        "bytes": len(raw),
        "canonical_exact_tuple_records": count,
        "input_pins": len(expected_inputs),
        "implementation_pins": len(expected_implementations),
        "inspection_math_providers_imported": False,
        "accepted_ancestors_replayed": False,
        "source_execution": False,
    }))


def produce() -> None:
    audit = preflight()
    setup = _setup()
    checkpoint = _compute_new(setup, audit)
    print(json.dumps({"checkpoint": str(CHECKPOINT), "sha256": sha(CHECKPOINT),
                      "kernel_sha256": checkpoint["kernel_sha256"]}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    actions = parser.add_mutually_exclusive_group(required=True)
    actions.add_argument("--compute-new", action="store_true")
    actions.add_argument("--write-report", action="store_true")
    actions.add_argument("--inspect-report", action="store_true")
    args = parser.parse_args()
    if args.compute_new:
        produce()
    elif args.write_report:
        package()
    else:
        inspect()
