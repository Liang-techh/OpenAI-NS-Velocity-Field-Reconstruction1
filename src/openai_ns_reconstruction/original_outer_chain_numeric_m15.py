"""Contract the six saved outer finite-N viscous joins.

The outer-chain adapter consumes the checked m15 packets and the checked
mixed-source graph archive.  It does not run an outer producer or replay a
constructor.  The returned packet removes the six consumed outer viscous
sectors from the correlated view and replaces the saved long--Rm external
term with one directed positive-union bound for the complete chain.
"""

from __future__ import annotations

import copy
import gzip
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from .original_core_numeric import PREFIX, ROOT
from .original_micro_viscous_m15 import contract_trace_chain


SAVED_Z = ("-.3", "0", ".3")
RAW_COUNT = 213
PREVIOUS_REDUCED_COUNT = 196
NEW_REDUCED_COUNT = 190
UNCHANGED_COUNT = 189

OLD_EXTERNAL_TERM = "long_Rm_Rh_finite_N_viscous_external_endpoints"
NEW_EXTERNAL_TERM = "long_outer_chain_finite_N_viscous_external_endpoints"

OUTER_M15_REPORT = PREFIX + "compliant_current_original_outer_m15_numeric.json"
AXIAL_M15_REPORT = PREFIX + "compliant_current_original_axial_m15_numeric.json"
POST_AXIAL_M15_REPORT = PREFIX + "compliant_current_original_post_axial_m15_numeric.json"
MIXED_SOURCE_ARCHIVE = PREFIX + "compliant_current_original_Rm_outer_mixed_source_joins.json.gz"
MIXED_SOURCE_CHECK = PREFIX + "compliant_current_original_Rm_outer_mixed_source_joins_check.json"
TRUE_OUTER_ARCHIVE = PREFIX + "compliant_current_original_true_outer_C6_C5_functions.json.gz"
TRUE_OUTER_CHECK = PREFIX + "compliant_current_original_true_outer_C6_C5_functions_check.json"

REPORT_SHA256 = {
    OUTER_M15_REPORT: "b8dd2f642c4c55ac183a087ed712e6a95a46aaaa7828cb4c70326b44deda8bc8",
    AXIAL_M15_REPORT: "57f5495394ab864ace99905fac07aab5bd315ff420329704121de121312b8d4a",
    POST_AXIAL_M15_REPORT: "48272c642ed8a28248aab35e7feda450f6a10aeff534f77b29f4b56a07f913f4",
    MIXED_SOURCE_ARCHIVE: "27510c39da6e84293c2eaa2355ce36d1c5622c39839c33a3b907739b85580a8c",
    MIXED_SOURCE_CHECK: "14d1cc6b641435c60a7e3662c996a260c07b328df5b8f486c4344f3d7eeb4867",
    TRUE_OUTER_ARCHIVE: "75895e1637411cfd2e490778f495f3adf904295a8840184cdc467e006b23305d",
    TRUE_OUTER_CHECK: "ce5669d01621683c3be857917f206b50af067b44feefe4b000c2ba99b2785bba",
}

JOIN_EXPECTED = {
    "patch_reference": {"left_chart": "Rm_patch", "right_chart": "Rh_reference", "left_query": 2948, "right_query": 3348},
    "reference_slope": {"left_chart": "Rh_reference", "right_chart": "O2_slope", "left_query": 0, "right_query": 0},
    "slope_axial": {"left_chart": "O2_slope", "right_chart": "O2_axial", "left_query": 1, "right_query": 0},
    "axial_buffer": {"left_chart": "O2_axial", "right_chart": "O2_buffer", "left_query": 1, "right_query": 0},
    "buffer_transition": {"left_chart": "O2_buffer", "right_chart": "O3_transition", "left_query": 4, "right_query": 0},
    "transition_power": {"left_chart": "O3_transition", "right_chart": "O3_power", "left_query": 1, "right_query": 0},
}

PACKET_SPECS = {
    "patch_reference": {
        "report": OUTER_M15_REPORT,
        "chart": "Rh_reference",
        "coordinate": "-5",
        "prefix_count": 168,
        "term_count": 172,
        "raw_index": 171,
        "term_name": "Rh_reference_finite_N_viscous_error",
        "aggregate_coordinate": "0",
    },
    "reference_slope": {
        "report": OUTER_M15_REPORT,
        "chart": "O2_slope",
        "coordinate": "0",
        "prefix_count": 172,
        "term_count": 176,
        "raw_index": 175,
        "term_name": "O2_slope_finite_N_viscous_error",
        "aggregate_coordinate": "1",
    },
    "slope_axial": {
        "report": AXIAL_M15_REPORT,
        "chart": "O2_axial",
        "coordinate": "0",
        "prefix_count": 176,
        "term_count": 195,
        "raw_index": 194,
        "term_name": "O2_axial_finite_N_viscous_error",
        "aggregate_coordinate": "1",
    },
    "axial_buffer": {
        "report": POST_AXIAL_M15_REPORT,
        "chart": "O2_buffer",
        "coordinate": "0",
        "prefix_count": 195,
        "term_count": 199,
        "raw_index": 198,
        "term_name": "O2_buffer_finite_N_viscous_error",
        "aggregate_coordinate": "11",
    },
    "buffer_transition": {
        "report": POST_AXIAL_M15_REPORT,
        "chart": "O3_transition",
        "coordinate": "0",
        "prefix_count": 199,
        "term_count": 203,
        "raw_index": 202,
        "term_name": "O3_transition_finite_N_viscous_error",
        "aggregate_coordinate": "1",
    },
    "transition_power": {
        "report": POST_AXIAL_M15_REPORT,
        "chart": "O3_power",
        "coordinate": "0",
        "prefix_count": 203,
        "term_count": 207,
        "raw_index": 206,
        "term_name": "O3_power_finite_N_viscous_error",
        "aggregate_coordinate": "2",
    },
}

REMOVED_OUTER_TERMS = tuple(spec["term_name"] for spec in PACKET_SPECS.values())


def _canonical(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def _digest(value: Any) -> str:
    return hashlib.sha256(_canonical(value).encode("utf-8")).hexdigest()


def _mapping(value: Any, label: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{label} must be a mapping")
    return value


def _find(rows: Any, *, label: str, **criteria: Any) -> Mapping[str, Any]:
    if not isinstance(rows, (list, tuple)):
        raise ValueError(f"Saved {label} rows are missing")
    for row in rows:
        if isinstance(row, Mapping) and all(row.get(key) == value for key, value in criteria.items()):
            return row
    raise ValueError(f"Saved {label} not found: {criteria}")


def _path(root: Any, relative: str) -> Path:
    candidates: list[Path] = []
    if root is not None:
        base = Path(root)
        candidates.extend((base / relative, base / "experiments" / "root_st073" / relative))
    candidates.extend((ROOT / relative, ROOT / "experiments" / "root_st073" / relative))
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    return candidates[0] if candidates else ROOT / relative


def _pin_file(owner: Any, relative: str, expected: str) -> dict[str, Any]:
    path = _path(getattr(owner, "directory", None), relative)
    actual = hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None
    input_hashes = _mapping(getattr(owner, "hashes", {}), "owner.hashes")
    code_hashes = _mapping(getattr(owner, "code_hashes", {}), "owner.code_hashes")
    basename = Path(relative).name
    for key, known in ((relative, code_hashes.get(relative)), (basename, input_hashes.get(basename))):
        if known is not None and known != actual:
            raise ValueError(f"Changed accepted source: {key}")
    if expected is not None and actual != expected:
        raise ValueError(f"Unexpected accepted source: {relative}")
    return {"path": relative, "exists": path.is_file(), "sha256": actual, "expected_sha256": expected, "passed": actual == expected}


def _load_json(owner: Any, relative: str, expected: str) -> Mapping[str, Any]:
    _pin_file(owner, relative, expected)
    path = _path(getattr(owner, "directory", None), relative)
    value = json.loads(path.read_text(encoding="utf-8"))
    return _mapping(value, relative)


def _load_gzip_json(owner: Any, relative: str, expected: str) -> Mapping[str, Any]:
    _pin_file(owner, relative, expected)
    path = _path(getattr(owner, "directory", None), relative)
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        value = json.load(handle)
    return _mapping(value, relative)


def _parent_report(owner: Any) -> Mapping[str, Any]:
    for candidate in (getattr(owner, "parent", None), getattr(owner, "parent_report", None), getattr(owner, "report", None)):
        if isinstance(candidate, Mapping) and isinstance(candidate.get("samples"), list):
            return candidate
        report = getattr(candidate, "report", None)
        if isinstance(report, Mapping) and isinstance(report.get("samples"), list):
            return report
    relative = PREFIX + "compliant_current_original_Rm_inlet_m15_numeric.json"
    path = _path(getattr(owner, "directory", None), relative)
    if not path.is_file():
        raise ValueError("Latest Rm inlet parent report required")
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, Mapping) or not isinstance(value.get("samples"), list):
        raise ValueError("Latest Rm inlet parent samples required")
    return value


def _selected_parent(parent_sample: Mapping[str, Any]) -> Mapping[str, Any]:
    selected = parent_sample.get("selected_N_evidence")
    if isinstance(selected, Mapping):
        return selected
    raise ValueError("Parent selected-N evidence required")


def _report_sample(report: Mapping[str, Any], Z: str, chart: str, coordinate: str) -> Mapping[str, Any]:
    return _find(report.get("samples"), label=f"{chart} m15 sample", Z=Z, chart=chart, native_coordinate=_coerce_coordinate(coordinate))


def _coerce_coordinate(value: Any) -> Any:
    # Reports serialize exact coordinates as strings.  Keep an exact string
    # when matching, while accepting an integer/float supplied by a caller.
    if isinstance(value, str):
        return value
    if isinstance(value, int):
        return str(value)
    return value


def _term_is_zero(term: Mapping[str, Any]) -> bool:
    coefficients = term.get("coefficient_C0_C1")
    if not isinstance(coefficients, list) or len(coefficients) != 2:
        return False
    for row in coefficients:
        if not isinstance(row, Mapping):
            return False
        if row.get("approximation") not in ("0.0", "0") or row.get("absolute_error_upper") not in ("0.0", "0"):
            return False
        if row.get("raw_directed_error_exact_mpf_tuple") != [0, 0, 0, 0]:
            return False
        if row.get("source_enclosure_exact_mpf_tuples") != [[0, 0, 0, 0], [0, 0, 0, 0]]:
            return False
    return True


def _term_names(rows: Any, label: str) -> list[str]:
    if not isinstance(rows, list):
        raise ValueError(f"{label} must be a list")
    names = [row.get("name") if isinstance(row, Mapping) else None for row in rows]
    if any(name is None for name in names):
        raise ValueError(f"{label} has an unnamed term")
    return names


def _archive_endpoint(value: Mapping[str, Any]) -> dict[str, Any]:
    physical = value.get("physical")
    if not isinstance(physical, Mapping):
        raise ValueError("Outer join endpoint physical graph rows required")
    intermediate = value.get("intermediate")
    compact_intermediate: dict[str, Any] = {}
    if isinstance(intermediate, Mapping):
        for key in ("coordinate", "domain", "radius", "fixed_pressure", "region", "complete_incoming_mass_not_reset"):
            if key in intermediate:
                compact_intermediate[key] = copy.deepcopy(intermediate[key])
        for key in ("fields", "cumulative"):
            if key in intermediate:
                compact_intermediate[key + "_sha256"] = _digest(intermediate[key])
    return {
        "native_query": copy.deepcopy(value.get("native_query")),
        "radius": copy.deepcopy(value.get("radius")),
        "physical": copy.deepcopy(dict(physical)),
        "physical_sha256": _digest(physical),
        "N_y": copy.deepcopy(value.get("N_y")),
        "intermediate": compact_intermediate,
        "own_sha256": _digest(value.get("own")) if isinstance(value.get("own"), Mapping) else None,
        "numerical_oracle_available": value.get("numerical_oracle_available"),
    }


def _load_join_archive(owner: Any) -> tuple[dict[str, Any], dict[str, Any]]:
    archive = _load_gzip_json(owner, MIXED_SOURCE_ARCHIVE, REPORT_SHA256[MIXED_SOURCE_ARCHIVE])
    binding = _mapping(archive.get("Rm_outer_mixed_source_join_bindings"), "mixed outer join bindings")
    charts = _mapping(binding.get("charts"), "mixed outer charts")
    joins = _mapping(binding.get("joins"), "mixed outer joins")
    compact_joins: dict[str, Any] = {}
    for name, expected in JOIN_EXPECTED.items():
        join = _mapping(joins.get(name), f"outer join {name}")
        if join.get("left_chart") != expected["left_chart"] or join.get("right_chart") != expected["right_chart"]:
            raise ValueError(f"Outer join chart identity changed: {name}")
        left = _mapping(join.get("left"), f"outer join {name} left")
        right = _mapping(join.get("right"), f"outer join {name} right")
        if left.get("native_query") != expected["left_query"] or right.get("native_query") != expected["right_query"]:
            raise ValueError(f"Outer join native query changed: {name}")
        compact_joins[name] = {
            "left_chart": expected["left_chart"],
            "right_chart": expected["right_chart"],
            "left": _archive_endpoint(left),
            "right": _archive_endpoint(right),
        }
    compact_charts: dict[str, Any] = {}
    for name in ("Rh_reference", "O2_slope", "O2_axial", "O2_buffer", "O3_transition", "O3_power"):
        chart = _mapping(charts.get(name), f"outer chart {name}")
        compact_charts[name] = {
            "domain": copy.deepcopy(chart.get("domain")),
            "coordinate": chart.get("coordinate"),
            "radius": chart.get("radius"),
            "jacobian": chart.get("jacobian"),
            "phase": chart.get("phase"),
            "physical": copy.deepcopy(chart.get("physical")),
            "physical_sha256": _digest(chart.get("physical")),
            "N_y": copy.deepcopy(chart.get("N_y")),
            "intermediate": {
                key: copy.deepcopy(chart["intermediate"].get(key))
                for key in ("coordinate", "domain", "radius", "fixed_pressure", "region", "complete_incoming_mass_not_reset")
                if isinstance(chart.get("intermediate"), Mapping) and key in chart["intermediate"]
            },
            "own_sha256": _digest(chart.get("own")) if isinstance(chart.get("own"), Mapping) else None,
        }
    source = {
        "source_family": copy.deepcopy(archive.get("source_family")),
        "definition": copy.deepcopy(binding.get("definition")),
        "charts": compact_charts,
        "joins": compact_joins,
        "flat_controls": copy.deepcopy(binding.get("flat_controls")),
        "contract": copy.deepcopy(binding.get("contract")),
        "archive_sha256": REPORT_SHA256[MIXED_SOURCE_ARCHIVE],
    }
    check = _load_json(owner, MIXED_SOURCE_CHECK, REPORT_SHA256[MIXED_SOURCE_CHECK])
    if check.get("all_passed") is not True:
        raise ValueError("Mixed outer source join receipt must pass")
    for name in JOIN_EXPECTED:
        flag = f"current_original_n1_{name}_full_mixed_source_join_installed"
        if check.get(flag) is not True:
            raise ValueError(f"Mixed outer source receipt missing {name}")
    source["receipt"] = {
        "all_passed": True,
        "source_family": copy.deepcopy(check.get("source_family")),
        "join_flags": {name: True for name in JOIN_EXPECTED},
        "receipt_sha256": REPORT_SHA256[MIXED_SOURCE_CHECK],
    }
    true_archive = _load_gzip_json(owner, TRUE_OUTER_ARCHIVE, REPORT_SHA256[TRUE_OUTER_ARCHIVE])
    true_check = _load_json(owner, TRUE_OUTER_CHECK, REPORT_SHA256[TRUE_OUTER_CHECK])
    if true_check.get("all_passed") is not True:
        raise ValueError("True outer C6/C5 receipt must pass")
    source["true_outer_archive"] = {
        "source_family": copy.deepcopy(true_archive.get("source_family")),
        "archive_sha256": REPORT_SHA256[TRUE_OUTER_ARCHIVE],
        "receipt_sha256": REPORT_SHA256[TRUE_OUTER_CHECK],
        "all_passed": True,
    }
    # Keep the returned binding compact.  The full gzip archives are pinned
    # and inspected above, but carrying their multi-megabyte graph payloads
    # into every sample would obscure the numeric packet and invite accidental
    # use as a producer input.
    return source, {
        "mixed_archive": {"path": MIXED_SOURCE_ARCHIVE, "sha256": REPORT_SHA256[MIXED_SOURCE_ARCHIVE], "passed": True},
        "mixed_check": {"path": MIXED_SOURCE_CHECK, "sha256": REPORT_SHA256[MIXED_SOURCE_CHECK], "passed": True},
        "true_archive": {"path": TRUE_OUTER_ARCHIVE, "sha256": REPORT_SHA256[TRUE_OUTER_ARCHIVE], "passed": True},
        "true_check": {"path": TRUE_OUTER_CHECK, "sha256": REPORT_SHA256[TRUE_OUTER_CHECK], "passed": True},
    }


def _proofs(seam_proofs: Mapping[str, Any]) -> dict[str, Mapping[str, Any]]:
    candidate = seam_proofs.get("joins") if isinstance(seam_proofs.get("joins"), Mapping) else seam_proofs
    result: dict[str, Mapping[str, Any]] = {}
    for name, expected in JOIN_EXPECTED.items():
        value = candidate.get(name) if isinstance(candidate, Mapping) else None
        if not isinstance(value, Mapping) or value.get("passed") is not True:
            raise ValueError(f"Accepted {name} source proof must pass")
        for key, expected_value in (("left_chart", expected["left_chart"]), ("right_chart", expected["right_chart"])):
            if key in value and value.get(key) != expected_value:
                raise ValueError(f"{name} source proof has the wrong {key}")
        result[name] = value
    return result


def _power_certificate(seam_proofs: Mapping[str, Any], owner: Any = None) -> tuple[bool, Mapping[str, Any] | None]:
    explicit = False
    for key in (
        "power_transport_certificate",
        "native_power_transport_certificate",
        "O3_power_transport_certificate",
        "power_transport",
        "O3_power_constant_D",
    ):
        value = seam_proofs.get(key)
        if value is not None:
            explicit = True
        if value is None and isinstance(seam_proofs.get("joins"), Mapping):
            value = seam_proofs["joins"].get(key)
            if value is not None:
                explicit = True
        if value is None:
            continue
        if not isinstance(value, Mapping) or value.get("passed") is not True:
            raise ValueError("Explicit O3 power transport certificate did not pass")
        flags = []
        for name, flag in value.items():
            normalized = str(name).lower().replace("_", "")
            if flag is True and any(token in normalized for token in ("constant", "equal", "identical", "lowerupper", "dpower")):
                flags.append(name)
        if not flags:
            raise ValueError("O3 power certificate must explicitly prove equal lower and upper D")
        return True, value
    # A native owner may expose the same accepted certificate separately.  A
    # failed native candidate is not a request to use the stronger bound, so
    # retain the safe O3-power aggregate unless the candidate itself passed.
    candidate = getattr(owner, "power_transport", None)
    if not explicit and isinstance(candidate, Mapping) and candidate.get("passed") is True:
        flags = [
            name for name, flag in candidate.items()
            if flag is True and any(token in str(name).lower().replace("_", "") for token in ("constant", "equal", "identical", "lowerupper", "dpower"))
        ]
        if flags:
            return True, candidate
        raise ValueError("O3 power certificate must explicitly prove equal lower and upper D")
    return False, None


def _validate_packets(
    owner: Any,
    parent: Mapping[str, Any],
    parent_sample: Mapping[str, Any],
) -> tuple[dict[str, Any], dict[str, Any]]:
    reports = {
        OUTER_M15_REPORT: _load_json(owner, OUTER_M15_REPORT, REPORT_SHA256[OUTER_M15_REPORT]),
        AXIAL_M15_REPORT: _load_json(owner, AXIAL_M15_REPORT, REPORT_SHA256[AXIAL_M15_REPORT]),
        POST_AXIAL_M15_REPORT: _load_json(owner, POST_AXIAL_M15_REPORT, REPORT_SHA256[POST_AXIAL_M15_REPORT]),
    }
    raw = parent_sample.get("raw_m15_C0_C1_terms")
    if not isinstance(raw, list) or len(raw) != RAW_COUNT:
        raise ValueError("Raw 213 m15 state must remain verbatim")
    _term_names(raw, "parent raw m15 terms")
    # The complete upper packet of one chart is the exact incoming prefix of
    # the next chart.  Keep the dictionaries, rather than only their names,
    # so serialized coefficients, source bindings, and log atoms are pinned
    # by the seam check as well.
    previous_upper_terms: list[Mapping[str, Any]] | None = None
    expected_family = getattr(owner, "family", None)
    selected_reference: Mapping[str, Any] | None = None
    packets: dict[str, Any] = {}
    for join_name, spec in PACKET_SPECS.items():
        report = reports[spec["report"]]
        if expected_family is not None and report.get("source_family") != expected_family:
            raise ValueError(f"{join_name} source family changed")
        endpoint = _report_sample(report, str(parent_sample.get("Z")), spec["chart"], spec["coordinate"])
        terms = endpoint.get("m15_C0_C1_terms")
        incoming = endpoint.get("incoming_m15_terms")
        if not isinstance(terms, list) or len(terms) != spec["term_count"]:
            raise ValueError(f"{join_name} complete endpoint packet count changed")
        if not isinstance(incoming, list) or len(incoming) != spec["prefix_count"]:
            raise ValueError(f"{join_name} incoming prefix count changed")
        names = _term_names(terms, f"{join_name} terms")
        incoming_names = _term_names(incoming, f"{join_name} incoming terms")
        parent_raw_prefix = raw[: spec["prefix_count"]]
        terms_prefix = terms[: spec["prefix_count"]]
        lower_terms_prefix_equals_incoming = terms_prefix == incoming
        lower_incoming_equals_parent_raw = incoming == parent_raw_prefix
        if not lower_terms_prefix_equals_incoming or not lower_incoming_equals_parent_raw:
            raise ValueError(f"{join_name} incoming m15 prefix is not the exact raw cumulative prefix")
        previous_upper_terms_equals_incoming = previous_upper_terms is None or previous_upper_terms == incoming
        if not previous_upper_terms_equals_incoming:
            raise ValueError(f"{join_name} incoming prefix does not equal the previous chart upper packet")
        if names[spec["raw_index"]] != spec["term_name"]:
            raise ValueError(f"{join_name} viscous raw index changed")
        if not endpoint.get("exact_empty_increment"):
            raise ValueError(f"{join_name} lower endpoint must carry an exact empty increment")
        tails = terms[spec["prefix_count"] :]
        if not tails or not all(_term_is_zero(term) for term in tails):
            raise ValueError(f"{join_name} lower endpoint tails must be exact zero packets")
        selected = endpoint.get("actual_selected_integer")
        if not isinstance(selected, Mapping):
            raise ValueError(f"{join_name} selected integer is missing")
        if selected_reference is None:
            selected_reference = selected
        elif selected != selected_reference:
            raise ValueError("Outer m15 packets do not use one selected N")
        aggregate_endpoint = _report_sample(report, str(parent_sample.get("Z")), spec["chart"], spec["aggregate_coordinate"])
        aggregate_terms = aggregate_endpoint.get("m15_C0_C1_terms")
        aggregate_incoming = aggregate_endpoint.get("incoming_m15_terms")
        if not isinstance(aggregate_terms, list) or len(aggregate_terms) != spec["term_count"]:
            raise ValueError(f"{join_name} upper aggregate packet count changed")
        if not isinstance(aggregate_incoming, list) or len(aggregate_incoming) != spec["prefix_count"]:
            raise ValueError(f"{join_name} upper aggregate incoming prefix count changed")
        aggregate_names = _term_names(aggregate_terms, f"{join_name} aggregate terms")
        _term_names(aggregate_incoming, f"{join_name} aggregate incoming terms")
        aggregate_terms_prefix = aggregate_terms[: spec["prefix_count"]]
        upper_terms_prefix_equals_incoming = aggregate_terms_prefix == aggregate_incoming
        upper_incoming_equals_parent_raw = aggregate_incoming == parent_raw_prefix
        if not upper_terms_prefix_equals_incoming or not upper_incoming_equals_parent_raw:
            raise ValueError(f"{join_name} upper aggregate incoming prefix changed")
        if aggregate_names[spec["raw_index"]] != spec["term_name"]:
            raise ValueError(f"{join_name} upper aggregate viscous raw index changed")
        aggregate = _find(aggregate_endpoint.get("m15_C0_C1_terms"), label=f"{join_name} viscous aggregate", name=spec["term_name"])
        logs = aggregate.get("log_scale_terms")
        if not isinstance(logs, list) or len(logs) != 1:
            raise ValueError(f"{join_name} complete aggregate must have one positive log atom")
        atom = logs[0]
        if atom.get("positive_error_bound_atom_only") is not True or atom.get("source_value_selected") is not False:
            raise ValueError(f"{join_name} aggregate must remain an error-only positive atom")
        if aggregate_endpoint.get("actual_selected_integer") != selected_reference:
            raise ValueError(f"{join_name} aggregate selected N changed")
        packets[join_name] = {
            "join": join_name,
            "chart": spec["chart"],
            "endpoint_coordinate": spec["coordinate"],
            "aggregate_coordinate": spec["aggregate_coordinate"],
            "prefix_count": spec["prefix_count"],
            "term_count": spec["term_count"],
            "raw_index": spec["raw_index"],
            "term_name": spec["term_name"],
            "endpoint_sample_sha256": _digest(endpoint),
            "aggregate_sample_sha256": _digest(aggregate_endpoint),
            "incoming_prefix_names_sha256": _digest(incoming_names),
            "exact_prefix_witness": {
                "passed": True,
                "parent_raw_prefix_count": len(parent_raw_prefix),
                "lower_incoming_count": len(incoming),
                "lower_terms_prefix_count": len(terms_prefix),
                "upper_incoming_count": len(aggregate_incoming),
                "upper_terms_prefix_count": len(aggregate_terms_prefix),
                "previous_upper_terms_count": None if previous_upper_terms is None else len(previous_upper_terms),
                "lower_incoming_equals_parent_raw": lower_incoming_equals_parent_raw,
                "lower_terms_prefix_equals_incoming": lower_terms_prefix_equals_incoming,
                "upper_incoming_equals_parent_raw": upper_incoming_equals_parent_raw,
                "upper_terms_prefix_equals_incoming": upper_terms_prefix_equals_incoming,
                "previous_upper_terms_equals_lower_incoming": previous_upper_terms_equals_incoming,
                "exact_dictionary_equality": True,
                "parent_raw_prefix_sha256": _digest(parent_raw_prefix),
                "lower_incoming_sha256": _digest(incoming),
                "lower_terms_prefix_sha256": _digest(terms_prefix),
                "upper_incoming_sha256": _digest(aggregate_incoming),
                "upper_terms_prefix_sha256": _digest(aggregate_terms_prefix),
                "previous_upper_terms_sha256": None if previous_upper_terms is None else _digest(previous_upper_terms),
            },
            "empty_tail_count": len(tails),
            "empty_tail_names": [term.get("name") for term in tails],
            "aggregate": copy.deepcopy(dict(aggregate)),
            "aggregate_log_atom": copy.deepcopy(dict(atom)),
            "selected_integer": copy.deepcopy(dict(selected)),
            "source_report": spec["report"],
        }
        # The next endpoint starts with this endpoint's complete upper packet;
        # the next iteration compares it by full dictionary equality.
        previous_upper_terms = copy.deepcopy(aggregate_terms)
    if selected_reference is None or selected_reference != _selected_parent(parent_sample):
        raise ValueError("Outer selected N must equal the parent selected-N evidence")
    return reports, packets


def _trace_chain() -> dict[str, Any]:
    pieces = [
        {"name": "long_reshape", "lower": "long_native0", "upper": "long_native1"},
        {"name": "reference", "lower": "reference_native0", "upper": "reference_native1"},
        {"name": "slope", "lower": "slope_native0", "upper": "slope_native1"},
        {"name": "axial", "lower": "axial_native0", "upper": "axial_native1"},
        {"name": "buffer", "lower": "buffer_native0", "upper": "buffer_native1"},
        {"name": "transition", "lower": "transition_native0", "upper": "transition_native1"},
        {"name": "power", "lower": "power_native0", "upper": "power_native1"},
    ]
    aliases = {
        "long_native1": "reference_native0",
        "reference_native1": "slope_native0",
        "slope_native1": "axial_native0",
        "axial_native1": "buffer_native0",
        "buffer_native1": "transition_native0",
        "transition_native1": "power_native0",
    }
    result = contract_trace_chain(pieces, aliases)
    result.update(
        {
            "external_endpoints": {"lower": "long native0", "upper": "O3 power/Rc"},
            "cancelled_interfaces_explicit": [
                "long native1 = reference native0",
                "reference native1 = slope native0",
                "slope native1 = axial native0",
                "axial native1 = buffer native0",
                "buffer native1 = transition native0",
                "transition native1 = power native0",
            ],
            "function": "D(long native0)-D(outer chain upper source)",
        }
    )
    return result


def _parent_long_witness(parent_sample: Mapping[str, Any]) -> Mapping[str, Any]:
    budget = _mapping(parent_sample.get("external_endpoint_budget_witness"), "parent external endpoint budget")
    witness = _mapping(budget.get("own_longreshape0"), "parent own long native0 witness")
    rows = witness.get("radius_weighted_C0_C1_log_bounds")
    if not isinstance(rows, list) or len(rows) != 2 or any(value is None for value in rows):
        raise ValueError("Parent own long native0 C0/Z1 positive bounds required")
    return witness


def _evaluate_sample(
    owner: Any,
    parent: Mapping[str, Any],
    source: Mapping[str, Any],
    seams: Mapping[str, Any],
    parent_sample: Mapping[str, Any],
    packets: Mapping[str, Any],
    power_certificate: Mapping[str, Any] | None,
    power_certificate_used: bool,
) -> dict[str, Any]:
    raw = parent_sample["raw_m15_C0_C1_terms"]
    previous = parent_sample.get("correlated_reduced_m15_C0_C1_terms")
    if not isinstance(previous, list) or len(previous) != PREVIOUS_REDUCED_COUNT:
        raise ValueError("Latest Rm inlet 196-term view required")
    old_external = [term for term in previous if term.get("name") == OLD_EXTERNAL_TERM]
    if len(old_external) != 1:
        raise ValueError("Latest Rm inlet external viscous replacement required")
    removed = [term for term in previous if term.get("name") in REMOVED_OUTER_TERMS]
    if len(removed) != len(REMOVED_OUTER_TERMS) or {term.get("name") for term in removed} != set(REMOVED_OUTER_TERMS):
        raise ValueError("All six outer viscous sectors must be present exactly once")

    aggregate_join = "buffer_transition" if power_certificate_used else "transition_power"
    aggregate_packet = packets[aggregate_join]
    aggregate_atom = aggregate_packet["aggregate_log_atom"]
    witness = _parent_long_witness(parent_sample)
    long_atoms = list(witness["radius_weighted_C0_C1_log_bounds"])
    positive_inputs = [owner.decode(atom) for atom in long_atoms]
    positive_inputs.append(owner.decode(aggregate_atom))
    bound = owner.positive_union(positive_inputs)

    replacement = copy.deepcopy(old_external[0])
    replacement["name"] = NEW_EXTERNAL_TERM
    replacement["log_scale_terms"] = [owner.record(bound)]
    replacement["actual_function_binding"] = {
        "error_enclosure_only": True,
        "function": "D(long native0)-D(outer chain upper source)",
        "delta_definition": "R*(Qy_modified-Qy_leading)",
        "ordinary_Z1": True,
        "physical_R_once": True,
        "external_endpoints_not_zeroed": True,
        "signed_centers_evaluated": False,
        "no_error_cap_subtraction": True,
        "no_interval_midpoint_or_overlap_closure": True,
        "directed_positive_union_long_native0_then_outer_aggregate": True,
        "parent_long_native0_budget_retained": True,
        "outer_aggregate_dominates_upper_chain": True,
        "outer_aggregate_chart": aggregate_packet["chart"],
        "outer_aggregate_raw_index": aggregate_packet["raw_index"],
        "outer_aggregate_coordinate": aggregate_packet["aggregate_coordinate"],
        "power_transport_certificate_used": power_certificate_used,
        "power_transport_certificate_required_for_transition_202": True,
        "source_functions_values_or_caps_not_selected": True,
        "six_source_joins_proved": True,
        "source_seam_proofs": copy.deepcopy(dict(seams)),
    }

    reduced: list[dict[str, Any]] = []
    removed_names: list[str] = []
    for term in previous:
        name = term.get("name")
        if name in REMOVED_OUTER_TERMS:
            removed_names.append(name)
            continue
        if name == OLD_EXTERNAL_TERM:
            reduced.append(copy.deepcopy(replacement))
        else:
            reduced.append(copy.deepcopy(term))
    if sorted(removed_names) != sorted(REMOVED_OUTER_TERMS) or len(reduced) != NEW_REDUCED_COUNT:
        raise ValueError("Exact new 190-term view required")
    unchanged_old = [term for term in previous if term.get("name") not in set(REMOVED_OUTER_TERMS) | {OLD_EXTERNAL_TERM}]
    unchanged_new = [term for term in reduced if term.get("name") != NEW_EXTERNAL_TERM]
    if len(unchanged_old) != UNCHANGED_COUNT or unchanged_new != unchanged_old:
        raise ValueError("Other 189 correlated sectors must remain verbatim")
    if sum(term.get("name") == NEW_EXTERNAL_TERM for term in reduced) != 1:
        raise ValueError("Exactly one new outer-chain replacement required")

    target = copy.deepcopy(parent_sample.get("independent_target"))
    if target is None:
        target = copy.deepcopy(parent.get("independent_target"))
    if target is None:
        raise ValueError("Independent target must remain available")
    removed_provenance = [copy.deepcopy(raw[spec["raw_index"]]) for spec in PACKET_SPECS.values()]
    return {
        "Z": parent_sample.get("Z"),
        "source_family": getattr(owner, "family", parent.get("source_family")),
        "raw_m15_C0_C1_terms": copy.deepcopy(raw),
        "correlated_reduced_m15_C0_C1_terms": reduced,
        "original_raw_terms": RAW_COUNT,
        "previous_reduced_terms": PREVIOUS_REDUCED_COUNT,
        "correlated_reduced_terms": NEW_REDUCED_COUNT,
        "unchanged_other_sectors": UNCHANGED_COUNT,
        "signed_trace_chain": _trace_chain(),
        "source_seam_proofs": copy.deepcopy(dict(seams)),
        "external_endpoint_budget_witness": {
            "own_longreshape0": copy.deepcopy(dict(witness)),
            "outer_chain_upper_aggregate": copy.deepcopy(aggregate_packet),
            "positive_union": owner.record(bound),
            "positive_union_inputs": copy.deepcopy(long_atoms) + [copy.deepcopy(aggregate_atom)],
            "source_binding": copy.deepcopy(dict(source)),
            "power_transport_certificate": copy.deepcopy(dict(power_certificate)) if power_certificate is not None else None,
            "power_transport_certificate_used": power_certificate_used,
            "no_cap_subtraction_value_midpoint_or_overlap": True,
        },
        "old_removed_outer_viscous_terms": removed_provenance,
        "previous_external_replacement_retained_as_provenance": copy.deepcopy(old_external[0]),
        "replacement_outer_endpoints": copy.deepcopy(replacement),
        "independent_target": target,
        "selected_N_evidence": copy.deepcopy(parent_sample.get("selected_N_evidence")),
        "largest_remaining_individual_log_caps": owner.largest_caps(reduced, target),
        "contract": {
            "raw_213_verbatim": True,
            "other_189_verbatim": True,
            "independent_target_unchanged": True,
            "reserved_retained_once": sum(term.get("name") == "reserved_n1_own_swirl_once" for term in reduced) == 1,
            "ordinary_Z1_and_physical_R_once": True,
            "leading_subtraction_retained": True,
            "all_six_outer_source_joins_proved": True,
            "six_signed_outer_viscous_sectors_removed": True,
            "outer_upper_aggregate_selected_as_error_only": True,
            "power_transport_certificate_used": power_certificate_used,
            "power_transport_certificate_required_for_202": True,
            "flat_U_does_not_zero_history_mean_or_Qy": True,
            "no_error_cap_as_value": True,
            "all_other_volume_width_numerical_errors_retained": True,
            "ancestors_not_reexecuted": True,
            "source_bound_signed_trace_contraction": True,
        },
        "limits": {
            "outer_chain_internal_trace_resolved": True,
            "inherited_pressure_bounds_independently_recertified": False,
            "external_signed_values_resolved": False,
            "power_transport_identity_assumed_without_certificate": False,
            "independent_numeric_Rp_seam": False,
            "global_five_moment_common_N_closure": False,
            "exact_heat_and_energy": False,
            "admissible_stress_flat_remainder": False,
            "temporal_n_recursion_and_summation": False,
            "independent_corrected_Cartesian_NS": False,
        },
    }


def build_outer_chain_numeric(owner: Any, seam_proofs: Mapping[str, Any]) -> dict[str, Any]:
    """Return the saved six-join outer-chain 190-term contraction packet."""
    if not isinstance(seam_proofs, Mapping):
        raise ValueError("Six outer source seam proofs are required")
    proofs = _proofs(seam_proofs)
    power_used, power_certificate = _power_certificate(seam_proofs, owner)
    parent = _parent_report(owner)
    source, archive_receipts = _load_join_archive(owner)
    source.update(
        {
            "passed": True,
            "function": "D=R*(Qy_modified-Qy_leading); long native0 minus outer-chain upper source",
            "physical_radius_power": 1,
            "ordinary_Z_orders": [0, 1],
            "physical_R_once": True,
            "leading_subtraction_retained": True,
            "actual_outer_join_graph_binding": source.pop("joins"),
            "actual_outer_chart_graph_binding": source.pop("charts"),
            "source_seam_proofs": copy.deepcopy(proofs),
            "source_proof_sha256": {name: _digest(proof) for name, proof in proofs.items()},
            "source_archive_receipts": archive_receipts,
            "outer_m15_packet_reports": {
                relative: {"path": relative, "sha256": REPORT_SHA256[relative]}
                for relative in (OUTER_M15_REPORT, AXIAL_M15_REPORT, POST_AXIAL_M15_REPORT)
            },
            "raw_outer_viscous_indices": {
                spec["term_name"]: spec["raw_index"] for spec in PACKET_SPECS.values()
            },
            "raw_outer_viscous_names": list(REMOVED_OUTER_TERMS),
            "numeric_endpoint_coordinates": {
                name: {"chart": spec["chart"], "lower_coordinate": spec["coordinate"], "aggregate_coordinate": spec["aggregate_coordinate"]}
                for name, spec in PACKET_SPECS.items()
            },
            "default_outer_aggregate_raw_index": 206,
            "default_outer_aggregate_chart": "O3_power",
            "default_outer_aggregate_coordinate": "2",
            "power_transport_certificate_required_for_raw_202": True,
            "power_transport_certificate_used": power_used,
            "power_transport_certificate": copy.deepcopy(dict(power_certificate)) if power_certificate is not None else None,
            "flat_U_does_not_zero_history_mean_or_Qy": True,
            "parent_long_native0_budget_retained": True,
            "full_original_O3_power_viscous_aggregate_used_by_default": not power_used,
            "full_original_O3_transition_viscous_aggregate_used_with_power_certificate": power_used,
            "outer_chain_upper_aggregate_dominates": True,
            "all_other_189_correlated_sectors_retained_verbatim": True,
            "independent_target_retained": True,
            "reserved_n1_own_swirl_retained_once": True,
            "source_functions_values_or_caps_not_selected": True,
            "no_error_cap_as_function_value": True,
            "accepted_ancestor_constructors_not_run": True,
            "parent_report_sha256": _digest(parent),
        }
    )
    samples: list[dict[str, Any]] = []
    for Z in SAVED_Z:
        parent_sample = _find(parent.get("samples"), label="latest Rm inlet parent sample", Z=Z)
        _reports, packets = _validate_packets(owner, parent, parent_sample)
        samples.append(_evaluate_sample(owner, parent, source, proofs, parent_sample, packets, power_certificate, power_used))
    selected = [sample["selected_N_evidence"] for sample in samples]
    if not selected or any(value != selected[0] for value in selected[1:]):
        raise ValueError("All outer-chain samples must use one selected N")
    return {
        "source_binding": source,
        "samples": samples,
        "source_seam_proofs": {
            "passed": True,
            "joins": copy.deepcopy(proofs),
            "power_transport_certificate": copy.deepcopy(dict(power_certificate)) if power_certificate is not None else None,
            "power_transport_certificate_used": power_used,
        },
        "contract": {
            "raw_213_retained": True,
            "previous_196_to_new_190": True,
            "unchanged_other_189": True,
            "saved_packets_only": True,
            "accepted_ancestor_constructors_not_run": True,
            "cap_is_not_signed_function_value": True,
            "six_outer_source_joins_required": True,
            "all_six_source_sectors_removed_once": True,
            "power_transport_certificate_required_for_202": True,
            "flat_U_does_not_zero_history_mean_or_Qy": True,
        },
        "limits": {
            "numeric_Rp_closure": False,
            "global_all_five_closure": False,
            "full_goal_complete": False,
        },
    }


__all__ = ["build_outer_chain_numeric"]
