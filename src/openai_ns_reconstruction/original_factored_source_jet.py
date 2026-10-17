"""Fail-closed provenance guard for a rejected R100-to-R0 source adapter.

The saved switch packet and the requested macro source use different chart
generators.  This module reads only pinned JSON/gzip metadata and source ASTs,
records that origin mismatch, and returns no source rows.  The historical
factored source algebra was deliberately removed after audit; the original
draft is preserved separately for review.
"""

from __future__ import annotations

import ast
import gzip
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping


_FILES = {
    "bridge_switch": (
        "lei_ren_part1_paper_compliant_actual_bridge_switch.json",
        "1f7a80703d66be2fb44749efedb6aa03e63537c2f56ba7be9c75cf9833b6388c",
    ),
    "bridge_integrals": (
        "lei_ren_part1_paper_compliant_actual_bridge_integrals.json",
        "7e25fa846ef77e408dea54acf331df081e96b45a702ff3c9bc87136221144455",
    ),
    "first_collar": (
        "lei_ren_part1_paper_compliant_current_original_first_collar_jets.json",
        "a1a2818d18bf734c59460701d9ec4b342d6a8392ab97f581915766e9e6701a2f",
    ),
    "high_archive": (
        "lei_ren_part1_paper_compliant_current_original_prescribed_high_radial_source_functions.json.gz",
        "5b54a210d48077a984b72b3c74a182c628bf00fb4c2b71f8f164845c228de991",
    ),
    "parent_report": (
        "lei_ren_part1_paper_compliant_current_original_leading_viscous_source_numeric.json",
        "d1aa3644ea8158fe78dccb210a0a75a4c81ad257b8504758cb5797845b8003d4",
    ),
}

_SOURCE_FILES = {
    # These source pins are deliberately independent of the serialized bridge
    # packets.  The public builder is allowed to inspect their ASTs, but it
    # must not infer that two different chart generators are interchangeable.
    "switch_signed_integrals_ast": (
        "experiments/root_st073/lei_ren_part1_paper_compliant_switch_signed_integrals.py",
        "713960de53db501905a9ad79c181755811d1f1a1d82da624ef615cd6db20a85a",
    ),
    "macro_signed_integrals_ast": (
        "experiments/root_st073/lei_ren_part1_paper_compliant_macro_signed_integrals.py",
        "046b7982464eb00eb16c41c3a652c46a6ff82a1c7386fc3283dc5dd654377913",
    ),
    "macro_source_ast": (
        "src/openai_ns_reconstruction/original_macro_pressure_rows.py",
        "6a96ec977ae8f4235bdb874ce9d2c7b5acc15b91f55710fef1ec040e025d6753",
    ),
    "macro_numeric_ast": (
        "src/openai_ns_reconstruction/original_macro_numeric.py",
        "899b7aec6844f97693c67f8562702888794a8ad3c089a484dc8ab82c4b11afb1",
    ),
    "second_micro_numeric_ast": (
        "src/openai_ns_reconstruction/original_second_micro_numeric.py",
        "1c7f0d61ec2c9bcbc84454c2e047dd6c4382c3fa789be46f2339ffbe3180af2a",
    ),
    "bridge_integrals_ast": (
        "experiments/root_st073/lei_ren_part1_paper_compliant_actual_bridge_integrals.py",
        "2d3f5453b2838353d6eb8a77f6d719bf51cfda48c9649d67604603e9aed5d8a3",
    ),
    "bridge_switch_ast": (
        "experiments/root_st073/lei_ren_part1_paper_compliant_actual_bridge_switch.py",
        "f23b8bcb770401a595563d0650760eef1e6e143c6bfed020b16386b78f05c2a4",
    ),
}


def _mapping(value: Any) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def _repo_root(root: Any, directory: Any) -> Path:
    for candidate in (root, directory):
        if candidate is None:
            continue
        try:
            path = Path(candidate)
        except TypeError:
            continue
        if path.is_file():
            path = path.parent
        if (path / "src").is_dir() and (path / "experiments").is_dir():
            return path
    return Path(__file__).resolve().parents[2]


def _relative(path: Path, repo: Path) -> str:
    try:
        return path.relative_to(repo).as_posix()
    except ValueError:
        return path.as_posix()


def _sha256(path: Path) -> str | None:
    try:
        return hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError:
        return None


def _load(value: Any, repo: Path, *, name: str, default: Path) -> tuple[Mapping[str, Any], dict[str, Any]]:
    if isinstance(value, Mapping):
        return value, {"name": name, "mode": "supplied_mapping", "passed": True}
    path = Path(value) if value is not None else default
    if not path.is_absolute():
        path = repo / path
    try:
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rt", encoding="utf-8") as stream:
            data = json.load(stream)
    except (OSError, ValueError, EOFError, json.JSONDecodeError) as exc:
        return {}, {
            "name": name,
            "path": _relative(path, repo),
            "passed": False,
            "reason": "serialized input is unreadable",
            "error": str(exc),
        }
    return (data if isinstance(data, Mapping) else {}), {
        "name": name,
        "path": _relative(path, repo),
        "sha256": _sha256(path),
        "passed": isinstance(data, Mapping),
    }


def _pin(path: Path, repo: Path, *, name: str, expected: str | None, parse_ast: bool = False) -> dict[str, Any]:
    actual = _sha256(path)
    parsed: bool | None = None
    error = None
    if parse_ast and actual is not None:
        try:
            ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
            parsed = True
        except (OSError, SyntaxError, UnicodeError) as exc:
            parsed = False
            error = str(exc)
    row: dict[str, Any] = {
        "name": name,
        "path": _relative(path, repo),
        "sha256": actual,
        "expected_sha256": expected,
        "passed": bool(actual is not None and (expected is None or actual == expected)),
    }
    if parsed is not None:
        row["ast_parsed_without_execution"] = parsed
        row["passed"] = bool(row["passed"] and parsed)
    if error:
        row["error"] = error
    return row


def _source_pins(repo: Path, directory: Any) -> dict[str, Any]:
    base = repo / "experiments" / "root_st073"
    pins: dict[str, Any] = {}
    for name, (basename, expected) in _FILES.items():
        pins[name] = _pin(base / basename, repo, name=name, expected=expected)
    for name, (relative, expected) in _SOURCE_FILES.items():
        pins[name] = _pin(repo / relative, repo, name=name, expected=expected, parse_ast=True)
    return pins


def _skip_large_archive(value: Any, repo: Path, *, name: str, default: Path) -> tuple[Mapping[str, Any], dict[str, Any]]:
    """Keep an optional supplied family mapping without hydrating the archive."""
    if isinstance(value, Mapping):
        return value, {
            "name": name,
            "mode": "supplied_mapping",
            "passed": True,
            "hydrated": False,
            "family_metadata_only": True,
        }
    path = Path(value) if value is not None else default
    if not path.is_absolute():
        path = repo / path
    return {}, {
        "name": name,
        "path": _relative(path, repo),
        "mode": "not_hydrated",
        "passed": True,
        "hydrated": False,
        "family_metadata_only": True,
    }


def _pin_summary(pins: Mapping[str, Mapping[str, Any]]) -> bool:
    return bool(pins) and all(row.get("passed") is True for row in pins.values())


_FAMILY_KEYS = (
    "actual_five_defect_family_sha256",
    "implicit_source_sha256",
    "datum_enclosure_sha256",
)

# This is the project family recorded by the pinned bridge packets.  Keeping
# the three values here makes the check an allow-list, rather than accepting a
# caller-supplied family label as proof of source identity.
_PROJECT_FAMILY_TRIPLE = {
    "actual_five_defect_family_sha256": "3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894",
    "implicit_source_sha256": "5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae",
    "datum_enclosure_sha256": "dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010",
}

# These are the saved parent-code/input-map edges that can be checked without
# executing any producer.  A missing or changed edge is a reason to remain
# blocked; it is never repaired from an arbitrary caller mapping.
_SOURCE_INPUT_BINDINGS = (
    (
        "switch_signed_integrals_ast",
        "bridge_switch",
        "input_hashes",
        "lei_ren_part1_paper_compliant_switch_signed_integrals.py",
    ),
    (
        "macro_signed_integrals_ast",
        "bridge_switch",
        "input_hashes",
        "lei_ren_part1_paper_compliant_macro_signed_integrals.py",
    ),
    (
        "macro_signed_integrals_ast",
        "bridge_integrals",
        "input_hashes",
        "lei_ren_part1_paper_compliant_macro_signed_integrals.py",
    ),
    (
        "bridge_integrals_ast",
        "bridge_integrals",
        "input_hashes",
        "lei_ren_part1_paper_compliant_actual_bridge_integrals.py",
    ),
    (
        "bridge_switch_ast",
        "bridge_switch",
        "input_hashes",
        "lei_ren_part1_paper_compliant_actual_bridge_switch.py",
    ),
    (
        "macro_source_ast",
        "parent_report",
        "implementation_hashes",
        "src/openai_ns_reconstruction/original_macro_pressure_rows.py",
    ),
    (
        "macro_numeric_ast",
        "parent_report",
        "implementation_hashes",
        "src/openai_ns_reconstruction/original_macro_numeric.py",
    ),
    (
        "second_micro_numeric_ast",
        "parent_report",
        "implementation_hashes",
        "src/openai_ns_reconstruction/original_second_micro_numeric.py",
    ),
)


_SOURCE_FUNCTION_SPECS = {
    "switch_signed_integrals_ast": {
        "class": "CompliantSwitchSignedIntegrals",
        "function": "packet",
        "required_fragments": (
            "comparison=self.bridge.comparison(z,self.bridge.r/100)",
            "modes=_mode_rows",
            "incoming=self.bridge.actual(z,self.bridge.r/100)",
            "quotient=actualphi/phi.truncate(5)",
            "ratio=quotient*ell.exp()",
        ),
    },
    "macro_signed_integrals_ast": {
        "class": "CompliantMacroSignedIntegrals",
        "function": "packet",
        "required_fragments": (
            "entries=self.comparison[\"comparison_point_packets\"][label][\"macro\"]",
            "modes=_mode_rows",
            "radius=self.bridge.r",
            "_radius_kernels(c,radius,self.Y,1)",
            "_weighted_modes(c,modes[\"drive_hydro\"],hydro_kernels)",
        ),
    },
    "bridge_integrals_ast": {
        "class": "CompliantActualBridgeIntegrals",
        "function": "prepare",
        "required_fragments": (
            "end=self.comparison.micro(Z,2,root)",
            "modes=_mode_rows",
            "R0=self.bridge.r*c.exp(c.mpf([0,endpoints(2*self.cap)[1]]))",
            "quotient=data['phi0'].truncate(5)/phi2.truncate(5)",
        ),
    },
    "macro_source_ast": {
        "class": "OriginalMacroPressureRows",
        "function": "source_Uy",
        "required_fragments": (
            "modes=[prepared['R0']**j*R**(power-j)forjinrange(3)]",
            "row=-prepared['quotient']*module.weighted_modes(prepared['modes']['drive_'+name],modes)",
            "nonlinear=norm*o._transfer(6*b.cap*b.ctx.exp(6*b.cap))",
        ),
    },
    "macro_numeric_ast": {
        "class": "OriginalMacroNumeric",
        "function": "_geometry",
        "required_fragments": (
            "length=self.Y-iv.mpf([0,_bounds(c,2*o.hcap)[1]])",
            "logrho=iv.ln(4)+base+iv.mpf([0,_bounds(c,tiny)[1]])",
        ),
    },
    "second_micro_numeric_ast": {
        "class": "OriginalSecondMicroNumeric",
        "function": "_background",
        "required_fragments": (
            "o.bridge.prepare(Z)",
            "_second_increments(o._prepared[Z],upper,order)",
        ),
    },
}


def _family_triple(value: Any) -> dict[str, Any] | None:
    """Return exactly the project family triple from a saved mapping."""
    mapping = _mapping(value)
    nested = _mapping(mapping.get("source_family"))
    candidates = (mapping, nested)
    for candidate in candidates:
        row = {key: candidate.get(key) for key in _FAMILY_KEYS}
        if all(isinstance(row[key], str) and row[key] for key in _FAMILY_KEYS):
            return row
    return None


def _family_triple_check(documents: Mapping[str, Any]) -> dict[str, Any]:
    rows: dict[str, Any] = {}
    passed = True
    for name, document in documents.items():
        actual = _family_triple(document)
        row_passed = actual == _PROJECT_FAMILY_TRIPLE
        rows[name] = {
            "keys": list(_FAMILY_KEYS),
            "expected": dict(_PROJECT_FAMILY_TRIPLE),
            "actual": actual,
            "passed": row_passed,
        }
        passed = passed and row_passed
    return {
        "project_only_three_keys": list(_FAMILY_KEYS),
        "rows": rows,
        "passed": bool(passed),
    }


def _ast_function_binding(repo: Path, pins: Mapping[str, Any], name: str) -> dict[str, Any]:
    """Inspect one pinned function body with ``ast`` and no module execution."""
    spec = _mapping(_SOURCE_FUNCTION_SPECS.get(name))
    pin = _mapping(pins.get(name))
    relative = _SOURCE_FILES.get(name, ("", None))[0]
    path = repo / relative
    row: dict[str, Any] = {
        "name": name,
        "path": relative,
        "sha256": pin.get("sha256"),
        "expected_sha256": pin.get("expected_sha256"),
        "class": spec.get("class"),
        "function": spec.get("function"),
        "ast_only": True,
        "passed": False,
    }
    try:
        source = path.read_text(encoding="utf-8")
        tree = ast.parse(source, filename=str(path))
        target = None
        for node in tree.body:
            if isinstance(node, ast.ClassDef) and node.name == spec.get("class"):
                for child in node.body:
                    if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef)) and child.name == spec.get("function"):
                        target = child
                        break
                break
        row["line"] = [target.lineno, target.end_lineno] if target is not None else None
        normalized = "".join(ast.get_source_segment(source, target).split()) if target is not None else ""
        fragments = [str(fragment) for fragment in spec.get("required_fragments", ())]
        row["required_fragments"] = fragments
        row["fragments_present"] = {fragment: "".join(fragment.split()) in normalized for fragment in fragments}
        row["passed"] = bool(
            pin.get("passed") is True
            and target is not None
            and all(row["fragments_present"].values())
        )
    except (OSError, SyntaxError, UnicodeError, AttributeError, TypeError) as exc:
        row["error"] = str(exc)
    return row


def _source_input_binding_check(pins: Mapping[str, Any], documents: Mapping[str, Any]) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    passed = True
    for source_pin, document_name, map_name, key in _SOURCE_INPUT_BINDINGS:
        pin = _mapping(pins.get(source_pin))
        document = _mapping(documents.get(document_name))
        saved_map = _mapping(document.get(map_name))
        expected = pin.get("sha256")
        actual = saved_map.get(key)
        row_passed = bool(expected and actual == expected)
        rows.append(
            {
                "source_pin": source_pin,
                "document": document_name,
                "map": map_name,
                "key": key,
                "expected_sha256": expected,
                "actual_sha256": actual,
                "passed": row_passed,
            }
        )
        passed = passed and row_passed
    return {"rows": rows, "passed": bool(passed)}


def _source_origin_guard(
    *,
    repo: Path,
    pins: Mapping[str, Any],
    packet: Any,
    documents: Mapping[str, Any],
) -> dict[str, Any]:
    """Block the historical draft until its chart generators are equivalent.

    The saved switch packet and the requested macro source are different
    generators.  This receipt records both origins explicitly.  There is no
    accepted equivalence proof in this draft, so this function intentionally
    returns ``passed=False`` even when all read-only pins and family labels
    match.  In particular, it never turns the R100 packet into a macro R0
    packet by changing a radial weight.
    """
    packet_mapping = _mapping(packet)
    family_check = _family_triple_check(documents)
    ast_bindings = {
        name: _ast_function_binding(repo, pins, name)
        for name in _SOURCE_FUNCTION_SPECS
    }
    input_bindings = _source_input_binding_check(pins, documents)
    packet_shape = {
        "comparison_direction_modes_present": isinstance(packet_mapping.get("comparison_direction_modes"), Mapping),
        "q_rows_present": isinstance(packet_mapping.get("actual_incoming_over_frozen_comparison_quotient_axial5"), list),
        "ell_rows_present": isinstance(packet_mapping.get("actual_first_chart_logF_prefix_axial5_enclosure"), list),
        "source_domain": packet_mapping.get("source_domain"),
        "newly_selected_point_values": packet_mapping.get("newly_selected_point_values"),
        "source_modes_use_complete_comparison_moments": packet_mapping.get(
            "source_modes_use_complete_comparison_moments"
        ),
    }
    requested_macro = {
        "generator": "src/openai_ns_reconstruction/original_macro_pressure_rows.py:OriginalMacroPressureRows.source_Uy",
        "modes": "prepared['modes'] from actual_bridge_integrals.prepare after comparison.micro(Z, 2, root) (second-micro exit)",
        "quotient": "prepared['quotient'] = data['phi0'] / phi2 at the second-micro exit",
        "radius": "prepared['R0'] = bridge.r * exp([0, 2*cap])",
        "ell_scope": "macro source_Uy exp(ell)-1 nonlinear remainder",
    }
    saved_switch = {
        "generator": "experiments/root_st073/lei_ren_part1_paper_compliant_switch_signed_integrals.py:CompliantSwitchSignedIntegrals.packet",
        "modes": "comparison(z, bridge.r / 100) and actual comparison moments at R100",
        "quotient": "actual(z, bridge.r / 100).F_actual_over_F0 / comparison phi at R100",
        "ell_scope": "actual first-chart logF prefix inside the R100-to-R110 switch packet",
        "radial_scope": packet_shape.get("source_domain"),
    }
    equivalence = {
        "mode_origin_equivalent": False,
        "quotient_origin_equivalent": False,
        "ell_scope_equivalent": False,
        "proof_present": False,
        "proof_source": None,
    }
    reasons = [
        "saved switch comparison modes and q originate at R100, while the requested macro generator consumes second-micro-exit modes and q at R0",
        "saved ell is a first-switch prefix, while the requested macro source uses its own exp(ell)-1 remainder scope",
        "no pinned same-source equivalence proof binds those packet fields to the macro generator",
    ]
    if not family_check["passed"]:
        reasons.append("the pinned inputs do not share the project family triple")
    if not input_bindings["passed"]:
        reasons.append("one or more pinned sourcepy parent-code/input-map edges are absent or changed")
    if not all(_mapping(row).get("passed") is True for row in ast_bindings.values()):
        reasons.append("one or more pinned sourcepy AST function bindings are absent or changed")
    return {
        "passed": False,
        "blocked_before_compute_source": True,
        "reason_code": "SOURCE_ORIGIN_EQUIVALENCE_UNPROVED",
        "reason": "; ".join(reasons),
        "requested_macro_origin": requested_macro,
        "saved_switch_origin": saved_switch,
        "packet_shape": packet_shape,
        "family_triple_check": family_check,
        "source_input_bindings": input_bindings,
        "ast_bindings": ast_bindings,
        "equivalence": equivalence,
        "required_repair": [
            "supply a pinned macro-generator packet whose modes, quotient, R0, and ell scope are produced by the same source family",
            "or provide a separately pinned proof equating the saved R100 switch packet with the requested second-micro-exit R0 generator",
            "retain the source as blocked until that proof exists; do not reinterpret R100 modes with macro R0 weights",
        ],
    }


def build_factored_source_jet(
    *,
    nodes: Any,
    high_archive: Any,
    parent_report: Any,
    directory: Any,
    root: Any,
) -> dict[str, Any]:
    """Return a provenance receipt while the rejected adapter remains blocked.

    The signature is retained for the draft runner.  ``nodes`` is deliberately
    unused, and the large high archive is not hydrated: a supplied mapping may
    contribute only its fixed family triple.  No source constructor, producer,
    checker, fixture, or numerical source kernel is invoked.
    """
    del nodes
    repo = _repo_root(root, directory)
    base = repo / "experiments" / "root_st073"
    high, high_meta = _skip_large_archive(
        high_archive,
        repo,
        name="high_archive",
        default=base / _FILES["high_archive"][0],
    )
    parent, parent_meta = _load(
        parent_report,
        repo,
        name="parent_report",
        default=base / _FILES["parent_report"][0],
    )
    switch, switch_meta = _load(
        None,
        repo,
        name="bridge_switch",
        default=base / _FILES["bridge_switch"][0],
    )
    bridge, bridge_meta = _load(
        None,
        repo,
        name="bridge_integrals",
        default=base / _FILES["bridge_integrals"][0],
    )
    first, first_meta = _load(
        None,
        repo,
        name="first_collar",
        default=base / _FILES["first_collar"][0],
    )

    pins = _source_pins(repo, directory)
    packet = _mapping(_mapping(switch.get("actual_R100_R110_packets")).get("whole_Z")).get(
        "signed_actual_switch_source_integrals"
    )
    family_documents = {
        "parent_report": parent,
        "bridge_switch": switch,
        "bridge_integrals": bridge,
        "first_collar": first,
    }
    if high:
        family_documents["high_archive_supplied"] = high
    origin_guard = _source_origin_guard(
        repo=repo,
        pins=pins,
        packet=packet,
        documents=family_documents,
    )
    # This rejected draft has no accepted proof path.  Keep the public entry
    # fail closed even if a caller monkeypatches the receipt during a smoke
    # test; adding a proof requires a new accepted implementation.
    if origin_guard.get("passed") is True:
        origin_guard = dict(origin_guard)
        origin_guard.update(
            passed=False,
            blocked_before_compute_source=True,
            reason_code="SOURCE_ORIGIN_GATE_DISABLED",
            reason="the rejected adapter has no public source-computation path",
        )

    source_pins_passed = _pin_summary(pins)
    blocked_receipt = {
        "name": "macro_source_Uy_and_first_switch_lower_source_Uy",
        "targets": [
            {
                "name": "macro_source_Uy",
                "chart": "frozen_macro",
                "coordinate": "1",
            },
            {
                "name": "first_switch_lower_source_Uy",
                "chart": "first_switch",
                "coordinate": "s=0",
                "label": "lower endpoint",
            },
        ],
        "missing": ["ordinary B/B_Z/B_ZZ direct enclosure"],
        "reason_code": origin_guard["reason_code"],
        "reason": origin_guard["reason"],
        "blocked_before_compute_source": True,
        "source_rows_evaluated": False,
        "approximation": None,
        "source_values_selected": False,
        "origin_receipt": origin_guard,
    }
    source_bindings = {
        "macro_source_Uy": {
            "chart": "frozen_macro",
            "coordinate": "1",
            "source_formula": "blocked: source-origin equivalence is unproved",
            "requested_generator": origin_guard["requested_macro_origin"],
            "saved_packet_generator": origin_guard["saved_switch_origin"],
            "passed": False,
        },
        "first_switch_lower_source_Uy": {
            "chart": "first_switch",
            "coordinate": "s=0",
            "label": "lower endpoint",
            "source_formula": "blocked: first-switch packet cannot be aliased to macro source",
            "passed": False,
        },
        "Taylor_convention": {
            "saved": "k=D_Z^k/k!",
            "ordinary": "factorial(k)*saved Taylor coefficient",
            "conversion_applied": False,
            "blocked_before_source_algebra": True,
        },
        "point_status": {
            "newly_selected_point_values": False,
            "approximation_available": False,
            "midpoint_used": False,
        },
        "origin_guard": origin_guard,
        "passed": False,
    }
    return {
        "passed": False,
        "partial": False,
        "source_origin_guard": origin_guard,
        "source_pins": pins,
        "source_bindings": source_bindings,
        "rows": [],
        "blocked_rows": [blocked_receipt],
        "source_input_loads": {
            "high_archive": high_meta,
            "parent_report": parent_meta,
            "bridge_switch": switch_meta,
            "bridge_integrals": bridge_meta,
            "first_collar": first_meta,
            "source_origin_guard": origin_guard,
        },
        "input_scope": {
            "accepted_ancestor_replayed": False,
            "accepted_producer_called": False,
            "accepted_checker_called": False,
            "saved_JSON_gzip_only": True,
            "high_archive_hydrated": False,
            "source_AST_read_without_execution": True,
            "raw46_or_C1_inversion": False,
            "midpoint_as_source": False,
            "zero_padding": False,
        },
        "contract": {
            "ordinary_B_ZZ_recovered": False,
            "source_enclosures_retained": False,
            "full_source_materialization": False,
            "formal_swirl_component_retained": False,
            "approximation_separate_and_unavailable": True,
            "epsilon_once": False,
            "physical_R_external_moment_once": False,
            "Pstar_once": False,
            "hbB_and_hbS_distinct": True,
            "swirl_scale_not_underflow_selected": True,
            "origin_guard_passed": False,
            "source_pins_passed": source_pins_passed,
        },
        "limits": {
            "point_modal_approximation": False,
            "whole_Z_enclosure_scope": False,
            "full_signed_source_integral": False,
            "accepted_ancestor_replayed": False,
            "source_origin_equivalence_proved": False,
        },
    }


__all__ = ["build_factored_source_jet"]
