"""Integrity audit for legacy logarithmic absolute-value caps.

The accepted ``original_native_log_jet4.abs_log`` helper converts an exact
MPF endpoint through the process-wide scalar precision before taking its
logarithm.  This module audits that arithmetic boundary on two fresh rational
witnesses and records the affected direct call sites.  It does not replay an
accepted producer, checker, builder, or self-check, and it does not rewrite
the historical packets that used the helper.

The audit therefore says only that legacy logarithmic certificates require
revalidation.  It does not claim that every historical value is false, that
all transitive callers were found, or that signed native interval data is
invalid.  A separate explicit receipt from the new directed primitive is
required before a consumer may treat the affected log caps as revalidated.
"""

from __future__ import annotations

import ast
import copy
import hashlib
import json
from collections.abc import Mapping
from pathlib import Path

from mpmath import make_mpf, mp
from mpmath.ctx_iv import MPIntervalContext
from mpmath.libmp import mpf_cmp

from .original_core_signed_direction import record, require, sha
from .original_directed_abs_log import abs_log_upper, fixture_bound_check


BASE = "experiments/root_st073/"
PARENT = BASE + "lei_ren_part1_paper_compliant_current_original_restoration_finite_N_I1.json"
PARENT_SHA256 = "fbede8e8c948a2025c65671bce6265af84e57fcb25d4a7e497bf456769d647e4"
PARENT_BYTES = 39009457
MODULE = "src/openai_ns_reconstruction/original_abs_log_integrity.py"
OLD_HELPER = "src/openai_ns_reconstruction/original_native_log_jet4.py"
NEW_PRIMITIVE = "src/openai_ns_reconstruction/original_directed_abs_log.py"
NEW_PRIMITIVE_SHA256 = "96126cf50406cfc87c8d9af854a42ac9b4ed9a0222b55587d5c913bfcd9ec767"
SUPPORTED_VERIFIED_REVALIDATION_METHODS = ()


def _require(value, message):
    require(bool(value), message)


def _record_interval(value):
    return record(value)


def _jsonable(value):
    if hasattr(value, "_mpi_"):
        return _record_interval(value)
    if isinstance(value, Mapping):
        return {key: _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    return value


def _compact(text):
    return "".join(text.split())


def _ast_guard(root, relative, function, snippets):
    path = root / relative
    text = path.read_text(encoding="utf-8")
    tree = ast.parse(text, filename=relative)
    candidates = []
    for node in ast.walk(tree):
        if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            continue
        if node.name != function:
            continue
        segment = ast.get_source_segment(text, node) or ""
        if all(_compact(snippet) in _compact(segment)
               for snippet in snippets):
            candidates.append(node)
    _require(len(candidates) == 1,
             "Source AST guard is ambiguous: " + relative + ":" + function)
    node = candidates[0]
    return dict(
        passed=True,
        path=relative,
        function=function,
        required_text=list(snippets),
        AST_sha256=hashlib.sha256(
            ast.dump(node, include_attributes=False).encode()
        ).hexdigest(),
        source_sha256=sha(path),
        source_executed=False,
    )


def _source_guards(root):
    guards = dict(
        legacy_raw_bounds=_ast_guard(
            root,
            OLD_HELPER,
            "_raw_bounds",
            (
                "rows = value._mpi_",
                "lower = _mpf_endpoint(rows[0])",
                "upper = _mpf_endpoint(rows[1])",
            ),
        ),
        legacy_mpf_endpoint=_ast_guard(
            root,
            OLD_HELPER,
            "_mpf_endpoint",
            ("return mp.make_mpf(tuple(value))",),
        ),
        legacy_abs_log=_ast_guard(
            root,
            OLD_HELPER,
            "abs_log",
            (
                "magnitude = max(abs(lower), abs(upper))",
                "return c.ln(_interval(c, magnitude, magnitude))",
            ),
        ),
        directed_abs_upper=_ast_guard(
            root,
            NEW_PRIMITIVE,
            "abs_upper",
            (
                "lo_abs, hi_abs = (0, *lo[1:]), (0, *hi[1:])",
                "return context.mpf([exact, exact])",
            ),
        ),
        directed_abs_log_upper=_ast_guard(
            root,
            NEW_PRIMITIVE,
            "abs_log_upper",
            (
                "magnitude = abs_upper(context, value)",
                "return context.ln(magnitude)",
            ),
        ),
    )
    _require(all(item.get("passed") is True for item in guards.values()),
             "Absolute-log source AST guard failed")
    return guards


def _pin(root, path, expected=None, *, self_pin=False):
    actual = sha(root / path)
    if expected is not None:
        _require(actual == expected, "Changed pinned source: " + path)
    result = dict(
        path=path,
        expected_sha256=actual if expected is None else expected,
        actual_sha256=actual,
        passed=True,
        source_executed=False,
    )
    if self_pin:
        result["self_pin"] = True
    return result


def _validate_parent(root, parent):
    path = root / PARENT
    _require(path.is_file(), "Canonical finite-N parent report is missing")
    raw = path.read_bytes()
    _require(len(raw) == PARENT_BYTES,
             "Canonical finite-N parent byte count changed")
    _require(hashlib.sha256(raw).hexdigest() == PARENT_SHA256,
             "Canonical finite-N parent hash changed")
    _require(json.loads(raw) == parent,
             "Build parent mapping is not the canonical report")
    _require(parent.get("checks", {}).get("passed") is True,
             "Canonical finite-N parent did not pass")
    _require(isinstance(parent.get("source_family"), Mapping),
             "Canonical finite-N source family is missing")
    return dict(
        passed=True,
        path=PARENT,
        byte_count=len(raw),
        expected_sha256=PARENT_SHA256,
        actual_sha256=PARENT_SHA256,
        source_family=copy.deepcopy(parent["source_family"]),
        historical_pass_flag_not_consumed_as_certification=True,
        parent_checks_passed_is_historical_metadata=True,
        source_executed=False,
    )


def _source_pins(root, parent):
    implementation_hashes = parent.get("implementation_hashes")
    _require(isinstance(implementation_hashes, Mapping),
             "Parent implementation_hashes mapping is required")
    old_expected = implementation_hashes.get(OLD_HELPER)
    _require(isinstance(old_expected, str) and len(old_expected) == 64,
             "Pinned legacy log helper hash is missing")
    pins = {
        PARENT: _pin(root, PARENT, PARENT_SHA256),
        OLD_HELPER: _pin(root, OLD_HELPER, old_expected),
        NEW_PRIMITIVE: _pin(root, NEW_PRIMITIVE, NEW_PRIMITIVE_SHA256),
        MODULE: _pin(root, MODULE, self_pin=True),
    }
    return pins, old_expected


def _old_log_witness(root):
    """Compare the pinned legacy helper against exact tuple-directed logs."""
    # This is the only old-helper call in this module, and it is restricted to
    # the two fresh rational witnesses below.  No old fixture or builder is
    # imported or called.
    from .original_native_log_jet4 import abs_log as legacy_abs_log

    previous = mp.prec
    rows = []
    restored = False
    try:
        mp.prec = 53
        context = MPIntervalContext()
        context.prec = 1600
        for numerator, denominator in ((7, 10), (1, 3)):
            value = context.mpf(numerator) / denominator
            legacy = legacy_abs_log(context, value)
            directed = abs_log_upper(context, value)
            rational_reference = context.ln(value)
            legacy_upper = make_mpf(legacy._mpi_[1])
            directed_lower = make_mpf(directed._mpi_[0])
            underbound_reference = mpf_cmp(
                legacy._mpi_[1], rational_reference._mpi_[0]
            ) < 0
            underbound_directed = mpf_cmp(
                legacy._mpi_[1], directed._mpi_[0]
            ) < 0
            _require(underbound_reference and underbound_directed,
                     "Legacy abs_log witness did not expose underbound")
            rows.append(dict(
                rational=f"{numerator}/{denominator}",
                scalar_precision_bits=53,
                interval_precision_bits=1600,
                input_interval=_record_interval(value),
                legacy_log_cap=_record_interval(legacy),
                rational_reference_log_cap=_record_interval(rational_reference),
                directed_reference_log_cap=_record_interval(directed),
                legacy_upper_exact_mpf_tuple=list(legacy._mpi_[1]),
                rational_reference_lower_exact_mpf_tuple=list(
                    rational_reference._mpi_[0]
                ),
                directed_lower_exact_mpf_tuple=list(directed._mpi_[0]),
                legacy_upper_lt_rational_reference_lower=True,
                legacy_upper_lt_directed_lower=True,
                exact_tuple_comparison=True,
            ))
    finally:
        mp.prec = previous
        restored = mp.prec == previous
    _require(restored, "Global scalar precision was not restored")
    return dict(
        passed=True,
        witness_count=len(rows),
        witnesses=rows,
        scalar_precision_bits=53,
        interval_precision_bits=1600,
        global_mp_precision_restored=True,
        old_helper_fixture_not_called=True,
        old_helper_called_only_on_new_witnesses=True,
        accepted_builders_or_checkers_called=False,
        source_executed=False,
    )


def _old_module_names(alias):
    return alias in {"original_native_log_jet4", OLD_HELPER.rsplit("/", 1)[-1][:-3]}


def _scan_direct_old_calls(root, parent):
    """Scan only parent-pinned Python paths for direct old abs_log calls."""
    implementation_hashes = parent["implementation_hashes"]
    paths = sorted(
        path for path in implementation_hashes
        if isinstance(path, str) and path.endswith(".py")
    )
    _require(paths, "Parent has no pinned Python implementation paths")
    affected = []
    checked = []
    for relative in paths:
        expected = implementation_hashes[relative]
        _require(isinstance(expected, str) and len(expected) == 64,
                 "Malformed parent implementation hash: " + relative)
        path = root / relative
        actual = sha(path)
        _require(actual == expected,
                 "Changed parent-pinned implementation: " + relative)
        text = path.read_text(encoding="utf-8")
        tree = ast.parse(text, filename=relative)
        module_aliases = set()
        function_aliases = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for imported in node.names:
                    if _old_module_names(imported.name.rsplit(".", 1)[-1]):
                        module_aliases.add(
                            imported.asname or imported.name.rsplit(".", 1)[-1]
                        )
            elif isinstance(node, ast.ImportFrom):
                module_name = node.module or ""
                if (module_name.endswith("original_native_log_jet4") or
                        (node.level and any(
                            alias.name == "original_native_log_jet4"
                            for alias in node.names))):
                    for imported in node.names:
                        if imported.name == "original_native_log_jet4":
                            module_aliases.add(
                                imported.asname or imported.name
                            )
                        elif imported.name == "abs_log":
                            function_aliases.add(
                                imported.asname or imported.name
                            )
        # Include straightforward local aliases such as ``old_abs_log =
        # jet4.abs_log`` without attempting whole-program dataflow.
        for node in ast.walk(tree):
            if not isinstance(node, (ast.Assign, ast.AnnAssign)):
                continue
            value = node.value
            if not (isinstance(value, ast.Attribute) and
                    value.attr == "abs_log" and
                    isinstance(value.value, ast.Name) and
                    value.value.id in module_aliases):
                continue
            targets = node.targets if isinstance(node, ast.Assign) else [node.target]
            for target in targets:
                if isinstance(target, ast.Name):
                    function_aliases.add(target.id)
        calls = []
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            func = node.func
            direct = (
                isinstance(func, ast.Attribute)
                and func.attr == "abs_log"
                and isinstance(func.value, ast.Name)
                and func.value.id in module_aliases
            ) or (
                isinstance(func, ast.Name)
                and func.id in function_aliases
            )
            if direct:
                calls.append(getattr(node, "lineno", None))
        checked.append(relative)
        if calls:
            affected.append(dict(
                path=relative,
                expected_sha256=expected,
                actual_sha256=actual,
                module_aliases=sorted(module_aliases),
                imported_abs_log_aliases=sorted(function_aliases),
                direct_call_count=len(calls),
                call_lines=sorted(calls),
                source_executed=False,
            ))
    return dict(
        passed=True,
        scope="parent.implementation_hashes entries whose path ends with .py",
        pinned_python_path_count=len(paths),
        pinned_python_paths_sha_checked=len(checked),
        affected_direct_modules=affected,
        affected_direct_module_count=len(affected),
        direct_old_abs_log_call_count=sum(
            item["direct_call_count"] for item in affected
        ),
        dynamic_aliases_not_scanned=True,
        transitive_calls_not_scanned=True,
        all_transitive_cases_known=False,
        historical_data_not_rewritten=True,
        source_executed=False,
    )


def _gate_fixture(new_sha):
    historical = dict(passed=True, revalidation_required=True)
    rejected = False
    try:
        require_revalidated_log_caps(historical)
    except ValueError:
        rejected = True
    _require(rejected,
             "Historical passed flag bypassed the revalidation gate")
    forged = dict(
        passed=True,
        verified_method="recompute",
        receipt_type="new_directed_log_revalidation",
        transport_or_recompute="recompute",
        revalidation_scope="all_required_log_caps",
        legacy_cap_source_replaced=True,
        directed_transport_verified=True,
        source_module=NEW_PRIMITIVE,
        source_sha256=new_sha,
        source_executed=False,
        accepted_ancestors_replayed=False,
    )
    forged_rejected = False
    try:
        require_revalidated_log_caps(dict(
            revalidation_required=True,
            log_cap_revalidation_receipt=forged,
        ))
    except ValueError:
        forged_rejected = True
    _require(forged_rejected,
             "Schema-only forged revalidation receipt was admitted")
    return dict(
        passed=True,
        historical_passed_flag_rejected=True,
        schema_only_forged_receipt_rejected=True,
        explicit_new_receipt_accepted=False,
        supported_verified_revalidation_methods=list(
            SUPPORTED_VERIFIED_REVALIDATION_METHODS
        ),
        no_revalidation_receipt_admitted_yet=True,
        source_executed=False,
    )


def require_revalidated_log_caps(packet):
    """Reject until an evidence-validating revalidation verifier exists.

    A historical packet's ``passed`` bit and a schema-only receipt are both
    insufficient.  No verifier is admitted in this audit revision, so every
    current input is rejected until a new method validates an actual directed
    transport or recomputation receipt.
    """
    _require(isinstance(packet, Mapping),
             "A revalidation packet mapping is required")
    receipt = packet.get("log_cap_revalidation_receipt")
    method = receipt.get("verified_method") if isinstance(receipt, Mapping) else None
    _require(method in SUPPORTED_VERIFIED_REVALIDATION_METHODS,
             "No evidence-validating log-cap revalidation method is admitted")
    raise ValueError(
        "No evidence-validating log-cap revalidation verifier is implemented"
    )


def build(root, parent):
    """Build the source-bound legacy absolute-log integrity audit."""
    root = Path(root)
    parent_scope = _validate_parent(root, parent)
    pins, old_helper_sha = _source_pins(root, parent)
    guards = _source_guards(root)
    primitive_fixture = fixture_bound_check()
    _require(primitive_fixture.get("passed") is True,
             "New directed absolute-log fixture did not pass")
    legacy_witness = _old_log_witness(root)
    scan = _scan_direct_old_calls(root, parent)
    new_sha = pins[NEW_PRIMITIVE]["actual_sha256"]
    gate_fixture = _gate_fixture(new_sha)
    self_hash = pins[MODULE]["actual_sha256"]
    result = dict(
        passed=True,
        stage="original_abs_log_integrity_audit",
        schema="original_abs_log_integrity_v1",
        passed_semantics="integrity_audit_valid_not_legacy_bounds_valid",
        parent_report=PARENT,
        parent_report_sha256=PARENT_SHA256,
        parent_scope=parent_scope,
        source_family=copy.deepcopy(parent["source_family"]),
        source_pins=pins,
        source_guards=guards,
        self_hash_fixture_receipt=dict(
            module=MODULE,
            self_sha256=self_hash,
            fixture=primitive_fixture,
            gate_fixture=gate_fixture,
            source_executed=False,
        ),
        legacy_helper_pin_sha256=old_helper_sha,
        new_directed_primitive_pin_sha256=NEW_PRIMITIVE_SHA256,
        legacy_abs_log_witness=legacy_witness,
        direct_old_abs_log_static_scan=scan,
        legacy_old_abs_log_underbound_witness=True,
        legacy_bounds_validity="not_assessed",
        not_every_legacy_value_declared_false=True,
        supported_verified_revalidation_methods=list(
            SUPPORTED_VERIFIED_REVALIDATION_METHODS
        ),
        no_revalidation_receipt_admitted_yet=True,
        inherited_log_cap_certification_verified=False,
        revalidation_required=True,
        old_caps_and_data_preserved=True,
        exact_source_identities_not_invalidated=True,
        signed_native_direct_interval_data_not_automatically_invalidated=True,
        no_new_Rm_certification=True,
        publication_scope=dict(
            root_final_publication_supersedes_old_strict_cap_claims=True,
            old_files_edited=False,
            accepted_old_producers_replayed=False,
            accepted_old_checkers_replayed=False,
            all_transitive_or_dynamic_old_abs_log_calls_known=False,
            historical_receipts_are_metadata_not_canonical_numerical_proof=True,
        ),
        contracts=dict(
            accepted_builders_called=False,
            accepted_producers_called=False,
            accepted_checkers_called=False,
            accepted_ancestors_replayed=False,
            signed_values_selected=False,
            source_executed=False,
            historical_pass_flag_not_a_revalidation=True,
        ),
        source_executed=False,
        accepted_ancestors_replayed=False,
        full_goal_complete=False,
    )
    json.dumps(_jsonable(result), allow_nan=False)
    return result


def self_check():
    # The named self-check is a source-only fixture audit; it never calls the
    # legacy helper or any accepted producer.
    return fixture_bound_check()


__all__ = (
    "PARENT",
    "PARENT_SHA256",
    "MODULE",
    "OLD_HELPER",
    "NEW_PRIMITIVE",
    "NEW_PRIMITIVE_SHA256",
    "build",
    "require_revalidated_log_caps",
    "self_check",
)
