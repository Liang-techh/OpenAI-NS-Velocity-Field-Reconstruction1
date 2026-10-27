"""Manager-facing, source-only parameter binding for a rebuilt W template.

The public entry point accepts the frozen trace mapping and returns plain JSON
data.  No source value, W algebra, proof, provider, or numerical profile is
constructed.  Internal provenance records are created only after pin validation.
"""

from dataclasses import dataclass
from typing import Any, Callable, Dict, Mapping, Optional, Tuple


INPUT_SCHEMA = "native_W_parameter_interface_input_v1"
RESULT_SCHEMA = "original_native_W_parameter_interface"
COMPONENT = "native_W_parameter_interface"
ALIAS_TRACE_SCHEMA = "wsector_parameter_alias_trace_static_v1"
ALIAS_TRACE_RAW_SHA256 = "b483e1bd0657fc2ad223f4153044875191b5beae7e2028e99f7a08b9d05ba112"
SOURCE_TRACE_PACKET_TYPE = "wsector_source_trace_v2"
SOURCE_TRACE_RAW_SHA256 = "eb0af8bce8ab9af4bb72f8bc8a3eb5e03f5a2a6db9688ab30f37104084003dcf"


@dataclass(frozen=True, slots=True)
class _SourceAnchor:
    owner_path: str
    source_sha256: str
    function: str
    function_ast_dump_sha256: str
    assignment_ast_dump_sha256: str
    line_span: Tuple[int, int]
    target: str
    assignment_source: str


@dataclass(frozen=True, slots=True)
class _SourceFunctionRef:
    kind: str
    domain: str
    owner_path: str
    source_sha256: str
    function: str
    function_ast_dump_sha256: str
    line_span: Tuple[int, int]
    formula_role: str
    data_only: bool


@dataclass(frozen=True, slots=True)
class _NativeHandle:
    name: str
    kind: str
    domain: str
    provenance_class: str
    anchor: _SourceAnchor
    upstream_anchor: Optional[_SourceAnchor]
    source_expression: str
    exact_native_handle: bool
    value_materialized: bool
    prohibited_promotions: Tuple[str, ...]


@dataclass(frozen=True, slots=True)
class _RebuiltOrigin:
    name: str
    kind: str
    domain: str
    provenance_class: str
    origin_anchor: _SourceAnchor
    owner_anchor: Optional[_SourceAnchor]
    source_expression: str
    exact_native_handle: bool
    value_materialized: bool
    prohibited_promotions: Tuple[str, ...]


@dataclass(frozen=True, slots=True)
class _ProviderRef:
    kind: str
    provider_id: str
    domain: str
    status: str
    source_only: bool
    provider_called: bool


@dataclass(frozen=True, slots=True)
class _DeferredSlot:
    slot_name: str
    kind: str
    domain: str
    provider: _ProviderRef
    status: str
    reason: str


@dataclass(frozen=True, slots=True)
class _TemplateSlot:
    slot_name: str
    binding_kind: str
    binding: Any
    status: str
    native_exact_required: bool


@dataclass(frozen=True, slots=True)
class _TemplateSignature:
    template_id: str
    kind: str
    status: str
    generic_w_source_formula: _SourceFunctionRef
    slots: Tuple[_TemplateSlot, ...]
    history_provider_prerequisite: _ProviderRef
    prohibited_operations: Tuple[str, ...]


@dataclass(frozen=True, slots=True)
class _ResolvedSourceReference:
    kind: str
    domain: str
    source_parameter: _NativeHandle
    value_materialized: bool
    exact_native_handle: bool


class _SourceParameterRejected(TypeError):
    def __init__(self, ref_kind: str, reason: str) -> None:
        self.ref_kind = ref_kind
        self.reason = reason
        super().__init__(reason)


def derive(
    source: Mapping[str, Any],
    diagnostic: Optional[Callable[[Dict[str, Any]], Any]] = None,
) -> Dict[str, Any]:
    """Bind source metadata only and return a complete plain-JSON result."""

    counts = _fresh_counts()
    counts["source_binding_calls"] += 1
    validation_rows = []
    native_rows = []
    rebuilt_rows = []
    admission_rows = []
    rejection_rows = []
    signature_rows = []

    if diagnostic is not None and not callable(diagnostic):
        return _failure("diagnostic must be callable or None", counts, validation_rows, native_rows, rebuilt_rows, admission_rows, rejection_rows, signature_rows)
    if not _validate_source_pin(source):
        return _failure("source schema, raw SHA, literal trace fields, or AST anchors failed validation", counts, validation_rows, native_rows, rebuilt_rows, admission_rows, rejection_rows, signature_rows)

    validation_event = _event("source_validation", True, {"validation": "frozen_literal_source_pin"})
    validation_rows.append(validation_event)
    if _emit_failed(diagnostic, validation_event):
        return _failure("diagnostic callback failed during source validation", counts, validation_rows, native_rows, rebuilt_rows, admission_rows, rejection_rows, signature_rows)

    try:
        # These records are deliberately constructed only after source-pin validation.
        pstar_anchor = _anchor(
            "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_core_inlet_axial6_generic_source.py",
            "29738d46c6eb5897f0fbeefaec741d0eda4b97baed57b43ce9a7890a7dc03890",
            "CurrentOriginalCoreInletAxial6GenericSource.__init__",
            "124e57dee4f6a07e24ca63283b090ba32ae119467a12bddd9733c60f1957d3cb",
            "8df5daa7005bb16cbd9adc38836c3899d54165b59103e4efac123617ea2363ce",
            (171, 171),
            "self.Pstar",
            "self.Pstar=g.unary('exp',self.core.parameters['logP'])",
        )
        logp_anchor = _anchor(
            "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions.py",
            "d71ff19a68cfed8ace9a6857ffeca860ecee90d8d89b35b2292a7310e180fb59",
            "CurrentOriginalN1CoreFixedPointFunctions.__init__",
            "6ae3e5de35c7e10868624e4a9330cc39e598cedd77a4aaa79c14384366545672",
            "409222d0c33bbcf89d8789583228bb8b51beec03f64c56488ab83615d6d31047",
            (82, 82),
            "logP",
            "logP=source.FunctionRef(g,Gamma.pulse.preRp.operator.background['parameters']['logP'])",
        )
        delta_anchor = _anchor(
            "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions.py",
            "d71ff19a68cfed8ace9a6857ffeca860ecee90d8d89b35b2292a7310e180fb59",
            "CurrentOriginalN1CoreFixedPointFunctions.__init__",
            "6ae3e5de35c7e10868624e4a9330cc39e598cedd77a4aaa79c14384366545672",
            "e5d4d9ba78458eb55560b259cfbaa25ef6b2d0f66cf0d154c5299b4bf08295c0",
            (85, 85),
            "self.parameters.update",
            "self.parameters.update(logP=logP,logLambda=logLambda,Lambda=Lambda,epsilon=epsilon,delta=Gamma.schedule['delta'])",
        )
        rebuilt_pstar_anchor = _anchor(
            "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rm_defect_patch_inverse.py",
            "c88ac849f55f64aa1b338b0eb91c55d561006e495e2256f848a27ade07e379df",
            "_ActualRmLeadingPatch.__init__",
            "2b5ca2ef0b98b831110d8d783f7df1d4b3bc7c729baeaa0c923531c329716c93",
            "c93f16bf4018351d309da435a6863a7fff8d5996dc21e98b12bf558f1b047f6a",
            (100, 100),
            "self.Pstar",
            "self.Pstar=f.factor((0,.5,0,0,0))",
        )
        rebuilt_delta_origin_anchor = _anchor(
            "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_long_reshape_endpoint.py",
            "9b41c9a587c52f3af5de35a39695c16e07de509155ce466680d984cabd54b231",
            "OriginalLongReshapeEndpoint.owner",
            "ea9e357bed63759f7221fa9b0550b5d35ece4ed9d3325f47fb27feb29f69adfb",
            "6c3805033142e3898ae5eb9141639a698b6b033e2bb9ede3e3f890ee4da9d594",
            (200, 200),
            "delta",
            "delta=c.exp(fields.previous.read_interval(c,records['pressure_source']['compliant_source']['parameter_bounds']['log_delta']))",
        )
        rebuilt_delta_owner_anchor = _anchor(
            "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_reference_restore_functions.py",
            "a86da4ef8791c65ce678a067a999f1f24945ac7c1ca0673f6bafc8c8969702b7",
            "ActualReferenceRestoreFunctions.__init__",
            "0219b5224273ec6ee3e22d8d524637cdabadf3e8ada9061c7bae7d85153b6380",
            "86d04ece4387d296d9968c780945f87690b11606b5df18eb24c401c8a63096a1",
            (91, 91),
            "self.delta",
            "self.delta=c.mpf(delta)",
        )
        rebuilt_radius_anchor = _anchor(
            "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rm_defect_patch_inverse.py",
            "c88ac849f55f64aa1b338b0eb91c55d561006e495e2256f848a27ade07e379df",
            "_ActualRmLeadingPatch.__init__",
            "2b5ca2ef0b98b831110d8d783f7df1d4b3bc7c729baeaa0c923531c329716c93",
            "7d351bf65e39f865ab5f8666d1ab86695548e80834bc37c9799efab983a461c4",
            (101, 101),
            "self.Rm_factor",
            "self.Rm_factor=f.factor((0,5,0,0,0),10*reference.logC+c.ln(110)-6)",
        )
        generic_w_ref = _SourceFunctionRef(
            "source_function_ref",
            "original_generic_W_source_formula",
            "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.py",
            "7398b8cbba75453c566b7b7281d98827e128e89441124ea67b96d9493ca092a7",
            "generic_functions",
            "cfc50a5e623c1d71b1c2551b6d8ed72c83565ff3c74329c58fa0e92a4c7be52f",
            (33, 53),
            "original generic W source formula; data-only provenance reference",
            True,
        )
        native_pstar = _NativeHandle("Pstar", "native_source_parameter", "native_graph_parameter", "exact_native_handle", pstar_anchor, logp_anchor, "graph exp(core.parameters['logP'])", True, False, ("enclosure", "interval", "bound", "midpoint", "numeric_value"))
        native_delta = _NativeHandle("delta", "native_source_parameter", "native_schedule_parameter", "exact_native_handle", delta_anchor, None, "Gamma.schedule['delta']", True, False, ("enclosure", "interval", "bound", "midpoint", "numeric_value"))
        rebuilt_pstar = _RebuiltOrigin("Pstar", "rebuilt_original_formal_factor", "rebuilt_original_formal_factor", "rebuilt_original_owner", rebuilt_pstar_anchor, None, "formal factor record only", False, False, ("native_source_parameter", "exact_native_handle", "numeric_value", "interval", "bound", "midpoint"))
        rebuilt_delta = _RebuiltOrigin("delta", "enclosure_origin", "rebuilt_original_interval_origin", "rebuilt_original_enclosure_owner", rebuilt_delta_origin_anchor, rebuilt_delta_owner_anchor, "exp(read_interval(parameter_bounds.log_delta)) then reference-owner record", False, False, ("native_source_parameter", "exact_native_handle", "numeric_value", "interval", "bound", "midpoint"))
        rebuilt_radius = _RebuiltOrigin("radius", "rebuilt_own_radius", "rebuilt_original_radius_factor", "rebuilt_original_owner", rebuilt_radius_anchor, None, "rebuiltOwnRadius", False, False, ("native_source_parameter", "native_radius", "physical_radius_equality", "numeric_value"))
        midpoint_delta = _RebuiltOrigin("delta", "midpoint", "interval_origin_candidate", "enclosure_derived_candidate", rebuilt_delta_origin_anchor, rebuilt_delta_owner_anchor, "midpoint candidate record only; no endpoint or midpoint value", False, False, ("native_source_parameter", "exact_native_handle", "numeric_value"))

        native_rows.extend((_native_to_json(native_pstar), _native_to_json(native_delta)))
        rebuilt_rows.extend((_rebuilt_to_json(rebuilt_pstar), _rebuilt_to_json(rebuilt_delta), _rebuilt_to_json(rebuilt_radius)))

        for handle in (native_pstar, native_delta):
            try:
                resolved = require_source_parameter(handle)
            except _SourceParameterRejected as rejected:
                return _failure("frozen native handle failed exact-native admission: " + rejected.reason, counts, validation_rows, native_rows, rebuilt_rows, admission_rows, rejection_rows, signature_rows)
            admission = {"parameter_name": handle.name, "ref_kind": handle.kind, "domain": resolved.domain, "exact_native_handle": True, "status": "admitted"}
            admission_rows.append(admission)
            counts["native_admissions"] += 1
            if _emit_failed(diagnostic, _event("native_handle_admission", True, admission)):
                return _failure("diagnostic callback failed during native admission", counts, validation_rows, native_rows, rebuilt_rows, admission_rows, rejection_rows, signature_rows)

        for record in (rebuilt_delta, midpoint_delta):
            try:
                require_source_parameter(record)
            except _SourceParameterRejected as rejected:
                rejection = {"parameter_name": record.name, "ref_kind": rejected.ref_kind, "domain": record.domain, "status": "rejected", "reject_reason": rejected.reason}
                rejection_rows.append(rejection)
                counts["rejections"] += 1
                if _emit_failed(diagnostic, _event("source_reference_rejection", True, rejection)):
                    return _failure("diagnostic callback failed during source-reference rejection", counts, validation_rows, native_rows, rebuilt_rows, admission_rows, rejection_rows, signature_rows)
            else:
                return _failure("a prohibited enclosure-like record was admitted", counts, validation_rows, native_rows, rebuilt_rows, admission_rows, rejection_rows, signature_rows)

        provider = _ProviderRef("deferred_source_history_provider", "source_native_history_provider", "source_native_history", "required_before_instantiation", True, False)
        signature = _TemplateSignature(
            "rebuilt_original_W_source_native_template_v1",
            "source_native_rebuilt_template_signature",
            "uninstantiated",
            generic_w_ref,
            (
                _TemplateSlot("S", native_pstar.kind, native_pstar, "bound_exact_native_handle", True),
                _TemplateSlot("delta", native_delta.kind, native_delta, "bound_exact_native_handle", True),
                _TemplateSlot("radius", rebuilt_radius.kind, rebuilt_radius, "rebuiltOwnRadius", False),
                _deferred_slot("V", provider),
                _deferred_slot("m", provider),
                _deferred_slot("e", provider),
                _deferred_slot("P0", provider),
                _deferred_slot("p", provider),
            ),
            provider,
            ("evaluate_W", "build_W_algebra", "promote_rebuilt_enclosure_to_native_parameter", "materialize_source_value", "instantiate_history_without_source_provider"),
        )
        signature_json = _template_to_json(signature)
        signature_rows.append(signature_json)
        counts["template_signatures"] += 1
        if _emit_failed(diagnostic, _event("template_signature", True, signature_json)):
            return _failure("diagnostic callback failed during template signature", counts, validation_rows, native_rows, rebuilt_rows, admission_rows, rejection_rows, signature_rows)
    except Exception:
        return _failure("source binding construction failed", counts, validation_rows, native_rows, rebuilt_rows, admission_rows, rejection_rows, signature_rows)

    return {
        "schema": RESULT_SCHEMA,
        "component": COMPONENT,
        "passed": True,
        "metadata_only": True,
        "source_only": True,
        "actual_counts": counts,
        "new_flags": {
            "native_W_parameter_source_handles_bound": True,
            "native_W_parameter_enclosure_parameters_separated": True,
            "native_W_parameter_uninstantiated_template_signature_defined": True,
        },
        "payload": {
            "source_input_provenance": {"input_schema": INPUT_SCHEMA, "alias_trace_schema": ALIAS_TRACE_SCHEMA, "alias_trace_raw_sha256": ALIAS_TRACE_RAW_SHA256, "source_trace_packet_type": SOURCE_TRACE_PACKET_TYPE, "source_trace_raw_sha256": SOURCE_TRACE_RAW_SHA256, "parent_data": "opaque_optional_not_consumed"},
            "native_source_handle_records": native_rows,
            "rebuilt_original_records": rebuilt_rows,
            "template_signature_records": signature_rows,
            "native_handle_admission_events": admission_rows,
            "source_reference_rejection_events": rejection_rows,
            "preserved_inherited_metadata": _preserved_metadata(),
        },
        "failure": None,
    }


def require_source_parameter(ref: object) -> _ResolvedSourceReference:
    """Resolve only a native handle carrying an exact frozen fingerprint."""

    if not isinstance(ref, _NativeHandle):
        ref_kind = str(getattr(ref, "kind", "unknown"))
        if ref_kind in ("enclosure", "enclosure_origin", "interval", "bound", "midpoint"):
            raise _SourceParameterRejected(ref_kind, "enclosure, interval, bound, or midpoint records are not source parameter values")
        raise _SourceParameterRejected(ref_kind, "only exact-native-handle provenance can resolve as a source parameter")
    if ref.kind != "native_source_parameter" or ref.provenance_class != "exact_native_handle":
        raise _SourceParameterRejected(ref.kind, "exact-native-handle provenance is required")
    if not ref.exact_native_handle or ref.value_materialized:
        raise _SourceParameterRejected(ref.kind, "materialized or non-exact source values are prohibited")
    if ref.name == "Pstar":
        if not _anchor_fingerprint(ref.anchor, "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_core_inlet_axial6_generic_source.py", "29738d46c6eb5897f0fbeefaec741d0eda4b97baed57b43ce9a7890a7dc03890", "8df5daa7005bb16cbd9adc38836c3899d54165b59103e4efac123617ea2363ce"):
            raise _SourceParameterRejected(ref.kind, "Pstar owner fingerprint is not an accepted native source anchor")
        if not _anchor_fingerprint(ref.upstream_anchor, "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions.py", "d71ff19a68cfed8ace9a6857ffeca860ecee90d8d89b35b2292a7310e180fb59", "409222d0c33bbcf89d8789583228bb8b51beec03f64c56488ab83615d6d31047"):
            raise _SourceParameterRejected(ref.kind, "Pstar upstream logP fingerprint is not an accepted native source anchor")
    elif ref.name == "delta":
        if not _anchor_fingerprint(ref.anchor, "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions.py", "d71ff19a68cfed8ace9a6857ffeca860ecee90d8d89b35b2292a7310e180fb59", "e5d4d9ba78458eb55560b259cfbaa25ef6b2d0f66cf0d154c5299b4bf08295c0"):
            raise _SourceParameterRejected(ref.kind, "delta owner fingerprint is not an accepted native source anchor")
        if ref.upstream_anchor is not None:
            raise _SourceParameterRejected(ref.kind, "delta has no accepted upstream alias anchor")
    else:
        raise _SourceParameterRejected(ref.kind, "parameter name is not an accepted native source parameter")
    return _ResolvedSourceReference("resolved_source_reference", ref.domain, ref, False, True)


def _fresh_counts() -> Dict[str, int]:
    return {"source_binding_calls": 0, "native_admissions": 0, "rejections": 0, "template_signatures": 0, "new_mathematical_proofs": 0, "new_algebra_family_calls": 0, "numeric_calls": 0, "source_provider_calls": 0, "old_math_calls": 0}


def _failure(reason: str, counts: Dict[str, int], validation_rows: list, native_rows: list, rebuilt_rows: list, admission_rows: list, rejection_rows: list, signature_rows: list) -> Dict[str, Any]:
    return {
        "schema": RESULT_SCHEMA,
        "component": COMPONENT,
        "passed": False,
        "metadata_only": True,
        "source_only": True,
        "actual_counts": dict(counts),
        "new_flags": {"native_W_parameter_source_handles_bound": False, "native_W_parameter_enclosure_parameters_separated": False, "native_W_parameter_uninstantiated_template_signature_defined": False},
        "payload": {"status": "failed_prefix", "source_input_provenance": {"input_schema": INPUT_SCHEMA, "alias_trace_schema": ALIAS_TRACE_SCHEMA, "alias_trace_raw_sha256": ALIAS_TRACE_RAW_SHA256, "source_trace_packet_type": SOURCE_TRACE_PACKET_TYPE, "source_trace_raw_sha256": SOURCE_TRACE_RAW_SHA256, "parent_data": "opaque_optional_not_consumed"}, "source_validation_events": list(validation_rows), "native_source_handle_records": list(native_rows), "rebuilt_original_records": list(rebuilt_rows), "template_signature_records": list(signature_rows), "native_handle_admission_events": list(admission_rows), "source_reference_rejection_events": list(rejection_rows), "preserved_inherited_metadata": _preserved_metadata()},
        "failure": {"kind": "source_binding_failure", "reject_reason": reason},
    }


def _validate_source_pin(source: object) -> bool:
    if not isinstance(source, Mapping) or source.get("schema") != INPUT_SCHEMA:
        return False
    if source.get("alias_trace_raw_sha256") != ALIAS_TRACE_RAW_SHA256 or source.get("source_trace_raw_sha256") != SOURCE_TRACE_RAW_SHA256:
        return False
    alias = source.get("alias_trace")
    trace = source.get("source_trace")
    if not isinstance(alias, Mapping) or not isinstance(trace, Mapping) or alias.get("schema") != ALIAS_TRACE_SCHEMA or trace.get("packet_type") != SOURCE_TRACE_PACKET_TYPE:
        return False
    checks = (
        (alias, ("source_files", "native_core_inlet", "path"), "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_core_inlet_axial6_generic_source.py"),
        (alias, ("source_files", "native_core_inlet", "source_sha256"), "29738d46c6eb5897f0fbeefaec741d0eda4b97baed57b43ce9a7890a7dc03890"),
        (alias, ("source_files", "native_core_parameter_initializer", "source_sha256"), "d71ff19a68cfed8ace9a6857ffeca860ecee90d8d89b35b2292a7310e180fb59"),
        (alias, ("source_files", "rebuilt_rm_patch_owner", "source_sha256"), "c88ac849f55f64aa1b338b0eb91c55d561006e495e2256f848a27ade07e379df"),
        (alias, ("source_files", "rebuilt_delta_source", "source_sha256"), "9b41c9a587c52f3af5de35a39695c16e07de509155ce466680d984cabd54b231"),
        (alias, ("source_files", "rebuilt_reference_delta_owner", "source_sha256"), "a86da4ef8791c65ce678a067a999f1f24945ac7c1ca0673f6bafc8c8969702b7"),
        (trace, ("source_functions", "native_generic_functions", "path"), "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.py"),
        (trace, ("source_functions", "native_generic_functions", "source_sha256"), "7398b8cbba75453c566b7b7281d98827e128e89441124ea67b96d9493ca092a7"),
        (trace, ("source_functions", "native_generic_functions", "function_ast_dump_sha256"), "cfc50a5e623c1d71b1c2551b6d8ed72c83565ff3c74329c58fa0e92a4c7be52f"),
    )
    for document, path, expected in checks:
        if _read_path(document, path) != expected:
            return False
    assignments = (
        (alias, ("source_files", "native_core_inlet", "assignments", 0, "ast_dump_sha256"), "8df5daa7005bb16cbd9adc38836c3899d54165b59103e4efac123617ea2363ce"),
        (alias, ("source_files", "native_core_parameter_initializer", "assignments_or_calls", 0, "ast_dump_sha256"), "409222d0c33bbcf89d8789583228bb8b51beec03f64c56488ab83615d6d31047"),
        (alias, ("source_files", "native_core_parameter_initializer", "assignments_or_calls", 1, "ast_dump_sha256"), "e5d4d9ba78458eb55560b259cfbaa25ef6b2d0f66cf0d154c5299b4bf08295c0"),
        (alias, ("source_files", "rebuilt_rm_patch_owner", "assignments", 2, "ast_dump_sha256"), "c93f16bf4018351d309da435a6863a7fff8d5996dc21e98b12bf558f1b047f6a"),
        (alias, ("source_files", "rebuilt_rm_patch_owner", "assignments", 3, "ast_dump_sha256"), "7d351bf65e39f865ab5f8666d1ab86695548e80834bc37c9799efab983a461c4"),
        (alias, ("source_files", "rebuilt_delta_source", "assignments", 0, "ast_dump_sha256"), "6c3805033142e3898ae5eb9141639a698b6b033e2bb9ede3e3f890ee4da9d594"),
        (alias, ("source_files", "rebuilt_reference_delta_owner", "assignment", "ast_dump_sha256"), "86d04ece4387d296d9968c780945f87690b11606b5df18eb24c401c8a63096a1"),
    )
    return all(_read_path(document, path) == expected for document, path, expected in assignments)


def _read_path(document: object, path: Tuple[Any, ...]) -> Any:
    current = document
    for part in path:
        if isinstance(part, int):
            if not isinstance(current, (list, tuple)) or part >= len(current):
                return None
            current = current[part]
        else:
            if not isinstance(current, Mapping):
                return None
            current = current.get(part)
    return current


def _anchor(*values: Any) -> _SourceAnchor:
    return _SourceAnchor(*values)


def _anchor_fingerprint(anchor: Optional[_SourceAnchor], owner_path: str, source_sha256: str, assignment_ast_dump_sha256: str) -> bool:
    return anchor is not None and anchor.owner_path == owner_path and anchor.source_sha256 == source_sha256 and anchor.assignment_ast_dump_sha256 == assignment_ast_dump_sha256


def _event(event_kind: str, passed: bool, event_payload: Mapping[str, Any]) -> Dict[str, Any]:
    return {"schema": RESULT_SCHEMA, "component": COMPONENT, "event_kind": event_kind, "passed": passed, "metadata_only": True, "source_only": True, "event_payload": dict(event_payload)}


def _emit_failed(diagnostic: Optional[Callable[[Dict[str, Any]], Any]], event: Dict[str, Any]) -> bool:
    if diagnostic is None:
        return False
    try:
        diagnostic(event)
    except Exception:
        return True
    return False


def _deferred_slot(slot_name: str, provider: _ProviderRef) -> _TemplateSlot:
    deferred = _DeferredSlot(slot_name, "deferred_source_slot", "source_native_history", provider, "deferred", "history provider required before field instantiation")
    return _TemplateSlot(slot_name, deferred.kind, deferred, deferred.status, False)


def _anchor_to_json(anchor: Optional[_SourceAnchor]) -> Optional[Dict[str, Any]]:
    if anchor is None:
        return None
    return {"owner_path": anchor.owner_path, "source_sha256": anchor.source_sha256, "function": anchor.function, "function_ast_dump_sha256": anchor.function_ast_dump_sha256, "assignment_ast_dump_sha256": anchor.assignment_ast_dump_sha256, "line_span": list(anchor.line_span), "target": anchor.target, "assignment_source": anchor.assignment_source}


def _native_to_json(handle: _NativeHandle) -> Dict[str, Any]:
    return {"name": handle.name, "kind": handle.kind, "domain": handle.domain, "provenance_class": handle.provenance_class, "anchor": _anchor_to_json(handle.anchor), "upstream_anchor": _anchor_to_json(handle.upstream_anchor), "source_expression": handle.source_expression, "exact_native_handle": handle.exact_native_handle, "value_materialized": handle.value_materialized, "prohibited_promotions": list(handle.prohibited_promotions)}


def _rebuilt_to_json(record: _RebuiltOrigin) -> Dict[str, Any]:
    return {"name": record.name, "kind": record.kind, "domain": record.domain, "provenance_class": record.provenance_class, "origin_anchor": _anchor_to_json(record.origin_anchor), "owner_anchor": _anchor_to_json(record.owner_anchor), "source_expression": record.source_expression, "exact_native_handle": record.exact_native_handle, "value_materialized": record.value_materialized, "prohibited_promotions": list(record.prohibited_promotions)}


def _source_ref_to_json(ref: _SourceFunctionRef) -> Dict[str, Any]:
    return {"kind": ref.kind, "domain": ref.domain, "owner_path": ref.owner_path, "source_sha256": ref.source_sha256, "function": ref.function, "function_ast_dump_sha256": ref.function_ast_dump_sha256, "line_span": list(ref.line_span), "formula_role": ref.formula_role, "data_only": ref.data_only}


def _provider_to_json(provider: _ProviderRef) -> Dict[str, Any]:
    return {"kind": provider.kind, "provider_id": provider.provider_id, "domain": provider.domain, "status": provider.status, "source_only": provider.source_only, "provider_called": provider.provider_called}


def _binding_to_json(binding: Any) -> Dict[str, Any]:
    if isinstance(binding, _NativeHandle):
        return _native_to_json(binding)
    if isinstance(binding, _RebuiltOrigin):
        return _rebuilt_to_json(binding)
    return {"slot_name": binding.slot_name, "kind": binding.kind, "domain": binding.domain, "provider": _provider_to_json(binding.provider), "status": binding.status, "reason": binding.reason}


def _template_to_json(signature: _TemplateSignature) -> Dict[str, Any]:
    return {"template_id": signature.template_id, "kind": signature.kind, "status": signature.status, "generic_w_source_formula": _source_ref_to_json(signature.generic_w_source_formula), "slots": [{"slot_name": slot.slot_name, "binding_kind": slot.binding_kind, "binding": _binding_to_json(slot.binding), "status": slot.status, "native_exact_required": slot.native_exact_required} for slot in signature.slots], "history_provider_prerequisite": _provider_to_json(signature.history_provider_prerequisite), "prohibited_operations": list(signature.prohibited_operations)}


def _preserved_metadata() -> Dict[str, Any]:
    return {"cumulative_terms": 212, "raw_signed_terms": 213, "reduced_terms": 183, "six_signed_channels_retained": True, "reserved_contribution_registered_once": True, "J_terms": 16, "K_terms": 17, "outer_B_terms": 2, "all_old_error_data_retained": True, "inherited_errors_deleted": False, "old_error_terms_removed": False, "old_errors_removed": False, "all_original_caps_retained": True, "cap_to_source_bridge_proved": False, "numeric_centers_used_as_source_values": False, "same_family_is_not_a_function_identity": True, "selected_N_metadata_matched": True}
