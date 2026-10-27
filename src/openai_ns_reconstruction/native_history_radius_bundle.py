"""Inert source-definition provider for the native history/radius bundle.

The module deliberately contains no import-time state and no executable entry
point.  ``derive`` consumes already parsed source packets and returns opaque
native graph references only.  It never imports a project module, opens an
archive, evaluates a graph, materializes ``N`` or calls a mathematical
provider.
"""


_HISTORY_SHA_V1 = "4c74e9ce767ed5d4e6d896d5bd411d5b137e8aa534058a8275cda95769e131ae"
_HISTORY_SHA_V2 = "4192b97d6d4135b3d5d3f3ff69c86982baa1b895b06412a188664a189420c114"
_HISTORY_SHA_V2_LEGACY = "d7ecc5982c0064e8e475df6be93c43c86edf9690f2e42e0b61954220a0a925bc"
_RADIUS_SHA = "1ecae6e29c77d95ae21545ff20bbfc6db69f3e5207e3fdf5cace1db023f400ee"
_RADIUS_SHA_LEGACY = "08b0fdbdef53b890e08fdb2c694be7e11ba51f62fd99cef66215e2676dd792ac"
_PARAMETER_SHA = "7d45286f9a4c7a619f2edf06f822ae709dcf4e2175ca001bae7bd7455d7d3877"
_PARENT_REPORT_SHA = "b817a0f0fca623bd40ec8ab074b1e7c851440f5c2b59a934c2ef5e2b3dce9c26"
_NATIVE_ARCHIVE_SHA = "75895e1637411cfd2e490778f495f3adf904295a8840184cdc467e006b23305d"
_INPUT_SCHEMA = "native_history_radius_bundle_input_v1"
_RESULT_SCHEMA = "original_native_history_radius_bundle"
_COMPONENT = "native_history_radius_bundle"
_PROVIDER_ID = "source_native_history_provider"


def _is_map(value):
    return hasattr(value, "get") and hasattr(value, "keys")


def _copy_plain(value):
    """Copy parsed JSON containers without copying or evaluating functions."""
    if isinstance(value, dict):
        return {key: _copy_plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_copy_plain(item) for item in value]
    return value


def _same_sequence(left, right):
    return isinstance(left, (list, tuple)) and list(left) == list(right)


def _expected_selectors():
    return {
        "background": {
            "V": [0, 0, 0, 0, 0],
            "m": [241647, 241648, 241649, 241650, 241651],
            "e": [241744, 241745, 241746, 241747, 241748],
            "P0": [32693, 32694, 32695, 32696, 66509],
            "p": [241772, 241773, 241774, 241775, 241776],
            "E": [241552, 241553, 241554, 241555, 241556],
            "h": [241687, 241688, 241689, 241690, 241691],
            "k": [241709, 241710, 241711, 241712, 241713],
            "pressure": [241779, 241780, 241781, 241782, 241783],
        },
        "modified": {
            "V": [244378, 244379, 244380, 244381, 244382],
            "m": [244710, 244711, 244712, 244713, 244714],
            "e": [244728, 244729, 244730, 244731, 244732],
            "P0": [32693, 32694, 32695, 32696, 66509],
            "p": [244734, 244735, 244736, 244737, 244738],
            "E": [244544, 244546, 244548, 244550, 244552],
            "h": [244716, 244717, 244718, 244719, 244720],
            "k": [244722, 244723, 244724, 244725, 244726],
            "pressure": [244740, 244741, 244742, 244743, 244744],
        },
    }


def _quantities():
    return ("V", "m", "e", "P0", "p", "E", "h", "k", "pressure")


def _zero_counts():
    return {
        "definition_provider_entries": 0,
        "native_branches_served": 0,
        "field_definition_queries": 0,
        "radius_definition_queries": 0,
        "native_bundle_manifests": 0,
        "mixed_branch_rejections": 0,
        "rebuilt_frame_rejections": 0,
        "duplicate_normalization_rejections": 0,
        "new_mathematical_proofs": 0,
        "new_algebra_family_calls": 0,
        "numeric_calls": 0,
        "source_value_provider_calls": 0,
        "old_math_calls": 0,
    }


def _false_flags():
    return {
        "native_history_radius_definitions_served": False,
        "native_branch_frame_separation_enforced": False,
        "native_ordinary_Z_projection_preserved": False,
    }


def _unwrap_packet(entry):
    if not _is_map(entry):
        raise ValueError("packet entry must be a parsed JSON mapping")
    raw_sha = None
    for key in ("raw_sha256", "raw_sha", "sha256"):
        if entry.get(key) is not None:
            raw_sha = entry.get(key)
            break
    packet = None
    for key in ("packet", "literal", "value", "data"):
        candidate = entry.get(key)
        if _is_map(candidate):
            packet = candidate
            break
    if packet is None:
        packet = {
            key: value
            for key, value in entry.items()
            if key not in ("raw_sha256", "raw_sha", "sha256")
        }
    return packet, raw_sha


def _get_input_packet(source, name):
    entry = source.get(name)
    packet, raw_sha = _unwrap_packet(entry)
    # A v2 literal must travel with the authenticated v1 literal.  Accept the
    # explicit sibling spelling used by the workflow as well as a combined
    # {base_packet, addendum} literal.
    if name == "history_trace" and _is_map(packet) and packet.get("packet_type") == "native_history_bundle_source_trace_v2_addendum":
        base_literal = None
        if _is_map(entry):
            for key in ("base_literal", "base_packet_literal", "base_packet_data", "base_v1", "base_packet"):
                candidate = entry.get(key)
                if _is_map(candidate) and candidate.get("packet_type") in ("native_history_bundle_source_trace", "native_history_bundle_source_trace_v1"):
                    base_literal = candidate
                    break
        if base_literal is None:
            for key in ("history_trace_base", "history_trace_v1", "history_trace_v1_packet"):
                sibling = source.get(key)
                if sibling is None:
                    continue
                candidate, sibling_sha = _unwrap_packet(sibling)
                if _is_map(candidate) and candidate.get("packet_type") in ("native_history_bundle_source_trace", "native_history_bundle_source_trace_v1"):
                    if sibling_sha is not None and str(sibling_sha).lower() != _HISTORY_SHA_V1:
                        raise ValueError("history v1 base raw SHA-256 mismatch")
                    base_literal = candidate
                    break
        if base_literal is not None:
            packet = {"base_packet": base_literal, "addendum": packet}
    if raw_sha is None:
        for key in (name + "_raw_sha256", name + "_raw_sha"):
            if source.get(key) is not None:
                raw_sha = source.get(key)
                break
    if not isinstance(raw_sha, str):
        raise ValueError(name + " raw SHA-256 is required")
    return packet, raw_sha.lower()


def _split_history_packet(packet, raw_sha):
    """Return base v1 data plus optional authoritative v2 addendum."""
    packet_type = packet.get("packet_type")
    if packet_type in ("native_history_bundle_source_trace", "native_history_bundle_source_trace_v1"):
        return packet, None
    if packet_type == "native_history_bundle_source_trace_v2_addendum":
        base = None
        for key in ("base_packet_literal", "base_packet_data", "base_packet_json"):
            candidate = packet.get(key)
            if _is_map(candidate) and candidate.get("packet_type"):
                base = candidate
                break
        if base is None:
            candidate = packet.get("base_packet")
            if _is_map(candidate) and candidate.get("packet_type"):
                base = candidate
        if base is None:
            raise ValueError("authoritative history v2 packet requires its authenticated v1 base literal")
        return base, packet
    if _is_map(packet.get("base_packet")) and _is_map(packet.get("addendum")):
        base = packet.get("base_packet")
        addendum = packet.get("addendum")
        if base.get("packet_type") in ("native_history_bundle_source_trace", "native_history_bundle_source_trace_v1"):
            return base, addendum
    if _is_map(packet.get("v1")) and _is_map(packet.get("v2")):
        return packet.get("v1"), packet.get("v2")
    # A root-frozen merged packet may retain v1 fields with a v2 packet type.
    if _is_map(packet.get("executable_native_definition_provider_plan")):
        return packet, packet
    raise ValueError("history_trace is not a recognized v1/v2 literal packet")


def _validate_selector_map(plan, expected):
    selectors = plan.get("branch_field_selectors") if _is_map(plan) else None
    if not _is_map(selectors):
        raise ValueError("history trace is missing typed branch selectors")
    for branch, field_map in expected.items():
        actual_fields = selectors.get(branch)
        if not _is_map(actual_fields):
            raise ValueError("history trace is missing branch " + branch)
        for quantity, expected_nodes in field_map.items():
            actual_nodes = actual_fields.get(quantity)
            if not _same_sequence(actual_nodes, expected_nodes):
                raise ValueError("history trace selector mismatch for " + branch + "." + quantity)
    return selectors


def _validate_p0(base, addendum):
    expected_nodes = [32693, 32694, 32695, 32696, 66509]
    expected_factors = [1, 1, 2, 6, 24]
    source = addendum.get("P0_factorial_once_guard") if _is_map(addendum) else None
    if _is_map(source):
        saved = source.get("saved_ordinary_nodes")
        if not _is_map(saved):
            raise ValueError("v2 P0 guard is missing saved ordinary nodes")
        for index, order in enumerate(range(5)):
            row = saved.get("order_" + str(order))
            if not _is_map(row):
                raise ValueError("v2 P0 guard is missing order " + str(order))
            if row.get("node") != expected_nodes[index] or row.get("factorial_already_applied") != expected_factors[index]:
                raise ValueError("v2 P0 ordinary projection mismatch")
        if source.get("double_factorial_guard") is not True:
            raise ValueError("v2 P0 double-factorial guard is required")
        return {"nodes": expected_nodes, "factors": expected_factors}
    source = base.get("P0_ordinary_projection") if _is_map(base) else None
    rows = source.get("ordinary_rows") if _is_map(source) else None
    if not isinstance(rows, list) or len(rows) < 5:
        raise ValueError("base P0 ordinary rows are missing")
    for index, row in enumerate(rows[:5]):
        if not _is_map(row) or row.get("order") != index or row.get("node") != expected_nodes[index] or row.get("factorial") != expected_factors[index]:
            raise ValueError("base P0 ordinary projection mismatch")
    if source.get("factorial_applied_once") is not True or source.get("do_not_multiply_saved_ordinary_P0_nodes_again") is not True:
        raise ValueError("base P0 factorial-once guard is required")
    return {"nodes": expected_nodes, "factors": expected_factors}


def _validate_history(packet, raw_sha):
    base, addendum = _split_history_packet(packet, raw_sha)
    packet_type = packet.get("packet_type")
    if packet_type is None and _is_map(addendum):
        packet_type = addendum.get("packet_type")
    if packet_type == "native_history_bundle_source_trace_v1":
        raise ValueError("authoritative v2 history addendum is required with the v1 base literal")
    elif packet_type == "native_history_bundle_source_trace_v2_addendum":
        if raw_sha == _HISTORY_SHA_V2_LEGACY:
            raise ValueError("legacy v2 history addendum is retained but not authoritative; use v4")
        if raw_sha != _HISTORY_SHA_V2:
            raise ValueError("history v2 raw SHA-256 mismatch")
        addendum_schema = addendum.get("schema") if _is_map(addendum) else packet.get("schema")
        if addendum_schema != "native_functionjet_row_convention_and_provider_refs_v2":
            raise ValueError("history v2 schema mismatch")
        base_pointer = addendum.get("base_packet") if _is_map(addendum) else packet.get("base_packet")
        if _is_map(base_pointer) and base_pointer.get("sha256") not in (None, _HISTORY_SHA_V1):
            raise ValueError("history v2 base packet SHA-256 mismatch")
        convention = addendum.get("declaration_source_trace") if _is_map(addendum) else None
        semantics = convention.get("ordinary_row_semantics") if _is_map(convention) else None
        if not _is_map(convention) or convention.get("row_convention_confidence") != "high" or not _is_map(semantics) or semantics.get("row_kind") != "ordinary Z derivatives" or semantics.get("shift_does_not_multiply_by_order") is not True:
            raise ValueError("v2 ordinary FunctionJet row convention is not authenticated")
    elif raw_sha != _HISTORY_SHA_V2:
        raise ValueError("history raw SHA-256 mismatch")

    if not _is_map(base):
        raise ValueError("history base packet must be a mapping")
    if base.get("packet_type") not in ("native_history_bundle_source_trace", "native_history_bundle_source_trace_v1"):
        raise ValueError("authenticated v1 history base literal is required")
    if base.get("schema") != "native_history_bundle_source_trace_v1":
        raise ValueError("authenticated v1 history base schema mismatch")
    if base.get("status") not in (None, "bounded_static_source_ast_literal_trace"):
        raise ValueError("history status is not a static trace")
    expected = _expected_selectors()
    plan = None
    if _is_map(addendum):
        plan = addendum.get("executable_native_definition_provider_plan")
    if not _is_map(plan):
        plan = base.get("frozen_provider_definition_plan")
    if not _is_map(plan):
        raise ValueError("history provider plan is missing")
    if plan.get("provider_id") != _PROVIDER_ID:
        raise ValueError("history provider id mismatch")
    if plan.get("explicit_branch_argument_required") is not True or plan.get("branch_mixing_forbidden") is not True:
        raise ValueError("history branch guard is missing")
    if plan.get("ordinary_orders") not in (None, [0, 1, 2, 3, 4]):
        raise ValueError("history ordinary order selector mismatch")
    selectors = _validate_selector_map(plan, expected)
    p0 = _validate_p0(base, addendum)
    p2 = base.get("P2_branch_and_inputs")
    if not _is_map(p2):
        raise ValueError("saved P2 route selectors are required")
    if p2.get("selected_p2_handles") != [242304, 242308, 242314, 242321] or p2.get("selected_orders") != [0, 1, 2, 3]:
        raise ValueError("saved background P2 route mismatch")
    if p2.get("separate_modified_generic_p2_handles") != [249570, 249574, 249580, 249587, 249595]:
        raise ValueError("saved modified P2 route mismatch")
    creators = base.get("source_creator_ast_records")
    if not isinstance(creators, list) or len(creators) < 1:
        raise ValueError("saved source creator AST records are required")
    if not _is_map(base.get("native_parameter_owners")):
        raise ValueError("saved native parameter owner chain is required")
    finite = base.get("finite_N_literal_source")
    if not _is_map(finite) or finite.get("N_node") != 18 or finite.get("N_definition") != "N=2^(2^J)" or finite.get("J_literal") != "1358356628656378313" or finite.get("N_materialized") is not False:
        raise ValueError("saved finite-N literal source guard is required")
    history_rows = base.get("history_normalization_and_Rbase")
    if not _is_map(history_rows) or history_rows.get("incoming_Rbase_selector") != "background_endpoints.O2_buffer.history":
        raise ValueError("saved history normalization selectors are required")

    # The v1 packet carries the chart and owner records; v2 carries exact
    # replacements.  Validate whichever records are present and retain them.
    chart = base.get("native_chart_signature")
    if _is_map(chart):
        if chart.get("chart") != "O3_transition" or chart.get("coordinate", {}).get("node") != 241536:
            raise ValueError("history chart selector mismatch")
        frame = chart.get("frame")
        if not _is_map(frame) or frame.get("radius_node") != 241787 or frame.get("offset_literal") != "17 + exp(40) + x":
            raise ValueError("history native frame mismatch")
    if not _is_map(addendum) or addendum.get("packet_type") != "native_history_bundle_source_trace_v2_addendum":
        raise ValueError("authoritative v2 history addendum is required")
    zero = addendum.get("native_zero_source") if _is_map(addendum) else None
    if _is_map(zero):
        if zero.get("field") != "V" or zero.get("graph_node", {}).get("node") != 0 or zero.get("graph_node", {}).get("not_padding") is not True:
            raise ValueError("history V zero source mismatch")
    radius_raw = addendum.get("native_radius_raw_record") if _is_map(addendum) else None
    if _is_map(radius_raw):
        if radius_raw.get("archive_path") != "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions.json.gz" or radius_raw.get("raw_archive_sha256") != _NATIVE_ARCHIVE_SHA or radius_raw.get("radius_node") != 241787 or radius_raw.get("offset_node") != 241542:
            raise ValueError("v2 native radius raw selector mismatch")
        raw_record = radius_raw.get("radius_raw_record")
        if not _is_map(raw_record) or raw_record.get("operation") != "product" or raw_record.get("arguments") != [85743, 241786]:
            raise ValueError("v2 native radius raw graph record mismatch")
    elif packet_type == "native_history_bundle_source_trace_v2_addendum":
        raise ValueError("v2 native radius raw record is required")
    p0_projection = addendum.get("P0_factorial_once_guard") if _is_map(addendum) else base.get("P0_ordinary_projection")
    if _is_map(p0_projection) and p0_projection.get("P0_numeric_values_absent") is False:
        raise ValueError("numeric P0 values are not accepted")
    return {
        "base": base,
        "addendum": addendum,
        "selectors": selectors,
        "p0": p0,
        "p2": p2,
        "packet_type": packet_type,
        "raw_sha256": raw_sha,
    }


def _validate_radius(packet, raw_sha):
    if raw_sha != _RADIUS_SHA:
        raise ValueError("radius trace raw SHA-256 mismatch")
    if packet.get("schema") != "native_radius_frame_source_trace_v1":
        raise ValueError("radius trace schema mismatch")
    if packet.get("status") != "passed_static_trace_with_rebuilt_radius_prerequisite_open":
        raise ValueError("radius trace status mismatch")
    boundary = packet.get("execution_boundary")
    if not _is_map(boundary):
        raise ValueError("radius execution boundary is missing")
    for key in (
        "source_only",
        "metadata_only",
        "forbidden_resume_or_switch_decode",
        "forbidden_huge_Rm_loop_density_materialization",
        "N_node18_materialized",
        "J_node18_materialized",
        "raw_gzip_hashed_before_decode",
        "raw_gzip_decoded",
    ):
        expected = key in ("source_only", "metadata_only", "raw_gzip_hashed_before_decode")
        if boundary.get(key) is not expected:
            raise ValueError("radius execution boundary mismatch for " + key)
    receipt = packet.get("native_archive_receipt")
    if not _is_map(receipt):
        raise ValueError("native archive receipt is missing")
    if receipt.get("raw_sha256_expected") != _NATIVE_ARCHIVE_SHA or receipt.get("raw_sha256_observed") != _NATIVE_ARCHIVE_SHA or receipt.get("sha_match") is not True:
        raise ValueError("native archive raw SHA-256 mismatch")
    if receipt.get("decode_performed_in_this_trace") is not False:
        raise ValueError("native archive decode is not allowed")
    frame = packet.get("native_chart_frame")
    if not _is_map(frame):
        raise ValueError("radius chart frame is missing")
    if (
        frame.get("chart") != "O3_transition"
        or frame.get("coordinate_node") != 241536
        or frame.get("O3_transition_domain") != [0, 1]
        or frame.get("O3_transition_log_offset") != "17 + exp(40) + x"
        or frame.get("O3_transition_physical_measure") != "1"
        or frame.get("radius_definition") != "R_native = self.Rm * exp(offset)"
    ):
        raise ValueError("radius chart/frame selector mismatch")
    function = packet.get("native_function_definition")
    if not _is_map(function) or function.get("function") != "generic_functions" or function.get("source_sha256") != "7398b8cbba75453c566b7b7281d98827e128e89441124ea67b96d9493ca092a7":
        raise ValueError("native generic function source selector mismatch")
    own = packet.get("rebuilt_own_radius_separation")
    if not _is_map(own) or own.get("native_alias_substitution_allowed") is not False:
        raise ValueError("rebuilt own-radius separation guard is missing")
    return packet


def _validate_parameter(packet, raw_sha):
    if raw_sha != _PARAMETER_SHA:
        raise ValueError("parameter result raw SHA-256 mismatch")
    if packet.get("schema") != "original_native_W_parameter_interface" or packet.get("component") != "native_W_parameter_interface":
        raise ValueError("parameter result schema mismatch")
    if packet.get("passed") is not True or packet.get("source_only") is not True or packet.get("metadata_only") is not True:
        raise ValueError("parameter result admission flags mismatch")
    payload = packet.get("payload")
    if not _is_map(payload):
        raise ValueError("parameter result payload is missing")
    records = payload.get("native_source_handle_records")
    if not isinstance(records, list) or len(records) != 2:
        raise ValueError("exact native parameter handle records are required")
    names = []
    for record in records:
        if not _is_map(record) or record.get("exact_native_handle") is not True or record.get("value_materialized") is not False or record.get("kind") != "native_source_parameter":
            raise ValueError("native parameter handle provenance mismatch")
        names.append(record.get("name"))
    if set(names) != {"Pstar", "delta"}:
        raise ValueError("Pstar and delta native handles are required")
    rebuilt = payload.get("rebuilt_original_records")
    template = payload.get("template_signature_records")
    if not isinstance(rebuilt, list) or not isinstance(template, list) or not template:
        raise ValueError("saved rebuilt/template records are required for refusal checks")
    if not any(_is_map(item) and item.get("kind") == "rebuilt_own_radius" for item in rebuilt):
        raise ValueError("saved rebuilt own-radius record is required")
    return packet


def _validate_source(source):
    if not _is_map(source):
        raise ValueError("source must be a parsed JSON mapping")
    history_packet, history_sha = _get_input_packet(source, "history_trace")
    radius_packet, radius_sha = _get_input_packet(source, "radius_trace")
    parameter_packet, parameter_sha = _get_input_packet(source, "parameter_result")
    history = _validate_history(history_packet, history_sha)
    radius = _validate_radius(radius_packet, radius_sha)
    parameter = _validate_parameter(parameter_packet, parameter_sha)
    return {
        "history": history,
        "radius": radius,
        "parameter": parameter,
        "parent_report_present": source.get("parent_report") is not None or source.get("parentreport") is not None,
        "parent_report_sha256": _PARENT_REPORT_SHA,
    }


class NativeDefinitionProvider:
    """Serve saved native FunctionJet definitions without evaluating them."""

    def __init__(self, validated_source):
        history = validated_source["history"]
        self._base = history["base"]
        self._addendum = history.get("addendum")
        self._selectors = _copy_plain(history["selectors"])
        self._p0 = _copy_plain(history["p0"])
        self._p2 = _copy_plain(history["p2"])
        self._radius_trace = validated_source["radius"]
        self._parameter = validated_source["parameter"]
        self._parameter_handles = {}
        for record in self._parameter["payload"]["native_source_handle_records"]:
            self._parameter_handles[record["name"]] = _copy_plain(record)
        self._frame = self._make_frame()
        self._radius_definition = None
        self._radius_queries = 0
        self._bundle_branches = []

    def _make_frame(self):
        frame = self._radius_trace["native_chart_frame"]
        radius_raw = self._addendum["native_radius_raw_record"]
        return {
            "chart": frame["chart"],
            "coordinate_symbol": frame.get("coordinate_symbol", "true_outer_O3_transition"),
            "coordinate_node": frame["coordinate_node"],
            "coordinate_role": frame.get("coordinate_role", "native radial/log-chart coordinate x; separate from ordinary axial Z"),
            "domain": _copy_plain(frame["O3_transition_domain"]),
            "offset": frame["O3_transition_log_offset"],
            "physical_measure": frame["O3_transition_physical_measure"],
            "radius_node": radius_raw["radius_node"],
            "radius_definition": frame["radius_definition"],
            "native_Rm_owner": frame.get("native_Rm_owner", "self.Rm = self.leading.inlet['radius']"),
            "value_evaluated": False,
        }

    def _field_provenance(self, branch, quantity, order, node):
        handles = self._base.get("native_field_function_handles")
        field = None
        if _is_map(handles) and _is_map(handles.get(branch)):
            fields = handles[branch].get("fields")
            if _is_map(fields):
                field = fields.get(quantity)
        record = {
            "source_chart": "O3_transition",
            "branch_selector": "background" if branch == "background" else "modified",
            "quantity": quantity,
            "ordinary_Z_order": order,
            "node": node,
            "source_row_convention": "saved FunctionJet row selector; ordinary shifted rows; no numeric evaluation",
            "generic_shift_output_evaluated": False,
            "source_row_selector_only": True,
            "source_creator_ast_records": [],
            "source_parameter_provenance": {},
        }
        creators = self._base.get("source_creator_ast_records")
        if isinstance(creators, list):
            for creator in creators:
                if not _is_map(creator):
                    continue
                selected = {"path": creator.get("path"), "source_sha256": creator.get("source_sha256")}
                if creator.get("selectors") is not None:
                    selected["selectors"] = _copy_plain(creator.get("selectors"))
                if creator.get("selector") is not None:
                    selected["selector"] = _copy_plain(creator.get("selector"))
                record["source_creator_ast_records"].append(selected)
        owners = self._base.get("native_parameter_owners")
        if _is_map(owners):
            for parameter_name in ("Pstar", "delta"):
                owner = owners.get(parameter_name)
                if _is_map(owner):
                    record["source_parameter_provenance"][parameter_name] = {
                        "exact_native_handle": owner.get("exact_native_handle"),
                        "value_materialized": owner.get("value_materialized"),
                        "outer_alias": owner.get("outer_alias"),
                    }
        if _is_map(field):
            for key in ("graph_source_pattern", "source_literal", "coordinate_independent", "source_function_not_range_endpoint"):
                if field.get(key) is not None:
                    record[key] = _copy_plain(field.get(key))
        if _is_map(handles) and _is_map(handles.get(branch)):
            auxiliary = handles[branch].get("auxiliary_required_by_generic_source")
            if _is_map(auxiliary):
                auxiliary_key = "E_N" if branch == "modified" and quantity == "E" else quantity
                if auxiliary.get(auxiliary_key) is not None:
                    record["auxiliary_row_nodes"] = _copy_plain(auxiliary.get(auxiliary_key))
        if quantity == "pressure":
            record["definition"] = "pressure=P0+p"
        if quantity == "P0":
            record["source_provider"] = "CurrentOriginalCoreInletAxial6GenericSource.P0_rows"
            record["ordinary_projection_provider"] = "CurrentOriginalCoreInletAxial6GenericSource.project_P0_function"
            record["factorial_once_guard"] = True
        if quantity == "V" and branch == "background":
            record["source_literal"] = "V=0"
            record["zero_source_node"] = 0
            record["zero_is_padding"] = False
        return record

    def _make_definition(self, quantity, branch, order):
        node = self._selectors[branch][quantity][order]
        definition = {
            "kind": "saved_native_graph_function_definition",
            "quantity": quantity,
            "branch": branch,
            "ordinary_Z_order": order,
            "archive_sha256": _NATIVE_ARCHIVE_SHA,
            "root_node": node,
            "row_order": order,
            "frame": _copy_plain(self._frame),
            "creator_provenance": self._field_provenance(branch, quantity, order, node),
            "source_selector": {
                "chart": "O3_transition",
                "branch": branch,
                "quantity": quantity,
                "ordinary_Z_order": order,
                "node": node,
            },
            "value_evaluated": False,
            "metadata_only": True,
            "source_only": True,
        }
        if quantity == "P0":
            definition["ordinary_projection"] = {
                "factorial_already_applied": self._p0["factors"][order],
                "do_not_multiply_again": True,
            }
        return definition

    def get_definition(self, quantity, *, branch, ordinary_Z_order):
        if branch not in ("background", "modified"):
            raise ValueError("branch must be explicit: background or modified")
        if quantity not in _quantities():
            raise KeyError("unknown native quantity: " + str(quantity))
        if type(ordinary_Z_order) is not int or ordinary_Z_order < 0 or ordinary_Z_order > 4:
            raise KeyError("ordinary_Z_order must be one of 0 through 4")
        if quantity not in self._selectors[branch]:
            raise KeyError("missing native quantity selector")
        return _copy_plain(self._make_definition(quantity, branch, ordinary_Z_order))

    def get_radius_definition(self):
        if self._radius_definition is None:
            self._radius_queries += 1
            raw_record = None
            if _is_map(self._addendum):
                raw_record = self._addendum.get("native_radius_raw_record")
            if not _is_map(raw_record):
                raise ValueError("validated v2 native radius raw record is unavailable")
            self._radius_definition = {
                "kind": "saved_native_radius_function_definition",
                "chart": "O3_transition",
                "branch": "shared_native_frame",
                "archive_sha256": _NATIVE_ARCHIVE_SHA,
                "root_node": 241787,
                "row_order": None,
                "frame": _copy_plain(self._frame),
                "creator_provenance": {
                    "owner": "CurrentOriginalTrueOuterC6C5Functions.__init__",
                    "radius_source_literal": "R=Rm*exp(offset)",
                    "radius_owner": "self.Rm=self.leading.inlet['radius']",
                    "raw_record": _copy_plain(raw_record.get("radius_raw_record")),
                },
                "value_evaluated": False,
                "metadata_only": True,
                "source_only": True,
            }
        return _copy_plain(self._radius_definition)

    def get_frame_context(self):
        return _copy_plain(self._frame)

    def _history_context(self):
        finite = self._base["finite_N_literal_source"]
        history = self._base.get("history_normalization_and_Rbase")
        context = {
            "provider_id": _PROVIDER_ID,
            "incoming_Rbase_selector": history.get("incoming_Rbase_selector"),
            "value_evaluated": False,
            "N": {
                "node": finite["N_node"],
                "definition": finite["N_definition"],
                "J_literal": finite["J_literal"],
                "materialized": finite["N_materialized"],
            },
        }
        if _is_map(history):
            for key in ("rates", "transport_rules", "background_history_rows", "finite_N_history_rows", "O2_buffer_endpoint_history_nodes"):
                if history.get(key) is not None:
                    context[key] = _copy_plain(history.get(key))
        finite = self._base.get("finite_N_literal_source")
        if _is_map(finite):
            context["invN_node"] = finite.get("invN_node")
        return context

    def _generic_definition(self):
        generic = self._base.get("native_generic_W_definition")
        if not _is_map(generic):
            generic = {}
        function = self._radius_trace.get("native_function_definition", {})
        return {
            "source_function": _copy_plain(generic),
            "function_definition": _copy_plain(function),
            "instantiated": False,
            "value_evaluated": False,
            "p2_is_quotient_selector": True,
            "p2_is_not_raw_W": True,
        }

    def build_native_bundle(self, branch):
        if branch not in ("background", "modified"):
            raise ValueError("branch must be explicit: background or modified")
        if branch in self._bundle_branches:
            raise ValueError("native bundle branch has already been served")
        fields = {}
        for quantity in _quantities():
            fields[quantity] = [self._make_definition(quantity, branch, order) for order in range(5)]
        p2 = self._base.get("P2_branch_and_inputs")
        if not _is_map(p2):
            raise ValueError("saved P2 route selectors are unavailable")
        p2_handles = p2.get("selected_p2_handles") if branch == "background" else p2.get("separate_modified_generic_p2_handles")
        p2_orders = p2.get("selected_orders") if branch == "background" else [0, 1, 2, 3, 4]
        manifest = {
            "kind": "native_history_radius_bundle_manifest",
            "manifest_id": _COMPONENT + "." + branch,
            "branch": branch,
            "chart": "O3_transition",
            "frame": _copy_plain(self._frame),
            "radius_definition": self.get_radius_definition(),
            "field_definitions": fields,
            "native_parameter_handles": {
                "Pstar": _copy_plain(self._parameter_handles["Pstar"]),
                "delta": _copy_plain(self._parameter_handles["delta"]),
            },
            "generic_definition": self._generic_definition(),
            "history_context": self._history_context(),
            "P2_definition": {
                "root_selectors": _copy_plain(p2_handles),
                "ordinary_Z_orders": _copy_plain(p2_orders),
                "not_raw_W": True,
                "not_a_P2_value": True,
                "route_distinction": _copy_plain(p2.get("route_distinction")) if _is_map(p2) and p2.get("route_distinction") is not None else "saved quotient-root selectors only; no route identity asserted",
            },
            "template_instantiated": False,
            "definition_rows_served": True,
            "template_history_slots_deferred": ["V", "m", "e", "P0", "p"],
            "rebuilt_template_status": "uninstantiated",
            "value_evaluated": False,
            "metadata_only": True,
            "source_only": True,
            "branch_frame_consistent": True,
        }
        self._bundle_branches.append(branch)
        return _copy_plain(manifest)

    def reject_mixed_branch_bundle(self, background_definition, modified_definition):
        if not _is_map(background_definition) or not _is_map(modified_definition):
            raise ValueError("mixed-branch refusal requires saved definition refs")
        if background_definition.get("branch") != "background" or modified_definition.get("branch") != "modified":
            raise ValueError("mixed-branch refusal requires one background and one modified ref")
        return {
            "kind": "mixed_branch_bundle",
            "status": "rejected",
            "reason": "background and modified field definitions cannot share one native branch bundle",
            "background_ref": {
                "quantity": background_definition.get("quantity"),
                "ordinary_Z_order": background_definition.get("ordinary_Z_order"),
                "root_node": background_definition.get("root_node"),
                "branch": background_definition.get("branch"),
            },
            "modified_ref": {
                "quantity": modified_definition.get("quantity"),
                "ordinary_Z_order": modified_definition.get("ordinary_Z_order"),
                "root_node": modified_definition.get("root_node"),
                "branch": modified_definition.get("branch"),
            },
            "value_evaluated": False,
        }

    def reject_rebuilt_frame_binding(self, candidate):
        if not _is_map(candidate):
            raise ValueError("rebuilt-frame refusal requires a saved parameter record")
        rebuilt = candidate.get("rebuilt_radius") if _is_map(candidate.get("rebuilt_radius")) else candidate
        template = candidate.get("template_signature")
        rebuilt_kind = rebuilt.get("kind") if _is_map(rebuilt) else None
        template_status = template.get("status") if _is_map(template) else None
        if rebuilt_kind != "rebuilt_own_radius" and template_status != "uninstantiated":
            raise ValueError("candidate is not the saved rebuilt radius/template record")
        return {
            "kind": "rebuilt_frame_binding",
            "status": "rejected",
            "reason": "rebuiltOwnRadius/template frame is not the native Rm*exp(offset) frame",
            "candidate_kind": rebuilt_kind,
            "template_status": template_status,
            "native_radius_root_node": 241787,
            "native_frame": _copy_plain(self._frame),
            "value_evaluated": False,
        }

    def reject_duplicate_normalization(self, definition, *, apply_factorial):
        if not _is_map(definition) or definition.get("quantity") != "P0" or definition.get("ordinary_Z_order") != 4:
            raise ValueError("duplicate-normalization refusal requires saved P0 order 4 ref")
        if apply_factorial is not True:
            raise ValueError("duplicate-normalization refusal must request the duplicate factor")
        if definition.get("root_node") != 66509:
            raise ValueError("saved P0 order 4 root node mismatch")
        return {
            "kind": "duplicate_normalization",
            "status": "rejected",
            "reason": "P0 ordinary order 4 already carries factorial 24; a second factorial is forbidden",
            "quantity": "P0",
            "ordinary_Z_order": 4,
            "root_node": 66509,
            "factorial_already_applied": 24,
            "value_evaluated": False,
        }


def _failure_result(counts, flags, payload, diagnostics, stage, reason):
    return {
        "schema": _RESULT_SCHEMA,
        "component": _COMPONENT,
        "passed": False,
        "runtime_model_metadata": None,
        "metadata_only": True,
        "source_only": True,
        "physical_layer_completed": False,
        "full_goal_complete": False,
        "actual_counts": _copy_plain(counts),
        "new_flags": _copy_plain(flags),
        "payload": _copy_plain(payload),
        "diagnostics": _copy_plain(diagnostics),
        "failure": {
            "stage": stage,
            "reason": reason,
            "finished_field_rows": len(payload.get("definitions", [])),
            "finished_radius_rows": 1 if payload.get("radius_definition") is not None else 0,
            "finished_bundle_manifests": len(payload.get("native_bundle_manifests", [])),
            "actual_counts": _copy_plain(counts),
            "new_flags": _copy_plain(flags),
        },
    }


def derive(source, diagnostic=None):
    """Return the frozen native definition bundle as plain JSON data."""
    counts = _zero_counts()
    flags = _false_flags()
    diagnostics = []
    definitions = []
    radius_definition = None
    manifests = []
    rejections = []
    payload = {
        "definitions": definitions,
        "radius_definition": radius_definition,
        "native_bundle_manifests": manifests,
        "rejections": rejections,
        "provider_entry": None,
        "claim_guards": {
            "W_evaluated": False,
            "P2_value_evaluated": False,
            "zero_identity_claimed": False,
            "native_rebuilt_equality_claimed": False,
            "error_cap_reduction_claimed": False,
            "physical_radius_equality_claimed": False,
            "global_claimed": False,
            "temporal_recursion_claimed": False,
        },
        "source_fingerprints": {
            "history_trace_raw_sha256": None,
            "radius_trace_raw_sha256": None,
            "parameter_result_raw_sha256": None,
            "parent_report_sha256": None,
            "parent_report_consumed": False,
        },
    }
    if diagnostic is not None and not callable(diagnostic):
        return _failure_result(counts, flags, payload, diagnostics, "input", "diagnostic must be callable or None")

    def emit(event, record):
        item = {"event": event}
        item.update(_copy_plain(record))
        diagnostics.append(item)
        if diagnostic is not None:
            try:
                diagnostic(item)
            except Exception as exc:
                item["callback_error"] = {"type": type(exc).__name__, "message": str(exc)}
                return False
        return True

    validation = None
    validation_error = None
    try:
        validation = _validate_source(source)
        payload["source_fingerprints"]["history_trace_raw_sha256"] = validation["history"]["raw_sha256"]
        payload["source_fingerprints"]["radius_trace_raw_sha256"] = _RADIUS_SHA
        payload["source_fingerprints"]["parameter_result_raw_sha256"] = _PARAMETER_SHA
        if validation["parent_report_present"]:
            payload["source_fingerprints"]["parent_report_sha256"] = _PARENT_REPORT_SHA
    except Exception as exc:
        validation_error = str(exc)
    validation_record = {
        "status": "passed" if validation_error is None else "rejected",
        "input_schema": _INPUT_SCHEMA,
        "history_trace_schema": validation["history"]["packet_type"] if validation is not None else None,
        "radius_trace_schema": validation["radius"].get("schema") if validation is not None else None,
        "parameter_result_schema": validation["parameter"].get("schema") if validation is not None else None,
        "raw_sha256_pins_checked": True,
        "typed_selectors_checked": validation_error is None,
    }
    if not emit("source_validation", validation_record):
        return _failure_result(counts, flags, payload, diagnostics, "source_validation", "diagnostic callback failed")
    if validation_error is not None:
        return _failure_result(counts, flags, payload, diagnostics, "source_validation", validation_error)

    counts["definition_provider_entries"] = 1
    payload["provider_entry"] = {
        "provider_id": _PROVIDER_ID,
        "input_schema": _INPUT_SCHEMA,
        "status": "registered_before_construction",
        "metadata_only": True,
        "source_only": True,
    }
    if not emit("definition_provider_entry", payload["provider_entry"]):
        counts["definition_provider_entries"] = 1
        return _failure_result(counts, flags, payload, diagnostics, "definition_provider_entry", "diagnostic callback failed")

    try:
        provider = NativeDefinitionProvider(validation)
    except Exception as exc:
        return _failure_result(counts, flags, payload, diagnostics, "provider_construction", str(exc))

    for branch in ("background", "modified"):
        for quantity in _quantities():
            for order in range(5):
                try:
                    definition = provider.get_definition(quantity, branch=branch, ordinary_Z_order=order)
                except Exception as exc:
                    return _failure_result(counts, flags, payload, diagnostics, "field_definition", str(exc))
                definitions.append(definition)
                counts["field_definition_queries"] += 1
                if not emit(
                    "field_definition_row",
                    {
                        "status": "served",
                        "branch": branch,
                        "quantity": quantity,
                        "ordinary_Z_order": order,
                        "definition": definition,
                    },
                ):
                    return _failure_result(counts, flags, payload, diagnostics, "field_definition_row", "diagnostic callback failed")

    try:
        radius_definition = provider.get_radius_definition()
    except Exception as exc:
        return _failure_result(counts, flags, payload, diagnostics, "radius_definition", str(exc))
    payload["radius_definition"] = radius_definition
    counts["radius_definition_queries"] = 1
    if not emit("radius_definition_row", {"status": "served", "definition": radius_definition}):
        return _failure_result(counts, flags, payload, diagnostics, "radius_definition_row", "diagnostic callback failed")

    for branch in ("background", "modified"):
        try:
            manifest = provider.build_native_bundle(branch)
        except Exception as exc:
            return _failure_result(counts, flags, payload, diagnostics, "native_bundle_manifest", str(exc))
        manifests.append(manifest)
        counts["native_branches_served"] += 1
        counts["native_bundle_manifests"] += 1
        if not emit("native_bundle_manifest", {"status": "served", "branch": branch, "manifest": manifest}):
            return _failure_result(counts, flags, payload, diagnostics, "native_bundle_manifest", "diagnostic callback failed")

    by_key = {(item["branch"], item["quantity"], item["ordinary_Z_order"]): item for item in definitions}
    try:
        mixed = provider.reject_mixed_branch_bundle(by_key[("background", "m", 0)], by_key[("modified", "V", 0)])
        rejections.append(mixed)
        counts["mixed_branch_rejections"] = 1
        if not emit("context_rejection", {"status": "rejected", "rejection": mixed}):
            return _failure_result(counts, flags, payload, diagnostics, "mixed_branch_rejection", "diagnostic callback failed")

        rebuilt_records = validation["parameter"]["payload"]["rebuilt_original_records"]
        rebuilt_radius = next(item for item in rebuilt_records if item.get("kind") == "rebuilt_own_radius")
        template = validation["parameter"]["payload"]["template_signature_records"][0]
        frame_rejection = provider.reject_rebuilt_frame_binding({"rebuilt_radius": rebuilt_radius, "template_signature": template})
        rejections.append(frame_rejection)
        counts["rebuilt_frame_rejections"] = 1
        if not emit("context_rejection", {"status": "rejected", "rejection": frame_rejection}):
            return _failure_result(counts, flags, payload, diagnostics, "rebuilt_frame_rejection", "diagnostic callback failed")

        duplicate = provider.reject_duplicate_normalization(by_key[("background", "P0", 4)], apply_factorial=True)
        rejections.append(duplicate)
        counts["duplicate_normalization_rejections"] = 1
        if not emit("context_rejection", {"status": "rejected", "rejection": duplicate}):
            return _failure_result(counts, flags, payload, diagnostics, "duplicate_normalization_rejection", "diagnostic callback failed")
    except Exception as exc:
        return _failure_result(counts, flags, payload, diagnostics, "context_refusal", str(exc))

    flags["native_history_radius_definitions_served"] = True
    flags["native_branch_frame_separation_enforced"] = True
    flags["native_ordinary_Z_projection_preserved"] = True
    return {
        "schema": _RESULT_SCHEMA,
        "component": _COMPONENT,
        "passed": True,
        "runtime_model_metadata": None,
        "metadata_only": True,
        "source_only": True,
        "physical_layer_completed": False,
        "full_goal_complete": False,
        "actual_counts": _copy_plain(counts),
        "new_flags": _copy_plain(flags),
        "payload": _copy_plain(payload),
        "diagnostics": _copy_plain(diagnostics),
        "failure": None,
    }
