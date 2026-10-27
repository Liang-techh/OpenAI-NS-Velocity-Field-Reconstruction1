"""Independent, unexecuted C3 endpoint-query projection kernel.

This module consumes authenticated result objects as data.  It resolves the
C3 right-source selectors from the authenticated parent result, reuses the
already-passed cached row records from the saved failed quiet result, and
composes provenance-only endpoint relations.  It does not import, replay, or
execute any ancestor quiet, native, algebra, provider, proof, checker, or
normalizer code.
"""

import hashlib


PARENT_RESULT_SHA256 = (
    "b91da155898c38a9ca4b26595983a38aa04b2965ae3db0a908fadcea97b34241"
)
FAILED_RESULT_SHA256 = (
    "e9b1711e64408b1d016d2a0e73adc55deef8db64e70c76acbc9f0b111b97790d"
)
C3_ARCHIVE_SHA256 = (
    "cf51e2131c3d1cb76f93a2deb76af12421e2deb75f7888bd310084b3b3eec193"
)
FIELDS = ("m", "h", "k", "e", "p")
ORDERS = (0, 1, 2, 3)
EXPECTED_TERMINAL_PORTS = {
    "m": [4920, 4921, 12042, 24177],
    "h": [4928, 4929, 12043, 24178],
    "k": [4930, 4931, 12044, 24179],
    "e": [4932, 4933, 12045, 24180],
    "p": [4888, 4889, 12041, 24176],
}
EXPECTED_MEMORY_NODES = {"m": 4919, "e": 4919, "h": 4927, "k": 4927, "p": 1}
EXPECTED_SOURCE_FAMILY = {
    "actual_five_defect_family_sha256": (
        "3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894"
    ),
    "datum_enclosure_sha256": (
        "dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010"
    ),
    "implicit_source_sha256": (
        "5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae"
    ),
}

INPUT_CONTRACT = {
    "source": {
        "parent_result": "authenticated parent result object",
        "failed_result": "saved failed quiet result object",
        "frozen_input_slice": "C3 literal slice object",
        "parent_result_raw_sha256": PARENT_RESULT_SHA256,
        "failed_result_raw_sha256": FAILED_RESULT_SHA256,
        "frozen_input_slice.archive_raw_sha256": C3_ARCHIVE_SHA256,
    },
    "selectors": {
        "parent_ports": "parent_result.C3_rebuilt_history_source_ports.terminal_ports",
        "parent_queries": "parent_result.signed_history_gap_function_nodes[*].right.source_node",
        "cached_terminal_rows": "failed_result.C3_source_binding.terminal_row_checks",
        "cached_memory_rows": "failed_result.C3_source_binding.memory_row_checks",
        "frozen_inlets": "frozen_input_slice.incoming_ports",
    },
    "scope": "cached-row reuse; no new density theorem",
    "ordinary_Z_orders": [0, 1, 2, 3],
    "Z_domain": ["-1", "1"],
    "function_not_samples": True,
    "numeric": False,
    "ancestor_execution": False,
}


def qep_journal(journal, event):
    if isinstance(journal, list):
        journal.append(event)


def _sha256_text(value):
    if isinstance(value, bytes):
        return hashlib.sha256(value).hexdigest()
    if isinstance(value, str):
        return hashlib.sha256(value.encode("utf-8")).hexdigest()
    return None


def _first_hash(source, names):
    if not isinstance(source, dict):
        return None
    for name in names:
        value = source.get(name)
        if isinstance(value, str):
            return value.lower()
    return None


def _object_hash(source, object_names):
    if not isinstance(source, dict):
        return None
    for name in object_names:
        value = source.get(name)
        if isinstance(value, dict):
            for key in ("raw_sha256", "_raw_sha256", "sha256"):
                declared = value.get(key)
                if isinstance(declared, str):
                    return declared.lower()
    return None


def _nested_named_hash(source, container_names, value_names):
    if not isinstance(source, dict):
        return None
    for container_name in container_names:
        container = source.get(container_name)
        if not isinstance(container, dict):
            continue
        for value_name in value_names:
            value = container.get(value_name)
            if isinstance(value, str):
                return value.lower()
    return None


def _select_mapping(source, names):
    if not isinstance(source, dict):
        return None
    for name in names:
        value = source.get(name)
        if isinstance(value, dict):
            return value
    return None


def _contains_forbidden_ref(value):
    if not isinstance(value, str):
        return False
    lowered = value.lower()
    return "native" in lowered or "primary" in lowered


def _validate_hash_pins(source, parent, failed, frozen, blockers):
    parent_sha = _first_hash(
        source,
        (
            "parent_result_raw_sha256",
            "parent_raw_sha256",
            "parent_result_sha256",
        ),
    ) or _nested_named_hash(
        source,
        ("raw_sha256", "raw_sha_pins", "hashes"),
        ("parent_result", "parent", "parent_result_raw_sha256"),
    ) or _object_hash(source, ("parent_result_metadata", "parent_metadata"))
    if parent_sha is None and isinstance(parent, dict):
        parent_sha = _first_hash(parent, ("_raw_sha256", "raw_sha256", "sha256"))
    if parent_sha != PARENT_RESULT_SHA256:
        blockers.append("parent result raw SHA pin is not the authenticated b91da result")

    failed_sha = _first_hash(
        source,
        (
            "failed_result_raw_sha256",
            "failed_saved_result_raw_sha256",
            "failed_saved_result_sha256",
            "saved_failed_result_raw_sha256",
            "failed_result_sha256",
        ),
    ) or _nested_named_hash(
        source,
        ("raw_sha256", "raw_sha_pins", "hashes"),
        (
            "failed_result",
            "failed_saved_result",
            "saved_failed_result",
            "failed_result_raw_sha256",
        ),
    ) or _object_hash(
        source,
        ("failed_result_metadata", "failed_saved_result_metadata", "failed_metadata"),
    )
    if failed_sha is None and isinstance(failed, dict):
        failed_sha = _first_hash(failed, ("_raw_sha256", "raw_sha256", "sha256"))
    if failed_sha != FAILED_RESULT_SHA256:
        blockers.append("failed result raw SHA pin is not the saved e9b171 result")

    if not isinstance(frozen, dict):
        blockers.append("frozen input slice is missing")
    elif str(frozen.get("archive_raw_sha256", "")).lower() != C3_ARCHIVE_SHA256:
        blockers.append("frozen C3 archive SHA pin is not authenticated")
    return {
        "parent_result_raw_sha256": parent_sha,
        "failed_result_raw_sha256": failed_sha,
        "frozen_archive_raw_sha256": (
            frozen.get("archive_raw_sha256") if isinstance(frozen, dict) else None
        ),
    }


def _validate_frozen_slice(frozen, blockers):
    if not isinstance(frozen, dict):
        return None, None
    if frozen.get("primary_alias_used") is not False:
        blockers.append("frozen C3 slice permits a primary alias")
    if frozen.get("primary_snapshot_compared") is not False:
        blockers.append("frozen C3 slice permits a primary snapshot comparison")
    family = frozen.get("source_family")
    if family != EXPECTED_SOURCE_FAMILY:
        blockers.append("frozen C3 source family is not the authenticated family")
    nodes = frozen.get("nodes")
    if not isinstance(nodes, dict):
        blockers.append("frozen C3 literal node map is missing")
    incoming = frozen.get("incoming_ports")
    if not isinstance(incoming, dict):
        blockers.append("frozen C3 incoming port map is missing")
        incoming = {}
    for field in FIELDS:
        values = incoming.get(field)
        if not isinstance(values, list) or len(values) != len(ORDERS):
            blockers.append("frozen C3 incoming port domain is incomplete: %s" % field)
            continue
        if len(set(values)) != len(values):
            blockers.append("frozen C3 incoming ports are not unique: %s" % field)
        if any(not isinstance(value, int) or isinstance(value, bool) for value in values):
            blockers.append("frozen C3 incoming ports are not integer selectors: %s" % field)
    terminal = frozen.get("terminal_ports")
    if not isinstance(terminal, dict):
        blockers.append("frozen C3 terminal port map is missing")
    return incoming, nodes


def _validate_parent_ports(parent, frozen_terminal, blockers):
    metadata = parent.get("C3_rebuilt_history_source_ports") if isinstance(parent, dict) else None
    if not isinstance(metadata, dict):
        blockers.append("parent C3 rebuilt history source port metadata is missing")
        return None
    if metadata.get("archive_sha256") != C3_ARCHIVE_SHA256:
        blockers.append("parent C3 source port archive SHA is not the C3 archive")
    if metadata.get("namespace") != "repair_C3_band":
        blockers.append("parent C3 source port namespace is not repair_C3_band")
    if metadata.get("ordinary_Z_orders") != list(ORDERS):
        blockers.append("parent C3 ordinary Z domain is not [0, 1, 2, 3]")
    if metadata.get("range_values_not_function_values") is not True:
        blockers.append("parent C3 port range/function distinction is missing")
    for key in ("namespace", "role", "source_role"):
        if _contains_forbidden_ref(metadata.get(key)):
            blockers.append("parent C3 source port metadata contains native/primary reference")
    ports = metadata.get("terminal_ports")
    if not isinstance(ports, dict):
        blockers.append("parent C3 terminal port selector map is missing")
        return None
    seen = set()
    for field in FIELDS:
        values = ports.get(field)
        if values != EXPECTED_TERMINAL_PORTS[field]:
            blockers.append("parent C3 terminal port selector mismatch: %s" % field)
            continue
        if any(not isinstance(value, int) or isinstance(value, bool) for value in values):
            blockers.append("parent C3 terminal ports are not integer selectors: %s" % field)
        if len(set(values)) != len(values):
            blockers.append("parent C3 terminal port selectors are not unique: %s" % field)
        for value in values:
            if value in seen:
                blockers.append("parent C3 terminal port selector is duplicated")
            seen.add(value)
    if isinstance(frozen_terminal, dict):
        for field in FIELDS:
            if frozen_terminal.get(field) != ports.get(field):
                blockers.append("frozen and parent C3 terminal selectors differ: %s" % field)
    return ports


def _project_parent_queries(parent, ports, blockers):
    rows = parent.get("signed_history_gap_function_nodes") if isinstance(parent, dict) else None
    if not isinstance(rows, list) or len(rows) != 20:
        blockers.append("parent C3 signed history query list is not exactly 20 rows")
        return []
    queries = []
    seen = set()
    for index, row in enumerate(rows):
        if not isinstance(row, dict):
            blockers.append("parent C3 signed history query row is not an object")
            continue
        field = row.get("field")
        order = row.get("ordinary_Z_order")
        right = row.get("right")
        if (
            field not in FIELDS
            or not isinstance(order, int)
            or isinstance(order, bool)
            or order not in ORDERS
            or not isinstance(right, dict)
        ):
            blockers.append("parent C3 query row has an invalid ordinary-Z role")
            continue
        key = (field, order)
        if key in seen:
            blockers.append("parent C3 query role is duplicated: %s/%s" % key)
            continue
        seen.add(key)
        source_node = right.get("source_node")
        if not isinstance(source_node, int) or isinstance(source_node, bool):
            blockers.append("parent C3 right source node is not an integer")
            continue
        port_values = ports.get(field) if isinstance(ports, dict) else None
        if not isinstance(port_values, list) or order >= len(port_values):
            blockers.append("parent C3 terminal port selector is unavailable for query role")
            continue
        if source_node != port_values[order]:
            blockers.append("parent C3 right source node disagrees with terminal_ports")
        namespace = right.get("namespace")
        role = right.get("role")
        source_role = right.get("source_role")
        if namespace != "repair_C3_band":
            blockers.append("parent C3 query right namespace is not repair_C3_band")
        if role != "C3_rebuilt_terminal_N_scaled_history":
            blockers.append("parent C3 query right role is not the C3 terminal role")
        if not isinstance(source_role, str) or "C3" not in source_role:
            blockers.append("parent C3 query right source role is not C3")
        if any(_contains_forbidden_ref(value) for value in (namespace, role, source_role)):
            blockers.append("parent C3 query right contains a native/primary reference")
        queries.append(
            {
                "field": field,
                "ordinary_Z_order": order,
                "terminal_node": source_node,
                "parent_query_index": index,
                "parent_terminal_selector": (
                    "parent.C3_rebuilt_history_source_ports.terminal_ports.%s[%s]"
                    % (field, order)
                ),
                "parent_query_selector": (
                    "parent.signed_history_gap_function_nodes[%s].right.source_node" % index
                ),
                "namespace": namespace,
                "role": role,
                "source_role": source_role,
            }
        )
    if seen != {(field, order) for field in FIELDS for order in ORDERS}:
        blockers.append("parent C3 query roles do not cover ordinary Z 0..3 for all five fields")
    return queries


def _cached_passed_rows(failed, blockers):
    binding = failed.get("C3_source_binding") if isinstance(failed, dict) else None
    if not isinstance(binding, dict):
        blockers.append("saved failed result C3_source_binding is missing")
        return None, None, None
    if failed.get("passed") is not False:
        blockers.append("saved failed result passed flag is not preserved as false")
    if binding.get("passed") is not False:
        blockers.append("saved failed C3_source_binding passed flag is not false")
    if binding.get("entry_gap_status", "unknown") != "unknown":
        blockers.append("saved failed C3 entry gap is not preserved as unknown")
    if binding.get("entry_gap_zero_or_nonzero_asserted") is not False:
        blockers.append("saved failed C3 entry gap has a zero/nonzero assertion")

    terminal_rows = binding.get("terminal_row_checks")
    if not isinstance(terminal_rows, list) or len(terminal_rows) != 20:
        blockers.append("saved failed C3 terminal cache is not exactly 20 rows")
        terminal_rows = []
    terminal_by_key = {}
    for row in terminal_rows:
        if not isinstance(row, dict) or row.get("passed") is not True:
            blockers.append("saved failed C3 terminal cache contains an unpassed row")
            continue
        field = row.get("field")
        node = row.get("node")
        key = (field, node)
        if field not in FIELDS or not isinstance(node, int) or key in terminal_by_key:
            blockers.append("saved failed C3 terminal cache has an invalid/duplicate selector")
            continue
        terminal_by_key[key] = row

    memory_rows = binding.get("memory_row_checks")
    if not isinstance(memory_rows, list) or len(memory_rows) != 6:
        blockers.append("saved failed C3 memory cache is not exactly six rows")
        memory_rows = []
    memory_by_node = {}
    for row in memory_rows:
        if not isinstance(row, dict) or row.get("passed") is not True:
            blockers.append("saved failed C3 memory cache contains an unpassed row")
            continue
        node = row.get("node")
        if not isinstance(node, int) or node in memory_by_node:
            blockers.append("saved failed C3 memory cache has an invalid/duplicate node")
            continue
        memory_by_node[node] = row
    return binding, terminal_by_key, memory_by_node


def _memory_sources(binding, memory_by_node, nodes, blockers):
    actual_memory = binding.get("actual_memory") if isinstance(binding, dict) else None
    if not isinstance(actual_memory, dict):
        blockers.append("saved failed C3 actual_memory metadata is missing")
        actual_memory = {}
    sources = {}
    for field in FIELDS:
        memory_node = EXPECTED_MEMORY_NODES[field]
        descriptor = actual_memory.get(field)
        if not isinstance(descriptor, dict):
            blockers.append("saved failed C3 memory descriptor is missing: %s" % field)
        elif descriptor.get("memory") != memory_node:
            blockers.append("saved failed C3 memory descriptor disagrees: %s" % field)
        cached = memory_by_node.get(memory_node)
        if cached is not None:
            sources[field] = {
                "memory_node": memory_node,
                "selector": (
                    "failed.C3_source_binding.memory_row_checks[node=%s]" % memory_node
                ),
                "cached_passed_row": cached,
                "descriptor": descriptor,
            }
            continue
        if field == "p" and memory_node == 1 and isinstance(nodes, dict):
            literal_row = nodes.get("1", nodes.get(1))
            if literal_row is not None:
                sources[field] = {
                    "memory_node": memory_node,
                    "selector": (
                        "failed.C3_source_binding.actual_memory.p + "
                        "frozen_input_slice.nodes['1']"
                    ),
                    "cached_passed_row": None,
                    "frozen_literal_row": literal_row,
                    "descriptor": descriptor,
                }
                continue
        blockers.append("cached C3 memory source is missing: %s" % field)
    return sources


def _cached_observed_arguments(cached_row, field, order, blockers):
    observed = cached_row.get("observed") if isinstance(cached_row, dict) else None
    arguments = observed.get("arguments") if isinstance(observed, dict) else None
    if not isinstance(arguments, list):
        blockers.append(
            "cached C3 observed arguments are missing for %s/%s" % (field, order)
        )
        return None
    if any(not isinstance(value, int) or isinstance(value, bool) for value in arguments):
        blockers.append(
            "cached C3 observed arguments are not integer selectors for %s/%s"
            % (field, order)
        )
        return None
    return arguments


def _compose_projection(queries, terminal_by_key, incoming, memory_sources, blockers):
    typed = []
    recipes = []
    for query in queries:
        field = query["field"]
        order = query["ordinary_Z_order"]
        key = (field, query["terminal_node"])
        cached_row = terminal_by_key.get(key)
        if cached_row is None:
            blockers.append("no passed cached terminal row for parent C3 query %s/%s" % (field, order))
            continue
        inlet_values = incoming.get(field) if isinstance(incoming, dict) else None
        inlet_node = (
            inlet_values[order]
            if isinstance(inlet_values, list) and order < len(inlet_values)
            else None
        )
        if not isinstance(inlet_node, int) or isinstance(inlet_node, bool):
            blockers.append("frozen C3 inlet selector is missing for %s/%s" % (field, order))
            continue
        memory = memory_sources.get(field)
        if memory is None:
            blockers.append("memory provenance is missing for %s" % field)
            continue
        observed_arguments = _cached_observed_arguments(
            cached_row, field, order, blockers
        )
        if observed_arguments is None:
            continue
        if field == "p":
            if query["terminal_node"] != inlet_node:
                blockers.append("p C3 endpoint is not an in/out alias for %s/%s" % (field, order))
                continue
            selected_inlet = query["terminal_node"]
            selected_memory_argument = None
            alias = True
        else:
            if len(observed_arguments) != 2:
                blockers.append(
                    "cached C3 product arguments are not a memory/inlet pair for %s/%s"
                    % (field, order)
                )
                continue
            selected_memory_argument = observed_arguments[0]
            selected_inlet = observed_arguments[1]
            if selected_memory_argument != memory["memory_node"]:
                blockers.append(
                    "cached C3 memory argument disagrees with descriptor for %s" % field
                )
                continue
            if selected_inlet != inlet_node:
                blockers.append(
                    "cached C3 inlet argument disagrees with frozen inlet for %s/%s"
                    % (field, order)
                )
                continue
            alias = False
        cache_arg_provenance = {
            "selector": (
                "failed.C3_source_binding.terminal_row_checks[field=%s,node=%s].observed.arguments"
                % (field, query["terminal_node"])
            ),
            "observed_arguments": observed_arguments,
            "selected_memory_argument": selected_memory_argument,
            "selected_inlet_argument": selected_inlet,
            "frozen_inlet_selector": (
                "frozen_input_slice.incoming_ports.%s[%s]" % (field, order)
            ),
            "frozen_inlet_node": inlet_node,
            "inlet_linked": selected_inlet == inlet_node,
            "p_in_out_alias": alias,
        }
        typed_ref = {
            "type": "C3_terminal_function_ref",
            "field": field,
            "ordinary_Z_order": order,
            "ordinary_Z_orders": list(ORDERS),
            "Z_domain": ["-1", "1"],
            "function_not_samples": True,
            "namespace": query["namespace"],
            "role": query["role"],
            "source_role": query["source_role"],
            "terminal_node": query["terminal_node"],
            "in_node": selected_inlet,
            "memory_node": memory["memory_node"],
            "parent_terminal_selector": query["parent_terminal_selector"],
            "parent_query_selector": query["parent_query_selector"],
            "cached_row_selector": (
                "failed.C3_source_binding.terminal_row_checks[field=%s,node=%s]"
                % (field, query["terminal_node"])
            ),
            "cached_passed_row": cached_row,
            "selected_cache_arguments": cache_arg_provenance,
        }
        typed.append(typed_ref)
        recipes.append(
            {
                "type": "C3_endpoint_memory_relation",
                "field": field,
                "ordinary_Z_order": order,
                "ordinary_Z_orders": list(ORDERS),
                "Z_domain": ["-1", "1"],
                "function_not_samples": True,
                "formula": "out_i^j = M_i * in_i^j",
                "out_node": query["terminal_node"],
                "in_node": selected_inlet,
                "memory_node": memory["memory_node"],
                "p_in_out_alias": alias,
                "out_selector": query["parent_query_selector"],
                "in_selector": (
                    "failed.C3_source_binding.terminal_row_checks[field=%s,node=%s].observed.arguments"
                    % (field, query["terminal_node"])
                ),
                "frozen_inlet_selector": (
                    "frozen_input_slice.incoming_ports.%s[%s]" % (field, order)
                ),
                "memory_selector": memory["selector"],
                "terminal_query_ref": typed_ref,
                "memory_provenance": memory,
                "selected_cache_arguments": cache_arg_provenance,
                "numeric": False,
                "source_only": True,
            }
        )
    return typed, recipes


def _result(
    source,
    blockers,
    pins=None,
    failed=None,
    failed_binding=None,
    typed=None,
    recipes=None,
    journal_events=None,
):
    typed = typed or []
    recipes = recipes or []
    failed_binding = failed_binding if isinstance(failed_binding, dict) else {}
    failed_passed = failed.get("passed") if isinstance(failed, dict) else None
    entry_status = failed_binding.get("entry_gap_status", "unknown")
    return {
        "passed": not blockers,
        "source_only": True,
        "numeric": False,
        "numeric_values_computed": False,
        "metadata_only": False,
        "C3_quiet_terminal_query_projection_bound": not blockers,
        "C3_quiet_endpoint_memory_relations_reused": not blockers,
        "typed_terminal_query_count": len(typed),
        "transport_recipe_count": len(recipes),
        "typed_terminal_queries": typed,
        "transport_recipes": recipes,
        "scope": "cached-row reuse; no new density theorem",
        "nonclaims": [
            "C3 quiet density zero",
            "native C3 execution",
            "common-kernel identity",
            "full quiet pair",
        ],
        "saved_failed_result_passed": failed_passed,
        "saved_failed_result_preserved_false": failed_passed is False,
        "saved_failed_C3_binding_passed": failed_binding.get("passed"),
        "entry_gap_status": entry_status,
        "entry_gap_zero_or_nonzero_asserted": failed_binding.get(
            "entry_gap_zero_or_nonzero_asserted", False
        ),
        "raw_sha_pins": pins or {},
        "input_contract": INPUT_CONTRACT,
        "journal_events": journal_events or [],
        "blockers": list(blockers),
    }


def derive(source, journal=None):
    """Project authenticated parent C3 endpoint queries onto cached passed rows.

    This function is intentionally independent of the entered quiet family. It
    consumes only result and literal-slice data under fixed raw SHA pins; it
    performs no provider, proof, algebra, density, numeric, or ancestor call.
    """
    blockers = []
    events = []
    if not isinstance(source, dict):
        return _result(source, ["source must be a mapping"], journal_events=events)

    parent = _select_mapping(source, ("parent_result", "parent"))
    failed = _select_mapping(
        source, ("failed_result", "failed_saved_result", "saved_failed_result", "failed")
    )
    frozen = _select_mapping(source, ("frozen_input_slice", "C3_literal_slice"))
    pins = _validate_hash_pins(source, parent, failed, frozen, blockers)
    qep_journal(journal, {"stage": "fixed_raw_sha_pins", "passed": not blockers})
    events.append({"stage": "fixed_raw_sha_pins", "passed": not blockers})
    if parent is None or failed is None or frozen is None:
        blockers.append("parent result, failed result, and frozen input slice are required")
        return _result(source, blockers, pins=pins, failed=failed, journal_events=events)

    incoming, nodes = _validate_frozen_slice(frozen, blockers)
    frozen_terminal = frozen.get("terminal_ports")
    ports = _validate_parent_ports(parent, frozen_terminal, blockers)
    qep_journal(
        journal,
        {
            "stage": "c3_selector_domains",
            "passed": not blockers,
            "selector": INPUT_CONTRACT["selectors"]["parent_queries"],
        },
    )
    events.append({"stage": "c3_selector_domains", "passed": not blockers})
    queries = _project_parent_queries(parent, ports or {}, blockers)
    binding, terminal_by_key, memory_by_node = _cached_passed_rows(failed, blockers)
    if binding is None:
        return _result(source, blockers, pins=pins, failed=failed, journal_events=events)
    memory_sources = _memory_sources(binding, memory_by_node, nodes, blockers)
    typed, recipes = _compose_projection(
        queries, terminal_by_key, incoming or {}, memory_sources, blockers
    )
    qep_journal(
        journal,
        {
            "stage": "cached_row_projection",
            "passed": not blockers,
            "typed_terminal_queries": len(typed),
            "transport_recipes": len(recipes),
        },
    )
    events.append(
        {
            "stage": "cached_row_projection",
            "passed": not blockers,
            "typed_terminal_queries": len(typed),
            "transport_recipes": len(recipes),
        }
    )
    return _result(
        source,
        blockers,
        pins=pins,
        failed=failed,
        failed_binding=binding,
        typed=typed,
        recipes=recipes,
        journal_events=events,
    )


__all__ = ["derive"]


