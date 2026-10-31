"""Independent synthetic-object fixture for the live Q4 source adapter.

The fixture does not import the adapter or any original producer.  It passes a
callable supplied by the owner, builds only tiny opaque graph/reference
objects, and returns JSON-safe checks and diagnostics.  Negative cases ensure
that missing history, a changed graph, a changed source-chart object, and
same-label/different-object data remain rejected or unresolved.
"""

from __future__ import annotations


SCHEMA = "original_q4_lower_live_bundle_fixture_v1"


class _Graph:
    pass


class _Ref:
    def __init__(self, graph, node):
        self.graph = graph
        self.node = node


class _Jet:
    def __init__(self, *rows):
        self.rows = tuple(rows)


class _Owner:
    pass


def _ref(graph, label):
    return _Ref(graph, label)


def _jet(graph, label):
    return _Jet(_ref(graph, label + ":Z0"), _ref(graph, label + ":Z1"),
                 _ref(graph, label + ":Z2"))


def _make_case(*, missing=None, swapped_graph=False, swapped_source=False,
               remove_history=False, same_label_source=False):
    graph = _Graph()
    other_graph = _Graph()
    background = _Owner()
    background.delta = _ref(graph, "delta")
    background.Pstar = _ref(graph, "Pstar")
    background.errors = {"background": "retained"}

    lower = _ref(graph, "first_micro:lower")
    upper = _ref(graph, "first_micro:upper")
    coordinate = _ref(graph, "coordinate")
    radius = _ref(graph, "radius")
    modified_fields = {
        "V": _jet(graph, "V"),
        "m": _jet(graph, "m"),
        "pressure": _jet(graph, "pressure"),
    }
    source_chart = {
        "coordinate": coordinate,
        "domain": [lower, upper],
        "radius": radius,
        "modified_fields": modified_fields,
    }
    loop_first = {
        "domain": [lower, upper],
        "source_chart": {"source": "bridge"},
        "incoming": {"m": _jet(graph, "incoming_m"), "p": _jet(graph, "incoming_p")},
        "history": {"m": _jet(graph, "history_m"), "p": _jet(graph, "history_p")},
    }
    if remove_history:
        del loop_first["history"]
    if swapped_source:
        # Deliberately retain equal graph labels in this replacement.  The
        # source-chart object identity must still be required.
        source_chart = {
            "coordinate": _ref(graph, "coordinate"),
            "domain": [_ref(graph, "first_micro:lower"), _ref(graph, "first_micro:upper")],
            "radius": _ref(graph, "radius"),
            "modified_fields": {"V": _jet(graph, "V")},
        }
    if same_label_source:
        source_chart = {
            "coordinate": _ref(graph, "coordinate"),
            "domain": [_ref(graph, "first_micro:lower"), _ref(graph, "first_micro:upper")],
            "radius": _ref(graph, "radius"),
            "modified_fields": modified_fields,
        }

    physical_graph = other_graph if swapped_graph else graph
    physical = {
        "U": _jet(physical_graph, "U"),
        "m": _jet(graph, "physical_m"),
        "U_y": _jet(graph, "U_y"),
        "Q_y": _jet(graph, "Q_y"),
    }
    finite_boundary = {"R_Q_y": _jet(graph, "finite_R_Q_y")}
    initial_R_Q_y = _Jet(_ref(graph, "initial:R_Q_y:Z0"),
                         _ref(graph, "initial:R_Q_y:Z1"),
                         _ref(graph, "initial:R_Q_y:Z2"))
    q4_source_chart = source_chart if (swapped_source or same_label_source) else loop_first
    first_micro = {
        "source_chart": q4_source_chart,
        "physical": physical,
        "finite_boundary": finite_boundary,
        "errors": {"chart": "retained"},
    }
    loop = _Owner()
    loop.g = graph
    loop.background = background
    loop.N = _ref(graph, "N")
    loop.charts = {"first_micro": loop_first}
    loop.errors = {"loop": "retained"}
    q4 = _Owner()
    q4.g = graph
    q4.loop = loop
    q4.background = background
    q4.charts = {"first_micro": first_micro}
    q4.initial = {"R_Q_y": initial_R_Q_y, "lower": lower}
    q4.initial_lower = lower
    q4.errors = {"q4": "retained"}

    if missing == "Q_y":
        del physical["Q_y"]
    elif missing == "R_Q_y":
        del finite_boundary["R_Q_y"]
    elif missing == "initial_Z1":
        initial_R_Q_y.rows = (initial_R_Q_y.rows[0],)
    elif missing == "source_domain":
        del source_chart["domain"]
    elif missing == "history":
        loop_first.pop("history", None)
    return q4, graph, other_graph, lower


def _check(checks, name, condition):
    checks[name] = bool(condition)


def _is_unresolved(value):
    return isinstance(value, dict) and value.get("status") == "unresolved"


def _case_diagnostic(case, passed, result=None):
    item = {"case": case, "passed": bool(passed)}
    if isinstance(result, dict):
        item["unresolved_count"] = len(result.get("unresolved", ()))
    return item


def run_fixture(extract):
    """Run source-only positive and negative identity checks."""

    checks = {}
    diagnostics = []
    if not callable(extract):
        return {"schema": SCHEMA, "passed": False,
                "checks": {"extract_callable": False},
                "diagnostics": [{"case": "input", "passed": False,
                                 "reason": "extract must be callable"}]}
    _check(checks, "extract_callable", True)

    q4, graph, other_graph, lower = _make_case()
    positive = extract(q4)
    _check(checks, "positive_bundle_mapping", isinstance(positive, dict))
    _check(checks, "positive_source_only", positive.get("source_only") is True)
    _check(checks, "positive_no_construction", positive.get("constructed") is False)
    _check(checks, "positive_no_numeric_evaluation",
           positive.get("numeric_values_computed") is False and
           positive.get("graph_handles_treated_as_values") is False)
    first = positive.get("first_micro", {})
    physical = first.get("physical_fields", {})
    _check(checks, "positive_physical_refs_preserved",
           physical.get("U") is q4.charts["first_micro"]["physical"]["U"] and
           physical.get("m") is q4.charts["first_micro"]["physical"]["m"] and
           physical.get("U_y") is q4.charts["first_micro"]["physical"]["U_y"] and
           physical.get("Q_y") is q4.charts["first_micro"]["physical"]["Q_y"])
    _check(checks, "positive_boundary_ref_preserved",
           first.get("finite_boundary_R_Q_y") is
           q4.charts["first_micro"]["finite_boundary"]["R_Q_y"])
    rows = first.get("initial_R_Q_y_rows", {})
    _check(checks, "positive_initial_Z0_preserved",
           rows.get("Z0") is q4.initial["R_Q_y"].rows[0])
    _check(checks, "positive_initial_Z1_preserved",
           rows.get("Z1") is q4.initial["R_Q_y"].rows[1])
    _check(checks, "positive_source_chart_ref_preserved",
           first.get("source_chart") is q4.loop.charts["first_micro"])
    # The production loop chart is the source chart for Q4; the synthetic
    # source fields are attached below so this check is explicit and identity-only.
    q4.loop.charts["first_micro"]["coordinate"] = q4.loop.charts["first_micro"].get("coordinate", _ref(graph, "coordinate"))
    q4.loop.charts["first_micro"]["radius"] = q4.loop.charts["first_micro"].get("radius", _ref(graph, "radius"))
    q4.loop.charts["first_micro"]["modified_fields"] = q4.loop.charts["first_micro"].get("modified_fields", {"V": _jet(graph, "V")})
    # Re-extract after attaching the exact source-chart fields, as a caller
    # would hand the same live chart object to Q4 and its loop.
    q4.charts["first_micro"]["source_chart"] = q4.loop.charts["first_micro"]
    positive = extract(q4)
    first = positive.get("first_micro", {})
    _check(checks, "positive_source_chart_fields_are_direct_refs",
           first.get("source_chart_fields", {}).get("coordinate") is
           q4.loop.charts["first_micro"]["coordinate"] and
           first.get("source_chart_fields", {}).get("radius") is
           q4.loop.charts["first_micro"]["radius"] and
           first.get("source_chart_fields", {}).get("modified_fields") is
           q4.loop.charts["first_micro"]["modified_fields"])
    _check(checks, "attached_lower_does_not_prove_substitution",
           positive.get("checks", {}).get("first_micro_initial_lower_is_source_chart_domain_lower") is False and
           positive.get("first_micro", {}).get("initial_lower_source") is
           q4.loop.charts["first_micro"]["domain"][0])
    _check(checks, "positive_error_state_refs_preserved",
           positive.get("errors", {}).get("q4") is q4.errors and
           positive.get("errors", {}).get("loop") is q4.loop.errors)
    diagnostics.append(_case_diagnostic("positive", True, positive))

    missing = extract(_make_case(missing="Q_y")[0])
    _check(checks, "missing_Q_y_is_unresolved",
           _is_unresolved(missing.get("first_micro", {}).get("physical_fields", {}).get("Q_y")))
    diagnostics.append(_case_diagnostic("missing_Q_y", checks["missing_Q_y_is_unresolved"], missing))

    swapped_graph = extract(_make_case(swapped_graph=True)[0])
    _check(checks, "changed_graph_fails_identity",
           swapped_graph.get("checks", {}).get("physical_U_graph_identity") is False)
    diagnostics.append(_case_diagnostic("changed_graph", checks["changed_graph_fails_identity"], swapped_graph))

    swapped_source = extract(_make_case(swapped_source=True)[0])
    _check(checks, "changed_source_chart_fails_identity",
           swapped_source.get("checks", {}).get("q4_first_micro_source_chart_is_loop_first_micro") is False)
    diagnostics.append(_case_diagnostic("changed_source_chart", checks["changed_source_chart_fails_identity"], swapped_source))

    no_history = extract(_make_case(remove_history=True)[0])
    _check(checks, "deleted_history_is_unresolved",
           _is_unresolved(no_history.get("loop_first_micro", {}).get("history")) and
           no_history.get("checks", {}).get("loop_first_micro_incoming_history_present") is False)
    diagnostics.append(_case_diagnostic("deleted_history", checks["deleted_history_is_unresolved"], no_history))

    same_labels = extract(_make_case(same_label_source=True)[0])
    _check(checks, "same_labels_do_not_prove_source_identity",
           same_labels.get("checks", {}).get("q4_first_micro_source_chart_is_loop_first_micro") is False)
    diagnostics.append(_case_diagnostic("same_numerical_labels", checks["same_labels_do_not_prove_source_identity"], same_labels))

    corrected_graph = _Graph()
    corrected_packet = {"U": _jet(corrected_graph, "U"), "m": _jet(corrected_graph, "m")}
    corrected_view = {"Q": _jet(corrected_graph, "Q")}
    optional = extract(q4, corrected_packet=corrected_packet, corrected_view=corrected_view)
    optional_checks = optional.get("optional_internal_views", {})
    _check(checks, "corrected_packet_internal_only",
           optional_checks.get("corrected_packet", {}).get("graph_identity") is True and
           optional_checks.get("corrected_packet", {}).get("cross_source_identity_checked") is False)
    _check(checks, "corrected_view_internal_only",
           optional_checks.get("corrected_view", {}).get("graph_identity") is True and
           optional_checks.get("corrected_view", {}).get("cross_source_identity_checked") is False)
    _check(checks, "no_cross_owner_claim",
           optional.get("cross_endpoint_identity_claimed") is False and
           optional.get("native_corrected_identity_proved") is False)
    diagnostics.append(_case_diagnostic("optional_corrected_views", True, optional))

    passed = all(value is True for value in checks.values())
    return {"schema": SCHEMA, "passed": passed, "checks": checks, "diagnostics": diagnostics}


run = run_fixture

__all__ = ["SCHEMA", "run_fixture", "run"]
