"""Source-bound lower O3 endpoint budget for the Rc F port.

This module reads the saved original post-axial packets.  It deliberately
does not construct an owner, replay a producer, or turn a finite-N cap into a
signed source value.  The only arithmetic performed here is directed interval
addition of the saved lower endpoint log enclosure and the seven saved
physical log-radius atoms at O3 native zero.
"""
import ast
import hashlib
import json
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext


ROOT = Path(__file__).resolve().parents[2]
POST_AXIAL = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_post_axial_m15_numeric.json"
PARENT = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_correlated_scaled_source.json"
COMPACT_RECEIPT = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_true_compact_C5_functions_check.json"
COMPACT_ARCHIVE = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_true_compact_C5_functions.json.gz"

EXPECTED_SHA256 = {
    POST_AXIAL: "48272c642ed8a28248aab35e7feda450f6a10aeff534f77b29f4b56a07f913f4",
    PARENT: "9beead66814da85325e86f4ee2d9c091da34c524f9512f8cbf4c77e18bf1d444",
    COMPACT_RECEIPT: "6c933d077f9c151eda246acd4d5f694d32dc8272fc3fb65a031b9b9cf28c8e7b",
    COMPACT_ARCHIVE: "3593dd1ba3e510081090dd315d8581c1253afd5192233a284079c0646a556269",
    "src/openai_ns_reconstruction/original_post_axial_m15.py": "1c7ff372dd75c4b69072ba06908ba1b3f881d82ba52fc21e46da621ccffb6a69",
    "src/openai_ns_reconstruction/original_axial_m15.py": "49f8d08c3f6418c172ba1339896cf631a6a54d10b886fc97e8aff234fccd0e50",
    "src/openai_ns_reconstruction/original_post_axial_pressure_rows.py": "ab37ca5a58935008ff2a8cc4717080ba40b3ec5268d38c924def02fc5517e812",
    "src/openai_ns_reconstruction/original_macro_m15.py": "5d2c162e6c61bee3d05955004114c9be73b44471555c0638a06b0bfb3d1febf5",
}

_METHOD_SHA256 = {
    ("src/openai_ns_reconstruction/original_axial_m15.py", "OriginalAxialM15", "finite_errors"):
        "5bafa555dab97b8ba234eeaa4f7cc10894445b3a5dcecac87aa59768c21fa499",
    ("src/openai_ns_reconstruction/original_macro_m15.py", "OriginalMacroM15", "correction_densities"):
        "3efeb14fd215e36a1b8127b1bf0eb6d55a017817594d6de47e6187bc3f812514",
}


def _sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _pin(root, relative):
    path = root / relative
    actual = _sha(path)
    expected = EXPECTED_SHA256[relative]
    if actual != expected:
        raise ValueError("Pinned source changed: " + relative)
    return dict(path=relative, sha256=actual, expected_sha256=expected, passed=True)


def _method_hash(root, relative, class_name, method_name):
    path = root / relative
    tree = ast.parse(path.read_text(encoding="utf-8"))
    for owner in ast.walk(tree):
        if isinstance(owner, ast.ClassDef) and owner.name == class_name:
            for node in owner.body:
                if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == method_name:
                    digest = hashlib.sha256(ast.dump(node, include_attributes=False).encode()).hexdigest()
                    expected = _METHOD_SHA256[(relative, class_name, method_name)]
                    if digest != expected:
                        raise ValueError("Pinned method changed: " + relative + ":" + method_name)
                    return dict(class_name=class_name, method=method_name, ast_sha256=digest,
                                expected_ast_sha256=expected, passed=True)
    raise ValueError("Pinned method missing: " + relative + ":" + method_name)


def _exact_interval(iv, pair):
    if not isinstance(pair, list) or len(pair) != 2:
        raise ValueError("Expected a directed two-endpoint mpf tuple pair")
    return iv.mpf([mp.make_mpf(tuple(pair[0])), mp.make_mpf(tuple(pair[1]))])


def _serialized(iv_value):
    return [list(endpoint) for endpoint in iv_value._mpi_]


def _read_json(root, relative):
    return json.loads((root / relative).read_text(encoding="utf-8"))


def _sample(samples, Z, coordinate, chart="O3_power"):
    matches = [row for row in samples if row.get("Z") == Z and row.get("chart") == chart
               and str(row.get("native_coordinate")) == coordinate]
    if len(matches) != 1:
        raise ValueError("Expected one saved sample for " + Z + " / " + coordinate)
    return matches[0]


def _compact_upper_certificate(root, parent):
    receipt = _read_json(root, COMPACT_RECEIPT)
    receipt_pin = _pin(root, COMPACT_RECEIPT)
    archive_pin = _pin(root, COMPACT_ARCHIVE)
    checks = receipt.get("checks", {})
    source_check = checks.get("actual_source", {})
    endpoint_check = checks.get("actual_inlet_and_endpoint", {})
    required = (
        receipt.get("all_passed") is True,
        receipt.get("current_original_true_Rc_to_2Rc_relative_five_C5_zero_identities_installed") is True,
        source_check.get("original_compact_log_transport_measure_and_relative_only_terminal_theorem") is True,
        endpoint_check.get("terminal_integral_identity_uses_pinned_original_exact_bump_theorem") is True,
    )
    if not all(required):
        raise ValueError("Saved compact upper-2 source certificate is incomplete")
    return dict(
        path=COMPACT_RECEIPT, sha256=receipt_pin["sha256"], passed=True, source_only=True,
        source_values_selected=False, numeric_caps_used=False,
        compact_archive=archive_pin, chart="compact", native_coordinate="2",
        endpoint="compact upper 2 Rc", relative_mean_C5_exact_zero=True,
        delta_F_zero_requires_actual_velocity_and_pressure_recovery_proof=True,
        source_identity="relative five C5 terminal correction rows are exact zero",
        terminal_nodes=[263339, 263340, 263341, 263342, 263343],
        parent_full_Rc_F_identity_consumes_certificate=True,
        upper_proof_repeated=False,
    )


def rebuild_Rc_endpoint_budget(repo_root=ROOT, *, precision=200):
    """Return exact directed C0/C1 log bounds for ``-R0² delta_F_lower``.

    The saved O3 native-2 lower endpoint row supplies the positive endpoint
    enclosure.  The saved O3 native-0 geometry supplies the seven physical
    log-radius atoms.  No cap, midpoint, aggregate subtraction, Qy term, or
    viscous term is selected as a value.
    """
    root = Path(repo_root)
    pins = {relative: _pin(root, relative) for relative in (
        POST_AXIAL, PARENT,
        "src/openai_ns_reconstruction/original_post_axial_m15.py",
        "src/openai_ns_reconstruction/original_axial_m15.py",
        "src/openai_ns_reconstruction/original_post_axial_pressure_rows.py",
        "src/openai_ns_reconstruction/original_macro_m15.py",
    )}
    methods = dict(
        finite_errors=_method_hash(root, "src/openai_ns_reconstruction/original_axial_m15.py", "OriginalAxialM15", "finite_errors"),
        correction_densities=_method_hash(root, "src/openai_ns_reconstruction/original_macro_m15.py", "OriginalMacroM15", "correction_densities"),
    )
    parent = _read_json(root, PARENT)
    identity = parent.get("compact_inlet_source_identity", {})
    if not (parent.get("contract", {}).get("incoming_pressure_and_all_finite_N_caps_retained") is True
            and identity.get("passed") is True
            and identity.get("endpoint_coordinate") == "compact x=1 equals outer O3_power upper endpoint at native x=2"):
        raise ValueError("Current Rc parent source identity changed")
    report = _read_json(root, POST_AXIAL)
    if report.get("checks", {}).get("passed") is not True:
        raise ValueError("Saved post-axial report is not passed")
    if report.get("source_family") != parent.get("source_family"):
        raise ValueError("Saved post-axial source family differs from current Rc parent")
    iv = MPIntervalContext(); iv.dps = int(precision)
    output_samples = []
    for Z in ("-.3", "0", ".3"):
        upper = _sample(report["samples"], Z, "2")
        origin = _sample(report["samples"], Z, "0")
        budget = upper.get("actual_lower_finite_N_endpoint_budget")
        if not isinstance(budget, dict) or budget.get("normalization") != "Rupper^2; physical dd once":
            raise ValueError("O3 lower endpoint budget normalization changed")
        endpoint_rows = budget.get("endpoint")
        atoms = origin.get("geometry", {}).get("atomic_log_radius_terms")
        if not isinstance(endpoint_rows, list) or len(endpoint_rows) != 2 or not isinstance(atoms, list) or len(atoms) != 7:
            raise ValueError("Expected two endpoint rows and seven native-zero radius atoms")
        radius_terms = [_exact_interval(iv, atom["source_enclosure_exact_mpf_tuples"]) for atom in atoms]
        radius_log = sum(radius_terms, iv.mpf(0)); radius_twice = 2 * radius_log
        radius_atom_records = [dict(index=index, log_mean_point=atom["approximation"],
            absolute_error_upper=atom["absolute_error_upper"],
            raw_directed_error_exact_mpf_tuple=atom["raw_directed_error_exact_mpf_tuple"],
            source_enclosure_exact_mpf_tuples=atom["source_enclosure_exact_mpf_tuples"],
            error_applies_to_serialized_decimal=atom["error_applies_to_serialized_decimal"],
            source_value_selected=False) for index, atom in enumerate(atoms)]
        components = []
        for index, row in enumerate(endpoint_rows):
            if row.get("exact_zero") is True:
                raise ValueError("Lower endpoint row unexpectedly exact zero")
            log_pair = row["log_absolute_upper"]
            log_bound = _exact_interval(iv, [log_pair["lower_exact_mpf_tuple"], log_pair["upper_exact_mpf_tuple"]])
            components.append(dict(component="C" + str(index), ordinary_Z_order=index,
                endpoint_log_absolute_upper_exact_mpf_tuples=_serialized(log_bound),
                weighted_log_bound_exact_mpf_tuples=_serialized(log_bound + radius_twice)))
        output_samples.append(dict(Z=Z, chart="O3_power", upper_native_coordinate="2",
            lower_native_coordinate="0", radius_atom_count=7,
            radius_log_exact_mpf_tuples=_serialized(radius_log),
            twice_radius_log_exact_mpf_tuples=_serialized(radius_twice), components=components,
            geometry_provenance=dict(
                physical_local_log_distance=origin["geometry"]["physical_local_log_distance"],
                native_coordinate_is_physical_log_offset=origin["geometry"]["native_coordinate_is_physical_log_offset"],
                physical_radius_materialized=origin["geometry"]["physical_radius_materialized"],
                exact_start=origin["geometry"]["exact_start"]),
            radius_atom_records=radius_atom_records,
            source_endpoint_budget=dict(endpoint_rows=[dict(exact_zero=row["exact_zero"],
                log_absolute_upper=row["log_absolute_upper"], source_value_selected=False)
                for row in endpoint_rows],
                normalization=budget["normalization"]),
            radius_atoms=[atom["source_enclosure_exact_mpf_tuples"] for atom in atoms]))
    return dict(schema="original_Rc_endpoint_budget_v1", passed=True,
        quantity="-R0^2*delta_F_lower", absolute_log_quantity="log|R0^2*delta_F_lower|",
        signed_orientation="negative", ordinary_Z_orders=[0, 1], samples=output_samples,
        source_family=parent.get("source_family"),
        source_pins=pins, source_methods=methods,
        compact_upper2_relative_C5_source_certificate=_compact_upper_certificate(root, parent),
        source_only=True, source_executed=False, source_values_selected=False,
        numeric_caps_used=False, contraction_applied=False, ancestor_replayed=False,
        no_aggregate_positive_bound_subtraction=True, no_midpoint_or_cap_as_value=True,
        qy_added=False, viscous_external_term_untouched=True,
        finite_errors_formula="lower endpoint log enclosure + 2 * native-zero physical log radius",
        remaining_term_provenance=dict(
            sign="-", radial_power=2, chart="O3_power", endpoint_native_coordinate="2",
            radius_anchor_native_coordinate="0", ordinary_components=["C0", "C1"],
            source_field="actual_lower_finite_N_endpoint_budget.endpoint",
            source_values_selected=False),
        F_coefficient_provenance="OriginalMacroM15.correction_densities -> endpoint rows.combine",
        parent_report=PARENT, post_axial_report=POST_AXIAL)


build = rebuild_Rc_endpoint_budget
