"""Bind a corrected physical preRp source to the independent pulse target.

Only --compute-new executes the new source-DAG constructor. Reports and
inspection consume immutable checkpoints; no ancestor mathematics is run.
"""
from __future__ import annotations

import argparse
import ast
import copy
from fractions import Fraction
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "experiments/root_st073"
PREFIX = "lei_ren_part1_paper_compliant_current_original_"
PARENT = PREFIX + "outer_pressure_enclosure.json"
PARENT_HASH = "e6c2725605ad8a6f76c1f07d204f13370e2fbf68456ed25044b72246e94f6dcb"
NATIVE = PREFIX + "native_directed_controls_v2.json"
NATIVE_HASH = "10abe68b98fdf87709a7a7c9f27103523701e746c357ae62386038700052490a"
PRE = PREFIX + "n1_preRp_moment_functions.json.gz"
FULL = PREFIX + "n1_full_defect_linear_functions.json.gz"
PULSE_BRIDGE = PREFIX + "n1_pulse_numeric_bridge.json.gz"
FULL_BRIDGE = PREFIX + "n1_full_five_numeric_bridge.json.gz"
ENDPOINT = PREFIX + "n1_preRc_endpoint_C1.json.gz"
KERNEL = "src/openai_ns_reconstruction/original_corrected_reserved_m15_source.py"
KERNEL_HASH = "0dad10c6721eb3456b186bc69eb7f25c1a7a5260b97ca119c606ff3af4607acd"
FIXTURE = "experiments/root_st073/original_corrected_reserved_m15_source_fixture.json"
FIXTURE_HASH = "3d1659ce68034752e8b124dc6ad4777c27be5eb7df4d8e346622b15a7fd15dcd"
NUMERIC_KERNEL = "src/openai_ns_reconstruction/original_corrected_reserved_m15_numeric.py"
NUMERIC_KERNEL_HASH = "5c1d17b998f475c505d50e3bac3a6fdaa7520085ead1a29701a1fcc55949a7bb"
NUMERIC_FIXTURE = "experiments/root_st073/original_corrected_reserved_m15_numeric_fixture.json"
NUMERIC_FIXTURE_HASH = "cd6b4d89a96229340fc014420f4b82ef00ace8bcf616e999fa267b3b2db036af"
RUNNER = "scripts/run_original_corrected_reserved_m15.py"
REPORT = BASE / (PREFIX + "corrected_reserved_m15_source.json")
CHECKPOINTS = ROOT / "tmp/original_corrected_reserved_m15"
SOURCE = CHECKPOINTS / "new_source.json"
CONSTRUCTED = CHECKPOINTS / "new_corrected_source.json"
NORMALIZED = CHECKPOINTS / "new_normalized_source.json"
NUMERIC = CHECKPOINTS / "new_corrected_numeric.json"


def require(value, message):
    if not value:
        raise ValueError(message)


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                   ensure_ascii=True, allow_nan=False).encode()).hexdigest()


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def path_for(name):
    p = Path(name)
    require(not p.is_absolute(), "absolute input path")
    if ".." in p.parts:
        require(name == "../../src/openai_ns_reconstruction/schedule_pressure.py", "unexpected traversal")
        return (BASE / p).resolve()
    return ROOT / p if p.parts[0] in ("experiments", "src", "scripts", "tmp") else BASE / p.name


def read_pinned(path, expected):
    raw = Path(path).read_bytes()
    require(hashlib.sha256(raw).hexdigest() == expected, "changed pinned input: " + str(path))
    return json.loads(gzip.decompress(raw) if str(path).endswith(".gz") else raw)


def source_preflight():
    for relative in (RUNNER, KERNEL, NUMERIC_KERNEL):
        raw = (ROOT / relative).read_bytes()
        require(b"\r" not in raw and raw.endswith(b"\n") and not raw.endswith(b"\n\n"), "LF/EOF: " + relative)
        tree = ast.parse(raw, filename=relative)
        names = [n.name for n in tree.body if isinstance(n, ast.FunctionDef)]
        require(len(names) == len(set(names)), "duplicate functions: " + relative)
    require(sha(ROOT / KERNEL) == KERNEL_HASH, "unreviewed new kernel")
    require(sha(ROOT / NUMERIC_KERNEL) == NUMERIC_KERNEL_HASH, "unreviewed numeric kernel")
    fixture = read_pinned(ROOT / FIXTURE, FIXTURE_HASH)
    require(fixture.get("passed") is True, "new synthetic fixture failed")
    require(read_pinned(ROOT / NUMERIC_FIXTURE, NUMERIC_FIXTURE_HASH)["passed"] is True,
            "new numeric synthetic fixture failed")


def setup():
    source_preflight()
    parent = read_pinned(BASE / PARENT, PARENT_HASH)
    native = read_pinned(BASE / NATIVE, NATIVE_HASH)
    for name, expected in {**parent["input_hashes"], **parent["implementation_hashes"]}.items():
        require(sha(path_for(name)) == expected, "changed inherited pin: " + name)
    pre = read_pinned(BASE / PRE, parent["input_hashes"][PRE])
    full = read_pinned(BASE / FULL, parent["input_hashes"][FULL])
    pulse_bridge = read_pinned(BASE / PULSE_BRIDGE, parent["input_hashes"][PULSE_BRIDGE])
    full_bridge = read_pinned(BASE / FULL_BRIDGE, parent["input_hashes"][FULL_BRIDGE])
    endpoint = read_pinned(BASE / ENDPOINT, parent["input_hashes"][ENDPOINT])
    require(pre["source_family"] == full["source_family"] == native["source_family"], "different source family")
    require(pulse_bridge["source_family"] == full_bridge["source_family"] == endpoint["source_family"] == native["source_family"],
            "normalization or endpoint source family changed")
    full_nodes = copy.deepcopy(full["exact_graph_nodes"])
    span, rows = pulse_bridge["normalization_graph_span"], pulse_bridge["normalization_graph_nodes"]
    require(span == [35558, 53438] and len(rows) == span[1] - span[0], "pulse normalization span changed")
    require(len(full_nodes) == 53398 and full_nodes[span[0]:] == rows[:53398-span[0]],
            "pulse normalization overlap is not exactly the existing full graph")
    full_nodes.extend(copy.deepcopy(rows[53398-span[0]:]))
    span2, rows2 = full_bridge["normalization_graph_span"], full_bridge["normalization_graph_nodes"]
    require(span2 == [53438, 53444] and len(rows2) == 6 and len(full_nodes) == span2[0],
            "full normalization span changed")
    full_nodes.extend(copy.deepcopy(rows2))
    require(full_nodes[53427] == dict(operation="analytic_unary", name="exp", argument=53421)
            and full_nodes[53430] == dict(operation="product", arguments=[53427,53427]),
            "common positive unit recipe changed")
    patch = native["actual_first_micro_seam"]["corrected_physical_endpoint_graph"]
    require(digest(patch) == "19322a78adef610dd8298a81c7157bce2ecd2743c71ef0311340a8f22f65b235", "changed corrected patch")
    seam = native["actual_first_micro_seam"]
    require(seam["actual_graph_scale_defect_checks"]["passed"] is True, "saved scale identity not admitted")
    binding = full["full_defect_linear_bindings"]
    pairs = {name: copy.deepcopy(binding["segments"][name]["m15"])
             for name in ("inner", "lower", "preRp", "pulse")}
    pairs["full"] = copy.deepcopy(binding["moments"]["m15"])
    old = native["actual_reserved_point_m15"]
    require(old["source_C1_pairs"] == pairs, "old point source ports changed")
    shared = old["shared_unit"]
    unit = dict(positive=old["source_admission"]["positive_reference_and_energy_epsilon_applied_once"],
                Z_independent=shared["U_Z_exact_zero"],
                ordinary_Z_derivative=shared["ordinary_C1_is_scaled_physical_derivative"],
                definition=shared["physical"], graph_namespace="complete_saved_physical_and_normalization_prefix",
                external_graph_ports=copy.deepcopy(shared),
                common_m15_unit_log=copy.deepcopy(old["common_m15_unit_log"]),
                physical_graph_operand=True, unit_node=53427, unit_square_node=53430,
                complete_prefix_node_count=len(full_nodes))
    end = endpoint["samples"]["('0', '.5')"]
    require(end["exact_Z_cell"] == ["0", ".5"] and end["proof"]["physical_common_Pstar_applied_once"] is True,
            "physical endpoint cell or scale identity changed")
    lower = end["lower"]
    require(lower["proof"]["actual_V_y_rows_0_through_3_factorial_converted_once"] is True
            and lower["proof"]["common_Pstar_applied_once"] is True,
            "ordinary physical lower Y source changed")
    for row in lower["source_rows"]["V_y"][:2]:
        require(row["encloses_original_source_function"] is True and row["point_value_selected"] is False,
                "Y rows must be source enclosures, not selected centers")
    def coordinate(bits):
        sign, mantissa, exponent, _ = bits
        return Fraction((-1 if sign else 1)*mantissa) * Fraction(2)**exponent
    point = old["selected_Z_interval"]
    zlo, zhi = (coordinate(point[k]) for k in ("lower_exact_mpf_tuple", "upper_exact_mpf_tuple"))
    require(0 < zlo <= Fraction(371,1000) <= zhi < Fraction(1,2), "numeric point outside saved endpoint cell")
    source = dict(schema="corrected_reserved_m15_inputs_v1", source_only=True,
                  parent_sha256=PARENT_HASH, native_sha256=NATIVE_HASH,
                  pre_graph_sha256=parent["input_hashes"][PRE], full_graph_sha256=parent["input_hashes"][FULL],
                  pre_node_count=len(pre["exact_graph_nodes"]), full_node_count=len(full_nodes),
                  pre_nodes_digest=digest(pre["exact_graph_nodes"]), full_nodes_digest=digest(full_nodes),
                  normalization_span_receipt=dict(full_original_count=53398, pulse_span=[35558,53438],
                      pulse_overlap_nodes=17840, pulse_appended_nodes=40, full_span=[53438,53444],
                      full_appended_nodes=6, exact_overlap_checked=True),
                  saved_patch=copy.deepcopy(patch), saved_patch_digest=digest(patch),
                  source_pairs=pairs, unit_receipt=unit,
                  physical_endpoint_source_rows=copy.deepcopy(lower["source_rows"]),
                  physical_endpoint_source_proof=copy.deepcopy(lower["proof"]),
                  numeric_point_Z=copy.deepcopy(old["selected_Z_interval"]),
                  numeric_point_exact_rational="371/1000", numeric_endpoint_cell=["0", ".5"],
                  delta=copy.deepcopy(native["actual_micro_Q_increment_bounds"]["actual_delta"]),
                  actual_logPstar_parameter=copy.deepcopy(native["actual_micro_background_I1_numeric"]["physical_scale"]["logPstar"]),
                  actual_logPstar_source="native.actual_micro_background_I1_numeric.physical_scale.logPstar; common original parameter, not endpoint magnitude cap",
                  old_discrepancy=[copy.deepcopy(old["ordinary_Z"][str(k)]["correlated_reserved_minus_target"]) for k in (0,1)],
                  independent_target=[copy.deepcopy(old["ordinary_Z"][str(k)]["independent_negative_complete_pulse_target"]) for k in (0,1)],
                  source_partition=copy.deepcopy(binding["partition"]),
                  old_graphs_unchanged=True, numerical_ancestor_executed=False)
    return parent, native, pre, full_nodes, source


def write_new(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("xb") as stream:
        stream.write((json.dumps(value, indent=2, ensure_ascii=True, allow_nan=False) + "\n").encode())


def read_checkpoint(path):
    return json.loads(path.read_bytes())


def compute_new():
    parent, native, pre, full_nodes, source = setup()
    require(not any(p.exists() for p in (SOURCE, CONSTRUCTED, NORMALIZED, NUMERIC)), "checkpoints exist; do not replay")
    write_new(SOURCE, source)
    spec = importlib.util.spec_from_file_location("new_corrected_reserved_source", ROOT / KERNEL)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    started = time.perf_counter()
    result = module.construct(pre["exact_graph_nodes"], full_nodes,
                              source["saved_patch"], source["source_pairs"],
                              expected_base_digest=source["pre_nodes_digest"],
                              expected_full_digest=source["full_nodes_digest"],
                              unit_receipt=source["unit_receipt"])
    result["execution_seconds"] = time.perf_counter() - started
    result["source_checkpoint_sha256"] = sha(SOURCE)
    write_new(CONSTRUCTED, result)
    result = read_checkpoint(CONSTRUCTED)
    normalized = normalize_new_source(result, source, full_nodes)
    write_new(NORMALIZED, normalized)
    numeric = compute_new_numeric(source)
    write_new(NUMERIC, numeric)
    return dict(passed=True, source_checkpoint_sha256=sha(SOURCE), constructed_sha256=sha(CONSTRUCTED),
                normalized_sha256=sha(NORMALIZED), numeric_sha256=sha(NUMERIC))


def normalize_new_source(constructed, source, full_nodes):
    """Append only quotients by the authenticated positive, frozen unit."""
    require(constructed["passed"] is True and source["unit_receipt"]["Z_independent"] is True,
            "source construction or frozen unit not admitted")
    nodes = full_nodes + constructed["delta_nodes"]
    require(digest(nodes) == constructed["graph_patch"]["full_canonical_nodes_sha256"], "constructed graph changed")
    delta, roots = [], {}
    started = time.perf_counter()
    for name, pair in constructed["physical_C1_new_roots"].items():
        roots[name] = []
        for numerator in pair:
            require(0 <= numerator < len(nodes), "physical source numerator outside prefix")
            roots[name].append(len(nodes) + len(delta))
            delta.append(dict(operation="positive_quotient", numerator=numerator, denominator=53427,
                              positivity_certificate="accepted positive reference*A*Rv^2; Z-independent"))
    return dict(passed=True, source_only=True, base_node_count=len(nodes), appended_nodes=delta,
                normalized_C1_roots=roots, common_unit_node=53427,
                ordinary_Z1_divided_by_same_unit=True, no_second_F0_derivative=True,
                unit_Z1_exact_zero=True, source_construction_sha256=sha(CONSTRUCTED),
                full_node_count=len(nodes)+len(delta), full_nodes_digest=digest(nodes+delta),
                execution_seconds=time.perf_counter()-started, numerical_value_installed=False,
                accepted_ancestors_replayed=False)


def compute_new_numeric(source):
    """Consume physical saved log bounds and the old fixed-unit point rows."""
    spec = importlib.util.spec_from_file_location("new_corrected_reserved_numeric", ROOT / NUMERIC_KERNEL)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    from mpmath.ctx_iv import MPIntervalContext
    c = MPIntervalContext()
    c.prec = 1800
    rows = source["physical_endpoint_source_rows"]
    started = time.perf_counter()
    result = module.enclose(c, source["numeric_point_Z"], source["delta"],
                            rows["R"]["log_absolute_upper"], source["actual_logPstar_parameter"],
                            [row["log_absolute_upper"] for row in rows["V_y"][:2]],
                            source["unit_receipt"]["common_m15_unit_log"],
                            source["old_discrepancy"], source["independent_target"],
                            radius_bits=4096, scale_Z_independent=True, same_source_repair=True)
    result = serialize(result)
    result.update(execution_seconds=time.perf_counter()-started,
                  normalized_source_checkpoint_sha256=sha(NORMALIZED),
                  source_checkpoint_sha256=sha(SOURCE), point_only=True,
                  exact_rational_Z="371/1000", context_precision_bits=1800)
    return result


def serialize(value):
    if hasattr(value, "_mpi_"):
        lo, hi = value._mpi_
        return dict(lower_exact_mpf_tuple=list(lo), upper_exact_mpf_tuple=list(hi), directed=True)
    if isinstance(value, dict):
        return {str(k): serialize(v) for k,v in value.items()}
    if isinstance(value, (list, tuple)):
        return [serialize(v) for v in value]
    return value


def write_report():
    parent, native, pre, full, expected_source = setup()
    source, result = read_checkpoint(SOURCE), read_checkpoint(CONSTRUCTED)
    require(source == expected_source and result["passed"] is True, "source checkpoint mismatch")
    normalized, numeric = read_checkpoint(NORMALIZED), read_checkpoint(NUMERIC)
    require(normalized["passed"] is True and numeric["passed"] is True, "new checkpoint failed")
    report = copy.deepcopy(parent)
    report.update(stage="corrected_physical_reserved_m15_source_function", parent_report=PARENT,
                  parent_report_sha256=PARENT_HASH, corrected_terminal_source_inputs=source,
                  corrected_reserved_m15_source=result,
                  corrected_reserved_normalized_source=normalized,
                  corrected_reserved_point_numeric=numeric,
                  execution_seconds=dict(new_source_construction=result["execution_seconds"],
                      new_normalization_source=normalized["execution_seconds"],
                      new_point_repair=numeric["execution_seconds"]),
                  corrected_reserved_numeric_enclosure=True, corrected_discrepancy_numeric_enclosure=True,
                  corrected_numeric_scope="point Z=371/1000 only; not axial functional closure",
                  legacy_reserved_numeric_boxes_not_transferred=True,
                  independent_negative_complete_pulse_target_unchanged=True,
                  old_reserved_point_packet_preserved=True,
                  source_function_view_only=False, existing_velocity_provider_rebound=False,
                  old_212_terms_and_all_errors_unchanged=True,
                  reserved_n1_own_swirl_registration_count=1,
                  extra_reserved_term_added=False,
                  corrected_source_checkpoint_sha256=sha(CONSTRUCTED),
                  normalized_source_checkpoint_sha256=sha(NORMALIZED), numeric_checkpoint_sha256=sha(NUMERIC))
    report["input_hashes"] = dict(parent["input_hashes"], **{
        PARENT: PARENT_HASH, NATIVE: NATIVE_HASH, FIXTURE: FIXTURE_HASH, NUMERIC_FIXTURE: NUMERIC_FIXTURE_HASH})
    report["implementation_hashes"] = dict(parent["implementation_hashes"], **{
        KERNEL: sha(ROOT / KERNEL), NUMERIC_KERNEL: sha(ROOT / NUMERIC_KERNEL), RUNNER: sha(ROOT / RUNNER)})
    report["scope"] = dict(parent["scope"], corrected_reserved_source_function=True,
                           independent_terminal_function_match=False, global_five_moment_closure=False,
                           full_goal_complete=False, numerical_corrected_point_target=True,
                           all_Z_functional_closure=False, existing_velocity_provider_rebound=False)
    write_new(REPORT, report)
    return dict(passed=True, bytes=REPORT.stat().st_size, sha256=sha(REPORT))


def inspect_report():
    parent, native, pre, full, expected_source = setup()
    report = read_checkpoint(REPORT)
    source, result = read_checkpoint(SOURCE), read_checkpoint(CONSTRUCTED)
    require(source == expected_source and report["corrected_terminal_source_inputs"] == source, "source changed")
    require(report["corrected_reserved_m15_source"] == result, "constructed checkpoint changed")
    require(report["corrected_source_checkpoint_sha256"] == sha(CONSTRUCTED), "constructed hash changed")
    normalized, numeric = read_checkpoint(NORMALIZED), read_checkpoint(NUMERIC)
    require(report["corrected_reserved_normalized_source"] == normalized
            and report["corrected_reserved_point_numeric"] == numeric, "new checkpoint/report mismatch")
    require(report["normalized_source_checkpoint_sha256"] == sha(NORMALIZED)
            and report["numeric_checkpoint_sha256"] == sha(NUMERIC), "new checkpoint hashes changed")
    require(numeric["source_checkpoint_sha256"] == sha(SOURCE)
            and numeric["normalized_source_checkpoint_sha256"] == sha(NORMALIZED), "numeric source links changed")
    require(report["scope"]["full_goal_complete"] is False, "invalid completion claim")
    require(report["implementation_hashes"][RUNNER] == sha(ROOT / RUNNER), "runner changed")
    for name, expected in {**report["input_hashes"], **report["implementation_hashes"]}.items():
        require(sha(path_for(name)) == expected, "report pin changed: " + name)
    require(normalized["base_node_count"] == result["graph_patch"]["full_node_count"], "normalization base changed")
    expected_pairs = result["physical_C1_new_roots"]
    cursor = 0
    for name, pair in expected_pairs.items():
        for order, numerator in enumerate(pair):
            row = normalized["appended_nodes"][cursor]
            require(row["operation"] == "positive_quotient" and row["numerator"] == numerator
                    and row["denominator"] == 53427, "new normalized source operand changed")
            require(normalized["normalized_C1_roots"][name][order] == normalized["base_node_count"]+cursor,
                    "new normalized root changed")
            cursor += 1
    require(cursor == len(normalized["appended_nodes"]), "extra normalization row")
    for key in ("actual_outer_pressure_enclosure", "actual_outer_signed_macro_scale", "actual_H_m15_Rm_transport",
                "actual_Rm_absolute_m15", "new_partial_rows", "new_absolute_rows"):
        require(report[key] == parent[key], "inherited packet changed: " + key)
    return dict(passed=True, serialized_only=True, bytes=REPORT.stat().st_size,
                sha256=sha(REPORT), source_checkpoint_sha256=sha(SOURCE),
                constructed_checkpoint_sha256=sha(CONSTRUCTED), old_mathematics_executed=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--compute-new", action="store_true")
    group.add_argument("--write-report", action="store_true")
    group.add_argument("--inspect-report", action="store_true")
    args = parser.parse_args()
    action = compute_new if args.compute_new else write_report if args.write_report else inspect_report
    print(json.dumps(action(), indent=2))
