"""Produce or inspect second-switch/post-power finite-N volume and pressure sources.

Only this stage's five new builders execute during production. Inspection reads
serialized data and byte pins. Accepted source functions are never replayed.
"""
from __future__ import annotations
import argparse
import copy
import gzip
import hashlib
import json
from pathlib import Path
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/"experiments/root_st073"
PARENT_NAME="lei_ren_part1_paper_compliant_current_original_first_switch_finite_N_I1.json"
PARENT_SHA256="2cc8b4ff9670858564297bc1c27d5322c9af7d2cb85381939dc3e37a3eea346a"
NAME="lei_ren_part1_paper_compliant_current_original_switch_post_finite_N_I1.json"
BACKGROUND_KEY="actual_switch_post_background_I1_numeric"
INPUT_KEY="actual_switch_post_loop_Z3_inputs"
LOOP_KEY="actual_switch_post_loop_Z3_bounds"
NUMERIC_KEY="actual_switch_post_finite_N_I1_numeric"
PRESSURE_KEY="actual_switch_post_pressure_source"
CHARTS=("second_switch","post_switch")
CODE=(
 "src/openai_ns_reconstruction/original_switch_post_background_I1_numeric.py",
 "src/openai_ns_reconstruction/original_switch_post_loop_Z3_inputs.py",
 "src/openai_ns_reconstruction/original_switch_post_loop_Z3_bounds.py",
 "src/openai_ns_reconstruction/original_switch_post_finite_N_I1_numeric.py",
 "src/openai_ns_reconstruction/original_switch_post_pressure_source.py",
 "scripts/run_original_switch_post_finite_N_I1.py",
)
REVIEWED={
 CODE[0]:"8bd1c99f39f82d1789566238ea3f6e699bc5f5f29bb32c69b2a47c662dad967f",
 CODE[1]:"dc101ab317dc3c78637419853799baa45558517c852c2d4ee7f97bd040bbfdd6",
 CODE[2]:"96d11eaad4db794d6e1d27d1fc9ac2bdcb9329c54a5489fc5e116054c742b1e4",
 CODE[3]:"c5bf5ad0ad21757105c940d8b734219405b691aa46aa51a179865bb6a1aec7d9",
 CODE[4]:"bafff5c2134ef6d7126f5e91219049d702bde153e98c4712a7c67d70334e4629",
}

def require(value, message):
    if not value:
        raise ValueError(message)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def load(path):
    raw = Path(path).read_bytes()
    if str(path).endswith(".gz"):
        raw = gzip.decompress(raw)
    return json.loads(raw)


def parent():
    path = BASE / PARENT_NAME
    require(path.is_file(), "Accepted macro finite-N parent is missing")
    require(sha(path) == PARENT_SHA256, "Accepted macro finite-N parent bytes changed")
    result = load(path)
    require(result.get("checks", {}).get("passed") is True,
            "Accepted macro finite-N parent did not pass")
    return result


def check_pins(report):
    count = 0
    for folder, group in ((BASE, report["input_hashes"]),
                          (ROOT, report["implementation_hashes"])):
        require(isinstance(group, dict), "Report source pin table is missing")
        for name, expected in group.items():
            require(isinstance(expected, str), "Source pin is not a SHA256: " + name)
            require(sha(folder / name) == expected,
                    "Source bytes changed: " + name)
            count += 1
    return count


def merge_pins(accepted, *packets):
    saved = dict(accepted["input_hashes"])
    implementations = dict(accepted["implementation_hashes"])
    saved[PARENT_NAME] = PARENT_SHA256
    for packet in packets:
        pins = packet.get("source_pins")
        require(isinstance(pins, dict) and pins,
                "Every first-switch packet must carry source_pins")
        for key, item in pins.items():
            path = Path(item.get("path", key) if isinstance(item, dict) else key)
            expected = item.get("expected_sha256") if isinstance(item, dict) else item
            require(isinstance(expected, str), "Missing source pin SHA256: " + str(path))
            if len(path.parts) == 1:
                actual = BASE / path
                name, target = path.name, saved
            elif path.as_posix().startswith("experiments/root_st073/"):
                actual = ROOT / path
                name, target = path.name, saved
            else:
                actual = ROOT / path
                name, target = path.as_posix(), implementations
            require(actual.is_file(), "Pinned source is missing: " + path.as_posix())
            require(sha(actual) == expected,
                    "New packet source bytes changed: " + path.as_posix())
            require(target.get(name, expected) == expected,
                    "Conflicting source pin: " + name)
            target[name] = expected
    for path in CODE:
        require((ROOT / path).is_file(), "New implementation is missing: " + path)
        implementations[path] = sha(ROOT / path)
    return saved, implementations


def check_records(value, path="", cached_source_evidence=False):
    """Check directed binary ordering without expanding large powers."""
    from mpmath.libmp import mpf_cmp

    if path.endswith((".charts.second_switch.source_proof.cache_source.geometry",
                      ".sourceproof.second_switch_cache.second_switch.geometry")):
        source=BASE/"lei_ren_part1_paper_compliant_current_original_second_switch_R110.json"
        require(sha(source)=="51b930434d4746c901eac1328eb0f1b8eae3ea2abb9fceca5cb8dc0b0ecc3751",
                "Legacy source geometry bytes changed")
        original=load(source)["second_switch_packets"]["0"][0]["function_evaluation"]["geometry"]
        require(value==original,"Copied legacy geometry differs from its exact saved source")
        # This is provenance data, not a new directed numerical enclosure.
        cached_source_evidence=True
    count = 0
    if isinstance(value, dict):
        if "lower_exact_mpf_tuple" in value and "upper_exact_mpf_tuple" in value:
            lo, hi = value["lower_exact_mpf_tuple"], value["upper_exact_mpf_tuple"]
            require(
                all(
                    isinstance(item, list)
                    and len(item) == 4
                    and all(type(number) is int for number in item)
                    and item[0] in (0, 1)
                    and item[1] >= 0
                    and item[3] >= 0
                    for item in (lo, hi)
                ),
                "Malformed or nonfinite directed interval",
            )
            directed = value.get("directed") is True
            legacy_exact_source = (
                path.endswith(".R100_geometry_certificate.actual_T_singleton_tuple")
                and "directed" not in value and lo == hi
            )
            require(mpf_cmp(tuple(lo), tuple(hi)) <= 0 and
                    (directed or legacy_exact_source or
                     (cached_source_evidence and "directed" not in value)),
                    "Reversed or nondirected interval: " + path)
            count += int(directed and not cached_source_evidence)
        count += sum(check_records(item, path + "." + str(key), cached_source_evidence)
                     for key, item in value.items())
    elif isinstance(value, (list, tuple)):
        count += sum(check_records(item, path + "." + str(key), cached_source_evidence)
                     for key, item in enumerate(value))
    return count


def _orders(rows, label, expected):
    require(isinstance(rows, list), label + " must be a list")
    require([row.get("ordinary_Z_order") for row in rows
             if isinstance(row, dict)] == list(expected),
            "Missing or unordered ordinary Z rows: " + label)
    return rows


def _record(row, field="value"):
    require(isinstance(row, dict), "Directed row is not a mapping")
    value = row.get(field)
    require(isinstance(value, dict), "Missing directed row field: " + field)
    require("lower_exact_mpf_tuple" in value and "upper_exact_mpf_tuple" in value,
            "Incomplete directed row field: " + field)
    return value


def _compare_mpf(left, right):
    from mpmath.libmp import mpf_cmp
    return mpf_cmp(tuple(left), tuple(right))


def _require_parent_packet(packet, label):
    require(isinstance(packet, dict) and packet.get("passed") is True,
            label + " packet did not pass")
    require(packet.get("parent_report_sha256") == PARENT_SHA256,
            label + " packet has the wrong parent SHA")
    require(isinstance(packet.get("source_pins"), dict)
            and packet["source_pins"],
            label + " source_pins are required")
    require(packet.get("physical_y_derivatives_claimed",
                       packet.get("limits", {}).get("physical_y_derivatives_claimed")) is False,
            label + " made a physical-y derivative claim")


def validate_background(packet,accepted):
    _require_parent_packet(packet,"Switch/post background")
    require(set(packet["charts"])==set(CHARTS),"Wrong background chart set")
    for name in CHARTS:
        chart=packet["charts"][name]
        require(chart["passed"] is True,"Unproved background chart: "+name)
        for field in ("phi","U","M","A","B","C"):
            _orders(chart["native_fields"][field],name+" native "+field,range(5))
            _orders(chart["incoming_fields"][field],name+" incoming "+field,range(5))
        _orders(chart["Q_b"],name+" Q",range(3))
        distinction=chart["source_proof"]["integral_vs_pointwise"]
        require(distinction["Q_rows_are_pointwise_decay_enclosures"] is True
                and distinction["integrated_q_rows_are_Iq_rows"] is True,
                "Pointwise Q must remain distinct from its weighted integral")
        geometry=chart["geometry"]
        require(geometry["physical_measure"]==
                ("R^2*hbS*dt / Wswitch" if name=="second_switch" else
                 "R^2*post_length*dp / Wmacro")
                and geometry["physical_J_name"]==
                ("hbS" if name=="second_switch" else "post_length"),
                "Native chart Jacobian or physical measure changed")
        require(geometry["normalization"]==("Wswitch" if name=="second_switch" else "Wmacro")
                and len(geometry["normalized_moments"])==3
                and geometry["mass"]==geometry["normalized_moments"][0],
                "Background analytic moments use the wrong physical mass")
        for value in geometry["normalized_moments"]:
            require(_compare_mpf(value["lower_exact_mpf_tuple"],[0,0,0,0])>0,
                    "Analytic physical moment must have a positive floor")
        I=chart["background_I1"]
        require(I["pointwise_Q_distinct_from_integrated_Iq"] is True
                and I["Q_rows_are_pointwise_decay_enclosures"] is True
                and I["integrated_q_formula"]=="I[q]: q=alpha*J0+beta*J1",
                "Integrated I[q] cannot provide a pointwise product bound")
        require(set(I["six_integrals"])==set(("q","q_Z","q_squared","U_q_Z","U_squared","U_U_Z")),
                "Incomplete analytic six-I1 representation")
        for key,rows in I["six_integrals"].items():
            _orders(rows,name+" background "+key,range(2))
        _orders(I["D_rows"],name+" background D",range(2))
        require("two_h_Z2_over_L2" in I["coefficientsC1rows"]
                and I["algebra_receipt"]["U_squared_coefficient"]=="2*h*Z^2/Lz^2",
                "D must explicitly use 2*h*Z^2/L^2")
        _orders(chart["background_B_rows"],name+" background B",range(2))
        _orders(chart["background_S1_rows"],name+" background S1",range(2))
        boundary=chart["background_pressure_algebra"]
        require(boundary["normalized_formula"]=="B/W=2*C0*J0+C1*J1"
                and boundary["C2_cancelled_before_interval_boxing"] is True
                and boundary["C2_not_used_in_B_rows"] is True
                and boundary["no_endpoint_subtraction"] is True
                and boundary["no_inverse_tiny_hb_boxing"] is True
                and boundary["S1_formula"]=="S1=B+D"
                and boundary["normalization"]==geometry["normalization"]
                and boundary["W_formula"]==
                ("10000*hbS" if name=="second_switch" else "10000")
                and boundary["normalized_moments"]["J0"]==geometry["normalized_moments"][0]
                and boundary["normalized_moments"]["J1"]==geometry["normalized_moments"][1],
                "Correlated background B algebra or physical units changed")
    def same_rows(left,right):
        return all(_record(l)[side]==_record(r)[side]
                   for l,r in zip(left,right)
                   for side in ("lower_exact_mpf_tuple","upper_exact_mpf_tuple"))
    first=accepted["actual_first_switch_background_I1_numeric"]["native_fields"]
    second=packet["charts"]["second_switch"]
    post=packet["charts"]["post_switch"]
    for field in ("phi","U","M","A","B","C"):
        require(same_rows(second["incoming_fields"][field],first[field])
                and same_rows(post["incoming_fields"][field],second["native_fields"][field]),
                "Native source inlet was reset or replaced: "+field)
    require(post["source_proof"]["incoming_M_endpoint"]["includes_second_endpoint"] is True
            and post["source_proof"]["no_reset"] is True,
            "Post mean must retain the actual second endpoint source")
    source=post["source_proof"]["post_switch"]
    required_strings={
        "analytic_M_formula":"M_R2=U+(M_R1-U)*exp(-hbS)",
        "memory_M_formula":"M_R1*exp(-hbS)+U*(1-exp(-hbS))",
        "B_source_convolution":"phi_R2^2*(5/6)*(exp(-4*L/5)-exp(-2*L))",
        "C_source_convolution":"phi_R2^2*5*(exp(-4*L/5)-exp(-L))",
        "native_M_hull":"hull(incoming_M,U)",
        "native_A_hull":"hull(incoming_A,U^2)",
        "native_B_hull":"hull(incoming_B,phi^2_uniform/2)",
        "native_C_hull":"hull(incoming_C,phi^2_uniform)",
    }
    require(all(source.get(key)==value for key,value in required_strings.items())
            and all(source.get(key) is True for key in (
                "native_histories_are_full_segment_hulls", "scalar_width_jet_is_constant",
                "endpoint_rows_retained_separately", "positive_source_kernel_factors"))
            and source["B_rate"]==2 and source["C_rate"]==1,
            "Post endpoint kernel or full-segment native hull proof changed")
    require(packet["geometry"]["weight_derivative_fixture"]["passed"] is True
            and packet["geometry"]["executable_weight_derivative_fixture"] is True,
            "Missing normalized-weight derivative certificate")
    geometry=packet["geometry"]["source"]
    require(geometry["second_formula"]=="R=100*exp(hbS*(1+t)), t in [0,1]"
            and geometry["post_formula"]==
            "y=post_length*p, R=100*exp(2*hbS+y), post_length=log(11/10)-2*hbS"
            and geometry["actual_weight_interval_not_selected"] is True,
            "Source geometry or signed weight interpretation changed")
    return dict(passed=True,analytic_signed_volume=True,chart_specific_native_rows=True,
                source_histories_not_reset=True,separate_physical_normalizations=True)


def validate_inputs(packet):
    _require_parent_packet(packet,"Switch/post loop input")
    for name in CHARTS:
        chart=packet["charts"][name]
        require(chart["passed"] is True and chart["ordinary_Z_orders"]==[0,1,2,3],
                "Incomplete source-related loop inputs")
        require(_compare_mpf(chart["kappa_log_upper"]["upper_exact_mpf_tuple"],[0,0,0,0])<0,
                "Current source kappa is not below one")
        for key in ("a","b","t0","E_over_frozen_E0","r","alpha","K"):
            _orders(chart["log_caps"][key],name+" "+key,range(4))
        for key in ("b","t0"):
            require(all(row.get("exact_zero") is True and row.get("log_absolute_upper_bound") is None
                        for row in chart["log_caps"][key]),
                    "Zero axial background control must retain exact-zero jets")
    return dict(passed=True,b_and_t0_exact_zero=True,kappa_cut_one=True,
                physical_y_derivatives_claimed=False)


def validate_loop(packet):
    _require_parent_packet(packet,"Switch/post loop bounds")
    for name in CHARTS:
        chart=packet["charts"][name]
        _require_parent_packet(chart,name+" loop")
        _orders(chart["B_over_Pstar_over_frozen_E0_log_caps"],name+" B",range(4))
        proof=chart["inverse_proof"]
        require(proof["passed"] is True and proof["ordinary_Z_derivatives"] is True
                and proof["Z_only"] is True and proof["physical_y_derivatives"] is False
                and chart["universal_B_over_E_value_cap_tightened"] is False,
                "Loop inverse proof changed derivative/conditioning scope")
        scale=chart["physical_scale_proof"]
        require(scale["B_quantity"]=="B/Pstar/frozen_E0"
                and scale["source_bound_to_current_switch_post_packet"] is True
                and scale["actual_N18"] is True and scale["phase_Z_independent"] is True,
                "Physical loop normalization or common N source changed")
    return dict(passed=True,ordinary_Z3=True,no_generic_C0_tightening=True)


def validate_numeric(packet):
    _require_parent_packet(packet,"Switch/post finite-N")
    require(packet["accepted_ancestors_replayed"] is False
            and packet["all_normalized_I1_errors_below_2_power_minus1024"] is True,
            "Finite-N numerical scope changed")
    identity=packet["pressure_B_plus_D_identity"]
    require(identity["passed"] is True and identity["symbolic_remainder"]=="0"
            and identity["finite_N_Uy_not_assumed_zero"] is True
            and identity["Q4_reduced_pressure_source_guarded"] is True,
            "S1=B+D requires the actual finite-N Q4 pressure identity")
    for name in CHARTS:
        chart=packet["charts"][name]
        require(chart["passed"] is True,"Incomplete finite-N chart")
        I,D=chart["I1_full"],chart["D_full"]
        switch=name=="second_switch"
        require(I["normalized_by_Wswitch"]==switch and I["normalized_by_Wmacro"]==(not switch)
                and D["computed_from_full_I1_once"] is True
                and D["separate_D_error_not_added"] is True,
                "Finite-N radius repeated or wrong chart physical units")
        require(chart["signed_I1_radius"]["lower_exact_mpf_tuple"]==[0,1,-1024,1]
                and chart["signed_I1_radius"]["upper_exact_mpf_tuple"]==[0,1,-1024,1],
                "Wrong current-N error radius")
        threshold=chart["I1_error_threshold"]["lower_exact_mpf_tuple"]
        for key,rows in chart["six_I1_difference_error_log_bounds"].items():
            _orders(rows,name+" error "+key,range(2))
            for row in rows:
                require(row.get("exact_zero") is True or
                    _compare_mpf(row["log_absolute_upper_bound"]["upper_exact_mpf_tuple"],threshold)<0,
                    "Normalized finite-N error exceeds 2^-1024")
        require(chart["source_contract"]["incoming_delta_U_not_added"] is True
                and chart["source_contract"]["nonreset_positive_kernel_mean"] is True,
                "Physical U or mean inlet was duplicated/reset")
        B,S1=chart["B_full"],chart["S1_full"]
        _orders(B["rows"],name+" B full",range(2))
        _orders(S1["rows"],name+" S1 full",range(2))
        require(B["normalized_by_Wswitch"]==switch
                and B["normalized_by_Wmacro"]==(not switch)
                and B["symmetric_radius_added_once"]==chart["signed_I1_radius"]
                and B["correlated_background_not_endpoint_interval_subtraction"] is True
                and B["endpoint_error_formula"]==
                "abs(deltaB)/W<=2*R_upper^2/W*sup(abs(deltaG))"
                and B["second_W_has_actual_hbS_in_denominator"]==switch
                and S1["computed_from_B_and_D_once"] is True
                and S1["separate_S1_error_not_added"] is True,
                "B/S1 correlation, normalization or error radius changed")
        _orders(B["finite_N_absolute_error_log_bounds"],name+" B error",range(2))
        for row in B["finite_N_absolute_error_log_bounds"]:
            require(row.get("exact_zero") is True or
                _compare_mpf(row["log_absolute_upper_bound"]["upper_exact_mpf_tuple"],threshold)<0,
                "Normalized finite-N B error exceeds 2^-1024")
    joined=packet["joined_macro_switch_post_volume"]
    require(joined["normalized_by_Wmacro"] is True
            and joined["complete_pressure_enclosed"] is False
            and joined["post_physical_weight_has_no_extra_hbS"] is True
            and joined["I1"]["macro_and_first_radii_not_repeated"] is True
            and joined["I1"]["post_finite_N_radius_not_repeated"] is True
            and joined["I1"]["signed_center_value_claimed"] is False,
            "Joined volume normalization or error accounting changed")
    require(packet["limits"]["complete_pressure_global_closure"] is False
            and packet["limits"]["outer_RQy_retained"] is True
            and packet["limits"]["chart_B_and_S1_numeric"] is True
            and packet["limits"]["joined_parent_B_numeric"] is False,
            "Volume enclosure was promoted to complete pressure")
    return dict(passed=True,finite_N_radius_once=True,joined_Wmacro=True,
                full_pressure_numeric=False,global_closure=False)


def validate_pressure(packet,accepted):
    _require_parent_packet(packet,"Switch/post pressure source")
    graph=packet["graph_patch"]
    require(graph["base_node_count"]==184898 and graph["base_prefix_unchanged"] is True
            and graph["new_node_start"]==184898
            and graph["appended_node_count"]==len(graph["appended_nodes"]),
            "Full accepted first-switch graph prefix changed")
    require(graph["base_prefix_canonical_nodes_sha256"]==
            "36ce7b65b98546cd14c644299d4b4527002dd45c1cd8230f33f503156ddabaf1",
            "Wrong full source prefix digest")
    limits=packet["limits"]
    require(limits["outer_RQy_retained"] is True
            and limits["outer_RQy_evaluated"] is False
            and limits["complete_pressure_global_closure"] is False
            and limits["all_caps_retained"] is True
            and limits["old_212_terms_and_errors_unchanged"] is True,
            "Pressure endpoints or global closure scope changed")
    seams=packet["endpoint_identity_receipts"]
    for name in ("first_upper_to_second_lower","second_upper_to_post_lower"):
        seam=seams[name]
        require(seam["status"]=="PROVED_FORMAL_CONGRUENCE"
                and seam["source_cache_semantic_rewrite"] is True
                and seam["arbitrary_symbolic_jet_proof"] is True
                and seam["arithmetic_fixture_not_proof"] is True
                and seam["finite_N_Qy_identity_proven"] is True,
                "Internal endpoint contraction lacks an actual source proof")
        proof=seam["formal_congruence"]
        require(proof["passed"] is True and proof["common_N"]==18
                and proof["R100_endpoint_bindings"]["incoming_actual_U_and_moments_bound"] is True
                and set(proof["Q4_port_bindings"])==
                {"first_switch","second_switch","post_power"}
                and proof["C5_total_radial_slow_over_J_and_Nfast_bound"] is True
                and proof["symbolic_endpoint_identity"]["all_symbolic_differences_zero"] is True
                and all(value=="0" for value in
                        proof["symbolic_endpoint_identity"]["differences"].values())
                and proof["symbolic_rational_congruence"]
                    ["Q4_difference_after_source_substitution"]=="0"
                and all(value=="0" for value in proof["symbolic_rational_congruence"]
                        ["C5_differences_after_source_substitution"])
                and all(value["passed"] is True and value["source_executed"] is False
                        for value in proof["source_formulas"].values()),
                "Actual source, cache or arbitrary-symbol endpoint proof changed")
        periodic=proof["periodic_phase_analytic_lemma"]
        require(periodic["passed"] is True and periodic["K_normalizes_Phi"] is True
                and periodic["B_phase_drift_after_t0"]=="0"
                and periodic["A_phase_drift"]=="0"
                and periodic["smooth_periodic_integer_extension"] is True,
                "Common fractional phase lacks its analytic periodic extension")
    joined=packet["joined_pressure_source"]
    require(joined["contraction_status"]=="PROVED_FORMAL_ENDPOINT_CONGRUENCES"
            and isinstance(joined["joined_contracted"],list)
            and isinstance(joined["joined_uncontracted"],list)
            and limits["first_switch_upper_standalone_retained"] is True
            and limits["second_switch_upper_standalone_retained"] is True
            and limits["post_power_upper_standalone_retained"] is True,
            "Source join lost uncontracted or standalone endpoint evidence")
    linear=joined["linearization_receipt"]["final"]
    require(linear["passed"] is True and linear["scope"]=="sum_and_negative_rows_only"
            and linear["targets"]==[148514,148515,154676,154677]
            and linear["removed_target_coefficients_zero"] is True
            and linear["surviving_target_dependencies_empty"] is True
            and joined["contracted_before_normalization_preserved"] is True,
            "Same-source endpoint identities were not cancelled before graph boxing")
    _orders(linear["rows"],"source exact signed linearization",range(2))
    for k,row in enumerate(linear["rows"]):
        require(row["exact_integer_arithmetic"] is True
                and all(value==0 for value in row["target_coefficients"].values())
                and not row["surviving_target_dependencies"]
                and row["input_handle"]==joined["contracted_before_linearization"][k]
                and row["reduced_handle"]==joined["joined_contracted"][k],
                "Internal RQy target survived exact signed graph normalization")
    return dict(passed=True,source_prefix_unchanged=True,
                complete_pressure_numeric=False,external_endpoints_retained=True)


def validate_report(report,accepted):
    mutable={"stage","parent_report","parent_report_sha256","input_hashes",
             "implementation_hashes","contract","checks","execution_seconds","execution_audit"}
    for key,value in accepted.items():
        if key not in mutable:
            require(report.get(key)==value,"Accepted mathematical payload changed: "+key)
    require(report["parent_report_sha256"]==PARENT_SHA256
            and report["parent_report"]==PARENT_NAME,"Wrong accepted predecessor")
    require(report["execution_audit"]["inherited_parent_execution_audit"]==accepted["execution_audit"]
            and report["execution_audit"]["current_turn_accepted_ancestors_replayed"] is False,
            "Current/inherited execution audit changed")
    counts={}
    for key in (BACKGROUND_KEY,INPUT_KEY,LOOP_KEY,NUMERIC_KEY,PRESSURE_KEY):
        _require_parent_packet(report[key],key)
        counts[key]=check_records(report[key])
        require(counts[key]>0,"Missing new directed records: "+key)
    return dict(background=validate_background(report[BACKGROUND_KEY],accepted),
        inputs=validate_inputs(report[INPUT_KEY]),loop=validate_loop(report[LOOP_KEY]),
        numeric=validate_numeric(report[NUMERIC_KEY]),
        pressure=validate_pressure(report[PRESSURE_KEY],accepted),directed_records=counts,
        accepted_math_payload_preserved=True)


def run():
    started=time.perf_counter()
    for path,expected in REVIEWED.items():
        require(expected!="REVIEW_PENDING" and sha(ROOT/path)==expected,
                "Independent review/source pin pending or changed: "+path)
    accepted=parent()
    check_pins(accepted)
    sys.path.insert(0,str(ROOT/"src"))
    from openai_ns_reconstruction import original_switch_post_background_I1_numeric as bg
    from openai_ns_reconstruction import original_switch_post_loop_Z3_inputs as inputs
    from openai_ns_reconstruction import original_switch_post_loop_Z3_bounds as loop
    from openai_ns_reconstruction import original_switch_post_finite_N_I1_numeric as numeric
    from openai_ns_reconstruction import original_switch_post_pressure_source as pressure
    background=bg.build(ROOT,accepted)
    source=inputs.build(ROOT,accepted,background)
    bounds=loop.build(ROOT,accepted,source)
    values=numeric.build(ROOT,accepted,background,bounds)
    pressure_source=pressure.build(ROOT,accepted,values)
    require(accepted==parent(),"New stage mutated accepted data")
    packets=(background,source,bounds,values,pressure_source)
    saved,implementations=merge_pins(accepted,*packets)
    report=copy.deepcopy(accepted)
    report.update(stage="current_N_second_switch_post_power_analytic_I1_and_correlated_pressure_source",
        parent_report=PARENT_NAME,parent_report_sha256=PARENT_SHA256,
        input_hashes=saved,implementation_hashes=implementations)
    for key,packet in zip((BACKGROUND_KEY,INPUT_KEY,LOOP_KEY,NUMERIC_KEY,PRESSURE_KEY),packets):
        report[key]=packet
    report["execution_audit"]=dict(current_turn_accepted_ancestors_replayed=False,
        current_turn_accepted_producer_owner_norm_recurrence_quadrature_replayed=False,
        new_stage_arithmetic_only=True,
        unaccepted_new_stage_integration_attempts=[dict(exit_code=1,
            phase="record-schema validation after the five new builders",
            reason="Copied legacy source geometry has no directed flag; now exact-source matched and excluded from new numerical counts",
            report_published=False)],
        inherited_parent_execution_audit=copy.deepcopy(accepted["execution_audit"]))
    report["contract"]=dict(second_switch_post_power_signed_I1_C1_enclosed=True,
        exact_constant_U_mean_formula=True,physical_normalization_once=True,
        complete_pressure_numeric=False,terminal_five_moment_closure=False,full_goal_complete=False)
    report["checks"]=dict(passed=True,packet_invariants=validate_report(report,accepted),
        independent_new_math_review_completed=True,
        independent_runtime_worker_model_metadata_available=False,
        reviewed_implementation_sha256=REVIEWED)
    report["execution_seconds"]=round(time.perf_counter()-started,3)
    json.dumps(report,allow_nan=False)
    return report


def inspect(path):
    report,accepted=load(path),parent()
    require(report["checks"]["passed"] is True
            and report["checks"]["reviewed_implementation_sha256"]==REVIEWED,
            "Saved report review bindings changed")
    for path,expected in REVIEWED.items():
        require(report["implementation_hashes"][path]==expected==sha(ROOT/path),
                "Reviewed source bytes changed")
    require(report["implementation_hashes"][CODE[-1]]==sha(ROOT/CODE[-1]),"Runner source changed")
    return dict(passed=True,pins=check_pins(report),packet_invariants=validate_report(report,accepted),
        serialized_only=True,accepted_builders_producers_or_checkers_replayed=False)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--write-report",action="store_true")
    group.add_argument("--inspect-report",action="store_true")
    args=parser.parse_args()
    path=BASE/NAME
    if args.write_report:
        report=run()
        tmp=path.with_suffix(".json.tmp")
        tmp.write_text(json.dumps(report,indent=2,ensure_ascii=False,allow_nan=False)+"\n",
                       encoding="utf-8",newline="\n")
        tmp.replace(path)
        print(json.dumps(dict(passed=True,execution_seconds=report["execution_seconds"],
                              bytes=path.stat().st_size,sha256=sha(path))))
    else:
        print(json.dumps(inspect(path),indent=2))
