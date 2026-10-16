"""Produce or inspect the source-bound reference native-derivative packet.

``--write-report`` is the only mode that constructs
``OriginalReferenceNativeM15`` and runs the focused synthetic fixture.  The
inspection path reads the saved packet and pinned parent phase report only;
it does not construct the owner or replay an accepted ancestor.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
import sys
import time


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
SCRIPTS = ROOT / "scripts"
DIRECTORY = ROOT / "experiments" / "root_st073"
for path in (str(ROOT), str(SRC), str(SCRIPTS)):
    if path not in sys.path:
        sys.path.insert(0, path)

from openai_ns_reconstruction.original_core_numeric import PREFIX  # noqa: E402
from openai_ns_reconstruction.original_reference_native_m15 import (  # noqa: E402
    GERM,
    GERM_SHA,
    PHASE,
    PHASE_SHA,
    OriginalReferenceNativeM15,
)


NAME = PREFIX + "compliant_current_original_reference_native_m15.json"
FIXTURE_NAME = "scripts/check_original_reference_native_m15.py"
RUNNER_NAME = "scripts/run_original_reference_native_m15.py"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical(value) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def _require(condition, message):
    if not condition:
        raise AssertionError(message)


def _target(base: Path, name: str) -> Path:
    path = Path(name)
    return path if path.is_absolute() else base / path


def _check_hash_group(base: Path, group, label: str) -> None:
    _require(isinstance(group, dict) and group, f"{label} hash map missing")
    for name, expected in group.items():
        _require(isinstance(name, str) and isinstance(expected, str) and len(expected) == 64,
                 (label, name, "malformed digest"))
        path = _target(base, name)
        _require(path.is_file(), (label, name, "missing pinned file"))
        _require(sha(path) == expected, (label, name, "hash changed"))


def _load_phase() -> dict:
    """Load the accepted phase packet without constructing any owner."""
    path = DIRECTORY / PHASE
    _require(path.is_file() and sha(path) == PHASE_SHA, "accepted phase report changed")
    parent = json.loads(path.read_text(encoding="utf-8"))
    _require(parent.get("checks", {}).get("passed") is True,
             "accepted phase report is not passed")
    _check_hash_group(DIRECTORY, parent.get("input_hashes"), "phase input")
    _check_hash_group(ROOT, parent.get("implementation_hashes"), "phase implementation")
    return parent


def _check_ledger(ledger, parent) -> dict:
    expected = parent.get("unchanged_m15_ledger")
    _require(isinstance(expected, list) and len(expected) == 3,
             "accepted phase ledger shape changed")
    _require(ledger == expected, "raw/reduced/target ledger changed")
    for row in ledger:
        _require(isinstance(row, dict), ("ledger row", row))
        _require(row.get("raw_terms") == 213 and row.get("reduced_terms") == 199,
                 ("ledger counts", row))
        for key in ("raw_ledger_sha256", "reduced_ledger_sha256", "independent_target_sha256"):
            _require(isinstance(row.get(key), str) and len(row[key]) == 64,
                     ("ledger digest", key, row))
    _require(parent.get("checks", {}).get("raw_213_reduced_199_and_target_unchanged") is True,
             "accepted phase raw/reduced contract changed")
    _require(parent.get("checks", {}).get("all_old_error_budgets_unchanged") is True,
             "accepted phase budget contract changed")
    return {"samples": len(ledger), "raw_terms": 213, "reduced_terms": 199,
            "independent_target_unchanged": True, "all_old_budgets_unchanged": True}


def _check_tangents(tangents) -> dict:
    _require(isinstance(tangents, dict) and tangents.get("passed") is True,
             "root tangent proof failed")
    expected = {"a_y", "b_y", "t0_y", "Delta_y", "q_y"}
    rows = tangents.get("rows")
    _require(isinstance(rows, list) and len(rows) == 5, "root tangent row count changed")
    _require({row.get("quantity") for row in rows} == expected,
             "root tangent quantities changed")
    _require(all(row.get("exact_zero") is True for row in rows),
             "root tangent is not exact zero")
    _require(tangents.get("Z_values_selected") is False,
             "root tangent selected an axial value")
    mixed=tangents.get('ordinary_Z1_rows')
    _require(isinstance(mixed,list) and len(mixed)==5 and all(r.get('exact_zero') is True and
        r.get('ordinary_Z_derivative_taken_before_source_substitutions') is True for r in mixed),'Mixed root-yZ source substitution missing')
    return {"passed": True, "rows": 5, "exact_zero_rows": 5, "mixed_yZ_rows":5}


def _check_response(response,phase) -> dict:
    _require(isinstance(response, dict), "common physical response missing")
    ports = response.get("ports")
    required = {"E", "Ez", "p", "pz", "py", "pyz", "q", "dstar"}
    _require(isinstance(ports, dict) and set(ports) == required,
             "common response ports changed")
    for name in ("E", "Ez", "p", "pz", "q"):
        port = ports[name]
        _require(port.get("kind") == "original_endpoint_function" and
                 isinstance(port.get("original_node"), int) and
                 port.get("native_endpoint") == 1 and
                 isinstance(port.get("native_coordinate"), int),
                 ("original response port", name, port))
    for name in ("py", "pyz"):
        port = ports[name]
        _require(port.get("kind") == "accepted_physical_recipe_function" and
                 port.get("parent_report") == GERM and
                 port.get("parent_report_sha256") == GERM_SHA and
                 isinstance(port.get("local_recipe_node"), int) and
                 port.get("not_an_original_graph_node") is True,
                 ("accepted local response port", name, port))
    _require(ports['q'].get('ordinary_Z1_zero_requires_actual_source_certificate') is True and
        isinstance(ports['q'].get('ordinary_Z1_node'),int),'Actual q_Z1 source port missing')
    _require(ports["dstar"].get("kind") == "original_source_parameter" and
             isinstance(ports["dstar"].get("original_node"), int),
             ("source parameter response port", ports["dstar"]))
    germ_path=DIRECTORY/GERM
    _require(sha(germ_path)==GERM_SHA,'Accepted physical recipe changed')
    germ=json.loads(germ_path.read_text(encoding='utf-8'))
    common={(r['group'],r['quantity'],r['ordinary_Z_order']):r['left_node'] for r in germ['common_source_identity_rows']}
    expected_original={'E':common['fields','E',0],'Ez':common['fields','E',1],
        'p':common['root_values','p2',0],'pz':common['root_values','p2',1],
        'q':phase['phase_lift_certificate']['graph_guards']['charts']['reference']['q_chain']['handles']['q0']}
    for name,node in expected_original.items():
        _require(ports[name]['original_node']==node,('wrong actual source port',name,node))
        _require(ports[name]['ordinary_Z_order']==(1 if name in ('Ez','pz') else 0),('wrong unscaled Z port',name))
    field_port=germ['physical_germ_recipes']['recipe_nodes'][germ['physical_germ_recipes']['common_fields']['y0']['E'][0]]
    _require(all(ports[name]['native_coordinate']==field_port['original_coordinate'] for name in expected_original),'Wrong native source coordinate')
    for k,name in enumerate(('py','pyz')):
        _require(ports[name]['local_recipe_node']==germ['physical_germ_recipes']['root_p2']['physical_y_Z0_Z1'][k] and
                 ports[name]['ordinary_Z_order']==k,('wrong actual physical p2 derivative port',name))
    _require(ports['dstar']['original_node']==germ['source_guards']['parameter_handles']['d_star'],'Wrong canonical dstar')
    definitions = response.get("definitions")
    _require(isinstance(definitions, list) and any("G=slow_B_y" in row for row in definitions) and
             any("G_Z=" in row for row in definitions), "common response definitions missing")
    _require(response.get("ordinary_Z_derivatives_unscaled") is True and
             response.get("phase_Z_exact_zero") is True and
             response.get("true_Z2_field_ports_consumed_by_pyz") is True and
             response.get("no_actual_phase_or_source_values_selected") is True,
             "common response source contract changed")
    _require("arbitrary actual fractional_part" in response.get("held_fast_phase", ""),
             "common response fast phase was specialized")
    return {"passed": True, "ports": sorted(required), "accepted_local_recipe_ports": ["py", "pyz"],
            "ordinary_Z1_retained": True, "fast_phase_arbitrary": True}


def _check_source_guards(guards) -> dict:
    _require(isinstance(guards, dict) and guards.get("passed") is True,
             "native derivative source guards failed")
    _require(guards.get("source_only") is True and guards.get("numeric_caps_used") is False and
             guards.get("ancestor_replayed") is False, "native guard provenance changed")
    checks = guards.get("checks")
    _require(isinstance(checks, dict) and checks, "native source guard checks missing")
    _require(all(isinstance(value, dict) and value.get("passed") is True
                 for value in checks.values()), "a native source guard check failed")
    return {"passed": True, "checks": sorted(checks), "numeric_caps_used": False}


def _check_history(history) -> dict:
    _require(isinstance(history, dict) and history.get("passed") is True,
             "corrected history seam failed")
    _require(history.get("source_only") is True and history.get("identity_row_count") == 138,
             "history seam row count/provenance changed")
    rows = history.get("identity_rows")
    _require(isinstance(rows, list) and len(rows) == 138, "history identity rows changed")
    required_orders = history.get("required_orders")
    _require(isinstance(required_orders, dict) and required_orders.get("ordinary_Z_orders_checked") == [0, 1, 2] and
             all(required_orders.get(name) is True for name in (
                 "finite_N_V_N", "finite_N_m_increment", "retained_history",
                 "modified_m_y", "modified_Q_y_Z2_inputs_available")),
             "history ordinary-Z2 source checks changed")
    claims = set(history.get("claims_not_made", ()))
    _require({"native_Jacobian_equality", "modified_U_y_equality", "modified_Q_y_equality",
              "D_native_B_zero"} <= claims, "history helper overclaimed an unproved derivative")
    return {"passed": True, "identity_rows": 138, "ordinary_Z_orders": [0, 1, 2],
            "native_derivative_aliases_not_claimed": True}


def _check_theorem(theorem) -> dict:
    _require(isinstance(theorem, dict) and theorem.get("passed") is True,
             "native derivative theorem failed")
    rows = theorem.get("rows")
    _require(isinstance(rows, list) and len(rows) == 2, "native theorem row count changed")
    _require([row.get("ordinary_Z_order") for row in rows] == [0, 1] and
             all(row.get("quantity") == "slow_B_y" and row.get("passed") is True
                 for row in rows), "native theorem ordinary-Z rows changed")
    _require([row.get("common_response") for row in rows] == ["G", "G_Z"],
             "native theorem response ports changed")
    _require(theorem.get('actual_root_parameter_Z1_certificate_consumed') is True and
             theorem.get('mixed_root_yZ_rows_consumed')==5,'Root-Z1 or mixed-yZ proof not consumed')
    _require(theorem.get("direct_full_modified_Uy_fingerprint") is False and
             theorem.get("arbitrary_fast_phase") is True and
             theorem.get("original_p2_y_not_replaced_by_zero") is True and
             theorem.get("no_error_caps_used") is True,
             "native theorem contract changed")
    _require("source-bound calculus deduction" in theorem.get("proof_method", ""),
             "native theorem method was not analytic source chain rule")
    return {"passed": True, "rows": 2, "ordinary_Z_orders": [0, 1],
            "direct_full_modified_Uy_fingerprint": False, "arbitrary_fast_phase": True}


def _check_trace_recovery(recovery, theorem) -> dict:
    _require(isinstance(recovery, dict) and recovery.get("passed") is True,
             "modified trace recovery failed")
    rows = recovery.get("rows")
    _require(isinstance(rows, list) and len(rows) == 17, "modified trace row count changed")
    ports = recovery.get("authorized_derivative_function_ports")
    _require(isinstance(ports, list) and len(ports) == 4 and len(set(ports)) == 4,
             "authorized derivative port count changed")
    expected = {row[side + "_node"] for row in theorem["rows"] for side in ("left", "right")}
    _require(set(ports) == expected, "authorized derivative ports do not match theorem rows")
    _require(recovery.get("actual_native_endpoint_scope_required") is True and
             recovery.get("original_graph_not_mutated") is True and
             recovery.get("proof_authorized_function_identity_not_unconditional_raw_fingerprint") is True and
             recovery.get("original_R_and_Pstar_source_factors_retained") is True and
             recovery.get("ordinary_Z2_U_mean_retained") is True and
             recovery.get("numeric_values_or_caps_used") is False,
             "trace recovery scope/provenance changed")
    return {"passed": True, "rows": 17, "authorized_derivative_function_ports": 4,
            "endpoint_scope_restricted": True, "graph_not_mutated": True}


def _check_contract(contract) -> dict:
    expected = {
        "accepted_ancestor_constructors_not_run": True,
        "raw_213_reduced_199_and_all_budgets_unchanged": True,
        "source_function_recipe_not_original_numeric_value": True,
        "no_error_cap_as_value": True,
        "actual_phase_not_selected_or_assumed_zero": True,
    }
    _require(isinstance(contract, dict), "native contract missing")
    for key, value in expected.items():
        _require(contract.get(key) is value, ("native contract", key, contract.get(key)))
    return {"raw_213_reduced_199_and_all_budgets_unchanged": True,
            "no_error_cap_subtraction": True, "actual_phase_not_assumed_zero": True}


def _check_limits(limits) -> dict:
    _require(isinstance(limits, dict), "native limits missing")
    expected_true = {
        "actual_native_B_y_naturality_proved",
        "actual_modified_Uy_Qy_RQy_seam_proved",
    }
    for key, value in limits.items():
        _require(value is (key in expected_true), ("native limit", key, value))
    _require(all(limits.get(key) is True for key in expected_true),
             "native theorem limits were not both proved")
    _require(len(limits) > len(expected_true), "native open-limit set missing")
    return {"actual_native_B_y_naturality_proved": True,
            "actual_modified_Uy_Qy_RQy_seam_proved": True,
            "all_other_limits_open": True}


def _check_root_Z(cert,response,phase):
    _require(isinstance(cert,dict) and cert.get('passed') is True,'Actual root-Z1 certificate failed')
    rows=cert.get('rows')
    _require(isinstance(rows,list) and len(rows)==18,'Actual root value/Z1 row coverage missing')
    _require(cert.get('source_values_selected') is False and cert.get('source_only') is True and
        cert.get('numeric_fixtures_used') is False,'Root-Z1 certificate selected source values')
    for name in ('phase_premise_receipt','germ_premise_receipt','parameter_premise_receipt'):
        _require(cert.get(name,{}).get('passed') is True,('Missing consumed root-Z1 premise',name))
    germ=json.loads((DIRECTORY/GERM).read_text(encoding='utf-8'))
    source={(r['quantity'],r['ordinary_Z_order']):r for r in germ['common_source_identity_rows'] if r['group']=='root_values'}
    expected={(name,key,k) for name in ('reference','restoration') for key in ('a','b','t0','Delta') for k in (0,1)}
    expected|={(name,'q',1) for name in ('reference','restoration')}
    seen=set();constants={'a':[4,5],'b':[0,1],'t0':[0,1],'Delta':[-6,5]}
    for row in rows:
        name,key,k=row['chart'],row['quantity'],row['ordinary_Z_order']
        token=(name,key,k);_require(token in expected and token not in seen,('Unexpected/duplicate root-Z1 row',token));seen.add(token)
        value=[0,1] if k==1 else constants[key]
        _require(row.get('exact_ring_row') is True and row.get('exact_value')==value and row.get('expected_value')==value,
                 ('Actual source root derivative not exact',token))
        if key=='q':
            actual=phase['phase_lift_certificate']['graph_guards']['charts'][name]['q_chain']['q_rows'][1]['handle']
        else:actual=source[key,k]['left_node' if name=='reference' else 'right_node']
        _require(row['original_node']==actual,('Wrong actual root-Z1 source node',token))
        if name=='reference' and key=='q':_require(response['ports']['q']['ordinary_Z1_node']==actual,'Response q_Z1 port differs')
    _require(seen==expected,'Incomplete actual root-Z1 grid')
    _require(cert['q_source_normalization']['q_Z1_direct_ring_proof'] is True,'q_Z1 was not proved by source arithmetic')
    _require(response['ports']['q']['ordinary_Z1_zero_requires_actual_source_certificate'] is True,'Response ignores root-Z1 premise')
    return dict(passed=True,actual_source_rows=18,no_source_values_selected=True,q_Z1_direct_source_ring=True)


def _check_packet(packet: dict, phase: dict) -> dict:
    _require(isinstance(packet, dict), "native packet missing")
    _require(packet.get("source_family") == phase.get("source_family"), "source family changed")
    _require(packet.get("parent_report") == PHASE and
             packet.get("parent_report_sha256") == PHASE_SHA, "accepted phase identity changed")
    _require(packet.get("input_hashes", {}).get(PHASE) == PHASE_SHA and
             packet.get("input_hashes", {}).get(GERM) == GERM_SHA,
             "accepted phase/germ pins missing")
    controls=packet.get('proof_negative_controls',{})
    _require(controls.get('passed') is True and len(controls)==8 and all(v is True for v in controls.values()),'Proof authorization or native geometry controls failed')
    ledger = _check_ledger(packet.get("unchanged_m15_ledger"), phase)
    source_guards = _check_source_guards(packet.get("native_derivative_source_guards"))
    history = _check_history(packet.get("corrected_history_seam"))
    tangents = _check_tangents(packet.get("root_tangent_proof"))
    response = _check_response(packet.get("common_physical_response"),phase)
    root_Z=_check_root_Z(packet.get('actual_root_Z1_certificate'),packet['common_physical_response'],phase)
    theorem = _check_theorem(packet.get("native_derivative_theorem"))
    recovery = _check_trace_recovery(packet.get("actual_modified_trace_recovery"),
                                     packet["native_derivative_theorem"])
    response_sha=hashlib.sha256(canonical(packet['common_physical_response']).encode('utf-8')).hexdigest()
    _require(packet['native_derivative_theorem'].get('response_definition_sha256')==response_sha,'Theorem uses wrong common response')
    for row in packet['native_derivative_theorem']['rows']:
        k=row['ordinary_Z_order']
        _require(row['common_definition_sha256']==response_sha,'Theorem row uses wrong response')
        for side,name in (('left','reference'),('right','restoration')):
            actual=packet['native_derivative_source_guards']['checks']['ordinary_Z1_rows']['charts'][name]['rows']['modified.radial.slow_B_y']['handles'][k]
            _require(row[side+'_node']==actual,('Theorem authorizes wrong original function',name,k))
    expected_trace={(name,k) for name,size in (('U',3),('m',3),('m_y',3),('U_y',2),('Q',2),('Q_y',2),('R_Q_y',2)) for k in range(size)}
    observed=[(r['quantity'],r['ordinary_Z_order']) for r in packet['actual_modified_trace_recovery']['rows']]
    _require(set(observed)==expected_trace and len(set(observed))==len(observed),'Incomplete or duplicated trace grid')
    contract = _check_contract(packet.get("contract"))
    limits = _check_limits(packet.get("limits"))
    return {
        "source_guards": source_guards,
        "history_helper": history,
        "root_tangent_proof": tangents,
        "actual_root_Z1_certificate": root_Z,
        "common_physical_response": response,
        "native_derivative_theorem": theorem,
        "actual_modified_trace_recovery": recovery,
        "ledger": ledger,
        "contract": contract,
        "limits": limits,
    }


def _check_fixture(fixture) -> dict:
    _require(isinstance(fixture, dict) and fixture.get("passed") is True,
             "independent native response fixture failed")
    expected_true = (
        "synthetic_parameters_only",
        "H_Hp_Hpp_independent_inverse_differentiation",
        "physical_y_and_ordinary_Z1_independent_differentiation",
        "unequal_J_naturality",
        "slow_derivative_nonzero",
        "omitted_p2_yZ_control_rejected",
        "wrong_J_control_rejected",
        "total_fast_phase_factor_N_control_rejected",
    )
    for key in expected_true:
        _require(fixture.get(key) is True, ("fixture evidence", key, fixture.get(key)))
    _require(fixture.get("original_source_values_evaluated") is False,
             "fixture evaluated original source values")
    _require(isinstance(fixture.get("root_tangents"), dict) and
             fixture["root_tangents"].get("passed") is True,
             "fixture root tangent evidence missing")
    return {"passed": True, "independent_synthetic_response": True,
            "ordinary_Z1": True, "fast_phase_factor_control": True,
            "original_source_values_evaluated": False}


def _proof_controls(owner,packet):
    from openai_ns_reconstruction.original_reference_native_m15 import ProvedNativeEndpointFingerprint
    from openai_ns_reconstruction.original_native_derivative_guards import _radius_log_jacobian_binding
    def rejected(theorem,response):
        try:ProvedNativeEndpointFingerprint(owner.nodes,owner.l,owner.m,theorem,response)
        except ValueError:return True
        return False
    theorem=copy.deepcopy(packet['native_derivative_theorem']);theorem['passed']=False
    _require(rejected(theorem,owner.response),'Unproved theorem was accepted')
    theorem=copy.deepcopy(packet['native_derivative_theorem']);theorem['actual_root_parameter_Z1_certificate_consumed']=False
    _require(rejected(theorem,owner.response),'Missing root-Z1 premise was accepted')
    theorem=copy.deepcopy(packet['native_derivative_theorem']);theorem['rows'][0]['left_node']+=1
    _require(rejected(theorem,owner.response),'Foreign function port was authorized')
    response=copy.deepcopy(owner.response);response['ports']['py']['local_recipe_node']+=1
    _require(rejected(packet['native_derivative_theorem'],response),'Wrong physical recipe was authorized')
    chart=owner.l['charts']['reference'];n=owner.l['N'];origin=owner.l['phase_origin']
    geo=owner.derivative_guards['checks']['held_actual_phase_and_geometry']['charts']['reference']['radius_log_jacobian']
    nodes=list(owner.nodes);i=geo['log_argument'];nodes[i]=dict(nodes[i],arguments=nodes[i]['arguments']+[chart['coordinate']])
    _require(not _radius_log_jacobian_binding(nodes,chart,n,origin)['passed'],'Extra native radius term accepted')
    nodes=list(owner.nodes);i=chart['global_offset'];nodes[i]=dict(nodes[i],arguments=nodes[i]['arguments']+[chart['coordinate']])
    _require(not _radius_log_jacobian_binding(nodes,chart,n,origin)['passed'],'Extra native phase term accepted')
    try:owner.evaluate('0.1')
    except ValueError:pass
    else:raise AssertionError('Source function owner silently evaluated an original numerical field')
    return dict(passed=True,unproved_theorem_rejected=True,missing_root_Z1_premise_rejected=True,foreign_native_port_rejected=True,
        wrong_physical_recipe_rejected=True,extra_native_radius_term_rejected=True,
        extra_native_phase_term_rejected=True,original_numeric_evaluation_rejected=True)


def run() -> dict:
    started = time.perf_counter()
    phase = _load_phase()
    from check_original_reference_native_m15 import run as run_fixture

    fixture = run_fixture()
    fixture_check = _check_fixture(fixture)
    owner = OriginalReferenceNativeM15(precision=100)
    packet = owner.build()
    packet['input_hashes']=dict(owner.hashes)
    packet['proof_negative_controls']=_proof_controls(owner,packet)
    packet_check = _check_packet(packet, phase)

    hashes = dict(owner.hashes)
    code = dict(owner.code_hashes)
    code[FIXTURE_NAME] = sha(ROOT / FIXTURE_NAME)
    code[RUNNER_NAME] = sha(ROOT / RUNNER_NAME)
    checks = {
        "passed": True,
        "independent_fixture": fixture,
        "fixture_invariants": fixture_check,
        "packet_invariants": packet_check,
        "source_theorem_ordinary_Z_rows": 2,
        "root_tangent_zero_rows": 5,
        "mixed_root_yZ_zero_rows": 5,
        "actual_root_Z1_source_rows": 18,
        "common_response_ports": 8,
        "history_helper_identity_rows": 138,
        "actual_trace_recovery_rows": 17,
        "authorized_derivative_function_ports": 4,
        "actual_native_B_y_naturality_proved": True,
        "actual_modified_Uy_Qy_RQy_seam_proved": True,
        "raw_213_reduced_199_ledger_unchanged": True,
        "all_old_budgets_unchanged": True,
        "independent_target_unchanged": True,
        "direct_full_modified_Uy_fingerprint": False,
        "no_error_cap_subtraction": True,
        "endpoint_scope_restricted_function_identity": True,
        "original_graph_not_mutated": True,
        "production_is_new_source_packet_only": True,
    }
    return dict(packet, input_hashes=hashes, implementation_hashes=code,
                checks=checks, execution_seconds=round(time.perf_counter() - started, 3))


def inspect(path: Path) -> dict:
    """Verify hashes and serialized invariants without constructing the owner."""
    data = json.loads(path.read_text(encoding="utf-8"))
    _require(data.get("checks", {}).get("passed") is True, "saved native packet is not passed")
    _check_hash_group(DIRECTORY, data.get("input_hashes"), "input")
    _check_hash_group(ROOT, data.get("implementation_hashes"), "implementation")
    _require(FIXTURE_NAME in data["implementation_hashes"] and
             RUNNER_NAME in data["implementation_hashes"],
             "runner/fixture implementation pins missing")
    phase = _load_phase()
    packet_check = _check_packet(data, phase)
    fixture_check = _check_fixture(data["checks"].get("independent_fixture"))
    checks = data["checks"]
    expected = {
        "source_theorem_ordinary_Z_rows": 2,
        "root_tangent_zero_rows": 5,
        "mixed_root_yZ_zero_rows": 5,
        "actual_root_Z1_source_rows": 18,
        "common_response_ports": 8,
        "history_helper_identity_rows": 138,
        "actual_trace_recovery_rows": 17,
        "authorized_derivative_function_ports": 4,
        "actual_native_B_y_naturality_proved": True,
        "actual_modified_Uy_Qy_RQy_seam_proved": True,
        "raw_213_reduced_199_ledger_unchanged": True,
        "all_old_budgets_unchanged": True,
        "independent_target_unchanged": True,
        "direct_full_modified_Uy_fingerprint": False,
        "no_error_cap_subtraction": True,
        "endpoint_scope_restricted_function_identity": True,
        "original_graph_not_mutated": True,
        "production_is_new_source_packet_only": True,
    }
    for key, value in expected.items():
        actual=checks.get(key)
        _require((actual is value) if isinstance(value,bool) else (type(actual) is type(value) and actual==value),
                 ("saved check changed", key, actual))
    return {"packet_invariants": packet_check, "fixture_invariants": fixture_check,
            "checks": checks, "execution_seconds": data.get("execution_seconds")}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-report", action="store_true")
    parser.add_argument("--inspect-report", action="store_true")
    args = parser.parse_args()
    path = DIRECTORY / NAME
    if args.inspect_report:
        data = inspect(path)
    else:
        data = run()
    if args.write_report:
        _require(not args.inspect_report, "inspection cannot rewrite an accepted report")
        temporary = path.with_suffix(".json.tmp")
        temporary.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n",
                             encoding="utf-8", newline="\n")
        temporary.replace(path)
    print(json.dumps({"checks": data.get("checks", {}),
                      "execution_seconds": data.get("execution_seconds")},
                     indent=2), flush=True)
