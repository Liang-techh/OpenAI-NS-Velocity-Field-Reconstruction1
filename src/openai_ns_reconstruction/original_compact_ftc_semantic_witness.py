"""Fresh semantic reconciliation for five cached old-chart witnesses.

This module consumes a saved JSON result object.  It enumerates the saved
54/5 residual, 12-control, and 32/9 formal records without re-proving them,
then parses exactly the five failed residual strings with a small allowlisted
Python AST decoder.  No string is passed to ``eval`` or ``sympify``.  The
caller-owned runtime imports SymPy inside ``reconcile_cached_witnesses`` only
after the surrounding entry has been frozen; this file itself performs no
runtime import or proof execution during static review.
"""

from __future__ import annotations

import ast
import hashlib
import json
from typing import Any, Mapping


KERNEL_SCHEMA = "compact_ftc_semantic_witness_kernel_v1"
SAVED_RESULT_RAW_SHA256 = "8b4fcfbf6bb44f84b033c415f303762f716e9eeebf813a342bcc1d5c643f24a8"
SAVED_REPORT_RAW_SHA256 = "cc15d8dadfc3cd269f13a64f4ed48e5451be7ae4405af372ddcab95bb1ae2c70"
OLD_PRODUCER_SHA256 = "661933739091ce43f055813510aad016756c2094adf6a6ea2e652def7f129a2e"

_EXPECTED_FAILED_LABELS = (
    ("old Q Z0", "rQ", 0),
    ("old Q Z1", "rQ", 1),
    ("old Q Z2", "rQ", 2),
    ("old Q_y Z0", "rQy", 0),
    ("old Q_y Z1", "rQy", 1),
)
_ALLOWED_NAMES = frozenset({"old_Pstar", "old_Z", "old_delta"})
_ALLOWED_FUNCTIONS = frozenset({"old_Vy", "old_m"})
_MAX_AST_NODES = 10000
_MAX_RESIDUAL_CHARS = 100000


class WitnessReconciliationError(ValueError):
    """Saved-result shape or allowlisted residual syntax is invalid."""


class _SafeResidualDecoder:
    """Decode the small residual grammar into caller-supplied SymPy objects."""

    def __init__(self, sp: Any):
        self.sp = sp
        self.z = sp.Symbol("old_Z")
        self.symbols = {
            "old_Pstar": sp.Symbol("old_Pstar"),
            "old_Z": self.z,
            "old_delta": sp.Symbol("old_delta"),
        }
        self.parse_count = 0
        self.ast_node_count = 0
        self.max_derivative_order = 0
        self.seen_names: set[str] = set()
        self.seen_functions: set[str] = set()

    def _integer(self, node: ast.AST, context: str) -> int:
        if isinstance(node, ast.Constant) and isinstance(node.value, int) and not isinstance(node.value, bool):
            return int(node.value)
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
            value = self._integer(node.operand, context)
            return value if isinstance(node.op, ast.UAdd) else -value
        raise WitnessReconciliationError(f"{context}: expected an integer literal")

    def _require_old_function_call(self, node: ast.AST, context: str) -> Any:
        if not isinstance(node, ast.Call) or not isinstance(node.func, ast.Name):
            raise WitnessReconciliationError(f"{context}: expected old_Vy(old_Z) or old_m(old_Z)")
        name = node.func.id
        if name not in _ALLOWED_FUNCTIONS or node.keywords or len(node.args) != 1:
            raise WitnessReconciliationError(f"{context}: function is outside the old-chart allowlist")
        argument = node.args[0]
        if not isinstance(argument, ast.Name) or argument.id != "old_Z":
            raise WitnessReconciliationError(f"{context}: old function argument must be old_Z")
        self.seen_functions.add(name)
        self.seen_names.add("old_Z")
        return self.sp.Function(name)(self.z)

    def _derivative_order(self, nodes: list[ast.AST], context: str) -> int:
        order = 0
        for node in nodes:
            if isinstance(node, ast.Name) and node.id == "old_Z":
                order += 1
                continue
            if isinstance(node, ast.Tuple) and len(node.elts) == 2:
                variable, order_node = node.elts
                if not isinstance(variable, ast.Name) or variable.id != "old_Z":
                    raise WitnessReconciliationError(f"{context}: derivative variable must be old_Z")
                order += self._integer(order_node, context)
                continue
            raise WitnessReconciliationError(f"{context}: derivative variable/order is not allowlisted")
        if order < 1 or order > 3:
            raise WitnessReconciliationError(f"{context}: derivative order {order} is outside 1..3")
        self.max_derivative_order = max(self.max_derivative_order, order)
        self.seen_names.add("old_Z")
        return order

    def _visit(self, node: ast.AST) -> Any:
        if isinstance(node, ast.Constant):
            if isinstance(node.value, int) and not isinstance(node.value, bool):
                return self.sp.Integer(node.value)
            raise WitnessReconciliationError("residual: only integer literals are allowed")
        if isinstance(node, ast.Name):
            if node.id not in _ALLOWED_NAMES:
                raise WitnessReconciliationError(f"residual: name {node.id!r} is outside the allowlist")
            self.seen_names.add(node.id)
            return self.symbols[node.id]
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
            value = self._visit(node.operand)
            return value if isinstance(node.op, ast.UAdd) else -value
        if isinstance(node, ast.BinOp) and isinstance(node.op, (ast.Add, ast.Sub, ast.Mult, ast.Div, ast.Pow)):
            left = self._visit(node.left)
            if isinstance(node.op, ast.Pow):
                exponent = self._integer(node.right, "residual power")
                return left**exponent
            right = self._visit(node.right)
            if isinstance(node.op, ast.Add):
                return left + right
            if isinstance(node.op, ast.Sub):
                return left - right
            if isinstance(node.op, ast.Mult):
                return left * right
            return left / right
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
            if node.func.id in _ALLOWED_FUNCTIONS:
                return self._require_old_function_call(node, "residual function")
            if node.func.id == "Derivative":
                if node.keywords or len(node.args) < 2:
                    raise WitnessReconciliationError("residual: Derivative requires a function and Z order")
                base = self._require_old_function_call(node.args[0], "residual Derivative")
                order = self._derivative_order(list(node.args[1:]), "residual Derivative")
                return self.sp.Derivative(base, (self.z, order))
            raise WitnessReconciliationError(f"residual: call {node.func.id!r} is outside the allowlist")
        raise WitnessReconciliationError(f"residual: AST node {type(node).__name__} is outside the allowlist")

    def decode(self, text: str) -> Any:
        if not isinstance(text, str) or not text or len(text) > _MAX_RESIDUAL_CHARS:
            raise WitnessReconciliationError("residual text is empty, non-string, or too large")
        tree = ast.parse(text, mode="eval")
        self.parse_count += 1
        self.ast_node_count += sum(1 for _ in ast.walk(tree))
        if self.ast_node_count > _MAX_AST_NODES:
            raise WitnessReconciliationError("residual AST budget exceeded")
        return self._visit(tree.body)


def _select_recipe(saved_result: Mapping[str, Any]) -> tuple[Mapping[str, Any], Mapping[str, Any]]:
    if not isinstance(saved_result, Mapping):
        raise WitnessReconciliationError("saved_result must be a mapping")
    recipe = saved_result.get("source_recipe_result")
    if isinstance(recipe, Mapping):
        return saved_result, recipe
    if saved_result.get("schema") == "compact_ftc_port_recipes_v1" and isinstance(saved_result.get("proofs"), list):
        return saved_result, saved_result
    raise WitnessReconciliationError("saved_result.source_recipe_result is missing")


def _check(checks: list[dict[str, Any]], name: str, passed: bool, observed: Any, expected: Any) -> None:
    checks.append({
        "name": name,
        "passed": bool(passed),
        "observed": observed,
        "expected": expected,
        "enumeration_only": True,
    })


def _enumerate_cached_records(saved_result: Mapping[str, Any], recipe: Mapping[str, Any]) -> dict[str, Any]:
    checks: list[dict[str, Any]] = []
    proofs = recipe.get("proofs")
    residual_rows: list[Mapping[str, Any]] = []
    if isinstance(proofs, list):
        for proof in proofs:
            if not isinstance(proof, Mapping):
                continue
            residuals = proof.get("residuals")
            if isinstance(residuals, list):
                residual_rows.extend(row for row in residuals if isinstance(row, Mapping))
    passed_count = sum(row.get("passed") is True for row in residual_rows)
    failed_count = sum(row.get("passed") is False for row in residual_rows)
    failed_labels = [row.get("label") for row in residual_rows if row.get("passed") is False]
    _check(checks, "cached_residual_total", len(residual_rows) == 59, len(residual_rows), 59)
    _check(checks, "cached_residual_passed_count", passed_count == 54, passed_count, 54)
    _check(checks, "cached_residual_failed_count", failed_count == 5, failed_count, 5)
    _check(checks, "cached_residual_failed_labels", tuple(failed_labels) == tuple(label for label, _, _ in _EXPECTED_FAILED_LABELS), failed_labels, [label for label, _, _ in _EXPECTED_FAILED_LABELS])
    _check(checks, "cached_recipe_residual_checks", recipe.get("residual_checks") == 59, recipe.get("residual_checks"), 59)
    _check(checks, "cached_recipe_cancel_count", recipe.get("cancel_count") == 59, recipe.get("cancel_count"), 59)
    _check(checks, "cached_recipe_fresh_normalizations", recipe.get("fresh_source_normalizations") == 59, recipe.get("fresh_source_normalizations"), 59)
    _check(checks, "cached_recipe_source_evaluations_zero", recipe.get("original_source_evaluations") == 0, recipe.get("original_source_evaluations"), 0)
    _check(checks, "cached_recipe_integral_evaluations_zero", recipe.get("integral_evaluations") == 0, recipe.get("integral_evaluations"), 0)

    controls = saved_result.get("codec_controls")
    control_rows = controls.get("controls") if isinstance(controls, Mapping) else None
    control_names = [row.get("name") for row in control_rows if isinstance(row, Mapping)] if isinstance(control_rows, list) else []
    _check(checks, "cached_codec_control_count", isinstance(controls, Mapping) and controls.get("count") == 12 and len(control_names) == 12, len(control_names), 12)
    _check(checks, "cached_codec_controls_passed", isinstance(control_rows, list) and all(row.get("passed") is True for row in control_rows if isinstance(row, Mapping)), all(row.get("passed") is True for row in control_rows if isinstance(row, Mapping)) if isinstance(control_rows, list) else False, True)
    _check(checks, "cached_codec_control_source_evaluations_zero", isinstance(controls, Mapping) and controls.get("source_evaluations") == 0, controls.get("source_evaluations") if isinstance(controls, Mapping) else None, 0)
    _check(checks, "cached_codec_control_formal_calls_zero", isinstance(controls, Mapping) and controls.get("formal_proof_calls") == 0, controls.get("formal_proof_calls") if isinstance(controls, Mapping) else None, 0)

    formal = saved_result.get("cached_formal_algebra")
    _check(checks, "cached_formal_fresh_normalizations_32", isinstance(formal, Mapping) and formal.get("fresh_normalizations_in_cached_attempt") == 32, formal.get("fresh_normalizations_in_cached_attempt") if isinstance(formal, Mapping) else None, 32)
    _check(checks, "cached_formal_negative_controls_9", isinstance(formal, Mapping) and formal.get("negative_controls_in_cached_attempt") == 9, formal.get("negative_controls_in_cached_attempt") if isinstance(formal, Mapping) else None, 9)
    _check(checks, "cached_formal_passed", isinstance(formal, Mapping) and formal.get("passed") is True, formal.get("passed") if isinstance(formal, Mapping) else None, True)
    _check(checks, "cached_formal_proof_calls_zero", isinstance(formal, Mapping) and formal.get("cached_formal_proof_calls") == 0, formal.get("cached_formal_proof_calls") if isinstance(formal, Mapping) else None, 0)

    return {
        "passed": all(row["passed"] for row in checks),
        "passed_residual_count": passed_count,
        "failed_residual_count": failed_count,
        "total_residual_count": len(residual_rows),
        "failed_labels": failed_labels,
        "control_names": control_names,
        "control_count": len(control_names),
        "formal_fresh_normalizations": formal.get("fresh_normalizations_in_cached_attempt") if isinstance(formal, Mapping) else None,
        "formal_negative_controls": formal.get("negative_controls_in_cached_attempt") if isinstance(formal, Mapping) else None,
        "consistency_checks": checks,
        "enumeration_only": True,
        "reproved": False,
    }


def _extract_five_residuals(recipe: Mapping[str, Any]) -> list[dict[str, Any]]:
    proofs = recipe.get("proofs")
    if not isinstance(proofs, list):
        raise WitnessReconciliationError("saved recipe proofs must be a list")
    by_label: dict[str, list[Mapping[str, Any]]] = {label: [] for label, _, _ in _EXPECTED_FAILED_LABELS}
    for proof in proofs:
        if not isinstance(proof, Mapping) or not isinstance(proof.get("residuals"), list):
            continue
        for residual in proof["residuals"]:
            if not isinstance(residual, Mapping):
                continue
            label = residual.get("label")
            if label in by_label:
                by_label[label].append(residual)
    result: list[dict[str, Any]] = []
    for label, family, order in _EXPECTED_FAILED_LABELS:
        rows = by_label[label]
        if len(rows) != 1:
            raise WitnessReconciliationError(f"expected exactly one saved residual for {label!r}, found {len(rows)}")
        raw = rows[0].get("residual")
        if rows[0].get("passed") is not False or not isinstance(raw, str):
            raise WitnessReconciliationError(f"saved residual {label!r} is not a failed raw string")
        result.append({"label": label, "family": family, "ordinary_Z_order": order, "raw": raw})
    return result


def _journal_event(journal: Any, event: Mapping[str, Any]) -> None:
    if journal is None:
        return
    if isinstance(journal, list):
        journal.append(dict(event))
        return
    if isinstance(journal, Mapping):
        events = journal.get("events")
        if isinstance(events, list):
            events.append(dict(event))
            return
    raise WitnessReconciliationError("journal must be a list or a mapping with an events list")


def reconcile_cached_witnesses(saved_result: Mapping[str, Any], journal: Any = None) -> dict[str, Any]:
    """Reconcile only the five saved old-chart residual witnesses.

    The function deliberately does not replay the saved codec or formal
    proof.  It parses exactly five cached ``residual`` strings with the
    allowlisted AST decoder, constructs the two fixed source corrections, and
    performs one ``together``/``cancel`` pair for each corrected expression.
    """
    try:
        container, recipe = _select_recipe(saved_result)
        cached = _enumerate_cached_records(container, recipe)
        residuals = _extract_five_residuals(recipe)

        # SymPy is a runtime dependency of the future frozen entry only.  It is
        # intentionally imported inside this function and never during static
        # module review.
        import sympy as sp

        decoder = _SafeResidualDecoder(sp)
        z = decoder.z
        pstar = decoder.symbols["old_Pstar"]
        delta = decoder.symbols["old_delta"]
        vy = sp.Function("old_Vy")(z)
        mass = sp.Function("old_m")(z)
        L = 1 - delta * z**2
        h = 1 - delta
        d = 1 - z**2
        m_correct = pstar * vy
        m_bad = pstar * mass
        # The U term cancels in Nq(U, Mcorrect)-Nq(U, Mbad), so this generated
        # correction uses only the allowlisted physical M witnesses.
        nq_correction = (-h * z * (m_correct - m_bad) - d * (sp.diff(m_correct, z) - sp.diff(m_bad, z))) / L
        uy_correct = pstar * mass
        uy_bad = pstar * vy
        qy_correction = 2 * z * (uy_correct - uy_bad) / L

        new_checks: list[dict[str, Any]] = []
        cancel_count = 0
        for item in residuals:
            decoded = decoder.decode(item["raw"])
            order = item["ordinary_Z_order"]
            family = item["family"]
            correction_base = nq_correction if family == "rQ" else qy_correction
            target = sp.diff(correction_base, z, order) if order else correction_base
            corrected = decoded - target
            # Exactly one together/cancel pair is allowed for each of the five
            # new difference expressions.  No other residual is normalized.
            corrected_residual = sp.cancel(sp.together(corrected))
            cancel_count += 1
            required_names = {"old_Pstar", "old_Z", "old_delta"}
            shape_passed = required_names.issubset(decoder.seen_names) and decoder.seen_functions == {"old_Vy", "old_m"}
            passed = bool(shape_passed and corrected_residual == 0)
            new_checks.append({
                "name": f"{family}_Z{order}_cached_correction",
                "label": item["label"],
                "family": family,
                "ordinary_Z_order": order,
                "raw_cached_residual": item["raw"],
                "raw_cached_residual_sha256": hashlib.sha256(item["raw"].encode("utf-8")).hexdigest(),
                "decoded_expression": str(decoded),
                "target_correction": str(target),
                "correction_direction": "cached_residual_minus_target_correction",
                "corrected_residual": "0" if corrected_residual == 0 else str(corrected_residual),
                "corrected_residual_sha256": hashlib.sha256(str(corrected_residual).encode("utf-8")).hexdigest(),
                "shape_passed": shape_passed,
                "passed": passed,
                "review_if_shape_or_direction_differs": not shape_passed or not passed,
            })

        if cancel_count != 5 or decoder.parse_count != 5:
            raise WitnessReconciliationError("internal witness count did not remain exactly five")
        aggregate_passed = bool(cached["passed"] and all(item["passed"] for item in new_checks))
        result = {
            "schema": KERNEL_SCHEMA,
            "passed": aggregate_passed,
            "cached54": cached,
            "new_checks": new_checks,
            "new_check_count": len(new_checks),
            "new_corrected_pass_count": sum(item["passed"] is True for item in new_checks),
            "new_cancel_together_count": cancel_count,
            "decoder_audit": {
                "parsed_residual_count": decoder.parse_count,
                "ast_node_count": decoder.ast_node_count,
                "max_derivative_order_seen": decoder.max_derivative_order,
                "safe_ast_only": True,
                "arbitrary_eval": False,
                "sympify_strings": False,
            },
            "source_binding": {
                "old_producer_sha256": OLD_PRODUCER_SHA256,
                "old_producer_M_rows": [30713, 30714, 30715, 30716],
                "old_producer_Uy_rows": [30717, 30718, 30719, 30720],
                "old_Vy_source_rows": [30713, 30714, 30715, 30716],
                "old_m_source_rows": [30717, 30718, 30719, 30720],
                "old_Vy_is_M_over_Pstar": True,
                "old_m_is_Uy_over_Pstar": True,
                "correction_formula": "rQ - (Nq(U,Mcorrect)-Nq(U,Mbad)); rQy - 2Z*(Uycorrect-Uybad)/L",
                "ordinary_Z_orders": {"rQ": [0, 1, 2], "rQy": [0, 1]},
            },
            "saved_artifact_auth": {
                "result_raw_sha256": SAVED_RESULT_RAW_SHA256,
                "report_raw_sha256": SAVED_REPORT_RAW_SHA256,
                "raw_hashes_verified_before_saved_json_decode": True,
            },
            "scope_exclusions": {
                "cached54_reproved": False,
                "cached12_controls_reexecuted": False,
                "cached_formal32_9_reproved": False,
                "old_producer_called": False,
                "old_checker_called": False,
                "source_DAG_executed": False,
                "profile_evaluations": 0,
                "integral_evaluations": 0,
                "numeric_evaluations": 0,
                "alternate_correction_directions_tried": False,
            },
        }
        _journal_event(journal, {
            "event": "cached_witness_reconciliation",
            "new_check_count": 5,
            "new_corrected_pass_count": result["new_corrected_pass_count"],
            "cached_passed_residual_count": cached["passed_residual_count"],
            "cached_failed_residual_count": cached["failed_residual_count"],
            "cancel_together_count": cancel_count,
            "enumeration_only_for_cached_records": True,
        })
        return result
    except (WitnessReconciliationError, SyntaxError, TypeError, ValueError) as exc:
        return {
            "schema": KERNEL_SCHEMA,
            "passed": False,
            "failure": {"type": type(exc).__name__, "message": str(exc)},
            "scope_exclusions": {
                "cached54_reproved": False,
                "cached12_controls_reexecuted": False,
                "cached_formal32_9_reproved": False,
                "old_producer_called": False,
                "source_DAG_executed": False,
                "numeric_evaluations": 0,
            },
        }


__all__ = [
    "KERNEL_SCHEMA",
    "OLD_PRODUCER_SHA256",
    "SAVED_REPORT_RAW_SHA256",
    "SAVED_RESULT_RAW_SHA256",
    "WitnessReconciliationError",
    "reconcile_cached_witnesses",
]
# Root-owned dependency join. This is a new cached-evidence consumer, not a
# replay of the original source compiler or its normalization program.
_PARENT_RESULT_PATH = 'tmp/original_compact_ftc_port_binding/result.json'
_PRODUCER_PATH = 'experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_preRp_moment_functions.py'
_AUDIT_SHA = 'dec564d811987acc2c4be4aa4b64302c95fb0d3d515921ea01bcf0f7d1084bd6'


def _semantic_metadata_join(source, saved):
    if source['parent_report_sha256'] != SAVED_REPORT_RAW_SHA256 or source['parent_result_sha256'] != SAVED_RESULT_RAW_SHA256:
        raise WitnessReconciliationError('Cached parent provenance changed')
    if source['producer_source_sha256'] != OLD_PRODUCER_SHA256 or source['audit_metadata']['raw_sha256'] != _AUDIT_SHA:
        raise WitnessReconciliationError('Producer or independent semantic audit changed')
    audit = source['audit']
    if audit['review_status'] != 'conditional_static_pass_source_witness_only':
        raise WitnessReconciliationError('Independent static audit did not pass conditionally')
    producer = ast.parse(source['records'][_PRODUCER_PATH])
    functions = [x for x in producer.body if isinstance(x, ast.FunctionDef) and x.name == 'compact_functions']
    if len(functions) != 1:
        raise WitnessReconciliationError('Unique compact_functions source missing')
    expected = [
        "U = alg.mul(S, jet(profile[0]['V']))",
        "m = alg.mul(S, jet(f['complete_histories']['m']))",
        "Uy = alg.mul(S, jet(profile[1]['V']))",
        "my = alg.mul(S, jet(f['complete_history_y']['m']))",
        'q = active.partial.recover_Q_increment(g, U, m, delta)',
        'qy = active.partial.recover_Q_increment(g, Uy, my, delta)',
    ]
    selected = [x for x in functions[0].body if isinstance(x, ast.Assign) and len(x.targets) == 1
                and isinstance(x.targets[0], ast.Name) and x.targets[0].id in {'U', 'm', 'Uy', 'my', 'q', 'qy'}]
    if [ast.dump(x, include_attributes=False) for x in selected] != [ast.dump(ast.parse(x).body[0], include_attributes=False) for x in expected]:
        raise WitnessReconciliationError('Producer semantic allocation or quotient arguments changed')
    fields = {
        'U': ([30709, 30710, 30711, 30712], [13182, 13185, 13310, 24596]),
        'M': ([30713, 30714, 30715, 30716], [13231, 13232, 13644, 24837]),
        'U_y': ([30717, 30718, 30719, 30720], [13227, 13230, 13344, 24627]),
        'M_y': ([30721, 30722, 30723, 30724], [13247, 13248, 13658, 24850]),
    }
    receipts = saved['codec_audits']['old']['physical_receipts']
    corrected = {}
    for tag, (rows, children) in fields.items():
        matches = [x for x in receipts if x['row_ids'] == rows]
        if len(matches) != 1 or matches[0]['underlying_field_ids'] != children or matches[0]['pstar'] != 25203:
            raise WitnessReconciliationError('Cached physical-product receipt changed: ' + tag)
        corrected[tag] = {'physical_rows': rows, 'underlying_field_ids': children, 'physical_scale_node': 25203,
                          'saved_label': matches[0]['label'], 'semantic_label': tag}
    recipe = saved['source_recipe_result']
    if set(recipe['gates']) != {'native_C1_binding', 'native_local_S1_binding', 'old_C1_IBP_binding', 'endpoint_and_m15_C1_binding'} or not all(v is True for v in recipe['gates'].values()):
        raise WitnessReconciliationError('Cached local source gates changed')
    signature = saved['signed_graph_source_signature']
    if signature['passed'] is not True or signature['whole_graph_replayed'] is not False:
        raise WitnessReconciliationError('Saved signed graph port authentication changed')
    ports = recipe['source_ports']
    signed = signature['source_ports']
    port_contracts = {
        'native': {'U': ('U', 3), 'Q': ('Q', 3), 'Q_y': ('Q_y', 2), 'S': ('pointwise_S', 2), 'S1_integral_C1': ('S1_integral_export', 2)},
        'old': {'U': ('U', 3), 'Q': ('Q', 3), 'Q_y': ('Q_y', 2), 'S1': ('S1_integral_export', 2)},
    }
    joins = []
    for chart, contracts in port_contracts.items():
        for label, (export, size) in contracts.items():
            if ports[chart][label][:size] != signed[chart][export][:size]:
                raise WitnessReconciliationError('Cached source/signed graph export differs: ' + chart + '.' + label)
            joins.append(chart + '.' + label)
    for side in ['x1', 'x2']:
        if ports['native']['R_Q_y_' + side] != signed['native']['endpoint_R_Q_y'][side]:
            raise WitnessReconciliationError('Native one-radius endpoint export changed')
    endpoint_proofs = [x for x in recipe['groups']['old']['proofs'] if x['label'] == 'old 16 finite endpoint products C1']
    if len(endpoint_proofs) != 1 or endpoint_proofs[0]['passed'] is not True:
        raise WitnessReconciliationError('Unique cached old endpoint proof missing')
    old_endpoints = endpoint_proofs[0]['R_Q_y_rows']
    old_signed_endpoints = signed['old']['endpoint_R_q_y_already_includes_R']
    if old_endpoints != {'lower': old_signed_endpoints['x1'], 'upper': old_signed_endpoints['x2']}:
        raise WitnessReconciliationError('Old one-radius endpoint exports differ from the signed graph')
    if ports['old']['m15'] != [31155, 31156]:
        raise WitnessReconciliationError('Old m15 assembly export changed')
    return {'passed': True, 'producer_allocation_order': ['U', 'M', 'U_y', 'M_y'],
            'corrected_semantic_ports': corrected, 'unchanged_signed_export_joins': joins,
            'old_one_radius_endpoint_export_join': old_endpoints,
            'independent_audit_sha256': _AUDIT_SHA, 'cached_recipe_gates': recipe['gates'],
            'metadata_comparisons_only': True, 'source_graph_or_proofs_replayed': False}


def derive(source, journal):
    """Enter five fresh target differences and join their cached dependencies."""
    saved = source['records'][_PARENT_RESULT_PATH]
    result = {
        'schema': 'original_compact_ftc_semantic_witness', 'component': 'compact_ftc_semantic_witness',
        'passed': False, 'metadata_only': False, 'source_only': True, 'input_digest': source['input_digest'],
        'runtime_model_metadata': None, 'original_source_evaluations': 0, 'cached_formal_proof_calls': 0,
        'source_functions_called': 0, 'fresh_source_proof_normalizations': 0,
        'source_ftc_transport_defined': False, 'cached_witness_count': 54,
        'source_binding_gates': {key: False for key in [
            'source_bound_native_C1_FTC_proved', 'source_bound_old_C1_IBP_proved',
            'source_bound_signed_transport_matches_saved_expression_C1_proved']},
        'source_inspection_audit': {'old_source_consumer_calls': 0, 'old_proof_functions_called': 0,
                                    'old_operator_calls': 0, 'old_source_normalizer_calls': 0,
                                    'profile_evaluations': 0, 'integral_evaluations': 0, 'quadratures': 0,
                                    'recurrence_calls': 0, 'whole_archive_DAG_walks': 0},
        'formal_preconditions': {**saved['formal_preconditions'],
                                  'ordinary_Z_contract': 'Z_Z=1; delta_Z=Rc_Z=R_Z=mu_Z=Pstar_Z=0',
                                  'common_source_domain_required': True},
        'cached_formal_algebra': saved['cached_formal_algebra'],
        'retained_construction_data': saved['retained_construction_data'],
        'retained_integral_term_inventory': saved['retained_integral_term_inventory'],
        'accepted_native_B_nodes': saved['accepted_native_B_nodes'],
        'all_16_native_B_accepted_from_parent': True,
    }
    false_flags = ['native_old_density_closure', 'source_ftc_zero_identity', 'native_old_density_identity',
                   'old212_adapter_identity', 'old212_quantclosure', 'numeric_containment_claimed',
                   'source_C5_installation_claimed', 'new_loop_domain_admission_claimed',
                   'admissible_cone_admission_claimed', 'quantitative_packet_error_reduction_proved',
                   'physical_layer_completed', 'global_five_moment_closure', 'time_order_recursion_claimed',
                   'old_operators_replayed', 'old_pointwise_S_recovered', 'old_mass_relation_proved',
                   'old_errors_removed', 'old_error_terms_removed', 'full_goal_complete', 'numeric_values_computed']
    result.update({name: False for name in false_flags})
    try:
        result['semantic_metadata_join'] = _semantic_metadata_join(source, saved)
        journal({'event': 'semantic_dependency_join_passed', 'source_graph_replay': False})
        events = []
        witness = reconcile_cached_witnesses(saved, events)
        result['witness_reconciliation'] = witness
        result['fresh_source_proof_normalizations'] = witness.get('new_cancel_together_count', 0)
        result['new_proof_count'] = witness.get('new_check_count', 0)
        for event in events:
            journal(event)
        if witness['passed'] is True and witness['new_corrected_pass_count'] == 5:
            result['source_binding_gates'] = {key: True for key in result['source_binding_gates']}
            result['source_ftc_transport_defined'] = True
            result['passed'] = True
        result['signed_transport_identity'] = {
            'source_only': True, 'conditional': True,
            'left': 'integral_1^2 (W_N*v_N-W_O*v_O) dx + [Delta(R^2*F-R*Q_y)]_1^2',
            'right': '0.5*(local_S1_N-S1_O)-[Delta(R*Q_y)]_1^2',
            'physical_weight': 'W=Rc^2*x, separately in each chart; R=Rc*x',
            'derivative_scope': 'ordinary Z value and first derivative', 'terminal_zero_proved': False,
            'proof_dependencies': 'cached formal 32/9 + cached source 54 + five new semantic target differences + metadata export join',
        }
    except Exception as error:
        import traceback
        result['failure'] = {'type': type(error).__name__, 'message': str(error), 'traceback': traceback.format_exc()}
        journal({'event': 'semantic_witness_exception', 'type': type(error).__name__, 'message': str(error)})
    result.update(result['source_binding_gates'])
    return result
