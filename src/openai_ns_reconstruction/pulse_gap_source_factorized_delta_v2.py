"""Conditional, source-only factorized P-minus-F pulse-gap packet.

This leaf is intentionally a boundary consumer rather than a numerical
adapter.  It authenticates every saved input on raw bytes, then reads only
the named direct fields needed by a later source-expression adapter:

* the published ``pulse_gap_actual_source_program_v1`` result;
* the actual physical pulse report;
* the saved original frame archive; and
* the saved exact n1 pulse-kernel archive and its check receipt.

The two gzip archives are not decompressed until *all* five raw SHA-256
checks have passed.  No old Python source is imported, no graph edge or
closure is followed, and no saved node id is used as a value.  The returned
trees are unevaluated expressions with opaque source references, namespace
labels, and explicit unresolved gates.

The formulas are conditional.  In particular, direct ``u``/``C``/kernel
references do not prove common native ownership, derivatives, source
semantics, ``mu``, geometry, or the historical F field radius.  The radius
factorization is retained even when a caller later supplies normalized
sources:

``R_P^2 DeltaS + (R_P-R_F)(R_P+R_F) S_F``

and

``R_P DeltaQy + (R_P-R_F) QyF``.

There is deliberately no ``__main__`` entry point.  The one-shot runner is
the only intended invocation boundary after the root agent freezes a
checkpoint and pending brief.
"""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Sequence


SCHEMA = "pulse_gap_source_factorized_delta_v2"
COMPONENT = "pulse_gap_source_factorized_delta_v2"

ACTUAL_RESULT = Path("tmp/pulse_gap_actual_source_program_v1/result.json")
PHYSICAL_REPORT = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_native_physical_pulse_source.json"
)
FRAME_ARCHIVE = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_Rp_native_C5_pulse_source.json.gz"
)
KERNEL_ARCHIVE = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_pulse_kernel_functions.json.gz"
)
KERNEL_CHECK = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_pulse_kernel_functions_check.json"
)

EXPECTED_RAW_SHA256 = {
    "published_actual_result": "e611cdc7e049ce40e6667275b4e9706c3478044880b92c93baebf8ab8d06ebf2",
    "physical_report": "8b21f8e868025d8cc95149ffab80b54603833f5f299ec9d100602af355cea82b",
    "frame_archive": "442df70e5ab4051066b241e943c8d0531067d057e774a04dcc341dadbd703ebe",
    "kernel_archive": "9584201b517fe3058ea6d1a169ddb7d319776eb97479b36024af4bffd4973b75",
    "kernel_check": "47343a839590c3cb08c8d726ef225415e15f4fb12a5f7cf888049ea1f20cffc0",
}

P_OWNER = "P"
F_OWNER = "F"
K_OWNER = "K"
P_NAMESPACE = "physical_current_native_report"
F_NAMESPACE = "frame_Rp_native_C5_archive"
K_NAMESPACE = "exact_n1_pulse_kernel_archive"

# Four rows are needed for the conditional m0--m3 rule.  We retain every
# saved row as well, without treating row count or ids as values.
MOMENT_ORDERS = (0, 1, 2, 3)
Z_ORDERS = (0, 1, 2)


class SourceFactorizedDeltaInputError(ValueError):
    """Raised only for an impossible direct-field access inside this leaf."""


def _sha256(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def _rel(path: Path) -> str:
    return path.as_posix()


def _expr(operator: str, *arguments: Any, **metadata: Any) -> dict[str, Any]:
    """Create a small named expression tree; never evaluate an argument."""

    value: dict[str, Any] = {
        "operator": operator,
        "arguments": list(arguments),
        "source_values_selected": False,
    }
    value.update(metadata)
    return value


def _literal(value: str) -> dict[str, Any]:
    return {"kind": "exact_scalar_literal", "value": value, "source_values_selected": False}


def _symbol(name: str, *, normalized: bool = False) -> dict[str, Any]:
    return {
        "kind": "formal_symbol",
        "name": name,
        "normalized_source": normalized,
        "source_values_selected": False,
    }


def _ref(name: str, *, normalized: bool = False, ordinary_Z_order: int | None = None) -> dict[str, Any]:
    """Reference a named DAG node without copying its expression tree."""

    value: dict[str, Any] = {
        "kind": "named_expression_ref",
        "name": name,
        "normalized_source": normalized,
        "source_values_selected": False,
    }
    if ordinary_Z_order is not None:
        value["ordinary_Z_order"] = ordinary_Z_order
    return value


def _named_definition(name: str, expression: Mapping[str, Any], **metadata: Any) -> dict[str, Any]:
    """Store one expression once; consumers use ``_ref(name)``."""

    value: dict[str, Any] = {
        "kind": "named_expression_definition",
        "name": name,
        "expression": dict(expression),
        "source_values_selected": False,
    }
    value.update(metadata)
    return value


def _direct(value: Any, path: Sequence[str], label: str) -> Any:
    """Read a literal mapping path only; never interpret graph nodes."""

    current = value
    for key in path:
        if not isinstance(current, Mapping) or key not in current:
            raise SourceFactorizedDeltaInputError(f"missing direct source field {label}")
        current = current[key]
    return current


def _optional(value: Any, path: Sequence[str]) -> Any | None:
    current = value
    for key in path:
        if not isinstance(current, Mapping) or key not in current:
            return None
        current = current[key]
    return current


def _raw_packet(
    root: Path,
    relative: Path,
    expected: str,
    label: str,
    *,
    gzip_encoded: bool,
) -> tuple[bytes | None, dict[str, Any], str | None]:
    """Read and authenticate raw bytes, but do not parse or decompress them."""

    path = root / relative
    meta: dict[str, Any] = {
        "path": _rel(relative),
        "expected_raw_sha256": expected,
        "gzip_encoded": gzip_encoded,
        "raw_hash_authenticated_before_parse": False,
    }
    if not path.is_file():
        meta["missing"] = True
        return None, meta, f"{label}:missing"
    raw = path.read_bytes()
    actual = _sha256(raw)
    meta["raw_sha256"] = actual
    meta["raw_hash_authenticated_before_parse"] = actual == expected
    if actual != expected:
        return None, meta, f"{label}:raw_sha256_changed"
    return raw, meta, None


def _parse_json(raw: bytes, label: str) -> Mapping[str, Any]:
    try:
        value = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise SourceFactorizedDeltaInputError(f"{label}:invalid_json") from exc
    if not isinstance(value, Mapping):
        raise SourceFactorizedDeltaInputError(f"{label}:top_level_not_object")
    return value


def _decode_after_all_raw_auth(raw: bytes, label: str) -> Mapping[str, Any]:
    """Decompress only after the caller authenticated every raw input."""

    try:
        decoded = gzip.decompress(raw)
        value = json.loads(decoded.decode("utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise SourceFactorizedDeltaInputError(f"{label}:invalid_gzip_json") from exc
    if not isinstance(value, Mapping):
        raise SourceFactorizedDeltaInputError(f"{label}:decoded_top_level_not_object")
    return value


def _opaque_ref(
    value: Any,
    pointer: str,
    owner: str,
    namespace: str,
    gates: list[str],
    label: str,
    *,
    ordinary_Z_order: int | None = None,
    role_kind: str = "direct_source_ref",
) -> dict[str, Any]:
    """Wrap a saved id as an opaque reference and never as a numeric value."""

    if not isinstance(value, int) or isinstance(value, bool):
        gates.append(f"{label}:direct_source_node_id_missing_or_noninteger")
        result: dict[str, Any] = {
            "status": "unresolved",
            "owner": owner,
            "namespace": namespace,
            "source_role": pointer,
            "node_id_is_not_a_numeric_value": True,
            "source_value_selected": False,
        }
    else:
        result = {
            "status": role_kind,
            "owner": owner,
            "namespace": namespace,
            "source_role": pointer,
            "source_node_id": value,
            "node_id_is_not_a_numeric_value": True,
            "source_value_selected": False,
        }
    if ordinary_Z_order is not None:
        result["ordinary_Z_order"] = ordinary_Z_order
    return result


def _rows(
    value: Any,
    pointer: str,
    owner: str,
    namespace: str,
    gates: list[str],
    label: str,
    *,
    minimum: int,
) -> dict[str, Any]:
    """Keep direct row references, with explicit ordinary-Z shape gates."""

    if not isinstance(value, list):
        gates.append(f"{label}:direct_rows_missing_or_not_list")
        return {
            "status": "unresolved",
            "owner": owner,
            "namespace": namespace,
            "source_role": pointer,
            "ordinary_Z_orders": list(range(minimum)),
            "all_saved_rows": [],
            "source_values_selected": False,
        }
    if len(value) < minimum:
        gates.append(f"{label}:fewer_than_{minimum}_saved_rows")
    refs = [
        _opaque_ref(
            row,
            f"{pointer}[{index}]",
            owner,
            namespace,
            gates,
            f"{label}[{index}]",
            ordinary_Z_order=index,
        )
        for index, row in enumerate(value)
    ]
    return {
        "status": "direct_source_rows" if len(value) >= minimum else "unresolved",
        "owner": owner,
        "namespace": namespace,
        "source_role": pointer,
        "ordinary_Z_orders": list(range(minimum)),
        "ordinary_Z0_to_required": refs[:minimum],
        "ordinary_Z0_to_Z3": refs[:4] if len(value) >= 4 else [],
        "all_saved_rows": refs,
        "saved_row_count": len(value),
        "source_values_selected": False,
    }


def _row_at(rows: Mapping[str, Any], order: int, label: str, gates: list[str]) -> dict[str, Any]:
    saved = rows.get("all_saved_rows")
    if not isinstance(saved, list) or order >= len(saved):
        gates.append(f"{label}:row_{order}_unavailable")
        return {
            "status": "unresolved",
            "source_role": f"{rows.get('source_role', label)}[{order}]",
            "source_value_selected": False,
        }
    return dict(saved[order])


def _source_family(document: Mapping[str, Any], label: str, gates: list[str]) -> dict[str, Any]:
    family = document.get("source_family")
    if not isinstance(family, Mapping) or not family:
        gates.append(f"{label}:source_family_missing")
        return {
            "status": "unresolved_missing_saved_source_family",
            "source_value_selected": False,
        }
    return {
        "status": "saved_source_family_metadata",
        "fields": dict(family),
        "source_role": f"{label}.source_family",
        "source_value_selected": False,
    }


def _ordinary_Z_jet(source: Mapping[str, Any], order: int) -> dict[str, Any]:
    return {
        "kind": "ordinary_Z_jet_source_handle",
        "source": dict(source),
        "ordinary_Z_order": order,
        "derivative_semantics_proved": False,
        "source_value_selected": False,
    }


def _jet_family(source: Mapping[str, Any], name: str) -> dict[str, Any]:
    return {
        f"{name}_Z{order}": _ordinary_Z_jet(source, order)
        for order in Z_ORDERS
    }


def _first_rows(
    binding: Mapping[str, Any],
    pointer: str,
    owner: str,
    namespace: str,
    gates: list[str],
    label: str,
) -> dict[str, Any]:
    value = binding.get("first_row")
    if not isinstance(value, list):
        gates.append(f"{label}:direct_first_row_missing_or_not_list")
        return {
            "status": "unresolved",
            "owner": owner,
            "namespace": namespace,
            "source_role": pointer,
            "all_saved_rows": [],
            "source_values_selected": False,
        }
    if len(value) < 2:
        gates.append(f"{label}:fewer_than_2_first_row_weights")
    refs = []
    for index, row in enumerate(value):
        ref = _opaque_ref(
            row,
            f"{pointer}[{index}]",
            owner,
            namespace,
            gates,
            f"{label}[{index}]",
            role_kind="first_row_weight_ref",
        )
        ref["index"] = index
        ref["first_row_weight"] = True
        refs.append(ref)
    return {
        "status": "direct_first_row" if len(value) >= 2 else "unresolved",
        "owner": owner,
        "namespace": namespace,
        "source_role": pointer,
        "first_row_weight_refs": refs,
        "all_saved_rows": refs,
        "saved_row_count": len(value),
        "source_values_selected": False,
    }


def _side_snapshot(
    binding: Mapping[str, Any],
    *,
    owner: str,
    namespace: str,
    binding_pointer: str,
    coordinate_key: str,
    gates: list[str],
) -> dict[str, Any]:
    current_frame = binding.get("current_frame")
    if not isinstance(current_frame, Mapping):
        gates.append(f"{owner}:current_frame_missing")
        current_frame = {}
    u_rows = _rows(
        current_frame.get("u"),
        f"{binding_pointer}.current_frame.u",
        owner,
        namespace,
        gates,
        f"{owner}.current_frame.u",
        minimum=4,
    )

    selected = binding.get("selected_pulse")
    if not isinstance(selected, Mapping):
        gates.append(f"{owner}:selected_pulse_missing")
        selected = {}
    coefficients = {
        name: _rows(
            selected.get(name),
            f"{binding_pointer}.selected_pulse.{name}",
            owner,
            namespace,
            gates,
            f"{owner}.selected_pulse.{name}",
            minimum=4,
        )
        for name in ("C1", "C2")
    }

    coordinates = binding.get(coordinate_key)
    if not isinstance(coordinates, Mapping):
        gates.append(f"{owner}:{coordinate_key}_missing")
        coordinates = {}
    gap = coordinates.get("gap")
    if not isinstance(gap, Mapping):
        gates.append(f"{owner}:{coordinate_key}.gap_missing")
        gap = {}
    radius_value = gap.get("radius")
    radius_ref = _opaque_ref(
        radius_value,
        f"{binding_pointer}.{coordinate_key}.gap.radius",
        owner,
        namespace,
        gates,
        f"{owner}.gap.radius",
        role_kind=(
            "conditional_field_radius_ref"
            if owner == P_OWNER
            else "conditional_current_radius_metadata_ref"
        ),
    )
    if owner == F_OWNER:
        gates.append("F:historical_field_radius_owner_unproved")
        radius_ref = {
            "status": "unresolved", "owner": F_OWNER,
            "source_role": "historical frame-only gap field radius",
            "current_radius_metadata_not_used_as_field_radius": radius_ref,
            "field_radius_identity_proved": False, "source_value_selected": False,
        }
    else:
        radius_ref["field_radius_identity_proved"] = False

    return {
        "owner": owner,
        "namespace": namespace,
        "binding_source_role": binding_pointer,
        "current_frame": {
            "u": u_rows,
            "required_ordinary_Z_orders": list(MOMENT_ORDERS),
            "four_rows_checked": u_rows.get("status") == "direct_source_rows",
            "source_values_selected": False,
        },
        "selected_pulse": coefficients,
        "radius": radius_ref,
        "source_values_selected": False,
        "source_graph_traversed": False,
    }


def _shared_u(p: Mapping[str, Any], f: Mapping[str, Any], order: int) -> dict[str, Any]:
    return {
        "kind": "conditional_shared_native_u_handle",
        "name": f"u_{order}_shared_native",
        "ordinary_Z_order": order,
        "side_refs": {"P": dict(p), "F": dict(f)},
        "same_native_owner_proved": False,
        "source_value_selected": False,
    }


def _build_q_dag(moment_rows: Mapping[str, Any]) -> dict[str, Any]:
    """Build the conditional q recursion from mF and delta_m by name only."""

    definitions: dict[str, dict[str, Any]] = {}
    symbols = {
        "delta": _symbol("delta"),
        "Z": _symbol("Z"),
        "T": {
            "kind": "conditional_common_field_scale",
            "name": "T",
            "definition": "exp(common original gap field logscale), including logPstar exactly once",
            "P_F_common_identity_proved": False,
            "scale_Z_independence_proved": False,
            "source_values_selected": False,
        },
    }

    def define(name: str, expression: Mapping[str, Any], **metadata: Any) -> None:
        definitions[name] = _named_definition(name, expression, **metadata)

    define(
        "one_minus_delta",
        _expr("subtract", _literal("1"), _ref("delta"), name="one_minus_delta"),
    )
    define(
        "one_plus_delta",
        _expr("add", _literal("1"), _ref("delta"), name="one_plus_delta"),
    )
    define(
        "one_minus_Z2",
        _expr("subtract", _literal("1"), _expr("power", _ref("Z"), _literal("2")), name="one_minus_Z2"),
    )
    define(
        "three_plus_delta",
        _expr("add", _literal("3"), _ref("delta"), name="three_plus_delta"),
    )
    define(
        "two_delta",
        _expr("multiply", _literal("2"), _ref("delta"), name="two_delta"),
    )
    define(
        "L",
        _expr(
            "subtract",
            _literal("1"),
            _expr("multiply", _ref("delta"), _expr("power", _ref("Z"), _literal("2"))),
            name="L",
            formula_text="L=1-delta*Z^2",
        ),
    )
    define(
        "L1",
        _expr("multiply", _literal("-2"), _ref("delta"), _ref("Z"), name="L1", formula_text="L1=-2*delta*Z"),
    )
    define(
        "L2",
        _expr("multiply", _literal("-2"), _ref("delta"), name="L2", formula_text="L2=-2*delta"),
    )

    def m_ref(prefix: str, order: int) -> dict[str, Any]:
        return _ref(f"{prefix}_m{order}", normalized=True, ordinary_Z_order=order)

    def define_D(prefix: str) -> None:
        d0 = _expr(
            "add",
            _expr("multiply", _ref("one_minus_delta"), _ref("Z"), m_ref(prefix, 0)),
            _expr("multiply", _ref("one_minus_Z2"), m_ref(prefix, 1)),
            name=f"D0_{prefix}",
            formula_text="D(m)=(1-delta)*Z*m0+(1-Z^2)*m1",
        )
        d1 = _expr(
            "add",
            _expr("multiply", _ref("one_minus_delta"), m_ref(prefix, 0)),
            _expr("negate", _expr("multiply", _ref("one_plus_delta"), _ref("Z"), m_ref(prefix, 1))),
            _expr("multiply", _ref("one_minus_Z2"), m_ref(prefix, 2)),
            name=f"D1_{prefix}",
            formula_text="D1=(1-delta)*m0-(1+delta)*Z*m1+(1-Z^2)*m2",
        )
        d2 = _expr(
            "add",
            _expr("negate", _expr("multiply", _ref("two_delta"), m_ref(prefix, 1))),
            _expr("negate", _expr("multiply", _ref("three_plus_delta"), _ref("Z"), m_ref(prefix, 2))),
            _expr("multiply", _ref("one_minus_Z2"), m_ref(prefix, 3)),
            name=f"D2_{prefix}",
            formula_text="D2=-2*delta*m1-(3+delta)*Z*m2+(1-Z^2)*m3",
        )
        define(f"D0_{prefix}", d0, ordinary_Z_order=0)
        define(f"D1_{prefix}", d1, ordinary_Z_order=1)
        define(f"D2_{prefix}", d2, ordinary_Z_order=2)

    def define_q_family(prefix: str, q_names: tuple[str, str, str]) -> None:
        define_D(prefix)
        q0, q1, q2 = q_names
        define(
            q0,
            _expr("divide", _expr("negate", _ref(f"D0_{prefix}")), _ref("L"), name=q0),
            ordinary_Z_order=0,
            normalized_source=True,
        )
        define(
            q1,
            _expr(
                "divide",
                _expr(
                    "add",
                    _expr("negate", _ref(f"D1_{prefix}")),
                    _expr("negate", _expr("multiply", _ref("L1"), _ref(q0))),
                ),
                _ref("L"),
                name=q1,
            ),
            ordinary_Z_order=1,
            normalized_source=True,
        )
        define(
            q2,
            _expr(
                "divide",
                _expr(
                    "add",
                    _expr("negate", _ref(f"D2_{prefix}")),
                    _expr("negate", _expr("multiply", _literal("2"), _ref("L1"), _ref(q1))),
                    _expr("negate", _expr("multiply", _ref("L2"), _ref(q0))),
                ),
                _ref("L"),
                name=q2,
            ),
            ordinary_Z_order=2,
            normalized_source=True,
        )

    define_q_family("mF", ("qF0", "qF1", "qF2"))
    define_q_family("delta_m", ("dq0", "dq1", "dq2"))

    define(
        "slinear_F",
        _expr(
            "divide",
            _expr("multiply", _ref("one_minus_delta"), _ref("Z"), _ref("qF1")),
            _expr("multiply", _literal("2"), _ref("L")),
            name="slinear_F",
            formula_text="slinear_F=(1-delta)*Z*qF1/(2*L)",
        ),
        normalized_source=True,
    )
    define(
        "delta_slinear",
        _expr(
            "divide",
            _expr("multiply", _ref("one_minus_delta"), _ref("Z"), _ref("dq1")),
            _expr("multiply", _literal("2"), _ref("L")),
            name="delta_slinear",
            formula_text="delta_slinear=(1-delta)*Z*dq1/(2*L)",
        ),
        normalized_source=True,
    )
    define(
        "slinear_F_Z1",
        _expr(
            "divide",
            _expr(
                "multiply",
                _ref("one_minus_delta"),
                _expr(
                    "add",
                    _ref("qF1"),
                    _expr("multiply", _ref("Z"), _ref("qF2")),
                    _expr("divide", _expr("multiply", _literal("2"), _ref("delta"), _expr("power", _ref("Z"), _literal("2")), _ref("qF1")), _ref("L")),
                ),
            ),
            _expr("multiply", _literal("2"), _ref("L")),
            name="slinear_F_Z1",
            ordinary_Z_order=1,
        ),
        normalized_source=True,
        conditions=["T_Z=0"],
    )
    define(
        "delta_slinear_Z1",
        _expr(
            "divide",
            _expr(
                "multiply",
                _ref("one_minus_delta"),
                _expr(
                    "add",
                    _ref("dq1"),
                    _expr("multiply", _ref("Z"), _ref("dq2")),
                    _expr("divide", _expr("multiply", _literal("2"), _ref("delta"), _expr("power", _ref("Z"), _literal("2")), _ref("dq1")), _ref("L")),
                ),
            ),
            _expr("multiply", _literal("2"), _ref("L")),
            name="delta_slinear_Z1",
            formula_text="(1-delta)/(2*L)*(dq1+Z*dq2+2*delta*Z^2*dq1/L)",
            ordinary_Z_order=1,
        ),
        normalized_source=True,
        conditions=["T_Z=0"],
    )
    define(
        "squad_F",
        _expr("negate", _expr("divide", _expr("power", _ref("qF0"), _literal("2")), _literal("2")), name="squad_F"),
        normalized_source=True,
    )
    define(
        "delta_squad",
        _expr(
            "negate",
            _expr("multiply", _expr("add", _ref("qF0"), _expr("divide", _ref("dq0"), _literal("2"))), _ref("dq0")),
            name="delta_squad",
            formula_text="delta_squad=-(qF0+dq0/2)*dq0",
        ),
        normalized_source=True,
    )
    define(
        "squad_F_Z1",
        _expr("negate", _expr("multiply", _ref("qF0"), _ref("qF1")), name="squad_F_Z1"),
        normalized_source=True,
        ordinary_Z_order=1,
        conditions=["T_Z=0"],
    )
    define(
        "delta_squad_Z1",
        _expr(
            "negate",
            _expr(
                "add",
                _expr("multiply", _expr("add", _ref("qF1"), _expr("divide", _ref("dq1"), _literal("2"))), _ref("dq0")),
                _expr("multiply", _expr("add", _ref("qF0"), _expr("divide", _ref("dq0"), _literal("2"))), _ref("dq1")),
            ),
            name="delta_squad_Z1",
        ),
        normalized_source=True,
        ordinary_Z_order=1,
        conditions=["T_Z=0"],
    )
    define(
        "S_F",
        _expr(
            "add",
            _expr("multiply", _ref("T"), _ref("slinear_F")),
            _expr("multiply", _expr("power", _ref("T"), _literal("2")), _ref("squad_F")),
            name="S_F",
            formula_text="S_F=T*slinear_F+T^2*squad_F",
        ),
        conditions=["T_Z=0"],
    )
    define(
        "DeltaS",
        _expr(
            "add",
            _expr("multiply", _ref("T"), _ref("delta_slinear")),
            _expr("multiply", _expr("power", _ref("T"), _literal("2")), _ref("delta_squad")),
            name="DeltaS",
            formula_text="DeltaS=T*delta_slinear+T^2*delta_squad",
        ),
        conditions=["T_Z=0"],
    )
    define(
        "S_F_Z1",
        _expr(
            "add",
            _expr("multiply", _ref("T"), _ref("slinear_F_Z1")),
            _expr("multiply", _expr("power", _ref("T"), _literal("2")), _ref("squad_F_Z1")),
            name="S_F_Z1",
            ordinary_Z_order=1,
        ),
        conditions=["T_Z=0"],
    )
    define(
        "DeltaS_Z1",
        _expr(
            "add",
            _expr("multiply", _ref("T"), _ref("delta_slinear_Z1")),
            _expr("multiply", _expr("power", _ref("T"), _literal("2")), _ref("delta_squad_Z1")),
            name="DeltaS_Z1",
            ordinary_Z_order=1,
        ),
        conditions=["T_Z=0"],
    )
    define(
        "QyF",
        _expr("negate", _expr("multiply", _ref("T"), _ref("qF0")), name="QyF", formula_text="QyF=-T*qF0"),
    )
    define(
        "DeltaQy",
        _expr("negate", _expr("multiply", _ref("T"), _ref("dq0")), name="DeltaQy", formula_text="DeltaQy=-T*dq0"),
    )
    define("QyF_Z1", _expr("negate", _expr("multiply", _ref("T"), _ref("qF1")), name="QyF_Z1"), conditions=["T_Z=0"])
    define("DeltaQy_Z1", _expr("negate", _expr("multiply", _ref("T"), _ref("dq1")), name="DeltaQy_Z1"), conditions=["T_Z=0"])

    normalized = {
        "qF": {f"Z{order}": _ref(f"qF{order}", normalized=True, ordinary_Z_order=order) for order in Z_ORDERS},
        "delta_q": {f"Z{order}": _ref(f"dq{order}", normalized=True, ordinary_Z_order=order) for order in Z_ORDERS},
        "normalization_contract": {
            "qF_is_normalized_source": True,
            "qDelta_is_normalized_source": True,
            "delta_q_definition": "dq0,dq1,dq2 from delta_m recursion (not qP-qF)",
            "source_semantics_proved": False,
            "source_values_selected": False,
        },
    }
    return {
        "symbols": symbols,
        "definitions": definitions,
        "normalized_sources": normalized,
        "source_values_selected": False,
    }


def _binomial(j: int, k: int) -> str:
    # The four rows are intentionally literal metadata, not evaluated source
    # data and not arithmetic on any saved node id.
    table = (
        (1,),
        (1, 1),
        (1, 2, 1),
        (1, 3, 3, 1),
    )
    return str(table[j][k])


def _moment_rows(
    p_u: Mapping[str, Any],
    f_u: Mapping[str, Any],
    cF_jets: Mapping[str, Any],
    delta_c_jets: Mapping[str, Any],
    gates: list[str],
) -> dict[str, Any]:
    def family(prefix: str, source_u: Mapping[str, Any], *, shared: bool) -> dict[str, Any]:
        rows: dict[str, Any] = {}
        for j in MOMENT_ORDERS:
            terms = []
            for k in range(j + 1):
                if shared:
                    u_k = _shared_u(
                        _row_at(p_u, k, "P.current_frame.u", gates),
                        _row_at(f_u, k, "F.current_frame.u", gates),
                        k,
                    )
                else:
                    u_k = _row_at(source_u, k, "F.current_frame.u", gates)
                c_name = f"{'cF' if prefix == 'mF' else 'delta_c'}_Z{j-k}"
                terms.append(
                    _expr(
                        "multiply",
                        _literal(_binomial(j, k)),
                        u_k,
                        _ref(c_name, ordinary_Z_order=j-k),
                    )
                )
            row_name = f"{prefix}_m{j}"
            rows[row_name] = {
                "ordinary_Z_order": j,
                "formula": _expr("negate", _expr("add", *terms), name=row_name),
                "formula_text": f"{row_name}=-sum(k=0..{j}) binom(j,k) u_k {prefix}_c_(j-k)",
                "conditional_on": [
                    "same_native_u_owner" if shared else "F_native_u_owner",
                    "same_first_row_kernel_owner",
                    "first_row_weights_Z_independent",
                    "source_semantics_and_normalization",
                ],
                "normalized_source": True,
                "source_values_selected": False,
            }
        return rows

    return {
        "mF": family("mF", f_u, shared=False),
        "delta_m": family("delta_m", f_u, shared=True),
        "ordinary_Z_orders": list(MOMENT_ORDERS),
        "source_values_selected": False,
    }


def _radius_factorized_expressions(
    p_radius: Mapping[str, Any],
    f_radius: Mapping[str, Any],
) -> dict[str, Any]:
    rp, rf = dict(p_radius), dict(f_radius)
    difference = _expr("subtract", rp, rf, name="R_P_minus_R_F")
    total = _expr("add", rp, rf, name="R_P_plus_R_F")
    return {
        "R_P": rp,
        "R_F": rf,
        "R_P_minus_R_F": difference,
        "R_P_plus_R_F": total,
        "source_values_selected": False,
    }


def _z_jets(source: Mapping[str, Any], name: str) -> dict[str, Any]:
    return _jet_family(source, name)


def _factorized_Z1_Z2(
    radius: Mapping[str, Any],
    dag: Mapping[str, Any],
) -> dict[str, Any]:
    """Expose named q/S/radius nodes without expanding a source expression tree."""

    local: dict[str, dict[str, Any]] = {}

    def define(name: str, expression: Mapping[str, Any], **metadata: Any) -> None:
        local[name] = _named_definition(name, expression, **metadata)

    rp = radius["R_P"]
    rf = radius["R_F"]
    dr = radius["R_P_minus_R_F"]
    rs = radius["R_P_plus_R_F"]
    rp_j = _z_jets(rp, "R_P")
    dr_j = _z_jets(dr, "R_P_minus_R_F")
    rs_j = _z_jets(rs, "R_P_plus_R_F")

    define("R_P_squared", _expr("power", rp, _literal("2"), name="R_P_squared"))
    define(
        "R_P_squared_Z1",
        _expr("multiply", _literal("2"), rp, rp_j["R_P_Z1"], name="R_P_squared_Z1"),
        ordinary_Z_order=1,
    )
    define(
        "radius_product",
        _expr("multiply", dr, rs, name="(R_P_minus_R_F)(R_P_plus_R_F)"),
    )
    define(
        "radius_product_Z1",
        _expr(
            "add",
            _expr("multiply", dr_j["R_P_minus_R_F_Z1"], rs),
            _expr("multiply", dr, rs_j["R_P_plus_R_F_Z1"]),
            name="radius_product_Z1",
            ordinary_Z_order=1,
        ),
    )
    define(
        "physical_radius_DeltaS",
        _expr(
            "add",
            _expr("multiply", _ref("R_P_squared"), _ref("DeltaS")),
            _expr("multiply", _ref("radius_product"), _ref("S_F")),
            name="physical_radius_DeltaS",
            formula_text="R_P^2 DeltaS +(R_P-R_F)(R_P+R_F) S_F",
        ),
        conditions=["T_Z=0"],
    )
    define(
        "physical_radius_DeltaS_Z1",
        _expr(
            "add",
            _expr("multiply", _ref("R_P_squared_Z1"), _ref("DeltaS")),
            _expr("multiply", _ref("R_P_squared"), _ref("DeltaS_Z1")),
            _expr("multiply", _ref("radius_product_Z1"), _ref("S_F")),
            _expr("multiply", _ref("radius_product"), _ref("S_F_Z1")),
            name="physical_radius_DeltaS_Z1",
            formula_text="Z1[R_P^2 DeltaS +(R_P-R_F)(R_P+R_F) S_F]",
            ordinary_Z_order=1,
        ),
        conditions=["T_Z=0", "R_P_Z1 and R_F_Z1 are conditional source jets"],
    )
    define(
        "physical_endpoint_DeltaQy",
        _expr(
            "add",
            _expr("multiply", rp, _ref("DeltaQy")),
            _expr("multiply", dr, _ref("QyF")),
            name="physical_endpoint_DeltaQy",
            formula_text="R_P DeltaQy +(R_P-R_F) QyF",
        ),
    )
    define(
        "physical_endpoint_DeltaQy_Z1",
        _expr(
            "add",
            _expr("multiply", rp_j["R_P_Z1"], _ref("DeltaQy")),
            _expr("multiply", rp, _ref("DeltaQy_Z1")),
            _expr("multiply", dr_j["R_P_minus_R_F_Z1"], _ref("QyF")),
            _expr("multiply", dr, _ref("QyF_Z1")),
            name="physical_endpoint_DeltaQy_Z1",
            formula_text="Z1[R_P DeltaQy +(R_P-R_F) QyF]",
            ordinary_Z_order=1,
        ),
        conditions=["T_Z=0", "R_P_Z1 and R_F_Z1 are conditional source jets"],
    )

    delta_Q = {
        f"Z{order}": _ref(f"dq{order}", normalized=True, ordinary_Z_order=order)
        for order in Z_ORDERS
    }
    source_Qy = {
        "Z0": _ref("DeltaQy"),
        "Z1": _ref("DeltaQy_Z1", ordinary_Z_order=1),
    }
    return {
        "expression_dag": {
            "definitions": local,
            "source_values_selected": False,
        },
        "delta_Q": {
            **delta_Q,
            "formula_text": "delta_Q=delta_q (normalized; dq0,dq1,dq2)",
            "source_values_selected": False,
        },
        "delta_Qy_source": source_Qy,
        "S_linear": {
            "F": _ref("slinear_F"),
            "delta": _ref("delta_slinear"),
            "source_values_selected": False,
        },
        "S_quadratic": {
            "F": _ref("squad_F"),
            "delta": _ref("delta_squad"),
            "source_values_selected": False,
        },
        "S_linear_Z1": {
            "expression": _ref("delta_slinear_Z1", ordinary_Z_order=1),
            "formula_text": "(1-delta)/(2*L)*(dq1+Z*dq2+2*delta*Z^2*dq1/L)",
            "conditions": ["T_Z=0"],
            "source_values_selected": False,
        },
        "S_quadratic_Z1": {
            "expression": _ref("delta_squad_Z1", ordinary_Z_order=1),
            "formula_text": "Z1[-(qF0+dq0/2)*dq0]",
            "conditions": ["T_Z=0"],
            "source_values_selected": False,
        },
        "S_F": {"expression": _ref("S_F"), "conditions": ["T_Z=0"], "source_values_selected": False},
        "DeltaS": {"expression": _ref("DeltaS"), "conditions": ["T_Z=0"], "source_values_selected": False},
        "S_F_Z1": {"expression": _ref("S_F_Z1", ordinary_Z_order=1), "conditions": ["T_Z=0"], "source_values_selected": False},
        "DeltaS_Z1": {"expression": _ref("DeltaS_Z1", ordinary_Z_order=1), "conditions": ["T_Z=0"], "source_values_selected": False},
        "actual_physical_radius_difference_endpoint_terms": {
            "linear_and_quadratic_density": {
                "Z0": _ref("physical_radius_DeltaS"),
                "Z1": _ref("physical_radius_DeltaS_Z1", ordinary_Z_order=1),
                "formula_text": "R_P^2 DeltaS +(R_P-R_F)(R_P+R_F) S_F",
                "m15_integration": "multiply by 1/2 and directed dxi/mu; measure Z-independence is unproved",
                "conditions": ["T_Z=0"],
                "source_values_selected": False,
            },
            "linear_Qy_endpoint": {
                "Z0": _ref("physical_endpoint_DeltaQy"),
                "Z1": _ref("physical_endpoint_DeltaQy_Z1", ordinary_Z_order=1),
                "formula_text": "R_P DeltaQy +(R_P-R_F) QyF",
                "m15_orientation": "minus (upper minus lower)",
                "conditions": ["T_Z=0"],
                "source_values_selected": False,
            },
        },
        "derivative_gates": {
            "coverage": "normalized Q through Z2; physical density and endpoint through Z1 only",
            "radius_Z_independence_proved": False,
            "normalized_source_derivative_semantics_proved": False,
            "ordinary_Z1_transport_proved": False,
            "ordinary_Z2_transport_proved": False,
            "T_Z_zero_assumption_proved": False,
        },
    }


def _provenance_meta(meta: Mapping[str, Any]) -> dict[str, Any]:
    return dict(meta)


def _failed_packet(provenance: Mapping[str, Any], gates: Sequence[str], diagnostic: str) -> dict[str, Any]:
    return {
        "schema": SCHEMA,
        "status": "failed_closed_source_boundary",
        "source_only": True,
        "admissible_for_numeric_use": False,
        "formal_source_expression_defined": False,
        "numeric_values_computed": False,
        "numeric_difference_evaluated": False,
        "source_values_selected": False,
        "diagnostic": diagnostic,
        "unresolved_gates": sorted(set(gates)),
        "provenance": dict(provenance),
        "scope_audit": {
            "gzip_archives_decoded": 0,
            "gzip_raw_hashes_authenticated_before_decode": True,
            "source_graph_traversed": False,
            "graph_edges_followed": False,
            "source_constructors_called": 0,
            "source_checkers_called": 0,
            "quadratures": 0,
            "numeric_values_computed": False,
            "numeric_difference_evaluated": False,
            "source_values_selected": False,
        },
    }


def _verify_dag_synthetic(dag: Mapping[str, Any]) -> dict[str, Any]:
    """Check emitted expression nodes against an independent exact polynomial model.

    Source handles are assigned synthetic polynomials by role and never read as
    numerical values.  This is a test of the emitted algebra, not source proof.
    """
    import re
    import sympy as sp
    z = sp.Symbol("z")
    zz, delta, T = sp.Rational(2, 7), sp.Rational(1, 5), sp.Rational(3, 2)
    a0, a1 = sp.Rational(2, 3), sp.Rational(3, 5)
    u = 1 + 2*z + 3*z**2 + z**3
    cf1, cf2 = 1-z+2*z**2, 2+z-z**2+z**3
    cp1, cp2 = 2+2*z-z**2, 1-2*z+z**2+2*z**3
    rp, rf = 2+z/5, 3-z/7
    definitions = dag["definitions"]
    memo = {}
    def convert(x):
        if not isinstance(x, Mapping):
            raise ValueError("synthetic DAG node must be a mapping")
        kind = x.get("kind")
        if kind == "exact_scalar_literal":
            return sp.Rational(x["value"])
        if kind == "named_expression_ref":
            name = x["name"]
            if name in {"Z", "delta", "T"}:
                return {"Z": z, "delta": delta, "T": T}[name]
            if name not in memo:
                memo[name] = convert(definitions[name]["expression"])
            return memo[name]
        if kind == "conditional_shared_native_u_handle":
            return sp.diff(u, z, x["ordinary_Z_order"])
        if kind == "ordinary_Z_jet_source_handle":
            return sp.diff(convert(x["source"]), z, x["ordinary_Z_order"])
        if x.get("first_row_weight") is True:
            return (a0, a1)[x["index"]]
        role = x.get("source_role", "")
        owner = x.get("owner")
        if role.endswith(".gap.radius") or role == "historical frame-only gap field radius":
            return rp if owner == "P" else rf
        m = re.search(r"\.current_frame\.u\[(\d+)\]$", role)
        if m:
            return sp.diff(u, z, int(m.group(1)))
        m = re.search(r"\.selected_pulse\.(C[12])\[(\d+)\]$", role)
        if m:
            polynomials = {("F", "C1"): cf1, ("F", "C2"): cf2,
                           ("P", "C1"): cp1, ("P", "C2"): cp2}
            return sp.diff(polynomials[(owner, m.group(1))], z, int(m.group(2)))
        op = x.get("operator")
        if op:
            args = [convert(v) for v in x["arguments"]]
            if op == "add": return sum(args)
            if op == "subtract": return args[0]-args[1]
            if op == "multiply": return sp.prod(args)
            if op == "divide": return args[0]/args[1]
            if op == "negate": return -args[0]
            if op == "power": return args[0]**args[1]
        raise ValueError("unsupported synthetic DAG role: " + str(role or kind or op))
    mf = -u*(a0*cf1+a1*cf2)
    mp = -u*(a0*cp1+a1*cp2)
    L = 1-delta*z*z
    def q(m):
        return -((1-delta)*z*m+(1-z*z)*sp.diff(m,z))/L
    qf, qp = q(mf), q(mp)
    sf = T*(1-delta)*z*sp.diff(qf,z)/(2*L)-T*T*qf*qf/2
    sp_ = T*(1-delta)*z*sp.diff(qp,z)/(2*L)-T*T*qp*qp/2
    targets = {"qF0":qf, "qF1":sp.diff(qf,z), "qF2":sp.diff(qf,z,2),
               "dq0":qp-qf, "dq1":sp.diff(qp-qf,z), "dq2":sp.diff(qp-qf,z,2),
               "DeltaS":sp_-sf, "DeltaS_Z1":sp.diff(sp_-sf,z),
               "physical_radius_DeltaS":rp**2*sp_-rf**2*sf,
               "physical_radius_DeltaS_Z1":sp.diff(rp**2*sp_-rf**2*sf,z),
               "physical_endpoint_DeltaQy":-rp*T*qp+rf*T*qf,
               "physical_endpoint_DeltaQy_Z1":sp.diff(-rp*T*qp+rf*T*qf,z)}
    for order in MOMENT_ORDERS:
        targets[f"mF_m{order}"] = sp.diff(mf,z,order)
        targets[f"delta_m_m{order}"] = sp.diff(mp-mf,z,order)
    checks = {name: sp.cancel(convert(_ref(name)).subs(z,zz)-expected.subs(z,zz)) == 0
              for name, expected in targets.items()}
    return {"passed": all(checks.values()), "checks": checks,
            "scope": "synthetic exact polynomial role substitution; no source function evaluated"}


def build_pulse_gap_source_factorized_delta(repo_root: str | Path) -> dict[str, Any]:
    """Build the unevaluated conditional P-minus-F source packet."""

    root = Path(repo_root)
    paths = {
        "published_actual_result": ACTUAL_RESULT,
        "physical_report": PHYSICAL_REPORT,
        "frame_archive": FRAME_ARCHIVE,
        "kernel_archive": KERNEL_ARCHIVE,
        "kernel_check": KERNEL_CHECK,
    }
    raw_values: dict[str, bytes] = {}
    input_meta: dict[str, Any] = {}
    gates: list[str] = []

    # This first pass is deliberately separate from parsing/decompression.
    # Every raw SHA, including both gzip inputs, is authenticated first.
    for key, relative in paths.items():
        raw, meta, error = _raw_packet(
            root,
            relative,
            EXPECTED_RAW_SHA256[key],
            key,
            gzip_encoded=key in {"frame_archive", "kernel_archive"},
        )
        input_meta[key] = meta
        if error is not None:
            gates.append(error)
        elif raw is not None:
            raw_values[key] = raw

    provenance = {
        "published_actual_result": _provenance_meta(input_meta["published_actual_result"]),
        "physical_report": _provenance_meta(input_meta["physical_report"]),
        "frame_archive": _provenance_meta(input_meta["frame_archive"]),
        "kernel_archive": _provenance_meta(input_meta["kernel_archive"]),
        "kernel_check": _provenance_meta(input_meta["kernel_check"]),
        "expected_raw_sha256": dict(EXPECTED_RAW_SHA256),
        "all_raw_sha256_authenticated_before_decode": not gates,
    }
    if gates:
        return _failed_packet(provenance, gates, "one or more frozen raw source hashes are unavailable or changed")

    try:
        actual_result = _parse_json(raw_values["published_actual_result"], "published_actual_result")
        physical = _parse_json(raw_values["physical_report"], "physical_report")
        kernel_check = _parse_json(raw_values["kernel_check"], "kernel_check")
        # No gzip operation occurs before this point and all five raw hashes
        # above have passed.
        frame = _decode_after_all_raw_auth(raw_values["frame_archive"], "frame_archive")
        kernel = _decode_after_all_raw_auth(raw_values["kernel_archive"], "kernel_archive")
    except SourceFactorizedDeltaInputError as exc:
        return _failed_packet(provenance, [str(exc)], "authenticated source bytes are not valid saved JSON")

    actual_gates: list[str] = []
    physical_gates: list[str] = []
    frame_gates: list[str] = []
    kernel_gates: list[str] = []

    actual_result_status = {
        "component": actual_result.get("component"),
        "status": actual_result.get("status"),
        "returned": actual_result.get("status") == "returned",
        "published_result_admissible_for_numeric_use": _optional(
            actual_result, ("result", "admissible_for_numeric_use")
        ),
        "published_result_source_only": _optional(actual_result, ("result", "source_only")),
        "published_result_source_graph_traversed": _optional(
            actual_result, ("result", "flags", "source_graph_traversed")
        ),
        "published_result_physical_frame_source_equality_proved": _optional(
            actual_result, ("result", "flags", "physical_frame_source_equality_proved")
        ),
        "source_values_selected": _optional(actual_result, ("result", "flags", "source_function_values_selected")),
    }
    if actual_result_status["component"] != "pulse_gap_actual_source_program_v1":
        actual_gates.append("published_actual_result:component_mismatch")
    if actual_result_status["status"] != "returned":
        actual_gates.append("published_actual_result:not_returned")
    if actual_result_status["published_result_admissible_for_numeric_use"] is not False:
        actual_gates.append("published_actual_result:admissible_or_missing")
    if actual_result_status["published_result_source_only"] is not True:
        actual_gates.append("published_actual_result:source_only_flag_missing")
    if actual_result_status["published_result_source_graph_traversed"] is not False:
        actual_gates.append("published_actual_result:source_graph_traversed_or_missing")
    if actual_result_status["published_result_physical_frame_source_equality_proved"] is not False:
        actual_gates.append("published_actual_result:physical_frame_equality_flag_changed")
    if actual_result_status["source_values_selected"] is not False:
        actual_gates.append("published_actual_result:source_values_selected_or_missing")

    physical_component = physical.get("component_result")
    if not isinstance(physical_component, Mapping):
        physical_gates.append("P:component_result_missing")
        physical_component = {}
    physical_binding = physical_component.get("physical_pulse_bindings")
    if not isinstance(physical_binding, Mapping):
        physical_gates.append("P:physical_pulse_bindings_missing")
        physical_binding = {}
    frame_binding = frame.get("Rp_native_C5_pulse_source_bindings")
    if not isinstance(frame_binding, Mapping):
        frame_gates.append("F:Rp_native_C5_pulse_source_bindings_missing")
        frame_binding = {}

    p_snapshot = _side_snapshot(
        physical_binding,
        owner=P_OWNER,
        namespace=P_NAMESPACE,
        binding_pointer="$.component_result.physical_pulse_bindings",
        coordinate_key="coordinates",
        gates=physical_gates,
    )
    f_snapshot = _side_snapshot(
        frame_binding,
        owner=F_OWNER,
        namespace=F_NAMESPACE,
        binding_pointer="$.Rp_native_C5_pulse_source_bindings",
        coordinate_key="pulse_coordinates",
        gates=frame_gates,
    )

    kernel_functions = kernel.get("exact_pulse_kernel_functions")
    if not isinstance(kernel_functions, Mapping):
        kernel_gates.append("K:exact_pulse_kernel_functions_missing")
        kernel_functions = {}
    first_row = _first_rows(
        kernel_functions,
        "$.exact_pulse_kernel_functions.first_row",
        K_OWNER,
        K_NAMESPACE,
        kernel_gates,
        "K.exact_pulse_kernel_functions.first_row",
    )
    first_row_refs = {
        "A0": _row_at(first_row, 0, "K.first_row", kernel_gates),
        "A1": _row_at(first_row, 1, "K.first_row", kernel_gates),
    }

    family_values = {
        "physical_report": _source_family(physical, "$.physical_report", physical_gates),
        "frame_archive": _source_family(frame, "$.frame_archive", frame_gates),
        "kernel_archive": _source_family(kernel, "$.kernel_archive", kernel_gates),
        "kernel_check": _source_family(kernel_check, "$.kernel_check", kernel_gates),
    }
    family_fields = [
        value.get("fields")
        for value in family_values.values()
        if value.get("status") == "saved_source_family_metadata"
    ]
    family_metadata_equal = bool(family_fields) and all(item == family_fields[0] for item in family_fields[1:])
    if len(family_fields) != 4:
        kernel_gates.append("saved_source_family:one_or_more_families_missing")
    elif not family_metadata_equal:
        kernel_gates.append("saved_source_family:metadata_disagrees")

    physical_prerequisites = {
        "component_result_passed": physical_component.get("passed") is True,
        "component_result_source_only": physical_component.get("source_only") is True,
        "component_result_full_goal_complete_false": physical_component.get("full_goal_complete") is False,
        "accepted_source_stage_report_present": isinstance(physical.get("accepted_source_stage_report"), str),
        "accepted_source_stage_sha256_present": isinstance(physical.get("accepted_source_stage_sha256"), str),
        "saved_source_family_present": family_values["physical_report"].get("status") == "saved_source_family_metadata",
    }
    for key, passed in physical_prerequisites.items():
        if not passed:
            physical_gates.append(f"P:accepted_prerequisite:{key}")

    kernel_open_keys = (
        "current_original_n1_full_pulse_C1_function_adapter_installed",
        "n1_complete_signed_defect_C1_function_adapter_installed",
        "n1_selected_linear_coefficient_functions_installed",
        "n1_repaired_own_M_V_P_installed",
        "positive_order_five_moment_repair_closed",
        "temporal_recursion",
        "full_corrected_Navier_Stokes_solution",
    )
    kernel_prerequisites = {
        "check_all_passed": kernel_check.get("all_passed") is True,
        "exact_kernel_installed": kernel_check.get("current_original_n1_exact_pulse_kernel_functions_installed") is True,
        "source_family_present": family_values["kernel_check"].get("status") == "saved_source_family_metadata",
        "all_open_future_adapters_false": all(kernel_check.get(key) is False for key in kernel_open_keys),
        "exact_pulse_kernel_functions_present": isinstance(kernel.get("exact_pulse_kernel_functions"), Mapping),
    }
    for key, passed in kernel_prerequisites.items():
        if not passed:
            kernel_gates.append(f"K:accepted_prerequisite:{key}")

    p_u = p_snapshot["current_frame"]["u"]
    f_u = f_snapshot["current_frame"]["u"]
    p_c1 = p_snapshot["selected_pulse"]["C1"]
    f_c1 = f_snapshot["selected_pulse"]["C1"]
    p_c2 = p_snapshot["selected_pulse"]["C2"]
    f_c2 = f_snapshot["selected_pulse"]["C2"]

    delta_C: dict[str, Any] = {}
    for name, p_rows, f_rows in (("C1", p_c1, f_c1), ("C2", p_c2, f_c2)):
        p_ref = _row_at(p_rows, 0, f"P.selected_pulse.{name}", physical_gates)
        f_ref = _row_at(f_rows, 0, f"F.selected_pulse.{name}", frame_gates)
        delta_C[name] = {
            "P": p_ref,
            "F": f_ref,
            "expression": _expr("subtract", p_ref, f_ref, name=f"delta_{name}"),
            "formula_text": f"delta_{name}=C_P{name[-1]}-C_F{name[-1]}",
            "ordinary_Z_rows": {
                str(order): _expr(
                    "subtract",
                    _row_at(p_rows, order, f"P.selected_pulse.{name}", physical_gates),
                    _row_at(f_rows, order, f"F.selected_pulse.{name}", frame_gates),
                    name=f"delta_{name}_Z{order}",
                )
                for order in MOMENT_ORDERS
            },
            "P_ordinary_Z_rows": {
                str(order): _row_at(p_rows, order, f"P.selected_pulse.{name}", physical_gates)
                for order in MOMENT_ORDERS
            },
            "F_ordinary_Z_rows": {
                str(order): _row_at(f_rows, order, f"F.selected_pulse.{name}", frame_gates)
                for order in MOMENT_ORDERS
            },
            "source_values_selected": False,
        }

    A0, A1 = first_row_refs["A0"], first_row_refs["A1"]
    cF_jets = {
        f"cF_Z{order}": _expr(
            "add",
            _expr("multiply", A0, delta_C["C1"]["F_ordinary_Z_rows"][str(order)]),
            _expr("multiply", A1, delta_C["C2"]["F_ordinary_Z_rows"][str(order)]),
            name=f"cF_Z{order}",
            ordinary_Z_order=order,
            formula_text=f"cF_Z{order}=A0 CF1_Z{order}+A1 CF2_Z{order}",
        )
        for order in MOMENT_ORDERS
    }
    delta_c_jets = {
        f"delta_c_Z{order}": _expr(
            "add",
            _expr(
                "multiply",
                A0,
                delta_C["C1"]["ordinary_Z_rows"][str(order)],
            ),
            _expr(
                "multiply",
                A1,
                delta_C["C2"]["ordinary_Z_rows"][str(order)],
            ),
            name=f"delta_c_Z{order}",
            ordinary_Z_order=order,
            formula_text=f"delta_c_Z{order}=A0 delta_C1_Z{order}+A1 delta_C2_Z{order}",
        )
        for order in MOMENT_ORDERS
    }
    cF_refs = {name: _ref(name, ordinary_Z_order=int(name.rsplit("Z", 1)[1])) for name in cF_jets}
    delta_c_refs = {
        name: _ref(name, ordinary_Z_order=int(name.rsplit("Z", 1)[1]))
        for name in delta_c_jets
    }
    cF = _ref("cF_Z0", ordinary_Z_order=0)
    delta_c = _ref("delta_c_Z0", ordinary_Z_order=0)
    moment_rows = _moment_rows(
        p_u,
        f_u,
        cF_refs,
        delta_c_refs,
        physical_gates + frame_gates + kernel_gates,
    )
    dag = _build_q_dag(moment_rows)
    normalized = dag["normalized_sources"]
    radius = _radius_factorized_expressions(p_snapshot["radius"], f_snapshot["radius"])
    z_formulas = _factorized_Z1_Z2(radius, dag)
    combined_definitions = dict(dag["definitions"])
    for name, expression in {**cF_jets, **delta_c_jets}.items():
        combined_definitions[name] = _named_definition(
            name,
            expression,
            ordinary_Z_order=int(name.rsplit("Z", 1)[1]),
        )
    for family_rows in (moment_rows["mF"], moment_rows["delta_m"]):
        for name, row in family_rows.items():
            combined_definitions[name] = _named_definition(
                name,
                row["formula"],
                ordinary_Z_order=row["ordinary_Z_order"],
                conditional_on=row["conditional_on"],
            )
    combined_definitions.update(z_formulas["expression_dag"]["definitions"])
    dag = dict(dag)
    dag["definitions"] = combined_definitions

    algebra_fixture = _verify_dag_synthetic(dag)
    if not algebra_fixture["passed"]:
        raise ValueError("synthetic DAG algebra verification failed")

    all_gates = sorted(set(actual_gates + physical_gates + frame_gates + kernel_gates + [
        "same_native_u_owner_unproved",
        "same_first_row_kernel_owner_unproved",
        "source_semantics_unproved",
        "derivative_semantics_unproved",
        "mu_common_owner_unproved",
        "delta_common_and_Z_independent_unproved",
        "common_field_scale_and_coordinate_Z_independence_unproved",
        "first_row_weights_Z_independence_unproved",
        "geometry_common_owner_unproved",
        "F_field_radius_owner_unproved",
        "P_F_W_owner_unproved",
        "P_F_incoming_owner_unproved",
        "P_F_future_owner_unproved",
    ]))
    return {
        "schema": SCHEMA,
        "component": COMPONENT,
        "status": "conditional_source_expression_not_admissible_for_numeric_use",
        "source_only": True,
        "formal_source_expression_defined": True,
        "admissible_for_numeric_use": False,
        "numeric_values_computed": False,
        "numeric_difference_evaluated": False,
        "source_values_selected": False,
        "inputs": {
            "published_actual_result": provenance["published_actual_result"],
            "physical_report": provenance["physical_report"],
            "frame_archive": provenance["frame_archive"],
            "kernel_archive": provenance["kernel_archive"],
            "kernel_check": provenance["kernel_check"],
            "all_raw_sha256_authenticated_before_decode": True,
            "gzip_decoded_after_all_raw_authentication": True,
        },
        "published_actual_result": actual_result_status,
        "source_family": {
            "metadata": family_values,
            "saved_metadata_equal": family_metadata_equal,
            "source_identity_proved": False,
            "source_values_selected": False,
        },
        "accepted_prerequisites": {
            "physical_report": physical_prerequisites,
            "kernel": kernel_prerequisites,
            "source_prerequisite_semantics_proved": False,
        },
        "sides": {"P": p_snapshot, "F": f_snapshot},
        "kernel": {
            "owner": K_OWNER,
            "namespace": K_NAMESPACE,
            "exact_pulse_kernel_functions": {
                "first_row": first_row,
                "A0": A0,
                "A1": A1,
                "source_values_selected": False,
            },
        },
        "coefficient_delta": {
            "delta_C1": delta_C["C1"],
            "delta_C2": delta_C["C2"],
            "cF": cF,
            "delta_c": delta_c,
            "cF_Z0_to_Z3": cF_jets,
            "delta_c_Z0_to_Z3": delta_c_jets,
            "formula_text": [
                "delta_Cj=C_Pj-C_Fj",
                "cF=cF_Z0=A0 CF1+A1 CF2",
                "cF_Zj=A0 CF1_Zj+A1 CF2_Zj",
                "delta_c=delta_c_Z0=A0 delta_C1+A1 delta_C2",
                "delta_c_Zj=A0 delta_C1_Zj+A1 delta_C2_Zj",
            ],
            "source_values_selected": False,
        },
        "moment_transport": {
            "mF": moment_rows["mF"],
            "delta_m": moment_rows["delta_m"],
            "ordinary_Z_orders": list(MOMENT_ORDERS),
            "source_values_selected": False,
        },
        "delta_m": {
            "ordinary_Z_orders": list(MOMENT_ORDERS),
            "rows": moment_rows["delta_m"],
            "conditional_formula": "delta_m_j=-sum(k=0..j) binom(j,k) u_k delta_c_(j-k)",
            "same_native_u_and_first_row_required": True,
            "source_values_selected": False,
        },
        "mF": {
            "ordinary_Z_orders": list(MOMENT_ORDERS),
            "rows": moment_rows["mF"],
            "conditional_formula": "mF_j=-sum(k=0..j) binom(j,k) uF_k cF_(j-k)",
            "source_values_selected": False,
        },
        "expression_dag": dag,
        "independent_algebra_fixture": algebra_fixture,
        "normalized_sources": normalized,
        "L": {
            "expression": _ref("L"),
            "formula_text": "L=1-delta*Z^2",
            "definition_name": "L",
            "delta_Z_source_semantics_proved": False,
            "source_values_selected": False,
        },
        "factorized_P_minus_F": {
            "delta_Q": z_formulas["delta_Q"],
            "delta_Qy_source": z_formulas["delta_Qy_source"],
            "S_linear": z_formulas["S_linear"],
            "S_quadratic": z_formulas["S_quadratic"],
            "S_linear_Z1": z_formulas["S_linear_Z1"],
            "S_quadratic_Z1": z_formulas["S_quadratic_Z1"],
            "S_F": z_formulas["S_F"],
            "DeltaS": z_formulas["DeltaS"],
            "S_F_Z1": z_formulas["S_F_Z1"],
            "DeltaS_Z1": z_formulas["DeltaS_Z1"],
            "T_Z_condition": "T_Z=0",
            "factorization_not_numeric_subtraction": True,
            "source_values_selected": False,
        },
        "actual_physical_radius_difference_endpoint_terms": z_formulas["actual_physical_radius_difference_endpoint_terms"],
        "conditional_proof_gates": {
            "all_proved": False,
            "source_family_metadata_equal_only": family_metadata_equal,
            "requirements": [
                "published_actual_result_is_returned_source_only_and_inadmissible",
                "physical_and_frame_current_frame_u_have_four_direct_rows",
                "physical_and_frame_selected_C1_C2_rows_are_present",
                "kernel_exact_pulse_kernel_functions_first_row_has_two_refs",
                "same_native_u_owner",
                "same_first_row_kernel_owner",
                "first_row_weights_ordinary_Z_independent",
                "common_mu_owner",
                "common_delta_and_ordinary_Z_independence",
                "common_gap_coordinate_and_field_scale_T_Z_independence",
                "common_geometry_owner",
                "common_source_semantics_and_normalization",
                "ordinary_Z1_and_Z2_derivative_transport",
                "historical_F_field_radius_owner",
                "P_F_W_incoming_future_common_owners",
                "common_coordinates_measures_units_and_endpoint_orientation",
                "finite_L_nonzero_domain",
            ],
            "physical_frame_source_equality_proved": False,
            "F_field_radius_identity_proved": False,
            "source_values_selected": False,
        },
        "retained_prior_errors": {
            "old_errors_removed": False,
            "numerical_caps_changed": False,
            "physical_pulse_corrections_removed": False,
            "independent_reference_offset_selected": False,
            "source_values_selected": False,
        },
        "unresolved_gates": all_gates,
        "scope_audit": {
            "direct_literal_field_retrieval_only": True,
            "source_graph_traversed": False,
            "graph_edges_followed": False,
            "graph_closure_followed": False,
            "source_constructors_called": 0,
            "source_checkers_called": 0,
            "quadratures": 0,
            "integrals": 0,
            "mathematical_evaluations": 0,
            "numeric_values_computed": False,
            "numeric_difference_evaluated": False,
            "source_values_selected": False,
            "gzip_archives_decoded": 2,
            "gzip_raw_hashes_authenticated_before_decode": True,
            "old_family_replayed": False,
            "old_errors_removed": False,
            "graph_id_arithmetic": False,
            "cross_namespace_id_comparison": False,
        },
    }


# A short alias makes the intended leaf callable without introducing a second
# implementation or importing an older family.
build_source_factorized_delta = build_pulse_gap_source_factorized_delta


__all__ = [
    "build_pulse_gap_source_factorized_delta",
    "build_source_factorized_delta",
    "SourceFactorizedDeltaInputError",
]
