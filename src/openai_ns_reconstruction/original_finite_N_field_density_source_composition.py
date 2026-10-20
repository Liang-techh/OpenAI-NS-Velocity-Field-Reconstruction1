"""Source-only finite-``N`` field and density composition.

The accepted R110--Rm projection supplies source functions and function ports;
it does not supply finite-``N`` fields.  This module composes those fields from
the authenticated C5 source recipe and checks the coefficientwise C2 identity
``C2 = N * C5``.  It never imports a producer, provider, graph evaluator,
transport operator, quadrature, recurrence or norm implementation.

``derive`` consumes a saved graph and authenticated source records as data.
The ordinary-jet helpers are deliberately independent of the reconstruction
graph.  A caller may supply symbolic scalar callbacks (for example, callbacks
that create source graph nodes) or ordinary Python scalar operations.  Jet
rows are ordinary derivatives, so the product rule contains the factor two in
the second row.  A callback result is never silently promoted to an exact
function value; the ``provenance`` field is retained on every composed jet.
"""

from __future__ import annotations

import ast
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
from fractions import Fraction
from typing import Any


SCHEMA = "original_finite_N_field_density_source_composition_v1"
EXPECTED_NODE_COUNT = 205300
EXPECTED_NODES_DIGEST = "a2dcc782ba3033ff6ca5f54a66632ac6ef007bb0c5b9960f10567d44ef2059a3"

PARENT_REPORT = "lei_ren_part1_paper_compliant_current_original_Rm_physical_Qy_Z_calculus.json"
PARENT_SHA256 = "8903b59267bae01d5182bf2dcc84ed3c02c44eb8c7d3dd9d7b14af67ac092746"
LATEST_REPORT = "lei_ren_part1_paper_compliant_current_original_R110_Rm_source_projection.json"
LATEST_SHA256 = "08ef44c6971d07a4c6583d481ef3a0ce1b37bf069995a732bdf2a9a7ccfa074e"

PREFIX = "lei_ren_part1_paper_compliant_current_original_"
SOURCE_FAMILY = {
    "actual_five_defect_family_sha256": "3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894",
    "implicit_source_sha256": "5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae",
    "datum_enclosure_sha256": "dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010",
}

FIRST_TEN_CHARTS = (
    "first_micro",
    "second_micro",
    "frozen_macro",
    "first_switch",
    "second_switch",
    "post_power",
    "long_reshape",
    "reference",
    "restoration",
    "postrestore",
)
OUTER_SIX_CHARTS = (
    "Rh_reference",
    "O2_slope",
    "O2_axial",
    "O2_buffer",
    "O3_transition",
    "O3_power",
)
ORDINARY_Z_ORDERS = (0, 1, 2)
SIGNED_DENSITIES = ("m", "h", "k", "e", "p")

# Plain source and receipt records are intentionally used here.  No gzip
# archive is decoded by this module.  A source record is a mapping containing
# ``path``, ``sha256`` and UTF-8 ``source_text``.
RECORD_HASHES = {
    PARENT_REPORT: PARENT_SHA256,
    LATEST_REPORT: LATEST_SHA256,
    PREFIX + "C2_phase_primitives.py": "84f5da3eb64096f82df88dd828b22b09bce18192a03925d91cfee0b193902d2c",
    PREFIX + "C2_phase_primitives_check.json": "76ed1af50ec799c49081ab67881b189b42e0d3137bcca770a67766d4111dff8b",
    PREFIX + "C2_density_transport.py": "8bfa32c87c73c5bedf3056c7171707ae934dafd950bd0df874dbd3b3450379d1",
    PREFIX + "C2_density_transport_check.json": "c4d2c62dd67e9dd7e171fabbdb9ba18308e420c5e0c98d01147c708840b8f583",
    PREFIX + "prescribed_loop_C5_functions.py": "69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d",
    PREFIX + "prescribed_loop_C5_functions_check.json": "2c77e80df04fae615f8ea909fc7b606bfccfa27659cb67744ae8d84020b1a9e7",
    PREFIX + "prescribed_bridge_C6_functions.py": "8ab2c26092db06453ed5a6047d8905e37fdad17eaf9c4947048498551425c605",
    PREFIX + "prescribed_bridge_C6_functions_check.json": "9a0969a8cd4a5989dbd6272d6dd223533400113824f9352b4e57a83196bff61e",
    PREFIX + "prescribed_R100_Rm_C6_functions.py": "70dd6d6022f68d601ac4f0fbaf38614bbaf215baab272f134f1c5a159af65d32",
    PREFIX + "prescribed_R100_Rm_C6_functions_check.json": "494caf482bffe5475c011a35a1b7d95fd02725b32eb55653f9ca25984dd98ad6",
    PREFIX + "whole_Z_bridge_source.py": "df8b204d023c1d8762da16d07190601ebd1e11447077b01a734bba8baa595a6a",
    PREFIX + "whole_Z_bridge_source_check.json": "91d9ba8119d9dba85322ea33791ac3126c22ad85dfd94bf9a961c49fc18268a3",
    PREFIX + "whole_Z_bridge_correlated_functions.py": "bb29100e671b6a1927f9897385056645a17d859108e537eec760a3904d66f832",
    PREFIX + "whole_Z_bridge_correlated_functions_check.json": "bbf33386a8a591eef946ca5f3f262dd92ac6dd1490bc8bfc7421ee1c6eca49da",
    PREFIX + "whole_Z_all_N_source_functions.py": "f0f6058bf77e52d9f5f04db3259ca9f737dc0fdb9281e2650b7a6ae5ca3f25fe",
    PREFIX + "whole_Z_all_N_source_functions_check.json": "b141175a0b476d2c7b325028ef6d45696ca5b6e52cdafcac775d2d98e2dd8ccb",
    PREFIX + "whole_Z_all_N_switch_functions.py": "8a3e9581a393fcbc59fdad7f24055b0176eb7edca8a76676a35c1b23b114922c",
    PREFIX + "whole_Z_all_N_switch_functions_check.json": "3b69567b7e51aef13be5f35492ebf8f4013b43b8d15c2d9f7d175bb4292c0497",
    PREFIX + "whole_Z_all_N_long_patch_functions.py": "847120150d0c8a83c2861f440c8f69514358b55633cfe8964f6de25bd5aa3ddc",
    PREFIX + "whole_Z_all_N_long_patch_functions_check.json": "2162a7b7fbf322a53065e9ba118e9466a7b17a3baff86d569ed35935b18e8a02",
    PREFIX + "true_Rm_patch_finite_N_C5_functions.py": "c4445ef8acfb33d29cf9022862f8b6396e6aa2f38bda7467409b5378a5deaf15",
    PREFIX + "true_Rm_patch_finite_N_C5_functions_check.json": "efd3b9a6fdfd4522ead1806d4ef4250864ef670909ca047d0fa85d9ce23bf207",
    PREFIX + "true_Rm_patch_C6_functions.py": "a6f1cb5deaf17e9482aa964aef7867aadbc815e4d703a415e4945cf665f14002",
    PREFIX + "true_Rm_patch_C6_functions_check.json": "0190fcf738fb8add5be85cf0767464a39cc0d18240e6b36ad4ec7930dce38aa7",
    "lei_ren_part1_paper_compliant_current_generic_loop_function_sources.py": "d694b8ee8a4bf1c10090d45a7f9d4d6de0792d02a1ebdd6b21872b8514a5b26f",
    "lei_ren_part1_paper_compliant_current_generic_loop_function_sources_check.json": "ac18d2621c67002769f8e0c76613ae91de54239699edd3801cccfaa052c19d07",
}

# The patch records are authenticated when present, but their absence keeps
# the patch row unresolved.  The first ten rows never depend on them.
OPTIONAL_PATCH_RECORDS = frozenset(
    {
        PREFIX + "true_Rm_patch_finite_N_C5_functions.py",
        PREFIX + "true_Rm_patch_finite_N_C5_functions_check.json",
        PREFIX + "true_Rm_patch_C6_functions.py",
        PREFIX + "true_Rm_patch_C6_functions_check.json",
    }
)


@dataclass(frozen=True)
class OrdinaryJet:
    """Rows ``(f, f_Z, f_ZZ, ...)`` of an ordinary source jet.

    ``provenance`` is intentionally explicit.  Values from an interval
    enclosure or an arbitrary callback remain ``unclassified`` unless the
    caller labels them.  The composer never changes ``unclassified`` into an
    exact function label.
    """

    rows: tuple[Any, ...]
    provenance: str = "unclassified"
    units: str | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.rows, tuple) or not self.rows:
            raise ValueError("ordinary jet rows must be a non-empty tuple")
        if self.provenance not in {"source_function", "interval_enclosure", "symbolic", "unclassified"}:
            raise ValueError("unknown ordinary jet provenance")

    def __getitem__(self, index: int) -> Any:
        return self.rows[index]

    def __len__(self) -> int:
        return len(self.rows)


@dataclass(frozen=True)
class ScalarOperations:
    """Scalar callbacks used by :func:`compose_finite_N_ordinary_jets`.

    The callbacks are binary for ``add``, ``sub`` and ``mul`` so graph
    builders with either binary or variadic primitives can adapt them with a
    small wrapper.  ``exprel_prime`` and ``exprel_second`` are required for a
    symbolic C2 composition unless the caller supplies ordinary numeric
    scalars, where stable entire formulas are used at ``A=0``.
    """

    add: Callable[[Any, Any], Any] = lambda a, b: a + b
    sub: Callable[[Any, Any], Any] = lambda a, b: a - b
    mul: Callable[[Any, Any], Any] = lambda a, b: a * b
    div: Callable[[Any, Any], Any] = lambda a, b: a / b
    neg: Callable[[Any], Any] = lambda a: -a
    exp: Callable[[Any], Any] = math.exp
    expm1: Callable[[Any], Any] = math.expm1
    exprel: Callable[[Any], Any] | None = None
    exprel_prime: Callable[[Any], Any] | None = None
    exprel_second: Callable[[Any], Any] | None = None
    constant: Callable[[Any], Any] = lambda value: value
    is_zero: Callable[[Any], bool] | None = None


def _ops(value: ScalarOperations | Mapping[str, Callable[..., Any]] | None) -> ScalarOperations:
    if value is None:
        return ScalarOperations()
    if isinstance(value, ScalarOperations):
        return value
    if isinstance(value, Mapping):
        fields = {name: value[name] for name in ScalarOperations.__dataclass_fields__ if name in value}
        return ScalarOperations(**fields)
    raise TypeError("ops must be ScalarOperations or a callback mapping")


def _scalar_is_zero(value: Any, ops: ScalarOperations) -> bool:
    """Check a source scalar is exactly zero without guessing for symbols."""

    if ops.is_zero is not None:
        try:
            return bool(ops.is_zero(value))
        except Exception:
            return False
    try:
        result = value == 0
    except Exception:
        return False
    return isinstance(result, bool) and result


def _jet(value: OrdinaryJet | Sequence[Any] | Any, *, provenance: str = "unclassified", units: str | None = None) -> OrdinaryJet:
    if isinstance(value, OrdinaryJet):
        return value
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return OrdinaryJet(tuple(value), provenance=provenance, units=units)
    return OrdinaryJet((value,), provenance=provenance, units=units)


def _kind(*jets: OrdinaryJet) -> str:
    names = {jet.provenance for jet in jets}
    if len(names) == 1:
        return next(iter(names))
    return "unclassified"


def _constant(ops: ScalarOperations, value: Any) -> Any:
    return ops.constant(value)


def _fold(fn: Callable[[Any, Any], Any], values: Sequence[Any], *, identity: Any | None = None) -> Any:
    if not values:
        if identity is None:
            raise ValueError("empty scalar fold")
        return identity
    result = values[0]
    for value in values[1:]:
        result = fn(result, value)
    return result


def _jet_binary(left: OrdinaryJet, right: OrdinaryJet, fn: Callable[[Any, Any], Any], *, units: str | None = None) -> OrdinaryJet:
    if len(left) != len(right):
        raise ValueError("ordinary jets must have equal order")
    return OrdinaryJet(tuple(fn(a, b) for a, b in zip(left.rows, right.rows)), provenance=_kind(left, right), units=units)


def _jet_add(left: OrdinaryJet, right: OrdinaryJet, ops: ScalarOperations, *, units: str | None = None) -> OrdinaryJet:
    return _jet_binary(left, right, ops.add, units=units)


def _jet_sub(left: OrdinaryJet, right: OrdinaryJet, ops: ScalarOperations, *, units: str | None = None) -> OrdinaryJet:
    return _jet_binary(left, right, ops.sub, units=units)


def _jet_mul(left: OrdinaryJet, right: OrdinaryJet, ops: ScalarOperations, *, units: str | None = None) -> OrdinaryJet:
    if len(left) != len(right):
        raise ValueError("ordinary jets must have equal order")
    rows: list[Any] = []
    for order in range(len(left)):
        if order == 0:
            rows.append(ops.mul(left[0], right[0]))
        elif order == 1:
            rows.append(ops.add(ops.mul(left[1], right[0]), ops.mul(left[0], right[1])))
        elif order == 2:
            cross = ops.mul(_constant(ops, 2), ops.mul(left[1], right[1]))
            rows.append(
                ops.add(
                    ops.add(ops.mul(left[2], right[0]), cross),
                    ops.mul(left[0], right[2]),
                )
            )
        else:
            raise ValueError("ordinary composition supports only Z orders 0, 1 and 2")
    return OrdinaryJet(tuple(rows), provenance=_kind(left, right), units=units)


def _jet_scale(jet: OrdinaryJet, scalar: Any, ops: ScalarOperations, *, units: str | None = None) -> OrdinaryJet:
    return OrdinaryJet(tuple(ops.mul(scalar, row) for row in jet.rows), provenance=jet.provenance, units=units or jet.units)


def _jet_exp(jet: OrdinaryJet, ops: ScalarOperations, *, units: str | None = None) -> OrdinaryJet:
    rows = [ops.exp(jet[0])]
    if len(jet) >= 2:
        rows.append(ops.mul(rows[0], jet[1]))
    if len(jet) >= 3:
        rows.append(ops.mul(rows[0], ops.add(jet[2], ops.mul(jet[1], jet[1]))))
    if len(jet) > 3:
        raise ValueError("ordinary composition supports only Z orders 0, 1 and 2")
    return OrdinaryJet(tuple(rows), provenance=jet.provenance, units=units)


def _jet_expm1(jet: OrdinaryJet, ops: ScalarOperations, *, units: str | None = None) -> OrdinaryJet:
    rows = [ops.expm1(jet[0])]
    if len(jet) >= 2:
        exponential = ops.exp(jet[0])
        rows.append(ops.mul(exponential, jet[1]))
    if len(jet) >= 3:
        exponential = ops.exp(jet[0])
        rows.append(ops.mul(exponential, ops.add(jet[2], ops.mul(jet[1], jet[1]))))
    if len(jet) > 3:
        raise ValueError("ordinary composition supports only Z orders 0, 1 and 2")
    return OrdinaryJet(tuple(rows), provenance=jet.provenance, units=units)


def _numeric_exprel(value: Any, ops: ScalarOperations) -> Any:
    if value == 0:
        return _constant(ops, 1)
    return ops.div(ops.expm1(value), value)


def _numeric_exprel_prime(value: Any, ops: ScalarOperations) -> Any:
    if value == 0:
        return _constant(ops, 0.5)
    numerator = ops.add(ops.mul(ops.sub(value, _constant(ops, 1)), ops.exp(value)), _constant(ops, 1))
    return ops.div(numerator, ops.mul(value, value))


def _numeric_exprel_second(value: Any, ops: ScalarOperations) -> Any:
    if value == 0:
        return _constant(ops, 1 / 3)
    value2 = ops.mul(value, value)
    numerator = ops.sub(ops.mul(ops.add(ops.sub(value2, ops.mul(_constant(ops, 2), value)), _constant(ops, 2)), ops.exp(value)), _constant(ops, 2))
    return ops.div(numerator, ops.mul(value2, value))


def _jet_exprel(jet: OrdinaryJet, ops: ScalarOperations, *, units: str | None = None) -> OrdinaryJet:
    exprel = ops.exprel or (lambda value: _numeric_exprel(value, ops))
    prime = ops.exprel_prime or (lambda value: _numeric_exprel_prime(value, ops))
    second = ops.exprel_second or (lambda value: _numeric_exprel_second(value, ops))
    rows = [exprel(jet[0])]
    if len(jet) >= 2:
        rows.append(ops.mul(prime(jet[0]), jet[1]))
    if len(jet) >= 3:
        first = ops.mul(second(jet[0]), ops.mul(jet[1], jet[1]))
        rows.append(ops.add(first, ops.mul(prime(jet[0]), jet[2])))
    if len(jet) > 3:
        raise ValueError("ordinary composition supports only Z orders 0, 1 and 2")
    return OrdinaryJet(tuple(rows), provenance=jet.provenance, units=units)


def _require_three(*jets: OrdinaryJet) -> None:
    if any(len(jet) != 3 for jet in jets):
        raise ValueError("finite-N composition requires ordinary Z orders 0, 1 and 2")


def _inverse_N_jet(N: OrdinaryJet, ops: ScalarOperations) -> tuple[OrdinaryJet, Any]:
    """Return ``(1/N, (1/N)_Z, (1/N)_ZZ)`` under the constant-N contract.

    The source contract makes ``N`` one positive Z-independent integer.  A
    symbolic callback must therefore supply a zero recognizer when its zero
    rows are graph objects; this helper never assumes that a name implies a
    constant.  A changing ``N`` is rejected before any field composition.
    """

    if not (_scalar_is_zero(N[1], ops) and _scalar_is_zero(N[2], ops)):
        raise ValueError("N must have exact zero ordinary Z and Z2 rows")
    one = _constant(ops, 1)
    zero = _constant(ops, 0)
    inverse_scalar = ops.div(one, N[0])
    inverse = _jet_scale(
        OrdinaryJet((one, zero, zero), provenance=N.provenance, units="1"),
        inverse_scalar,
        ops,
        units="1/N",
    )
    return inverse, inverse_scalar


def compose_finite_N_ordinary_jets(
    E: OrdinaryJet | Sequence[Any] | Any,
    V: OrdinaryJet | Sequence[Any] | Any,
    A: OrdinaryJet | Sequence[Any] | Any,
    B_over_Pstar: OrdinaryJet | Sequence[Any] | Any,
    N: OrdinaryJet | Sequence[Any] | Any,
    *,
    ops: ScalarOperations | Mapping[str, Callable[..., Any]] | None = None,
) -> dict[str, Any]:
    """Compose the prescribed C5 finite-``N`` fields and five densities.

    Inputs are ordinary ``Z`` jets.  ``N`` must be a constant jet in a source
    graph; the function does not infer that property from a name.  The caller
    can provide callbacks that make source nodes instead of Python values.
    """

    algebra = _ops(ops)
    E, V, A, B_over_Pstar, N = (_jet(v) for v in (E, V, A, B_over_Pstar, N))
    _require_three(E, V, A, B_over_Pstar, N)
    invN, invN_scalar = _inverse_N_jet(N, algebra)
    x = _jet_mul(A, invN, algebra, units="A/N")
    exp_x = _jet_exp(x, algebra, units="1")
    expm1_x = _jet_expm1(x, algebra, units="1")
    r = _jet_mul(E, expm1_x, algebra, units="E")
    s = _jet_scale(B_over_Pstar, invN_scalar, algebra, units="B_over_Pstar/N")
    E_N = _jet_mul(E, exp_x, algebra, units="E")
    V_N = _jet_add(V, s, algebra, units="V")
    half = algebra.div(_constant(algebra, 1), _constant(algebra, 2))
    r2_half = _jet_scale(_jet_mul(r, r, algebra), half, algebra, units="E^2")
    s2 = _jet_mul(s, s, algebra, units="B_over_Pstar^2/N^2")
    inc = {
        "m": s,
        "h": r,
        "k": _jet_add(_jet_add(_jet_mul(V, r, algebra), _jet_mul(E, s, algebra), algebra), _jet_mul(r, s, algebra), algebra),
        "e": _jet_sub(_jet_add(_jet_scale(_jet_mul(V, s, algebra), _constant(algebra, 2), algebra), s2, algebra), _jet_add(_jet_mul(E, r, algebra), r2_half, algebra), algebra),
        "p": _jet_add(_jet_mul(E, r, algebra), r2_half, algebra),
    }
    return {
        "E_N": E_N,
        "V_N": V_N,
        "deltaE": r,
        "deltaV": s,
        "r": r,
        "s": s,
        "A_over_N": x,
        "five_signed_density_increments": inc,
        "units": {
            "E": "normalized E",
            "V": "normalized V",
            "B_over_Pstar": "normalized B_over_Pstar",
            "N": "same positive Z-independent integer",
            "ordinary_Z": "derivatives in global fixed phase Z",
        },
    }


def compose_c2_density_ordinary_jets(
    E: OrdinaryJet | Sequence[Any] | Any,
    V: OrdinaryJet | Sequence[Any] | Any,
    A: OrdinaryJet | Sequence[Any] | Any,
    B: OrdinaryJet | Sequence[Any] | Any,
    N: OrdinaryJet | Sequence[Any] | Any,
    *,
    ops: ScalarOperations | Mapping[str, Callable[..., Any]] | None = None,
) -> dict[str, Any]:
    """Compose C2 signed densities using full exponential derivatives."""

    algebra = _ops(ops)
    E, V, A, B, N = (_jet(v) for v in (E, V, A, B, N))
    _require_three(E, V, A, B, N)
    one = _constant(algebra, 1)
    invN, invN_scalar = _inverse_N_jet(N, algebra)
    argument = _jet_mul(A, invN, algebra, units="A/N")
    exp_x = _jet_exp(argument, algebra, units="1")
    expm1_x = _jet_expm1(argument, algebra, units="1")
    # G = A*exprel(A/N) = N*expm1(A/N) is entire, including A = 0.
    # Writing its first two ordinary rows directly avoids cancellation in
    # numerical exprel derivatives while preserving the full product rule.
    G = OrdinaryJet(
        (
            algebra.mul(N[0], expm1_x[0]),
            algebra.mul(exp_x[0], A[1]),
            algebra.mul(
                exp_x[0],
                algebra.add(A[2], algebra.mul(A[1], algebra.div(A[1], N[0]))),
            ),
        ),
        provenance=_kind(A, N),
        units="A",
    )
    F = _jet_mul(E, G, algebra, units="E*A")
    two = _constant(algebra, 2)
    half = algebra.div(one, _constant(algebra, 2))
    first = {
        "m": B,
        "h": F,
        "k": _jet_add(_jet_mul(V, F, algebra), _jet_mul(E, B, algebra), algebra),
        "e": _jet_sub(_jet_scale(_jet_mul(V, B, algebra), two, algebra), _jet_mul(E, F, algebra), algebra),
        "p": _jet_mul(E, F, algebra),
    }
    zero = OrdinaryJet((_constant(algebra, 0),) * 3, provenance=_kind(B, F))
    F2 = _jet_mul(F, F, algebra)
    B2 = _jet_mul(B, B, algebra)
    second = {
        "m": zero,
        "h": zero,
        "k": _jet_mul(F, B, algebra),
        "e": _jet_sub(B2, _jet_scale(F2, half, algebra), algebra),
        "p": _jet_scale(F2, half, algebra),
    }
    density = {name: _jet_add(first[name], _jet_scale(second[name], invN_scalar, algebra), algebra) for name in SIGNED_DENSITIES}
    return {
        "F": F,
        "first": first,
        "second": second,
        "density": density,
        "units": {
            "E": "normalized E",
            "V": "normalized V",
            "B": "normalized B; equals B_over_Pstar only after source-unit admission",
            "N": "same positive Z-independent integer",
            "ordinary_Z": "derivatives in global fixed phase Z",
        },
    }


def finite_N_source_recipe() -> dict[str, Any]:
    """Return the source-bound formulas used by later field evaluators."""

    return {
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "fields": {
            "r": "E*expm1(A/N)",
            "s": "B_over_Pstar/N",
            "E_N": "E*exp(A/N)",
            "V_N": "V+s",
            "deltaE": "r",
            "deltaV": "s",
        },
        "C5_density": {
            "m": "s",
            "h": "r",
            "k": "V*r+E*s+r*s",
            "e": "2*V*s+s^2-E*r-r^2/2",
            "p": "E*r+r^2/2",
        },
        "C2_density": {
            "F": "E*A*exprel(A/N)=N*E*expm1(A/N)",
            "first": {
                "m": "B",
                "h": "F",
                "k": "V*F+E*B",
                "e": "2*V*B-E*F",
                "p": "E*F",
            },
            "second": {
                "m": "0",
                "h": "0",
                "k": "F*B",
                "e": "B^2-F^2/2",
                "p": "F^2/2",
            },
            "density": "first+second/N",
        },
        "identity": {
            "relation": "C2=N*C5",
            "conditions": [
                "B_C2=B_over_Pstar_C5 in the same normalized Pstar units",
                "F=N*r as an entire identity, including A=0",
                "same ordinary Z jet product rule through order 2",
            ],
            "coefficientwise_fields": list(SIGNED_DENSITIES),
        },
        "source_contract": {
            "E": "normalized E",
            "V": "normalized V",
            "B_over_Pstar": "normalized loop B divided by the accepted physical Pstar",
            "N": "one positive Z-independent integer shared by all admitted charts",
            "phase": "fixed global phase while taking ordinary Z derivatives",
        },
    }


def source_composition_evidence() -> dict[str, Any]:
    """Return immutable recipe/equality evidence without evaluating a field."""

    recipe = finite_N_source_recipe()
    return {
        "source_only": True,
        "numeric_values_computed": False,
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "recipe": recipe,
        "equality": recipe["identity"],
        "function_definitions_distinct_from_values": True,
        "interval_enclosure_rows_are_not_promoted": True,
    }


class _ExactLaurent:
    """Small exact Laurent-polynomial ring used only by the internal check.

    The formal ``N`` variable may occur with exponent ``-1``.  This lets the
    check compare the C2 ``second/N`` expression with ``N*C5`` without
    assigning a numerical N or importing a symbolic algebra package.
    """

    __slots__ = ("terms",)
    _width = 32

    def __init__(self, terms: Mapping[tuple[int, ...], Fraction] | None = None) -> None:
        clean: dict[tuple[int, ...], Fraction] = {}
        for monomial, coefficient in (terms or {}).items():
            if coefficient:
                clean[tuple(monomial)] = Fraction(coefficient)
        self.terms = clean

    @classmethod
    def constant(cls, value: Any) -> "_ExactLaurent":
        return cls({(0,) * cls._width: Fraction(value)})

    @classmethod
    def variable(cls, index: int) -> "_ExactLaurent":
        powers = [0] * cls._width
        powers[index] = 1
        return cls({tuple(powers): Fraction(1)})

    def inverse(self) -> "_ExactLaurent":
        if len(self.terms) != 1:
            raise ValueError("exact Laurent check only inverts a monomial")
        monomial, coefficient = next(iter(self.terms.items()))
        if not coefficient:
            raise ZeroDivisionError("exact Laurent zero inverse")
        return _ExactLaurent({tuple(-power for power in monomial): Fraction(1, 1) / coefficient})

    def __add__(self, other: Any) -> "_ExactLaurent":
        other = _exact_coerce(other)
        result = dict(self.terms)
        for monomial, coefficient in other.terms.items():
            result[monomial] = result.get(monomial, Fraction(0)) + coefficient
        return _ExactLaurent(result)

    def __radd__(self, other: Any) -> "_ExactLaurent":
        return self + other

    def __neg__(self) -> "_ExactLaurent":
        return _ExactLaurent({monomial: -coefficient for monomial, coefficient in self.terms.items()})

    def __sub__(self, other: Any) -> "_ExactLaurent":
        return self + (-_exact_coerce(other))

    def __rsub__(self, other: Any) -> "_ExactLaurent":
        return _exact_coerce(other) - self

    def __mul__(self, other: Any) -> "_ExactLaurent":
        other = _exact_coerce(other)
        result: dict[tuple[int, ...], Fraction] = {}
        for left_monomial, left_coefficient in self.terms.items():
            for right_monomial, right_coefficient in other.terms.items():
                monomial = tuple(a + b for a, b in zip(left_monomial, right_monomial))
                result[monomial] = result.get(monomial, Fraction(0)) + left_coefficient * right_coefficient
        return _ExactLaurent(result)

    def __rmul__(self, other: Any) -> "_ExactLaurent":
        return self * other

    def __truediv__(self, other: Any) -> "_ExactLaurent":
        return self * _exact_coerce(other).inverse()

    def __eq__(self, other: Any) -> bool:
        return isinstance(other, _ExactLaurent) and self.terms == other.terms


def _exact_coerce(value: Any) -> _ExactLaurent:
    return value if isinstance(value, _ExactLaurent) else _ExactLaurent.constant(value)


def _generic_density_identity_check() -> dict[str, Any]:
    """Verify C2 = N*C5 in an exact Laurent-polynomial ordinary-jet ring."""

    exact = ScalarOperations(
        add=lambda a, b: a + b,
        sub=lambda a, b: a - b,
        mul=lambda a, b: a * b,
        div=lambda a, b: a / b,
        neg=lambda a: -a,
        constant=_ExactLaurent.constant,
        is_zero=lambda value: isinstance(value, _ExactLaurent) and value == _ExactLaurent.constant(0),
    )

    def jet(name: str, offset: int) -> OrdinaryJet:
        return OrdinaryJet(
            tuple(_ExactLaurent.variable(offset + index) for index in range(3)),
            provenance="symbolic",
            units=name,
        )

    E = jet("E", 0)
    V = jet("V", 3)
    r = jet("r", 6)
    s = jet("s", 9)
    N = OrdinaryJet(
        (_ExactLaurent.variable(12), _ExactLaurent.constant(0), _ExactLaurent.constant(0)),
        provenance="symbolic",
        units="N",
    )
    F = _jet_mul(N, r, exact)
    B = _jet_mul(N, s, exact)
    two = exact.constant(2)
    half = exact.div(exact.constant(1), exact.constant(2))
    first = {
        "m": B,
        "h": F,
        "k": _jet_add(_jet_mul(V, F, exact), _jet_mul(E, B, exact), exact),
        "e": _jet_sub(_jet_scale(_jet_mul(V, B, exact), two, exact), _jet_mul(E, F, exact), exact),
        "p": _jet_mul(E, F, exact),
    }
    zero = OrdinaryJet((_ExactLaurent.constant(0),) * 3, provenance="symbolic")
    second = {
        "m": zero,
        "h": zero,
        "k": _jet_mul(F, B, exact),
        "e": _jet_sub(_jet_mul(B, B, exact), _jet_scale(_jet_mul(F, F, exact), half, exact), exact),
        "p": _jet_scale(_jet_mul(F, F, exact), half, exact),
    }
    invN = exact.div(exact.constant(1), N[0])
    c2 = {name: _jet_add(first[name], _jet_scale(second[name], invN, exact), exact) for name in SIGNED_DENSITIES}
    c5 = {
        "m": s,
        "h": r,
        "k": _jet_add(_jet_add(_jet_mul(V, r, exact), _jet_mul(E, s, exact), exact), _jet_mul(r, s, exact), exact),
        "e": _jet_sub(
            _jet_add(_jet_scale(_jet_mul(V, s, exact), two, exact), _jet_mul(s, s, exact), exact),
            _jet_add(_jet_mul(E, r, exact), _jet_scale(_jet_mul(r, r, exact), half, exact), exact),
            exact,
        ),
        "p": _jet_add(_jet_mul(E, r, exact), _jet_scale(_jet_mul(r, r, exact), half, exact), exact),
    }
    scaled = {name: _jet_scale(c5[name], N[0], exact) for name in SIGNED_DENSITIES}
    rows = {
        name: all(c2[name][order] == scaled[name][order] for order in ORDINARY_Z_ORDERS)
        for name in SIGNED_DENSITIES
    }
    api_check = _generic_density_public_api_check()
    return {
        "proved": bool(all(rows.values()) and api_check.get("proved")),
        "signed_densities": rows,
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "method": "exact Laurent-polynomial ordinary-jet algebra through public composers",
        "assumptions": ["F=N*r", "B=N*s", "N_Z=N_ZZ=0"],
        "public_api": api_check,
    }


def _generic_density_public_api_check() -> dict[str, Any]:
    """Exercise both public composers in the exact formal Laurent ring."""

    exp_symbols: dict[tuple[Any, ...], int] = {}

    def exact_exp(value: Any) -> _ExactLaurent:
        if not isinstance(value, _ExactLaurent):
            raise TypeError("formal exponential expects an exact Laurent value")
        if value == _ExactLaurent.constant(0):
            return _ExactLaurent.constant(1)
        key = tuple(sorted(value.terms.items()))
        index = exp_symbols.setdefault(key, 13 + len(exp_symbols))
        return _ExactLaurent.variable(index)

    one = _ExactLaurent.constant(1)
    api_ops = ScalarOperations(
        add=lambda a, b: a + b,
        sub=lambda a, b: a - b,
        mul=lambda a, b: a * b,
        div=lambda a, b: a / b,
        neg=lambda a: -a,
        exp=exact_exp,
        expm1=lambda value: exact_exp(value) - one,
        constant=_ExactLaurent.constant,
        is_zero=lambda value: isinstance(value, _ExactLaurent) and value == _ExactLaurent.constant(0),
    )

    def jet(name: str, offset: int) -> OrdinaryJet:
        return OrdinaryJet(
            tuple(_ExactLaurent.variable(offset + index) for index in range(3)),
            provenance="symbolic",
            units=name,
        )

    E = jet("E", 0)
    V = jet("V", 3)
    A = jet("A", 6)
    s = jet("s", 9)
    N = OrdinaryJet(
        (_ExactLaurent.variable(12), _ExactLaurent.constant(0), _ExactLaurent.constant(0)),
        provenance="symbolic",
        units="N",
    )
    B = _jet_mul(N, s, api_ops)
    c5 = compose_finite_N_ordinary_jets(E, V, A, B, N, ops=api_ops)
    c2 = compose_c2_density_ordinary_jets(E, V, A, B, N, ops=api_ops)
    invN, _ = _inverse_N_jet(N, api_ops)
    # Independently written scalar derivative equations exercise the API;
    # they do not reuse the jet product, inverse or exponential helpers.
    n = N[0]
    x_expected = OrdinaryJet((A[0]/n, A[1]/n, A[2]/n), provenance="symbolic")
    X = exact_exp(x_expected[0])
    XZ = X*A[1]/n
    XZZ = X*(A[2]/n + A[1]*A[1]/(n*n))
    exp_expected = OrdinaryJet((X, XZ, XZZ), provenance="symbolic")
    expm1_expected = OrdinaryJet((X-one, XZ, XZZ), provenance="symbolic")
    r_expected = OrdinaryJet((E[0]*(X-one),
        E[1]*(X-one)+E[0]*XZ,
        E[2]*(X-one)+2*E[1]*XZ+E[0]*XZZ), provenance="symbolic")
    E_N_expected = OrdinaryJet((E[0]*X,
        E[1]*X+E[0]*XZ,
        E[2]*X+2*E[1]*XZ+E[0]*XZZ), provenance="symbolic")
    checks: dict[str, bool] = {
        "A_over_N_rows": all(c5["A_over_N"][order] == x_expected[order] for order in ORDINARY_Z_ORDERS),
        "r_rows": all(c5["r"][order] == r_expected[order] for order in ORDINARY_Z_ORDERS),
        "E_N_rows": all(c5["E_N"][order] == E_N_expected[order] for order in ORDINARY_Z_ORDERS),
        "F_equals_N_r": all(c2["F"][order] == _jet_mul(N, c5["r"], api_ops)[order] for order in ORDINARY_Z_ORDERS),
    }
    for name in SIGNED_DENSITIES:
        scaled = _jet_scale(c5["five_signed_density_increments"][name], N[0], api_ops)
        checks["density_" + name] = all(c2["density"][name][order] == scaled[order] for order in ORDINARY_Z_ORDERS)
    r2 = r_expected
    checks["factor_two_second_row"] = r2[2] == api_ops.add(
        api_ops.add(api_ops.mul(E[2], expm1_expected[0]), api_ops.mul(api_ops.constant(2), api_ops.mul(E[1], expm1_expected[1]))),
        api_ops.mul(E[0], expm1_expected[2]),
    )
    zero_A = OrdinaryJet((_ExactLaurent.constant(0), A[1], A[2]), provenance="symbolic")
    seam5 = compose_finite_N_ordinary_jets(E, V, zero_A, B, N, ops=api_ops)
    seam2 = compose_c2_density_ordinary_jets(E, V, zero_A, B, N, ops=api_ops)
    checks["A_zero_entire_extension"] = all(
        seam2["density"][name][order] == N[0]*seam5["five_signed_density_increments"][name][order]
        for name in SIGNED_DENSITIES for order in ORDINARY_Z_ORDERS)
    checks["constant_inverse_N_rows"] = invN.rows == (one/N[0], _ExactLaurent.constant(0), _ExactLaurent.constant(0))
    bad_N = OrdinaryJet(
        (N[0], _ExactLaurent.variable(14), _ExactLaurent.constant(0)),
        provenance="symbolic",
        units="N",
    )
    try:
        compose_finite_N_ordinary_jets(E, V, A, B, bad_N, ops=api_ops)
    except ValueError as error:
        checks["reject_nonconstant_N"] = str(error) == "N must have exact zero ordinary Z and Z2 rows"
    except Exception:
        checks["reject_nonconstant_N"] = False
    else:
        checks["reject_nonconstant_N"] = False
    return {
        "proved": bool(all(checks.values())),
        "checks": checks,
        "public_composers_called": ["compose_finite_N_ordinary_jets", "compose_c2_density_ordinary_jets"],
        "formal_exp_relation": "expm1(x)=exp(x)-1",
    }


def _digest(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False).encode("utf-8")
    ).hexdigest()


def canonical_digest(value: Any) -> str:
    """Canonical graph digest helper exposed for runners and static checks."""

    return _digest(value)


def _module_hash() -> str | None:
    try:
        return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    except OSError:
        return None


def _gate(gates: list[str], name: str, *, chart: str | None = None) -> None:
    value = name if chart is None else f"{name}:{chart}"
    if value not in gates:
        gates.append(value)


def _record(source_records: Mapping[str, Any], name: str, gates: list[str], *, optional: bool = False) -> str | None:
    record = source_records.get(name) if isinstance(source_records, Mapping) else None
    if not isinstance(record, Mapping):
        if not optional:
            _gate(gates, "source_record_missing")
        return None
    expected_path = "experiments/root_st073/" + name
    text = record.get("source_text")
    if record.get("path") != expected_path:
        _gate(gates, "source_path_mismatch")
        return None
    if not isinstance(text, str):
        _gate(gates, "source_text_missing")
        return None
    expected_hash = RECORD_HASHES[name]
    if record.get("sha256") != expected_hash or hashlib.sha256(text.encode("utf-8")).hexdigest() != expected_hash:
        _gate(gates, "source_hash_mismatch")
        return None
    return text


def _parse_json(text: str | None, gates: list[str]) -> Mapping[str, Any]:
    if not isinstance(text, str):
        return {}
    try:
        value = json.loads(text)
    except (TypeError, ValueError):
        _gate(gates, "source_json_invalid")
        return {}
    return value if isinstance(value, Mapping) else {}


def _tree(text: str | None, gates: list[str]) -> ast.AST:
    if not isinstance(text, str):
        return ast.Module(body=[], type_ignores=[])
    try:
        return ast.parse(text)
    except (SyntaxError, ValueError):
        _gate(gates, "source_ast_invalid")
        return ast.Module(body=[], type_ignores=[])


def _fn(tree: ast.AST, name: str, cls: str | None = None) -> ast.AST | None:
    if cls is not None:
        for node in ast.walk(tree):
            if isinstance(node, ast.ClassDef) and node.name == cls:
                return next((child for child in node.body if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef)) and child.name == name), None)
        return None
    return next((node for node in ast.walk(tree) if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == name), None)


def _norm(text: str) -> str:
    return "".join(text.split())


def _has(node: ast.AST | None, *needles: str) -> bool:
    if node is None:
        return False
    try:
        body = _norm(ast.unparse(node))
    except (AttributeError, TypeError):
        return False
    return all(_norm(needle) in body for needle in needles)


def _source_has(text: str | None, *needles: str) -> bool:
    return isinstance(text, str) and all(_norm(needle) in _norm(text) for needle in needles)


def _target_text(target: ast.AST) -> str:
    try:
        return _norm(ast.unparse(target))
    except (AttributeError, TypeError):
        return ""


def _assign_values(node: ast.AST | None, target: str) -> list[ast.AST]:
    if node is None:
        return []
    wanted = _norm(target)
    values: list[ast.AST] = []
    for item in ast.walk(node):
        if isinstance(item, ast.Assign) and any(_target_text(t) == wanted for t in item.targets):
            values.append(item.value)
        elif isinstance(item, ast.AnnAssign) and _target_text(item.target) == wanted:
            values.append(item.value)
    return values


def _same_ast(value: ast.AST | None, expression: str) -> bool:
    if value is None:
        return False
    try:
        wanted = ast.parse(expression, mode="eval").body
    except (SyntaxError, ValueError):
        return False
    return ast.dump(value, include_attributes=False) == ast.dump(wanted, include_attributes=False)


def _assign(node: ast.AST | None, target: str, expression: str) -> bool:
    return any(_same_ast(value, expression) for value in _assign_values(node, target))


def _dict_field(node: ast.AST | None, variable: str, field: str) -> ast.AST | None:
    for value in _assign_values(node, variable):
        if isinstance(value, ast.Dict):
            for key, result in zip(value.keys, value.values):
                if isinstance(key, ast.Constant) and key.value == field:
                    return result
        if isinstance(value, ast.Call) and isinstance(value.func, ast.Name) and value.func.id == "dict":
            for keyword in value.keywords:
                if keyword.arg == field:
                    return keyword.value
    return None


def _receipt(text: str | None, required: Sequence[Sequence[str]], gates: list[str]) -> bool:
    value = _parse_json(text, gates)
    ok = value.get("all_passed") is True
    for path in required:
        current: Any = value
        for key in path:
            current = current.get(key) if isinstance(current, Mapping) else None
        if current is not True:
            ok = False
    family = value.get("source_family")
    if isinstance(family, Mapping):
        ok = ok and all(family.get(k) == v for k, v in SOURCE_FAMILY.items())
    elif any(k in value for k in SOURCE_FAMILY):
        ok = ok and all(value.get(k) == v for k, v in SOURCE_FAMILY.items())
    return bool(ok)


def _report_witness(texts: Mapping[str, str], gates: list[str]) -> dict[str, Any]:
    parent = _parse_json(texts.get(PARENT_REPORT), gates)
    latest = _parse_json(texts.get(LATEST_REPORT), gates)
    parent_ok = bool(parent) and parent.get("source_only") is True and parent.get("numeric_values_computed") is False
    component = latest.get("component_result") if isinstance(latest, Mapping) else None
    if not isinstance(component, Mapping):
        component = {}
    latest_ok = bool(component) and all(
        (
            component.get("source_only") is True,
            component.get("passed") is True,
            component.get("full_R110_Rm_background_projection_proved") is True,
            component.get("full_R110_Rm_own_history_projection_proved") is True,
            component.get("ordinary_Z_orders") == [0, 1, 2],
            component.get("endpoint_ports", {}).get("endpoint_rows_are_values") is False,
            component.get("numeric_values_computed") is False,
            component.get("old_operators_replayed") is False,
        )
    )
    if not parent_ok:
        _gate(gates, "accepted_parent_source_projection_unproved")
    if not latest_ok:
        _gate(gates, "accepted_R110_Rm_source_projection_unproved")
    return {"parent": parent, "latest": latest, "component": component, "passed": parent_ok and latest_ok}


def _formula_witness(trees: Mapping[str, ast.AST], texts: Mapping[str, str], gates: list[str]) -> dict[str, Any]:
    c2_density = trees.get(PREFIX + "C2_density_transport.py", ast.Module(body=[], type_ignores=[]))
    c2_phase = trees.get(PREFIX + "C2_phase_primitives.py", ast.Module(body=[], type_ignores=[]))
    c5 = trees.get(PREFIX + "prescribed_loop_C5_functions.py", ast.Module(body=[], type_ignores=[]))
    c2_fn = _fn(c2_density, "coefficient_C2")
    c2_build = _fn(c2_density, "build")
    c2_primitive = _fn(c2_phase, "primitive_C2")
    c5_finite = _fn(c5, "finite_N_functions")
    c5_roots = _fn(c5, "generic_roots")
    c5_loop = _fn(c5, "loop_functions")

    c2_first = {
        "m": _target_text(_dict_field(c2_fn, "first", "m")) == "B",
        "h": _target_text(_dict_field(c2_fn, "first", "h")) == "F",
        "k": _same_ast(_dict_field(c2_fn, "first", "k"), "alg.add(VF, EB)"),
        "e": _same_ast(_dict_field(c2_fn, "first", "e"), "alg.add(alg.scale(VB, 2), alg.neg(EF))"),
        "p": _same_ast(_dict_field(c2_fn, "first", "p"), "EF"),
    }
    c2_second = {
        "m": _target_text(_dict_field(c2_fn, "second", "m")) == "zero",
        "h": _target_text(_dict_field(c2_fn, "second", "h")) == "zero",
        "k": _same_ast(_dict_field(c2_fn, "second", "k"), "alg.mul(F, B)"),
        "e": _same_ast(_dict_field(c2_fn, "second", "e"), "alg.add(B2, alg.neg(alg.scale(F2, '1/2')))"),
        "p": _same_ast(_dict_field(c2_fn, "second", "p"), "alg.scale(F2, '1/2')"),
    }
    c2_ok = all(c2_first.values()) and all(c2_second.values()) and _has(
        c2_fn,
        "arg=g.mul(A.value,invN)",
        "exprel=g.unary('exprel',arg)",
        "exponential=g.unary('exp',arg)",
        "FZZ=g.add",
    ) and _has(c2_build, "coefficients[-1][key]", "coefficients[-2][key]", "density={key:alg.add")

    c5_expected = {
        "invN": _assign(c5_finite, "invN", "g.quotient(g.one, N, 'same exact positive integer N on every source chart')"),
        "s": _assign(c5_finite, "s", "a.scale(loop['B_over_Pstar'], invN)"),
        "r": _assign(c5_finite, "r", "a.mul(roots['E'], exm1)"),
        "E_N": _assign(c5_finite, "E_N", "a.mul(E, ex)"),
        "V_N": _assign(c5_finite, "V_N", "a.add(V, s)"),
        "increments": _assign(
            c5_finite,
            "inc",
            "dict(m=s, h=r, k=a.add(a.mul(V, r), a.mul(E, s), a.mul(r, s)), e=a.sub(a.add(a.scale(a.mul(V, s), 2), a.mul(s, s)), a.add(a.mul(E, r), a.scale(a.mul(r, r), '1/2'))), p=a.add(a.mul(E, r), a.scale(a.mul(r, r), '1/2')))",
        ),
        "full_exp_rows": _has(c5_finite, "ex=a.exp(x)", "FunctionJet((g.unary('expm1',x.value),)+rows(ex)[1:])"),
    }
    c5_root_expected = {
        "C": _assign(c5_roots, "C", "a.sub(E, a.scale(a.trim(radial['E_y']), 2))"),
        "B": _assign(c5_roots, "B", "a.scale(a.trim(radial['V_y']), 2)"),
        "aa_bb": _assign(c5_roots, "(aa, bb)", "(divide(C), divide(B))"),
        "t0": _assign(c5_roots, "t0", "a.neg(a.div(B, C, 'same original a>0 and C=a*E'))"),
        "p2": _has(c5_roots, "p2=generic['p2']"),
    }
    c5_loop_expected = {
        "B": _assign(c5_loop, "B", "lazy(a.scale(a.mul(roots['E'],M),'1/2'))"),
        "A": _has(c5_loop, "A=lazy(a.scale(a.mul(roots['a'],chi),'1/2'))"),
        "inverse": _has(c5_loop, "original_monotone_phase_inverse", "true_prescribed_source_C5_implicit_inverse=True", "phase_endpoints=[0,1]"),
    }
    c5_ok = all(c5_expected.values()) and all(c5_root_expected.values()) and all(c5_loop_expected.values())
    c2_phase_ok = _has(
        c2_primitive,
        "same_original_inverse_not_a_new_branch=True",
        "phase_Z_and_ZZ_exact_zero=True",
        "BZZ=flat(g.mul(g.constant('1/2'),g.add(g.mul(E.ZZ,M),g.mul(g.constant(2),E.Z,MZ),g.mul(E.value,MZZ))))",
    )
    if not c2_ok:
        _gate(gates, "c2_full_density_formula_ast_unproved")
    if not c5_ok:
        _gate(gates, "c5_finite_N_field_formula_ast_unproved")
    if not c2_phase_ok:
        _gate(gates, "c2_common_phase_inverse_ast_unproved")
    return {
        "c2": {"passed": c2_ok, "first": c2_first, "second": c2_second, "density": "first+second/N"},
        "c5": {"passed": c5_ok, "finite_N": c5_expected, "roots": c5_root_expected, "loop": c5_loop_expected},
        "c2_phase": {"passed": c2_phase_ok},
        "passed": c2_ok and c5_ok and c2_phase_ok,
    }


def _owner_witness(trees: Mapping[str, ast.AST], texts: Mapping[str, str], gates: list[str]) -> dict[str, Any]:
    c5 = texts.get(PREFIX + "prescribed_loop_C5_functions.py", "")
    bridge = texts.get(PREFIX + "prescribed_bridge_C6_functions.py", "")
    r100 = texts.get(PREFIX + "prescribed_R100_Rm_C6_functions.py", "")
    c2_phase = texts.get(PREFIX + "C2_phase_primitives.py", "")
    c5_tokens = (
        "type(background) is not previous.CurrentOriginalPrescribedR100RmC6Functions",
        "self.background=background",
        "self.Pstar=background.Pstar;self.P0=background.P0",
        "self.source_charts={",
        "'first_micro':bridge.micro",
        "'second_micro':bridge.micro",
        "'frozen_macro':bridge.macro",
        "'first_switch':background.charts['first_switch']",
        "'second_switch':background.charts['second_switch']",
        "'post_power':background.charts['post_switch']",
        "'long_reshape':background.charts['long_reshape']",
        "'reference':background.charts['reference']",
        "'restoration':background.charts['restoration']",
        "'postrestore':background.charts['post_restore']",
        "roots=generic_roots(g,chart,self.Pstar)",
        "phase=g.unary('fractional_part',g.mul(self.N,offset))",
    )
    bridge_tokens = (
        "class CurrentOriginalPrescribedBridgeC6Functions",
        "self.P0=self.pressure7()",
        "fields,radial,generic=generic_fields",
        "fields=dict(E=",
        "m=a.scale(J['M'],inverse)",
        "P0=a.trim(P0)",
        "true_prescribed_background_before_finite_N_modification=True",
    )
    r100_tokens = (
        "type(bridge_owner) is not previous.CurrentOriginalPrescribedBridgeC6Functions",
        "self.bridge=bridge_owner",
        "self.P0=bridge_owner.P0",
        "self.Pstar=Pstar=bridge_owner.inlet.Pstar",
        "self.add_chart('long_reshape'",
        "self.add_chart('reference'",
        "self.add_chart('restoration'",
        "self.add_chart('post_restore'",
        "original_independent_P0_and_physical_Pstar_F0_factors_once=True",
    )
    source_packet_tokens = (
        "packet['source_identity']!=self.identity",
        "packet['exact_common_P0_axial5'] is not op.P0",
        "source_radius_and_phase_Z_ZZ_independent=True",
        "same_current_Z_Taylor_coordinate=True",
    )
    passed = all(_source_has(text, *tokens) for text, tokens in ((c5, c5_tokens), (bridge, bridge_tokens), (r100, r100_tokens), (c2_phase, source_packet_tokens)))
    if not passed:
        _gate(gates, "actual_C2_provider_C5_C6_owner_chain_unproved")
    return {
        "passed": passed,
        "c2_provider": "CurrentC2PhasePrimitives.source_packet",
        "c5_provider": "CurrentOriginalPrescribedLoopC5Functions.source_charts -> generic_roots -> finite_N_functions",
        "c6_provider": "CurrentOriginalPrescribedR100RmC6Functions -> CurrentOriginalPrescribedBridgeC6Functions",
        "same_Pstar_source": _source_has(c5, "self.Pstar=background.Pstar;self.P0=background.P0") and _source_has(r100, "self.Pstar=Pstar=bridge_owner.inlet.Pstar"),
        "same_phase_source": _source_has(c5, "phase=g.unary('fractional_part',g.mul(self.N,offset))") and _source_has(c2_phase, "same_current_Z_Taylor_coordinate=True"),
        "same_N_source": _source_has(c5, "one_common_N_for_all_charts_and_control_denominators=True"),
        "source_packet_identity": _source_has(c2_phase, "packet['source_identity']!=self.identity", "packet['exact_common_P0_axial5'] is not op.P0"),
    }








def _provider_admission(source_records: Mapping[str, Any], gates: list[str]) -> dict[str, Any]:
    """Keep live C2/C5/C6 admission closed until a source theorem binds handles.

    Caller-supplied booleans, names, hashes or enclosure-overlap records are
    descriptive only at this stage.  They cannot establish equality of the
    actual function handles, so this function intentionally has no admitting
    path.  The exact recipe and independent generic algebra check remain
    available to the later theorem.
    """

    _gate(gates, "actual_C2_provider_C5_C6_field_function_identity_unbound")
    return {"present": False, "key": "deferred_source_theorem", "rows": {}, "passed": False}




def _receipt_witness(texts: Mapping[str, str], gates: list[str]) -> dict[str, bool]:
    requirements = {
        PREFIX + "C2_phase_primitives_check.json": (
            ("current_original_17_chart_source_Z2_rows_installed",),
            ("current_original_same_inverse_C2_phase_primitive_functions_installed",),
            ("upstream_six_chart_C2_provider_installed",),
        ),
        PREFIX + "C2_density_transport_check.json": (
            ("current_original_17_chart_C2_signed_density_functions_installed",),
            ("current_original_17_chart_C2_own_rate_transport_and_five_target_functions_installed",),
            ("same_actual_terminal_amplitude_P0_source_and_global_phase_retained",),
        ),
        PREFIX + "prescribed_loop_C5_functions_check.json": (
            ("checks", "actual_source", "passed"),
            ("checks", "actual_source", "genuine_C6_inputs_and_C5_generic_roots"),
            ("checks", "actual_source", "same_original_N_eta_dstar_and_global_phase"),
            ("checks", "actual_source", "all_ten_correction_predecessors_preserved"),
            ("current_original_true_prescribed_ten_chart_signed_loop_C5_functions_installed",),
            ("current_original_true_prescribed_ten_chart_finite_N_density_C5_functions_installed",),
        ),
        PREFIX + "prescribed_bridge_C6_functions_check.json": (
            ("current_original_true_comparison_bridge_C7_functions_installed",),
            ("current_original_true_prescribed_micro_macro_C6_functions_installed",),
            ("current_original_true_prescribed_bridge_Q_inertial_p2_C5_functions_installed",),
        ),
        PREFIX + "prescribed_R100_Rm_C6_functions_check.json": (
            ("current_original_true_prescribed_R100_R110_C6_functions_installed",),
            ("current_original_true_prescribed_finite_long_reference_restore_Rm_C6_functions_installed",),
            ("current_original_true_prescribed_R100_Rm_Q_inertial_p2_C5_functions_installed",),
        ),
        PREFIX + "whole_Z_bridge_source_check.json": (("current_original_whole_Z_live_leading_bridge_source_functions_installed",),),
        PREFIX + "whole_Z_bridge_correlated_functions_check.json": (("current_original_whole_Z_correlated_bridge_C1_functions_and_Rc_retransport_installed",),),
        PREFIX + "whole_Z_all_N_source_functions_check.json": (("current_original_whole_Z_bridge_full_phase_all_N_Nscaled_C1_functions_and_R100_transport_installed",),),
        PREFIX + "whole_Z_all_N_switch_functions_check.json": (("current_original_whole_Z_full_phase_all_N_switch_C1_functions_and_normalized_R110_transport_installed",),),
        PREFIX + "whole_Z_all_N_long_patch_functions_check.json": (("current_original_whole_Z_full_phase_all_N_long_patch_C1_functions_and_normalized_Rh_transport_installed",),),
    }
    result = {name: _receipt(texts.get(name), paths, gates) for name, paths in requirements.items()}
    return result


def _source_witness(source_records: Mapping[str, Any], gates: list[str]) -> dict[str, Any]:
    texts: dict[str, str] = {}
    for name in RECORD_HASHES:
        text = _record(source_records, name, gates, optional=name in OPTIONAL_PATCH_RECORDS)
        if text is not None:
            texts[name] = text
    trees = {name: _tree(text, gates) for name, text in texts.items() if name.endswith(".py")}
    reports = _report_witness(texts, gates)
    receipts = _receipt_witness(texts, gates)
    formulas = _formula_witness(trees, texts, gates)
    owners = _owner_witness(trees, texts, gates)
    receipt_ok = bool(receipts) and all(receipts.values())
    if not receipt_ok:
        _gate(gates, "accepted_source_receipts_unproved")
    return {
        "texts": texts,
        "trees": trees,
        "reports": reports,
        "receipts": receipts,
        "formulas": formulas,
        "owners": owners,
        "receipt_ok": receipt_ok,
        "passed": reports["passed"] and formulas["passed"] and owners["passed"] and receipt_ok,
    }


def _endpoint_ports(component: Mapping[str, Any]) -> dict[str, Any]:
    ports = component.get("endpoint_ports") if isinstance(component, Mapping) else None
    if not isinstance(ports, Mapping):
        return {"endpoint_rows_are_values": False}
    value = dict(ports)
    value["endpoint_rows_are_values"] = False
    return value


def _chart_row(chart: str, admitted: bool, *, provider: Any, reason: str | None = None) -> dict[str, Any]:
    row = {
        "chart": chart,
        "status": "proved" if admitted else "unresolved",
        "field_source_composition_proved": admitted,
        "density_C2_NC5_identity_proved": admitted,
        "same_live_C2_C5_C6_provider": admitted,
        "source_provider": provider,
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "signed_densities": list(SIGNED_DENSITIES),
        "source_only": True,
        "numeric_values_computed": False,
        "endpoint_rows_are_values": False,
        "history_identity_proved": False,
    }
    if reason is not None:
        row["unresolved_gates"] = [reason]
    return row


def _safe_digest(nodes: Any, gates: list[str]) -> str:
    try:
        return _digest(nodes)
    except (TypeError, ValueError, OverflowError):
        _gate(gates, "nodes_digest_error")
        return ""


def derive(nodes: Sequence[Any], source_records: Mapping[str, Any], *, expected_nodes_digest: str) -> dict[str, Any]:
    """Return a fail-closed source composition result.

    The graph is hashed but no graph node is evaluated.  Missing records,
    malformed maps and missing provider fields produce unresolved chart rows;
    they never raise a ``TypeError`` or silently widen admission.
    """

    gates: list[str] = []
    actual_digest = _safe_digest(nodes, gates)
    if expected_nodes_digest != EXPECTED_NODES_DIGEST:
        _gate(gates, "expected_nodes_digest_mismatch")
    if actual_digest != EXPECTED_NODES_DIGEST:
        _gate(gates, "nodes_digest_mismatch")
    if not isinstance(nodes, Sequence) or isinstance(nodes, (str, bytes, bytearray)) or len(nodes) != EXPECTED_NODE_COUNT:
        _gate(gates, "node_count_mismatch")
    records = source_records if isinstance(source_records, Mapping) else {}
    if not isinstance(source_records, Mapping):
        _gate(gates, "source_records_not_mapping")
    try:
        witness = _source_witness(records, gates)
    except Exception as error:  # fail closed without exposing a provider TypeError
        _gate(gates, "source_inspection_error:" + type(error).__name__)
        witness = {"passed": False, "reports": {"component": {}}, "formulas": {}, "owners": {}, "receipts": {}}
    try:
        generic_identity = _generic_density_identity_check()
    except Exception as error:  # exact algebra is an independent hard gate
        _gate(gates, "generic_density_identity_check_error:" + type(error).__name__)
        generic_identity = {
            "proved": False,
            "signed_densities": {},
            "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
            "method": "exact Laurent-polynomial ordinary-jet algebra",
        }
    providers = _provider_admission(records, gates)
    base_ok = bool(
        witness.get("passed")
        and generic_identity.get("proved")
        and providers.get("present")
        and providers.get("passed")
    )
    rows: list[dict[str, Any]] = []
    for chart in FIRST_TEN_CHARTS:
        admitted = bool(base_ok and providers.get("rows", {}).get(chart))
        rows.append(
            _chart_row(
                chart,
                admitted,
                provider={
                    "map_key": providers.get("key"),
                    "binding_present": providers.get("rows", {}).get(chart, False),
                },
                reason=None if admitted else "first_ten_live_C2_C5_C6_provider_binding",
            )
        )
    patch_admitted = bool(base_ok and providers.get("rows", {}).get("actual_patch"))
    if patch_admitted:
        rows.append(_chart_row("actual_patch", True, provider={"map_key": providers.get("key"), "binding_present": True}))
    unresolved_outer = [chart for chart in OUTER_SIX_CHARTS if chart not in {row["chart"] for row in rows}]
    first_ten_field = bool(rows) and all(row["field_source_composition_proved"] for row in rows if row["chart"] in FIRST_TEN_CHARTS)
    first_ten_density = bool(rows) and all(row["density_C2_NC5_identity_proved"] for row in rows if row["chart"] in FIRST_TEN_CHARTS)
    if not first_ten_field:
        _gate(gates, "first_ten_finite_N_field_source_composition_unproved")
    if not first_ten_density:
        _gate(gates, "first_ten_density_C2_NC5_identity_unproved")
    component = witness.get("reports", {}).get("component", {}) if isinstance(witness.get("reports"), Mapping) else {}
    endpoint_ports = _endpoint_ports(component)
    recipe = finite_N_source_recipe()
    result = {
        "schema": SCHEMA,
        "source_only": True,
        "passed": bool(first_ten_field and first_ten_density and not gates),
        "nodes_digest": actual_digest,
        "expected_nodes_digest": EXPECTED_NODES_DIGEST,
        "node_count": len(nodes) if isinstance(nodes, Sequence) and not isinstance(nodes, (str, bytes, bytearray)) else None,
        "covered_charts": [row["chart"] for row in rows if row["status"] == "proved"],
        "rows": rows,
        "unresolved_rows": [row for row in rows if row["status"] != "proved"],
        "first_ten_finite_N_field_source_composition_proved": first_ten_field,
        "first_ten_density_C2_NC5_identity_proved": first_ten_density,
        "full_finite_N_field_source_composition_proved": False,
        "generic_C2_NC5_coefficientwise_identity": {
            "relation": "C2=N*C5",
            "signed_densities": list(SIGNED_DENSITIES),
            "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
            "proved_as_formula": bool(witness.get("formulas", {}).get("passed")),
            "generic_density_identity_proved": bool(generic_identity.get("proved")),
            "requires_same_B_provider": True,
            "requires_same_Pstar_units": True,
            "F_equals_N_r_entire_at_A_zero": True,
        },
        "generic_density_identity_proved": bool(generic_identity.get("proved")),
        "generic_density_identity_check": generic_identity,
        "patch_source_composition_proved": False,
        "finite_N_source_recipe": recipe,
        "source_composition_evidence": source_composition_evidence(),
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "ordinary_Z_not_temporal_n": True,
        "endpoint_ports": endpoint_ports,
        "numeric_values_computed": False,
        "old_operators_replayed": False,
        "mathematical_operator_replayed": False,
        "finite_N_history_identity": False,
        "same_density_certificate": False,
        "physical_Q_y_Z1": False,
        "global_five_moment_closure": False,
        "time_order_recursion_claimed": False,
        "timeorder_recursion_claimed": False,
        "full_goal_complete": False,
        "ast_witnesses": {
            "accepted_parent_and_R110_Rm_projection": witness.get("reports", {}).get("passed") is True,
            "finite_N_and_C2_density_formulas": witness.get("formulas", {}).get("passed") is True,
            "same_C2_phase_inverse": witness.get("formulas", {}).get("c2_phase", {}).get("passed") is True,
            "actual_C2_provider_C5_C6_owner_chain": witness.get("owners", {}).get("passed") is True,
            "accepted_source_receipts": witness.get("receipt_ok") is True,
            "provider_map_rows": providers.get("passed") is True,
        },
        "accepted_source_reports": {
            "parent": PARENT_REPORT,
            "parent_sha256": PARENT_SHA256,
            "latest": LATEST_REPORT,
            "latest_sha256": LATEST_SHA256,
        },
        "source_provider_witness": {
            "map_key": providers.get("key"),
            "present": providers.get("present", False),
            "first_ten_rows": providers.get("rows", {}),
            "raw_map": providers.get("raw"),
            "owner_chain": witness.get("owners", {}),
        },
        "source_inspection_audit": {
            "source_records_only": True,
            "static_source_AST_graph_and_saved_bound_inspection_only": True,
            "nodes_hashed_only": True,
            "new_source_consumer_calls": 1,
            "accepted_ancestor_operator_calls": 0,
            "producer_provider_builder_checker_replayed": False,
            "quadrature_or_norm_calls": 0,
            "numeric_field_evaluated": False,
            "recurrence_calls": 0,
            "target_calls": 0,
        },
        "record_hashes": dict(RECORD_HASHES),
        "module_sha256": _module_hash(),
        "runtime_model_metadata": None,
        "unresolved_gates": sorted(set(gates)),
        "scope_limits": {
            "outer_six_unresolved": unresolved_outer,
            "global_same_density_certificate": False,
            "history_requires_same_rate_J_partitions_and_nonreset_incoming": True,
            "endpoint_rows_are_values": False,
        },
    }
    return result


__all__ = [
    "EXPECTED_NODE_COUNT",
    "EXPECTED_NODES_DIGEST",
    "RECORD_HASHES",
    "OrdinaryJet",
    "ScalarOperations",
    "canonical_digest",
    "compose_c2_density_ordinary_jets",
    "compose_finite_N_ordinary_jets",
    "derive",
    "finite_N_source_recipe",
    "source_composition_evidence",
]
