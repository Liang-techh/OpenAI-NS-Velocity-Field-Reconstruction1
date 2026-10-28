"""Consumer-only direct interval matrix enclosure draft, v2.

The root agent must construct one frozen JSON bundle before calling
``consume_frozen_source_bundle``.  This module does not read a result file,
import an old operator or producer, run a checker, select a midpoint, or
materialize ``F0``.  The bundle contains only the selected actual interval
records as exact mpmath endpoint tuples; complete upstream records and their
old errors are carried through an opaque ``upstream`` field.

The two source leaves and exact extraction paths are deliberately explicit:

* A: ``tmp/native_O3_eta_n1_core_source/result.json`` / ``result.payload``;
  ``parameters.delta``, ``parameters.epsilon_core``, ``parameters.j``;
  ``mixed_jet_table.DU.uniform_eta.uniform_rho`` keys ``m0_k0``,
  ``m0_k1``, ``m1_k0``; and ``gauge_table.uniform_eta.uniform_rho`` keys
  ``z`` and ``F_jets.m0_k0/m0_k1/m1_k0``.
* B: ``tmp/native_O3_anchor_second_Z_eta_chart/result.json`` /
  ``result.payload.compatible_pressure_Q.uniform_eta.rows.uniform_rho.Q_chart``.

The output is a conditional source enclosure.  Hashes and common source ids
are recorded and cross-checked as caller-owned identity premises, never
treated as a proof.  The A exact-view registration and B direct eta-chart
proof are retained as conditional premises.  No claim is made about
source-function point values,
matrix C1, W or old-212 errors, or global closure.

The root-frozen identity contract additionally requires
``A.source_ids.eta_result_sha256 == B.result_sha256`` and the literal
``e77879ef7c7175ca2d042719ba875c3125c508c952c08366f98e7c16b86c6f76``.
Selected interval records contain only exact endpoint tuples; decimal
``lower``/``upper``/``width`` fields and midpoint/value fields are rejected.
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from typing import Any, Mapping

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from mpmath.libmp import finf, fnan, fninf, mpf_cmp


SCHEMA = "native_O3_eta_n1_matrix_enclosure_draft_v2"
RUNTIME_MODEL_METADATA = None
CONDITIONAL = True
HASH_IS_PROOF = False
F0_ANCHOR = "F0_anchor"
POINT_VALUE = False
INTERVAL_ENGINE = "mpmath.iv (MPIntervalContext)"

SOURCE_PATHS = {
    "A_result": (
        "C:/Users/z5242/Documents/Codex/2026-09-13/"
        "before-doing-any-project-work-verify/work/st073-transition-next/"
        "tmp/native_O3_eta_n1_core_source/result.json"
    ),
    "B_result": (
        "C:/Users/z5242/Documents/Codex/2026-09-13/"
        "before-doing-any-project-work-verify/work/st073-transition-next/"
        "tmp/native_O3_anchor_second_Z_eta_chart/result.json"
    ),
}
COMPONENTS = {
    "A": "native_O3_eta_n1_core_source",
    "B": "native_O3_anchor_second_Z_eta_chart",
}
EXPECTED_ETA_RESULT_SHA256 = (
    "e77879ef7c7175ca2d042719ba875c3125c508c952c08366f98e7c16b86c6f76"
)
ZERO_MPF = (0, 0, 0, 0)
NONFINITE_MPF_TUPLES = frozenset((fnan, finf, fninf))

FIELDS = ("F1", "Uz1", "K1", "P1", "d_x_F1", "d_x_Uz1")
DIAGONAL = (0, 0, 2, 0, 3, 1)
B0_KINEMATIC = {"1,5": 1, "2,6": 1, "3,6": -1}
F_DEPENDENT_B0 = frozenset(("4,1", "5,2", "5,3", "6,1"))
F_DEPENDENT_B1 = frozenset(("5,2", "5,3"))

# This is the raw six-row matrix, not a call to the saved operator.  The
# expressions are transcribed as algebra below and retained in metadata.
RAW_FORMULAS = {
    "equation": "d_x W + D W/x = B0 W + B1 d_Z W + g",
    "x": "sqrt(rho)",
    "H": "Z*(1-delta)/2 + d*U",
    "B": "Q + (1-2*Z*U)/L",
    "E": "Z*U - 1/2",
    "Cf": "F + rho*F_r",
    "Zf": "(F*(-2-delta)*Z + d*F_Z - 2*Z*rho*F_r)/L",
    "Zw": "(U*(-1-delta)*Z + d*U_Z - 2*Z*rho*U_r)/L",
    "B0_4_1": "epsilon_core*(4*x*F)",
    "B0_5_1": "epsilon_core*2*(Q + E*(-2+delta)/L)",
    "B0_5_2": "epsilon_core*(2*Zf - 2*Cf*Z*(-1+delta)/L)",
    "B0_5_3": "epsilon_core*(-2*Cf*Z*(1+delta)/L)",
    "B0_5_5": "epsilon_core*B*x",
    "B0_6_1": "epsilon_core**2*(-8*Z*rho*F/L)",
    "B0_6_2": "epsilon_core*2*(E*(-1+delta)/L + Zw - U_r*Z*(-1+delta)*rho/L)",
    "B0_6_3": "epsilon_core*(-2*U_r*Z*(1+delta)*rho/L)",
    "B0_6_4": "epsilon_core*(-4*Z/L)",
    "B0_6_6": "epsilon_core*B*x",
    "B1_5_1": "epsilon_core*(2*H/L)",
    "B1_5_2": "epsilon_core*(-2*Cf*d/L)",
    "B1_5_3": "epsilon_core*(-2*Cf*d/L)",
    "B1_6_2": "epsilon_core*2*(H-U_r*d*rho)/L",
    "B1_6_3": "epsilon_core*(-2*U_r*d*rho/L)",
    "B1_6_4": "epsilon_core*(2*d/L)",
}


class MatrixEnclosureInputError(ValueError):
    """Raised when the root-frozen source contract is absent or malformed."""


@dataclass(frozen=True)
class FrozenInterval:
    """An input interval represented only by exact mpmath endpoint tuples."""

    lower_exact_mpf_tuple: tuple[int, int, int, int]
    upper_exact_mpf_tuple: tuple[int, int, int, int]
    point_value: bool = False

    @classmethod
    def from_json(cls, record: Mapping[str, Any], label: str) -> "FrozenInterval":
        if not isinstance(record, Mapping):
            raise MatrixEnclosureInputError(label + " must be a JSON object")
        if record.get("point_value") is not False:
            raise MatrixEnclosureInputError(label + " must have point_value=false")
        for forbidden in ("lower", "upper", "width", "value", "midpoint"):
            if forbidden in record:
                raise MatrixEnclosureInputError(
                    label + " may carry only exact endpoint tuples, not " + forbidden
                )
        lower = _mpf_tuple(record.get("lower_exact_mpf_tuple"), label + ".lower")
        upper = _mpf_tuple(record.get("upper_exact_mpf_tuple"), label + ".upper")
        return cls(lower, upper, False)

    def interval(self, iv: MPIntervalContext, scalar: Any) -> Any:
        # ``iv.make_mpf`` accepts the exact backend endpoint tuples directly.
        # Converting through a scalar context first can round an endpoint
        # inward, so the scalar argument is intentionally unused here.
        del scalar
        if mpf_cmp(self.lower_exact_mpf_tuple, self.upper_exact_mpf_tuple) > 0:
            raise MatrixEnclosureInputError("reversed exact endpoints")
        return iv.make_mpf((self.lower_exact_mpf_tuple, self.upper_exact_mpf_tuple))


def _mpf_tuple(value: Any, label: str) -> tuple[int, int, int, int]:
    if not isinstance(value, (list, tuple)) or len(value) != 4:
        raise MatrixEnclosureInputError(label + " must be a four-integer mpf tuple")
    if any(type(item) is not int for item in value):
        raise MatrixEnclosureInputError(label + " may not contain decimal strings or floats")
    result = tuple(int(item) for item in value)
    if result in NONFINITE_MPF_TUPLES:
        raise MatrixEnclosureInputError(label + " must be finite")
    return result


def _finite_tuple(value: Any, label: str) -> tuple[int, int, int, int]:
    result = tuple(value)
    if len(result) != 4 or any(type(item) is not int for item in result):
        raise MatrixEnclosureInputError(label + " is not an exact four-integer mpf tuple")
    if result in NONFINITE_MPF_TUPLES:
        raise MatrixEnclosureInputError(label + " is non-finite")
    return result


def _mapping(value: Any, label: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise MatrixEnclosureInputError(label + " must be a JSON object")
    return value


def _path(value: Any, label: str) -> str:
    if not isinstance(value, str) or not value:
        raise MatrixEnclosureInputError(label + " must be a nonempty source path")
    return value.replace("\\", "/")


def _sha(value: Any, label: str) -> str:
    if not isinstance(value, str) or len(value) != 64:
        raise MatrixEnclosureInputError(label + " must be a 64-character SHA-256 string")
    try:
        int(value, 16)
    except ValueError as exc:
        raise MatrixEnclosureInputError(label + " is not hexadecimal SHA-256") from exc
    return value.lower()


def _validate_source_identity(bundle: Mapping[str, Any]) -> dict[str, Any]:
    if bundle.get("entered") is not True or bundle.get("source_frozen") is not True:
        raise MatrixEnclosureInputError("root must freeze source and entry before consumer call")
    A = _mapping(bundle.get("A"), "bundle.A")
    B = _mapping(bundle.get("B"), "bundle.B")
    _mapping(A.get("upstream"), "A.upstream")
    _mapping(B.get("upstream"), "B.upstream")
    if A.get("component") != COMPONENTS["A"] or B.get("component") != COMPONENTS["B"]:
        raise MatrixEnclosureInputError("source components do not match A/B capsule paths")
    if _path(A.get("result_path"), "A.result_path") != SOURCE_PATHS["A_result"]:
        raise MatrixEnclosureInputError("A.result_path is not the frozen core-source result")
    if _path(B.get("result_path"), "B.result_path") != SOURCE_PATHS["B_result"]:
        raise MatrixEnclosureInputError("B.result_path is not the frozen Q-chart result")
    A_hash = _sha(A.get("result_sha256"), "A.result_sha256")
    B_hash = _sha(B.get("result_sha256"), "B.result_sha256")
    A_ids = _mapping(A.get("source_ids"), "A.source_ids")
    B_ids = _mapping(B.get("source_ids"), "B.source_ids")
    eta_hash = _sha(A_ids.get("eta_result_sha256"), "A.source_ids.eta_result_sha256")
    if eta_hash != B_hash:
        raise MatrixEnclosureInputError(
            "A.source_ids.eta_result_sha256 must equal B.result_sha256"
        )
    if eta_hash != EXPECTED_ETA_RESULT_SHA256:
        raise MatrixEnclosureInputError(
            "A/B eta source identity does not match the frozen expected SHA-256"
        )
    common_A = _sha(A_ids.get("common_capsule_sha256"), "A common capsule id")
    common_B = _sha(B_ids.get("common_capsule_sha256"), "B common capsule id")
    if common_A != common_B:
        raise MatrixEnclosureInputError("A/B common capsule identity mismatch")
    axis_A = _sha(A_ids.get("direct_axis_result_sha256"), "A direct-axis id")
    axis_B = _sha(B_ids.get("anchor_axis_result_sha256"), "B anchor-axis id")
    if axis_A != axis_B:
        raise MatrixEnclosureInputError("A/B direct-axis identity mismatch")
    exact_view = _mapping(
        A.get("exact_view_registration"), "A.exact_view_registration"
    )
    if exact_view.get("conditional_leading_identity_registered") is not True:
        raise MatrixEnclosureInputError(
            "A exact_view_registration lacks conditional leading identity"
        )
    direct_eta_proof = B.get("direct_eta_chart_proof")
    if direct_eta_proof is None:
        # The source leaf records this proof as ``proof_packet``.  A root
        # freeze may expose it under the more explicit consumer key above.
        direct_eta_proof = B.get("proof_packet")
    direct_eta_proof = _mapping(
        direct_eta_proof, "B.direct_eta_chart_proof"
    )
    if B.get("conditional") is not True:
        raise MatrixEnclosureInputError("B direct eta-chart proof is not conditional")
    if not (
        direct_eta_proof.get("direct_eta_chart")
        or direct_eta_proof.get("eta_chart")
        or direct_eta_proof.get("kind") == "native_O3_second_Z_eta_chart_proof_v1"
    ):
        raise MatrixEnclosureInputError(
            "B direct eta-chart proof marker is missing"
        )
    for label, item in (("A", A), ("B", B)):
        if item.get("conditional") is not True:
            raise MatrixEnclosureInputError(label + " source must remain conditional")
        if item.get("absolute_F0_exp_materialized") is True:
            raise MatrixEnclosureInputError(label + " materialized absolute F0")
    return {
        "A_result_sha256": A_hash,
        "B_result_sha256": B_hash,
        "eta_result_sha256": eta_hash,
        "common_capsule_sha256": common_A,
        "direct_axis_result_sha256": axis_A,
        "A_source_ids": deepcopy(dict(A_ids)),
        "B_source_ids": deepcopy(dict(B_ids)),
        "conditional_premises": {
            "A_exact_view_registration": {
                "schema": exact_view.get("schema"),
                "conditional_leading_identity_registered": True,
                "regular_n1_view": exact_view.get("regular_n1_view"),
            },
            "B_direct_eta_chart_proof": {
                "kind": direct_eta_proof.get("kind"),
                "eta_chart": direct_eta_proof.get("eta_chart"),
                "conditional": B.get("conditional"),
            },
        },
        "hash_is_proof": HASH_IS_PROOF,
        "conditional": CONDITIONAL,
    }


def _selected_record(container: Mapping[str, Any], key: str, label: str) -> FrozenInterval:
    return FrozenInterval.from_json(container.get(key), label + "." + key)


def _contains_zero(value: Any, scalar: Any = None) -> bool:
    del scalar
    lower, upper = value._mpi_
    lower = _finite_tuple(lower, "interval.lower")
    upper = _finite_tuple(upper, "interval.upper")
    return mpf_cmp(lower, ZERO_MPF) <= 0 and mpf_cmp(upper, ZERO_MPF) >= 0


def _endpoint_record(value: Any, label: str = "interval") -> dict[str, Any]:
    lower, upper = value._mpi_
    lower = _finite_tuple(lower, label + ".lower")
    upper = _finite_tuple(upper, label + ".upper")
    if mpf_cmp(lower, upper) > 0:
        raise MatrixEnclosureInputError(label + " has reversed interval endpoints")
    return {
        "lower_exact_mpf_tuple": [int(item) for item in lower],
        "upper_exact_mpf_tuple": [int(item) for item in upper],
        "point_value": POINT_VALUE,
    }


def _entry(value: Any, factor: Any) -> dict[str, Any]:
    return {
        "factor": factor,
        "C0": _endpoint_record(value, "matrix.C0"),
        "point_value": POINT_VALUE,
    }


def _checked(value: Any, label: str) -> Any:
    """Reject a non-finite interval before it participates in later algebra."""

    _endpoint_record(value, label)
    return value


def _validate_factor_layout(
    B0: Mapping[str, Mapping[str, Any]], B1: Mapping[str, Mapping[str, Any]]
) -> dict[str, Any]:
    audits: dict[str, Any] = {}
    for name, table, dependent in (
        ("B0", B0, F_DEPENDENT_B0),
        ("B1", B1, F_DEPENDENT_B1),
    ):
        seen_dependent: list[str] = []
        seen_unit: list[str] = []
        for key, entry in table.items():
            expected = F0_ANCHOR if key in dependent else 1
            if entry.get("factor") != expected:
                raise MatrixEnclosureInputError(
                    f"{name}[{key}] has an unexpected formal F factor"
                )
            (seen_dependent if key in dependent else seen_unit).append(key)
        if set(seen_dependent) != set(dependent):
            raise MatrixEnclosureInputError(
                f"{name} formal F-dependent factor list is incomplete"
            )
        audits[name] = {
            "F_dependent": sorted(seen_dependent),
            "unit_factor": sorted(seen_unit),
        }
    return audits


def _compute_interval_matrix(bundle: Mapping[str, Any], dps: int) -> dict[str, Any]:
    A = _mapping(bundle["A"], "bundle.A")
    B = _mapping(bundle["B"], "bundle.B")
    parameters = _mapping(A.get("parameters"), "A.parameters")
    DU = _mapping(A.get("DU"), "A.DU")
    F_jets = _mapping(A.get("F_jets"), "A.F_jets")

    scalar = mp.mp.clone()
    scalar.dps = dps
    iv = MPIntervalContext()
    iv.dps = dps

    def load(record: Any, label: str) -> Any:
        frozen = FrozenInterval.from_json(record, label)
        return _checked(frozen.interval(iv, scalar), label)

    delta = load(parameters.get("delta"), "A.parameters.delta")
    epsilon = load(parameters.get("epsilon_core"), "A.parameters.epsilon_core")
    j = load(parameters.get("j"), "A.parameters.j")
    z = load(A.get("z"), "A.gauge_table.uniform_eta.uniform_rho.z")
    DU_eta = _mapping(DU.get("uniform_eta"), "A.DU.uniform_eta")
    DU_rho = _mapping(
        DU_eta.get("uniform_rho"), "A.mixed_jet_table.DU.uniform_eta.uniform_rho"
    )
    du00 = _checked(
        _selected_record(
            DU_rho, "m0_k0", "A.mixed_jet_table.DU.uniform_eta.uniform_rho"
        ).interval(iv, scalar),
        "A.mixed_jet_table.DU.uniform_eta.uniform_rho.m0_k0",
    )
    du10 = _checked(
        _selected_record(
            DU_rho, "m0_k1", "A.mixed_jet_table.DU.uniform_eta.uniform_rho"
        ).interval(iv, scalar),
        "A.mixed_jet_table.DU.uniform_eta.uniform_rho.m0_k1",
    )
    du01 = _checked(
        _selected_record(
            DU_rho, "m1_k0", "A.mixed_jet_table.DU.uniform_eta.uniform_rho"
        ).interval(iv, scalar),
        "A.mixed_jet_table.DU.uniform_eta.uniform_rho.m1_k0",
    )
    f = load(
        F_jets.get("m0_k0"),
        "A.gauge_table.uniform_eta.uniform_rho.F_jets.m0_k0",
    )
    fr = load(
        F_jets.get("m0_k1"),
        "A.gauge_table.uniform_eta.uniform_rho.F_jets.m0_k1",
    )
    fz = load(
        F_jets.get("m1_k0"),
        "A.gauge_table.uniform_eta.uniform_rho.F_jets.m1_k0",
    )
    Q = load(
        _mapping(B.get("Q_chart"), "B.compatible_pressure_Q.uniform_eta.rows.uniform_rho.Q_chart"),
        "B.compatible_pressure_Q.uniform_eta.rows.uniform_rho.Q_chart",
    )

    # These are the only newly introduced coordinate intervals.  Construct
    # the rho endpoints directly in the interval context; no hand-written
    # endpoint tuple is allowed here.
    rho = _checked(iv.mpf([0, 4]), "rho")
    x = _checked(iv.sqrt(rho), "x")
    U = _checked(4 * z + j + du00, "U")
    U_r = _checked(du10, "U_r")
    U_Z = _checked(4 + du01, "U_Z")
    d = _checked(1 - z * z, "d")
    L = _checked(1 - delta * z * z, "L")
    if _contains_zero(L, scalar):
        raise MatrixEnclosureInputError("L interval contains zero")
    H = _checked(z * ((1 - delta) / 2) + d * U, "H")
    Bcoef = _checked(Q + (1 - 2 * z * U) / L, "B")
    Ecoef = _checked(z * U - iv.mpf(1) / 2, "E")
    Cf = _checked(f + fr * rho, "Cf")
    Zf = _checked(
        (f * (z * (-2 - delta)) + fz * d - fr * (z * (2 * rho))) / L,
        "Zf",
    )
    Zw = _checked(
        (U * (z * (-1 - delta)) + U_Z * d - U_r * (z * (2 * rho))) / L,
        "Zw",
    )
    eps = lambda value: epsilon * value

    B0: dict[str, dict[str, Any]] = {
        key: _entry(iv.mpf(value), 1) for key, value in B0_KINEMATIC.items()
    }
    B0.update(
        {
            "4,1": _entry(eps(f * (4 * x)), F0_ANCHOR),
            "5,1": _entry(eps((Q + Ecoef * (-2 + delta) / L) * 2), 1),
            "5,2": _entry(eps(Zf * 2 - Cf * (z * (2 * (-1 + delta))) / L), F0_ANCHOR),
            "5,3": _entry(eps(-Cf * (z * (2 * (1 + delta))) / L), F0_ANCHOR),
            "5,5": _entry(eps(Bcoef * x), 1),
            "6,1": _entry(eps(eps(-f * (z * (8 * rho)) / L)), F0_ANCHOR),
            "6,2": _entry(eps((Ecoef * (-1 + delta) / L + Zw - U_r * (z * ((-1 + delta) * rho)) / L) * 2), 1),
            "6,3": _entry(eps(-U_r * (z * (2 * (1 + delta) * rho)) / L), 1),
            "6,4": _entry(eps(z * (-4) / L), 1),
            "6,6": _entry(eps(Bcoef * x), 1),
        }
    )
    B1: dict[str, dict[str, Any]] = {
        "5,1": _entry(eps(H * 2 / L), 1),
        "5,2": _entry(eps(-Cf * d * 2 / L), F0_ANCHOR),
        "5,3": _entry(eps(-Cf * d * 2 / L), F0_ANCHOR),
        "6,2": _entry(eps((H - U_r * d * rho) * 2 / L), 1),
        "6,3": _entry(eps(-U_r * d * (2 * rho) / L), 1),
        "6,4": _entry(eps(d * 2 / L), 1),
    }
    factor_audit = _validate_factor_layout(B0, B1)
    return {
        "rho": _endpoint_record(rho, "rho"),
        "x": _endpoint_record(x, "x"),
        "z": _endpoint_record(z, "z"),
        "delta": _endpoint_record(delta, "delta"),
        "epsilon_core": _endpoint_record(epsilon, "epsilon_core"),
        "j": _endpoint_record(j, "j"),
        "U": _endpoint_record(U, "U"),
        "U_r": _endpoint_record(U_r, "U_r"),
        "U_Z": _endpoint_record(U_Z, "U_Z"),
        "Q": _endpoint_record(Q, "Q"),
        "d": _endpoint_record(d, "d"),
        "L": _endpoint_record(L, "L"),
        "B0": B0,
        "B1": B1,
        "factor_audit": factor_audit,
        "F_inputs": {
            "factor": F0_ANCHOR,
            "F_relative": _endpoint_record(f, "F_relative"),
            "F_r": _endpoint_record(fr, "F_r"),
            "F_Z": _endpoint_record(fz, "F_Z"),
            "F0_materialized": False,
        },
        "point_value": POINT_VALUE,
        "all_computed_outputs_finite": True,
        "interval_dps": dps,
    }


def consume_frozen_source_bundle(
    bundle: Mapping[str, Any], *, interval_dps: int = 280
) -> dict[str, Any]:
    """Consume one root-frozen A/B bundle and return new B0/B1 C0 enclosures.

    ``bundle.A`` and ``bundle.B`` contain already-extracted selected records,
    not paths to be opened.  The caller must retain full source result records
    (including old ``error`` objects) under ``A['upstream']`` and
    ``B['upstream']``; this function copies them unchanged into the result.
    """

    if type(interval_dps) is not int or interval_dps < 32:
        raise MatrixEnclosureInputError("interval_dps must be an integer >= 32")
    bundle = _mapping(bundle, "bundle")
    identity = _validate_source_identity(bundle)
    A = _mapping(bundle.get("A"), "bundle.A")
    B = _mapping(bundle.get("B"), "bundle.B")
    matrix = _compute_interval_matrix(bundle, interval_dps)
    return {
        "schema": SCHEMA,
        "runtime_model_metadata": RUNTIME_MODEL_METADATA,
        "conditional": CONDITIONAL,
        "hash_is_proof": HASH_IS_PROOF,
        "consumer_call_count": 1,
        "root_freeze_required": True,
        "source_identity": identity,
        "source_paths": dict(SOURCE_PATHS),
        "source_components": dict(COMPONENTS),
        "matrix": matrix,
        "raw_formulas": dict(RAW_FORMULAS),
        "interval_engine": INTERVAL_ENGINE,
        "fields": list(FIELDS),
        "diagonal": list(DIAGONAL),
        "formal_F0_anchor": F0_ANCHOR,
        "F0_materialized": False,
        "point_value": POINT_VALUE,
        "all_computed_outputs_finite": matrix["all_computed_outputs_finite"],
        "matrix_C1": False,
        "source_function_value_claim": False,
        "W_value_claim": False,
        "old212_error_reduction_claim": False,
        "global_closure_claim": False,
        "selected_inputs": {
            "A": {
                "parameters": deepcopy(A.get("parameters")),
                "DU": deepcopy(A.get("DU")),
                "F_jets": deepcopy(A.get("F_jets")),
                "z": deepcopy(A.get("z")),
            },
            "B": {"Q_chart": deepcopy(B.get("Q_chart"))},
        },
        "upstream": {
            "A": deepcopy(A.get("upstream")),
            "B": deepcopy(B.get("upstream")),
        },
    }


__all__ = [
    "SCHEMA",
    "RUNTIME_MODEL_METADATA",
    "SOURCE_PATHS",
    "COMPONENTS",
    "EXPECTED_ETA_RESULT_SHA256",
    "INTERVAL_ENGINE",
    "RAW_FORMULAS",
    "FrozenInterval",
    "MatrixEnclosureInputError",
    "consume_frozen_source_bundle",
]
