"""Fail-closed physical-``Rp`` support witness for reserved m15.

This module is a source-port consumer, not a producer.  It does not import a
field, graph, archive, interval context, or numerical helper.  The caller
must supply two live source-port records (reserved bump and quiet ``Rp``),
the same *object identities* for the live ``Rc`` owner and parameter owner,
and one shared opaque proof-token object.  Identity checks catch accidental
mixing of source ports; they cannot authenticate an arbitrary object supplied
by a caller.  Consequently this v1 bridge always reports
``source_owner_unverified=True`` and the resulting full-support statement is
conditional until an independently trusted owner layer authenticates those
ports.

The geometry part is exact at the symbolic/log-coordinate level supplied by
the handoff.  Geometry alone does **not** imply the inverse response row.
The caller must also provide an explicit, separately unverified premise that
the complete three-bump m15 moment row equals the third row of the same
invertible angular matrix, and that this equality is a full ordinary-Z C1
function identity.  The bridge records those premises and keeps them
conditional; it does not inspect row entries, derivatives, or a determinant.

The geometric facts are:

* reserved angular channels are ``(1, 2, 3)``;
* the five source bump center labels are ``(12, 13, 14, 15, 16)``;
* the largest bump support is ``(151/400) * log(2)`` above ``Rc``;
* ``2 Rc`` is at ``1 * log(2)`` above ``Rc``, leaving the exact margin
  ``(249/400) * log(2)``;
* ``Rp = Rc * exp(Tw - 2)`` and the supplied quiet witness is
  ``Tw - 2 - log(2) > 0``.

Thus the bridge conditionally records that the source-template support is
before ``2 Rc`` and, conditionally on the supplied physical mapping witness,
before the actual ``Rp``.  It never evaluates ``log`` or ``exp`` and never
treats a name, boolean, or equal value as an owner identity.

Minimal use (illustrative only; the result remains conditional)::

    token = object()
    rc_owner = object()
    parameter_owner = object()
    ports = (
        LiveSourceOwnerPort("reserved_bump", rc_owner, parameter_owner, token),
        LiveSourceOwnerPort("quiet_rp", rc_owner, parameter_owner, token),
    )
    mapping = PhysicalRpMapping(
        quiet_length_positive=PositiveSourceWitness("Tw-2-log2>0", token),
        rc_positive=PositiveSourceWitness("Rc>0", token),
        exponential_order=PositiveSourceWitness(
            "exp_strictly_increasing_on_real", token
        ),
    )
    matrix = FullSupportMatrixPremise(
        complete_m15_moment_row=object(),
        angular_matrix_third_row=object(),
        row_identity_witness=SourcePremiseWitness(
            "same_three_bump_complete_m15_moment_row_equals_angular_matrix_third_row",
            token,
        ),
        row_C1_identity_witness=SourcePremiseWitness(
            "same_three_bump_complete_m15_moment_row_equals_angular_matrix_third_row_as_full_Z_C1_function",
            token,
        ),
        invertibility_witness=SourcePremiseWitness(
            "same_three_bump_angular_matrix_invertible", token
        ),
    )
    witness = prove_reserved_rp_full_support(
        owner_ports=ports, mapping=mapping, matrix_premise=matrix
    )
    # ``caller_iv`` is the caller's mpmath interval context.
    row_inputs = witness.conditional_variable_z_response_inputs(
        caller_iv, acknowledge_unverified_conditions=True
    )
    # row_inputs has only ``r_i`` and ``r_i_Z``; J/P remain caller inputs.
    # row == (0, 0, 1), but witness.exact_matrix_row_certified is False.

The example deliberately uses opaque objects rather than source names.  It is
not a claim that those objects are the repository's live owners.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from fractions import Fraction
from typing import Any, Final


SCHEMA: Final[str] = "reserved_rp_support_bridge_v1"

RESERVED_ANGULAR: Final[tuple[int, int, int]] = (1, 2, 3)
RESERVED_BUMP_CENTERS: Final[tuple[int, int, int, int, int]] = (12, 13, 14, 15, 16)
BUMP_SUPPORT_UPPER_LOG2: Final[Fraction] = Fraction(151, 400)
TWO_RC_LOG2: Final[Fraction] = Fraction(1, 1)
FULL_SUPPORT_MARGIN_LOG2: Final[Fraction] = TWO_RC_LOG2 - BUMP_SUPPORT_UPPER_LOG2
EXACT_RESPONSE_ROW: Final[tuple[int, int, int]] = (0, 0, 1)
EXACT_RESPONSE_ROW_Z1: Final[tuple[int, int, int]] = (0, 0, 0)

EXPECTED_RESERVED_ROLE: Final[str] = "reserved_bump"
EXPECTED_QUIET_ROLE: Final[str] = "quiet_rp"
EXPECTED_RP_FORMULA: Final[str] = "Rp=Rc*exp(Tw-2)"
EXPECTED_TWO_RC_FORMULA: Final[str] = "2Rc=Rc*exp(log2)"
EXPECTED_QUIET_LENGTH_FORMULA: Final[str] = "quiet_length=Tw-2-log2"
EXPECTED_QUIET_LENGTH_PROOF: Final[str] = "Tw-2-log2>0"
EXPECTED_RC_POSITIVE_PROOF: Final[str] = "Rc>0"
EXPECTED_EXP_ORDER_PROOF: Final[str] = "exp_strictly_increasing_on_real"
EXPECTED_ROW_IDENTITY_PROOF: Final[str] = (
    "same_three_bump_complete_m15_moment_row_equals_angular_matrix_third_row"
)
EXPECTED_ROW_C1_IDENTITY_PROOF: Final[str] = (
    "same_three_bump_complete_m15_moment_row_equals_angular_matrix_third_row_as_full_Z_C1_function"
)
EXPECTED_MATRIX_INVERTIBLE_PROOF: Final[str] = "same_three_bump_angular_matrix_invertible"


class SupportBridgeError(ValueError):
    """Raised when a required source or mapping witness is not fail-closed."""


@dataclass(frozen=True, slots=True)
class LiveSourceOwnerPort:
    """A caller-owned live source port, never a name or graph-node number.

    ``proof_token`` is checked by object identity with the other ports and
    mapping witnesses.  This validates the submitted wiring only; v1 has no
    authority to certify that an arbitrary token came from the real producer.
    """

    role: str
    rc_owner: object
    parameter_owner: object
    proof_token: object


@dataclass(frozen=True, slots=True)
class PositiveSourceWitness:
    """Opaque caller evidence for one exact positive/mapping premise.

    The expression is restricted to the literals used by this bridge.  The
    token is linked by identity to the source-port token, but its provenance
    remains unverified by this module.
    """

    expression: str
    proof_token: object


@dataclass(frozen=True, slots=True)
class SourcePremiseWitness:
    """Opaque caller evidence for a non-numeric source identity premise."""

    statement: str
    proof_token: object


@dataclass(frozen=True, slots=True)
class PhysicalRpMapping:
    """Explicit assumptions needed to transport compact support to physical Rp."""

    quiet_length_positive: PositiveSourceWitness
    rc_positive: PositiveSourceWitness
    exponential_order: PositiveSourceWitness
    rp_formula: str = EXPECTED_RP_FORMULA
    two_rc_formula: str = EXPECTED_TWO_RC_FORMULA
    quiet_length_formula: str = EXPECTED_QUIET_LENGTH_FORMULA


@dataclass(frozen=True, slots=True)
class FullSupportMatrixPremise:
    """Unverified source premise needed in addition to geometric support.

    The two row objects are intentionally opaque.  This module does not
    inspect or numerically compare them.  Their pointwise relation, their
    full ordinary-Z C1 function relation, and invertibility arrive as
    caller-owned witnesses tied by object identity to the source proof token;
    v1 records the premises as unverified rather than promoting them to a
    matrix theorem.
    """

    complete_m15_moment_row: object
    angular_matrix_third_row: object
    row_identity_witness: SourcePremiseWitness
    row_C1_identity_witness: SourcePremiseWitness
    invertibility_witness: SourcePremiseWitness


@dataclass(frozen=True, slots=True)
class OwnerBindingReport:
    """What was checked about caller-supplied source-port wiring."""

    roles: tuple[str, str]
    live_rc_identity_checked: bool
    parameter_identity_checked: bool
    proof_token_identity_checked: bool
    source_owner_authenticated: bool
    source_owner_unverified: bool
    status: str
    limitation: str


@dataclass(frozen=True, slots=True)
class MatrixPremiseReport:
    """Validation status for the row identity and angular invertibility premise."""

    complete_m15_row_supplied: bool
    angular_third_row_supplied: bool
    row_identity_witness_checked: bool
    row_C1_identity_witness_checked: bool
    invertibility_witness_checked: bool
    same_source_token_checked: bool
    matrix_premise_authenticated: bool
    matrix_premise_unverified: bool
    row_C1_identity_premise_unverified: bool
    status: str
    limitation: str


@dataclass(frozen=True, slots=True)
class OrdinaryZ1Contract:
    """Derivative bookkeeping retained after selecting the exact response row."""

    exact_response_row: tuple[int, int, int]
    exact_response_row_Z1: tuple[int, int, int]
    response_row_constant_in_Z_conditionally: bool
    response_row_constant_in_Z_certified: bool
    row_C1_identity_premise_unverified: bool
    source_function_derivatives_retained: bool
    retained_source_Z1_fields: tuple[str, ...]
    statement: str


@dataclass(frozen=True, slots=True)
class ReservedRpSupportWitness:
    """Conditional full-support result returned by the pure bridge."""

    schema: str
    angular_channels: tuple[int, int, int]
    bump_center_labels: tuple[int, int, int, int, int]
    bump_support_upper_log2_coefficient: Fraction
    two_rc_log2_coefficient: Fraction
    exact_margin_log2_coefficient: Fraction
    exact_margin_expression: str
    support_before_two_rc_conditionally: bool
    support_before_two_rc_certified: bool
    rp_beyond_two_rc_conditionally: bool
    support_before_rp_conditionally: bool
    full_support_status: str
    owner_binding: OwnerBindingReport
    matrix_premise: MatrixPremiseReport
    physical_mapping_assumptions: tuple[str, ...]
    exact_matrix_row: tuple[int, int, int]
    exact_matrix_row_certified: bool
    exact_matrix_row_conditionally_available: bool
    ordinary_z1: OrdinaryZ1Contract
    limitations: tuple[str, ...]

    def __post_init__(self) -> None:
        _validate_witness_constants(self)

    def conditional_exact_response_row(
        self,
        *,
        acknowledge_unverified_conditions: bool = False,
    ) -> tuple[int, int, int]:
        """Return the synthetic/conditional row only after explicit acknowledgment.

        This is not a production authentication path.  A real consumer must
        use a separately implemented certified interface; v1's certified
        method always fails closed.
        """

        _validate_witness_constants(self)
        if acknowledge_unverified_conditions is not True:
            raise SupportBridgeError(
                "conditional row requires acknowledge_unverified_conditions=True"
            )
        if not self.exact_matrix_row_conditionally_available:
            raise SupportBridgeError("physical full-support premise is unresolved")
        return self.exact_matrix_row

    def certified_exact_response_row(self) -> tuple[int, int, int]:
        """Always reject: v1 has no mechanism for authenticating physical owners."""

        raise SupportBridgeError(
            "v1 has no certified owner/matrix interface; certified physical row unavailable"
        )

    def certified_variable_z_response_inputs(self, iv: Any) -> dict[str, tuple[Any, Any, Any]]:
        """Always reject: v1 cannot provide certified variable-Z inputs."""

        del iv
        raise SupportBridgeError(
            "v1 has no certified owner/matrix interface; certified variable-Z inputs unavailable"
        )

    def conditional_variable_z_response_inputs(
        self,
        iv: Any,
        *,
        acknowledge_unverified_conditions: bool = False,
    ) -> dict[str, tuple[Any, Any, Any]]:
        """Return synthetic/conditional interval rows after explicit acknowledgment.

        ``iv`` must be the caller's mpmath interval context.  The returned
        ``r_i``/``r_i_Z`` contain interval values produced by ``iv.mpf`` and
        can be passed as the only two bridge-owned keyword arguments to a
        caller-owned ``evaluate_partial_reserved_minus_pulse`` invocation.
        ``J``, ``P`` and their ordinary-Z1 rows remain live inputs; this method
        deliberately supplies no zero for any source moment or pulse term.
        A real authenticated consumer must use a separate certified interface;
        v1's certified method always fails closed.
        """

        _validate_witness_constants(self)
        if acknowledge_unverified_conditions is not True:
            raise SupportBridgeError(
                "conditional variable-Z inputs require acknowledge_unverified_conditions=True"
            )
        if not self.exact_matrix_row_conditionally_available:
            raise SupportBridgeError("physical full-support premise is unresolved")
        if iv is None or not callable(getattr(iv, "mpf", None)):
            raise SupportBridgeError("iv must be the caller-provided interval context")
        try:
            zero = iv.mpf("0")
            one = iv.mpf("1")
        except (AttributeError, TypeError, ValueError) as exc:
            raise SupportBridgeError("iv.mpf must construct interval constants") from exc
        if not _exact_interval_constant(zero, 0) or not _exact_interval_constant(one, 1):
            raise SupportBridgeError(
                "iv.mpf('0') and iv.mpf('1') must be trusted exact intervals with matching endpoints"
            )
        return {
            "r_i": (zero, zero, one),
            "r_i_Z": (zero, zero, zero),
        }

    def report(self) -> dict[str, Any]:
        """Return a serialization-safe report without leaking owner objects/tokens."""

        return {
            "schema": self.schema,
            "angular_channels": self.angular_channels,
            "bump_center_labels": self.bump_center_labels,
            "support": {
                "upper_log2_coefficient": f"{self.bump_support_upper_log2_coefficient.numerator}/{self.bump_support_upper_log2_coefficient.denominator}",
                "two_rc_log2_coefficient": f"{self.two_rc_log2_coefficient.numerator}/{self.two_rc_log2_coefficient.denominator}",
                "exact_margin": self.exact_margin_expression,
                "support_before_two_rc_conditionally": self.support_before_two_rc_conditionally,
                "support_before_two_rc_certified": self.support_before_two_rc_certified,
                "geometry_source_template_only": True,
                "support_before_rp_conditionally": self.support_before_rp_conditionally,
            },
            "full_support_status": self.full_support_status,
            "owner_binding": {
                "roles": self.owner_binding.roles,
                "live_rc_identity_checked": self.owner_binding.live_rc_identity_checked,
                "parameter_identity_checked": self.owner_binding.parameter_identity_checked,
                "proof_token_identity_checked": self.owner_binding.proof_token_identity_checked,
                "source_owner_authenticated": self.owner_binding.source_owner_authenticated,
                "source_owner_unverified": self.owner_binding.source_owner_unverified,
                "status": self.owner_binding.status,
                "limitation": self.owner_binding.limitation,
            },
            "matrix_row_premise": {
                "complete_m15_row_supplied": self.matrix_premise.complete_m15_row_supplied,
                "angular_third_row_supplied": self.matrix_premise.angular_third_row_supplied,
                "row_identity_witness_checked": self.matrix_premise.row_identity_witness_checked,
                "row_C1_identity_witness_checked": self.matrix_premise.row_C1_identity_witness_checked,
                "invertibility_witness_checked": self.matrix_premise.invertibility_witness_checked,
                "same_source_token_checked": self.matrix_premise.same_source_token_checked,
                "matrix_premise_authenticated": self.matrix_premise.matrix_premise_authenticated,
                "matrix_premise_unverified": self.matrix_premise.matrix_premise_unverified,
                "row_C1_identity_premise_unverified": self.matrix_premise.row_C1_identity_premise_unverified,
                "status": self.matrix_premise.status,
                "limitation": self.matrix_premise.limitation,
            },
            "physical_mapping_assumptions": self.physical_mapping_assumptions,
            "exact_matrix_row": self.exact_matrix_row,
            "exact_matrix_row_certified": self.exact_matrix_row_certified,
            "exact_matrix_row_conditionally_available": self.exact_matrix_row_conditionally_available,
            "variable_z_response_metadata": {
                "conditional_method_name": "conditional_variable_z_response_inputs",
                "certified_method_name": "certified_variable_z_response_inputs",
                "consumer_keyword_names": ("r_i", "r_i_Z"),
                "interval_context_required": True,
                "exact_interval_endpoints_required": {"zero": 0, "one": 1},
                "acknowledge_unverified_conditions_required": True,
                "response_row_is_conditional": not self.exact_matrix_row_certified,
                "response_row_Z1_is_conditional": not self.ordinary_z1.response_row_constant_in_Z_certified,
                "full_Z_C1_row_identity_premise_unverified": self.ordinary_z1.row_C1_identity_premise_unverified,
                "J_must_remain_live": True,
                "P_must_remain_live": True,
                "J_Z_must_remain_live": True,
                "P_Z_must_remain_live": True,
                "zero_substitution_forbidden": ("J", "P", "J_Z", "P_Z"),
            },
            "ordinary_Z1": {
                "exact_response_row": self.ordinary_z1.exact_response_row,
                "exact_response_row_Z1": self.ordinary_z1.exact_response_row_Z1,
                "response_row_constant_in_Z_conditionally": self.ordinary_z1.response_row_constant_in_Z_conditionally,
                "response_row_constant_in_Z_certified": self.ordinary_z1.response_row_constant_in_Z_certified,
                "row_C1_identity_premise_unverified": self.ordinary_z1.row_C1_identity_premise_unverified,
                "source_function_derivatives_retained": self.ordinary_z1.source_function_derivatives_retained,
                "retained_source_Z1_fields": self.ordinary_z1.retained_source_Z1_fields,
                "statement": self.ordinary_z1.statement,
            },
            "limitations": self.limitations,
        }


def _validate_witness_constants(witness: ReservedRpSupportWitness) -> None:
    """Reject forged fixed rows/geometry even for hand-built dataclasses."""

    if not isinstance(witness, ReservedRpSupportWitness):
        raise SupportBridgeError("reserved Rp witness type required")
    if witness.schema != SCHEMA:
        raise SupportBridgeError("reserved Rp witness schema is not v1")
    if witness.angular_channels != RESERVED_ANGULAR:
        raise SupportBridgeError("reserved angular channels are not the fixed (1,2,3) template")
    if witness.bump_center_labels != RESERVED_BUMP_CENTERS:
        raise SupportBridgeError("reserved bump center labels are not the fixed source template")
    if witness.bump_support_upper_log2_coefficient != BUMP_SUPPORT_UPPER_LOG2:
        raise SupportBridgeError("reserved bump support upper coefficient was changed")
    if witness.two_rc_log2_coefficient != TWO_RC_LOG2:
        raise SupportBridgeError("2Rc log2 coefficient was changed")
    if witness.exact_margin_log2_coefficient != Fraction(249, 400):
        raise SupportBridgeError("reserved support margin is not 249/400*log2")
    if witness.exact_margin_expression != "249/400*log2":
        raise SupportBridgeError("reserved support margin expression was changed")
    if witness.support_before_two_rc_conditionally is not True:
        raise SupportBridgeError("conditional source-template support premise is missing")
    if witness.support_before_two_rc_certified is not False:
        raise SupportBridgeError("source-template support cannot be certified by v1")
    if witness.rp_beyond_two_rc_conditionally is not True:
        raise SupportBridgeError("conditional Rp mapping premise is missing")
    if witness.support_before_rp_conditionally is not True:
        raise SupportBridgeError("conditional physical Rp support premise is missing")
    if witness.full_support_status != "conditional_support_and_matrix_source_unverified":
        raise SupportBridgeError("full-support status must remain conditional in v1")
    if witness.exact_matrix_row != EXACT_RESPONSE_ROW:
        raise SupportBridgeError("conditional response row is not the fixed (0,0,1) row")
    if witness.exact_matrix_row_certified is not False:
        raise SupportBridgeError("v1 cannot certify the exact response row")
    if witness.exact_matrix_row_conditionally_available is not True:
        raise SupportBridgeError("conditional response row premise is missing")
    if not isinstance(witness.owner_binding, OwnerBindingReport):
        raise SupportBridgeError("owner binding report is missing")
    if witness.owner_binding.roles != (EXPECTED_RESERVED_ROLE, EXPECTED_QUIET_ROLE):
        raise SupportBridgeError("owner binding roles were changed")
    if witness.owner_binding.live_rc_identity_checked is not True:
        raise SupportBridgeError("live Rc identity check is missing")
    if witness.owner_binding.parameter_identity_checked is not True:
        raise SupportBridgeError("parameter identity check is missing")
    if witness.owner_binding.proof_token_identity_checked is not True:
        raise SupportBridgeError("proof-token identity check is missing")
    if witness.owner_binding.source_owner_authenticated is not False:
        raise SupportBridgeError("v1 cannot authenticate source owners")
    if witness.owner_binding.source_owner_unverified is not True:
        raise SupportBridgeError("source owner must remain explicitly unverified")
    if witness.owner_binding.status != "conditional_source_owner_unverified":
        raise SupportBridgeError("owner binding status must remain conditional")
    if not isinstance(witness.matrix_premise, MatrixPremiseReport):
        raise SupportBridgeError("matrix premise report is missing")
    if witness.matrix_premise.complete_m15_row_supplied is not True:
        raise SupportBridgeError("complete m15 row premise is missing")
    if witness.matrix_premise.angular_third_row_supplied is not True:
        raise SupportBridgeError("angular third-row premise is missing")
    if witness.matrix_premise.row_identity_witness_checked is not True:
        raise SupportBridgeError("matrix row identity witness is missing")
    if witness.matrix_premise.row_C1_identity_witness_checked is not True:
        raise SupportBridgeError("full-Z C1 row identity witness is missing")
    if witness.matrix_premise.invertibility_witness_checked is not True:
        raise SupportBridgeError("matrix invertibility witness is missing")
    if witness.matrix_premise.same_source_token_checked is not True:
        raise SupportBridgeError("matrix/source token binding is missing")
    if witness.matrix_premise.matrix_premise_authenticated is not False:
        raise SupportBridgeError("v1 cannot authenticate the matrix premise")
    if witness.matrix_premise.matrix_premise_unverified is not True:
        raise SupportBridgeError("matrix premise must remain explicitly unverified")
    if witness.matrix_premise.row_C1_identity_premise_unverified is not True:
        raise SupportBridgeError("full-Z C1 row premise must remain unverified")
    if witness.matrix_premise.status != "conditional_matrix_identity_C1_and_invertibility_unverified":
        raise SupportBridgeError("matrix premise status must remain conditional")
    if not isinstance(witness.ordinary_z1, OrdinaryZ1Contract):
        raise SupportBridgeError("ordinary-Z1 contract is missing")
    if witness.ordinary_z1.exact_response_row != EXACT_RESPONSE_ROW:
        raise SupportBridgeError("ordinary-Z1 response row was changed")
    if witness.ordinary_z1.exact_response_row_Z1 != EXACT_RESPONSE_ROW_Z1:
        raise SupportBridgeError("ordinary-Z1 response derivative row was changed")
    if witness.ordinary_z1.response_row_constant_in_Z_conditionally is not True:
        raise SupportBridgeError("conditional constant-row premise is missing")
    if witness.ordinary_z1.response_row_constant_in_Z_certified is not False:
        raise SupportBridgeError("constant-row C1 identity cannot be certified by v1")
    if witness.ordinary_z1.row_C1_identity_premise_unverified is not True:
        raise SupportBridgeError("full-Z C1 row identity must remain unverified")
    if witness.ordinary_z1.source_function_derivatives_retained is not True:
        raise SupportBridgeError("source-function ordinary-Z1 rows must remain live")


_NON_OWNER_VALUES = (type(None), bool, int, float, complex, str, bytes, bytearray, list, tuple, dict, set, frozenset)


def _require_live_identity(value: object, label: str) -> None:
    """Reject scalar/name/container stand-ins before any identity comparison."""

    if isinstance(value, _NON_OWNER_VALUES):
        raise SupportBridgeError(f"{label} must be a live opaque object, not a name/value")


def _exact_interval_constant(value: object, target: int) -> bool:
    """Recognize a genuine exact interval constant without coercion."""

    raw = getattr(value, "_mpi_", None)
    if not isinstance(raw, tuple) or len(raw) != 2:
        return False
    if not hasattr(value, "a") or not hasattr(value, "b"):
        return False
    try:
        lower_equal = value.a == target
        upper_equal = value.b == target
    except (AttributeError, TypeError, ValueError):
        return False
    return lower_equal is True and upper_equal is True


def _require_mapping_witness(value: object, expected: str, label: str, token: object) -> None:
    if not isinstance(value, PositiveSourceWitness):
        raise SupportBridgeError(f"{label} must be a PositiveSourceWitness")
    if value.expression != expected:
        raise SupportBridgeError(
            f"{label} has unsupported expression {value.expression!r}; expected {expected!r}"
        )
    _require_live_identity(value.proof_token, f"{label}.proof_token")
    if value.proof_token is not token:
        raise SupportBridgeError(f"{label}.proof_token is not the live source proof token")


def _require_source_premise(value: object, expected: str, label: str, token: object) -> None:
    if not isinstance(value, SourcePremiseWitness):
        raise SupportBridgeError(f"{label} must be a SourcePremiseWitness")
    if value.statement != expected:
        raise SupportBridgeError(
            f"{label} has unsupported statement {value.statement!r}; expected {expected!r}"
        )
    _require_live_identity(value.proof_token, f"{label}.proof_token")
    if value.proof_token is not token:
        raise SupportBridgeError(f"{label}.proof_token is not the live source proof token")


def _validated_matrix_premise(
    premise: FullSupportMatrixPremise,
    token: object,
) -> MatrixPremiseReport:
    """Validate the submitted row/invertibility premise without evaluating it."""

    if not isinstance(premise, FullSupportMatrixPremise):
        raise SupportBridgeError("matrix_premise must be a FullSupportMatrixPremise")
    _require_live_identity(premise.complete_m15_moment_row, "matrix_premise.complete_m15_moment_row")
    _require_live_identity(premise.angular_matrix_third_row, "matrix_premise.angular_matrix_third_row")
    _require_source_premise(
        premise.row_identity_witness,
        EXPECTED_ROW_IDENTITY_PROOF,
        "matrix_premise.row_identity_witness",
        token,
    )
    _require_source_premise(
        premise.row_C1_identity_witness,
        EXPECTED_ROW_C1_IDENTITY_PROOF,
        "matrix_premise.row_C1_identity_witness",
        token,
    )
    _require_source_premise(
        premise.invertibility_witness,
        EXPECTED_MATRIX_INVERTIBLE_PROOF,
        "matrix_premise.invertibility_witness",
        token,
    )
    return MatrixPremiseReport(
        complete_m15_row_supplied=True,
        angular_third_row_supplied=True,
        row_identity_witness_checked=True,
        row_C1_identity_witness_checked=True,
        invertibility_witness_checked=True,
        same_source_token_checked=True,
        matrix_premise_authenticated=False,
        matrix_premise_unverified=True,
        row_C1_identity_premise_unverified=True,
        status="conditional_matrix_identity_C1_and_invertibility_unverified",
        limitation=(
            "v1 validates the explicit witness statements and shared token wiring, "
            "but does not inspect the complete three-bump m15 row, compare matrix "
            "entries, or certify nonzero determinant"
        ),
    )


def _validated_owner_ports(
    owner_ports: Sequence[LiveSourceOwnerPort],
) -> tuple[LiveSourceOwnerPort, LiveSourceOwnerPort, OwnerBindingReport, object]:
    if isinstance(owner_ports, (str, bytes, bytearray)) or not isinstance(owner_ports, Sequence):
        raise SupportBridgeError("owner_ports must be a sequence of exactly two source ports")
    ports = tuple(owner_ports)
    if len(ports) != 2:
        raise SupportBridgeError("owner_ports must contain reserved_bump and quiet_rp exactly once")
    if any(not isinstance(port, LiveSourceOwnerPort) for port in ports):
        raise SupportBridgeError("owner_ports entries must be LiveSourceOwnerPort values")
    by_role = {port.role: port for port in ports}
    if set(by_role) != {EXPECTED_RESERVED_ROLE, EXPECTED_QUIET_ROLE}:
        raise SupportBridgeError(
            "owner_ports roles must be exactly reserved_bump and quiet_rp"
        )
    reserved = by_role[EXPECTED_RESERVED_ROLE]
    quiet = by_role[EXPECTED_QUIET_ROLE]
    for port in (reserved, quiet):
        _require_live_identity(port.rc_owner, f"{port.role}.rc_owner")
        _require_live_identity(port.parameter_owner, f"{port.role}.parameter_owner")
        _require_live_identity(port.proof_token, f"{port.role}.proof_token")
    rc_same = reserved.rc_owner is quiet.rc_owner
    parameter_same = reserved.parameter_owner is quiet.parameter_owner
    token_same = reserved.proof_token is quiet.proof_token
    if not rc_same:
        raise SupportBridgeError("reserved and quiet ports do not share one Rc owner object")
    if not parameter_same:
        raise SupportBridgeError("reserved and quiet ports do not share one parameter owner object")
    if not token_same:
        raise SupportBridgeError("reserved and quiet ports do not share one proof-token object")
    report = OwnerBindingReport(
        roles=(EXPECTED_RESERVED_ROLE, EXPECTED_QUIET_ROLE),
        live_rc_identity_checked=True,
        parameter_identity_checked=True,
        proof_token_identity_checked=True,
        source_owner_authenticated=False,
        source_owner_unverified=True,
        status="conditional_source_owner_unverified",
        limitation=(
            "v1 checks object-identity wiring only; it cannot authenticate that an "
            "arbitrary caller token was issued by the live source owner"
        ),
    )
    return reserved, quiet, report, reserved.proof_token


def prove_reserved_rp_full_support(
    *,
    owner_ports: Sequence[LiveSourceOwnerPort],
    mapping: PhysicalRpMapping,
    matrix_premise: FullSupportMatrixPremise,
) -> ReservedRpSupportWitness:
    """Build a conditional physical ``Rp`` full-support witness.

    The function is pure with respect to its inputs: it performs no I/O,
    imports no producer, evaluates no source function, and does not calculate
    a floating logarithm or exponential.  Invalid/mixed owner or mapping
    records raise ``SupportBridgeError``.  Valid records still produce a
    conditional result because source provenance and the submitted three-bump
    matrix premise cannot be authenticated by this standalone module.
    ``matrix_premise`` is mandatory: support ordering by itself never yields
    ``(0, 0, 1)``.
    """

    _, _, owner_report, token = _validated_owner_ports(owner_ports)
    matrix_report = _validated_matrix_premise(matrix_premise, token)
    if not isinstance(mapping, PhysicalRpMapping):
        raise SupportBridgeError("mapping must be a PhysicalRpMapping")
    _require_mapping_witness(
        mapping.quiet_length_positive,
        EXPECTED_QUIET_LENGTH_PROOF,
        "mapping.quiet_length_positive",
        token,
    )
    _require_mapping_witness(mapping.rc_positive, EXPECTED_RC_POSITIVE_PROOF, "mapping.rc_positive", token)
    _require_mapping_witness(
        mapping.exponential_order,
        EXPECTED_EXP_ORDER_PROOF,
        "mapping.exponential_order",
        token,
    )
    if mapping.rp_formula != EXPECTED_RP_FORMULA:
        raise SupportBridgeError("mapping.rp_formula is not the fixed physical Rp recipe")
    if mapping.two_rc_formula != EXPECTED_TWO_RC_FORMULA:
        raise SupportBridgeError("mapping.two_rc_formula is not the fixed 2Rc recipe")
    if mapping.quiet_length_formula != EXPECTED_QUIET_LENGTH_FORMULA:
        raise SupportBridgeError("mapping.quiet_length_formula is not the fixed quiet-length recipe")
    if FULL_SUPPORT_MARGIN_LOG2 != Fraction(249, 400):
        raise SupportBridgeError("internal support margin constant changed")

    ordinary_z1 = OrdinaryZ1Contract(
        exact_response_row=EXACT_RESPONSE_ROW,
        exact_response_row_Z1=EXACT_RESPONSE_ROW_Z1,
        response_row_constant_in_Z_conditionally=True,
        response_row_constant_in_Z_certified=False,
        row_C1_identity_premise_unverified=True,
        source_function_derivatives_retained=True,
        retained_source_Z1_fields=(
            "theta0_Z",
            "theta1_Z",
            "theta2_Z",
            "amplitude_Z",
            "m15_pressure_Z",
            "mu_Z",
        ),
        statement=(
            "Only after the explicit same-three-bump m15-row/third-matrix-row "
            "identity, its full-Z C1 function identity, and invertibility "
            "premises does the selected inverse row become the constant "
            "(0,0,1), whose own ordinary-Z1 row is (0,0,0); "
            "source-function/product-rule Z1 rows remain required and are not "
            "erased by the row selection."
        ),
    )
    assumptions = (
        EXPECTED_TWO_RC_FORMULA,
        EXPECTED_RP_FORMULA,
        EXPECTED_QUIET_LENGTH_FORMULA + " with " + EXPECTED_QUIET_LENGTH_PROOF,
        EXPECTED_RC_POSITIVE_PROOF,
        EXPECTED_EXP_ORDER_PROOF,
    )
    limitations = (
        "The 249/400*log2 margin is an exact symbolic separation, not a numerical radius value.",
        "Physical support at Rp is conditional on the caller's quiet-length and mapping witnesses.",
        "Source-owner authenticity, source-function equality, and archive/graph identity are unverified.",
        "Support ordering alone does not imply (0,0,1); the same-three-bump complete m15 row identity, its full-Z C1 function identity, and angular-matrix invertibility are separate unverified premises.",
        "The result does not evaluate ordinary-Z1 source jets or an m15 coefficient.",
        "A certified exact row requires a separate trusted owner layer; v1 never promotes the conditional row.",
    )
    return ReservedRpSupportWitness(
        schema=SCHEMA,
        angular_channels=RESERVED_ANGULAR,
        bump_center_labels=RESERVED_BUMP_CENTERS,
        bump_support_upper_log2_coefficient=BUMP_SUPPORT_UPPER_LOG2,
        two_rc_log2_coefficient=TWO_RC_LOG2,
        exact_margin_log2_coefficient=FULL_SUPPORT_MARGIN_LOG2,
        exact_margin_expression="249/400*log2",
        support_before_two_rc_conditionally=True,
        support_before_two_rc_certified=False,
        rp_beyond_two_rc_conditionally=True,
        support_before_rp_conditionally=True,
        full_support_status="conditional_support_and_matrix_source_unverified",
        owner_binding=owner_report,
        matrix_premise=matrix_report,
        physical_mapping_assumptions=assumptions,
        exact_matrix_row=EXACT_RESPONSE_ROW,
        exact_matrix_row_certified=False,
        exact_matrix_row_conditionally_available=True,
        ordinary_z1=ordinary_z1,
        limitations=limitations,
    )


__all__ = [
    "BUMP_SUPPORT_UPPER_LOG2",
    "EXACT_RESPONSE_ROW",
    "EXACT_RESPONSE_ROW_Z1",
    "FULL_SUPPORT_MARGIN_LOG2",
    "FullSupportMatrixPremise",
    "LiveSourceOwnerPort",
    "MatrixPremiseReport",
    "OrdinaryZ1Contract",
    "PhysicalRpMapping",
    "PositiveSourceWitness",
    "RESERVED_ANGULAR",
    "RESERVED_BUMP_CENTERS",
    "ReservedRpSupportWitness",
    "SCHEMA",
    "SourcePremiseWitness",
    "SupportBridgeError",
    "prove_reserved_rp_full_support",
]
