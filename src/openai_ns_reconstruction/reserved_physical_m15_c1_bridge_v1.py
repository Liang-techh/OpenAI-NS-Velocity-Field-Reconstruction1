"""Fail-closed source bridge for the accepted reserved physical ``m15`` row.

This module is deliberately a *consumer* of already authenticated, parsed
archives.  It does not open paths, decompress gzip members, import the old
producers, execute a checker, evaluate an integral, or select an interval
endpoint.  The caller supplies :class:`AuthenticatedArchive` objects whose
raw SHA-256 values were checked before this function is called.

The accepted terminal proof already contains the exact change-of-variables
proof for 25 physical basis columns and 25 ordinary-Z1 matrix rows.  The
bridge below only verifies that saved proof and returns the source graph
handles.  In particular, it never reconstructs the response matrix and it
never turns the unrelated forward-212 packets into a zero residual.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any, Final


SCHEMA: Final[str] = "reserved_physical_m15_c1_bridge_v1"

TERMINAL_ARCHIVE_NAME: Final[str] = (
    "lei_ren_part1_paper_compliant_current_original_n1_terminal_integral_identities.json.gz"
)
TERMINAL_CHECK_NAME: Final[str] = (
    "lei_ren_part1_paper_compliant_current_original_n1_terminal_integral_identities_check.json"
)
OWN_ARCHIVE_NAME: Final[str] = (
    "lei_ren_part1_paper_compliant_current_original_n1_reserved_own_fields.json.gz"
)
FULL_ARCHIVE_NAME: Final[str] = (
    "lei_ren_part1_paper_compliant_current_original_n1_full_defect_linear_functions.json.gz"
)

EXPECTED_SHA256: Final[dict[str, str]] = {
    TERMINAL_ARCHIVE_NAME: "58d90a1d34bbbece45798b3fbd7db75ec751a321eb1740b26f6ac434b76c3b81",
    TERMINAL_CHECK_NAME: "a5c0ca5789b00fb2a3cdeffd07d3f4dff77bef2d169f5b86a8048151afe5efe",
    OWN_ARCHIVE_NAME: "5d43944d68e9e1c85b52493d1e1a4b86a4928209dd6fb97cba9e1d4c63bdad32",
    FULL_ARCHIVE_NAME: "d6b193bf527d4e8cd4afe415b1e941c928f32288b63362783074381ef7480fc1",
}

EXPECTED_SOURCE_FAMILY: Final[dict[str, str]] = {
    "actual_five_defect_family_sha256": "3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894",
    "implicit_source_sha256": "5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae",
    "datum_enclosure_sha256": "dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010",
}


class M15BridgeError(ValueError):
    """Raised whenever the supplied archive set is not the accepted proof."""


@dataclass(frozen=True, slots=True)
class AuthenticatedArchive:
    """A parsed archive plus an external raw-hash attestation.

    ``payload`` must already be parsed by the caller.  This class intentionally
    has no path-opening or decompression method.  ``hash_authenticated`` is an
    assertion made by the caller's archive reader; the bridge additionally
    checks the exact accepted SHA-256 constant and the terminal check's raw
    hash pins.
    """

    name: str
    sha256: str
    payload: Mapping[str, Any]
    hash_authenticated: bool


@dataclass(frozen=True, slots=True)
class GraphHandlePair:
    """Two handles for one genuine C1 source function (value and ordinary Z1)."""

    archive: str
    archive_sha256: str
    role: str
    value: int
    Z1: int


@dataclass(frozen=True, slots=True)
class SignedM15Patch:
    """Actual reserved patch handles in the own-archive prefix.

    The own graph is checked as an exact prefix of the terminal graph, so the
    same canonical integer coordinates are usable in the terminal namespace;
    the archive field remains explicit to prevent silently crossing owners.
    """

    per_dR: GraphHandlePair
    per_ds: GraphHandlePair
    integrand: GraphHandlePair
    history: GraphHandlePair
    endpoint: GraphHandlePair
    terminal: GraphHandlePair


@dataclass(frozen=True, slots=True)
class ReducedM15Sources:
    """Reduced canonical handles consumed by the accepted terminal proof."""

    patch: GraphHandlePair
    terminal: GraphHandlePair
    full_moment: GraphHandlePair


@dataclass(frozen=True, slots=True)
class M15ResponseRow:
    """The five signed physical response columns and their angular factors."""

    columns: tuple[GraphHandlePair, ...]
    angular_third_row: tuple[int, int, int]
    scale: GraphHandlePair
    full_moment: GraphHandlePair


@dataclass(frozen=True, slots=True)
class PhysicalM15Binding:
    """Sign, physical scale and measure retained by the source graph."""

    sign: int
    sign_handle: int
    per_dR_formula: str
    per_ds_formula: str
    measure: str
    radius_handle: int
    log_radius_handle: int
    response_scale_formula: str
    patch_power_formula: str
    scale_value_handle: int
    scale_Z1_handle: int


@dataclass(frozen=True, slots=True)
class ReusedProof:
    """Saved proof facts consumed without redoing their algebra.

    The terminal handles returned by this module are trusted only because the
    authenticated graph prefixes preserve their integer indices and the saved
    terminal receipt certifies the corresponding rows.  This record does not
    claim that the bridge re-evaluated any integral or re-derived a handle.
    """

    basis_rows: int
    basis_Z1_rows: int
    m15_basis_rows: int
    m15_Z1_rows: int
    transported_terminal_rows: tuple[str, ...]
    selected_inverse_rows: int
    terminal_value_and_Z1_zero_certified: bool
    recomputed: bool


@dataclass(frozen=True, slots=True)
class ReservedPhysicalM15C1Bridge:
    """Source-bound result returned by :func:`bind_reserved_physical_m15_c1`."""

    schema: str
    source_family: tuple[tuple[str, str], ...]
    archive_sha256: tuple[tuple[str, str], ...]
    graph_prefix_lengths: tuple[tuple[str, int], ...]
    signed_physical_m15_patch: SignedM15Patch
    reduced_sources: ReducedM15Sources
    response_row: M15ResponseRow
    physical: PhysicalM15Binding
    reused_proof: ReusedProof
    new_forward212_bound: bool
    residual_zero: bool
    limitations: tuple[str, ...]


def _fail(message: str) -> None:
    raise M15BridgeError(message)


def _mapping(value: Any, label: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        _fail(f"{label} must be a parsed mapping")
    return value


def _sequence(value: Any, label: str) -> Sequence[Any]:
    if isinstance(value, (str, bytes, bytearray)) or not isinstance(value, Sequence):
        _fail(f"{label} must be a parsed sequence")
    return value


def _required(mapping: Mapping[str, Any], key: str, label: str) -> Any:
    if key not in mapping:
        _fail(f"missing {label}.{key}")
    return mapping[key]


def _true(mapping: Mapping[str, Any], key: str, label: str) -> None:
    if _required(mapping, key, label) is not True:
        _fail(f"{label}.{key} is not certified true")


def _handle(value: Any, nodes: Sequence[Any], label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        _fail(f"{label} is not a canonical non-negative integer graph handle")
    index = value
    if index < 0 or index >= len(nodes):
        _fail(f"{label} is outside its authenticated graph")
    if not isinstance(nodes[index], Mapping):
        _fail(f"{label} does not address a graph row")
    return index


def _pair(value: Any, nodes: Sequence[Any], label: str, *, archive: AuthenticatedArchive) -> GraphHandlePair:
    pair = _sequence(value, label)
    if len(pair) != 2:
        _fail(f"{label} must contain exactly value and Z1 handles")
    return GraphHandlePair(
        archive=archive.name,
        archive_sha256=archive.sha256,
        role=label,
        value=_handle(pair[0], nodes, f"{label}.value"),
        Z1=_handle(pair[1], nodes, f"{label}.Z1"),
    )


def _row(nodes: Sequence[Any], handle: int, label: str) -> Mapping[str, Any]:
    return _mapping(nodes[_handle(handle, nodes, label)], label)


def _prefix(shorter: Sequence[Any], longer: Sequence[Any], label: str) -> None:
    if len(shorter) >= len(longer) or list(shorter) != list(longer[: len(shorter)]):
        _fail(f"{label} is not an exact graph prefix")


def _auth(archive: AuthenticatedArchive, expected_name: str) -> None:
    if archive.name != expected_name:
        _fail(f"unexpected archive name: {archive.name!r}")
    if archive.hash_authenticated is not True:
        _fail(f"external raw hash attestation missing for {archive.name}")
    expected = EXPECTED_SHA256[expected_name]
    if archive.sha256 != expected:
        _fail(f"raw SHA-256 mismatch for {archive.name}")
    if not isinstance(archive.payload, Mapping):
        _fail(f"{archive.name} was not supplied as a parsed archive")


def _source_family(archive: AuthenticatedArchive) -> tuple[tuple[str, str], ...]:
    family = _mapping(_required(archive.payload, "source_family", archive.name), archive.name + ".source_family")
    if dict(family) != EXPECTED_SOURCE_FAMILY:
        _fail(f"unexpected source family in {archive.name}")
    return tuple((key, family[key]) for key in EXPECTED_SOURCE_FAMILY)


def _check_raw_pins(check: Mapping[str, Any], archives: Sequence[AuthenticatedArchive]) -> None:
    pins = _mapping(_required(check, "input_hashes", "terminal_check"), "terminal_check.input_hashes")
    for archive in archives:
        if pins.get(archive.name) != archive.sha256:
            _fail(f"terminal check does not pin raw hash for {archive.name}")


def _check_terminal_proof(check: Mapping[str, Any]) -> ReusedProof:
    if check.get("all_passed") is not True:
        _fail("terminal check is not an accepted all_passed receipt")
    _true(check, "current_original_n1_actual_terminal_integrals_reduce_to_reserved_response", "terminal_check")
    _true(check, "n1_repaired_own_terminal_five_function_zeros", "terminal_check")
    _true(check, "actual_terminal_integral_functions_zero_not_only_linear_response", "terminal_check")

    checks = _mapping(_required(check, "checks", "terminal_check"), "terminal_check.checks")
    actual = _mapping(_required(checks, "actual_integrals", "terminal_check.checks"), "actual_integrals")
    _true(actual, "passed", "actual_integrals")
    if _required(actual, "basis_rows", "actual_integrals") != 25:
        _fail("accepted terminal proof does not contain 25 physical basis rows")
    basis = _sequence(_required(actual, "actual_translated_physical_basis_integrals", "actual_integrals"), "basis rows")
    if len(basis) != 25 or any(not isinstance(row, Mapping) or row.get("passed") is not True for row in basis):
        _fail("25 physical basis rows are not all saved passed rows")
    m15_basis = [row for row in basis if row.get("moment") == "m15"]
    if len(m15_basis) != 5:
        _fail("saved basis proof does not contain five m15 rows")
    z1_rows = _sequence(_required(actual, "actual_C1_matrix_amplitude_derivative_rows", "actual_integrals"), "matrix Z1 rows")
    if len(z1_rows) != 25 or any(not isinstance(row, str) or not row.endswith("_Z1") for row in z1_rows):
        _fail("saved matrix Z1 proof is not a complete 25-row list")
    if len([row for row in z1_rows if row.startswith("m15_column")]) != 5:
        _fail("saved matrix Z1 proof does not contain five m15 rows")
    _true(actual, "all_five_value_and_Z1_patch_integral_response_identities", "actual_integrals")
    _true(actual, "actual_beta_normalizer_and_support_domain_bound", "actual_integrals")
    _true(actual, "no_numeric_weights_or_enclosure_centers_consumed", "actual_integrals")

    partition = _mapping(_required(checks, "transported_partition", "terminal_check.checks"), "transported_partition")
    _true(partition, "passed", "transported_partition")
    terminal_rows = tuple(_sequence(_required(partition, "actual_transported_terminal_rows", "transported_partition"), "terminal rows"))
    if "m15_value" not in terminal_rows or "m15_Z1" not in terminal_rows:
        _fail("saved transported partition omits m15 value or Z1")
    _true(partition, "m15_negative_weight_not_negated_again", "transported_partition")

    inverse = _mapping(_required(checks, "selected_inverse", "terminal_check.checks"), "selected_inverse")
    _true(inverse, "passed", "selected_inverse")
    _true(inverse, "actual_selected_graph_expressions", "selected_inverse")
    _true(inverse, "accepted_exact_determinant_identities_consumed", "selected_inverse")
    _true(inverse, "original_mu_Rc_Astar_factors_retained", "selected_inverse")
    _true(inverse, "Astar_Z_not_omitted", "selected_inverse")
    inverse_rows = _sequence(
        _required(inverse, "arbitrary_defect_and_amplitude_C1_identity_rows", "selected_inverse"),
        "selected inverse rows",
    )
    if len(inverse_rows) != 10:
        _fail("saved selected inverse proof does not contain ten value/Z1 rows")

    return ReusedProof(
        basis_rows=25,
        basis_Z1_rows=25,
        m15_basis_rows=5,
        m15_Z1_rows=5,
        transported_terminal_rows=terminal_rows,
        selected_inverse_rows=10,
        terminal_value_and_Z1_zero_certified=True,
        recomputed=False,
    )


def _signed_patch(
    terminal: AuthenticatedArchive,
    own: AuthenticatedArchive,
    full: AuthenticatedArchive,
) -> tuple[SignedM15Patch, PhysicalM15Binding, GraphHandlePair]:
    own_nodes = _sequence(_required(own.payload, "exact_graph_nodes", own.name), own.name + ".exact_graph_nodes")
    full_nodes = _sequence(_required(full.payload, "exact_graph_nodes", full.name), full.name + ".exact_graph_nodes")
    own_bindings = _mapping(_required(own.payload, "own_field_bindings", own.name), own.name + ".own_field_bindings")
    patch = _mapping(_required(own_bindings, "patch_bodies", "own_field_bindings"), "patch_bodies")
    per_dr = _mapping(_required(patch, "per_dR", "patch_bodies"), "patch_bodies.per_dR")
    per_ds = _mapping(_required(patch, "per_ds", "patch_bodies"), "patch_bodies.per_ds")
    bodies = _mapping(_required(patch, "bodies", "patch_bodies"), "patch_bodies.bodies")
    histories = _mapping(_required(own_bindings, "patch_histories", own.name), own.name + ".patch_histories")
    endpoint_map = _mapping(_required(own_bindings, "band_endpoint", own.name), own.name + ".band_endpoint")
    terminal_map = _mapping(_required(own_bindings, "terminal", own.name), own.name + ".terminal")

    d_r = _pair(_required(per_dr, "m15", "patch_bodies.per_dR"), own_nodes, "reserved.m15.per_dR", archive=own)
    d_s = _pair(_required(per_ds, "m15", "patch_bodies.per_ds"), own_nodes, "reserved.m15.per_ds", archive=own)
    body = _pair(_required(bodies, "m15", "patch_bodies.bodies"), own_nodes, "reserved.m15.integrand", archive=own)
    history = _pair(_required(histories, "m15", own.name), own_nodes, "reserved.m15.history", archive=own)
    endpoint = _pair(_required(endpoint_map, "m15", own.name), own_nodes, "reserved.m15.endpoint", archive=own)
    terminal_source = _pair(_required(terminal_map, "m15", own.name), own_nodes, "reserved.m15.terminal", archive=own)

    # The graph itself, not a formula string, carries the physical sign and
    # measure.  Both ordinary rows must have the exact -1 factor.
    sign_handles: list[int] = []
    for pair in (d_r,):
        for order, handle in enumerate((pair.value, pair.Z1)):
            row = _row(own_nodes, handle, f"{pair.role}.{order}")
            args = _sequence(row.get("arguments"), f"{pair.role}.{order}.arguments")
            candidates = tuple(
                arg
                for arg in args
                if isinstance(arg, int) and not isinstance(arg, bool)
                and _row(own_nodes, arg, f"{pair.role}.{order}.sign").get("operation") == "exact_rational"
                and _row(own_nodes, arg, f"{pair.role}.{order}.sign").get("numerator") == -1
                and _row(own_nodes, arg, f"{pair.role}.{order}.sign").get("denominator") == 1
            )
            if row.get("operation") != "product" or len(candidates) != 1:
                _fail("reserved m15 per_dR graph does not retain an exact negative sign")
            sign_handles.append(candidates[0])
    if len(set(sign_handles)) != 1:
        _fail("reserved m15 value and Z1 rows do not share the same signed factor")

    full_bindings = _mapping(_required(full.payload, "full_defect_linear_bindings", full.name), full.name + ".full_defect_linear_bindings")
    repair = _mapping(_required(full_bindings, "repair_velocity", "full_defect_linear_bindings"), "repair_velocity")
    radius = _handle(_required(repair, "physical_radius", "repair_velocity"), full_nodes, "physical radius")
    log_radius = _handle(_required(repair, "log_R_over_Rc", "repair_velocity"), full_nodes, "log R/Rc")
    radius_row = _row(full_nodes, radius, "physical radius")
    if radius_row.get("operation") != "analytic_unary" or radius_row.get("name") != "exp":
        _fail("physical radius handle is not the accepted exp graph")

    for order, handle in enumerate((d_s.value, d_s.Z1)):
        row = _row(own_nodes, handle, f"reserved.m15.per_ds.{order}")
        args = _sequence(row.get("arguments"), f"reserved.m15.per_ds.{order}.arguments")
        d_r_handle = (d_r.value, d_r.Z1)[order]
        if row.get("operation") != "product" or radius not in args or d_r_handle not in args:
            _fail("reserved m15 per_ds does not retain R*dR density")

    variable = _handle(_required(patch, "variable", "patch_bodies"), own_nodes, "history variable")
    variable_row = _row(own_nodes, variable, "history variable")
    if variable_row.get("operation") != "bound_variable":
        _fail("reserved history variable is not a bound graph variable")
    for order, handle in enumerate((body.value, body.Z1)):
        row = _row(own_nodes, handle, f"reserved.m15.integrand.{order}")
        if row.get("operation") != "function_substitution" or row.get("expression") != (d_s.value, d_s.Z1)[order]:
            _fail("reserved m15 integrand is not the per_ds substitution")
        if row.get("value") != variable or row.get("Z_independent_substitution") is not True:
            _fail("reserved m15 substitution variable is not source-bound")

    history_variable = variable_row.get("name")
    for order, handle in enumerate((history.value, history.Z1)):
        row = _row(own_nodes, handle, f"reserved.m15.history.{order}")
        if (
            row.get("operation") != "definite_integral"
            or row.get("integrand") != (body.value, body.Z1)[order]
            or row.get("variable") != history_variable
            or row.get("endpoint_independent_of_Z") is not True
            or row.get("exact_function_integral") is not True
            or row.get("numerical_value_not_installed") is not True
            or row.get("ordinary_slow_Z_derivative_order") != order
            or row.get("physical_measure") != "dR=R ds"
            or row.get("paper_moment") != "m15"
        ):
            _fail("reserved m15 history is not the saved physical C1 integral pair")

    response = _mapping(_required(full_bindings, "response_matrix", "full_defect_linear_bindings"), "response_matrix")
    response_pairs = _sequence(_required(response, "m15", "response_matrix"), "response_matrix.m15")
    if len(response_pairs) != 5:
        _fail("saved physical m15 response row does not have five columns")
    columns = tuple(
        _pair(pair, full_nodes, f"response.m15.column{index}", archive=full)
        for index, pair in enumerate(response_pairs)
    )
    zero_nodes = {0}
    if any(columns[index].value not in zero_nodes or columns[index].Z1 not in zero_nodes for index in (0, 1)):
        _fail("m15 axial columns are not the exact inactive zero columns")

    factors: list[int] = []
    value_scales: list[int] = []
    z_scales: list[int] = []
    for index in (2, 3, 4):
        value_row = _row(full_nodes, columns[index].value, f"response.m15.column{index}.value")
        z_row = _row(full_nodes, columns[index].Z1, f"response.m15.column{index}.Z1")
        va = _sequence(value_row.get("arguments"), f"response.m15.column{index}.value.arguments")
        za = _sequence(z_row.get("arguments"), f"response.m15.column{index}.Z1.arguments")
        if value_row.get("operation") != "product" or z_row.get("operation") != "product":
            _fail("m15 angular response columns are not graph products")
        common = tuple(arg for arg in va if arg in za)
        if len(common) != 1:
            _fail("m15 response value/Z1 pair does not share one angular row handle")
        angular = common[0]
        if not isinstance(angular, int) or isinstance(angular, bool):
            _fail("m15 angular row handle is not a graph handle")
        factors.append(_handle(angular, full_nodes, f"response.m15.angular{index - 2}"))
        value_scales.append(_handle(next(arg for arg in va if arg != angular), full_nodes, f"response.m15.scale{index}.value"))
        z_scales.append(_handle(next(arg for arg in za if arg != angular), full_nodes, f"response.m15.scale{index}.Z1"))
    if len(set(factors)) != 3 or len(set(value_scales)) != 1 or len(set(z_scales)) != 1:
        _fail("m15 response row does not contain three angular columns and one common value/Z1 scale")

    scale_value = _row(full_nodes, value_scales[0], "m15 physical scale")
    scale_z = _row(full_nodes, z_scales[0], "m15 physical scale Z1")
    if scale_value.get("operation") != "product" or scale_z.get("operation") != "product":
        _fail("m15 scale handles are not actual product graph expressions")
    if not any(_row(full_nodes, arg, "m15 scale sign").get("operation") == "negative" for arg in scale_value.get("arguments", ()) if isinstance(arg, int) and not isinstance(arg, bool)):
        _fail("m15 physical response scale lost its negative sign")
    if not any(_row(full_nodes, arg, "m15 scale Z1 sign").get("operation") == "negative" for arg in scale_z.get("arguments", ()) if isinstance(arg, int) and not isinstance(arg, bool)):
        _fail("m15 physical response scale Z1 lost its negative sign")

    full_moment = _pair(_required(full_bindings, "moments", "full_defect_linear_bindings")["m15"], full_nodes, "full.moment.m15", archive=full)
    return (
        SignedM15Patch(per_dR=d_r, per_ds=d_s, integrand=body, history=history, endpoint=endpoint, terminal=terminal_source),
        PhysicalM15Binding(
            sign=-1,
            sign_handle=sign_handles[0],
            per_dR_formula="-theta0*theta",
            per_ds_formula="R*(-theta0*theta)",
            measure="dR=R ds",
            radius_handle=radius,
            log_radius_handle=log_radius,
            response_scale_formula="-Astar*exp(-mu*log(2Rc))",
            patch_power_formula="theta0=Astar*exp(-(1/2+mu)*(log(2)+log(Rc)+s))",
            scale_value_handle=value_scales[0],
            scale_Z1_handle=z_scales[0],
        ),
        full_moment,
    )


def bind_reserved_physical_m15_c1(
    *,
    terminal_archive: AuthenticatedArchive,
    terminal_check: AuthenticatedArchive,
    reserved_own_archive: AuthenticatedArchive,
    full_defect_archive: AuthenticatedArchive,
) -> ReservedPhysicalM15C1Bridge:
    """Bind the accepted reserved physical m15 C1 source row.

    The function consumes saved graph/check facts only.  It has no ``forward``
    argument by design: ``new_forward212_bound`` and ``residual_zero`` are
    always false until a later, separately source-proved comparison adapter is
    added.
    """

    archives = (terminal_archive, terminal_check, reserved_own_archive, full_defect_archive)
    for archive, name in zip(archives, (TERMINAL_ARCHIVE_NAME, TERMINAL_CHECK_NAME, OWN_ARCHIVE_NAME, FULL_ARCHIVE_NAME)):
        _auth(archive, name)

    family = _source_family(terminal_archive)
    for archive in (terminal_check, reserved_own_archive, full_defect_archive):
        if _source_family(archive) != family:
            _fail(f"source family mismatch in {archive.name}")

    check = terminal_check.payload
    _check_raw_pins(check, (terminal_archive, reserved_own_archive, full_defect_archive))
    terminal_nodes = _sequence(_required(terminal_archive.payload, "exact_graph_nodes", terminal_archive.name), "terminal graph")
    own_nodes = _sequence(_required(reserved_own_archive.payload, "exact_graph_nodes", reserved_own_archive.name), "reserved own graph")
    full_nodes = _sequence(_required(full_defect_archive.payload, "exact_graph_nodes", full_defect_archive.name), "full defect graph")
    _prefix(full_nodes, own_nodes, "full defect graph -> reserved own graph")
    _prefix(own_nodes, terminal_nodes, "reserved own graph -> terminal graph")

    # Prefix equality is the only handle-identity relation established here:
    # the same non-negative integer points at the same source row in the
    # larger graph.  The terminal receipt is consumed as saved proof below;
    # no integral, derivative, or handle is recomputed by this bridge.

    proof = _check_terminal_proof(check)
    patch, physical, full_moment = _signed_patch(terminal_archive, reserved_own_archive, full_defect_archive)
    terminal_bindings = _mapping(_required(terminal_archive.payload, "terminal_identity_bindings", terminal_archive.name), "terminal_identity_bindings")
    own_bindings = _mapping(_required(reserved_own_archive.payload, "own_field_bindings", reserved_own_archive.name), "own_field_bindings")
    full_bindings = _mapping(_required(full_defect_archive.payload, "full_defect_linear_bindings", full_defect_archive.name), "full_defect_linear_bindings")
    reduced = _mapping(_required(terminal_bindings, "reduced_patch", "terminal_identity_bindings"), "reduced_patch")
    reduced_terminal = _mapping(_required(terminal_bindings, "reduced_terminal", "terminal_identity_bindings"), "reduced_terminal")
    full_moments = _mapping(_required(full_bindings, "moments", "full_defect_linear_bindings"), "moments")
    reduced_sources = ReducedM15Sources(
        patch=_pair(_required(reduced, "m15", "reduced_patch"), terminal_nodes, "reduced.patch.m15", archive=terminal_archive),
        terminal=_pair(_required(reduced_terminal, "m15", "reduced_terminal"), terminal_nodes, "reduced.terminal.m15", archive=terminal_archive),
        full_moment=_pair(_required(full_moments, "m15", "full moments"), full_nodes, "full.moment.m15", archive=full_defect_archive),
    )
    response_matrix = _mapping(_required(full_bindings, "response_matrix", "full_defect_linear_bindings"), "response_matrix")
    response_pairs = _sequence(_required(response_matrix, "m15", "response_matrix"), "response_matrix.m15")
    response_columns = tuple(
        _pair(pair, full_nodes, f"response.m15.column{index}", archive=full_defect_archive)
        for index, pair in enumerate(response_pairs)
    )
    factors = []
    for index in (2, 3, 4):
        row = _row(full_nodes, response_columns[index].value, f"response.m15.column{index}")
        args = _sequence(row.get("arguments"), f"response.m15.column{index}.arguments")
        zrow = _row(full_nodes, response_columns[index].Z1, f"response.m15.column{index}.Z1")
        zargs = _sequence(zrow.get("arguments"), f"response.m15.column{index}.Z1.arguments")
        shared = tuple(arg for arg in args if arg in zargs)
        if len(shared) != 1:
            _fail("cannot identify the saved angular third-row factors")
        factors.append(_handle(shared[0], full_nodes, f"angular third row {index - 2}"))
    scale_pair = GraphHandlePair(
        archive=full_defect_archive.name,
        archive_sha256=full_defect_archive.sha256,
        role="response.m15.common physical scale",
        value=physical.scale_value_handle,
        Z1=physical.scale_Z1_handle,
    )
    return ReservedPhysicalM15C1Bridge(
        schema=SCHEMA,
        source_family=family,
        archive_sha256=tuple((archive.name, archive.sha256) for archive in archives),
        graph_prefix_lengths=((FULL_ARCHIVE_NAME, len(full_nodes)), (OWN_ARCHIVE_NAME, len(own_nodes)), (TERMINAL_ARCHIVE_NAME, len(terminal_nodes))),
        signed_physical_m15_patch=patch,
        reduced_sources=reduced_sources,
        response_row=M15ResponseRow(columns=response_columns, angular_third_row=tuple(factors), scale=scale_pair, full_moment=full_moment),
        physical=physical,
        reused_proof=proof,
        new_forward212_bound=False,
        residual_zero=False,
        limitations=(
            "This result is only the accepted reserved source graph; forward212 has no binding here.",
            "No interval, name, range, metadata, or source-family label is used as function identity.",
            "Terminal handle identity is established only by the verified exact graph prefix and saved proof; no handle or integral is recomputed here.",
            "A later forward comparison must provide genuine source C1 handles, physical measure, endpoints, scale, sign, and ordinary Z1 equality.",
        ),
    )


__all__ = [
    "AuthenticatedArchive",
    "GraphHandlePair",
    "M15BridgeError",
    "M15ResponseRow",
    "PhysicalM15Binding",
    "ReducedM15Sources",
    "ReservedPhysicalM15C1Bridge",
    "SignedM15Patch",
    "bind_reserved_physical_m15_c1",
]
