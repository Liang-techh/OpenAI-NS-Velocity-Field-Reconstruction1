"""Inert conditional cylindrical-chart metric/coordinate-operator kernel.

This draft is deliberately inert at import time.  It consumes only parsed,
hash-pinned source packets supplied by the caller and never opens source
files, imports project modules, calls a provider, materializes ``N``, or
replays a native field definition.  The only executable symbolic work is the
new universal cylindrical-chart algebra in :func:`derive`, after source and
branch validation has completed.

The native packets are provenance data.  They bind two branch-local frames
to the same universal geometry definition; they do not instantiate a native
field pullback or any phase identity.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
from typing import Any

__all__ = ["derive"]


INPUT_SCHEMA = "native_O3_cylindrical_pullback_input_v1"
RESULT_SCHEMA = "original_native_O3_cylindrical_pullback"
COMPONENT = "native_O3_cylindrical_pullback"

# These are source-data pins only.  No path is opened by this module.
EXPECTED_SOURCE_SHA256 = {
    "accepted_result": "b541342e47bed12ed9b7394facef236d03b4828b6436776739245715b93e0df9",
    "rm_trace": "901ed6eca9505e20eb7b3819ddc6d87d8e162974f37132655c05679b4c028f3f",
    "phase_trace": "8255f24dc44c6d0a47aad3c854b7a197cbf2fec809cbd965e715df3f5f3c49fb",
    "overlap_audit": "d163d3997ad3b7b85d5dc5b4eba066a0d7b4ab398d3cc608b5c3377cf2ea7b13",
}
PARENT_REPORT_SHA256 = "7786f54ac15bf400cba235613ee5c2c89fd4a2556cb36eda3d1f1cfa58b43a93"

EXPECTED_BRANCHES = ("background", "modified")
EXPECTED_COORDINATE_NODE = 241536
EXPECTED_RADIUS_NODE = 241787
EXPECTED_CHART = "O3_transition"

_FLAG_NAMES = (
    "conditional_native_chart_metric_identity_proved",
    "native_branch_total_derivative_policy_enforced",
    "native_cylindrical_volume_density_defined",
    "actual_native_field_pullback_instantiated",
    "native_Rm_Z_zero_proved",
    "actual_physical_time_map_installed",
    "phase_identity_proved",
    "physical_layer_completed",
    "full_goal_complete",
    "caps_deleted",
    "numeric_closure",
    "temporal_recursion",
)

_COUNT_NAMES = (
    "branch_frame_bindings",
    "new_algebra_family_calls",
    "new_mathematical_proofs",
    "jacobian_identity_rows",
    "inverse_coordinate_rows",
    "cartesian_commutator_rows",
    "scalar_laplacian_rows",
    "branch_operator_manifests",
    "policy_rejections",
    "numeric_calls",
    "source_value_provider_calls",
    "old_math_calls",
    "field_definition_replays",
    "N_materializations",
)


class _InputFailure(Exception):
    """Internal validation failure; never exposed as a traceback."""


class _DiagnosticFailure(Exception):
    """Internal callback failure used to stop all later work."""

    def __init__(self, event: str, cause: BaseException) -> None:
        super().__init__(str(cause))
        self.event = event
        self.cause = cause


class _PolicyRefusal(Exception):
    """A source-backed forbidden promotion was actually refused."""

    def __init__(self, row: Mapping[str, Any]) -> None:
        super().__init__(str(row.get("policy", "policy refusal")))
        self.row = dict(row)


def _new_counts() -> dict[str, int]:
    return {name: 0 for name in _COUNT_NAMES}


def _false_flags() -> dict[str, bool]:
    return {name: False for name in _FLAG_NAMES}


def _source_fingerprint_view(saved_sources: Mapping[str, Mapping[str, Any]]) -> dict[str, dict[str, Any]]:
    """Return provenance without echoing the potentially huge literals."""

    return {
        name: {
            "raw_sha256": str(wrapper["raw_sha256"]).lower(),
            "literal_consumed": True,
            "data_only": True,
        }
        for name, wrapper in saved_sources.items()
    }


def _failure_result(
    *,
    counts: Mapping[str, int],
    payload: Mapping[str, Any],
    events: list[dict[str, Any]],
    failure: Mapping[str, Any],
) -> dict[str, Any]:
    """Build the stable failure contract without flattening counters."""

    return {
        "schema": RESULT_SCHEMA,
        "component": COMPONENT,
        "passed": False,
        "runtime_model_metadata": None,
        "metadata_only": False,
        "source_only": True,
        "conditional": True,
        "physical_layer_completed": False,
        "full_goal_complete": False,
        "actual_counts": dict(counts),
        "new_flags": _false_flags(),
        "payload": dict(payload),
        "diagnostics": {
            "events": list(events),
            "event_count": len(events),
            "callback_failure_stops_later_work": failure.get("kind") == "diagnostic_callback_failure",
        },
        "failure": dict(failure),
    }


def _success_result(
    *,
    counts: Mapping[str, int],
    payload: Mapping[str, Any],
    events: list[dict[str, Any]],
) -> dict[str, Any]:
    flags = _false_flags()
    flags["conditional_native_chart_metric_identity_proved"] = True
    flags["native_branch_total_derivative_policy_enforced"] = True
    flags["native_cylindrical_volume_density_defined"] = True
    return {
        "schema": RESULT_SCHEMA,
        "component": COMPONENT,
        "passed": True,
        "runtime_model_metadata": None,
        "metadata_only": False,
        "source_only": True,
        "conditional": True,
        "physical_layer_completed": False,
        "full_goal_complete": False,
        "actual_counts": dict(counts),
        "new_flags": flags,
        "payload": dict(payload),
        "diagnostics": {
            "events": list(events),
            "event_count": len(events),
            "callback_failure_stops_later_work": True,
        },
        "failure": None,
    }


def derive(source: dict, diagnostic: Callable[[dict], Any] | None = None) -> dict:
    """Derive the new conditional cylindrical coordinate/operator identities.

    Parameters
    ----------
    source:
        A parsed JSON-like mapping following
        ``native_O3_cylindrical_pullback_input_v1``.  Four saved packets are
        supplied as ``{literal, raw_sha256}`` wrappers; the parent report is
        a separate hash-only provenance field.  Literals are retained as
        opaque data; this function never reads their paths.
    diagnostic:
        Optional callback.  It receives one plain JSON-like mapping for each
        completed row.  A callback exception returns the completed prefix and
        stops all later work.
    """

    counts = _new_counts()
    events: list[dict[str, Any]] = []
    payload: dict[str, Any] = {
        "source_fingerprints": {},
        "branch_frame_bindings": [],
        "branch_operator_manifests": [],
        "policy_rejections": [],
        "witness_rows": {
            "jacobian_identity": [],
            "inverse_coordinate": [],
            "cartesian_commutator": [],
            "scalar_laplacian": [],
        },
    }

    if diagnostic is not None and not callable(diagnostic):
        return _failure_result(
            counts=counts,
            payload=payload,
            events=events,
            failure={
                "kind": "invalid_diagnostic",
                "message": "diagnostic must be callable or None",
                "stopped": True,
            },
        )

    def emit(event: str, row: Mapping[str, Any]) -> None:
        record = {"event": event, **dict(row)}
        events.append(record)
        if diagnostic is None:
            return
        try:
            diagnostic(record)
        except BaseException as exc:  # callback boundary is intentionally broad
            raise _DiagnosticFailure(event, exc) from exc

    try:
        if not isinstance(source, Mapping):
            raise _InputFailure("source must be a mapping")
        if source.get("schema") != INPUT_SCHEMA:
            raise _InputFailure(f"source.schema must be {INPUT_SCHEMA!r}")

        saved_sources = source.get("saved_sources")
        if not isinstance(saved_sources, Mapping):
            raise _InputFailure("source.saved_sources must be a mapping of literal wrappers")
        required_source_names = tuple(EXPECTED_SOURCE_SHA256)
        if set(saved_sources) != set(required_source_names):
            raise _InputFailure(
                "source.saved_sources must contain exactly " + ", ".join(required_source_names)
            )
        checked_sources: dict[str, Mapping[str, Any]] = {}
        for name in required_source_names:
            wrapper = saved_sources.get(name)
            if not isinstance(wrapper, Mapping):
                raise _InputFailure(f"saved_sources.{name} must be a {{literal, raw_sha256}} wrapper")
            if "literal" not in wrapper or "raw_sha256" not in wrapper:
                raise _InputFailure(f"saved_sources.{name} must expose literal and raw_sha256")
            observed = str(wrapper["raw_sha256"]).lower()
            expected = EXPECTED_SOURCE_SHA256[name]
            if observed != expected:
                raise _InputFailure(f"saved_sources.{name}.raw_sha256 does not match its accepted pin")
            checked_sources[name] = wrapper

        parent_report_raw_sha256 = str(source.get("parent_report_raw_sha256", "")).lower()
        if parent_report_raw_sha256 != PARENT_REPORT_SHA256:
            raise _InputFailure("parent_report_raw_sha256 must match the accepted report pin")
        accepted_literal = checked_sources["accepted_result"]["literal"]
        if not isinstance(accepted_literal, Mapping):
            raise _InputFailure("saved_sources.accepted_result.literal must be parsed result data")
        accepted_payload = accepted_literal.get("payload")
        if not isinstance(accepted_payload, Mapping):
            raise _InputFailure("accepted_result.literal.payload is required for manifest binding")
        saved_manifests = accepted_payload.get("native_bundle_manifests")
        if not isinstance(saved_manifests, (list, tuple)) or len(saved_manifests) != 2:
            raise _InputFailure("accepted_result.literal.payload.native_bundle_manifests must contain two rows")
        saved_by_id: dict[str, Mapping[str, Any]] = {}
        for saved_manifest in saved_manifests:
            if not isinstance(saved_manifest, Mapping):
                raise _InputFailure("saved native bundle manifests must be mappings")
            saved_id = saved_manifest.get("manifest_id")
            if not isinstance(saved_id, str) or saved_id in saved_by_id:
                raise _InputFailure("saved native bundle manifest ids must be unique strings")
            saved_by_id[saved_id] = saved_manifest

        branches = source.get("branch_manifests")
        if not isinstance(branches, (list, tuple)) or len(branches) != 2:
            raise _InputFailure("exactly two branch_manifests are required")
        validated_manifests: list[dict[str, Any]] = []
        seen_branches: set[str] = set()
        seen_ids: set[str] = set()
        for manifest in branches:
            if not isinstance(manifest, Mapping):
                raise _InputFailure("each branch manifest must be a mapping")
            manifest_id = manifest.get("manifest_id")
            branch = manifest.get("branch")
            frame = manifest.get("frame")
            radius_node = manifest.get("radius_node")
            radius_definition = manifest.get("radius_definition")
            if not isinstance(manifest_id, str) or not manifest_id:
                raise _InputFailure("each branch manifest needs a nonempty manifest_id")
            if branch not in EXPECTED_BRANCHES:
                raise _InputFailure("branch manifests must be exactly background and modified")
            if branch in seen_branches or manifest_id in seen_ids:
                raise _InputFailure("branch manifest ids and branch names must be unique")
            if not isinstance(frame, Mapping):
                raise _InputFailure(f"{branch} manifest.frame must be a mapping")
            if not isinstance(radius_definition, Mapping):
                raise _InputFailure(f"{branch} manifest.radius_definition must be the saved source mapping")
            saved_manifest = saved_by_id.get(manifest_id)
            if saved_manifest is None:
                raise _InputFailure(f"{branch} manifest_id is absent from the accepted saved result")
            saved_branch = saved_manifest.get("branch")
            saved_frame = saved_manifest.get("frame")
            saved_radius_definition = saved_manifest.get("radius_definition")
            if saved_branch != branch or not isinstance(saved_frame, Mapping):
                raise _InputFailure(f"{branch} manifest branch/frame does not match the accepted saved result")
            if not isinstance(saved_radius_definition, Mapping):
                raise _InputFailure(f"{branch} saved radius_definition source is missing")
            if dict(frame) != dict(saved_frame):
                raise _InputFailure(f"{branch} frame must match the accepted saved result literally")
            if dict(radius_definition) != dict(saved_radius_definition):
                raise _InputFailure(f"{branch} radius_definition must match the accepted source literally")
            if frame.get("chart") != EXPECTED_CHART:
                raise _InputFailure(f"{branch} frame must identify chart {EXPECTED_CHART!r}")
            if frame.get("coordinate_node") != EXPECTED_COORDINATE_NODE:
                raise _InputFailure(f"{branch} frame must bind coordinate node {EXPECTED_COORDINATE_NODE}")
            if radius_node != frame.get("radius_node") or radius_node != EXPECTED_RADIUS_NODE:
                raise _InputFailure(f"{branch} manifest must bind saved radius node {EXPECTED_RADIUS_NODE}")
            if saved_radius_definition.get("root_node") != EXPECTED_RADIUS_NODE:
                raise _InputFailure(f"{branch} saved radius_definition root node is not {EXPECTED_RADIUS_NODE}")
            if bool(manifest.get("rebuilt")) or bool(frame.get("rebuilt")):
                raise _InputFailure("rebuilt frames cannot bind the native branch geometry")
            field_branches = manifest.get("field_branches")
            if field_branches is not None:
                try:
                    field_branch_set = set(field_branches)
                except TypeError as exc:
                    raise _InputFailure("field_branches must be an iterable branch selector") from exc
                if field_branch_set != {branch}:
                    raise _InputFailure("a branch manifest cannot blend background and modified field references")
            seen_branches.add(branch)
            seen_ids.add(manifest_id)
            validated_manifests.append({
                "manifest_id": manifest_id,
                "branch": branch,
                "frame": dict(saved_frame),
                "radius_node": radius_node,
                "radius_definition": dict(saved_radius_definition),
                "radius_handle": {
                    "owner": saved_frame.get(
                        "native_Rm_owner", "self.Rm = self.leading.inlet['radius']"
                    )
                },
            })
        if tuple(sorted(seen_branches)) != tuple(sorted(EXPECTED_BRANCHES)):
            raise _InputFailure("branch manifests must contain one background and one modified manifest")

        requested_context = source.get("requested_context", {})
        if not isinstance(requested_context, Mapping):
            raise _InputFailure("requested_context must be a mapping when supplied")

        payload["source_fingerprints"] = _source_fingerprint_view(checked_sources)
        payload["source_fingerprints"]["parent_report"] = {
            "raw_sha256": parent_report_raw_sha256,
            "literal_consumed": False,
            "data_only": True,
        }
        payload["source_contract"] = {
            "schema": INPUT_SCHEMA,
            "literal_wrappers": {
                name: {"keys": ["literal", "raw_sha256"], "data_only": True}
                for name in required_source_names
            },
            "source_files_opened": False,
            "project_modules_imported": False,
            "old_math_called": False,
            "source_values_consumed": False,
            "N_materialized": False,
        }
        emit(
            "source_validation",
            {
                "schema": INPUT_SCHEMA,
                "saved_source_names": list(required_source_names),
                "saved_source_hashes": {
                    name: checked_sources[name]["raw_sha256"] for name in required_source_names
                },
                "parent_report_raw_sha256": parent_report_raw_sha256,
                "parent_report_consumed": False,
                "literal_wrappers_validated": True,
                "data_only": True,
                "branch_manifest_count": 2,
            },
        )

        # Bind the two saved frames before entering the new universal algebra.
        for manifest in validated_manifests:
            branch = manifest["branch"]
            binding = {
                "manifest_id": manifest["manifest_id"],
                "branch": branch,
                "frame": manifest["frame"],
                "radius_node": manifest["radius_node"],
                "radius_definition": manifest["radius_definition"],
                "native_Rm_maps_to": "rho(Z)",
                "native_Rm_owner": manifest["radius_handle"].get(
                    "owner", "self.Rm = self.leading.inlet['radius']"
                ),
                "rho_Z_status": "unknown; no zero specialization accepted",
                "ordinary_Z_ports": "referenced saved selectors only; not rebuilt",
                "field_values": False,
                "mixed_derivatives_available": False,
            }
            payload["branch_frame_bindings"].append(binding)
            counts["branch_frame_bindings"] += 1
            emit("branch_frame_binding", binding)

        # The generic symbolic engine is intentionally imported only now,
        # after the root entry has frozen source and branch provenance.
        counts["new_algebra_family_calls"] += 1
        emit(
            "algebra_family_entry",
            {
                "family": "conditional_native_cylindrical_chart_metric",
                "coordinates": ["x", "alpha", "Z"],
                "physical_coordinates": ["X", "Y", "Zphysical"],
                "rho_assumption": "rho(Z) > 0 and C2",
                "c_assumption": "c real constant",
                "F_assumption": "F(x, alpha, Z) arbitrary C2",
                "alpha_role": "physical azimuth; distinct from loop angle psi and inertial theta",
                "native_Rm_mapping": "Rm handle maps to rho(Z), with rho_Z left unknown",
                "source_only": True,
            },
        )
        import sympy as sp  # noqa: PLC0415 - guarded by the source/root entry above
        payload["symbolic_backend"] = {
            "library": "sympy",
            "version": str(getattr(sp, "__version__", "unknown")),
            "imported_after_source_validation": True,
        }

        symbolic_counts = {
            "diff_api_calls": 0,
            "determinant_api_calls": 0,
            "reduction_wrapper_calls": 0,
            "planned_diff_api_calls": 70,
            "planned_determinant_api_calls": 1,
            "planned_reduction_wrapper_calls": 14,
        }

        def diff(expr: Any, *variables: Any) -> Any:
            symbolic_counts["diff_api_calls"] += 1
            return sp.diff(expr, *variables)

        def determinant(matrix: Any) -> Any:
            symbolic_counts["determinant_api_calls"] += 1
            return matrix.det()

        def reduce_expr(expr: Any) -> Any:
            symbolic_counts["reduction_wrapper_calls"] += 1
            return sp.simplify(sp.trigsimp(expr))

        x, Z, alpha, c = sp.symbols("x Z alpha c", real=True)
        rho_function = sp.Function("rho")
        rho = rho_function(Z)
        F = sp.Function("F")(x, alpha, Z)
        R = rho * sp.exp(c + x)

        payload["symbolic_call_counts"] = symbolic_counts
        rho_Z = diff(rho, Z)
        a = rho_Z / rho
        coordinates = (x, alpha, Z)
        embedding = (
            R * sp.cos(alpha),
            R * sp.sin(alpha),
            Z,
        )
        J = sp.Matrix([[diff(component, coordinate) for coordinate in coordinates] for component in embedding])
        B = sp.Matrix(
            [
                [sp.cos(alpha) / R, sp.sin(alpha) / R, -a],
                [-sp.sin(alpha) / R, sp.cos(alpha) / R, sp.Integer(0)],
                [sp.Integer(0), sp.Integer(0), sp.Integer(1)],
            ]
        )

        determinant_value = determinant(J)
        determinant_expected = R**2
        determinant_residual = reduce_expr(determinant_value - determinant_expected)
        determinant_row = {
            "row_id": "jacobian_determinant",
            "embedding": ["R*cos(alpha)", "R*sin(alpha)", "Zphysical=Z"],
            "determinant": str(sp.sstr(determinant_value)),
            "expected": str(sp.sstr(determinant_expected)),
            "residual": str(sp.sstr(determinant_residual)),
            "proved": determinant_residual == sp.Integer(0),
            "not_R_cubed": True,
            "numeric": False,
        }
        payload["determinant"] = determinant_row
        payload["witness_rows"]["jacobian_identity"].append(determinant_row)
        counts["jacobian_identity_rows"] += 1
        emit("jacobian_identity_row", determinant_row)

        inverse_product = J * B
        inverse_matrix = inverse_product - sp.eye(3)
        inverse_rows: list[dict[str, Any]] = []
        for row_index in range(3):
            for column_index in range(3):
                residual = reduce_expr(inverse_matrix[row_index, column_index])
                row = {
                    "row_id": f"J_times_B_minus_I_{row_index}_{column_index}",
                    "row": row_index,
                    "column": column_index,
                    "J_row": [str(sp.sstr(J[row_index, k])) for k in range(3)],
                    "B_column": [str(sp.sstr(B[k, column_index])) for k in range(3)],
                    "lhs_J_times_B": str(sp.sstr(inverse_product[row_index, column_index])),
                    "expected": "1" if row_index == column_index else "0",
                    "unreduced_residual": str(sp.sstr(inverse_matrix[row_index, column_index])),
                    "residual": str(sp.sstr(residual)),
                    "proved": residual == sp.Integer(0),
                    "numeric": False,
                }
                inverse_rows.append(row)
                payload["witness_rows"]["inverse_coordinate"].append(row)
                counts["inverse_coordinate_rows"] += 1
                emit("inverse_coordinate_row", row)

        def P(index: int, expression: Any) -> Any:
            return sum(
                (B[coordinate_index, index] * diff(expression, coordinates[coordinate_index])
                 for coordinate_index in range(3)),
                sp.Integer(0),
            )

        operator_labels = ("P_X", "P_Y", "P_Z")
        operator_coefficients = {
            label: {
                "physical_coordinate": physical,
                "definition": "sum_k B[k,i] D_k F",
                "coefficients": {
                    coordinate_name: str(sp.sstr(B[coordinate_index, operator_index]))
                    for coordinate_index, coordinate_name in enumerate(("x", "alpha", "Z"))
                },
                "function": "F(x, alpha, Z)",
            }
            for operator_index, (label, physical) in enumerate(zip(operator_labels, ("X", "Y", "Z")))
        }
        payload["full_operator_coefficients"] = {
            "coordinate_order": ["x", "alpha", "Z"],
            "physical_order": ["X", "Y", "Z"],
            "R": str(sp.sstr(R)),
            "a": str(sp.sstr(a)),
            "B_rows": [[str(sp.sstr(B[row, column])) for column in range(3)] for row in range(3)],
            "operators": operator_coefficients,
        }

        commutator_rows: list[dict[str, Any]] = []
        for left_index, right_index in ((0, 1), (0, 2), (1, 2)):
            left = P(left_index, P(right_index, F))
            right = P(right_index, P(left_index, F))
            residual = reduce_expr(left - right)
            row = {
                "row_id": f"[{operator_labels[left_index]},{operator_labels[right_index]}]F",
                "operators": [operator_labels[left_index], operator_labels[right_index]],
                "left": str(sp.sstr(left)),
                "right": str(sp.sstr(right)),
                "residual": str(sp.sstr(residual)),
                "proved": residual == sp.Integer(0),
                "function": "arbitrary C2 F",
                "numeric": False,
            }
            commutator_rows.append(row)
            payload["witness_rows"]["cartesian_commutator"].append(row)
            counts["cartesian_commutator_rows"] += 1
            emit("cartesian_commutator_row", row)

        laplace_value = sum((P(index, P(index, F)) for index in range(3)), sp.Integer(0))
        F_xx = diff(F, x, 2)
        F_alphaalpha = diff(F, alpha, 2)
        F_ZZ = diff(F, Z, 2)
        F_xZ = diff(F, x, Z)
        F_x = diff(F, x)
        a_Z = diff(a, Z)
        laplace_expected = (
            (F_xx + F_alphaalpha) / R**2
            + F_ZZ
            - 2 * a * F_xZ
            + a**2 * F_xx
            - a_Z * F_x
        )
        laplace_residual = reduce_expr(laplace_value - laplace_expected)
        laplace_row = {
            "row_id": "cartesian_scalar_laplacian",
            "operator": "P_X(P_X F) + P_Y(P_Y F) + P_Z(P_Z F)",
            "computed": str(sp.sstr(laplace_value)),
            "expected": str(sp.sstr(laplace_expected)),
            "residual": str(sp.sstr(laplace_residual)),
            "proved": laplace_residual == sp.Integer(0),
            "a": str(sp.sstr(a)),
            "a_Z": str(sp.sstr(a_Z)),
            "function": "arbitrary C2 F",
            "numeric": False,
        }
        payload["witness_rows"]["scalar_laplacian"].append(laplace_row)
        counts["scalar_laplacian_rows"] += 1
        emit("scalar_laplacian_row", laplace_row)

        all_rows = [
            payload["determinant"]["proved"],
            *(row["proved"] for row in payload["witness_rows"]["inverse_coordinate"]),
            *(row["proved"] for row in payload["witness_rows"]["cartesian_commutator"]),
            payload["witness_rows"]["scalar_laplacian"][0]["proved"],
        ]
        if not all(all_rows):
            return _failure_result(
                counts=counts,
                payload=payload,
                events=events,
                failure={
                    "kind": "symbolic_identity_not_proved",
                    "message": "At least one of the 14 required witness rows was not zero.",
                    "stopped": True,
                },
            )
        counts["new_mathematical_proofs"] += 1

        payload["premises"] = {
            "rho": {
                "expression": "rho(Z)",
                "positive": True,
                "smoothness": "C2",
                "rho_Z": "symbolic and unresolved",
                "rho_Z_zero_proved": False,
            },
            "c": {"real_constant": True},
            "R": "rho(Z)*exp(c + x)",
            "alpha": {
                "role": "physical azimuth",
                "distinct_from": ["loop_angle psi", "inertial theta"],
            },
            "F": {"expression": "F(x, alpha, Z)", "smoothness": "arbitrary C2"},
            "native_field_pullback": {
                "instantiated": False,
                "mixed_derivatives_claimed": False,
                "source_function_values": False,
            },
            "native_Rm": {
                "owner": "self.Rm = self.leading.inlet['radius']",
                "positive": "source premise only; no value evaluated",
                "Z_derivative": "unknown",
                "maps_to": "rho(Z)",
            },
            "phase_context": {
                "opaque_source_recipe": True,
                "fractional_part_differentiated": False,
                "phase_identity_proved": False,
            },
        }
        payload["volume_density"] = {
            "density": str(sp.sstr(R**2)),
            "measure": "R^2 dx dalpha dZ",
            "conditional_on": "native cylindrical embedding X=R*cos(alpha), Y=R*sin(alpha), Zphysical=Z",
            "actual_physical_time_volume": False,
            "L2_certification": False,
        }

        # Both branch manifests refer to the one universal algebra only as a
        # definition.  No branch field row is requested or replayed.
        for binding in payload["branch_frame_bindings"]:
            branch_manifest = {
                "manifest_id": binding["manifest_id"],
                "branch": binding["branch"],
                "frame": binding["frame"],
                "radius_node": binding["radius_node"],
                "radius_definition": binding["radius_definition"],
                "universal_geometry": "conditional_native_cylindrical_chart_metric",
                "universal_definition_only": True,
                "native_Rm_maps_to": "rho(Z)",
                "same_universal_geometry_as": [item["manifest_id"] for item in payload["branch_frame_bindings"]],
                "ordinary_Z_ports": "saved selectors referenced; not rebuilt",
                "mixed_derivatives_available": False,
                "actual_field_pullback_values": False,
                "phase_context": "opaque source recipe; no phase proof",
            }
            payload["branch_operator_manifests"].append(branch_manifest)
            counts["branch_operator_manifests"] += 1
            emit("branch_operator_manifest", branch_manifest)

        def context_guard(candidate: Mapping[str, Any]) -> None:
            """Raise only when the actual saved source supports a refusal."""

            policy = candidate["policy"]
            if policy == "held_phase_NativeDerivative_not_total_x":
                literal = candidate["source_literal"]
                definitions = literal.get("derivative_operator_definitions_already_declared")
                held_rows = [
                    row for row in (definitions if isinstance(definitions, (list, tuple)) else [])
                    if isinstance(row, Mapping)
                    and row.get("name") == "NativeDerivative"
                    and "held" in str(row.get("phase_policy", "")).lower()
                ]
                if candidate["attempted"] and held_rows:
                    raise _PolicyRefusal({
                        "policy": policy,
                        "status": "rejected",
                        "requested": candidate["requested"],
                        "guard_attempted": True,
                        "source_basis": held_rows[0].get("phase_policy"),
                        "source_ref": "phase_trace.derivative_operator_definitions_already_declared",
                        "source_raw_sha256": candidate["source_raw_sha256"],
                        "accepted_action": "keep held-phase rows opaque; do not differentiate fractional_part",
                    })
            elif policy == "source_physical_measure1_not_3D_volume_weight":
                frame = candidate["source_frame"]
                if candidate["attempted"] and frame.get("physical_measure") in (1, "1"):
                    raise _PolicyRefusal({
                        "policy": policy,
                        "status": "rejected",
                        "requested": candidate["requested"],
                        "guard_attempted": True,
                        "source_basis": "physical_measure=1 is dlogR/dx only",
                        "source_ref": "accepted_result.literal.payload.native_bundle_manifests[*].frame.physical_measure",
                        "source_raw_sha256": candidate["source_raw_sha256"],
                        "accepted_action": "do not promote source measure one to a 3D volume weight",
                    })
            elif policy == "rho_Z_zero_specialization_without_accepted_certificate":
                literal = candidate["source_literal"]
                conclusion = literal.get("dependency_conclusion", {})
                ordinary = conclusion.get("ordinary_Z_dependency_of_Rm", {})
                has_zero_row = bool(ordinary.get("saved_Rm_DZ_zero_row"))
                has_positive_certificate = bool(ordinary.get("saved_Rm_positive_quotient_certificate"))
                if candidate["attempted"] and not has_zero_row:
                    raise _PolicyRefusal({
                        "policy": policy,
                        "status": "rejected",
                        "requested": candidate["requested"],
                        "guard_attempted": True,
                        "source_basis": ordinary.get(
                            "status", "no direct saved Rm ordinary-Z derivative or positive certificate"
                        ),
                        "positive_certificate_present": has_positive_certificate,
                        "source_ref": "rm_trace.dependency_conclusion.ordinary_Z_dependency_of_Rm",
                        "source_raw_sha256": candidate["source_raw_sha256"],
                        "accepted_action": "retain rho_Z symbolic and unresolved",
                    })

        policy_candidates = (
            {
                "policy": "held_phase_NativeDerivative_not_total_x",
                "requested": bool(requested_context.get("held_phase_rows_as_total_x", False)),
                "attempted": bool(requested_context.get("held_phase_rows_as_total_x", False)),
                "source_literal": checked_sources["phase_trace"]["literal"],
                "source_raw_sha256": EXPECTED_SOURCE_SHA256["phase_trace"],
            },
            {
                "policy": "source_physical_measure1_not_3D_volume_weight",
                "requested": bool(requested_context.get("physical_measure_as_3d_volume", False)),
                "attempted": bool(requested_context.get("physical_measure_as_3d_volume", False)),
                "source_frame": validated_manifests[0]["frame"],
                "source_raw_sha256": EXPECTED_SOURCE_SHA256["accepted_result"],
            },
            {
                "policy": "rho_Z_zero_specialization_without_accepted_certificate",
                "requested": bool(requested_context.get("rho_Z_zero_specialization", False)),
                "attempted": bool(requested_context.get("rho_Z_zero_specialization", False)),
                "source_literal": checked_sources["rm_trace"]["literal"],
                "source_raw_sha256": EXPECTED_SOURCE_SHA256["rm_trace"],
            },
        )
        for candidate in policy_candidates:
            try:
                context_guard(candidate)
            except _PolicyRefusal as refusal:
                policy_row = refusal.row
                payload["policy_rejections"].append(policy_row)
                counts["policy_rejections"] += 1
                emit("policy_rejection", policy_row)
            else:
                return _failure_result(
                    counts=counts,
                    payload=payload,
                    events=events,
                    failure={
                        "kind": "context_policy_not_refused",
                        "policy": candidate["policy"],
                        "message": "A source-backed forbidden context was accepted; later work stopped.",
                        "stopped": True,
                    },
                )

        payload["native_nonclaims"] = {
            "actual_native_field_pullback_instantiated": False,
            "native_Rm_Z_zero_proved": False,
            "actual_physical_time_map_installed": False,
            "phase_identity_proved": False,
            "physical_layer_completed": False,
            "full_goal_complete": False,
            "caps_deleted": False,
            "numeric_closure": False,
            "temporal_recursion": False,
        }

        return _success_result(
            counts=counts,
            payload=payload,
            events=events,
        )

    except _DiagnosticFailure as exc:
        return _failure_result(
            counts=counts,
            payload=payload,
            events=events,
            failure={
                "kind": "diagnostic_callback_failure",
                "event": exc.event,
                "message": str(exc.cause),
                "stopped": True,
            },
        )
    except _InputFailure as exc:
        return _failure_result(
            counts=counts,
            payload=payload,
            events=events,
            failure={
                "kind": "source_validation_failure",
                "message": str(exc),
                "stopped": True,
            },
        )
    except Exception as exc:
        return _failure_result(
            counts=counts,
            payload=payload,
            events=events,
            failure={
                "kind": "unexpected_symbolic_exception",
                "message": str(exc),
                "exception_type": type(exc).__name__,
                "stopped": True,
            },
        )

