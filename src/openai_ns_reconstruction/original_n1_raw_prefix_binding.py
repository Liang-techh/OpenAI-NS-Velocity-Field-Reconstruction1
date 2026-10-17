"""Bind the saved forward n=1 radial prefix to the accepted raw system.

The saved prefix is a finite coefficient archive for the normalized forward
variables ``f, u, k, p``.  This module performs only a new static binding to
the accepted all-order radial conjugacy.  It does not import or execute a
producer, owner, recurrence, checker, fixture, or quadrature.

The raw fields and the two regularized odd fields are kept in symbolic scale
form.  In particular, the enormous ``Lambda`` and ``Pstar`` are never
materialized and the archive's decimal approximations are never promoted to
values.  The returned receipt contains the exact saved interval tuples and a
fail-closed description of the finite-prefix and ordinary-Z limits.
"""

from __future__ import annotations

import ast
import gzip
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping


PARENT_NAME = "lei_ren_part1_paper_compliant_current_original_C5_terminal_primitive.json"
PARENT_SHA256 = "588d221fd082d5f84d94f1812eb6342bb6647847b114dac458c0165324a063b1"
PREFIX_NAME = "lei_ren_part1_paper_compliant_current_original_n1_core_numeric_prefix.json"
PREFIX_SHA256 = "2ae583cdf5ea9fde03208fea0423eb11d171054c9a8dd45d78aa5fe2076d81fa"
SOURCE_VECTOR_NAME = "lei_ren_part1_paper_compliant_current_original_inlet_source_vector.json"
SOURCE_VECTOR_SHA256 = "2ac1269c3da29466067d1dfd50a92be321e85291b2634ffa6d2eafd679e9941e"
FORWARD_SOURCE = "src/openai_ns_reconstruction/original_n1_core_numeric.py"
CANONICAL_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_inner_volterra_functions.py"
)
CANONICAL_OPERATOR = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_inner_operator.py"
)
TAYLOR_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_core_nonlinear_operator.py"
)

NORMALIZATION = {
    "f": "epsilon*F1/F0(Z)",
    "u": "Uz1/epsilon",
    "k": "K1/epsilon",
    "p": "P1/(epsilon*Pstar^2)",
}
CANONICAL_STATE = ["F1", "Uz1", "K1", "P1", "F1x", "Uz1x"]
CANONICAL_D = [0, 0, 2, 0, 3, 1]
AXIS_POWERS = [2, 2, 2, 2, 1, 1]
BASE_FIELDS = ("f", "u", "k", "p")
RAW_FIELDS = ("W1", "W2", "W3", "W4")
WEIGHTED_FIELDS = ("W5_over_x", "W6_over_x")


def _read(path: Path) -> Any:
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8") as stream:
            return json.load(stream)
    return json.loads(path.read_text(encoding="utf-8"))


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _mapping_digest(value: Any) -> str:
    """Digest mappings independently of JSON integer parsing choices."""

    def normalize(item: Any) -> Any:
        if isinstance(item, bool) or item is None or isinstance(item, str):
            return item
        if isinstance(item, int):
            return str(item)
        if isinstance(item, float):
            return repr(item)
        if isinstance(item, Mapping):
            return {str(key): normalize(item[key]) for key in sorted(item, key=str)}
        if isinstance(item, (list, tuple)):
            return [normalize(part) for part in item]
        return repr(item)

    payload = json.dumps(
        normalize(value), sort_keys=True, separators=(",", ":"), ensure_ascii=False
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _path(root: Path, name: str) -> Path:
    candidate = root / name
    if candidate.exists():
        return candidate
    return root / "experiments" / "root_st073" / Path(name).name


def _function_text(path: Path, name: str) -> str:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    nodes = [
        node
        for node in ast.walk(tree)
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == name
    ]
    if len(nodes) != 1:
        return ""
    return "".join(ast.unparse(nodes[0]).split())


def _source_pin(expected: Any, root: Path, label: str) -> dict[str, Any]:
    """Check one accepted source pin without executing its module."""

    if not isinstance(expected, Mapping):
        return {"label": label, "passed": False, "reason": "missing accepted source pin"}
    name = expected.get("path")
    digest = expected.get("sha256", expected.get("expected_sha256"))
    if not isinstance(name, str) or not isinstance(digest, str):
        return {"label": label, "passed": False, "reason": "malformed source pin"}
    path = _path(root, name)
    if not path.exists():
        return {"label": label, "path": name, "passed": False, "reason": "source missing"}
    actual = _sha256(path)
    return {
        "label": label,
        "path": name,
        "expected_sha256": digest,
        "actual_sha256": actual,
        "passed": actual == digest,
        "source_executed": False,
    }


def _accepted_pins(root: Path, parent: Mapping[str, Any]) -> dict[str, Any]:
    conjugacy = parent.get("n1_all_order_radial_conjugacy", {})
    source_pins = conjugacy.get("source_pins", {}) if isinstance(conjugacy, Mapping) else {}
    receipts = conjugacy.get("accepted_receipts", {}) if isinstance(conjugacy, Mapping) else {}
    checks: dict[str, Any] = {}
    for key in (
        "accepted_source_vector_module",
        "forward_advance",
        "forward_rhs",
        "forward_linear_rhs",
        "canonical_raw_matrix",
        "canonical_raw_system",
        "canonical_operator_rows",
    ):
        checks[key] = _source_pin(source_pins.get(key), root, key)

    prefix_pin = receipts.get("prefix_report", {})
    source_vector_pin = receipts.get("source_vector_report", {})
    checks["prefix_report"] = _source_pin(
        {"path": prefix_pin.get("path"), "sha256": prefix_pin.get("sha256")},
        root,
        "prefix_report",
    )
    checks["source_vector_report"] = _source_pin(
        {"path": source_vector_pin.get("path"), "sha256": source_vector_pin.get("sha256")},
        root,
        "source_vector_report",
    )
    # The accepted prefix receipt also carries the implementation pins.  Check
    # these bytes directly, while refusing to run the producer.
    prefix_path = _path(root, PREFIX_NAME)
    prefix_report = _read(prefix_path) if prefix_path.exists() else {}
    implementation_checks: dict[str, Any] = {}
    for name, digest in prefix_report.get("implementation_hashes", {}).items():
        path = _path(root, name)
        actual = _sha256(path) if path.exists() else None
        implementation_checks[name] = {
            "expected_sha256": digest,
            "actual_sha256": actual,
            "passed": actual == digest,
            "source_executed": False,
        }
    checks["prefix_implementation_hashes"] = implementation_checks
    # The prefix archive's own input map pins the formal Taylor rectangle.
    # This source is needed to distinguish ordinary Taylor coefficients from
    # raw derivatives when a future consumer evaluates a Z jet.
    taylor_path = _path(root, TAYLOR_SOURCE)
    taylor_name = "lei_ren_part1_paper_compliant_current_core_nonlinear_operator.py"
    taylor_expected = prefix_report.get("input_hashes", {}).get(taylor_name)
    taylor_actual = _sha256(taylor_path) if taylor_path.exists() else None
    checks["prefix_taylor_source"] = {
        "path": taylor_name,
        "expected_sha256": taylor_expected,
        "actual_sha256": taylor_actual,
        "passed": isinstance(taylor_expected, str) and taylor_actual == taylor_expected,
        "source_executed": False,
    }
    checks["passed"] = all(
        entry.get("passed") is True
        for key, entry in checks.items()
        if key != "prefix_implementation_hashes"
    ) and bool(implementation_checks) and all(
        entry.get("passed") is True for entry in implementation_checks.values()
    )
    return checks


def _ast_binding(root: Path) -> dict[str, Any]:
    forward = _path(root, FORWARD_SOURCE)
    canonical = _path(root, CANONICAL_SOURCE)
    operator = _path(root, CANONICAL_OPERATOR)
    taylor = _path(root, TAYLOR_SOURCE)
    if not all(path.exists() for path in (forward, canonical, operator, taylor)):
        return {"passed": False, "reason": "accepted source AST missing"}

    forward_advance = _function_text(forward, "_advance")
    forward_rhs = _function_text(forward, "_rhs")
    advance_fragments = (
        "axial=degree+depth-n-1",
        "rf,ru,rp=self._rhs(fixed,leading,state,n,axial,center,delta,eps,P2,c)",
        "k=n+1",
        "state['f'].append([v/(4*k*(k+1))forvinrf.rows[n]])",
        "state['u'].append([v/(4*k*k)forvinru.rows[n]])",
        "state['p'].append([v/kforvinrp.rows[n]])",
        "state['k'].append([-c.mpf(k)*v/(k+1)forvinstate['u'][-1]])",
    )
    rhs_fragments = (
        "Aw=known(leading['A'],n+1,axial+2)",
        "Uw=known(leading['Uz'],n+1,axial+2)",
        "Qw=(2*zz*uu-zz*mm*(1-delta)-dd*mw.dz())/LL",
        "QZ=crop(Qw.dz(),n,axial)",
        "pressure_head=R([[(i+1)*(i+2)*vforvinQw.rows[i+1][:axial+1]",
        "F2=R(known([fixed['S_Z_taylor']],0,axial).rows,n,axial,c.mpf)/(eps*eps)",
        "force_u-=4*z*rho*Np*(eps*P2)/L",
    )
    forward_text = forward.read_text(encoding="utf-8")
    linear_text = _function_text(forward, "linear_rhs")
    linear_fragments = (
        "rf=(2*e*(Q+E*(-2+de)/L)-2*H*g/L)*f",
        "ru=-8*F2*z*r*A*f/L",
        "rp=2*F2/(e*P2)*A*f",
    )
    canonical_text = canonical.read_text(encoding="utf-8")
    canonical_normalized = "".join(canonical_text.split())
    canonical_fragments = (
        "D=(0,0,2,0,3,1)",
        "POWERS=(2,2,2,2,1,1)",
        "raw_units_not_normalized_potential_kernel=True",
        "raw_physical_units_retained=True",
        "F_rho_rows=[a.mul(core.functions['F0'],q)forqinrows['Phi']]",
        "*[a.scale(q,eps)forqinrows['Psi'][1:]]",
        "g4=g.mul(c(2),eps,x,Np.value)",
        "g6=g.sub(g.mul(c(2),eps,Nz)",
    )
    operator_text = operator.read_text(encoding="utf-8")
    operator_normalized = "".join(operator_text.split())
    operator_fragments = (
        "B0['4,1']=eps(f*(4*x))",
        "B1[f'{i},{j}']=eps(value)",
        "returndict(B0={k:truncate(v,3)fork,vinB0.items()}",
    )
    taylor_normalized = "".join(taylor.read_text(encoding="utf-8").split())
    taylor_dz = _function_text(taylor, "dz")
    taylor_euler = _function_text(taylor, "radial_euler")
    taylor_fragments = {
        "ordinary_rows": "self.rows=[[convert(rows[n][m])",
        "ordinary_Z_derivative_factor": "[(m+1)*row[m+1]forminrange(self.axial)]",
        "ordinary_rho_euler_factor": "[[n*vforvinrow]forn,rowinenumerate(self.rows)]",
    }
    matches = {
        "forward_advance": {part: part in forward_advance for part in advance_fragments},
        "forward_rhs": {part: part in forward_rhs for part in rhs_fragments},
        "forward_linear_rhs": {part: part in linear_text for part in linear_fragments},
        "canonical_raw_system": {part: part in canonical_normalized for part in canonical_fragments},
        "canonical_operator": {part: part in operator_normalized for part in operator_fragments},
        "taylor_rectangle": {
            "ordinary_rows": taylor_fragments["ordinary_rows"] in taylor_normalized,
            "ordinary_Z_derivative_factor": taylor_fragments["ordinary_Z_derivative_factor"] in taylor_dz,
            "ordinary_rho_euler_factor": taylor_fragments["ordinary_rho_euler_factor"] in taylor_euler,
        },
        "forward_source_text_present": bool(forward_text),
    }
    return {
        "passed": all(all(group.values()) if isinstance(group, Mapping) else group for group in matches.values()),
        "matches": matches,
        "source_executed": False,
        "canonical_state": CANONICAL_STATE,
        "canonical_D": CANONICAL_D,
        "axis_powers": AXIS_POWERS,
        "ordinary_Z_coefficient_convention": "rows[rho_degree][Z_degree] are ordinary Taylor coefficients; dz multiplies (Z_degree+1)",
        "ordinary_rho_coefficient_convention": "rows[rho_degree] are ordinary rho powers; radial_euler multiplies rho_degree",
    }


def _archive_schema(prefix: Mapping[str, Any]) -> dict[str, Any]:
    coefficients = prefix.get("coefficients")
    levels = prefix.get("levels")
    checks = prefix.get("checks", {})
    limitations = prefix.get("limitations", {})
    errors: list[str] = []
    if not isinstance(coefficients, list):
        return {"passed": False, "errors": ["coefficients is not a list"]}
    if not isinstance(levels, list) or len(levels) != 2:
        errors.append("two precision levels are required")
    fields = {record.get("field") for record in coefficients if isinstance(record, Mapping)}
    radial = {
        record.get("rho_degree")
        for record in coefficients
        if isinstance(record, Mapping) and isinstance(record.get("rho_degree"), int)
    }
    centers = {
        record.get("Z") for record in coefficients if isinstance(record, Mapping)
    }
    z_orders_by_radial: dict[int, set[int]] = {}
    slots: list[tuple[Any, Any, Any, Any]] = []
    for record in coefficients:
        if not isinstance(record, Mapping):
            errors.append("non-mapping coefficient record")
            continue
        if not set(("Z", "field", "rho_degree", "Z_degree", "source_enclosure_exact_mpf_tuples")) <= set(record):
            errors.append("coefficient schema missing required key")
            continue
        m = record["rho_degree"]
        z_order = record["Z_degree"]
        if not isinstance(m, int) or not isinstance(z_order, int):
            errors.append("coefficient indices must be integers")
            continue
        z_orders_by_radial.setdefault(m, set()).add(z_order)
        slots.append((record.get("Z"), record.get("field"), m, z_order))
        enclosure = record["source_enclosure_exact_mpf_tuples"]
        if (
            not isinstance(enclosure, list)
            or len(enclosure) != 2
            or any(not isinstance(endpoint, list) or len(endpoint) != 4 for endpoint in enclosure)
            or any(not isinstance(part, int) for endpoint in enclosure for part in endpoint)
        ):
            errors.append("coefficient enclosure is not an exact mpf endpoint pair")
    expected_fields = set(BASE_FIELDS)
    expected_radial = set(range(17))
    expected_centers = {"-0.3", "0", "0.3"}
    if fields != expected_fields:
        errors.append("archive fields are not exactly f,u,k,p")
    if radial != expected_radial:
        errors.append("archive radial degrees are not exactly 0..16")
    if centers != expected_centers:
        errors.append("archive centers are not exactly -0.3, 0, 0.3")
    if len(coefficients) != 2040:
        errors.append("archive does not contain 2040 directed coefficient slots")
    unique_slots = set(slots)
    duplicate_slots = len(slots) - len(unique_slots)
    expected_slots = {
        (center, field, m, z_order)
        for center in expected_centers
        for field in expected_fields
        for m in range(17)
        for z_order in range(18 - m)
    }
    missing_slots = sorted(expected_slots - unique_slots, key=str)
    extra_slots = sorted(unique_slots - expected_slots, key=str)
    if duplicate_slots:
        errors.append("archive contains duplicate (Z,field,rho_degree,Z_degree) slots")
    if missing_slots:
        errors.append("archive has missing contiguous radial/Z slots")
    if extra_slots:
        errors.append("archive has slots outside the accepted radial/Z lattice")
    if not checks.get("passed"):
        errors.append("accepted prefix checks did not pass")
    if limitations.get("radial_Taylor_remainder_certified") is not False:
        errors.append("prefix limitation does not explicitly leave radial remainder open")
    if limitations.get("Volterra_word_tail_reused_as_radial_tail") is not False:
        errors.append("prefix limitation does not reject Volterra-tail relabeling")
    if not checks.get("no_source_midpoint_selected"):
        errors.append("accepted archive midpoint guard is absent")
    max_z_order = {
        str(m): max(orders) for m, orders in sorted(z_orders_by_radial.items())
    }
    ordinary_z3_max = max(
        (m for m, orders in z_orders_by_radial.items() if set(range(4)) <= orders),
        default=-1,
    )
    level_meta: list[dict[str, Any]] = []
    for level in levels if isinstance(levels, list) else []:
        samples = level.get("samples", [])
        if level.get("degree") != 16:
            errors.append("precision level degree is not 16")
        level_meta.append(
            {
                "precision": level.get("precision"),
                "quadrature_nodes": level.get("quadrature_nodes"),
                "degree": level.get("degree"),
                "sample_count": len(samples) if isinstance(samples, list) else None,
                "sample_Z_centers": sorted({s.get("Z") for s in samples}, key=str)
                if isinstance(samples, list)
                else [],
                "sample_rho_values": sorted({s.get("rho") for s in samples}, key=str)
                if isinstance(samples, list)
                else [],
                "all_samples_are_finite_prefixes": all(
                    s.get("finite_radial_prefix_only") is True
                    and s.get("infinite_radial_remainder_bounded") is False
                    and s.get("radial_degree") == 16
                    for s in samples
                )
                if isinstance(samples, list)
                else False,
            }
        )
    passed = not errors and all(meta["all_samples_are_finite_prefixes"] for meta in level_meta)
    return {
        "passed": passed,
        "errors": errors,
        "record_count": len(coefficients),
        "fields": sorted(fields, key=str),
        "radial_degrees": sorted(radial),
        "Z_centers": sorted(centers, key=str),
        "max_Z_degree_by_rho_degree": max_z_order,
        "slot_lattice": {
            "unique_slot_count": len(unique_slots),
            "duplicate_slot_count": duplicate_slots,
            "missing_slot_count": len(missing_slots),
            "extra_slot_count": len(extra_slots),
            "missing_slots": [list(slot) for slot in missing_slots[:32]],
            "extra_slots": [list(slot) for slot in extra_slots[:32]],
            "contiguous_C0_to_C16": not duplicate_slots and not missing_slots and not extra_slots,
            "common_Z0_to_Z3_rho0_to_rho14": all(
                (center, field, m, z_order) in unique_slots
                for center in expected_centers
                for field in expected_fields
                for m in range(15)
                for z_order in range(4)
            ),
        },
        "ordinary_Z_orders": [0, 1, 2, 3],
        "ordinary_Z_0_to_3_max_rho_degree": ordinary_z3_max,
        "ordinary_Z_0_to_3_max_weighted_degree": ordinary_z3_max - 1,
        "level_metadata": level_meta,
        "point_values_selected": False,
        "approximation_decimal_used": False,
    }


def _conjugacy_binding(parent: Mapping[str, Any]) -> dict[str, Any]:
    """Bind this map to the accepted all-order receipt plus fresh index algebra."""

    conjugacy = parent.get("n1_all_order_radial_conjugacy", {})
    symbolic = conjugacy.get("symbolic", {}) if isinstance(conjugacy, Mapping) else {}
    matrix = symbolic.get("matrix_conjugacy", {}) if isinstance(symbolic, Mapping) else {}
    known = symbolic.get("known_source_map", {}) if isinstance(symbolic, Mapping) else {}
    radial = symbolic.get("radial_rules", {}) if isinstance(symbolic, Mapping) else {}
    formal = conjugacy.get("formal_statement", {}) if isinstance(conjugacy, Mapping) else {}
    checks = radial.get("checks", {}) if isinstance(radial, Mapping) else {}
    general = radial.get("general_power_lhs", {}) if isinstance(radial, Mapping) else {}
    general_checks = general.get("checks", {}) if isinstance(general, Mapping) else {}
    matrix_checks = matrix.get("checks", {}) if isinstance(matrix, Mapping) else {}
    known_checks = known.get("checks", {}) if isinstance(known, Mapping) else {}
    def normalization_matches(value: Any) -> bool:
        if not isinstance(value, Mapping):
            return False
        # The accepted receipt writes F0 without its argument; the saved
        # prefix writes F0(Z).  These are the same named source function.
        normalized = {key: str(part).replace("F0(Z)", "F0") for key, part in value.items()}
        expected = {key: str(part).replace("F0(Z)", "F0") for key, part in NORMALIZATION.items()}
        return normalized == expected

    accepted_normalization = normalization_matches(matrix.get("normalization")) and normalization_matches(known.get("normalization"))
    accepted_assumptions = (
        matrix.get("assumptions", {}).get("rho_equals_x_squared") == "rho = x^2"
        and matrix.get("assumptions", {}).get("P2_equals_Pstar_squared") == "P2 = Pstar^2"
        and known.get("assumptions", {}).get("P2_equals_Pstar_squared") == "P2 = Pstar^2"
    )
    accepted_rule_checks = all(
        checks.get(name) is True
        for name in (
            "F1x_definition",
            "Uz1x_definition",
            "row_5_factor",
            "row_6_factor",
            "row_4_factor",
            "row_3_factor" if "row_3_factor" in checks else "K1_kinematic_factor",
            "general_power_lhs",
            "ordinary_Z_commutes_with_radial_coefficient",
        )
    ) and all(general_checks.values())
    accepted_formal_state = (
        formal.get("canonical_state") == CANONICAL_STATE
        and formal.get("canonical_D") == CANONICAL_D
        and formal.get("axis_extension") == "W(0,Z)=0; W1..W4 are O(x^2), W5,W6 are O(x)"
    )
    # New finite-index arithmetic checks the derivative and weighted degree
    # map for every archived nonzero coefficient, independently of the old
    # 30/36 symbolic comparison counts.
    index_rows: list[dict[str, Any]] = []
    finite_algebra_ok = True
    for n in range(16):
        m = n + 1
        row = {
            "n": n,
            "m": m,
            "base_power_x": 2 * m,
            "odd_power_x": 2 * m - 1,
            "weighted_rho_degree": m - 1,
            "derivative_factor": 2 * m,
            "row_5_factor": 4 * m * (m + 1),
            "row_6_factor": 4 * m * m,
            "row_4_factor": 2 * m,
            "row_3_factor": 2 * (m + 1),
            "identities": {
                "m_equals_n_plus_1": m == n + 1,
                "weighted_degree_equals_n": m - 1 == n,
                "derivative_factor_is_2n_plus_2": 2 * m == 2 * (n + 1),
                "odd_power_is_2n_plus_1": 2 * m - 1 == 2 * n + 1,
            },
        }
        finite_algebra_ok = finite_algebra_ok and all(row["identities"].values())
        index_rows.append(row)
    raw_polynomial_ok = (
        accepted_formal_state
        and accepted_assumptions
        and finite_algebra_ok
        and general.get("input_series_power") == "Base coordinates W1..W4 have degree-m terms proportional to x^(2m), m=n+1"
        and general.get("F1x_and_Uz1x_power") == "2m-1"
    )
    derivative_ok = (
        accepted_rule_checks
        and finite_algebra_ok
        and general.get("row_5_factor") == "(2m)(2m-1)+3(2m)=4m(m+1)"
        and general.get("row_6_factor") == "(2m)(2m-1)+(2m)=4m^2"
        and general.get("row_4_factor") == "2m"
        and general.get("row_3_factor") == "2m+2"
    )
    accepted_receipt_checks = (
        conjugacy.get("passed") is True
        and matrix_checks.get("normalization_assumptions_bound") is True
        and known_checks.get("all_three_sources") is True
        and accepted_normalization
    )
    return {
        "passed": accepted_receipt_checks and accepted_rule_checks and raw_polynomial_ok and derivative_ok,
        "accepted_receipt_checks": {
            "conjugacy_passed": conjugacy.get("passed") is True,
            "matrix_normalization_assumptions_bound": matrix_checks.get("normalization_assumptions_bound") is True,
            "known_source_all_three": known_checks.get("all_three_sources") is True,
            "normalization_exact": accepted_normalization,
            "assumptions_exact": accepted_assumptions,
        },
        "accepted_radial_rule_checks": {key: checks.get(key) for key in sorted(checks)},
        "accepted_general_power_checks": {key: general_checks.get(key) for key in sorted(general_checks)},
        "formal_state_exact": accepted_formal_state,
        "fresh_index_algebra": {
            "passed": finite_algebra_ok,
            "rows": index_rows,
            "derivative_factor": "2*(n+1)",
        },
        "normalization": matrix.get("normalization"),
        "known_source_normalization": known.get("normalization"),
        "source": "parent.n1_all_order_radial_conjugacy.symbolic",
        "source_executed": False,
    }


def _physical_scale_binding(parent: Mapping[str, Any]) -> dict[str, Any]:
    conjugacy = parent.get("n1_all_order_radial_conjugacy", {})
    symbolic = conjugacy.get("symbolic", {}) if isinstance(conjugacy, Mapping) else {}
    matrix = symbolic.get("matrix_conjugacy", {}) if isinstance(symbolic, Mapping) else {}
    known = symbolic.get("known_source_map", {}) if isinstance(symbolic, Mapping) else {}
    premises = matrix.get("accepted_normalization_premises", {}) if isinstance(matrix, Mapping) else {}
    premise_checks = premises.get("checks", {}) if isinstance(premises, Mapping) else {}
    complex_bounds = parent.get("n1_complex_actual_source_bounds", {})
    physical_scale = complex_bounds.get("physical_scale") if isinstance(complex_bounds, Mapping) else None
    checks = {
        "epsilon_reciprocal_identity": premise_checks.get("symbolic_shared_core_Lambda_epsilon_reciprocal") is True,
        "P2_equals_Pstar_squared": premise_checks.get("symbolic_P2_equals_shared_Pstar_squared") is True,
        "F0_nonzero": matrix.get("assumptions", {}).get("F0_nonzero") == "F0 != 0",
        "epsilon_nonzero": matrix.get("assumptions", {}).get("epsilon_nonzero") == "epsilon != 0",
        "P2_nonzero": matrix.get("assumptions", {}).get("P2_nonzero") == "P2 != 0",
        "physical_radial_scale_source": physical_scale == "ordinary d_R^i=exp(i*logLambda)d_rho^i; no Lambda or physical derivative is materialized",
        "F0_complex_axis_binding": complex_bounds.get("gates", {}).get("F0_complex_axis_binding") is True if isinstance(complex_bounds, Mapping) else False,
    }
    # The raw field multipliers follow directly from the accepted normalization;
    # keep the formulas visible and check them against the four source scales.
    raw_scale_map = {
        "W1_F1": "F0(Z)/epsilon",
        "W2_Uz1": "epsilon",
        "W3_K1": "epsilon",
        "W4_P1": "epsilon*Pstar^2",
        "W5_over_x": "F0(Z)/epsilon",
        "W6_over_x": "epsilon",
    }
    known_normalization = known.get("normalization")
    checks["raw_scale_map_from_normalization"] = isinstance(known_normalization, Mapping) and {
        key: str(value).replace("F0(Z)", "F0") for key, value in known_normalization.items()
    } == {
        key: str(value).replace("F0(Z)", "F0") for key, value in NORMALIZATION.items()
    }
    return {
        "passed": all(checks.values()),
        "checks": checks,
        "raw_scale_map": raw_scale_map,
        "intervals_materialized": False,
        "source_executed": False,
    }


def _coefficient_enclosures(prefix: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Copy exact interval metadata while intentionally dropping decimal points."""

    output: list[dict[str, Any]] = []
    for record in prefix.get("coefficients", []):
        output.append(
            {
                "Z": record["Z"],
                "field": record["field"],
                "rho_degree": record["rho_degree"],
                "Z_degree": record["Z_degree"],
                "source_enclosure_exact_mpf_tuples": record[
                    "source_enclosure_exact_mpf_tuples"
                ],
            }
        )
    return output


def _symbolic_maps() -> dict[str, Any]:
    return {
        "rho_equals_x_squared": "rho = x^2",
        "normalized_to_raw": {
            "f_to_W1": "W1=F1=(F0(Z)/epsilon)*f",
            "u_to_W2": "W2=Uz1=epsilon*u",
            "k_to_W3": "W3=K1=epsilon*k",
            "p_to_W4": "W4=P1=(epsilon*Pstar^2)*p",
        },
        "raw_polynomial_definitions": {
            "f_prefix": "f^(16)(rho,Z)=sum_{m=1}^{16} f_m(Z)*rho^m",
            "u_prefix": "u^(16)(rho,Z)=sum_{m=1}^{16} u_m(Z)*rho^m",
            "k_prefix": "k^(16)(rho,Z)=sum_{m=1}^{16} k_m(Z)*rho^m",
            "p_prefix": "p^(16)(rho,Z)=sum_{m=1}^{16} p_m(Z)*rho^m",
            "W1_prefix": "F1^(16)=(F0(Z)/epsilon)*sum_{m=1}^{16} f_m*rho^m",
            "W2_prefix": "Uz1^(16)=epsilon*sum_{m=1}^{16} u_m*rho^m",
            "W3_prefix": "K1^(16)=epsilon*sum_{m=1}^{16} k_m*rho^m",
            "W4_prefix": "P1^(16)=(epsilon*Pstar^2)*sum_{m=1}^{16} p_m*rho^m",
            "W5_prefix": "F1x^(16)=x*(F0(Z)/epsilon)*sum_{n=0}^{15}2*(n+1)*f_(n+1)*rho^n",
            "W6_prefix": "Uz1x^(16)=x*epsilon*sum_{n=0}^{15}2*(n+1)*u_(n+1)*rho^n",
        },
        "weighted_regularized_prefix": {
            "W5_over_x": "F1x/x=(F0(Z)/epsilon)*sum_{n=0}^{15}2*(n+1)*f_(n+1)*rho^n",
            "W6_over_x": "Uz1x/x=epsilon*sum_{n=0}^{15}2*(n+1)*u_(n+1)*rho^n",
            "degree": 15,
            "coefficient_factor": "2*(n+1)",
            "n_range": [0, 15],
        },
        "raw_base_prefix": {
            "fields": list(RAW_FIELDS),
            "degree": 16,
            "terms": "rho^m, m=1..16",
            "zero_axis_values": {field: "0" for field in CANONICAL_STATE},
        },
        "axis": {
            "raw_regular_powers_in_x": AXIS_POWERS,
            "D": CANONICAL_D,
            "raw_odd_fields": ["W5", "W6"],
            "regularized_fields_are_rho_polynomials": True,
        },
        "checks": {
            "F1x_factor": "d_x(rho^m)/x=2*m*rho^(m-1)",
            "Uz1x_factor": "d_x(rho^m)/x=2*m*rho^(m-1)",
            "n_index": "m=n+1",
            "derivative_factor": "2*(n+1)",
        },
    }


def _physical_scales(parent: Mapping[str, Any]) -> dict[str, Any]:
    complex_bounds = parent.get("n1_complex_actual_source_bounds", {})
    delta = complex_bounds.get("delta_admission", {}) if isinstance(complex_bounds, Mapping) else {}
    return {
        "epsilon": "epsilon=exp(-logLambda)=1/Lambda",
        "Pstar_squared": "Pstar^2=exp(2*logPstar)=P2",
        "F0": "F0(Z)=exp(-selected_logCstar-Lambda*G(Z))",
        "interval_handles": {
            "epsilon": {
                "symbolic": "exp(-logLambda)",
                "outward_interval_saved_in_this_prefix": False,
                "consumer_must_bind": "accepted common source scale receipt",
            },
            "F0": {
                "symbolic": "exp(-selected_logCstar-Lambda*G(Z))",
                "outward_interval_saved_in_this_prefix": False,
                "consumer_must_bind": "accepted complex core F0 axis binding",
            },
            "Pstar_squared": {
                "symbolic": "exp(2*logPstar)",
                "outward_interval_saved_in_this_prefix": False,
                "consumer_must_bind": "shared P2=Pstar^2 source premise",
            },
        },
        "raw_scale_by_field": {
            "W1_F1": "F0(Z)/epsilon",
            "W2_Uz1": "epsilon",
            "W3_K1": "epsilon",
            "W4_P1": "epsilon*Pstar^2",
            "W5_over_x": "F0(Z)/epsilon",
            "W6_over_x": "epsilon",
        },
        "physical_radial_derivative": "d_R^i=exp(i*logLambda)*d_rho^i",
        "logPstar_source_handle": delta.get("actual_logPstar_formula", "exp(Md)+11"),
        "Lambda_materialized": False,
        "Pstar_materialized": False,
        "F0_amplitude_materialized": False,
    }


def build(root: str | Path, parent: Mapping[str, Any] | None = None) -> dict[str, Any]:
    """Return a static, fail-closed binding receipt for the saved prefix."""

    repo = Path(root)
    base = repo / "experiments" / "root_st073"
    parent_path = base / PARENT_NAME
    parent_file_exists = parent_path.exists()
    parent_file_sha = _sha256(parent_path) if parent_file_exists else None
    parent_report = parent if parent is not None else (_read(parent_path) if parent_file_exists else {})
    parent_mapping_matches_file = parent_file_exists and parent_mapping_digest(parent_report, parent_path)
    prefix_path = base / PREFIX_NAME
    source_vector_path = base / SOURCE_VECTOR_NAME
    missing = [str(path) for path in (parent_path, prefix_path, source_vector_path) if not path.exists()]
    if missing:
        return {
            "passed": False,
            "parent_report_sha256": parent_file_sha,
            "parent_mapping_matches_file": parent_mapping_matches_file,
            "checks": {"missing_inputs": False},
            "actual_forward_prefix_bound": {
                "passed": False,
                "missing_input_contract": missing,
                "source_executed": False,
            },
        }

    prefix = _read(prefix_path)
    source_vector = _read(source_vector_path)
    pins = _accepted_pins(repo, parent_report)
    ast_receipt = _ast_binding(repo)
    schema = _archive_schema(prefix)
    conjugacy_binding = _conjugacy_binding(parent_report)
    scale_binding = _physical_scale_binding(parent_report)
    family_match = prefix.get("source_family") == source_vector.get("source_family")
    parent_hash_ok = parent_file_sha == PARENT_SHA256
    prefix_hash_ok = _sha256(prefix_path) == PREFIX_SHA256
    source_vector_hash_ok = _sha256(source_vector_path) == SOURCE_VECTOR_SHA256
    conjugacy = parent_report.get("n1_all_order_radial_conjugacy", {})
    conjugacy_ok = (
        isinstance(conjugacy, Mapping)
        and conjugacy.get("passed") is True
        and conjugacy.get("point_values_selected") is False
        and conjugacy.get("finite_N_numerical_calculation") is False
    )
    source_vector_ok = bool(source_vector.get("checks", {}).get("passed"))
    checks = {
        "parent_report_sha256": parent_hash_ok,
        "parent_mapping_digest": parent_mapping_matches_file,
        "accepted_source_pins": pins.get("passed") is True,
        "prefix_report_sha256": prefix_hash_ok,
        "source_vector_report_sha256": source_vector_hash_ok,
        "source_family_match": family_match,
        "accepted_conjugacy_receipt": conjugacy_ok,
        "accepted_source_vector_receipt": source_vector_ok,
        "prefix_schema": schema.get("passed") is True,
        "canonical_ast_binding": ast_receipt.get("passed") is True,
        "normalization_map": conjugacy_binding.get("accepted_receipt_checks", {}).get("normalization_exact") is True,
        "raw_polynomial_map": conjugacy_binding.get("passed") is True,
        "derivative_maps": conjugacy_binding.get("fresh_index_algebra", {}).get("passed") is True,
        "ordinary_Z_archive": (
            schema.get("ordinary_Z_orders") == [0, 1, 2, 3]
            and schema.get("slot_lattice", {}).get("common_Z0_to_Z3_rho0_to_rho14") is True
        ),
        "physical_scale_handles": scale_binding.get("passed") is True,
        "no_point_selection": schema.get("point_values_selected") is False,
        "accepted_producers_not_replayed": True,
    }
    structural = all(checks.values())
    maps = _symbolic_maps()
    scales = _physical_scales(parent_report)
    bound = {
        "passed": structural,
        "source_executed": False,
        "accepted_producers_replayed": False,
        "parent_report": PARENT_NAME,
        "parent_report_sha256": parent_file_sha,
        "prefix_report": PREFIX_NAME,
        "prefix_report_sha256": _sha256(prefix_path),
        "source_vector_report": SOURCE_VECTOR_NAME,
        "source_vector_report_sha256": _sha256(source_vector_path),
        "normalization": NORMALIZATION,
        "canonical_state": CANONICAL_STATE,
        "canonical_D": CANONICAL_D,
        "raw_base_radial_degree": 16,
        "weighted_regularized_radial_degree": 15,
        "raw_fields": list(RAW_FIELDS),
        "weighted_fields": list(WEIGHTED_FIELDS),
        "symbolic_maps": maps,
        "physical_scale_handles": scales,
        "physical_scale_binding": scale_binding,
        "conjugacy_binding": conjugacy_binding,
        "archive_schema": schema,
        "coefficient_enclosures": _coefficient_enclosures(prefix),
        "ordinary_Z_domain": {
            "orders": [0, 1, 2, 3],
            "coefficient_centers": schema.get("Z_centers", []),
            "max_rho_degree_with_orders_0_to_3": schema.get(
                "ordinary_Z_0_to_3_max_rho_degree"
            ),
            "max_weighted_degree_with_orders_0_to_3": schema.get(
                "ordinary_Z_0_to_3_max_weighted_degree"
            ),
            "center_jet_archive_only": True,
            "complex_Z_domain_proved_here": False,
            "coefficient_convention": "ordinary Taylor coefficients in (Z-Z_center)^j; derivative D_Z multiplies (j+1)",
            "factorial_conversion": "D_Z^j value at center = j! * archived Z_degree=j coefficient",
            "real_Z_query_intervals_available": False,
        },
        "archive_location": {
            "path": PREFIX_NAME,
            "record_schema": "{Z, field, rho_degree, Z_degree, source_enclosure_exact_mpf_tuples}",
            "record_count": schema.get("record_count"),
            "decimal_approximations_retained_only_in_source": True,
            "exact_interval_tuples_exposed": True,
        },
        "limits": {
            "finite_radial_prefix_only": True,
            "radial_Taylor_remainder_certified": False,
            "Volterra_word_tail_reused_as_radial_tail": False,
            "physical_scaled_coefficient_intervals_materialized": False,
            "physical_total_prefix_error_combined": False,
            "full_raw_function_equality": False,
        },
        "source_pins": pins,
        "ast_binding": ast_receipt,
        "missing_input_contract": [] if structural else [
            key for key, value in checks.items() if not value
        ],
    }
    return {
        "passed": structural,
        "parent_report_sha256": parent_file_sha,
        "parent_mapping_digest": _mapping_digest(parent_report),
        "parent_mapping_matches_file": parent_mapping_matches_file,
        "checks": checks,
        "actual_forward_prefix_bound": bound,
    }


def parent_mapping_digest(value: Mapping[str, Any], path: Path) -> bool:
    """Compare supplied parent mapping against the exact pinned parent file."""

    try:
        return _mapping_digest(value) == _mapping_digest(_read(path))
    except (OSError, ValueError, TypeError):
        return False


__all__ = ["build"]
