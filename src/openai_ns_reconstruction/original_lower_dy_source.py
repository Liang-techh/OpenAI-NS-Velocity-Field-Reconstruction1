"""Recover the saved second-micro native-2 lower ``d_y`` source.

This adapter reads serialized source rows and graph bindings only.  It does
not construct an original owner, replay a producer, or call an ancestor
calculation.  The lower row is the ordinary-Z jet

``actual.Phi / comparison.Phi[0] * comparison.Psi[1] * rho``.

The interval rows are kept separate from their point approximations.  In
particular, the positive Cauchy continuation estimates widen comparison
rows; they are never used as signed values.  The physical epsilon and hb
factors stay in the source metadata and are not folded into this normalized
jet.
"""

from __future__ import annotations

import ast
import gzip
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from mpmath.libmp import mpf_abs, mpf_cmp

from .original_leading_viscous_kernel import decode_pair, encode_pair
from .original_source_jet_kernel import transported_lower_dy


_FAMILY = {
    "actual_five_defect_family_sha256": "3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894",
    "implicit_source_sha256": "5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae",
    "datum_enclosure_sha256": "dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010",
}

_HIGH = (
    "lei_ren_part1_paper_compliant_current_original_prescribed_high_radial_source_functions.json.gz",
    "5b54a210d48077a984b72b3c74a182c628bf00fb4c2b71f8f164845c228de991",
)
_PARENT = (
    "lei_ren_part1_paper_compliant_current_original_leading_viscous_source_numeric.json",
    "d1aa3644ea8158fe78dccb210a0a75a4c81ad257b8504758cb5797845b8003d4",
)
_INPUTS = {
    "second_micro_numeric": (
        "lei_ren_part1_paper_compliant_current_original_second_micro_numeric.json",
        "ffe82eda3eb0292eb628f2325bfd7eebeb0045911375c622da8609ec5c082156",
    ),
    "first_collar_jets": (
        "lei_ren_part1_paper_compliant_current_original_first_collar_jets.json",
        "a1a2818d18bf734c59460701d9ec4b342d6a8392ab97f581915766e9e6701a2f",
    ),
    "bridge_integrals": (
        "lei_ren_part1_paper_compliant_actual_bridge_integrals.json",
        "7e25fa846ef77e408dea54acf331df081e96b45a702ff3c9bc87136221144455",
    ),
    "analytic_domain": (
        "lei_ren_part1_paper_compliant_current_original_n1_inner_analytic_domain.json.gz",
        "c9905f7a9d13d4f651d1f55036932ea6c69db4fc5fa890e67a9839e611af4db7",
    ),
}

_SOURCE_FILES = {
    "second_micro_m15_source": (
        "src/openai_ns_reconstruction/original_second_micro_m15.py",
        "458be364bcc329e103e9be5d98fa14f9c1cfab39b74032d1642d864da6a83a6b",
    ),
    "second_micro_numeric_source": (
        "src/openai_ns_reconstruction/original_second_micro_numeric.py",
        "1c7f0d61ec2c9bcbc84454c2e047dd6c4382c3fa789be46f2339ffbe3180af2a",
    ),
    "first_collar_jets_source": (
        "src/openai_ns_reconstruction/original_first_collar_jets.py",
        "71c424d11b6bd06f26e86c01fca2640ecdf13fcacec69253b297dc7a31972e5c",
    ),
    "first_collar_transport_source": (
        "src/openai_ns_reconstruction/original_first_collar_transport.py",
        "10e4552fa09c9392e4e42a1fc1a76930905121f8c6bc5ac9d96e4d9418bd4305",
    ),
    "leading_kernel_source": (
        "src/openai_ns_reconstruction/original_leading_viscous_kernel.py",
        "fbd4b151f1e681eff6be9e4c8e30820416ee75fda4ac8d546a7cb335bd396ef6",
    ),
    "macro_m15_source": (
        "src/openai_ns_reconstruction/original_macro_m15.py",
        "5d2c162e6c61bee3d05955004114c9be73b44471555c0638a06b0bfb3d1febf5",
    ),
    "source_jet_kernel_source": (
        "src/openai_ns_reconstruction/original_source_jet_kernel.py",
        None,
    ),
}


def _mapping(value: Any) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def _repo_root(root: Any, directory: Any) -> Path:
    for candidate in (root, directory):
        if candidate is None:
            continue
        try:
            path = Path(candidate)
        except TypeError:
            continue
        if path.is_file():
            path = path.parent
        if (path / "src").is_dir() and (path / "experiments").is_dir():
            return path
    return Path(__file__).resolve().parents[2]


def _relative(path: Path, repo: Path) -> str:
    try:
        return path.relative_to(repo).as_posix()
    except ValueError:
        return path.as_posix()


def _sha256(path: Path) -> str | None:
    try:
        return hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError:
        return None


def _load(value: Any, repo: Path, *, name: str, default: Path | None = None) -> tuple[Mapping[str, Any], dict[str, Any]]:
    if isinstance(value, Mapping):
        return value, {"name": name, "mode": "supplied_mapping", "passed": True}
    path = Path(value) if value is not None else default
    if path is None:
        return {}, {"name": name, "passed": False, "reason": "serialized input was not supplied"}
    if not path.is_absolute():
        path = repo / path
    try:
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rt", encoding="utf-8") as stream:
            data = json.load(stream)
    except (OSError, ValueError, EOFError, json.JSONDecodeError) as exc:
        return {}, {"name": name, "path": _relative(path, repo), "passed": False, "reason": "serialized input is unreadable", "error": str(exc)}
    return (data if isinstance(data, Mapping) else {}), {
        "name": name,
        "path": _relative(path, repo),
        "sha256": _sha256(path),
        "passed": isinstance(data, Mapping),
    }


def _pin(path: Path, repo: Path, *, name: str, expected: str | None, parse_ast: bool = False) -> dict[str, Any]:
    actual = _sha256(path)
    parsed = None
    error = None
    if parse_ast and actual is not None:
        try:
            ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
            parsed = True
        except (OSError, SyntaxError, UnicodeError) as exc:
            parsed = False
            error = str(exc)
    passed = actual is not None and (expected is None or actual == expected)
    if parsed is False:
        passed = False
    row: dict[str, Any] = {
        "name": name,
        "path": _relative(path, repo),
        "sha256": actual,
        "expected_sha256": expected,
        "passed": bool(passed),
    }
    if parsed is not None:
        row["ast_parsed_without_execution"] = parsed
    if error:
        row["error"] = error
    return row


def _ast_expression_guard(path: Path, repo: Path, *, name: str, function: str, fragments: list[str]) -> dict[str, Any]:
    """Parse source and require the exact source expressions used here."""
    try:
        source = path.read_text(encoding="utf-8")
        tree = ast.parse(source, filename=str(path))
    except (OSError, SyntaxError, UnicodeError) as exc:
        return {"name": name, "path": _relative(path, repo), "function": function, "passed": False, "reason": str(exc)}
    compact = "".join(source.split())
    function_present = any(
        isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == function
        for node in ast.walk(tree)
    )
    checks = {fragment: fragment.replace(" ", "") in compact for fragment in fragments}
    return {
        "name": name,
        "path": _relative(path, repo),
        "function": function,
        "function_present": function_present,
        "required_fragments": checks,
        "passed": bool(function_present and all(checks.values())),
    }


def _ast_expression_guards(repo: Path) -> dict[str, dict[str, Any]]:
    return {
        "second_micro_inputs": _ast_expression_guard(
            repo / "src/openai_ns_reconstruction/original_second_micro_m15.py",
            repo,
            name="second_micro_inputs_formula",
            function="inputs",
            fragments=[
                "actual={name:a[name]+o._norm_error_jet(o._width_error(v),order) for name,v in delta.items()}",
                "comparison={name:[inlet['truth']['radial'][name][i]+o._norm_error_jet(",
                "rho=J.constant(iv,iv.mpf([4,_bounds(self.audit,iv.mpf(4)+",
            ],
        ),
        "second_micro_numeric_envelope": _ast_expression_guard(
            repo / "src/openai_ns_reconstruction/original_second_micro_numeric.py",
            repo,
            name="second_micro_numeric_finite_N_envelope",
            function="_numeric",
            fragments=[
                "delta_phi=bt['Phi']*o.Jet(iv,[o._symmetric(v) for v in exp_error])",
                "point[name]=bp[name];truth[name]=bt[name]+change",
            ],
        ),
        "macro_lower_dy": _ast_expression_guard(
            repo / "src/openai_ns_reconstruction/original_macro_m15.py",
            repo,
            name="macro_lower_dy_formula",
            function="leading",
            fragments=[
                "lower_dy = (lower['actual']['Phi']/lower['comparison']['Phi'][0]",
                "lower_dy=J(iv,[o._transfer(v) for v in lower_dy.coefficients])",
            ],
        ),
        "source_jet_transport": _ast_expression_guard(
            repo / "src/openai_ns_reconstruction/original_source_jet_kernel.py",
            repo,
            name="ordinary_transport_kernel",
            function="transported_lower_dy",
            fragments=[
                "ratio = ordinary_quotient(actual_phi,comparison_phi)",
                "return ordinary_product(ordinary_product(ratio,comparison_psi_radial),rho)",
            ],
        ),
    }


def _path_for(directory: Any, repo: Path, basename: str) -> Path:
    candidates: list[Path] = []
    if directory is not None:
        try:
            base = Path(directory)
            candidates.append(base if base.is_file() else base / basename)
        except TypeError:
            pass
    candidates.extend((repo / "experiments" / "root_st073" / basename, repo / basename))
    for path in candidates:
        if path.is_file():
            return path
    return candidates[0]


def _meta_path(meta: Mapping[str, Any], repo: Path, fallback: Path) -> Path:
    path_text = meta.get("path")
    if not isinstance(path_text, str):
        return fallback
    path = Path(path_text)
    return path if path.is_absolute() else repo / path


def _exact(ctx: Any, iv: Any, row: Mapping[str, Any]) -> Any:
    lower = row.get("lower_exact_mpf_tuple")
    upper = row.get("upper_exact_mpf_tuple")
    if isinstance(lower, (list, tuple)) and isinstance(upper, (list, tuple)):
        return iv.mpf([ctx.make_mpf(tuple(lower)), ctx.make_mpf(tuple(upper))])
    endpoints = row.get("source_enclosure_exact_mpf_tuples")
    if isinstance(endpoints, (list, tuple)) and len(endpoints) == 2:
        return iv.mpf([ctx.make_mpf(tuple(endpoints[0])), ctx.make_mpf(tuple(endpoints[1]))])
    raise ValueError("exact directed interval endpoints are required")


def _upper_abs(ctx: Any, value: Any) -> Any:
    if hasattr(value, "_mpi_"):
        lower, upper = (mpf_abs(tuple(v)) for v in value._mpi_)
    else:
        lower = upper = mpf_abs(value._mpf_)
    return ctx.make_mpf(lower if mpf_cmp(lower, upper) >= 0 else upper)


def _symmetric(iv: Any, ctx: Any, value: Any) -> Any:
    bound = _upper_abs(ctx, value)
    # Negation and endpoint magnitudes must not round in the lower-precision
    # approximation context before constructing the directed interval.
    return iv.mpf([-1, 1]) * iv.mpf(bound)


def _compact(row: Mapping[str, Any], *, kind: str, signed: bool = True) -> dict[str, Any]:
    endpoints = row.get("source_enclosure_exact_mpf_tuples")
    valid = isinstance(endpoints, (list, tuple)) and len(endpoints) == 2 and all(
        isinstance(item, (list, tuple)) and len(item) == 4 for item in endpoints
    )
    return {
        "approximation": row.get("approximation"),
        "absolute_error_upper": row.get("absolute_error_upper"),
        "raw_directed_error_exact_mpf_tuple": row.get("raw_directed_error_exact_mpf_tuple"),
        "source_enclosure_exact_mpf_tuples": [list(item) for item in endpoints] if valid else endpoints,
        "error_applies_to_serialized_decimal": row.get("error_applies_to_serialized_decimal") is True,
        "source_kind": kind,
        "signed_source_row": signed,
        "approximation_is_not_exact_source_value": True,
        "source_enclosure_is_not_a_selected_value": True,
        "exact_directed_endpoints_present": valid,
    }


def _ordinary_pairs(ctx: Any, iv: Any, rows: list[Mapping[str, Any]]) -> list[tuple[Any, Any]]:
    """Decode saved Taylor coefficients and convert them to ordinary rows."""
    pairs: list[tuple[Any, Any]] = []
    for index, row in enumerate(rows):
        point, truth = decode_pair(ctx, iv, row)
        factor = math.factorial(index)
        pairs.append((point * factor, truth * factor))
    return pairs


def _samples(report: Mapping[str, Any]) -> list[Mapping[str, Any]]:
    return [sample for sample in report.get("samples", ()) if isinstance(sample, Mapping)]


def _family(data: Mapping[str, Any]) -> Mapping[str, Any]:
    nested = data.get("source_family")
    source = nested if isinstance(nested, Mapping) else data
    return {key: source.get(key) for key in _FAMILY}


def _z_targets(parent: Mapping[str, Any]) -> list[str]:
    packet = _mapping(parent.get("inherited_m15_packet"))
    values: list[str] = []
    for sample in _samples(packet):
        z = str(sample.get("Z"))
        if z not in values:
            values.append(z)
    if not values:
        values = [str(sample.get("Z")) for sample in _samples(parent) if sample.get("Z") is not None]
    return values or ["-.3", "0", ".3"]


def _row_sample(report: Mapping[str, Any], z: str, coordinate: str | None = None, phase: str | None = None) -> Mapping[str, Any] | None:
    for sample in _samples(report):
        if str(sample.get("Z")) != z:
            continue
        if coordinate is not None:
            item = _mapping(sample.get("coordinate"))
            if str(item.get("value", sample.get("coordinate"))) != coordinate:
                continue
        if phase is not None and str(sample.get("phase")) != phase:
            continue
        return sample
    return None


def _coefficient_convention_guard(numeric: Mapping[str, Any], first: Mapping[str, Any], targets: list[str]) -> dict[str, Any]:
    numeric_rows = []
    first_rows = []
    for z in targets:
        numeric_sample = _row_sample(numeric, z, coordinate="2")
        first_sample = _row_sample(first, z, phase="0")
        numeric_rows.append(
            isinstance(numeric_sample, Mapping)
            and "k=D_Z^k/k!" in str(numeric_sample.get("coefficient_convention", ""))
        )
        first_rows.append(
            isinstance(first_sample, Mapping)
            and "native fixed-phase Z derivative / k!" in str(first_sample.get("units", ""))
        )
    return {
        "actual_numeric_Taylor_convention": "k=D_Z^k/k!",
        "comparison_core_Taylor_convention": "coefficient k = native fixed-phase Z derivative / k!",
        "actual_rows_passed": all(numeric_rows),
        "comparison_rows_passed": all(first_rows),
        "passed": bool(numeric_rows and first_rows and all(numeric_rows) and all(first_rows)),
    }


def _node(nodes: Any, number: int) -> Any:
    if isinstance(nodes, Mapping):
        return nodes.get(number, nodes.get(str(number)))
    try:
        return nodes[number]
    except (IndexError, KeyError, TypeError):
        return None


def _graph_binding(nodes: Any, high: Mapping[str, Any]) -> dict[str, Any]:
    binding = _mapping(high.get("prescribed_high_radial_source_bindings"))
    chart = _mapping(_mapping(binding.get("charts")).get("second_micro"))
    background = _mapping(_mapping(binding.get("backgrounds")).get("second_micro"))
    fields = _mapping(background.get("fields"))
    radial = _mapping(background.get("radial"))
    handles = {
        "V": list(fields.get("V", ()))[:3],
        "m": list(fields.get("m", ()))[:3],
        "V_y": list(radial.get("V_y", ()))[:3],
    }
    expected = {
        "coordinate": 66988,
        "radius": 68644,
        "jacobian": 8,
        "phase": 97301,
        "domain": [1, 48],
    }
    ports = {key: chart.get(key) for key in expected}
    port_passed = all(ports.get(key) == value for key, value in expected.items())
    handle_passed = handles == {
        "V": [69298, 69299, 69300],
        "m": [69305, 69306, 69307],
        "V_y": [69540, 69541, 69542],
    }
    source_nodes_present = all(_node(nodes, value) is not None for values in handles.values() for value in values)
    inverse_pstar_node = 69245
    query = _node(nodes, 48)
    query_passed = query == {"operation": "exact_rational", "numerator": 2, "denominator": 1}
    # A caller may supply only the selected nodes.  The archive itself is the
    # fallback witness for the normalization node; it is still recorded.
    normalization_node = _node(nodes, inverse_pstar_node)
    if normalization_node is None:
        normalization_node = _node(_mapping(high).get("exact_graph_nodes"), inverse_pstar_node)
    pstar = _mapping(binding.get("downstream")).get("Pstar")
    return {
        "chart": "second_micro",
        "ports": ports,
        "expected_ports": expected,
        "ports_passed": port_passed,
        "source_handles": handles,
        "source_handles_expected": {
            "V": [69298, 69299, 69300],
            "m": [69305, 69306, 69307],
            "V_y": [69540, 69541, 69542],
        },
        "source_handles_passed": handle_passed,
        "source_nodes_present": source_nodes_present,
        "normalization_invPstar_node": inverse_pstar_node,
        "normalization_invPstar_present": normalization_node is not None,
        "Pstar_node": pstar,
        "query_node": 48,
        "query_value": "2",
        "query_passed": query_passed,
        "macro_J_node": 66987,
        "macro_J_applied_during_copy": False,
        "passed": bool(port_passed and handle_passed and source_nodes_present and query_passed),
    }


def _width_cap_summary(ctx: Any, cap: Any) -> dict[str, Any]:
    value = ctx.make_mpf(tuple(cap._mpi_[1])) if hasattr(cap, "_mpi_") else cap
    return {
        "upper_exact_mpf_tuple": list(value._mpf_),
        "positive_remainder_only": True,
        "used_to_widen_source_enclosure": True,
        "source_value_selected": False,
    }


def _comparison_rows(
    *,
    ctx: Any,
    iv: Any,
    first: Mapping[str, Any],
    bridge: Mapping[str, Any],
    domain: Mapping[str, Any],
    name: str,
    radial_index: int,
    order: int = 5,
) -> tuple[list[Any], list[Any], list[dict[str, Any]]]:
    raw = _mapping(first.get("core_inlet_radial_derivatives")).get(name)
    if not isinstance(raw, list) or len(raw) <= radial_index or len(raw[radial_index]) < 3:
        raise ValueError(f"saved core comparison {name}[{radial_index}] Z0/Z1/Z2 rows are missing")
    # ``core_inlet_radial_derivatives`` stores the axial Taylor jet.  The
    # radial index selects the fixed radial derivative in the Cauchy bound;
    # it is not the ordinary-Z index below.
    base = [decode_pair(ctx, iv, row) for row in raw[radial_index][:3]]
    hcap = _exact(ctx, iv, _mapping(bridge.get("width_enclosure_cap")))
    weight = _exact(ctx, iv, _mapping(bridge.get("axial_jet_weight")))
    d = _mapping(domain.get("joint_core_domain"))
    majorants = _mapping(d.get("Xh_norm_majorants"))
    bound = _exact(ctx, iv, _mapping(majorants.get(name))) * iv.mpf(d.get("core_value_factor", "4"))
    rr = iv.mpf(d.get("Cauchy_rho_radius", ".5"))
    rz = _exact(ctx, iv, _mapping(d.get("Cauchy_Z_radius")))
    change = 4 * iv.mpf(2) * iv.exp(hcap * 2) * bound * math.factorial(radial_index + 1) / rr ** (radial_index + 1)
    change = sum((change * weight ** k / rz ** k for k in range(order + 1)), iv.mpf(0))
    cap = hcap * change
    widened: list[tuple[Any, Any]] = []
    caps: list[dict[str, Any]] = []
    for ordinary_index, (point, truth) in enumerate(base):
        # _norm_error_jet stores cap/weight**k as a Taylor coefficient.  The
        # ordinary derivative row therefore carries k! as well.
        taylor_cap = cap / weight ** ordinary_index
        ordinary_factor = math.factorial(ordinary_index)
        ordinary_cap = taylor_cap * ordinary_factor
        widened.append((
            point * ordinary_factor,
            truth * ordinary_factor + _symmetric(iv, ctx, ordinary_cap),
        ))
        caps.append({
            "radial_index": radial_index,
            "ordinary_Z_index": ordinary_index,
            "order": order,
            "formula": "factorial(k)*hb_cap*sum_j[0..order](4*s*exp(hb_cap*s)*Xh_norm_majorant*factorial(radial_index+1)/rho_radius^(radial_index+1)*(weight/Z_radius)^j)/weight^k",
            "cap": _width_cap_summary(ctx, ordinary_cap),
            "positive_remainder_only": True,
        })
    return [pair[0] for pair in widened], [pair[1] for pair in widened], caps


def _actual_rows(ctx: Any, iv: Any, numeric: Mapping[str, Any], z: str) -> tuple[list[Any], list[Any]]:
    sample = _row_sample(numeric, z, coordinate="2")
    if sample is None:
        raise ValueError(f"saved second-micro coordinate 2 sample is missing at Z={z}")
    leading = _mapping(sample.get("normalized_modified_leading"))
    raw = leading.get("Phi")
    if not isinstance(raw, list) or len(raw) < 3:
        raise ValueError(f"saved actual Phi Z0/Z1/Z2 rows are missing at Z={z}")
    pairs = _ordinary_pairs(ctx, iv, raw[:3])
    return [pair[0] for pair in pairs], [pair[1] for pair in pairs]


def _rho_rows(ctx: Any, iv: Any, bridge: Mapping[str, Any]) -> tuple[list[Any], list[Any], dict[str, Any]]:
    hcap = _exact(ctx, iv, _mapping(bridge.get("width_enclosure_cap")))
    cap = hcap * 8 * iv.exp(hcap * 2)
    point = [ctx.mpf(4), ctx.mpf(0), ctx.mpf(0)]
    upper = _upper_abs(ctx, cap)
    # Keep the upper sum in the interval context.  ``4 + upper`` in the
    # point context would round the directed endpoint before reboxing it.
    truth = [iv.mpf([4, 4]) + iv.mpf([0, upper]), iv.mpf(0), iv.mpf(0)]
    return point, truth, {
        "formula": "rho=4*exp(hb*s), evaluated at s=2; width cap hb*(4*s*exp(hb*s))",
        "ordinary_Z_derivatives": ["0", "0"],
        "width_cap": _width_cap_summary(ctx, cap),
        "epsilon_and_hb_factored": True,
    }


def _one_row(
    *,
    ctx: Any,
    iv: Any,
    z: str,
    numeric: Mapping[str, Any],
    first: Mapping[str, Any],
    bridge: Mapping[str, Any],
    domain: Mapping[str, Any],
) -> dict[str, Any]:
    numeric_sample = _row_sample(numeric, z, coordinate="2")
    first_sample = _row_sample(first, z, phase="0")
    if numeric_sample is None:
        raise ValueError(f"saved second-micro coordinate 2 sample is missing at Z={z}")
    if first_sample is None:
        raise ValueError(f"saved first-collar phase 0 sample is missing at Z={z}")
    actual_point, actual_truth = _actual_rows(ctx, iv, numeric, z)
    comparison_phi_point, comparison_phi_truth, phi_caps = _comparison_rows(
        ctx=ctx, iv=iv, first=first_sample, bridge=bridge, domain=domain, name="Phi", radial_index=0
    )
    comparison_psi_point, comparison_psi_truth, psi_caps = _comparison_rows(
        ctx=ctx, iv=iv, first=first_sample, bridge=bridge, domain=domain, name="Psi", radial_index=1
    )
    rho_point, rho_truth, rho_meta = _rho_rows(ctx, iv, bridge)
    point = transported_lower_dy(actual_point, comparison_phi_point, comparison_psi_point, rho_point)
    truth = transported_lower_dy(actual_truth, comparison_phi_truth, comparison_psi_truth, rho_truth)
    return {
        "name": "second_micro_lower_dy",
        "Z": z,
        "coordinate": {"kind": "exact_decimal_phase", "value": "2"},
        "ordinary_Z_orders": [0, 1, 2],
        "formula": "actual.Phi / comparison.Phi[0] * comparison.Psi[1] * rho",
        "coefficient_convention": {
            "saved_rows": "Taylor coefficients k=D_Z^k/k!",
            "ordinary_kernel_rows": "factorial(k)*saved Taylor coefficient",
            "factorial_conversion_applied": True,
        },
        "actual": {
            "field": "normalized_modified_leading.Phi",
            "source_rows": [_compact(row, kind="saved_second_micro_actual_Phi") for row in _mapping(numeric_sample.get("normalized_modified_leading")).get("Phi", [])[:3]],
            "inherited_conservative_envelope": True,
            "finite_N_symmetric_widening_retained": True,
            "signed_correction_selected": False,
        },
        "comparison": {
            "Phi_radial_index": 0,
            "Psi_radial_index": 1,
            "Phi_source_rows": [_compact(row, kind="saved_first_collar_core_comparison") for row in _mapping(first_sample.get("core_inlet_radial_derivatives")).get("Phi", [])[0][:3]],
            "Psi_source_rows": [_compact(row, kind="saved_first_collar_core_comparison") for row in _mapping(first_sample.get("core_inlet_radial_derivatives")).get("Psi", [])[1][:3]],
            "Phi_continuation_caps": phi_caps,
            "Psi_continuation_caps": psi_caps,
            "positive_caps_only_widen_enclosure": True,
        },
        "rho": rho_meta,
        "result": [encode_pair(ctx, iv, p, t) for p, t in zip(point, truth)],
        "point_remainder_rows": [encode_pair(ctx, iv, p, t) for p, t in zip(point, truth)],
        "transport": {
            "exact_interval_endpoints_copied_between_contexts": True,
            "macro_J66987_applied": False,
            "native_query_node": 48,
            "native_query_value": "2",
        },
        "physical_factors": {
            "epsilon": "Lambda^-1 factored",
            "hb": "hb factored through dy=hb ds",
            "physical_R_once": True,
            "nonzero_history_retained": True,
            "source_values_selected": False,
            "approximation_is_not_exact_source_value": True,
        },
        "passed": True,
    }


def build_lower_dy_source(
    *,
    nodes: Any,
    high_archive: Any,
    parent_report: Any,
    directory: Any,
    root: Any,
) -> dict[str, Any]:
    """Load saved rows and recover ordinary native-2 ``lower_dy``."""
    repo = _repo_root(root, directory)
    high_default = repo / "experiments" / "root_st073" / _HIGH[0]
    high, high_meta = _load(high_archive, repo, name="high_archive", default=high_default)
    parent_default = repo / "experiments" / "root_st073" / _PARENT[0]
    parent, parent_meta = _load(parent_report, repo, name="parent_report", default=parent_default)
    inputs: dict[str, Mapping[str, Any]] = {}
    input_meta: dict[str, dict[str, Any]] = {}
    for name, (basename, _) in _INPUTS.items():
        data, meta = _load(None, repo, name=name, default=_path_for(directory, repo, basename))
        inputs[name], input_meta[name] = data, meta

    pins: dict[str, dict[str, Any]] = {
        "high_archive": _pin(
            _meta_path(high_meta, repo, high_default),
            repo,
            name="high_archive",
            expected=_HIGH[1],
        ),
        "parent_report": _pin(
            _meta_path(parent_meta, repo, parent_default),
            repo,
            name="parent_report",
            expected=_PARENT[1],
        ),
    }
    for name, (relative, expected) in _SOURCE_FILES.items():
        pins[name] = _pin(repo / relative, repo, name=name, expected=expected, parse_ast=True)
    for name, (_, expected) in _INPUTS.items():
        path = repo / "experiments" / "root_st073" / _INPUTS[name][0]
        pins[name] = _pin(path, repo, name=name, expected=expected)

    family_checks = {
        name: {
            "family": dict(_family(data)),
            "passed": all(_family(data).get(key) == value for key, value in _FAMILY.items()),
        }
        for name, data in {"high_archive": high, "parent_report": parent, **inputs}.items()
    }
    graph = _graph_binding(nodes if nodes is not None else high.get("exact_graph_nodes"), high)
    ast_guards = _ast_expression_guards(repo)
    convention_guard = _coefficient_convention_guard(
        inputs["second_micro_numeric"], inputs["first_collar_jets"], _z_targets(parent)
    )
    source_bindings = {
        "role": "second_micro_transported_lower_dy",
        "source_chart": "second_micro",
        "native_coordinate_node": 66988,
        "native_radius_node": 68644,
        "native_query_node": 48,
        "native_query_value": "2",
        "formula": "lower.actual.Phi / lower.comparison.Phi[0] * lower.comparison.Psi[1] * lower.rho",
        "transport_formula": "copy exact interval endpoints between point and interval contexts; no macro J=66987 factor is applied during that copy",
        "graph": graph,
        "ast_expression_guards": ast_guards,
        "coefficient_convention": convention_guard,
        "modified_Qy_not_used": True,
        "source_caps_not_signed_values": True,
        "epsilon": "Lambda^-1 factored",
        "hb": "hb factored",
        "nonzero_histories_retained": True,
    }
    rows: list[dict[str, Any]] = []
    blocked: list[dict[str, Any]] = []
    ctx = mp.mp.clone()
    ctx.dps = 110
    iv = MPIntervalContext()
    iv.dps = 180
    for z in _z_targets(parent):
        try:
            rows.append(_one_row(
                ctx=ctx,
                iv=iv,
                z=z,
                numeric=inputs["second_micro_numeric"],
                first=inputs["first_collar_jets"],
                bridge=inputs["bridge_integrals"],
                domain=inputs["analytic_domain"],
            ))
        except (KeyError, TypeError, ValueError, ArithmeticError) as exc:
            blocked.append({
                "name": "second_micro_lower_dy",
                "Z": z,
                "coordinate": {"kind": "exact_decimal_phase", "value": "2"},
                "missing": ["saved signed actual/comparison/rho row needed by formula"],
                "reason": str(exc),
                "source_values_selected": False,
            })
    source_pins_passed = all(row.get("passed") is True for row in pins.values())
    family_passed = all(row.get("passed") is True for row in family_checks.values())
    rows_passed = bool(rows) and not blocked and all(row.get("passed") is True for row in rows)
    source_bindings["passed"] = bool(
        graph.get("passed")
        and family_passed
        and convention_guard.get("passed") is True
        and all(row.get("passed") is True for row in ast_guards.values())
        and rows_passed
    )
    limits = {
        "lower_dy_numeric": bool(rows),
        "lower_dy_C0_C1": bool(rows),
        "lower_dy_C2": bool(rows) and all(len(row.get("result", ())) >= 3 for row in rows),
        "macro_J66987_applied": False,
        "source_values_selected": False,
        "numeric_caps_used": False,
        "ancestor_replayed": False,
        "epsilon_factored": True,
        "hb_factored": True,
        "complete": False,
    }
    return {
        "passed": bool(source_pins_passed and family_passed and source_bindings["passed"]),
        "partial": bool(rows),
        "source_family": dict(_FAMILY),
        "parent_report": _PARENT[0],
        "parent_report_sha256": _PARENT[1],
        "source_pins": pins,
        "family_checks": family_checks,
        "source_bindings": source_bindings,
        "rows": rows,
        "point_remainder_rows": rows,
        "blocked_rows": blocked,
        "limits": limits,
        "contract": {
            "source_only": True,
            "actual_and_comparison_rows_kept_separate": True,
            "point_approximations_not_exact_source_values": True,
            "positive_remainders_widen_only": True,
            "caps_not_promoted": True,
            "modified_Qy_not_used": True,
            "macro_J66987_not_applied": True,
            "ancestor_constructors_producers_checkers_not_run": True,
            "epsilon_and_hb_factored": True,
            "nonzero_history_retained": True,
        },
    }


__all__ = ["build_lower_dy_source"]
