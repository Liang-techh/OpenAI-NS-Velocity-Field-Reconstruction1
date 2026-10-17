"""New same-generator modal source point algebra.

This module deliberately reads the saved first-collar and comparison leaves.
It does not call a core producer, a bridge producer, or an old source
checker.  The pressure sector is evaluated through ordinary Z order two from
the saved normalized pressure Taylor jet.  The other sectors remain explicit
formal blocks until their same-generator modal transport is implemented.
"""

from __future__ import annotations

import ast
import gzip
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from .original_leading_viscous_kernel import decode_pair, encode_pair


_FIRST_COLLAR = "lei_ren_part1_paper_compliant_current_original_first_collar_jets.json"
_BRIDGE = "lei_ren_part1_paper_compliant_actual_bridge_integrals.json"
_CORE = "lei_ren_part1_paper_compliant_current_original_core_numeric_approximation.json"
_PARENT = "lei_ren_part1_paper_compliant_current_original_source_jet_mass_numeric.json"
_PARENT_SHA256 = "c3b2e2a974ebca402e2d2762ab2c508150fdd72f80aa2882152e944bdad84efe"
_Z_TEXT = ("-.3", "0", ".3")


def _sha256(path: Path) -> str | None:
    try:
        return hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError:
        return None


def _load(path: Path) -> Mapping[str, Any]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as stream:
        value = json.load(stream)
    if not isinstance(value, Mapping):
        raise ValueError(f"saved mapping required: {path}")
    return value


def _repo(root: Any, directory: Any) -> Path:
    for candidate in (root, directory):
        if candidate is None:
            continue
        path = Path(candidate)
        if path.is_file():
            path = path.parent
        if (path / "src").is_dir() and (path / "experiments").is_dir():
            return path
    return Path(__file__).resolve().parents[2]


def _point_context() -> tuple[Any, MPIntervalContext]:
    ctx = mp.mp.clone()
    ctx.dps = 180
    iv = MPIntervalContext()
    iv.dps = 220
    return ctx, iv


def _pin(path: Path, repo: Path, *, name: str, parse_ast: bool = False) -> dict[str, Any]:
    digest = _sha256(path)
    parsed = True
    error = None
    if parse_ast and digest is not None:
        try:
            ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        except (OSError, SyntaxError, UnicodeError) as exc:
            parsed = False
            error = str(exc)
    try:
        relative = path.relative_to(repo).as_posix()
    except ValueError:
        relative = path.as_posix()
    row: dict[str, Any] = {
        "name": name,
        "path": relative,
        "sha256": digest,
        "passed": bool(digest is not None and parsed),
    }
    if parse_ast:
        row["ast_parsed_without_execution"] = parsed
    if error:
        row["error"] = error
    return row


def _source_pins(repo: Path, directory: Path) -> dict[str, Any]:
    base = directory
    entries = (
        ("parent_report", base / _PARENT, False),
        ("first_collar_C5", base / _FIRST_COLLAR, False),
        ("actual_bridge_C5", base / _BRIDGE, False),
        ("core_inlet_saved", base / _CORE, False),
        (
            "macro_signed_integrals_AST",
            repo / "experiments/root_st073/lei_ren_part1_paper_compliant_macro_signed_integrals.py",
            True,
        ),
        (
            "macro_source_AST",
            repo / "src/openai_ns_reconstruction/original_macro_pressure_rows.py",
            True,
        ),
        (
            "inner_bridge_profiles_AST",
            repo / "experiments/root_st073/lei_ren_part1_paper_compliant_inner_bridge_profiles.py",
            True,
        ),
        (
            "actual_bridge_integrals_AST",
            repo / "experiments/root_st073/lei_ren_part1_paper_compliant_actual_bridge_integrals.py",
            True,
        ),
        (
            "original_core_numeric_AST",
            repo / "src/openai_ns_reconstruction/original_core_numeric.py",
            True,
        ),
        (
            "original_first_collar_jets_AST",
            repo / "src/openai_ns_reconstruction/original_first_collar_jets.py",
            True,
        ),
        ("pressure_datum_AST",base/"lei_ren_part1_paper_logarithmic_pressure_datum.py",True),
        ("compliant_pressure_AST",base/"lei_ren_part1_paper_compliant_pressure_source.py",True),
        ("comparison_continuation_AST",repo/"src/openai_ns_reconstruction/original_first_collar_transport.py",True),
        ("saved_comparison_export_AST",repo/"src/openai_ns_reconstruction/original_lower_dy_source.py",True),
        ("second_micro_pressure_AST",repo/"src/openai_ns_reconstruction/original_second_micro_m15.py",True),
        ("saved_scale_AST",repo/"src/openai_ns_reconstruction/original_first_collar_numeric.py",True),
    )
    return {name: _pin(path, repo, name=name, parse_ast=parse_ast) for name, path, parse_ast in entries}


def _sample(samples: list[Mapping[str, Any]], ztext: str, *, phase: str = "0") -> Mapping[str, Any]:
    target = mp.mpf(ztext)
    for sample in samples:
        if str(sample.get("phase")) != phase:
            continue
        try:
            if mp.mpf(str(sample.get("Z"))) == target:
                return sample
        except (TypeError, ValueError):
            continue
    raise KeyError(f"first_collar_jets sample Z={ztext}, phase={phase}")


def _factor_pair(ctx: Any, iv: Any, row: Mapping[str, Any], order: int) -> tuple[Any, Any]:
    point, truth = decode_pair(ctx, iv, dict(row))
    if order < 0:
        raise ValueError("Taylor order must be nonnegative")
    factor = ctx.mpf((1, 1, 2, 6)[order] if order < 4 else mp.factorial(order))
    return point * factor, truth * factor


def _cap_interval(ctx: Any, iv: Any, row: Mapping[str, Any]) -> Any:
    cap = row.get("cap", row)
    raw = cap.get("upper_exact_mpf_tuple")
    if not isinstance(raw, (list, tuple)) or len(raw) != 4:
        raise ValueError("missing exact continuation cap upper endpoint")
    upper = ctx.make_mpf(tuple(raw))
    return iv.mpf([-1,1])*iv.mpf(upper)


def _positive_width(ctx: Any, iv: Any, bridge: Mapping[str, Any]) -> Any:
    row = bridge.get("width_enclosure_cap")
    if not isinstance(row, Mapping):
        raise ValueError("missing actual_bridge_integrals.width_enclosure_cap")
    lower_raw = row.get("lower_exact_mpf_tuple")
    upper_raw = row.get("upper_exact_mpf_tuple")
    if not all(isinstance(raw, (list, tuple)) and len(raw) == 4 for raw in (lower_raw, upper_raw)):
        raise ValueError("missing exact positive width cap endpoints")
    return iv.mpf([ctx.make_mpf(tuple(lower_raw)), ctx.make_mpf(tuple(upper_raw))])


def _log_h_interval(ctx: Any, iv: Any, bridge: Mapping[str, Any]) -> Any:
    row = bridge.get("source_log_hb_enclosure")
    if not isinstance(row, Mapping):
        raise ValueError("missing actual_bridge_integrals.source_log_hb_enclosure")
    lower_raw = row.get("lower_exact_mpf_tuple")
    upper_raw = row.get("upper_exact_mpf_tuple")
    if not all(isinstance(raw, (list, tuple)) and len(raw) == 4 for raw in (lower_raw, upper_raw)):
        raise ValueError("missing exact source log-h endpoints")
    return iv.mpf([ctx.make_mpf(tuple(lower_raw)), ctx.make_mpf(tuple(upper_raw))])


def _function_ast_guard(path, function, fragments, *, class_name=None):
    tree=ast.parse(path.read_text(encoding="utf-8"))
    scope=next(node for node in tree.body if isinstance(node,ast.ClassDef) and node.name==class_name) if class_name else tree
    functions=[node for node in ast.walk(scope) if isinstance(node,ast.FunctionDef) and node.name==function]
    if len(functions)!=1:
        raise ValueError("Unique pinned source function required: "+str(path)+":"+function)
    node=functions[0];nodes={ast.dump(part,include_attributes=False) for part in ast.walk(node)}
    checks={}
    for fragment in fragments:
        wanted=ast.parse(fragment).body[0]
        if isinstance(wanted,ast.Expr):wanted=wanted.value
        checks[fragment]=ast.dump(wanted,include_attributes=False) in nodes
    return {"function":function,"class":class_name,"required_AST":checks,
            "function_AST_sha256":hashlib.sha256(ast.dump(node,include_attributes=False).encode()).hexdigest(),
            "source_executed":False,"passed":all(checks.values())}


def _source_guards(repo,parent,first,bridge):
    base=repo/"experiments/root_st073";src=repo/"src/openai_ns_reconstruction"
    specs=(
        ("pressure_mode_identity",base/"lei_ren_part1_paper_compliant_macro_signed_integrals.py","_mode_rows",None,
         ["[base['drive_pressure'], _zero(base['drive_pressure']), _zero(base['drive_pressure'])]"]),
        ("pressure_direction",base/"lei_ren_part1_paper_compliant_inner_bridge_profiles.py","direction",None,
         ["((z*pressure)*(2*(1+delta))-d*derivative(pressure))/(2*L)"]),
        ("bridge_p0",base/"lei_ren_part1_paper_compliant_inner_bridge_profiles.py","inputs","CompliantInnerBridgeProfiles",
         ["datum=self.core.datum.normalized_jets(endpoints(Z),6)",
          "pressure=IntervalTaylor(c,[c.mpf(endpoints(value)) for value in datum['normalized_pressure_coefficients']])"]),
        ("macro_quotient",base/"lei_ren_part1_paper_compliant_actual_bridge_integrals.py","prepare","CompliantActualBridgeIntegrals",
         ["data=self.comparison.inlet(Z,root)","end=self.comparison.micro(Z,2,root)",
          "phi2=value_enclosure(end['phi'])","data['phi0'].truncate(5)/phi2.truncate(5)"]),
        ("macro_source",src/"original_macro_pressure_rows.py","source_Uy","OriginalMacroPressureRows",
         ["-prepared['quotient']*module.weighted_modes(prepared['modes']['drive_'+name], modes)",
          "norm*o._transfer(6*b.cap*b.ctx.exp(6*b.cap))", "total+o._norm_error_jet(nonlinear,order)"]),
        ("saved_p0_normalization",src/"original_first_collar_jets.py","inlet_jets","OriginalFirstCollarJets",
         ["c.exp(2*(self.core.lp_iv if directed else self.core.logP)-(self.core.ll_iv if directed else self.core.logLambda))",
          "pressure=J(c, fixed['P0_Z_taylor'][:order+1])/pscale"]),
        ("constant_source_scales",src/"original_core_numeric.py","__init__","OriginalCoreNumeric",
         ["self.logP=c.exp(40)+11", "self.logLambda=4*self.logP+1000",
          "self.delta=c.exp(-4*self.logP-30)", "self.epsilon=c.exp(-self.logLambda)"]),
        ("saved_epsilon_log",src/"original_first_collar_numeric.py","evaluate","OriginalFirstCollarNumeric",
         ["self._summary(-self.core.logLambda,-self.core.ll_iv)"]),
        ("core_p0",src/"original_core_numeric.py","_seed","OriginalCoreNumeric",
         ["M2=_box(c,self.pressure['fixed_beta2_mass_normalized'])",
          "M0=_box(c,self.pressure['fixed_beta0_mass_normalized'])", "qs=self.q_jets(c,z0,n)",
          "-(M2*q+(M0 if k == 0 else 0)+remainder)", "[c.exp(2*lp-ll)*v for v in pressure]"]),
        ("datum_p0",base/"lei_ren_part1_paper_logarithmic_pressure_datum.py","normalized_jets","LogarithmicPressureDatum",
         ["q_jets(c,z,order)","-(self.m2*a+(self.m0 if n==0 else 0)+remainder)",
          "remainder=self.stages['z_flatten']['mass']", "self.flatten_complex_upper/self.rho**n"]),
        ("compliant_datum",base/"lei_ren_part1_paper_compliant_pressure_source.py","__init__","CompliantPressureDatum",
         ["self.rho=c.mpf('.25')", "sum((self.stages[n]['mass'] for n in BETA2),c.mpf(0))",
          "sum((self.stages[n]['mass'] for n in BETA0),c.mpf(0))"]),
        ("linear_continuation_cap",src/"original_lower_dy_source.py","_comparison_rows",None,
         ["cap=hcap*change", "taylor_cap=cap/weight**ordinary_index", "ordinary_cap=taylor_cap*ordinary_factor"]),
        ("second_micro_comparison",src/"original_second_micro_m15.py","inputs","OriginalSecondMicroPressureRows",
         ["self.t._core_comparison_changes(name,i,upper,order)","inlet['truth']['radial'][name][i]"]),
    )
    guards={name:_function_ast_guard(path,function,fragments,class_name=cls) for name,path,function,cls,fragments in specs}
    family=parent["source_family"]
    family_checks={name:all((data.get('source_family',data)).get(k)==v for k,v in family.items())
                   for name,data in (("first_collar",first),("bridge",bridge))}
    lower=parent["transported_lower_dy"]
    inherited=lower["passed"] is True and lower["source_bindings"]["passed"] is True
    convention=all("k!" in _sample(first["samples"],z)["units"] and
                   _sample(first["samples"],z)["axial_Taylor_order"]>=3 for z in _Z_TEXT)
    passed=all(row["passed"] for row in guards.values()) and all(family_checks.values()) and inherited and convention
    return {"passed":passed,"AST":guards,"family_checks":family_checks,
            "inherited_native2_comparison_source_passed":inherited,"saved_Taylor_Z3_available":convention,
            "same_analytic_core_source_family_and_datum":dict(family),
            "shared_phi_identity":"prepared phi2 = core Phi0 + native2 comparison continuation Delta",
            "linear_cap_identity":"saved ordinary cap=hcap*change(hcap)*k!/weight^k; actual Delta/h bounded by same monotone Cauchy coefficient",
            "pressure_mode_j1_j2_exact_zero":guards["pressure_mode_identity"]["passed"]}


def _q_jet(
    ctx: Any,
    iv: Any,
    parent: Mapping[str, Any],
    ztext: str,
    bridge: Mapping[str, Any],
    guards: Mapping[str,Any],
) -> tuple[list[Any], list[Any], Mapping[str, Any]]:
    """Return q point proxy and a saved comparison quotient enclosure.

    The selected production point ratio is deliberately the documented q=1
    proxy.  The truth ratio is formed from the saved core/comparison Phi
    leaves and the positive continuation caps; no midpoint is promoted.
    """
    rows = parent.get("transported_lower_dy", {}).get("rows", ())
    selected = None
    for row in rows:
        if str(row.get("Z")) == ztext or mp.mpf(str(row.get("Z"))) == mp.mpf(ztext):
            selected = row
            break
    if not isinstance(selected, Mapping):
        raise KeyError(f"parent transported comparison row Z={ztext}")
    comparison = selected.get("comparison")
    if not isinstance(comparison, Mapping):
        raise KeyError(f"parent comparison block Z={ztext}")
    phi_rows = comparison.get("Phi_source_rows")
    caps = comparison.get("Phi_continuation_caps")
    if not isinstance(phi_rows, (list, tuple)) or len(phi_rows) < 3:
        raise KeyError(f"parent comparison Phi_source_rows[0:3] Z={ztext}")
    if not isinstance(caps, (list, tuple)) or len(caps) < 3:
        raise KeyError(f"parent Phi_continuation_caps[0:3] Z={ztext}")
    if (comparison.get("Phi_radial_index")!=0 or selected["transport"].get("native_query_value")!="2"
            or selected["coefficient_convention"].get("factorial_conversion_applied") is not True):
        raise ValueError("Same native2 comparison Phi and ordinary-cap convention required")

    phi: list[Any] = []
    continuation: list[Any] = []
    for order in range(3):
        if caps[order].get("radial_index")!=0 or caps[order].get("ordinary_Z_index")!=order:
            raise ValueError("Saved comparison cap radial/axial indices changed")
        point, truth = _factor_pair(ctx, iv, phi_rows[order], order)
        # The source rows are saved first-collar Taylor coefficients.  The
        # continuation cap is already an ordinary derivative cap.
        phi.append(truth)
        continuation.append(_cap_interval(ctx, iv, caps[order]))
    from .original_signed_C1_source import correlated_phi_ratio

    hcap=_positive_width(ctx,iv,bridge);log_h=_log_h_interval(ctx,iv,bridge)
    from mpmath.libmp import mpf_cmp
    if mpf_cmp(log_h._mpi_[1],iv.log(hcap)._mpi_[0])>=0:
        raise ValueError("Original positive width must admit the saved continuation cap")
    qpoint, qtruth, proof = correlated_phi_ratio(
        ctx,
        iv,
        phi,
        continuation,
        hcap,
        log_h,
        shared_phi_identity_proved=guards["passed"],
        cap_linear_width_proved=guards["AST"]["linear_continuation_cap"]["passed"],
    )
    return qpoint, qtruth, proof


def _pressure_pair(
    ctx: Any,
    iv: Any,
    sample: Mapping[str, Any],
    qpoint: list[Any],
    qtruth: list[Any],
    bridge: Mapping[str,Any],
    guards: Mapping[str,Any],
) -> tuple[list[Any], list[Any],Mapping[str,Any]]:
    pressure_rows = sample.get("pressure", {}).get("P0_over_Pstar_squared")
    if not isinstance(pressure_rows, (list, tuple)) or len(pressure_rows) < 4:
        raise KeyError("first_collar pressure.P0_over_Pstar_squared[0:4]")
    ppoint: list[Any] = []
    ptruth: list[Any] = []
    for order in range(4):
        point, truth = _factor_pair(ctx, iv, pressure_rows[order], order)
        ppoint.append(point)
        ptruth.append(truth)

    scales = sample.get("source_scales")
    if not isinstance(scales, Mapping) or not isinstance(scales.get("log_epsilon"), Mapping):
        raise KeyError("first_collar source_scales.log_epsilon")
    logeps_point, logeps_truth = decode_pair(ctx, iv, dict(scales["log_epsilon"]))
    logp_point = (-logeps_point - ctx.mpf(1000)) / ctx.mpf(4)
    logp_truth = (-logeps_truth - iv.mpf(1000)) / iv.mpf(4)
    logscale_point = ctx.log(100) + 2 * logp_point - logeps_point
    logscale_truth = iv.log(100) + 2 * logp_truth - logeps_truth
    scale_point = ctx.exp(logscale_point)
    scale_truth = iv.exp(logscale_truth)
    delta_point = ctx.exp(logeps_point + ctx.mpf(970))
    delta_truth = iv.exp(logeps_truth + iv.mpf(970))

    from .original_signed_C1_source import macro_nonlinear_enclosure, pressure_source_ordinary

    source_point = pressure_source_ordinary(
        ctx.mpf(sample["Z"]), delta_point, ppoint, qpoint, scale_point
    )
    source_truth = pressure_source_ordinary(
        iv.mpf(str(sample["Z"])), delta_truth, ptruth, qtruth, scale_truth
    )

    raw_weight=bridge['axial_jet_weight']
    weight=iv.mpf([ctx.make_mpf(tuple(raw_weight[k+'_exact_mpf_tuple'])) for k in ('lower','upper')])
    source_truth,nonlinear=macro_nonlinear_enclosure(ctx,iv,source_truth,weight,
        _positive_width(ctx,iv,bridge),macro_norm_bound_proved=guards['AST']['macro_source']['passed'])
    return source_point, source_truth,nonlinear


def _blocked() -> list[dict[str, Any]]:
    return [
        {
            "sector": "hydro",
            "status": "not_implemented",
            "reason": "same-generator hydro modal transport and its comparison-moment remainder are not implemented",
            "missing_leaf": "factored finite-width comparison V/M/A jets; saved core C5 inlets alone are not the full transported source",
        },
        {
            "sector": "swirl",
            "status": "not_implemented",
            "reason": "dressed relative F0^2 jets and their formal absolute source scale are not yet propagated",
            "missing_leaf": "same-generator swirl modal jets with a combined-log physical scale proof",
            "formal_log_enclosure_retained": True,
        },
    ]


def build_modal_source_point(*, parent_report: Any, directory: Any, root: Any) -> dict[str, Any]:
    """Build three signed pressure-sector B/B_Z/B_ZZ source rows.

    A failed leaf lookup returns an explicit blocked result with no padded
    source rows.  All numerical values are new arithmetic over saved leaves.
    """
    repo = _repo(root, directory)
    base = Path(directory) if directory is not None else repo / "experiments/root_st073"
    if base.is_file():
        base = base.parent
    pins = _source_pins(repo, base)
    if not all(row.get("passed") is True for row in pins.values()):
        missing = [row.get("name") for row in pins.values() if row.get("passed") is not True]
        return {
            "passed": False,
            "source_family": None,
            "rows": [],
            "source_scope": "blocked_saved_leaf",
            "full_source_scope": False,
            "blocked": [{"reason": "unreadable or unparsable saved input", "missing_leaf": name} for name in missing],
            "source_pins": pins,
        }

    parent = _load(Path(parent_report)) if isinstance(parent_report, (str, Path)) else dict(parent_report)
    if _sha256(base / _PARENT) != _PARENT_SHA256:
        raise ValueError("Pinned accepted parent source-jet report required")
    if parent.get("checks", {}).get("passed") is not True:
        raise ValueError("Passed accepted parent report required")
    for pin in pins.values():
        relative=pin['path']
        expected=(_PARENT_SHA256 if relative=='experiments/root_st073/'+_PARENT else
                  parent['input_hashes'].get(Path(relative).name) if relative.startswith('experiments/root_st073/') else
                  parent['implementation_hashes'].get(relative))
        if not expected or pin['sha256']!=expected:
            raise ValueError('Source pin is not bound to the accepted parent: '+relative)
        pin.update(expected_sha256=expected,bound_to_accepted_parent=True)
    first = _load(base / _FIRST_COLLAR)
    bridge = _load(base / _BRIDGE)
    family = parent.get("source_family")
    if not isinstance(family, Mapping):
        raise ValueError("parent source family is missing")
    guards=_source_guards(repo,parent,first,bridge)
    if guards['passed'] is not True:
        raise ValueError('Specific same-generator source AST/convention guard failed: '+str(guards))

    ctx, iv = _point_context()
    samples = first.get("samples")
    if not isinstance(samples, list):
        raise ValueError("first_collar_jets.samples is missing")
    rows = []
    qproofs = []
    for ztext in _Z_TEXT:
        sample = _sample(samples, ztext)
        qpoint, qtruth, qproof = _q_jet(ctx, iv, parent, ztext, bridge,guards)
        qproofs.append(qproof)
        point, truth,nonlinear = _pressure_pair(ctx, iv, sample, qpoint, qtruth, bridge,guards)
        rows.append(
            {
                "Z": ztext,
                "ordinary_Z_orders": [0, 1, 2],
                "source_rows": [encode_pair(ctx, iv, p, t) for p, t in zip(point, truth)],
                "source_binding": "original_macro_pressure_rows.source_Uy at q1",
                "generator": "same_original_core_generator",
                "pressure_sector_only": True,
                "approximation_is_not_exact_source_value": True,
                "source_enclosure_is_not_a_selected_value": True,
                "macro_nonlinear_remainder":nonlinear,
            }
        )

    return {
        "passed": True,
        "source_family": dict(family),
        "rows": rows,
        "row_count": len(rows),
        "source_scope": "same_generator_pressure_sector_q1_factored_macro",
        "full_source_scope": False,
        "blocked": _blocked(),
        "source_pins": pins,
        "source_bindings":guards,
        "comparison_q": {
            "point_proxy": "q=1 with ordinary Z derivatives zero",
            "truth_source": "saved core Phi rows plus symmetric ordinary continuation caps",
            "point_not_selected_from_finite_width": True,
            "correlated_shared_phi_proof": qproofs,
        },
        "finite_width": {
            "width_cap_used_as_remainder_only": True,
            "expell_relative_remainder": "6*cap*exp(6*cap)",
            "saved_width_enclosure": "actual_bridge_integrals.width_enclosure_cap",
        },
        "contract": {
            "accepted_ancestor_replayed": False,
            "midpoint_as_function_value": False,
            "arbitrary_width_as_point": False,
            "raw46_or_switch_alias_used": False,
            "Taylor_factorials_applied": True,
            "source_rows_are_ordinary_derivatives": True,
        },
    }


__all__ = ["build_modal_source_point"]
