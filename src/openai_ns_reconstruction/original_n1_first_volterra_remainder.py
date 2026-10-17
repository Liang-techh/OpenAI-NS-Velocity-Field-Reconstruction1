"""Source-bound finite-word residual bridge for the actual first n1 chart.

The saved first-chart graph contains the actual six-row operator and its
complete regular Volterra solution.  This module only binds that source and
provides the elementary finite-word bookkeeping needed by a later directed
enclosure.  It deliberately does not turn the normalized Green majorant in
the older reports into a raw first-chart value or tail.

For a natural word depth ``n`` let

    W_0 = (I + J) G g,
    W_(k+1) = (I + J) G (C_0 W_k + B_1 d_Z W_k),
    S_n = sum_(k=0)^n W_k.

The accepted raw equation gives

    L S_n - g = -(C_0 + B_1 d_Z) W_n.

Given source-bound raw first-chart constants ``M0``, ``q`` and ``r`` with
``||W_0|| <= M0``, ``||H|| <= q < 1`` and
``||C_0 + B_1 d_Z|| <= r``, the bridge returns
``M0*q**(n+1)/(1-q)`` for the omitted words and
``r*M0*q**n`` for the residual.  No such raw constants are silently taken
from the old normalized/potential-kernel reports.
"""

from __future__ import annotations

import ast
import gzip
import hashlib
import json
import re
from fractions import Fraction
from pathlib import Path
from typing import Any, Mapping, Sequence


SOURCE_FAMILY = {
    "actual_five_defect_family_sha256": "3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894",
    "implicit_source_sha256": "5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae",
    "datum_enclosure_sha256": "dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010",
}

PARENT_REPORT = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_native_collar_transport.json"
PARENT_REPORT_SHA256 = "7c05f7438a27d9caa4378cb89b9ab6d33b44138086e86292b73f2fd8e3f4bf96"

INNER_SOURCE = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_inner_volterra_functions.py"
INNER_ARCHIVE = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_inner_volterra_functions.json.gz"
INNER_CHECK = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_inner_volterra_functions_check.json"
INNER_SOURCE_SHA256 = "ff82f4edb03922a551cba3e2df20217a55dc80ed42bd39ccb8b77b83f63b0208"
INNER_ARCHIVE_SHA256 = "1b44be8c6cbfffb692c7a0430ecbc8a590f46001a233c99c59ca8c5e2c44a6fe"
INNER_CHECK_SHA256 = "3cfdc3ce6a7e8c2a98a31fae985c02a2acb4e28025d313ee35db4cffc6e7b47c"

REGULAR_SOURCE = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_regular_series_function.py"
REGULAR_ARCHIVE = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_regular_series_function.json.gz"
REGULAR_CHECK = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_regular_series_function_check.json"
REGULAR_SOURCE_SHA256 = "b64b7cae8529cfc5caa158fe06ac5660ec1d2217fc4d845d02002dac3a195555"
REGULAR_ARCHIVE_SHA256 = "25486283c104fc3679337ad8da686ae4879c122a3c7bd3d7f3ba275999209860"
REGULAR_CHECK_SHA256 = "1071e24e5fc37fb413d95190573839583f2fb9d4b6de180a3e5c9769a959d53c"

OPERATOR_SOURCE = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_inner_operator.py"
OPERATOR_ARCHIVE = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_inner_operator.json.gz"
OPERATOR_CHECK = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_inner_operator_check.json"
OPERATOR_SOURCE_SHA256 = "268217c6b33bc717472742fa64d76a853c6760a956abfad92ac45bcf385df014"
OPERATOR_ARCHIVE_SHA256 = "63b181010f1715d3b20589b06c1c2139ef63bcc73d7b0f3c5b7c43d9c5e03f52"
OPERATOR_CHECK_SHA256 = "846cfedcbd0355a4a8e6a5ed445558d367cd55c7459085dd9c5d74821d969bd4"

COUPLED_SOURCE = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_coupled_first_step.py"
COUPLED_ARCHIVE = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_coupled_first_step.json.gz"
COUPLED_CHECK = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_coupled_first_step_check.json"
COUPLED_SOURCE_SHA256 = "97913b82b8570b220921cd89918569a4019ff3495bdbea0548b3c53e9d8c93e7"
COUPLED_ARCHIVE_SHA256 = "30b807af4147658ec0978ebd253cc1b2b31bed67bfed430a965e66461450add1"
COUPLED_CHECK_SHA256 = "6a58caf1db41697620820830cbe7f76620874a8b429b351ad3d58f637e4e1416"

ANALYTIC_DOMAIN_SOURCE = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_inner_analytic_domain.py"
ANALYTIC_DOMAIN_ARCHIVE = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_inner_analytic_domain.json.gz"
ANALYTIC_DOMAIN_CHECK = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_inner_analytic_domain_check.json"
ANALYTIC_DOMAIN_SOURCE_SHA256 = "56507e29d1c069c6ca71c13722717dfaeef662265401652afb3b62232bfe73f3"
ANALYTIC_DOMAIN_ARCHIVE_SHA256 = "c9905f7a9d13d4f651d1f55036932ea6c69db4fc5fa890e67a9839e611af4db7"
ANALYTIC_DOMAIN_CHECK_SHA256 = "01a4a5f497a224d3dd6ebcb7b44809df5b3702e8d32028170a764903238e0ae2"

OLD_NUMERIC_SOURCE = "src/openai_ns_reconstruction/original_n1_core_volterra_numeric.py"
OLD_CUTOFF_SOURCE = "src/openai_ns_reconstruction/original_n1_cutoff_numeric.py"
OLD_M15_SOURCE = "src/openai_ns_reconstruction/original_n1_inner_m15_numeric.py"
OLD_NUMERIC_SOURCE_SHA256 = "19317da946c95430b6fa9c7915b2d063e5b186c87dd77b96c4d9df3c14691adc"
OLD_CUTOFF_SOURCE_SHA256 = "7a2004b378ee495c6c053d04a8394406937f3a6fca235a721138f1ce4661c731"
OLD_M15_SOURCE_SHA256 = "004341ab4eb1f4693995bc86a748eae4c20ff2b5406fb0d784c4416963f5988b"

EXACT_INNER_WORD_DEFINITION = "W=sum_k>=0 H^k W0; W0=(I+J)Gg; H=(I+J)G(C0+B1*dZ)"
EXACT_REGULAR_WORD_DEFINITION = "W=sum_k>=0 V^k W0; W0=Gprincipal*g; V f=Gprincipal*(C0*f+B1*dZ(f))"
EXACT_RESIDUAL_IDENTITY = "L(sum_0^N Wk)-g=-(C0*WN+B1*dZ(WN))"
FIELDS = ("F1", "Uz1", "K1", "P1", "d_x_F1", "d_x_Uz1")
ORDINARY_Z_ORDERS = (0, 1, 2, 3)


def _repo_root(root: str | Path | None) -> Path:
    if root is None:
        return Path(__file__).resolve().parents[2]
    return Path(root).resolve()


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _mapping_digest(value: Any) -> str:
    """Digest a JSON mapping without treating its upstream parent field as identity."""
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("utf-8")
    ).hexdigest()


def _load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _load_gzip_json(path: Path) -> Any:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return json.load(handle)


def _pin(root: Path, relative: str, expected: str | None = None) -> dict[str, Any]:
    path = root / relative
    result: dict[str, Any] = {
        "path": relative,
        "exists": path.is_file(),
        "sha256": None,
        "expected_sha256": expected,
        "sha256_match": False,
    }
    if path.is_file():
        result["sha256"] = _sha256(path)
        result["sha256_match"] = expected is None or result["sha256"] == expected
    return result


def _ast_inventory(path: Path) -> dict[str, Any]:
    if not path.is_file():
        return {"exists": False, "functions": [], "classes": []}
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    return {
        "exists": True,
        "functions": sorted(
            node.name
            for node in ast.walk(tree)
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        ),
        "classes": sorted(node.name for node in ast.walk(tree) if isinstance(node, ast.ClassDef)),
    }


def _contains_all(text: str, fragments: Sequence[str]) -> bool:
    return all(fragment in text for fragment in fragments)


def _fraction(value: Any, *, label: str) -> Fraction:
    """Parse an exact integer/decimal input; reject binary floats."""
    if isinstance(value, bool):
        raise TypeError(f"{label}: boolean is not a bound")
    if isinstance(value, Fraction):
        return value
    if isinstance(value, int):
        return Fraction(value, 1)
    if isinstance(value, (tuple, list)):
        # The saved source reports use mpmath's exact tuple
        # ``(sign, mantissa, exponent, bitcount)``.  Accept only bounded
        # tuples so the ordinary huge N/Mbase scales stay lazy/source-only.
        if len(value) != 4 or any(isinstance(item, bool) or not isinstance(item, int) for item in value):
            raise TypeError(f"{label}: invalid exact mpf tuple")
        sign, mantissa, exponent, bitcount = value
        if sign not in (0, 1) or mantissa < 0 or bitcount < 0:
            raise ValueError(f"{label}: invalid exact mpf tuple signs")
        if mantissa.bit_length() > 10000 or abs(exponent) > 100000:
            raise ValueError(f"{label}: exact mpf tuple too large for bounded arithmetic")
        result = Fraction(mantissa, 1)
        if exponent >= 0:
            result *= 2 ** exponent
        else:
            result /= 2 ** (-exponent)
        return -result if sign else result
    if isinstance(value, str):
        text = value.strip()
        if not text:
            raise ValueError(f"{label}: empty bound")
        # Saved interval reports intentionally retain huge exponential scales.
        # Keep them as source records; do not materialize their enormous exact
        # integers merely because a caller requested the generic arithmetic API.
        match = re.search(r"[eE]([+-]?\d+)$", text)
        if match and abs(int(match.group(1))) > 10000:
            raise ValueError(f"{label}: decimal exponent too large for exact arithmetic")
        return Fraction(text)
    raise TypeError(f"{label}: exact integer, Fraction, or decimal string required")


def _upper(value: Any, *, label: str) -> Fraction:
    """Read an upper endpoint from a source-bound record or exact scalar."""
    if isinstance(value, Mapping):
        for key in ("upper_exact_mpf_tuple", "exact_upper_tuple", "upper", "upper_decimal", "value", "exact_upper"):
            if key in value:
                return _fraction(value[key], label=label)
        raise ValueError(f"{label}: source record has no upper endpoint")
    return _fraction(value, label=label)


def _lower(value: Any, *, label: str) -> Fraction:
    """Read a positive lower endpoint for a denominator or Cauchy gap."""
    if isinstance(value, Mapping):
        for key in ("lower_exact_mpf_tuple", "exact_lower_tuple", "lower", "lower_decimal", "value", "exact_lower"):
            if key in value:
                return _fraction(value[key], label=label)
        raise ValueError(f"{label}: source record has no lower endpoint")
    return _fraction(value, label=label)


def source_bound(
    upper: int | str | Fraction,
    *,
    path: str,
    sha256: str,
    ordinary_Z_order: int,
    chart: str = "first",
    units: str = "raw_first_chart",
) -> dict[str, Any]:
    """Construct an explicit source-bound input record for the formula API.

    The record is metadata, not an evaluated field value.  ``build`` never
    creates one from a normalized tail report.
    """
    if not isinstance(ordinary_Z_order, int) or ordinary_Z_order < 0:
        raise ValueError("ordinary_Z_order must be a nonnegative integer")
    if not path or not sha256:
        raise ValueError("source path and digest are required")
    _fraction(upper, label="upper")
    return {
        "upper": str(upper),
        "source_path": path,
        "source_sha256": sha256,
        "ordinary_Z_order": ordinary_Z_order,
        "chart": chart,
        "units": units,
        "source_bound": True,
        "point_value": False,
    }


def finite_word_residual_formula(
    natural_word_depth: int,
    *,
    source_norm: Any = None,
    operator_ratio: Any = None,
    residual_multiplier: Any = None,
    ordinary_Z_order: int = 0,
    legacy_polynomial_degree: int | None = None,
) -> dict[str, Any]:
    """Return a deterministic finite-word tail/residual formula.

    ``natural_word_depth`` indexes the Volterra words.  It is intentionally
    independent of a legacy radial polynomial degree.  If raw source-bound
    constants are absent, symbolic formulas and the required inputs are
    returned.  Numeric arithmetic is performed only when all three explicit
    constants are supplied as exact scalars or source-bound records.
    """
    if isinstance(natural_word_depth, bool) or not isinstance(natural_word_depth, int):
        raise TypeError("natural_word_depth must be an integer")
    if natural_word_depth < 0:
        raise ValueError("natural_word_depth must be nonnegative")
    if isinstance(ordinary_Z_order, bool) or not isinstance(ordinary_Z_order, int):
        raise TypeError("ordinary_Z_order must be an integer")
    if ordinary_Z_order < 0 or ordinary_Z_order > 3:
        raise ValueError("ordinary_Z_order must lie in 0..3")
    if legacy_polynomial_degree is not None:
        if isinstance(legacy_polynomial_degree, bool) or not isinstance(legacy_polynomial_degree, int):
            raise TypeError("legacy_polynomial_degree must be an integer when supplied")
        if legacy_polynomial_degree < 0:
            raise ValueError("legacy_polynomial_degree must be nonnegative")

    symbolic = {
        "word_term": f"M0_Z{ordinary_Z_order}*q_Z{ordinary_Z_order}^{natural_word_depth}",
        "omitted_words": (
            f"M0_Z{ordinary_Z_order}*q_Z{ordinary_Z_order}^{natural_word_depth + 1}/"
            f"(1-q_Z{ordinary_Z_order})"
        ),
        "finite_residual": f"r_Z{ordinary_Z_order}*M0_Z{ordinary_Z_order}*q_Z{ordinary_Z_order}^{natural_word_depth}",
    }
    result: dict[str, Any] = {
        "natural_word_depth": natural_word_depth,
        "ordinary_Z_order": ordinary_Z_order,
        "legacy_polynomial_degree": legacy_polynomial_degree,
        "word_recurrence": "W_(k+1)=(I+J)G(C0*W_k+B1*dZ(W_k))",
        "finite_sum": f"S_{natural_word_depth}=sum_(k=0)^{natural_word_depth} W_k",
        "residual_identity": EXACT_RESIDUAL_IDENTITY,
        "symbolic": symbolic,
        "raw_first_chart_inputs": {
            "source_norm": source_norm,
            "operator_ratio": operator_ratio,
            "residual_multiplier": residual_multiplier,
        },
        "numeric_evaluation": False,
        "tail_certified": False,
        "required_source_inputs": [
            f"raw_first_chart M0 bound for ordinary Z order {ordinary_Z_order}",
            f"raw_first_chart H ratio q < 1 for ordinary Z order {ordinary_Z_order}",
            f"raw_first_chart residual multiplier r for ordinary Z order {ordinary_Z_order}",
        ],
    }
    if source_norm is None or operator_ratio is None or residual_multiplier is None:
        return result

    m0 = _upper(source_norm, label="source_norm")
    q = _upper(operator_ratio, label="operator_ratio")
    residual = _upper(residual_multiplier, label="residual_multiplier")
    if m0 < 0 or q < 0 or residual < 0:
        raise ValueError("nonnegative source/operator/residual bounds required")
    if q >= 1:
        raise ValueError("operator_ratio upper endpoint must be below one")
    word = m0 * q**natural_word_depth
    omitted = m0 * q ** (natural_word_depth + 1) / (1 - q)
    remainder = residual * word
    result.update(
        {
            "numeric_evaluation": True,
            "tail_certified": True,
            "word_term_exact_fraction": {"numerator": word.numerator, "denominator": word.denominator},
            "omitted_words_exact_fraction": {"numerator": omitted.numerator, "denominator": omitted.denominator},
            "finite_residual_exact_fraction": {"numerator": remainder.numerator, "denominator": remainder.denominator},
            "numeric_note": "This is a generic evaluation from caller-supplied raw source bounds; build supplies none.",
        }
    )
    return result


def _source_interval_summary(value: Any, *, include_exact: bool = False) -> Any:
    """Keep saved endpoints, optionally retaining bounded exact MPF tuples."""
    if not isinstance(value, Mapping):
        return value
    result = {
        key: value[key]
        for key in ("lower", "upper", "width")
        if key in value
    }
    if include_exact:
        for key in ("lower_exact_mpf_tuple", "upper_exact_mpf_tuple"):
            if key in value:
                result[key] = list(value[key]) if isinstance(value[key], (tuple, list)) else value[key]
    return result


def _saved_upper_below_one(value: Any) -> bool:
    """Check ``upper < 1`` without materializing a huge exponential bound."""
    if isinstance(value, Mapping):
        value = value.get("upper")
    if not isinstance(value, str):
        try:
            return _fraction(value, label="saved upper") < 1
        except (TypeError, ValueError, OverflowError):
            return False
    text = value.strip()
    match = re.search(r"^([+-]?(?:\d+(?:\.\d*)?|\.\d+))[eE]([+-]?\d+)$", text)
    if match:
        mantissa = match.group(1)
        exponent = int(match.group(2))
        try:
            if Fraction(mantissa) <= 0:
                return False
        except (TypeError, ValueError, OverflowError):
            return False
        return exponent < 0
    try:
        return _fraction(text, label="saved upper") < 1
    except (TypeError, ValueError, OverflowError):
        return False


def _saved_lower_positive(value: Any) -> bool:
    """Check a saved lower endpoint is positive without expanding its scale."""
    if isinstance(value, Mapping):
        value = value.get("lower")
    if not isinstance(value, str):
        try:
            return _fraction(value, label="saved lower") > 0
        except (TypeError, ValueError, OverflowError):
            return False
    text = value.strip()
    match = re.search(r"^([+-]?(?:\d+(?:\.\d*)?|\.\d+))[eE]([+-]?\d+)$", text)
    if match:
        try:
            return Fraction(match.group(1)) > 0
        except (TypeError, ValueError, OverflowError):
            return False
    try:
        return _fraction(text, label="saved lower") > 0
    except (TypeError, ValueError, OverflowError):
        return False


def preconditioned_finite_word_tail_formula(
    natural_word_depth: int,
    *,
    source_tail: Mapping[str, Any] | None = None,
    base_bound: Any = None,
    kinematic_bound: Any = None,
    off_kinematic_bound: Any = None,
    derivative_bound: Any = None,
    z_loss: Any = None,
    ratio_upper: Any = None,
    reference_ratio_upper: Any = None,
    omitted_upper: Any = None,
    residual_multiplier: Any = None,
    ordinary_Z_gap: Any = None,
    ordinary_Z_gap_reference: str | None = None,
    ordinary_Z_order: int = 0,
) -> dict[str, Any]:
    """Describe the raw ``(I+J)G`` preconditioned word majorant.

    The coupled first-step source resums the actual raw Akin feedback.  Its
    ``Mbase``, ``a``, ``Cp``, ``C1`` and ``loss`` therefore have different
    provenance from the older ``Gprincipal``/Kchi regular-series report.  A
    caller may provide the saved source tail record; the record is retained as
    an enclosure certificate, while its very large exponential endpoints are
    never converted into Python integers.  Small exact caller inputs may still
    evaluate the finite word term and, when an exact ratio is supplied, the
    geometric remainder.
    """
    if isinstance(natural_word_depth, bool) or not isinstance(natural_word_depth, int):
        raise TypeError("natural_word_depth must be an integer")
    if natural_word_depth < 0:
        raise ValueError("natural_word_depth must be nonnegative")
    if isinstance(ordinary_Z_order, bool) or not isinstance(ordinary_Z_order, int):
        raise TypeError("ordinary_Z_order must be an integer")
    if ordinary_Z_order not in ORDINARY_Z_ORDERS:
        raise ValueError("ordinary_Z_order must lie in 0..3")

    symbolic = {
        "word_term": "Mbase*[a*(1+a)*(Cp+C1*k/loss)]^k/k!",
        "word_term_at_depth": (
            f"Mbase*[a*(1+a)*(Cp+C1*{natural_word_depth}/loss)]^{natural_word_depth}/"
            f"{natural_word_depth}!"
        ),
        "ratio_at_depth": (
            f"beta_{natural_word_depth}=e*a*(1+a)*(Cp/(%d)+C1/loss)"
            % (natural_word_depth + 1)
        ),
        "omitted_words": (
            f"Mbase*beta_{natural_word_depth}^{natural_word_depth + 1}/"
            f"(1-beta_{natural_word_depth})"
        ),
        "ordinary_Z_derivative_factor": (
            "1" if ordinary_Z_order == 0 else
            f"{ordinary_Z_order}!/(Z_gap)^{ordinary_Z_order}"
        ),
    }
    result: dict[str, Any] = {
        "natural_word_depth": natural_word_depth,
        "ordinary_Z_order": ordinary_Z_order,
        "raw_units": "actual_first_chart_preconditioned_Akin",
        "raw_preconditioned_Akin_resummation": True,
        "normalized_Kchi": False,
        "source_formula": symbolic,
        "source_inputs": {
            "Mbase": base_bound,
            "a": kinematic_bound,
            "Cp": off_kinematic_bound,
            "C1": derivative_bound,
            "loss": z_loss,
            "ratio_upper": ratio_upper,
            "reference_ratio_upper": reference_ratio_upper,
            "omitted_upper": omitted_upper,
            "residual_multiplier": residual_multiplier,
        },
        "numeric_evaluation": False,
        "word_term_numeric": False,
        "geometric_remainder_numeric": False,
        "tail_certified": False,
        "source_tail_record_bound": False,
        "required_source_inputs": [
            "initial_Gg_common_tube_bound from the actual inner analytic domain",
            "sqrt(real_rho_upper_bound) for the actual common radial domain",
            "Cp and C1 from perturbation_bound applied to actual core plus first-collar rows",
            "positive common Z loss and the source beta_N admission",
        ],
    }
    result["ordinary_Z_derivative_rule"] = (
        "coupled derivative_tail.Cauchy_gap=solution_Z_radius/2; use p!/(Cauchy_gap)^p"
    )
    result["ordinary_Z_gap_input"] = ordinary_Z_gap
    result["ordinary_Z_gap_source_reference"] = ordinary_Z_gap_reference
    result["ordinary_Z_tail_bound_symbolic"] = (
        f"({symbolic['ordinary_Z_derivative_factor']})*"
        f"({symbolic['omitted_words']})"
    )
    result["uniform_N_ge_one_tail_symbolic"] = (
        "Mbase*beta_1^(N+1)/(1-beta_1), for N>=1 and beta_N<=beta_1"
    )

    source_tail_contract = (
        isinstance(source_tail, Mapping)
        and source_tail.get("retained_preconditioned_indices")
        and all(isinstance(index, int) and not isinstance(index, bool)
                for index in source_tail.get("retained_preconditioned_indices", ()))
        and natural_word_depth >= max(source_tail.get("retained_preconditioned_indices", ()))
        and source_tail.get("initial_base_common_tube_bound") is not None
        and source_tail.get("kinematic_J_norm_upper") is not None
        and source_tail.get("I_plus_J_norm_upper") is not None
        and source_tail.get("off_kinematic_B0_norm_upper") is not None
        and source_tail.get("B1_norm_upper") is not None
        and source_tail.get("total_available_Z_loss") is not None
        and source_tail.get("geometric_tail_ratio_upper") is not None
        and _saved_upper_below_one(source_tail.get("geometric_tail_ratio_upper"))
        and source_tail.get("omitted_sum_bound") is not None
        and source_tail.get("equal_nested_radius_loss_for_length_k") == "loss/k"
        and source_tail.get("arbitrary_k_term_formula") == symbolic["word_term"]
        and source_tail.get("factorial_to_geometric_formula") == "k! >= (k/e)^k"
        and source_tail.get("all_k_derivative_words_allowed") is True
        and source_tail.get("half_derivative_block_reduction_used") is False
        and source_tail.get("old_unpreconditioned_tail_reused") is False
    )
    if source_tail_contract:
        result["source_tail_record_bound"] = True
        result["tail_certified"] = True
        retained_max = max(source_tail.get("retained_preconditioned_indices", ()))
        result["tail_bound_mode"] = (
            "saved_beta_1_uniform_for_deeper_N" if natural_word_depth > retained_max else "saved_retained_index"
        )
        result["source_tail_record"] = {
            key: _source_interval_summary(source_tail[key])
            for key in (
                "retained_preconditioned_indices",
                "initial_base_common_tube_bound",
                "kinematic_J_norm_upper",
                "I_plus_J_norm_upper",
                "off_kinematic_B0_norm_upper",
                "B1_norm_upper",
                "total_available_Z_loss",
                "geometric_tail_ratio_upper",
                "omitted_sum_bound",
            )
            if key in source_tail
        }
        result["source_tail_scope"] = source_tail.get("scope")

    # Exact finite-term arithmetic is useful for moderate caller-provided
    # bounds.  Saved huge decimal exponent endpoints are deliberately kept as
    # records, so a ValueError here means "symbolic/source-bound", not failure
    # of the underlying source certificate.
    scalar_inputs = (base_bound, kinematic_bound, off_kinematic_bound, derivative_bound, z_loss)
    if all(value is not None for value in scalar_inputs):
        try:
            mbase = _upper(base_bound, label="base_bound")
            a = _upper(kinematic_bound, label="kinematic_bound")
            cp = _upper(off_kinematic_bound, label="off_kinematic_bound")
            c1 = _upper(derivative_bound, label="derivative_bound")
            loss = _lower(z_loss, label="z_loss")
            if min(mbase, a, cp, c1) < 0 or loss <= 0:
                raise ValueError("raw preconditioned majorants must be nonnegative and loss positive")
            import math

            term = mbase * (a * (1 + a) * (cp + c1 * natural_word_depth / loss)) ** natural_word_depth
            term /= math.factorial(natural_word_depth)
            result.update(
                {
                    "word_term_numeric": True,
                    "word_term_exact_fraction": {"numerator": term.numerator, "denominator": term.denominator},
                }
            )
            if residual_multiplier is not None:
                residual = _upper(residual_multiplier, label="residual_multiplier")
                if residual < 0:
                    raise ValueError("residual_multiplier must be nonnegative")
                residual_term = residual * term
                result["finite_residual_exact_fraction"] = {
                    "numerator": residual_term.numerator,
                    "denominator": residual_term.denominator,
                }
        except (TypeError, ValueError, OverflowError) as exc:
            result["numeric_input_error"] = str(exc)

    if ratio_upper is not None and base_bound is not None:
        try:
            mbase = _upper(base_bound, label="base_bound")
            beta = _upper(ratio_upper, label="ratio_upper")
            if mbase < 0 or beta < 0 or beta >= 1:
                raise ValueError("geometric ratio must satisfy 0 <= beta < 1")
            omitted = mbase * beta ** (natural_word_depth + 1) / (1 - beta)
            result.update(
                {
                    "geometric_remainder_numeric": True,
                    "geometric_remainder_exact_fraction": {
                        "numerator": omitted.numerator,
                        "denominator": omitted.denominator,
                    },
                    "tail_certified": True,
                }
            )
        except (TypeError, ValueError, OverflowError) as exc:
            result["numeric_ratio_input_error"] = str(exc)

    # The saved N=1 admission also certifies every deeper N by beta_N <= beta_1.
    # Keep this as a separate coarse bound; it does not pretend that beta_N was
    # numerically materialized for the huge source exponent.
    reference = reference_ratio_upper
    if reference is None and source_tail_contract and natural_word_depth >= 1:
        reference = source_tail.get("geometric_tail_ratio_upper")
        result["uniform_reference_ratio_source"] = "saved N=1 geometric_tail_ratio_upper"
    if reference is not None and base_bound is not None and natural_word_depth >= 1:
        try:
            mbase = _upper(base_bound, label="base_bound")
            beta_ref = _upper(reference, label="reference_ratio_upper")
            if mbase < 0 or beta_ref < 0 or beta_ref >= 1:
                raise ValueError("reference geometric ratio must satisfy 0 <= beta_1 < 1")
            uniform_omitted = mbase * beta_ref ** (natural_word_depth + 1) / (1 - beta_ref)
            result.update(
                {
                    "uniform_geometric_remainder_numeric": True,
                    "uniform_geometric_remainder_exact_fraction": {
                        "numerator": uniform_omitted.numerator,
                        "denominator": uniform_omitted.denominator,
                    },
                }
            )
        except (TypeError, ValueError, OverflowError) as exc:
            result["uniform_ratio_input_error"] = str(exc)

    remainder_record = result.get("geometric_remainder_exact_fraction")
    if remainder_record is None:
        remainder_record = result.get("uniform_geometric_remainder_exact_fraction")
    if ordinary_Z_order and ordinary_Z_gap is not None and remainder_record is not None:
        try:
            import math

            gap = _lower(ordinary_Z_gap, label="ordinary_Z_gap")
            if gap <= 0:
                raise ValueError("ordinary_Z_gap must be positive")
            factor = Fraction(math.factorial(ordinary_Z_order), 1) / gap ** ordinary_Z_order
            z_tail = factor * Fraction(remainder_record["numerator"], remainder_record["denominator"])
            result["ordinary_Z_tail_exact_fraction"] = {
                "numerator": z_tail.numerator,
                "denominator": z_tail.denominator,
            }
            result["ordinary_Z_tail_basis"] = (
                "beta_N" if "geometric_remainder_exact_fraction" in result else "beta_1"
            )
        except (TypeError, ValueError, OverflowError) as exc:
            result["ordinary_Z_gap_input_error"] = str(exc)

    if omitted_upper is not None:
        result["caller_supplied_omitted_upper"] = _source_interval_summary(omitted_upper)
    return result


def finite_word_bridge(
    natural_word_depth: int,
    source_inputs: Mapping[int, Mapping[str, Any]] | None = None,
    *,
    ordinary_Z_orders: Sequence[int] = (0, 1),
) -> dict[str, Any]:
    """Apply the generic formula independently to requested ordinary Z orders."""
    source_inputs = source_inputs or {}
    orders = tuple(ordinary_Z_orders)
    if not orders or any(order not in ORDINARY_Z_ORDERS for order in orders):
        raise ValueError("ordinary_Z_orders must be a nonempty subset of 0..3")
    rows = {}
    for order in orders:
        inputs = source_inputs.get(order, {})
        rows[str(order)] = finite_word_residual_formula(
            natural_word_depth,
            source_norm=inputs.get("source_norm"),
            operator_ratio=inputs.get("operator_ratio"),
            residual_multiplier=inputs.get("residual_multiplier"),
            ordinary_Z_order=order,
            legacy_polynomial_degree=inputs.get("legacy_polynomial_degree"),
        )
    complete = all(row["numeric_evaluation"] for row in rows.values())
    return {
        "natural_word_depth": natural_word_depth,
        "ordinary_Z_orders": list(orders),
        "rows": rows,
        "all_raw_inputs_supplied": complete,
        "tail_certified": complete and all(row["tail_certified"] for row in rows.values()),
        "legacy_polynomial_is_not_word_depth": True,
    }


def _graph_checks(graph: list[Any], bindings: Mapping[str, Any], family: Mapping[str, str]) -> dict[str, Any]:
    function_checks = []
    for field in FIELDS:
        rows = bindings.get("functions", {}).get("first", {}).get(field, [])
        field_ok = len(rows) == 4
        for order, node_id in enumerate(rows):
            if not isinstance(node_id, int) or not (0 <= node_id < len(graph)):
                field_ok = False
                continue
            node = graph[node_id]
            field_ok = field_ok and node.get("operation") == "current_original_unique_regular_n1_component"
            field_ok = field_ok and node.get("chart") == "first" and node.get("quantity") == field
            field_ok = field_ok and node.get("coordinate_domain_lower") == 0
            field_ok = field_ok and node.get("coordinate_domain_upper") == 35824
            field_ok = field_ok and node.get("ordinary_slow_Z_derivative_order") == order
            field_ok = field_ok and node.get("source_family") == dict(family)
            field_ok = field_ok and node.get("exact_function_definition") == EXACT_INNER_WORD_DEFINITION
            field_ok = field_ok and node.get("all_words_retained_in_analytic_limit") is True
            field_ok = field_ok and node.get("raw_units_not_normalized_potential_kernel") is True
            field_ok = field_ok and node.get("source_parameter_bounds_not_selected_as_values") is True
        function_checks.append({"field": field, "rows": rows, "passed": field_ok})

    first_system = bindings.get("systems", {}).get("first", {})
    required_b0 = {"1,5", "2,6", "3,6", "4,1", "5,1", "5,2", "5,3", "5,5", "6,1", "6,2", "6,3", "6,4", "6,6"}
    required_b1 = {"5,1", "5,2", "5,3", "6,2", "6,3", "6,4"}
    operator_rows = (
        set(first_system.get("B0", {})) >= required_b0
        and set(first_system.get("B1", {})) >= required_b1
        and set(first_system.get("g", {})) == {"4", "5", "6"}
        and first_system.get("D") == [0, 0, 2, 0, 3, 1]
    )
    definition = bindings.get("definition", {})
    semantic_flags = {
        "raw_C0_is_B0_minus_Akin": definition.get("raw_C0_is_B0_minus_Akin") is True,
        "one_core_to_collar_path_no_seam_reset": definition.get("one_core_to_collar_path_no_seam_reset") is True,
        "first_domain_upper_is_sc_quarter": definition.get("common_first_phase_upper") == 35824,
        "raw_units_distinct_from_normalized_kernel": definition.get("normalized_Kchi_not_used_with_raw_matrix") is True,
        "complete_limit_distinct_from_finite_words": definition.get("full_regular_solution_distinct_from_finite_words_and_tail_enclosures") is True,
        "numeric_point_accuracy_not_claimed": definition.get("numerical_original_n1_point_accuracy_claimed") is False,
    }
    return {
        "first_function_rows": function_checks,
        "first_function_rows_all_orders_0_to_3": all(row["passed"] for row in function_checks),
        "first_operator_rows_present": operator_rows,
        "first_operator_row_counts": {"B0": len(first_system.get("B0", {})), "B1": len(first_system.get("B1", {})), "g": len(first_system.get("g", {}))},
        "definition_flags": semantic_flags,
        "definition_flags_all": all(semantic_flags.values()),
        "parameters": bindings.get("parameters", {}),
        "coordinates": bindings.get("coordinates", {}),
    }


def build(root: str | Path | None = None, parent: Any = None) -> dict[str, Any]:
    """Read saved source artifacts and return a bounded receipt.

    No source owner, numerical producer, recurrence, fixture, or quadrature
    is imported or executed.  ``parent`` may be a path or a mapping carrying a
    parent digest; the pinned native-collar report remains authoritative.
    """
    repo = _repo_root(root)
    pins = {
        "native_collar_parent": _pin(repo, PARENT_REPORT, PARENT_REPORT_SHA256),
        "inner_volterra_source": _pin(repo, INNER_SOURCE, INNER_SOURCE_SHA256),
        "inner_volterra_archive": _pin(repo, INNER_ARCHIVE, INNER_ARCHIVE_SHA256),
        "inner_volterra_check": _pin(repo, INNER_CHECK, INNER_CHECK_SHA256),
        "regular_series_source": _pin(repo, REGULAR_SOURCE, REGULAR_SOURCE_SHA256),
        "regular_series_archive": _pin(repo, REGULAR_ARCHIVE, REGULAR_ARCHIVE_SHA256),
        "regular_series_check": _pin(repo, REGULAR_CHECK, REGULAR_CHECK_SHA256),
        "inner_operator_source": _pin(repo, OPERATOR_SOURCE, OPERATOR_SOURCE_SHA256),
        "inner_operator_archive": _pin(repo, OPERATOR_ARCHIVE, OPERATOR_ARCHIVE_SHA256),
        "inner_operator_check": _pin(repo, OPERATOR_CHECK, OPERATOR_CHECK_SHA256),
        "coupled_tail_source": _pin(repo, COUPLED_SOURCE, COUPLED_SOURCE_SHA256),
        "coupled_tail_archive": _pin(repo, COUPLED_ARCHIVE, COUPLED_ARCHIVE_SHA256),
        "coupled_tail_check": _pin(repo, COUPLED_CHECK, COUPLED_CHECK_SHA256),
        "analytic_domain_source": _pin(repo, ANALYTIC_DOMAIN_SOURCE, ANALYTIC_DOMAIN_SOURCE_SHA256),
        "analytic_domain_archive": _pin(repo, ANALYTIC_DOMAIN_ARCHIVE, ANALYTIC_DOMAIN_ARCHIVE_SHA256),
        "analytic_domain_check": _pin(repo, ANALYTIC_DOMAIN_CHECK, ANALYTIC_DOMAIN_CHECK_SHA256),
        "legacy_numeric_polynomial_source": _pin(repo, OLD_NUMERIC_SOURCE, OLD_NUMERIC_SOURCE_SHA256),
        "legacy_cutoff_source": _pin(repo, OLD_CUTOFF_SOURCE, OLD_CUTOFF_SOURCE_SHA256),
        "legacy_m15_source": _pin(repo, OLD_M15_SOURCE, OLD_M15_SOURCE_SHA256),
    }
    pins_ok = all(item["exists"] and item["sha256_match"] for item in pins.values())

    supplied_parent_sha = None
    supplied_parent_mapping_digest = None
    pinned_parent_mapping_digest = None
    if isinstance(parent, (str, Path)):
        parent_path = Path(parent)
        if not parent_path.is_absolute():
            parent_path = repo / parent_path
        if parent_path.is_file():
            supplied_parent_sha = _sha256(parent_path)
    elif isinstance(parent, Mapping):
        # The mapping's own content is the parent identity.  Its
        # ``parent_report_sha256`` field records an upstream predecessor and
        # must not be compared with this module's pinned report digest.
        supplied_parent_mapping_digest = _mapping_digest(parent)
        pinned_path = repo / PARENT_REPORT
        if pinned_path.is_file():
            try:
                pinned_parent_mapping_digest = _mapping_digest(_load_json(pinned_path))
            except (OSError, ValueError, json.JSONDecodeError):
                pinned_parent_mapping_digest = None
        parent_argument_ok = (
            supplied_parent_mapping_digest is not None
            and pinned_parent_mapping_digest is not None
            and supplied_parent_mapping_digest == pinned_parent_mapping_digest
        )
    else:
        parent_argument_ok = parent is None
    if not isinstance(parent, Mapping) and "parent_argument_ok" not in locals():
        parent_argument_ok = supplied_parent_sha is None or supplied_parent_sha == PARENT_REPORT_SHA256

    checks: dict[str, Any] = {
        "source_pins_present_and_unchanged": pins_ok,
        "parent_argument_matches_pinned_report": parent_argument_ok,
        "accepted_parent_sha256": PARENT_REPORT_SHA256,
        "supplied_parent_mapping_digest": supplied_parent_mapping_digest,
        "pinned_parent_mapping_digest": pinned_parent_mapping_digest,
        "source_family_literal": dict(SOURCE_FAMILY),
    }

    missing: list[str] = []
    try:
        inner = _load_gzip_json(repo / INNER_ARCHIVE)
        inner_check = _load_json(repo / INNER_CHECK)
        regular = _load_gzip_json(repo / REGULAR_ARCHIVE)
        regular_check = _load_json(repo / REGULAR_CHECK)
        coupled = _load_gzip_json(repo / COUPLED_ARCHIVE)
        coupled_check = _load_json(repo / COUPLED_CHECK)
        analytic_domain = _load_gzip_json(repo / ANALYTIC_DOMAIN_ARCHIVE)
        analytic_domain_check = _load_json(repo / ANALYTIC_DOMAIN_CHECK)
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        inner = inner_check = regular = regular_check = None
        coupled = coupled_check = analytic_domain = analytic_domain_check = None
        missing.append(f"saved source archive could not be read: {exc}")

    if inner is not None:
        graph = inner.get("exact_graph_nodes", [])
        bindings = inner.get("inner_Volterra_bindings", {})
        graph_check = _graph_checks(graph, bindings, SOURCE_FAMILY)
        checks.update(
            {
                "inner_graph_node_count_49024": len(graph) == 49024,
                "inner_graph_source_family": inner.get("source_family") == SOURCE_FAMILY,
                "inner_graph_actual_raw_operator_complete_limit": inner.get("actual_raw_six_row_operator_and_complete_regular_limit") is True,
                "inner_graph_checks": graph_check,
                "inner_archive_does_not_claim_numeric_point_accuracy": (
                    inner.get("numerical_original_n1_point_accuracy_claimed") is False
                    or inner.get("inner_Volterra_bindings", {}).get("definition", {}).get(
                        "numerical_original_n1_point_accuracy_claimed"
                    ) is False
                ),
            }
        )
        if not checks["inner_graph_node_count_49024"]:
            missing.append("saved inner Volterra graph node count 49,024")
    else:
        checks["inner_graph_checks"] = {}

    if inner_check is not None:
        checks.update(
            {
                "inner_check_all_passed": inner_check.get("all_passed") is True,
                "inner_check_function_gate": inner_check.get("current_original_n1_true_inner_Volterra_function_adapter_installed") is True,
                "inner_check_numeric_oracle_not_executed": inner_check.get("original_numeric_n1_solution_oracle_executed") is False,
                "inner_check_numeric_point_accuracy_false": inner_check.get("numerical_original_n1_point_accuracy_claimed") is False,
            }
        )
    else:
        checks.update({"inner_check_all_passed": False, "inner_check_function_gate": False})

    if regular is not None:
        contract = regular.get("contract", {})
        checks.update(
            {
                "regular_contract_family": regular.get("source_family") == SOURCE_FAMILY,
                "regular_contract_word_definition": contract.get("exact_function_definition") == EXACT_REGULAR_WORD_DEFINITION,
                "regular_contract_residual_identity": contract.get("finite_equation_residual_identity") == EXACT_RESIDUAL_IDENTITY,
                "regular_contract_zeroth_index_only": contract.get("finite_word_index") == 0,
                "regular_contract_keeps_omitted_words": contract.get("all_remaining_words_enclosed_not_set_zero") is True,
                "regular_contract_separates_centers": contract.get("numeric_centers_compute_only_zeroth_source") is True,
                "regular_contract_tail_is_normalized_reference": contract.get("normalization") == "Wphysical=Snew*W; original F0 Z derivative applied once",
            }
        )
    else:
        checks["regular_contract_word_definition"] = False
        checks["regular_contract_residual_identity"] = False

    # The accepted coupled tail is a raw Akin-preconditioned majorant.  Bind
    # its five inputs to the saved analytic-domain report and the source call
    # chain before exposing it to a later evaluator.  The older regular-series
    # report remains a separate normalized Kchi reference above.
    analytic_tail = None
    derivative_cauchy_gap = None
    if analytic_domain is not None:
        analytic_domain_check_ok = (
            analytic_domain.get("source_family") == SOURCE_FAMILY
            and analytic_domain.get("initial_Gg_common_tube_bound") is not None
            and analytic_domain.get("actual_inner_geometry", {}).get("real_rho_upper_bound") is not None
            and analytic_domain.get("common_tubes", {}).get("outer_Z_radius") is not None
            and analytic_domain.get("common_tubes", {}).get("matrix_Z_radius") is not None
            and analytic_domain.get("common_tubes", {}).get("Gg_Z_radius") is not None
            and analytic_domain.get("common_tubes", {}).get("solution_Z_radius") is not None
            and analytic_domain.get("denominator_bounds", {}).get("L_modulus_lower") is not None
        )
        checks.update(
            {
                "analytic_domain_source_family": analytic_domain.get("source_family") == SOURCE_FAMILY,
                "analytic_domain_required_bound_records": analytic_domain_check_ok,
                "analytic_domain_check_all_passed": analytic_domain_check.get("all_passed") is True
                if analytic_domain_check is not None
                else False,
                "analytic_domain_common_tube_guard": analytic_domain_check.get("common_tube_nozero_guards_passed") is True
                if analytic_domain_check is not None
                else False,
                "analytic_domain_no_computed_iteration_terms": analytic_domain_check.get("bounds_do_not_supply_computed_iteration_terms") is True
                if analytic_domain_check is not None
                else False,
                "analytic_domain_width_not_materialized": analytic_domain_check.get("source_width_is_not_a_materialized_parameter") is True
                if analytic_domain_check is not None
                else False,
                "analytic_domain_ordinary_Z_derivative_rule": analytic_domain.get("ordinary_Z_derivative_rule")
                == "for fields on the Gg tube use k!/(eta/4)^k times corresponding middle-tube bound; solution tube reserves another eta/8",
            }
        )
    else:
        checks.update(
            {
                "analytic_domain_source_family": False,
                "analytic_domain_required_bound_records": False,
                "analytic_domain_check_all_passed": False,
                "analytic_domain_common_tube_guard": False,
                "analytic_domain_no_computed_iteration_terms": False,
                "analytic_domain_width_not_materialized": False,
                "analytic_domain_ordinary_Z_derivative_rule": False,
            }
        )

    if coupled is not None:
        coupled_samples = coupled.get("samples", {})
        first_sample = coupled_samples.get("actual_first_keep")
        if not isinstance(first_sample, Mapping):
            first_sample = None
        analytic_tail = first_sample.get("analytic_preconditioned_tail") if first_sample else None
        derivative_tail = first_sample.get("derivative_tail") if first_sample else None
        derivative_cauchy_gap = derivative_tail.get("Cauchy_gap") if isinstance(derivative_tail, Mapping) else None
        tail_shape_ok = isinstance(analytic_tail, Mapping) and all(
            analytic_tail.get(key) is not None
            for key in (
                "initial_base_common_tube_bound",
                "kinematic_J_norm_upper",
                "I_plus_J_norm_upper",
                "off_kinematic_B0_norm_upper",
                "B1_norm_upper",
                "total_available_Z_loss",
                "geometric_tail_ratio_upper",
                "omitted_sum_bound",
            )
        )
        tail_formula_ok = isinstance(analytic_tail, Mapping) and (
            analytic_tail.get("equal_nested_radius_loss_for_length_k") == "loss/k"
            and analytic_tail.get("arbitrary_k_term_formula") == "Mbase*[a*(1+a)*(Cp+C1*k/loss)]^k/k!"
            and analytic_tail.get("factorial_to_geometric_formula") == "k! >= (k/e)^k"
            and analytic_tail.get("all_k_derivative_words_allowed") is True
            and analytic_tail.get("half_derivative_block_reduction_used") is False
            and analytic_tail.get("old_unpreconditioned_tail_reused") is False
        )
        tail_ratio_admitted = (
            isinstance(analytic_tail, Mapping)
            and _saved_upper_below_one(analytic_tail.get("geometric_tail_ratio_upper"))
            and "if ep(beta)[1]>=1:raise ValueError('Preconditioned geometric tail not admitted')" in (
                (repo / COUPLED_SOURCE).read_text(encoding="utf-8")
                if (repo / COUPLED_SOURCE).is_file()
                else ""
            )
        )
        derivative_gap_record = (
            isinstance(derivative_cauchy_gap, Mapping)
            and _saved_lower_positive(derivative_cauchy_gap)
            and _saved_lower_positive(derivative_cauchy_gap.get("upper"))
        )
        checks.update(
            {
                "coupled_tail_source_family": coupled.get("source_family") == SOURCE_FAMILY,
                "coupled_tail_check_all_passed": coupled_check.get("all_passed") is True
                if coupled_check is not None
                else False,
                "coupled_tail_first_iteration_gate": coupled_check.get("current_original_n1_first_coupled_regular_iteration_computed") is True
                if coupled_check is not None
                else False,
                "coupled_tail_actual_first_keep_present": first_sample is not None,
                "coupled_tail_record_shape": tail_shape_ok,
                "coupled_tail_formula_contract": tail_formula_ok,
                "coupled_tail_ratio_admitted": tail_ratio_admitted,
                "coupled_tail_derivative_cauchy_gap_record": derivative_gap_record,
                "coupled_tail_derivative_gap_source": _contains_all(
                    (repo / COUPLED_SOURCE).read_text(encoding="utf-8")
                    if (repo / COUPLED_SOURCE).is_file()
                    else "",
                    (
                        "derivative_margin=tail['solution_Z_radius']/2",
                        "derivative_tail=dict(value_norm_on_Z_radius=tail['solution_Z_radius']",
                        "Cauchy_gap=derivative_margin",
                    ),
                ),
                "coupled_tail_raw_equation_source": _contains_all(
                    (repo / COUPLED_SOURCE).read_text(encoding="utf-8")
                    if (repo / COUPLED_SOURCE).is_file()
                    else "",
                    (
                        "Wbase=(I+J)Gg",
                        "H=(I+J)G(B0-Akin+B1*d_Z)",
                        "def preconditioned_tail",
                    ),
                ),
            }
        )
    else:
        checks.update(
            {
                "coupled_tail_source_family": False,
                "coupled_tail_check_all_passed": False,
                "coupled_tail_first_iteration_gate": False,
                "coupled_tail_actual_first_keep_present": False,
                "coupled_tail_record_shape": False,
                "coupled_tail_formula_contract": False,
                "coupled_tail_ratio_admitted": False,
                "coupled_tail_derivative_cauchy_gap_record": False,
                "coupled_tail_derivative_gap_source": False,
                "coupled_tail_raw_equation_source": False,
            }
        )

    inner_text = (repo / INNER_SOURCE).read_text(encoding="utf-8") if (repo / INNER_SOURCE).is_file() else ""
    regular_text = (repo / REGULAR_SOURCE).read_text(encoding="utf-8") if (repo / REGULAR_SOURCE).is_file() else ""
    operator_text = (repo / OPERATOR_SOURCE).read_text(encoding="utf-8") if (repo / OPERATOR_SOURCE).is_file() else ""
    coupled_text = (repo / COUPLED_SOURCE).read_text(encoding="utf-8") if (repo / COUPLED_SOURCE).is_file() else ""
    old_numeric_text = (repo / OLD_NUMERIC_SOURCE).read_text(encoding="utf-8") if (repo / OLD_NUMERIC_SOURCE).is_file() else ""
    old_cutoff_text = (repo / OLD_CUTOFF_SOURCE).read_text(encoding="utf-8") if (repo / OLD_CUTOFF_SOURCE).is_file() else ""
    checks["source_ast"] = {
        "inner": _ast_inventory(repo / INNER_SOURCE),
        "regular": _ast_inventory(repo / REGULAR_SOURCE),
        "operator": _ast_inventory(repo / OPERATOR_SOURCE),
        "coupled": _ast_inventory(repo / COUPLED_SOURCE),
        "analytic_domain": _ast_inventory(repo / ANALYTIC_DOMAIN_SOURCE),
    }
    checks["inner_source_definitions"] = _contains_all(
        inner_text,
        ("def raw_matrix_expressions", "def raw_system", "def green_apply", "def feedback_vectors", "class CurrentOriginalN1InnerVolterraFunctions"),
    )
    checks["operator_source_definitions"] = _contains_all(
        operator_text,
        ("def known_forcing", "def matrix", "def first_rows", "equation='d_x W + D W/x = B0 W + B1 d_Z W + g'"),
    )
    checks["coupled_source_tail_provider_is_only_pinned"] = _contains_all(
        coupled_text,
        ("def preconditioned_tail", "all_k_derivative_words_allowed=True", "old_unpreconditioned_tail_reused=False"),
    )
    checks["coupled_tail_raw_input_binding_source"] = _contains_all(
        coupled_text,
        (
            "view['initial_Gg_common_tube_bound']",
            "view['actual_inner_geometry']['real_rho_upper_bound']",
            "rows=perturbation_bound(c,together",
            "rows['Cp'],rows['C1']",
            "view['common_tubes']['outer_Z_radius']",
            "tail=preconditioned_tail(c,view['initial_Gg_common_tube_bound'],c.sqrt(rho_max),rows['Cp'],rows['C1'],eta/8,1)",
            "view['actual_core_radial_bounds']",
            "view['actual_collar_radial_bounds']",
        ),
    )
    checks["raw_operator_is_distinct_from_normalized_Kchi"] = (
        _contains_all(
            coupled_text,
            (
                "Wbase=(I+J)Gg",
                "H=(I+J)G(B0-Akin+B1*d_Z)",
                "old_unpreconditioned_tail_reused=False",
            ),
        )
        and checks.get("inner_graph_checks", {}).get("definition_flags", {}).get(
            "raw_units_distinct_from_normalized_kernel"
        ) is True
        and checks.get("regular_contract_tail_is_normalized_reference") is True
    )
    checks["analytic_domain_source_bound_chain"] = _contains_all(
        (repo / ANALYTIC_DOMAIN_SOURCE).read_text(encoding="utf-8")
        if (repo / ANALYTIC_DOMAIN_SOURCE).is_file()
        else "",
        (
            "def hierarchy_bounds",
            "def neumann_tail",
            "initial_Gg_common_tube_bound=M0",
            "regular_B1_infinity_norm_upper=C1",
            "common_tubes=dict",
            "source_measure_identity='dx=(sqrt(rho)/2)*hb*ds; first-collar phase length 1/2'",
        ),
    )
    checks["legacy_numeric_is_polynomial_not_word_prefix"] = _contains_all(
        old_numeric_text,
        ("Finite forcing polynomials", "not interpreted as", "def _forcing_polynomial", "def _principal", "def _source_tail"),
    ) and _contains_all(old_cutoff_text, ("def _raw_uniform", "def _cutoff", "def _mass", "def _evaluate"))
    checks["legacy_numeric_not_executed"] = True

    if not checks.get("source_pins_present_and_unchanged"):
        missing.append("one or more saved source pins changed or are unavailable")
    if not checks.get("parent_argument_matches_pinned_report"):
        missing.append("parent report argument does not match the pinned native-collar report")
    for key in (
        "inner_graph_source_family",
        "inner_graph_actual_raw_operator_complete_limit",
        "inner_graph_checks",
        "inner_check_all_passed",
        "inner_check_function_gate",
        "regular_contract_family",
        "regular_contract_word_definition",
        "regular_contract_residual_identity",
        "inner_source_definitions",
        "operator_source_definitions",
        "analytic_domain_source_family",
        "analytic_domain_required_bound_records",
        "analytic_domain_check_all_passed",
        "analytic_domain_common_tube_guard",
        "analytic_domain_no_computed_iteration_terms",
        "analytic_domain_width_not_materialized",
        "coupled_tail_source_family",
        "coupled_tail_check_all_passed",
        "coupled_tail_first_iteration_gate",
        "coupled_tail_actual_first_keep_present",
        "coupled_tail_record_shape",
        "coupled_tail_formula_contract",
        "coupled_tail_ratio_admitted",
        "coupled_tail_derivative_cauchy_gap_record",
        "coupled_tail_derivative_gap_source",
        "coupled_tail_raw_equation_source",
        "coupled_tail_raw_input_binding_source",
        "raw_operator_is_distinct_from_normalized_Kchi",
        "analytic_domain_source_bound_chain",
        "legacy_numeric_is_polynomial_not_word_prefix",
    ):
        if checks.get(key) is False:
            missing.append(key)
    graph_checks = checks.get("inner_graph_checks", {})
    if isinstance(graph_checks, dict) and (
        graph_checks.get("first_function_rows_all_orders_0_to_3") is False
        or graph_checks.get("first_operator_rows_present") is False
        or graph_checks.get("definition_flags_all") is False
    ):
        missing.append("actual first-chart graph rows/operator/definition guards")

    ast_inventory = checks["source_ast"]
    checks["source_ast_function_names_bound"] = (
        {"green_apply", "feedback_vectors", "raw_system"}.issubset(ast_inventory["inner"]["functions"])
        and {"known_forcing", "matrix", "first_rows"}.issubset(ast_inventory["operator"]["functions"])
        and {"preconditioned_tail"}.issubset(ast_inventory["coupled"]["functions"])
        and {"regular_series_bound"}.issubset(ast_inventory["regular"]["functions"])
        and {"hierarchy_bounds", "neumann_tail"}.issubset(ast_inventory["analytic_domain"]["functions"])
    )
    if not checks["source_ast_function_names_bound"]:
        missing.append("required first-chart source functions in AST")

    preconditioned_tail_source_bound = all(
        checks.get(key) is True
        for key in (
            "analytic_domain_source_family",
            "analytic_domain_required_bound_records",
            "analytic_domain_check_all_passed",
            "analytic_domain_common_tube_guard",
            "analytic_domain_no_computed_iteration_terms",
            "analytic_domain_width_not_materialized",
            "coupled_tail_source_family",
            "coupled_tail_check_all_passed",
            "coupled_tail_first_iteration_gate",
            "coupled_tail_actual_first_keep_present",
            "coupled_tail_record_shape",
            "coupled_tail_formula_contract",
            "coupled_tail_ratio_admitted",
            "coupled_tail_derivative_cauchy_gap_record",
            "coupled_tail_derivative_gap_source",
            "coupled_tail_raw_equation_source",
            "coupled_tail_raw_input_binding_source",
            "raw_operator_is_distinct_from_normalized_Kchi",
            "analytic_domain_source_bound_chain",
        )
    )
    preconditioned_tail_receipts = {
        str(order): preconditioned_finite_word_tail_formula(
            1,
            source_tail=analytic_tail if isinstance(analytic_tail, Mapping) else None,
            ordinary_Z_gap=(
                _source_interval_summary(derivative_cauchy_gap, include_exact=True)
                if isinstance(derivative_cauchy_gap, Mapping)
                else None
            ),
            ordinary_Z_gap_reference=(
                "coupled.samples.actual_first_keep.derivative_tail.Cauchy_gap"
                if isinstance(derivative_cauchy_gap, Mapping)
                else None
            ),
            ordinary_Z_order=order,
        )
        for order in ORDINARY_Z_ORDERS
    }
    preconditioned_tail_receipt = preconditioned_tail_receipts["0"]
    source_missing = list(missing)
    if preconditioned_tail_source_bound:
        source_missing.extend(
            [
                "ordinary Z order-specific residual multiplier r remains source-bound input for the exact PDE residual formula",
                "signed/numerical first-chart field values and a pointwise finite-word evaluation remain unavailable",
            ]
        )
    else:
        source_missing.extend(
            [
                "raw first-chart preconditioned Mbase, Akin factor, Cp, C1, and loss source receipt",
                "source-bound unit/conjugacy receipt connecting those constants to the actual raw operator",
            ]
        )

    structural_passed = not missing
    return {
        "module": "original_n1_first_volterra_remainder",
        "passed": structural_passed,
        "source_only": True,
        "source_executed": False,
        "natural_word_definition_bound": structural_passed,
        "finite_word_residual_formula_available": True,
        "actual_first_chart_operator_source_bound": structural_passed,
        "ordinary_Z_orders_source_bound": [0, 1, 2, 3],
        "legacy_polynomial_not_promoted_to_word_prefix": checks.get("legacy_numeric_is_polynomial_not_word_prefix", False),
        "raw_first_chart_tail_inputs_available": preconditioned_tail_source_bound,
        "raw_preconditioned_tail_majorant_source_bound": preconditioned_tail_source_bound,
        "certified_raw_preconditioned_tail_majorant": preconditioned_tail_source_bound,
        "actual_all_word_cauchy_majorant_available": preconditioned_tail_source_bound,
        "actual_tail_does_not_use_fixed_H_ratio": True,
        "certified_first_chart_tail": False,
        "numeric_actual_first_chart_field_evaluated": False,
        "normalized_tail_not_used_as_raw_tail": True,
        "old_numeric_producer_executed": False,
        "source_family": dict(SOURCE_FAMILY),
        "parent_report": PARENT_REPORT,
        "parent_report_sha256": PARENT_REPORT_SHA256,
        "source_pins": pins,
        "checks": checks,
        "formula_contract": {
            "word_definition": EXACT_INNER_WORD_DEFINITION,
            "residual_identity": EXACT_RESIDUAL_IDENTITY,
            "omitted_words_bound": "M0*q^(n+1)/(1-q)",
            "finite_residual_bound": "r*M0*q^n",
            "q_requirement": "0 <= q < 1",
            "raw_units_required": True,
            "ordinary_Z_orders_supported": [0, 1, 2, 3],
            "actual_preconditioned_word_term": "Mbase*[a*(1+a)*(Cp+C1*k/loss)]^k/k!",
            "actual_preconditioned_ratio": "beta_N=e*a*(1+a)*(Cp/(N+1)+C1/loss)",
            "actual_preconditioned_omitted_words": "Mbase*beta_N^(N+1)/(1-beta_N)",
            "deeper_N_coarse_tail": "Mbase*beta_1^(N+1)/(1-beta_1), N>=1, because beta_N<=beta_1",
            "ordinary_Z_derivative_factor": "p!/(coupled derivative_tail.Cauchy_gap)^p, Cauchy_gap=solution_Z_radius/2",
        },
        "preconditioned_tail_receipt": preconditioned_tail_receipt,
        "preconditioned_tail_receipts_by_ordinary_Z_order": preconditioned_tail_receipts,
        "missing_inputs": source_missing,
        "limitations": [
            "The saved coupled preconditioned-tail constants are consumed only as the source-bound raw Akin majorant; no normalized Kchi quantity is substituted.",
            "The saved regular-series report retains omitted words but does not supply a raw first-chart point oracle.",
            "The raw preconditioned omitted-word majorant is certified from saved source records, while no signed or numerical actual first-chart field or finite-word value is claimed.",
        ],
    }


__all__ = [
    "SOURCE_FAMILY",
    "PARENT_REPORT",
    "PARENT_REPORT_SHA256",
    "source_bound",
    "finite_word_residual_formula",
    "preconditioned_finite_word_tail_formula",
    "finite_word_bridge",
    "build",
]
