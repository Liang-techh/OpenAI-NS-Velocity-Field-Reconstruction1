"""Bounded, read-only C5 mixed-phase topology checker.

This helper reads the frozen C5 JSON graph as data.  It does not import or
execute any source module.  The local expression algebra mirrors only the
ordinary six-row algebra and the normalized MixedTaylor recurrence needed by
the two new actual phase charts.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable


CHARTS = ("restoration", "postrestore")
ORDER = 5
INDICES = tuple((k, j) for k in range(ORDER + 1) for j in range(ORDER + 1 - k))
ZERO = ("rat", 0, 1)
ONE = ("rat", 1, 1)


class _ActualGraph:
    """Canonical view of raw graph rows, retaining source operation shape."""

    def __init__(self, nodes: list[dict[str, Any]]):
        self.nodes = nodes
        self._memo: dict[int, tuple] = {}
        self.opaque_inputs: set[int] = set()

    def pin_inputs(self, values: Iterable[int]) -> None:
        """Keep saved q/p2/parameter rows as typed source atoms."""
        self.opaque_inputs.update(values)

    def key(self, index: int) -> tuple:
        if type(index) is not int or not 0 <= index < len(self.nodes):
            raise ValueError("Invalid source graph handle")
        if index in self._memo:
            return self._memo[index]
        if index in self.opaque_inputs:
            out = ("raw", index)
            self._memo[index] = out
            return out
        row = self.nodes[index]
        op = row.get("operation")
        if op == "exact_rational":
            out = ("rat", int(row["numerator"]), int(row["denominator"]))
        elif op in ("sum", "product"):
            out = (op, tuple(self.key(child) for child in row["arguments"]))
        elif op == "negative":
            out = ("negative", self.key(row["argument"]))
        elif op in ("positive_quotient", "positive_function_quotient"):
            out = (
                op,
                self.key(row["numerator"]),
                self.key(row["denominator"]),
                row.get("source_positive_certificate", row.get("positive_certificate")),
            )
        elif op == "analytic_unary":
            out = ("analytic_unary", row.get("name"), self.key(row["argument"]))
        elif op == "bound_variable":
            out = ("bound_variable", row.get("name"))
        elif op == "mathematical_pi":
            out = ("mathematical_pi",)
        else:
            # Inputs q, p2, eta, d_star, and t0 are already pinned source
            # rows.  Treat an unrecognized source row as an opaque atom; this
            # keeps the check local to the phase recurrence and avoids
            # silently widening into a general graph verifier.
            out = ("raw", index)
        self._memo[index] = out
        return out


def _rat(value: Any) -> tuple:
    value = Fraction(value)
    return ("rat", value.numerator, value.denominator)


def _add(*values: tuple) -> tuple:
    values = tuple(value for value in values if value != ZERO)
    return ZERO if not values else values[0] if len(values) == 1 else ("sum", values)


def _mul(*values: tuple) -> tuple:
    if any(value == ZERO for value in values):
        return ZERO
    values = tuple(value for value in values if value != ONE)
    return ONE if not values else values[0] if len(values) == 1 else ("product", values)


def _neg(value: tuple) -> tuple:
    return ZERO if value == ZERO else ("negative", value)


def _sub(left: tuple, right: tuple) -> tuple:
    return _add(left, _neg(right))


def _quotient(numerator: tuple, denominator: tuple, certificate: str) -> tuple:
    return ZERO if numerator == ZERO else ("positive_quotient", numerator, denominator, certificate)


def _unary(name: str, argument: tuple) -> tuple:
    return ("analytic_unary", name, argument)


def _fixed(value: tuple) -> list[tuple]:
    return [value] + [ZERO] * ORDER


def _jet_add(*jets: list[tuple]) -> list[tuple]:
    return [_add(*(jet[k] for jet in jets)) for k in range(ORDER + 1)]


def _jet_neg(jet: list[tuple]) -> list[tuple]:
    return [_neg(value) for value in jet]


def _jet_sub(left: list[tuple], right: list[tuple]) -> list[tuple]:
    return _jet_add(left, _jet_neg(right))


def _jet_mul(left: list[tuple], right: list[tuple]) -> list[tuple]:
    return [
        _add(
            *(
                _mul(_rat(math.comb(k, j)), left[j], right[k - j])
                for j in range(k + 1)
            )
        )
        for k in range(ORDER + 1)
    ]


def _jet_scale(jet: list[tuple], value: Any) -> list[tuple]:
    return _jet_mul(jet, _fixed(_rat(value) if not isinstance(value, tuple) else value))


def _jet_div(left: list[tuple], right: list[tuple], certificate: str) -> list[tuple]:
    out: list[tuple] = []
    for k in range(ORDER + 1):
        rest = _add(
            *(
                _mul(_rat(math.comb(k, j)), right[j], out[k - j])
                for j in range(1, k + 1)
            )
        )
        out.append(_quotient(_sub(left[k], rest), right[0], certificate))
    return out


def _jet_sqrt(jet: list[tuple]) -> list[tuple]:
    out = [_unary("positive_sqrt", jet[0])]
    for k in range(1, ORDER + 1):
        rest = _add(
            *(
                _mul(_rat(math.comb(k, j)), out[j], out[k - j])
                for j in range(1, k)
            )
        )
        out.append(
            _quotient(
                _sub(jet[k], rest),
                _mul(_rat(2), out[0]),
                "original active positive square root; lazy flat extension outside active domain",
            )
        )
    return out


def _mx_constant(value: tuple) -> dict[tuple[int, int], tuple]:
    return {(0, 0): value}


def _mx_axial(jet: list[tuple]) -> dict[tuple[int, int], tuple]:
    return {
        (k, 0): _mul(_rat(Fraction(1, math.factorial(k))), jet[k])
        for k in range(ORDER + 1)
    }


def _mx_add(*values: dict[tuple[int, int], tuple]) -> dict[tuple[int, int], tuple]:
    return {ij: _add(*(value.get(ij, ZERO) for value in values)) for ij in INDICES}


def _mx_scale(value: dict[tuple[int, int], tuple], ref: tuple) -> dict[tuple[int, int], tuple]:
    return {ij: _mul(ref, entry) for ij, entry in value.items()}


def _mx_mul(left: dict[tuple[int, int], tuple], right: dict[tuple[int, int], tuple]) -> dict[tuple[int, int], tuple]:
    return {
        (k, j): _add(
            *(
                _mul(left.get((a, b), ZERO), right.get((k - a, j - b), ZERO))
                for a in range(k + 1)
                for b in range(j + 1)
            )
        )
        for k, j in INDICES
    }


def _mx_div(left: dict[tuple[int, int], tuple], right: dict[tuple[int, int], tuple]) -> dict[tuple[int, int], tuple]:
    out: dict[tuple[int, int], tuple] = {}
    for degree in range(ORDER + 1):
        for k in range(degree + 1):
            j = degree - k
            rest = _add(
                *(
                    _mul(right.get((a, b), ZERO), out.get((k - a, j - b), ZERO))
                    for a in range(k + 1)
                    for b in range(j + 1)
                    if a + b > 0
                )
            )
            out[k, j] = _quotient(
                _sub(left.get((k, j), ZERO), rest),
                right[0, 0],
                "1-2*r*cos(psi)+r^2=(cos(psi)-r)^2+sin(psi)^2>0 since |r|<1",
            )
    return out


def _mx_trig(angle: tuple) -> dict[tuple[int, int], tuple]:
    cosine = _unary("cos", angle)
    sine = _unary("sin", angle)
    cycle = (cosine, _neg(sine), _neg(cosine), sine)
    return {
        (0, j): _mul(_rat(Fraction(1, math.factorial(j))), cycle[j % 4])
        for j in range(ORDER + 1)
    }


def _phase_expected(actual: _ActualGraph, data: dict[str, Any], chart: str) -> dict[tuple[int, int], tuple]:
    bindings = data["prescribed_loop_C5_bindings"]
    chart_data = bindings["charts"][chart]
    loop = chart_data["loop"]
    roots = chart_data["roots"]
    atom = actual.key
    eta, d_star = atom(bindings["eta"]), atom(bindings["dstar"])
    q = [atom(value) for value in loop["q"]]
    p2 = [atom(value) for value in roots["p2"]]
    t0 = [atom(value) for value in roots["t0"]]

    u = _jet_scale(_jet_mul(p2, q), _quotient(ONE, d_star, "same positive canonical d_star"))
    hinv = _jet_div(
        _fixed(ONE),
        _jet_sqrt(_jet_add(_fixed(ONE), _jet_mul(u, u))),
        "positive sqrt(1+u^2)",
    )
    radius = _jet_mul(u, hinv)
    alpha = _jet_scale(_jet_mul(q, hinv), 2)

    angle = atom(loop["angle"])
    rr = _mx_axial(radius)
    cc = _mx_trig(angle)
    numerator = _mx_add(cc, _mx_scale(rr, _rat(-1)))
    denominator = _mx_add(
        _mx_constant(_rat(1)),
        _mx_scale(_mx_mul(rr, cc), _rat(-2)),
        _mx_mul(rr, rr),
    )
    t = _mx_add(
        _mx_axial(t0),
        _mx_mul(_mx_axial(alpha), _mx_div(numerator, denominator)),
    )
    t_squared = _mx_mul(t, t)

    ratio = _jet_add(
        _fixed(ONE),
        _jet_mul(t0, t0),
        _jet_scale(_jet_mul(q, q), 2),
    )
    two_pi = _mul(_rat(2), ("mathematical_pi",))
    K = _jet_div(
        _fixed(ONE),
        _jet_scale(ratio, two_pi),
        "same positive 2*pi*(1+t0^2+2*q^2)",
    )
    lam = _mx_mul(_mx_axial(K), _mx_add(_mx_constant(_rat(1)), t_squared))

    expected: dict[tuple[int, int], tuple] = {}
    phi = [atom(value) for value in loop["Phi"]]
    for k in range(ORDER + 1):
        expected[k, 0] = _mul(_rat(Fraction(1, math.factorial(k))), phi[k])
    for k, j in INDICES:
        if j:
            expected[k, j] = _mul(_rat(Fraction(1, j)), lam[k, j - 1])
    return expected


def verify_phase_mixed(
    cache_path: str | Path,
    charts: Iterable[str] = CHARTS,
) -> dict[str, Any]:
    """Return bounded row-level results for the two new actual phase charts."""

    cache_path = Path(cache_path)
    charts = tuple(charts)
    source_path = cache_path.with_name(cache_path.name.removesuffix(".json.gz") + ".py")
    with gzip.open(cache_path, "rt", encoding="utf-8") as handle:
        data = json.load(handle)
    actual = _ActualGraph(data["exact_graph_nodes"])
    bindings = data["prescribed_loop_C5_bindings"]
    pinned = {bindings["eta"], bindings["dstar"]}
    for chart in charts:
        chart_data = bindings["charts"][chart]
        pinned.update(chart_data["loop"]["q"])
        pinned.update(chart_data["roots"]["p2"])
        pinned.update(chart_data["roots"]["t0"])
    pinned = {index for index in pinned
              if data["exact_graph_nodes"][index]["operation"] != "exact_rational"}
    actual.pin_inputs(pinned)
    out: dict[str, Any] = {
        "cache": str(cache_path),
        "cache_sha256": hashlib.sha256(cache_path.read_bytes()).hexdigest(),
        "source": str(source_path),
        "source_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
        "source_executed": False,
        "charts": {},
        "total_rows": 0,
        "matched_rows": 0,
        "passed": True,
        "first_mismatch": None,
    }
    for chart in charts:
        loop = bindings["charts"][chart]["loop"]
        expected = _phase_expected(actual, data, chart)
        raw_rows = loop["phase_partial_coefficients"]
        if (len(raw_rows) != len(INDICES)
                or any(type(row.get(key)) is not int for row in raw_rows
                       for key in ("Z", "angle", "node"))
                or {(row["Z"], row["angle"]) for row in raw_rows} != set(INDICES)
                or len({row["node"] for row in raw_rows}) != len(INDICES)):
            raise ValueError("Actual mixed phase table has repeated or malformed rows")
        actual_rows = {
            (row["Z"], row["angle"]): actual.key(row["node"])
            for row in raw_rows
        }
        mismatches = []
        for index in INDICES:
            if index not in actual_rows:
                mismatches.append({"index": list(index), "reason": "missing actual row"})
            elif expected[index] != actual_rows[index]:
                mismatches.append({"index": list(index), "reason": "topology mismatch"})
        matched = len(INDICES) - len(mismatches)
        lambda_indices = tuple(index for index in INDICES if index[1] > 0)
        lambda_mismatches = [
            index for index in lambda_indices
            if index not in actual_rows or expected[index] != actual_rows[index]
        ]
        out["total_rows"] += len(INDICES)
        out["matched_rows"] += matched
        out["charts"][chart] = {
            "actual_rows": len(actual_rows),
            "expected_rows": len(INDICES),
            "matched_rows": matched,
            "raw_lambda_rows": len(lambda_indices),
            "raw_lambda_rows_matched": len(lambda_indices) - len(lambda_mismatches),
            "raw_lambda_convolution_bound": not lambda_mismatches,
            "mismatches": mismatches,
            "rows": [dict(Z=row["Z"], angle=row["angle"], node=row["node"],
                          passed=expected[row["Z"], row["angle"]]
                          == actual_rows[row["Z"], row["angle"]])
                     for row in raw_rows],
        }
        if mismatches and out["first_mismatch"] is None:
            out["first_mismatch"] = {"chart": chart, **mismatches[0]}
    out["passed"] = out["matched_rows"] == out["total_rows"]
    return out


BASE = "experiments/root_st073/"
PARENT = BASE + "lei_ren_part1_paper_compliant_current_original_long_reference_modified_seams.json"
PARENT_SHA256 = "204d8b45aff263a3df1eae319dd7534e304162794ac375a1c0ec2cfe9beafc35"
SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py"
CACHE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.json.gz"
MODULE = "src/openai_ns_reconstruction/original_restoration_phase_mixed.py"
SOURCE_GUARDS = (
    "self.indices=tuple((k,j) for k in range(order+1) for j in range(order+1-k))",
    "lam=m.mul(m.axial(K),m.add(m.constant(1),t2))",
    "P=a.add(a.fixed(angle),T2);Phi=a.mul(K,P)",
    "true_loop_C5_original_integrand=True,angle_endpoints_independent_of_Z=True)",
    "partials={(k,0):g.mul(c(Fraction(1,math.factorial(k))),rows(Phi)[k]) for k in range(6)}",
    "partials[k,j]=g.mul(c(Fraction(1,j)),lam[k,j-1])",
    "psi,delta,at=implicit_inverse(g,partials,inverse,angle)",
)


def build(root: str | Path, parent: dict[str, Any]) -> dict[str, Any]:
    """Check only the 42 newly consumed restoration phase coefficients."""
    root = Path(root)
    parent_bytes = (root / PARENT).read_bytes()
    if (hashlib.sha256(parent_bytes).hexdigest() != PARENT_SHA256
            or json.loads(parent_bytes) != parent or parent["checks"]["passed"] is not True):
        raise ValueError("Exact accepted modified-seams predecessor required")
    pins = {}
    for path in (PARENT, SOURCE, CACHE, MODULE):
        actual = hashlib.sha256((root / path).read_bytes()).hexdigest()
        expected = (PARENT_SHA256 if path == PARENT else actual if path == MODULE
                    else parent["input_hashes"][Path(path).name])
        if actual != expected:
            raise ValueError("Changed source: " + path)
        pins[path] = dict(path=path, expected_sha256=expected, actual_sha256=actual,
                          passed=True, source_executed=False)
    source = (root / SOURCE).read_text(encoding="utf-8")
    guards = [dict(token=token, passed=token in source) for token in SOURCE_GUARDS]
    if not all(guard["passed"] for guard in guards):
        raise ValueError("Actual mixed coefficient generator changed")
    receipt = verify_phase_mixed(root / CACHE, CHARTS)
    if receipt["passed"] is not True or receipt["total_rows"] != 42:
        raise ValueError("Actual restoration phase coefficient mismatch: "
                         + str(receipt["first_mismatch"]))
    return dict(passed=True, parent_report=PARENT, parent_report_sha256=PARENT_SHA256,
                source_family=parent["source_family"], source_pins=pins,
                source_guards=guards, sourceguards_tuple=True, source_executed=False,
                accepted_chart_roots_recomputed=False, accepted_builders_executed=False,
                charts=list(CHARTS), phase_coefficients=receipt,
                normalized_mixed_orders="Z+angle<=5; 1/k! axial and 1/j angular",
                limits=dict(physical_seed_identity_discharged=False,
                            positive_loop_IFT_discharged=False,
                            modified_Qy_seams_proved=False,
                            global_five_moment_closure=False, full_goal_complete=False))
