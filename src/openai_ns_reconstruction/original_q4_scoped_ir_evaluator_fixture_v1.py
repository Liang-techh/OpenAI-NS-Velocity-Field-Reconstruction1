"""Small synthetic-only checks for ``original_q4_scoped_ir_evaluator_v1``.

The fixture deliberately has no project imports and no graph/archive access.
It accepts the evaluator callable so a caller can run it once after freezing
the implementation.  Every numeric value is an exact ``Fraction`` and the
synthetic oracle exposes all source and integration calls explicitly.
"""

from fractions import Fraction


def _rational(numerator, denominator=1):
    return {
        "op": "exact_rational",
        "numerator": numerator,
        "denominator": denominator,
    }


def _binder(token, name="t"):
    return {
        "kind": "fresh_binder",
        "token": token,
        "source_name": name,
        "source_ref": "synthetic-binder-" + token,
    }


def _bound(name, source_ref):
    return {
        "op": "bound_variable",
        "name": name,
        "resolved": False,
        "source_node_reference": source_ref,
        "source_row": {
            "op": "bound_variable",
            "name": name,
        },
    }


def _integral(source_ref, token, integrand, lower=None, upper=None):
    return {
        "op": "definite_integral",
        "source_ref": source_ref,
        "source_row": {
            "op": "definite_integral",
            "source_ref": source_ref,
        },
        "lower": _rational(0) if lower is None else lower,
        "upper": _rational(1) if upper is None else upper,
        "binder": _binder(token),
        "integrand": integrand,
    }


class _ExactBackend:
    """Fraction backend with no float conversion."""

    def rational(self, numerator, denominator):
        return Fraction(numerator, denominator)

    def unary(self, name, value):
        if name.lower() == "abs":
            return abs(value)
        raise ValueError("unexpected unary")


class _SyntheticOracle:
    mode = "synthetic_reference"

    def __init__(self, parameter_value=Fraction(3), output=None):
        self.parameter_value = parameter_value
        self.output = output
        self.calls = []

    def parameter(self, row):
        self.calls.append(("parameter", row))
        if self.output is not None:
            return self.output
        if row.get("name") != "synthetic_parameter":
            raise ValueError("unknown synthetic parameter")
        return self.parameter_value

    def component(self, row, rho, z_value):
        self.calls.append(("component", row, rho, z_value))
        if self.output is not None:
            return self.output
        if row.get("quantity") != "Phi":
            raise ValueError("unknown synthetic component")
        if row.get("provider_refs", {}).get("radial") != "fake-radial-provider":
            raise ValueError("component provider was not passed opaquely")
        return rho + z_value

    def integrate(self, function, lower, upper):
        self.calls.append(("integrate", lower, upper))
        # Simpson's three exact nodes is exact for the linear/quadratic
        # polynomial cases used below, including the nested integral.
        midpoint = (lower + upper) / 2
        return (upper - lower) * (
            function(lower) + 4 * function(midpoint) + function(upper)
        ) / 6


class _OriginalOracle(_SyntheticOracle):
    mode = "original_source"

    def __init__(self, family, graph_hash, certificate=None):
        super().__init__()
        self.source_family = family
        self.graph_hash = graph_hash
        if certificate is not None:
            self.admission_certificate = certificate


def _evaluated(result, expected):
    return (
        isinstance(result, dict)
        and result.get("status") == "evaluated"
        and result.get("value") == expected
        and result.get("certification") is False
        and result.get("numeric_errors_verified") is False
        and result.get("provenance", {}).get("caller_provenance_independently_verified")
        is False
    )


def _rejected(result):
    return (
        isinstance(result, dict)
        and result.get("status") == "rejected"
        and result.get("value") is None
        and result.get("certification") is False
    )


def run_fixture(evaluate):
    """Run bounded synthetic checks against an injected ``evaluate`` API.

    The returned receipt contains no claim about an archived source or a
    numerical error bound.  It only records whether the caller-supplied
    evaluator obeyed the synthetic contract.
    """

    cases = {}
    if not callable(evaluate):
        return {
            "status": "failed",
            "passed": 0,
            "failed": 1,
            "cases": {"callable_evaluate": {"passed": False}},
        }

    backend = _ExactBackend()

    def run(name, expression, oracle, free_variables=None, expected=None):
        result = evaluate(
            expression,
            oracle=oracle,
            free_variables={} if free_variables is None else free_variables,
            backend=backend,
        )
        passed = _evaluated(result, expected)
        cases[name] = {
            "passed": passed,
            "status": result.get("status") if isinstance(result, dict) else None,
            "value": result.get("value") if isinstance(result, dict) and passed else None,
        }
        return result

    # Nested integrals deliberately reuse the source binder name while using
    # different fresh tokens.  The inner lower limit reads the outer token,
    # so bounds must be evaluated in the outer environment.
    outer_token = "outer_t"
    inner_token = "inner_t"
    nested = _integral(
        "outer-integral",
        outer_token,
        _integral(
            "inner-integral",
            inner_token,
            {
                "op": "sum",
                "terms": [_binder(outer_token), _binder(inner_token)],
            },
            lower=_binder(outer_token),
            upper=_rational(1),
        ),
    )
    nested_result = run(
        "nested_integral_same_name_distinct_tokens",
        nested,
        _SyntheticOracle(),
        expected=Fraction(1, 2),
    )

    # A non-zero free variable is resolved by an explicit source-aware key;
    # the integer-looking source reference is never used as a node index.
    free_expression = {
        "op": "sum",
        "terms": [_bound("x", "free-x-source"), _rational(1)],
    }
    free_result = run(
        "nonzero_source_aware_free_variable",
        free_expression,
        _SyntheticOracle(),
        free_variables={("x", "free-x-source"): Fraction(7)},
        expected=Fraction(8),
    )

    # Name-only lookup is allowed only with an explicit lexical certificate.
    lexical_result = run(
        "lexical_free_variable_with_certificate",
        free_expression,
        _SyntheticOracle(),
        free_variables={
            "values": {"x": Fraction(5)},
            "lexical_certificate": {"names": ["x"]},
        },
        expected=Fraction(6),
    )

    source_expression = {
        "op": "product",
        "factors": [
            {
                "op": "current_original_source_parameter",
                "source_ref": "parameter-source",
                "source_row": {
                    "op": "current_original_source_parameter",
                    "name": "synthetic_parameter",
                    "Z_independent": True,
                },
            },
            {
                "op": "current_original_unique_analytic_core_component",
                "source_ref": "component-source",
                "source_row": {
                    "op": "current_original_unique_analytic_core_component",
                    "quantity": "Phi",
                },
                "provider": "fake-component-provider",
                "provider_refs": {"radial": "fake-radial-provider"},
                "quantity": "Phi",
                "arguments": {
                    "rho": _rational(2),
                    "Z_variable": _rational(3),
                },
            },
        ],
    }
    source_result = run(
        "opaque_parameter_and_component_oracle",
        source_expression,
        _SyntheticOracle(),
        expected=Fraction(15),
    )

    unary_result = run(
        "restricted_backend_unary_and_exact_arithmetic",
        {
            "op": "analytic_unary",
            "function": "abs",
            "argument": {
                "op": "negative",
                "argument": _rational(3),
            },
        },
        _SyntheticOracle(),
        expected=Fraction(3),
    )

    # The evaluator must reject an operation that v4 should have lowered.
    unknown = evaluate(
        {"op": "not_an_operation"},
        oracle=_SyntheticOracle(),
        free_variables={},
        backend=backend,
    )
    cases["reject_unknown_op"] = {"passed": _rejected(unknown), "status": unknown.get("status")}

    missing = evaluate(
        _rational(1),
        oracle=None,
        free_variables={},
        backend=backend,
    )
    cases["reject_missing_oracle"] = {"passed": _rejected(missing), "status": missing.get("status")}

    dictionary_output = evaluate(
        {
            "op": "current_original_source_parameter",
            "source_row": {"op": "current_original_source_parameter", "name": "synthetic_parameter"},
        },
        oracle=_SyntheticOracle(output={"saved_bound": Fraction(1)}),
        free_variables={},
        backend=backend,
    )
    cases["reject_dict_or_saved_bound_output"] = {
        "passed": _rejected(dictionary_output),
        "status": dictionary_output.get("status"),
    }

    list_output = evaluate(
        {
            "op": "current_original_source_parameter",
            "source_row": {"op": "current_original_source_parameter", "name": "synthetic_parameter"},
        },
        oracle=_SyntheticOracle(output=[Fraction(1)]),
        free_variables={},
        backend=backend,
    )
    cases["reject_list_output"] = {
        "passed": _rejected(list_output),
        "status": list_output.get("status"),
    }

    # Original-source admission includes every integral source_ref/token/name
    # receipt and an explicit (empty here) free-variable namespace receipt.
    family = {"name": "synthetic-q4-family"}
    graph_hash = "a" * 64
    original_ir = _integral(
        "original-integral-source",
        "original_t",
        _binder("original_t"),
    )
    certificate = {
        "source_family": family,
        "graph_hash": graph_hash,
        "capture_safe": True,
        "binder_identity": {"original_t": "opaque-integral-identity"},
        "binder_receipts": [
            {
                "source_ref": "original-integral-source",
                "token": "original_t",
                "source_name": "t",
            }
        ],
        "free_variable_namespace": [],
    }
    original_wrapper = {
        "ir": original_ir,
        "source_family": family,
        "graph_hash": graph_hash,
    }
    original_result = run(
        "original_source_complete_admission",
        original_wrapper,
        _OriginalOracle(family, graph_hash, certificate),
        expected=Fraction(1, 2),
    )

    mismatch_oracle = _OriginalOracle({"name": "other-family"}, graph_hash, certificate)
    mismatch = evaluate(
        original_wrapper,
        oracle=mismatch_oracle,
        free_variables={},
        backend=backend,
    )
    cases["reject_original_family_hash_mismatch"] = {
        "passed": _rejected(mismatch),
        "status": mismatch.get("status"),
    }

    missing_certificate = evaluate(
        original_wrapper,
        oracle=_OriginalOracle(family, graph_hash),
        free_variables={},
        backend=backend,
    )
    cases["reject_missing_original_admission_certificate"] = {
        "passed": _rejected(missing_certificate),
        "status": missing_certificate.get("status"),
    }

    # A receipt that omits the integral source relation must fail closed even
    # when capture_safe is true and the token map looks plausible.
    incomplete_certificate = dict(certificate)
    incomplete_certificate["binder_receipts"] = []
    incomplete = evaluate(
        original_wrapper,
        oracle=_OriginalOracle(family, graph_hash, incomplete_certificate),
        free_variables={},
        backend=backend,
    )
    cases["reject_incomplete_binder_receipt"] = {
        "passed": _rejected(incomplete),
        "status": incomplete.get("status"),
    }

    passed = sum(1 for item in cases.values() if item.get("passed") is True)
    failed = len(cases) - passed
    return {
        "status": "passed" if failed == 0 else "failed",
        "passed": passed,
        "failed": failed,
        "cases": cases,
        "limitations": [
            "Synthetic Fraction checks do not authenticate any original graph.",
            "The evaluator must keep certification false and numeric errors unverified.",
        ],
    }


__all__ = ["run_fixture"]
