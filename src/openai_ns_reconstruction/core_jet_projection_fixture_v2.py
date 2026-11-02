"""Small synthetic source strings for ``core_jet_projection_v2``.

The fixtures are source-only.  They are intentionally short but retain the
implementation anchors used by the real source leaves: factorial source
normalisation, Taylor convolution/reciprocal/differentiation, the relative-F0
recurrence, the low Dual product/reciprocal, and the explicit low ``fz``
formula.  The mutation cases exercise a signature omission, a denominator
exponent, a recurrence factor, and the low ``fz`` factor.
"""


SYNTHETIC_SOURCE_HIGH = r'''
class IntervalTaylorJet:
    def __init__(self, context, coefficients):
        self.context, self.coefficients = context, coefficients

    def __mul__(self, other):
        other = self._coerce(other)
        order = min(self.order, other.order)
        coefficients = []
        for index in range(order + 1):
            value = self.context.mpf(0)
            for left in range(index + 1):
                value += self[left] * other[index - left]
            coefficients.append(value)
        return IntervalTaylorJet(self.context, coefficients)

    def differentiate(self):
        return IntervalTaylorJet(self.context, [(index + 1) * self.coefficients[index + 1] for index in range(self.order)])

    def reciprocal(self):
        _require_nonzero(self[0], "Taylor denominator")
        coefficients = [self.context.mpf(1) / self[0]]
        for index in range(1, self.order + 1):
            value = self.context.mpf(0)
            for left in range(1, index + 1):
                value += self[left] * coefficients[index - left]
            coefficients.append(-value / self[0])
        return IntervalTaylorJet(self.context, coefficients)


def _source_jet(context, grid, radial_order, axial_order, label):
    coefficients = []
    for order in range(axial_order + 1):
        value = _source_value(grid, radial_order, order, label)
        coefficients.append(value / context.mpf(math.factorial(order)))
    return IntervalTaylorJet(context, coefficients)


def _relative_f0_jet(context, z, delta, Lambda, required_j, order):
    one = IntervalTaylorJet.constant(context, 1, order)
    H = z * (one - delta) / 3
    sigma = required_j / context.mpf(700)
    denominator = _square(H) + sigma ** 2
    _require_nonzero(denominator[0], "H^2+(j/500)^2")
    ell = -(Lambda * H) / denominator
    coefficients = [context.mpf(1)]
    for n in range(1, order + 1):
        value = context.mpf(0)
        for k in range(n):
            value += ell[k] * coefficients[n - 1 - k]
        coefficients.append(value / n)
    return IntervalTaylorJet(context, coefficients)
'''


SYNTHETIC_SOURCE_LOW = r'''
class Dual:
    def __init__(self, context, value, derivative=0):
        self.context, self.value, self.derivative = context, value, derivative

    def __mul__(self, other):
        other = self._coerce(other)
        return Dual(
            self.context,
            self.value * other.value,
            self.derivative * other.value + self.value * other.derivative,
        )

    def reciprocal(self):
        _require_nonzero(self.value, "dual denominator")
        inverse = self.context.mpf(1) / self.value
        return Dual(self.context, inverse, -self.derivative / (self.value ** 2))


def build(c, rho, z, delta, Lambda, j, Phi_grid, U_grid, M_grid):
    ell = Dual(c, ell0, ell1)
    relative_f0 = Dual(c, c.mpf(1), ell.value)
    phi0 = Dual(c, phi0_0, phi0_1)
    f = phi0 * relative_f0
    fz = Dual(
        c,
        f.derivative,
        phi0_2 + 2 * ell.value * phi0_1 + (ell.derivative + ell.value ** 2) * phi0_0,
    )
    return fz
'''


SOURCE_HIGH = SYNTHETIC_SOURCE_HIGH
SOURCE_LOW = SYNTHETIC_SOURCE_LOW
BASE_INPUTS = (SOURCE_HIGH, SOURCE_LOW)

# Mutations are made against unique literal fragments so that each failure is
# attributable to one static anchor rather than to a large source diff.
SOURCE_HIGH_BAD_SIGNATURE = SOURCE_HIGH.replace(
    "def _source_jet(context, grid, radial_order, axial_order, label):",
    "def _source_jet(context, grid, radial_order, axial_order, label, extra):",
)
SOURCE_HIGH_BAD_DENOMINATOR = SOURCE_HIGH.replace("sigma ** 2", "sigma ** 3")
SOURCE_HIGH_BAD_FACTOR = SOURCE_HIGH.replace("value / n", "value / (n + 1)")
SOURCE_HIGH_OMITTED_SOURCE_JET = SOURCE_HIGH.replace("def _source_jet", "def source_jet")
SOURCE_LOW_BAD_SIGNATURE = SOURCE_LOW.replace(
    "def build(c, rho, z, delta, Lambda, j, Phi_grid, U_grid, M_grid):",
    "def build(c, rho, z, delta, Lambda, j, Phi_grid, U_grid, M_grid, extra):",
)
SOURCE_LOW_BAD_DENOMINATOR = SOURCE_LOW.replace("self.value ** 2", "self.value ** 3")
SOURCE_LOW_BAD_FACTOR = SOURCE_LOW.replace("phi0_2 + 2 * ell.value", "phi0_2 + 3 * ell.value")

BAD_SIGNATURE_INPUTS = (SOURCE_HIGH_BAD_SIGNATURE, SOURCE_LOW_BAD_SIGNATURE)
BAD_DENOMINATOR_INPUTS = (SOURCE_HIGH_BAD_DENOMINATOR, SOURCE_LOW_BAD_DENOMINATOR)
BAD_FACTOR_INPUTS = (SOURCE_HIGH_BAD_FACTOR, SOURCE_LOW_BAD_FACTOR)
OMISSION_INPUTS = (SOURCE_HIGH_OMITTED_SOURCE_JET, SOURCE_LOW)


def run_fixture(leaf=None):
    """Run source-string mutation checks through the new proof leaf."""

    if leaf is None:
        from . import core_jet_projection_v2 as leaf

    cases = {
        "base": leaf.build(*BASE_INPUTS),
        "bad_signature": leaf.build(*BAD_SIGNATURE_INPUTS),
        "bad_denominator": leaf.build(*BAD_DENOMINATOR_INPUTS),
        "bad_factor": leaf.build(*BAD_FACTOR_INPUTS),
        "omission": leaf.build(*OMISSION_INPUTS),
        "syntax": leaf.build("def invalid(:", SOURCE_LOW),
    }
    checks = {
        "base_proved": cases["base"].get("status") == "proved" and cases["base"].get("passed") is True,
        "default_signature_anchored": cases["base"].get("anchors", {}).get("low_fz", {}).get("matched") is True,
        "exact_theorem_present": cases["base"].get("theorem", {}).get("exact") is True,
        "bad_signature_rejected": cases["bad_signature"].get("passed") is False and bool(cases["bad_signature"].get("differences")),
        "bad_denominator_rejected": cases["bad_denominator"].get("passed") is False and bool(cases["bad_denominator"].get("differences")),
        "bad_factor_rejected": cases["bad_factor"].get("passed") is False and bool(cases["bad_factor"].get("differences")),
        "omission_rejected": cases["omission"].get("passed") is False and bool(cases["omission"].get("differences")),
        "syntax_rejected": cases["syntax"].get("passed") is False and bool(cases["syntax"].get("differences")),
        "no_physical_sum_claim": cases["base"].get("physical_field") is False and cases["base"].get("sum_claim") is False,
    }
    return {
        "schema": "core_jet_projection_fixture_v2",
        "synthetic_only": True,
        "passed": all(checks.values()),
        "check_count": len(checks),
        "checks": checks,
        "reports": cases,
    }


__all__ = [
    "BASE_INPUTS",
    "BAD_DENOMINATOR_INPUTS",
    "BAD_FACTOR_INPUTS",
    "BAD_SIGNATURE_INPUTS",
    "OMISSION_INPUTS",
    "SOURCE_HIGH",
    "SOURCE_HIGH_BAD_DENOMINATOR",
    "SOURCE_HIGH_BAD_FACTOR",
    "SOURCE_HIGH_BAD_SIGNATURE",
    "SOURCE_HIGH_OMITTED_SOURCE_JET",
    "SOURCE_LOW",
    "SOURCE_LOW_BAD_DENOMINATOR",
    "SOURCE_LOW_BAD_FACTOR",
    "SOURCE_LOW_BAD_SIGNATURE",
    "SYNTHETIC_SOURCE_HIGH",
    "SYNTHETIC_SOURCE_LOW",
    "run_fixture",
]
