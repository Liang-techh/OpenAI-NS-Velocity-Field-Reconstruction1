from decimal import Decimal, localcontext
from fractions import Fraction
import math
from dataclasses import dataclass, replace
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from openai_ns_reconstruction.axis_coefficient_data import (
    actual_schedule_axis_coefficient_data,
)
from openai_ns_reconstruction.axis_coefficient_operators import (
    actual_schedule_coefficient_operators,
)
from openai_ns_reconstruction.axis_coefficient_reference_state import (
    AxisCoefficientJetState,
)
from openai_ns_reconstruction.axis_coefficient_wide_first_picard import (
    actual_schedule_wide_first_picard_state,
)
from openai_ns_reconstruction.axis_coefficient_wide_first_picard_slow2_axial_quadratic import (
    MixedScaleFirstPicardSlow2AxialQuadraticCoefficientJet,
)
from openai_ns_reconstruction.axis_coefficient_wide_first_picard_slow2_average_dot import (
    MixedScaleFirstPicardSlow2AverageDotCoefficientJet,
)
from openai_ns_reconstruction.axis_coefficient_wide_first_picard_slow2_average_mixed import (
    MixedScaleFirstPicardSlow2AverageMixedCoefficientJet,
)
from openai_ns_reconstruction.axis_coefficient_wide_first_picard_slow2_param import (
    MixedScaleFirstPicardSlow2ParamCoefficientJet,
)
from openai_ns_reconstruction.axis_coefficient_wide_first_picard_slow2 import (
    ActualScheduleWideFirstPicardSlow2State,
    MixedScaleFirstPicardSlow2CoefficientJet,
    wide_first_picard_slow2_state,
)
from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData


PRECISION = 96
FIELDS = (
    "ordinary_reference",
    "ordinary_inverse_lambda_numerator",
    "ordinary_inverse_lambda_squared_numerator",
    "ordinary_inverse_lambda_cubed_numerator",
    "ordinary_inverse_lambda_fourth_numerator",
    "pressure_linear_inverse_lambda_numerator",
    "pressure_linear_inverse_lambda_squared_numerator",
    "pressure_linear_inverse_lambda_cubed_numerator",
    "pressure_square_inverse_lambda_squared_numerator",
)

# The four input channels of the first-Picard axial state are
#
#   u = u_0 + u_1/Lambda + u_2/Lambda^2 + a^2 p_1/Lambda.
#
# The oracle below keeps those channels separate and forms the nine output
# channels by direct radial convolution.  It never calls any of the four
# slow2 branch modules.  The pressure-normalized jets already include the
# Bell factors from the actual, eta-dependent amplitude; the test also checks
# that this amplitude has a nonzero eta derivative.
INPUT_CHANNELS = 4
OUTPUT_CHANNELS = 9


def _schedule_data() -> TailData:
    return TailData(
        OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0),
        h=0.01,
    )


@pytest.fixture(scope="module")
def x1():
    return actual_schedule_wide_first_picard_state(_schedule_data(), 0.05)


@pytest.fixture(scope="module")
def complete_state(x1):
    return wide_first_picard_slow2_state(x1)


def _zero(size: int) -> list[Decimal]:
    return [Decimal(0) for _ in range(size)]


def _eta_decimal(value) -> Decimal:
    return value if isinstance(value, Decimal) else Decimal.from_float(float(value))


def _add_in_place(target: list[Decimal], source: list[Decimal], weight=Decimal(1)) -> None:
    for index, value in enumerate(source):
        target[index] += weight * value


def _u_channels(x1, n: int, m: int, eta: float) -> tuple[Decimal, ...]:
    """Return the four theorem-scale input channels at one coefficient jet."""

    _, axial = x1.jet_pair(n, m, eta)
    with localcontext() as ctx:
        ctx.prec = PRECISION
        return (
            +axial.reference,
            +axial.inverse_lambda_numerator,
            +axial.inverse_lambda_squared_numerator,
            +(
                x1.remainder.axial.pressure_normalized_factor(n, m, eta)
                / Decimal(2)
            ),
        )


def _bilinear(left: tuple[Decimal, ...], right: tuple[Decimal, ...]) -> list[Decimal]:
    """Multiply two four-channel split jets, collecting all nine channels."""

    assert len(left) == INPUT_CHANNELS
    assert len(right) == INPUT_CHANNELS
    result = _zero(OUTPUT_CHANNELS)
    with localcontext() as ctx:
        ctx.prec = PRECISION
        # Ordinary times ordinary gives Lambda^-(i+j).
        for left_power in range(3):
            for right_power in range(3):
                result[left_power + right_power] += (
                    left[left_power] * right[right_power]
                )
        # Exactly one pressure channel gives a^2 Lambda^-(power+1).
        for ordinary_power in range(3):
            result[5 + ordinary_power] += left[ordinary_power] * right[3]
            result[5 + ordinary_power] += left[3] * right[ordinary_power]
        # Two pressure channels give a^4 Lambda^-2.
        result[8] += left[3] * right[3]
        return [+value for value in result]


def _sum_output(*values: list[Decimal]) -> list[Decimal]:
    result = _zero(OUTPUT_CHANNELS)
    with localcontext() as ctx:
        ctx.prec = PRECISION
        for value in values:
            _add_in_place(result, value)
        return [+value for value in result]


def _u_family(x1):
    return lambda n, m, eta: _u_channels(x1, n, m, eta)


def _average_input(family):
    def averaged(n: int, m: int, eta: float) -> tuple[Decimal, ...]:
        with localcontext() as ctx:
            ctx.prec = PRECISION
            return tuple(value / Decimal(n + 1) for value in family(n, m, eta))

    return averaged


def _eta_scaled_input(family, scalar: Decimal):
    def scaled(n: int, m: int, eta: float) -> tuple[Decimal, ...]:
        with localcontext() as ctx:
            ctx.prec = PRECISION
            value = [_eta_decimal(eta) * item for item in family(n, m, eta)]
            if m:
                _add_in_place(value, list(family(n, m - 1, eta)), Decimal(m))
            return tuple(scalar * item for item in value)

    return scaled


def _eta_scaled_output(family, scalar: Decimal):
    def scaled(n: int, m: int, eta: float) -> list[Decimal]:
        with localcontext() as ctx:
            ctx.prec = PRECISION
            result = [_eta_decimal(eta) * item for item in family(n, m, eta)]
            if m:
                _add_in_place(result, family(n, m - 1, eta), Decimal(m))
            return [scalar * item for item in result]

    return scaled


def _d_times_input(family):
    """Multiply an input family by d(eta)=1-eta^2 using exact jets."""

    def transformed(n: int, m: int, eta: float) -> tuple[Decimal, ...]:
        result = _zero(INPUT_CHANNELS)
        with localcontext() as ctx:
            ctx.prec = PRECISION
            eta_decimal = _eta_decimal(eta)
            d_jets = (
                Decimal(1) - eta_decimal * eta_decimal,
                -Decimal(2) * eta_decimal,
                Decimal(-2),
            )
            for order, d_jet in enumerate(d_jets):
                if order <= m:
                    _add_in_place(
                        result,
                        list(family(n, m - order, eta)),
                        Decimal(math.comb(m, order)) * d_jet,
                    )
            return tuple(+value for value in result)

    return transformed


def _d_times_output(family):
    """Multiply a nine-channel family by d(eta)=1-eta^2."""

    def transformed(n: int, m: int, eta: float) -> list[Decimal]:
        result = _zero(OUTPUT_CHANNELS)
        with localcontext() as ctx:
            ctx.prec = PRECISION
            eta_decimal = _eta_decimal(eta)
            d_jets = (
                Decimal(1) - eta_decimal * eta_decimal,
                -Decimal(2) * eta_decimal,
                Decimal(-2),
            )
            for order, d_jet in enumerate(d_jets):
                if order <= m:
                    _add_in_place(
                        result,
                        family(n, m - order, eta),
                        Decimal(math.comb(m, order)) * d_jet,
                    )
            return [+value for value in result]

    return transformed


def _j1(source):
    def inverse(n: int, m: int, eta: float) -> list[Decimal]:
        if n == 0:
            # Keep the source validation live on the exact zero row.
            source(0, m, eta)
            return _zero(OUTPUT_CHANNELS)
        with localcontext() as ctx:
            ctx.prec = PRECISION
            return [+(value / Decimal(n * n)) for value in source(n - 1, m, eta)]

    return inverse


def _radial_product(left, right, *, right_euler: bool, left_eta_shift: bool):
    """Independent coefficient convolution for dot, mixed, and param terms."""

    def product(n: int, m: int, eta: float) -> list[Decimal]:
        if n == 0:
            return _zero(OUTPUT_CHANNELS)
        result = _zero(OUTPUT_CHANNELS)
        with localcontext() as ctx:
            ctx.prec = PRECISION
            divisor = Decimal(n * n)
            for i in range(n):
                j = n - 1 - i
                right_weight = Decimal(j) if right_euler else Decimal(1)
                for k in range(m + 1):
                    left_order = k + 1 if left_eta_shift else k
                    term = _bilinear(
                        left(i, left_order, eta),
                        right(j, m - k, eta),
                    )
                    _add_in_place(
                        result,
                        term,
                        Decimal(math.comb(m, k))
                        * right_weight
                        / divisor,
                    )
            return [+value for value in result]

    return product


def _slow2_from_family(
    U,
    A: Decimal,
    D: Decimal,
    n: int,
    m: int,
    eta: float,
) -> tuple[Decimal, ...]:
    """Evaluate the four-term slow2 formula from a split input family.

    This is deliberately a test-only differential/convolution implementation;
    it does not import or call any of the four production slow2 branches.
    """

    # J1[2 A eta * u^2].
    square = _j1(
        _eta_scaled_output(
            lambda i, q, z: _product_input(U, U, i, q, z),
            Decimal(2) * A,
        )
    )

    # J1[(2 D eta) average(u) * (Y d_Y u)].
    averaged = _average_input(U)
    dot_left = _eta_scaled_input(averaged, Decimal(2) * D)
    dot = _radial_product(dot_left, U, right_euler=True, left_eta_shift=False)

    # d * J1[(d_eta average(u)) * (Y d_Y u)].
    mixed = _radial_product(
        averaged,
        U,
        right_euler=True,
        left_eta_shift=True,
    )
    mixed = _d_times_output(mixed)

    # -J1[(d_eta u) * d*u].
    param_source = _radial_product(
        U,
        _d_times_input(U),
        right_euler=False,
        left_eta_shift=True,
    )

    def param(n: int, m: int, eta: float) -> list[Decimal]:
        return [-value for value in param_source(n, m, eta)]

    # The square branch currently uses a helper whose input product is written
    # below so the oracle's channel classification stays visible here.
    return tuple(
        _sum_output(
            square(n, m, eta),
            dot(n, m, eta),
            mixed(n, m, eta),
            param(n, m, eta),
        )
    )


def _independent_slow2_oracle(x1, n: int, m: int, eta: float) -> tuple[Decimal, ...]:
    data = actual_schedule_axis_coefficient_data(x1.reference)
    return _slow2_from_family(
        _u_family(x1),
        Decimal.from_float(data.A),
        Decimal.from_float(data.D),
        n,
        m,
        eta,
    )


@dataclass(frozen=True)
class _SyntheticScaleFixture:
    """Synthetic finite polynomial input; never treated as paper evidence."""

    # Each ordinary entry is (radial degree -> eta polynomial coefficients).
    ordinary: tuple[tuple[tuple[Decimal, ...], ...], ...]
    pressure: tuple[tuple[Decimal, ...], ...]
    log_amplitude_slope: Decimal

    def family(self):
        def evaluate(n: int, m: int, eta: float) -> tuple[Decimal, ...]:
            eta_decimal = Decimal(str(eta))
            if n >= len(self.pressure):
                return (Decimal(0),) * INPUT_CHANNELS

            def derivative(coefficients: tuple[Decimal, ...], order: int) -> Decimal:
                return sum(
                    Decimal(math.factorial(power))
                    / Decimal(math.factorial(power - order))
                    * coefficient
                    * eta_decimal ** (power - order)
                    for power, coefficient in enumerate(coefficients)
                    if power >= order
                )

            ordinary = [
                derivative(channel[n], m) if n < len(channel) else Decimal(0)
                for channel in self.ordinary
            ]
            # The pressure channel is the normalized derivative of a^2 p.  For
            # a(eta)=exp(k eta), its Bell factor is (2 k)^q.
            pressure = sum(
                Decimal(math.comb(m, q))
                * (Decimal(2) * self.log_amplitude_slope) ** q
                * derivative(self.pressure[n], m - q)
                for q in range(m + 1)
            )
            return tuple(ordinary + [pressure])

        return evaluate


_SYNTHETIC_FIXTURE = _SyntheticScaleFixture(
    ordinary=(
        (
            (Decimal(1), Decimal(2)),
            (Decimal(2), Decimal(-1), Decimal(1)),
            (Decimal(3), Decimal(1)),
            (Decimal(-1), Decimal(2), Decimal(1)),
        ),
        (
            (Decimal(2), Decimal(-1)),
            (Decimal(1), Decimal(1), Decimal(1)),
            (Decimal(-2), Decimal(2)),
            (Decimal(1), Decimal(-2), Decimal(1)),
        ),
        (
            (Decimal(-1), Decimal(1)),
            (Decimal(3), Decimal(1)),
            (Decimal(1), Decimal(-1), Decimal(1)),
            (Decimal(2), Decimal(1)),
        ),
    ),
    pressure=(
        (Decimal(1), Decimal(-1), Decimal(1)),
        (Decimal(2), Decimal(1)),
        (Decimal(-1), Decimal(2), Decimal(1)),
        (Decimal(1), Decimal(1), Decimal(-1)),
    ),
    log_amplitude_slope=Decimal("1.5"),
)


def _product_input(left, right, n: int, m: int, eta: float) -> list[Decimal]:
    result = _zero(OUTPUT_CHANNELS)
    with localcontext() as ctx:
        ctx.prec = PRECISION
        for i in range(n + 1):
            j = n - i
            for k in range(m + 1):
                _add_in_place(
                    result,
                    _bilinear(left(i, k, eta), right(j, m - k, eta)),
                    Decimal(math.comb(m, k)),
                )
        return [+value for value in result]


_INJECTED_BRANCH_ROWS = {
    "axial_quadratic": tuple(
        Fraction(value)
        for value in (
            "25/3",
            "6/5",
            "109/9",
            "8/3",
            "217/45",
            "34976/1125",
            "8521/375",
            "224639/5625",
            "178741/5625",
        )
    ),
    "average_dot": tuple(
        Fraction(value)
        for value in (
            "157/30",
            "-1/25",
            "959/150",
            "-46/75",
            "269/150",
            "24389/1875",
            "9429/1250",
            "111214/9375",
            "60469/6250",
        )
    ),
    "average_mixed": tuple(
        Fraction(value)
        for value in (
            "-776/225",
            "487/225",
            "-37/25",
            "-119/225",
            "-6/25",
            "43139/1875",
            "61277/9375",
            "768967/46875",
            "5477024/140625",
        )
    ),
    # This row is already the signed -param1 contribution.  The aggregate
    # must add it as supplied rather than negating it a second time.
    "param": tuple(
        Fraction(value)
        for value in (
            "194/75",
            "32/75",
            "32/225",
            "26/15",
            "46/225",
            "-918572/28125",
            "-85354/3125",
            "-3318488/140625",
            "-11445124/140625",
        )
    ),
}
_INJECTED_TOTALS = tuple(
    Fraction(value)
    for value in (
        "5717/450",
        "844/225",
        "103/6",
        "733/225",
        "329/50",
        "322916/9375",
        "11861/1250",
        "2090866/46875",
        "-55609/56250",
    )
)


def _fraction_decimal(value: Fraction) -> Decimal:
    with localcontext() as ctx:
        ctx.prec = PRECISION
        return +(Decimal(value.numerator) / Decimal(value.denominator))


def _injected_branch(
    x1,
    *,
    key: str,
    row: tuple[Fraction, ...],
    amplitude_delta: Decimal = Decimal(0),
):
    jet_types = {
        "axial_quadratic": MixedScaleFirstPicardSlow2AxialQuadraticCoefficientJet,
        "average_dot": MixedScaleFirstPicardSlow2AverageDotCoefficientJet,
        "average_mixed": MixedScaleFirstPicardSlow2AverageMixedCoefficientJet,
        "param": MixedScaleFirstPicardSlow2ParamCoefficientJet,
    }
    materialized_flags = {
        "axial_quadratic": "slow2_axial_quadratic_branch_materialized",
        "average_dot": "slow2_average_dot_branch_materialized",
        "average_mixed": "slow2_average_mixed_branch_materialized",
        "param": "slow2_param_branch_materialized",
    }

    def jet(_n: int, _m: int, eta: float):
        with localcontext() as ctx:
            ctx.prec = PRECISION
            amplitude_log_source = (
                x1.remainder.axial.wide_pressure.amplitude.log_amplitude_source(eta)
            )
            if amplitude_delta:
                amplitude_log_source = replace(
                    amplitude_log_source,
                    midpoint=+(amplitude_log_source.midpoint + amplitude_delta),
                )
            amplitude_log = +amplitude_log_source.midpoint
        values = {
            name: _fraction_decimal(value)
            for name, value in zip(FIELDS, row)
        }
        return jet_types[key](
            **values,
            Lambda=x1.Lambda,
            amplitude_log=amplitude_log,
            amplitude_log_source=amplitude_log_source,
        )

    return SimpleNamespace(
        x1=x1,
        epsilon=x1.epsilon,
        Lambda=x1.Lambda,
        **{materialized_flags[key]: True, "jet": jet},
    )


def _inject_branch_rows(state, x1, rows=None, *, malformed_key=None):
    rows = _INJECTED_BRANCH_ROWS if rows is None else rows
    for key, attribute in (
        ("axial_quadratic", "axial_quadratic"),
        ("average_dot", "average_dot"),
        ("average_mixed", "average_mixed"),
        ("param", "param"),
    ):
        object.__setattr__(
            state,
            attribute,
            _injected_branch(
                x1,
                key=key,
                row=rows[key],
                amplitude_delta=(
                    Decimal("1e700") if key == malformed_key else Decimal(0)
                ),
            ),
        )


def test_production_aggregate_sums_injected_nonzero_channels_and_preserves_scales(
    x1,
) -> None:
    """Exercise the production sum with independent nonzero branch jets.

    The injected rows are a finite polynomial regression fixture only.  They
    replace each branch's ``jet`` after the actual aggregate has been built,
    so this test checks the aggregate's production path without reusing any
    branch implementation to produce its expected values.
    """

    state = wide_first_picard_slow2_state(x1)
    _inject_branch_rows(state, x1)
    eta = 0.2
    actual = state.jet(3, 2, eta)

    with localcontext() as ctx:
        ctx.prec = PRECISION
        expected = tuple(_fraction_decimal(value) for value in _INJECTED_TOTALS)
        for name, value in zip(FIELDS, expected):
            assert abs(getattr(actual, name) - value) < Decimal("1e-90")

        amplitude_log = x1.remainder.axial.wide_pressure.amplitude.log_amplitude(eta)
        for name, term in zip(FIELDS[5:8], actual.pressure_linear_terms_log()):
            numerator = getattr(actual, name)
            assert numerator != 0
            assert term.sign == (1 if numerator > 0 else -1)
            assert term.log_scale == Decimal(2) * amplitude_log
            power = FIELDS[5:8].index(name) + 1
            assert term.log_factor == +(
                abs(numerator).ln() - Decimal(power) * actual.Lambda.ln()
            )

        square = actual.pressure_square_term_log()
        assert actual.pressure_square_inverse_lambda_squared_numerator < 0
        assert square.sign == -1
        assert square.log_scale == Decimal(4) * amplitude_log
        assert square.log_factor == +(
            abs(actual.pressure_square_inverse_lambda_squared_numerator).ln()
            - Decimal(2) * actual.Lambda.ln()
        )


def test_production_aggregate_cancellation_and_bad_amplitude_metadata_fail_closed(
    x1,
) -> None:
    """Check signed cancellation and per-jet amplitude identity guards."""

    cancelling_rows = {
        key: list(row) for key, row in _INJECTED_BRANCH_ROWS.items()
    }
    # These terminating values make the first ordinary and pressure-square
    # totals cancel exactly while keeping the parameter contribution negative.
    cancelling_rows["axial_quadratic"][0] = Fraction(1)
    cancelling_rows["average_dot"][0] = Fraction(2)
    cancelling_rows["average_mixed"][0] = Fraction(4)
    cancelling_rows["param"][0] = Fraction(-7)
    cancelling_rows["axial_quadratic"][8] = Fraction(1)
    cancelling_rows["average_dot"][8] = Fraction(2)
    cancelling_rows["average_mixed"][8] = Fraction(4)
    cancelling_rows["param"][8] = Fraction(-7)
    cancelling_rows = {
        key: tuple(row) for key, row in cancelling_rows.items()
    }

    state = wide_first_picard_slow2_state(x1)
    _inject_branch_rows(state, x1, cancelling_rows)
    actual = state.jet(3, 2, 0.2)
    assert actual.ordinary_reference == Decimal(0)
    assert actual.pressure_square_inverse_lambda_squared_numerator == Decimal(0)
    assert actual.pressure_square_term_log().sign == 0
    assert actual.pressure_square_term_log().log_scale is None
    assert actual.pressure_square_term_log().log_factor is None

    malformed = wide_first_picard_slow2_state(x1)
    _inject_branch_rows(malformed, x1, malformed_key="average_dot")
    with pytest.raises(ValueError, match="amplitude log mismatch"):
        malformed.jet(3, 2, 0.2)


def test_complete_state_binds_one_x1_and_preserves_truth_flags(complete_state, x1) -> None:
    assert isinstance(complete_state, ActualScheduleWideFirstPicardSlow2State)
    assert complete_state.x1 is x1
    assert complete_state.Lambda == x1.Lambda
    assert complete_state.epsilon == x1.epsilon
    assert complete_state.slow2_axial_quadratic_branch_materialized is True
    assert complete_state.slow2_average_dot_branch_materialized is True
    assert complete_state.slow2_average_mixed_branch_materialized is True
    assert complete_state.slow2_param_branch_materialized is True
    assert complete_state.slow2_all_constituent_branches_materialized is True
    assert complete_state.slow2_complete is True
    assert complete_state.natural_remainder_x1_materialized is False
    assert complete_state.fixed_point_materialized is False
    assert complete_state.fixed_point_convergence_certified is False
    assert complete_state.paper_exact is False
    assert all(branch.x1 is x1 for branch in (
        complete_state.axial_quadratic,
        complete_state.average_dot,
        complete_state.average_mixed,
        complete_state.param,
    ))


def test_ordinary_reference_matches_independent_pinned_operator_composition(
    x1,
    complete_state,
) -> None:
    """Cross-check the Lambda^0 channel through the landed operator record."""

    import math

    data = actual_schedule_axis_coefficient_data(x1.reference)
    operators = actual_schedule_coefficient_operators(x1.reference)
    U = x1.reference.u

    def scaled(state, factor):
        return AxisCoefficientJetState(
            epsilon=state.epsilon,
            origin="test-only scalar multiple",
            _jet_provider=lambda n, m, eta: factor * state.jet(n, m, eta),
        )

    two_A_eta = scaled(data.eta, 2.0 * data.A)
    two_D_eta = scaled(data.eta, 2.0 * data.D)
    quadratic = operators.j1(operators.product(two_A_eta, operators.product(U, U)))
    dot = operators.dot1(operators.product(two_D_eta, operators.average(U)), U)
    mixed = operators.product(
        data.d,
        operators.mixed1(operators.average(U), U),
    )
    param = operators.param1(U, operators.product(data.d, U))

    terms = (quadratic, dot, mixed, scaled(param, -1.0))
    expected = AxisCoefficientJetState(
        epsilon=x1.epsilon,
        origin="test-only ordinary slow2 operator sum",
        _jet_provider=lambda n, m, eta: math.fsum(
            term.jet(n, m, eta) for term in terms
        ),
    )
    expected_value = expected.jet(2, 1, -0.11)
    actual_value = complete_state.jet(2, 1, -0.11).ordinary_reference
    assert math.isclose(float(actual_value), expected_value, rel_tol=3e-13, abs_tol=3e-13)


def test_all_nine_channels_match_independent_differential_oracle(x1, complete_state) -> None:
    """Check every split channel from raw input jets, including pressure scales."""

    n, m, eta = 2, 1, 0.03
    actual = complete_state.jet(n, m, eta)
    expected = _independent_slow2_oracle(x1, n, m, eta)
    assert isinstance(actual, MixedScaleFirstPicardSlow2CoefficientJet)
    for name, value in zip(FIELDS, expected):
        assert getattr(actual, name) == value

    amplitude = x1.remainder.axial.wide_pressure.amplitude
    data = actual_schedule_axis_coefficient_data(x1.reference)
    amplitude_log = amplitude.log_amplitude(eta)
    with localcontext() as ctx:
        ctx.prec = PRECISION
        # The actual schedule has a nonconstant amplitude.  Its normalized
        # pressure jets, used by the oracle above, include this Bell factor.
        # The first two actual radial rows can have a zero pressure source;
        # the synthetic finite polynomial regression below supplies explicit
        # nonzero pressure-linear and pressure-square channels.
        log_gradient = amplitude.Lambda * Decimal.from_float(
            data.normalizedGradient.jet(0, 0, eta)
        )
        assert log_gradient != 0
        assert amplitude.log_amplitude(eta - 0.02) != amplitude_log
        assert amplitude.log_amplitude(eta + 0.02) != amplitude_log

        assert actual.amplitude_log == amplitude_log
        for power, term in enumerate(actual.pressure_linear_terms_log(), start=1):
            numerator = getattr(actual, FIELDS[4 + power], None)
            if numerator == 0:
                assert term.sign == 0
            else:
                assert term.log_scale == Decimal(2) * amplitude_log
        square_term = actual.pressure_square_term_log()
        if actual.pressure_square_inverse_lambda_squared_numerator == 0:
            assert square_term.sign == 0
        else:
            assert square_term.log_scale == Decimal(4) * amplitude_log


def test_synthetic_finite_polynomial_oracle_is_explicitly_noncertifying() -> None:
    """Exercise higher jets and nonconstant amplitude on a test-only fixture.

    This fixture is an analytic polynomial regression for the channel algebra;
    it is intentionally not a paper-data certificate and is not passed to the
    actual-schedule constructor.
    """

    eta = Decimal("0.2")
    result = _slow2_from_family(
        _SYNTHETIC_FIXTURE.family(),
        Decimal("1.25"),
        Decimal("0.75"),
        3,
        2,
        eta,
    )
    # These are the independently expanded finite-polynomial values at
    # (n,m,eta)=(3,2,1/5), retained here so this test does not merely compare
    # two calls to the production branch sum.  Convert exact rational values
    # under the same 96-digit context used by the oracle.
    expected_fractions = (
        Fraction(5717, 450),
        Fraction(844, 225),
        Fraction(103, 6),
        Fraction(733, 225),
        Fraction(329, 50),
        Fraction(322916, 9375),
        Fraction(11861, 1250),
        Fraction(2090866, 46875),
        Fraction(-55609, 56250),
    )
    with localcontext() as ctx:
        ctx.prec = PRECISION
        expected = tuple(
            +(Decimal(value.numerator) / Decimal(value.denominator))
            for value in expected_fractions
        )
    assert len(result) == len(FIELDS)
    assert all(value != 0 for value in result)
    for actual, reference in zip(result, expected):
        # The direct convolution rounds each intermediate at 96 digits;
        # compare against the exact rational expansion before that rounding.
        assert abs(actual - reference) < Decimal("1e-24")

    family = _SYNTHETIC_FIXTURE.family()
    pressure_jet_1 = family(0, 1, eta)[3]
    base_pressure_derivative = Decimal(-1) + Decimal(2) * Decimal(1) * Decimal(str(eta))
    assert pressure_jet_1 != base_pressure_derivative


def test_zero_row_is_exact_zero_but_keeps_scale_metadata(complete_state, x1) -> None:
    actual = complete_state.jet(0, 2, 0.07)
    for name in FIELDS:
        assert getattr(actual, name) == 0
    assert actual.Lambda == x1.Lambda
    assert actual.amplitude_log == (
        x1.remainder.axial.wide_pressure.amplitude.log_amplitude(0.07)
    )


def test_scale_helpers_keep_four_ordinary_and_three_pressure_powers(complete_state) -> None:
    actual = complete_state.jet(2, 1, 0.03)
    ordinary = actual.ordinary_correction_terms_decimal()
    assert len(ordinary) == 4
    assert all(isinstance(value, Decimal) for value in ordinary)

    linear = actual.pressure_linear_terms_log()
    assert len(linear) == 3
    for name, term in zip(FIELDS[5:8], linear):
        if getattr(actual, name) == 0:
            assert term.sign == 0
        else:
            assert term.log_scale == Decimal(2) * actual.amplitude_log
    square = actual.pressure_square_term_log()
    if actual.pressure_square_inverse_lambda_squared_numerator == 0:
        assert square.sign == 0
    else:
        assert square.log_scale == Decimal(4) * actual.amplitude_log


def test_guards_reject_bad_indices_eta_and_x1(complete_state, x1) -> None:
    with pytest.raises(ValueError, match="n must be a nonnegative integer"):
        complete_state.jet(-1, 0, 0.0)
    with pytest.raises(ValueError, match="m must be a nonnegative integer"):
        complete_state.jet(0, True, 0.0)
    with pytest.raises(ValueError, match="pinned window"):
        complete_state.jet(1, 0, 2.0)
    with pytest.raises(TypeError, match="x1 must be ActualScheduleWideFirstPicardState"):
        wide_first_picard_slow2_state(object())


def test_provenance_manifest_keeps_complete_slow2_boundary_fail_closed(x1) -> None:
    root = Path(__file__).resolve().parents[1]
    manifest = json.loads(
        (
            root
            / "references"
            / "provenance_manifest_addendum_axis_coefficient_wide_first_picard_slow2.json"
        ).read_text(encoding="utf-8")
    )
    assert manifest["full_reconstruction"] is False
    assert manifest["paper_exact_velocity_available"] is False
    layer = manifest["layer"]
    assert layer["id"] == "stage-1-leading-profile"
    assert layer["status"] == "formal-structure"
    truth = layer["truth_boundary"]
    assert truth["slow2_all_constituent_branches_materialized"] is True
    assert truth["slow2_complete"] is True
    assert truth["natural_remainder_x1_materialized"] is False
    assert truth["picard_x2_materialized"] is False
    assert truth["fixed_point_materialized"] is False
    assert (root / layer["provenance"]).is_file()
    for artifact in layer["artifacts"]:
        assert (root / artifact).is_file()

    mismatched = wide_first_picard_slow2_state(x1)
    object.__setattr__(
        mismatched,
        "average_dot",
        SimpleNamespace(
            x1=x1,
            epsilon=x1.epsilon + 1.0,
            Lambda=x1.Lambda,
            slow2_average_dot_branch_materialized=True,
        ),
    )
    with pytest.raises(ValueError, match="branch epsilon mismatch"):
        mismatched.jet(1, 0, 0.0)

    mismatched_lambda = wide_first_picard_slow2_state(x1)
    object.__setattr__(
        mismatched_lambda,
        "average_dot",
        SimpleNamespace(
            x1=x1,
            epsilon=x1.epsilon,
            Lambda=x1.Lambda + Decimal(1),
            slow2_average_dot_branch_materialized=True,
        ),
    )
    with pytest.raises(ValueError, match="branch Lambda mismatch"):
        mismatched_lambda.jet(1, 0, 0.0)

    mismatched_x1 = wide_first_picard_slow2_state(x1)
    object.__setattr__(
        mismatched_x1,
        "average_dot",
        SimpleNamespace(
            x1=object(),
            epsilon=x1.epsilon,
            Lambda=x1.Lambda,
            slow2_average_dot_branch_materialized=True,
        ),
    )
    with pytest.raises(ValueError, match="exact same x1 object"):
        mismatched_x1.jet(1, 0, 0.0)
