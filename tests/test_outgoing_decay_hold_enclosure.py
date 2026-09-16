from fractions import Fraction

import pytest

from openai_ns_reconstruction.outgoing_decay_hold_enclosure import (
    validated_decay_hold_enclosure,
)
from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData


def _actual_data() -> TailData:
    return TailData(
        OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0),
        h=0.01,
    )


def test_actual_decay_hold_interval_contains_legacy_value_and_is_positive() -> None:
    data = _actual_data()
    tolerance = Fraction(1, 64)

    result = validated_decay_hold_enclosure(data, absolute_tolerance=tolerance)

    assert result.h == Fraction.from_float(data.h)
    assert result.lam == Fraction.from_float(data.core.lam)
    assert result.width <= tolerance
    assert result.lag.lower > result.debt.upper > 0
    assert result.log_ratio.lower > 0
    assert result.hold.lower > 0
    assert result.lag_tolerance != result.debt_tolerance
    assert result.paper_exact is False
    assert result.global_axis_norm_certified is False
    assert result.full_reconstruction is False

    legacy_value = Fraction.from_float(data.decay_hold)
    assert result.hold.lower <= legacy_value <= result.hold.upper


def test_decay_hold_refinement_tightens_and_overlaps() -> None:
    data = _actual_data()
    loose = validated_decay_hold_enclosure(
        data,
        absolute_tolerance=Fraction(1, 32),
    )
    tight = validated_decay_hold_enclosure(
        data,
        absolute_tolerance=Fraction(1, 64),
    )

    assert loose.width <= Fraction(1, 32)
    assert tight.width <= Fraction(1, 64)
    assert tight.width < loose.width
    assert max(loose.hold.lower, tight.hold.lower) <= min(loose.hold.upper, tight.hold.upper)
    assert tight.lag_tolerance < loose.lag_tolerance
    assert tight.debt_tolerance <= loose.debt_tolerance
    assert tight.log_tolerance < loose.log_tolerance


def test_supplied_sources_must_match_same_actual_parameters() -> None:
    data = _actual_data()
    result = validated_decay_hold_enclosure(
        data,
        absolute_tolerance=Fraction(1, 32),
    )
    other = TailData(
        OutgoingCoreParameters(P=2.0, m=1.0, lam=0.04, wait=30.0),
        h=0.01,
    )

    with pytest.raises(ValueError, match="release lag parameters do not match TailData"):
        validated_decay_hold_enclosure(
            other,
            absolute_tolerance=Fraction(1, 32),
            release_lag=result.release_lag_source,
            tail_debt=result.tail_debt_source,
        )

    with pytest.raises(ValueError, match="must be supplied together"):
        validated_decay_hold_enclosure(
            data,
            release_lag=result.release_lag_source,
        )


def test_invalid_controls_and_computation_cap_fail_closed() -> None:
    data = _actual_data()

    with pytest.raises(TypeError, match="data must be TailData"):
        validated_decay_hold_enclosure(object())
    with pytest.raises(ValueError, match="absolute_tolerance must be positive"):
        validated_decay_hold_enclosure(data, absolute_tolerance=Fraction(0))
    with pytest.raises(ValueError, match="max_cells must be a positive integer"):
        validated_decay_hold_enclosure(data, max_cells=True)
    with pytest.raises(ValueError, match="max_terms must be a positive integer"):
        validated_decay_hold_enclosure(data, max_terms=0)
    with pytest.raises(ValueError, match="max_squarings must be a positive integer"):
        validated_decay_hold_enclosure(data, max_squarings=False)
    with pytest.raises(ValueError, match="max_refinements must be a positive integer"):
        validated_decay_hold_enclosure(data, max_refinements=0)
    with pytest.raises(ArithmeticError, match="max_cells"):
        validated_decay_hold_enclosure(
            data,
            absolute_tolerance=Fraction(1, 64),
            max_cells=1,
        )
