import json

import pytest

from openai_ns_reconstruction.constrained_eq45_supported_phi10_temporal_replay import (
    DEFAULT_SLOPES,
    replay_supported_phi10_temporal_slopes,
)


def test_phi10_temporal_replay_calibration():
    report = replay_supported_phi10_temporal_slopes()
    selected = report["selected_trial"]
    assert selected is not None
    summary = {
        "baseline_early": report["baseline_early_radial_outer_fraction"],
        "baseline_late": report["baseline_late_radial_outer_fraction"],
        "late_guard": report["late_guard_limit"],
        "trials": [
            {
                "slope": row["slope"],
                "early": row["early_radial_outer_fraction"],
                "early_relative_change": row["early_relative_change"],
                "late": row["late_radial_outer_fraction"],
                "coefficient_start": row["coefficient_start"],
                "coefficient_end": row["coefficient_end"],
                "feasible": row["feasible"],
            }
            for row in report["trials"]
        ],
        "selected_slope": selected["slope"],
    }
    pytest.fail("CALIBRATION " + json.dumps(summary, sort_keys=True))


def test_phi10_replay_rejects_invalid_slope_grid():
    with pytest.raises(ValueError, match="positive finite"):
        replay_supported_phi10_temporal_slopes(slopes=(0.8, 0.0))
    with pytest.raises(ValueError, match="unique"):
        replay_supported_phi10_temporal_slopes(slopes=(0.8, 0.8))


def test_default_slope_grid_stays_inside_screened_range():
    assert DEFAULT_SLOPES == (0.8, 1.1, 1.4, 1.7, 2.0)
