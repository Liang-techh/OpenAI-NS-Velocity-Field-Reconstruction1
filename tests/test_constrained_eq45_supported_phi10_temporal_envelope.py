import json

import numpy as np

from openai_ns_reconstruction.constrained_eq45_supported_phi10_temporal_envelope import (
    DEFAULT_RESOLUTIONS,
    DEFAULT_SLOPE,
    DEFAULT_TIMES,
    audit_supported_phi10_temporal_envelope,
)


def test_supported_phi10_temporal_envelope_calibration():
    report = audit_supported_phi10_temporal_envelope()

    assert report["schema"] == "eq45_supported_phi10_temporal_envelope_audit_v1"
    assert report["slope"] == DEFAULT_SLOPE
    assert report["resolutions"] == list(DEFAULT_RESOLUTIONS)
    assert report["times"] == list(DEFAULT_TIMES)
    assert len(report["rows"]) == len(DEFAULT_RESOLUTIONS) * len(DEFAULT_TIMES)

    # Calibration-only sentinel: expose deterministic hosted summaries once,
    # then replace this assertion with fail-closed numerical regressions.
    raise AssertionError(json.dumps(report["time_summaries"], sort_keys=True))
