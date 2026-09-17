import json

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

    fine = [row for row in report["rows"] if row["resolution"] == DEFAULT_RESOLUTIONS[-1]]
    raise AssertionError(
        json.dumps(
            {"time_summaries": report["time_summaries"], "fine_rows": fine},
            sort_keys=True,
        )
    )
