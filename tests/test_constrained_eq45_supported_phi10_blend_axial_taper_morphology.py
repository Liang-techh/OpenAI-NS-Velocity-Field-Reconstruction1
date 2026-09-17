import json

from openai_ns_reconstruction.constrained_eq45_supported_phi10_blend_axial_taper_morphology import (
    audit_axial_taper_vorticity_morphology,
)


def test_axial_taper_vorticity_morphology_calibration_receipt():
    report = audit_axial_taper_vorticity_morphology()
    raise AssertionError(json.dumps(report["finest_summaries"], sort_keys=True))
