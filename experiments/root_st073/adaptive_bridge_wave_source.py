"""Build a local Kelvin/stress source for the three-knot moment-closed bridge."""

import json

from adaptive_bridge_recursive_defect import build_fields
from radial_peak_cone import run as cone_run
from radial_continuation import ROOT
from separated_moment_modes import SeparatedMomentModes, RADIAL_WINDOWS_THREE


def run(k=11.0, eta=-0.2, y=0.05,
        output_name="adaptive_bridge_wave_source.json"):
    inner, fields = build_fields()
    report = json.loads((ROOT / "separated_moment_three_knots.json").read_text())
    field = SeparatedMomentModes(
        fields["two_sided_cone"], report["constrained_l4"]["amplitudes"],
        windows=RADIAL_WINDOWS_THREE, knots=(11.0, 15.0, 19.0))
    tau = 0.5 * 2.0**-k
    X = inner.p.X_max * (1.0 + 15.0 * y)**2
    point = inner.from_similarity([X], [eta], tau)[0]
    cone_run(radius=float(point[0]), z=float(point[2]),
             tau=tau, field=field, field_id="three_knot_moment_closed",
             output_name=output_name)


if __name__ == "__main__":
    run()
