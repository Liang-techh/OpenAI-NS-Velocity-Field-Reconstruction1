"""Matched fixed-slope comparator using the feedback run's initial control."""
import json
import numpy as np
from outer_feedback_evolution import build_current, P_WINDOWS, P_BREAKS
from outer_swirl_slope import OuterSwirlSlope
from outer_pressure_modes import OuterPressure
from midplane_integrated_moment_balance import evaluate
from radial_continuation import ROOT


def run():
    report = json.loads((ROOT/'outer_feedback_evolution.json').read_text(encoding='utf-8'))
    inner, _, current = build_current()
    control = np.array(report['nodes'][0]['control'])
    field = OuterSwirlSlope(OuterPressure(current, control[:9], windows=P_WINDOWS),
                            control[9:], 11., P_WINDOWS)
    values = evaluate(field, inner, 11.0005, 96, .002, radial_breaks=P_BREAKS)
    result = dict(k=11.0005, order=96, moments=values.tolist(),
                  moment_max=float(max(abs(values))), initial_control=control.tolist(),
                  scope='Fixed slope comparator from exactly the feedback trajectory initial state and control; same integrated moment replay',
                  accepted=False)
    (ROOT/'outer_feedback_fixed_comparison.json').write_bytes(
        (json.dumps(result, indent=2)+'\n').encode())
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    run()
