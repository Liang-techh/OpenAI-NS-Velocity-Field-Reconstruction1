"""Apply the same local cone map to the grid-aware mean correction."""

import json

from adaptive_bridge_cone_overlap import run
from radial_continuation import ROOT


if __name__ == "__main__":
    report = json.loads((ROOT/"adaptive_bridge_cone_grid_fit.json").read_text())
    selected = next(row for row in report["candidates"]
                    if row["cone_weight"] == report["selected_cone_weight"])
    run(coefficients=selected["amplitudes"],
        output_name="adaptive_bridge_cone_grid_overlap.json")
