"""Fit cone geometry on both axial sides of the bridge hotspots."""

from adaptive_bridge_cone_grid_fit import run


if __name__ == "__main__":
    run(hotspot_indices=(0, 10),
        output_name="adaptive_bridge_cone_biside_fit.json")
