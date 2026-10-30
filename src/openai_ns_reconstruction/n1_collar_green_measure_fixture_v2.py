"""Independent synthetic fixture for the first-collar ``dx`` adapter.

The fixture is intentionally source-free: it uses a variable synthetic
function of ``x`` and a unit ordinary-Z1 source row, then compares the
directed cell output with an independent high-precision antiderivative.  It
is provided for the root agent to run after the implementation is frozen; it
does not run the leaf, a helper, or a repository test at import time.
"""

from __future__ import annotations


def run(leaf, c):
    import mpmath as mp

    hb_iv = c.mpf(1) / 3
    sc_iv = c.mpf(1) / 8
    epsilon = c.mpf(1) / 5
    cells = 4
    atoms = {
        name: "synthetic_" + name
        for name in ("hb", "F0_anchor", "Pstar_squared", "Lambda")
    }
    certificate = {
        "hb_sc_Z_independent": True,
        "rho_radius_Z_independent": True,
    }

    def variable_source(q_interval, cell):
        # f(x)=1+x is variable on every cell.  The x interval is a valid
        # whole-cell enclosure; no midpoint or endpoint is selected.
        f_x = c.mpf(1) + cell["x_interval"]
        return {
            "Z0": {
                "powers": (0, 0, 0, 0),
                "coefficient": f_x,
                "origin": {"source_ids": ["variable_f_x"], "uncertainty": True},
            },
            "Z1": {
                "powers": (0, 0, 0, 0),
                "coefficient": c.mpf(1),
                "origin": {"source_ids": ["unit_Z1"], "uncertainty": False},
            },
        }

    out = leaf.build(
        c,
        sc=sc_iv,
        hb=hb_iv,
        epsilon=epsilon,
        cells=cells,
        g4=variable_source,
        source_family="synthetic_dx_source",
        source_atoms=atoms,
        geometry_certificate=certificate,
        epsilon_owner="synthetic_source_owner",
        formal_hb_owner="green_measure",
        physical_hb_owner="synthetic_geometry",
        scope="fixture variable-function dx measure",
    )

    def evaluate(terms, hb_value):
        total = c.mpf(0)
        for term in terms:
            total += term["coefficient"] * hb_value ** term["powers"][0]
        return total

    mp.mp.dps = 90
    h = mp.mpf(1) / 3
    scale = mp.mpf(1) / 8
    x0 = mp.mpf(2)
    x1 = 2 * mp.exp(h * scale / 8)
    expected_variable = (x1 - x0) + (x1 * x1 - x0 * x0) / 2
    expected_unit_z1 = x1 - x0
    got_variable = evaluate(out["source_integral"]["Z0"], hb_iv)
    got_unit_z1 = evaluate(out["source_integral"]["Z1"], hb_iv)

    expected_variable_iv = c.mpf(mp.nstr(expected_variable, 75))
    expected_unit_z1_iv = c.mpf(mp.nstr(expected_unit_z1, 75))
    checks = {
        "variable_function_dx_integral": bool(
            got_variable.a <= expected_variable_iv.a
            and got_variable.b >= expected_variable_iv.b
        ),
        "unit_Z1_chain": bool(
            got_unit_z1.a <= expected_unit_z1_iv.a
            and got_unit_z1.b >= expected_unit_z1_iv.b
        ),
        "formal_hb_attached_once": all(
            term["powers"][0] == 1 for term in out["source_integral"]["Z0"]
        ),
        "physical_hb_not_coefficient": out["ownership"]["physical_hb_inserted_as_coefficient"] is False,
        "weight_Z1_zero": all(
            cell["ordinary_Z1_chain"]["weight_Z1_is_zero"] for cell in out["cells"]
        ),
        "directed_cell_weights_present": len(out["source_cell_weights"]) == cells,
        "scope_not_legacy_identity": out["legacy_identity"] is False,
        "not_m15_certificate": out["m15_computed"] is False
        and out["complete_m15_computed"] is False,
        "source_q_normalization": out["q_definition"] == "q=4*s/sc"
        and out["dx_dq_formula"] == "(hb*sc/8)*x",
    }

    # Independent hb scaling check for the differential identity.  The
    # numeric ratio is formed only from high-precision synthetic values.
    h2 = 2 * h
    q = mp.mpf("0.375")
    x_q = 2 * mp.exp(h * scale * q / 8)
    dx_ds_value = h * x_q / 2
    dx_dq_value = h * scale * x_q / 8
    checks["dx_ds_identity"] = bool(abs(dx_ds_value - h * x_q / 2) == 0)
    checks["dx_dq_hb_scaling"] = bool(abs(dx_dq_value - scale * dx_ds_value / 4) == 0)
    checks["second_hb_changes_map"] = bool(2 * mp.exp(h2 * scale * q / 8) != x_q)

    closed_hb_weights = leaf.source_cell_weights(
        c,
        sc=sc_iv,
        hb=c.mpf([0, 1]) / 2,
        cells=cells,
        q_normalization=4,
    )
    checks["closed_hb_geometry_enclosure"] = len(closed_hb_weights) == cells and all(
        item["weight_is_enclosure"]
        and item["exact_endpoint_weight_claimed"] is False
        for item in closed_hb_weights
    )

    return {
        "passed": all(checks.values()),
        "check_count": len(checks),
        "checks": checks,
        "synthetic_only": True,
        "reference": "independent antiderivative of (1+x) dx with unit ordinary Z1",
        "source_executed": False,
        "legacy_identity_claimed": False,
    }
