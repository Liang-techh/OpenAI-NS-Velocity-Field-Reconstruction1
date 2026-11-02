"""Independent source-text fixtures for ``core_operator_formula_identity_v1``.

These strings are intentionally self-contained and do not import or read the
real producer modules.  They contain the same sixteen matrix rows and seven
intermediate formula shapes under different names.  The two final high-side
fixtures change exactly one column expression or one intermediate expression;
the identity leaf should retain the evidence and report a non-proved result.
"""


SYNTHETIC_SOURCE_HIGH = r'''
def synthetic_high(context, rho_iv):
    degree = 5
    one = IntervalTaylorJet.constant(context, 1, degree)
    Q = (2 * z_j * U0 - (1 - delta_j) * z_j * M0 - d_j * M0_Z) / L_j
    H = z_j * (1 - delta_j) / 2 + d_j * w
    B = Q + (1 - 2 * z_j * w) / L_j
    E = z_j * w - context.mpf("0.5")
    Cf = f + rho_j * fr
    Zf = ((-2 - delta_j) * z_j * f + d_j * fz - 2 * z_j * rho_j * fr) / L_j
    Zw = ((-1 - delta_j) * z_j * w + d_j * wz - 2 * z_j * rho_j * wr) / L_j
    c0_jets = {
        "4,1": (1, -1, f * (4 * context.sqrt(rho_iv))),
        "5,1": (0, -1, 2 * (Q + E * (-2 + delta_j) / L_j)),
        "5,2": (1, -1, 2 * Zf - 2 * (-1 + delta_j) * z_j * Cf / L_j),
        "5,3": (1, -1, -2 * (1 + delta_j) * z_j * Cf / L_j),
        "5,5": (0, -1, B * context.sqrt(rho_iv)),
        "6,1": (1, -2, -8 * z_j * rho_j * f / L_j),
        "6,2": (0, -1, 2 * (E * (-1 + delta_j) / L_j + Zw - wr * z_j * (-1 + delta_j) * rho_j / L_j)),
        "6,3": (0, -1, -2 * (1 + delta_j) * z_j * rho_j * wr / L_j),
        "6,4": (0, -1, -4 * z_j / L_j),
        "6,6": (0, -1, B * context.sqrt(rho_iv)),
    }
    b1_jets = {
        "5,1": (0, -1, 2 * H / L_j),
        "5,2": (1, -1, -2 * d_j * Cf / L_j),
        "5,3": (1, -1, -2 * d_j * Cf / L_j),
        "6,2": (0, -1, 2 * (H - wr * d_j * rho_j) / L_j),
        "6,3": (0, -1, -2 * d_j * rho_j * wr / L_j),
        "6,4": (0, -1, 2 * d_j / L_j),
    }
    return c0_jets, b1_jets
'''


SYNTHETIC_SOURCE_LOW = r'''
def synthetic_low(c, rho_value):
    one = Dual.constant(c, 1)
    x = c.sqrt(rho_value)
    Q = (2 * z_dual * w - (one - delta_dual) * z_dual * m - d * mz) / L
    H = z_dual * (one - delta_dual) / 2 + d * w
    B = Q + (one - 2 * z_dual * w) / L
    E = z_dual * w - c.mpf("0.5")
    Cf = f + rho_dual * fr
    Zf = ((-2 - delta_dual) * z_dual * f + d * fz - 2 * z_dual * rho_dual * fr) / L
    Zw = ((-1 - delta_dual) * z_dual * w + d * wz - 2 * z_dual * rho_dual * wr) / L
    c0_jets = {
        "4,1": (1, -1, f * (4 * x)),
        "5,1": (0, -1, 2 * (Q + E * (-2 + delta_dual) / L)),
        "5,2": (1, -1, 2 * Zf - 2 * (-1 + delta_dual) * z_dual * Cf / L),
        "5,3": (1, -1, -2 * (1 + delta_dual) * z_dual * Cf / L),
        "5,5": (0, -1, B * x),
        "6,1": (1, -2, -8 * z_dual * rho_dual * f / L),
        "6,2": (0, -1, 2 * (E * (-1 + delta_dual) / L + Zw - wr * z_dual * (-1 + delta_dual) * rho_dual / L)),
        "6,3": (0, -1, -2 * (1 + delta_dual) * z_dual * rho_dual * wr / L),
        "6,4": (0, -1, -4 * z_dual / L),
        "6,6": (0, -1, B * x),
    }
    b1_jets = {
        "5,1": (0, -1, 2 * H / L),
        "5,2": (1, -1, -2 * d * Cf / L),
        "5,3": (1, -1, -2 * d * Cf / L),
        "6,2": (0, -1, 2 * (H - wr * d * rho_dual) / L),
        "6,3": (0, -1, -2 * d * rho_dual * wr / L),
        "6,4": (0, -1, 2 * d / L),
    }
    return c0_jets, b1_jets
'''


SYNTHETIC_SOURCE_HIGH_COLUMN_CHANGED = r'''
def synthetic_high_column_changed(context, rho_iv):
    degree = 5
    one = IntervalTaylorJet.constant(context, 1, degree)
    Q = (2 * z_j * U0 - (1 - delta_j) * z_j * M0 - d_j * M0_Z) / L_j
    H = z_j * (1 - delta_j) / 2 + d_j * w
    B = Q + (1 - 2 * z_j * w) / L_j
    E = z_j * w - context.mpf("0.5")
    Cf = f + rho_j * fr
    Zf = ((-2 - delta_j) * z_j * f + d_j * fz - 2 * z_j * rho_j * fr) / L_j
    Zw = ((-1 - delta_j) * z_j * w + d_j * wz - 2 * z_j * rho_j * wr) / L_j
    c0_jets = {
        "4,1": (1, -1, f * (4 * context.sqrt(rho_iv))),
        "5,1": (0, -1, 2 * (Q + E * (-2 + delta_j) / L_j)),
        "5,2": (1, -1, 2 * Zf - 2 * (-1 + delta_j) * z_j * Cf / L_j),
        "5,3": (1, -1, -2 * (1 + delta_j) * z_j * Cf / L_j),
        "5,5": (0, -1, B * context.sqrt(rho_iv)),
        "6,1": (1, -2, -8 * z_j * rho_j * f / L_j),
        "6,2": (0, -1, 2 * (E * (-1 + delta_j) / L_j + Zw - wr * z_j * (-1 + delta_j) * rho_j / L_j)),
        "6,3": (0, -1, -2 * (1 + delta_j) * z_j * rho_j * wr / L_j),
        "6,4": (0, -1, -4 * z_j / L_j),
        "6,6": (0, -1, B * (2 * context.sqrt(rho_iv))),
    }
    b1_jets = {
        "5,1": (0, -1, 2 * H / L_j),
        "5,2": (1, -1, -2 * d_j * Cf / L_j),
        "5,3": (1, -1, -2 * d_j * Cf / L_j),
        "6,2": (0, -1, 2 * (H - wr * d_j * rho_j) / L_j),
        "6,3": (0, -1, -2 * d_j * rho_j * wr / L_j),
        "6,4": (0, -1, 2 * d_j / L_j),
    }
    return c0_jets, b1_jets
'''


SYNTHETIC_SOURCE_HIGH_INTERMEDIATE_CHANGED = r'''
def synthetic_high_intermediate_changed(context, rho_iv):
    degree = 5
    one = IntervalTaylorJet.constant(context, 1, degree)
    Q = (2 * z_j * U0 + (1 - delta_j) * z_j * M0 - d_j * M0_Z) / L_j
    H = z_j * (1 - delta_j) / 2 + d_j * w
    B = Q + (1 - 2 * z_j * w) / L_j
    E = z_j * w - context.mpf("0.5")
    Cf = f + rho_j * fr
    Zf = ((-2 - delta_j) * z_j * f + d_j * fz - 2 * z_j * rho_j * fr) / L_j
    Zw = ((-1 - delta_j) * z_j * w + d_j * wz - 2 * z_j * rho_j * wr) / L_j
    c0_jets = {
        "4,1": (1, -1, f * (4 * context.sqrt(rho_iv))),
        "5,1": (0, -1, 2 * (Q + E * (-2 + delta_j) / L_j)),
        "5,2": (1, -1, 2 * Zf - 2 * (-1 + delta_j) * z_j * Cf / L_j),
        "5,3": (1, -1, -2 * (1 + delta_j) * z_j * Cf / L_j),
        "5,5": (0, -1, B * context.sqrt(rho_iv)),
        "6,1": (1, -2, -8 * z_j * rho_j * f / L_j),
        "6,2": (0, -1, 2 * (E * (-1 + delta_j) / L_j + Zw - wr * z_j * (-1 + delta_j) * rho_j / L_j)),
        "6,3": (0, -1, -2 * (1 + delta_j) * z_j * rho_j * wr / L_j),
        "6,4": (0, -1, -4 * z_j / L_j),
        "6,6": (0, -1, B * context.sqrt(rho_iv)),
    }
    b1_jets = {
        "5,1": (0, -1, 2 * H / L_j),
        "5,2": (1, -1, -2 * d_j * Cf / L_j),
        "5,3": (1, -1, -2 * d_j * Cf / L_j),
        "6,2": (0, -1, 2 * (H - wr * d_j * rho_j) / L_j),
        "6,3": (0, -1, -2 * d_j * rho_j * wr / L_j),
        "6,4": (0, -1, 2 * d_j / L_j),
    }
    return c0_jets, b1_jets
'''


# Short aliases make the intended runner inputs discoverable without reading
# any external source.  All four values above remain independent literals.
BASE_INPUTS = (SYNTHETIC_SOURCE_HIGH, SYNTHETIC_SOURCE_LOW)
COLUMN_CHANGED_INPUTS = (SYNTHETIC_SOURCE_HIGH_COLUMN_CHANGED, SYNTHETIC_SOURCE_LOW)
INTERMEDIATE_CHANGED_INPUTS = (
    SYNTHETIC_SOURCE_HIGH_INTERMEDIATE_CHANGED,
    SYNTHETIC_SOURCE_LOW,
)

SOURCE_HIGH = SYNTHETIC_SOURCE_HIGH
SOURCE_LOW = SYNTHETIC_SOURCE_LOW
SOURCE_HIGH_COLUMN_CHANGED = SYNTHETIC_SOURCE_HIGH_COLUMN_CHANGED
SOURCE_HIGH_INTERMEDIATE_CHANGED = SYNTHETIC_SOURCE_HIGH_INTERMEDIATE_CHANGED


def run_fixture(leaf):
    """Check the comparator on independent source strings and mutations."""
    import traceback

    try:
        base = leaf.build(*BASE_INPUTS)
        changed_column = leaf.build(*COLUMN_CHANGED_INPUTS)
        changed_intermediate = leaf.build(*INTERMEDIATE_CHANGED_INPUTS)
        broken = leaf.build("def invalid(:", SYNTHETIC_SOURCE_LOW)
        checks = {
            "base_proved": base.get("status") == "proved" and base.get("passed") is True,
            "all_16_terms": base.get("expected_matrix_item_count") == 16,
            "changed_column_detected": changed_column.get("passed") is False and bool(changed_column.get("differences")),
            "changed_intermediate_detected": changed_intermediate.get("passed") is False and bool(changed_intermediate.get("differences")),
            "invalid_source_detected": broken.get("passed") is False and bool(broken.get("differences")),
            "conditional_scope": base.get("conditional") is True and bool(base.get("unproved")),
        }
        return {"schema": "core_operator_formula_identity_fixture_v1", "passed": all(checks.values()), "check_count": len(checks), "checks": checks, "synthetic_only": True, "reports": {"base": base, "changed_column": changed_column, "changed_intermediate": changed_intermediate, "broken": broken}}
    except Exception as exc:
        return {"schema": "core_operator_formula_identity_fixture_v1", "passed": False, "error": str(exc), "traceback": traceback.format_exc(), "synthetic_only": True}

__all__ = [
    "BASE_INPUTS",
    "COLUMN_CHANGED_INPUTS",
    "INTERMEDIATE_CHANGED_INPUTS",
    "SOURCE_HIGH",
    "SOURCE_LOW",
    "SOURCE_HIGH_COLUMN_CHANGED",
    "SOURCE_HIGH_INTERMEDIATE_CHANGED",
    "SYNTHETIC_SOURCE_HIGH",
    "SYNTHETIC_SOURCE_LOW",
    "SYNTHETIC_SOURCE_HIGH_COLUMN_CHANGED",
    "SYNTHETIC_SOURCE_HIGH_INTERMEDIATE_CHANGED",
]
