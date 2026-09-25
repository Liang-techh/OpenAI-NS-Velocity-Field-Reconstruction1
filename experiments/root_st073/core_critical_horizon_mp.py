"""Independent mpmath audit of the ST073-V critical-horizon core.

The production ``full_radial.py`` evaluator stores its coefficient jets in
NumPy (and uses the host ``longdouble`` implementation).  This file keeps the
same radial recurrence, but implements the axis jets, two-variable products,
derivatives, and cylindrical momentum algebra with ``mpmath`` scalars.  The
single requested point is ``X=1/64, eta=.3``; no floating-point coefficient
generation is used by the high-precision path.

This is a diagnostic of the unchanged local core.  It does not extend the
field's registered time domain into a global or matched solution.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import mpmath as mp


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "core_critical_horizon_mp.json"
ORDERS = (8, 10, 12, 14)
K_VALUES = (6, 14, 18, 22)
DPS_VALUES = (50, 70)

# ST073-V data are decimal/rational constants.  Module imports normally use
# mpmath's 15-digit default, so construct and retain the constants at higher
# precision before any per-case workdps context is entered.
with mp.workdps(100):
    H = mp.mpf("0.005")
    NU = mp.mpf("0.01")
    AXIS_SWIRL = mp.mpf("4")
    AXIAL_SLOPE = mp.mpf("4")
    AXIAL_BIAS = mp.mpf("0.02")
    PRESSURE = mp.mpf("1")
    X_POINT = mp.mpf(1) / 64
    ETA_POINT = mp.mpf(3) / 10
AXIS_TERMS = 120
JET_EXTRA = 3


def zeros(rows: int, cols: int) -> list[list[mp.mpf]]:
    return [[mp.mpf("0") for _ in range(cols)] for _ in range(rows)]


def binom_jet(beta: mp.mpf, n: int) -> list[mp.mpf]:
    out = [mp.mpf("1")]
    for i in range(1, n):
        out.append(out[-1] * (beta - i + 1) / i)
    return out


def q_power_jet(beta: mp.mpf, eta: mp.mpf, shape: tuple[int, int], terms: int = AXIS_TERMS):
    """Lagrange-inversion axis jet from full_radial.q_power_jet, in mp arithmetic."""
    rows, cols = shape
    beta = mp.mpf(beta)
    eps = 2 * H
    d = 1 - eta * eta
    out = zeros(rows, cols)
    for m in range(terms):
        if m == 0:
            cm = mp.mpf("1")
        else:
            cm = beta / m
            for j in range(1, m):
                cm *= (beta + eps * m - j) / j
        gamma = beta + (eps - 1) * m
        zlen = min(2 * m + 1, rows)
        z = [mp.mpf("0") for _ in range(zlen)]
        zbinom = binom_jet(mp.mpf(2 * m), zlen)
        for i in range(zlen):
            z[i] = zbinom[i] * eta ** (2 * m - i)
        tbinom = binom_jet(gamma, cols)
        for i in range(zlen):
            for j in range(cols):
                out[i][j] += cm * z[i] * tbinom[j] * d ** (gamma - j)
    return out


def add_inplace(dst, src, factor=mp.mpf("1")) -> None:
    for i in range(len(dst)):
        di, si = dst[i], src[i]
        for j in range(len(di)):
            di[j] += factor * si[j]


def derivative(a, axis: int, scale: mp.mpf):
    rows, cols = len(a), len(a[0])
    out = zeros(rows, cols)
    if axis == 0:
        for i in range(rows - 1):
            for j in range(cols):
                out[i][j] = (i + 1) * a[i + 1][j] / scale
    else:
        for i in range(rows):
            for j in range(cols - 1):
                out[i][j] = (j + 1) * a[i][j + 1] / scale
    return out


def multiply(a, b, shape: tuple[int, int]):
    """Truncated bivariate Cauchy product, with row/column order preserved."""
    rows, cols = shape
    out = zeros(rows, cols)
    for i in range(rows):
        imax = min(i, len(a) - 1)
        for j in range(cols):
            jmax = min(j, len(a[0]) - 1)
            total = mp.mpf("0")
            for ia in range(imax + 1):
                if ia >= len(b):
                    break
                ib = i - ia
                if ib >= len(b):
                    continue
                arow, brow = a[ia], b[ib]
                for ja in range(jmax + 1):
                    jb = j - ja
                    if ja < len(arow) and jb < len(brow):
                        total += arow[ja] * brow[jb]
            out[i][j] = total
    return out


def recurrence(b0, c0, p0, order: int, S, wz, wt):
    """Return A/B/C/P coefficient jet lists using the full ST073 recurrence."""
    shape = (len(b0), len(b0[0]))
    B, C, P = [b0], [c0], [p0]
    A = [None]
    # Incompressibility is A_0=-d_z C_0/2; derivative() uses the requested
    # physical source scale, so form it directly before the q-dependent loop.
    A[0] = derivative(C[0], 0, wz)
    for i in range(shape[0]):
        for j in range(shape[1]):
            A[0][i][j] *= -mp.mpf("0.5")
    for n in range(order):
        dz_b = derivative(B[n], 0, wz)
        dz_c = derivative(C[n], 0, wz)
        dz_a = derivative(A[n], 0, wz)
        rb = derivative(B[n], 1, wt)
        for i in range(shape[0]):
            for j in range(shape[1]):
                rb[i][j] *= -1
        add_inplace(rb, derivative(dz_b, 0, wz), factor=-1)
        rc = derivative(C[n], 1, wt)
        for i in range(shape[0]):
            for j in range(shape[1]):
                rc[i][j] *= -1
        add_inplace(rc, derivative(P[n], 0, wz))
        add_inplace(rc, derivative(dz_c, 0, wz), factor=-1)
        ra = derivative(A[n], 1, wt)
        for i in range(shape[0]):
            for j in range(shape[1]):
                ra[i][j] *= -1
        add_inplace(ra, derivative(dz_a, 0, wz), factor=-1)
        for i in range(n + 1):
            j = n - i
            add_inplace(rb, multiply(A[i], B[j], shape), factor=2 * (j + 1))
            add_inplace(rb, multiply(C[i], derivative(B[j], 0, wz), shape))
            add_inplace(rc, multiply(A[i], C[j], shape), factor=2 * j)
            add_inplace(rc, multiply(C[i], derivative(C[j], 0, wz), shape))
            add_inplace(ra, multiply(A[i], A[j], shape), factor=1 + 2 * j)
            add_inplace(ra, multiply(C[i], derivative(A[j], 0, wz), shape))
            add_inplace(ra, multiply(B[i], B[j], shape), factor=-1)
        Bn = rb
        Cn = rc
        for i in range(shape[0]):
            for j in range(shape[1]):
                Bn[i][j] *= S / (4 * (n + 1) * (n + 2))
                Cn[i][j] *= S / (4 * (n + 1) ** 2)
        B.append(Bn)
        C.append(Cn)
        An = derivative(Cn, 0, wz)
        for i in range(shape[0]):
            for j in range(shape[1]):
                An[i][j] *= -mp.mpf(1) / (2 * (n + 2))
        A.append(An)
        Pn = zeros(*shape)
        for i in range(shape[0]):
            for j in range(shape[1]):
                Pn[i][j] = 2 * (n + 2) * An[i][j] - S * ra[i][j] / (2 * (n + 1))
        P.append(Pn)
    return A, B, C, P


def polyval(coeffs, x):
    value = mp.mpf("0")
    for coeff in reversed(coeffs):
        value = value * x + coeff
    return value


def polyder(coeffs, count=1):
    out = list(coeffs)
    for _ in range(count):
        out = [i * out[i] for i in range(1, len(out))]
    return out


class HighPrecisionRadial:
    """One-point ST073-V evaluator with a max-order recurrence cache."""

    def __init__(self, dps: int, max_order: int = max(ORDERS)):
        self.dps = dps
        self.max_order = max_order
        rows = 2 * max_order + JET_EXTRA + 2
        cols = max_order + JET_EXTRA + 1
        shape = (rows, cols)
        z = zeros(rows, cols)
        z[0][0] = ETA_POINT
        z[1][0] = mp.mpf("1")
        powers = {
            -1 - H: q_power_jet(-1 - H, ETA_POINT, shape),
            -1: q_power_jet(-1, ETA_POINT, shape),
            -(mp.mpf("0.5") + H): q_power_jet(-(mp.mpf("0.5") + H), ETA_POINT, shape),
            -2 * (mp.mpf("0.5") + H): q_power_jet(-2 * (mp.mpf("0.5") + H), ETA_POINT, shape),
            -2: q_power_jet(-2, ETA_POINT, shape),
        }
        self.axis_b = [[AXIS_SWIRL * v for v in row] for row in powers[-1 - H]]
        self.axis_c = multiply(z, powers[-1], shape)
        add_inplace(self.axis_c, powers[-(mp.mpf("0.5") + H)], factor=AXIAL_BIAS / AXIAL_SLOPE)
        for i in range(rows):
            for j in range(cols):
                self.axis_c[i][j] *= AXIAL_SLOPE
        self.axis_p = [[-v for v in row] for row in powers[-2 * (mp.mpf("0.5") + H)]]
        p2 = multiply(z, z, shape)
        p2 = multiply(p2, powers[-2], shape)
        add_inplace(self.axis_p, p2, factor=mp.mpf("0.5"))
        self._coeff_cache = {}

    def coefficients(self, k: int):
        if k in self._coeff_cache:
            return self._coeff_cache[k]
        tau = mp.mpf(1) / (mp.mpf(2) ** (k + 1))
        q = tau / (1 - ETA_POINT * ETA_POINT)
        Aexp = mp.mpf("0.5") + H
        Dexp = mp.mpf("0.5") - H
        S, wz, wt = 2 * q, q ** Dexp, q
        B0 = [[q ** (-1 - H) * v for v in row] for row in self.axis_b]
        C0 = [[q ** (-Aexp) * v for v in row] for row in self.axis_c]
        P0 = [[q ** (-2 * Aexp) * v for v in row] for row in self.axis_p]
        raw = recurrence(B0, C0, P0, self.max_order, S, wz, wt)
        self._coeff_cache[k] = (tau, q, raw, S, wz, wt)
        return self._coeff_cache[k]

    def evaluate(self, k: int, order: int):
        tau, q, raw, S, wz, wt = self.coefficients(k)
        X = X_POINT
        values = []
        for seq in raw:
            coeffs = []
            for n in range(order + 1):
                row = seq[n]
                coeffs.append([
                    row[0][0],
                    row[1][0] / wz,
                    2 * row[2][0] / (wz * wz),
                    -row[0][1] / wt,
                ])
            values.append(coeffs)
        v = [[polyval([values[c][n][j] for n in range(order + 1)], X) for j in range(4)] for c in range(4)]
        vr = [[polyval(polyder([values[c][n][0] for n in range(order + 1)], j), X) / (S ** j) for j in (1, 2)] for c in range(4)]
        a, b, c, p = [v[i][0] for i in range(4)]
        az, bz, cz, pz = [v[i][1] for i in range(4)]
        azz, bzz, czz, pzz = [v[i][2] for i in range(4)]
        at, bt, ct, pt = [v[i][3] for i in range(4)]
        ass, bs, cs, ps = [vr[i][0] for i in range(4)]
        ass2, bss, css, pss = [vr[i][1] for i in range(4)]
        r = mp.sqrt(S * X)
        ut = [r * at, r * bt, ct]
        adv = [r * (a * a + 2 * S * X * a * ass + c * az - b * b),
               r * (2 * a * b + 2 * S * X * a * bs + c * bz),
               2 * S * X * a * cs + c * cz]
        grad = [2 * r * ps, mp.mpf("0"), pz]
        lap = [r * (8 * ass + 4 * S * X * ass2 + azz),
               r * (8 * bs + 4 * S * X * bss + bzz),
               4 * cs + 4 * S * X * css + czz]
        res = [ut[i] + adv[i] + grad[i] - lap[i] for i in range(3)]
        div = 2 * a + 2 * S * X * ass + cz
        rn = mp.sqrt(NU)
        physical = [rn * x for x in res]
        normalized = [res[0] * q ** mp.mpf("1.5"),
                      res[1] * q ** (mp.mpf("1.5") + H),
                      res[2] * q ** (mp.mpf("1.5") + H)]
        term_norms = {name: mp.sqrt(sum(x * x for x in term)) for name, term in (
            ("time_derivative", [rn * x for x in ut]),
            ("advection", [rn * x for x in adv]),
            ("pressure_gradient", [rn * x for x in grad]),
            ("viscous_term", [rn * x for x in lap]),
        )}
        return {
            "k": k,
            "order": order,
            "tau": tau,
            "q": q,
            "X": X,
            "eta": ETA_POINT,
            "residual_cylindrical": physical,
            "momentum_norm": mp.sqrt(sum(x * x for x in physical)),
            "normalized_residual": normalized,
            "normalized_momentum_norm": mp.sqrt(sum(x * x for x in normalized)),
            "divergence": div,
            "term_norms": term_norms,
        }


def text_number(value, digits=None):
    if digits is None:
        digits = min(36, max(10, mp.mp.dps - 4))
    return mp.nstr(value, digits)


def serial_row(row):
    out = {}
    for key, value in row.items():
        if isinstance(value, dict):
            out[key] = serial_row(value)
        elif isinstance(value, list):
            out[key] = [text_number(v) if isinstance(v, mp.mpf) else v for v in value]
        elif isinstance(value, mp.mpf):
            out[key] = text_number(value)
        else:
            out[key] = value
    return out


def run(output: Path = OUT):
    started = time.time()
    raw_by_dps = {}
    rows_by_dps = {}
    for dps in DPS_VALUES:
        with mp.workdps(dps):
            evaluator = HighPrecisionRadial(dps)
            raw = {}
            rows = []
            for k in K_VALUES:
                for order in ORDERS:
                    row = evaluator.evaluate(k, order)
                    raw[(order, k)] = row
                    rows.append(serial_row(row))
            raw_by_dps[dps] = raw
            rows_by_dps[str(dps)] = rows
    mismatch = []
    with mp.workdps(max(DPS_VALUES) + 10):
        for order in ORDERS:
            for k in K_VALUES:
                a = raw_by_dps[50][(order, k)]
                b = raw_by_dps[70][(order, k)]
                da = abs(a["momentum_norm"] - b["momentum_norm"])
                mismatch.append({
                    "order": order,
                    "k": k,
                    "dps50_vs_dps70_momentum_abs": text_number(da),
                    "dps50_vs_dps70_momentum_rel": text_number(da / max(abs(b["momentum_norm"]), mp.mpf("1e-200"))),
                    "dps50_vs_dps70_divergence_abs": text_number(abs(a["divergence"] - b["divergence"])),
                })
    report = {
        "schema": "st073_core_critical_horizon_mp_v1",
        "source": "NS_ST073_Full_Local_Recurrence/data/ST073-V.json",
        "independent_path": "mpmath axis Lagrange jets + bivariate recurrence + complete cylindrical momentum",
        "point": {"X": "1/64", "eta": "0.3", "angle": "0"},
        "orders": list(ORDERS),
        "k_values": list(K_VALUES),
        "dps_values": list(DPS_VALUES),
        "rows_by_dps": rows_by_dps,
        "precision_mismatch": mismatch,
        "runtime_seconds": time.time() - started,
        "scope": "Single local similarity point; k beyond registered k<=6 is diagnostic only. No outer matching, forcing, global energy or scale recursion claim.",
        "pde_validated": False,
        "scale_recursion_established": False,
    }
    output.write_bytes((json.dumps(report, indent=2) + "\n").encode())
    print(json.dumps({"output": str(output), "runtime_seconds": report["runtime_seconds"], "precision_mismatch": mismatch}, indent=2))
    return report


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--order", type=int, default=10)
    parser.add_argument("--k", type=int, default=18)
    parser.add_argument("--dps", type=int, default=30)
    parser.add_argument("--full", action="store_true", help="Run the expensive 4-order, 4-scale, 2-precision matrix")
    args = parser.parse_args()
    if args.full:
        run()
    else:
        started = time.time()
        with mp.workdps(args.dps):
            evaluator = HighPrecisionRadial(args.dps, max_order=args.order)
            row = serial_row(evaluator.evaluate(args.k, args.order))
        print(json.dumps(dict(order=args.order, k=args.k, dps=args.dps,
                              row=row, runtime_seconds=time.time() - started), indent=2))
