"""Ask whether pressure can finish the eight-mode mean field's cone window.

Uses the cached full nonlinear mean response. Pressure leaves velocity and
shear fixed, so free-axial stress feasibility is an exact pointwise
necessary condition before fitting any pressure modes.
"""

import json

import numpy as np
from scipy.optimize import linprog

from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from delayed_multimode_cone_fit import BASE_NAME, load_cache
from delayed_pressure_cone_geometry import line_interval
from delayed_pressure_cone_window import pressure_stress_columns
from delayed_similarity_pressure_screen import SimilarityPressurePatch
from radial_continuation import ROOT


def geometry(row, coefficients):
    u = row['u0']+row['U']@coefficients
    grad = row['g0']+np.einsum('abi,i->ab', row['G'], coefficients)
    target = (row['T0']+row['Tlin']@coefficients
              +np.einsum('aij,i,j->a', row['Tquad'],
                         coefficients, coefficients))
    F = u[1]/row['radius']
    shear = np.array([grad[1, 0]-F, grad[2, 0]])
    N = shear/np.linalg.norm(shear)
    K = np.array([-N[1], N[0]])
    lam2 = -2*F*N[0]*(2*F*N[0]+np.linalg.norm(shear))
    multiplier = (float(np.sqrt(lam2)/(2*F*N[0]))
                  if lam2 > 0 and abs(2*F*N[0]) > 1e-14 else None)
    return target, N, K, multiplier, float(lam2)


def run():
    source = json.loads((ROOT/'delayed_multimode_cone_fit.json').read_text())
    coefficients = np.asarray(source['coefficients'])
    rows, _ = load_cache(ROOT/source['response_cache'])
    base = CoupledMomentPhysicalLift(slice_filename=BASE_NAME)
    pressure = SimilarityPressurePatch(base, np.zeros(6))
    tau = source['tau']
    results, A, b = [], [], []
    for row in rows:
        target, N, K, multiplier, lam2 = geometry(row, coefficients)
        if multiplier is None:
            results.append(dict(X=row['X'], eta=row['eta'],
                                lambda_squared=lam2,
                                free_axial_080=False,
                                free_axial_0999=False))
            continue
        strict = [.8*N+sign*multiplier*K for sign in (-1., 1.)]
        relaxed = [.999*N+sign*multiplier*K for sign in (-1., 1.)]
        feasible_080 = line_interval(strict, target, .1, 1) is not None
        feasible_0999 = line_interval(relaxed, target, 0., 1) is not None
        results.append(dict(X=row['X'], eta=row['eta'],
                            lambda_squared=lam2,
                            target=target.tolist(),
                            free_axial_080=feasible_080,
                            free_axial_0999=feasible_0999))
        point = base.compact.joined.inner.from_similarity(
            np.array([row['X']]), np.array([row['eta']]), tau)[0]
        columns = pressure_stress_columns(pressure, point, tau)
        for direction in strict:
            A.append(direction@columns)
            b.append(-.1-direction@target)
    linear_feasible = False
    if len(A) == 2*len(rows):
        A, b = np.asarray(A), np.asarray(b)
        scale = 1/np.maximum(np.max(np.abs(A), axis=0), 1e-12)
        linear_feasible = bool(linprog(
            np.zeros(6), A_ub=A*scale, b_ub=b,
            bounds=[(None, None)]*6, method='highs').success)
    report = dict(source='delayed_multimode_cone_fit.json',
                  tau=tau, rows=results,
                  pressure_linear_feasible=linear_feasible,
                  scope='Pointwise free-axial geometry and six compact '
                        'pressure-mode linear feasibility for the '
                        'selected eight-mode mean field. No wave or '
                        'continuous cone acceptance.',
                  accepted=False)
    (ROOT/'delayed_multimode_pressure_admission.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps(dict(
        free_axial_080=sum(r['free_axial_080'] for r in results),
        free_axial_0999=sum(r['free_axial_0999'] for r in results),
        total=len(results),
        pressure_linear_feasible=linear_feasible,
        failures=[(r['X'], r['eta']) for r in results
                  if not r['free_axial_0999']])), flush=True)


if __name__ == '__main__':
    run()
