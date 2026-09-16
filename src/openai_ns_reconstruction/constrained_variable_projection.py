"""Bounded variable projection of pressure and restricted force coefficients."""
from dataclasses import replace
import numpy as np
from scipy.optimize import lsq_linear
from .constrained_force import RestrictedForce


def fit_linear_coefficients(candidate, points, times, nu, step=0.001):
    from .constrained_optimize import training_residual
    names=('pressure_constant','pressure_radial','pressure_axial')
    c=replace(candidate,**dict.fromkeys(names,0.0))
    base=training_residual(c,RestrictedForce(0,0),points,times,nu,step)
    columns=[]
    for name in names:
        basis=replace(c,**{name:1.0})
        grad=np.zeros_like(base)
        for j in range(3):
            d=np.eye(3)[j]*step
            grad[:,j]=(basis.pressure(points+d,times)-basis.pressure(points-d,times))/(2*step)
        columns.append(grad.ravel())
    columns += [-RestrictedForce(1,0)(points,times).ravel(),
                -RestrictedForce(0,1)(points,times).ravel()]
    matrix=np.column_stack(columns)
    fit=lsq_linear(matrix,-base.ravel(),bounds=([-100]*3+[0,0],[100]*3+[10,10]),tol=1e-10,max_iter=100)
    if not fit.success:
        raise RuntimeError('bounded pressure/force solve failed: '+fit.message)
    return replace(c,**dict(zip(names,fit.x[:3]))),RestrictedForce(*fit.x[3:])
