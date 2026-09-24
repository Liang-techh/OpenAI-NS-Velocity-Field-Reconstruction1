"""Appendix B.3 scalar comparison profile, NOT a nonlinear axis-core solution.
Source: OpenAI-hosted navier-stokes.pdf, printed page 146, (B.11)-(B.13).
The argument is Y*chi, not the old physical X. Constants remain uncalibrated.
"""
import numpy as np
from scipy.special import hyp0f1

def comparison(z):
    z=np.asarray(z,float)
    if np.any(~np.isfinite(z)) or np.any(z<0) or np.any(z>4.1):raise ValueError('Registered scalar comparison interval is [0,4.1]')
    f=hyp0f1(2,-z/2);fp=-.25*hyp0f1(3,-z/2)
    return dict(f=f,derivative=fp,logarithmic_shear=-2*z*fp/f)
