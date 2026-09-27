"""Reconstruct the integrated diagnostic velocity field from its saved maps."""
import json
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.integrate import solve_ivp
from scipy.optimize import root
from adaptive_bridge_recursive_defect import build_fields
from midplane_outer_moment_dae import FREE
from midplane_resolved_feasibility import ZeroBackground,_evaluate,_jacobian
from separated_moment_modes import SeparatedMomentModes,RADIAL_WINDOWS_THREE
from radial_continuation import ROOT

class SavedOuterDAEField:
    def __init__(self):
        self.report=json.loads((ROOT/'midplane_outer_dae_interval.json').read_text())
        inner,fields=build_fields();self.base=fields['two_sided_cone']
        a=np.array(json.loads((ROOT/'midplane_outer_axial_interscale.json').read_text())['amplitudes'])
        self.current=SeparatedMomentModes(self.base,a,windows=RADIAL_WINDOWS_THREE,knots=(11.,15.,19.),axial_powers=(0,1,2))
        for name in ('inner','nu','join_X','ratio'):setattr(self,name,getattr(self.current,name))
        samples=self.report['coefficient_samples'];self.nodes=np.array([x['k'] for x in samples])
        self.curves=[CubicSpline(self.nodes,np.array([x[key] for x in samples]),axis=0) for key in ('offset','linear','quadratic')]
        self.response=CubicSpline(self.nodes,np.array([x['slope_response'] for x in samples]),axis=0)
        self.seed=np.array(json.loads((ROOT/'midplane_outer_moment_dae_k13.json').read_text())['coefficient_values'])
        self.solution=solve_ivp(self.rhs,(self.nodes[0],self.nodes[-1]),self.seed[:6],rtol=1e-9,atol=1e-11,dense_output=True,max_step=.01)
        if not self.solution.success:raise RuntimeError(self.solution.message)
        np.testing.assert_allclose(self.solution.sol(self.report['integration_nodes']).T,self.report['swirl_values'],rtol=1e-9,atol=1e-11)
    def coefficients(self,k):return tuple(c(k) for c in self.curves)
    def project(self,k,swirl):
        values=self.seed.copy();values[:6]=swirl;c=self.coefficients(k)
        def fun(p):
            values[6:8]=p;return _evaluate(c,values)[[1,3]]
        def jac(p):
            values[6:8]=p;return _jacobian(c,values)[[1,3],6:8]
        sol=root(fun,self.seed[6:8],jac=jac,tol=1e-10)
        if max(abs(fun(sol.x)))>1e-7:raise RuntimeError('Axial projection failed')
        values[6:8]=sol.x;return values.copy()
    def rhs(self,k,swirl):
        v=self.project(k,swirl);m=_evaluate(self.coefficients(k),v)
        return np.linalg.lstsq(self.response(k)[[0,2],:6],-m[[0,2]],rcond=None)[0]
    def fields(self,points,tau):
        k=-np.log2(2*float(np.asarray(tau).ravel()[0]))
        if k<self.nodes[0] or k>self.nodes[-1]:raise ValueError('Outside integrated interval')
        b=np.zeros(18);b[FREE]=self.project(k,self.solution.sol(k))
        correction=SeparatedMomentModes(ZeroBackground(self.base),np.tile(b,3),windows=RADIAL_WINDOWS_THREE,knots=(11.,15.,19.),axial_powers=(0,1,2))
        u,p=self.current.fields(points,tau);v,q=correction.fields(points,tau)
        return u+v,p+q
