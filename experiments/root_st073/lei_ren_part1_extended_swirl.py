"""Actual angular-matched extended swirl candidate, before pressure iteration.

Keeps the saved finite core near the axis, transports its negative shear
into a weak power-law anchor, and joins the actual heat collar with a flat
logarithmic blend. A scalar blend timing is solved separately for each Z
from the ACTUAL inward angular target. This is not the paper's complete
outer construction: axial moments, common core-pressure iteration and the
whole stress cone still have to be constructed.
"""
from __future__ import annotations

from functools import lru_cache
import math

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.special import expit

from lei_ren_part1_exterior_targets import exterior_targets
from lei_ren_part1_heat_collar import HeatCollar
from lei_ren_part1_pressure_core import load_core


def _step(x, shift=0.0):
    x = np.asarray(x, dtype=float)
    result = np.zeros_like(x)
    result[x >= 1] = 1
    mask = (x > 0) & (x < 1)
    result[mask] = expit(-1/x[mask] + 1/(1-x[mask]) + shift)
    return result


class ExtendedSwirl:
    def __init__(self, core=None, collar=None, *, R_core=.005,
                 R_anchor=.01, anchor_power=.005, n=128):
        self.core = load_core()[0] if core is None else core
        self.collar = HeatCollar() if collar is None else collar
        self.R_core, self.R_anchor = float(R_core), float(R_anchor)
        self.anchor_power, self.n = float(anchor_power), int(n)
        if not 0 < self.R_core < self.R_anchor < self.collar.collar_inner_radius:
            raise ValueError("invalid layer radii")
        if not 0 < self.anchor_power < .1 or self.n < 32:
            raise ValueError("invalid anchor power or quadrature")

    @lru_cache(maxsize=1024)
    def _local(self, z):
        r0, r1, p = self.R_core, self.R_anchor, self.anchor_power
        def rhs(r, value):
            chi = 1-float(_step((r-r0)/(r1-r0)))
            return [chi*float(self.core.F_radial_derivative(r,z))
                    -(1-chi)*p*value[0]/r]
        sol = solve_ivp(rhs, (r0,r1), [float(self.core.F(r0,z))],
                        rtol=2e-11, atol=1e-15, dense_output=True,
                        max_step=(r1-r0)/16)
        if not sol.success or sol.y[0,-1] <= 0:
            raise ArithmeticError("local negative-shear transport failed")
        return sol.sol, float(sol.y[0,-1])

    def _anchor(self, r, z):
        return self._local(z)[1]*(np.asarray(r)/self.R_anchor)**(-self.anchor_power)

    @lru_cache(maxsize=65536)
    def _outer(self, r, z):
        # This analytic continuation is only a blend target below the collar;
        # its stress is NOT inherited from the genuine collar.
        y = math.log(self.collar.R_b/r)
        flat = math.exp(-4/y**2) if y > 0 else 0.0
        return (1-self.collar.epsilon*flat)*self.collar.F_heat(r,z)

    @lru_cache(maxsize=1024)
    def _cross(self,z):
        ra = self.collar.collar_inner_radius
        def gap(logr):
            r=math.exp(logr)
            return float(self._anchor(r,z))-self._outer(r,z)
        a,b = math.log(self.R_anchor),math.log(ra)
        if gap(b) <= 0:
            raise ArithmeticError("anchor below actual collar endpoint")
        return self.R_anchor if gap(a) >= 0 else math.exp(brentq(gap,a,b,xtol=1e-13))

    def _integral(self, z, shift, density, n):
        nodes,weights=leggauss(n)
        r0,r1,rc,ra=self.R_core,self.R_anchor,self._cross(z),self.collar.collar_inner_radius
        value=0.0
        for a,b in ((0,r0),(r0,r1)):
            rs=a+(b-a)*(nodes+1)/2
            fs=np.array([self._prefix(float(r),z) for r in rs])
            value += float(np.dot(weights*(b-a)/2,density(rs,fs)))
        for a,b in ((r1,rc),(rc,ra)):
            if b <= a: continue
            logs=math.log(a)+(math.log(b/a))*(nodes+1)/2
            rs=np.exp(logs)
            anchors=self._anchor(rs,z)
            if a == rc:
                x=(logs-math.log(rc))/math.log(ra/rc)
                s=_step(x,shift)
                outer=np.array([self._outer(float(r),z) for r in rs])
                fs=(1-s)*anchors+s*outer
            else: fs=anchors
            value += float(np.dot(weights*math.log(b/a)/2*rs,density(rs,fs)))
        return value

    def _prefix(self,r,z):
        if r <= self.R_core: return float(self.core.F(r,z))
        if r < self.R_anchor: return float(self._local(z)[0](r)[0])
        return float(self._anchor(r,z))

    @lru_cache(maxsize=1024)
    def timing(self,z):
        z=float(z)
        target=exterior_targets(self.collar.F,self.collar.heat,
                               self.collar.collar_inner_radius,self.collar.R_b,z)
        wanted=target["angular_target"]
        def defect(shift):
            return self._integral(z,shift,lambda r,f:2*r*f,self.n)-wanted
        low,high=defect(-30),defect(30)
        if low*high > 0:
            raise ArithmeticError(f"angular target unbracketed at Z={z}: {low}, {high}")
        shift=brentq(defect,-30,30,xtol=2e-11)
        return {"Z":z,"shift":shift,"cross_radius":self._cross(z),
                "angular_target":wanted,"bracket_defects":[low,high],
                "angular_solve_defect":defect(shift)}

    def F(self,r,z):
        r,z=float(r),float(z)
        if not math.isfinite(r) or not math.isfinite(z) or r < 0 or abs(z)>1:
            raise ValueError("require R>=0 and |Z|<=1")
        if r >= self.collar.R_b: return self.collar.F_heat(r,z)
        if r >= self.collar.collar_inner_radius: return self.collar.F(r,z)
        rc=self._cross(z)
        if r <= rc: return self._prefix(r,z)
        x=math.log(r/rc)/math.log(self.collar.collar_inner_radius/rc)
        s=float(_step(x,self.timing(z)["shift"]))
        return (1-s)*float(self._anchor(r,z))+s*self._outer(r,z)

    def F_R(self,r,z):
        r,z=float(r),float(z)
        if r < 0 or not math.isfinite(r) or not math.isfinite(z) or abs(z)>1:
            raise ValueError("require finite R>=0 and |Z|<=1")
        if r <= self.R_core: return float(self.core.F_radial_derivative(r,z))
        if r < self.R_anchor:
            chi=1-float(_step((r-self.R_core)/(self.R_anchor-self.R_core)))
            return chi*float(self.core.F_radial_derivative(r,z))-(1-chi)*self.anchor_power*self.F(r,z)/r
        if r >= self.collar.R_b: return self.collar.F_heat_R(r,z)
        if r >= self.collar.collar_inner_radius: return self.collar.F_R(r,z)
        anchor=float(self._anchor(r,z))
        da=-self.anchor_power*anchor/r
        rc=self._cross(z)
        if r <= rc: return da
        length=math.log(self.collar.collar_inner_radius/rc)
        x=math.log(r/rc)/length
        s=float(_step(x,self.timing(z)["shift"]))
        ds=s*(1-s)*(1/x**2+1/(1-x)**2)/(r*length)
        y=math.log(self.collar.R_b/r)
        flat=math.exp(-4/y**2)
        factor=1-self.collar.epsilon*flat
        db=(factor*self.collar.F_heat_R(r,z)
            +8*self.collar.epsilon*flat/(r*y**3)*self.collar.F_heat(r,z))
        return (1-s)*da+s*db+ds*(self._outer(r,z)-anchor)

    def moments(self,z,*,n=256):
        shift=self.timing(float(z))["shift"]
        angular=self._integral(float(z),shift,lambda r,f:2*r*f,n)
        pressure=self._integral(float(z),shift,lambda r,f:f*f,n)
        swirl_energy=self._integral(float(z),shift,lambda r,f:r*f*f,n)
        targets=exterior_targets(self.collar.F,self.collar.heat,
                                 self.collar.collar_inner_radius,self.collar.R_b,float(z))
        p0=targets["pressure_at_inner"]-pressure
        return {**self.timing(float(z)),"angular_holdout_defect":angular-targets["angular_target"],
                "pressure_integral_inner":pressure,"new_axis_pressure":p0,
                "old_axis_pressure":float(self.core.Pi(0,float(z))),
                "required_inner_axial_energy":targets["quadratic_target"]+swirl_energy,
                "common_pressure_core_rebuilt":False,"whole_cone_validated":False}

    def axis_pressure(self,z,*,n=128):
        z=float(z)
        shift=self.timing(z)["shift"]
        return self.collar.P(self.collar.collar_inner_radius,z)-self._integral(
            z,shift,lambda r,f:f*f,n)
