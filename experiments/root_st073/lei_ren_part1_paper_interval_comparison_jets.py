"""Explicit-center-family comparison from a resumable candidate core.

Production loading requires degree 124 and the analytic tail gate. Diagnostic
loading is explicit and does not promote a partial core. RK4 interval arithmetic
does not enclose its discretization error or propagate the analytic core tail.
"""
import hashlib
import json
import operator
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from lei_ren_part1_paper_candidate_comparison_jets import CandidateComparisonJets
from lei_ren_part1_paper_candidate_gauge_core import _unpack_fixed, _unpack_rows
from lei_ren_part1_paper_candidate_interval_core import STEM, tail_budget
from lei_ren_part1_paper_candidate_general_center_factory import PARAMETERS
from lei_ren_part1_paper_candidate_pressure_function import ACCEPTED_SHA
from lei_ren_part1_paper_candidate_shared_inlet import normalized_inlet, radial_product
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE = Path(__file__).resolve().parent


class IntervalComparisonJets(CandidateComparisonJets):
    def __init__(self, steps=32, stem=STEM, diagnostic_partial=False):
        self.steps = operator.index(steps)
        if self.steps < 4:
            raise ValueError('Need at least four comparison integration steps')
        raw = (HERE/(stem+'_state.json')).read_bytes()
        self.state_hash = hashlib.sha256(raw).hexdigest()
        state = json.loads(raw)
        identity = state['identity']
        if identity['parameters'] != PARAMETERS or identity['accepted_schedule_sha256'] != ACCEPTED_SHA:
            raise ValueError('Candidate parameters or accepted pressure datum changed')
        for name, digest in state['input_hashes'].items():
            if hashlib.sha256((HERE/name).read_bytes()).hexdigest() != digest:
                raise ValueError('Interval core dependency changed: '+name)
        self.degree = operator.index(state['completed_radial_order'])
        target = operator.index(identity['target_radial_degree'])
        self.ctx = MPIntervalContext()
        self.ctx.dps = identity['precision']
        self.precision = self.ctx.dps
        self.diagnostic_partial = bool(diagnostic_partial)
        self.center_family = identity['center']
        with mp.workdps(self.precision+40):
            c = self.ctx
            self.tail = tail_budget(c, self.degree)
            if not diagnostic_partial and (target != 124 or self.degree != target or not self.tail['target_met']):
                raise ValueError('Production comparison requires degree 124 and passing analytic core tail')
            fixed = _unpack_fixed(c, state['fixed'])
            rows = {name: _unpack_rows(c, state['rows'][name]) for name in ('A','Uz','P')}
            for rowset in rows.values():
                if len(rowset) != self.degree+1 or any(len(row) != target+4-n for n,row in enumerate(rowset)):
                    raise ValueError('Invalid radial row count or retained axial depth')
            self.lam = c.mpf(PARAMETERS['Lambda'])
            self.eps = 1/self.lam
            self.delta = c.mpf(PARAMETERS['delta'])
            self.z = IntervalTaylor.variable(c, c.mpf(self.center_family), 3)
            def jet(row):
                if len(row) < 4:
                    raise ValueError('Three axial derivative orders required')
                return IntervalTaylor(c, row[:4])
            self.phi = [jet(row)/self.lam**n for n,row in enumerate(rows['A'])]
            self.u = [jet(row)/self.lam**n for n,row in enumerate(rows['Uz'])]
            self.ell = jet(fixed['ell_Z_taylor'])
            self.S = jet(fixed['S_Z_taylor'])
            self.p0 = jet(fixed['P0_Z_taylor'])
            self.F0 = fixed['F0_interval']
            if endpoints(self.F0)[0] <= 0 or endpoints((1-self.z*self.z*self.delta)[0])[0] <= 0:
                raise ValueError('Amplitude or coordinate denominator positivity lost')
            self.hb = mp.mpf('.005')
            self.initial_phi = self.radial(self.phi, c.mpf(4))
            if endpoints(self.initial_phi[0])[0] <= 0:
                raise ValueError('Finite interval core does not establish positive inlet phi; refine axial family')
            self.phi_squared = radial_product(self.phi, self.phi)
            self.phi_u = radial_product(self.phi, self.u)
            self.u_squared = radial_product(self.u, self.u)
        self.cache = {}
        self.radial_cache = {}

    def rhs(self, y, state):
        s = 4*self.ctx.exp(self.ctx.mpf(y))
        if endpoints(self.radial(self.phi,s)[0])[0] <= 0:
            raise ValueError('Radial comparison slope denominator crosses zero; refine axial family')
        return super().rhs(y, state)

    def evaluate(self, y, Z=None):
        """Return bounds over the saved center family, never a scalar slice."""
        if Z is not None:
            raise ValueError('Family evaluator takes no scalar Z; regenerate a scalar core for scalar results')
        with mp.workdps(self.precision+40):
            yy = mp.mpf(str(y))
            if yy < 0:
                raise ValueError('Comparison requires y >= 0')
            key = mp.nstr(yy, self.precision)
            if key in self.cache:
                return self.cache[key]
            c = self.ctx
            s = 4*c.exp(c.mpf(yy))
            if yy > 2*self.hb:
                endpoint = self.evaluate(2*self.hb)
                state = dict(endpoint['state'])
                sb = endpoint['scaled_radius']
                ds, ds2 = s-sb, s*s-sb*sb
                f, u = state['phi'], state['U']
                for name, increment in dict(theta=f*ds2, z=u*ds,
                    theta_z=f*u*ds2, p=f*f*ds, u_squared=u*u*ds,
                    weighted_phi_squared=f*f*(ds2/2)).items():
                    state[name] = state[name]+increment
            elif yy <= self.hb:
                state = self.core_state(s)
            else:
                state = self.core_state(4*c.exp(c.mpf(self.hb)))
                h = (yy-self.hb)/self.steps
                for n in range(self.steps):
                    at = self.hb+n*h
                    k1 = self.rhs(at, state)
                    k2 = self.rhs(at+h/2, self.add(state,k1,c.mpf(h/2)))
                    k3 = self.rhs(at+h/2, self.add(state,k2,c.mpf(h/2)))
                    k4 = self.rhs(at+h, self.add(state,k3,c.mpf(h)))
                    state = {k:state[k]+(k1[k]+k2[k]*2+k3[k]*2+k4[k])*c.mpf(h/6) for k in state}
            if endpoints(state['phi'][0])[0] <= 0:
                raise ValueError('Comparison phi positivity lost; refine axial family or bound propagation')
            zero = state['phi']*0
            inlet = normalized_inlet(self.phi,self.u,state,self.S,self.ell,self.p0,
                self.lam,s,self.z,self.delta,
                endpoints_override=dict(phi_exit=state['phi'],u_exit=state['U'],phi_s=zero,u_s=zero))
            packet = dict(y=yy, scaled_radius=s, state=state, center_family=self.center_family,
                D=inlet['ratio'], I_z=inlet['iz']*c.sqrt(self.eps), pressure=inlet['pressure'],
                Fbar_over_F0=state['phi'], Fa_over_Fbar=self.initial_phi/state['phi'],
                diagnostic_partial=self.diagnostic_partial, analytic_tail_propagated=False,
                ODE_discretization_error_enclosed=False,
                outside_analytic_core_radial_domain=bool(endpoints(s)[1] > mp.mpf('4.1')),
                comparison_driver_only=True, stress_derivatives_supplied=False,
                adapter_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
            self.cache[key] = packet
            return packet

    def exit_driver_jets(self, y, Z, chi):
        if Z is not None:
            raise ValueError('Pass Z=None for an interval-family driver')
        return super().exit_driver_jets(y, None, chi)
