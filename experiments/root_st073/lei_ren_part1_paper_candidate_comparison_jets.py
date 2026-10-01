"""Actual Lambda120 center core propagated through the Section 9.23 ODE.

Axial derivatives use Taylor algebra, with the amplitude factored exactly.
Intervals enclose input/rounding errors of the finite RK4 calculation, NOT
its discretization error or the omitted analytic radial core. Center only.
"""
import hashlib
import json
import operator
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from lei_ren_part1_paper_candidate_shared_inlet import integral, radial_product, normalized_inlet
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE = Path(__file__).resolve().parent
STATE = "lei_ren_part1_paper_candidate_gauge_core_Lambda120_state.json"


def restore(ctx, row):
    return ctx.mpf([mp.make_mpf(tuple(row['lower_exact_mpf_tuple'])),
                    mp.make_mpf(tuple(row['upper_exact_mpf_tuple']))])


class CandidateComparisonJets:
    def __init__(self, steps=32):
        self.steps = operator.index(steps)
        if self.steps < 4:
            raise ValueError('Need at least four comparison integration steps')
        raw = (HERE/STATE).read_bytes()
        self.state_hash = hashlib.sha256(raw).hexdigest()
        state = json.loads(raw)
        target = state['target']
        if target['Lambda'] != '1e120' or state['completed_radial_order'] != 124:
            raise ValueError('Require completed Lambda120 candidate')
        for name, digest in state['source_hashes'].items():
            name = {'driver': target['driver_file'],
                    'pressure_input': target['pressure_file']}.get(name, name)
            if hashlib.sha256((HERE/name).read_bytes()).hexdigest() != digest:
                raise ValueError('Core dependency changed: '+name)
        self.ctx = MPIntervalContext()
        self.ctx.dps = target['precision']
        self.precision = self.ctx.dps
        with mp.workdps(self.precision+40):
            c = self.ctx
            self.center = mp.mpf(target['Z'])
            self.lam = c.mpf(target['Lambda'])
            self.eps = 1/self.lam
            self.delta = c.mpf(target['delta'])
            self.z = IntervalTaylor.variable(c, c.mpf(target['Z']), 3)
            def jet(row):
                if len(row) < 4:
                    raise ValueError('Three axial orders required')
                return IntervalTaylor(c, [restore(c, v) for v in row[:4]])
            self.phi = [jet(row)/self.lam**n for n,row in enumerate(state['A_rows'])]
            self.u = [jet(row)/self.lam**n for n,row in enumerate(state['Uz_rows'])]
            fixed = state['fixed_jets']
            self.ell = jet(fixed['ell_Z_taylor'])
            self.S = jet(fixed['S_Z_taylor'])
            self.p0 = jet(fixed['P0_Z_taylor'])
            self.F0 = restore(c, fixed['F0_interval'])
            self.hb = mp.mpf('.005')
            self.initial_phi = self.radial(self.phi, c.mpf(4))
            self.phi_squared = radial_product(self.phi,self.phi)
            self.phi_u = radial_product(self.phi,self.u)
            self.u_squared = radial_product(self.u,self.u)
        self.cache = {}
        self.radial_cache = {}

    @staticmethod
    def radial(rows, s, derivative=0):
        out = rows[0]*0
        for n,row in enumerate(rows):
            if n >= derivative:
                out += row*s**(n-derivative)*(n if derivative else 1)
        return out

    def core_state(self, s):
        return dict(phi=self.radial(self.phi,s), U=self.radial(self.u,s),
                    theta=integral(self.phi,s,1)*2,z=integral(self.u,s),
                    theta_z=integral(self.phi_u,s,1)*2,p=integral(self.phi_squared,s),
                    u_squared=integral(self.u_squared,s),
                    weighted_phi_squared=integral(self.phi_squared,s,1))

    def rhs(self, y, state):
        c = self.ctx
        s = 4*c.exp(c.mpf(y))
        key = mp.nstr(y,self.precision)
        if key not in self.radial_cache:
            phi = self.radial(self.phi,s)
            self.radial_cache[key] = (self.radial(self.phi,s,1)*s/phi,
                                      self.radial(self.u,s,1)*s)
        slope,u_slope = self.radial_cache[key]
        if y <= self.hb:
            alpha = c.mpf(1)
        elif y >= 2*self.hb:
            alpha = c.mpf(0)
        else:
            x = c.mpf((y-self.hb)/self.hb)
            left, right = c.exp(-1/(x*x)), c.exp(-1/((1-x)**2))
            alpha = right/(left+right)
        f, u = state['phi'], state['U']
        return dict(phi=f*slope*alpha, U=u_slope*alpha,
                    theta=f*(2*s*s), z=u*s, theta_z=f*u*(2*s*s),
                    p=f*f*s, u_squared=u*u*s, weighted_phi_squared=f*f*s*s)

    @staticmethod
    def add(a,b,factor):
        return {k:a[k]+b[k]*factor for k in a}

    def evaluate(self,y,Z='.3'):
        with mp.workdps(self.precision+40):
            yy, zz = mp.mpf(str(y)), mp.mpf(str(Z))
            if abs(zz-self.center) > mp.mpf(10)**(-self.precision+5):
                raise ValueError('This completed tensor supplies only Z=.3 jets; no extrapolation')
            if yy < 0:
                raise ValueError('Comparison requires y >= 0')
            key = mp.nstr(yy,self.precision)
            if key in self.cache:
                return self.cache[key]
            c = self.ctx
            s = 4*c.exp(c.mpf(yy))
            if yy > 2*self.hb:
                endpoint = self.evaluate(2*self.hb,Z)
                state = dict(endpoint['state'])
                sb = endpoint['scaled_radius']
                ds, ds2 = s-sb, s*s-sb*sb
                f,u = state['phi'],state['U']
                for name,increment in dict(theta=f*ds2,z=u*ds,
                    theta_z=f*u*ds2,p=f*f*ds,u_squared=u*u*ds,
                    weighted_phi_squared=f*f*(ds2/2)).items():
                    state[name] = state[name]+increment
            elif yy <= self.hb:
                state = self.core_state(s)
            else:
                state = self.core_state(4*c.exp(c.mpf(self.hb)))
                h = (yy-self.hb)/self.steps
                for n in range(self.steps):
                    at = self.hb+n*h
                    k1 = self.rhs(at,state)
                    k2 = self.rhs(at+h/2,self.add(state,k1,c.mpf(h/2)))
                    k3 = self.rhs(at+h/2,self.add(state,k2,c.mpf(h/2)))
                    k4 = self.rhs(at+h,self.add(state,k3,c.mpf(h)))
                    state = {k:state[k]+(k1[k]+k2[k]*2+k3[k]*2+k4[k])*c.mpf(h/6)
                             for k in state}
            zero = state['phi']*0
            inlet = normalized_inlet(self.phi,self.u,state,self.S,self.ell,self.p0,
                self.lam,s,self.z,self.delta,
                endpoints_override=dict(phi_exit=state['phi'],u_exit=state['U'],
                                         phi_s=zero,u_s=zero))
            # I_z scales by sqrt(eps); D=I_theta/F has its factor cancelled.
            packet = dict(y=yy, scaled_radius=s, state=state,
                          D=inlet['ratio'], I_z=inlet['iz']*c.sqrt(self.eps),
                          pressure=inlet['pressure'],
                          Fbar_over_F0=state['phi'], Fa_over_Fbar=self.initial_phi/state['phi'])
            self.cache[key] = packet
            return packet

    def exit_driver_jets(self,y,Z,chi):
        """A/B and their analytic Z derivatives; chi depends only on y."""
        p = self.evaluate(y,Z)
        c = self.ctx
        A = p['D']*(-c.mpf(chi)/2)
        B = p['Fa_over_Fbar']*p['I_z']*(-c.mpf(chi)*c.sqrt(p['scaled_radius']*self.eps/2))
        return dict(A=A, B=B, A_Z=A[1], B_Z=B[1])


def run():
    calc = CandidateComparisonJets(32)
    with mp.workdps(calc.precision+40):
        start = calc.evaluate('0')
        end = calc.evaluate('.01')
        driver = calc.exit_driver_jets('.01','.3','1')
        finer = CandidateComparisonJets(64).evaluate('.01')
        changes = {name:[abs(end[name][k]-finer[name][k]) for k in range(2)]
                   for name in ('D','I_z')}
        reference_name = 'lei_ren_part1_paper_candidate_shared_inlet_Lambda120.json'
        reference = json.loads((HERE/reference_name).read_text(encoding='utf-8'))
        if reference['state_sha256'] != calc.state_hash:
            raise AssertionError('Reference inlet comes from a different core')
        for k,row in enumerate(reference['analytic_core_inlet_enclosures']['ratio']):
            lo,hi = endpoints(restore(calc.ctx,row))
            fl,fh = endpoints(start['D'][k])
            if not lo <= fl <= fh <= hi:
                raise AssertionError(('Comparison inlet disagrees with analytic core',k))
        frozen = calc.evaluate(mp.log(110*mp.mpf('1e120')/4))
        if frozen['state']['phi'].coefficients != end['state']['phi'].coefficients:
            raise AssertionError('Frozen comparison angular field changed')
        # Outside the transition the comparison fields are frozen in R.
        # These checks compare actual core-derived fields; no free fit.
        if endpoints(end['Fbar_over_F0'][0])[0] <= 0:
            raise AssertionError('Comparison angular amplitude lost positivity')
        output = dict(state_sha256=calc.state_hash, steps=32, refinement_steps=64,
                      actual_candidate_core_used=True, axial_finite_difference_used=False,
                      amplitude_cancelled_algebraically_not_reset=True,
                      axial_center='.3', comparison_y_interval=['0','.01'],
                      analytic_core_inlet_consistency_checks=3,
                      frozen_comparison_propagated_to_physical_R='110',
                      frozen_R110_D=list(frozen['D'].coefficients),
                      frozen_R110_I_z=list(frozen['I_z'].coefficients),
                      endpoint_D=list(end['D'].coefficients),
                      endpoint_I_z=list(end['I_z'].coefficients),
                      endpoint_phi=list(end['Fbar_over_F0'].coefficients),
                      endpoint_U=list(end['state']['U'].coefficients),
                      normalized_exit_A=list(driver['A'].coefficients),
                      normalized_exit_B=list(driver['B'].coefficients),
                      step_refinement_absolute_changes=changes,
                      radial_infinite_tail_propagated=False,
                      ODE_discretization_error_enclosed=False,
                      whole_axis_transition_generated=False,
                      terminal_five_moment_closure=False,
                      input_hashes={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                          for name in (STATE,reference_name,
                              'lei_ren_part1_paper_candidate_shared_inlet.py',
                              'lei_ren_part1_paper_interval_taylor.py')},
                      source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
        Path(__file__).with_suffix('.json').write_text(json.dumps(encode(output),indent=2)+'\n',encoding='utf-8')
        print('Actual Lambda120 comparison and analytic driver jets generated; center-only RK4, not certified ODE closure',flush=True)
        return output


if __name__ == '__main__':
    run()
