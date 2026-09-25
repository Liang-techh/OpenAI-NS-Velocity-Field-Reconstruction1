"""Exact-curl wave with an advected local cylindrical phase and support."""
import json

import numpy as np

from curl_wave_prototype import (LocalizedCurlWave, WavePerturbedField,
                                 bump, cylindrical_residual)
from joined_field import ROOT
from joint_collar_fit import kinematics
from local_poloidal_basis_screen import load_robust_candidate
from transported_phase_screen import flow


class TransportedCurlWave(LocalizedCurlWave):
    def __init__(self, source, base, radial_halfwidth=.00275,
                 axial_halfwidth=.000075, time_halfwidth=.00015):
        super().__init__(source, radial_halfwidth, axial_halfwidth,
                         time_halfwidth=time_halfwidth)
        self.trajectories = []
        self.state_cache = {}
        for wave in self.waves:
            kr, _, kz = wave['normal']
            state = np.array([self.radius, self.zcenter, kr, kz, 0.])
            self.trajectories.append(flow(base, wave['m'], self.tau0,
                                          state,
                                          self.tau0-time_halfwidth,
                                          self.tau0+time_halfwidth))

    def mode_state(self, index, tau):
        key = (index, float(tau))
        if key not in self.state_cache:
            state = self.trajectories[index](tau)
            rc, _, kr, kz, _ = state
            wave = self.waves[index]
            normal = np.array([kr, wave['m']/rc, kz])
            amp = wave['amplitude']
            amp = amp-normal*np.dot(normal, amp)/np.dot(normal, normal)
            potential = np.cross(normal, amp)/np.dot(normal, normal)
            self.state_cache[key] = (state, potential)
        return self.state_cache[key]

    def fields(self, points, tau):
        pts = np.asarray(points, float)
        ts = np.broadcast_to(tau, (len(pts),))
        velocity = np.zeros_like(pts)
        for i, (point, t) in enumerate(zip(pts, ts)):
            t = float(t)
            if not self.tau0-self.time_halfwidth < t < self.tau0+self.time_halfwidth:
                continue
            time_cut = bump(t, self.tau0, self.time_halfwidth)[0]
            x, y, z = point
            r = np.hypot(x, y)
            if r == 0:
                continue
            theta = np.arctan2(y, x)
            wr = wt = wz = 0.
            for index, (weight, wave) in enumerate(zip(self.weights, self.waves)):
                state, potential = self.mode_state(index, t)
                rc, zc, kr, kz, phi0 = state
                br, br_r = bump(r, rc, self.radial_halfwidth)
                bz, bz_z = bump(z, zc, self.axial_halfwidth)
                if br == 0 or bz == 0:
                    continue
                envelope = br*bz*time_cut
                er = br_r*bz*time_cut
                ez = br*bz_z*time_cut
                cr, ct, cz = potential
                m = wave['m']
                phase = m*theta+kr*(r-rc)+kz*(z-zc)+phi0
                sn, cs = np.sin(phase), np.cos(phase)
                factor = np.sqrt(weight)
                wr += factor*((-m*cz/r+kz*ct)*envelope*sn-ct*ez*cs)
                wt += factor*((-kz*cr+kr*cz)*envelope*sn
                              +(cr*ez-cz*er)*cs)
                wz += factor*((-kr*ct+m*cr/r)*envelope*sn
                              +ct*(er+envelope/r)*cs)
            ca, sa = x/r, y/r
            velocity[i] = [wr*ca-wt*sa, wr*sa+wt*ca, wz]
        return velocity, np.zeros(len(pts))


def run():
    source = json.loads((ROOT/'compact_potential'/'local_poloidal_10pct_source.json').read_text())
    base = load_robust_candidate(.1)
    frozen = LocalizedCurlWave(source, .00275, .000075,
                               time_halfwidth=.00015)
    transported = TransportedCurlWave(source, base)
    radius, _, z = source['point']
    angles = np.arange(16)*2*np.pi/16
    points = np.column_stack((radius*np.cos(angles), radius*np.sin(angles),
                              np.full(len(angles), z)))
    rows = []
    for tau in (source['tau'], source['tau']+.00006):
        hs = .0005*np.sqrt(base.nu*tau)
        ht = .0001*tau
        def measure(field):
            u, J, part = kinematics(field, points, tau, hs, ht)
            residual = part+np.einsum('nij,nj->ni', J, u)
            cylindrical = cylindrical_residual(residual, points)
            return {'max_residual': float(np.max(np.linalg.norm(residual, axis=1))),
                    'mean_residual_cylindrical': cylindrical.mean(axis=0).tolist(),
                    'max_fd_divergence': float(np.max(np.abs(np.trace(J, axis1=1, axis2=2))))}
        row = {'tau': tau, 'base': measure(base),
               'frozen': measure(WavePerturbedField(base, frozen)),
               'transported': measure(WavePerturbedField(base, transported))}
        rows.append(row)
        print(json.dumps(row), flush=True)
    report = {'source': 'local_poloidal_10pct_source.json',
              'rows': rows,
              'selected_pulse_indices': transported.pulse_indices,
              'scope': 'Two compact exact-curl harmonic waves on a 16-angle ring at two times. Transported phase/center and time-varying spatially constant potentials; no amplitude or pressure evolution equation, global stress closure, or full-domain momentum certificate.',
              'accepted': False}
    path = ROOT/'compact_potential'/'transported_curl_wave.json'
    path.write_bytes((json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
