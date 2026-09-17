"""Periodic-box spectral diagnostics for a compactly supported candidate."""
import json
from pathlib import Path
import numpy as np
from .constrained_candidate import CompactCandidate


def run():
    folder=Path('artifacts/constrained/adaptive_v4')
    c=CompactCandidate.load(folder/'candidate.json');rows=[]
    for n in (32,48,64):
        axis=-2+4*np.arange(n)/n
        points=np.stack(np.meshgrid(axis,axis,axis,indexing='ij'),axis=-1)
        modes=np.fft.fftfreq(n)*n
        kk=np.stack(np.meshgrid(modes,modes,modes,indexing='ij'),axis=-1)
        radius=np.floor(np.linalg.norm(kk,axis=-1)).astype(int)
        tail=np.max(np.abs(kk),axis=-1)>=.8*(n/2)
        for t in (.25,.5,.75):
            u=c.velocity(points,t);hat=np.fft.fftn(u,axes=(0,1,2))/n**3
            power=.5*64*np.sum(np.abs(hat)**2,axis=-1)
            energy=float(power.sum());direct=float(.5*64*np.mean(np.sum(u*u,axis=-1)))
            assert abs(energy-direct)<1e-10
            shells=np.bincount(radius.ravel(),weights=power.ravel())
            rows.append({'n':n,'time':t,'energy':energy,'tail_energy_fraction':float(power[tail].sum()/energy),
                'shell_energy':shells.tolist()})
    report={'candidate':'adaptive_v4/candidate.json','box':[-2,2],
        'shell_coordinate':'integer Fourier mode radius; physical wavenumber = pi/2 times mode',
        'scope':'finite-window spectral diagnostic, no blow-up inference; no spectral acceptance threshold was preregistered',
        'rows':rows}
    (folder/'spectrum.json').write_text(json.dumps(report,indent=2)+'\n')
    print([(r['n'],r['time'],r['tail_energy_fraction']) for r in rows])

if __name__=='__main__':run()
