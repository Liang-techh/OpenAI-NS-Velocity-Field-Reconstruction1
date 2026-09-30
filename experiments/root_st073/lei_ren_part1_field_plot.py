"""Standalone scientific plot of the measured background scale diagnostics."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def run():
    path=Path(__file__).with_name('lei_ren_part1_field_diagnostics.json')
    report=json.loads(path.read_text(encoding='utf-8'))
    rows=report['rows']; tau=np.array([r['tau'] for r in rows])
    fig,axes=plt.subplots(2,2,figsize=(10,7),layout='constrained')
    def values(key): return np.array([r[key] for r in rows])
    ax=axes[0,0]
    for key,label in [('radial_length','Radial scale'),('axial_length','Axial scale')]:
        y=values(key); ax.loglog(tau,y/y[0],'o-',label=label)
    ax.set_ylabel('Length / initial length'); ax.legend(); ax.grid(alpha=.3)
    ax=axes[0,1]
    y=values('axial_to_radial_aspect')
    ax.semilogx(tau,y/y[0],'o-')
    ax.set_ylabel('Axial / radial aspect, normalized')
    ax.set_title('Weak relative elongation: h=0.001'); ax.grid(alpha=.3)
    ax=axes[1,0]
    for key,label in [('swirl','Swirl speed'),('axial_velocity','Axial speed'),('axial_vorticity','Axial vorticity')]:
        y=values(key); ax.loglog(tau,y/y[0],'o-',label=label)
    ax.set_ylabel('Magnitude / initial magnitude'); ax.legend(); ax.grid(alpha=.3)
    ax=axes[1,1]
    y=np.array([r['global_energy']['total'] for r in rows])
    ax.semilogx(tau,y,'o-',label='Full energy including radial heat tail')
    ax.set_ylabel('Kinetic energy'); ax.legend(fontsize=8); ax.grid(alpha=.3)
    for ax in axes.flat:
        ax.invert_xaxis(); ax.set_xlabel('tau = T - t')
    fig.suptitle('Self-similar representation measurements; NS stress/recursion closure remains open')
    destination=path.with_suffix('.png')
    fig.savefig(destination,dpi=160); plt.close(fig)
    print(destination)


if __name__=='__main__': run()
