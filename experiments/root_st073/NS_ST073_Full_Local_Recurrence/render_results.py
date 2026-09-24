"""Plot already executed calibration; does not run optimization or a MATLAB UI."""
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parent
rows=json.loads((root/'evidence/full_order_calibration.json').read_text())
fig,ax=plt.subplots(figsize=(8.5,5))
ax.semilogy([r['order'] for r in rows],[r['max'] for r in rows],'-o',label='Full physical residual: calibration maximum')
ax.axhline(1e-3,linestyle='--',label='0.001 reference threshold')
ax.set_xticks([r['order'] for r in rows]);ax.set_xlabel('Radial Taylor order (not time scale)');ax.set_ylabel('Full-vector physical residual')
ax.set_title('ST073-F: local full-momentum recurrence\nFixed calibration points; not global NS acceptance')
ax.grid(True,which='both',alpha=.3);ax.legend();fig.tight_layout()
(root/'figures').mkdir(exist_ok=True);fig.savefig(root/'figures/full_momentum_radial_convergence.png',dpi=180);plt.close(fig)
print('Created actual calibration diagnostic; no global or native MATLAB result.')
