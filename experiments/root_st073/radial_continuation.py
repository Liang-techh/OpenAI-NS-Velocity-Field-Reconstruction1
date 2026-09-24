"""Explicit radial continuation experiment; frozen ST073 data remains unchanged."""
from dataclasses import replace
from pathlib import Path
import sys,json
import numpy as np
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'NS_ST073_Full_Local_Recurrence'))
from full_radial import FullRadialField

def run():
    source=FullRadialField.load(ROOT/'NS_ST073_Full_Local_Recurrence/data/ST073-V.json')
    fields={n:FullRadialField(replace(source.p,order=n,X_max=.125)) for n in (8,10)}
    eta=np.linspace(-.5,.5,25);rows=[]
    for k in (0,3,6):
        tau=.5*2**(-k)
        for x in (1/64,.0234375,.03125,.046875,.0625,.09375,.125):
            results={n:f.evaluate_similarity(np.full(len(eta),x),eta,tau) for n,f in fields.items()}
            rows.append(dict(k=k,X=x,radius_ratio_to_frozen=float(np.sqrt(x/(1/64))),residual_max_by_order={n:float(np.max(np.linalg.norm(d['residual'],axis=1))) for n,d in results.items()},velocity_order_difference=float(np.max(np.linalg.norm(results[10]['velocity']-results[8]['velocity'],axis=1))),pressure_order_difference=float(np.max(np.abs(results[10]['pressure']-results[8]['pressure'])))))
    report=dict(rows=rows,scope='Explicit autonomous extension beyond frozen Xmax; same axis data, nu=.01, unforced. Boundary samples at 25 eta nodes only. No whole-annulus L2, finite-total-energy or exterior matching acceptance.',global_field_ready=False)
    out=ROOT/'radial_continuation';out.mkdir(exist_ok=True);(out/'screen.json').write_bytes((json.dumps(report,indent=2)+'\n').encode());print(json.dumps([r for r in rows if r['k']==6],indent=2))
if __name__=='__main__':run()
