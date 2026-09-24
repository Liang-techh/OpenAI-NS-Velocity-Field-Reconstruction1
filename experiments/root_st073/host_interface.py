"""Regenerate portable interface traces from frozen ST073 JSON, without NPZ dtype reinterpretation."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parent
BUNDLE=ROOT/'NS_ST073_Full_Local_Recurrence'
sys.path.insert(0,str(BUNDLE))
from full_radial import FullRadialField

def run():
    f=FullRadialField.load(BUNDLE/'data/ST073-V.json')
    out=ROOT/'host_interface';out.mkdir(exist_ok=True)
    reference=json.loads((BUNDLE/'evidence/interface_audit.json').read_text())['rows']
    rows=[]
    for k in (0,3,6):
        tau=.5*2**(-k);eta=np.linspace(-.5,.5,17);X=np.full(17,1/64)
        d=f.evaluate_similarity(X,eta,tau)
        portable={key:np.asarray(v,dtype=np.float64) for key,v in d.items()}
        np.savez_compressed(out/f'V_k{k}_float64.npz',eta=eta,X=X,tau=tau,**portable)
        previous=next(r for r in reference if r['id']=='ST073-V' and r['k']==k)
        pressure_error=float(np.max(np.abs(portable['pressure_gradient']-previous['physical_pressure_gradient'])))
        rows.append(dict(k=k,tau=tau,pressure_gradient_reference_max_difference=pressure_error,fields={key:list(v.shape) for key,v in portable.items()},maximum_momentum_residual=float(np.max(np.linalg.norm(portable['residual'],axis=-1)))))
    report=dict(model_sha256=hashlib.sha256((BUNDLE/'data/ST073-V.json').read_bytes()).hexdigest(),numpy_version=np.__version__,host_longdouble_epsilon=float(np.finfo(np.longdouble).eps),output_dtype='float64',original_npz_issue='Host NumPy cannot load archived <f16 arrays; originals preserved byte-for-byte. Values regenerated from model, not reinterpreted.',rows=rows,global_field_ready=False,scope='17 interface nodes at three registered scales only; not exterior matching or finite-total-energy validation')
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':run()
