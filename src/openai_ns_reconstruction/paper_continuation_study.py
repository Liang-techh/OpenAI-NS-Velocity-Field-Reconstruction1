"""Save successful and failed direct-continuation runs incrementally."""
import json
from pathlib import Path
import numpy as np
from .paper_core_continuation import PaperCoreContinuation
from .paper_profile_residual import residual


def run():
    out=Path('artifacts/function_first/core_continuation')
    out.mkdir(parents=True,exist_ok=True)
    runs=[]
    for nodes in [129,257,513]:
        entry=dict(eta_nodes=nodes,sigma=.3,start=.03,end=.41,rtol=1e-10)
        try:
            f=PaperCoreContinuation(eta_nodes=nodes)
            entry['metadata']=f.metadata();entry['rows']=[]
            for X in [.03,.1,.2,.3,.4]:
                values=[];errors=[]
                for eta in np.linspace(-.95,.95,39):
                    try:
                        values.append(dict(eta=float(eta),residual=residual(
                            f.profile,X,float(eta),step=5e-5).tolist()))
                    except (ValueError,ArithmeticError) as err:
                        errors.append(dict(eta=float(eta),error=str(err)))
                entry['rows'].append(dict(X=X,values=values,errors=errors,
                    maximum=np.max(np.abs([v['residual'] for v in values]),axis=0).tolist() if values else None))
        except (ValueError,ArithmeticError,RuntimeError) as err:
            entry['failure']=str(err)
        runs.append(entry)
        (out/'resolution.json').write_text(json.dumps(dict(
            status='failed to stabilize extension; not selected as final field',runs=runs),indent=2)+'\n')
        print(nodes,entry.get('failure','completed; inspect residuals'),flush=True)
    return runs


if __name__=='__main__':
    run()
