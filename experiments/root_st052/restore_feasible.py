"""Pre-holdout feasibility repair, with unchanged physical/structure constraints.
Only the auxiliary epigraph number is freely recomputed. Actual field coefficients
are changed solely by bounded, saved minimum-norm constraint-restoration steps.
"""
from pathlib import Path
import json,time
import numpy as np
from scipy.optimize import minimize
from minimax_exchange import ROOT,EdgeModel,AlignedObjective,atomic_json,scan,pool_points

def main():
    run=ROOT/'evidence_st052/run_M';dest=ROOT/'evidence_st052/feasibility_repair'
    if dest.exists():raise ValueError('Refuse overwriting repair')
    dest.mkdir();start=time.monotonic()
    m=EdgeModel(ROOT/'artifacts/research/ST051-B/candidate.json',edge_modes=True)
    o=AlignedObjective(m,shear_reference=ROOT/'artifacts/research/ST048-S/candidate.json',robust_grid=True,poisson_weight=.003,morph_ratio=.9999,axis_ratio=1.75,shear_ratio=.999,edge_weight=.08,pressure_target=0.,profile_ratio=1.,pressure_ratio=1.,acceleration_ratio=1.,axis_weight=.02,space_order=(32,48),anchor_tolerance=.08,time_cap=1.,peak_weight=0.,softmax_weight=.03)
    T=np.load(run/'coordinate_map.npy');lo,hi=np.array(m.bounds).T;bs=1/np.maximum(hi-lo,1e-4)
    c=np.load(run/'round2/modifiers.npy');hist=[]
    def evaluate(cc):
        v,J=o.constraints(cc)
        return np.r_[v,(cc-lo)*bs,(hi-cc)*bs],np.vstack((J@T,T*bs[:,None],-T*bs[:,None]))
    for it in range(8):
        v,K=evaluate(c);feas=float(v.min());row=dict(iteration=it,min_constraint=feas,elapsed=time.monotonic()-start)
        hist.append(row);atomic_json(dest/'history.json',hist)
        m.f.save(m.candidate(c),dest/'checkpoint.json',row);np.save(dest/'checkpoint.npy',c)
        print('REPAIR',json.dumps(row),flush=True)
        if feas>=-1e-7 or time.monotonic()-start>240:break
        scale=np.clip(1/np.maximum(np.linalg.norm(K,axis=1),1e-6),.001,1e4)
        A=K*scale[:,None];b=(v+1e-10)*scale
        r=minimize(lambda d:(.5*d@d,d),np.zeros(m.dim),jac=True,method='SLSQP',bounds=[(-.15,.15)]*m.dim,constraints=[dict(type='ineq',fun=lambda d:b+A@d,jac=lambda d:A)],options=dict(maxiter=70,ftol=1e-12,disp=False))
        row['qp_success']=bool(r.success);row['qp_message']=str(r.message);row['qp_iterations']=int(r.nit)
        improved=False
        for alpha in (1.,.5,.25,.125,.0625):
            new=c+alpha*(T@r.x)
            try:m.candidate(new)
            except ValueError:continue
            vv,_=evaluate(new)
            if float(vv.min())>feas:
                c=new;row['step_fraction']=alpha;improved=True;break
        if not improved:break
    pool,ng=pool_points();selection=[]
    for name,cc in [('round1',np.load(run/'round1/modifiers.npy')),('round2',np.load(run/'round2/modifiers.npy')),('repaired',c),('original_selection',np.load(run/'modifiers.npy'))]:
        v,K=evaluate(cc);good=float(v.min())>=-1e-7
        if good:
            try:raw=m.candidate(cc)
            except ValueError:good=False
        row=dict(id=name,feasible=good,min_constraint=float(v.min()))
        if good:
            peak=float(scan(m,cc,pool).max());row['training_pool_max']=peak
            if not selection or peak<min(x[0] for x in selection):selection.append((peak,name,cc.copy()))
        hist.append({'selection_check':row});print('SELECT',json.dumps(row),flush=True)
    if not selection:raise RuntimeError('No feasible field after bounded repair')
    peak,name,c=min(selection,key=lambda x:x[0]);old=json.loads((run/'summary.json').read_text())
    atomic_json(dest/'original_minimax_summary.json',old)
    meta=dict(selection=name,training_pool_max=peak,min_constraint=float(evaluate(c)[0].min()),elapsed=time.monotonic()-start,scope='Training-only selection with original constraint tolerance; recomputed epigraph value, no held-out data')
    m.f.save(m.candidate(c),run/'candidate.json',meta);np.save(run/'modifiers.npy',c)
    old['post_fit_feasibility_repair']=meta;old['selected_pool_peak']=peak;old['new_improvement_found']=peak<old['parent_pool_peak']
    import hashlib
    old['raw_sha256']=hashlib.sha256((run/'candidate.json').read_bytes()).hexdigest()
    atomic_json(run/'summary.json',old);atomic_json(dest/'summary.json',meta);atomic_json(dest/'history.json',hist)
    print('FEASIBILITY_REPAIR_DONE',json.dumps(meta),flush=True)
if __name__=='__main__':main()
