"""Actual finite-window particle trajectories; not a source-identity certificate."""
from pathlib import Path
import argparse,hashlib
import numpy as np
import minimax_exchange
from spacetime import Family
from structure_and_particles import integrate
from minimax_exchange import atomic_json

def compute(candidate,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);f,raw=Family.load(candidate)
    R,Z=np.meshgrid([.04,.09,.14,.19],[-.21,-.14,-.075,-.035,0.,.035,.075,.14,.21],indexing='ij')
    initial=np.c_[R.ravel()*np.sqrt(.75),np.zeros(R.size),Z.ravel()*.75**.495]
    field=lambda x,t:f.fields(raw,x,t)[0]
    t,p,n=integrate(field,initial,1e-8,.02);_,fine,nfine=integrate(field,initial,1e-10,.01)
    rf=np.linalg.norm(fine[-1,:,:2],axis=1);r0=np.linalg.norm(initial[:,:2],axis=1)
    phi=np.unwrap(np.arctan2(fine[:,:,1],fine[:,:,0]),axis=0);turn=phi[-1]-phi[0]
    mask=initial[:,2]!=0;outward=np.sign(initial[mask,2])*(fine[-1,mask,2]-initial[mask,2])
    result=dict(candidate_sha256=hashlib.sha256(Path(candidate).read_bytes()).hexdigest(),seeds=len(initial),radial_contraction_fraction=float(np.mean(rf<r0)),positive_turn_fraction=float(np.mean(turn>0)),bipolar_outward_displacement_fraction=float(np.mean(outward>0)),midplane_upward_fraction=float(np.mean(fine[-1,~mask,2]>0)),radial_ratio_range=[float((rf/r0).min()),float((rf/r0).max())],turn_radians_range=[float(turn.min()),float(turn.max())],integration_difference=float(np.max(abs(p-fine))),function_evaluations=[n,nfine],scope='36 finite-window nonautonomous trajectories. Two DOP853 configurations, not independent integrators or a proof for all trajectories. No steady-streamline/source-identity assertion.')
    atomic_json(out/'summary.json',result);np.savez_compressed(out/'paths.npz',time=t,initial=initial,positions=fine);print(result,flush=True);return result
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('candidate');p.add_argument('--out',required=True);a=p.parse_args();compute(a.candidate,a.out)
