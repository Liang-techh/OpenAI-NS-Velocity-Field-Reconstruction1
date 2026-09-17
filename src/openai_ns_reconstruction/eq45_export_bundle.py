"""Portable candidate JSON plus matching Python/MATLAB velocity samples."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.io import savemat
from .eq45_supported_delivery import default_field


def export_bundle(candidate,output,grid_size=17):
    if not isinstance(grid_size,int) or not 2<=grid_size<=65:
        raise ValueError('grid_size must be an integer in [2,65]')
    out=Path(output);out.mkdir(parents=True,exist_ok=True)
    path=out/'candidate.json'
    candidate.save_json(path)
    replay=type(candidate).load_json(path)
    if replay.sha256!=candidate.sha256:
        raise RuntimeError('candidate identity changed on reload')
    x=np.linspace(-2,2,grid_size);times=np.linspace(.25,.75,9)
    values=replay.grid(x,x,x,times)
    if not np.isfinite(values).all() or not np.any(values):
        raise RuntimeError('export must be finite and nonzero')
    data=dict(x=x,y=x,z=x,times=times,u=values[...,0],v=values[...,1],w=values[...,2])
    np.savez_compressed(out/'velocity.npz',**data)
    savemat(out/'velocity.mat',data,do_compression=True,oned_as='column')
    (out/'evaluate_velocity.m').write_text("""% Interpolated sampled candidate; not exact off-grid evaluation.
S = load('velocity.mat');
Fu = griddedInterpolant({S.times,S.x,S.y,S.z},S.u,'linear','none');
Fv = griddedInterpolant({S.times,S.x,S.y,S.z},S.v,'linear','none');
Fw = griddedInterpolant({S.times,S.x,S.y,S.z},S.w,'linear','none');
velocity = @(x,y,z,t) [Fu(t,x,y,z),Fv(t,x,y,z),Fw(t,x,y,z)];
disp(velocity(0.1,0,0.1,0.5));
""",encoding='utf-8')
    point=[.1,0,.1];time=.5
    manifest=dict(candidate_sha256=replay.sha256,family=replay.to_dict()['schema'],
                  shape=list(values.shape),layout=['time','x','y','z','component'],
                  components=['u','v','w'],time_interval=[.25,.75],box=[[-2,2]]*3,
                  exact_python_sample=dict(point=point,time=time,velocity=replay.at_points(np.array(point),time).tolist()),
                  matlab='linear interpolation of saved grid; differs from exact candidate off-grid',
                  matlab_runtime_tested=False,
                  truth_boundary=replay.to_dict()['truth_boundary'])
    manifest['files']={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                       for p in [path,out/'velocity.npz',out/'velocity.mat',out/'evaluate_velocity.m']}
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    return manifest


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--output',default='artifacts/delivery/eq45_supported')
    p.add_argument('--grid-size',type=int,default=17)
    args=p.parse_args()
    print(json.dumps(export_bundle(default_field().candidate,args.output,args.grid_size),indent=2))


if __name__=='__main__':
    main()
