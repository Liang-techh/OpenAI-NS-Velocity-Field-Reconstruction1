"""Export the immutable ST054 mathematical coefficients for native MATLAB.
Usage: python export_st054.py --source-root SOURCE_CHECKOUT --out data/st054_models.mat
No fitting. Normalization and original QR transforms are baked into the coefficient
matrices; MATLAB uses the same separable bump/Legendre basis and continuous time.
"""
from __future__ import annotations
import argparse, hashlib, json, sys
from pathlib import Path
import numpy as np
from scipy.io import savemat, loadmat
SOURCE_SHA='c77492a48e9c0f13d4d51244987c57c28519409b'
RAW_SHA={'ST054-M3':'2b2b571986297b52bc884a2ee4808a7ade1ade86f38f938604f2a64118882df3',
         'ST054-Q2':'d772949621e0f7703a9d6b36f28828532ba3e5ec678dec14db5ce14eb8b851d5'}

def export(source_root: Path,out: Path):
    root=source_root.resolve()
    sys.path.insert(0,str(root/'experiments/root_st054'))
    sys.path.insert(0,str(root/'experiments/root_st030'))
    from replay_st054 import reconstruct
    from spacetime import Family
    rng=np.random.default_rng(9186801)
    x=rng.uniform(-2.05,2.05,(96,3))
    x[:8]=[[0,0,0],[0,0,.16],[.2,0,0],[2,0,0],[0,0,2],[0,0,-2],[.03,0,1.93],[.5,.5,0]]
    times=np.array([.25,.33173,.5,.63791,.75])
    models=[];report=[]
    for ident in ['ST054-Q2','ST054-M3']:
        f,raw=reconstruct(ident)
        rawfile=root/f'artifacts/research/{ident}/candidate.json'
        if rawfile.exists():
            assert hashlib.sha256(rawfile.read_bytes()).hexdigest()==RAW_SHA[ident]
            f0,original=Family.load(rawfile)
            np.testing.assert_allclose(raw,original,rtol=0,atol=1e-11)
        a,b,p,fc,_=f.coefficients(raw)
        F=(f.Tp@a.reshape(f.ns,f.nt)).reshape(f.nr,f.nz,f.nt)
        G=(f.Tw@b.reshape(f.ns,f.nt)).reshape(f.nr,f.nz,f.nt)
        P=(f.Tq@p.reshape(f.ns,f.nt)).reshape(f.nr,f.nz,f.nt)
        vel=[];pres=[];res=[]
        for t in times:
            u,q=f.fields(raw,x,t)
            rr=f.analytic_residual(raw,x,t)
            vel.append(u);pres.append(q);res.append(rr)
        # Small reference data are independent of MATLAB implementation.
        model=dict(id=ident,mode='spectral',basis_kind=f.basis_kind,F=F,G=G,P=P,
                   force=fc,nu=.01,tmin=.25,tmax=.75,support_radius=2.,support_z=2.,
                   source_commit=SOURCE_SHA,original_raw_sha256=RAW_SHA[ident],
                   raw_coefficients=raw,pde_validated=np.uint8(0),
                   ref_points=x,ref_times=times,ref_velocity=np.array(vel),
                   ref_pressure=np.array(pres),ref_residual=np.array(res))
        models.append(model)
        report.append({'id':ident,'raw_array_sha256':hashlib.sha256(raw.astype('<f8').tobytes()).hexdigest(),
                       'original_raw_sha256':RAW_SHA[ident],'reference_points':len(x),'times':times.tolist()})
    result=root/'experiments/root_st054/results.json'
    evidence=json.loads(result.read_text())
    out.parent.mkdir(parents=True,exist_ok=True)
    if out.exists():raise FileExistsError('Refusing to overwrite '+str(out))
    savemat(out,dict(schema='ns_matlab_spectral_v1',models=np.array(models,dtype=object),
                    source_commit=SOURCE_SHA,evidence_json=json.dumps(evidence),
                    pde_validated=np.uint8(0)),do_compression=True,oned_as='row')
    loaded=loadmat(out,simplify_cells=True)
    for a,b in zip(models,loaded['models']):
        for key in ['F','G','P','raw_coefficients','ref_velocity','ref_pressure','ref_residual']:
            assert np.array_equal(a[key],b[key]),key
    receipt={'source_commit':SOURCE_SHA,'schema':'ns_matlab_spectral_v1',
             'mat_sha256':hashlib.sha256(out.read_bytes()).hexdigest(),
             'models':report,'roundtrip_exact':True,'native_matlab_executed':False,
             'scope':'Coefficient transfer, not a new fit or PDE validation. MATLAB tests run separately.'}
    out.with_suffix('.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))
    return receipt
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source-root',type=Path,required=True)
    p.add_argument('--out',type=Path,default=Path(__file__).parent/'data/st054_models.mat')
    a=p.parse_args();export(a.source_root,a.out)
