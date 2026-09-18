"""Publication tests: immutable evidence and API, not NS scientific acceptance."""
import hashlib,json,subprocess,sys
from pathlib import Path
import numpy as np
import pytest
from research_baseline import ROOT, load_best, verify_integrity
from research_baseline.visualization import (
    load_visualization_bundle,
    sample_reference_grid,
    write_visualization_bundle,
)


def test_integrity():
    r=verify_integrity();assert r['artifact_integrity'] and r['pde_validated'] is False


def test_evaluation_matches_original_runtime(tmp_path):
    f=load_best();rng=np.random.default_rng(9172999)
    x=rng.uniform(-2,2,(73,3));x[:3]=[[0,0,0],[0,0,.4],[0,0,-.4]];t=rng.uniform(.25,.75,len(x))
    out=tmp_path/'original.npz';data=tmp_path/'points.npz';np.savez(data,x=x,t=t)
    code="from spacetime import Family; import numpy as np; from pathlib import Path; f,r=Family.load(%r); d=np.load(%r);u,p=f.fields(r,d['x'],d['t']);np.savez(%r,u=u,p=p)"%(str(ROOT/'artifacts/research/ST006/candidate.json'),str(data),str(out))
    subprocess.run([sys.executable,'-c',code],cwd=ROOT/'experiments/root_st030',check=True)
    ref=np.load(out);u,p=f.fields(x,t)
    np.testing.assert_allclose(u,ref['u'],rtol=1e-11,atol=1e-11);np.testing.assert_allclose(p,ref['p'],rtol=1e-11,atol=1e-11)


def test_shapes_support_and_rotation():
    f=load_best();p=np.array([[2.,0,0],[0,0,2.],[3,3,3]])
    u,pr=f.fields(p,.5);assert not u.any() and not pr.any()
    x=np.array([[.1,.2,.3],[.7,0.,-.2]]);a=.43;Q=np.array([[np.cos(a),-np.sin(a),0],[np.sin(a),np.cos(a),0],[0,0,1]])
    np.testing.assert_allclose(f.at_points(x@Q.T,.5),f.at_points(x,.5)@Q.T,atol=1e-11)
    assert f.velocity(.1,0.,.1,np.array([.25,.5,.75])).shape==(3,3)
    assert not f._raw.flags.writeable


def test_visualization_bundle_roundtrip_and_tamper(tmp_path):
    f=load_best();grid=sample_reference_grid(5)
    assert grid['velocity'].shape==(3,5,5,5,3) and grid['speed'].shape==(3,5,5,5)
    assert grid['candidate_sha256']==f.sha256 and grid['pde_validated'] is False
    direct=f.at_points([grid['x'][3],grid['y'][2],grid['z'][1]],.5)
    # Grid evaluation is batched; the same point evaluated alone can differ at roundoff level.
    np.testing.assert_allclose(grid['velocity'][1,3,2,1],direct,rtol=2e-13,atol=2e-14)
    path=tmp_path/'st006-vis.npz';meta=write_visualization_bundle(path,5)
    loaded=load_visualization_bundle(path)
    assert meta['grid_sha256']==loaded['grid_sha256'] and loaded['visualization_ready'] is False
    np.testing.assert_array_equal(loaded['velocity'],grid['velocity'])
    assert not loaded['velocity'].flags.writeable and not loaded['times'].flags.writeable
    with np.load(path,allow_pickle=False) as data:
        payload={name:np.array(data[name],copy=True) for name in data.files}
    payload['velocity'][0,2,2,2,0]+=1e-3
    np.savez_compressed(path,**payload)
    with pytest.raises(ValueError,match='checksum'):
        load_visualization_bundle(path)


@pytest.mark.parametrize('time',[.249,.751,float('nan')])
def test_bad_time(time):
    with pytest.raises(ValueError):load_best().at_points([.1,0.,.1],time)


def test_tamper_rejected(tmp_path):
    import shutil
    for rel in ['research_baseline','artifacts/research/ST006']:shutil.copytree(ROOT/rel,tmp_path/rel)
    candidate=tmp_path/'artifacts/research/ST006/candidate.json';candidate.write_bytes(candidate.read_bytes()+b' ')
    with pytest.raises(ValueError,match='Checksum'):verify_integrity(tmp_path)


def test_original_evidence_stays_failed():
    f=ROOT/'artifacts/research/ST006/evidence/round2/ST006_validation.json';r=json.loads(f.read_text())
    assert r['pde_validated'] is False and r['gates']['momentum_max'] is False and r['gates']['momentum_L2'] is False
    assert r['gates']['divergence_max'] is False
    run=subprocess.run([sys.executable,str(ROOT/'experiments/root_st030/check_acceptance.py'),str(f)],capture_output=True,text=True)
    assert run.returncode==1 and 'momentum_max' in run.stdout
