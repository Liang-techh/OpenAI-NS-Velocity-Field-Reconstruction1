import numpy as np
from scipy.io import loadmat
from openai_ns_reconstruction.eq45_export_bundle import export_bundle
from openai_ns_reconstruction.eq45_supported_delivery import default_field


def test_matlab_python_layout_and_exact_grid_agree(tmp_path):
    f=default_field().candidate
    m=export_bundle(f,tmp_path,5)
    a=np.load(tmp_path/'velocity.npz');b=loadmat(tmp_path/'velocity.mat')
    for key in ['u','v','w']:
        np.testing.assert_array_equal(a[key],b[key])
    index=(3,3,2,3)
    p=np.array([a['x'][index[1]],a['y'][index[2]],a['z'][index[3]]])
    expected=f.at_points(p,float(a['times'][index[0]]))
    np.testing.assert_allclose([a[key][index] for key in ['u','v','w']],expected,atol=1e-14)
    assert m['candidate_sha256']==f.sha256
    assert not m['truth_boundary']['pde_validated']
