"""Complete mean forcing from a compact exact-curl Fourier wave.

Physical Cartesian wave advection is averaged in the cylindrical frame and
compared with cylindrical divergence of the complete covariance tensor.
No principal-wave or radial-flux-only approximation is used.
"""
import json
from pathlib import Path
import numpy as np
from supported_fourier_basis import basis_data
from fourier_patch_evolution import basis_jets


def angular_points(rz, count):
    rz = np.asarray(rz, float)
    theta = .173 + np.arange(count)*2*np.pi/count
    c, s = np.cos(theta), np.sin(theta)
    p = np.stack((rz[:, 0, None]*c, rz[:, 0, None]*s,
                  np.broadcast_to(rz[:, 1, None], (len(rz), count))), axis=-1)
    return p.reshape(-1, 3), c, s


def cylindrical(values, c, s):
    v = np.asarray(values).reshape(-1, len(c), 3)
    return np.stack((c*v[..., 0]+s*v[..., 1],
                     -s*v[..., 0]+c*v[..., 1], v[..., 2]), axis=-1)


def covariance(rz, center, widths, mode, degree, carrier, coefficients, angles=16):
    points, c, s = angular_points(rz, angles)
    V = basis_data(points, center, widths, mode, degree, carrier)[0]
    w = cylindrical(np.einsum('niq,q->ni', V, coefficients).real, c, s)
    return np.einsum('nai,naj->nij', w, w)/angles


def mean_force(rz, center, widths, mode, degree, carrier, coefficients, nu, h, angles=16):
    points, c, s = angular_points(rz, angles)
    V, J, _, _, _ = basis_jets(points, center, widths, mode, degree, carrier, nu, h)
    w = np.einsum('niq,q->ni', V, coefficients).real
    grad = np.einsum('nijq,q->nij', J, coefficients).real
    adv = np.einsum('nij,nj->ni', grad, w)
    return cylindrical(adv, c, s).mean(axis=1)


def single_mode_force(rz, center, widths, mode, degree, carrier, coefficients, nu, h):
    """Analytic angular average for one nonzero harmonic, including full curl.

At theta=0 the cylindrical frame equals Cartesian. Rotational covariance of
the complete velocity/Jacobian yields <Re(J e^imtheta) Re(w e^imtheta)>.
Different harmonics or mode zero require their own cross-term treatment.
"""
    if int(mode) != mode or mode <= 0:
        raise ValueError('Single-mode average requires a positive integer harmonic')
    rz=np.asarray(rz,float)
    points=np.column_stack((rz[:,0],np.zeros(len(rz)),rz[:,1]))
    V,J,_,_,_=basis_jets(points,center,widths,mode,degree,carrier,nu,h)
    w=np.einsum('niq,q->ni',V,coefficients)
    grad=np.einsum('nijq,q->nij',J,coefficients)
    return .5*np.einsum('nij,nj->ni',grad,w.conj()).real


def covariance_divergence(rz, center, widths, mode, degree, carrier, coefficients, h, angles=16):
    args = (center, widths, mode, degree, carrier, coefficients, angles)
    W = covariance(rz, *args)
    derivatives = []
    for axis in np.eye(2):
        wm2, wm, wp, wp2 = [covariance(rz+j*h*axis, *args) for j in (-2,-1,1,2)]
        derivatives.append((wm2-8*wm+8*wp-wp2)/(12*h))
    dr, dz = derivatives; r = np.asarray(rz)[:, 0]
    return np.column_stack((dr[:, 0, 0]+dz[:, 2, 0]+(W[:, 0, 0]-W[:, 1, 1])/r,
                            dr[:, 0, 1]+dz[:, 2, 1]+2*W[:, 0, 1]/r,
                            dr[:, 0, 2]+dz[:, 2, 2]+W[:, 0, 2]/r))


def run():
    root = Path(__file__).resolve().parent
    saved = json.loads((root/'broad_shear_growth.json').read_text())
    raw = np.asarray(saved['selected_candidate']['potential_coefficients'])
    coefficient = raw[..., 0]+1j*raw[..., 1]
    center = np.asarray(saved['center']); widths = np.asarray(saved['widths'])
    mode = int(saved['selected_candidate']['mode']); carrier = saved['carriers'][str(mode)]
    rz = center + widths*np.array([[0.,0.],[-.35,.15],[.4,-.25],[.8,.6]])
    tau = .5*2.**(-saved['initial_k']); nu=.01; h=5e-4*np.sqrt(nu*tau)
    args=(center,widths,mode,saved['degree'],carrier,coefficient)
    direct=mean_force(rz,*args,nu,h)
    analytic=single_mode_force(rz,*args,nu,h)
    tensor=covariance(rz,*args)
    projected=covariance_divergence(rz,*args,h)
    coarse=covariance(rz,*args,angles=8)
    report=dict(accepted=False,pde_validated=False,scale_recursion_established=False,
        source='broad_shear_growth.json', mode=mode, nu=nu, tau=tau, h=h,
        coefficient_normalization='Unchanged saved selected potential coefficient; no renormalization.',
        rz=rz.tolist(), covariance=tensor.tolist(), mean_advection=direct.tolist(),
        covariance_divergence=projected.tolist(),
        maximum_absolute_identity_error=float(np.max(abs(direct-projected))),
        single_mode_average_max_absolute_error=float(np.max(abs(direct-analytic))),
        maximum_scaled_identity_error=float(np.max(np.linalg.norm(direct-projected,axis=1)/np.maximum(np.linalg.norm(direct,axis=1),1.))),
        angular_covariance_8_vs_16_max=float(np.max(abs(coarse-tensor))),
        scope='Four-point complete actual-wave mean-advection identity. Every covariance component retained, including curl/cutoff terms. Not a matched mean-wave solution or time integration.')
    root.joinpath('wave_mean_flux.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('rz','covariance','mean_advection','covariance_divergence')}),flush=True)


if __name__=='__main__':run()
