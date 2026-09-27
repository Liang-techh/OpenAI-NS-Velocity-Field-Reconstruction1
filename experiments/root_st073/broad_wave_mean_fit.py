"""Fit nonnegative energy of the saved growing wave to complete mean force.

Only the angular mean is optimized; oscillatory momentum remains uncorrected.
The result is a candidate for joint mean-wave work, never a full PDE solution.
"""
import json
from pathlib import Path
import time
import numpy as np
from numpy.polynomial.legendre import leggauss
from affine_momentum import jets, momentum
from joined_field import coordinates
from broad_shear_dynamic_control import load_saved_field
from grouped_joined_field import install_in_field
from meridional_state_cache import StatefulMean, value_and_pressure_modes
from wave_mean_flux import mean_force


def norms(residual, weights):
    return dict(sampled_max=float(np.linalg.norm(residual,axis=1).max()),
                volume_L2=float(np.sqrt(np.sum(weights[:,None]*residual**2))),
                volume_RMS=float(np.sqrt(np.sum(weights[:,None]*residual**2)/sum(weights))))


def patch_nodes(field, center, widths, tau, nz, nr, breaks):
    gz,wz=leggauss(nz);gr,wr=leggauss(nr);points=[];weights=[]
    rlo,rhi=center[0]+np.array([-1,1])*widths[0]
    for z,zw in zip(center[1]+widths[1]*gz,widths[1]*wz):
        q=float(coordinates(0.,z/np.sqrt(field.nu),tau,field.inner.h)['q'])
        ri=np.sqrt(2*field.nu*q*field.join_X)
        cuts=sorted(set([rlo,rhi]+[float(r) for r in ri*(1+(field.ratio-1)*np.asarray(breaks)) if rlo<r<rhi]))
        for lo,hi in zip(cuts[:-1],cuts[1:]):
            for r,rw in zip((lo+hi)/2+(hi-lo)*gr/2,(hi-lo)*wr/2):
                points.append([r,0.,z]);weights.append(2*np.pi*r*rw*zw)
    return np.asarray(points),np.asarray(weights),dict(z_order=nz,panel_order=nr,
        radial_breaks=list(breaks),point_count=len(points))


def run():
    root=Path(__file__).resolve().parent
    wave=json.loads((root/'broad_shear_growth.json').read_text())
    evolution=json.loads((root/'meridional_constrained_evolution.json').read_text())
    dynamic,source=load_saved_field(); install_in_field(dynamic)
    base,values,pressures,_,_=value_and_pressure_modes(dynamic,source)
    fit=evolution['initial_fit']
    field=StatefulMean(base,values,pressures,fit['state'],fit['slope'],fit['pressure'],fit['k'])
    center=np.asarray(wave['center']); widths=np.asarray(wave['widths'])
    mode=int(wave['selected_candidate']['mode']);carrier=wave['carriers'][str(mode)]
    raw=np.asarray(wave['selected_candidate']['potential_coefficients'])
    coefficient=raw[...,0]+1j*raw[...,1]
    tau=.5*2.**(-fit['k']);h=5e-4*np.sqrt(field.nu*tau)
    volume=np.pi*((center[0]+widths[0])**2-(center[0]-widths[0])**2)*2*widths[1]
    rows=[]; energy=None; started=time.perf_counter()
    saved=json.loads((root/'broad_meridional_constrained.json').read_text())
    breaks=sorted(set(saved['quadrature']['radial_split_breaks']+[.62]))
    for label,nz,nr in [('fit',6,6),('holdout',9,9)]:
        points,weights,geometry=patch_nodes(field,center,widths,tau,nz,nr,breaks)
        np.testing.assert_allclose(sum(weights),volume,rtol=1e-12)
        residual=momentum(jets(field,points,tau,h,1e-4*tau))
        force=mean_force(points[:,[0,2]],center,widths,mode,wave['degree'],carrier,
                         coefficient,field.nu,h,angles=16)
        linear=float(np.sum(weights[:,None]*residual*force))
        quadratic=float(np.sum(weights[:,None]*force*force))
        if energy is None:
            energy=max(0.,-linear/max(quadratic,1e-300))
        rows.append(dict(label=label,geometry=geometry,domain_volume=float(volume),
            base=norms(residual,weights),with_wave_mean=norms(residual+energy*force,weights),
            energy_fit=energy,amplitude_fit=float(np.sqrt(energy)),
            residual_force_inner_product=linear,unit_force_squared_norm=quadratic,
            unconstrained_energy_optimum=-linear/max(quadratic,1e-300)))
    result=dict(accepted=False,pde_validated=False,scale_recursion_established=False,
        wave_integrated=False,oscillatory_residual_fitted=False,
        mean_source='meridional_constrained_evolution.json:initial_fit',
        wave_source='broad_shear_growth.json:selected_candidate',
        physical_domain=dict(radial_bounds=(center[0]+np.array([-1,1])*widths[0]).tolist(),
                             axial_bounds=(center[1]+np.array([-1,1])*widths[1]).tolist(),
                             angular_bounds=[0.,2*np.pi]),
        k=fit['k'],tau=tau,rows=rows,elapsed_seconds=time.perf_counter()-started,
        scope='Nonnegative scalar energy fit of the complete exact-curl wave mean force on its fixed physical patch. Independent spatial holdout; no oscillatory residual, joint mean refit, PDE or recursion acceptance.')
    root.joinpath('broad_wave_mean_fit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)


if __name__=='__main__':run()
