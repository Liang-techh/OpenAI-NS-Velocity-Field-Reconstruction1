"""Physical mean-plus-wave replay with local oscillatory derivative fitting.

An affine physical-time wave coefficient is a local tangent, not a globally
integrated trajectory. Pressure harmonics and wave time derivative have zero
angular mean, so this step cannot conceal an uncorrected mean residual.
"""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from affine_momentum import jets, momentum
from supported_fourier_basis import basis_data
from broad_wave_covariance import load_saved_wave, ExactCurlWave
from broad_meridional_constrained import load_saved_field as load_constrained_mean
from broad_wave_mean_fit import patch_nodes, norms
from broad_shear_dynamic_control import load_saved_field as load_saved_field_mean
from grouped_joined_field import install_in_field
from meridional_state_cache import StatefulMean, value_and_pressure_modes


class MeanWaveField:
    def __init__(self, mean, wave, amplitude, rate, pressure, tau0):
        self.mean,self.wave=mean,wave
        self.amplitude,self.rate=float(amplitude),complex(rate)
        self.pressure=np.asarray(pressure,complex)
        self.tau0=float(tau0);self.nu=mean.nu

    def fields(self, points, tau):
        tau=float(np.asarray(tau).ravel()[0])
        u,p=self.mean.fields(points,tau)
        factor=self.amplitude+(self.tau0-tau)*self.rate
        u=u+(factor*self.wave.complex_velocity(points)).real
        p=p.copy()
        for m,coefficient in zip((1,2),self.pressure):
            P=basis_data(points,self.wave.center,self.wave.widths,m,self.wave.degree,
                         np.asarray(self.wave.carrier)*m)[1]
            p+=(P@coefficient).real
        return u,p


def grid(field,wave,tau,breaks,order,angles,shift):
    meridional,weights,metadata=patch_nodes(field,np.asarray(wave.center),
        np.asarray(wave.widths),tau,order,order,breaks)
    theta=shift+np.arange(angles)*2*np.pi/angles
    points=np.stack((meridional[:,0,None]*np.cos(theta),
                     meridional[:,0,None]*np.sin(theta),
                     np.broadcast_to(meridional[:,2,None],(len(meridional),angles))),axis=-1)
    return points.reshape(-1,3),np.repeat(weights/angles,angles),metadata


def columns(wave,points):
    w=wave.complex_velocity(points)
    blocks=[w.real[:,:,None],-w.imag[:,:,None]]
    for m in (1,2):
        G=basis_data(points,wave.center,wave.widths,m,wave.degree,
                     np.asarray(wave.carrier)*m)[2]
        blocks.extend((G.real,-G.imag))
    return np.concatenate(blocks,axis=2)


def load_saved_field(path=None):
    root=Path(__file__).resolve().parent
    report=json.loads(Path(path or root/'mean_wave_replay.json').read_text())
    source=json.loads((root/'meridional_constrained_evolution.json').read_text())
    dynamic,original=load_saved_field_mean();install_in_field(dynamic)
    base,values,pressures,_,_=value_and_pressure_modes(dynamic,original)
    initial=source['initial_fit']
    mean=StatefulMean(base,values,pressures,initial['state'],initial['slope'],
                      initial['pressure'],initial['k'])
    wave,_=load_saved_wave()
    if report.get('wave_report'):
        wave=codesigned_wave(root/report['wave_report'])
        mean,_=load_constrained_mean();install_in_field(mean)
    raw=np.asarray(report['pressure_coefficients']);pressure=raw[...,0]+1j*raw[...,1]
    raw=np.asarray(report['physical_time_rate']);rate=raw[0]+1j*raw[1]
    return MeanWaveField(mean,wave,report['amplitude'],rate,pressure,report['tau']),report


def codesigned_wave(path):
    data=json.loads(Path(path).read_text())
    raw=np.asarray(data['selected']['coefficients_original'])
    return ExactCurlWave(data['center'],data['widths'],data['degree'],
                         data['carrier'],raw[...,0]+1j*raw[...,1])


def run(fit_order=6,holdout_order=9,output_name='mean_wave_replay.json',wave_report=None):
    root=Path(__file__).resolve().parent
    source=json.loads((root/'meridional_constrained_evolution.json').read_text())
    scalar=json.loads((root/'broad_wave_mean_fit.json').read_text())
    saved=json.loads((root/'broad_meridional_constrained.json').read_text())
    dynamic,original=load_saved_field_mean();install_in_field(dynamic)
    base,values,pressures,_,_=value_and_pressure_modes(dynamic,original)
    initial=source['initial_fit'];tau=.5*2.**(-initial['k'])
    mean=StatefulMean(base,values,pressures,initial['state'],initial['slope'],
                      initial['pressure'],initial['k'])
    wave,_=load_saved_wave();q=(wave.degree+1)**2
    amplitude=scalar['rows'][0]['amplitude_fit']
    if wave_report:
        wave=codesigned_wave(root/wave_report)
        mean,_=load_constrained_mean();install_in_field(mean)
        amplitude=1.
    frozen=MeanWaveField(mean,wave,amplitude,0,np.zeros((2,q)),tau)
    breaks=sorted(set(saved['quadrature']['radial_split_breaks']+[.62]))
    h=5e-4*np.sqrt(mean.nu*tau);ht=1e-4*tau
    started=time.perf_counter()
    pts,weights,_=grid(mean,wave,tau,breaks,fit_order,8,.13)
    residual=momentum(jets(frozen,pts,tau,h,ht))
    C=columns(wave,pts);A=(C*np.sqrt(weights[:,None,None])).reshape(-1,C.shape[-1])
    b=-(residual*np.sqrt(weights[:,None])).reshape(-1)
    scales=np.maximum(np.linalg.norm(A,axis=0),1e-30)
    z,_,rank,s=np.linalg.lstsq(A/scales,b,rcond=1e-10);control=z/scales
    pressure=np.stack((control[2:2+q]+1j*control[2+q:2+2*q],
                       control[2+2*q:2+3*q]+1j*control[2+3*q:2+4*q]))
    field=MeanWaveField(mean,wave,amplitude,control[0]+1j*control[1],pressure,tau)
    pack=lambda x:np.stack((np.asarray(x).real,np.asarray(x).imag),axis=-1).tolist()
    result=dict(status='fit_complete',accepted=False,pde_validated=False,scale_recursion_established=False,
        amplitude=amplitude,physical_time_rate=pack(field.rate),pressure_coefficients=pack(pressure),
        wave_report=wave_report,
        mean_report=('broad_meridional_constrained.json' if wave_report else 'meridional_constrained_evolution.json:initial_fit'),
        k=initial['k'],tau=tau,rank=int(rank),column_count=C.shape[-1],
        fit_order=fit_order,holdout_order=holdout_order,
        velocity_pressure_normalized_inner_product_max=float(np.max(abs(
            (A[:,:2]/scales[:2]).T@(A[:,2:]/scales[2:])))),
        fit_points=len(pts),frozen_fit=norms(residual,weights),
        fitted_fit=norms(residual+np.einsum('nip,p->ni',C,control),weights),
        scope='Actual physical mean plus exact-curl mode1 wave; local affine physical-time coefficient, fitted zero-mean mode1/2 pressures. Independent Cartesian full-momentum holdout on fixed patch. Not time integration, whole-support or recursion acceptance.')
    output=root/output_name;output.write_text(json.dumps(result,indent=2)+'\n')
    hp,hw,metadata=grid(mean,wave,tau,breaks,holdout_order,12,.31)
    result['holdout_geometry']=metadata;result['holdout_points']=len(hp)
    result['mean_holdout']=norms(momentum(jets(mean,hp,tau,h,ht)),hw)
    result['wave_holdout']=norms(momentum(jets(field,hp,tau,h,ht)),hw)
    result['status']='completed';result['elapsed_seconds']=time.perf_counter()-started
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--fit-order',type=int,default=6)
    parser.add_argument('--holdout-order',type=int,default=9)
    parser.add_argument('--output',default='mean_wave_replay.json')
    parser.add_argument('--wave-report',default=None)
    args=parser.parse_args()
    run(args.fit_order,args.holdout_order,args.output,args.wave_report)
