"""Actual-field replay for the 236-control, degree-3 mode-2 candidate."""
import argparse
import hashlib
import json
import time
from pathlib import Path
import numpy as np
from affine_momentum import jets, momentum
from broad_meridional_constrained import load_saved_field
from broad_shear_dynamic_control import load_saved_field as load_dynamic
from full_wave_tangent import LocalPotentialField, _unpack_full
from grouped_joined_field import install_in_field
from meridional_state_cache import value_and_pressure_modes
from vortex_state_observables import fixed_cylinder, observe
from wave_higher_harmonic_tangent import Mode2TangentCorrection, _decode_mode_block
from wave_residual_harmonics import budget
from supported_fourier_basis import basis_data

ROOT = Path(__file__).resolve().parent


class Mode0TangentCorrection:
    def __init__(self,base,geometry,control,tau0):
        self.base=base
        self.geometry=geometry
        self.control=np.asarray(control,float)
        self.tau0=tau0
        if len(self.control)!=64:
            raise ValueError('Degree-3 mode0 needs 48 velocity and 16 pressure controls')

    def fields(self,points,tau):
        velocity,pressure=self.base.fields(points,tau)
        V,P,_=basis_data(points,self.geometry['center'],self.geometry['widths'],0,3,np.zeros(2))
        delta=self.tau0-float(np.asarray(tau).ravel()[0])
        return (velocity+delta*np.einsum('niq,q->ni',V,self.control[:48]).real,
                pressure+np.einsum('nq,q->n',P,self.control[48:]).real)


def build_field(candidate):
    snapshot = json.loads((ROOT/'full_wave_frozen_cache.json').read_text())
    geometry = snapshot['inputs']['wave']
    control = np.array(candidate['selected']['tangent_coefficients'])
    if len(control) not in (236,264) or geometry['degree']!=2:
        raise ValueError('Expected 236 or 264 controls with degree-2 initial wave')
    mean, report = load_saved_field()
    install_in_field(mean)
    if report['coefficients']!=snapshot['inputs']['mean']['coefficients']:
        raise ValueError('Frozen mean changed')
    c = np.array(candidate['selected']['coefficients_original'])
    wave = c[:,0]+1j*c[:,1]
    legacy = (np.r_[control[:108],np.zeros(72)] if len(control)==236
              else np.r_[np.zeros(36),control[64:136],np.zeros(72)])
    derivatives, pressures = _unpack_full(legacy,9)
    carrier = np.array(geometry['carrier'])
    tau = snapshot['inputs']['mean']['tau']
    base = LocalPotentialField(mean,geometry['center'],geometry['widths'],2,
        {0:np.zeros(2),1:carrier,2:2*carrier},
        (np.zeros(27,complex),wave,np.zeros(27,complex)),derivatives,pressures,tau)
    if len(control)==264:
        base=Mode0TangentCorrection(base,geometry,control[:64],tau)
    dv,dp = _decode_mode_block(control[-128:],3)
    field = Mode2TangentCorrection(base,geometry['center'],geometry['widths'],carrier,3,dv,dp,tau)
    field.nu = mean.nu
    return field,snapshot


def run(source,output,mode,delta_k=1e-6):
    started = time.perf_counter()
    raw = Path(source).read_bytes()
    candidate = json.loads(raw)
    if candidate['status']!='completed' or not candidate['assembled_feasible']:
        raise ValueError('Finish a feasible frozen candidate before replay')
    field,snapshot = build_field(candidate)
    tau = snapshot['inputs']['mean']['tau']
    report = dict(status='running',accepted=False,pde_validated=False,
        scale_recursion_established=False,source=Path(source).name,
        source_sha256=hashlib.sha256(raw).hexdigest(),mode=mode,
        tangent_control_count=len(candidate['selected']['tangent_coefficients']),scope='Sampled actual Cartesian field. No full-domain, time-evolution or scale-recursion acceptance.')
    def save():
        Path(output).write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    save()
    if mode=='momentum':
        dense = json.loads((ROOT/'full_wave_dense_tangent.json').read_text())
        data = dense['new_frozen_cache']
        points,weights = np.array(data['points']),np.array(data['weights'])
        residual = momentum(jets(field,points,tau,snapshot['timesteps']['hspace'],snapshot['timesteps']['htime']))
        report['momentum'] = budget(points,weights,residual,12)
    else:
        dynamic,dynamic_report = load_dynamic()
        install_in_field(dynamic)
        base,*_ = value_and_pressure_modes(dynamic,dynamic_report)
        k0 = dynamic_report['k']
        points,weights,domain = fixed_cylinder(base.inner,k0,angles=12)
        report.update(domain=domain,delta_k=delta_k,rows=[])
        for k in (k0,k0+delta_k):
            row = observe(field,points,weights,k)
            report['rows'].append(row)
            save()
            print(json.dumps(row),flush=True)
        left,right = report['rows']
        report['forward_interval_direction_flags'] = dict(
            radial_contraction=right['enstrophy_radial_rms']<left['enstrophy_radial_rms'],
            relative_axial_elongation=right['enstrophy_aspect_ratio']>left['enstrophy_aspect_ratio'],
            angular_speed_magnitude_increase=abs(right['enstrophy_weighted_angular_speed'])>abs(left['enstrophy_weighted_angular_speed']))
    report.update(status='completed',elapsed_seconds=time.perf_counter()-started)
    save()
    print(json.dumps(report.get('momentum',report.get('forward_interval_direction_flags'))),flush=True)
    return report


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--source',type=Path,default=ROOT/'enriched_shape_tangent.json')
    parser.add_argument('--mode',choices=['momentum','shape'],default='momentum')
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    run(args.source,args.output or ROOT/f'enriched_shape_{args.mode}_replay.json',args.mode)
