"""Callable frozen reference for the bounded original-wave amplitude fit."""
import hashlib
import json
from pathlib import Path
import numpy as np
from global_two_patch_candidate import load as load_global
from supported_fourier_basis import basis_data

ROOT=Path(__file__).resolve().parent


class WaveReference:
    def __init__(self,base,snapshot,drift,amplitude,pressure_coefficients):
        self.base=base
        self.snapshot=snapshot
        self.drift=drift
        self.amplitude=float(amplitude)
        self.tau0=float(snapshot['inputs']['mean']['tau'])
        self.h=float(snapshot['inputs']['mean'].get('h',.005))
        self.nu=float(snapshot['inputs']['mean']['nu'])
        parent=json.loads((ROOT/'enriched_mean_endpoint_tangent.json').read_text())
        c=np.asarray(parent['selected']['coefficients_original'])
        self.wave_coefficients=c[:,0]+1j*c[:,1]
        self.pressure_coefficients=np.asarray(pressure_coefficients).reshape(2,45)

    def _increments(self,points):
        points=np.asarray(points,float)
        g=self.snapshot['inputs']['wave']
        carrier=np.asarray(g['carrier'])
        wave,_,_=basis_data(points,g['center'],g['widths'],1,2,carrier)
        du=(self.amplitude-1)*np.einsum('niq,q->ni',wave,self.wave_coefficients).real
        dp=np.zeros(len(points))
        for patch,name in enumerate(('first','second')):
            c=self.pressure_coefficients[patch]
            mode_coeff=(c[:9],c[9:27:2]+1j*c[10:27:2],c[27:45:2]+1j*c[28:45:2])
            for mode in (0,1,2):
                _,p,_=basis_data(points,self.drift['inputs'][name+'_patch_center'],
                    self.drift['inputs'][name+'_patch_widths'],mode,2,mode*carrier)
                dp+=np.einsum('nq,q->n',p,mode_coeff[mode]).real
        return du,dp

    def fields(self,points):
        u,p=self.base.fields(points,self.tau0)
        du,dp=self._increments(points)
        return u+du,p+dp

    def velocity(self,points):return self.fields(points)[0]
    def pressure(self,points):return self.fields(points)[1]


def load_reference(source=ROOT/'scale_wave_amplitude_fit.json'):
    report=json.loads(Path(source).read_text())
    if report['status']!='completed':raise ValueError('Complete amplitude fit first')
    for name,expected in report['sources'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=expected:
            raise ValueError('Changed amplitude-fit source: '+name)
    field,_,snapshot,drift,_=load_global(ROOT/'localized_drift_1800_fit.json')
    result=WaveReference(field,snapshot,drift,report['wave_amplitude'],report['pressure_coefficients'])
    projection=json.loads((ROOT/'scale_generator_pressure_projection.json').read_text())
    result.base_reference=WaveReference(field,snapshot,drift,1,projection['fit']['coefficients'])
    return result
