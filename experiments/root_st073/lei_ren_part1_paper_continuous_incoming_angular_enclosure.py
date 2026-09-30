"""Incoming mixed angular row bounds with an enclosed sigma primitive."""
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_continuous_incoming_enclosure import ContinuousIncomingEnclosure


class ContinuousIncomingAngularEnclosure(ContinuousIncomingEnclosure):
    def unit_atoms(self,*,panels=4096):
        if not isinstance(panels,int) or panels<16:raise ValueError('Integer panels >=16 required')
        v=self.iv;step=v.mpf(1)/panels
        sigma=[self.sigma_node(j,panels) for j in range(panels+1)]
        primitive=v.mpf(0);total=v.mpf(0)
        for j in range(panels):
            after=primitive+step*self._hull(sigma[j],sigma[j+1])
            J=self._hull(primitive,after)
            y=self._hull(v.mpf(j)/panels,v.mpf(j+1)/panels)
            total+=step*v.exp(v.mpf('1.6')*y-v.mpf('.6')*J)
            primitive=after
        return dict(unit=total,terminal_J=primitive)

    def mixed_atoms(self,logPstar,Z,Md='.5',*,panels=4096):
        v=self.iv;z=v.mpf(Z)
        unit=self.unit_atoms(panels=panels)
        incoming=self.weighted_atoms(Md,panels=panels)
        factor=4*v.exp(v.mpf(logPstar))*(v.mpf(1)/v.mpf('1.6')+
            unit['unit']+v.exp(v.mpf('.3'))*incoming['transition'][0])
        return dict(factor=factor,I_theta_z=factor*z/(1+z**2),
            I_theta_z_Z=factor*(1-z**2)/(1+z**2)**2,
            unit=unit,incoming=incoming)

    def report(self,logPstar,Z,Md='.5',*,panels=4096,nominal_factor=None):
        atoms=self.mixed_atoms(logPstar,Z,Md,panels=panels)
        result=dict(logPstar=str(logPstar),Z=str(Z),Md=str(Md),panels=panels,
            interval_precision=self.precision,
            mixed_factor=self.describe(atoms['factor']),
            I_theta_z=self.describe(atoms['I_theta_z']),
            I_theta_z_Z=self.describe(atoms['I_theta_z_Z']),
            sigma_primitive_terminal=self.describe(atoms['unit']['terminal_J']),
            mixed_reference_quadrature_enclosed=True,
            parameter_scope='Declared logPstar, Z and Md; ideal continuous incoming angular source only',
            complete_angular_velocity_installed=False,inner_offsets_enclosed=False,
            global_mean_closed=False,finite_energy_certified=False)
        with mp.workdps(self.precision+20):
            lo,hi=self.bounds(atoms['factor'])
            result['relative_width']=mp.nstr((hi-lo)/lo,30)
            jlo,jhi=self.bounds(atoms['unit']['terminal_J'])
            result['exact_half_primitive_contained']=bool(jlo<=mp.mpf('.5')<=jhi)
            if nominal_factor is not None:
                result['nominal_factor_contained']=bool(lo<=mp.mpf(nominal_factor)<=hi)
        return result


def run():
    folder=Path(__file__).parent
    source=json.loads((folder/'lei_ren_part1_paper_shared_candidate_1.json').read_text())
    logP=source['shared_parameters']['logPstar']
    incoming=json.loads((folder/'lei_ren_part1_paper_continuous_incoming_outer.json').read_text())['incoming']
    nominal=incoming['continuous_incoming']['reference_linear_factors']['I_theta_z_times_1plusZ2_over_Z']
    enclosure=ContinuousIncomingAngularEnclosure(precision=80)
    reports=[]
    for panels in (1024,4096):
        r=enclosure.report(logP,str(incoming['Z']),panels=panels,nominal_factor=nominal)
        if not r['nominal_factor_contained'] or not r['exact_half_primitive_contained']:
            raise AssertionError('Mixed angular reference containment failed')
        reports.append(r)
    result=dict(reports=reports,full_angular_field_certified=False)
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps([dict(panels=r['panels'],relative_width=r['relative_width'],contained=r['nominal_factor_contained']) for r in reports]),flush=True)
    return result


if __name__=='__main__':run()
