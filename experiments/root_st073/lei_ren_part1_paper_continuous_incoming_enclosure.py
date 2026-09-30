"""Outward bounds for the complete continuous incoming axial rows.

The finite transition y=exp(Md*t), t in [0,1], includes all cutoff support.
Monotone switch node bounds avoid singular interval endpoint divisions.
Angular and inner-seed errors are outside this reference axial calculation.
"""
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_continuous_pulse_energy_enclosure import ContinuousPulseEnergyEnclosure


class ContinuousIncomingEnclosure(ContinuousPulseEnergyEnclosure):
    def weighted_atoms(self,Md='.5',*,power=1,panels=4096):
        if not isinstance(panels,int) or panels<16:
            raise ValueError('Integer panels >=16 required')
        v=self.iv;Md=v.mpf(Md);power=v.mpf(power)
        if self.bounds(Md)[0]<=0:raise ValueError('Positive Md required')
        if self.bounds(power)[0]<=0:raise ValueError('Positive weight power required')
        step=v.mpf(1)/panels
        sigma=[self.sigma_node(j,panels) for j in range(panels+1)]
        transition=[v.mpf(0),v.mpf(0)]
        for j in range(panels):
            t=self._hull(v.mpf(j)/panels,v.mpf(j+1)/panels)
            y=v.exp(Md*t)
            c=self._hull(sigma[panels-j-1],sigma[panels-j])
            weight=Md*y*v.exp(power*y)*step
            transition[0]+=weight*c
            transition[1]+=weight*c**2
        plateau=(v.exp(power)-1)/power
        return dict(plateau=plateau,transition=transition,
                    weighted=[plateau+t for t in transition])

    def row_atoms(self,Z,Md='.5',*,panels=4096):
        v=self.iv;z=v.mpf(Z)
        atoms=self.weighted_atoms(Md,panels=panels)
        mass_factor=4*(1+atoms['weighted'][0])
        squared_factor=16*(1+atoms['weighted'][1])
        return dict(I_z=z*mass_factor,I_uz2=z**2*squared_factor,
                    I_z_Z=mass_factor,I_uz2_Z=2*z*squared_factor,
                    weighted_atoms=atoms)

    def report(self,Z,Md='.5',*,panels=4096,nominal=None):
        atoms=self.row_atoms(Z,Md,panels=panels)
        keys=('I_z','I_uz2','I_z_Z','I_uz2_Z')
        report=dict(Z=str(Z),Md=str(Md),panels=panels,
            interval_precision=self.precision,
            rows={key:self.describe(atoms[key]) for key in keys},
            reference_axial_quadrature_enclosed=True,
            full_support_covered=True,omitted_support_bound='0',
            parameter_scope='Declared real decimal Z and Md; angular, source schedule and inner-seed uncertainty are not enclosed',
            mixed_angular_row_enclosed=False,inner_offsets_enclosed=False,
            global_mean_closed=False,finite_energy_certified=False,
            scale_recursion_certified=False)
        with mp.workdps(self.precision+20):
            widths={}
            for key in ('I_z','I_uz2'):
                lo,hi=self.bounds(atoms[key]);scale=min(abs(lo),abs(hi))
                widths[key]=mp.nstr((hi-lo)/scale,30) if scale else None
            report['relative_widths']=widths
            if nominal is not None:
                values=nominal.full_incoming_rows(Z)
                report['nominal_containment']={key:bool(self.bounds(atoms[key])[0]<=values[key]<=self.bounds(atoms[key])[1]) for key in values}
        return report


def run():
    from lei_ren_part1_paper_continuous_incoming import ContinuousIncomingAxial
    enclosure=ContinuousIncomingEnclosure(precision=80)
    nominal=ContinuousIncomingAxial('.5',precision=100)
    reports=[]
    for panels in (1024,4096):
        report=enclosure.report('.3',panels=panels,nominal=nominal)
        if not all(report['nominal_containment'].values()):
            raise AssertionError('Installed reference incoming rows outside bounds')
        reports.append(report)
    result=dict(method='Analytic plateau and outward transformed rectangles with monotone cutoff bounds',reports=reports)
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps([dict(panels=r['panels'],widths=r['relative_widths'],containment=r['nominal_containment']) for r in reports]),flush=True)
    return result


if __name__=='__main__':run()
