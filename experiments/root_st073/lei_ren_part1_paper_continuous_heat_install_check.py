"""Shared-field heat jet and separately retained tiny atom checks."""
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_joined_outer import build_joined_field
from lei_ren_part1_paper_continuous_incoming_outer import ContinuousIncomingOuterField
from lei_ren_part1_paper_axial_correction import signed_log


def run():
    folder=Path(__file__).parent
    print('building shared field for heat installation',flush=True)
    source=build_joined_field()
    seeded=json.loads((folder/'lei_ren_part1_paper_seeded_shared_candidate.json').read_text())
    atoms=json.loads((folder/'lei_ren_part1_paper_continuous_axial_solve.json').read_text())
    field=ContinuousIncomingOuterField(source,prepared={.3:(seeded,atoms)})
    p=field.outer;s=p.schedule
    with mp.workdps(field.precision):
        z=mp.mpf('.3');step=mp.mpf('1e-4');samples=[]
        for local in ('.6','4','30'):
            center=p.log_at(s.logR_tail,local)
            row=s.at_log_radius(center,z);actual=p.values_with_jets(center,z)
            if not(row['heat_deficit']>0 and actual['angular_heat_relative_correction']<0):
                raise ArithmeticError('Shared heat deficit disappeared')
            neighbors={i:s.at_log_radius(p.log_at(center,i*step),z) for i in (-2,-1,1,2)}
            derivative=sum(w*neighbors[i]['heat_log_amplitude_correction'] for i,w in ((-2,1),(-1,-8),(1,8),(2,-1)))/(12*step)
            slope_error=abs(derivative/row['heat_logarithmic_slope_correction']-1)
            neighborsZ={i:s.at_log_radius(center,z+i*step) for i in (-2,-1,1,2)}
            derivativeZ=sum(w*neighborsZ[i]['heat_log_amplitude_correction'] for i,w in ((-2,1),(-1,-8),(1,8),(2,-1)))/(12*step)
            jet_error=abs(derivativeZ/mp.mpf(str(row['dlogU_dZ']))-1)
            field_error=abs(actual['Utheta_Z']/(actual['Utheta']*mp.mpf(str(row['dlogU_dZ'])))-1)
            if max(slope_error,jet_error,field_error)>mp.mpf('1e-8'):
                raise ArithmeticError('Installed heat field jet mismatch')
            samples.append(dict(t=local,deficit=signed_log(row['heat_deficit'],50),
                log_correction=signed_log(row['heat_log_amplitude_correction'],50),
                slope_relative_error=mp.nstr(slope_error,40),Z_relative_error=mp.nstr(jet_error,40),
                field_Z_relative_error=mp.nstr(field_error,40)))
        report=dict(samples=samples,continuous_heat_point_kernel_installed=True,
            separate_tiny_heat_atoms_retained=True,heat_integral_targets_regenerated=False,
            heat_integral_targets_inherited_Taylor=True,complete_heat_moments=False,
            arithmetic_error_enclosed=False,finite_energy_certified=False,scale_recursion_certified=False)
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'samples':[{k:v for k,v in r.items() if 'error' in k or k=='t'} for r in samples]}),flush=True)
    return report


if __name__=='__main__':run()
