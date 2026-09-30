"""Resolved source-lineage and independent radial-density reference fixture."""
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_axial_dual import AxialDual
from lei_ren_part1_paper_five_moment_reference_background import ReferenceDefectBackground
from lei_ren_part1_paper_five_bump_map import FiveBumpMomentMap
from lei_ren_part1_paper_five_bump_field import evaluate_correction,join


class Source:
    precision=100
    work_precision=100
    def __init__(self):
        self.Rm=mp.mpf(16);self.logPstar=mp.mpf('1.4');self.delta=mp.mpf('.001');self.calls=0
    def evaluate(self,Z):
        self.calls+=1;z=mp.mpf(Z)
        dual=lambda a,b=0:AxialDual(a,b,pressure_order=0,width_order=0)
        d={i:dual(mp.mpf(i)*mp.mpf('1e-7')*(1+z*z),2*z*mp.mpf(i)*mp.mpf('1e-7')) for i in range(1,6)}
        p0=dual(mp.mpf('-.2')+mp.mpf('.04')*z,mp.mpf('.04'))
        return dict(d_dual=d,d={i:v.value for i,v in d.items()},d_Z={i:v.tangent for i,v in d.items()},
          parts_dual={str(i):{i:v} for i,v in d.items()},P0=p0.value,P0_Z=p0.tangent)


def run():
    with mp.workdps(100):
        source=Source();reference=ReferenceDefectBackground(source)
        assert reference.delta==source.delta
        z=mp.mpf('.3');x=mp.mpf('1.25');R=source.Rm*x
        result=reference.evaluate_x(x,z)
        scalar=lambda v:v.component(0,0)
        reconstruction=max(abs(sum(scalar(rows[k]) for rows in result['moment_parts'].values())-scalar(result['moments'][k])) for k in result['moments'])
        assert reconstruction<mp.mpf('1e-90')
        assert 'reference_power' in result['moment_parts']
        assert result['P0']==reference.defect_data(z)['P0']
        step=mp.mpf('1e-6');samples={j:reference.evaluate_x((R+j*step)/source.Rm,z) for j in (-2,-1,1,2)}
        u=scalar(result['Utheta']);v=scalar(result['Uz'])
        densities=dict(theta=mp.sqrt(2*R)*u,z=v,theta_z=mp.sqrt(2*R)*u*v,z_theta=v*v-u*u/2,p=u*u/(2*R))
        errors={}
        for k,expected in densities.items():
            actual=(-scalar(samples[2]['moments'][k])+8*scalar(samples[1]['moments'][k])-8*scalar(samples[-1]['moments'][k])+scalar(samples[-2]['moments'][k]))/(12*step)
            errors[k]=abs((actual-expected)/expected)
        assert max(errors.values())<mp.mpf('1e-25')
        terminal=reference.evaluate_phase(3,z)
        labels=[k for k in result['moment_parts'] if k!='reference_power']
        constant=max(abs(scalar(result['moment_parts'][label][k])-scalar(terminal['moment_parts'][label][k])) for label in labels for k in densities)
        assert constant==0
        assert source.calls==1
        assert result['raw_quadratic_integrals'] is None
        m=FiveBumpMomentMap(precision=100,order=48)
        data=reference.defect_data(z);h=[-v for v in m.linear_inverse([data['d_dual'][i] for i in range(1,6)])]
        amplitude=mp.exp(source.logPstar-mp.mpf('.6'))/(1+z*z)
        amp=AxialDual(amplitude,-2*z*amplitude/(1+z*z),pressure_order=0,width_order=0)
        correction=evaluate_correction(m,h,amp,source.Rm,z,x,delta=source.delta)
        corrected=join(result,correction,delta=source.delta)
        assert corrected['P0']==result['P0']
        assert corrected['raw_quadratic_integrals'] is None
        assert not corrected['raw_quadratic_parts_complete']
        assert 'five_bump_correction' in corrected['moment_parts']
        report=dict(radial_density_relative_errors={k:mp.nstr(v,30) for k,v in errors.items()},
          moment_part_reconstruction_absolute_error=mp.nstr(reconstruction,30),
          constant_defect_parts=True,same_source_cached_once=True,P0_preserved=True,optional_raw_reference_join_pass=True,
          raw_individual_energies_not_inferred=True,actual_source_scope=False,functional_closure=False,quadrature_enclosed=False)
        Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps(report,indent=2))
        return report

if __name__=='__main__':run()
