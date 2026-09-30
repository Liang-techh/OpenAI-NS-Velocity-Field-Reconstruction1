"""Resolved axial-stencil checks of complete finite C2 corrected field."""
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_axial_second_jet import AxialSecondJet
from lei_ren_part1_paper_five_bump_map import FiveBumpMomentMap
from lei_ren_part1_paper_five_bump_inverse import iterate_inverse
from lei_ren_part1_paper_second_five_bump_field import evaluate_second_bump_field


def run():
    with mp.workdps(120):
        moment_map=FiveBumpMomentMap(precision=120,order=96)
        def field(z):
            zd=AxialSecondJet(z,1,0);am=20/(1+zd*zd)
            d=[mp.mpf('1e-10')*(i+1)*(1+zd*zd) for i in range(5)]
            datum=-3+zd*zd
            source=dict(Z=z,Am=am.value,Am_Z=am.tangent,Am_ZZ=am.second,
                        P0_dual=datum,parts_second_jet={f'analytic_{i+1}':{i+1:v} for i,v in enumerate(d)})
            inverse=iterate_inverse(moment_map,d,am,steps=16)
            return evaluate_second_bump_field(moment_map,inverse,source,Rm=2,Z=z,x='1.25',delta='.01')
        z=mp.mpf('.3');step=mp.mpf('1e-8')
        samples={k:field(z+k*step) for k in (-2,-1,0,1,2)}
        errors={}
        for name in ('F','Utheta','Uz','P','P0'):
            v={k:r['second_jet_fields'][name].value for k,r in samples.items()}
            expected=(-v[2]+16*v[1]-30*v[0]+16*v[-1]-v[-2])/(12*step*step)
            actual=samples[0]['second_jet_fields'][name].second
            errors[name]=abs((actual-expected).evaluate())/max(1,abs(expected.evaluate()))
        for name in samples[0]['second_jet_moments']:
            v={k:r['second_jet_moments'][name].value for k,r in samples.items()}
            expected=(-v[2]+16*v[1]-30*v[0]+16*v[-1]-v[-2])/(12*step*step)
            actual=samples[0]['second_jet_moments'][name].second
            errors['moment_'+name]=abs((actual-expected).evaluate())/max(1,abs(expected.evaluate()))
        ur={k:v['Ur'] for k,v in samples.items()}
        expected=(ur[-2]-8*ur[-1]+8*ur[1]-ur[2])/(12*step)
        errors['Ur_Z']=abs((samples[0]['Ur_Z']-expected).evaluate())/max(1,abs(expected.evaluate()))
        maximum=max(errors.values())
        assert maximum<mp.mpf('1e-24'),errors
        report=dict(maximum_scaled_error=mp.nstr(maximum,60),errors={k:mp.nstr(v,60) for k,v in errors.items()},
                    second_Z_stencil_step='1e-8',scope='resolved analytic input family and finite inverse/field',
                    source_fixture_parts=5,actual_source_certified=False,uniform_Z_certified=False)
        Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
        print('resolved C2 corrected field axial-stencil error',mp.nstr(maximum,20),flush=True)


if __name__=='__main__':run()
