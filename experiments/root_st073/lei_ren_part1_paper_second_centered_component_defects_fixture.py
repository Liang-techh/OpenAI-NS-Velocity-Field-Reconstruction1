"""Second-Z check against direct scalar density integration and axial stencil."""
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_axial_second_jet import AxialSecondJet
from lei_ren_part1_paper_pressure_width_long_reshape_fixture import _ResolvedR110Provider,_build_comparison
from lei_ren_part1_paper_centered_component_defects_fixture import _direct_rows,LOG_C,LOG_PSTAR,A
from lei_ren_part1_paper_second_centered_component_defects import SecondCenteredComponentDefects


class ResolvedSecondSource(_ResolvedR110Provider):
    """The fixture's F,Uz,P and five moments are explicitly affine in Z."""
    pressure_order=3
    width_order=2
    def evaluate_R(self,R,Z):
        result=super().evaluate_R(R,Z);zero=self.comparison._base_jet(0)
        result['second_jet_fields']={name:AxialSecondJet(result[name],result[name+'_Z'],zero)
                                     for name in ('F','Uz','P')}
        result['second_jet_moments']={name:AxialSecondJet(v,result['moments_Z'][name],zero)
                                      for name,v in result['moments'].items()}
        result['P0_ZZ']=zero
        return result


def run():
    with mp.workdps(100):
        z=mp.mpf('.3');step=mp.mpf('1e-4')
        source=ResolvedSecondSource(_build_comparison(100))
        defects=SecondCenteredComponentDefects(source,A=A,logC=LOG_C,logPstar=LOG_PSTAR,
                                               order=32,window=24,restore_order=64)
        result=defects.evaluate(z)
        samples={k:_direct_rows(source,z+k*step) for k in (-2,-1,0,1,2)}
        errors={};replay_errors={}
        for row in range(1,6):
            expected=(-samples[2][row]+16*samples[1][row]-30*samples[0][row]
                      +16*samples[-1][row]-samples[-2][row])/(12*step*step)
            actual=result['d_ZZ'][row].evaluate()
            if row==1:
                # This resolved family's first row is exactly affine. A
                # relative error to its zero second derivative is undefined.
                assert not result['d_ZZ'][row].atoms
                errors[str(row)]=abs(expected)/max(abs(samples[0][row]),mp.mpf('1e-100'))
            else:
                errors[str(row)]=abs(actual-expected)/max(abs(actual),abs(expected)) if actual or expected else mp.mpf(0)
            with mp.workdps(defects.work_precision):
                total=sum((rows.get(row,0) for rows in result['parts_ZZ'].values()),defects._jet(0))
                replay_errors[str(row)]=max((abs(v) for v in (total-result['d_ZZ'][row]).atoms.values()),default=mp.mpf(0))
        maximum=max(errors.values())
        print('second centered relative errors',{k:mp.nstr(v,15) for k,v in errors.items()},flush=True)
        assert maximum<mp.mpf('1e-8'),errors
        assert max(replay_errors.values())==0
        assert not (result['P0_ZZ']-source.evaluate_R(110,z)['P0_ZZ']).atoms
        report=dict(second_Z_relative_errors={k:mp.nstr(v,60) for k,v in errors.items()},
                    row1_metric='exact zero second derivative; stencil absolute error divided by input value',
                    maximum_relative_error=mp.nstr(maximum,60),part_second_replay_exact=True,
                    part_count=len(result['parts_ZZ']),P0_ZZ_preserved=True,
                    scalar_direct_density_oracle=True,second_Z_stencil_step='1e-4',
                    uniform_Z_certified=False,quadrature_enclosed=False,actual_source_certified=False)
        Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':run()
