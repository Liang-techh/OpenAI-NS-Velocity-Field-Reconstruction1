"""Actual serialized pressure9/width2 B atoms through flat defect composition."""
import ast
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_pressure_width_jet import PressureWidthJet
from lei_ren_part1_paper_axial_correction import signed_log


def run():
    from lei_ren_part1_paper_flat_shape_component import compose_flat_shape_defect
    source=Path(__file__).with_name('lei_ren_part1_paper_pressure_width_switches_check.json')
    inputs=json.loads(source.read_text(encoding='utf-8'))['local_reshape_input_components']
    with mp.workdps(500):
        def decode(name):
            return PressureWidthJet({ast.literal_eval(k):mp.mpf(v['arbitrary_exponent_value'])
                for k,v in inputs[name].items()},pressure_order=9,width_order=2)
        B=decode('reshape_B');BZ=decode('reshape_B_Z');T=mp.mpf('4e152')
        encode=lambda a:{str(k):signed_log(v,120) for k,v in a.atoms.items()}
        rows=[]
        for name,k,m in [('angular','1.6',1),('pressure','.2',2),('angular_energy','1.2',2)]:
            print('actual component flat defect '+name,flush=True)
            v=compose_flat_shape_defect(k,m,B,BZ,T,precision=500,order=32,window=24)
            assert v['value'].component(0,0)!=0
            assert v['value'].component(1,0)!=0
            assert v['value'].component(0,1)!=0
            assert v['tangent'].component(1,0)!=0
            assert v['value_terms'][2].component(2,0)!=0
            rows.append(dict(name=name,value=encode(v['value']),tangent=encode(v['tangent']),
                value_terms={str(n):encode(a) for n,a in v['value_terms'].items()},
                tangent_terms={str(n):encode(a) for n,a in v['tangent_terms'].items()},
                metadata=v.get('metadata',{})))
        report=dict(source=source.name,source_scope='Serialized actual B/B_Z pressure9/width2 atoms; source precision inherited',
            pressure_order=9,width_order=2,T=mp.nstr(T,80),rows=rows,
            per_derivative_order_terms_preserved=True,all_five_functional_defects_installed=False,
            quadrature_error_enclosed=False,finite_ring_remainder_enclosed=False,
            functional_terminal_moments_closed=False,temporal_recursion_certified=False)
        Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
        print('actual component flat sources saved',flush=True)
        return report

if __name__=='__main__':
    run()
