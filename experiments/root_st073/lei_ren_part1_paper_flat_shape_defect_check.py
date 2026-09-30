"""Signed-log flat-shape sources from the committed actual leading atom."""
import json
from pathlib import Path
import mpmath as mp


def run():
    from lei_ren_part1_paper_flat_shape_defect import evaluate_flat_shape_defect
    source=Path(__file__).with_name('lei_ren_part1_paper_pressure_width_switches_check.json')
    data=json.loads(source.read_text(encoding='utf-8'))
    atoms=data['local_reshape_input_components']
    with mp.workdps(300):
        B=mp.mpf(atoms['reshape_B']['(0, 0)']['arbitrary_exponent_value'])
        BZ=mp.mpf(atoms['reshape_B_Z']['(0, 0)']['arbitrary_exponent_value'])
        T=mp.mpf('4e152'); logalpha=T-10*(mp.mpf('5e151')+14)+6
        rows=[]
        for name,k,m,power,half,sign in [('angular','1.6',1,'1.6',False,1),('pressure','.2',2,'.2',True,1),('angular_energy','1.2',2,'1.2',True,-1)]:
            print('actual flat defect '+name,flush=True)
            v=evaluate_flat_shape_defect(k,m,B,T,precision=300,order=32,window=24)
            pref=mp.mpf(power)*logalpha-(mp.log(2) if half else 0)
            kk=mp.mpf(k); mode0=(2*T*T/kk)**(mp.mpf(1)/3)
            leading=-mp.mpf('1.5')*kk*mode0+1+mp.log(abs(m*B))+mp.log(2*mp.pi*mode0/(3*kk))/2
            asymptotic_error=abs(v['log_abs']-leading)
            assert asymptotic_error<mp.mpf('1e-8'),asymptotic_error
            rows.append(dict(name=name,source_sign=sign*int(v['sign']),source_log_abs=mp.nstr(pref+v['log_abs'],100),
                source_Z_sign=sign*int(v['derivative_sign'])*int(mp.sign(BZ)),
                source_Z_log_abs=mp.nstr(pref+v['derivative_log_abs']+mp.log(abs(BZ)),100),
                leading_saddle_log_discrepancy=mp.nstr(asymptotic_error,60),
                kernel_sign=int(v['sign']),kernel_log_abs=mp.nstr(v['log_abs'],100),
                derivative_log_abs=mp.nstr(v['derivative_log_abs'],100),
                mode=mp.nstr(v['mode'],100),width=mp.nstr(v['width'],100),
                old_cutoff='1842.0680743952365',mode_beyond_old_cutoff=bool(v['mode']>1843)))
            assert v['sign']!=0 and v['mode']>1843
        report=dict(source=str(source.name),source_scope='Serialized actual leading (0,0) B/B_Z atom; not full pressure-width ring',
            T=mp.nstr(T,80),B=mp.nstr(B,80),BZ=mp.nstr(BZ,80),rows=rows,
            defect_kernel_nonzero=True,all_five_functional_defects_installed=False,
            finite_window_error_enclosed=False,quadrature_error_enclosed=False,
            functional_terminal_moments_closed=False,cone_certified=False,temporal_recursion_certified=False)
        Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
        return report

if __name__=='__main__':
    run()
