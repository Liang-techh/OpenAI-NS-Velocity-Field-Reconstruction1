"""Independent quadrature fixtures for family defect kernel/transport bounds."""
import json
import hashlib
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_interval_taylor import IntervalTaylor,constant
import lei_ren_part1_paper_interval_functional_defects as target
from lei_ren_part1_paper_axial_primitive import _sigma_mp
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints


def contained(box,reference):
    lo,hi=endpoints(box)
    if not lo<=reference<=hi:raise AssertionError('independent reference outside interval')


def run():
    c=MPIntervalContext();c.dps=90
    with mp.workdps(130):
        B=IntervalTaylor(c,[c.mpf(['-.4','-.2']),c.mpf(['-.3','.5'])])
        kernels=0
        for k,m in (('1.6','1'),('1.2','2'),('.2','2')):
            box=target.flat_kernel(c,c.mpf(k),c.mpf(m),B,c.mpf(16))['jet']
            for b in ('-.4','-.2'):
                bv=mp.mpf(b);kv=mp.mpf(k);mv=mp.mpf(m)
                val=mp.quad(lambda l:mp.exp(-kv*l)*mp.expm1(mv*bv*_sigma_mp(l/16)),[0,8,16])
                der=mp.quad(lambda l:mp.exp(-kv*l)*mv*_sigma_mp(l/16)*mp.exp(mv*bv*_sigma_mp(l/16)),[0,8,16])
                contained(box[0],val);kernels+=1
                for bz in ('-.3','.5'):
                    contained(box[1],der*mp.mpf(bz));kernels+=1
        K=target.restore_constants(c,128)
        for key,power,m in (('K1','1',1),('K16','1.6',1),('K2','1',2)):
            ref=mp.quad(lambda t:mp.exp(mp.mpf(power)*t)*(1-_sigma_mp(t))**m,[0,mp.mpf('.5'),1])
            contained(K[key],ref)
        z=IntervalTaylor.variable(c,c.mpf('.5'),2)
        calc=SimpleNamespace(ctx=c,eps=c.mpf('.1'),u=[4*z+c.mpf('1e-14'),constant(c,'.001',2),constant(c,'-.0001',2)])
        zero=z*0
        bar=dict(phi=2+z/20,U=4*z)
        bar.update({key:zero for key in ('theta','z','theta_z','p','u_squared','weighted_phi_squared')})
        continuation=dict(endpoint_comparison_state=bar,endpoint_scaled_radius=c.mpf('.1'),
            exit_shear_epsilon=c.mpf('.01'),normalized_exit_state=dict(phi=(1+z/10).truncate(1)),
            delta_u_target=zero.truncate(1))
        bridge=dict(U_increment=(constant(c,'.00001',2)+z*c.mpf('.000002')).truncate(1))
        def driver(*args):
            return dict(D=(z/100+c.mpf('.2')).truncate(1),I_z=(z/50-c.mpf('.1')).truncate(1))
        with patch.object(target,'tail_jet',return_value=zero),patch.object(target,'_driver',side_effect=driver):
            g,parts=target.centered_axial_exit(calc,bridge,continuation)
        def direct(Z):
            phi=1+Z/10;bp=2+Z/20;D=mp.mpf('.2')+Z/100;Iz=-mp.mpf('.1')+Z/50
            hb=mp.mpf('.005');eps=mp.mpf('.01')
            change=mp.quad(lambda x:-phi/bp*mp.exp(-eps*D*x/2)*Iz*eps
                *mp.sqrt(100*mp.exp(x)/2)*(1-_sigma_mp(x/hb)),[0,hb/2,hb])
            return mp.mpf('1e-14')+mp.mpf('.001')*4-mp.mpf('.0001')*16+mp.mpf('.00001')+Z*mp.mpf('.000002')+change
        contained(g[0],direct(mp.mpf('.5')))
        contained(g[1],mp.diff(direct,mp.mpf('.5')))
        if not parts['baseline_4Z_cancelled_algebraically']:raise AssertionError('baseline not cancelled')
    here=Path(__file__).parent
    report=dict(passed=True,fixture_only=True,flat_kernel_value_and_axial_coefficients_contained=kernels,
        restoration_quadratures_contained=3,centered_axial_transport_coefficients_contained=2,
        five_row_assembly_independent_validation=False,
        input_hashes={name:hashlib.sha256((here/name).read_bytes()).hexdigest() for name in
            (target.__file__.split('\\')[-1].split('/')[-1],Path(__file__).name)})
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print('Functional defect fixtures:18 kernel,3 restoration,2 centered transport coefficients contained',flush=True)
    return report

if __name__=='__main__':run()
