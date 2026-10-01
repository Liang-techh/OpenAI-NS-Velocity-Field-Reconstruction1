"""Fresh resumable local-family radial core for the Md1.1 pressure datum.

Old finite tensors and old analytic-tail certificates are not reused. This is
radial coefficient generation, not temporal n-dependent background recursion.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_Md11_pressure_datum import pressure_jets
from lei_ren_part1_paper_candidate_general_center_factory import PARAMETERS
from lei_ren_part1_paper_candidate_exact_amplitude import exact_amplitude,_encode
from lei_ren_part1_paper_uniform_axis_jets import uniform_axis_jets
from lei_ren_part1_paper_factored_core_positivity import squared_axis_rows
from lei_ren_part1_paper_functional_core_recursion import coupled_rows,_source_hashes
from lei_ren_part1_paper_functional_core_step import advance_one
from lei_ren_part1_paper_candidate_gauge_core import _pack_fixed,_unpack_fixed,_pack_rows,_unpack_rows,_atomic_write_json
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
DATUM='lei_ren_part1_paper_Md11_pressure_datum.json'
STEM='lei_ren_part1_paper_Md11_interval_core_Z049_Z051'


class Md11CoreFactory:
    def __init__(self,precision=260):
        self.precision=precision;self.ctx=MPIntervalContext();self.ctx.dps=precision
        raw=(HERE/DATUM).read_bytes();self.datum=json.loads(raw)
        for name,digest in self.datum['input_hashes'].items():
            if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:
                raise ValueError('new pressure dependency changed:'+name)
        if not self.datum['explicit_Md_greater_than_one'] or not self.datum['all_14_true_pressure_stages_included']:
            raise ValueError('new fourteen-stage Md>1 pressure datum required')
        self.hashes=dict(self.datum['input_hashes']);self.hashes[DATUM]=hashlib.sha256(raw).hexdigest()
        self.hashes.update(_source_hashes())
        for name in (Path(__file__).name,'lei_ren_part1_paper_Md11_pressure_datum.py',
                     'lei_ren_part1_paper_candidate_exact_amplitude.py','lei_ren_part1_paper_uniform_axis_jets.py',
                     'lei_ren_part1_paper_factored_core_positivity.py','lei_ren_part1_paper_functional_core_step.py',
                     'lei_ren_part1_paper_candidate_gauge_core.py','lei_ren_part1_paper_candidate_general_center_factory.py'):
            self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()

    def seed(self,center,length):
        c=self.ctx
        with mp.workdps(self.precision+40):
            z=c.mpf(center)
            amplitude=exact_amplitude(c,Lambda=PARAMETERS['Lambda'],Z=z)
            axis=uniform_axis_jets(c,radius=1,j=PARAMETERS['j'],Lambda=PARAMETERS['Lambda'],
                logC=PARAMETERS['logC'],delta=PARAMETERS['delta'],length=length,axial_interval=z)
            pressure=pressure_jets(c,z,length-1,self.datum)
            gradient=list(axis['gradient_coefficients']);f0=amplitude['F0_interval']
            a,b=endpoints(axis['F0_interval']);d,e=endpoints(f0)
            if not a<=d<=e<=b:raise ArithmeticError('amplitude outside analytic enclosure')
            lam=c.mpf(PARAMETERS['Lambda'])
            fixed=dict(ell_Z_taylor=[-lam*g for g in gradient],
                S_Z_taylor=squared_axis_rows(c,gradient,lam,f0,length),U0_Z_taylor=list(axis['U0']),
                P0_Z_taylor=list(pressure['physical_pressure_coefficients']),
                gradient_coefficients=gradient,F0_interval=f0)
            rows=coupled_rows(c,fixed,z,0,PARAMETERS['delta'],required_depth=length-1)
            return fixed,{name:rows[name] for name in ('A','Uz','P')}


def run(seconds=40,max_steps=124,target_degree=124,center=('.49','.51'),stem=STEM):
    if seconds<0 or max_steps<0 or target_degree<3:raise ValueError('invalid compute budget or target')
    factory=Md11CoreFactory();c=factory.ctx;center=list(center)
    identity=dict(center=center,target_radial_degree=target_degree,required_axial_depth=3,
        initial_axis_length=target_degree+4,parameters=PARAMETERS,precision=factory.precision,
        new_schedule_sha256=factory.datum['new_schedule_sha256'])
    state_path=HERE/(stem+'_state.json');receipt_path=HERE/(stem+'.json')
    with mp.workdps(factory.precision+40):
        if state_path.exists():
            state=json.loads(state_path.read_bytes())
            if state['identity']!=identity or state['input_hashes']!=factory.hashes:
                raise ValueError('resume source or identity changed')
            fixed=_unpack_fixed(c,state['fixed'])
            rows={name:_unpack_rows(c,state['rows'][name]) for name in ('A','Uz','P')}
            completed=state['completed_radial_order']
            if not 0<=completed<=target_degree or any(len(rs)!=completed+1 or any(
                len(row)!=target_degree+4-n for n,row in enumerate(rs)) for rs in rows.values()):
                raise ValueError('resume row/depth mismatch')
        else:
            fixed,rows=factory.seed(center,target_degree+4);completed=0
            state=dict(identity=identity,input_hashes=factory.hashes,completed_radial_order=0,
                fixed=_pack_fixed(fixed),rows={n:_pack_rows(rs) for n,rs in rows.items()},
                step_timings_seconds=[],elapsed_compute_seconds=0,temporal_recursion=False,
                old_core_tensor_used=False,analytic_tail_certified=False)
            _atomic_write_json(state_path,state)
        started=time.monotonic();updates=0
        print('Fresh Md1.1 interval radial core:',completed,'of',target_degree,flush=True)
        while completed<target_degree and updates<max_steps and time.monotonic()-started<seconds:
            before=time.monotonic()
            advance_one(c,fixed,rows,completed,c.mpf(center),c.mpf(PARAMETERS['delta']))
            duration=time.monotonic()-before;completed+=1;updates+=1
            state['completed_radial_order']=completed
            state['rows']={name:_pack_rows(rs) for name,rs in rows.items()}
            state['step_timings_seconds'].append(duration);state['elapsed_compute_seconds']+=duration
            _atomic_write_json(state_path,state)
            print('Fresh Md1.1 radial order',completed,'of',target_degree,'seconds',round(duration,3),flush=True)
        receipt=dict(identity=identity,completed_radial_order=completed,completed_target=completed==target_degree,
            state_file=state_path.name,state_sha256=hashlib.sha256(state_path.read_bytes()).hexdigest(),
            input_hashes=factory.hashes,exact_interval_checkpoint=True,retained_axial_depth=3,
            new_fourteen_stage_pressure_used=True,old_core_tensor_used=False,analytic_tail_certified=False,
            production_inlet_certified=False,new_moment_inverse_generated=False,
            temporal_recursion=False,last_run_updates=updates,
            elapsed_compute_seconds=state['elapsed_compute_seconds'])
        receipt_path.write_text(json.dumps(_encode(receipt),indent=2)+'\n',encoding='utf-8')
        return receipt

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--seconds',type=float,default=40)
    parser.add_argument('--max-steps',type=int,default=124);args=parser.parse_args()
    run(seconds=args.seconds,max_steps=args.max_steps)
