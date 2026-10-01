"""Same-source Md40 normalized finite core enclosures, with resumable rows.

Axis F0 is implicit and positive. Its very small squared-amplitude jets are
retained as nonzero analytic Cauchy enclosures, not replaced by zero. This
is radial profile generation; temporal n-dependent recursion is separate.
"""
import argparse
import hashlib
import json
import time
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_logarithmic_pressure_datum import LogarithmicPressureDatum
from lei_ren_part1_paper_interval_taylor import IntervalTaylor, constant
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_functional_core_step import initial_rows
from lei_ren_part1_paper_logarithmic_core_step import advance_scaled_one
from lei_ren_part1_paper_candidate_gauge_core import (
    _pack_vector,_unpack_vector,_pack_rows,_unpack_rows,_atomic_write_json)

HERE=Path(__file__).parent
STEM='lei_ren_part1_paper_logarithmic_interval_core_Z049_Z051'

class LogarithmicCoreFactory:
    def __init__(self):
        self.datum=LogarithmicPressureDatum('40',precision=160)
        self.ctx=c=self.datum.ctx;self.precision=c.dps
        self.majorant=json.loads((HERE/'lei_ren_part1_paper_logarithmic_core_majorant.json').read_bytes())
        self.tail=json.loads((HERE/'lei_ren_part1_paper_logarithmic_core_tail_admission.json').read_bytes())
        self.binding=json.loads((HERE/'lei_ren_part1_paper_logarithmic_source_binding.json').read_bytes())
        for record in (self.majorant,self.tail,self.binding):
            for name,digest in record['input_hashes'].items():
                if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:
                    raise ValueError('dependency changed: '+name)
            if record['implicit_source_sha256']!=self.datum.source_sha or record['datum_enclosure_sha256']!=self.datum.datum_sha:
                raise ValueError('source or datum identity mismatch')
        if not self.majorant['contraction_proved'] or not self.binding['core_analytic_contraction_and_scaled_tail_admitted']:
            raise ValueError('analytic admission required')
        self.degree=self.binding['selected_radial_degree']
        self.hashes=dict(self.binding['input_hashes'])
        names=(Path(__file__).name,'lei_ren_part1_paper_logarithmic_source_binding.json',
            'lei_ren_part1_paper_logarithmic_core_step.py','lei_ren_part1_paper_functional_core_step.py',
            'lei_ren_part1_paper_candidate_gauge_core.py','lei_ren_part1_paper_interval_taylor.py',
            'lei_ren_part1_paper_uniform_fixed_beta_error.py','lei_ren_part1_paper_schedule_endpoint_enclosures.py')
        self.hashes.update({n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names})
        with mp.workdps(c.dps+40):
            self.eps=read_interval(c,self.majorant['epsilon'])
            self.logLambda=read_interval(c,self.majorant['logLambda'])
            self.eta=read_interval(c,json.loads((HERE/'lei_ren_part1_paper_analytic_radial_tail.json').read_bytes())['complex_tube_radius'])
            self.delta=self.datum.parameters.delta
            if not 0<endpoints(self.delta)[0]<=endpoints(self.delta)[1]<=mp.mpf('1e-200'):
                raise ValueError('delta outside uniform complex tube bound')

    def seed(self,center,length):
        c=self.ctx
        with mp.workdps(c.dps+40):
            zc=c.mpf(center)
            if endpoints(zc)[0]<-1 or endpoints(zc)[1]>1:raise ValueError('centers outside real axis')
            z=IntervalTaylor.variable(c,zc,length-1);one=constant(c,1,length-1)
            j=c.mpf('1e-14');sigma=j/500;dt=self.delta
            u=4*z+j;L=one-z*z*dt;H=z*((1-dt)/2)+(one-z*z)*u
            den=H*H+sigma*sigma
            co=list(den.coefficients);co[0]=H.coefficients[0]**2+sigma*sigma
            g=(L*H)/IntervalTaylor(c,co)
            # F0=exp(-logC-Lambda G). On the certified complex tube,
            # logC=Lambda Gbar+2logLambda+1000, |G|<=Gbar.
            # Thus epsilon^2 |F0|^2 <=exp(-6logLambda-2000).
            bound=c.exp(-6*self.logLambda-2000);rho=self.eta/2
            S=[]
            for k in range(length):
                b=bound/rho**k
                S.append(c.mpf([0,endpoints(b)[1]]) if k==0 else
                    c.mpf([endpoints(-b)[0],endpoints(b)[1]]))
            if endpoints(S[0])[1]<=0:raise ArithmeticError('swirl bound lost')
            pressure=self.datum.normalized_jets(zc,length-1)
            factor=c.exp(2*self.datum.parameters.logPstar-self.logLambda)
            fixed=dict(ell_Z_taylor=[-x for x in g.coefficients],S_Z_taylor=S,
                U0_Z_taylor=list(u.coefficients),
                P0_Z_taylor=[factor*x for x in pressure['normalized_pressure_coefficients']])
            metadata=dict(implicit_positive_F0_definition='exp(-logC-Lambda G), G primitive anchored at H root in [-j,0]',
                logC_definition=self.majorant['logC_definition'],
                S_scaled_complex_supremum=_pack_vector([bound])[0],
                S_Cauchy_radius=_pack_vector([rho])[0],
                actual_positive_F0_not_numerically_selected=True,
                tiny_swirl_not_replaced_by_zero=True,
                amplitude_lower_bound_not_materialized=True,
                seed_is_same_source_analytic_enclosure=True)
            return fixed,initial_rows(c,fixed,zc,self.degree,required_depth=3),metadata

def run(seconds=30,max_steps=110):
    if seconds<0 or max_steps<0:raise ValueError('negative budget')
    f=LogarithmicCoreFactory();c=f.ctx;center=['.49','.51'];target=f.degree
    identity=dict(center_Z=center,target_radial_degree=target,initial_axis_length=target+4,
        retained_axial_depth=3,precision=c.dps,implicit_source_sha256=f.datum.source_sha,
        datum_enclosure_sha256=f.datum.datum_sha,
        row_units=dict(A='physical_A[n]/Lambda**n',Uz='physical_Uz[n]/Lambda**n',
            P='epsilon*physical_P[n]/Lambda**n'),
        epsilon=_pack_vector([f.eps])[0],delta=_pack_vector([f.delta])[0])
    state_path=HERE/(STEM+'_state.json');receipt_path=HERE/(STEM+'.json')
    with mp.workdps(c.dps+40):
        if state_path.exists():
            state=json.loads(state_path.read_bytes())
            if state['identity']!=identity or state['input_hashes']!=f.hashes:
                raise ValueError('resume identity changed')
            fixed={name:_unpack_vector(c,row) for name,row in state['fixed'].items()}
            rows={name:_unpack_rows(c,row) for name,row in state['rows'].items()}
            completed=state['completed_radial_order']
            if not 0<=completed<=target or any(len(rs)!=completed+1 or any(
                len(row)!=target+4-n for n,row in enumerate(rs)) for rs in rows.values()):
                raise ValueError('resume shapes invalid')
        else:
            fixed,rows,seed_metadata=f.seed(center,target+4);completed=0
            state=dict(identity=identity,input_hashes=f.hashes,seed_metadata=seed_metadata,
                completed_radial_order=0,fixed={k:_pack_vector(v) for k,v in fixed.items()},
                rows={k:_pack_rows(v) for k,v in rows.items()},elapsed_compute_seconds=0,
                step_timings_seconds=[],temporal_recursion=False,old_core_tensor_used=False)
            _atomic_write_json(state_path,state)
        updates=0;start=time.monotonic()
        print('Md40 direct scaled radial core',completed,'of',target,flush=True)
        while completed<target and updates<max_steps and time.monotonic()-start<seconds:
            step=time.monotonic()
            advance_scaled_one(c,fixed,rows,completed,c.mpf(center),f.delta,f.eps)
            duration=time.monotonic()-step;completed+=1;updates+=1
            state['completed_radial_order']=completed;state['elapsed_compute_seconds']+=duration
            state['step_timings_seconds'].append(duration)
            state['rows']={name:_pack_rows(rs) for name,rs in rows.items()}
            _atomic_write_json(state_path,state)
            print('Md40 radial order',completed,'of',target,'seconds',round(duration,3),flush=True)
        receipt=dict(identity=identity,input_hashes=f.hashes,completed_radial_order=completed,
            completed_target=completed==target,state_file=state_path.name,
            state_sha256=hashlib.sha256(state_path.read_bytes()).hexdigest(),
            fresh_same_source_seed_generated=True,finite_radial_rows_generated=completed>0,
            seed_is_analytic_enclosure_not_point_amplitude=True,
            tiny_swirl_not_replaced_by_zero=True,old_core_tensor_used=False,
            analytic_tail_admission_available=True,finite_plus_tail_exit_validation_completed=False,
            reference_inlet_built=False,whole_axis_core_generated=False,
            temporal_recursion=False,background_assembled=False,
            last_run_updates=updates,elapsed_compute_seconds=state['elapsed_compute_seconds'])
        _atomic_write_json(receipt_path,receipt)
        return receipt

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--seconds',type=float,default=30)
    p.add_argument('--max-steps',type=int,default=110);a=p.parse_args()
    run(a.seconds,a.max_steps)
