"""Fourteen-stage preheat pressure envelope for the logarithmic source.

All post-yd atoms remain positive intervals. Their shared analytic tail
bound encloses the complete integral regardless of enormous later radii.
The datum is an enclosure of an implicitly defined continuous source;
it is not a fitted pressure and does not instantiate its corrected core.
"""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_logarithmic_outer_parameters import LogarithmicOuterParameters
from lei_ren_part1_paper_candidate_pressure_function import q_jets
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_interval_exit_continuation_enclosure import _pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
BETA2=('reference_extension','slope_transition_ref','axial_turnoff',
    'slope_transition_mu','power_buffer','pulse_reserved')
BETA0=('power_buffer_rel','steep_transition_in','steep_power','steep_transition_out',
    'waiting','heat_collar','exterior_power_tail')


class LogarithmicPressureDatum:
    def __init__(self,Md='40',precision=160):
        self.parameters=LogarithmicOuterParameters(Md,precision=precision)
        self.ctx=c=self.parameters.ctx;p=self.parameters
        with mp.workdps(precision+40):
            refinement_name='lei_ren_part1_paper_candidate_pressure_mass_refinement.json'
            refined=json.loads((HERE/refinement_name).read_bytes())
            expected='736bbadbde99bc2f3d098d279d61ef4cb64418368263a4aba7b275e7f8892de4'
            if refined['accepted_schedule_sha256']!=expected or 'slope_transition_ref' not in refined['stages']:
                raise ValueError('unexpected refinement source or atom identity')
            for name,digest in refined['input_hashes'].items():
                if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:
                    raise ValueError('refinement dependency changed: '+name)
            self.input_hashes={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
                (Path(__file__).name,refinement_name,'lei_ren_part1_paper_logarithmic_outer_parameters.py',
                 'lei_ren_part1_paper_candidate_pressure_function.py','lei_ren_part1_paper_outer.py')}
            self.input_hashes.update(refined['input_hashes'])
            self.tail_upper=c.exp(c.mpf('.6')-p.yd)/(2*(1-p.epsilon)**2)
            box=c.mpf([0,endpoints(self.tail_upper)[1]])
            self.stages={name:dict(beta=2,mass=box) for name in BETA2}
            self.stages.update({name:dict(beta=0,mass=box) for name in BETA0})
            self.stages['reference_extension']['mass']=c.mpf('2.5')
            self.stages['slope_transition_ref']['mass']=read_interval(c,refined['stages']['slope_transition_ref']['refined_mass'])
            self.stages['axial_turnoff']['mass']=c.exp(c.mpf('.6'))*(c.exp(-1)-c.exp(-p.yd))/2
            self.stages['z_flatten']=dict(beta='variable [0,2]',mass=box)
            self.m2=sum((self.stages[n]['mass'] for n in BETA2),c.mpf(0))
            self.m0=sum((self.stages[n]['mass'] for n in BETA0),c.mpf(0))
            self.rho=c.mpf('.25')
            # For |Im Z|<=rho and real theta in [0,1],
            # |(1+Z^2)^(-theta)| <= (1-rho^2)^(-1).
            self.flatten_complex_upper=self.tail_upper/(1-self.rho**2)**2
            definition=dict(Md=str(Md),logPstar='exp(Md)+11',c_mu='.001',c_delta='.001',c_epsilon='.01',
                delta='min(1e-200,exp(-4logPstar-30))',
                waiting='unique positive root of the continuous raw preheat waiting equation',
                angular_profile='Section 6.1 reference-plus-outer ansatz, H replaced by 1',
                Tw='-60log(mu)',Ts='4log(2/delta)',Tf=100,
                cutoff_and_schedule_python_sha256=hashlib.sha256((HERE/'lei_ren_part1_paper_outer.py').read_bytes()).hexdigest())
            self.definition=definition
            self.source_sha=hashlib.sha256(json.dumps(definition,sort_keys=True).encode()).hexdigest()
            enclosure=dict(source_sha256=self.source_sha,input_hashes=self.input_hashes,
                stages=encode(_pack(self.stages)),complex_flatten_upper=encode(self.flatten_complex_upper),
                normalized_m2=encode(self.m2),normalized_m0=encode(self.m0))
            self.datum_sha=hashlib.sha256(json.dumps(enclosure,sort_keys=True).encode()).hexdigest()

    def normalized_jets(self,z,order):
        c=self.ctx
        if not isinstance(order,int) or isinstance(order,bool) or order<0:
            raise ValueError('nonnegative Taylor order required')
        z=c.mpf(z)
        if endpoints(z)[0]<-1 or endpoints(z)[1]>1:raise ValueError('real centers in [-1,1] required')
        with mp.workdps(c.dps+40):
            rows=[]
            for n,a in enumerate(q_jets(c,z,order)):
                if n==0:remainder=self.stages['z_flatten']['mass']
                elif n%2 and endpoints(z)==(mp.mpf(0),mp.mpf(0)):
                    remainder=c.mpf(0)
                else:
                    bound=self.flatten_complex_upper/self.rho**n
                    remainder=c.mpf([endpoints(-bound)[0],endpoints(bound)[1]])
                rows.append(-(self.m2*a+(self.m0 if n==0 else 0)+remainder))
            return dict(center_Z=z,ordinary_Taylor_order=order,normalized_pressure_coefficients=rows,
                pressure_units='P/Pstar^2',physical_log_scale=2*self.parameters.logPstar,
                exact_even_pressure_parity_preserved=True,
                implicit_source_sha256=self.source_sha,all_14_pressure_atoms_included=True,
                datum_enclosure_sha256=self.datum_sha,datum_input_hashes=self.input_hashes,
                no_post_yd_mass_replaced_by_zero=True,finite_quadrature_fit_used=False)

    def report(self):
        c=self.ctx
        with mp.workdps(c.dps+40):
            if len(self.stages)!=14:raise ValueError('fourteen atoms required')
            if endpoints(self.tail_upper)[0]<=0:raise ValueError('tiny tail must remain positive')
            return dict(implicit_source_definition=self.definition,implicit_source_sha256=self.source_sha,
                datum_enclosure_sha256=self.datum_sha,input_hashes=self.input_hashes,
                slope_atom_transfer_invariance='Utheta/Pstar=(1+Z^2)^-1 exp(y/10-.6J(y)), y in [0,1], independent of Md/mu/delta/waiting',
                transferred_atom_units='normalized Utheta^2/(2Pstar^2) dy with (1+Z^2)^-2 factored out',
                parameter_bounds=self.parameters.report(),stage_count_total=14,stages=self.stages,
                shared_post_yd_mass_upper=self.tail_upper,
                shared_tail_bound_is_joint_not_additive=True,
                individual_mass_boxes_used_for_safe_beta_sums=True,
                fixed_beta2_mass_normalized=self.m2,fixed_beta0_mass_normalized=self.m0,
                flatten_complex_mass_upper=self.flatten_complex_upper,complex_strip_half_width=self.rho,
                local_C3_normalized_pressure=self.normalized_jets(c.mpf(['.49','.51']),3),
                global_real_axis_jets_supported=True,complete_preheat_pressure_enclosed=True,
                implicit_continuous_waiting_root_enclosed=True,
                separate_late_atoms_integrated_sharply=False,exact_waiting_root_selected=False,
                corrected_heat_pressure_restoration_completed=False,new_core_generated=False,
                admissible_background_assembled=False,temporal_recursion=False)


def run():
    datum=LogarithmicPressureDatum();result=datum.report()
    names=(Path(__file__).name,'lei_ren_part1_paper_logarithmic_outer_parameters.py',
        'lei_ren_part1_paper_candidate_pressure_mass_refinement.json',
        'lei_ren_part1_paper_candidate_pressure_function.py','lei_ren_part1_paper_outer.py')
    result['input_hashes'].update({n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names})
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(_pack(result)),indent=2)+'\n',encoding='utf-8')
    print('Md40 fourteen-stage normalized pressure envelope complete; positive tail retained',flush=True)
    return result


if __name__=='__main__':run()
