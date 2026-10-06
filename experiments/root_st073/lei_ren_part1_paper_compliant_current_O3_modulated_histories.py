"""Modified local five histories, common-axis pressure and radial recovery."""
import gzip
import hashlib
import json
import operator
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_compliant_current_O3_finite_frequency_profiles import (
    CurrentO3FiniteFrequencyProfiles,CurrentO3TransitionDirection,
    HERE,PREFIX,sha,pack,encode,endpoints,source_precision,_verify_hashes,
    NAME as PROFILE_NAME,RECEIPT as PROFILE_RECEIPT)
from lei_ren_part1_paper_compliant_current_O3_modulated_histories_operator import (
    SHAPES,exact_modulated_history_theorem,cumulative_scalar_enclosures,
    frequency_uniform_bounds,increment_rows)
from lei_ren_part1_paper_compliant_current_pre_pulse_stress_operator import raw_pre_velocity_rows,copy_jet
from lei_ren_part1_paper_compliant_collar_stress_C3 import shifted_rows
from lei_ren_part1_paper_interval_taylor import IntervalTaylor

NAME=PREFIX+'current_O3_modulated_histories.json'
RECEIPT=PREFIX+'current_O3_modulated_histories_check.json'
VIEWS_NAME=PREFIX+'current_O3_modulated_histories_views.json.gz'
GATES=('current_modulated_five_cumulative_history_source_functions_available',
       'current_modulated_common_axis_absolute_pressure_source_available',
       'current_modulated_radial_velocity_from_own_moments_available',
       'current_modulated_frequency_uniform_axial_O1_over_N_bounds_available')
OPEN=('current_O3_transition_finite_N_modified_profiles_installed',
      'current_O3_transition_modified_five_moment_repair_certified',
      'current_whole_O3_transition_signed_two_vector_cone_certified',
      'current_O2_buffer_modified_whole_cone_certified',
      'current_modified_complete_tensor_and_physical_source_joins_certified',
      'current_modified_analytic_preheat_exterior_compatible_after_repair',
      'current_modified_global_admissible_tensor_and_full_NS_certified')

def valid_integer(value,label):
    if isinstance(value,bool):raise ValueError(label+' requires an integer')
    try:value=operator.index(value)
    except TypeError as exc:raise ValueError(label+' requires an integer') from exc
    if value<1:raise ValueError(label+' must be positive')
    return value

class CurrentO3ModulatedHistories:
    @source_precision
    def __init__(self,candidate=None,N=10**12,cells=8,require_checked=True):
        self.candidate=candidate if candidate is not None else CurrentO3FiniteFrequencyProfiles(CurrentO3TransitionDirection())
        if type(self.candidate) is not CurrentO3FiniteFrequencyProfiles:raise ValueError('Actual supported profile source required')
        self.candidate.direction.assert_graph();self.ctx=c=self.candidate.ctx
        self.pre=self.candidate.pre;self.mu=self.candidate.mu
        self.family=self.candidate.family;self.source=self.candidate.source;self.datum_sha=self.candidate.datum_sha
        self.N=valid_integer(N,'Finite frequency');self.cells=valid_integer(cells,'Signed range subdivisions')
        if self.cells<4 or self.cells>256 or self.cells&(self.cells-1):
            raise ValueError('A power-of-two subdivision count in[4,256] required')
        receipt=json.loads((HERE/PROFILE_RECEIPT).read_bytes());_verify_hashes(receipt)
        if not receipt['all_passed'] or tuple(receipt[k] for k in
          ('actual_five_defect_family_sha256','implicit_source_sha256','datum_enclosure_sha256'))!=(self.family,self.source,self.datum_sha):
            raise ValueError('Checked same-source finite-profile construction receipt required')
        profile=json.loads((HERE/PROFILE_NAME).read_bytes())
        if encode(pack(self.candidate.theorem))!=profile['exact_finite_frequency_profile_theorem']:
            raise ValueError('Actual local profile theorem differs')
        self.Ua=c.mpf(endpoints(self.candidate.direction.bounds['actual_Ua']))
        self.delta=c.mpf(self.pre.delta);self.theorem=exact_modulated_history_theorem()
        self.bounds=frequency_uniform_bounds(c,self.mu);self.cache={}
        self.hashes={**self.candidate.hashes,**receipt['input_hashes'],**self.theorem['input_hashes'],
          PROFILE_RECEIPT:sha(PROFILE_RECEIPT)}
        for stem in ('current_O3_modulated_histories_operator','current_O3_modulated_histories'):
            name=PREFIX+stem+'.py';self.hashes[name]=sha(name)
        definition=dict(parent_family=self.family,parent_source=self.source,parent_datum=self.datum_sha,
          finite_integer_N=self.N,signed_range_subdivisions=self.cells,
          profile_program_sha256=sha(PREFIX+'current_O3_finite_frequency_profiles.py'),
          history_program_sha256=sha(PREFIX+'current_O3_modulated_histories.py'),
          history_operator_sha256=sha(PREFIX+'current_O3_modulated_histories_operator.py'))
        self.modified_source_definition_sha256=hashlib.sha256(json.dumps(definition,sort_keys=True).encode()).hexdigest()
        self.acceptance_loaded=False
        if require_checked:
            accepted=json.loads((HERE/RECEIPT).read_bytes());_verify_hashes(accepted)
            if not accepted['all_passed'] or not all(accepted[k] for k in GATES) or any(accepted[k] for k in OPEN):
                raise ValueError('Modified histories receipt exceeds source scope')
            if tuple(accepted[k] for k in ('actual_five_defect_family_sha256','implicit_source_sha256','datum_enclosure_sha256')) \
                    !=(self.family,self.source,self.datum_sha):raise ValueError('Foreign modified history receipt')
            if accepted['modified_source_definition_sha256']!=self.modified_source_definition_sha256:
                raise ValueError('Different modified profile/history source definition')
            if (accepted['finite_integer_N'],accepted['signed_range_subdivisions_per_fixed_piece'])!=(self.N,self.cells):
                raise ValueError('Exact worked frequency/range source differs')
            raw=json.loads((HERE/NAME).read_bytes())
            if raw['frequency_uniform_O1_over_N_constants']!=encode(pack(self.bounds)):raise ValueError('Uniform source bounds differ')
            if raw['exact_modulated_history_pressure_radial_theorem']!=encode(pack(self.theorem)):
                raise ValueError('Actual modified history source theorem differs')
            self.hashes.update(accepted['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.acceptance_loaded=True

    @source_precision
    def scalars(self,t):
        t=self.ctx.mpf(t);lo,hi=endpoints(t)
        if not all(mp.isfinite(x) for x in (lo,hi)) or lo<-11 or hi>3:
            raise ValueError('This local history API supports log(R/Rd) in[-11,3]')
        key=t._mpi_
        if key not in self.cache:self.cache[key]=cumulative_scalar_enclosures(self.ctx,self.mu,self.N,t,self.cells)
        return dict(self.cache[key])

    @source_precision
    def history(self,region,Z,coordinate):
        self.candidate.direction.assert_graph();c=self.ctx;v=c.mpf(coordinate);zv=c.mpf(Z)
        if not all(mp.isfinite(x) for x in endpoints(v)+endpoints(zv)) or endpoints(zv)[0]<-1 or endpoints(zv)[1]>1:
            raise ValueError('Finite original coordinates and Z[-1,1] required')
        lo,hi=endpoints(v)
        if region=='O2_buffer' and 0<=lo<=hi<=11:
            pre=self.pre.axial(zv,buffer_offset=v);t=v-11
            local=self.candidate.profile(region,zv,v,self.N)
        elif region=='O3_slope_mu' and 0<=lo<=hi<=1:
            pre=self.pre.slope_mu(zv,v);t=v
            local=self.candidate.profile(region,zv,v,self.N)
        elif region=='quiet_O3_power' and 0<=lo<=hi<=2:
            Tw=c.mpf(self.pre.params.Tw)
            if endpoints(Tw)[0]<=2:raise ValueError('Quiet power repair interval does not fit source')
            pre=self.pre.power(zv,v/Tw);t=1+v;local=None
        else:raise ValueError('Original buffer/transition or quiet power offset[0,2] required')
        original=raw_pre_velocity_rows(c,pre)
        if any(any(endpoints(x)!=(0,0) for x in row.coefficients) for row in original['axial']):
            raise ValueError('The source before modulation must have exactly zero axial velocity')
        if local is None:
            Unew=original['theta'];Vnew=[row*0 for row in Unew];du=[row*0 for row in Unew]
        else:
            Unew=local['modified_theta_over_Pstar_ordinary_logR_axial5']
            Vnew=local['modified_axial_over_Pstar_ordinary_logR_axial5']
            du=local['swirl_increment_over_Pstar_ordinary_logR_axial5']
        scalars=self.scalars(t);z=IntervalTaylor.variable(c,zv,5)
        delta_rows,Q=increment_rows(c,z,self.delta,self.Ua,scalars,Unew,original['theta'],Vnew,du)
        old={name:[copy_jet(c,row) for row in pre['actual_normalized_primitive_y_derivative_axial5'][name]] for name in SHAPES}
        datum=IntervalTaylor(c,[c.mpf(endpoints(x)) for x in pre['original_P0_axial5_coefficients']])
        original_P=[old['p'][0]+datum]+old['p'][1:]
        pressure=[a+b for a,b in zip(original_P,delta_rows['p'])]
        radial_increment=shifted_rows(Q,c.mpf('.5'),4)
        histories={}
        for name in SHAPES:
            if name in ('m','k'):
                histories[name]=[dict(Pstar_power=0,ordinary_logR_rows=old[name]),
                  dict(Pstar_power=1,ordinary_logR_rows=delta_rows[name])]
            else:histories[name]=[dict(Pstar_power=0,ordinary_logR_rows=[a+b for a,b in zip(old[name],delta_rows[name])])]
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
          datum_enclosure_sha256=self.datum_sha,modified_source_definition_sha256=self.modified_source_definition_sha256,
          region=region,coordinate=v,Z=zv,actual_shared_logR_offset=t,
          finite_integer_N=self.N,current_original_pre_source=pre,
          exact_source_Pstar_log=c.mpf(self.pre.params.logPstar),
          signed_actual_cumulative_scalar_enclosures=scalars,
          exact_cumulative_integral_definitions=self.theorem['scalar_integral_definitions'],
          modified_five_histories_in_original_normalized_units_Pstar_sectors=histories,
          normalized_units=dict(m='Mz/R',h='Mtheta/(sqrt2*R^1.5*Pstar)',
            k='Mtheta_z/(sqrt2*R^1.5*Pstar)',e='Mztheta/(R*Pstar^2)',p='Mp/Pstar^2'),
          exact_modified_history_increment_rows=delta_rows,
          same_original_axis_pressure_datum_axial5=datum,
          original_absolute_pressure_over_Pstar2_ordinary_logR_rows=original_P,
          modified_absolute_pressure_over_Pstar2_ordinary_logR_rows=pressure,
          pressure_correction_from_own_cumulative_Cp_ordinary_logR_rows=delta_rows['p'],
          modified_cylindrical_velocity_source_log_sectors=dict(
            theta=[dict(Pstar_power=1,ordinary_logR_rows=Unew)],
            axial=[dict(Pstar_power=1,ordinary_logR_rows=Vnew)],
            radial=[dict(Pstar_power=0,ordinary_logR_rows=original['radial']),
                    dict(Pstar_power=1,ordinary_logR_rows=radial_increment)]),
          radial_common_prefactor='sqrt(current_R/2)',original_incoming_meridional_history_retained=True,
          full_variable_profile_jets_and_new_own_moment_derivatives_retained=True,
          same_original_recovery_equation_preserves_exact_divergence_as_function=True,
          same_axis_pressure_datum_preserved=True,
          exact_zero_increment_before_original_support=endpoints(t)[1]<=-2,
          after_support_profile_returns_but_cumulative_histories_still_transport=endpoints(t)[0]>=.5,
          terminal_five_moment_restore_not_performed=True,
          source_enclosures_not_resolved_point_values=True,
          **{k:self.acceptance_loaded for k in GATES},**dict.fromkeys(OPEN,False))

    def manifest(self):
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
          datum_enclosure_sha256=self.datum_sha,modified_source_definition_sha256=self.modified_source_definition_sha256,
          finite_integer_N=self.N,
          signed_range_subdivisions_per_fixed_piece=self.cells,
          exact_modulated_history_pressure_radial_theorem=self.theorem,
          frequency_uniform_O1_over_N_constants=self.bounds,
          complete_new_history_source_views=VIEWS_NAME,
          current_strict_original_nonzero_whole_regions=15,remaining_original_whole_regions_without_cone=17,
          scope='Actual changed cumulative five histories as source functions on buffer/transition/quiet power, full-Z shapes and ordinary logR/Z jets, shared-axis absolute pressure and radial velocity from own M. Correct relative Pstar sectors preserved. Signed cell range-integration and uniform O(1/N) axial bounds are enclosures of the exact oscillatory integrals, not cycle-resolving samples. Quiet power still carries transported defects. Independent five-bump restoration, admitted uniform N, modified complete tensor/joins/cone, preheat exterior compatibility, actual recursion/waves, flatness and corrected NS/energy remain open.',
          input_hashes=self.hashes,**dict.fromkeys(GATES+OPEN,False))

@source_precision
def run(field=None):
    field=field if field is not None else CurrentO3ModulatedHistories(require_checked=False)
    result=field.manifest()
    views={key:field.history(region,(-1,1),coordinate) for key,region,coordinate in (
      ('buffer_whole_support','O2_buffer',(9,11)),
      ('transition_whole','O3_slope_mu',(0,1)),
      ('quiet_power_repair_preview','quiet_O3_power',(1,2)))}
    payload=(json.dumps(encode(pack(views)),separators=(',',':'))+'\n').encode()
    (HERE/VIEWS_NAME).write_bytes(gzip.compress(payload,compresslevel=6,mtime=0))
    result['input_hashes'][VIEWS_NAME]=sha(VIEWS_NAME)
    (HERE/NAME).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Actual changed five histories, same-axis pressure and own radial recovery constructed; repair/cone remain open',flush=True)
    return result
