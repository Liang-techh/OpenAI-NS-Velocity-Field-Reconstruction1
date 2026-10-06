"""Callable own repaired five histories, pressure and radial source at O3 exit."""
import gzip
import hashlib
import json
from lei_ren_part1_paper_compliant_current_O3_independent_repair import (
    CurrentO3IndependentRepair,HERE,PREFIX,NAME as REPAIR_NAME,RECEIPT as REPAIR_RECEIPT,
    sha,pack,encode,endpoints,source_precision,_verify_hashes)
from lei_ren_part1_paper_compliant_current_O3_repaired_histories_operator import (
    actual_repaired_history_theorem,bump_rows,partial_weights,partial_repair_primitives)
from lei_ren_part1_paper_compliant_current_O3_modulated_histories_operator import SHAPES,increment_rows
from lei_ren_part1_paper_compliant_current_pre_pulse_stress_operator import raw_pre_velocity_rows,copy_jet
from lei_ren_part1_paper_compliant_collar_stress_C3 import shifted_rows
from lei_ren_part1_paper_compliant_frozen_comparison_field import square
from lei_ren_part1_paper_interval_taylor import IntervalTaylor

NAME=PREFIX+'current_O3_repaired_histories.json'
RECEIPT=PREFIX+'current_O3_repaired_histories_check.json'
VIEWS_NAME=PREFIX+'current_O3_repaired_histories_views.json.gz'
GATES=('current_modified_actual_bump_logR4_axial5_source_available',
       'current_modified_repaired_partial_histories_pressure_radial_source_installed',
       'current_O3_transition_modified_five_moment_repair_certified',
       'current_modified_quiet_power_functional_exit_join_to_original_source_certified')
OPEN=('current_O3_transition_finite_N_modified_profiles_installed',
      'current_whole_O3_transition_signed_two_vector_cone_certified',
      'current_O2_buffer_modified_whole_cone_certified',
      'current_modified_complete_tensor_and_physical_source_joins_certified',
      'current_modified_analytic_preheat_exterior_compatible_after_repair',
      'current_modified_global_admissible_tensor_and_full_NS_certified')

class CurrentO3RepairedHistories:
    @source_precision
    def __init__(self,repair=None,partial_cells=128,require_checked=True):
        self.repair=repair if repair is not None else CurrentO3IndependentRepair()
        if type(self.repair) is not CurrentO3IndependentRepair or not self.repair.acceptance_loaded:
            raise ValueError('Checked same-source independent implicit repair required')
        self.histories=self.repair.histories;self.histories.candidate.direction.assert_graph()
        self.ctx=c=self.repair.ctx;self.pre=self.histories.pre;self.mu=self.repair.mu;self.N=self.repair.N
        if isinstance(partial_cells,bool) or not isinstance(partial_cells,int) or partial_cells<16 or partial_cells>1024 or partial_cells&(partial_cells-1):
            raise ValueError('Power-of-two partial subdivisions in[16,1024] required')
        self.partial_cells=partial_cells;self.f2=c.exp(-1-c.mpf('1.5')*self.mu)
        self.theorem=actual_repaired_history_theorem();self.cache={}
        source_key=self.theorem['checked_parent_D_source_identity_required']
        if not self.repair.theorem['passed'] or self.repair.theorem['identities'].get(source_key) is not True:
            raise ValueError('Checked actual inlet J/M correlation is required for functional exit')
        if not self.repair.certificate['unique_exact_implicit_controls_exist'] or not self.repair.certificate['no_source_defect_midpoint_or_zero_substitution']:
            raise ValueError('Exit requires the same exact implicit root of the actual source defects')
        self.hashes={**self.repair.hashes,**self.theorem['input_hashes'],
          REPAIR_NAME:sha(REPAIR_NAME),REPAIR_RECEIPT:sha(REPAIR_RECEIPT)}
        for stem in ('current_O3_repaired_histories_operator','current_O3_repaired_histories'):
            name=PREFIX+stem+'.py';self.hashes[name]=sha(name)
        definition=dict(parent_repair_definition=self.repair.repair_definition_sha256,
          finite_integer_N=self.N,partial_cells=self.partial_cells,
          operator_sha256=sha(PREFIX+'current_O3_repaired_histories_operator.py'),
          source_sha256=sha(PREFIX+'current_O3_repaired_histories.py'))
        self.modified_source_definition_sha256=hashlib.sha256(json.dumps(definition,sort_keys=True).encode()).hexdigest()
        self.acceptance_loaded=False
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes());_verify_hashes(receipt)
            if not receipt['all_passed'] or not all(receipt[k] for k in GATES) or any(receipt[k] for k in OPEN):
                raise ValueError('Repaired source receipt exceeds its scoped field construction')
            if receipt['modified_source_definition_sha256']!=self.modified_source_definition_sha256:
                raise ValueError('Foreign repaired source or implicit control vector')
            raw=json.loads((HERE/NAME).read_bytes())
            expected=self.manifest();expected['input_hashes'][VIEWS_NAME]=sha(VIEWS_NAME)
            if encode(pack(expected))!=raw:raise ValueError('Actual repaired field manifest differs')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.acceptance_loaded=True

    @source_precision
    def partial(self,y):
        y=self.ctx.mpf(y);key=y._mpi_
        if key not in self.cache:
            W=partial_weights(self.repair,y,self.partial_cells)
            values=partial_repair_primitives(self.ctx,self.mu,self.N,self.repair.controls,W)
            self.cache[key]=(W,values)
        return self.cache[key]

    @source_precision
    def history(self,region,Z,coordinate):
        parent=self.histories.history(region,Z,coordinate);c=self.ctx
        result=dict(parent)
        result.update(parent_modified_source_definition_sha256=self.histories.modified_source_definition_sha256,
          modified_source_definition_sha256=self.modified_source_definition_sha256,
          same_unique_implicit_repair_definition_sha256=self.repair.repair_definition_sha256,
          same_exact_control_vector_for_all_queries=True,
          terminal_five_moment_restore_not_performed=False,
          **{k:self.acceptance_loaded for k in GATES},**dict.fromkeys(OPEN,False))
        if region!='quiet_O3_power':return result
        q=c.mpf(coordinate);y=q-1;x=c.exp(y);z=IntervalTaylor.variable(c,c.mpf(Z),5)
        original=raw_pre_velocity_rows(c,parent['current_original_pre_source'])
        Ahat=(1+square(z)).reciprocal()*(self.histories.Ua*self.f2)
        rows=[bump_rows(c,y,center,self.repair.weights['radius'],self.repair.normalization)
          for center in self.repair.weights['centers']]
        h=self.repair.controls;sa=c.sqrt(self.mu)/self.N;se=self.mu/self.N
        du=[Ahat*(se*sum((h[i+2]*rows[i][n] for i in range(3)),c.mpf(0))) for n in range(5)]
        Vnew=[Ahat*(sa*(h[0]*rows[0][n]+h[1]*rows[2][n])) for n in range(5)]
        Unew=[a+b for a,b in zip(original['theta'],du)]
        W,primitives=self.partial(y)
        previous=parent['signed_actual_cumulative_scalar_enclosures']
        scalar=dict(m=previous['m']+self.f2*primitives['M']/x,
          h=previous['h']+self.f2*primitives['I']/x**c.mpf('1.5'),
          k=previous['k']+self.f2**2*primitives['J']/x**c.mpf('1.5'),
          e=previous['e']+self.f2**2*primitives['S']/x,
          p=previous['p']+self.f2**2*primitives['Cp'])
        exit_identity=endpoints(q)[0]>=2
        if exit_identity:
            # This is a source-level algebraic reduction of the very same
            # implicit map, only after all supports have ended. The signed
            # enclosure sums above are preserved as diagnostic evidence.
            if not self.theorem['exact_exit_reduction_requires_all_full_support_weights_and_checked_defining_equation']:
                raise ValueError('Exact terminal source identity is required')
            if any(any(endpoints(v)!=(0,0) for v in row.coefficients) for row in du+Vnew):
                raise ValueError('Exact exit cannot precede flat profile support end')
            if not self.repair.acceptance_loaded:raise ValueError('Unbound implicit control vector at source exit')
            scalar=dict.fromkeys(SHAPES,c.mpf(0))
        delta_rows,Q=increment_rows(c,z,self.histories.delta,self.histories.Ua,scalar,Unew,original['theta'],Vnew,du)
        old=parent['current_original_pre_source']['actual_normalized_primitive_y_derivative_axial5']
        new_histories={}
        for name in SHAPES:
            if name in ('m','k'):
                new_histories[name]=[dict(Pstar_power=0,ordinary_logR_rows=[copy_jet(c,row) for row in old[name]]),
                  dict(Pstar_power=1,ordinary_logR_rows=delta_rows[name])]
            else:new_histories[name]=[dict(Pstar_power=0,ordinary_logR_rows=[copy_jet(c,a)+b for a,b in zip(old[name],delta_rows[name])])]
        original_pressure=parent['original_absolute_pressure_over_Pstar2_ordinary_logR_rows']
        result.update(actual_repair_local_logx=y,actual_repair_radius_over_R0=x,
          current_independent_bump_logR4_rows=rows,
          current_same_source_partial_bump_weights=W,
          current_same_source_partial_repair_primitive_changes=primitives,
          original_unrepaired_modulation_scalar_enclosures=previous,
          unreduced_interval_sum_before_exact_exit_rewrite=dict(m=previous['m']+self.f2*primitives['M']/x,
            h=previous['h']+self.f2*primitives['I']/x**c.mpf('1.5'),
            k=previous['k']+self.f2**2*primitives['J']/x**c.mpf('1.5'),
            e=previous['e']+self.f2**2*primitives['S']/x,p=previous['p']+self.f2**2*primitives['Cp']),
          signed_actual_cumulative_scalar_enclosures=scalar,
          actual_positive_modulation_kinetic_history_retained=previous['e_kinetic'],
          actual_nonnegative_repair_kinetic_primitive=sa**2*(h[0]**2*W['energy'][0]+h[1]**2*W['energy'][2]),
          signed_energy_exit_cancellation_does_not_remove_kinetic_source=True,
          modified_five_histories_in_original_normalized_units_Pstar_sectors=new_histories,
          exact_modified_history_increment_rows=delta_rows,
          modified_absolute_pressure_over_Pstar2_ordinary_logR_rows=[a+b for a,b in zip(original_pressure,delta_rows['p'])],
          pressure_correction_from_own_cumulative_Cp_ordinary_logR_rows=delta_rows['p'],
          modified_cylindrical_velocity_source_log_sectors=dict(
            theta=[dict(Pstar_power=1,ordinary_logR_rows=Unew)],
            axial=[dict(Pstar_power=1,ordinary_logR_rows=Vnew)],
            radial=[dict(Pstar_power=0,ordinary_logR_rows=original['radial']),
                    dict(Pstar_power=1,ordinary_logR_rows=shifted_rows(Q,c.mpf('.5'),4))]),
          same_original_axis_datum_retained_in_repaired_pressure=True,
          new_radial_velocity_uses_own_repaired_moments=True,
          exact_functional_exit_reduction_applied=exit_identity,
          all_five_deltas_pressure_and_radial_rows_equal_original_at_exit=exit_identity,
          completed_modified_tensor_join_and_common_cone_N_remain_open=True)
        return result

    def manifest(self):
        return dict(actual_five_defect_family_sha256=self.histories.family,
          implicit_source_sha256=self.histories.source,datum_enclosure_sha256=self.histories.datum_sha,
          parent_repair_definition_sha256=self.repair.repair_definition_sha256,
          modified_source_definition_sha256=self.modified_source_definition_sha256,
          finite_integer_N=self.N,partial_range_subdivisions=self.partial_cells,
          exact_repaired_partial_history_and_exit_source_theorem=self.theorem,
          complete_repaired_source_views=VIEWS_NAME,
          current_strict_original_nonzero_whole_regions=15,remaining_original_whole_regions_without_cone=17,
          scope='Actual independent fixed-bump logR4/Z5 profiles and continuous partial weights, own five cumulative moments with all physical/Pstar factors, same-axis pressure and own-moment radial source. Source-bound exact implicit controls give functional q=2 exit equality to original field/histories through retained derivative orders. Source values remain enclosures; no control midpoints or zeroed incoming moments. Modified complete tensor/physical joins, analytic preheat exterior after repair, common cone N, global stress, coefficient recursion/waves and full corrected NS/energy remain open.',
          input_hashes=self.hashes,**dict.fromkeys(GATES+OPEN,False))

@source_precision
def run(field=None):
    field=field if field is not None else CurrentO3RepairedHistories(require_checked=False)
    result=field.manifest()
    views={key:field.history('quiet_O3_power',(-1,1),q) for key,q in (
      ('whole_repair_band',(1,2)),('first_overlap_interior','1.201337'),('middle_swirl_interior','1.501337'),
      ('last_overlap_interior','1.801337'),('exact_exit',2))}
    payload=(json.dumps(encode(pack(views)),separators=(',',':'))+'\n').encode()
    (HERE/VIEWS_NAME).write_bytes(gzip.compress(payload,compresslevel=6,mtime=0))
    result['input_hashes'][VIEWS_NAME]=sha(VIEWS_NAME)
    (HERE/NAME).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Actual repaired profiles/partial histories/pressure/radial source and exact five-moment exit installed; modified cone open',flush=True)
    return result

if __name__=='__main__':run()
