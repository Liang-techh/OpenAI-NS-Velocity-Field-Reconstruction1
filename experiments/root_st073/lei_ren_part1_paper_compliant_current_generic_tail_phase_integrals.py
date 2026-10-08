"""Original chart-generic true-phase density oracle and directed C1 integrals.

Uses actual N*y, conditioned inverse and original first jets, without whole
period primitive caps. The producer evaluates one actual buffer period at two
resolutions and records full-domain source branch refusals explicitly.
"""
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import time
import lei_ren_part1_paper_compliant_current_weighted_O2_to_Rc_tail as tail

common=tail.common;rc=tail.rc;HERE,PREFIX,sha=tail.HERE,tail.PREFIX,tail.sha
ep,iv,encode=tail.ep,tail.iv,tail.encode;N=tail.N;KEYS=tail.KEYS
NAME=PREFIX+'current_generic_tail_phase_integrals.json'
RECEIPT=PREFIX+'current_generic_tail_phase_integrals_check.json'
GATE='original_generic_tail_true_phase_density_oracle_and_actual_buffer_period_C1_integrals_executed'
CHARTS=('O2_axial','O2_buffer','O3_slope_mu','O3_power')


class GenericTailPhaseIntegrals:
    def __init__(self,owner):
        if type(owner) is not rc.NativeRcC1Histories:raise ValueError('Original Rc source owner required')
        self.owner=owner;self.transfer=owner.transfer;self.coordinates=owner.coordinates
        self.density_owner=self.transfer.owner.owner
        if type(self.density_owner) is not rc.density.NativeDensityC1LocalIntegrals:
            raise ValueError('Original source/phase/first-jet density owner required')
        self.c=owner.ctx;self.family=owner.family

    def geometry(self,chart,left,right):
        if chart not in CHARTS:raise ValueError('Original tail chart required')
        if chart=='O3_power':
            return self.transfer.geometry.cell(chart,dict(original_power_offset=str(left)),dict(original_power_offset=str(right)))
        return self.transfer.geometry.cell(chart,str(left),str(right))

    def query(self,chart,left,right,N):
        N=rc.density.spatial.density.candidate_integer(N)
        if N<160:raise ValueError('Original candidate N>=160 required')
        geometry=self.geometry(chart,left,right)
        got=self.density_owner.spatial_query(chart,self.Z,geometry['coordinate'],N)
        statuses=[cell['record']['status'] for cell in got['cells']]
        enclosed=all(cell['values'] is not None for cell in got['cells'])
        record=dict(chart=chart,exact_native_endpoints=[str(left),str(right)],candidate_N=N,
            source_family=self.family,exact_Z_range=list(self.Z),actual_true_geometry=geometry['record'],
            actual_original_phase_inverse_and_density_source=got['record'],
            phase_cell_statuses=statuses,status='enclosed' if enclosed else 'requires_original_source_phase_subdivision',
            actual_N_log_radius_phase_not_independent_phase_samples=True,
            whole_period_primitive_caps_not_used_as_A_B=True,
            original_ordinary_Z_density_derivatives_not_cap_derivatives=True,
            original_source_query_encloses_functions_not_selected_midpoints=True)
        if not enclosed:return dict(record=record,geometry=geometry,values=None)
        kernels={k:rc.transfer.local.same_source_union([cell['values']['kernels'][k] for cell in got['cells']]) for k in KEYS}
        jets={k:rc.transfer.local.same_source_union([cell['values']['Z_derivatives'][k] for cell in got['cells']]) for k in KEYS}
        factors={k:rc.transfer.true_width_kernel(self.coordinates,geometry,rate) for k,rate in rc.RATES.items()}
        values={k:self.coordinates.rebase(kernels[k],self.family)*factors[k]['mass'] for k in KEYS}
        derivatives={k:self.coordinates.rebase(jets[k],self.family)*factors[k]['mass'] for k in KEYS}
        record.update(original_native_source_log_bases=got['source']['roots']['E'][(0,0)].scale.bases,
            full_phase_union_signed_C0_density={k:v.record() for k,v in kernels.items()},
            full_phase_union_signed_Z_density={k:v.record() for k,v in jets.items()},
            true_positive_kernel_factors={k:dict(branch=f['branch'],mass=f['mass'].record(),decay=f['decay'].record()) for k,f in factors.items()},
            original_local_C0_integral={k:v.record() for k,v in values.items()},
            original_local_Z_integral={k:v.record() for k,v in derivatives.items()},
            true_log_radius_Jacobian_integrated_once=True,
            pressure_P0_not_added_to_density_or_integral=True)
        return dict(record=record,geometry=geometry,values=values,Z_derivatives=derivatives)

    def integrate(self,chart,left,right,cells,N):
        lo,hi=Fraction(left),Fraction(right)
        if lo>=hi or type(cells) is not int or cells<1:raise ValueError('Ordered exact nonempty partition required')
        operator=rc.history.C1DuhamelOperator(self.coordinates);rows=[]
        for index in range(cells):
            a=lo+(hi-lo)*index/cells;b=lo+(hi-lo)*(index+1)/cells
            got=self.query(chart,a,b,N);rows.append(got['record'])
            if got['values'] is None:
                return dict(record=dict(status='requires_original_source_phase_subdivision',chart=chart,
                    exact_native_endpoints=[str(lo),str(hi)],attempted_partition_cells=cells,
                    first_unresolved_cell=index,queried_original_cells=rows,
                    no_partial_integral_substituted_for_complete_window=True),values=None)
            rc.transfer.append_true_cell(operator,got['geometry'],got['values'],got['Z_derivatives'],self.family)
        return dict(record=dict(status='enclosed',source_family=self.family,exact_Z_range=list(self.Z),chart=chart,
            exact_native_endpoints=[str(lo),str(hi)],candidate_N=N,ordered_source_cells=cells,
            queried_original_cells=rows,original_C1_integral_operator=operator.record(),
            five_original_C0_integrals={k:v.record() for k,v in operator.increments.items()},
            five_original_Z_integrals={k:v.record() for k,v in operator.Z_increments.items()},
            original_interval_value_and_Z_functions_evaluated_not_caps_selected=True,
            zero_operator_additive_identity_not_zero_actual_incoming=True,
            actual_source_phase_radius_Jacobian_and_density_graph_retained=True,
            numerical_local_integral_enclosure_not_exact_field_point_or_whole_chart_closure=True),
            values=operator.increments,Z_derivatives=operator.Z_increments)


@rc.native.inlet.source_precision
def run():
    begin=time.monotonic();prior=json.loads((HERE/tail.NAME).read_bytes());hashes={}
    common.attach_receipt(hashes,tail,prior['source_family'])
    bridge,seed=rc.native.inlet.native_bridge_owner()
    with rc.native.inlet.CheckedSourceRuntime():
        owner=rc.NativeRcC1Histories(rc.preceding.NativeO2C1Histories(tail.baseline.build_owner(bridge)))
        oracle=GenericTailPhaseIntegrals(owner);archives=[];summaries=[]
        for tag,Z in (('positive',('.36','.38')),('negative',('-.38','-.36'))):
            oracle.Z=Z;queries=[];integrals=[]
            for label,chart,left,right in (('axial','O2_axial',0,1),('buffer_0_9','O2_buffer',0,9),
                ('buffer_9_11','O2_buffer',9,11),('transition','O3_slope_mu',0,1),
                ('power_to_r_plus','O3_power',0,1),('power_to_Rc','O3_power',1,2)):
                got=oracle.query(chart,Fraction(left),Fraction(right),N)
                got['record']['original_tail_label']=label;queries.append(got['record'])
                print('Genuine original tail phase/density query:',tag,label,got['record']['status'],flush=True)
            for count in (8,32):
                got=oracle.integrate('O2_buffer',Fraction(5),Fraction(5)+Fraction(1,N),count,N)
                if got['values'] is None:raise ArithmeticError('Refine actual buffer period source before claiming integral')
                integrals.append(got['record']);print('Genuine actual buffer-period C1 integral:',tag,count,flush=True)
            payload=dict(source_family=owner.family,candidate_N=N,exact_Z_range=list(Z),
                full_original_tail_domain_phase_density_queries=queries,actual_buffer_period_integral_refinements=integrals,
                original_period_window='buffer selector[5,5+1/N], true log-radius length1/N; actual start phase retained',
                actual_integer_N_phase_full_period_not_fixed_free_phase=True,
                full_axial_transition_or_complete_tail_numerical_integrals_admitted=False)
            raw=json.dumps(encode(payload),indent=2).encode()+b'\n';compressed=gzip.compress(raw,compresslevel=9,mtime=0)
            filename=PREFIX+'current_generic_tail_phase_integrals_'+tag+'.json.gz'
            if len(compressed)>=100*1024*1024:raise ValueError('Split large lossless numerical evidence')
            (HERE/filename).write_bytes(compressed);hashes[filename]=sha(filename)
            archives.append(dict(filename=filename,exact_Z_range=list(Z),compressed_bytes=len(compressed),uncompressed_bytes=len(raw),
                lossless_original_json_sha256=hashlib.sha256(raw).hexdigest()))
            summaries.append(dict(exact_Z_range=list(Z),full_domain_status={q['original_tail_label']:q['status'] for q in queries},
                actual_buffer_period_counts=[8,32],actual_buffer_period_C0_Z_integrals=[
                    dict(ordered_source_cells=q['ordered_source_cells'],C0=q['five_original_C0_integrals'],Z=q['five_original_Z_integrals']) for q in integrals]))
        for name,digest in owner.service.hashes.items():
            if name in hashes and hashes[name]!=digest:raise ValueError('Original source closures disagree: '+name)
            hashes[name]=digest
        hashes[Path(__file__).name]=sha(Path(__file__).name)
        result=dict(**{GATE:True},source_family=owner.family,candidate_N=N,
            actual_original_tail_phase_density_source_archives=archives,actual_tail_oracle_and_period_summaries=summaries,
            source_seed_evidence=seed,input_hashes=hashes,execution_seconds=time.monotonic()-begin,
            same_original_generic_chart_phase_inverse_and_density_value_Z_interface_executed=True,
            whole_period_primitive_caps_not_used_as_numerical_values=True,
            numerical_complete_tail_integrals_or_Rc_targets_or_closure_admitted=False,
            full_Rc_closure_or_global_N_or_stress_or_recursion_claimed=False,**dict.fromkeys(common.current.FLAGS,False),
            scope='Original chart-generic spatial source/phase/first-jet density evaluation on both strict-sign tiles; true buffer period[5,5+1/N]8/32-cell C0/Z directed integrals atN1024. Unresolved whole source branches explicit. Not complete tail integrals, Rc defects/controls/global N/stress/recursion/full NS.')
        (HERE/NAME).write_text(json.dumps(encode(result),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()
