"""Actual N1024 upstream corrections feed the complete original O2 driver.

The five original C0/Z histories use fresh full-predicate direct integrals.
A second complete history enclosure uses the accepted all-N spatial bound.
Original P0 stays separate and rate-zero pressure retains upstream memory.
"""
import gzip
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_O2_all_N_spatial_envelope as spatial
import lei_ren_part1_paper_compliant_current_original_O2_source_incoming_common_N as upstream

integrals=spatial.integrals;middle=upstream.middle;functions=upstream.functions
base,prior,ep=spatial.base,spatial.prior,spatial.ep
HERE,PREFIX,sha=spatial.HERE,spatial.PREFIX,spatial.sha
KEYS=integrals.KEYS;N=1024
NAME=PREFIX+'current_original_O2_full_predicate_incoming_driver.json.gz'
RECEIPT=PREFIX+'current_original_O2_full_predicate_incoming_driver_check.json'
GATE='original_N1024_actual_upstream_full_Z_correction_bound_to_full_predicate_O2_C0_Z_driver'


def accepted(hashes,module,family):
    receipt=json.loads((HERE/module.RECEIPT).read_bytes())
    data=(HERE/module.NAME).read_bytes()
    raw=gzip.decompress(data) if module.NAME.endswith('.gz') else data
    report=json.loads(raw)
    if not receipt.get('all_passed') or not receipt.get(module.GATE) or report.get('source_family')!=family:
        raise ValueError('Accepted same original source required: '+module.RECEIPT)
    if module.NAME.endswith('.gz') and hashlib.sha256(raw).hexdigest()!=receipt['compressed_producer_report']['lossless_original_json_sha256']:
        raise ValueError('Accepted source archive changed: '+module.NAME)
    closure={**receipt['input_hashes'],module.NAME:sha(module.NAME),module.RECEIPT:sha(module.RECEIPT)}
    for name,digest in closure.items():
        if sha(name)!=digest:raise ValueError('Original driver dependency changed: '+name)
        if name in hashes and hashes[name]!=digest:raise ValueError('Original driver source closures disagree: '+name)
        hashes[name]=digest
    return report


class OriginalO2FullPredicateIncomingDriver:
    def __init__(self):
        self.direct=integrals.OriginalO2PredicateFiveIntegrals()
        self.spatial=spatial.OriginalO2AllNSpatialEnvelope()
        self.ctx=c=self.direct.ctx;self.out=self.direct.out;self.family=self.direct.family
        if self.spatial.family!=self.family:raise ValueError('Same original O2 source family required')
        self.hashes=dict(self.direct.hashes)
        for name,digest in self.spatial.hashes.items():
            if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Original O2 owners disagree')
            self.hashes[name]=digest
        accepted(self.hashes,integrals,self.family);accepted(self.hashes,spatial,self.family)
        report=accepted(self.hashes,middle,self.family)
        self.exact_binding=upstream.source_function_prefix_binding(accepted(self.hashes,functions,self.family))
        if report['candidate_N']!=N or self.exact_binding['common_N_lower']>N:
            raise ValueError('Accepted actual upstream candidate N1024 required')
        self.incoming=report['actual_inlet_to_O2_inlet_C1_history_records']['whole_Z']
        incoming=self.incoming
        if incoming['candidate_N']!=N or incoming['source_family']!=self.family:
            raise ValueError('Same source-owned N1024 upstream correction required')
        if ep(integrals.interval(c,incoming['Z_box']))!=(-1,1) or not incoming['no_original_interval_skipped_between_true_inlet_and_O2_inlet']:
            raise ValueError('Complete actual upstream continuous whole-Z route required')
        background=incoming['actual_O2_inlet_original_background_and_separate_P0_Z']
        if background['original_units']!=integrals.five.UNITS or not background['P0_not_merged_into_pressure_history']:
            raise ValueError('Same five normalized units and separate P0 required')
        common=background['common_directed_coordinate_theorem']
        if common['source_family']!=self.family or not common['common_basis_is_arithmetic_not_a_new_source_function']:
            raise ValueError('Original shared directed coordinate contract required')
        old=common['common_log_bases']
        self.coordinates=middle.history.CommonSourceCoordinates(c,integrals.interval(c,old[1]),self.family)
        if any(ep(integrals.interval(c,old[i]))!=(0,0) for i in (0,2,3,4)):
            raise ValueError('Original common Pstar-squared basis required')
        if s.simplify(self.out.frame.definitions['logPstar']-(s.exp(40)+11))!=0:
            raise ValueError('Exact original Pstar definition required')
        with mp.workdps(c.dps+40):
            nativeP2=2*self.out.bases[0]
            oldP2=self.coordinates.logP_squared
            if max(ep(nativeP2)[0],ep(oldP2)[0])>min(ep(nativeP2)[1],ep(oldP2)[1]):
                raise ValueError('Original shared Pstar definitions disagree')
        self.require_output_atlas(self.out)
        self.require_output_atlas(self.spatial.out)
        self.raw_incoming=[{key:upstream.restore_common_source(value,self.coordinates) for key,value in incoming[name].items()}
            for name in ('actual_original_inlet_to_O2_inlet_correction_C0','actual_original_inlet_to_O2_inlet_correction_Z')]
        if any(set(group)!=set(KEYS) for group in self.raw_incoming):raise ValueError('All five source-owned C0/Z incoming corrections required')
        self.native_incoming=[{key:self.common_to_native(value) for key,value in group.items()} for group in self.raw_incoming]
        self.P0=self.common_to_native(upstream.restore_common_source(background['original_separate_P0_over_Pstar_squared'],self.coordinates))
        self.P0_Z=self.common_to_native(upstream.restore_common_source(background['original_separate_P0_Z_over_Pstar_squared'],self.coordinates))
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)

    def require_output_atlas(self,atlas):
        record=atlas.record()
        if record['source_family']!=self.family or record['exact_Z_bounds']!=['-1','1']:
            raise ValueError('Same whole-Z original output atlas required')
        if record['source_basis_order']!=['logPstar','selected_logCstar','logL','original_logq','zero']:
            raise ValueError('Original output basis semantics required')
        if ep(atlas.bases[3])!=(0,0) or ep(atlas.bases[4])!=(0,0):
            raise ValueError('Only collected whole-Z output values can be rebound')
        if atlas.frame.selected_logCstar_mpf_tuple!=self.out.frame.selected_logCstar_mpf_tuple or any(
            atlas.frame.definitions[key]!=self.out.frame.definitions[key] for key in ('logPstar','log_delta','logRref')):
            raise ValueError('Same exact original selected parameters required')
        if any(ep(left)!=ep(right) for left,right in zip(atlas.bases,self.out.bases,strict=True)):
            raise ValueError('Same directed original output bases required')

    def common_to_native(self,value):
        """Exact same-source Pstar^(2k) rebase, including full directed offsets."""
        if type(value) is not prior.ScaledEnclosure or value.scale.bases is not self.coordinates.bases or value.ledger is not self.coordinates.ledger:
            raise ValueError('Live accepted common-coordinate source cover required')
        powers=value.scale.powers
        if any(powers[i] for i in (0,2,3,4)):raise ValueError('Collected common Pstar-squared factor required')
        c=self.ctx
        return prior.ScaledEnclosure(prior.FormalScale(self.out.bases,(2*powers[1],0,0,0,0),c.mpf(ep(value.scale.offset))),
            c.mpf(ep(value.coefficient)),self.out.ledger)

    def spatial_to_native(self,value):
        if type(value) is not prior.ScaledEnclosure or value.scale.bases is not self.spatial.out.bases or value.ledger is not self.spatial.out.ledger:
            raise ValueError('Issued same-source spatial output cover required')
        if value.scale.powers not in (spatial.slow.UNITS[spatial.Y],spatial.slow.UNITS[spatial.YZ]):
            raise ValueError('Original full-spatial dominating units required')
        c=self.ctx
        return prior.ScaledEnclosure(prior.FormalScale(self.out.bases,value.scale.powers,c.mpf(ep(value.scale.offset))),
            c.mpf(ep(value.coefficient)),self.out.ledger)

    def add_native(self,left,right):
        """Collect incoming offsets as well as O2 powers before exponentials.

        The ceiling is an arithmetic normalization, never a source point or
        a Z-dependent field to differentiate. Both whole signed covers remain.
        """
        for value in (left,right):
            if type(value) is not prior.ScaledEnclosure or value.scale.bases is not self.out.bases or value.ledger is not self.out.ledger:
                raise ValueError('Same rebound original output basis and ledger required')
            if value.scale.powers[4]!=0:raise ValueError('Original actual radius must already be collected')
        if left.zero:return right
        if right.zero:return left
        c=self.ctx;lp,rp=left.scale.powers,right.scale.powers
        powers=(max(lp[0],rp[0]),max(lp[1],rp[1]),min(lp[2],rp[2]),min(lp[3],rp[3]),0)
        offset=c.mpf(max(ep(left.scale.offset)[1],ep(right.scale.offset)[1]))
        anchor=prior.FormalScale(self.out.bases,powers,offset)
        normalized=[value.coefficient*value.bounded_exp((value.scale-anchor).evaluate()) for value in (left,right)]
        self.out.ledger['incoming_directed_offset_and_power_collections']=self.out.ledger.get('incoming_directed_offset_and_power_collections',0)+1
        return prior.ScaledEnclosure(anchor,normalized[0]+normalized[1],self.out.ledger)

    def sum_native(self,values):
        result=self.out.scalar(0)
        for value in values:result=self.add_native(result,value)
        return result

    def compose(self,direct_result,spatial_result,*,N=1024):
        if type(N) is not int or N!=self.incoming['candidate_N']:
            raise ValueError('Do not reuse N1024 incoming for another candidate')
        if self.direct.issued_integrals.get(id(direct_result)) is not direct_result:
            raise ValueError('Issued full-predicate actual direct integral required')
        direct_record=direct_result['report']
        if direct_record['explicit_candidate_N']!=N or direct_record['source_family']!=self.family:
            raise ValueError('Same actual source candidate phase required')
        if direct_record['exact_y_window']!=['0','1'] or direct_record['exact_Z_range']!=['-1','1']:
            raise ValueError('Complete actual original O2 y/Z window required')
        candidate=self.spatial.at_candidate(spatial_result,N=N)
        c=self.ctx;alternate=[];correction=[];direct_correction=[];spatial_native=[];primary=[]
        with mp.workdps(c.dps+40):
            for j,(original,changed,values) in enumerate(zip((direct_result['original'],direct_result['originalZ']),
                (direct_result['changed'],direct_result['changedZ']),(candidate['C0'],candidate['ordinary_Z']),strict=True)):
                local={key:self.spatial_to_native(value) for key,value in values.items()};spatial_native.append(local)
                transported={key:self.native_incoming[j][key]*c.exp(-c.mpf(integrals.five.RATES[key])) for key in KEYS}
                corrected={key:self.add_native(transported[key],local[key]) for key in KEYS};correction.append(corrected)
                direct_correction.append({key:self.add_native(transported[key],changed[key]) for key in KEYS})
                alternate.append({key:self.add_native(original[key],corrected[key]) for key in KEYS})
                primary.append({key:self.sum_native((original[key],changed[key],transported[key])) for key in KEYS})
            pressures=[self.add_native(self.P0,primary[0]['p']),self.add_native(self.P0_Z,primary[1]['p'])]
            alternate_pressure=[self.add_native(self.P0,alternate[0]['p']),self.add_native(self.P0_Z,alternate[1]['p'])]
        encode=lambda group:{key:value.record() for key,value in group.items()}
        record=dict(source_family=self.family,explicit_candidate_N=N,exact_y_window=['0','1'],exact_Z_range=['-1','1'],
            normalized_own_units=integrals.five.UNITS,own_rates=integrals.five.RATES,
            original_P0_datum_sha256=self.family['datum_enclosure_sha256'],
            actual_upstream_source_binding=dict(manifest=middle.NAME,manifest_sha256=sha(middle.NAME),receipt=middle.RECEIPT,
                receipt_sha256=sha(middle.RECEIPT),source_record_key='whole_Z',candidate_N=N,
                exact_13_cell_upstream_function_binding=self.exact_binding,
                complete_original_C0_Z_correction_covers_not_original_plus_correction=True,
                same_source_Pstar_power_identity='(0,k,0,0,0) in common log(Pstar^2) -> (2k,0,0,0,0) in original O2 log(Pstar)',
                original_Pstar_definition='logPstar=exp(40)+11',common_coordinates=self.coordinates.record(),
                output_original_source_atlas=self.out.record(),full_signed_coefficients_and_directed_offsets_retained=True,
                no_native_log_radius_or_Pstar_exponential_materialized=True),
            actual_upstream_C0=encode(self.raw_incoming[0]),
            actual_upstream_Z=encode(self.raw_incoming[1]),
            rebound_source_owned_C0_Z_corrections=[encode(group) for group in self.native_incoming],
            full_predicate_direct_C0_Z_changed_contributions=[encode(direct_result[name]) for name in ('changed','changedZ')],
            source_defined_original_background_C0_Z=[encode(direct_result[name]) for name in ('original','originalZ')],
            propagated_actual_direct_C0_Z_corrections=[encode(group) for group in direct_correction],
            actual_complete_direct_C0_Z_histories=[encode(group) for group in primary],
            original_source_all_N_spatial_C0_Z_alternative_contributions=[encode(group) for group in spatial_native],
            propagated_actual_spatial_C0_Z_alternative_corrections=[encode(group) for group in correction],
            actual_complete_spatial_C0_Z_alternative_histories=[encode(group) for group in alternate],
            separate_analytic_P0=self.P0.record(),separate_analytic_P0_Z=self.P0_Z.record(),
            actual_absolute_direct_pressure_C0_Z=[value.record() for value in pressures],
            actual_absolute_spatial_alternative_pressure_C0_Z=[value.record() for value in alternate_pressure],
            actual_direct_radius_phase_endpoints=direct_record['actual_original_radius_phase_endpoints'],
            actual_spatial_candidate_radius_phase_endpoints=candidate['record']['actual_original_radius_common_N_phase_endpoints'],
            source_incoming_driver_full_Z_including_axis_and_both_signs=True,
            downstream_fresh_original_common_N_phase_not_old_N7_archive=True,
            upstream_error_covers_not_sharpened_by_local_O2_frequency_bound=True,
            alternative_envelopes_not_summed_as_two_distinct_physical_contributions=True,
            ordinary_Z_rows_translated_without_differentiating_normalization_or_caps=True,
            incoming_positive_offsets_collected_in_arithmetic_normalization_not_exponentiated=True,
            normalized_offset_ceiling_is_not_a_source_field_point_or_differentiated_function=True,
            original_background_added_once_and_analytic_P0_separate=True,
            exact_rate_zero_pressure_keeps_all_upstream_memory=True,
            O2_to_Rc_route_still_open=True,actual_all_route_incoming_histories_installed=False,
            current_whole_N_selected=False,functional_terminal_identity_solved=False)
        return dict(C0=primary[0],ordinary_Z=primary[1],record=record,
            correction=dict(values=direct_correction[0],Z_derivatives=direct_correction[1],source_family=self.family,candidate_N=N),
            spatial_alternative_correction=dict(values=correction[0],Z_derivatives=correction[1],source_family=self.family,candidate_N=N))

    def integrate(self,*,N=1024):
        if type(N) is not int or N!=1024:raise ValueError('Only accepted actual upstream N1024 is bound')
        direct=self.direct.integrate(N=N);spatial_result=self.spatial.integrate()
        return self.compose(direct,spatial_result,N=N),direct,spatial_result


def run():
    began=time.monotonic();owner=OriginalO2FullPredicateIncomingDriver();result,direct,bound=owner.integrate()
    report=dict(**{GATE:True},source_family=owner.family,candidate_N=N,
        actual_full_predicate_original_O2_source_incoming_driver=result['record'],
        actual_fresh_original_N1024_direct_integrals=direct['report'],
        full_spatial_all_N_coefficients=bound['report']['normalized_full_spatial_C0_Z_Nminus2_coefficient_caps'],
        actual_full_Z_source_owned_O2_incoming_driver_installed=True,
        actual_all_route_incoming_histories_installed=False,current_whole_N_selected=False,
        functional_terminal_identity_solved=False,all_17_chart_or_24_cell_oracle_installed=False,
        actual_five_controls_installed=False,actual_phase_primitive_point_evaluator_installed=False,
        **dict.fromkeys(spatial.source.ordered.base.point.source.inertial.profiles.loop.OPEN,False),
        input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
        scope='Actual same-source N1024 inlet-to-O2 corrections feed fresh complete original O2 y[0,1]/Z[-1,1] full-predicate C0/Z direct integrals and alternative full all-N spatial contribution bounds. Exact source basis adapter, original background once, separate analytic P0 and rate-zero memory. Upstream covers remain broad. Not O2-to-Rc/all-route closure, terminal controls/global N, scale recursion or corrected NS.')
    (HERE/NAME).write_bytes(gzip.compress(json.dumps(base.encoded(report),indent=2).encode()+b'\n',compresslevel=9,mtime=0))
    print('Actual N1024 upstream to full-predicate O2 C0/Z driver generated',flush=True)
    return report


if __name__=='__main__':run()
