"""Actual source-bound Rm incoming consumed by the complete local Rm patch.

Hydrate accepted finite-N prefix and complete all-u/weighted source rows.
Rebind them to one live source algebra before own-rate transport to Rh.
No old source constructor, integration producer, or whole checker is run.
"""
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_micro_finite_N as prefix
import lei_ren_part1_paper_compliant_current_original_Rm_weighted_averaging as weighted

fields,base,ep=prefix.fields,prefix.macro.base,prefix.micro.ep
switch=prefix.switch
HERE,PREFIX,sha=prefix.HERE,prefix.PREFIX,prefix.sha
RATES=prefix.RATES
NAME=PREFIX+'current_original_Rm_finite_N_prefix.json.gz'
RECEIPT=PREFIX+'current_original_Rm_finite_N_prefix_check.json'
GATE='current_original_source_bound_finite_N_prefix_consumed_through_Rh'
PHASE_RECIPE='frac(N*(log(110/4)+14logPstar+10logCstar+994+logx-hb*s_c/2))'


def canonical_saved(row):
    return {key:row[key] for key in ('coefficient_interval','formal_positive_scale','exact_zero')}


def row_packet(packet,key):
    rows=packet[key]
    if set(rows)!=set(RATES) or any(len(pair)!=2 for pair in rows.values()):
        raise ValueError('All five signed C0 and ordinary Z source rows required')
    if any(row.get('point_value_selected') is not False or row.get('encloses_original_source_function') is not True
           for pair in rows.values() for row in pair):
        raise ValueError('Actual source enclosures, not selected values, required')
    return rows


def guard_prefix(family,label,N,P0,packet):
    switch.guard_saved_source(family,label,N,P0,packet)
    if packet.get('genuine_current_finite_N_R0_R100_R110_Rm_boundary_enclosures_supplied') is not True:
        raise ValueError('Genuine current finite-N source prefix required')
    if packet.get('local_drivers_composed_with_real_current_prefix_not_arbitrary_zero') is not True:
        raise ValueError('Source-bound incoming, not arbitrary zero, required')
    row_packet(packet,'actual_current_Rm_incoming_correction_C0_Z')


def guard_local(family,label,N,P0,Rm_factor,packet,phase_points):
    switch.guard_saved_source(family,label,N,P0,packet)
    expected=base.encoded(weighted.PARTITION)
    if packet.get('exact_partition')!=expected:
        raise ValueError('Complete accepted original Rm-to-Rh partition required')
    row_packet(packet,'actual_averaged_whole_patch_local_integral_C0_Z')
    cells=packet['actual_source_weighted_cells']
    if len(cells)!=len(weighted.PARTITION)-1:
        raise ValueError('Every original weighted source cell required')
    radius=switch.canonical_expression([Rm_factor])[0]
    pressure=switch.canonical_expression(P0)
    for cell,left,right in zip(cells,weighted.PARTITION,weighted.PARTITION[1:]):
        source=cell['source']
        if source['source_family']!=family or source['candidate_N']!=N:
            raise ValueError('Same actual weighted source family and finite N required')
        if canonical_saved(source['exact_source_Rm_factor'])!=radius:
            raise ValueError('Same canonical original Rm radius factor required')
        if [canonical_saved(row) for row in source['exact_common_P0_axial5']]!=pressure:
            raise ValueError('Same canonical analytic P0 required in every source cell')
        if source['source_geometry']!=weighted.terminal.exact_geometry(left,right):
            raise ValueError('Same exact actual cell geometry required')
        for which,coordinate in (('left',left),('right',right)):
            saved=source['actual_endpoint_phases'][which]
            if saved!=phase_points[coordinate] or saved['exact_original_radius_phase']!=PHASE_RECIPE:
                raise ValueError('Same true positive-origin current radius phase required')
        if source.get('original_full_density_mean_not_zeroed') is not True:
            raise ValueError('Original nonlinear mean must be retained')
    if packet.get('ordinary_Z_Nminus2_averaging_installed') is not False:
        raise ValueError('Ordinary Z must retain its accepted direct all-u bound')


class OriginalRmFiniteNPrefix:
    mode='current_original_actual_Rm_incoming_plus_complete_weighted_local_source_to_Rh'
    def __init__(self,dps=500,require_checked=True):
        # This is the hydrate-only current owner chain, not contribution().
        self.upstream=prefix.OriginalMicroFiniteN(dps)
        self.c=self.upstream.c;self.family=self.upstream.family
        self.hashes=dict(self.upstream.hashes);self.saved={};self.cache={}
        all_u=weighted.previous.previous
        for module in (prefix,weighted,all_u):
            checked=json.loads((HERE/module.RECEIPT).read_bytes())
            if not checked.get('all_passed') or not checked.get(module.GATE) or checked['source_family']!=self.family:
                raise ValueError('Accepted same-family prefix/local source receipt required')
            report=json.loads(gzip.decompress((HERE/module.NAME).read_bytes()))
            if not report.get(module.GATE) or report['source_family']!=self.family or report['candidate_N']!=257:
                raise ValueError('Accepted same-family current N257 report required')
            for source in (checked,report):
                for name,digest in source['input_hashes'].items():fields.previous.bind(self.hashes,name,digest)
            for name in (module.NAME,module.RECEIPT,Path(module.__file__).name):
                fields.previous.bind(self.hashes,name,sha(name))
            self.saved[module.NAME]=report
        if RATES!=weighted.RATES:raise ValueError('Same five original own rates required')
        fields.previous.bind(self.hashes,Path(__file__).name,sha(Path(__file__).name))
        if require_checked:
            checked=json.loads((HERE/RECEIPT).read_bytes())
            if not checked.get('all_passed') or not checked.get(GATE) or checked['source_family']!=self.family:
                raise ValueError('Accepted actual Rm incoming consumer receipt required')
            for name,digest in checked['input_hashes'].items():fields.previous.bind(self.hashes,name,digest)
            fields.previous.bind(self.hashes,RECEIPT,sha(RECEIPT))

    def owner(self,label):return self.upstream.parameters.upstream.upstream.owner(label).op

    def contribution(self,label,N=257):
        N=prefix.phase.candidate_N(N)
        if N!=257:raise ValueError('Same accepted prefix and local candidate N257 required')
        if label not in ('0','.5'):raise ValueError('Only actual current source frames0,.5 supplied')
        if (label,N) in self.cache:return self.cache[(label,N)]
        op=self.owner(label);f,c=op.flow,op.c
        actual=self.saved[prefix.NAME]['frames'][label]
        local=self.saved[weighted.NAME]['frames'][label]
        baseline=self.saved[weighted.previous.previous.NAME]['frames'][label]
        with mp.workdps(c.dps+40):
            guard_prefix(self.family,label,N,op.P0,actual)
            switch.guard_saved_source(self.family,label,N,op.P0,baseline)
            if (baseline['complete_local_whole_patch_C0_Z_integral_bounds_installed'] is not True
                    or baseline['unknown_local_source_integral_count']!=0):
                raise ValueError('Complete all-u local source baseline required for this frame')
            parameter_family=self.upstream.parameters.upstream.upstream.upstream.repair
            mapper=weighted.previous.phase.RmRadiusPhase(op,self.family,parameter_family)
            if mapper.sc._mpi_!=self.upstream.macro.sc._mpi_:
                raise ValueError('Same actual positive micro phase origin required')
            phases={x:base.encoded(fields.serialized(mapper.point(x,N))) for x in weighted.PARTITION}
            guard_local(self.family,label,N,op.P0,op.Rm_factor,local,phases)
            restore=lambda row:switch.first.endpoint.restore_row(f,row)
            hydrate=lambda packet,key:{name:[restore(row) for row in values]
                for name,values in row_packet(packet,key).items()}
            incoming=hydrate(actual,'actual_current_Rm_incoming_correction_C0_Z')
            drive=hydrate(local,'actual_averaged_whole_patch_local_integral_C0_Z')
            # Original interval x in[1,e], y=log(x), width=1 exactly.
            width=c.mpf(1)
            if ep(fields.previous.read_interval(c,baseline['full_original_logarithmic_interval_width']))!=(1,1):
                raise ValueError('Original complete Rm-to-Rh logarithmic width1 required')
            outgoing=weighted.terminal.affine_transport(f,incoming,drive,width)
            memory={name:c.mpf(1) if not rate else c.exp(-c.mpf(rate.numerator)/rate.denominator)
                for name,rate in RATES.items()}
            result=dict(source_family=self.family,source_frame=label,candidate_N=N,
                exact_common_P0_axial5=op.P0,exact_source_Rm_factor=op.Rm_factor,
                exact_source_Rh_factor=op.Rm_factor*c.exp(1),exact_radial_partition=weighted.PARTITION,
                full_original_logarithmic_interval_width=width,original_incoming_Rm_to_Rh_decays=memory,
                actual_current_Rm_incoming_correction_C0_Z=incoming,
                accepted_complete_local_Rm_Rh_driver_C0_Z=drive,
                actual_current_Rh_correction_C0_Z=outgoing,
                bound_original_endpoint_phase_source=PHASE_RECIPE,
                bound_complete_local_source_cell_count=len(local['actual_source_weighted_cells']),
                real_finite_N_Rm_incoming_correction_is_supplied=True,
                genuine_current_finite_N_prefix_through_Rh_boundary_enclosures_supplied=True,
                original_pressure_incoming_decay_exactly_one=True,
                real_incoming_memory_preserved_not_reset=True,
                accepted_ordinary_Z_direct_local_source_bitwise_preserved=True,
                C0_weighted_local_source_caps_not_selected_values=True,
                no_old_incomplete_whole_atlas_transport_executed=True,
                older_Z0_whole_atlas_four_unknown_cells_not_promoted=True,
                source_bound_two_frame_enclosures_not_whole_Z_or_five_moment_closure=True,
                ordinary_Z_Nminus2_averaging_installed=False,
                no_ancestor_producer_checker_or_source_integration_rerun=True,
                **dict.fromkeys(fields.previous.OPEN,False))
        self.cache[(label,N)]=result;return result


def run():
    began=time.monotonic();owner=OriginalRmFiniteNPrefix(require_checked=False)
    frames={label:owner.contribution(label) for label in ('0','.5')}
    report=dict(mode=owner.mode,source_family=owner.family,**{GATE:True},candidate_N=257,
        frames=fields.serialized(frames),input_hashes=owner.hashes,
        execution_seconds=time.monotonic()-began,**dict.fromkeys(fields.previous.OPEN,False))
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(base.encoded(report),separators=(',',':'))+'\n').encode(),compresslevel=6,mtime=0))
    print('Actual source-bound finite-N prefix consumed through Rh',flush=True);return report


if __name__=='__main__':run()
