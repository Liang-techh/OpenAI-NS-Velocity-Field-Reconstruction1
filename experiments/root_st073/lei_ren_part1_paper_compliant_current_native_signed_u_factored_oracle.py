"""Connect the accepted signed-u cover to original role-bound graph sources.

Declared interior queries exercise all seventeen chart providers. Each
unresolved source stays explicit; a point query is never whole-chart admission.
"""
import json
from pathlib import Path
import time
import lei_ren_part1_paper_compliant_current_native_signed_u_phase_cover as cover
import lei_ren_part1_paper_compliant_current_native_Rc_factored_source_oracle as accepted

HERE,PREFIX,sha=cover.HERE,cover.PREFIX,cover.sha
NAME=PREFIX+'current_native_signed_u_factored_oracle.json'
RECEIPT=PREFIX+'current_native_signed_u_factored_oracle_check.json'
GATE='current_original_signed_u_role_bound_range_oracle_all17_declared_queries_executed'
POINTS=dict(bridge_first='.1337',bridge_second='1.831',bridge_macro='.537',
    switch_first='.537',switch_second='1.337',switch_power='.537',reshape='.537',
    inner_reference='.537',axial_restore='.537',restore_buffer='-6.337',
    actual_patch='1.337',Rh_reference='-2.337',O2_slope='.537',O2_axial='.1337',
    O2_buffer='5.337',O3_slope_mu='.537',O3_power={'original_power_offset':'.537'})


class SignedUCoverFactoredOracle(accepted.NativeFactoredSourceOracle):
    def __init__(self,role_owner,built=None):
        super().__init__(role_owner,built)
        receipt=json.loads((HERE/cover.RECEIPT).read_bytes())
        if not receipt.get('all_passed') or not receipt.get(cover.GATE) or receipt['source_family']!=self.source_family:
            raise ValueError('Accepted same-original-family signed-u range cover required')
        self.original_density_owner=self.owner
        self.owner=cover.NativeSignedUDensityCover(self.original_density_owner)
        self.hashes={**self.hashes,**receipt['input_hashes'],cover.RECEIPT:sha(cover.RECEIPT),
            Path(__file__).name:sha(Path(__file__).name)}
        self.service.bind_hashes(self.hashes)

    @cover.native.inlet.source_precision
    def density_frame(self,*,chart,Z,coordinate,N):
        frame=super().density_frame(chart=chart,Z=Z,coordinate=coordinate,N=N)
        frame.record['phase_solver_backend']='accepted original first-jet body with overlapping signed-u range cover'
        frame.record['original_signed_u_cover_receipt']=cover.RECEIPT
        frame.record['original_role_namespace_and_graph_hash_dispatch_unchanged']=True
        frame.record['original_nonlinear_density_precedes_overlapping_branch_hull']=True
        frame.record['selected_declared_source_query_not_whole_chart_admission']=True
        return frame


@cover.native.inlet.source_precision
def run(role_owner):
    began=time.monotonic();built=role_owner.build();oracle=SignedUCoverFactoredOracle(role_owner,built)
    if set(POINTS)!=set(cover.native.DOMAINS):raise ValueError('All seventeen original chart names required')
    rows={(row.get('chart'),row.get('function_role')):row for row in built['graph'].nodes if row['operation']=='original_function_graph'}
    records={};enclosed=0;dispatches=0
    for chart,coordinate in POINTS.items():
        try:frame=oracle.density_frame(chart=chart,Z=(-1,1),coordinate=coordinate,N=2048)
        except ArithmeticError as error:
            records[chart]=dict(status='requires_original_source_range_refinement',chart=chart,coordinate=coordinate,
                Z_box=[-1,1],candidate_N=2048,arithmetic_obstruction=str(error),zero_or_midpoint_fallback_used=False)
        else:
            records[chart]=frame.record
            if frame.values is not None:
                enclosed+=1
                for key in cover.density.RATES:
                    for order in ('C0','Z'):
                        row=rows.get((chart,'density_'+key+'_'+order))
                        if row is not None:oracle.dispatch_function_range(row,frame);dispatches+=1
        print('Original signed-u chart provider:',chart,records[chart]['status'],flush=True)
    c=oracle.ctx;box=c.mpf((cover.ep(c.mpf('.12'))[0],cover.ep(c.mpf('.15'))[1]))
    broad=oracle.density_frame(chart='O2_slope',Z=(-1,1),coordinate=box,N=2048)
    if broad.values is None:raise ArithmeticError('Previously accepted broad signed-crossing O2 must resolve')
    for key in cover.density.RATES:
        for order in ('C0','Z'):oracle.dispatch_function_range(rows['O2_slope','density_'+key+'_'+order],broad);dispatches+=1
    amplitude=oracle.amplitude_frame(Z=(-1,1),N=2048)
    for role in ('Rc_E_C0','Rc_E_Z'):oracle.dispatch_function_range(rows['O3_power',role],amplitude);dispatches+=1
    result=dict(source_family=oracle.source_family,**{GATE:True},declared_original_chart_query_count=17,
        enclosed_original_chart_queries=enclosed,unresolved_original_chart_queries=17-enclosed,
        original_chart_query_records=records,original_broad_signed_u_O2_frame=broad.record,
        separate_original_Rc_amplitude_frame=amplitude.record,actual_original_role_dispatches=dispatches,
        declared_exact_interior_coordinates=POINTS,full_Z_interval=[-1,1],candidate_N=2048,
        all17_original_providers_exercised=True,all17_whole_chart_or_continuous_route_ranges_admitted=False,
        full_factored_function_graph_evaluator_installed=False,actual_original_full_route_numerical_integrals_evaluated=False,
        actual_five_controls_installed=False,actual_terminal_Z_function_closure_installed=False,current_whole_N_selected=False,
        **dict.fromkeys(cover.packets.OPEN,False),execution_seconds=time.monotonic()-began,input_hashes=oracle.hashes,
        scope='Original role-bound factored range dispatch with accepted signed-u cover at seventeen declared full-Z interior queries, broad O2 cell and original Rc E/E_Z. Unresolved provider records remain explicit; sampled provider coverage is not whole-chart/full-route integration or control/closure/recursion admission.')
    (HERE/NAME).write_text(json.dumps(cover.packets.encode(result),indent=2)+'\n',encoding='utf8')
    return result
