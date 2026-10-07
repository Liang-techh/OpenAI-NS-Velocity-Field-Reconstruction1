"""Fixed-N original whole-cell C1 Duhamel range transport.

This evaluates conservative source-function integral enclosures, not scalar
quadrature or uniform-in-N coefficients. Microscopic widths remain formal.
Unresolved source cells prevent downstream history and target admission.
"""
import json
from pathlib import Path
import time
import lei_ren_part1_paper_compliant_current_native_cutoff_density_oracle as density

transport=density.original.accepted.transport;current=transport.current
widths=current.transfer;packets=density.packets;native=density.native
HERE,PREFIX,sha=density.HERE,density.PREFIX,density.sha;ep=density.cutoff.ep
NAME=PREFIX+'current_native_cutoff_range_transport.json'
RECEIPT=PREFIX+'current_native_cutoff_range_transport_check.json'
GATE='current_original_cutoff_local_fixed_N_24_cell_C1_range_transport_executed'
RATES=transport.RATES


def records(rows):return {key:value.record() for key,value in rows.items()}


def fixed_N_target_rows(history,jets,A,AZ,logA,mu,logmu):
    """Normalize evaluated fixed-N histories; no all-N coefficient slots."""
    ratio=AZ.positive_divide(A,logA);den=A*A;values={};derivatives={}
    for out,key,degree in (('M','m',1),('I','h',1),('S','e',2),('Cp','p',2)):
        divisor=A if degree==1 else den
        value=history[key].positive_divide(divisor,degree*logA)
        values[out]=value;derivatives[out]=jets[key].positive_divide(divisor,degree*logA)-value*ratio*degree
    joint=history['k']-A*history['m'];jointZ=jets['k']-AZ*history['m']-A*jets['m']
    divided=current.repair.ROWS[1];divisor=den*mu
    values[divided]=joint.positive_divide(divisor,2*logA+logmu)
    derivatives[divided]=(jointZ-joint*ratio*2).positive_divide(divisor,2*logA+logmu)
    return dict(values=values,Z_derivatives=derivatives,joint_numerator=joint,joint_numerator_Z=jointZ)


class NativeCutoffRangeTransport:
    def __init__(self,role_owner):
        self.built=role_owner.build();self.oracle=density.NativeCutoffFactoredOracle(role_owner,self.built)
        self.target=self.oracle.target_owner;self.coordinates=self.target.coordinates
        self.ctx=self.oracle.ctx;self.family=self.oracle.source_family;self.service=self.oracle.service
        checked=json.loads((HERE/density.RECEIPT).read_bytes())
        if not checked.get('all_passed') or not checked.get(density.GATE) or checked['source_family']!=self.family:
            raise ValueError('Checked original cutoff density provider required')
        self.hashes={**self.oracle.hashes,**checked['input_hashes'],density.RECEIPT:sha(density.RECEIPT),
            Path(__file__).name:sha(Path(__file__).name)}
        self.coordinates.require_family(self.family);self.service.bind_hashes(self.hashes)
        if len(self.built['cells'])!=24 or [(q['label'],q['chart']) for q in self.built['cells']]!=[(q[0],q[1]) for q in current.ROUTE]:
            raise ValueError('Same exact original24-cell function route required')

    @native.inlet.source_precision
    def cell(self,label,chart,left,right,*,Z,N):
        N=density.density.density.candidate_integer(N)
        if N<160:raise ValueError('Same common original integer N>=160 required, including quiet cells')
        geom=self.target.geometry(label,chart,left,right);coords=self.coordinates
        factors={key:widths.true_width_kernel(coords,geom,rate) for key,rate in RATES.items()}
        for factor in factors.values():
            if ep(factor['mass'].coefficient)[0]<=0:raise ArithmeticError('Strict positive original true-width mass required')
            for value in (factor['mass'],factor['decay']):
                coords.scalar(0).coerce(value)
                if value.ctx is not self.ctx or value.ledger is not coords.ledger:
                    raise ValueError('Same common original mass/decay context and ledger required')
        frame=None
        if label=='initial_flat_collar':
            root_owner=self.target.q_owner.owner.owner
            eta=packets.interval(self.ctx,root_owner.scales['selected_positive_eta_log'])
            excess=packets.interval(self.ctx,self.target.owner.domain['left_collar_kappa_minus2_lower'])
            if ep(eta)[1]>ep(self.ctx.ln(excess)-self.ctx.ln(2))[0]:
                raise ValueError('Original entire initial collar exact-flat cutoff theorem required')
            values={key:coords.scalar(0) for key in RATES};jets=dict(values)
            source=dict(status='enclosed',original_exact_flat_initial_collar=True,
                original_flat_source_receipt=widths.RECEIPT,original_nonzero_E_V_P0_not_zeroed=True)
        else:
            frame=self.oracle.density_frame(chart=chart,Z=Z,coordinate=geom['coordinate'],N=N)
            source=frame.record
            if frame.values is None:
                return dict(record=dict(label=label,chart=chart,status='requires_original_whole_cell_source_refinement',
                    original_geometry=geom['record'],original_whole_cell_source=source),geometry=geom,factors=factors,values=None,Z_derivatives=None)
            if chart=='O3_power' and any(not value.zero for group in frame.values.values() for value in group.values()):
                raise ArithmeticError('Original O3 quiet cell source must be exactly flat')
            if chart=='O3_power':
                qsource=source['actual_original_spatial_source']['original_q_slow_jet_source']
                if not qsource['source_q_and_jet_ZERO_same_object']:
                    raise ArithmeticError('Original quiet q and ordinary ZERO jet must be the same exact source')
                branches=qsource['conditional_branches']
                if {branch['conditional_branch'] for branch in branches}!={'flat'}:
                    raise ArithmeticError('Original O3 quiet whole-cell cutoff branch must be proved flat')
            values={key:coords.rebase(frame.values['C0'][key],self.family)*factors[key]['mass'] for key in RATES}
            jets={key:coords.rebase(frame.values['Z'][key],self.family)*factors[key]['mass'] for key in RATES}
        record=dict(label=label,chart=chart,status='enclosed',candidate_N=N,Z_box=self.ctx.mpf(Z),
            original_geometry=geom['record'],original_whole_cell_source=source,
            true_log_radius_kernel_factors={key:dict(branch=item['branch'],mass=item['mass'].record(),decay=item['decay'].record())
                for key,item in factors.items()},C0_contributions=records(values),Z_contributions=records(jets),
            source_covers_entire_native_coordinate_and_Z_cell=True,
            cutoff_signed_u_and_phase_unions_before_one_positive_mass=True,
            original_true_width_collected_before_endpoint_subtraction=True,
            logarithmic_radius_measure_no_second_native_Jacobian=True,
            original_Z_independent_weights_and_endpoints_derivative_under_integral=True,
            original_P0_and_P0_Z_separate=True,no_scalar_quadrature_or_midpoint_field_used=True,
            **dict.fromkeys(packets.OPEN,False))
        return dict(record=record,geometry=geom,factors=factors,values=values,Z_derivatives=jets,frame=frame)

    @native.inlet.source_precision
    def route(self,*,Z=(-1,1),N=2048):
        N=density.density.density.candidate_integer(N)
        if N<160:raise ValueError('Same common original integer N>=160 required')
        coords=self.coordinates;c=self.ctx
        inlet=self.target.transfer.owner.inlet(Z,N=N)
        original_inlet=self.target.transfer.owner.native.left_inlet(Z)
        inlet_flags={key:original_inlet[key] for key in ('left_loop_q_A_B_exact_zero_by_existing_checked_collar',
            'original_P0_and_all_five_incoming_histories_retained')}
        if not all(inlet_flags.values()):raise ValueError('Original exact-flat inlet and full background/P0 memory theorem required')
        if any(not value.zero for group in ('initial_defect','initial_defect_Z') for value in inlet[group].values()):
            raise ValueError('Original exact zero correction inlet required')
        history={key:coords.scalar(0) for key in RATES};jets=dict(history);cells=[];unresolved=[]
        for label,chart,left,right in current.ROUTE:
            try:got=self.cell(label,chart,left,right,Z=Z,N=N)
            except ArithmeticError as error:
                got=dict(record=dict(label=label,chart=chart,status='requires_original_whole_cell_source_refinement',
                    arithmetic_obstruction=str(error),zero_or_midpoint_fallback_used=False),values=None,Z_derivatives=None)
            incoming,incomingZ=history,jets
            if got['values'] is None:
                unresolved.append(label);history=jets=None
            elif history is not None:
                history={key:got['factors'][key]['decay']*incoming[key]+got['values'][key] for key in RATES}
                jets={key:got['factors'][key]['decay']*incomingZ[key]+got['Z_derivatives'][key] for key in RATES}
            got['record'].update(incoming_C0=None if incoming is None else records(incoming),
                incoming_Z=None if incomingZ is None else records(incomingZ),
                outgoing_C0=None if history is None else records(history),outgoing_Z=None if jets is None else records(jets),
                incoming_memory_never_reset=True,quiet_source_does_not_zero_inherited_history=True,
                pressure_rate_zero_memory_retained=True)
            cells.append(got);print('Original whole-cell range transport:',label,got['record']['status'],flush=True)
        target=None;amplitude=None
        if history is not None:
            amplitude=self.oracle.amplitude_frame(Z=Z,N=N)
            logA=packets.interval(c,self.target.owner.reservation['positive_Ac_over_S_log_lower'])
            A=coords.rebase(amplitude.values['Rc_E_C0'],self.family).positive_intersection(logA)
            AZ=coords.rebase(amplitude.values['Rc_E_Z'],self.family)
            mu=packets.interval(c,self.target.owner.domain['right_collar_mu']);logmu=c.mpf(ep(c.ln(mu))[0])
            if mu._mpi_!=self.target.repair_mu._mpi_:raise ValueError('Same exact original target/repair mu required')
            # Fixed-N histories are not the -1/-2 all-N coefficient rows.
            target=dict(**fixed_N_target_rows(history,jets,A,AZ,logA,coords.scalar(mu),logmu),
                A=A,AZ=AZ,mu=mu,logA=logA,logmu=logmu)
        record=dict(source_family=self.family,candidate_N=N,Z_box=c.mpf(Z),original_exact_function_route_cells=24,
            original_whole_cell_queries=len(cells),enclosed_original_whole_cells=24-len(unresolved),
            unresolved_original_whole_cell_labels=unresolved,original_serial_cells=[q['record'] for q in cells],
            original_zero_inlet=inlet['record'],original_checked_inlet_theorem_flags=inlet_flags,original_common_coordinates=coords.record(),
            full24_original_C1_integral_range_transport_enclosed=history is not None,
            original_Rc_C0_history_ranges=None if history is None else records(history),
            original_Rc_Z_history_ranges=None if jets is None else records(jets),
            original_Rc_amplitude=None if amplitude is None else amplitude.record,
            original_Rc_target_C0_ranges=None if target is None else records(target['values']),
            original_Rc_target_Z_ranges=None if target is None else records(target['Z_derivatives']),
            original_joint_numerator=None if target is None else target['joint_numerator'].record(),
            original_joint_numerator_Z=None if target is None else target['joint_numerator_Z'].record(),
            fixed_N_range_integrals_not_scalar_values_or_uniform_coefficients=True,
            independent_hulls_do_not_prove_joint_target_cancellation=True,
            original_P0_P0_Z_and_incoming_memories_preserved=True,
            actual_five_controls_installed=False,actual_terminal_Z_function_closure_installed=False,current_whole_N_selected=False,
            **dict.fromkeys(packets.OPEN,False))
        return dict(record=record,cells=cells,history=history,Z_derivatives=jets,target=target,inlet=inlet)


@native.inlet.source_precision
def run(role_owner,*,return_live=False):
    began=time.monotonic();owner=NativeCutoffRangeTransport(role_owner);got=owner.route()
    result=dict(**got['record'],**{GATE:True},execution_seconds=time.monotonic()-began,input_hashes=owner.hashes,
        scope='Original fixed-N whole-cell C1 Duhamel range transport with true microscopic widths, full source/phase/cutoff covers and exact inlet. Unresolved cells block history/targets; scalar quadrature, quantitative joint cancellation, controls/global N/terminal closure and recursion remain separate.')
    (HERE/NAME).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    return (result,owner,got) if return_live else result
