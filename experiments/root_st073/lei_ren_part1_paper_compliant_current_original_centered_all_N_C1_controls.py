"""Fresh original24 centered first-Z coefficient bounds for every N>=160.

Exact N-dependent functions are inherited unchanged. Only branchwise
ordinary-Z covers are selected; no fixed-N integrated row is rescaled.
"""
import gzip
import json
from pathlib import Path
import time
import lei_ren_part1_paper_compliant_current_native_Rc_all_N_function_controls as current
import lei_ren_part1_paper_compliant_current_original_centered_first_bridge_C1_controls as centered

controls=current.controls;paired=current.paired;target_module=current.target_module
repair=current.repair;packets=current.packets;base=centered.base
HERE,PREFIX,sha=current.HERE,current.PREFIX,current.sha
ZERO,DZ,ORDERS,RATES,ep=current.ZERO,current.DZ,current.ORDERS,current.RATES,current.ep
encode=centered.encode
NAME=PREFIX+'current_original_centered_all_N_C1_controls.json.gz'
RECEIPT=PREFIX+'current_original_centered_all_N_C1_controls_check.json'
GATE='actual_original24_centered_C1_all_N_source_coefficients_and_conditional_repair_frequency_bound_installed'


class OriginalCenteredAllNC1Controls(current.NativeRcAllNFunctionControls):
    def __init__(self,functional_owner,range_owner):
        super().__init__(functional_owner,range_owner)
        self.previous=base.driver.accepted(self.hashes,current,self.family)
        base.driver.accepted(self.hashes,centered,self.family)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.service.bind_hashes(self.hashes)
        self.issued=None

    def primitive(self,conditional,eta,logamin,dstar):
        roots=conditional['source']['roots'];a=roots['a'][ZERO]
        if conditional['record']['source_family']!=self.family or a.ctx is not self.ctx:
            raise ValueError('Same live original source family/context required')
        for name in ('a','b','t0','kappa_minus2','E','p2'):
            for order in (ZERO,DZ):
                row=roots[name][order]
                if row.ctx is not self.ctx or row.scale.bases is not a.scale.bases or row.ledger is not a.ledger:
                    raise ValueError('Same original branch native basis/ledger required')
        # Retain the old complete C0 objects and same original q/q2 provenance.
        old=current.whole_period_primitives(roots,conditional['rows'],conditional['q2_rows'],logamin,dstar)
        # This repeats the original conditional branch construction, including
        # theta, correlated Delta_Z and direct q2_Z; aggregate q boxes are not
        # used as inputs to the centered branch theorem.
        alternative=centered.whole_centered_C1(roots,eta,logamin,dstar)
        values=dict(old['values']);decisions={}
        for key in ('A_Z','B_Z_over_Pstar'):
            values[key],decisions[key]=centered.accepted.smaller_cover(values[key],alternative['values'][key])
        return dict(values=values,old=old,alternative=alternative,decisions=decisions,
            record=dict(original_previous_primitive_C0_Z=base.records(old['values']),
                original_complete_centered_first_Z=alternative['record'],
                selected_original_primitive_C0_Z=base.records(values),primitive_Z_comparisons=decisions,
                original_C0_functions_and_full_period_covers_unchanged=True,
                same_source_family=self.family,first_Z_only=True,uniform_over_every_fractional_phase=True))

    @paired.native.inlet.source_precision
    def route(self,Z=(-1,1)):
        if tuple(Z)!=(-1,1):raise ValueError('Accepted complete original Z[-1,1] required')
        c=self.ctx;coords=self.target.coordinates;qcover=self.range_owner.oracle.owner.qcover
        root_owner=self.target.q_owner.owner.owner;signed=self.target.transfer.owner.signed_owner
        eta=packets.interval(c,root_owner.scales['selected_positive_eta_log'])
        dstar=packets.interval(c,root_owner.scales['logarithmic_selected_positive_lower_constants']['d_star'])
        if ep(eta)[1]>ep(-c.ln(2))[0]:raise ValueError('Original constant eta<=1/2 required')
        zeros=lambda:{key:{p:coords.scalar(0) for p in ORDERS} for key in RATES}
        history,jets=zeros(),zeros();cells=[];source_branches=0
        inlet=self.target.transfer.owner.inlet(Z,N=160)
        if any(not q.zero for q in [*inlet['initial_defect'].values(),*inlet['initial_defect_Z'].values()]):
            raise ValueError('Original exact zero correction inlet required')
        for label,chart,left,right in target_module.ROUTE:
            geometry=self.target.geometry(label,chart,left,right)
            factors={key:target_module.transfer.true_width_kernel(coords,geometry,rate) for key,rate in RATES.items()}
            branch_results=[];query=None
            if label=='initial_flat_collar':
                values,densityZ=zeros(),zeros()
                record=dict(original_zero_initial_collar_proof=target_module.transfer.RECEIPT,all_N_zero_support=True)
            else:
                query=qcover.query(chart,Z,geometry['coordinate'])
                positive=root_owner.decode(root_owner.inventory[chart]['actual_positive_denominator_theorem'])
                positivity='whole_actual_source_positive_not_inferred_from_saved_denominator_box' if chart.startswith('O3_') else 'source_function_positivity_not_inferred_from_saved_box'
                if positive.get(positivity) is not True:raise ValueError('Whole original chart source positivity required')
                for branch in query['branches']:
                    conditional=branch['query'];roots=conditional['source']['roots'];packet=conditional['source']['packet']
                    primitive=self.primitive(conditional,eta,positive['log_actual_a_positive_lower'],dstar)
                    def axial(k):
                        row=target_module.prior.signed.ordinary_axial_coefficient(packet.velocity['axial'][0],k)
                        return signed.leaf(target_module.prior.signed.expressions.RadiusPolynomial(packet.algebra,{0:row}),
                            roots['E'][ZERO].scale.bases,roots['E'][ZERO].ledger)
                    axial_values=[axial(0),axial(1)]
                    got=current.uniform_coefficients(roots['E'][ZERO],roots['E'][DZ],*axial_values,primitive['values'])
                    branch_results.append(dict(conditional=conditional,primitives=primitive,density=got,axial=axial_values,
                        record=dict(original_cutoff_branch=branch['record'],original_full_period_primitive_C1=primitive['record'],
                            original_all_N_density_coefficients=got['record'],original_positive_a_theorem=positive,
                            original_E_EZ_V_VZ=base.records(dict(E=roots['E'][ZERO],EZ=roots['E'][DZ],V=axial_values[0],VZ=axial_values[1])))))
                source_branches+=len(branch_results)
                union=paired.density.density.local.same_source_union
                values={key:{p:coords.rebase(union([b['density']['values'][key][p] for b in branch_results]),self.family)*factors[key]['mass']
                    for p in ORDERS} for key in RATES}
                densityZ={key:{p:coords.rebase(union([b['density']['Z_derivatives'][key][p] for b in branch_results]),self.family)*factors[key]['mass']
                    for p in ORDERS} for key in RATES}
                record=dict(original_source_query=query['record'],original_conditional_branch_results=[b['record'] for b in branch_results],
                    nonlinear_coefficient_ranges_before_overlapping_branch_hull=True,actual_N_or_selected_phase_not_used=True)
            incoming,incomingZ=history,jets
            history={key:{p:factors[key]['decay']*incoming[key][p]+values[key][p] for p in ORDERS} for key in RATES}
            jets={key:{p:factors[key]['decay']*incomingZ[key][p]+densityZ[key][p] for p in ORDERS} for key in RATES}
            record.update(label=label,chart=chart,original_geometry=geometry['record'],
                coefficient_contributions=target_module.records(values),coefficient_Z_contributions=target_module.records(densityZ),
                coefficient_outgoing=target_module.records(history),coefficient_Z_outgoing=target_module.records(jets),
                original_true_mass_applied_once=True,quiet_and_pressure_memories_not_reset=True)
            cells.append(dict(record=record,geometry=geometry,query=query,branches=branch_results,factors=factors,
                values=values,Z_derivatives=densityZ,incoming=incoming,incomingZ=incomingZ,outgoing=history,outgoingZ=jets))
            print('Original centered all-N coefficient C1 transport:',label,'branches',len(branch_results),flush=True)
        endpoint=2/self.target.transfer.geometry.binder.fixed['Tw']
        amplitude=self.target.q_owner.query('O3_power',Z,endpoint)['source']['roots']
        logA=packets.interval(c,self.target.owner.reservation['positive_Ac_over_S_log_lower'])
        A=coords.rebase(amplitude['E'][ZERO],self.family).positive_intersection(logA);AZ=coords.rebase(amplitude['E'][DZ],self.family)
        mu=packets.interval(c,self.target.owner.domain['right_collar_mu']);logmu=c.mpf(ep(c.ln(mu))[0])
        if mu._mpi_!=self.target.repair_mu._mpi_ or ep(mu)[0]<=0:raise ValueError('Same actual positive source/repair mu required')
        targets=target_module.target_rows(history,jets,A,AZ,logA,coords.scalar(mu),logmu)
        caps=target_module.target_caps(targets)
        W=repair.fresh_weights(c,mu,cells=512);matrix=repair.fresh_linear_inverse(c,mu,W)
        conditions=repair.contraction_log_conditions(c,caps,matrix,W,logmu,c.ln(target_module.MIN_N))
        self.hashes.update(self.target.service.hashes)
        self.issued=dict(cells=cells,history=history,Z_derivatives=jets,targets=targets,caps=caps,conditions=conditions,
            amplitude=A,amplitude_Z=AZ,mu=mu,logmu=logmu,logA=logA,weights=W,matrix=matrix,
            Z_box=c.mpf(Z),source_branch_count=source_branches,original_inlet=inlet)
        return self.issued

    def frequency_conditions(self,live):
        if live is not self.issued or live is None:raise ValueError('Issued fresh original all-N source route required')
        return current.combined_log_conditions(self,live)


def owner_factory(bridge):
    downstream=base.joint.downstream
    shell=downstream.NativeRcC1Histories(downstream.preceding.NativeO2C1Histories(downstream.preceding.make_middle_owner(bridge)))
    target=target_module.NativeRcParameterTargets(shell)
    role=controls.roles.RoleBoundNativeRcFunctions(controls.source.NativeRcFunctionTransport(target))
    return OriginalCenteredAllNC1Controls(controls.NativeRcFunctionalControls(role),paired.NativePairedC1Transport(role))


def comparison(owner,live):
    old=owner.previous['actual_uniform_N_scaled_target_C1_bounds']['transformed_N_scaled_target_C1_caps'];rows={}
    for key,new in live['caps']['transformed_N_scaled_target_C1_caps'].items():
        before=packets.interval(owner.ctx,old[key]['log_absolute_upper']);after=new['log_absolute_upper']
        rows[key]=dict(previous_log_upper=before,current_log_upper=after,strict_log_upper_reduction=ep(after)[1]<ep(before)[0])
    return rows


@paired.native.inlet.source_precision
def run(*,return_live=False,write=True):
    began=time.monotonic();bridge,construction=paired.native.inlet.native_bridge_owner()
    with paired.native.inlet.CheckedSourceRuntime() as runtime:
        owner=owner_factory(bridge);built=owner.build();live=owner.route();conditions=owner.frequency_conditions(live)
    old=owner.previous
    source=dict(source_family=owner.family,**{GATE:True},exact_function_graph_nodes=built['graph'].nodes,
        exact_function_source_graph_sha256=built['source_graph_sha256'],
        exact_coefficient_history_roots=current.encode_orders(built['coefficient_history']),
        exact_coefficient_target_roots=current.encode_orders(built['coefficient_target_orders']),
        exact_reconstructed_N_scaled_target_roots=controls.pair_roots(built['N_scaled_targets'].values(),built['N_scaled_targets'].keys()),
        exact_finite_Picard_function_roots=[controls.pair_roots(row,controls.CONTROLS) for row in built['finite_picard_sequence']],
        exact_all_N_source_reference_nodes=built['all_N_source_reference_nodes'],
        actual_all_N_continuous_cell_range_records=[cell['record'] for cell in live['cells']],
        actual_uniform_coefficient_history_C0_ranges=target_module.records(live['history']),
        actual_uniform_coefficient_history_Z_ranges=target_module.records(live['Z_derivatives']),
        actual_uniform_target_order_C0_ranges=target_module.records(live['targets']['values']),
        actual_uniform_target_order_Z_ranges=target_module.records(live['targets']['Z_derivatives']),
        actual_uniform_N_scaled_target_C1_bounds=live['caps'],actual_uniform_repair_C1_log_N_conditions=live['conditions'],
        original_source_and_repair_combined_log_conditions=conditions,
        previous_all_N_target_bound_comparison=comparison(owner,live),
        previous_all_N_log_conditions=old['original_source_and_repair_combined_log_conditions'],
        original_uniform_integer_N_lower=160,whole_Z_box=live['Z_box'],original_continuous_source_cells=24,
        original_conditional_coefficient_branch_count=live['source_branch_count'],
        primitive_strict_derivative_reductions=sum(r['strict_absolute_upper_reduction'] for cell in live['cells'] for b in cell['branches'] for r in b['primitives']['decisions'].values()),
        original_same_C0_and_exact_N_dependent_source_functions=True,
        exact_actual_phase='fractional_part(N*original_log_radius_minus_inlet)',
        exact_F_N='E*A*exprel(A/N)',exact_F_N_Z='E_Z*A*exprel(A/N)+E*A_Z*exp(A/N)',
        coefficient_functions_N_independent=False,actual_fixed_N_range_rescaled=False,
        original_P0_P0_Z_and_quiet_pressure_memories_preserved=True,
        original_bridge_construction=construction,original_source_runtime=runtime.record(),
        centered_theorem_receipt=dict(filename=centered.RECEIPT,sha256=sha(centered.RECEIPT)),
        previous_all_N_report=dict(filename=current.NAME,sha256=sha(current.NAME)),
        actual_five_controls_installed=False,certified_actual_fixed_point_tail_installed=False,
        actual_terminal_Z_function_closure_installed=False,functional_terminal_identity_solved=False,current_whole_N_selected=False,
        **dict.fromkeys(packets.OPEN,False),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
        scope='Fresh original24 all-phase/all-N centered ordinary-Z source bounds and conditional source/repair frequency threshold. Exact N-phase/coefficient functions, C0, geometry, P0 and pressure memory preserved. No selected global N, evaluated controls, functional terminal closure, heat/stress/recursion/pulses/full corrected NS admission.')
    if write:(HERE/NAME).write_bytes(gzip.compress(json.dumps(encode(source),indent=2).encode()+b'\n',compresslevel=9,mtime=0))
    print('Fresh original24 centered all-N C1 source bound; global N remains open',flush=True)
    return (source,owner,built,live) if return_live else source


if __name__=='__main__':run()
