"""Exact N-dependent coefficient functions and whole-source C1 repair contract.

The original phase frac(N*y) and exprel(A/N) stay in every exact source
coefficient. Uniform ranges come from new original full-period queries,
not by rescaling a fixed-N result. A repair-only log-N threshold does not
select a globally compatible N or install numerical fixed-point controls.
"""
import json
from pathlib import Path
import time
import sympy as sy
import lei_ren_part1_paper_compliant_current_native_Rc_functional_controls as controls

source=controls.source;target_module=source.current;paired=controls.paired
repair=controls.repair;packets=controls.packets
HERE,PREFIX,sha=controls.HERE,controls.PREFIX,controls.sha
ROWS,ORDERS,RATES=controls.ROWS,(-1,-2),source.RATES
ZERO,DZ=paired.ZERO,paired.DZ;ep=controls.ep
NAME=PREFIX+'current_native_Rc_all_N_function_controls.json'
RECEIPT=PREFIX+'current_native_Rc_all_N_function_controls_check.json'
GATE='current_original_N_dependent_coefficient_functions_and_uniform_C1_repair_contract_bound'


def normalize_history(g,history,A,mu):
    den=g.mul(A.value,A.value)
    ratio=g.quotient(A.Z,A.value,'same original positive Rc amplitude A')
    result={}
    for row,key,degree in (('M','m',1),('I','h',1),('S','e',2),('Cp','p',2)):
        divisor=A.value if degree==1 else den
        value=g.quotient(history[key].value,divisor,'same original positive Rc amplitude power')
        jet=g.sub(g.quotient(history[key].Z,divisor,'same original positive Rc amplitude power'),
            g.mul(g.constant(degree),ratio,value))
        result[row]=source.C1Function(value,jet)
    numerator=source.C1Function(g.sub(history['k'].value,g.mul(A.value,history['m'].value)),
        g.sub(history['k'].Z,g.add(g.mul(A.Z,history['m'].value),g.mul(A.value,history['m'].Z))))
    divisor=g.mul(mu,den)
    value=g.quotient(numerator.value,divisor,'same original positive mu and Rc amplitude squared')
    jet=g.quotient(g.sub(numerator.Z,g.mul(g.constant(2),ratio,numerator.value)),divisor,
        'same original positive mu and Rc amplitude squared')
    result[ROWS[1]]=source.C1Function(value,jet)
    return {key:result[key] for key in ROWS},numerator


def coefficient_pairs(g,E,V,A,B,N):
    """Exact C1 coefficients of the original full signed density increments."""
    arg=g.quotient(A.value,N,'one common original positive integer N')
    F=source.C1Function(g.mul(E.value,A.value,g.unary('exprel',arg)),
        g.add(g.mul(E.Z,A.value,g.unary('exprel',arg)),g.mul(E.value,A.Z,g.unary('exp',arg))))
    mul=g.c1mul;add=g.c1add;scale=g.c1scale
    negative=lambda pair:scale(g.constant(-1),pair)
    half=lambda pair:scale(g.constant('1/2'),pair)
    zero=source.C1Function(g.zero,g.zero)
    EF=mul(E,F);VF=mul(V,F);EB=mul(E,B);VB=mul(V,B)
    F2=mul(F,F);B2=mul(B,B)
    first=dict(m=B,h=F,k=add(VF,EB),e=add(scale(g.constant(2),VB),negative(EF)),p=EF)
    second=dict(m=zero,h=zero,k=mul(F,B),e=add(B2,negative(half(F2))),p=half(F2))
    return {-1:first,-2:second},F


def exact_all_N_graph(role_owner,*,iterations=3):
    built=role_owner.build();g=built['graph'];N=built['N'];parameters=built['parameters']
    maps,symbols,x=source.exact_radius_maps()
    history={p:{key:source.C1Function(g.zero,g.zero) for key in RATES} for p in ORDERS}
    cells=[];source_refs=[]
    for i,original in enumerate(built['cells']):
        chart=original['chart'];quiet=original['source_flat_exact_zero'];view=role_owner.owner.views[chart]
        t=g.symbol('coordinate_'+str(i));offset=source.expression(g,maps[chart],parameters,t)
        phase=g.unary('fractional_part',g.mul(N,offset))
        jacobian=source.expression(g,sy.diff(maps[chart],x),parameters,t)
        def ref(node,namespace,role):
            value=g.node('original_function_graph',chart=chart,graph_file=source.sources.VIEWS,
                graph_sha256=built['source_graph_sha256'],source_node=node,coordinate=t.node,
                phase=phase.node,shared_N=N.node,Z_variable='Z',source_graph_namespace=namespace,function_role=role,
                source_coefficients_are_function_recipes_not_cover_values=True)
            source_refs.append(value.node);return value
        coefficients={p:{key:source.C1Function(g.zero,g.zero) for key in RATES} for p in ORDERS};F=None;bound_inputs=None
        if not quiet:
            roots=view['roots'];original_roots=view['original_signed_input_graph']['jet_expression_dag']['roots']
            E=source.C1Function(*[ref(original_roots['E'][order],
                'original_signed_input_graph.jet_expression_dag','all_N_original_E_'+suffix)
                for order,suffix in (('y0_Z0','C0'),('y0_Z1','Z'))])
            V_nodes=[]
            for order in ('y0_Z0','y0_Z1'):
                matches=[j for j,row in enumerate(view['function_graph_nodes']) if row['operation']=='source_derivative' and row['name']=='V_'+order]
                if len(matches)!=1:raise ValueError('Unique original axial source derivative required')
                V_nodes.append(matches[0])
            V=source.C1Function(*[ref(j,'function_graph_nodes','all_N_original_V_'+suffix) for j,suffix in zip(V_nodes,('C0','Z'))])
            A=source.C1Function(ref(roots['A'],'function_graph_nodes','all_N_periodic_A_C0'),
                ref(roots['A_Z_slow'],'function_graph_nodes','all_N_periodic_A_Z'))
            B=source.C1Function(ref(roots['B_over_Pstar'],'function_graph_nodes','all_N_periodic_B_C0'),
                ref(roots['B_Z_slow'],'function_graph_nodes','all_N_periodic_B_Z'))
            coefficients,F=coefficient_pairs(g,E,V,A,B,N);bound_inputs=dict(E=E,V=V,A=A,B=B)
        incoming=history;contributions={p:{} for p in ORDERS};outgoing={p:{} for p in ORDERS}
        for p in ORDERS:
            for key,rate in RATES.items():
                decay=g.unary('exp',g.neg(g.mul(g.constant(rate),source.FunctionRef(g,original['width']))))
                kernel=g.unary('exp',g.neg(g.mul(g.constant(rate),g.sub(source.FunctionRef(g,original['right_radius_offset']),offset))))
                pair=coefficients[p][key];integrals=[]
                for handle in (pair.value,pair.Z):
                    if handle==g.zero:integrals.append(g.zero);continue
                    integrals.append(controls.integral(g,g.mul(kernel,handle,jacobian),'coordinate_'+str(i),
                        source.FunctionRef(g,original['lower']),source.FunctionRef(g,original['upper']),
                        measure='native coordinate; original dy/dcoordinate applied exactly once',
                        extracted_N_power=p,coefficient_still_depends_on_N_and_actual_phase=True))
                contribution=source.C1Function(*integrals);contributions[p][key]=contribution
                outgoing[p][key]=g.c1add(g.c1scale(decay,incoming[p][key]),contribution)
        history=outgoing
        cells.append(dict(original_cell=original,coefficient_density_pairs=coefficients,
            coefficient_contributions=contributions,incoming=incoming,outgoing=outgoing,F_N=F,source_inputs=bound_inputs))
    orders={};joint={}
    for p in ORDERS:orders[p],joint[p]=normalize_history(g,history[p],built['amplitude'],parameters['mu'])
    invN=g.quotient(g.one,N,'one common original positive integer N')
    scaled={key:g.c1add(orders[-1][key],g.c1scale(invN,orders[-2][key])) for key in ROWS}
    total_targets={key:g.c1scale(invN,pair) for key,pair in scaled.items()}
    total_history={key:g.c1add(g.c1scale(invN,history[-1][key]),
        g.c1scale(g.mul(invN,invN),history[-2][key])) for key in RATES}
    combined=dict(**built)
    combined.update(original_full_density_history=built['history'],original_full_density_targets=built['targets'],
        history=total_history,targets=total_targets,N_scaled_targets=scaled,
        all_N_cells=cells,coefficient_history=history,coefficient_target_orders=orders,
        coefficient_joint_numerators=joint,all_N_source_reference_nodes=source_refs)
    return controls.exact_control_graph(combined,iterations=iterations)


def uniform_coefficients(E,EZ,V,VZ,primitive):
    """Uniform for all N>=160 and all phases; coefficients themselves vary with N."""
    c=E.ctx;A,AZ,B,BZ=[primitive[key] for key in ('A','A_Z','B_over_Pstar','B_Z_over_Pstar')]
    zero=E.scalar(0)
    argument=c.mpf(0) if A.zero else c.mpf(('-1.25','1.25'))/target_module.MIN_N
    factor=c.exp(argument);F=E*A*factor;FZ=EZ*A*factor+E*AZ*factor
    values=dict(m={-1:B,-2:zero},h={-1:F,-2:zero},k={-1:V*F+E*B,-2:F*B},
        e={-1:2*V*B-E*F,-2:B*B-F*F*c.mpf('.5')},p={-1:E*F,-2:F*F*c.mpf('.5')})
    jets=dict(m={-1:BZ,-2:zero},h={-1:FZ,-2:zero},
        k={-1:VZ*F+V*FZ+EZ*B+E*BZ,-2:FZ*B+F*BZ},
        e={-1:2*VZ*B+2*V*BZ-EZ*F-E*FZ,-2:2*B*BZ-F*FZ},p={-1:EZ*F+E*FZ,-2:F*FZ})
    return dict(values=values,Z_derivatives=jets,record=dict(N_lower=160,fractional_phase_cover=[0,1],
        original_exprel_and_exp_argument_uniform_cover=argument,coefficients_still_depend_on_N=True,
        analytic_original_A_cap_1p25_bounds_exponential_factors_only=True,actual_fixed_N_result_rescaled=False))


def whole_period_primitives(roots,qrows,q2rows,log_a,dstar):
    """Reuse original paired derivative proofs and derive C0 caps before density."""
    got=paired.paired_C1_support(roots,qrows,q2rows,log_a,dstar)
    c=roots['a'][ZERO].ctx;scalar=roots['a'][ZERO].scalar
    if qrows[ZERO].zero and qrows[DZ].zero:
        A_cap=B_cap=scalar(0)
    else:
        restrict=paired.accepted.accepted.accepted.signed_C0_support_range
        a,unused=restrict(roots['a'][ZERO],scalar(c.mpf((0,'2.5'))))
        b,unused=restrict(roots['b'][ZERO],scalar(c.mpf(('-1.5','1.5'))))
        minimum=paired.accepted.serial.minimum_upper;absolute=paired.accepted.serial.absolute_upper
        T=minimum((-b).positive_divide(a,log_a),roots['t0'][ZERO])
        Q,R,AU=[absolute(v) for v in (qrows[ZERO],q2rows[ZERO],a)]
        A_cap=minimum(scalar(c.mpf('1.25')),AU*(T*Q*4+R*101))
        B_cap=minimum(scalar(c.mpf('1.5')),T*A_cap+AU*Q*2)
    values=dict(A=paired.accepted.serial.symmetric_bound(A_cap),
        B_over_Pstar=roots['E'][ZERO]*paired.accepted.serial.symmetric_bound(B_cap),**got['values'])
    return dict(values=values,record=dict(**got['record'],original_C0_caps=dict(A=A_cap.record(),B_over_E=B_cap.record()),
        C0_and_Z_covers_uniform_over_all_fractional_phases=True))


class NativeRcAllNFunctionControls:
    def __init__(self,functional_owner,range_owner):
        if type(functional_owner) is not controls.NativeRcFunctionalControls or type(range_owner) is not paired.NativePairedC1Transport:
            raise TypeError('Existing accepted original functional and paired range owners required')
        if functional_owner.target is not range_owner.target:raise ValueError('Same actual source target owner required')
        self.functional_owner,self.range_owner=functional_owner,range_owner
        self.target=functional_owner.target;self.role_owner=functional_owner.role_owner
        self.ctx,self.family,self.service=self.target.ctx,self.target.family,self.target.service
        checked=json.loads((HERE/controls.RECEIPT).read_bytes())
        if not checked.get('all_passed') or not checked.get(controls.GATE) or checked['source_family']!=self.family:
            raise ValueError('Checked same-original functional control map required')
        self.hashes={**functional_owner.hashes,**range_owner.hashes,**checked['input_hashes'],
            controls.NAME:sha(controls.NAME),controls.RECEIPT:sha(controls.RECEIPT),Path(__file__).name:sha(Path(__file__).name)}
        self.service.bind_hashes(self.hashes)

    def build(self,*,iterations=3):return exact_all_N_graph(self.role_owner,iterations=iterations)

    @paired.native.inlet.source_precision
    def route(self,Z=(-1,1)):
        c=self.ctx;coords=self.target.coordinates;qcover=self.range_owner.oracle.owner.qcover
        root_owner=self.target.q_owner.owner.owner;signed=self.target.transfer.owner.signed_owner
        dstar=packets.interval(c,root_owner.scales['logarithmic_selected_positive_lower_constants']['d_star'])
        zeros=lambda:{key:{p:coords.scalar(0) for p in ORDERS} for key in RATES}
        history,jets=zeros(),zeros();cells=[];source_branches=0
        # Proven original zero inlet for every N; the geometry adapter's
        # N argument here is not a repair selection or a source evaluation.
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
                for branch in query['branches']:
                    conditional=branch['query'];roots=conditional['source']['roots'];packet=conditional['source']['packet']
                    primitive=whole_period_primitives(roots,conditional['rows'],conditional['q2_rows'],
                        positive['log_actual_a_positive_lower'],dstar)
                    def axial(k):
                        row=target_module.prior.signed.ordinary_axial_coefficient(packet.velocity['axial'][0],k)
                        return signed.leaf(target_module.prior.signed.expressions.RadiusPolynomial(packet.algebra,{0:row}),
                            roots['E'][ZERO].scale.bases,roots['E'][ZERO].ledger)
                    got=uniform_coefficients(roots['E'][ZERO],roots['E'][DZ],axial(0),axial(1),primitive['values'])
                    branch_results.append(dict(conditional=conditional,primitives=primitive,density=got,
                        record=dict(original_cutoff_branch=branch['record'],original_full_period_primitive_C1=primitive['record'],
                            original_all_N_density_coefficients=got['record'])))
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
            print('Original all-N coefficient C1 transport:',label,'branches',len(branch_results),flush=True)
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
        return dict(cells=cells,history=history,Z_derivatives=jets,targets=targets,caps=caps,conditions=conditions,
            amplitude=A,amplitude_Z=AZ,mu=mu,logmu=logmu,logA=logA,weights=W,matrix=matrix,
            Z_box=c.mpf(Z),source_branch_count=source_branches,original_inlet=inlet)


def encode_orders(rows):
    return {str(p):controls.pair_roots(pairs.values(),pairs.keys()) for p,pairs in rows.items()}


def combined_log_conditions(owner,live):
    """Include original generic source requirements without claiming global N."""
    c=owner.ctx;views=owner.role_owner.owner.views
    source_rows={chart:packets.interval(c,view['actual_source_frequency_majorants']['required_positive_log_N_lower'])
        for chart,view in views.items()}
    source_lower=c.mpf(max(ep(value)[1] for value in source_rows.values()))
    repair_lower=live['conditions']['repair_sufficient_common_log_N_lower']
    combined=c.mpf(max(ep(c.ln(160))[1],ep(source_lower)[1],ep(repair_lower)[1]))
    return dict(original_generic_source_log_N_requirements=source_rows,
        original_generic_source_log_N_lower=source_lower,uniform_repair_and_positivity_log_N_lower=repair_lower,
        source_and_repair_sufficient_log_N_lower=combined,
        finite_integer_existence_recipe='Any integer N>=ceil(exp(the recorded directed log-N upper threshold))',
        source_and_repair_requirements_all_lower_bounds=True,integer_or_exponential_not_materialized=True,
        higher_y_N_power0_terms_retained_in_original_sidecars=True,
        remaining_global_cone_higher_jet_join_N_conditions_not_certified=True,current_whole_N_selected=False)


@paired.native.inlet.source_precision
def record(owner,built,live,*,execution_seconds=None,warm_original_route_reused=False):
    old=json.loads((HERE/target_module.NAME).read_bytes());comparison={}
    for key,row in live['caps']['transformed_N_scaled_target_C1_caps'].items():
        before=packets.interval(owner.ctx,old['actual_original_Rc_parameter_target_records']['whole_Z']
            ['actual_uniform_N_scaled_repair_C1_caps']['transformed_N_scaled_target_C1_caps'][key]['log_absolute_upper'])
        after=row['log_absolute_upper'];comparison[key]=dict(previous_log_upper=before,current_log_upper=after,
            strict_log_upper_reduction=ep(after)[1]<ep(before)[0])
    result=dict(source_family=owner.family,**{GATE:True},exact_function_graph_nodes=built['graph'].nodes,
        exact_original_full_density_history_roots=controls.pair_roots(built['original_full_density_history'].values(),built['original_full_density_history'].keys()),
        exact_coefficient_history_roots=encode_orders(built['coefficient_history']),
        exact_coefficient_target_roots=encode_orders(built['coefficient_target_orders']),
        exact_reconstructed_history_roots=controls.pair_roots(built['history'].values(),built['history'].keys()),
        exact_reconstructed_N_scaled_target_roots=controls.pair_roots(built['N_scaled_targets'].values(),built['N_scaled_targets'].keys()),
        exact_finite_Picard_function_roots=[controls.pair_roots(row,controls.CONTROLS) for row in built['finite_picard_sequence']],
        exact_all_N_source_reference_nodes=built['all_N_source_reference_nodes'],
        coefficient_functions_N_independent=False,exact_actual_phase='fractional_part(N*original_log_radius_minus_inlet)',
        exact_F_N='E*mathcal_A*exprel(mathcal_A/N)',exact_F_N_Z='E_Z*mathcal_A*exprel(mathcal_A/N)+E*mathcal_A_Z*exp(mathcal_A/N)',
        original_continuous_source_cells=24,whole_Z_box=live['Z_box'],uniform_integer_N_lower=160,
        original_conditional_coefficient_branch_count=live['source_branch_count'],
        actual_all_N_continuous_cell_range_records=[row['record'] for row in live['cells']],
        actual_uniform_coefficient_history_C0_ranges=target_module.records(live['history']),
        actual_uniform_coefficient_history_Z_ranges=target_module.records(live['Z_derivatives']),
        actual_uniform_target_order_C0_ranges=target_module.records(live['targets']['values']),
        actual_uniform_target_order_Z_ranges=target_module.records(live['targets']['Z_derivatives']),
        actual_uniform_N_scaled_target_C1_bounds=live['caps'],actual_uniform_repair_C1_log_N_conditions=live['conditions'],
        original_source_and_repair_combined_log_conditions=combined_log_conditions(owner,live),
        previous_all_N_target_bound_comparison=comparison,
        exact_original_Banach_family_contract=dict(source_graph_bound=True,domain='Z in[-1,1]',
            frequency='any integer N>=160 satisfying the recorded repair/source/positivity log bounds',
            target='exact reconstructed N*r(N,Z), not enclosure endpoints',
            existence_and_uniqueness_conditional_on_compatible_frequency=True,
            tail='L^n/(1-L)*norm_C1(h1-h0) after the same-source ball and contraction conditions'),
        analytic_source_expm1_argument_safety=dict(N_lower=160,primitive_absolute_upper='5/4',
            argument_absolute_upper='1/128',current_C0_Z_source_exp_and_exprel_factors_uniform=True,
            raw_higher_y_majorant_or_other_global_N_conditions_not_admitted=True),
        original_P0_P0_Z_and_all_histories_preserved=True,actual_fixed_N_range_rescaled=False,
        source_caps_or_midpoints_chosen_as_functions=False,actual_source_ancestor_constructors_called=False,
        actual_five_controls_installed=False,certified_actual_fixed_point_tail_installed=False,
        actual_terminal_Z_function_closure_installed=False,current_whole_N_selected=False,
        **dict.fromkeys(packets.OPEN,False),input_hashes=owner.hashes,execution_seconds=execution_seconds,
        warm_original_route_reused=warm_original_route_reused,
        scope='Exact N-dependent 1/N and 1/N^2 coefficient functions with actual spatial phase, all24 original continuous/full-Z uniform C1 source/history/target ranges, and source-bound repair-only log-N conditions. No fixed-N rescaling, chosen global N, evaluated fixed point/terminal closure or recursion.')
    (HERE/NAME).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    print('Exact all-N coefficient functions and uniform original C1 repair contract bound; global N remains open',flush=True)
    return result


@paired.native.inlet.source_precision
def run(functional_owner,range_owner,*,return_live=False):
    began=time.monotonic();owner=NativeRcAllNFunctionControls(functional_owner,range_owner)
    built=owner.build();live=owner.route()
    result=record(owner,built,live,execution_seconds=time.monotonic()-began)
    return (result,owner,built,live) if return_live else result
