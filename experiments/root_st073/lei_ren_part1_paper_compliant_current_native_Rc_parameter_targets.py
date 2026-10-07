"""Source-backed all-N C0/Z Rc correction orders and original repair targets.

The coefficients below are uniform covers of N-dependent source functions,
not selected values or exact N-independent polynomial coefficients. Every
original cell uses new periodic bounds; no N1024 output is rescaled.
"""
import json
from pathlib import Path
import time
import sympy as sy
import lei_ren_part1_paper_compliant_current_native_Rc_C1_histories as preceding
import lei_ren_part1_paper_compliant_current_generic_moment_repair_operator as repair

serial=preceding.serial;transfer=preceding.transfer;history=preceding.history
prior=preceding.prior;native=preceding.native;packets=preceding.packets;density=preceding.density
HERE,PREFIX,sha=preceding.HERE,preceding.PREFIX,preceding.sha;ep=preceding.ep;RATES=preceding.RATES
ZERO,DZ=serial.ZERO,serial.DZ;ORDERS=(-1,-2);MIN_N=160
NAME=PREFIX+'current_native_Rc_parameter_targets.json';RECEIPT=PREFIX+'current_native_Rc_parameter_targets_check.json'
GATE='current_actual_original_Rc_all_N_C1_correction_orders_and_repair_targets_executed'
MID=preceding.preceding.preceding;SWITCH=MID.preceding
ROUTE=(('initial_flat_collar','bridge_first',None,None),('active_first_bridge','bridge_first',None,1),
    ('second_bridge','bridge_second',1,2))+tuple((chart,chart,left,right) for chart,left,right in SWITCH.ROUTE)+\
    tuple((chart,chart,left,right) for chart,left,right in MID.ROUTE)+preceding.preceding.ROUTE+preceding.ROUTE


def exact_density_order_theorem():
    N,E,V,F,B=sy.symbols('N E V F_N B',nonzero=True)
    EZ,VZ,FZ,BZ=sy.symbols('E_Z V_Z F_N_Z B_Z')
    dE,dV=F/N,B/N
    exact=dict(m=dV,h=dE,k=V*dE+E*dV+dE*dV,
        e=2*V*dV+dV*dV-E*dE-dE*dE/2,p=E*dE+dE*dE/2)
    coefficients=dict(m={-1:B,-2:0},h={-1:F,-2:0},
        k={-1:V*F+E*B,-2:F*B},e={-1:2*V*B-E*F,-2:B*B-F*F/2},
        p={-1:E*F,-2:F*F/2})
    for key,value in exact.items():
        if sy.expand(value-sum(coefficient*N**power for power,coefficient in coefficients[key].items()))!=0:
            raise ArithmeticError('Original all-N signed density order identity failed: '+key)
        derivative=sum(sy.diff(value,variable)*jet for variable,jet in ((E,EZ),(V,VZ),(F,FZ),(B,BZ)))
        reconstructed=sum(sum(sy.diff(coefficient,variable)*jet for variable,jet in
            ((E,EZ),(V,VZ),(F,FZ),(B,BZ)))*N**power for power,coefficient in coefficients[key].items())
        if sy.expand(derivative-reconstructed)!=0:raise ArithmeticError('Original first-Z order identity failed: '+key)
    return dict(passed=True,exact_five_signed_density_and_first_Z_order_identities=True,
        original_phase_primitive_name='mathcal_A (distinct from repair amplitude A_c/S)',
        g_N='integral_0^1 exp(t*mathcal_A/N)dt',F_N='E*mathcal_A*g_N',
        F_N_Z='E_Z*mathcal_A*g_N+E*mathcal_A_Z*exp(mathcal_A/N)',
        F_N_derivative_identity='d_mathcal_A(mathcal_A*g_N)=exp(mathcal_A/N)',
        F_N_and_coefficient_functions_still_depend_on_N=True,
        coefficient_enclosures_uniform_for_every_integer_N_ge=MIN_N,
        original_global_phase='fractional_part(N*log(R/r_minus)); same original inverse; y_Z=0',
        all_N_cover_not_rescaled_from_fixed_N_results=True,higher_jets_not_claimed=True)


def uniform_density_orders(E,EZ,V,VZ,primitive):
    c=E.ctx;A,AZ,B,BZ=[primitive[key] for key in ('A','A_Z','B_over_Pstar','B_Z_over_Pstar')]
    argument=density.phase.bounded_value(A*(c.mpf(1)/MIN_N))
    lo,hi=ep(argument);argument=c.mpf((min(0,lo),max(0,hi)))
    # The analytic source theorem |mathcal_A|<=159 is stronger than the
    # outward log/exp arithmetic hull at its last rounding ulp. Intersect
    # bound coordinates with that theorem, without clipping field values.
    universal=c.mpf((-159,159))/MIN_N
    argument=c.mpf((max(ep(argument)[0],ep(universal)[0]),min(ep(argument)[1],ep(universal)[1])))
    expcover=c.exp(argument);gcover=expcover
    F=E*A*gcover;FZ=EZ*A*gcover+E*AZ*expcover;zero=E.scalar(0)
    values=dict(m={-1:B,-2:zero},h={-1:F,-2:zero},k={-1:V*F+E*B,-2:F*B},
        e={-1:2*V*B-E*F,-2:B*B-F*F*c.mpf('.5')},p={-1:E*F,-2:F*F*c.mpf('.5')})
    jets=dict(m={-1:BZ,-2:zero},h={-1:FZ,-2:zero},
        k={-1:VZ*F+V*FZ+EZ*B+E*BZ,-2:FZ*B+F*BZ},
        e={-1:2*VZ*B+2*V*BZ-EZ*F-E*FZ,-2:2*B*BZ-F*FZ},
        p={-1:EZ*F+E*FZ,-2:F*FZ})
    return dict(values=values,Z_derivatives=jets,record=dict(primitive_over_N_union=argument,
        g_N_uniform_positive_cover=gcover,exp_primitive_over_N_uniform_cover=expcover,
        original_universal_primitive_cap159_intersected_as_function_theorem=True,
        uniform_N_lower=MIN_N,coefficients_remain_N_dependent_functions=True,
        every_original_axial_cross_and_quadratic_term_retained=True))


def records(rows):return {key:{str(power):value.record() for power,value in row.items()} for key,row in rows.items()}


def evaluate(rows,N):
    N=density.spatial.density.candidate_integer(N)
    if N<MIN_N:raise ValueError('All-N source lift requires integer N>=160')
    return {key:sum((value*(value.ctx.mpf(1)/N)**(-power) for power,value in row.items()),next(iter(row.values())).scalar(0))
        for key,row in rows.items()}


def magnitude(value):
    upper=serial.absolute_upper(value)
    return repair.LogUpper(value.ctx,None if upper.zero else upper.scale.evaluate())


def target_rows(values,jets,A,AZ,logA,mu,logmu):
    """Signed joint numerator is formed before normalization/division by mu."""
    ratio=AZ.positive_divide(A,logA);den=A*A;denmu=den*mu
    target={key:{} for key in repair.ROWS};targetZ={key:{} for key in repair.ROWS};joint={};jointZ={}
    for power in ORDERS:
        for out,key,degree in (('M','m',1),('I','h',1),('S','e',2),('Cp','p',2)):
            divisor=A if degree==1 else den
            value=values[key][power].positive_divide(divisor,degree*logA)
            jet=jets[key][power].positive_divide(divisor,degree*logA)-value*ratio*degree
            target[out][power]=value;targetZ[out][power]=jet
        num=values['k'][power]-A*values['m'][power]
        numZ=jets['k'][power]-AZ*values['m'][power]-A*jets['m'][power]
        joint[power]=num;jointZ[power]=numZ
        target[repair.ROWS[1]][power]=num.positive_divide(denmu,2*logA+logmu)
        targetZ[repair.ROWS[1]][power]=(numZ-num*ratio*2).positive_divide(denmu,2*logA+logmu)
    return dict(values=target,Z_derivatives=targetZ,joint_numerator=joint,joint_numerator_Z=jointZ)


def target_caps(target):
    c=next(iter(target['values']['M'].values())).ctx;caps={}
    for key in repair.ROWS:
        terms=[]
        for power in ORDERS:
            cap=repair.LogUpper.add(c,[magnitude(target['values'][key][power]),magnitude(target['Z_derivatives'][key][power])])
            if power==-2:cap=cap*repair.LogUpper.constant(c,c.mpf(1)/MIN_N)
            terms.append(cap)
        caps[key]=repair.LogUpper.add(c,terms)
    nonzero=[cap for cap in caps.values() if cap.log is not None]
    maximum=repair.LogUpper(c,c.mpf(max(ep(cap.log)[1] for cap in nonzero))) if nonzero else repair.LogUpper(c,None)
    return dict(transformed_N_scaled_target_C1_caps={key:cap.record() for key,cap in caps.items()},
        whole_target_C1_cap=maximum.record(),norm='max_i(sup_Z|N*r_i|+sup_Z|N*r_i_Z|)',
        uniform_for_integer_N_ge=MIN_N,joint_divided_row_cancellation_not_proved_by_interval_hulls=True)


def decode(c,value):
    if isinstance(value,dict):
        if 'lower' in value and 'upper' in value:return packets.interval(c,value)
        return {key:decode(c,item) for key,item in value.items()}
    if isinstance(value,list):return [decode(c,item) for item in value]
    return value


class NativeRcParameterTargets:
    def __init__(self,owner):
        if type(owner) is not preceding.NativeRcC1Histories:raise ValueError('Same accepted original Rc owner required')
        receipt=json.loads((HERE/preceding.RECEIPT).read_bytes())
        if not receipt['all_passed'] or not receipt[preceding.GATE] or receipt['source_family']!=owner.family:
            raise ValueError('Checked actual original no-gap Rc route required')
        self.owner=owner;self.ctx=owner.ctx;self.family=owner.family;self.coordinates=owner.coordinates
        self.transfer=owner.transfer;self.service=owner.service;self.q_owner=owner.q_owner
        self.service.bind_hashes(receipt['input_hashes'])
        self.service.bind_hashes({name:sha(name) for name in (preceding.NAME,preceding.RECEIPT,Path(__file__).name,
            repair.NAME,repair.RECEIPT)})
        accepted=json.loads((HERE/repair.RECEIPT).read_bytes());manifest=json.loads((HERE/repair.NAME).read_bytes())
        if not accepted['all_passed'] or not accepted[repair.GATE] or manifest['source_family']!=self.family:
            raise ValueError('Same original checked exact five-bump inverse required')
        self.service.bind_hashes(accepted['input_hashes']);self.repair_manifest=decode(self.ctx,manifest)
        name=PREFIX+'current_generic_shear_O3_sources.json'
        O3=json.loads((HERE/name).read_bytes())
        if O3['source_family']!=self.family:raise ValueError('Original repair mu source family required')
        self.repair_mu=packets.interval(self.ctx,O3['original_O3_quotient_log_norms']['O3_power']['actual_positive_denominator_theorem']['actual_positive_mu'])
        self.service.bind_hashes({name:sha(name)})

    def geometry(self,label,chart,left,right):
        c=self.ctx;geometry=self.transfer.geometry;sc=geometry.binder.sc
        if label=='initial_flat_collar':return geometry.cell(chart,{'selected_sc_multiple':'1/2'},{'selected_sc_multiple':'3/4'})
        if label=='active_first_bridge':
            a=sc*3/4;return geometry.build(chart,0,{'bridge':1-a},c.mpf((ep(a)[0],1)),[{'selected_sc_multiple':'3/4'},'1'])
        if chart=='actual_patch':return geometry.build(chart,1,{},c.mpf((1,ep(c.exp(1))[1])),['1','original exp(1)'])
        if chart=='O3_power':return geometry.cell(chart,{'original_power_offset':left},{'original_power_offset':right})
        return geometry.cell(chart,left,right)

    @native.inlet.source_precision
    def route(self,Z=(-1,1)):
        c=self.ctx;coords=self.coordinates;root_owner=self.q_owner.owner.owner;signed=self.transfer.owner.signed_owner
        eta=packets.interval(c,root_owner.scales['selected_positive_eta_log'])
        dstar=packets.interval(c,root_owner.scales['logarithmic_selected_positive_lower_constants']['d_star'])
        # The original inlet is exactly unmodified for every N. N here is
        # only the old geometry API argument; no evaluated defect is lifted.
        inlet=self.transfer.owner.inlet(Z,N=1024)
        if any(not v.zero for v in [*inlet['initial_defect'].values(),*inlet['initial_defect_Z'].values()]):
            raise ValueError('Original zero initial correction source required')
        zeros=lambda:{key:{power:coords.scalar(0) for power in ORDERS} for key in RATES}
        cumulative,cumulativeZ=zeros(),zeros();cells=[];previous='bridge_first';seams=[]
        excess=packets.interval(c,self.owner.domain['left_collar_kappa_minus2_lower'])
        if ep(eta)[1]>ep(c.ln(excess)-c.ln(2))[0]:raise ValueError('Original whole initial collar must be flat')
        for label,chart,left,right in ROUTE:
            geom=self.geometry(label,chart,left,right)
            if chart!=previous:
                seam=previous+' -> '+chart
                if seam not in geom['record']['exact_original_radius_Jacobian_identities']['original_same_radius_periodic_phase_seam_identities']:
                    raise ValueError('Exact original route seam missing: '+seam)
                seams.append(seam)
            factors={key:transfer.true_width_kernel(coords,geom,rate) for key,rate in RATES.items()}
            if label=='initial_flat_collar':
                values,jets=zeros(),zeros();primitive=None
                source_record=dict(original_full_initial_collar_flat_receipt=transfer.RECEIPT,
                    original_inlet_and_initial_collar_source_bound=True,original_cutoff_and_primitive_C0_Z_exact_zero=True)
                got=dict(values=zeros(),Z_derivatives=zeros(),record=dict(original_initial_flat_support=True))
            else:
                source=self.q_owner.query(chart,Z,geom['coordinate']);roots=source['source']['roots'];packet=source['source']['packet']
                positive=root_owner.decode(root_owner.inventory[chart]['actual_positive_denominator_theorem'])
                key='whole_actual_source_positive_not_inferred_from_saved_denominator_box' if chart.startswith('O3_') else 'source_function_positivity_not_inferred_from_saved_box'
                if positive.get(key) is not True:raise ValueError('Original chart-uniform source positivity required: '+chart)
                primitive=serial.whole_period_C1(roots,eta,positive['log_actual_a_positive_lower'],dstar)
                def axial(k):
                    row=prior.signed.ordinary_axial_coefficient(packet.velocity['axial'][0],k)
                    return signed.leaf(prior.signed.expressions.RadiusPolynomial(packet.algebra,{0:row}),roots['E'][ZERO].scale.bases,roots['E'][ZERO].ledger)
                got=uniform_density_orders(roots['E'][ZERO],roots['E'][DZ],axial(0),axial(1),primitive['values'])
                values={key:{power:coords.rebase(got['values'][key][power],self.family)*factors[key]['mass'] for power in ORDERS} for key in RATES}
                jets={key:{power:coords.rebase(got['Z_derivatives'][key][power],self.family)*factors[key]['mass'] for power in ORDERS} for key in RATES}
                source_record=dict(source_family=self.family,chart=chart,source_provenance=packet.provenance,
                    original_chart_uniform_positive_theorem=positive,
                    original_source_root_C0_Z={key:{str(order):row[order].record() for order in (ZERO,DZ)} for key,row in roots.items() if key in ('E','a','b','p2','kappa_minus2')},
                    original_new_periodic_C1_cover=primitive['record'],
                    original_old_q_status=source['record'].get('status'),unresolved_old_q_rows_not_used=True)
            incoming,incomingZ=cumulative,cumulativeZ
            cumulative={key:{power:factors[key]['decay']*incoming[key][power]+values[key][power] for power in ORDERS} for key in RATES}
            cumulativeZ={key:{power:factors[key]['decay']*incomingZ[key][power]+jets[key][power] for power in ORDERS} for key in RATES}
            record=dict(label=label,chart=chart,original_geometry=geom['record'],original_signed_source_C1=source_record,
                original_all_N_density_order_proof=got['record'],
                actual_true_width_density_C0_order_covers=records(values),actual_true_width_density_Z_order_covers=records(jets),
                actual_inherited_C0_order_covers=records(incoming),actual_inherited_Z_order_covers=records(incomingZ),
                actual_right_correction_C0_order_covers=records(cumulative),actual_right_correction_Z_order_covers=records(cumulativeZ),
                all_N_covers_rebuilt_from_original_source_and_new_primitive_bounds=True,
                true_log_radius_mass_applied_once=True,no_incoming_reset=True)
            cells.append(dict(record=record,geometry=geom,factors=factors,values=values,Z_derivatives=jets,
                density=got,primitives=primitive,incoming=incoming,incoming_Z=incomingZ,
                cumulative=cumulative,cumulative_Z=cumulativeZ));previous=chart
        endpoint=2/self.transfer.geometry.binder.fixed['Tw'];amplitude_source=self.q_owner.query('O3_power',Z,endpoint)
        roots=amplitude_source['source']['roots']
        logA=packets.interval(c,self.owner.reservation['positive_Ac_over_S_log_lower'])
        A=coords.rebase(roots['E'][ZERO],self.family).positive_intersection(logA)
        AZ=coords.rebase(roots['E'][DZ],self.family)
        mu=c.mpf(packets.interval(c,self.owner.domain['right_collar_mu']));logmu=c.mpf(ep(c.ln(mu))[0]);mu_source=coords.scalar(mu)
        if mu._mpi_!=self.repair_mu._mpi_:raise ValueError('Native target and checked exact repair must use the same original mu')
        cover=self.repair_manifest['new_repair_geometry']['actual_mu_cover_for_constants_only']
        if ep(cover)[0]>ep(mu)[0] or ep(cover)[1]<ep(mu)[1]:raise ValueError('Original mu outside checked inverse bounds')
        target=target_rows(cumulative,cumulativeZ,A,AZ,logA,mu_source,logmu);caps=target_caps(target)
        W=self.repair_manifest['fresh_exact_weight_definitions_and_enclosures'];matrix=self.repair_manifest['fresh_divided_linear_inverse_and_enclosures']
        conditions=repair.contraction_log_conditions(c,caps,matrix,W,logmu,c.ln(MIN_N))
        # Identify a dominant contribution after the actual remaining quiet
        # and active widths. These are diagnostics of bounds, not field values.
        suffix={key:coords.scalar(1) for key in RATES};dominants={};weighted=[]
        for cell in reversed(cells):
            v={key:{power:cell['values'][key][power]*suffix[key] for power in ORDERS} for key in RATES}
            j={key:{power:cell['Z_derivatives'][key][power]*suffix[key] for power in ORDERS} for key in RATES}
            part=target_caps(target_rows(v,j,A,AZ,logA,mu_source,logmu))
            weighted.append(dict(label=cell['record']['label'],transformed_N_scaled_target_C1_caps=part['transformed_N_scaled_target_C1_caps']))
            for key,row in part['transformed_N_scaled_target_C1_caps'].items():
                if row['exact_zero']:continue
                upper=ep(packets.interval(c,row['log_absolute_upper']))[1]
                if key not in dominants or upper>dominants[key][0]:dominants[key]=(upper,cell['record']['label'])
            suffix={key:cell['factors'][key]['decay']*suffix[key] for key in RATES}
        contract=dict(normalized_incoming_rows='r=(delta_m/A,(delta_k-A*delta_m)/(mu*A^2),delta_h/A,delta_e/A^2,delta_p/A^2)',
            amplitude='A=Ac/S=original Rc E; A_Z from same source jet',
            joint_numerator='C(N,Z)=delta_k(N,Z)-A(Z)*delta_m(N,Z)',
            joint_numerator_Z='C_Z=delta_k_Z-A_Z*delta_m-A*delta_m_Z',
            divided_row_Z='(C_Z-2*(A_Z/A)*C)/(mu*A^2)',
            shared_original_signed_density_functions_and_Duhamel_tree_retained=True,
            common_basis_does_not_restore_lost_interval_correlations=True,
            joint_divided_row_quantitative_cancellation_not_established=True,
            repair_increment_target='-r',exact_original_implicit_equation='B_exact(mu)*h+N*r(N,Z)+Q_exact(mu,h)/N=0',
            P0_and_P0_Z_unchanged=True,actual_control_functions_installed=False,
            actual_terminal_Z_function_closure_installed=False)
        record=dict(source_family=self.family,Z_box=c.mpf(Z),all_integer_N_lower=MIN_N,
            actual_original_route=[dict(label=label,chart=chart,left=left,right=right) for label,chart,left,right in ROUTE],
            actual_original_initial_inlet= inlet['record'],actual_original_all_N_serial_cells=[cell['record'] for cell in cells],
            actual_Rc_correction_C0_order_covers=records(cumulative),actual_Rc_correction_Z_order_covers=records(cumulativeZ),
            original_Rc_amplitude=dict(actual_endpoint_phase=endpoint,original_positive_log_lower=logA,
                A=A.record(),A_Z=AZ.record(),mu=mu,log_mu_lower=logmu,
                checked_native_and_exact_repair_original_mu_equal=True,
                same_original_source_provenance=amplitude_source['source']['packet'].provenance),
            actual_original_signed_repair_target_C0_order_covers=records(target['values']),
            actual_original_signed_repair_target_Z_order_covers=records(target['Z_derivatives']),
            original_joint_divided_row_numerator_C0_orders={str(power):value.record() for power,value in target['joint_numerator'].items()},
            original_joint_divided_row_numerator_Z_orders={str(power):value.record() for power,value in target['joint_numerator_Z'].items()},
            actual_uniform_N_scaled_repair_C1_caps=caps,actual_source_target_repair_log_conditions=conditions,
            original_joint_target_function_contract=contract,
            original_seams=seams,complete_original_chart_count=len(set(chart for _,chart,_,_ in ROUTE)),
            original_true_cell_count=len(cells),uniform_new_periodic_bounds_applied_to_entire_route=True,
            actual_weighted_cell_target_diagnostics=list(reversed(weighted)),
            dominant_target_bound_cell_labels={key:label for key,(_,label) in dominants.items()},
            no_log_N_exponential_or_huge_integer_materialized=True,
            one_global_finite_N_or_control_field_admitted=False,global_inlet_to_Rc_histories_admitted=False,
            **dict.fromkeys(packets.OPEN,False))
        return dict(record=record,cells=cells,values=cumulative,Z_derivatives=cumulativeZ,amplitude=A,amplitude_Z=AZ,
            logA=logA,mu=mu_source,logmu=logmu,targets=target,caps=caps,conditions=conditions)


@native.inlet.source_precision
def run(bridge=None):
    began=time.monotonic()
    if bridge is None:bridge,_=native.inlet.native_bridge_owner()
    with native.inlet.CheckedSourceRuntime():
        owner=NativeRcParameterTargets(preceding.NativeRcC1Histories(preceding.preceding.NativeO2C1Histories(
            preceding.preceding.make_middle_owner(bridge))));out={}
        for name,Z in (('whole_Z',(-1,1)),('Z_interval',('.49','.51'))):
            out[name]=owner.route(Z)['record'];print('Actual all-N Rc correction/repair target orders:',name,flush=True)
    result=dict(source_family=owner.family,**{GATE:True},native_Z_query_count=2,all_integer_N_lower=MIN_N,
        actual_original_Rc_parameter_target_records=out,exact_density_order_theorem=exact_density_order_theorem(),
        native_all_N_source_family_rebuilt=True,original_new_periodic_bounds_used_on_entire_route=True,
        original_generic_repair_inverse_conditions_connected_to_native_targets=True,
        actual_control_functions_or_terminal_closure_installed=False,
        one_global_finite_N_or_control_field_admitted=False,global_inlet_to_Rc_histories_admitted=False,
        **dict.fromkeys(packets.OPEN,False),execution_seconds=time.monotonic()-began,input_hashes=dict(owner.service.hashes),
        scope='Same-source all-N>=160 signed N^-1/N^-2 C0/Z coefficient covers over all24 true inlet-to-Rc cells and17 original charts; new periodic bounds used everywhere. Actual original Rc amplitude and correction targets connect to exact generic repair inverse/log conditions. Covers do not define coefficients, restore cross-moment cancellation, select global N or install controls/terminal/higher jets/heat/energy/recursion/full NS.')
    (HERE/NAME).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()
