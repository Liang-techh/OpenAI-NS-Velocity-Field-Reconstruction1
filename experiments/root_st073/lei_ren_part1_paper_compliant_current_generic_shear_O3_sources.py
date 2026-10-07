"""Original O3 source extension, logarithmic input bounds and repair reservation.

Rebuild only the checked raw row arithmetic from saved original source jets.
The power tensor wrapper has different amplitude units and is not the raw
source. No ancestor graph, physical radius, source amplitude or inverse
microscopic width is evaluated. These exports are covers, not point fields.
"""
import gzip
import json
import math
from pathlib import Path
import mpmath as mp
import lei_ren_part1_paper_compliant_current_generic_shear_source_bounds as bounds
from lei_ren_part1_paper_compliant_long_reshape_mixed_C4 import exponential_derivatives
from lei_ren_part1_paper_compliant_collar_stress_C3 import shifted_rows

packets=bounds.packets
HERE,PREFIX,sha=packets.HERE,packets.PREFIX,packets.sha
NAME=PREFIX+'current_generic_shear_O3_sources.json'
RECEIPT=PREFIX+'current_generic_shear_O3_sources_check.json'
VIEWS=PREFIX+'current_generic_shear_O3_sources_views.json.gz'
GATE='current_original_O3_common_unit_packets_and_quotient_log_bounds_certified'
RESERVE_GATE='current_original_O3_right_edge_and_reserved_repair_geometry_certified'
OPEN=packets.OPEN
CHARTS=('O3_slope_mu','O3_power')
SPEC={
 'O3_slope_mu':('current_O3_transition_background_tensor',
     ('current_actual_O3_transition_tensor_views','whole_transition'),
     'actual_upstream_original_pre_transition_source',
     ('current_actual_O3_transition_full_tensor_available',
      'current_actual_O3_transition_power_completed_tensor_join_certified')),
 'O3_power':('current_O3_power_cone',
     ('whole_current_O3_power','original_complete_view'),
     'actual_upstream_original_pre_power_source',
     ('current_whole_O3_power_signed_two_vector_cone_certified',
      'current_O3_power_actual_M_K_X_full_energy_pressure_correlation_bound',
      'current_O3_power_transition_and_entrance_function_joins_consumed'))}


def raw_rows(c,original):
    """Same raw_pre_velocity_rows equations on exact decoded source coefficients."""
    read=lambda v:packets.jet(c,v)
    u0=read(original['Utheta_over_Pstar_axial5_coefficients'])
    logjet=[read(v) for v in original['log_Utheta_ordinary_y_derivatives']]
    velocity=dict(theta=[u0*v for v in exponential_derivatives(logjet)],
        axial=[read(v) for v in original['Uz_ordinary_y_derivative_axial5']],
        radial=shifted_rows([read(v) for v in original['actual_Q_y_derivative_axial4']],c.mpf('.5'),4))
    histories={k:[read(v) for v in values] for k,values in
        original['actual_normalized_primitive_y_derivative_axial5'].items()}
    datum=read(original['original_P0_axial5_coefficients'])
    pressure=[histories['p'][0]+datum]+histories['p'][1:]
    return velocity,histories,pressure,datum


def original_row_theorem():
    asts=packets.recovery.numeric.transport.SourceAST()
    for target,wanted in (
        ('u', '[u0*value for value in exponential_derivatives(logjet)]'),
        ('V', "[copy_jet(c,value) for value in pre['Uz_ordinary_y_derivative_axial5']]"),
        ('Q', "[copy_jet(c,value) for value in pre['actual_Q_y_derivative_axial4']]"),
        ('radial', "shifted_rows(Q,c.mpf('.5'),4)")):
        asts.expression('current_pre_pulse_stress_operator','raw_pre_velocity_rows',target,wanted=wanted)
    asts.expression('pre_pulse_mixed_C4','packet','data',
        wanted='physical_mixed(c,Z,self.delta,u,logU,V,history,p0,self.invP2)')
    asts.expression('pre_pulse_mixed_C4','slope_mu','factor',wanted="c.exp(-t/2-mu*K['J'])")
    asts.expression('pre_pulse_mixed_C4','power','parent',wanted='self.slope_mu(Z,1)')
    asts.expression('pre_pulse_mixed_C4','power','t',wanted='self.params.Tw*phase')
    asts.expression('pre_pulse_mixed_C4','power','f',wanted="c.exp((-c.mpf('.5')-mu)*t)")
    asts.method('pre_pulse_mixed_C4','physical_mixed')
    return dict(passed=True,original_source_AST_bindings=asts.bindings,
        source_histories_are_original_m_h_k_e_p_not_power_specific_m_n_e_P_X=True,
        absolute_pressure_is_original_p_plus_separate_original_P0=True,
        original_radial_physical_half_shift_applied_once=True,
        power_wrapper_amplitude_units_not_used_as_raw_source=True,input_hashes=asts.hashes)


class CurrentO3Sources:
    def __init__(self,service=None):
        self.service=service if service is not None else packets.CurrentSourcePackets()
        self.ctx=self.service.ctx;self.family=self.service.family;self.cache={}
        self.theorem=original_row_theorem();self.service.bind_hashes(self.theorem['input_hashes'])
        record=json.loads((HERE/bounds.RECEIPT).read_bytes())
        if not record.get('all_passed') or not record.get(bounds.GATE) or record['source_family']!=self.family:
            raise ValueError('Same checked current inner-to-O2 source bounds required')
        self.service.bind_hashes(record['input_hashes'])
        self.service.bind_hashes({bounds.RECEIPT:sha(bounds.RECEIPT),Path(__file__).name:sha(Path(__file__).name)})
        self.previous=json.loads((HERE/bounds.NAME).read_bytes())
        self.previous_lower=self.previous['positive_noncore_source_denominator_theorem']['source_charts']['O2_buffer']

    def saved_view(self,chart):
        if chart not in SPEC:raise ValueError('Original O3 source chart required')
        if chart in self.cache:return self.cache[chart]
        stem,path,original_key,gates=SPEC[chart]
        receipt=PREFIX+stem+'_check.json'
        admitted=json.loads((HERE/receipt).read_bytes())
        if not admitted.get('all_passed') or not all(admitted.get(g) for g in gates):
            raise ValueError('Checked original whole O3 source/cone receipt required')
        if {k:admitted[k] for k in packets.FAMILY_KEYS}!=self.family or any(admitted.get(k) for k in OPEN):
            raise ValueError('Original O3 source family or regional scope differs')
        self.service.bind_hashes(admitted['input_hashes']);self.service.bind_hashes({receipt:sha(receipt)})
        if chart=='O3_power':
            producer=PREFIX+stem+'.json';viewfile=PREFIX+stem+'_views.json.gz'
            data=json.loads((HERE/producer).read_bytes())
            views=json.loads(gzip.decompress((HERE/viewfile).read_bytes()))
        else:
            producer=PREFIX+stem+'.json.gz';viewfile=producer
            data=json.loads(gzip.decompress((HERE/producer).read_bytes()));views=data
        if {k:data[k] for k in packets.FAMILY_KEYS}!=self.family:
            raise ValueError('Saved original O3 producer has another family')
        self.service.bind_hashes({producer:sha(producer),viewfile:sha(viewfile)})
        # Both files must already be covered by the checked upstream receipt.
        if any(admitted['input_hashes'].get(n)!=sha(n) for n in {producer,viewfile}):
            raise ValueError('Unbound original O3 source view')
        view=packets.path_get(views,path)
        if view['chart']!=chart or not view['source_bounds_not_resolved_physical_point_values']:
            raise ValueError('Original whole source cover required')
        proof=data['current_actual_O3_transition_tensor_and_power_join_theorem'] if chart=='O3_slope_mu' else data['current_actual_O3_power_full_source_composition']
        self.cache[chart]=(view,dict(mode='saved_original_O3_cover',producer=producer,
            producer_sha256=sha(producer),viewfile=viewfile,view_path=path,receipt=receipt,
            receipt_source_function_theorem=proof,original_source_path=original_key,
            cache_cover=True,arbitrary_coordinates_evaluated=False))
        return self.cache[chart]

    def saved(self,chart):
        view,origin=self.saved_view(chart)
        return self.adapt(chart,view,origin)

    def adapt(self,chart,view,origin):
        if chart not in SPEC or view['chart']!=chart or not view.get('source_bounds_not_resolved_physical_point_values'):
            raise ValueError('Original O3 cover with chart provenance required')
        c=self.ctx;original=view[SPEC[chart][2]]
        physical=view['actual_upstream_physical_spatial4_time1_packet']
        if {k:physical[k] for k in packets.FAMILY_KEYS}!=self.family:
            raise ValueError('Original O3 source/axis datum differs')
        source=view['current_actual_source_stress_packet'];logP=packets.interval(c,source['exact_pulse_reference_logB_parts']['logPstar'])
        lo,hi=packets.recovery.endpoints(logP);el,eh=packets.recovery.endpoints(self.service.data['logP'])
        if hi<el or eh<lo:raise ValueError('Original constant Pstar unit differs')
        velocity,histories,pressure,datum=raw_rows(c,original)
        if set(histories)!=set(packets.recovery.RATES):raise ValueError('All five original histories required')
        algebra=packets.FactoredAlgebra(c,(c.mpf(0),2*logP,c.mpf(0),c.mpf(0)),[])
        def rows(values):
            if len(values)!=5 or any(v.order<4 for v in values):
                raise ValueError('Five complete ordinary-y source rows required')
            return tuple(algebra.lift(v) for v in values)
        native_v={k:rows(v) for k,v in velocity.items()};native_m={k:rows(v) for k,v in histories.items()}
        common_v={k:tuple(algebra.shift(v,packets.INVERSE_S) if k in ('axial','radial') else v for v in values) for k,values in native_v.items()}
        common_m={k:tuple(algebra.shift(v,packets.INVERSE_S) if k in ('m','k') else v for v in values) for k,values in native_m.items()}
        Z=packets.IntervalTaylor.variable(c,packets.interval(c,source['Z']),5)
        coordinate=packets.interval(c,view['coverage_coordinate'])
        endpoints=packets.recovery.endpoints
        zl,zh=endpoints(Z[0]);cl,ch=endpoints(coordinate)
        if not all(mp.isfinite(v) for v in (zl,zh,cl,ch)) or zl<-1 or zh>1 or cl<0 or ch>1:
            raise ValueError('Whole O3 axial/coordinate cover must lie in checked domain')
        provenance=dict(origin,source_family=self.family,chart=chart,Z_box=Z[0],coordinate_box=coordinate,
            coordinate_contract='offset=t=log(R/Rd)' if chart=='O3_slope_mu' else 'phase=t/Tw; y rows remain ordinary logR',
            derivative_coordinate='ordinary y=log R',original_radius_source=original['exact_radius_source'],
            logR_cover=packets.interval(c,source['exact_logR']),requested_log_tau=view['requested_log_tau'],
            physical_viscosity=view['physical_viscosity'],raw_record_form='converted_plain_raw_cover',
            original_raw_rows_replayed_from_saved_source_jets=True,
            source_factor_resolution_performed=False,raw_coordinate_conversion_not_reapplied=True,
            P0_read_from_original_source_not_subtracted=True,per_view_bounds_not_a_seam_identity=True,
            power_normalized_tensor_wrapper_not_used_as_raw_velocity=True)
        return packets.CurrentSourcePacket(chart,dict(self.family),algebra,Z,common_v,common_m,
            rows(pressure),algebra.lift(datum),native_v,native_m,provenance)

    def denominator(self,chart):
        c=self.ctx;read=lambda v:packets.interval(c,v);lower=lambda v:c.mpf(packets.recovery.endpoints(v)[0])
        view,_=self.saved_view(chart);original=view[SPEC[chart][2]];mu=read(original['selected_positive_mu'])
        if not 0<packets.recovery.endpoints(mu)[0]<=packets.recovery.endpoints(mu)[1]<c.mpf('1/6'):
            raise ValueError('Same original positive small mu required')
        parent=read(self.previous_lower['log_E_positive_lower'])
        ell=parent-c.mpf('.5')-mu
        if chart=='O3_power':ell-=(c.mpf('.5')+mu)*read(original['total_log_length'])
        return dict(log_E_positive_lower=lower(ell),log_C_positive_lower=lower(ell+c.ln(2)),
            log_actual_a_positive_lower=lower(c.ln(2)),actual_positive_mu=mu,
            transition_J_bound='0<=integral_0^t sigma(s)ds<=1 on offset[0,1]',
            original_parent='same checked O2 buffer at offset11; no new amplitude seed',
            whole_actual_source_positive_not_inferred_from_saved_denominator_box=True,
            exact_a_minus2_source='2*mu*sigma(t)' if chart=='O3_slope_mu' else '2*mu',
            transition_left_a_minus2_can_be_zero=chart=='O3_slope_mu')

    def log_bounds(self,packet,field):
        c=self.ctx;E,V=packet.velocity['theta'],packet.velocity['axial']
        C=tuple(E[j]-2*E[j+1] for j in range(3));B=tuple(2*V[j+1] for j in range(3))
        stress=field['full_signed_stress_ordinary_y_rows']
        def numerator(component):
            physical=[left+packet.algebra.shift(right,(0,.5,0,0)) for left,right in
                zip(stress['inertial_'+component+'_linear'],stress['inertial_'+component+'_quadratic'])]
            return tuple(sum((physical[i]*(math.comb(j,i)*c.mpf('.5')**(j-i)) for i in range(j+1)),
                packet.algebra.lift(0)) for j in range(3))
        tables=dict(E=bounds.row_table(packet,E),V=bounds.row_table(packet,V),
            C=bounds.row_table(packet,C),B=bounds.row_table(packet,B),
            nt=bounds.row_table(packet,numerator('theta'),radius_power=1),
            nz=bounds.row_table(packet,numerator('axial'),radius_power=1))
        positive=self.denominator(packet.chart);quotients={}
        for key,num,den,ell in (('a','C','E','log_E_positive_lower'),('b','B','E','log_E_positive_lower'),
            ('p1','nt','E','log_E_positive_lower'),('p2','nz','E','log_E_positive_lower'),('t0','B','C','log_C_positive_lower')):
            quotients[key]=bounds.quotient_table(c,tables[num],tables[den],positive[ell])
        return dict(actual_positive_denominator_theorem=positive,
            ordinary_mixed_source_log_norms={k:bounds.encode_table(v) for k,v in tables.items()},
            original_five_history_log_norms={k:bounds.encode_table(bounds.row_table(packet,v)) for k,v in packet.histories.items()},
            absolute_pressure_log_norms=bounds.encode_table(bounds.row_table(packet,packet.absolute_pressure)),
            separate_P0_log_norms=bounds.encode_table({(0,k):bounds.modal_partial_bound(packet,packet.P0,k) for k in range(2)}),
            admitted_original_quotient_log_norms={k:bounds.encode_table(v) for k,v in quotients.items()},
            ordinary_y_Z_orders=[list(v) for v in bounds.ORDERS],
            full_signed_inertial_pressure_energy_radial_source_retained=True,
            source_factors_radius_width_or_inverse_denominator_not_materialized=True)

    def right_reservation(self):
        c=self.ctx;view,_=self.saved_view('O3_power');original=view[SPEC['O3_power'][2]]
        Tw=packets.interval(c,original['total_log_length']);mu=packets.interval(c,original['selected_positive_mu'])
        endpoint=2+c.ln(2);margin=Tw-endpoint
        if packets.recovery.endpoints(margin)[0]<=0:raise ValueError('Reserved repair interval must lie strictly inside original power chart')
        return dict(log_radius_offsets_from_same_Rw=dict(r_plus=c.mpf(1),Rc=c.mpf(2),twice_Rc=endpoint,Rp=Tw),
            original_power_phase_locations=dict(r_plus=1/Tw,Rc=2/Tw,twice_Rc=endpoint/Tw),
            positive_log_radius_gaps=dict(r_plus_to_Rc=c.mpf(1),Rc_to_twice_Rc=c.ln(2),twice_Rc_to_Rp=margin),
            original_reserved_profile='Utheta=Ac_theta(Z)*(R/Rc)^(-1/2-mu), Uz=0',
            original_Ac_definition='Ac_theta(Z)=Utheta(Rw,Z)*exp(-2*(1/2+mu))',
            positive_Ac_over_S_log_lower=self.denominator('O3_slope_mu')['log_E_positive_lower']-2*(c.mpf('.5')+mu),
            canonical_a='2+2*mu',canonical_b_exact_zero=True,strict_kappa_minus2_source=2*mu,
            strict_excess_not_computed_by_subtracting_rounded_a_minus2=True,
            original_whole_power_cone_consumed_on_r_plus_through_twice_Rc=True,
            original_whole_power_cone_receipt=PREFIX+'current_O3_power_cone_check.json',
            right_edge_and_repair_have_strict_original_power_collars=True,
            new_repair_controls_or_terminal_closure_constructed=False,
            left_r_minus_and_admissibility_through_Rb_still_required=True,
            complete_Section11_loop_domain_certified=False)

    def run(self):
        covers={};norms={}
        for chart in CHARTS:
            packet=self.saved(chart);field=packet.recover_original(self.service.data['delta'])
            covers[chart]=dict(packet=packet.record(),original_recovered_field=field)
            norms[chart]=self.log_bounds(packet,field)
            print('Original O3 raw packet and complete signed quotient bounds: '+chart,flush=True)
        reservation=self.right_reservation()
        (HERE/VIEWS).write_bytes(gzip.compress((json.dumps(packets.encode(covers),separators=(',',':'))+'\n').encode(),compresslevel=6,mtime=0))
        self.service.bind_hashes({VIEWS:sha(VIEWS)})
        result=dict(source_family=self.family,original_raw_source_row_theorem=self.theorem,
            common_unit_theorem=packets.unit_theorem(),original_O3_quotient_log_norms=norms,
            right_edge_and_new_repair_reservation=reservation,normalized_original_O3_views=VIEWS,
            original_extension_chart_count=2,quotient_mixed_derivative_bound_count=60,
            combined_saved_original_source_cover_chart_count=18,
            combined_positive_quotient_bound_chart_count=17,
            **{GATE:True,RESERVE_GATE:True},**dict.fromkeys(OPEN,False),
            whole_generic_scales_instantiated=False,phase_held_loop_primitive_derivative_bounds_certified=False,
            whole_upstream_source_derivative_norms_certified=False,
            source_graph_ancestor_constructors_called=False,
            scope='Two original whole O3 source covers and full signed y2/Z1 quotient bounds; strict right edge and new repair geometry inside original power. No whole loop domain, new modified histories/repair, common N, recursion or corrected NS admission.',
            input_hashes=self.service.hashes)
        (HERE/NAME).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
        return result


def run():return CurrentO3Sources().run()


if __name__=='__main__':run()
