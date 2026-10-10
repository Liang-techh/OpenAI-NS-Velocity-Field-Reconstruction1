"""Check actual endpoints, radial units, target substitution and C1 FTC.

Checks this new source adapter without replaying the admitted reconstruction
fixtures. Exact source functions and shared atoms, not interval overlap,
identify the complete future integrals.
"""
import json
import time
from pathlib import Path
import sympy as s
import lei_ren_part1_paper_compliant_current_limit_absolute_future_integrals as source


def exact_zero(jet):
    return all(source.joined.history.endpoints(row)==(0,0) for row in jet.coefficients)


@source.joined.source_precision
def run(field=None):
    began=time.monotonic();raw=json.loads((source.HERE/source.NAME).read_bytes())
    assert raw['all_passed'] and raw['candidate_absolute_future_source_constructed']
    assert not any(raw[key] for key in source.GATES+source.OPEN)
    hashes=dict(raw['input_hashes']);family=raw['source_family']
    for name,digest in hashes.items():assert source.sha(name)==digest,name
    field=field if field is not None else source.CurrentLimitAbsoluteFutureIntegrals(require_checked=False)
    assert field.family_record==family and not field.acceptance_loaded
    assert field.pulse is field.owner.history.selected.pulse
    assert field.buffer is field.pulse.pulse.buffer
    assert field.heat is field.owner.history.heat
    assert field.params is field.owner.history.selected.future.params
    assert all(field.owner.object_graph().values())
    encoded=source.joined.serialized
    assert encoded(field.FTC)==raw['actual_current_cumulative_FTC_source_proof']
    assert encoded(field.heat_proof)==raw['actual_current_final_heat_reference_and_C1_terminal_limits']
    assert field.FTC['passed'] and field.heat_proof['passed']
    assert all(field.FTC['identities'].values())
    assert field.FTC['actual_raw_m_k_divided_by_Pstar_exactly_once']
    assert field.FTC['same_native_Rw_Rc_mu_Tw_definitions_consumed']
    assert not field.FTC['interval_overlap_used_as_function_identity']
    # Rebind to the actual production methods, rather than trusting a
    # saved proof's flags or another branch's interval rows.
    frame=json.loads((source.HERE/source.joined.selected.source.inlet.identity.NAME).read_bytes())
    pulse=json.loads((source.HERE/(source.PREFIX+'pulse_interface_certificate.json')).read_bytes())
    assert source.source_FTC_proof(field.owner,frame,field.power_report,pulse)==field.FTC
    assert source.heat_reference_and_limits_proof(field.owner)==field.heat_proof
    graph=field.installed_target_graph
    assert encoded(graph)==raw['actual_final_reference_target_graph']
    original=field.band_report['exact_graph_nodes'];old=graph['replaced_c_infinity_symbol']
    assert original[old]==dict(operation='bound_variable',name='future_source_c_infinity')
    assert graph['original_graph_prefix_length']==len(original)
    recipe=graph['appended_graph_nodes'][0]
    assert graph['actual_c_infinity_recipe']==len(original)
    assert recipe['operation']=='current_selected_final_heat_reference_function_recipe'
    assert recipe['exact_definition']=='Ev0*theta_base*Rtail^((1+delta)/2)'
    assert recipe['source_family']==family and recipe['strictly_positive'] and recipe['Z_independent']
    assert recipe['defining_quantity_not_a_range_value']
    assert recipe['source_AST']==field.heat_proof['actual_amplitude_and_radius_AST']
    for key,rows in graph['actual_source_target_handles'].items():
        for n,node in enumerate(rows):
            row=graph['appended_graph_nodes'][node-len(original)]
            assert row==dict(operation='simultaneous_function_substitution',
                expression=field.required_target_nodes[key][n],variables=[old],
                values=[graph['actual_c_infinity_recipe']],simultaneous=True,source_functions_held_fixed=True)
    # Independent physical scale derivation at Rp/Rv/Rtail. Original
    # paper Jacobians are R, sqrt(2) R^(3/2), R and pressure/Pstar^2.
    L,B,mu,a=s.symbols('quiet_length tail_finite mu a',positive=True)
    radius_logs=dict(Rp=L,Rv=L+13/mu,Rt=B+13/mu)
    theta_logs=dict(Rp=s.Integer(0),Rv=-13/(2*mu)-13,Rt=-13/(2*mu)-13)
    expected={};scale_proofs={}
    for stage in radius_logs:
        rl,ul=radius_logs[stage],theta_logs[stage]
        expected[stage]=dict(M=rl+ul,J=s.Rational(3,2)*rl+2*ul,
            I=s.Rational(3,2)*rl+ul,S_energy=rl+2*ul,pressure_tail=s.Integer(0))
    assert s.expand((1-a)*L-13*a/mu+13/mu-(1-a)*radius_logs['Rv'])==0
    assert s.expand((1-a)*B-13*a/mu+13/mu-(1-a)*radius_logs['Rt'])==0
    original=source.buffer.SharedOuterBuffer.power
    assert source.buffer.fraction_box is original.__globals__['fraction_box']
    assert source.buffer.fraction_box is not source.exact_phase_box
    views={}
    for name,Z in (('fresh','.519'),('axis',0)):
        packet=field.future_integrals(Z);assert encoded(packet)==raw['source_views'][name]
        views[name]=packet;ledger=packet['actual_four_segment_source_ledger'];atoms=ledger['atoms']
        assert set(ledger['segments'])==set(source.SEGMENTS)
        assert set(ledger['normalized_endpoints'])=={'R2','Rp','Rv','Rt'}
        assert len(atoms)==24
        assert packet['exact_entire_future_telescope']['passed']
        assert packet['native_power_code_unchanged_only_real_phase_acquisition_adapted']
        assert source.telescope_proof(ledger)==packet['exact_entire_future_telescope']
        c=field.ctx;length=c.mpf(source.joined.history.endpoints(field.params.Tw))-2-c.ln(2)
        tail_offset=100+field.owner.history.outer.Lrel+2+field.owner.history.steep.Ts+field.owner.history.steep.wait
        Bvalue=length+tail_offset
        for stage in ('Rp','Rv','Rt'):
            for key in ('M','J','I','S_energy','pressure_tail'):
                entry=atoms[stage+'_'+key];parts=entry['exact_scale_log']
                assert entry['source_family']==family
                assert entry['source_provider']['source_sha256'] in hashes.values()
                assert entry['numerical_rows_are_enclosures_not_defining_values']
                assert entry['exact_positive_scale_never_replaced_with_cap']
                if stage=='Rt' and key in ('M','J'):
                    assert exact_zero(entry['actual_C1_coefficient_enclosure']);continue
                exact=s.expand(expected[stage][key])
                inv=s.simplify(s.diff(exact,mu)*(-mu**2))
                finite=s.simplify(exact-inv/mu)
                expected_finite=c.mpf(str(finite)) if not finite.free_symbols else (
                    (c.mpf('1.5') if finite.coeff(L)==s.Rational(3,2) else c.mpf(1))*length
                    +c.mpf(str(finite.subs(L,0))) if L in finite.free_symbols else
                    (c.mpf('1.5') if finite.coeff(B)==s.Rational(3,2) else c.mpf(1))*Bvalue
                    +c.mpf(str(finite.subs(B,0))))
                assert source.joined.history.endpoints(parts['finite'])==source.joined.history.endpoints(expected_finite),(stage,key)
                expected_inv=c.mpf(int(s.numer(inv)))/int(s.denom(inv))
                assert source.joined.history.endpoints(parts['inverse_mu_coefficient'])==source.joined.history.endpoints(expected_inv)
                scale_proofs[stage+'_'+key]=True
            refparts=atoms[stage+'_Ipow']['exact_scale_log']
            ref_finite=field.heat.k*(length if stage!='Rt' else Bvalue)
            if stage!='Rp':ref_finite-=13*field.heat.a/field.heat.mu
            assert source.joined.history.endpoints(refparts['finite'])==source.joined.history.endpoints(ref_finite)
            expected_inverse=(0,0) if stage=='Rp' else (13,13)
            assert source.joined.history.endpoints(refparts['inverse_mu_coefficient'])==expected_inverse
        for key in source.KEYS:
            assert encoded(packet['absolute_future_integrals'][key])==encoded(packet['required_future_integral_C1_in_declared_units'][key])
            assert packet['absolute_future_integrals'][key].order==1
        assert not any(packet[key] for key in source.GATES+source.OPEN)
        for stage in ('Rv','Rt'):
            assert exact_zero(atoms[stage+'_M']['actual_C1_coefficient_enclosure'])
            assert exact_zero(atoms[stage+'_J']['actual_C1_coefficient_enclosure'])
        # The R2 absolute histories retain incoming source functions.
        if name=='fresh':
            assert not exact_zero(packet['current_2Rc_common_own_C1']['m'])
            assert not exact_zero(packet['current_2Rc_common_own_C1']['k'])
        assert not exact_zero(packet['independent_original_P0_C1'])
        assert exact_zero(source.IntervalTaylor(c,[packet['actual_final_reference_primitive_in_R2_units'][1]]))
        assert ledger['same_internal_source_atom_reused_by_adjacent_segments']
    # True axis derivative survives even where odd absolute moment values
    # vanish. This would fail if C1 were replaced by a scalar snapshot.
    assert source.joined.history.endpoints(views['axis']['current_2Rc_common_own_C1']['m'][0])==(0,0)
    assert source.joined.history.endpoints(views['axis']['current_2Rc_common_own_C1']['m'][1])!=(0,0)
    rejected=[]
    for Z in (-2,2):
        try:field.future_integrals(Z)
        except ValueError:rejected.append(Z);continue
        raise AssertionError('Original Z domain not enforced')
    bindings=source.joined.selected.source.inlet.identity.ast_assignments(
        'current_limit_absolute_future_integrals','CurrentLimitAbsoluteFutureIntegrals','__init__',{
            'self.pulse':'self.owner.history.selected.pulse','self.buffer':'self.pulse.pulse.buffer',
            'self.heat':'self.owner.history.heat','self.params':'self.owner.history.selected.future.params'})
    hashes[source.NAME]=source.sha(source.NAME);hashes[Path(__file__).name]=source.sha(Path(__file__).name)
    receipt=dict(source_family=family,all_passed=True,**dict.fromkeys(source.GATES,True),
        **dict.fromkeys(source.OPEN,False),actual_current_owner_AST=bindings,
        actual_final_reference_target_graph_checked=True,
        exact_real_2Rc_phase_uses_native_power_code_without_rational_point_selection=True,
        final_reference_Rv_Rtail_full_inverse_mu_radius_scale_checked=True,
        five_actual_absolute_endpoint_functions_and_four_segments_checked=True,
        independent_paper_radius_and_swirl_scale_log_identities=scale_proofs,
        native_2Rc_backward_unique_C1_source_and_one_Pstar_conversion_checked=True,
        nonzero_incoming_M_J_and_independent_P0_retained=True,
        selected_terminal_M_J_zero_only_after_actual_linear_repair=True,
        actual_complete_full_Gamma_physical_C1_majorants_checked=True,
        entire_four_segment_future_targets_equal_as_functions_and_true_first_Z=True,
        fresh_and_axis_actual_source_C1_views_checked=True,invalid_Z_rejected=rejected,
        numeric_overlap_or_caps_not_used_for_function_equalities=True,
        source_scope='absolute five post-2Rc similarity-source integrals; C1 only',
        patched_Rh_global_physical_C4_energy_cone_and_time_recursion_still_open=True,
        input_hashes=hashes,execution_seconds=time.monotonic()-began)
    (source.HERE/source.RECEIPT).write_text(json.dumps(encoded(receipt),indent=2)+'\n',encoding='utf8')
    print('PASS_CURRENT_LIMIT_ABSOLUTE_FUTURE: five absolute C1 source integrals, actual heat reference, four source segments',flush=True)
    return receipt


if __name__=='__main__':run()
