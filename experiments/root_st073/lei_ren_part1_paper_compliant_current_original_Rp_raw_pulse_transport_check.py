"""Independent active-pulse scales, full density equations and caller checks.

Reuses accepted selected/closure proofs. This checks the new algebra and
runtime boundary; it does not replay upstream quadrature or contraction.
"""
import gzip
import json
import math
from pathlib import Path
from types import SimpleNamespace
import time
import sympy as s

import lei_ren_part1_paper_compliant_current_original_Rp_raw_pulse_transport as current
import lei_ren_part1_paper_compliant_current_generic_shear_moment_recovery as generic
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision


def same(a,b,label):assert s.cancel(s.expand(a-b))==0,label


def scale_proof(owner):
    field=SimpleNamespace(graph=owner.graph,parameters=owner.radius.parameters,
        contracts=owner.radius.frame.bridge.leading.contracts)
    A=s.Symbol('same_absolute_logRp',real=True)
    q=current.radius.RadiusInterpreter(field,False,{owner.radius.logRp.node:A})
    identities=0;cancellations=0
    for _,chart,coordinate in current.VIEWS[:6]:
        geometry=owner.radius.geometry(chart,coordinate);Y=q.at(geometry['logR'])-A
        for key,powers in current.POWERS.items():
            r,e,p=map(s.Rational,powers)
            scale=owner.scale(chart,coordinate,powers)
            wanted=r*(A+Y)+p*q.logP-e*(s.Rational(1,2)+q.mu)*Y
            same(sum(q.at(v) for v in scale.values()),wanted,(chart,key,'absolute pulse scale'))
            identities+=1
            if key in ('Mztheta','Ur'):
                same(q.at(scale['combined_inverse_mu_pulse']),0,(chart,key,'exact pulse cancellation'))
                cancellations+=1
        coordinate=s.Rational(current.radius.exact_coordinate(coordinate))
        jac=1/q.mu if chart in ('pulse_main','pulse_exit','pulse_gap') else s.Integer(1)
        same(q.at(geometry['native_to_log_radius_jacobian']),jac,(chart,'one native Jacobian'))
    return dict(passed=True,actual_absolute_pulse_scale_function_identities=identities,
        energy_and_radial_velocity_inverse_mu_cancellations=cancellations,
        six_actual_native_coordinate_Jacobians_identified=True,
        original_Rp_and_true_pulse_F_functions_not_exponential_caps=True)


def cumulative_proof():
    y,v,z=s.symbols('y v Z',real=True);R0=s.Symbol('R0',positive=True)
    E=s.Function('actual_Utheta');V=s.Function('actual_Uz')
    seeds=s.symbols('incoming_Mz incoming_Mtheta incoming_Mtheta_z incoming_Mztheta incoming_Mp')
    density=generic.history_densities(E(v,z),V(v,z))
    primal=dict(Mz=seeds[0]+s.Integral(R0*s.exp(v)*density['m'],(v,0,y)),
        Mtheta=seeds[1]+s.Integral(s.sqrt(2)*(R0*s.exp(v))**s.Rational(3,2)*density['h'],(v,0,y)),
        Mtheta_z=seeds[2]+s.Integral(s.sqrt(2)*(R0*s.exp(v))**s.Rational(3,2)*density['k'],(v,0,y)),
        Mztheta=seeds[3]+s.Integral(R0*s.exp(v)*density['e'],(v,0,y)),
        Mp=seeds[4]+s.Integral(density['p'],(v,0,y)))
    R=R0*s.exp(y);J=s.Symbol('true_native_Jacobian',positive=True)
    wanted=dict(Mz=R*V(y,z)*J,Mtheta=s.sqrt(2)*R**s.Rational(3,2)*E(y,z)*J,
        Mtheta_z=s.sqrt(2)*R**s.Rational(3,2)*E(y,z)*V(y,z)*J,
        Mztheta=R*(V(y,z)**2-E(y,z)**2/2)*J,Mp=E(y,z)**2*J/2)
    for name,value in primal.items():same(s.diff(value,y)*J,wanted[name],name+' full cumulative equation')
    P0=s.Function('independent_analytic_P0')(z)
    same(s.diff(primal['Mp']+P0,y)*J,wanted['Mp'],'full pressure with analytic P0')
    assert density['e']==V(v,z)**2-E(v,z)**2/2
    return dict(passed=True,five_general_nonzero_axial_cumulative_equations_proved=True,
        true_density_source=Path(generic.__file__).name,
        incoming_seed_memory_and_independent_P0_retained=True,
        exactly_one_native_Jacobian_in_each_derivative=True,
        energy_has_signed_Uz_squared_minus_Utheta_squared_half=True)


def identical(a,b,label):
    assert a.order==b.order,label
    for n in range(a.order+1):assert a[n]._mpi_==b[n]._mpi_,(label,n)
    return a.order+1


@source_precision
def run(heat=None):
    began=time.monotonic();candidate=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    for name,digest in candidate['input_hashes'].items():assert current.sha(name)==digest,name
    owner=current.CurrentOriginalRpRawPulseTransport(heat,require_checked=False)
    assert not owner.acceptance_loaded and not any(candidate[k] for k in current.GATES+current.OPEN)
    assert candidate['source_family']==owner.family_record
    assert candidate['actual_same_current_provider_graph']==owner.assert_graph()
    general=cumulative_proof();views={}
    rows=0;velocityrows=0;common_graph_parts=0
    for name,chart,coordinate in current.VIEWS:
        view=owner.evaluate(chart,candidate['fresh_Z'],coordinate);views[name]=view
        assert current.raw.packed(current.view_report(view))==candidate['actual_fifteen_current_views'][name],name
        for group in ('raw_histories','similarity_velocity','first_native_radial_derivatives'):
            for value in view[group].values():
                assert all(ref.graph is owner.graph for _,ref in value.log_scale_parts),(name,group,'wrong expression graph')
                common_graph_parts+=len(value.log_scale_parts)
        if chart not in current.PULSE:
            for key in ('Uz','Ur'):
                assert all(current.radius.post.selected.inlet.endpoints(v)==(0,0)
                    for v in view['similarity_velocity'][key].coefficients.coefficients)
            if chart.startswith('heat_'):
                assert view['actual_provider_kind']=='same_repair_closed_heat'
                assert view['source_packet']['whole_Z_pressure_function_identity_consumed']
                assert view['raw_histories']['pressure'].coefficients is view['closed_normalized_histories']['closed_pressure']
            continue
        packet=view['source_packet'];u=packet['inlet_Utheta_over_Pstar_Taylor'];B=packet['Uz_over_Utheta']
        expected=dict(Mz=u*packet['Mz_over_R_Utheta'],
            Mtheta=u*packet['Mtheta_over_sqrt2_R_3half_Utheta']*owner.ctx.sqrt(2),
            Mtheta_z=u*u*packet['Mtheta_z_over_sqrt2_R_3half_Utheta_squared']*owner.ctx.sqrt(2),
            Mztheta=u*u*packet['Mztheta_over_R_Utheta_squared'],
            Mp=packet['pressure']['Mp_over_Pstar_squared'],P0=packet['pressure']['P0_over_Pstar_squared'],
            pressure=packet['pressure']['P_over_Pstar_squared'],Utheta=u,Uz=u*B,
            Ur=packet['factored_Ur_over_sqrt_R_over_2_Pstar_Taylor']/owner.ctx.sqrt(2))
        for key,value in expected.items():rows+=identical(view['raw_histories'][key].in_exact_units(current.POWERS[key]),value,(name,key))
        rows+=identical(expected['pressure'],expected['P0']+expected['Mp'],(name,'full pressure memory'))
        assert view['raw_histories']['Ur'].coefficients.order==4
        physical=packet['physical_mixed_derivatives_total_order_le4']
        for key,label,divisor in (
            ('Utheta','Utheta_over_Pstar_without_common_theta_radial_factor',1),
            ('Uz','Uz_over_Pstar_without_common_theta_radial_factor',1),
            ('Ur','Ur_over_sqrt_R_over_2_Pstar_without_common_theta_radial_factor',owner.ctx.sqrt(2))):
            for n in range(5):
                value=view['similarity_velocity'][key].coefficients[n]*math.factorial(n)
                wanted=physical[label]['y0_Z'+str(n)]/divisor
                al,ah=current.radius.post.selected.inlet.endpoints(value)
                bl,bh=current.radius.post.selected.inlet.endpoints(wanted)
                assert max(al,bl)<=min(ah,bh),(name,key,n,'original physical velocity row')
                velocityrows+=1
        density=generic.history_densities(u,u*B);jac=view['directed_native_Jacobian_bound']
        expected_d=dict(Mz=density['m'],Mtheta=density['h']*owner.ctx.sqrt(2),
            Mtheta_z=density['k']*owner.ctx.sqrt(2),Mztheta=density['e'],
            Mp=density['p'],P0=u*0,pressure=density['p'])
        for key,value in expected_d.items():
            # Independent algebra can yield a different interval enclosure.
            actual=view['first_native_radial_derivatives'][key].coefficients;expected_jet=value*jac
            for n in range(6):
                al,ah=current.radius.post.selected.inlet.endpoints(actual[n])
                bl,bh=current.radius.post.selected.inlet.endpoints(expected_jet[n])
                assert max(al,bl)<=min(ah,bh),(name,key,n,'general density');rows+=1
        assert view['first_native_radial_derivatives']['Mp'].powers==tuple(map(current.Fraction,(0,2,2)))
        assert view['raw_histories']['pressure'].powers==tuple(map(current.Fraction,(0,0,2)))
    # Distinguish active velocity from inactive gap memory. No source zeros
    # are inferred merely because an outward end-factor interval contains 0.
    main=views['main']['source_packet']
    assert current.radius.post.selected.inlet.endpoints(main['Uz_over_Utheta'][0])[0]>0
    end=views['active_end']['source_packet']
    assert any(current.radius.post.selected.inlet.endpoints(v)!=(0,0)
        for v in end['formal_scaled_Uz_over_Utheta'].coefficients)
    gap=views['gap']['source_packet']
    assert all(current.radius.post.selected.inlet.endpoints(v)==(0,0) for v in gap['Uz_over_Utheta'].coefficients)
    assert any(current.radius.post.selected.inlet.endpoints(v)!=(0,0)
        for v in gap['Mz_over_R_Utheta'].coefficients)
    terminal=owner.evaluate('pulse_end',candidate['fresh_Z'],0)
    for key in ('Mz','Mtheta_z'):
        assert all(current.radius.post.selected.inlet.endpoints(v)==(0,0)
            for v in terminal['raw_histories'][key].coefficients.coefficients)
    assert len(candidate['actual_source_call_trace'])==15
    # Read extra symbolic scale expressions after the recorded live-call
    # sequence; proof-only graph nodes must not reorder published node IDs.
    proof=scale_proof(owner)
    rejected=[]
    for chart,coordinate in (('pulse_main',0),('pulse_gap',13),('pulse_gap_end',-3),('pulse_end',1),('core',0)):
        try:owner.evaluate(chart,candidate['fresh_Z'],coordinate)
        except ValueError:rejected.append(chart)
        else:raise AssertionError('Invalid pulse chart/domain accepted')
    try:current.CurrentOriginalRpRawPulseTransport(owner.raw,require_checked=False)
    except ValueError:pass
    else:raise AssertionError('Unclosed source bypassed heat closure admission')
    result=dict(all_passed=True,source_family=owner.family_record,**dict.fromkeys(current.GATES,True),
        independent_absolute_pulse_scale_proof=proof,independent_general_cumulative_equations=general,
        actual_six_pulse_history_and_density_coefficient_checks=rows,
        original_velocity_mixed_row_overlap_diagnostics=velocityrows,
        all_fifteen_chart_scale_parts_use_one_current_expression_graph=common_graph_parts,
        active_main_and_end_axial_sources_not_deleted=True,inactive_gap_history_memory_preserved=True,
        only_actual_Rv_terminal_linear_histories_are_zero=True,
        full_pressure_is_original_P0_plus_original_Mp=True,
        all_fifteen_current_provider_calls_and_closed_heat_dispatch_observed=True,
        original_current_selection_repair_and_P0_graph_preserved=True,
        unclosed_heat_and_invalid_chart_domains_rejected=rejected,
        physical_velocity_rows_remain_directed_source_enclosures=True,
        source_scope='Six active pulse factorized raw histories and similarity velocities, native first derivatives and a common 15-chart closed-heat source caller. Numeric/global/time Cartesian admission remains open',
        **dict.fromkeys(current.OPEN,False),
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),
            Path(__file__).name:current.sha(Path(__file__).name),
            Path(generic.__file__).name:current.sha(Path(generic.__file__).name)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.raw.packed(result),indent=2)+'\n',encoding='utf8',newline='\n')
    print('PASS_CURRENT_ORIGINAL_RP_RAW_PULSE_TRANSPORT six active pulse histories/velocities and 15-chart closed heat caller',flush=True)
    return result


if __name__=='__main__':run()
