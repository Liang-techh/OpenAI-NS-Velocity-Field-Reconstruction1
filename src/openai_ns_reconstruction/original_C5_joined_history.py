"""Contract the O2/O3 signed IBP seam only after a same-source certificate.

The accepted 197-node history graph remains immutable. References below name
its nested expression namespace explicitly. Neither signed values nor old
positive error caps are evaluated or subtracted.
"""
from __future__ import annotations

from pathlib import Path

from .original_finite_N_signed_source import require,sha,digest,load
from .original_signed_density_expression import SignedDensityExpressionGraph,ScalarExpression

PARENT='lei_ren_part1_paper_compliant_current_original_radial_IBP.json'
PARENT_SHA256='67945b4e9a8f2b1cdaa9a255b42a3bd9e5377f514cb12409ba41e9f17f957af0'
GRAPH_PATH=['actual_C5_signed_history_IBP','exact_local_expression_nodes']


def contracted_history(H0,M1,M2,C_in,C_out,I1,I2,inverse_N):
    """The two rate-one charts with their common internal C endpoint removed."""
    return M2*M1*H0+inverse_N*(C_out-M2*M1*C_in-M2*I1-I2)


def exact_controls():
    import sympy as s
    z=s.symbols('Z');M1,M2,inv=s.symbols('M1 M2 inverse_N')
    H,Cin,Cross,Cout,I1,I2=(s.Function(name)(z) for name in ('H0','Cin','Cross','Cout','I1','I2'))
    first=M1*H+inv*(Cross-M1*Cin-I1)
    nested=M2*first+inv*(Cout-M2*Cross-I2)
    result=contracted_history(H,M1,M2,Cin,Cout,I1,I2,inv)
    row_checks=[s.expand(s.diff(nested-result,z,k))==0 for k in range(4)]
    X1,X2=s.symbols('Cross1 Cross2')
    unequal=M2*(M1*H+inv*(X1-M1*Cin-I1))+inv*(Cout-M2*X2-I2)
    residual=s.expand(unequal-result)
    checks=dict(same_source_internal_endpoint_cancels=all(row_checks),
        unequal_endpoints_leave_exact_signed_difference=s.expand(residual-inv*M2*(X1-X2))==0,
        incoming_history_coefficient=s.expand(s.diff(result,H)-M1*M2)==0,
        both_slow_integrals_remain=s.diff(result,I1)==-inv*M2 and s.diff(result,I2)==-inv,
        only_one_inverse_N_in_joined_bracket=s.diff(result,inv,2)==0)
    require(all(checks.values()),'New signed-history seam algebra failed')
    return dict(passed=True,checks=checks,ordinary_Z0_Z3= row_checks,
        unequal_seam_residual='M2*(Cbuffer11-Ctransition0)/N',
        source_values_or_error_caps_used=False)


def _require_seam(seam):
    require(seam.get('passed') is True and seam.get('source_only') is True and
        seam.get('accepted_parent_sha256')==PARENT_SHA256,'Passed same-parent seam certificate required')
    require(seam.get('ordinary_Z_orders')==[0,1,2,3] and
        seam.get('arbitrary_issued_phase_B_identity') is True and
        seam.get('anchored_primitive_identity') is True and
        seam.get('common_unwrapped_phase_identity') is True,
        'An actual-phase sample or unmatched primitive cannot remove the seam')


def build(root,parent,seam):
    root=Path(root);directory=root/'experiments/root_st073'
    require(sha(directory/PARENT)==PARENT_SHA256,'Accepted source history changed')
    require(digest(parent)==digest(load(directory/PARENT)),'Supplied parent is not the accepted report')
    _require_seam(seam)
    accepted=parent['actual_C5_signed_history_IBP'];nodes=accepted['exact_local_expression_nodes']
    require(accepted['passed'] is True and accepted['local_nodes']==197 and
        accepted['numeric_value_available'] is False,'Actual accepted signed-history functions required')
    ports=accepted['source_ports'];b,t=(ports[n] for n in ('O2_buffer','O3_transition'))
    require([r['incoming_H'] for r in t['rows']]==b['H_out_endpoint_rows'],
        'Transition does not retain the actual buffer endpoint')
    require(len(b['rows'])==len(t['rows'])==4,'All four ordinary-Z rows required')
    g=SignedDensityExpressionGraph()
    def ref(index):
        require(type(index) is int and 0<=index<len(nodes),'Missing accepted expression')
        operation=nodes[index]['operation']
        require(not any(s in operation.lower() for s in ('cap','enclosure','majorant','bound_only')),
            'A bound cannot become a signed history function')
        return g.insert(dict(operation='accepted_signed_expression_reference',report=PARENT,
            report_sha256=PARENT_SHA256,graph_path=GRAPH_PATH,expression_node=index,
            expression_operation=operation,reference_is_function_not_numeric_enclosure=True))
    def at(value,chart,endpoint):
        return g.insert(dict(operation='function_substitution',expression=value.node,
            variable=ref(chart['coordinate']).node,value=ref(endpoint).node,
            substitution_through_reference=True,capture_avoiding=True,
            referenced_graph_path=GRAPH_PATH,source_native_substitution=True))
    M1=at(ref(b['rows'][0]['memory']),b,b['domain_upper_node'])
    M2=ref(t['rows'][0]['memory']);inverse_N=ref(t['inverse_N'])
    rows=[];moving=[];endpoint=[]
    for k,(br,tr) in enumerate(zip(b['rows'],t['rows'])):
        require(br['ordinary_Z_order']==tr['ordinary_Z_order']==k,'Wrong ordinary derivative order')
        require(br['memory']==b['rows'][0]['memory'] and tr['memory']==t['rows'][0]['memory'],
            'Source memory must be Z-independent in all rows')
        H0=ref(br['incoming_H']);Cin=ref(br['C_in']);Cout=ref(tr['C_out'])
        I1=at(ref(br['IBP_integral']),b,b['domain_upper_node']);I2=ref(tr['IBP_integral'])
        joined=contracted_history(H0,M1,M2,Cin,Cout,I1,I2,inverse_N)
        end=at(joined,t,t['domain_upper_node'])
        moving.append(joined.node);endpoint.append(end.node)
        rows.append(dict(ordinary_Z_order=k,H0=H0.node,C_in=Cin.node,C_out=Cout.node,
            buffer_slow_integral_at11=I1.node,transition_slow_integral=I2.node,
            original_buffer_C_at11=at(ref(br['C_out']),b,b['domain_upper_node']).node,
            original_transition_C_at0=ref(tr['C_in']).node,
            original_nested_moving_H=ref(tr['H_out']).node,
            original_nested_endpoint_H=ref(t['H_out_endpoint_rows'][k]).node,
            contracted_moving_H=joined.node,contracted_endpoint_H=end.node))
    controls=exact_controls()
    # This is a new application of the accepted pure pressure kernel, with
    # the newly joined source functions as inputs. Its old controls/builders
    # and its coefficient interval evaluation are not replayed.
    from .original_O3_inlet_functional import quiet_Nscaled_sectors
    from .original_finite_N_density_kernel import ordinary_jet
    inlet=parent['O3_Nscaled_inlet_functional'];inlet_nodes=inlet['exact_local_expression_nodes']
    require(inlet['passed'] is True and inlet['physical_background_Q0']['ordinary_Z1_Z2_exact_zero'] is True,
        'Accepted constant physical background Q0 required')
    require(accepted['H']['source_archive_outgoing_anchor_Z0_to_Z3']==
        inlet['normalized_history_source']['ordinary_m_rows'] and
        inlet['normalized_history_source']['N']==18,
        'Joined transition endpoint is not the actual quiet-functional inlet')
    def inlet_ref(index):
        return g.insert(dict(operation='accepted_signed_expression_reference',report=PARENT,
            report_sha256=PARENT_SHA256,
            graph_path=['O3_Nscaled_inlet_functional','exact_local_expression_nodes'],
            expression_node=index,expression_operation=inlet_nodes[index]['operation'],
            reference_is_function_not_numeric_enclosure=True))
    require(inlet_nodes[3]['operation']=='original_function_reference' and
        inlet_nodes[3]['original_node']==250215 and inlet_nodes[3]['original_operation']=='bound_variable' and
        inlet_nodes[3]['report']==inlet['normalized_history_source']['report'] and
        inlet_nodes[3]['report_sha256']==inlet['normalized_history_source']['report_sha256'] and
        inlet_nodes[12]['original_node']==20 and inlet_nodes[18]['original_node']==25220 and
        inlet['physical_background_Q0']['physical_Q0_local_node']==16 and
        inlet_nodes[17]['operation']=='bounded_expm1' and inlet_nodes[17]['argument']==3,
        'Actual quiet functional ports changed')
    quiet=quiet_Nscaled_sectors([ScalarExpression(g,n) for n in endpoint],
        ordinary_jet(inlet_ref(12),g.one,g.zero),inlet_ref(18),
        ordinary_jet(inlet_ref(16),g.zero,g.zero),inlet_ref(3),inlet_ref(17))
    quiet_ports=dict(H_input_rows=endpoint,qhat_C0_C1_C2=[quiet['qhat'].value.node,
        quiet['qhat'].dz.node,quiet['qhat'].dzz.node],
        linear_C0_C1=[quiet['linear'].value.node,quiet['linear'].dz.node],
        quadratic_C0_C1=[quiet['quadratic'].value.node,quiet['quadratic'].dz.node],
        recombination='R0^2*((Pstar/N)*linear+(Pstar/N)^2*quadratic)',
        H_already_N_scaled=True,physical_Q0='-4*mstar*exp(-12)',
        ordinary_Z1_and_H_Z3_retained=True,new_source_kernel_application=True,
        signed_H_or_pressure_increment_numerically_evaluated=False)
    return dict(passed=True,source_only=True,accepted_parent=PARENT,accepted_parent_sha256=PARENT_SHA256,
        accepted_seam_certificate_sha256=digest(seam),accepted_graph_path=GRAPH_PATH,
        ordinary_Z_orders=[0,1,2,3],M1=M1.node,M2=M2.node,inverse_N=inverse_N.node,
        rows=rows,H_moving_Z0_Z3=moving,H_endpoint_Z0_Z3=endpoint,O3_quiet_functional=quiet_ports,
        formula='M2*M1*H0+(Cout-M2*M1*Cin-M2*I1-I2)/N',
        only_same_source_internal_C_pair_contracted=True,
        original_graph_and_all_error_views_unchanged=True,
        predecessor_history_and_outer_endpoints_retained=True,
        both_signed_slow_integrals_retained=True,exact_controls=controls,
        numeric_signed_H_values_installed=False,physical_N_or_scale_materialized=False,
        old_positive_caps_subtracted=False,global_five_moment_closure=False,
        accepted_owner_fixture_recurrence_or_quadrature_replayed=False,**g.packet())
