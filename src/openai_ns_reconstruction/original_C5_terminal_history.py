"""Apply a proved outgoing flat-loop identity to the terminal signed history."""
from __future__ import annotations

from pathlib import Path

from .original_C5_directed_primitive import PARENT,PARENT_SHA256
from .original_finite_N_signed_source import require,sha,digest,load
from .original_signed_density_expression import SignedDensityExpressionGraph


def terminal_history(H0,M1,M2,Cin,I1,I2,inverse_N):
    return M2*M1*H0-inverse_N*(M2*M1*Cin+M2*I1+I2)


def controls():
    import sympy as s
    z=s.symbols('Z');M1,M2,inv=s.symbols('M1 M2 inverse_N')
    H,Cin,I1,I2=(s.Function(n)(z) for n in ('H0','Cin','I1','I2'))
    before=M2*M1*H+inv*(0-M2*M1*Cin-M2*I1-I2)
    after=terminal_history(H,M1,M2,Cin,I1,I2,inv)
    rows=[s.expand(s.diff(before-after,z,k))==0 for k in range(4)]
    checks=dict(all_four_ordinary_rows_equal=all(rows),
        predecessor_retained=s.diff(after,H)==M1*M2,
        initial_primitive_retained=s.diff(after,Cin)==-M1*M2*inv,
        both_signed_slow_integrals_retained=s.diff(after,I1)==-M2*inv and s.diff(after,I2)==-inv)
    require(all(checks.values()),'Terminal signed-history contraction failed')
    return dict(passed=True,checks=checks,ordinary_Z0_Z3=rows,error_cap_values_used=False)


def build(root,parent,flat):
    root=Path(root);directory=root/'experiments/root_st073'
    require(sha(directory/PARENT)==PARENT_SHA256 and digest(parent)==digest(load(directory/PARENT)),
        'Accepted joined-history source required')
    require(flat.get('passed') is True and flat.get('outgoing_C_Z0_Z3_exact_zero') is True and
        flat.get('accepted_parent_sha256')==PARENT_SHA256,
        'Actual terminal flat-loop C identity required')
    joined=parent['joined_C5_signed_history'];ibp=parent['actual_C5_signed_history_IBP']
    require(joined['passed'] is True and joined['numeric_signed_H_values_installed'] is False,
        'Source-bound joined history required')
    g=SignedDensityExpressionGraph()
    paths={'joined':['joined_C5_signed_history','exact_local_expression_nodes'],
           'ibp':['actual_C5_signed_history_IBP','exact_local_expression_nodes'],
           'quiet':['O3_Nscaled_inlet_functional','exact_local_expression_nodes']}
    def ref(kind,index):
        nodes=parent
        for key in paths[kind]:nodes=nodes[key]
        require(type(index) is int and 0<=index<len(nodes),'Missing accepted function expression')
        return g.insert(dict(operation='accepted_signed_expression_reference',report=PARENT,
            report_sha256=PARENT_SHA256,graph_path=paths[kind],expression_node=index,
            expression_operation=nodes[index]['operation'],reference_is_function_not_numeric_enclosure=True))
    transition=ibp['source_ports']['O3_transition']
    variable=ref('ibp',transition['coordinate']);upper=ref('ibp',transition['domain_upper_node'])
    def at(value):
        return g.insert(dict(operation='function_substitution',expression=value.node,
            variable=variable.node,value=upper.node,substitution_through_reference=True,
            capture_avoiding=True,referenced_graph_path=paths['ibp'],source_native_substitution=True))
    M1=ref('joined',joined['M1']);M2=at(ref('joined',joined['M2']))
    inverse_N=ref('joined',joined['inverse_N']);rows=[];endpoint=[]
    for row in joined['rows']:
        H0=ref('joined',row['H0']);Cin=ref('joined',row['C_in'])
        I1=ref('joined',row['buffer_slow_integral_at11'])
        I2=at(ref('joined',row['transition_slow_integral']))
        H=terminal_history(H0,M1,M2,Cin,I1,I2,inverse_N)
        endpoint.append(H.node)
        rows.append(dict(ordinary_Z_order=row['ordinary_Z_order'],H0=H0.node,Cin=Cin.node,
            I1=I1.node,I2_endpoint=I2.node,terminal_H=H.node,
            original_joined_endpoint_H=ref('joined',row['contracted_endpoint_H']).node,
            removed_outgoing_C_endpoint='exact_zero_by_actual_flat_loop_certificate'))
    # New source application only; no accepted coefficient computations are replayed.
    from .original_O3_inlet_functional import quiet_Nscaled_sectors
    from .original_finite_N_density_kernel import ordinary_jet
    inlet=parent['O3_Nscaled_inlet_functional'];ns=inlet['exact_local_expression_nodes']
    require(inlet['physical_background_Q0']['ordinary_Z1_Z2_exact_zero'] is True and
        joined['O3_quiet_functional']['H_input_rows']==joined['H_endpoint_Z0_Z3'],
        'Original quiet pressure source/input correspondence changed')
    require(ns[12]['original_node']==20 and ns[18]['original_node']==25220 and
        ns[3]['original_node']==250215 and ns[17]['operation']=='bounded_expm1' and
        ns[17]['argument']==3 and inlet['physical_background_Q0']['physical_Q0_local_node']==16,
        'Actual quiet kernel source ports changed')
    H=[type(g.one)(g,n) for n in endpoint]
    quiet=quiet_Nscaled_sectors(H,ordinary_jet(ref('quiet',12),g.one,g.zero),ref('quiet',18),
        ordinary_jet(ref('quiet',16),g.zero,g.zero),ref('quiet',3),ref('quiet',17))
    quiet_ports=dict(H_input_rows=endpoint,linear_C0_C1=[quiet['linear'].value.node,quiet['linear'].dz.node],
        quadratic_C0_C1=[quiet['quadratic'].value.node,quiet['quadratic'].dz.node],
        recombination='R0^2*((Pstar/N)*linear+(Pstar/N)^2*quadratic)',
        ordinary_Z1_and_H_Z3_retained=True,signed_H_or_pressure_increment_numerically_evaluated=False)
    return dict(passed=True,source_only=True,parent_report=PARENT,parent_sha256=PARENT_SHA256,
        actual_flat_certificate_sha256=digest(flat),formula='M2*M1*H0-(M2*M1*Cin+M2*I1+I2_at1)/N',
        ordinary_Z_orders=[0,1,2,3],M1=M1.node,M2_endpoint=M2.node,inverse_N=inverse_N.node,
        rows=rows,H_terminal_Z0_Z3=endpoint,O3_quiet_functional=quiet_ports,
        contraction_only_at_actual_transition_endpoint1=True,moving_history_unchanged=True,
        predecessor_initial_primitive_and_both_slow_integrals_retained=True,
        old_positive_caps_subtracted=False,old_graphs_or_error_views_modified=False,
        numeric_signed_H_values_installed=False,actual_phase_chosen_or_materialized=False,
        accepted_ancestor_or_quadrature_replayed=False,exact_controls=controls(),**g.packet())
