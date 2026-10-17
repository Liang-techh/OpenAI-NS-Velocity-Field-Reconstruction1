"""Complete same-source pressure-boundary congruence at Rc.

This binds exact endpoint inputs before constructing the common algebraic
boundary. It does not identify independent enclosures or select a cap value.
"""
from __future__ import annotations

from pathlib import Path

from .original_finite_N_signed_source import load,sha,digest,require
from .original_compact_inlet_identity import _Reducer
from .original_inner_collar_source import _ast_guard
from .original_finite_N_density_kernel import ordinary_jet
from .original_scaled_signed_density import correlated_scaled_density
from .original_signed_density_expression import SignedDensityExpressionGraph

PARENT='lei_ren_part1_paper_compliant_current_original_correlated_scaled_source.json'
PARENT_SHA256='9beead66814da85325e86f4ee2d9c091da34c524f9512f8cbf4c77e18bf1d444'
PREFIX='lei_ren_part1_paper_compliant_current_original_'
OUTER=PREFIX+'true_outer_C6_C5_functions.json.gz'
COMPACT=PREFIX+'true_compact_C5_functions.json.gz'


def quiet_pressure_sectors(q0,dq,Z,delta,length,expm1_length):
    """Exact R0^2*(Pstar*linear+Pstar^2*quadratic) for a quiet chart.

    q0/dq are normalized inlet ordinary jets through Z2. The caller proves
    U=0 and Q_y=-Q from the actual source. Input length is independent of Z.
    """
    require(q0.dzz is not None and dq.dzz is not None,'Inlet Q_ZZ is required for ordinary C1')
    L=1-Z*Z*delta
    dqZ=ordinary_jet(dq.dz,dq.dzz)
    return dict(linear=Z*(1-delta)*dqZ/(4*L)*expm1_length,
        quadratic=-(2*q0*dq+dq*dq)/4*length)


def quiet_exact_controls():
    """New algebraic integration, no original datum or integral execution."""
    import sympy as s
    from fractions import Fraction as F
    z,de,x,T=s.symbols('Z delta x T',real=True)
    q0=1+z+z*z;dq=2-z+z**3
    row=lambda q:ordinary_jet(q,s.diff(q,z),s.diff(q,z,2))
    result=quiet_pressure_sectors(row(q0),row(dq),ordinary_jet(z,1,0),de,T,s.exp(T)-1)
    L=1-de*z*z
    # The original source domain has L>0. Integrate only x here so SymPy
    # does not add a spurious L=0 branch outside that admitted domain.
    direct_linear=z*(1-de)*s.diff(dq,z)/(4*L)*s.integrate(s.exp(x),(x,0,T))
    direct_quadratic=s.integrate(-((q0+dq)**2-q0**2)/4,(x,0,T))
    checks=dict(linear_integral=s.simplify(result['linear'].value-direct_linear)==0,
        quadratic_integral=s.simplify(result['quadratic'].value-direct_quadratic)==0,
        ordinary_C1_linear=s.simplify(result['linear'].dz-s.diff(direct_linear,z))==0,
        ordinary_C1_quadratic=s.simplify(result['quadratic'].dz-s.diff(direct_quadratic,z))==0)
    P,R0=s.symbols('Pstar R0',positive=True)
    Q=P*(q0+dq)*s.exp(-x);Qb=P*q0*s.exp(-x);R=R0*s.exp(x)
    source_half_S=lambda q:((q+z*(1-de)*s.diff(q,z)/2-q)/L+q*(q/2-q))/2
    source_difference=s.expand(R*R*(source_half_S(Q)-source_half_S(Qb)))
    primitive=R0*R0*(P*result['linear'].value+P*P*result['quadratic'].value)
    checks['actual_reduced_S_integrand_matches_primitive_derivative']=s.simplify(
        source_difference-s.diff(primitive,T).subs(T,x))==0
    checks['viscous_R_Qy_endpoint_difference_zero']=s.simplify(s.diff(-R*(Q-Qb),x))==0
    v=quiet_pressure_sectors(ordinary_jet(F(1),0,0),ordinary_jet(F(1,10),0,0),
        ordinary_jet(F(0),1,0),F(1,5),F(2),F(6))
    checks['zero_linear_does_not_zero_quadratic']=v['linear'].value==0 and v['quadratic'].value==F(-21,200)
    require(all(checks.values()),'New quiet exact integration control failed')
    return dict(passed=True,checks=checks,original_source_functions_or_quadratures_executed=False,
        synthetic_polynomial_controls_only=True)


def build_quiet_primitive(root,parent,profile,join):
    """Replace the quiet O3 source quadrature by an exact signed primitive."""
    root=Path(root);work=root/'experiments/root_st073'
    require(profile['passed'] is True and profile['whole_O3_power_quiet'] is True and
        profile['parent_report_sha256']==PARENT_SHA256,'Actual whole-chart quiet branch proof required')
    require(join['passed'] is True and join['required_profile_proof_digest']==digest(profile) and
        join['AST']['physical_recovery']['passed'] is True,'Same actual pressure recovery proof required')
    for name in (OUTER,COMPACT):require(sha(work/name)==parent['input_hashes'][name],'Changed quiet inlet source')
    raw=load(work/OUTER);ob=raw['true_outer_C6_C5_bindings'];chart=ob['charts']['O3_power'];bg=ob['background']['O3_power']
    nodes=raw['exact_graph_nodes'];g=SignedDensityExpressionGraph()
    from .original_Rc_profile_join import _GraphProof,_zeros
    proof=_GraphProof(nodes,'outer')
    _zeros(proof,chart['finite_N']['five_signed_density_increments']['m'])
    _zeros(proof,bg['fields']['V'])
    local_rows=[]
    for label,rows in (('correction',chart['contributions']['m']),('background',bg['history_contributions']['m'])):
        for k,index in enumerate(rows[:4]):
            row=nodes[index]
            require(row['operation']=='definite_integral' and row['original_rate']=='1' and
                row['original_chart']=='O3_power' and row['lower']==0 and row['upper']==chart['coordinate'] and
                row['physical_Jacobian_applied_once'] is True,'Changed actual mean transport integral')
            _zeros(proof,[row['integrand']])
            local_rows.append(dict(passed=True,sector=label,ordinary_Z_order=k,node=index,
                integrand=row['integrand'],integrand_identically_zero=True))
    require(chart['memory']['m']==bg['memory']['m'] and chart['physical_measure']==1,
        'Changed common quiet mean memory or measure')
    outer_source='experiments/root_st073/'+PREFIX+'true_outer_C6_C5_functions.py'
    guards={
        'mean_transport':_ast_guard(root/outer_source,root,'transport_geometry',text=(
            'width=g.sub(offset,origin);suffix=g.sub(offset,at_t)',
            'decay=g.unary(\'exp\',g.mul(g.constant(-rate),width))',
            'out[key]=a.add(a.scale(incoming[key],decay),local[key])')),
        'quiet_geometry':_ast_guard(root/outer_source,root,'geometry',text=(
            "'O3_power':g.add(c(18),Y,x)","else g.one","'O3_power':(g.zero,c(2))")),
    }
    for key in ('reduced_S','compact_FTC'):
        receipt=parent['finite_N_signed_source']['AST'][key]
        require(receipt['passed'] is True and receipt['source_executed'] is False and
            sha(root/receipt['path'])==parent['input_hashes'][Path(receipt['path']).name],
            'Accepted actual reduced S/physical m15 source changed')
        guards[key]=dict(receipt,inherited_saved_only=True)
    ref=lambda n:g.reference(OUTER,parent['input_hashes'][OUTER],n,source_nodes=nodes)
    z=ordinary_jet(ref(20),g.one,g.zero);de=ref(25220);D=1-z*z;L=1-z*z*de
    m=[ref(n) for n in bg['incoming']['m'][:4]];dm=[ref(n) for n in chart['incoming']['m'][:4]]
    pressure_from_m=lambda values:(-z*(1-de)*ordinary_jet(*values[:3])-D*ordinary_jet(*values[1:4]))/L
    q0,dq=pressure_from_m(m),pressure_from_m(dm)
    x=ref(chart['coordinate']);exp_minus_one=g.insert(dict(operation='bounded_expm1',argument=x.node,
        source_domain='0<=x<=2',exact_function=True,physical_scale_materialized=False))
    sectors=quiet_pressure_sectors(q0,dq,z,de,x,exp_minus_one)
    # R0 is an exact source substitution, not a numerical radius.
    r0=g.insert(dict(operation='function_substitution',expression=ref(chart['radius']).node,
        variable=x.node,value=g.zero.node,endpoint_independent_of_Z=True,source_endpoint='O3_power native0'))
    P=ref(33572);physical=(sectors['linear']*P+sectors['quadratic']*(P*P))*(r0*r0)
    lower_F=-(dq/(2*L)*P+(2*q0*dq+dq*dq)/4*(P*P))*(r0*r0)
    pair=lambda q:[q.value.node,q.dz.node]
    return dict(passed=True,source_executed=False,required_profile_proof_digest=digest(profile),
        report=OUTER,report_sha256=parent['input_hashes'][OUTER],native_domain=[0,2],AST=guards,
        actual_zero_mean_integrands=local_rows,common_mean_decay_node=chart['memory']['m'],
        normalized_quiet_increment_reduction='S/2=Z*(1-delta)*Q_Z/(4L)-Q^2/4',
        normalized_inlet_background_m=bg['incoming']['m'][:4],normalized_inlet_correction_m=chart['incoming']['m'][:4],
        quiet_source_equations=['U=U_y=0','m_y=-m','Q_y=-Q','Q=Pstar*q_in(Z)*exp(-x)'],
        normalized_sector_C0_C1={key:pair(value) for key,value in sectors.items()},
        physical_increment_C0_C1=pair(physical),negative_remaining_lower_R0_squared_deltaF_C0_C1=pair(lower_F),
        inlet_R0=r0.node,**g.packet(),
        exact_physical_primitive='R0^2*(Pstar*Z*(1-delta)*deltaq_in_Z*expm1(x)/(4L)-Pstar^2*(2q0_in*deltaq_in+deltaq_in^2)*x/4)',
        all_quadratic_and_inlet_history_terms_retained=True,viscous_endpoint_difference_exact_zero=True,
        physical_volume_measure='R^2 dx',incoming_pressure_not_reset=True,
        original_signed_inlet_values_numerically_evaluated=False,
        exact_source_quadrature_eliminated=True,existing_aggregate_caps_retained=True)


def _substitution(nodes,index,body,variable,value):
    row=nodes[index]
    require(row.get('operation')=='function_substitution' and row.get('expression')==body and
        row.get('variable')==variable and row.get('value')==value,'Changed actual endpoint substitution')
    return dict(passed=True,node=index,expression=body,variable=variable,value=value)


def build(root,parent,profile):
    root=Path(root);work=root/'experiments/root_st073'
    require(sha(work/PARENT)==PARENT_SHA256 and digest(load(work/PARENT))==digest(parent),'Accepted source parent required')
    reports={}
    for name in (OUTER,COMPACT):
        require(sha(work/name)==parent['input_hashes'][name],'Changed actual Rc source graph')
        reports[name]=load(work/name)
    ob=reports[OUTER]['true_outer_C6_C5_bindings'];cb=reports[COMPACT]['true_compact_C5_bindings']
    nodes=reports[COMPACT]['exact_graph_nodes'];prefix=reports[OUTER]['exact_graph_nodes']
    require(nodes[:len(prefix)]==prefix,'Different source prefix')
    X=260181
    require(X>=len(prefix) and nodes[X].get('name')=='true_compact_repair_x','Wrong compact variable')
    inherited=parent['compact_inlet_source_identity']
    require(inherited['passed'] is True and inherited['unresolved_dependencies']==[] and
        inherited['source_graph']['same_actual_source_digest']==digest(prefix),'Accepted 39 source identities missing')
    # The separately proved profile theorem is a required premise. No endpoint
    # may be silently assigned the common zero without this source proof.
    require(profile['passed'] is True and profile['parent_report_sha256']==PARENT_SHA256,
        'Actual Rc profile continuity unproved or different source parent')
    require(profile['full_V_and_V_y_exact_zero'] is True and profile['baseline_V_and_V_y_exact_zero'] is True,
        'Actual velocity/radial source inputs are not zero at the Rc join')
    require(profile['input_source_prefix_sha256']==digest(prefix) and
        all(v['passed'] is True for v in profile['AST'].values()),'Different or unguarded profile source proof')
    ib=cb['inlet'];oe=ob['endpoints']['O3_power'];bg=ob['background_endpoints']['O3_power'];chart=ob['charts']['O3_power']
    require(cb['Rc']==oe['radius'],'Different physical Rc anchors')
    reducer=_Reducer(nodes,len(prefix));background=[];endpoint_maps=[];full_mean=[]
    for k in range(4):
        before=len(reducer.unknown);value=reducer.rewrite(ib['background_fields']['m'][k])
        expected=('atom',bg['fields']['m'][k])
        require(value==expected and len(reducer.unknown)==before,'Background mean inlet differs')
        background.append(dict(passed=True,ordinary_Z_order=k,compact_node=ib['background_fields']['m'][k],
            outer_node=bg['fields']['m'][k],same_actual_source_atom=True))
        for label,b,mean,correction,rawbg,variable,at in (
            ('outer',chart,oe['modified_fields']['m'][k],oe['correction']['m'][k],
                ob['background']['O3_power']['fields']['m'][k],chart['coordinate'],chart['domain'][1]),
            ('compact',cb,ib['modified_fields']['m'][k],ib['correction']['m'][k],
                cb['background']['fields']['m'][k],X,1)):
            fields=b['modified_fields'];history=b['history']
            row=nodes[fields['m'][k]]
            require(row.get('operation')=='sum' and sorted(row.get('arguments',()))==sorted([rawbg,history['m'][k]]),
                'Full mean is not actual background plus signed history')
            endpoint_maps.append(_substitution(nodes,mean,fields['m'][k],variable,at))
            endpoint_maps.append(_substitution(nodes,correction,history['m'][k],variable,at))
        full_mean.append(dict(passed=True,ordinary_Z_order=k,outer_full_mean=oe['modified_fields']['m'][k],
            compact_full_mean=ib['modified_fields']['m'][k],
            common_background=bg['fields']['m'][k],common_signed_correction=oe['correction']['m'][k],
            theorem='substitution distributes over the verified source sum; both summands have the same source'))
    radial='experiments/root_st073/'+PREFIX+'modified_Q4_pressure_C3_functions.py'
    outer_source='experiments/root_st073/'+PREFIX+'true_outer_C6_C5_functions.py'
    compact_source='experiments/root_st073/'+PREFIX+'true_compact_C5_functions.py'
    for rel in (radial,outer_source,compact_source):
        require(sha(root/rel)==parent['input_hashes'][Path(rel).name],'Changed endpoint recovery source')
    guards={
        'physical_recovery':_ast_guard(root/radial,root,'physical_pressure_fields',text=(
            "U=a.scale(fields['V'],Pstar);Uy=a.scale(radial['V_y'],Pstar)",
            "m=a.scale(fields['m'],Pstar);my=a.sub(U,m)",
            'Q=div(numerator(U,m));Qy=div(numerator(Uy,my))')),
        'outer_physical_call':_ast_guard(root/outer_source,root,'__init__',cls='CurrentOriginalTrueOuterC6C5Functions',
            text=('physical=radial.physical_pressure_fields(g,modified,total,self.delta,self.Pstar)',
                'boundary,pressure_local,cumulative,pressure=pressure_geometry(g,S,physical,R,x,domain[0],jac,integrals,self.initial,name)')),
        'compact_physical_call':_ast_guard(root/compact_source,root,'__init__',cls='CurrentOriginalTrueCompactC5Functions',
            text=('self.N,self.mu,self.Pstar,self.P0,self.delta=outer.N,outer.mu,outer.Pstar,outer.P0,outer.delta',
                'self.physical=radial.physical_pressure_fields(g,self.modified_fields,self.radial,self.delta,self.Pstar)',
                "self.boundary,self.pressure_local,self.integrals,self.pressure=previous.pressure_geometry(g,self.S,self.physical,R,")),
    }
    for name in ('normalized_baseline_Q','half_S1_coefficients'):
        receipt=parent['finite_N_signed_source']['AST'][name];rel=receipt['path']
        expected=(parent['implementation_hashes'].get(rel) if rel.startswith('src/') else
            parent['input_hashes'].get(Path(rel).name))
        require(receipt['passed'] is True and sha(root/rel)==expected,'Changed accepted baseline/F definition')
        guards[name]=dict(receipt,inherited_saved_only=True)
    require(all(v['passed'] is True for v in guards.values()),'Actual pressure endpoint call/recovery unbound')
    from .original_Rc_pressure_ports import build_pressure_ports
    ports=build_pressure_ports(dict(outer=prefix,compact=nodes),ob,cb,
        dict(passed=True,source_only=True,background_mean_identities=background,full_mean_identities=full_mean),profile)
    require(ports['passed'] is True and all(ports['checks'].values()),'Actual pressure port normalization failed')
    g=SignedDensityExpressionGraph();ref=lambda n:g.reference(COMPACT,parent['input_hashes'][COMPACT],n,source_nodes=nodes)
    m=[ref(n) for n in bg['fields']['m'][:4]];dm=[ref(n) for n in oe['correction']['m'][:4]]
    z=ordinary_jet(ref(20),g.one,g.zero);de=ref(cb['delta']);D=1-z*z;L=1-z*z*de
    mjet=ordinary_jet(*m[:3]);mz=ordinary_jet(*m[1:4])
    q0=(-z*(1-de)*mjet-D*mz)/L
    common=correlated_scaled_density([g.zero]*3,[q0.value,q0.dz,q0.dzz],
        [g.zero]*3,dm,[g.zero]*2,z,de)
    R=ref(cb['Rc']);P=ref(cb['Pstar']);F=common['F1']*P+common['F2']*(P*P)
    Qy=common['delta_q_y']*P;B=F*(R*R)-Qy*R
    # Endpoint congruence is justified above by actual source input identities.
    # The common B is not zeroed; only its oppositely signed occurrence cancels.
    gap=B-B
    require(gap.value.node==g.zero.node and gap.dz.node==g.zero.node,'Common source boundary pair did not contract')
    pair=lambda q:[q.value.node,q.dz.node]
    return dict(passed=True,source_executed=False,parent_report=PARENT,parent_report_sha256=PARENT_SHA256,
        input_source_prefix_sha256=digest(prefix),inherited_39_identity_receipt_digest=digest(inherited),
        required_profile_proof_digest=digest(profile),background_mean_identities=background,
        actual_physical_pressure_ports=ports,
        full_mean_identities=full_mean,actual_endpoint_maps=endpoint_maps,AST=guards,
        endpoint_interfaces=dict(outer=dict(physical=oe['physical'],boundary=oe['boundary'],radius=oe['radius']),
            compact=dict(generic=ib['generic'],radial=ib['radial'],boundary=ib['boundary'],radius=ib['radius'])),
        deterministic_recovery_congruence=dict(passed=True,mean_Z0_through_Z3=True,velocity_and_radial_same_source=True,
            Q_Qy_ordinary_C1_equal=True,F_ordinary_C1_equal=True,radius_same_source=True,
            meridional_domain='original common positive-L C5 source domain; no extension'),
        normalized_roots={key:pair(common[key]) for key in ('delta_q','delta_q_y','F1','F2')},
        common_boundary_C0_C1=pair(B),oppositely_signed_pair_C0_C1=pair(gap),**g.packet(),
        same_source_Rc_boundary_pair_exact_zero=True,common_boundary_individually_zero_claim=False,
        common_boundary_numeric_value_computed=False,all_aggregate_caps_retained=True,
        original_archives_or_discrepancy_views_changed=False,global_terminal_closure=False)
