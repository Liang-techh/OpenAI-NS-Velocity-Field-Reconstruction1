"""Actual O3 inlet pressure functional with Pstar/N scale sectors.

The input H=N*delta-m remains the actual signed history function. No norm
bound or zero center is admitted as H. Physical Q0 is the background inlet
pressure, whose true source is independent of Z.
"""
from __future__ import annotations

from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from .original_finite_N_signed_source import load,sha,digest,require,ast_guard
from .original_finite_N_density_kernel import ordinary_jet,OrdinaryJet
from .original_signed_density_expression import SignedDensityExpressionGraph

PARENT='lei_ren_part1_paper_compliant_current_original_Rc_pressure_join.json'
PARENT_SHA256='5f160cfa1f251253155bdad0466c462c74c5d7349bdf5cbcb9b3b7209822e9c4'
POST='lei_ren_part1_paper_compliant_current_original_post_axial_m15_numeric.json'
OUTER='lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions.json.gz'


def _exact_zero(value):
    if hasattr(value,'graph') and hasattr(value,'node'):
        node=value.graph.nodes[value.node]
        return node['operation']=='exact_rational' and node['numerator']==0
    if hasattr(value,'_mpi_'):return all(tuple(v)==(0,0,0,0) for v in value._mpi_)
    return value==0


def quiet_Nscaled_sectors(H,Z,delta,Q0,x,expm1x):
    """Return C1 sectors for R0^2*(Pstar/N*A+(Pstar/N)^2*B).

    H has precisely four ordinary Z rows. Q0 has its actual ordinary rows
    (including two exact zero derivatives); no missing row is inferred.
    """
    require(isinstance(H,(list,tuple)) and len(H)==4,'Actual H ordinary Z0..3 rows required')
    require(isinstance(Z,OrdinaryJet) and Z.dzz is not None,'Explicit Z jet required')
    require(isinstance(Q0,OrdinaryJet) and Q0.dzz is not None,'Explicit background Q0 jet required')
    require(_exact_zero(Q0.dz) and _exact_zero(Q0.dzz),'Source-proved Z-independent Q0 required')
    require(not any(isinstance(v,OrdinaryJet) for v in (delta,x,expm1x)),
        'delta and native coordinate factors must be Z-independent source scalars')
    D=1-Z*Z;L=1-Z*Z*delta
    qhat=-(Z*(1-delta)*ordinary_jet(*H[:3])+D*ordinary_jet(*H[1:]))/L
    qhatZ=ordinary_jet(qhat.dz,qhat.dzz)
    linear=Z*(1-delta)*qhatZ*expm1x/(4*L)-Q0*qhat*x/2
    quadratic=-qhat*qhat*x/4
    return dict(qhat=qhat,linear=linear,quadratic=quadratic)


def functional_coefficients(Z,delta,Q0,x,expm1x):
    """Actual linear and symmetric quadratic coefficient forms on H0..H3.

    The returned pairs are coefficients for the value and ordinary Z1 of
    the output. H is never assigned numerical basis values or set to zero.
    """
    require(_exact_zero(Q0.dz) and _exact_zero(Q0.dzz),'Source-proved Z-independent background Q0 required')
    zero=Z.value*0;L=1-Z*Z*delta
    c0=-Z*(1-delta)/L;c1=-(1-Z*Z)/L
    c=[c0.value,c1.value,zero,zero]
    dc=[c0.dz,c0.value+c1.dz,c1.value,zero]
    ddc=[c0.dzz,2*c0.dz+c1.dzz,c0.value+2*c1.dz,c1.value]
    alpha=Z*(1-delta)*expm1x/(4*L)
    A0=[alpha.value*dc[j]-Q0.value*x*c[j]/2 for j in range(4)]
    A1=[alpha.dz*dc[j]+alpha.value*ddc[j]-Q0.value*x*dc[j]/2 for j in range(4)]
    B0=[[-x*c[i]*c[j]/4 for j in range(4)] for i in range(4)]
    B1=[[-x*(dc[i]*c[j]+c[i]*dc[j])/4 for j in range(4)] for i in range(4)]
    return dict(linear=[A0,A1],quadratic=[B0,B1])


def exact_controls():
    import sympy as s
    z,de,t,lam,Q0=s.symbols('Z delta x lambda Q0',real=True)
    H=s.Function('H')(z);rows=[s.diff(H,z,k) for k in range(4)]
    result=quiet_Nscaled_sectors(rows,ordinary_jet(z,1,0),de,ordinary_jet(Q0,0,0),t,s.exp(t)-1)
    qhat=-(z*(1-de)*H+(1-z*z)*s.diff(H,z))/(1-de*z*z)
    primitive=lam*result['linear'].value+lam*lam*result['quadratic'].value
    Qb=Q0*s.exp(-t);Q=(Q0+lam*qhat)*s.exp(-t);L=1-de*z*z
    halfS=lambda q:z*(1-de)*s.diff(q,z)/(4*L)-q*q/4
    derivative=s.simplify(s.diff(primitive,t)-s.exp(2*t)*(halfS(Q)-halfS(Qb)))==0
    C1=s.simplify(lam*result['linear'].dz+lam*lam*result['quadratic'].dz-s.diff(primitive,z))==0
    inlet=s.simplify(primitive.subs(t,0))==0
    viscous=s.simplify(s.diff(-s.exp(t)*(Q-Qb),t))==0
    # The background m=4Z*M/Pstar identity is universal in delta and Z.
    M,P=s.symbols('M Pstar',nonzero=True)
    recovered=-P*(z*(1-de)*(4*z*M/P)+(1-z*z)*(4*M/P))/L
    bg=s.simplify(recovered+4*M)==0
    require(derivative and C1 and inlet and bg and viscous,'New N-scaled inlet functional control failed')
    try:quiet_Nscaled_sectors(rows[:3],ordinary_jet(z,1,0),de,ordinary_jet(Q0,0,0),t,s.exp(t)-1)
    except ValueError:missing=True
    else:missing=False
    require(missing,'Missing ordinary derivative was padded')
    try:quiet_Nscaled_sectors(rows,ordinary_jet(z,1,0),de,ordinary_jet(Q0,1,0),t,s.exp(t)-1)
    except ValueError:dependent=True
    else:dependent=False
    require(dependent,'Z-dependent background pressure was admitted as a constant')
    symbols=s.symbols('H0:4');Qjet=ordinary_jet(Q0,0,0);zjet=ordinary_jet(z,1,0)
    coeff=functional_coefficients(zjet,de,Qjet,t,s.exp(t)-1)
    direct=quiet_Nscaled_sectors(list(symbols),zjet,de,Qjet,t,s.exp(t)-1)
    cchecks=[]
    for k in range(2):
        linear=sum(a*h for a,h in zip(coeff['linear'][k],symbols))
        quadratic=sum(coeff['quadratic'][k][i][j]*symbols[i]*symbols[j] for i in range(4) for j in range(4))
        cchecks.extend((s.simplify(linear-[direct['linear'].value,direct['linear'].dz][k])==0,
            s.simplify(quadratic-[direct['quadratic'].value,direct['quadratic'].dz][k])==0))
    require(all(cchecks),'Coefficient tensors do not recover the exact signed functional')
    return dict(passed=True,full_reduced_S_derivative=derivative,ordinary_C1=C1,zero_inlet_increment=inlet,
        universal_background_pressure_identity=bg,missing_H_Z3_rejected=missing,
        source_R_delta_Qy_endpoint_cancellation=viscous,Z_dependent_Q0_rejected=dependent,
        linear_quadratic_coefficient_reconstruction=all(cchecks),coefficient_identity_count=len(cchecks),
        original_owners_fixtures_or_quadratures_executed=False)


def _exact(iv,record):
    return iv.mpf([mp.make_mpf(tuple(v)) for v in record['source_enclosure_exact_mpf_tuples']])


def _record(value):
    return dict(source_enclosure_exact_mpf_tuples=[list(v) for v in value._mpi_],
        signed_enclosure_only=True,midpoint_selected_as_value=False)


def coefficient_enclosures(nodes,Q0,iv):
    """New bounded coefficient evaluation on two entire rational Z cells."""
    memo={}
    def scalar(index):
        if index in memo:return memo[index]
        node=nodes[index];op=node['operation']
        if op=='exact_rational':value=iv.mpf(node['numerator'])/node['denominator']
        elif op=='negative':value=-scalar(node['argument'])
        elif op=='sum':value=sum((scalar(i) for i in node['arguments']),iv.mpf(0))
        elif op=='product':
            value=iv.mpf(1)
            for i in node['arguments']:value*=scalar(i)
        elif op=='analytic_unary' and node['name']=='exp':value=iv.exp(scalar(node['argument']))
        else:raise ValueError('Unexpected scalar operation in actual delta source')
        memo[index]=value;return value
    delta=scalar(25220);require(delta.a>0 and delta.b<iv.mpf('1e-200'),'Changed original small positive delta')
    results=[]
    for a,b in ((-3,0),(0,3)):
        Z=iv.mpf([a,b])/10;z=ordinary_jet(Z,iv.mpf(1),iv.mpf(0))
        require((1-Z*Z*delta).a>0,'Positive original meridional denominator required')
        coeff=functional_coefficients(z,delta,ordinary_jet(Q0,iv.mpf(0),iv.mpf(0)),iv.mpf(2),iv.expm1(2))
        results.append(dict(exact_Z_cell=[[a,10],[b,10]],native_endpoint=2,
            linear_C0_C1=[[_record(v) for v in row] for row in coeff['linear']],
            quadratic_C0_C1=[[[_record(v) for v in row] for row in matrix] for matrix in coeff['quadratic']],
            full_cell_directed_enclosures=True,actual_H_values_evaluated=False))
    return dict(passed=True,precision=iv.dps,actual_delta_source_node=25220,delta_enclosure=_record(delta),
        evaluated_scalar_source_nodes=sorted(memo),cells=results,
        quantity='operator coefficients on H0..H3, not H or pressure increment values',
        basis_H_values_not_selected=True,full_Z_domain_or_terminal_closure_claim=False)


def build(root,parent):
    """Bind actual H ports and recover the known signed background coefficient."""
    root=Path(root);work=root/'experiments/root_st073'
    require(sha(work/PARENT)==PARENT_SHA256 and digest(load(work/PARENT))==digest(parent),'Accepted Rc parent required')
    for name in (POST,OUTER):require(sha(work/name)==parent['input_hashes'][name],'Changed actual inlet source')
    require(parent['quiet_O3_signed_primitive']['passed'] is True and
        parent['Rc_profile_join']['whole_O3_power_quiet'] is True,'Accepted exact quiet chart source required')
    post=load(work/POST);outer=load(work/OUTER);ob=outer['true_outer_C6_C5_bindings']
    require(post['source_family']==outer['source_family']==parent['source_family'],'Different construction data')
    axial='src/openai_ns_reconstruction/original_axial_pressure_rows.py'
    post_source='src/openai_ns_reconstruction/original_post_axial_pressure_rows.py'
    loop_source='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py'
    for rel in (axial,post_source):require(sha(root/rel)==parent['implementation_hashes'][rel],'Changed analytic background source')
    require(sha(root/loop_source)==parent['input_hashes'][Path(loop_source).name],'Changed original finite N density source')
    guards={
        'exact_mean_recovery':ast_guard(root,axial,'fields',(
            'U=4*z*beta;M=4*z*mean;Uy=4*z*beta_y',
            'Q=(z*z*(8*beta)/A-4*mean).truncate(2)')),
        'same_global_mean':ast_guard(root,post_source,'scalars',(
            "mean=self.mean_truth*self.iv.exp(-offset-t)",)),
        'actual_m_density_units':ast_guard(root,loop_source,'finite_N_functions',(
            "invN=g.quotient(g.one,N,'same exact positive integer N on every source chart')",
            "s=a.scale(loop['B_over_Pstar'],invN)","inc=dict(m=s,h=r")),
    }
    # All three saved axial diagnostics name the same scalar enclosure. The
    # analytic recovery proves its Z independence; agreement is not the proof.
    samples=[v for v in post['pressure_row_samples'] if v['chart']=='O3_power' and v['native_coordinate']=='0']
    require(len(samples)==3 and len({v['Z'] for v in samples})==3,'Three saved inlet diagnostics required')
    base=samples[0];mean_source=base['positive_original_axial_mean_source'];mean_inlet=base['leading_mean']
    require(all(v['positive_original_axial_mean_source']==mean_source and v['leading_mean']==mean_inlet for v in samples),
        'Actual background scalar differs across the saved sources')
    accepted_mean=parent['quiet_scalar_source']
    require(accepted_mean['passed'] is True and accepted_mean['source_mean_equality']['passed'] is True and
        accepted_mean['mean_source_summary']['summary']==mean_source and
        ob['seed_identity']['propagated_original_histories_identified_by_linear_ODE_uniqueness'] is True,
        'Same exact global mean source/transport identity required')
    iv=MPIntervalContext();iv.dps=160
    Q0=-4*_exact(iv,mean_inlet)
    require(Q0.b<0,'The actual signed background pressure must be negative')
    for sample in samples:
        qrows=sample['physical_pressure_rows']['Q']['signed_background_ordinary_Z']
        require(len(qrows)==3,'Three background Q ordinary rows required')
        require(_exact(iv,qrows[0])._mpi_==Q0._mpi_ and
            all(_exact(iv,v)._mpi_==iv.mpf(0)._mpi_ for v in qrows[1:]),'Saved background pressure does not match the exact source identity')
    nodes=outer['exact_graph_nodes'];g=SignedDensityExpressionGraph()
    coefficients=coefficient_enclosures(nodes,Q0,iv)
    ref=lambda n:g.reference(OUTER,parent['input_hashes'][OUTER],n,source_nodes=nodes)
    chart=ob['charts']['O3_power'];N=ref(ob['N']);x=ref(chart['coordinate'])
    H=[N*ref(n) for n in chart['incoming']['m'][:4]]
    z=ordinary_jet(ref(20),g.one,g.zero)
    # Use the actual physical background pressure at native zero. The signed
    # scalar enclosure stays separate and never replaces this function node.
    background=ob['background']['O3_power']
    Q0ref=g.insert(dict(operation='function_substitution',
        expression=(ref(33572)*ref(background['generic']['Q'][0])).node,
        variable=x.node,value=g.zero.node,endpoint_independent_of_Z=True,
        theorem='same actual background source is -4*mstar*exp(-12), independent of Z'))
    E=g.insert(dict(operation='bounded_expm1',argument=x.node,source_domain='0<=x<=2',exact_function=True))
    result=quiet_Nscaled_sectors(H,z,ref(25220),ordinary_jet(Q0ref,g.zero,g.zero),x,E)
    pair=lambda q:[q.value.node,q.dz.node]
    graph_packet=g.packet();graph_packet.pop('original_sources_numerically_executed')
    return dict(passed=True,source_executed=False,parent_report=PARENT,parent_report_sha256=PARENT_SHA256,AST=guards,
        source_family=parent['source_family'],normalized_history_source=dict(report=OUTER,
            report_sha256=parent['input_hashes'][OUTER],ordinary_m_rows=chart['incoming']['m'][:4],N=ob['N'],
            new_H_rows=[v.node for v in H],definition='H=N*delta-m_in; actual signed source, no numerical value selected'),
        physical_background_Q0=dict(exact_function='-4*mstar*exp(-12)',ordinary_Z1_Z2_exact_zero=True,
            actual_normalized_Q_source=background['generic']['Q'][:3],physical_Q0_local_node=Q0ref.node,
            inherited_common_mean_identity_sha256=digest(accepted_mean['source_mean_equality']),
            inlet_mean_source_enclosure=mean_inlet,mstar_source_enclosure=mean_source,
            signed_Q0_enclosure=_record(Q0),strictly_negative=True,
            proof='Universal rational mean recovery plus actual global rate-one background source',
            independent_diagnostic_count=3,diagnostic_agreement_is_not_function_identity_proof=True),
        Nscaled_sector_roots={key:pair(result[key]) for key in ('linear','quadratic')},
        actual_directed_coefficient_enclosures=coefficients,
        qhat_C0_C1_C2=[result['qhat'].value.node,result['qhat'].dz.node,result['qhat'].dzz.node],
        physical_recombination='R0^2*((Pstar/N)*linear+(Pstar/N)^2*quadratic)',
        physical_linear_sector='Z*(1-delta)*qhat_Z*expm1(x)/(4L)-Q0*qhat*x/2',
        physical_quadratic_sector='-qhat^2*x/4',
        N_Pstar_and_R0_not_materialized=True,unknown_H_not_set_to_zero=True,
        N_scaled_history_values_or_signed_centers_evaluated=False,
        finite_N_error_logs_not_used_as_H_values=True,
        existing_Rc178_views_and_all_error_terms_unchanged=True,
        archived_source_owners_numerically_executed=False,actual_delta_source_interval_evaluated=True,
        referenced_history_functions_numerically_evaluated=False,
        source_execution_scope='source_executed=False means no archived owner/producer; the bounded delta parameter graph is interval evaluated',
        **graph_packet)
