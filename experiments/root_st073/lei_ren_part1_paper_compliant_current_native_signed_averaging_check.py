"""Independent split/IBP references and actual source averaging route checks."""
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as sy
import lei_ren_part1_paper_compliant_current_native_signed_averaging as current
import lei_ren_part1_paper_compliant_current_native_Rc_parameter_targets_check as previous

HERE,PREFIX,sha=current.HERE,current.PREFIX,current.sha
packets,ep=current.packets,current.ep
require,same=previous.require,previous.same


def exact_mixed_source_theorem():
    y,Z,w=sy.symbols('y Z psi');psi=sy.Function('psi')(y,Z)
    # Independent quadratic free-angle functions realize every partial
    # needed by the general implicit second chain rule.
    coefficients=sy.symbols('c0:10')
    polynomial=sum(c*t for c,t in zip(coefficients,(1,w,y,Z,w*w,w*y,w*Z,y*y,y*Z,Z*Z)))
    substitute=lambda f:f.subs(w,psi)
    ly,lZ=sy.diff(psi,y),sy.diff(psi,Z);lyZ=sy.diff(psi,y,Z)
    lam=sy.diff(polynomial,w)
    chain=substitute(sy.diff(polynomial,y,Z))+substitute(sy.diff(lam,y))*lZ+substitute(sy.diff(lam,Z))*ly+substitute(sy.diff(lam,w))*ly*lZ+substitute(lam)*lyZ
    require(sy.expand(sy.diff(substitute(polynomial),y,Z)-chain)==0,'Mixed original implicit inverse chain changed')
    a,b,E=(sy.Function(k)(y,Z) for k in ('a','b','E'));phi=sy.Symbol('fixed_phi')
    T=substitute(polynomial);chi=phi-psi/(2*sy.pi);A=a*chi/2;M=-a*T/(2*sy.pi)-b*phi;B=E*M/2
    mixedA=sy.diff(a,y,Z)*chi/2-(sy.diff(a,y)*lZ+sy.diff(a,Z)*ly+a*lyZ)/(4*sy.pi)
    mixedB=(sy.diff(E,y,Z)*M+sy.diff(E,y)*sy.diff(M,Z)+sy.diff(E,Z)*sy.diff(M,y)+E*sy.diff(M,y,Z))/2
    require(sy.expand(sy.diff(A,y,Z)-mixedA)==0 and sy.expand(sy.diff(B,y,Z)-mixedB)==0,
        'Explicit mixed original A/B chain changed')
    Tmixed=substitute(sy.diff(polynomial,y,Z))+substitute(sy.diff(polynomial,w,y))*lZ+substitute(sy.diff(polynomial,w,Z))*ly+substitute(sy.diff(polynomial,w,2))*ly*lZ+substitute(sy.diff(polynomial,w))*lyZ
    require(sy.expand(sy.diff(T,y,Z)-Tmixed)==0,'Composed original T1 mixed chain changed')
    return dict(passed=True,independent_general_implicit_inverse_yZ_chain=True,
        independent_original_A_B_and_T1_mixed_product_chain_identities=3,
        fixed_fractional_phi_not_total_y=True,all_partial_function_operators_hold_other_free_variables_fixed=True)


def independent_references(c):
    p=mp.mp.clone();p.dps=110;C=lambda v:current.LogUpper.constant(c,v)
    # All bounds are of modest independent analytic fixtures, not values
    # substituted for the native source functions.
    jet=lambda values:current.SlowJet(c,{k:C(v) for k,v in zip(current.ORDERS,values)})
    E,V,A,B=jet(('1.4','.2','.15','.05')),jet(('.7','.1','.1','.02')),jet(('3','1','.6','.2')),jet(('2','.4','.3','.1'))
    leading=dict(m=B,h=E*A,k=E*(V*A+B),e=current.SlowJet.constant(c,2)*V*B+E*E*A,p=E*E*A)
    bounds=current.remainder_caps(dict(E=E,V=V,A=A,B=B,leading=leading))
    def fields(y,z,phi):
        return (p.mpf('1.3')+p.mpf('.1')*y+p.mpf('.1')*z+p.mpf('.02')*y*z,
            p.mpf('-.4')+p.mpf('.05')*y-p.mpf('.05')*z+p.mpf('.01')*y*z,
            (2+p.mpf('.5')*y+p.mpf('.3')*z+p.mpf('.1')*y*z)*p.sin(2*p.pi*phi),
            (p.mpf('1.2')+p.mpf('.2')*y-p.mpf('.1')*z+p.mpf('.04')*y*z)*p.sin(4*p.pi*phi))
    def linear(y,z,phi):
        e,v,a,b=fields(y,z,phi)
        return dict(m=b,h=e*a,k=e*(v*a+b),e=2*v*b-e*e*a,p=e*e*a)
    def delta(y,z,phi,N):
        e,v,a,b=fields(y,z,phi);changed_e=e*p.exp(a/N);changed_v=v+b/N
        return dict(m=changed_v-v,h=changed_e-e,k=changed_e*changed_v-e*v,
            e=changed_v**2-changed_e**2/2-(v*v-e*e/2),p=(changed_e**2-e*e)/2)
    checked=0
    for N in (160,1024,10**20):
        for phi in ('0','.137','.51','.87','1'):
            ph=p.mpf(phi);y=p.mpf('.13');z=p.mpf('.07')
            for key in current.RATES:
                remainder=(delta(y,z,ph,N)[key]-linear(y,z,ph)[key]/N)*N*N
                remZ=p.diff(lambda zz:(delta(y,zz,ph,N)[key]-linear(y,zz,ph)[key]/N)*N*N,z)
                for value,cap in ((remainder,bounds['values'][key]),(remZ,bounds['Z_derivatives'][key])):
                    require(abs(value)<p.mpf('1e-70') if cap.log is None else p.log(abs(value)+p.mpf('1e-100'))<=p.mpf(str(ep(cap.log)[1])),
                        'Independent exact original-density remainder outside cap: '+key)
                    checked+=1
    # Finite intervals deliberately end at nonintegral phases, so endpoint
    # terms are nonzero. Quadrature subtracts the original density directly.
    ibp_count=0;endpoint_nonzero=0
    for rate in (0,1,p.mpf('1.5')):
        aa,bb=p.mpf('.013'),p.mpf('.061');N=17;z=p.mpf('.07');phase0=p.mpf('.23')
        amplitude=lambda y,Z:1+p.mpf('.3')*y+p.mpf('.2')*Z+p.mpf('.1')*y*Z
        f=lambda y,Z,phi:amplitude(y,Z)*p.sin(2*p.pi*phi)
        G=lambda y,Z,phi:amplitude(y,Z)*(1-p.cos(2*p.pi*phi))/(2*p.pi)
        phase=lambda y:N*y+phase0;K=lambda y:p.exp(-rate*(bb-y))
        for differentiated in (False,True):
            ff=lambda y:(p.diff(lambda zz:f(y,zz,phase(y)),z) if differentiated else f(y,z,phase(y)))
            gg=lambda y:(p.diff(lambda zz:G(y,zz,phase(y)),z) if differentiated else G(y,z,phase(y)))
            gy=lambda y:(p.diff(lambda zz:p.diff(lambda yy:G(yy,zz,phase(y)),y),z) if differentiated
                else p.diff(lambda yy:G(yy,z,phase(y)),y))
            endpoint=K(bb)*gg(bb)-K(aa)*gg(aa)
            left=p.quad(lambda y:K(y)*ff(y)/N,[aa,bb])
            right=(endpoint-p.quad(lambda y:K(y)*(gy(y)+rate*gg(y)),[aa,bb]))/(N*N)
            require(abs(left-right)<p.mpf('1e-95'),'Independent C0/Z signed IBP value or kernel sign failed')
            endpoint_nonzero+=abs(endpoint)>p.mpf('1e-10');ibp_count+=1
    require(endpoint_nonzero==6,'Independent fixture must detect omitted endpoint terms')
    return dict(passed=True,independent_original_density_second_remainder_C0_Z_comparisons=checked,
        independent_finite_interval_signed_IBP_C0_Z_comparisons=ibp_count,
        nonzero_endpoint_term_cases=endpoint_nonzero,N_remainder_cases=[160,1024,'10^20'],
        fixtures_do_not_define_actual_source_values=True)


@current.native.inlet.source_precision
def run(bridge=None):
    began=time.monotonic();saved=json.loads((HERE/current.NAME).read_bytes())
    require(saved[current.GATE] and saved['source_graph_chart_count']==17 and saved['all_integer_N_lower']==160,
        'Same original all-N averaged source stage required')
    for name,digest in saved['input_hashes'].items():require(sha(name)==digest,'Changed original signed averaging source: '+name)
    require(not saved['actual_control_functions_or_terminal_closure_installed'] and not any(saved.get(k) for k in packets.OPEN),
        'Averaged C0/Z covers cannot install terminal controls or higher/global admission')
    if bridge is None:bridge,_=current.native.inlet.native_bridge_owner()
    with current.native.inlet.CheckedSourceRuntime():
        owner=current.NativeSignedAveraging(current.current.NativeRcParameterTargets(current.current.preceding.NativeRcC1Histories(
            current.current.preceding.preceding.NativeO2C1Histories(current.current.preceding.preceding.make_middle_owner(bridge)))))
        require(saved['exact_source_averaging_theorem']==packets.encode(owner.theorem),'Exact original averaging theorem changed')
        mixed=exact_mixed_source_theorem();independent=independent_references(owner.ctx)
        graphs=json.loads(gzip.decompress((HERE/current.VIEWS).read_bytes()))
        require(graphs==packets.encode(owner.graphs),'Exact mixed source graph adapter changed')
        for chart,graph in graphs.items():
            require(graph['partial_derivative_nodes_have_function_semantics'] and graph['derivative_caps_not_function_values']
                and not graph['actual_point_inverse_or_mixed_jet_evaluator_installed'],
                'Mixed derivative operator definitions cannot become point values')
            require(set(graph['leading_signed_density_slow_roots'])==set(current.RATES), 'Five signed leading functions required')
            for rows in graph['leading_signed_density_slow_roots'].values():
                require(set(rows)=={'y0_Z0','y1_Z0','y0_Z1','y1_Z1'},'Four fixed-phi leading density derivatives required')
        rows=inherited=quiet=pressure=targets=0;regions={};comparisons={}
        for name,record in saved['actual_original_native_signed_averaging_records'].items():
            got=owner.route(packets.interval(owner.ctx,record['Z_box']))
            require(packets.encode(got['record'])==record,'Actual endpoint-retaining averaged route changed: '+name)
            incoming={key:owner.coordinates.scalar(0) for key in current.RATES};incomingZ=dict(incoming)
            for cell in got['cells']:
                label=cell['record']['label'];chart=cell['record']['chart']
                for key in current.RATES:
                    decay=cell['factors'][key]['decay']
                    same(cell['incoming'][key],incoming[key],'Averaging reset inherited C0 memory: '+label)
                    same(cell['incoming_Z'][key],incomingZ[key],'Averaging reset inherited Z memory: '+label)
                    same(cell['cumulative'][key],decay*incoming[key]+cell['values'][key],'Averaged C0 semigroup mismatch')
                    same(cell['cumulative_Z'][key],decay*incomingZ[key]+cell['Z_derivatives'][key],'Averaged Z semigroup mismatch')
                    for v in (cell['values'][key],cell['Z_derivatives'][key],cell['cumulative'][key],cell['cumulative_Z'][key]):
                        require(v.scale.bases is owner.coordinates.bases and v.ledger is owner.coordinates.ledger,
                            'Averaged covers lost same original common-unit source context');rows+=1
                    inherited+=2
                    if chart=='O3_power':
                        require(cell['values'][key].zero and cell['Z_derivatives'][key].zero,'Averaging created nonzero quiet source');quiet+=2
                        if key=='p':
                            same(cell['cumulative'][key],incoming[key],'Quiet averaging discarded C0 pressure memory')
                            same(cell['cumulative_Z'][key],incomingZ[key],'Quiet averaging discarded Z pressure memory');pressure+=2
                    if not cell['flat']:
                        proof=cell['record']['per_density_IBP_bounds'][key]
                        require(proof['C0']['all_each_cell_endpoint_terms_kept'] and proof['Z']['all_each_cell_endpoint_terms_kept'],
                            'Zero mean cannot discard finite-cell endpoints')
                incoming,incomingZ=cell['cumulative'],cell['cumulative_Z']
            for key in current.RATES:
                require(got['values'][key][-1].zero and got['Z_derivatives'][key][-1].zero,
                    'Original zero inlet and IBP must remove N^-1 route certificate term')
            for key in current.repair.ROWS:
                require(got['targets']['values'][key][-1].zero and got['targets']['Z_derivatives'][key][-1].zero,
                    'Normalized original averaged target retains an N^-1 bound term');targets+=4
            comparisons[name]=record['same_source_direct_and_averaged_bound_comparison']
            require(not record['tightened_first_bridge_yZ_caps_or_actual_controls_installed']
                and record['cell_endpoint_cancellation_or_selected_phase_not_assumed'],
                'Old mixed caps or exact zero means cannot establish tight controls')
            regions[name]=dict(original_true_cells=24,original_charts=17,
                all_N_C0_Z_Nminus2_source_integral_bounds=True,
                actual_controls_or_terminal_Z_closure=False)
            print('Original signed averaging checked:',name,flush=True)
    hashes=dict(owner.service.hashes);hashes[current.NAME]=sha(current.NAME);hashes[Path(__file__).name]=sha(Path(__file__).name)
    result=dict(all_passed=True,source_family=owner.family,**{current.GATE:True},native_Z_queries_checked=2,
        original_true_cells_checked=24,original_charts_checked=17,source_graph_mixed_leading_rows_checked=17*5*4,
        actual_C0_Z_contribution_and_transport_rows_checked=rows,inherited_C0_Z_memory_rows_checked=inherited,
        exact_quiet_C0_Z_rows_checked=quiet,preserved_pressure_memory_rows_checked=pressure,normalized_target_order_rows_checked=targets,
        independent_exact_mixed_source_theorem=mixed,independent_reference_checks=independent,
        same_source_direct_and_averaged_bound_comparisons=comparisons,regions=regions,
        exact_source_averaging_theorem=owner.theorem,
        actual_controls_or_terminal_closure_installed=False,
        one_global_finite_N_or_control_field_admitted=False,global_inlet_to_Rc_histories_admitted=False,
        **dict.fromkeys(packets.OPEN,False),execution_seconds=time.monotonic()-began,input_hashes=hashes,
        scope=saved['scope'])
    (HERE/current.RECEIPT).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    print('Native signed zero-mean C0/Z averaging PASS',rows,inherited,quiet,pressure,targets,flush=True)
    return result


if __name__=='__main__':run()
