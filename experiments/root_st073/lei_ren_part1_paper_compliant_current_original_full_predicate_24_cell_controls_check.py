"""Actual original24 affine memory, joint target and finite repair diagnostics."""
import gzip
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_full_predicate_24_cell_controls as current

ep=current.ep;iv=current.packets.interval;KEYS=current.KEYS;controls=current.controls


def same(record,value):
    scale=record['formal_positive_scale'];c=value.ctx
    assert tuple(scale['source_exponents'])+(scale['radius_power'],)==value.scale.powers
    assert ep(iv(c,scale['additional_log_interval']))==ep(value.scale.offset)
    assert ep(iv(c,record['coefficient_interval']))==ep(value.coefficient)
    assert record['exact_zero']==value.zero and not record['point_value_selected'] and record['encloses_original_source_function']


def overlaps(a,b):
    assert a.ctx is b.ctx and a.scale.bases is b.scale.bases and a.ledger is b.ledger
    anchor=a.ctx.mpf(max(ep(a.scale.evaluate())[1],ep(b.scale.evaluate())[1]))
    left=a.coefficient*a.bounded_exp(a.scale.evaluate()-anchor)
    right=b.coefficient*b.bounded_exp(b.scale.evaluate()-anchor)
    return max(ep(left)[0],ep(right)[0])<=min(ep(left)[1],ep(right)[1])


def cap(c,values,jets):
    return {key:current.parameters.repair.LogUpper.add(c,[current.parameters.magnitude(values[key]),
        current.parameters.magnitude(jets[key])]).record() for key in values}


def first_difference(a,b,path=''):
    if type(a) is not type(b):return path,type(a).__name__,type(b).__name__
    if isinstance(a,dict):
        if a.keys()!=b.keys():return path,list(a),list(b)
        for key in a:
            got=first_difference(a[key],b[key],path+'/'+str(key))
            if got:return got
    elif isinstance(a,list):
        if len(a)!=len(b):return path,len(a),len(b)
        for i,(x,y) in enumerate(zip(a,b,strict=True)):
            got=first_difference(x,y,path+'/'+str(i))
            if got:return got
    elif a!=b:return path,str(a)[:160],str(b)[:160]


def source_chain(owner,live,saved):
    assert current.encode(live['record'])==saved
    assert saved['actual_source_prefix_slope_downstream_counts']==[13,5,6]
    c=owner.ctx;coords=owner.coordinates;incoming=[{key:coords.scalar(0) for key in KEYS} for unused in range(2)]
    manifests={};affine=0;background=0;geometry=0;local=0;oldoverlaps=0;quiet=0;budgets=[]
    composite=current.history.C1DuhamelOperator(coords)
    for index,(cell,spec) in enumerate(zip(live['cells'],current.parameters.ROUTE,strict=True)):
        row=cell['record'];label,chart,left,right=spec;binding=row['actual_original_source_cell_binding']
        assert (row['label'],row['chart'])==(label,chart) and row['candidate_N']==1024
        assert row['source_family']==owner.family and ep(row['Z_box'])==(-1,1)
        name=binding['manifest'];assert current.sha(name)==binding['manifest_sha256']
        if name not in manifests:
            data=(current.HERE/name).read_bytes();manifests[name]=json.loads(gzip.decompress(data) if name.endswith('.gz') else data)
        raw=manifests[name]
        for part in binding['record_path']:raw=raw[part]
        assert raw['source_family']==owner.family and raw['candidate_N']==1024 and ep(iv(c,raw['Z_box']))==(-1,1)
        fresh=owner.target_owner.geometry(label,chart,left,right)
        assert current.encode(fresh['record'])==row['actual_original_geometry'];geometry+=1
        for key,rate in current.history.RATES.items():
            decay=current.preceding.downstream.transfer.true_width_kernel(coords,fresh,rate)['decay']
            assert overlaps(decay,cell['operator'].coefficients[key])
        op=cell['operator']
        for key in KEYS:
            composite.increments[key]=op.coefficients[key]*composite.increments[key]+op.increments[key]
            composite.Z_increments[key]=op.coefficients[key]*composite.Z_increments[key]+op.Z_increments[key]
            composite.coefficients[key]=op.coefficients[key]*composite.coefficients[key]
        composite.steps+=1
        for j,(field,suffix) in enumerate((('values','C0'),('Z_derivatives','Z'))):
            increments=cell[field]
            for key in KEYS:
                same(row['actual_inherited_correction_'+suffix][key],incoming[j][key])
                same(row['actual_local_'+suffix+'_contributions'][key],increments[key])
                assert overlaps(increments[key],(op.increments if j==0 else op.Z_increments)[key]);local+=1
                # Independent affine formula; integrated source increments are not multiplied by a second mass.
                expected=incoming[j][key]*op.coefficients[key]+increments[key]
                same(row['actual_right_correction_'+suffix][key],cell['correction'][field][key])
                assert overlaps(expected,cell['correction'][field][key]);affine+=1
                bg=cell['background'][field][key];same(row['actual_right_background_'+suffix][key],bg)
                same(row['actual_right_own_history_'+suffix][key],bg+cell['correction'][field][key]);background+=1
                assert overlaps(cell['correction'][field][key],owner.restore(cell['saved_ancestor_out'][j][key]));oldoverlaps+=1
                if index>=22:
                    assert increments[key].zero and overlaps(incoming[j][key],cell['correction'][field][key]);quiet+=1
            p0=cell['background']['P0' if j==0 else 'P0_Z']
            same(row['original_separate_P0' if j==0 else 'original_separate_P0_Z'],p0)
            same(row['absolute_pressure_'+suffix],p0+(cell['background'][field]['p']+cell['correction'][field]['p']))
            incoming[j]=cell['correction'][field]
        budgets.append(dict(label=label,local_C1_caps=cap(c,cell['values'],cell['Z_derivatives']),
            outgoing_correction_C1_caps=cap(c,cell['correction']['values'],cell['correction']['Z_derivatives'])))
        if index==12:
            for j in range(2):
                for key in KEYS:assert overlaps(incoming[j][key],owner.restore(owner.upstream_saved[j][key]))
    composed=composite.apply({key:coords.scalar(0) for key in KEYS},{key:coords.scalar(0) for key in KEYS},owner.family)
    for j,field in enumerate(('values','Z_derivatives')):
        for key in KEYS:
            same(saved['actual_Rc_correction_C0_Z'][j][key],incoming[j][key])
            assert overlaps(composed[field][key],incoming[j][key])
    assert not incoming[0]['p'].zero and not incoming[1]['p'].zero
    assert live['amplitude']['record']['original_Rc_power_offset']==2
    assert saved['spatial_O2_global_alternative_not_summed_into_direct24_history']
    return dict(passed=True,actual_original24_source_manifest_geometry_bindings=geometry,
        independent_original_own_rate_decays=120,integrated_local_C0_Z_increment_comparisons=local,
        independent_actual_serial_C0_Z_affine_memory_comparisons=affine,
        actual_background_once_C0_Z_own_history_checks=background,
        inherited_accepted_source_output_C0_Z_overlaps=oldoverlaps,quiet_power_zero_local_nonzero_inherited_checks=quiet,
        independent_serial_composite_C0_Z_overlaps=10,
        source_cell_C1_error_budget_inventory=budgets,
        original_P0_separate_actual_pressure_memory_retained=True,ancestor_source_suites_not_replayed=True)


def exact_function_binding(owner,source,result):
    transport=controls.source;built=result['live']['built'];g=built['graph'];binding_owner=result['role_owner'].owner
    maps,symbols,x=transport.exact_radius_maps();symbolic_route=transport.exact_route_endpoints(current.parameters.ROUTE,symbols)
    import sympy as sy
    views=binding_owner.views;viewhash=built['source_graph_sha256']
    assert current.sha(transport.sources.VIEWS)==viewhash
    parameter=built['parameters'];N=built['N'];source_bindings=[];geometry=0;original_nodes=len(g.nodes)
    for index,(cell,spec,typed) in enumerate(zip(built['cells'],symbolic_route,source['cells'],strict=True)):
        label,chart,left,right=spec;t=g.symbol('coordinate_'+str(index))
        make=lambda expr:transport.expression(g,expr,parameter,t)
        lo,hi=make(left),make(right);offset=make(maps[chart]);yleft=make(sy.simplify(maps[chart].subs(x,left)))
        yright=make(sy.simplify(maps[chart].subs(x,right)))
        width=make(sy.simplify(maps[chart].subs(x,right)-maps[chart].subs(x,left)))
        jac=make(sy.diff(maps[chart],x))
        phase=None if cell['source_flat_exact_zero'] else g.unary('fractional_part',g.mul(N,offset))
        assert (cell['lower'],cell['upper'],cell['left_radius_offset'],cell['right_radius_offset'],cell['width'])==(
            lo.node,hi.node,yleft.node,yright.node,width.node);geometry+=1
        assert cell['all_incoming_histories_retained'] and cell['original_P0_and_P0_Z_unchanged']
        if index:
            assert cell['incoming']==built['cells'][index-1]['outgoing']
        for key,rate in transport.RATES.items():
            pair=cell['contributions'][key]
            if cell['source_flat_exact_zero']:
                assert pair==dict(value=g.zero.node,Z=g.zero.node)
                assert typed['values'][key].zero and typed['Z_derivatives'][key].zero
                continue
            kernel=g.unary('exp',g.neg(g.mul(g.constant(rate),g.sub(yright,offset))))
            for order,rootkey,role in (('value','five_signed_increment_rate_roots','density_'+key+'_C0'),
                ('Z','five_signed_increment_rate_first_derivatives','density_'+key+'_Z')):
                integral=g.nodes[pair[order]];assert integral['operation']=='definite_integral'
                assert (integral['lower'],integral['upper'],integral['variable'])==(lo.node,hi.node,'coordinate_'+str(index))
                assert integral['measure']=='native coordinate; original dy/dcoordinate applied exactly once'
                args=g.nodes[integral['integrand']]['arguments']
                refs=[i for i in args if g.nodes[i]['operation']=='original_function_graph'];assert len(refs)==1
                node=refs[0];ref=g.nodes[node]
                expectedroot=views[chart][rootkey][key] if order=='value' else views[chart][rootkey]['Z'][key]
                assert ref['source_node']==expectedroot and ref['chart']==chart
                assert ref['graph_file']==transport.sources.VIEWS and ref['graph_sha256']==viewhash
                assert (ref['coordinate'],ref['phase'],ref['shared_N'],ref['Z_variable'])==(t.node,phase.node,N.node,'Z')
                assert ref['source_coefficients_are_function_recipes_not_cover_values']
                assert built['source_roles'][node]==('function_graph_nodes',role)
                assert (ref['source_graph_namespace'],ref['function_role'])==('function_graph_nodes',role)
                expected=g.mul(kernel,transport.FunctionRef(g,node),jac)
                assert expected.node==integral['integrand']
                source_bindings.append(dict(label=label,rate=key,order=order,source_node=expectedroot,
                    exact_graph_contribution_node=pair[order],source_graph_sha256=viewhash,
                    accepted_actual_numerical_source_binding=typed['record']['actual_original_source_cell_binding']))
    assert len(source_bindings)==210 and len(built['source_roles'])==212
    endpoint=make(2/symbols['Tw']);rc_offset=make(maps['O3_power'].subs(x,2/symbols['Tw']))
    rcphase=g.unary('fractional_part',g.mul(N,rc_offset))
    roots=views['O3_power']['original_signed_input_graph']['jet_expression_dag']['roots']['E']
    for q,key,role in ((built['amplitude'].value,'y0_Z0','Rc_E_C0'),(built['amplitude'].Z,'y0_Z1','Rc_E_Z')):
        ref=g.nodes[q.node]
        assert ref['source_node']==roots[key] and ref['graph_sha256']==viewhash
        assert (ref['coordinate'],ref['phase'],ref['shared_N'])==(endpoint.node,rcphase.node,N.node)
        assert built['source_roles'][q.node]==('original_signed_input_graph.jet_expression_dag',role)
    assert len(g.nodes)==original_nodes
    return dict(passed=True,original24_exact_endpoint_radius_width_bindings=geometry,
        original_C0_Z_density_graph_root_phase_integrand_Jacobian_bindings=source_bindings,
        original_Rc_amplitude_C0_Z_graph_bindings=2,source_graph_sha256=viewhash,
        exact_global_phase_is_fractional_part_N_times_original_radius_offset=True,
        Jacobian_applied_once_quiet_flat_cells_zero_and_all_incoming_graph_functions_retained=True)


def equivalent_density_definition(result):
    """Bind equivalent numerical density recipes to the original view AST.

    The archived numerical providers do not use the view's node numbering.
    Compare their defining signed moment differences and ordinary-Z product
    rules symbolically; no source function is replaced by a cover endpoint.
    """
    import sympy as s
    z=s.Symbol('Z');N=s.Symbol('N',positive=True)
    E,V,A,B=[s.Function(key)(z) for key in ('E','V','A','B')]
    EN=E*s.exp(A/N);VN=V+B/N
    expected=dict(m=VN-V,h=EN-E,k=EN*VN-E*V,
        e=VN**2-V**2-(EN**2-E**2)/2,p=(EN**2-E**2)/2)
    comparisons=0
    for chart,view in result['role_owner'].owner.views.items():
        nodes=view['function_graph_nodes'];roots=view['roots'];inputroots=view['original_signed_input_graph']['jet_expression_dag']['roots']
        replacements={inputroots['E']['y0_Z0']:E,inputroots['E']['y0_Z1']:s.diff(E,z),
            roots['A']:A,roots['A_Z_slow']:s.diff(A,z),roots['B_over_Pstar']:B,roots['B_Z_slow']:s.diff(B,z)}
        for i,node in enumerate(nodes):
            if node['operation']=='source_derivative' and node.get('name') in ('V_y0_Z0','V_y0_Z1'):
                replacements[i]=V if node['name']=='V_y0_Z0' else s.diff(V,z)
            if node['operation']=='shared_positive_integer_parameter' and node.get('name')=='common_N':replacements[i]=N
        cache={}
        def walk(i):
            if i in replacements:return replacements[i]
            if i in cache:return cache[i]
            row=nodes[i];op=row['operation']
            if op=='constant':value=s.Integer(row['value'])
            elif op=='sum':value=s.Add(*(walk(j) for j in row['arguments']))
            elif op=='product':value=s.Mul(*(walk(j) for j in row['arguments']))
            elif op=='negative':value=-walk(row['argument'])
            elif op=='positive_function_quotient':value=walk(row['numerator'])/walk(row['denominator'])
            elif op=='analytic_unary' and row['name'] in ('exp','expm1'):
                value=s.exp(walk(row['argument']))-(1 if row['name']=='expm1' else 0)
            else:raise AssertionError('Unexpected unbound density recipe node: '+str((chart,i,op)))
            cache[i]=value;return value
        for key,value in expected.items():
            assert s.simplify(s.expand(walk(view['five_signed_increment_rate_roots'][key])-value))==0
            assert s.simplify(s.expand(walk(view['five_signed_increment_rate_first_derivatives']['Z'][key])-s.diff(value,z)))==0
            comparisons+=2
        assert view['source_family']==result['role_owner'].family
        phase=view['common_phase_binding']
        assert phase['phi']=='fractional_part(N*log(R/r_minus))' and phase['total_Z_phase_derivative_exact_zero']
    return dict(passed=True,independent_original_changed_moment_difference_and_Z_AST_identities=comparisons,
        original_charts=17,exact_density_formula='deltaE=E*expm1(A/N),deltaV=B_over_Pstar/N; original signed five moment differences and complete first-Z product rules',
        numerical_archives_use_source_rooted_equivalent_recipes_not_generic_view_node_numbers=True,
        source_graph_sha256=result['live']['built']['source_graph_sha256'],
        no_field_value_or_derivative_selected_from_exported_covers=True)


def joint_target(owner,live,saved):
    H,J=live['history'],live['Z_derivatives'];a=live['amplitude'];target=live['target'];c=owner.ctx
    A,AZ,mu=a['A'],a['AZ'],a['mu_source'];logA,logmu=a['logA'],a['logmu']
    assert ep(a['mu'])==ep(owner.target_owner.repair_mu) and ep(a['mu'])[0]>0
    C=H['k']-A*H['m'];CZ=J['k']-AZ*H['m']-A*J['m']
    same(saved['actual_Rc_joint_numerator_C0_Z'][0],C);same(saved['actual_Rc_joint_numerator_C0_Z'][1],CZ)
    expected={};expectedZ={}
    for out,key,degree in (('M','m',1),('I','h',1),('S','e',2),('Cp','p',2)):
        den=A if degree==1 else A*A
        expected[out]=H[key].positive_divide(den,degree*logA)
        expectedZ[out]=(J[key]*A-degree*H[key]*AZ).positive_divide(den*A,(degree+1)*logA)
    divided=controls.ROWS[1]
    expected[divided]=C.positive_divide(A*A*mu,2*logA+logmu)
    expectedZ[divided]=(CZ*A-2*C*AZ).positive_divide(A*A*A*mu,3*logA+logmu)
    for j,field in enumerate(('values','Z_derivatives')):
        for key in controls.ROWS:
            same(saved['actual_Rc_joint_target_C0_Z'][j][key],target[field][key])
            assert overlaps((expected if j==0 else expectedZ)[key],target[field][key])
    return dict(passed=True,actual_new_Rc_joint_numerator_C0_Z_comparisons=2,
        independently_rearranged_new_Rc_quotient_C0_Z_overlaps=10,
        actual_new_N_scaled_target_C1_caps=cap(c,{key:target['values'][key]*1024 for key in controls.ROWS},
            {key:target['Z_derivatives'][key]*1024 for key in controls.ROWS}),
        same_actual_positive_mu_and_fresh_original_Rc_amplitude=True,old_target_or_P0_not_substituted=True)


def independent_Q(c,target,matrix,row):
    a,b,e,f,k=row;C,E,P=[matrix[key] for key in ('cross_weights','energy_weights','pressure_weights')]
    zero=a.value.scalar(0);scalar=zero.scalar;mu=scalar(target['mu'])
    div=lambda q:q.positive_divide(mu,target['logmu'])
    # Group derivatives by controls rather than copying range_quadratic's product ordering.
    values=[zero,div(k.value*b.value*C[1]+e.value*a.value*C[0]),zero,
        (b.value*b.value-k.value*k.value*c.mpf('.5'))*E[2]+(a.value*a.value-e.value*e.value*c.mpf('.5'))*E[0]-f.value*f.value*E[1]*c.mpf('.5'),
        sum((q.value*q.value*weight for weight,q in zip(P,(e,f,k))),zero)*c.mpf('.5')]
    jets=[zero,div((k.value*b.Z+k.Z*b.value)*C[1]+(e.value*a.Z+e.Z*a.value)*C[0]),zero,
        (b.Z*b.value*2-k.Z*k.value)*E[2]+(a.Z*a.value*2-e.Z*e.value)*E[0]-f.Z*f.value*E[1],
        sum((q.Z*q.value*weight for weight,q in zip(P,(e,f,k))),zero)]
    return [controls.C1Enclosure(v,j) for v,j in zip(values,jets,strict=True)]


def finite_repair(owner,source,result,saved):
    # JSON object keys are strings; source-role node IDs are integer keys live.
    serialized=json.loads(json.dumps(current.encode(result['report'])))
    assert serialized==saved,first_difference(serialized,saved)
    live=result['live'];c=owner.ctx;matrix=live['exact_matrix_enclosures'];B=matrix['linear_enclosure'];inv=matrix['inverse_enclosure']
    assert result['range_owner'].target is result['control_owner'].target is owner.target_owner
    built=live['built'];assert built['source_family']==owner.family
    assert [(row['label'],row['chart']) for row in built['cells']]==[(row[0],row[1]) for row in current.parameters.ROUTE]
    assert len(live['ranges'])==saved['finite_picard_depth']+1==4
    assert all(q.value.zero and q.Z.zero for q in live['ranges'][0])
    # Independent block inverse action must include the exact identity matrix.
    for i in range(5):
        for j in range(5):
            box=sum((B[i][k]*inv[k][j] for k in range(5)),c.mpf(0));v=int(i==j)
            assert ep(box)[0]<=v<=ep(box)[1]
    residuals=0;updates=0;quadratics=0
    d=live['N_scaled_target_ranges'];zero=d[0].value.scalar(0)
    for depth,(row,record) in enumerate(zip(live['ranges'],saved['actual_finite_control_and_residual_ranges'],strict=True)):
        Q=independent_Q(c,source['target'],matrix,row)
        productionQ=controls.range_quadratic(c,source['target']['mu'],source['target']['logmu'],matrix,row)
        for j,label in enumerate(controls.ROWS):
            for field,suffix in (('value','C0'),('Z','Z')):
                assert overlaps(getattr(Q[j],field),getattr(productionQ[j],field));quadratics+=1
                # B acts in the residual; the inverse acts only in Picard updates.
                value=getattr(d[j],field)+getattr(Q[j],field)*(c.mpf(1)/1024)
                value+=sum((getattr(row[k],field)*B[j][k] for k in reversed(range(5))),zero)
                assert overlaps(value,owner.restore(record['actual_repair_equation_residual_C0_Z_ranges'][label][suffix]));residuals+=1
            if depth<len(live['ranges'])-1:
                for field in ('value','Z'):
                    rhs=[getattr(d[k],field)+getattr(Q[k],field)*(c.mpf(1)/1024) for k in range(5)]
                    nxt=-sum((rhs[k]*inv[j][k] for k in reversed(range(5))),zero)
                    assert overlaps(nxt,getattr(live['ranges'][depth+1][j],field));updates+=1
    diagnostic=controls.fixed_N_contraction_diagnostic(live)
    assert current.encode(diagnostic)==saved['fixed_N_sufficient_contraction_diagnostic']
    assert diagnostic['half_contraction_sufficient_bound_passed'] is False
    assert diagnostic['failed_sufficient_bound_does_not_prove_actual_map_diverges']
    assert not diagnostic['globally_compatible_N_admitted']
    assert not saved['certified_fixed_point_tail_installed']
    return dict(passed=True,exact_original24_integral_graph_order_and_same_live_owner=True,
        independent_matrix_inverse_identity_enclosures=25,independent_quadratic_C0_Z_overlaps=quadratics,
        independent_finite_Picard_C0_Z_update_overlaps=updates,actual_numerical_linear_B_not_inverse_residual_C0_Z_overlaps=residuals,
        finite_depth=3,sufficient_half_contraction_bound_passed=False,
        failed_outer_bound_not_actual_divergence_or_solved_controls=True)


def guards(owner,source,result):
    count=0;control=result['control_owner'];range_owner=result['range_owner']
    wrongsource=dict(source,record=dict(source['record'],source_family=dict(source['record']['source_family'],source='foreign')))
    missing=dict(source,cells=source['cells'][:-1]);unclosed=dict(source,record=dict(source['record'],full24_original_C1_integral_range_transport_enclosed=False))
    attempts=(lambda:owner.assemble(N=True),lambda:owner.assemble(N=2048),
        lambda:owner.finite_controls(dict(source)),lambda:owner.finite_controls(None),
        lambda:owner.finite_controls(source,iterations=0),lambda:owner.finite_controls(source,iterations=17),
        lambda:control.controls((-1,1),True,range_owner=range_owner,source_ranges=source),
        lambda:control.controls((-1,1),2048,range_owner=range_owner,source_ranges=source),
        lambda:control.controls((-1,0),1024,range_owner=range_owner,source_ranges=source),
        lambda:control.controls((-1,1),1024,range_owner=range_owner,source_ranges=wrongsource),
        lambda:control.controls((-1,1),1024,range_owner=range_owner,source_ranges=missing),
        lambda:control.controls((-1,1),1024,range_owner=range_owner,source_ranges=unclosed),
        lambda:control.controls((-1,1),1024,range_owner=object(),source_ranges=source))
    for callback in attempts:
        try:callback()
        except (ValueError,TypeError):count+=1
        else:raise AssertionError('Wrong frequency/domain/source/owner or unsupplied original24 admitted')
    return dict(passed=True,candidate_domain_source_issued_object_full24_live_owner_depth_guards=count)


@current.native.inlet.source_precision
def run():
    began=time.monotonic();path=current.HERE/current.NAME;raw=gzip.decompress(path.read_bytes());manifest=json.loads(raw)
    assert manifest[current.GATE] and manifest['candidate_N']==1024
    assert manifest['full24_original_C1_integral_range_transport_enclosed']
    assert manifest['actual_finite_picard_and_residual_C0_Z_ranges_installed']
    flags=('actual_five_controls_installed','functional_terminal_identity_solved','current_whole_N_selected',*current.packets.OPEN)
    assert all(manifest[flag] is False for flag in flags)
    for name,digest in manifest['input_hashes'].items():assert current.sha(name)==digest,name
    bridge,_=current.native.inlet.native_bridge_owner();checks={}
    with current.native.inlet.CheckedSourceRuntime():
        owner=current.OriginalFullPredicate24CellControls(bridge);live=owner.assemble()
        with mp.workdps(owner.ctx.dps+40):
            checks['actual_original24_source_affine_memory_and_budget_inventory']=source_chain(owner,live,manifest['actual_original24_full_predicate_source_range_bridge'])
            print('Actual original24 source affine memory and budget inventory PASS',flush=True)
            checks['actual_new_Rc_joint_quotient_target']=joint_target(owner,live,manifest['actual_original24_full_predicate_source_range_bridge'])
            print('Actual newly issued Rc joint target PASS',flush=True)
        result=owner.finite_controls(live)
        with mp.workdps(owner.ctx.dps+40):
            checks['exact_original24_source_integrand_radius_phase_binding']=exact_function_binding(owner,live,result)
            checks['equivalent_original_density_C0_Z_definition']=equivalent_density_definition(result)
            print('Exact original24 source graph radius phase integrands PASS',flush=True)
            checks['actual_finite_control_and_numerical_residual']=finite_repair(owner,live,result,manifest['actual_original24_finite_control_diagnostics'])
            checks['actual_original_source_control_guards']=guards(owner,live,result)
            print('Actual finite control numerical residual and guards PASS',flush=True)
    receipt=dict(all_passed=True,**{current.GATE:True},source_family=owner.family,candidate_N=1024,**checks,
        compressed_producer_report=dict(filename=current.NAME,lossless_original_json_sha256=hashlib.sha256(raw).hexdigest(),
            compressed_bytes=path.stat().st_size,uncompressed_bytes=len(raw)),
        full24_original_C1_integral_range_transport_enclosed=True,actual_finite_picard_and_residual_C0_Z_ranges_installed=True,
        **dict.fromkeys(flags,False),input_hashes={**manifest['input_hashes'],current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began,
        scope='Actual original13+five exact full-predicate O2+six downstream affine source records and new joint Rc target bind the exact original24 control functions. Finite C0/Z iterations and numerical residuals checked, outer contraction insufficient. No fixed point, terminal/global N, actual recursion or full corrected NS admission.')
    (current.HERE/current.RECEIPT).write_bytes(json.dumps(current.encode(receipt),indent=2).encode()+b'\n')
    print('Actual original full-predicate24 finite control residual ranges PASS',flush=True);return receipt


if __name__=='__main__':run()
