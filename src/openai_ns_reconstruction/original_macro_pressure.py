"""Apply actual pressure to new micro/macro source seams with explicit queries.

This generic kernel receives explicit chart/background registries. It never
constructs an ancestor or reruns an accepted seam; ordinary pressure algebra
is applied only to the new input identity witnesses supplied by its caller.
"""
import ast
import hashlib

from .original_core_numeric import ROOT,PREFIX
from .original_upstream_m15 import digest,pressure_congruence

def build_macro_pressure_traces(owner):
    """Recover physical Q/Qy from source-proved U, Uy and mean germs."""
    self=owner
    source_path='experiments/root_st073/'+PREFIX+'compliant_current_original_modified_Q4_pressure_C3_functions.py'
    raw=(ROOT/source_path).read_bytes();sha=hashlib.sha256(raw).hexdigest()
    if self.hashes.get(source_path.rsplit('/',1)[-1])!=sha:raise ValueError('Pinned physical pressure formula required')
    text=raw.decode('utf-8');module=ast.parse(text)
    function=next(r for r in module.body if isinstance(r,ast.FunctionDef) and r.name=='physical_pressure_fields')
    body=ast.get_source_segment(text,function)
    if not all(s in body for s in ('my=a.sub(U,m)','Q=div(numerator(U,m));Qy=div(numerator(Uy,my))',
            'b.mul(d,b.shift(mean))',"m=a.scale(fields['m'],Pstar)")):
        raise ValueError('Original pressure formula changed')
    formula=dict(passed=True,source=source_path,sha256=sha,
        function_AST_sha256=hashlib.sha256(ast.dump(function,include_attributes=False).encode()).hexdigest(),
        Q='[2Z U-(1-delta)Z m-(1-Z^2)mz]/(1-delta Z^2)',
        Qy='[2Z Uy-(1-delta)Z(U-m)-(1-Z^2)(Uz-mz)]/(1-delta Z^2)',
        actual_m_y_equals_U_minus_m=True,ordinary_Q_Z2_requires_m_Z3=True,
        ordinary_Qy_Z1_requires_U_Z2_m_Z2_Uy_Z1=True,source_executed=False)
    pstar=self.pstar;factor_rows=[];denominators=set();specs={}
    for name in self.chart_names:
        loop=self.l['charts'][name];modified=self.m['charts'][name]
        background=self.background_charts[name]
        if loop['coordinate']!=background['phase'] or loop['radius']!=background['radius']:
            raise ValueError('Independent original leading geometry mismatch')
        specs[name]=dict(chart=name,coordinate=loop['coordinate'],radius=loop['radius'],
            physical=modified['physical'],boundary=modified['finite_boundary'],
            fields=loop['modified_fields'],radial=modified['radial'],
            endpoint=self.m['endpoints'][name],upper_query=self.chart_upper_queries[name],
            leading=dict(V=background['fields']['V'],m=background['fields']['m'],Vy=background['radial']['V_y']))
        for quantity,inputs,count in (('U',loop['modified_fields']['V'],3),('U_y',modified['radial']['V_y'],2),('m',loop['modified_fields']['m'],4)):
            for k in range(count):
                node=modified['physical'][quantity][k];value=inputs[k];r=self.nodes[node]
                zero=lambda i:self.nodes[i].get('operation')=='exact_rational' and self.nodes[i].get('numerator')==0
                valid=(zero(node) and zero(value)) or (r.get('operation')=='product' and sorted(r['arguments'])==sorted([pstar,value]))
                if not valid:raise ValueError('Original physical factor mismatch '+str((name,quantity,k)))
                factor_rows.append(dict(chart=name,quantity=quantity,ordinary_Z_order=k,physical_node=node,
                    actual_input_node=value,Pstar_node=pstar,passed=True))
        for key in ('Q','Q_y'):
            r=self.nodes[modified['physical'][key][0]]
            if r.get('operation')!='positive_quotient' or r.get('source_positive_certificate')!='same original L=1-delta*Z^2':
                raise ValueError('Original positive pressure denominator required')
            denominators.add(r['denominator'])
        if not self.native['pressure_ports'][name]['m_y_equals_U_minus_m']:
            raise ValueError('Original physical mass ODE required')
    if len(denominators)!=1:raise ValueError('Same original pressure denominator required')
    formula['actual_common_positive_L_node']=next(iter(denominators))
    # Keep native variables and substitutions scoped while proving R and
    # Pstar independent of the ordinary axial variable Z.
    seen=set();coordinates={spec['coordinate'] for spec in specs.values()}
    def fixed(i,bound):
        key=i,tuple(sorted(bound))
        if key in seen:return
        seen.add(key);row=self.nodes[i];op=row['operation']
        if op in ('sum','product'):
            for child in row['arguments']:fixed(child,bound)
        elif op in ('negative','analytic_unary'):fixed(row['argument'],bound)
        elif op=='positive_quotient':fixed(row['numerator'],bound);fixed(row['denominator'],bound)
        elif op=='function_substitution':
            if not row['Z_independent_substitution'] or row['variable']==20:raise ValueError('Moving axial radius substitution')
            fixed(row['value'],bound);fixed(row['expression'],bound|{row['variable']})
        elif op=='exact_rational':pass
        elif op=='bound_variable':
            if i not in bound:raise ValueError('Free variable in physical radius or parameter')
        elif op in ('current_original_source_parameter','original_source_parameter','exact_positive_integer_expression'):
            if not row.get('Z_independent'):raise ValueError('Axially varying physical parameter')
        else:raise ValueError('Unproved radius operation '+op)
    for spec in specs.values():fixed(spec['radius'],coordinates)
    fixed(pstar,set())
    L=formula['actual_common_positive_L_node'];lr=self.nodes[L]
    negatives=[i for i in lr.get('arguments',[]) if self.nodes[i].get('operation')=='negative']
    if lr.get('operation')!='sum' or len(negatives)!=1 or sorted(lr['arguments'])!=sorted([1,negatives[0]]):
        raise ValueError('Actual L=1-delta*Z^2 source expression required')
    product=self.nodes[self.nodes[negatives[0]]['argument']]
    square=[i for i in product.get('arguments',[]) if self.nodes[i].get('operation')=='product' and self.nodes[i].get('arguments')==[20,20]]
    if product.get('operation')!='product' or len(product['arguments'])!=2 or len(square)!=1 or self.nodes[20]!={'operation':'bound_variable','name':'Z'}:
        raise ValueError('Actual ordinary Z pressure denominator source required')
    delta=next(i for i in product['arguments'] if i!=square[0]);fixed(delta,set())
    formula['actual_delta_node']=delta;formula['actual_Z_node']=20
    formula['delta_Z_independent']=True
    radius_guard=dict(passed=True,checked_scoped_nodes=len(seen),ordinary_R_Z_derivatives_zero=True,
        Pstar_Z_independent=True,delta_node=delta,delta_Z_independent=True,parameter_values_not_selected=True)

    def port(spec,key,k,point,*,boundary=False):
        expression=spec['boundary' if boundary else 'physical'][key][k]
        recipe=dict(operation='local_endpoint_recipe',expression=expression,
            variable=spec['coordinate'],value=point,Z_independent_substitution=True)
        out=dict(recipe=recipe,recipe_is_archived_node=False,full_chart_function_node=expression,passed=True)
        # Existing upper boundary wrappers are bound separately; physical
        # U/Uy/m/Q/Qy wrappers absent from these archives remain recipes.
        if boundary and point==spec['upper_query']:
            node=spec['endpoint']['boundary'][key][k]
            expected=dict(operation='function_substitution',expression=expression,variable=spec['coordinate'],
                value=point,Z_independent_substitution=True)
            if self.nodes[node]!=expected:raise ValueError('Foreign upper endpoint boundary wrapper')
            out.update(node=node,actual_wrapper_bound=True)
        return out

    def source_port(row,side,spec,key,k,point,*,leading=False):
        expression=(spec['leading'] if leading else spec['physical'])[key][k]
        actual=row[side+'_node']
        if actual==expression:return
        expected=dict(operation='function_substitution',expression=expression,variable=spec['coordinate'],
            value=point,Z_independent_substitution=True)
        if self.nodes[actual]!=expected:raise ValueError('Foreign source germ port '+str((side,key,k,actual,expression)))

    native_wiring_sha=digest(self.derivative['charts'])
    def issued_input(row,quantity,source,*,leading=False):
        identity=row.get('actual_function_identity',{})
        namespace='leading' if leading else 'modified'
        if (identity.get('passed') is not True or identity.get('quantity')!=quantity
                or identity.get('namespace')!=namespace):
            raise ValueError('Foreign modified/independent-leading input identity namespace')
        if quantity in ('Uy','V_y') and identity.get('issued_derivative_wiring_sha256')!=native_wiring_sha:
            raise ValueError('Actual physical derivative wiring missing for '+quantity)
        if identity.get('exact_endpoint_fingerprints') is True:
            group=({'V':'leading_V','m':'leading_m','V_y':'leading_Vy'} if leading else
                   {'U':'modified_U','Uy':'modified_Uy','m':'modified_m'})[quantity]
            direct=source.get('direct_endpoint_fingerprints',{})
            rows=direct.get('groups',{}).get(group,{}).get('rows',[])
            order=row['ordinary_Z_order']
            if (identity.get('direct_group')!=group or len(rows)<=order
                    or rows[order].get('passed') is not True):
                raise ValueError('Direct fingerprint has a foreign source group')
            witness=rows[order]
            for side in ('left','right'):
                if (witness.get(side+'_node')!=row.get(side+'_node')
                        or witness.get(side+'_endpoint_sha256')!=identity.get(side+'_endpoint_sha256')
                        or direct.get(side+'_endpoint_scope')!=source[side+'_scope']):
                    raise ValueError('Direct fingerprint has a foreign actual node or endpoint scope')

    result={}
    for name,(left_name,right_name) in self.seam_pairs.items():
        source=self.source['seams'][name];left=specs[left_name];right=specs[right_name]
        left_query,right_query=self.seam_queries[name]
        for query in (left_query,right_query):
            row=self.nodes[query]
            if row.get('operation')!='exact_rational' or row.get('denominator')!=1:
                raise ValueError('Exact archived native endpoint query required')
        if not source['passed']:raise ValueError('Actual upstream field source proof required')
        for side,spec,point in (('left',left,left_query),('right',right,right_query)):
            scope=source[side+'_scope']
            if scope['coordinate_node']!=spec['coordinate'] or scope['native_query_node']!=point:
                raise ValueError('Actual source theorem has foreign endpoint scope')
        geometry=source.get('radius_phase_identity',{})
        if (geometry.get('passed') is not True
                or geometry.get('common_radius_source_identity_proved') is not True
                or geometry.get('global_phase_transfer') is not True
                or geometry.get('source_values_selected') is not False
                or geometry.get('numeric_caps_used') is not False):
            raise ValueError('Issued common physical-radius/source-phase identity required')
        geometry_rows=geometry.get('rows',[])
        if len(geometry_rows)!=3 or {r.get('role') for r in geometry_rows}!={'radius','global_offset','actual_phase'}:
            raise ValueError('Actual source geometry role coverage required')
        for row in geometry_rows:
            role=row['role']
            for side,chart in (('left',left_name),('right',right_name)):
                if (row.get(side+'_node')!=self.l['charts'][chart][role]
                        or row.get(side+'_scope')!=source[side+'_scope']):
                    raise ValueError('Foreign physical-radius or global-phase source port')
        receipt=geometry.get('source_receipt_path','')
        if (not receipt.startswith('experiments/root_st073/')
                or self.hashes.get(receipt.rsplit('/',1)[-1])!=geometry.get('source_receipt_sha256')):
            raise ValueError('Actual source geometry receipt is not pinned')
        geometry_sha=digest(geometry)
        traces=[];bindings=[];modified_rows=[]
        for quantity,key,count,source_key in (('U','U',3,'modified_U_rows'),('Uy','U_y',2,'modified_Uy_rows'),('m','m',4,'modified_m_rows')):
            rows=source[source_key]
            if len(rows)!=count:raise ValueError('Complete original modified source germ required')
            for k,row in enumerate(rows):
                if not row['passed'] or row['ordinary_Z_order']!=k:raise ValueError('Unproved ordinary source derivative')
                issued_input(row,quantity,source)
                ports={}
                for side,spec,point in (('left',left,left_query),('right',right,right_query)):
                    source_port(row,side,spec,key,k,point)
                    ports[side]=port(spec,key,k,point)
                    bindings.append(dict(side=side,quantity=quantity,ordinary_Z_order=k,**ports[side]))
                modified_rows.append(dict(quantity=quantity,ordinary_Z_order=k,source_row_sha256=digest(row),
                    left_node=row['left_node'],right_node=row['right_node'],passed=True))
                if quantity!='m' or k<3:
                    traces.append(dict(quantity=quantity,ordinary_Z_order=k,**ports,
                        common_function_sha256=digest([name,quantity,k,row]),passed=True))
        leading=source['leading_Qy_inputs']
        if not leading['passed'] or len(leading['rows'])!=8:raise ValueError('Independent leading source inputs required')
        leading_rows=[]
        for row in leading['rows']:
            if not row['passed']:raise ValueError('Unproved original leading derivative')
            issued_input(row,{'Vy':'V_y'}.get(row['quantity'],row['quantity']),source,leading=True)
            key={'V_y':'Vy','Vy':'Vy'}.get(row['quantity'],row['quantity']);k=row['ordinary_Z_order']
            for side,spec,point in (('left',left,left_query),('right',right,right_query)):
                source_port(row,side,spec,key,k,point,leading=True)
            leading_rows.append(dict(quantity=key,ordinary_Z_order=k,source_row_sha256=digest(row),
                left_node=row['left_node'],right_node=row['right_node'],passed=True))
        for quantity,key,count in (('Q','Q',3),('Qy','Q_y',2)):
            for k in range(count):
                a=port(left,key,k,left_query);b=port(right,key,k,right_query)
                traces.append(dict(quantity=quantity,ordinary_Z_order=k,left=a,right=b,
                    physical_formula=formula,source_field_identity_sha256=digest(source),
                    theorem_application_not_raw_fingerprint=True,passed=True))
        for k in (0,1):
            ports={}
            for side,spec,point in (('left',left,left_query),('right',right,right_query)):
                product=spec['boundary']['R_Q_y'][k];r=self.nodes[product]
                if r.get('operation')!='product' or sorted(r['arguments'])!=sorted([spec['radius'],spec['physical']['Q_y'][k]]):
                    raise ValueError('Original radius factor must occur once')
                ports[side]=port(spec,'R_Q_y',k,point,boundary=True)
            traces.append(dict(quantity='RQy',ordinary_Z_order=k,**ports,common_radius_and_Qy_source_deduction=True,
                radius_source_identity_sha256=geometry_sha,passed=True))
        algebra=pressure_congruence(source,formula,loop_identity=self.loop_identity,derivative_identity=self.derivative)
        result[name]=dict(passed=True,left_chart=left_name,right_chart=right_name,left_endpoint=self.nodes[left_query]['numerator'],right_endpoint=self.nodes[right_query]['numerator'],
            left_endpoint_node=left_query,right_endpoint_node=right_query,left_coordinate=left['coordinate'],right_coordinate=right['coordinate'],
            left_scope=source['left_scope'],right_scope=source['right_scope'],
            trace_rows=traces,modified_source_rows=modified_rows,leading_source_rows=leading_rows,
            pressure_algebra_witness=algebra,
            leading_Qy_common=dict(passed=True,input_rows=leading['rows'],original_leading_namespace_bound=True,
                physical_formula=formula,source_values_selected=False,leading_output_alias_assigned=False),
            endpoint_binding_rows=bindings,endpoint_wrapper_binding_passed=True,
            missing_physical_wrappers_are_local_recipes=True,extra_mean_Z3_input_retained=True,
            radius_phase_identity=geometry,radius_source_identity_sha256=geometry_sha,
            D_equality_rows=[dict(ordinary_Z_order=k,passed=True,physical_R_once=True,
                radius_source_identity_sha256=geometry_sha,
                recipe='R*(Qy_modified-Qy_original_leading)',signed_endpoint_value_evaluated=False) for k in (0,1)],
            incoming_history_retained=True,source_values_selected=False,numeric_caps_used=False)
    return dict(passed=True,joins=result,physical_pressure_formula=formula,physical_factor_rows=factor_rows,
        radius_parameter_guard=radius_guard,trace_row_count=sum(len(r['trace_rows']) for r in result.values()),
        source_values_selected=False,numeric_caps_used=False,ancestor_replayed=False,
        original_graph_not_mutated=True,endpoint_recipes_not_registered_as_graph_nodes=True)
