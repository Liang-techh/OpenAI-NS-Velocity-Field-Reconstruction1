"""Independent declarative bound tables for compact physical pressure.

This successor reads saved data. It imports no mathematical implementation
from the failed pressure entry or any ancestor. No source integral is run.
"""
from collections import Counter
from itertools import product
from math import comb, factorial


def derive(source, journal):
    if not __debug__:
        raise RuntimeError('Certificate guards require unoptimized Python')
    packet=source['input'];data=packet['data'];geo=data['receipts']['geometry']
    assert packet['schema']=='compact_prefix_bound_table_data_v1' and packet['ancestor_math_calls']==0
    assert data['schema']=='compact_prefix_pressure_certificate_input_v1' and data['source_evaluations_performed']==0
    lineage=data['receipts']['logC'];memory=data['receipts']['memory']
    upstream=packet['cached_upstream_result']
    assert upstream==source['cached_upstream_proof']
    assert packet['cached_upstream_result_raw_sha256']=='417636494e5346a42261a2c9c37a69b6ef3e3e8b03afbbad5f4a8350580adb42'
    assert geo['schema']=='compact_native_numeric_volume_geometry_trace_v2'
    assert lineage['schema']=='compact_logC_parameter_lineage_join/v1'
    assert upstream['native_BG_exact_numeric_baseline_C3_identity_proved'] is True
    assert upstream['relative_mass_C3_source_identity_proved'] is True
    assert upstream['passed'] is True

    # Literal source equations, not a generic source DAG evaluator. Each
    # row is a defining equation; substitution is restricted to this table.
    nodes=dict(geo['native_Rc_source_graph']['graph_nodes'])
    nodes.update(data['native_geometry_literal_supplement']['graph_nodes'])
    constants={1:1,2:40,4:11,48:2,2083:10,79968:110,84529:-7,205156:18}
    for index,value in constants.items():
        row=nodes[str(index)]
        assert (row['operation'],row['numerator'],row['denominator'])==('exact_rational',value,1)
    equations={
        5:('sum',[3,4]),81513:('sum',[35560,5]),81514:('product',[2083,81513]),
        84530:('sum',[84529,84528]),84531:('sum',[81514,84530]),
        84533:('product',[79968,84532]),250222:('sum',[205156,3,250215]),
        250457:('product',[85743,250456]),260862:('product',[251066,260181])}
    for index,(operation,arguments) in equations.items():
        assert nodes[str(index)]['operation']==operation
        assert nodes[str(index)]['arguments']==arguments
    for index,expression,variable,value in ((85743,84533,84528,1),(251066,250457,250215,48)):
        row=nodes[str(index)]
        assert row['operation']=='function_substitution'
        assert (row['expression'],row['variable'],row['value'])==(expression,variable,value)
        assert row['Z_independent_substitution'] is True
    for index,argument in ((3,2),(84532,84531),(250456,250222)):
        row=nodes[str(index)]
        assert row['operation']=='analytic_unary' and row['name']=='exp' and row['argument']==argument
    # Native log Rc is log110+10(logC+logP)-7+1+18+exp40+2.
    # The six columns are log110,logC,logP,exp40,constant,logx.
    native_terms=[(1,0,0,0,0,0),(0,10,10,0,0,0),
                  (0,0,0,0,-7,0),(0,0,0,0,1,0),
                  (0,0,0,0,18,0),(0,0,0,1,0,0),
                  (0,0,0,0,2,0),(0,0,0,0,0,1)]
    numeric_terms=[(1,0,0,0,0,0),(0,10,0,0,0,0),(0,0,10,0,0,0),
                   (0,0,0,1,0,0),(0,0,0,0,14,0),(0,0,0,0,0,1)]
    native_columns=[sum(t[j] for t in native_terms) for j in range(6)]
    numeric_columns=[sum(t[j] for t in numeric_terms) for j in range(6)]
    assert native_columns==numeric_columns==[1,10,10,1,14,1]
    assert geo['numeric_cap_geometry']['exact_source_formula']=='log(Rupper)=log(110)+10*logC+10*logP+exp(40)+14+log(x)'
    assert geo['source_formula_chain']['source_formula_match'] is True
    assert geo['source_formula_chain']['logP_source']=='original_core_numeric.logP=exp(40)+11'
    leaf=nodes['35560'];assert leaf['source_path']=='common.core.logC'
    assert leaf==lineage['native_logC_target']['node']
    assert leaf['exact_fixed_source_parameter'] and leaf['parameter_enclosure_not_selected_value'] and leaf['Z_independent']
    bg=geo['native_Rc_source_graph']['compact_binding']['background']
    assert (bg['radius'],bg['phase'],bg['domain'])==(260862,260181,[1,48])
    assert bg['fields']['V']==[0]*7 and bg['radial']['V_y']==[0]*7
    assert upstream['exact_function_results']['physical_BG_Q_C2'][1:]==['0','0']
    measure=data['cached_C2_measure_normalization']
    assert (measure['upper_R2_mass'],measure['lower_Rc2_mass'])==('3/8','3/2')
    assert measure['exact_rational_check_passed'] and measure['new_exact_rational_mass_checks']==0
    assert measure['cached_source_sha256']=='1072a7882c1b7f5ce62ddde08650acd48e675e890d484350f8bb99c3a176812d'
    assert lineage['join_assessment']['native_numeric_logC_source_identity_join'] is True
    selected=lineage['numeric_logC_target']['selected_logCstar_exact_record']
    assert selected['lower_exact_mpf_tuple']==selected['upper_exact_mpf_tuple']
    assert lineage['numeric_logC_target']['parameter_is_midpoint_or_range_endpoint'] is False
    fixed=lineage['fixed_leading_N0_identity']
    js=[geo['numeric_cap_geometry']['selected_N_exact']['J'],
        fixed['native_branch_selected_integer']['J'],fixed['C3_repair_selected_integer']['J']]
    assert all(type(v) is int for v in js) and js==[1358356628656378313]*3
    geometry={'passed':True,'native_equations':len(equations)+len(constants)+5,
              'columns':['log110','logC','logP','exp40','constant','logx'],
              'native_terms':native_terms,'numeric_terms':numeric_terms,
              'native_sum':native_columns,'numeric_sum':numeric_columns,
              'residual':[a-b for a,b in zip(native_columns,numeric_columns)],
              'scalar_logC_conditional_on_cached_source_provider':True,
              'selected_J':js[0],'physical_R_evaluations':0,
              'coordinate':'x=R/Rc in [1,2]; dy=dx/x exactly once',
              'cached_measure':data['cached_C2_measure_normalization']}
    journal({'event':'new_declarative_geometry_table_closed','columns':6})

    # One fixed scalar, not a theorem for arbitrary s. The cached directed
    # construction certificate and its homogeneous term are explicit premises.
    rows=memory['selected_pressure_rows'];assert len(rows)==3
    assert memory['schema']=='compact_exact_memory_enclosure_join_trace_v1'
    assert memory['authenticated_artifacts']['raw_report']['sha256']=='57f5495394ab864ace99905fac07aab5bd315ff420329704121de121312b8d4a'
    lower=rows[0]['source_enclosure_lower'];upper=rows[0]['source_enclosure_upper']
    for raw_tuple in (lower,upper,rows[0]['raw_directed_error_tuple']):
        assert len(raw_tuple)==4 and all(type(v) is str and str(int(v))==v for v in raw_tuple)
    lo=[int(v) for v in lower];hi=[int(v) for v in upper]
    assert lo[0]==hi[0]==0 and lo[1]>0 and hi[1]>0
    assert lo[3]==lo[1].bit_length() and hi[3]==hi[1].bit_length()
    assert hi[2]-lo[2]==1 and lo[1]<=2*hi[1]
    for row,index,z in zip(rows,[5,11,17],['-.3','0','.3']):
        assert row['zero_based_pressure_row_index']==str(index) and row['Z']==z
        assert row['chart']=='O2_axial' and row['native_coordinate']=='1'
        assert row['selector']=='$.pressure_row_samples[%d]'%index
        assert row['source_family']==leaf['source_family']
        assert row['mstar_field_selector']=='$.pressure_row_samples[%d].leading_scalars.mean'%index
        assert row['source_enclosure_lower']==lower and row['source_enclosure_upper']==upper
        assert row['raw_directed_error_tuple']==rows[0]['raw_directed_error_tuple']
        cert=row['positive_receipt']
        assert cert['selector']=='$.pressure_row_samples[%d].positive_mean_memory_receipt'%index
        assert (cert['power'],cert['kernel_rate'])==('1','1')
        assert cert['positive_actual_memory_lower']==lower
        assert cert['method']=='peak-normalized Taylor2/fourth remainder with tangent tails'
        assert cert['finite_endpoint_rounding_bound']=={'lower':['0']*4,'upper':['0']*4}
        assert cert['local_axial_zero_does_not_zero_memory'] is True
    target=memory['mstar_source_contract']['target_call']
    assert target['class']=='OriginalAxialRateMemory' and target['evaluate_method']=='evaluate'
    assert target['evaluate_arguments']=={'s':'1','power':'1'} and target['constructor_kernel_rate']=='1'
    recipe=memory['mstar_source_contract']['source_recipe']
    cached_recipe=upstream['upstream_mean_identity']['axial_transport']['mstar_binding']['recipe']
    assert recipe['ordinary_Z_independent'] is True and recipe['homogeneous_term']=='exp(1-Y)'
    assert recipe['exact_mstar'].replace(' ','')==cached_recipe['definition'].split(';')[0].replace(' ','')
    assert cached_recipe['ordinary_Z_independent'] and cached_recipe['initial_memory']==1 and cached_recipe['kernel_rate']==1
    assert memory['certificate_assessment']['raw_report_checks']['checks_passed'] is True
    fixed_memory={'passed':True,'conditional_on_saved_directed_certificate':True,
                  'scalar_domain':{'s':'1','power':'1','rate':'1','ordinary_Z_independent':True},
                  'lower_exact_tuple':lower,'upper_exact_tuple':upper,
                  'directed_error_tuple':rows[0]['raw_directed_error_tuple'],
                  'source_recipe':recipe,'homogeneous_term_retained':True,
                  'universal_in_s_theorem':False,'independent_spot_checks_used_as_proof':False,
                  'midpoint_or_decimal_selected':False,'new_quadratures':0}
    journal({'event':'new_fixed_scalar_table_bound','cached_certificate_premise':True})

    # Independent integer coefficient tables. A monomial is a seven-column
    # exponent vector, rather than an expression tree or old normalizer.
    # Variables B,e0,e1,e2,u0,u1,u2; b is signed in the signed tables below.
    width=7;unit=[tuple(int(i==j) for i in range(width)) for j in range(width)]
    zero=(0,)*width
    def multiply(left,right):
        out=Counter()
        for (a,ca),(b,cb) in product(left.items(),right.items()):
            out[tuple(x+y for x,y in zip(a,b))]+=ca*cb
        return +out if all(v>=0 for v in out.values()) else Counter({k:v for k,v in out.items() if v})
    def add(*tables):
        out=Counter()
        for table in tables:
            for monomial,value in table.items():out[monomial]+=value
        return Counter({k:v for k,v in out.items() if v})
    def derivative(table):
        out=Counter()
        for exponents,c in table.items():
            for old,new in ((1,2),(2,3),(4,5),(5,6)):
                if exponents[old]:
                    q=list(exponents);q[old]-=1;q[new]+=1
                    out[tuple(q)]+=c*exponents[old]
        return Counter({k:v for k,v in out.items() if v})
    B,e0,e1,e2,u0,u1,u2=[Counter({m:1}) for m in unit]
    blank=Counter();negative_B2={m:-v for m,v in multiply(B,B).items()}
    signed=[e0,e1,add(multiply(add(B,e0),add(B,e0)),negative_B2),
            multiply(u0,e1),multiply(u0,u0),multiply(u0,u1)]
    bgQ=[B,blank,blank];bgU=[blank,blank,blank]
    errQ=[e0,e1,e2];errU=[u0,u1,u2]
    def product_error(a,b,ea,eb,k):
        return add(*(multiply(p[j],q[k-j]) for p,q in ((a,eb),(b,ea),(ea,eb)) for j in range(k+1)))
    comparison=[errQ[:2],errQ[1:],
                [product_error(bgQ,bgQ,errQ,errQ,k) for k in range(2)],
                [product_error(bgU,bgQ[1:],errU,errQ[1:],k) for k in range(2)],
                [product_error(bgU,bgU,errU,errU,k) for k in range(2)],
                [product_error(bgU,bgU[1:],errU,errU[1:],k) for k in range(2)]]
    table_rows=[]
    for name,signed0,expected in zip(['q','qZ','q2','UqZ','U2','UUZ'],signed,comparison):
        pair=[signed0,derivative(signed0)]
        for k,polynomial in enumerate(pair):
            absolute={m:abs(c) for m,c in polynomial.items()}
            assert absolute==dict(expected[k]),(name,k)
            table_rows.append({'channel':name,'ordinary_Z_order':k,
                 'signed_coefficients':[[list(m),c] for m,c in sorted(polynomial.items())],
                 'absolute_coefficients':[[list(m),abs(c)] for m,c in sorted(polynomial.items())],
                 'archived_three_term_bound_equal':True})
    pressure=[]
    for k in range(3):
        for j in range(k+1):
            taylor_multiplier=factorial(j)*comb(k,j)
            ordinary_multiplier=comb(k,j)
            assert taylor_multiplier==ordinary_multiplier*factorial(j)
            pressure.append({'target_order':k,'coefficient_order':j,
                             'Taylor_coefficient_multiplier':taylor_multiplier,
                             'ordinary_coefficient_multiplier':ordinary_multiplier,
                             'ordinary_from_Taylor_factor':factorial(j)})
    assert pressure[-1]['Taylor_coefficient_multiplier']==2
    volume=[]
    for k in range(2):
        # Independent coefficient A[channel,j] is tagged separately; equality
        # follows term by term from the twelve verified coefficient tables.
        summands=[{'coefficient_abs_port':'A_%d_%d'%(i,j),'channel':rowname,
                   'channel_order':k-j,'binomial':comb(k,j)}
                  for i,rowname in enumerate(['q','qZ','q2','UqZ','U2','UUZ']) for j in range(k+1)]
        volume.append({'ordinary_Z_order':k,'summands':summands,'termwise_bound_proved':True})
    dominance={'passed':True,'signed_variables':['b','deltaQ0','deltaQ1','deltaQ2','deltaU0','deltaU1','deltaU2'],
               'envelope_variables':['B','E0','E1','E2','U0','U1','U2'],
               'channel_tables':table_rows,'channel_order_checks':12,
               'volume_bound_tables':volume,'volume_order_checks':2,
               'ordinary_pressure_conversion':pressure,
               'hypotheses':['abs(Q_BG)<=B, (Q_BG)_Z=(Q_BG)_ZZ=0',
                             'U_BG=0 through ordinary Z2',
                             'abs(d_Z^k deltaQ)<=e_k; abs(d_Z^k deltaU)<=u_k, k<=2',
                             'A_i,j=abs(d_Z^j a_i), j<=1'],
               'conclusion':'abs(d_Z^k(v_N-v_BG)) bounded by the displayed coefficient sum, k=0,1',
               'source_inequality_hypotheses_discharged':False}
    journal({'event':'new_integer_bound_tables_proved','channels':6,'channel_order_checks':12,'volume_checks':2})

    prefix={'conditional_inequality_proved':True,'ordinary_Z_orders':[0,1,2,3],
            'defined_history':'H_k(s)=exp(-s) I_k + integral_0^s exp(-(s-u)) rho_k(u)du',
            'premises':['abs(I_k)<=Icap_k','abs(rho_k)<=Dcap_k*indicator_S','0<=u<=s<=log(2)',
                        'support and s are ordinary-Z independent'],
            'argument':['exp(-s)<=1','0<exp(-(s-u))<=1','integrate the absolute inequality over S intersect [0,s]'],
            'conclusion':'sup_0<=v<=s abs(H_k(v))<=Icap_k+Dcap_k*measure(S intersect [0,s])',
            'incoming_retained':True,'extra_N_divisions':0,
            'native_C3_cap_to_C5_source_ports_proved':False,
            'remaining_obligations':[
                'C3 terminal_N_scaled_histories/N = native incoming 258646..258649 through ordinary Z3',
                'Same exact normalized targets imply the C3 controls and C5 system are the same unique branch through Z3',
                'C3 profile/density source rows = native deltaV 260375/260378/260382/260387 and Duhamel integrands',
                'Predecessor and band incoming minimum bounds this same native source; unioned U/U_y errors retained']}
    result={'schema':'original_compact_prefix_bound_table_certificate',
            'component':'compact_prefix_bound_table_certificate','passed':True,'metadata_only':False,
            'source_only':True,'input_digest':source['input_digest'],'runtime_model_metadata':None,
            'bound_table_defined':True,'signed_error_rows_defined':True,
            'signed_error_rows_scope':'Formal signed difference polynomials, not evaluated finite-N source values',
            'original_source_evaluations':0,'source_packet_defined':True,'source_ftc_transport_defined':True,
            'native_BG_exact_numeric_baseline_C3_identity_proved':True,'relative_mass_C3_source_identity_proved':True,
            'compact_prefix_geometry_logC_binding_proved':True,'compact_prefix_fixed_mstar_enclosure_proved':True,
            'compact_prefix_12channel_2volume_dominance_proved':True,
            'compact_prefix_source_binding_proved':False,'compact_volume_cap_same_source_bound_proved':False,
            'numeric_mean_enclosure_transferred':False,
            'prefix_pressure_certificate':{'geometry':geometry,'fixed_mstar':fixed_memory,
                                            'C1_absolute_dominance':dominance,'conditional_C3_prefix':prefix},
            'source_binding':{'passed':True,'conditional_on_saved_source_and_interval_contracts':True,
                              'actual_C3_cap_to_C5_port_bound':False,
                              'same_family_metadata_used_as_function_identity':False},
            'coordinate_mapping':{'x':'R/Rc','dy':'dx/x once','logR':'log110+10logC+10logP+exp40+14+logx'},
            'analytic_integral_definition':{'signed_Evol':'integral_Rc^2Rc R^2*(v_N-v_BG)dy',
                                           'value_computed':False,'independent_old_reference_offset_retained':True},
            'source_inspection_audit':{'original_source_evaluations':0,'old_source_consumer_calls':0,
                 'old_proof_functions_called':0,'old_operator_calls':0,'old_source_normalizer_calls':0,
                 'quadratures':0,'numeric_field_evaluations':0,'temporal_recurrence_calls':0,
                 'new_channel_order_checks':12,'new_volume_order_checks':2,
                 'new_pressure_conversion_checks':6},
            'remaining_scope':prefix['remaining_obligations']+['Signed Evol and the independent old reference offset remain unevaluated; global construction incomplete.']}
    for name in ('retained_cap_provenance','retained_construction_data','retained_integral_term_inventory'):
        result[name]=upstream[name]
    result['bound_table_certificate']={'schema':'compact_prefix_bound_table_manifest_v1',
         'details_selector':'$.prefix_pressure_certificate','geometry_source_columns':6,
         'fixed_scalar_intervals_bound':1,'channel_order_checks':12,'volume_order_checks':2,
         'ordinary_pressure_conversion_checks':6,'actual_C3_cap_to_C5_port_bound':False}
    for name in source['inherited_false_flags']:result[name]=False
    journal({'event':'new_bound_table_certificate_complete','actual_volume_cap_transfer':False})
    return result
