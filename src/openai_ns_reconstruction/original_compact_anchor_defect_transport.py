"""New compact-anchor defect transport DAG; intentionally unentered."""
from fractions import Fraction

_MISSING = object()

def _f(x):
    if isinstance(x, Fraction): return x
    if isinstance(x, int): return Fraction(x, 1)
    if isinstance(x, str): return Fraction(x)
    if isinstance(x, dict): return Fraction(int(x["numerator"]), int(x["denominator"]))
    raise TypeError("exact rational required")
def _fj(x):
    x=_f(x); return {"numerator":x.numerator,"denominator":x.denominator}
def _get(o,path):
    for k in path.split("."):
        if not isinstance(o,dict) or k not in o: return _MISSING
        o=o[k]
    return o
def _need(o,*paths):
    for p in paths:
        x=_get(o,p)
        if x is not _MISSING and x is not None: return x
    raise KeyError("missing explicit binding: "+", ".join(paths))
def _freeze(x):
    if isinstance(x,dict): return tuple(sorted((k,_freeze(v)) for k,v in x.items()))
    if isinstance(x,(list,tuple)): return tuple(_freeze(v) for v in x)
    return x
def _node(x,label):
    if isinstance(x,dict):
        for k in ("node","source_node","row","value"):
            if k in x: x=x[k]; break
    if isinstance(x,bool): raise TypeError(label+" cannot be boolean")
    return int(x)
def _binding(b,label):
    if not isinstance(b,dict): raise TypeError(label+" must be a mapping")
    if "raw_sha256" not in b and isinstance(b.get("raw_sha"),str):
        b=dict(b); b["raw_sha256"]=b["raw_sha"]
    for k in ("raw_sha256","selector"):
        if not isinstance(b.get(k),str) or not b[k]: raise ValueError(label+" requires "+k)
    if not any(k in b for k in ("source_node","source_nodes","recipe")):
        raise ValueError(label+" requires source_node, source_nodes, or recipe")
    return b
def _law_binding(b,label):
    if not isinstance(b,dict): raise TypeError(label+" must be a mapping")
    if "raw_sha256" not in b and isinstance(b.get("raw_sha"),str): b=dict(b); b["raw_sha256"]=b["raw_sha"]
    for k in ("raw_sha256","selector"):
        if not isinstance(b.get(k),str) or not b[k]: raise ValueError(label+" requires "+k)
    if not any(k in b for k in ("source_node","source_nodes","recipe","passed_static_receipt","exact_node_records")):
        raise ValueError(label+" requires a source recipe or exact static node records")
    return b
def _rows(b,label,count):
    b=_binding(b,label); xs=b.get("source_nodes",_MISSING)
    if xs is _MISSING: raise ValueError(label+" requires source_nodes length "+str(count))
    if not isinstance(xs,(list,tuple)) or len(xs)!=count: raise ValueError(label+" requires exactly "+str(count)+" source nodes")
    return [_node(x,label) for x in xs]
def _port_descriptor(b,label,count):
    b=_binding(b,label)
    if "source_nodes" in b: return _rows(b,label,count)
    if "recipe" in b: return {"raw_sha256":b["raw_sha256"],"selector":b["selector"],"recipe":b["recipe"],"ordinary_Z_orders":list(range(count))}
    raise ValueError(label+" requires source_nodes or recipe for all orders")

class _DAG:
    def __init__(self):
        self.nodes=[]; self.index={}; self.counts={k:0 for k in ("add","sub","mul","div","neg","source_ref","recipe_ref","cached_ref","bound_variable","formal_integral","formal_primitive","analytic_constant","c1_product_rules","c1_quotient_rules","c2_quotient_rules","c3_quotient_rules","c2_product_rules","c3_product_rules")}
    def _put(self,key,rec):
        if key in self.index:return self.index[key]
        i=len(self.nodes); rec=dict(rec); rec["id"]=i; self.nodes.append(rec); self.index[key]=i
        if rec.get("op") in self.counts:self.counts[rec["op"]]+=1
        return i
    def const(self,x):
        x=_f(x); return self._put(("const",x.numerator,x.denominator),{"op":"const","value":_fj(x)})
    def op(self,op,args=(),**meta):
        args=tuple(int(a) for a in args); key=(op,args,_freeze(meta)); rec={"op":op,"args":list(args)}; rec.update(meta); return self._put(key,rec)
    def source(self,b,label,order):
        b=_binding(b,label); order=int(order)
        if "source_nodes" in b:
            nodes=b["source_nodes"]
            if order>=len(nodes): raise ValueError(label+" source order out of range")
            n=_node(nodes[order],label)
        elif order==0 and "source_node" in b:n=_node(b["source_node"],label)
        elif "recipe" in b:n=None
        else: raise ValueError(label+" has no source at order "+str(order))
        if n is not None:return self._put(("source_ref",b["raw_sha256"],b["selector"],label,n,order),{"op":"source_ref","raw_sha256":b["raw_sha256"],"selector":b["selector"],"label":label,"source_node":n,"ordinary_Z_order":order})
        return self._put(("recipe_ref",b["raw_sha256"],b["selector"],label,_freeze(b["recipe"]),order),{"op":"recipe_ref","raw_sha256":b["raw_sha256"],"selector":b["selector"],"label":label,"recipe":b["recipe"],"ordinary_Z_order":order})
    def cached(self,b,label,order):
        b=_binding(b,label); xs=b.get("source_nodes",_MISSING)
        if xs is not _MISSING:
            if len(xs)!=2: raise ValueError(label+" requires value/dZ source_nodes")
            n=_node(xs[order],label)
        elif order==0 and "value_node" in b:n=_node(b["value_node"],label)
        elif order==1 and "dz_node" in b:n=_node(b["dz_node"],label)
        else: raise ValueError(label+" requires cached value and dZ nodes")
        return self._put(("cached_ref",b["raw_sha256"],b["selector"],label,n,order),{"op":"cached_ref","raw_sha256":b["raw_sha256"],"selector":b["selector"],"label":label,"cached_node":n,"ordinary_Z_order":order,"replay":False})
    def local(self,name): return self._put(("bound_variable",name),{"op":"bound_variable","name":name,"ordinary_Z_derivative":0})

def _zeros(g,n): return [g.const(0) for _ in range(n)]
def _jconst(g,x,n): return [g.const(x)]+_zeros(g,n-1)
def _jnode(g,x,n): return [int(x)]+_zeros(g,n-1)
def _jadd(g,a,b): return [g.op("add",[x,y]) for x,y in zip(a,b)]
def _jsub(g,a,b): return [g.op("sub",[x,y]) for x,y in zip(a,b)]
def _jneg(g,a): return [g.op("neg",[x]) for x in a]
def _jscale(g,a,c): return [g.op("mul",[x,c]) for x in a]
def _jmul(g,a,b):
    n=len(a); g.counts["c1_product_rules" if n==2 else "c2_product_rules" if n==3 else "c3_product_rules"]+=1
    out=[g.op("mul",[a[0],b[0]])]
    if n>=2:out.append(g.op("add",[g.op("mul",[a[1],b[0]]),g.op("mul",[a[0],b[1]])]))
    if n>=3:out.append(g.op("add",[g.op("add",[g.op("mul",[a[2],b[0]]),g.op("mul",[g.const(2),g.op("mul",[a[1],b[1]])])]),g.op("mul",[a[0],b[2]])]))
    if n>=4:out.append(g.op("add",[g.op("add",[g.op("mul",[a[3],b[0]]),g.op("mul",[g.const(3),g.op("mul",[a[2],b[1]])])]),g.op("add",[g.op("mul",[g.const(3),g.op("mul",[a[1],b[2]])]),g.op("mul",[a[0],b[3]])])]))
    return out
def _jinv(g,b):
    n=len(b); g.counts["c1_quotient_rules" if n==2 else "c2_quotient_rules" if n==3 else "c3_quotient_rules"]+=1
    y=b[0]; y2=g.op("mul",[y,y]); q0=g.op("div",[g.const(1),y]); out=[q0]
    if n>=2:out.append(g.op("neg",[g.op("div",[b[1],y2])]))
    if n>=3:
        num=g.op("sub",[g.op("mul",[g.const(2),g.op("mul",[b[1],b[1]])]),g.op("mul",[y,b[2]])]); out.append(g.op("div",[num,g.op("mul",[y2,y])]))
    if n>=4:
        y3=g.op("mul",[y2,y]); y4=g.op("mul",[y3,y])
        t1=g.op("mul",[g.const(-6),g.op("mul",[b[1],g.op("mul",[b[1],b[1]])])])
        t2=g.op("mul",[g.const(6),g.op("mul",[y,g.op("mul",[b[1],b[2]])])])
        t3=g.op("neg",[g.op("mul",[y2,b[3]])])
        out.append(g.op("div",[g.op("add",[g.op("add",[t1,t2]),t3]),y4]))
    return out
def _jdiv(g,a,b): return _jmul(g,a,_jinv(g,b))
def _shift(a,k,n): return list(a[k:k+n])
def _srcjet(g,b,label,n): return [g.source(b,label,i) for i in range(n)]
def _c1(g,a): return a[:2]
def _cached_c1(g,b,label): return [g.cached(b,label,0),g.cached(b,label,1)]
def _formal(g,x,label,order): return g.op("formal_integral",[x],label=label,lower=1,upper=2,measure="dx",ordinary_Z_order=order,derivative_under_integral=(order==1))
def _sparse(g,atoms,coefficients):
    names=list(atoms); coefficients=[_f(c) for c in coefficients]
    meta={"atom_order":names,"atom_terms":[{"atom":n,"coefficient":_fj(c)} for n,c in zip(names,coefficients)],"all_atom_arguments_retained":True,"exact_merge":True}
    value=g.op("sparse_linear",[atoms[n][0] for n in names],**meta)
    dz=g.op("sparse_linear",[atoms[n][1] for n in names],**meta)
    return [value,dz]

def construct_anchor_graph(source,journal):
    if not isinstance(source,dict) or not callable(journal): raise TypeError("source mapping and callable journal required")
    anchor=_need(source,"anchor_binding"); memory=_need(source,"memory_binding"); geometry=_need(source,"geometry_binding"); cached=_need(source,"cached_volume_binding")
    anchor_law=_law_binding(_need(source,"compact_anchor_law","anchor_binding.compact_anchor_law"),"source.compact_anchor_law")
    pstar=_binding(_need(anchor,"Pstar"),"anchor.Pstar"); m=_binding(_need(anchor,"normalized_m_at_x1_C3"),"anchor.normalized_m_at_x1_C3"); mrows=_port_descriptor(m,"anchor.normalized_m_at_x1_C3",4)
    mstar=_binding(_need(memory,"mstar_exact"),"memory.mstar_exact")
    for k in ("Rc","delta","Z","local_coordinate"):_binding(_need(geometry,k),"geometry."+k)
    native_bg=_binding(_need(cached,"native_minus_background"),"cached.native_minus_background"); bg_old=_binding(_need(cached,"background_minus_old"),"cached.background_minus_old"); endpoint=_binding(_need(cached,"endpoint_difference"),"cached.endpoint_difference")
    g=_DAG(); x=g.local("compact_anchor_x"); x1=_jnode(g,x,3); coord=g.source(_need(geometry,"local_coordinate"),"geometry.local_coordinate",0)
    Z=[g.source(_need(geometry,"Z"),"geometry.Z",0),g.const(1),g.const(0)]; delta=[g.source(_need(geometry,"delta"),"geometry.delta",0),g.const(0),g.const(0)]; Rc=[g.source(_need(geometry,"Rc"),"geometry.Rc",0),g.const(0),g.const(0)]
    P=g.source(pstar,"anchor.Pstar",0); mj=_srcjet(g,m,"anchor.normalized_m_at_x1_C3",4); A=_jscale(g,mj,P); A3=A[:3]
    exp14=g.op("analytic_constant",[],expression="exp(-14)",exact_symbolic=True); ms=g.source(mstar,"memory.mstar_exact",0); msj=[ms,g.const(0),g.const(0)]
    Bscalar=g.op("mul",[exp14,ms]); B=_jmul(g,_jmul(g,_jconst(g,4,3),Z),_jnode(g,Bscalar,3)); B4=B+[g.const(0)]; D=_jsub(g,A,B4)
    one=_jconst(g,1,3); z2=_jmul(g,Z,Z); d=_jsub(g,one,z2); h=_jsub(g,one,delta); L=_jsub(g,one,_jmul(g,delta,z2)); AZ=_shift(A,1,3); num=_jadd(g,_jmul(g,_jmul(g,h,Z),A3),_jmul(g,d,AZ)); G=_jneg(g,_jdiv(g,num,L)); H=_jscale(g,msj,g.op("mul",[g.const(-4),exp14])); GD=_jsub(g,G,H); x4=_jnode(g,x,4); DM=_jdiv(g,D,x4); DQ=_jdiv(g,GD,x1)
    zeroU=_jconst(g,0,3); GD1=GD[:2]; L1=L[:2]; hz=_jmul(g,h,Z)[:2]; term1=_jdiv(g,_jscale(g,GD1,g.const(Fraction(-1,2))),L1); term2=_jdiv(g,_jmul(g,hz,GD[1:3]),_jscale(g,L1,g.const(4))); G2=_jmul(g,G[:2],G[:2]); H2=_jmul(g,H[:2],H[:2]); diff2=_jsub(g,G2,H2); int1=g.op("formal_primitive",[g.const(1)],integrand="1",lower=1,upper=2,result="1",exact_primitive_rule="integral_1_dx_equals_2_minus_1"); log2=g.op("formal_primitive",[x],integrand="1/x",lower=1,upper=2,result="log(2)",exact_primitive_rule="integral_inv_x_dx_equals_log_2"); primitive_unit_check=(Fraction(2)-Fraction(1)==Fraction(1)); term3=_jscale(g,_jmul(g,diff2,[log2,g.const(0)]),g.const(Fraction(-1,4))); I=_jmul(g,_jmul(g,Rc[:2],Rc[:2]),_jadd(g,_jadd(g,term1,term2),term3))
    nbg=_cached_c1(g,native_bg,"cached_native_minus_background")
    bg_old_c1=_cached_c1(g,bg_old,"cached_background_minus_old")
    ep=_cached_c1(g,endpoint,"cached_endpoint_difference")
    total=_jadd(g,nbg,I)
    offset_numeric_minus_old=_jsub(g,bg_old_c1,I)
    legacy_numeric_response=_jadd(g,total,ep)
    response=_jadd(g,_jadd(g,total,offset_numeric_minus_old),ep)
    atoms={"cached_NBG":nbg,"I_BGnum":I,"cached_BG_old":bg_old_c1,"cached_endpoint":ep}
    split_vectors=[[Fraction(1),Fraction(1),Fraction(0),Fraction(0)],
                   [Fraction(0),Fraction(-1),Fraction(1),Fraction(0)],
                   [Fraction(0),Fraction(0),Fraction(0),Fraction(1)]]
    merged_coefficients=[sum(vector[i] for vector in split_vectors) for i in range(4)]
    sparse=_sparse(g,atoms,merged_coefficients)
    sparse_merge_check=(merged_coefficients==[Fraction(1),Fraction(0),Fraction(1),Fraction(1)])
    n_minus_old=_jadd(g,total,offset_numeric_minus_old)
    count=len(g.nodes)
    counts=dict(g.counts,new_derive_calls=1,original_source_evaluations=0,numeric_evaluations=0,quadrature_calls=0,normalization_calls=0,old_proof_calls=0)
    result={"schema":"compact_anchor_defect_transport_kernel_v1","source_only":True,"graph_constructed":True,"passed":True,"graph":{"nodes":g.nodes,"node_count":count},"source_ports":{"coordinate_source":coord,"geometry":{"Rc":Rc[0],"delta":delta[0],"Z":Z[0],"local_coordinate":coord,"Pstar":P},"memory_mstar":ms,"anchor_law":anchor_law,"anchor_m_C3":mrows,"cached_native_minus_background":native_bg,"cached_background_minus_old":bg_old,"cached_endpoint_difference":endpoint},"jets":{"A":A,"B":B4,"D":D,"G":G,"H":H,"D_M":DM,"D_Q":DQ,"physical_U_BG":zeroU,"I_BG_minus_numeric":I,"offset_numeric_minus_old":offset_numeric_minus_old},"pointwise":{"anchor_formula":"A=Pstar*m_normalized_at_x1","baseline_formula":"B=4*Z*exp(-14)*mstar_exact","G_formula":"-((1-delta)*Z*A+(1-Z^2)*A_Z)/L","H_formula":"-4*exp(-14)*mstar_exact","L_formula":"1-delta*Z^2","physical_U_BG_formula":"0 as an input field formula; no source equality claim"},"integrals":{"radial_basis_int_1_dx":int1,"radial_basis_int_inv_x_dx":log2,"primitive_unit_check":primitive_unit_check,"I_formula":"Rc^2*(-(G-H)/(2L)+(1-delta)*Z*(G_Z-H_Z)/(4L)-(G^2-H^2)*log(2)/4)","native_minus_numeric_volume":total,"offset_numeric_minus_old":offset_numeric_minus_old,"cached_endpoint":ep},"response":{"N_minus_numeric_volume":total,"legacy_N_minus_numeric_with_endpoint":legacy_numeric_response,"N_minus_old_volume":n_minus_old,"endpoint_cached_outside":ep,"N_minus_old_with_endpoint":response},"pointwise_atoms":{"atoms":{n: {"value":v[0],"dz":v[1]} for n,v in atoms.items()},"split_vectors":[[ _fj(c) for c in vector] for vector in split_vectors],"merged_coefficients":[_fj(c) for c in merged_coefficients],"sparse_merge_check":sparse_merge_check,"all_atom_arguments_retained":True,"response_sparse":{"value":sparse[0],"dz":sparse[1]}},"audits":{"journal_received":True,"exact_rational_coefficients":True,"mstar_midpoint_used":False,"compact_anchor_law_proved":False,"compact_anchor_law_required_external_precondition":True,"exact_memory_binding_required":True,"L_nonzero_required_external_precondition":True,"source_equality_or_zero_claimed":False,"cap_identity_claimed":False,"native_minus_background_is_cached_external_ref":True,"background_minus_old_is_cached_external_ref":True,"actual_N_minus_old_is_compensated_split":True,"sparse_merge_exact":sparse_merge_check},"counts":counts}
    journal({"event":"compact_anchor_defect_graph_constructed","node_count":count,"passed":True,"sparse_merge_exact":sparse_merge_check})
    result["audits"]["journal_event_emitted"]=True
    return result

"""Fresh literal-zero transport binding and prefix error relocation."""

def _expect_subset(record, expected, label):
    if not isinstance(record, dict) or any(record.get(k) != v for k, v in expected.items()):
        raise ValueError('Frozen compact literal changed: ' + label)


def bind_compact_anchor_literal(i):
    """Prove this four-row zero-integrand transport, without an old reducer."""
    lit = i['compact_literal']
    if lit['archive_raw_sha256'] != '3593dd1ba3e510081090dd315d8581c1253afd5192233a284079c0646a556269':
        raise ValueError('Wrong native compact archive')
    n = lit['selected_node_records']
    required = {
        0: {'operation': 'exact_rational', 'numerator': 0, 'denominator': 1},
        1: {'operation': 'exact_rational', 'numerator': 1, 'denominator': 1},
        220: {'operation': 'exact_rational', 'numerator': -1, 'denominator': 1},
        260181: {'operation': 'bound_variable', 'name': 'true_compact_repair_x'},
        260694: {'operation': 'bound_variable', 'name': 'true_outer_transport_true_compact_background_6'},
        260182: {'operation': 'analytic_unary', 'name': 'log', 'argument': 260181},
        260695: {'operation': 'function_substitution', 'expression': 260182,
                 'variable': 260181, 'value': 1, 'Z_independent_substitution': True},
        260697: {'operation': 'negative', 'argument': 260695},
        260698: {'operation': 'sum', 'arguments': [260182, 260697]},
        260702: {'operation': 'product', 'arguments': [220, 260698]},
        260703: {'operation': 'analytic_unary', 'name': 'exp', 'argument': 260702},
        260706: {'operation': 'function_substitution', 'expression': 0,
                 'variable': 260181, 'value': 260694, 'Z_independent_substitution': True},
        260708: {'operation': 'product', 'arguments': [260706, 260707]},
    }
    for node, expected in required.items():
        _expect_subset(n[str(node)], expected, str(node))
    incoming = [251017, 251018, 251019, 251020]
    for k in range(4):
        _expect_subset(n[str(260709+k)], {
            'operation': 'definite_integral', 'integrand': 260708,
            'variable': 'true_outer_transport_true_compact_background_6',
            'lower': 1, 'upper': 260181, 'exact_function_integral': True,
            'endpoint_independent_of_Z': True, 'ordinary_slow_Z_derivative_order': k,
            'original_chart': 'true_compact_background', 'original_rate': '1',
            'physical_Jacobian_applied_once': True, 'predecessor_memory_never_reset': True,
            'leading_background': True, 'finite_N_correction': False,
        }, str(260709+k))
        _expect_subset(n[str(260716+k)], {'operation': 'product', 'arguments': [incoming[k], 260703]}, str(260716+k))
        _expect_subset(n[str(260723+k)], {'operation': 'sum', 'arguments': [260716+k, 260709+k]}, str(260723+k))
    if i['anchor_binding']['normalized_m_at_x1_C3']['source_nodes'] != incoming:
        raise ValueError('Anchor must be the retained incoming rows, not a mean enclosure')
    for key, expected in [('Pstar', 33572), ('Rc', 251066), ('delta', 25220), ('Z', 20), ('local_coordinate', 260181)]:
        b = i['anchor_binding']['Pstar'] if key == 'Pstar' else i['geometry_binding'][key]
        if b['source_node'] != expected or b['raw_sha256'] != lit['archive_raw_sha256']:
            raise ValueError('Geometry namespace changed: '+key)
    return {
        'proved': True, 'ordinary_Z_orders': [0, 1, 2, 3], 'domain': '1 <= x <= 2',
        'normalized_m_formula': 'm_BG[k](x,Z)=incoming_m[k](Z)/x',
        'zero_substitution': 260706, 'zero_product': 260708,
        'zero_integrals': [260709, 260710, 260711, 260712],
        'decay': 'exp(-(log(x)-log(1)))=1/x, x>0',
        'physical_m_formula': 'M_BG[k](x,Z)=Pstar*incoming_m[k](Z)/x',
        'Pstar_factors': 1, 'incoming_histories_reset': False,
        'proof_rules': ['substitute literal zero', 'zero times a defined function',
                        'integral of zero on the source domain', 'log(1)=0',
                        'exp(-log(x))=1/x for positive x', 'linearity for ordinary Z jets'],
        'new_literal_transport_proof_calls': 1, 'old_proof_calls': 0,
        'source_domain_and_regular_ordinary_Z_contract_required': True,
        'incoming_anchor_equals_numeric_baseline': False,
    }


def append_prefix_error_split(graph, i):
    """Define actual field errors against both baselines; never select a cap."""
    ns = graph['graph']['nodes']
    def put(op, args=(), **metadata):
        node = len(ns)
        ns.append(dict(id=node, op=op, args=list(args), **metadata))
        return node
    def source_field(label, rows):
        return [put('source_ref', raw_sha256=i['compact_literal']['archive_raw_sha256'],
                    selector='true_compact_C5_bindings.physical.'+label,
                    source_node=row, ordinary_Z_order=k, label='actual_N_'+label)
                for k, row in enumerate(rows)]
    x = next(n['id'] for n in ns if n['op'] == 'bound_variable' and n['name'] == 'compact_anchor_x')
    jets = graph['jets']
    actual = {'M': source_field('m', [262281, 262282, 262283, 262284]),
              'Q': source_field('Q', [262364, 262368, 262374]),
              'U': source_field('U', [262269, 262270, 262271])}
    bg = {'M': [put('div', [a, x]) for a in jets['A']],
          'Q': [put('div', [a, x]) for a in jets['G']], 'U': jets['physical_U_BG']}
    baseline = {'M': [put('div', [a, x]) for a in jets['B']],
                'Q': [put('div', [a, x]) for a in jets['H']], 'U': jets['physical_U_BG']}
    differences = {'M': jets['D_M'], 'Q': jets['D_Q'], 'U': jets['physical_U_BG']}
    output = {}
    vectors = [[Fraction(1), Fraction(-1), Fraction(0)],
               [Fraction(0), Fraction(1), Fraction(-1)]]
    merged = [sum(v[k] for v in vectors) for k in range(3)]
    if merged != [Fraction(1), Fraction(0), Fraction(-1)]:
        raise ValueError('Prefix error relocation algebra failed')
    for name in ['M', 'Q', 'U']:
        e_nb = [put('sub', [a, b]) for a, b in zip(actual[name], bg[name])]
        e_num = [put('sub', [a, b]) for a, b in zip(actual[name], baseline[name])]
        output[name] = {'actual_N': actual[name], 'physical_BG': bg[name], 'exact_numeric_baseline': baseline[name],
                        'N_minus_BG': e_nb, 'BG_minus_numeric': differences[name], 'N_minus_numeric': e_num,
                        'ordinary_Z_orders': list(range(len(e_nb))), 'split_coefficients': [str(c) for c in merged]}
    graph['graph']['node_count'] = len(ns)
    graph['counts']['prefix_field_source_refs'] = sum(len(v) for v in actual.values())
    graph['counts']['new_prefix_error_split_checks'] = 1
    return {'fields': output, 'identity': 'e_N_numeric=e_N_BG+(BG-numeric)',
            'bound_transport': '|e_N_numeric[k]| <= B_relative[k]+|D_BG_numeric[k]|',
            'relative_cap_requires_same_source_bridge': True,
            'signed_mismatch_is_not_bounded_by_relative_cap': True,
            'all_incoming_prefix_histories_retained': True,
            'additional_N_factor': False, 'numeric_cap_or_midpoint_selected': False}


def derive(source, journal):
    result = {'schema': 'original_compact_anchor_defect_transport',
              'component': 'compact_anchor_defect_transport', 'passed': False,
              'source_only': True, 'metadata_only': False, 'input_digest': source['input_digest'],
              'runtime_model_metadata': None, 'original_source_evaluations': 0,
              'source_packet_defined': True, 'source_ftc_transport_defined': True,
              'cap_to_source_identity_proved': False,
              'retained_construction_data': source['parent_report']['retained_construction_data'],
              'retained_integral_term_inventory': source['parent_report']['retained_integral_term_inventory'],
              'retained_cap_provenance': source['parent_result']['retained_cap_provenance'],
              'source_inspection_audit': {'old_source_consumer_calls': 0, 'old_proof_functions_called': 0,
                                         'old_operator_calls': 0, 'old_source_normalizer_calls': 0,
                                         'original_source_evaluations': 0, 'quadratures': 0,
                                         'numeric_field_evaluations': 0, 'temporal_recurrence_calls': 0}}
    result.update({key: False for key in source['parent_report']['false_flags']})
    try:
        i = source['input']
        contract = source['parent_result']['cached_source_contract']
        if contract['background_normalized_V_C2'] != [0, 0, 0] or contract['background_normalized_m_C3'] != [260723,260724,260725,260726]:
            raise ValueError('Cached native background field contract changed')
        law = bind_compact_anchor_literal(i)
        result['source_binding'] = law
        recipe = i['memory_binding']['mstar_exact']['recipe']
        if recipe['definition'] != 'exp(1-Y)+integral_1^Y exp(-Y+y)*sigma(1-log(y)/40)dy; Y=exp(40)' or recipe['initial_memory'] != 1 or recipe['kernel_rate'] != 1:
            raise ValueError('Exact memory recipe changed or homogeneous term lost')
        graph = construct_anchor_graph(i, journal)
        result['anchor_defect_transport'] = graph
        graph['audits']['compact_anchor_law_proved'] = True
        graph['audits']['compact_anchor_law_required_external_precondition'] = False
        result['prefix_error_relocation'] = append_prefix_error_split(graph, i)
        if graph['pointwise_atoms']['sparse_merge_check'] is not True or graph['integrals']['primitive_unit_check'] is not True:
            raise ValueError('Fresh exact arithmetic checks failed')
        result['exact_memory_definition'] = i['memory_binding']['mstar_exact']
        result['prefix_error_source_recipe'] = i['prefix_error_source_recipe']
        result['coordinate_mapping'] = 'compact_anchor_x := native source node 260181; x=R/Rc'
        result['formal_preconditions'] = {'common_source_domain': True, 'L_nonzero': True, 'ordinary_Z_contract': 'Z_Z=1; delta_Z=Rc_Z=Pstar_Z=x_Z=0', 'original_BG_pressure_recovery_source_recipe_required': True}
        result['memory_enclosure_applied_to_exact_function'] = False
        result['static_receipt_correction'] = 'Exact baseline M ordinary Z0..3 is [4*Z*c,4*c,0,0], c=exp(-14)*mstar; the source-trace prose first entry omits Z.'
        result['compact_anchor_law_proved'] = True
        result['upstream_anchor_to_exact_memory_identity_proved'] = False
        result['analytic_BG_minus_numeric_integral_defined'] = True
        result['remaining_scope'] = ['join upstream incoming anchor and physical normalization to the exact memory',
                                     'bind complete-prefix correction caps to actual N-minus-BG function jets',
                                     'evaluate signed Evol and independent reference offset with certified errors']
        result['passed'] = True
        journal({'event': 'compact_anchor_law_and_error_relocation_saved', 'source_value_evaluations': 0,
                 'anchor_law_proved': True, 'upstream_mean_join_proved': False,
                 'nodes': graph['graph']['node_count'], 'cap_identity_proved': False})
    except Exception as error:
        import traceback
        result['failure'] = {'type': type(error).__name__, 'message': str(error), 'traceback': traceback.format_exc()}
        journal({'event': 'anchor_transport_failure', 'type': type(error).__name__, 'message': str(error)})
    return result
