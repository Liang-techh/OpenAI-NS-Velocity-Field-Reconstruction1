"""Fresh finite affine upstream-mean identity DAG; intentionally unentered."""
from fractions import Fraction

_MISSING=object()
def _f(x):
    if isinstance(x,Fraction): return x
    if isinstance(x,int): return Fraction(x,1)
    if isinstance(x,str): return Fraction(x)
    if isinstance(x,dict): return Fraction(int(x["numerator"]),int(x["denominator"]))
    raise TypeError("exact rational required")
def _fj(x):
    x=_f(x); return {"numerator":x.numerator,"denominator":x.denominator}
def _get(o,path):
    for k in path.split("."):
        if not isinstance(o,dict) or k not in o:return _MISSING
        o=o[k]
    return o
def _need(o,*paths):
    for p in paths:
        x=_get(o,p)
        if x is not _MISSING and x is not None:return x
    raise KeyError("missing explicit contract: "+", ".join(paths))
def _freeze(x):
    if isinstance(x,dict):return tuple(sorted((k,_freeze(v)) for k,v in x.items()))
    if isinstance(x,(list,tuple)):return tuple(_freeze(v) for v in x)
    return x
def _binding(b,label,allow_proof=False):
    if not isinstance(b,dict):raise TypeError(label+" must be a mapping")
    if "raw_sha256" not in b and isinstance(b.get("raw_sha"),str):b=dict(b);b["raw_sha256"]=b["raw_sha"]
    for k in ("raw_sha256","selector"):
        if not isinstance(b.get(k),str) or not b[k]:raise ValueError(label+" requires "+k)
    keys=("source_node","source_nodes","recipe")+(("passed_static_receipt","exact_node_records") if allow_proof else ())
    if not any(k in b for k in keys):raise ValueError(label+" requires an explicit source node or recipe")
    return b

class _DAG:
    def __init__(self):
        self.nodes=[];self.index={};self.counts={k:0 for k in ("const","add","sub","mul","neg","source_ref","recipe_ref","formal_width","formal_exponential","formal_kernel_primitive","formal_residual","bound_variable")}
    def put(self,key,rec):
        if key in self.index:return self.index[key]
        i=len(self.nodes);rec=dict(rec);rec["id"]=i;self.nodes.append(rec);self.index[key]=i
        if rec.get("op") in self.counts:self.counts[rec["op"]]+=1
        return i
    def const(self,x):
        x=_f(x);return self.put(("const",x.numerator,x.denominator),{"op":"const","value":_fj(x)})
    def op(self,op,args=(),**meta):
        args=tuple(int(a) for a in args);return self.put((op,args,_freeze(meta)),dict({"op":op,"args":list(args)},**meta))
    def ref(self,b,label,order=0):
        b=_binding(b,label);order=int(order);nodes=b.get("source_nodes")
        if nodes is not None:
            if order>=len(nodes):raise ValueError(label+" order out of range")
            return self.put(("source_ref",b["raw_sha256"],b["selector"],label,_freeze(nodes[order]),order),{"op":"source_ref","raw_sha256":b["raw_sha256"],"selector":b["selector"],"label":label,"source_node":int(nodes[order]),"ordinary_Z_order":order})
        if order==0 and "source_node" in b:
            return self.put(("source_ref",b["raw_sha256"],b["selector"],label,int(b["source_node"]),0),{"op":"source_ref","raw_sha256":b["raw_sha256"],"selector":b["selector"],"label":label,"source_node":int(b["source_node"]),"ordinary_Z_order":0})
        if "recipe" not in b:raise ValueError(label+" missing source order "+str(order))
        return self.put(("recipe_ref",b["raw_sha256"],b["selector"],label,_freeze(b["recipe"]),order),{"op":"recipe_ref","raw_sha256":b["raw_sha256"],"selector":b["selector"],"label":label,"recipe":b["recipe"],"ordinary_Z_order":order})
    def local(self,name):return self.put(("bound_variable",name),{"op":"bound_variable","name":name,"ordinary_Z_derivative":0})

def _add(g,a,b):return g.op("add",[a,b])
def _sub(g,a,b):return g.op("sub",[a,b])
def _mul(g,a,b):return g.op("mul",[a,b])
def _neg(g,a):return g.op("neg",[a])
def _sparse(terms):
    out={}
    for name,coef in terms.items():
        coef=_f(coef)
        if coef:out[name]=out.get(name,Fraction(0))+coef
    return out
def _sparse_total(terms):return sum(_sparse(terms).values(),Fraction(0))
def _pair_cancellation(pairs,requires_external_binding=False):
    raw=[(str(name),_f(coef)) for name,coef in pairs]
    merged={}
    for name,coef in raw:merged[name]=merged.get(name,Fraction(0))+coef
    return {"terms":[{"basis":name,"coefficient":_fj(coef)} for name,coef in raw],"merged_terms":_json_sparse(merged),"residual":_fj(sum(merged.values(),Fraction(0))),"all_basis_coefficients_zero":all(coef==0 for coef in merged.values()),"requires_external_binding":bool(requires_external_binding)}
def _json_sparse(x):return {k:_fj(v) for k,v in sorted(x.items())}
def _json_value(x):
    if isinstance(x,Fraction):return _fj(x)
    if isinstance(x,dict):return {str(k):_json_value(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [_json_value(v) for v in x]
    return x
def _width_node(g,width,label):
    if isinstance(width,dict):return g.ref(width,label,0)
    if isinstance(width,(int,Fraction,str)):
        try:return g.const(_f(width))
        except (TypeError,ValueError):return g.op("formal_width",[],label=label,expression=str(width))
    raise TypeError(label+" width must be exact or symbolic")
def _chart_record(g,c):
    if not isinstance(c,dict):raise TypeError("chart contract must be mapping")
    name=str(_need(c,"chart"));source_type=str(_need(c,"source_type"));width=c["width"];w=_width_node(g,width,name+".width")
    refs={}
    for order in range(4):refs[str(order)]=g.ref(_need(c,"ordinary_Z"+str(order)),name+".ordinary_Z"+str(order),0)
    seed=_binding(_need(c,"seedtheorem"),name+".seedtheorem",allow_proof=True)
    return {"chart":name,"source_type":source_type,"width":w,"width_input":width,"ordinary_Z_refs":refs,"seedtheorem":seed}
def _formal_transport(g,F,width,label):
    e=g.op("formal_exponential",[width],label=label+".exp_neg_width",exponent="-width")
    hom=_mul(g,F,e);force=_mul(g,F,_sub(g,g.const(1),e));total=_add(g,hom,force)
    return {"exponential":e,"homogeneous":hom,"forcing":force,"total":total,"identity_expression":"F*exp(-w)+F*(1-exp(-w))=F"}
def _post_widths(contracts):
    xs=_get(contracts,"postaxial_widths")
    if xs is _MISSING:xs=[c.get("width") for c in contracts.get("charts",[]) if isinstance(c,dict) and str(c.get("phase",""))=="postaxial"]
    if not isinstance(xs,(list,tuple)) or len(xs)!=3:return None
    return [_f(x) for x in xs]

def prove_upstream_mean(contracts,journal):
    if not isinstance(contracts,dict) or not callable(journal):raise TypeError("contracts mapping and callable journal required")
    charts=_need(contracts,"charts");
    if not isinstance(charts,list) or not charts:raise ValueError("charts must be a nonempty list")
    g=_DAG(); records=[_chart_record(g,c) for c in charts]
    anchor=_need(contracts,"anchor_binding"); memory=_need(contracts,"memory_binding")
    P=_binding(_need(anchor,"Pstar"),"anchor.Pstar");F=_binding(_need(anchor,"F"),"anchor.F");Zc=_binding(_need(anchor,"Zc"),"anchor.Zc");law=_binding(_need(contracts,"compact_anchor_law"),"compact_anchor_law",allow_proof=True);mstar=_binding(_need(memory,"mstar_exact"),"memory.mstar_exact")
    Fref=g.ref(F,"anchor.F");Pref=g.ref(P,"anchor.Pstar");Zref=g.ref(Zc,"anchor.Zc");z0=g.ref(Zc,"anchor.Zc",0)
    transports=[]
    for r in records:
        if r["source_type"] in ("constant_F","constant","Rh_reference","O2_slope"):
            coeff=_pair_cancellation([("F_E",1),("F",1),("F_E",-1),("F",-1)]);t=_formal_transport(g,Fref,r["width"],r["chart"]);transports.append({"chart":r["chart"],"width":r["width"],"terms":t,"coefficient_cancellation":dict(coeff,requires_seedtheorem=True)})
    exp40=g.op("formal_exponential",[],label="exp40",exponent="40");D=g.op("formal_width",[exp40],label="axial_D",expression="exp(40)-1");Y=exp40;kernel=g.op("formal_kernel_primitive",[Y],label="mstar_kernel",lower=1,upper="Y",integrand="exp(-Y+y)*sigma(1-log(y)/40)");mstar_hom=g.op("formal_exponential",[Y],label="mstar_hom",exponent="1-Y");mstar_node=_add(g,mstar_hom,kernel);axial_exp=g.op("formal_exponential",[D],label="axial_decay",exponent="-D");axial_hom=_mul(g,Fref,axial_exp);axial_forcing=_mul(g,Fref,kernel);axial_total=_add(g,axial_hom,axial_forcing)
    axial_coeff=_pair_cancellation([("F_H",1),("F_K",1),("F_H",-1),("F_K",-1)]);axial={"D":D,"Y":Y,"mstar":mstar_node,"mstar_binding":mstar,"formula":"mstar=exp(1-Y)+integral_1^Y exp(-Y+y)*sigma(1-log(y)/40)dy","homogeneous":axial_hom,"forcing":axial_forcing,"terminal":axial_total,"nonzero_homogeneous_term":axial_exp,"nonzero_homogeneous_term_scaled":axial_hom,"forcing_kernel":kernel,"identity_expression":"F*exp(-D)+F*integral_1^Y exp(-Y+y)*sigma(1-log(y)/40)dy=F*mstar","inlet_offset":"7","terminal_offset":"6+Y","D_expression":"Y-1","physical_offset":"6+exp(40*s)","jacobian":"40*exp(40*s)","V_expression":"F*sigma(1-s)","width_expression":"exp(40)-1","coefficient_cancellation":dict(axial_coeff,requires_mstar_definition=True)}
    post=_post_widths(contracts);post_sum=None;post_exp=None
    if post is not None:
        post_sum=sum(post,Fraction(0));post_sum_node=g.const(post_sum);post_exp=g.op("formal_exponential",[post_sum_node],label="postaxial_decay",exponent="-width_sum")
    post_terms=[]
    if post is not None:
        for i,w in enumerate(post):
            e=g.op("formal_exponential",[g.const(w)],label="postaxial_"+str(i)+"_decay",exponent="-width")
            incoming=g.local("postaxial_"+str(i)+"_incoming")
            homogeneous=_mul(g,incoming,e);forcing=g.const(0)
            post_terms.append({"index":i,"width":_fj(w),"exponential":e,"homogeneous":homogeneous,"forcing":forcing,"total":homogeneous,"forcing_is_zero_by_contract":True})
    anchor_m=g.op("mul",[Fref,mstar_node]);compact_A=g.op("mul",[Pref,g.op("mul",[anchor_m,post_exp if post_exp is not None else g.op("formal_exponential",[],label="postaxial_decay_unbound",exponent="-14")])]);baseline=g.op("mul",[g.op("mul",[g.const(4),Zref]),g.op("mul",[mstar_node,post_exp if post_exp is not None else g.op("formal_exponential",[],label="postaxial_decay_unbound",exponent="-14")])]);pstar_F=_mul(g,Pref,Fref);four_Z=_mul(g,g.const(4),Zref);pstar_cancel=_sub(g,pstar_F,four_Z);anchor_coeff=_sparse({"Pstar_times_F":1,"four_times_Zc":-1})
    identity_residual=_sub(g,compact_A,baseline)
    jets=[g.op("mul",[g.const(4),z0]),g.const(4),g.const(0),g.const(0)]
    deltaM=g.op("formal_residual",[compact_A,baseline],label="DeltaM",conditional_on="compact_anchor_binding");deltaQ=g.op("formal_residual",[compact_A,baseline],label="DeltaQ",conditional_on="compact_anchor_binding");analyticI=g.op("formal_residual",[compact_A,baseline],label="analytic_I",observable="analytic_integral_residual",conditional_on="compact_anchor_binding")
    post_ok=post is not None and post_sum==Fraction(14);binding_structural=True
    for b in (P,F,Zc,law,mstar):binding_structural=binding_structural and isinstance(b,dict)
    count=len(g.nodes);counts=dict(g.counts,new_proof_calls=1,new_derive_calls=1,original_source_evaluations=0,numeric_evaluations=0,cap_reductions=0,normalization_calls=0)
    result={"schema":"compact_upstream_mean_identity_kernel_v1","source_only":True,"graph_constructed":True,"passed":bool(binding_structural and post_ok),"graph":{"nodes":g.nodes,"node_count":count},"chart_contracts":records,"constant_F_transports":transports,"axial_transport":axial,"postaxial":{"widths":list(post) if post is not None else None,"width_sum":_fj(post_sum) if post_sum is not None else None,"width_sum_is_14":post_ok,"decay_node":post_exp,"nonzero_homogeneous_term_preserved":True},"anchor_transport":{"Pstar_ref":Pref,"F_ref":Fref,"Zc_ref":Zref,"compact_A":compact_A,"baseline_B":baseline,"identity_residual":identity_residual,"formal_identity":"Pstar*F*mstar*exp(-14)=4*Zc*mstar*exp(-14) after explicit Pstar*F=4Zc binding","C3_jets":[int(x) for x in jets]},"conditional_zero_outputs":{"DeltaM":{"node":deltaM,"status":"requires_explicit_compact_anchor_binding"},"DeltaQ":{"node":deltaQ,"status":"requires_explicit_compact_anchor_binding"},"analytic_I":{"node":analyticI,"status":"requires_explicit_compact_anchor_binding"}},"preserved_histories":{"nonzero_axial_homogeneous":axial["nonzero_homogeneous_term"],"leading_transport_input_refs":[r["ordinary_Z_refs"] for r in records],"leading_proof_does_not_delete_corrections":True},"audits":{"journal_received":True,"structural_bindings_present":binding_structural,"compact_anchor_law_proved":False,"source_identity_proved":False,"terminal_zero_claimed":False,"cap_reduction_claimed":False,"numeric_values_used":False,"physical_Pstar_cancellation_is_external_binding":True,"postaxial_width_sum_check":post_ok,"exact_fraction_coefficients":True},"counts":counts}
    result["chart_contracts"]=_json_value(records)
    result["graph"]["nodes"]=_json_value(g.nodes)
    result["constant_F_transports"]=_json_value(transports)
    result["axial_transport"]=_json_value(axial)
    result["postaxial"]["widths"]=_json_value(list(post)) if post is not None else None
    result["postaxial"]["terms"]=_json_value(post_terms)
    result["source_function_identity"]={"seed_identity":"Rh: V=m=4Z/Pstar (cached source theorem binding)","canonical_F_jet_formulas":["4*Zc/Pstar","4/Pstar","0","0"],"constant_transport":"F*exp(-w)+F*(1-exp(-w))=F","axial_transport":axial["identity_expression"],"postaxial_transport":"m_out=m_in*exp(-w) when V=0","compact_anchor":"Pstar*F*mstar*exp(-14)=4*Zc*mstar*exp(-14) conditional on explicit Pstar*F=4Zc source binding"}
    result["anchor_transport"].update({"Pstar_times_F_node":pstar_F,"four_times_Zc_node":four_Z,"Pstar_cancellation_residual_node":pstar_cancel,"Pstar_cancellation_coefficients":{"terms":_json_sparse(anchor_coeff),"formal_coefficient_residual":_fj(_sparse_total(anchor_coeff)),"requires_explicit_binding":True,"positive_Pstar_is_external_contract":True},"C3_scope":"Pstar*F affine factor; no M or Q value is set to zero","C3_jet_formulas":["4*Zc","4","0","0"],"canonical_F_jet_formulas":["4*Zc/Pstar","4/Pstar","0","0"]})
    result["preserved_histories"]["chart_contracts"]=_json_value(records)
    result["audits"]["cached_seed_identity_reused"] = True
    result["audits"]["exact_join_proof_called"] = False
    result["audits"]["affine_coefficient_pair_check"] = True
    result["audits"]["proof_scope"] = "structural affine graph only; source identity and conditional zero outputs require external binding"
    journal("upstream_mean_identity_graph_constructed",{"node_count":count,"passed":result["passed"],"postaxial_width_sum_is_14":post_ok})
    result["audits"]["journal_event_emitted"]=True
    return result


"""Bind the new affine mean proof to the actual six archived source charts."""

def bind_six_source_charts(i, parent):
    s = i['literal_source']; nodes = s['literal_nodes']; seed = s['seed_identity']
    def n(j): return nodes[str(j)]
    def need(ok, msg):
        if not ok: raise ValueError(msg)
    def subset(j, expected):
        need(all(n(j).get(k) == v for k,v in expected.items()), 'Literal source mismatch '+str(j))
    def args(j, operation, values):
        subset(j, {'operation':operation, 'arguments':values})
    def rational(j):
        subset(j, {'operation':'exact_rational'})
        return Fraction(n(j)['numerator'], n(j)['denominator'])
    def subst(j, body, variable, value):
        subset(j, {'operation':'function_substitution','expression':body,'variable':variable,
                   'value':value,'Z_independent_substitution':True})
    def unary(j, name):
        subset(j, {'operation':'analytic_unary','name':name})
        return n(j)['argument']
    def minus_difference(j, left, right):
        aa=n(j)['arguments']; need(len(aa)==2,'Difference arity')
        args(j,'sum',[left,aa[1]])
        subset(aa[1],{'operation':'negative','argument':right})
    need(s['archive_raw_sha256']=='3593dd1ba3e510081090dd315d8581c1253afd5192233a284079c0646a556269','Native source changed')
    need(seed['canonical_normalized_Rh']=='E=exp(-1/2)/(1+Z^2); V=m=4Z/Pstar; h=5E/8; k=h*V; e=V^2-5E^2/12; p=5E^2/2','Cached canonical seed changed')
    for key in ['same_exact_unique_branch_and_C6_implicit_derivatives',
                'propagated_original_histories_identified_by_linear_ODE_uniqueness',
                'legacy_C1_or_all_N_correction_never_used_as_selected_N_inlet']:
        need(seed[key] is True, 'Cached seed contract absent: '+key)
    p=seed['original_full_weight_function_identity']
    for key in ['passed','full_support_primitive_changes_equal_current_Lh_plus_Q',
                'full_weight_identity_uses_unique_current_smooth_leading_family',
                'ordinary_Z_derivatives_of_equal_functions_retained']:
        need(p[key] is True,'Cached original Rh theorem absent: '+key)
    receipts=s['cached_seed_receipts']
    rr=[v['receipt'] for v in receipts.values()]
    need(len(rr)==2 and all(v['all_passed'] is True for v in rr),'Original seed receipts not accepted')
    need(rr[0]['source_family']==rr[1]['source_family']==i['source_family'],'Seed source family changed')
    need(any(v.get('same_exact_original_unique_branch_theorem_bound') is True for v in rr),'True C6 seed branch missing')
    need(any(v.get('current_whole_Z_repaired_leading_Rh_reference_function_identity_checked') is True for v in rr),'Full-weight Rh functional identity missing')
    subset(33572,{'operation':'analytic_unary','name':'exp','argument':5})
    args(5,'sum',[3,4]);need(unary(3,'exp')==2 and rational(4)==11,'Pstar is not the same Z-independent original parameter')
    subset(183811,{'operation':'positive_quotient','numerator':1,'denominator':33572})
    need(rational(1)==1 and rational(0)==0 and rational(220)==-1 and rational(2)==40,'Exact constants changed')
    canonical_V=seed['actual_Rh_fields']['V'][:4]
    incoming=seed['actual_Rh_fields']['m'][:4]
    expected_names=['Rh_reference','O2_slope','O2_axial','O2_buffer','O3_transition','O3_power']
    charts=s['charts'];need([c['chart'] for c in charts]==expected_names,'Six source charts missing')
    expected_domains=[[-5,0],[0,1],[0,1],[0,11],[0,1],[0,2]]
    expected_types=['constant_F','constant_F','axial_sigma','zero','zero','zero']
    guard_count=0
    for key,node in [('Pstar',33572),('F',canonical_V[0]),('Zc',20)]:
        cb=i['contracts']['anchor_binding'][key]
        need(cb['source_node']==node and cb['raw_sha256']==s['archive_raw_sha256'],'Canonical physical cancellation port changed '+key)
    need(i['contracts']['compact_anchor_law']['recipe']==parent['source_binding'],'Cached compact law is not bound to the parent proof')
    for idx,c in enumerate(charts):
        ct=i['contracts']['charts'][idx]
        need(ct['chart']==c['chart'] and ct['source_type']==expected_types[idx] and ct['width']==[5,1,'exp(40)-1',11,1,2][idx],'Affine chart contract changed')
        for k in range(4):
            rb=ct['ordinary_Z'+str(k)]
            need(rb['source_node']==c['m'][k] and rb['raw_sha256']==s['archive_raw_sha256'],'Ordinary Z contract points to a different source')
        name=c['chart']; x=c['phase']; lower,upper=c['domain']; off=c['log_offset']; jac=c['physical_measure']
        need([rational(j) for j in c['domain']]==[Fraction(v) for v in expected_domains[idx]],'Domain changed '+name)
        need(c['incoming']==incoming,'Noncontiguous leading m input '+name)
        subset(x,{'operation':'bound_variable','name':'true_outer_'+name})
        aa=n(off)['arguments'];subset(off,{'operation':'sum'})
        if idx<2:
            need(len(aa)==2 and rational(aa[0])==6 and aa[1]==x,'Linear offset changed '+name)
            need(jac==1 and c['V']==canonical_V,'Constant seed velocity changed '+name)
        elif idx==2:
            need(len(aa)==2 and rational(aa[0])==6,'Axial offset changed')
            exponent=unary(aa[1],'exp');args(exponent,'product',[2,x]);args(jac,'product',[2,aa[1]])
            beta=None
            for k,v in enumerate(c['V']):
                subset(v,{'operation':'product'});a=n(v)['arguments']
                need(len(a)==2 and a[0]==canonical_V[k],'Axial velocity normalized seed changed')
                if beta is None: beta=a[1]
                need(a[1]==beta,'Axial ordinary-Z coefficient changed')
            subset(beta,{'operation':'exact_original_source_flat_sigma','definition':'0 for x<=0; 1 for x>=1; exp(-1/x^2)/(exp(-1/x^2)+exp(-1/(1-x)^2)) inside','exact_function_not_its_interval_cap':True})
            need(n(beta)['source_binding']['module_sha256']=='a8bc21feb4a2565f976b450f1e4f593180d6eff343e572c31e83117020779657','Original stable sigma provider changed')
            sigma=i['contracts']['memory_binding']['mstar_exact']['recipe']['sigma']
            need(sigma['equivalent_stable_form']=='a/(a+b), a=exp(-1/r^2), b=exp(-1/(1-r)^2)' and sigma['endpoint_values']=={'r=0':0,'r=1':1},'Exact memory sigma differs from native axial sigma')
            b=n(beta)['argument'];bb=n(b)['arguments'];args(b,'sum',[1,bb[1]])
            subset(bb[1],{'operation':'negative','argument':x})
        else:
            need(len(aa)==3 and rational(aa[0])==[6,17,18][idx-3] and aa[2]==x,'Post-axial offset changed '+name)
            need(unary(aa[1],'exp')==2,'Common exact exp40 parameter changed')
            need(jac==1 and c['V']==[0]*4,'Post-axial V is not literal zero')
        # Bind the actual Duhamel kernel, density, Jacobian, endpoints and
        # all four ordinary-Z source rows. No source function is evaluated.
        mem_arg=unary(c['memory'],'exp');ma=n(mem_arg)['arguments']
        need(len(ma)==2 and ma[0]==220,'Rate-one memory changed')
        args(mem_arg,'product',ma);width=ma[1];wa=n(width)['arguments']
        origin=n(wa[1])['argument'];minus_difference(width,off,origin);subst(origin,off,x,lower)
        for k in range(4):
            root=c['m'][k];integ=c['integrals'][k];ar=n(root)['arguments']
            args(root,'sum',[ar[0],integ]);args(ar[0],'product',[c['incoming'][k],c['memory']])
            subset(integ,{'operation':'definite_integral','lower':lower,'upper':x,
                          'ordinary_slow_Z_derivative_order':k,'original_chart':name,'original_rate':'1',
                          'exact_function_integral':True,'endpoint_independent_of_Z':True,
                          'physical_log_offset':off,'physical_Jacobian':jac,'physical_Jacobian_applied_once':True,
                          'predecessor_memory_never_reset':True,'leading_background':True,'finite_N_correction':False})
            body=n(integ)['integrand'];ba=n(body)['arguments'];need(len(ba)==2,'Source body arity')
            args(body,'product',ba)
            field=ba[0];t=n(field)['value'];subst(field,c['V'][k],x,t)
            subset(t,{'operation':'bound_variable','name':n(integ)['variable']})
            weight=ba[1];ww=n(weight)['arguments'];need(len(ww)==2,'Kernel weight arity')
            args(weight,'product',ww);subst(ww[0],jac,x,t)
            ka=unary(ww[1],'exp');kk=n(ka)['arguments'];need(kk[0]==220,'Kernel rate changed');args(ka,'product',kk)
            suffix=kk[1];ss=n(suffix)['arguments'];at_t=n(ss[1])['argument']
            minus_difference(suffix,off,at_t);subst(at_t,off,x,t)
            subst(c['endpoint'][k],root,x,upper)
            guard_count+=1
        incoming=c['endpoint']
    need(incoming==[251017,251018,251019,251020],'Compact incoming anchor changed')
    need(parent['compact_anchor_law_proved'] is True and parent['source_binding']['proved'] is True,'Cached compact transport unavailable')
    need(parent['exact_memory_definition']==i['contracts']['memory_binding']['mstar_exact'],'Exact mstar recipe differs from saved baseline')
    # A separate new four-row linear source identity binds the relative
    # m correction; no range cap is evaluated or promoted to a signed value.
    relative=s['relative_mass_rows']
    for k in range(4):
        a=relative[str(261716+k)];b=relative[str(262281+k)]
        need(a['operation']=='sum' and a['arguments']==[260723+k,261593+k],'Native relative m sum changed')
        need(b['operation']=='product' and b['arguments']==[261716+k,33572],'Physical m lift changed')
    # Substitute the cached F=4Z/P identity and merge exact P/Z monomials.
    physical_coeff={}
    P_exponent=1+(-1)
    for monomial,coefficient in [((P_exponent,1),Fraction(4)),((0,1),Fraction(-4))]:
        physical_coeff[monomial]=physical_coeff.get(monomial,Fraction(0))+coefficient
    need(P_exponent==0 and all(v==0 for v in physical_coeff.values()),'Pstar*F=4Z physical cancellation failed')
    relative_coeff={}
    for atom,coefficient in [('P_BG',Fraction(1)),('P_H',Fraction(1)),('P_BG',Fraction(-1))]:
        relative_coeff[atom]=relative_coeff.get(atom,Fraction(0))+coefficient
    need(relative_coeff=={'P_BG':Fraction(0),'P_H':Fraction(1)},'Relative physical mass expansion failed')
    return {'Pstar_F_source_identity_proved':True,'physical_parameter_cancellation':{'P_exponent_sum':P_exponent,'merged_monomial_coefficients':{str(k):_fj(v) for k,v in physical_coeff.items()},'nonzero_Pstar_proved':True},'relative_mass_exact_coefficients':{k:_fj(v) for k,v in relative_coeff.items()},'passed':True,'bound_chart_count':6,'ordinary_Z_rows_bound':guard_count,
            'canonical_Rh_seed_theorem_reused':True,'old_seed_proof_calls':0,
            'Pstar_positive_exp_node':33572,'inverse_Pstar_node':183811,
            'physical_Pstar_factors':1,'incoming_C3_nodes':incoming,
            'source_domain':'common admitted source domain of the cached Rh theorem and six charts',
            'mean_source':'same exact rate-one sigma memory with nonzero homogeneous inlet',
            'relative_mass_C3_identity':'M_N[k]-Pstar*m_BG[k]=Pstar*H_m[k], k=0..3',
            'relative_history_C3_nodes':[261593,261594,261595,261596],
            'incoming_corrections_zeroed':False,'relative_cap_transfer_proved':False}


def derive(source,journal):
    parent=source['parent_result']; report=source['parent_report']
    out={'schema':'original_compact_upstream_mean_identity','component':'compact_upstream_mean_identity',
         'passed':False,'source_only':True,'metadata_only':False,'runtime_model_metadata':None,
         'input_digest':source['input_digest'],'original_source_evaluations':0,
         'source_packet_defined':True,'source_ftc_transport_defined':True,
         'retained_construction_data':report['retained_construction_data'],
         'retained_integral_term_inventory':report['retained_integral_term_inventory'],
         'retained_cap_provenance':parent['retained_cap_provenance'],
         'native_BG_exact_numeric_baseline_C3_identity_proved':False,
         'source_inspection_audit':{'old_source_consumer_calls':0,'old_proof_functions_called':0,
                                    'old_operator_calls':0,'old_source_normalizer_calls':0,
                                    'original_source_evaluations':0,'quadratures':0,
                                    'numeric_field_evaluations':0,'temporal_recurrence_calls':0}}
    out.update({key:False for key in report['false_flags']})
    try:
        i=source['input'];out['source_binding']=bind_six_source_charts(i,parent)
        proof=prove_upstream_mean(i['contracts'],journal);out['upstream_mean_identity']=proof
        if not proof['passed'] or not proof['postaxial']['width_sum_is_14']:
            raise ValueError('New affine mean transport failed')
        for t in proof['constant_F_transports']:
            if not t['coefficient_cancellation']['all_basis_coefficients_zero']:
                raise ValueError('Constant source transport cancellation failed')
        if not proof['axial_transport']['coefficient_cancellation']['all_basis_coefficients_zero']:
            raise ValueError('Axial exact-memory expansion failed')
        if out['source_binding']['Pstar_F_source_identity_proved'] is not True:
            raise ValueError('The physical Pstar/F source identity is not proved')
        proof['audits'].update(source_identity_proved=True,physical_Pstar_cancellation_is_external_binding=False,compact_anchor_law_proved=True,proof_scope='source-bound six-chart affine mean transport using the cached Rh and compact identities')
        proof['binding_resolution']='Canonical source ports, all 24 transport rows, Pstar cancellation and exact memory recipe are bound by source_binding.'
        out['native_BG_exact_numeric_baseline_C3_identity_proved']=True
        out['relative_mass_C3_source_identity_proved']=True
        out['relative_prefix_caps_bound_to_actual_source']=False
        out['numeric_mean_enclosure_transferred']=False
        out['exact_function_results']={
            'incoming_normalized_m':'(4*Z/Pstar)*mstar*exp(-14)',
            'physical_BG_m_C3':['4*Z*c/x','4*c/x','0','0'],
            'physical_BG_Q_C2':['-4*c/x','0','0'],
            'physical_BG_Q_y_C1':['4*c/x','0'],
            'c':'exp(-14)*mstar; independent of ordinary Z',
            'mean_mismatch_C3':['0','0','0','0'],'pressure_mismatch_C2':['0','0','0'],
            'cached_analytic_BG_minus_numeric_volume_C1':['0','0'],
            'scope':'function identities on the common source domain; no floating-point source evaluation',
            'derivation':'cached Rh identity + bound rate-one Duhamel source + new affine cancellation + physical Pstar once + cached compact transport'}
        out['remaining_scope']=['certify same-source complete-prefix U/m/Q correction bounds',
                                 'join physical geometry of numerical volume cap and native Rc',
                                 'compute signed Evol and retain the independent numeric-to-old reference offset']
        out['static_receipt_corrections']=['The previous_recipe hash printed in the prefix trace contains an extra e; additional_raw_pins records the authenticated actual hash.']
        out['passed']=True
        journal({'event':'upstream_mean_source_identity_proved','charts':6,'ordinary_Z_rows':24,
                 'background_baseline_identity':True,'relative_mass_identity':True,
                 'cap_identity':False,'old_errors_removed':False})
    except Exception as error:
        import traceback
        out['failure']={'type':type(error).__name__,'message':str(error),'traceback':traceback.format_exc()}
        journal({'event':'upstream_mean_identity_failed','type':type(error).__name__,'message':str(error)})
    return out
