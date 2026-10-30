"""Frozen one-entry same-source Kchi finite feedback and m15 subterms."""
from __future__ import annotations
import gzip,hashlib,importlib.util,json,sys,traceback
from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext
ROOT=Path(__file__).resolve().parents[1];BRIEF=ROOT.parent.parent/'RESUME_BRIEF_2026_10_16.json';COMPONENT='n1_core_kchi_feedback_v2';DIR=ROOT/'tmp'/COMPONENT
SOURCE,ENTRY,RESULT=(DIR/n for n in ('source.json','entered.json','result.json'))
KCHI='experiments/root_st073/n1_core_kchi_first_v1_receipt.json.gz';LOW='experiments/root_st073/n1_core_low_word_v3_receipt.json.gz';PREV='experiments/root_st073/n1_core_w0_m15_v1_receipt.json.gz';CORE='experiments/root_st073/lei_ren_part1_paper_compliant_core_physical_field.json';SWIRL='experiments/root_st073/n1_first_word_swirl_m15_v1_receipt.json.gz';CORRECTION='experiments/root_st073/n1_first_word_swirl_m15_v1_provenance_correction.json'
INNER='src/openai_ns_reconstruction/n1_kchi_inner_feedback_v2.py';IFIX='src/openai_ns_reconstruction/n1_kchi_inner_feedback_fixture_v2.py';OUTER='src/openai_ns_reconstruction/n1_kchi_outer_feedback_v2.py';OFIX='src/openai_ns_reconstruction/n1_kchi_outer_feedback_fixture_v2.py';M15='src/openai_ns_reconstruction/n1_kchi_feedback_m15_v2.py';MFIX='src/openai_ns_reconstruction/n1_kchi_feedback_m15_fixture_v2.py';RUNNER='scripts/run_n1_core_kchi_feedback_v2_once.py'
FAILURE='experiments/root_st073/n1_core_kchi_feedback_v1_failure.json.gz'
REQUIRED=(KCHI,LOW,PREV,CORE,SWIRL,CORRECTION,INNER,IFIX,OUTER,OFIX,M15,MFIX,RUNNER,FAILURE)
FIXED={FAILURE:'7e98bb70f811151bce271000ca2297b37cc2cf4da1ec5cd067ff97a8f6bd65e6',KCHI:'b0d654203e8c4acabdd0b61ca7b79c1b7050314b0c003a0c2f8757bde94973a0',LOW:'5a5fe7c57a2efe10e58089ced608e22937c54fcf54c5f224a42abe54e450c6d6',PREV:'157f5677742b71a469c67e0efd21de0d94d230b4d3c2968d82acd2972413bda7',CORE:'ad000166f6fce1ea94a7beeccf3324536fbffca073487b873d5460133f380096',SWIRL:'a626e71fa2b0ddf40b79ce992016da2af8ac9955db3f17337029f3c86c2513b4',CORRECTION:'f2133bf5824b2f9636b430a33a2ba1d19fbf391a7f90a81de0c8269d1fe638c4'}
def sha(b):return hashlib.sha256(b).hexdigest()
def require(v,msg):
    if not v:raise ValueError(msg)
def save(path,value):
    b=(json.dumps(value,sort_keys=True,separators=(',',':'))+'\n').encode()
    with path.open('xb') as f:f.write(b)
    return sha(b)
def preflight():
    require(not ENTRY.exists() and not RESULT.exists(),'already entered');b=BRIEF.read_bytes();d=json.loads(b);src=SOURCE.read_bytes();ck=d['source_checkpoint']
    require(d['component']==COMPONENT and d['phase']==COMPONENT+'_pending' and d['goal_active'] is True and d['full_goal_complete'] is False and d['one_entry_only'] is True and d['do_not_reenter'] is True,'brief guard')
    require(Path(ck['path']).resolve()==SOURCE.resolve() and ck['sha256']==sha(src),'source checkpoint');m=json.loads(src);require(m['component']==COMPONENT and m['one_entry_only'] is True,'manifest');h={}
    for row in m['files']:
        p=row['path'];require(p in REQUIRED and p not in h,'unexpected input');v=sha((ROOT/p).read_bytes());require(v==row['sha256'] and (p not in FIXED or v==FIXED[p]),'changed '+p);h[p]=v
    require(set(h)==set(REQUIRED),'missing input');return h,sha(src),sha(b)
def iv(c,r):
    lo=tuple(r['lower_exact_mpf_tuple']);hi=tuple(r['upper_exact_mpf_tuple']);require(len(lo)==len(hi)==4,'interval endpoints');x=c.make_mpf((lo,hi));require(x.a<=x.b,'reversed interval');return x
def conv(c,x):
    if isinstance(x,dict):
        if 'lower_exact_mpf_tuple' in x and 'upper_exact_mpf_tuple' in x:return iv(c,x)
        return {k:conv(c,v) for k,v in x.items()}
    if isinstance(x,list):return [conv(c,v) for v in x]
    return x
def enc(x):
    if hasattr(x,'_mpi_'):return {'lower_exact_mpf_tuple':list(x._mpi_[0]),'upper_exact_mpf_tuple':list(x._mpi_[1]),'directed':True}
    if isinstance(x,dict):return {str(k):enc(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [enc(v) for v in x]
    if isinstance(x,(str,int,float,bool)) or x is None:return x
    raise TypeError(type(x).__name__)
def module(path,name):
    spec=importlib.util.spec_from_file_location(name,ROOT/path);m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m
def consume(hashes):
    raw={}
    for p in REQUIRED:
        b=(ROOT/p).read_bytes();require(sha(b)==hashes[p],'changed after entry '+p);raw[p]=b
    kchi=json.loads(gzip.decompress(raw[KCHI]));low=json.loads(gzip.decompress(raw[LOW]));prev=json.loads(gzip.decompress(raw[PREV]));swirl=json.loads(gzip.decompress(raw[SWIRL]));core=json.loads(raw[CORE]);corr=json.loads(raw[CORRECTION])
    for receipt in (kchi,low,prev,swirl):require(receipt['result']['status']=='returned','input status')
    k=kchi['result']['output'];a=low['result']['output'];p=prev['result']['output'];s=swirl['result']['output'];family=a['source_family'];atoms=a['source_logs']
    require(k['source_family']==p['source_family']==s['source_family']==family and k['source_logs']==p['source_logs']==s['source_logs']==atoms,'source family')
    pins={v['path']:v['sha256'] for v in kchi['source']['files']}; oldpins={v['path']:v['sha256'] for v in prev['source']['files']}
    require(pins[LOW]==oldpins[LOW]==hashes[LOW] and pins[CORE]==oldpins[CORE]==hashes[CORE] and oldpins[SWIRL]==hashes[SWIRL] and oldpins[CORRECTION]==hashes[CORRECTION] and corr['receipt_sha256']==hashes[SWIRL],'authenticated source chain')
    require(k['feedback_cross_terms_added'] is False and k['appended_to_cumulative_m15'] is False and k['Kchi_tail_is_uncertainty_not_signed_field'] is True,'input scope')
    c=MPIntervalContext();c.prec=1800
    inner=module(INNER,'n1_kchi_inner_feedback_v2');ifix=module(IFIX,'n1_kchi_inner_feedback_fixture_v2');outer=module(OUTER,'n1_kchi_outer_feedback_v2');ofix=module(OFIX,'n1_kchi_outer_feedback_fixture_v2');m15=module(M15,'n1_kchi_feedback_m15_v2');mfix=module(MFIX,'n1_kchi_feedback_m15_fixture_v2')
    proofs={'inner':ifix.run(inner.build,c),'outer':ofix.run(outer.build,c),'m15':mfix.run(m15,c)}
    for name,v in proofs.items():require(v['passed'],name+' fixture failed '+str([key for key,value in v['checks'].items() if not value]))
    index={};matrices=conv(c,a['matrix_cells']);w0=conv(c,a['result']['outer_cells']);delta0=k['first_correction']['target_cells'];tail0=k['W0_separate_symmetric_kernel_tail'];require(len(matrices)==len(w0)==len(delta0)==8,'eight source cells')
    for cell in range(8):
        for name in ('C0','B1'):
            for order,entries in matrices[cell][name].items():
                for key,terms in entries.items():
                    for i,t in enumerate(terms):
                        sid=f'matrix:{cell}:{name}:{order}:{key}:{i}';index[sid]={'receipt':LOW,'path':['matrix_cells',cell,name,order,key,i]};t['source_id']=sid;t['origin']={'source_ids':[sid],'error_classes':['whole_core_matrix_interval'],'uncertainty':True}
        for order,cols in w0[cell]['W0_ordinary_Z'].items():
            if order not in ('0','1','2'):continue
            for col,terms in cols.items():
                for i,t in enumerate(terms):
                    sid=f'W0:{cell}:{order}:{col}:{i}';index[sid]={'receipt':LOW,'path':['result','outer_cells',cell,'W0_ordinary_Z',order,col,i]};t['source_id']=sid;t['origin']={'source_ids':[sid],'error_classes':['whole_core_W0_interval'],'uncertainty':True}
    delta=[];tailcells=[];delta_packets=[];tail_packets=[]
    for cell,rawcell in enumerate(delta0):
        x=iv(c,rawcell['x']);d={'x':x};e={'x':x};dp={'x':x,'W0_ordinary_Z':{str(n):{str(col):[] for col in range(1,7)} for n in range(3)}};ep={'x':x,'W0_ordinary_Z':{str(n):{str(col):[] for col in range(1,7)} for n in range(3)}}
        for field,col in (('F1','1'),('d_x_F1','5')):
            d[field]={};e[field]={}
            for order in (0,1,2):
                terms=conv(c,rawcell[field]['Z'+str(order)]);groups={}
                for pos,t in enumerate(terms):
                    key=(tuple(t['powers']),t['chi_order'],t['source_order']);g=groups.get(key)
                    if g is None:g={'coefficient':t['coefficient'],'ids':set(),'classes':set(),'parents':[]};groups[key]=g
                    else:g['coefficient']+=t['coefficient']
                    g['ids'].update(t['source_ids']);g['classes'].update(t['error_classes']);g['parents'].append(pos)
                vals=[]
                for number,(key,g) in enumerate(groups.items()):
                    sid=f'delta:{cell}:{field}:{order}:{number}';index[sid]={'receipt':KCHI,'path':['first_correction','target_cells',cell,field,'Z'+str(order)],'term_indices':g['parents']}
                    vals.append({'powers':key[0],'coefficient':g['coefficient'],'source_id':sid,'origin':{'source_ids':sorted(g['ids'])+[sid],'error_classes':sorted(g['classes']),'uncertainty':True,'aggregation':'outward sum only'}})
                ev=conv(c,tail0[field][str(order)])
                for n,t in enumerate(ev):
                    sid=f'W0tail:{field}:{order}:{n}';index[sid]={'receipt':KCHI,'path':['W0_separate_symmetric_kernel_tail',field,str(order),n]};t['source_id']=sid;t['origin']={'source_ids':t['source_ids']+[sid],'error_classes':t['error_classes'],'uncertainty':True,'error_only':True};t['error_only']=True
                d[field]['Z'+str(order)]=vals;e[field]['Z'+str(order)]=ev;dp['W0_ordinary_Z'][str(order)][col]=vals;ep['W0_ordinary_Z'][str(order)][col]=ev
        delta.append(d);tailcells.append(e);delta_packets.append(dp);tail_packets.append(ep)
    cj=conv(c,k['chi']['chi']);chi={n:[{'powers':[0,0,0,0],'coefficient':cj[n],'source_id':f'chi:Z{n}','origin':{'source_ids':[f'chi:Z{n}'],'error_classes':['whole_Z_chi_interval'],'uncertainty':True}}] for n in (0,1)}
    for n in (0,1):index[f'chi:Z{n}']={'receipt':KCHI,'path':['chi','chi',n]}
    iout=inner.build(c,delta,matrices,uncertainty_cells=tailcells)
    oout=outer.build(c,w0,matrices,chi,delta_w0_cells=delta_packets,include_d_x_F1=False)
    otail=outer.build(c,tail_packets,matrices,chi,branch='delta_W0',include_d_x_F1=False)
    kt=conv(c,k['kernel_tail_bounds']['F1'])
    def outer_remainder(branch,target,order,kind):
        pieces=[];qmap=branch['feedback_q5_products_by_source_cell']
        for source in range(target+1):
            for qorder in range(order+1):
                for term in qmap[str(source)]['Z'+str(qorder)]:
                    upper=c.mpf(1)/4*kt[order-qorder]*abs(term['coefficient']);sym=c.make_mpf(((-upper)._mpi_[0],upper._mpi_[1]));ids=list(term.get('source_ids',[]))+['chi:Z0','chi:Z1']
                    pieces.append({'powers':term['powers'],'coefficient':sym,'source_ids':ids,'origin':{'source_ids':ids,'error_classes':list(term.get('error_classes',[]))+[kind,'outer_Kchi_kernel_tail'],'uncertainty':True,'error_only':True},'error_only':True,'source_cell':source,'target_cell':target})
        return pieces
    whole=core['whole_core'];require(iv(c,whole['Z']).a==-1 and iv(c,whole['Z']).b==1 and iv(c,atoms['F0_anchor'])._mpi_==iv(c,whole['F0_exact_positive_log_enclosure'])._mpi_,'same leading whole core')
    U={};F={};phi=[]
    for order in (0,1):
        for field,out,key,powers in (('U',U,'Uz_derivatives',[0,0,0,0]),('Phi',F,'Phi_derivatives',[0,1,0,0])):
            sid=f'leading:{field}:{order}';index[sid]={'receipt':CORE,'path':['whole_core',key,'rho0_Z'+str(order)]};value=iv(c,whole[key]['rho0_Z'+str(order)]);out[order]=[{'powers':powers,'coefficient':value,'origin':{'source_ids':[sid],'error_classes':['whole_core_source_interval'],'uncertainty':True}}]
            if field=='Phi':phi.append(value)
    ell=iv(c,s['F0_connection']['ell_over_lambda']);index['F0_connection']={'receipt':SWIRL,'path':['F0_connection','ell_over_lambda'],'path_correction':CORRECTION}
    F[1].append({'powers':[0,1,0,1],'coefficient':phi[0]*ell,'origin':{'source_ids':['leading:Phi:0','F0_connection'],'error_classes':['whole_core_source_interval','F0_connection_interval'],'uncertainty':True}})
    obase=oout['branches']['K0_W0'];odelta=oout['branches']['delta_W0'];otailbranch=otail['branches']['delta_W0'];profile=[]
    for cell in range(8):
        t=w0[cell]['t'];bdata={name:{} for name in ('W0_chi','E01','E10','E11')}
        for order in (0,1):
            zo='Z'+str(order);inner_f=iout['signed_target_cells'][cell]['G0_'+zo];inner_e=iout['uncertainty_target_cells'][cell]['G0_'+zo];fbase=obase['target_cells'][cell]['F1'][zo];fdelta=odelta['target_cells'][cell]['F1'][zo];ftail=otailbranch['target_cells'][cell]['F1'][zo]
            def packed(field,error):return {'field':field,'symmetric_error':error}
            bdata['W0_chi'][zo]={'F1':packed(delta[cell]['F1'][zo],tailcells[cell]['F1'][zo])}
            bdata['E01'][zo]={fld:packed(inner_f[fld],inner_e[fld]) for fld in ('F1','Uz1')}
            bdata['E10'][zo]={'F1':packed(fbase,outer_remainder(obase,cell,order,'outer_tail_on_K0_W0'))}
            e11error=ftail+outer_remainder(odelta,cell,order,'outer_tail_on_delta_W0')+outer_remainder(otailbranch,cell,order,'outer_tail_on_W0_tail')
            bdata['E11'][zo]={'F1':packed(fdelta,e11error)}
        profile.append({'t':t,'rho':w0[cell]['rho'],'branches':bdata})
    contribution=m15.build(c,profile=profile,leading_U=U,leading_F=F,source_family=family,source_atoms=atoms,certificate={'geometry_epsilon_Z_independent':True})
    return {'source_family':family,'source_logs':atoms,'domain':'x in [0,2], real Z in [-1,1]','fixtures':proofs,'m15':enc(contribution),'input_index':index,
      'finite_branch_formulas':{'W0_chi':'Dg','E01':'G0 Q Dg','E10':'D Q G0g','E11':'D Q Dg'},'tail_error_terms':'G0 Q Rg; D Q Rg; R Q G0g; R Q Dg; R Q Rg',
      'first_correction_aggregated_only_by_formal_powers_and_derivative_pair':True,'aggregation_is_outward_only':True,'original_terms_retained_by_receipt_reference':True,
      'kernel_tails_kept_as_uncertainty':True,'feedback_word_tail_included':False,'appended_to_cumulative_m15':False,'complete_W1_computed':False,'complete_m15_computed':False,'terminal_closure_proved':False}
def main():
    h,s,b=preflight();entry=save(ENTRY,{'component':COMPONENT,'source_sha256':s,'brief_sha256':b,'input_hashes':h,'entered':True});out={'component':COMPONENT,'source_sha256':s,'entry_sha256':entry,'status':'failed'}
    try:out.update(status='returned',output=consume(h),error=None)
    except BaseException as exc:out['error']={'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()}
    digest=save(RESULT,out);print(json.dumps({'component':COMPONENT,'status':out['status'],'result_sha256':digest}));return 0 if out['status']=='returned' else 1
if __name__=='__main__':raise SystemExit(main())
