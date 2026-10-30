"""One-entry Kchi compensation, fixed inlet subterm, and collar dx g4 audit."""
from __future__ import annotations
import gzip,hashlib,importlib.util,json,sys,traceback
from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext
ROOT=Path(__file__).resolve().parents[1];BRIEF=ROOT.parent.parent/'RESUME_BRIEF_2026_10_16.json';COMPONENT='n1_core_kchi_compensation_v2';DIR=ROOT/'tmp'/COMPONENT
SOURCE,ENTRY,RESULT=(DIR/n for n in ('source.json','entered.json','result.json'))
LOW='experiments/root_st073/n1_core_word1_z2_v2_receipt.json.gz';LOW3='experiments/root_st073/n1_core_low_word_v3_receipt.json.gz';KCHI='experiments/root_st073/n1_core_kchi_first_v1_receipt.json.gz';FEEDBACK='experiments/root_st073/n1_core_kchi_feedback_v2_receipt.json.gz';FORCING='experiments/root_st073/n1_actual_collar_forcing_v1_receipt.json.gz';PRESS='experiments/root_st073/n1_first_collar_g4_pressure_v1_receipt.json.gz';SELECT='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_selected_words.json';CORE='experiments/root_st073/lei_ren_part1_paper_compliant_core_physical_field.json';SWIRL='experiments/root_st073/n1_first_word_swirl_m15_v1_receipt.json.gz';CORR='experiments/root_st073/n1_first_word_swirl_m15_v1_provenance_correction.json'
WEIGHT='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_weighted_propagator.py';INITIAL='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_initial_integral.py';OP='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_inner_operator.py'
EXIT='src/openai_ns_reconstruction/n1_kchi_core_exit_v2.py';EF='src/openai_ns_reconstruction/n1_kchi_core_exit_fixture_v2.py';COMP='src/openai_ns_reconstruction/n1_kchi_physical_compensation_v2.py';CF='src/openai_ns_reconstruction/n1_kchi_physical_compensation_fixture_v2.py';MEASURE='src/openai_ns_reconstruction/n1_collar_green_measure_v2.py';MF='src/openai_ns_reconstruction/n1_collar_green_measure_fixture_v2.py';RUNNER='scripts/run_n1_core_kchi_compensation_v2_once.py'
FAILURE='experiments/root_st073/n1_core_kchi_compensation_v1_failure.json.gz'
REQUIRED=(LOW,LOW3,KCHI,FEEDBACK,FORCING,PRESS,SELECT,CORE,SWIRL,CORR,WEIGHT,INITIAL,OP,EXIT,EF,COMP,CF,MEASURE,MF,RUNNER,FAILURE)
FIXED={FAILURE:'33babbe901df97c2f73c00de748707cfa7f33dfd88d7bcbd8de43da2cf9b842b',LOW:'c53ff32ad3ce5288279953e84bc7417b71477c82b00d94466d0d2880adbc7ac9',LOW3:'5a5fe7c57a2efe10e58089ced608e22937c54fcf54c5f224a42abe54e450c6d6',KCHI:'b0d654203e8c4acabdd0b61ca7b79c1b7050314b0c003a0c2f8757bde94973a0',FEEDBACK:'0ff1992f08754f2e15c8ec8e65f3ed396b38dc3101821c35eb90635664a753aa',FORCING:'96e8a3b64140fb2c111e8117951e96674dfa3c2bc7ff41ebe34448dccbbd8be2',PRESS:'c229cf188d306413f6edd0641a729321ad4567c2587ebf4701ada47720853660',SELECT:'0e04c87a011bfec74e0ceb0f822774a03495a08a839aaf4e378d21949736cf1e',CORE:'ad000166f6fce1ea94a7beeccf3324536fbffca073487b873d5460133f380096',SWIRL:'a626e71fa2b0ddf40b79ce992016da2af8ac9955db3f17337029f3c86c2513b4',CORR:'f2133bf5824b2f9636b430a33a2ba1d19fbf391a7f90a81de0c8269d1fe638c4',WEIGHT:'5219a968018601d8c445815058f7c1c832398dc417fcb51c689688d03011ba62',OP:'268217c6b33bc717472742fa64d76a853c6760a956abfad92ac45bcf385df014'}
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
    require(Path(ck['path']).resolve()==SOURCE.resolve() and ck['sha256']==sha(src),'checkpoint');m=json.loads(src);require(m['component']==COMPONENT and m['one_entry_only'] is True,'manifest');h={}
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
    for path in REQUIRED:
        b=(ROOT/path).read_bytes();require(sha(b)==hashes[path],'changed after entry '+path);raw[path]=b
    low=json.loads(gzip.decompress(raw[LOW]));low3=json.loads(gzip.decompress(raw[LOW3]));kc=json.loads(gzip.decompress(raw[KCHI]));fb=json.loads(gzip.decompress(raw[FEEDBACK]));forcing=json.loads(gzip.decompress(raw[FORCING]));pressure=json.loads(gzip.decompress(raw[PRESS]));swirl=json.loads(gzip.decompress(raw[SWIRL]));core=json.loads(raw[CORE]);selected=json.loads(raw[SELECT]);correction=json.loads(raw[CORR])
    for r in (low,low3,kc,fb,forcing,pressure,swirl):require(r['result']['status']=='returned','receipt not returned')
    a=low['result']['output'];b=low3['result']['output'];k=kc['result']['output'];f=forcing['result']['output'];p=pressure['result']['output'];s=swirl['result']['output'];family=b['source_family'];atoms=b['source_logs']
    require(all(v['source_family']==family for v in (a,k,fb['result']['output'],f,p,s)) and selected['source_family']==family and all(v['source_logs']==atoms for v in (k,fb['result']['output'],f,p,s)),'same source family and atoms')
    require(b['input_receipts'][LOW]==hashes[LOW] and correction['receipt_sha256']==hashes[SWIRL] and p['core_m15_inlet_added'] is False,'prior source linkage')
    ks={v['path']:v['sha256'] for v in kc['source']['files']};ps={v['path']:v['sha256'] for v in pressure['source']['files']};require(ks[LOW]==hashes[LOW] and ks[LOW3]==hashes[LOW3] and ps[FORCING]==hashes[FORCING],'same source receipt pins')
    weighted=raw[WEIGHT].decode('utf8');initial=raw[INITIAL].decode('utf8');op=raw[OP].decode('utf8')
    require('WEIGHTS=((0,1,0,1)' in weighted and 'S^-1 B0 S+S^-1 B1 S_Z' in weighted and 'potential=algebra.lift(-2*chi)' in weighted and "residual['5,1']=residual['5,1']-potential" in weighted and 'if collar:value=algebra.width(value)' in initial and "B0['5,5']=eps(B*x)" in op,'source normalization and formal Jacobian evidence')
    measure_contract=selected['actual_cutoff_complete_support_measure'];require(measure_contract['passed'] and measure_contract['source_family']==family and measure_contract['measure_is_Z_independent'] and '0<hb<1/2' in measure_contract['inequality_proof'],'geometry certificate')
    c=MPIntervalContext();c.prec=1800;sc=iv(c,p['sc']);bits=tuple(measure_contract['source_sc_exact_mpf_tuple']);require(sc._mpi_==(bits,bits),'same exact sc')
    ex=module(EXIT,'n1_kchi_core_exit_v2');ef=module(EF,'n1_kchi_core_exit_fixture_v2');comp=module(COMP,'n1_kchi_physical_compensation_v2');cf=module(CF,'n1_kchi_physical_compensation_fixture_v2');measure=module(MEASURE,'n1_collar_green_measure_v2');mf=module(MF,'n1_collar_green_measure_fixture_v2')
    proofs={'exit':ef.run_fixture(ex.build,c),'compensation':cf.run_fixture(comp.build_compensation,c),'measure':mf.run(measure,c)}
    for name,v in proofs.items():require(v['passed'],name+' fixture failed '+str([q for q,x in v['checks'].items() if not x]))
    index={};w0=[];delta=[];g5cells=[];chi={}
    for n in (0,1):
        sid=f'chi:Z{n}';index[sid]={'receipt':KCHI,'path':['chi','chi',n]};chi[n]=[{'powers':[0,0,0,0],'coefficient':iv(c,k['chi']['chi'][n]),'origin':{'source_ids':[sid],'error_classes':['whole_Z_chi_interval'],'uncertainty':True}}]
    for cell in range(8):
        old=b['result']['outer_cells'][cell];dc=k['first_correction']['target_cells'][cell];t=iv(c,old['t']);x=iv(c,dc['x']);require((2*t)._mpi_==x._mpi_,'same x/t source cell');w={'t':t,'F1':{}};d={'t':t,'F1':{}};g={'t':t,'F1':{}}
        for n in (0,1):
            o=[]
            for i,term in enumerate(old['W0_ordinary_Z'][str(n)]['1']):
                sid=f'W0:{cell}:{n}:{i}';index[sid]={'receipt':LOW3,'path':['result','outer_cells',cell,'W0_ordinary_Z',str(n),'1',i]};o.append({'powers':[0,term['F0_power'],0,term['Lambda_power']],'coefficient':iv(c,term['coefficient']),'origin':{'source_ids':[sid],'error_classes':['whole_core_W0_interval'],'uncertainty':True}})
            w['F1']['Z'+str(n)]=o;groups={}
            for i,term in enumerate(dc['F1']['Z'+str(n)]):
                key=(tuple(term['powers']),term['chi_order'],term['source_order']);item=groups.get(key)
                if item is None:item={'coefficient':iv(c,term['coefficient']),'ids':set(),'classes':set(),'parents':[]};groups[key]=item
                else:item['coefficient']+=iv(c,term['coefficient'])
                item['ids'].update(term['source_ids']);item['classes'].update(term['error_classes']);item['parents'].append(i)
            signed=[]
            for i,(key,item) in enumerate(groups.items()):
                sid=f'delta:{cell}:{n}:{i}';index[sid]={'receipt':KCHI,'path':['first_correction','target_cells',cell,'F1','Z'+str(n)],'term_indices':item['parents']};signed.append({'powers':key[0],'coefficient':item['coefficient'],'origin':{'source_ids':sorted(item['ids'])+[sid],'error_classes':sorted(item['classes']),'uncertainty':True}})
            d['F1']['Z'+str(n)]=signed;source=[]
            for i,term in enumerate(a['global_low_source'][str(n)]['5']):
                sid=f'g5:{n}:{i}';index[sid]={'receipt':LOW,'path':['global_low_source',str(n),'5',i]};source.append({'powers':[0,term['F0_power'],0,term['Lambda_power']],'coefficient':iv(c,term['coefficient']),'origin':{'source_ids':[sid],'error_classes':['whole_core_source_interval'],'uncertainty':True}})
            g['F1']['Z'+str(n)]=source
        w0.append(w);delta.append(d);g5cells.append(g)
    exit_result=ex.build(c,g5cells,chi);comp_result=comp.build_compensation(c,w0,delta,chi,include_d_x_F1=True)
    whole=core['whole_core'];require(iv(c,whole['Z']).a==-1 and iv(c,whole['Z']).b==1 and iv(c,atoms['F0_anchor'])._mpi_==iv(c,whole['F0_exact_positive_log_enclosure'])._mpi_,'same leading field')
    phi0=iv(c,whole['Phi_derivatives']['rho0_Z0']);phi1=iv(c,whole['Phi_derivatives']['rho0_Z1']);ell=iv(c,s['F0_connection']['ell_over_lambda'])
    leading={0:[{'powers':(0,1,0,0),'coefficient':phi0,'ids':['leading:Phi:0']}],1:[{'powers':(0,1,0,0),'coefficient':phi1,'ids':['leading:Phi:1']},{'powers':(0,1,0,1),'coefficient':phi0*ell,'ids':['leading:Phi:0','F0_connection']}]}
    index['leading:Phi:0']={'receipt':CORE,'path':['whole_core','Phi_derivatives','rho0_Z0']};index['leading:Phi:1']={'receipt':CORE,'path':['whole_core','Phi_derivatives','rho0_Z1']};index['F0_connection']={'receipt':SWIRL,'path':['F0_connection','ell_over_lambda'],'path_correction':CORR}
    m15={}
    for branch,packet in comp_result['branches'].items():
        m15[branch]={}
        for order in (0,1):
            cells={}
            for cell,target in enumerate(packet['target_cells']):
                lo=c.mpf(cell)/8;hi=c.mpf(cell+1)/8;mass=-16*(hi**4-lo**4);terms=[];rules=((0,0),) if order==0 else ((1,0),(0,1))
                for lead_order,word_order in rules:
                    for lv in leading[lead_order]:
                        for v in target['fields']['Z'+str(word_order)]['F1']:
                            powers=[lv['powers'][j]+v['powers'][j]+(0,0,0,-2)[j] for j in range(4)];ids=lv['ids']+list(v['source_ids'])
                            terms.append({'powers':powers,'coefficient':mass*lv['coefficient']*v['coefficient'],'origin':{'source_ids':ids,'error_classes':list(v.get('error_classes',[]))+['whole_core_source_interval'],'uncertainty':True,'target_cell':cell,'lead_order':lead_order,'word_order':word_order,'symbolic_cancellation_claimed':False}})
                cells[cell]=terms
            m15[branch][order]=cells
    g4={}
    for n in (0,1):
        g4[n]=[]
        for i,t in enumerate(f['g'][str(n)]['4']):
            sid=f'collar_g4:{n}:{i}';index[sid]={'receipt':FORCING,'path':['g',str(n),'4',i],'parent_source_ids':t['origin']['source_ids']};origin=dict(t['origin']);origin['source_ids']=list(origin['source_ids'])+[sid];g4[n].append({'powers':t['powers'],'coefficient':iv(c,t['coefficient']),'origin':origin})
    geom=measure.build(c,sc=sc,hb=c.mpf([0,'0.5']),cell_count=8,g4=g4,source_family=family,source_atoms=atoms,geometry_certificate={'hb_sc_Z_independent':True,'rho_radius_Z_independent':True},formal_hb_owner='single_dx_jacobian',physical_hb_owner='0<hb<1/2 geometry enclosure',geometry_owner='same_original_chart')
    for n in (0,1):
        require(len(geom['source_integral']['Z'+str(n)])==8*len(g4[n]),'collar source missing or duplicated')
        for term in geom['source_integral']['Z'+str(n)]:require(term['powers'][0]>=0,'uncanceled inverse width')
    compact={name:{'exit_aggregate':packet['exit_aggregate'],'exit_by_source_cell':packet['exit_by_source_cell']} for name,packet in comp_result['branches'].items()}
    return {'source_family':family,'source_logs':atoms,'fixtures':proofs,'input_index':index,'signed_delta_W0_fixed_exit':enc(exit_result['signed_delta_g5_exit']),'saved_W0_symmetric_kernel_tail':k['W0_separate_symmetric_kernel_tail'],
      'potential_compensation_exit':enc(compact),'potential_compensation_m15':enc(m15),'collar_local_g4_dx':enc({'source_integral':geom['source_integral'],'source_cell_weights':geom['source_cell_weights'],'ownership':geom['ownership'],'q_definition':geom['q_definition'],'dx_dq_formula':geom['dx_dq_formula']}),
      'normalization_identity_scope':'S-conjugated Achi extraction implies physical +2chi E51; remaining gauge/source identity remains conditional',
      'W0_kernel_tail_not_reintegrated':True,'potential_feedback_tail_included':False,'collar_g4_is_alternative_P1_local_source_not_m15':True,'collar_Kchi_homogeneous_propagator_included':False,'no_reset_collar_join_computed':False,'appended_to_cumulative_m15':False,'complete_W1_computed':False,'complete_m15_computed':False,'terminal_closure_proved':False}
def main():
    h,s,b=preflight();entry=save(ENTRY,{'component':COMPONENT,'source_sha256':s,'brief_sha256':b,'input_hashes':h,'entered':True});out={'component':COMPONENT,'source_sha256':s,'entry_sha256':entry,'status':'failed'}
    try:out.update(status='returned',output=consume(h),error=None)
    except BaseException as exc:out['error']={'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()}
    digest=save(RESULT,out);print(json.dumps({'component':COMPONENT,'status':out['status'],'result_sha256':digest}));return 0 if out['status']=='returned' else 1
if __name__=='__main__':raise SystemExit(main())
