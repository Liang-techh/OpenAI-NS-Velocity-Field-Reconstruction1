"""One-entry formal collar homogeneous j>=2 tail for source-bound inlet groups."""
from __future__ import annotations
import gzip,hashlib,importlib.util,json,sys,traceback
from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext
ROOT=Path(__file__).resolve().parents[1];BRIEF=ROOT.parent.parent/'RESUME_BRIEF_2026_10_16.json';COMPONENT='n1_collar_inlet_j2_tail_v3';DIR=ROOT/'tmp'/COMPONENT
SOURCE,ENTRY,RESULT=(DIR/n for n in ('source.json','entered.json','result.json'))
BASE='experiments/root_st073/n1_core_six_inlet_v2_receipt.json.gz';KCHI='experiments/root_st073/n1_core_kchi_first_v1_receipt.json.gz';COMP='experiments/root_st073/n1_core_kchi_compensation_v2_receipt.json.gz';DIRECT='experiments/root_st073/n1_kchi_inlet_tail_bind_v1_receipt.json.gz';BIND='experiments/root_st073/n1_collar_axis_potential_binding_v1_receipt.json.gz';FIRST='experiments/root_st073/n1_collar_core_inlet_potential_v5_receipt.json.gz';LOCAL='experiments/root_st073/n1_collar_kchi_source_first_v3_receipt.json.gz';PRESS='experiments/root_st073/n1_first_collar_g4_pressure_v1_receipt.json.gz';SELECT='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_selected_words.json'
LEAF='src/openai_ns_reconstruction/n1_collar_inlet_tail_v3.py';FIXTURE='src/openai_ns_reconstruction/n1_collar_inlet_tail_fixture_v3.py';RUNNER='scripts/run_n1_collar_inlet_j2_tail_v3_once.py'
FAILURE='experiments/root_st073/n1_collar_inlet_j2_tail_v1_failure.json.gz'
FAILURE2='experiments/root_st073/n1_collar_inlet_j2_tail_v2_failure.json.gz'
REQUIRED=(BASE,KCHI,COMP,DIRECT,BIND,FIRST,LOCAL,PRESS,SELECT,LEAF,FIXTURE,RUNNER,FAILURE,FAILURE2)
FIXED={FAILURE2:'0396a90dc86e86c8d47b040e7e8400ef853f56e9ed284e7d9d35fd849e8e1220',FAILURE:'4da65b5a4d2c8ae7f6304baae8bd1d954ef8cd0f8afcbc89a542d7bb8723f8d5',BASE:'59c0ff88652d83d7b92ba9c3a564ee90eac5ae794220b8b858c112936b76342d',KCHI:'b0d654203e8c4acabdd0b61ca7b79c1b7050314b0c003a0c2f8757bde94973a0',COMP:'c862f0560e8cb4b053abbbfcc982e4b4c1b602a806628e2d709399878e1c3aaf',DIRECT:'5ee828e3ab4ea5ecc9d1e3eeb283077796b51a48eaa3a2b0bc056d5c5775bfd9',BIND:'5859d5cb4bb82d372d8ac2f3ea1629e07c9787ffbaca149c99d1aba5edb8c6e8',FIRST:'1c644fb0f41a8ae407c09f56b00e7aee1720f91d08163a0c7483b4b1b6c61e7a',LOCAL:'bdfd705e43db41badd69b737a014f9ff288eaa58a18060fe736aebb196f18262',PRESS:'c229cf188d306413f6edd0641a729321ad4567c2587ebf4701ada47720853660',SELECT:'0e04c87a011bfec74e0ceb0f822774a03495a08a839aaf4e378d21949736cf1e'}
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
        v=(ROOT/path).read_bytes();require(sha(v)==hashes[path],'changed after entry '+path);raw[path]=v
    c=MPIntervalContext();c.prec=1800;leaf=module(LEAF,'n1_collar_inlet_tail_v3');fixture=module(FIXTURE,'n1_collar_inlet_tail_fixture_v3');proof=fixture.run_fixture(leaf,c);require(proof['passed'],'fixture failed '+str([k for k,v in proof['checks'].items() if not v]))
    packets={path:json.loads(gzip.decompress(raw[path])) for path in (BASE,KCHI,COMP,DIRECT,BIND,FIRST,LOCAL,PRESS)};selected=json.loads(raw[SELECT]);require(all(x['result']['status']=='returned' for x in packets.values()),'accepted receipt status')
    b,k,comp,d,bind,first,local,p=(packets[path]['result']['output'] for path in (BASE,KCHI,COMP,DIRECT,BIND,FIRST,LOCAL,PRESS));family=b['source_family'];atoms=b['source_logs'];require(all(v['source_family']==family and v['source_logs']==atoms for v in (k,comp,d,bind,first,local,p)) and selected['source_family']==family,'same source family and formal atoms')
    fp={v['path']:v['sha256'] for v in packets[FIRST]['source']['files']};require(all(fp[path]==hashes[path] for path in (BASE,KCHI,COMP,DIRECT,BIND,LOCAL,PRESS,SELECT)),'authenticated first-word source chain')
    contract=selected['actual_cutoff_complete_support_measure'];require(contract['passed'] and contract['source_family']==family and contract['measure_is_Z_independent'] and '0<hb<1/2' in contract['inequality_proof'],'certified geometry');sc=iv(c,p['sc']);bits=tuple(contract['source_sc_exact_mpf_tuple']);require(sc._mpi_==(bits,bits) and sc.a>0,'same exact sc')
    require(bind['source_checks']['one_query_jet_shared_in_integration'] is True and d['direct_tail_is_symmetric_error_only'] is True and first['finite_input_branches_kept_distinct'] is True,'source scope')
    sources={'G0g':(BASE,b['baseline_G0g_fixed_exit']['G0g']['aggregate'],['result','output','baseline_G0g_fixed_exit','G0g','aggregate']),
       'core_Dg':(COMP,comp['signed_delta_W0_fixed_exit']['exit_aggregate'],['result','output','signed_delta_W0_fixed_exit','exit_aggregate'])}
    index={};inlets={}
    for name,(receipt,data,path) in sources.items():
        packet={'Fc':{},'Fxc':{}}
        for field,component in (('F1','Fc'),('d_x_F1','Fxc')):
            for order in (0,1):
                key='Z'+str(order);terms=[]
                for i,item in enumerate(data[key][field]):
                    t=conv(c,item);sid=f'{name}:{field}:{key}:{i}';origin=t.get('origin');om=origin if isinstance(origin,dict) else {};parents=list(t.get('source_ids',[]))+list(om.get('source_ids',[]));classes=list(t.get('error_classes',[]))+list(om.get('error_classes',[]));powers=t.get('powers',[0,t.get('F0_power',0),0,t.get('Lambda_power',0)]);index[sid]={'receipt':receipt,'path':path+[key,field,i],'parent_source_ids':parents,'parent_origin':item.get('origin')};terms.append({'powers':powers,'coefficient':t['coefficient'],'origin':{'source_ids':parents+[sid],'error_classes':classes,'uncertainty':True}})
                require(len(terms)>0,'nonempty inlet '+name+field+key);packet[component][key]=terms
        inlets[name]=packet
    tail={'Fc':{},'Fxc':{}}
    for field,component in (('F1','Fc'),('d_x_F1','Fxc')):
        for order in (0,1):
            key='Z'+str(order);terms=[]
            for i,item in enumerate(d['direct_W0_kernel_tail'][field][str(order)]):
                t=conv(c,item);sid=f'core_tail:{field}:{key}:{i}';parents=list(t['source_ids']);classes=list(t['error_classes']);index[sid]={'receipt':DIRECT,'path':['result','output','direct_W0_kernel_tail',field,str(order),i],'parent_source_ids':parents};terms.append({'powers':t['powers'],'coefficient':t['coefficient'],'origin':{'source_ids':parents+[sid],'error_classes':classes,'uncertainty':True,'error_only':True,'symmetric_error':True},'error_only':True})
            require(len(terms)>0,'nonempty direct core tail');tail[component][key]=terms
    inlets['core_symmetric_tail']=tail;chi={}
    for order in (0,1):
        sid=f'chi:Z{order}';index[sid]={'receipt':KCHI,'path':['result','output','chi','chi',order]};chi['Z'+str(order)]=[{'powers':[0,0,0,0],'coefficient':iv(c,k['chi']['chi'][order]),'origin':{'source_ids':[sid],'error_classes':['whole_Z_axis_potential_interval'],'uncertainty':True}}]
    cert={'strict_hb_admitted':True,'hb_sc_Z_independent':True,'rho_radius_Z_independent':True,'source':'selected actual_cutoff_complete_support_measure','source_sha256':hashes[SELECT]};outputs=leaf.build(c,inlet_groups=inlets,chi=chi,hb=c.mpf([0,'0.5']),sc=sc,source_family=family,source_atoms=['hb','F0_anchor','Pstar_squared','Lambda'],geometry_certificate=cert)
    return {'source_family':family,'source_logs':atoms,'source_index':index,'fixture':proof,'inlet_j2_tail':enc(outputs),'first_word_reference':{'receipt':FIRST,'path':['result','output','inlet_branches']},'input_receipt_sha256':{path:hashes[path] for path in (BASE,KCHI,COMP,DIRECT,BIND,FIRST,LOCAL,PRESS,SELECT)},'partial_inlet_groups_only':True,'collar_kernel_tail_distinct_from_feedback_tail':True,'formal_atoms_not_materialized':True,'first_word_not_recomputed':True,'complete_collar_field':False,'feedback_word_tail_included':False,'m15_computed':False,'appended_to_ledger':False,'terminal_closure_proved':False}
def main():
    h,s,b=preflight();entry=save(ENTRY,{'component':COMPONENT,'source_sha256':s,'brief_sha256':b,'input_hashes':h,'entered':True});out={'component':COMPONENT,'source_sha256':s,'entry_sha256':entry,'status':'failed'}
    try:out.update(status='returned',output=consume(h),error=None)
    except BaseException as exc:out['error']={'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()}
    digest=save(RESULT,out);print(json.dumps({'component':COMPONENT,'status':out['status'],'result_sha256':digest}));return 0 if out['status']=='returned' else 1
if __name__=='__main__':raise SystemExit(main())
