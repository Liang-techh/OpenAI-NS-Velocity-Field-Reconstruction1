"""One-entry source-bound first potential correction to collar row-5 forcing."""
from __future__ import annotations
import gzip,hashlib,importlib.util,json,sys,traceback
from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext
ROOT=Path(__file__).resolve().parents[1];BRIEF=ROOT.parent.parent/'RESUME_BRIEF_2026_10_16.json';COMPONENT='n1_collar_local_g5_j2_tail_v1';DIR=ROOT/'tmp'/COMPONENT
SOURCE,ENTRY,RESULT=(DIR/n for n in ('source.json','entered.json','result.json'))
FORCING='experiments/root_st073/n1_actual_collar_forcing_v1_receipt.json.gz';KCHI='experiments/root_st073/n1_core_kchi_first_v1_receipt.json.gz';PRESS='experiments/root_st073/n1_first_collar_g4_pressure_v1_receipt.json.gz';COMP='experiments/root_st073/n1_core_kchi_compensation_v2_receipt.json.gz';BIND='experiments/root_st073/n1_collar_axis_potential_binding_v1_receipt.json.gz';SELECT='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_selected_words.json'
FSRC='src/openai_ns_reconstruction/n1_actual_collar_forcing_compact.py';ORIGINAL='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_inner_volterra_functions.py';LEAF='src/openai_ns_reconstruction/n1_collar_local_g5_j2_tail_v1.py';FIXTURE='src/openai_ns_reconstruction/n1_collar_local_g5_j2_tail_fixture_v1.py';RUNNER='scripts/run_n1_collar_local_g5_j2_tail_v1_once.py'
PRIOR='experiments/root_st073/n1_collar_kchi_source_first_v3_receipt.json.gz'
FAILURE='experiments/root_st073/n1_collar_kchi_source_first_v1_failure.json.gz'
FAILURE2='experiments/root_st073/n1_collar_kchi_source_first_v2_failure.json.gz'
REQUIRED=(FORCING,KCHI,PRESS,COMP,BIND,SELECT,FSRC,ORIGINAL,LEAF,FIXTURE,RUNNER,PRIOR,FAILURE,FAILURE2)
FIXED={PRIOR:'bdfd705e43db41badd69b737a014f9ff288eaa58a18060fe736aebb196f18262',FAILURE2:'d1504b4a68d94f4fef60509bfc38279d07d5fd41b9b42bb9ae79ca467857341b',FAILURE:'276d7e98d57bb56a5e283a7ec46008f648a779499f3ccd0ecc5444d6e7ae7b9a',ORIGINAL:'ff82f4edb03922a551cba3e2df20217a55dc80ed42bd39ccb8b77b83f63b0208',FORCING:'96e8a3b64140fb2c111e8117951e96674dfa3c2bc7ff41ebe34448dccbbd8be2',KCHI:'b0d654203e8c4acabdd0b61ca7b79c1b7050314b0c003a0c2f8757bde94973a0',PRESS:'c229cf188d306413f6edd0641a729321ad4567c2587ebf4701ada47720853660',COMP:'c862f0560e8cb4b053abbbfcc982e4b4c1b602a806628e2d709399878e1c3aaf',BIND:'5859d5cb4bb82d372d8ac2f3ea1629e07c9787ffbaca149c99d1aba5edb8c6e8',SELECT:'0e04c87a011bfec74e0ceb0f822774a03495a08a839aaf4e378d21949736cf1e'}
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
    c=MPIntervalContext();c.prec=1800;leaf=module(LEAF,'n1_collar_local_g5_j2_tail_v1');fixture=module(FIXTURE,'n1_collar_local_g5_j2_tail_fixture_v1');proof=fixture.run_fixture(leaf,c);require(proof['passed'],'fixture failed '+str([k for k,v in proof['checks'].items() if not v]))
    forcing=json.loads(gzip.decompress(raw[FORCING]));kchi=json.loads(gzip.decompress(raw[KCHI]));pressure=json.loads(gzip.decompress(raw[PRESS]));comp=json.loads(gzip.decompress(raw[COMP]));binding=json.loads(gzip.decompress(raw[BIND]));prior=json.loads(gzip.decompress(raw[PRIOR]));selected=json.loads(raw[SELECT]);require(all(x['result']['status']=='returned' for x in (forcing,kchi,pressure,comp,binding,prior)),'accepted receipt status')
    fp={v['path']:v['sha256'] for v in forcing['source']['files']};pp={v['path']:v['sha256'] for v in pressure['source']['files']};cp={v['path']:v['sha256'] for v in comp['source']['files']};bp={v['path']:v['sha256'] for v in binding['source']['files']};prp={v['path']:v['sha256'] for v in prior['source']['files']}
    require(fp[FSRC]==hashes[FSRC] and fp[ORIGINAL]==hashes[ORIGINAL] and pp[FORCING]==hashes[FORCING] and all(cp[p]==hashes[p] for p in (FORCING,KCHI,PRESS,SELECT)) and bp[KCHI]==hashes[KCHI] and all(prp[p]==hashes[p] for p in (FORCING,KCHI,PRESS,COMP,BIND,SELECT)),'authenticated source chain')
    f=forcing['result']['output'];k=kchi['result']['output'];p=pressure['result']['output'];old=comp['result']['output'];bind=binding['result']['output'];prev=prior['result']['output'];family=f['source_family'];atoms=f['source_logs'];require(all(x['source_family']==family and x['source_logs']==atoms for x in (k,p,old,bind,prev)) and selected['source_family']==family,'same source family and formal atoms')
    contract=selected['actual_cutoff_complete_support_measure'];require(contract['passed'] and contract['source_family']==family and contract['measure_is_Z_independent'] and '0<hb<1/2' in contract['inequality_proof'],'certified physical geometry')
    sc=iv(c,p['sc']);scbits=tuple(contract['source_sc_exact_mpf_tuple']);require(sc._mpi_==(scbits,scbits) and sc.a>0,'same exact sc')
    fs=raw[FSRC].decode('utf8');source=raw[ORIGINAL].decode('utf8');require('F0_ratio, F1_ratio, F2_ratio = (Phi0 * relative_f0)' in fs and 'Nt_over_F0_anchor = -_axial_second' in fs and '"5": ((2 * Nt_over_F0_anchor, (0, 1, 0, -1)' in fs and '"F0_relative_derivative_applied_once": True' in fs and 'F(Z+h)/F0(Z) jet at each current Z; no extra D_F' in fs and 'def raw_system' in source and 'Nt=g.neg(axial_second(g,F,r,de,g.sub(c(-2),de)))' in source and 'g5=g.mul(c(2),eps,Nt)' in source,'raw g5 and F0 derivative source convention')
    require(bind['phase_shear_cutoff_identified_with_axis_potential'] is False and bind['source_checks']['one_query_jet_shared_in_integration'] is True,'axis potential binding')
    index={};g5={}
    for n in (0,1):
        terms=[]
        for i,t in enumerate(f['g'][str(n)]['5']):
            sid=f'collar_g5:Z{n}:{i}';origin=dict(t['origin']);parents=list(origin.get('source_ids',[]));origin['source_ids']=parents+[sid];index[sid]={'receipt':FORCING,'path':['result','output','g',str(n),'5',i],'parent_source_ids':parents};terms.append({'powers':t['powers'],'coefficient':iv(c,t['coefficient']),'origin':origin})
        require(len(terms)>0,'nonempty raw g5');g5['Z'+str(n)]=terms
    chi={}
    for n in (0,1):
        sid=f'chi:Z{n}';index[sid]={'receipt':KCHI,'path':['result','output','chi','chi',n]};chi['Z'+str(n)]=[{'powers':[0,0,0,0],'coefficient':iv(c,k['chi']['chi'][n]),'origin':{'source_ids':[sid],'error_classes':['whole_Z_axis_potential_interval'],'uncertainty':True}}]
    cert={'strict_hb_admitted':True,'hb_sc_Z_independent':True,'rho_radius_Z_independent':True,'source':'selected actual_cutoff_complete_support_measure','source_sha256':hashes[SELECT]}
    result=leaf.build(c,hb=c.mpf([0,'0.5']),sc=sc,g5=g5,chi=chi,source_family=family,source_atoms=['hb','F0_anchor','Pstar_squared','Lambda'],geometry_certificate=cert)
    return {'source_family':family,'source_logs':atoms,'source_index':index,'fixture':proof,'actual_local_j2_tail':enc(result),'first_potential_reference':{'receipt':PRIOR,'sha256':hashes[PRIOR]},'input_receipt_sha256':{path:hashes[path] for path in (FORCING,KCHI,PRESS,COMP,BIND,SELECT)},'raw_g5_frame_proven_from_pinned_source':True,'F0_derivative_not_reapplied':True,'source_enclosures_not_point_values':True,'hb_not_materialized':True,'baseline_K0_not_recomputed':True,'first_potential_recomputed':False,'local_j2_tail_only':True,'complete_collar_field':False,'feedback_word_tail_included':False,'m15_computed':False,'appended_to_ledger':False,'terminal_closure_proved':False}
def main():
    h,s,b=preflight();entry=save(ENTRY,{'component':COMPONENT,'source_sha256':s,'brief_sha256':b,'input_hashes':h,'entered':True});out={'component':COMPONENT,'source_sha256':s,'entry_sha256':entry,'status':'failed'}
    try:out.update(status='returned',output=consume(h),error=None)
    except BaseException as exc:out['error']={'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()}
    digest=save(RESULT,out);print(json.dumps({'component':COMPONENT,'status':out['status'],'result_sha256':digest}));return 0 if out['status']=='returned' else 1
if __name__=='__main__':raise SystemExit(main())
