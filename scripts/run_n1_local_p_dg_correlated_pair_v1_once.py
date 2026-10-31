"""One-entry conditional local P/Dg correlated m15 source pair."""
from __future__ import annotations
import gzip,hashlib,importlib.util,json,sys,traceback
from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext
ROOT=Path(__file__).resolve().parents[1];BRIEF=ROOT.parent.parent/'RESUME_BRIEF_2026_10_16.json';COMPONENT='n1_local_p_dg_correlated_pair_v1';DIR=ROOT/'tmp'/COMPONENT
SOURCE,ENTRY,RESULT=(DIR/n for n in ('source.json','entered.json','result.json'))
PHYS='experiments/root_st073/n1_local_g5_potential_physical_m15_v1_receipt.json.gz';FIRST='experiments/root_st073/n1_collar_kchi_source_first_v3_receipt.json.gz';W0='experiments/root_st073/n1_w0_collar_sectors_v1_receipt.json.gz';FORCING='experiments/root_st073/n1_actual_collar_forcing_v1_receipt.json.gz';KCHI='experiments/root_st073/n1_core_kchi_first_v1_receipt.json.gz';BIND='experiments/root_st073/n1_collar_axis_potential_binding_v1_receipt.json.gz';COMP='experiments/root_st073/n1_core_kchi_compensation_v2_receipt.json.gz';SELECT='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_selected_words.json';PRESS='experiments/root_st073/n1_first_collar_g4_pressure_v1_receipt.json.gz'
W0LEAF='src/openai_ns_reconstruction/n1_w0_collar_sectors_v1.py';FIRSTLEAF='src/openai_ns_reconstruction/n1_collar_sector_green_v3.py';COMPLEAF='src/openai_ns_reconstruction/n1_kchi_physical_compensation_v2.py';LEAF='src/openai_ns_reconstruction/n1_local_p_dg_correlated_pair_v1.py';FIXTURE='src/openai_ns_reconstruction/n1_local_p_dg_correlated_pair_fixture_v1.py';RUNNER='scripts/run_n1_local_p_dg_correlated_pair_v1_once.py'
REQUIRED=(PHYS,FIRST,W0,FORCING,KCHI,BIND,COMP,SELECT,PRESS,W0LEAF,FIRSTLEAF,COMPLEAF,LEAF,FIXTURE,RUNNER)
FIXED={PHYS:'55f25fa95c455f03f3e2af99645da3f29590040b1732f0b266275efb35dc10ba',FIRST:'bdfd705e43db41badd69b737a014f9ff288eaa58a18060fe736aebb196f18262',W0:'4229c29a2a21207f8b8da3ba53d70bbf6f50e09eb67110bd3b7cb30fa6577c06',FORCING:'96e8a3b64140fb2c111e8117951e96674dfa3c2bc7ff41ebe34448dccbbd8be2',KCHI:'b0d654203e8c4acabdd0b61ca7b79c1b7050314b0c003a0c2f8757bde94973a0',BIND:'5859d5cb4bb82d372d8ac2f3ea1629e07c9787ffbaca149c99d1aba5edb8c6e8',COMP:'c862f0560e8cb4b053abbbfcc982e4b4c1b602a806628e2d709399878e1c3aaf',SELECT:'0e04c87a011bfec74e0ceb0f822774a03495a08a839aaf4e378d21949736cf1e',PRESS:'c229cf188d306413f6edd0641a729321ad4567c2587ebf4701ada47720853660'}
def sha(b):return hashlib.sha256(b).hexdigest()
def require(v,msg):
    if not v:raise ValueError(msg)
def save(path,value):
    b=(json.dumps(value,sort_keys=True,separators=(',',':'))+'\n').encode()
    with path.open('xb') as f:f.write(b)
    return sha(b)
def preflight():
    require(not ENTRY.exists() and not RESULT.exists(),'already entered');b=BRIEF.read_bytes();d=json.loads(b);src=SOURCE.read_bytes();ck=d['source_checkpoint'];require(d['component']==COMPONENT and d['phase']==COMPONENT+'_pending' and d['goal_active'] is True and d['full_goal_complete'] is False and d['one_entry_only'] is True and d['do_not_reenter'] is True,'brief guard');require(Path(ck['path']).resolve()==SOURCE.resolve() and ck['sha256']==sha(src),'checkpoint');m=json.loads(src);require(m['component']==COMPONENT and m['one_entry_only'] is True,'manifest');h={}
    for row in m['files']:
        p=row['path'];require(p in REQUIRED and p not in h,'unexpected input');v=sha((ROOT/p).read_bytes());require(v==row['sha256'] and (p not in FIXED or v==FIXED[p]),'changed '+p);h[p]=v
    require(set(h)==set(REQUIRED),'missing input');return h,sha(src),sha(b)
def iv(c,r):
    lo=tuple(r['lower_exact_mpf_tuple']);hi=tuple(r['upper_exact_mpf_tuple']);require(len(lo)==len(hi)==4,'interval endpoints');x=c.make_mpf((lo,hi));require(x.a<=x.b,'reversed interval');return x
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
    c=MPIntervalContext();c.prec=1800;leaf=module(LEAF,'n1_local_p_dg_correlated_pair_v1');fixture=module(FIXTURE,'n1_local_p_dg_correlated_pair_fixture_v1');proof=fixture.run_fixture(leaf.build,c);require(proof['passed'],'fixture failed '+json.dumps(proof,sort_keys=True,separators=(',',':')))
    receipts={p:json.loads(gzip.decompress(raw[p])) for p in (PHYS,FIRST,W0,FORCING,KCHI,BIND,COMP,PRESS)};require(all(r['result']['status']=='returned' for r in receipts.values()),'receipt status');phys,first,w0,forcing,kchi,bind,comp,press=(receipts[k]['result']['output'] for k in (PHYS,FIRST,W0,FORCING,KCHI,BIND,COMP,PRESS));selected=json.loads(raw[SELECT]);family=phys['source_family'];atoms=phys['source_logs'];require(all(x['source_family']==family and x['source_logs']==atoms for x in (first,w0,forcing,kchi,bind,comp,press)) and selected['source_family']==family,'common source family and logs')
    pins=lambda r:{x['path']:x['sha256'] for x in r['source']['files']};ph=pins(receipts[PHYS]);fi=pins(receipts[FIRST]);wo=pins(receipts[W0]);cp=pins(receipts[COMP]);bi=pins(receipts[BIND]);require(ph[FIRST]==hashes[FIRST] and ph[SELECT]==hashes[SELECT] and fi[FORCING]==wo[FORCING]==cp[FORCING]==hashes[FORCING] and fi[KCHI]==cp[KCHI]==bi[KCHI]==hashes[KCHI] and fi[BIND]==hashes[BIND] and fi[COMP]==hashes[COMP] and fi[SELECT]==cp[SELECT]==hashes[SELECT] and wo[W0LEAF]==hashes[W0LEAF] and fi[FIRSTLEAF]==hashes[FIRSTLEAF] and cp[COMPLEAF]==hashes[COMPLEAF],'authenticated common owner/code chain')
    require(w0['source_forcing_receipt']==FORCING and w0['result']['g_uniform_on_entire_first_collar_required'] and w0['result']['local_added_once_in_definition'] and w0['allword_tail_added'] is False,'local W0 domain');aj=first['actual_first_correction'];require(aj['whole_domain_source_box'] and aj['kernel']=='Dg=(-2*chi) K0 o K0 g5' and aj['baseline_K0_reentered'] is False,'local Dg source/word');require(bind['source_checks']['one_query_jet_shared_in_integration'] and bind['phase_shear_cutoff_identified_with_axis_potential'] is False,'same axis chi, not phase cutoff');require(comp['normalization_identity_scope']=='S-conjugated Achi extraction implies physical +2chi E51; remaining gauge/source identity remains conditional','physical P contract');contract=selected['actual_cutoff_complete_support_measure'];sc=iv(c,press['sc']);bits=tuple(contract['source_sc_exact_mpf_tuple']);require(sc._mpi_==(bits,bits) and sc.a>0 and w0['sc']==press['sc'] and aj['geometry']['sc_interval']==press['sc'] and contract['passed'] and contract['measure_is_Z_independent'] and '0<hb<1/2' in contract['inequality_proof'],'same chart and ordinary Z')
    w0src=raw[W0LEAF].decode('utf-8');firstsrc=raw[FIRSTLEAF].decode('utf-8');compsrc=raw[COMPLEAF].decode('utf-8');require("ROWS=('5','6','6','4','5','6')" in w0src and "'local_green_cell'" in w0src and 'Dg=(-2*chi) K0 o K0 g5' in firstsrc and 'P = +2 chi E51' in compsrc,'saved operator definitions')
    index_w0=w0['input_index'];index_first=first['source_index'];source_matches=[];source_count={}
    for order in (0,1):
        base=forcing['g'][str(order)]['5'];source_count[str(order)]=len(base);require(len(base)==16,'raw g5 grouping count')
        for i,source in enumerate(base):
            w_id=f'w0_g:{order}:5:{i}';f_id=f'collar_g5:Z{order}:{i}';wr=index_w0[w_id];fr=index_first[f_id];require(wr['receipt']==fr['receipt']==FORCING and ['result','output']+wr['path']==fr['path'] and wr['parent_source_ids']==fr['parent_source_ids']==source['origin']['source_ids'],'same g5 source pointer and parents');source_matches.append({'order':order,'index':i,'w0_id':w_id,'first_id':f_id,'source_path':fr['path']})
            for cell in range(4):
                term=w0['result']['local'][str(order)][str(cell)]['F1'][i];origin=term['origin'];require(w_id in origin['source_ids'] and origin['source_row']=='5' and origin['phase_cell']==cell and origin['operation']=='local_green_cell','local F1 row owner');require(term['powers'][0]==source['powers'][0]+2 and term['powers'][1:]==source['powers'][1:],'local source physical dx power')
    projection=phys['projection']['J1'];packet={'Z0':[],'Z1':[]};owners={};parent_index={};physical_index=phys['source_index']
    for order in (0,1):
        z='Z'+str(order);rows=projection[z];require(len(rows)==(32 if order==0 else 128),'physical Dg row count')
        for i,t in enumerate(rows):
            pid=f'physical_Dg:{z}:{i}';ids=t['source_ids'];local_ids=[x for x in ids if str(x).startswith('local_J1:')];lead_ids=[x for x in ids if str(x).startswith('lead_F:')];require(len(local_ids)==len(lead_ids)==1 and all(x in physical_index for x in ids),'indexed Dg and leading owner');direct=physical_index[local_ids[0]];require(direct['receipt']==FIRST and direct['path'][:4]==['result','output','actual_first_correction','response'] and any(str(x).startswith('collar_g5:') for x in direct['parent_source_ids']),'saved Dg parent path');parent_index[pid]={'receipt':PHYS,'path':['result','output','projection','J1',z,i],'parent_source_ids':ids,'Dg_source':direct};owners[pid]={'parent_id':pid,'owner':'Dg_physical','passed':True,'proved':True};packet[z].append({'parent_id':pid,'coefficient':iv(c,t['coefficient']),'powers':t['powers'],'source_ids':ids,'source_paths':['source_index:'+x for x in ids],'error_classes':t['error_classes'],'uncertainty':t['uncertainty'],'error_only':t['error_only'],'symmetric_error':t['symmetric_error']})
    cert={'passed':True,'owner_by_parent':owners,'same_g5':True,'same_K0':True,'same_chi':True,'P_sign':2,'physical_measure_same':True,'Z_independent_ordinary_derivative':True};result=leaf.build(c,packet,cert);require(result['schema']=='n1_local_p_dg_correlated_pair_v1' and result['passed'] and result['ledger_written'] is False,'pair scope');pairs=result['pairs'];counts={}
    for z in ('Z0','Z1'):
        counts[z]=len(pairs[z]);require(counts[z]==len(packet[z]),'unpaired parent')
        for pair in pairs[z]:
            a=pair['Dg']['coefficient'];b=pair['P']['coefficient'];require(a.a==-b.b and a.b==-b.a and pair['Dg']['powers']==pair['P']['powers'] and pair['symbolic_sum']['exact_zero'] and pair['symbolic_sum']['independent_boxed_interval_sum_performed'] is False,'lost exact correlated sign')
    return {'source_family':family,'source_logs':atoms,'fixture':proof,'source_matches':source_matches,'source_group_counts':source_count,'parent_index':parent_index,'input_source_index':physical_index,'pair_counts':counts,'pairs':enc(pairs),'operator_trees':{'W0_local_F':['K0','g5'],'Dg_F':['-2','K0','chi',['K0','g5']],'P_W0_local_F':['+2','K0','chi',['K0','g5']]},'source_paths_aligned':True,'operator_identity_scope':'conditional same underlying grouped g5 and axis chi, same K0 and physical chart','old_errors_removed':False,'independent_interval_sum_performed':False,'appended_to_ledger':False,'all_physical_compensation_complete':False,'full_m15_computed':False,'terminal_closure_proved':False}
def main():
    h,s,b=preflight();entry=save(ENTRY,{'component':COMPONENT,'source_sha256':s,'brief_sha256':b,'input_hashes':h,'entered':True});out={'component':COMPONENT,'source_sha256':s,'entry_sha256':entry,'status':'failed'}
    try:out.update(status='returned',output=consume(h),error=None)
    except BaseException as exc:out['error']={'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()}
    digest=save(RESULT,out);print(json.dumps({'component':COMPONENT,'status':out['status'],'result_sha256':digest}));return 0 if out['status']=='returned' else 1
if __name__=='__main__':raise SystemExit(main())
