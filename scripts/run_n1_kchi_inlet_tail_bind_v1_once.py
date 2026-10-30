"""Bind the saved direct W0 Kchi tail to its same-source six-field baseline."""
from __future__ import annotations
import gzip,hashlib,json,traceback
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];BRIEF=ROOT.parent.parent/'RESUME_BRIEF_2026_10_16.json';COMPONENT='n1_kchi_inlet_tail_bind_v1';DIR=ROOT/'tmp'/COMPONENT
SOURCE,ENTRY,RESULT=(DIR/n for n in ('source.json','entered.json','result.json'))
BASE='experiments/root_st073/n1_core_six_inlet_v2_receipt.json.gz';GAUGE='experiments/root_st073/n1_core_gauge_exit_v2_receipt.json.gz';KCHI='experiments/root_st073/n1_core_kchi_first_v1_receipt.json.gz';RUNNER='scripts/run_n1_kchi_inlet_tail_bind_v1_once.py'
REQUIRED=(BASE,GAUGE,KCHI,RUNNER)
FIXED={BASE:'59c0ff88652d83d7b92ba9c3a564ee90eac5ae794220b8b858c112936b76342d',GAUGE:'4c1e49a83dbc7c13e5d9f09f4601ce34b26b487405a2d14a2f066ee2dfe27921',KCHI:'b0d654203e8c4acabdd0b61ca7b79c1b7050314b0c003a0c2f8757bde94973a0'}
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
def consume(hashes):
    raw={}
    for p in REQUIRED:
        b=(ROOT/p).read_bytes();require(sha(b)==hashes[p],'changed after entry '+p);raw[p]=b
    base=json.loads(gzip.decompress(raw[BASE]));gauge=json.loads(gzip.decompress(raw[GAUGE]));kchi=json.loads(gzip.decompress(raw[KCHI]));require(all(r['result']['status']=='returned' for r in (base,gauge,kchi)),'input status')
    bp={v['path']:v['sha256'] for v in base['source']['files']};gp={v['path']:v['sha256'] for v in gauge['source']['files']};require(bp[GAUGE]==hashes[GAUGE] and gp[KCHI]==hashes[KCHI],'authenticated receipt chain')
    b=base['result']['output'];g=gauge['result']['output'];k=kchi['result']['output'];require(b['source_family']==g['source_family']==k['source_family'] and b['source_logs']==g['source_logs']==k['source_logs'],'same source family and formal logs')
    require(k['Kchi_tail_is_uncertainty_not_signed_field'] is True and k['tail_includes_feedback'] is False and k['source_boxes_not_point_values'] is True,'tail scope')
    require(b['shared_global_source_reused_across_cells'] is True and b['feedback_word_tail_included'] is False and b['complete_six_component_inlet'] is False,'baseline scope')
    direct={};counts={}
    for field in ('F1','d_x_F1'):
        direct[field]={};counts[field]={}
        for order in ('0','1'):
            terms=k['W0_separate_symmetric_kernel_tail'][field][order];require(isinstance(terms,list) and len(terms)>0,'nonempty tail')
            out=[]
            for i,term in enumerate(terms):
                iv=term['coefficient'];lo=iv['lower_exact_mpf_tuple'];hi=iv['upper_exact_mpf_tuple'];require(iv['directed'] is True and len(lo)==len(hi)==4 and lo[1:]==hi[1:] and ((lo[0],hi[0])==(1,0) or (lo[1]==0 and hi[1]==0)),'directed symmetric tail')
                require(len(term['powers'])==4 and bool(term['source_ids']) and 'Kchi_j_ge_2_kernel_tail' in term['error_classes'] and term['uncertainty'] is True,'tail provenance')
                copy=dict(term);copy.update(error_only=True,source_record_path=['result','output','W0_separate_symmetric_kernel_tail',field,order,i],source_receipt=KCHI,shared_global_source=True);out.append(copy)
            direct[field][order]=out;counts[field][order]=len(out)
    return {'source_family':b['source_family'],'source_logs':b['source_logs'],'direct_W0_kernel_tail':direct,'tail_term_counts':counts,
      'baseline_reference':{'receipt':BASE,'path':['result','output','baseline_G0g_fixed_exit']},'first_Dg_reference':{'receipt':GAUGE,'path':['result','output','saved_W0_first_correction_fixed_exit']},
      'first_feedback_kernel_tail_references':{'receipt':GAUGE,'paths':[['result','output','inner_W0_kernel_tail_fixed_exit'],['result','output','outer_kernel_tail_fixed_exit'],['result','output','potential_compensation_kernel_tail_fixed_exit']]},
      'direct_tail_not_reintegrated':True,'direct_tail_is_symmetric_error_only':True,'formal_source_atoms_materialized':False,'merged_with_signed_terms':False,'feedback_word_tail_included':False,'complete_six_component_inlet':False,'appended_to_m15_ledger':False,'terminal_closure_proved':False}
def main():
    h,s,b=preflight();entry=save(ENTRY,{'component':COMPONENT,'source_sha256':s,'brief_sha256':b,'input_hashes':h,'entered':True});out={'component':COMPONENT,'source_sha256':s,'entry_sha256':entry,'status':'failed'}
    try:out.update(status='returned',output=consume(h),error=None)
    except BaseException as exc:out['error']={'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()}
    digest=save(RESULT,out);print(json.dumps({'component':COMPONENT,'status':out['status'],'result_sha256':digest}));return 0 if out['status']=='returned' else 1
if __name__=='__main__':raise SystemExit(main())
