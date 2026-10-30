"""One-entry finite W0 axial and swirl m15 enclosures by phase cell."""
from __future__ import annotations
import gzip,hashlib,importlib.util,json,sys,traceback
from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext
ROOT=Path(__file__).resolve().parents[1];BRIEF=ROOT.parent.parent/'RESUME_BRIEF_2026_10_16.json';COMPONENT='n1_w0_physical_m15_v1';DIR=ROOT/'tmp'/COMPONENT
SOURCE,ENTRY,RESULT=(DIR/n for n in ('source.json','entered.json','result.json'))
W0='experiments/root_st073/n1_w0_collar_sectors_v1_receipt.json.gz';RADIAL='experiments/root_st073/n1_actual_collar_radial_v1_receipt.json.gz';SWIRL='experiments/root_st073/n1_first_word_swirl_m15_v1_receipt.json.gz';AXIAL='experiments/root_st073/n1_first_word_axial_m15_v1_receipt.json.gz';CORRECTION='experiments/root_st073/n1_first_word_swirl_m15_v1_provenance_correction.json';TAIL='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_selected_words.json'
CUT_SOURCE='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_inner_cutoff_functions.py';GAMMA='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_Gamma_future_functions.py';CUTOFF='src/openai_ns_reconstruction/original_n1_cutoff_numeric.py'
LEAF='src/openai_ns_reconstruction/n1_w0_physical_m15_v1.py';FIXTURE='src/openai_ns_reconstruction/n1_w0_physical_m15_fixture_v1.py';RUNNER='scripts/run_n1_w0_physical_m15_v1_once.py'
REQUIRED=(W0,RADIAL,SWIRL,AXIAL,CORRECTION,TAIL,CUT_SOURCE,GAMMA,CUTOFF,LEAF,FIXTURE,RUNNER)
FIXED={W0:'4229c29a2a21207f8b8da3ba53d70bbf6f50e09eb67110bd3b7cb30fa6577c06',RADIAL:'27d0e2a97568c50ff62f373cbdf04dea262b4980392d7ac4b13b9c7db21752cf',SWIRL:'a626e71fa2b0ddf40b79ce992016da2af8ac9955db3f17337029f3c86c2513b4',AXIAL:'02acdcb22d9b0f03b4201f502388814d54d49937124ebb2288e66b9706770ab5',CORRECTION:'f2133bf5824b2f9636b430a33a2ba1d19fbf391a7f90a81de0c8269d1fe638c4',TAIL:'0e04c87a011bfec74e0ceb0f822774a03495a08a839aaf4e378d21949736cf1e',CUT_SOURCE:'b9f8b854b51b9eb0eda338503c8d25e6a1b69307da54b75bf40c1bd212ff1b29',GAMMA:'38d3a8e5abd964c42b135410485f6ad58dfe9f7656035b4b5c7b1ed48c2ebb01',CUTOFF:'7a2004b378ee495c6c053d04a8394406937f3a6fca235a721138f1ce4661c731'}
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
    require(Path(ck['path']).resolve()==SOURCE.resolve() and ck['sha256']==sha(src),'source checkpoint');m=json.loads(src);require(m['component']==COMPONENT and m['one_entry_only'] is True,'manifest');hashes={}
    for row in m['files']:
        p=row['path'];require(p in REQUIRED and p not in hashes,'unexpected input');h=sha((ROOT/p).read_bytes());require(h==row['sha256'] and (p not in FIXED or h==FIXED[p]),'changed '+p);hashes[p]=h
    require(set(hashes)==set(REQUIRED),'missing input');return hashes,sha(src),sha(b)
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
    for p in REQUIRED:
        b=(ROOT/p).read_bytes();require(sha(b)==hashes[p],'changed after entry '+p);raw[p]=b
    w0=json.loads(gzip.decompress(raw[W0]));radial=json.loads(gzip.decompress(raw[RADIAL]));swirl=json.loads(gzip.decompress(raw[SWIRL]));axial=json.loads(gzip.decompress(raw[AXIAL]));correction=json.loads(raw[CORRECTION]);selected=json.loads(raw[TAIL])
    for receipt in (w0,radial,swirl,axial):require(receipt['result']['status']=='returned','accepted inputs')
    a=w0['result']['output'];r=radial['result']['output'];f=swirl['result']['output'];u=axial['result']['output'];family=a['source_family'];atoms=a['source_logs']
    require(r['source_family']==f['source_family']==u['source_family']==selected['source_family']==family and r['source_logs']==f['source_logs']==u['source_logs']==atoms and a['sc']==f['sc']==u['sc'],'same family, atoms and phase')
    require(a['prior_axial_receipt']==AXIAL and a['prior_swirl_receipt']==SWIRL and a['allword_tail_receipt']==TAIL and f['first_word_W1_only'] and u['first_word_W1_only'] and f['appended_to_cumulative_m15'] is False and u['appended_to_cumulative_m15'] is False,'independent old contributions')
    require(correction['receipt_sha256']==hashes[SWIRL] and correction['numerical_result_changed'] is False and correction['input_id']=='F0_connection','F0 source path correction')
    measure=selected['actual_cutoff_complete_support_measure'];require(measure['passed'] and measure['source_family']==family and measure['measure_is_Z_independent'] and '0<hb<1/2' in measure['inequality_proof'],'geometry source')
    src=raw[CUT_SOURCE].decode('utf8');gamma=raw[GAMMA].decode('utf8');cutoff=raw[CUTOFF].decode('utf8')
    require("pressure=a.add(a.scale(a.mul(f,F),2),system['N1p'])" in src and 'a.neg(a.scale(pressure,g.mul(eps2,r)))' in src and "definition='1-sigma(4q-2)',only_f_and_u_multiplied=True" in cutoff and 'exp(-1/x^2)/(exp(-1/x^2)+exp(-1/(1-x)^2)) inside' in gamma,'original density and exact velocity cutoff')
    c=MPIntervalContext();c.prec=1800;sc=iv(c,a['sc']);bits=tuple(measure['source_sc_exact_mpf_tuple']);require(sc._mpi_==(bits,bits),'same exact sc')
    markers=(('Cauchy_phase','symmetric_certified_error_not_field_value'),('Cauchy_phase','cauchy_phase_variation_error'),('mean_phase','U_phase_error_history_integral'),('cutoff_value','minus_sigma_factor'),('cutoff_value','hb_sigma_factor'),('cutoff_derivative','sigma_prime_factor'))
    index={};leading_U={};leading_F={};word={}
    for order in (0,1):
        key='(0, '+str(order)+')';leading_U[order]=[];leading_F[order]=[];word[order]={}
        for i,t in enumerate(r['actual_radial_sector_enclosures']['U'][key]):
            sid=f'w0_U:{order}:{i}';origin_json=json.dumps(t['origin'],sort_keys=True,separators=(',',':'));classes=sorted({label for label,needle in markers if needle in origin_json});index[sid]={'receipt':RADIAL,'path':['actual_radial_sector_enclosures','U',key,i],'origin_sha256':sha(origin_json.encode())}
            leading_U[order].append({'powers':t['powers'],'coefficient':iv(c,t['coefficient']),'origin':{'source_ids':[sid],'uncertainty':bool(classes),'error_classes':classes}})
        for i,t in enumerate(f['F0_connection']['leading'][str(order)]):
            sid=f'w0_F:{order}:{i}';index[sid]={'receipt':SWIRL,'path':['F0_connection','leading',str(order),i],'parent_source_ids':t['origin']['source_ids'],'connection_path_correction':CORRECTION}
            leading_F[order].append({'powers':t['powers'],'coefficient':iv(c,t['coefficient']),'origin':{'source_ids':[sid],'uncertainty':t['origin']['uncertainty'],'error_classes':t['origin']['error_classes']}})
        for field in ('Uz1','F1'):
            word[order][field]={'core':[],'local':{}}
            for i,t in enumerate(a['result']['core_transport'][str(order)][field]):
                sid=f'w0_word:{order}:{field}:core:{i}';index[sid]={'receipt':W0,'path':['result','core_transport',str(order),field,i],'parent_source_ids':t['origin']['source_ids']}
                word[order][field]['core'].append({'powers':t['powers'],'coefficient':iv(c,t['coefficient']),'origin':{'source_ids':[sid],'uncertainty':t['origin']['uncertainty'],'error_classes':t['origin']['error_classes']}})
            for cell in range(4):
                word[order][field]['local'][cell]=[]
                for i,t in enumerate(a['result']['local'][str(order)][str(cell)][field]):
                    sid=f'w0_word:{order}:{field}:{cell}:{i}';index[sid]={'receipt':W0,'path':['result','local',str(order),str(cell),field,i],'parent_source_ids':t['origin']['source_ids'],'zero_velocity_cutoff':cell==3}
                    word[order][field]['local'][cell].append({'powers':t['powers'],'coefficient':iv(c,t['coefficient']),'origin':{'source_ids':[sid],'uncertainty':t['origin']['uncertainty'],'error_classes':t['origin']['error_classes']}})
    leaf=module(LEAF,'n1_w0_physical_m15_v1');fixture=module(FIXTURE,'n1_w0_physical_m15_fixture_v1');proof=fixture.run(leaf,c);require(proof['passed'],'fixture failed '+str([k for k,v in proof['checks'].items() if not v]))
    cert={'hb_strictly_between_zero_and_half':True,'sc_phase_Z_independent':True,'hb_epsilon_Z_independent':True,'original_cutoff_mass_5_over_8':True,'original_cutoff_mass_by_cell':True}
    result=leaf.build(c,sc=sc,leading_U=leading_U,leading_F=leading_F,w0=word,source_family=family,source_atoms=atoms,certificate=cert)
    seen=set();classes=set()
    for channel in ('axial','swirl'):
        for branch in result[channel].values():
            collections=[branch['core']]+list(branch['local'].values());require(set(branch['local'])=={0,1,2},'zero-cutoff cell handling')
            for terms in collections:
                for term in terms:
                    o=term['origin'];require(all(k in index for k in o['source_ids']),'unmapped product source');seen.update(o['source_ids']);classes.update(o['error_classes'])
    expected={sid for sid,item in index.items() if item.get('zero_velocity_cutoff') is not True};require(seen==expected and classes,'lost input provenance')
    return {'source_family':family,'source_logs':atoms,'sc':a['sc'],'domain':a['domain'],'fixture':proof,'input_index':index,'referenced_input_ids':len(seen),'output_error_classes':sorted(classes),'result':enc(result),
      'exact_eta_cell_masses':['1/4','1/4','1/8','0'],'zero_cell_applies_only_to_cut_velocities':True,'F0_connection_reused_once':True,'W0_only':True,'W1_components_separate':[AXIAL,SWIRL],
      'whole_cell_uniform_boxes_not_point_values':True,'appended_to_cumulative_m15':False,'allword_tail_added':False,'complete_m15_computed':False,'terminal_closure_proved':False,'allword_tail_receipt':TAIL}
def main():
    h,s,b=preflight();entry=save(ENTRY,{'component':COMPONENT,'source_sha256':s,'brief_sha256':b,'input_hashes':h,'entered':True});out={'component':COMPONENT,'source_sha256':s,'entry_sha256':entry,'status':'failed'}
    try:out.update(status='returned',output=consume(h),error=None)
    except BaseException as exc:out['error']={'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()}
    digest=save(RESULT,out);print(json.dumps({'component':COMPONENT,'status':out['status'],'result_sha256':digest}));return 0 if out['status']=='returned' else 1
if __name__=='__main__':raise SystemExit(main())
