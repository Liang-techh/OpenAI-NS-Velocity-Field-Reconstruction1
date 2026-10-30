"""One-entry leading-pressure increment and finite velocity-prefix binding."""
from __future__ import annotations
import gzip,hashlib,importlib.util,json,sys,traceback
from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext
ROOT=Path(__file__).resolve().parents[1];BRIEF=ROOT.parent.parent/'RESUME_BRIEF_2026_10_16.json';COMPONENT='n1_first_collar_g4_pressure_v1';DIR=ROOT/'tmp'/COMPONENT
SOURCE,ENTRY,RESULT=(DIR/n for n in ('source.json','entered.json','result.json'))
W0='experiments/root_st073/n1_w0_physical_m15_v1_receipt.json.gz';W0SECT='experiments/root_st073/n1_w0_collar_sectors_v1_receipt.json.gz';W1U='experiments/root_st073/n1_first_word_axial_m15_v1_receipt.json.gz';W1F='experiments/root_st073/n1_first_word_swirl_m15_v1_receipt.json.gz'
FORCING='experiments/root_st073/n1_actual_collar_forcing_v1_receipt.json.gz';RADIAL='experiments/root_st073/n1_actual_collar_radial_v1_receipt.json.gz';TAIL='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_selected_words.json'
OP='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_inner_operator.py';CUT='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_inner_cutoff_functions.py';RADIAL_SOURCE='src/openai_ns_reconstruction/n1_actual_collar_radial.py';FORCING_SOURCE='src/openai_ns_reconstruction/n1_actual_collar_forcing_compact.py'
LEAF='src/openai_ns_reconstruction/n1_first_collar_g4_pressure_v1.py';FIXTURE='src/openai_ns_reconstruction/n1_first_collar_g4_pressure_fixture_v1.py';RUNNER='scripts/run_n1_first_collar_g4_pressure_v1_once.py'
REQUIRED=(W0,W0SECT,W1U,W1F,FORCING,RADIAL,TAIL,OP,CUT,RADIAL_SOURCE,FORCING_SOURCE,LEAF,FIXTURE,RUNNER)
FIXED={W0:'f827579eb137c936eb0a2c00fbba873cccf6e8204d00e03d4c478531d0b770ab',W0SECT:'4229c29a2a21207f8b8da3ba53d70bbf6f50e09eb67110bd3b7cb30fa6577c06',W1U:'02acdcb22d9b0f03b4201f502388814d54d49937124ebb2288e66b9706770ab5',W1F:'a626e71fa2b0ddf40b79ce992016da2af8ac9955db3f17337029f3c86c2513b4',FORCING:'96e8a3b64140fb2c111e8117951e96674dfa3c2bc7ff41ebe34448dccbbd8be2',RADIAL:'27d0e2a97568c50ff62f373cbdf04dea262b4980392d7ac4b13b9c7db21752cf',TAIL:'0e04c87a011bfec74e0ceb0f822774a03495a08a839aaf4e378d21949736cf1e',OP:'268217c6b33bc717472742fa64d76a853c6760a956abfad92ac45bcf385df014',CUT:'b9f8b854b51b9eb0eda338503c8d25e6a1b69307da54b75bf40c1bd212ff1b29'}
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
    w0=json.loads(gzip.decompress(raw[W0]));w0sect=json.loads(gzip.decompress(raw[W0SECT]));w1u=json.loads(gzip.decompress(raw[W1U]));w1f=json.loads(gzip.decompress(raw[W1F]));forcing=json.loads(gzip.decompress(raw[FORCING]));radial=json.loads(gzip.decompress(raw[RADIAL]));selected=json.loads(raw[TAIL])
    for r in (w0,w0sect,w1u,w1f,forcing,radial):require(r['result']['status']=='returned','accepted receipt')
    a=w0['result']['output'];s=w0sect['result']['output'];u=w1u['result']['output'];v=w1f['result']['output'];f=forcing['result']['output'];r=radial['result']['output'];family=a['source_family'];atoms=a['source_logs']
    require(all(x['source_family']==family and x['source_logs']==atoms for x in (s,u,v,f,r)) and selected['source_family']==family,'same source family and atom records')
    require(a['sc']==s['sc']==u['sc']==v['sc'] and a['W0_only'] and a['W1_components_separate']==[W1U,W1F] and all(x['appended_to_cumulative_m15'] is False for x in (a,u,v)),'same geometry and nonduplicating finite words')
    require({row['path']:row['sha256'] for row in w0['source']['files']}[W0SECT]==hashes[W0SECT],'same W0 sector receipt')
    require(s['source_forcing_receipt']==FORCING and f['source_index_receipt_sha256']==hashes[RADIAL] and f['prior_radial_result_sha256']==radial['result_sha256'] and f['same_source_radial_and_core_receipts_verified'] and r['inherited_interval_errors_preserved_in_sector_origins'],'same full radial forcing chain')
    meta=r['actual_radial_metadata'];require('chi_r*U_core_rho' in meta['radial_formulas']['U_rhorho'] and 'sigma_prime' in meta['chi_r_formula'],'uncut leading cutoff radial derivative retained')
    op=raw[OP].decode('utf8');cut=raw[CUT].decode('utf8');rs=raw[RADIAL_SOURCE].decode('utf8');fs=raw[FORCING_SOURCE].decode('utf8')
    require("g={'4':eps(forcing['N1p']*(2*x))" in op and "pressure=a.add(a.scale(a.mul(f,F),2),system['N1p'])" in cut and 'a.neg(a.scale(pressure,g.mul(eps2,r)))' in cut and 'chi_r*U_core_rho' in rs and 'g4_power0, g4_power1 = 2 * x_formal * P, -x_formal * S' in fs,'original g4 identity, physical sign and full source')
    measure=selected['actual_cutoff_complete_support_measure'];require(measure['passed'] and measure['source_family']==family and measure['measure_is_Z_independent'] and '0<hb<1/2' in measure['inequality_proof'],'geometry source')
    c=MPIntervalContext();c.prec=1800;sc=iv(c,a['sc']);bits=tuple(measure['source_sc_exact_mpf_tuple']);require(sc._mpi_==(bits,bits),'same exact sc')
    index={};g4={}
    for order in (0,1):
        g4[order]=[]
        for i,t in enumerate(f['g'][str(order)]['4']):
            sid=f'pressure_g4:{order}:{i}';origin=t['origin'];index[sid]={'receipt':FORCING,'path':['g',str(order),'4',i],'parent_source_ids':origin['source_ids']}
            g4[order].append({'powers':t['powers'],'coefficient':iv(c,t['coefficient']),'origin':{'source_ids':[sid],'uncertainty':origin['uncertainty'],'error_classes':origin['error_classes']}})
    leaf=module(LEAF,'n1_first_collar_g4_pressure_v1');fixture=module(FIXTURE,'n1_first_collar_g4_pressure_fixture_v1');proof=fixture.run(leaf,c);require(proof['passed'],'fixture failed '+str([k for k,v in proof['checks'].items() if not v]))
    cert={'hb_strictly_between_zero_and_half':True,'sc_phase_Z_independent':True,'hb_epsilon_Z_independent':True,'g4_equals_two_epsilon_x_N1p':True}
    pressure=leaf.build(c,sc=sc,g4=g4,source_family=family,source_atoms=atoms,certificate=cert)
    for order in (0,1):
        require(len(pressure['pressure'][order])==len(g4[order]),'g4 pressure term duplicated or lost')
        for original,term in zip(g4[order],pressure['pressure'][order]):
            powers=list(original['powers']);powers[0]+=1;powers[3]-=1
            require(list(term['powers'])==powers and term['origin']['source_ids']==original['origin']['source_ids'] and term['origin']['error_classes']==original['origin']['error_classes'] and term['origin']['uncertainty']==original['origin']['uncertainty'],'g4 source or error changed')
    require(pressure['velocity_cutoff_applied'] is False and pressure['explicit_pressure_endpoint_added'] is False and pressure['allword_tail_added'] is False,'uncut pressure density contract')
    prefix={};counts={}
    for order in (0,1):
        prefix[order]={}
        for channel,word1,word1path in (('axial',u,W1U),('swirl',v,W1F)):
            terms=[];paths=[(W0,[channel,str(order),'core'],a['result'][channel][str(order)]['core'])]
            paths.extend((W0,[channel,str(order),'local',str(cell)],a['result'][channel][str(order)]['local'][str(cell)]) for cell in range(3))
            paths.extend((word1path,[channel,str(order),branch],word1['result'][channel][str(order)][branch]) for branch in ('core','local'))
            for pathnum,(receipt,path,items) in enumerate(paths):
                for i,t in enumerate(items):
                    sid=f'prefix:{order}:{channel}:{pathnum}:{i}';index[sid]={'receipt':receipt,'path':['result']+path+[i],'parent_source_ids':t['origin']['source_ids']}
                    o=dict(t['origin']);o['source_ids']=list(o['source_ids'])+[sid];terms.append({'powers':t['powers'],'coefficient':iv(c,t['coefficient']),'origin':o})
            prefix[order][channel]=terms;counts[str(order)+':'+channel]=len(terms)
        prefix[order]['pressure']=pressure['pressure'][order]
    return {'source_family':family,'source_logs':atoms,'sc':a['sc'],'domain':a['domain'],'fixture':proof,'input_index':index,'pressure':enc(pressure),'finite_prefix_channels':enc(prefix),'velocity_group_counts':counts,
      'g4_full_radial_derivative_including_cutoff_prime':True,'pressure_uncut_by_eta':True,'pressure_direct_density_no_endpoint_added':True,
      'velocity_prefix_is_W0_plus_W1_once':True,'finite_prefix_is_unmerged_sector_lists':True,'source_intervals_not_point_values':True,
      'core_m15_inlet_added':False,'allword_tail_added':False,'appended_to_cumulative_m15':False,'complete_m15_computed':False,'terminal_closure_proved':False,
      'W0_receipt':W0,'W1_axial_receipt':W1U,'W1_swirl_receipt':W1F,'forcing_receipt':FORCING,'allword_tail_receipt':TAIL}
def main():
    h,s,b=preflight();entry=save(ENTRY,{'component':COMPONENT,'source_sha256':s,'brief_sha256':b,'input_hashes':h,'entered':True});out={'component':COMPONENT,'source_sha256':s,'entry_sha256':entry,'status':'failed'}
    try:out.update(status='returned',output=consume(h),error=None)
    except BaseException as exc:out['error']={'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()}
    digest=save(RESULT,out);print(json.dumps({'component':COMPONENT,'status':out['status'],'result_sha256':digest}));return 0 if out['status']=='returned' else 1
if __name__=='__main__':raise SystemExit(main())
