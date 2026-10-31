"""One-entry fixed E00 physical velocity-m15 projection by q cell."""
from __future__ import annotations
import gzip,hashlib,importlib.util,json,sys,traceback
from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext
ROOT=Path(__file__).resolve().parents[1];BRIEF=ROOT.parent.parent/'RESUME_BRIEF_2026_10_16.json';COMPONENT='n1_e00_physical_m15_cells_v2';DIR=ROOT/'tmp'/COMPONENT
SOURCE,ENTRY,RESULT=(DIR/n for n in ('source.json','entered.json','result.json'))
SEAM='experiments/root_st073/n1_e00_collar_fixed_seam_v4_receipt.json.gz';TAIL='experiments/root_st073/n1_e00_collar_j2_volterra_v2_receipt.json.gz';RADIAL='experiments/root_st073/n1_actual_collar_radial_v1_receipt.json.gz';SWIRL='experiments/root_st073/n1_first_word_swirl_m15_v1_receipt.json.gz';CORRECTION='experiments/root_st073/n1_first_word_swirl_m15_v1_provenance_correction.json';SELECT='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_selected_words.json';PREVIOUS='experiments/root_st073/n1_collar_qp_k0_m15_cells_v1_receipt.json.gz'
LEAF='src/openai_ns_reconstruction/n1_e00_physical_m15_cells_v2.py';FIXTURE='src/openai_ns_reconstruction/n1_e00_physical_m15_cells_fixture_v2.py';RUNNER='scripts/run_n1_e00_physical_m15_cells_v2_once.py'
FAILURE='experiments/root_st073/n1_e00_physical_m15_cells_v1_failure.json.gz'
REQUIRED=(SEAM,TAIL,RADIAL,SWIRL,CORRECTION,SELECT,PREVIOUS,FAILURE,LEAF,FIXTURE,RUNNER)
FIXED={FAILURE:'48f0382f6a259b5756d325546a2c2674f82a3fc09934e2aa453a258ea9816659',SEAM:'6f8cc4c4fb4baedb9d14274c2f8ddc92f4ae536766d3d911a807d6f626539766',TAIL:'94ca27da585d7e3d252ba94bac1453f037f9601e5fd3d14a9b1b30d2de4f11c5',RADIAL:'27d0e2a97568c50ff62f373cbdf04dea262b4980392d7ac4b13b9c7db21752cf',SWIRL:'a626e71fa2b0ddf40b79ce992016da2af8ac9955db3f17337029f3c86c2513b4',CORRECTION:'f2133bf5824b2f9636b430a33a2ba1d19fbf391a7f90a81de0c8269d1fe638c4',SELECT:'0e04c87a011bfec74e0ceb0f822774a03495a08a839aaf4e378d21949736cf1e',PREVIOUS:'48409ed314e95fe0f3117f6c27ed74b5aa8b56822c810ba5050d2b000b08e5c2'}
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
    c=MPIntervalContext();c.prec=1800;leaf=module(LEAF,'n1_e00_physical_m15_cells_v2');fixture=module(FIXTURE,'n1_e00_physical_m15_cells_fixture_v2');proof=fixture.run_fixture(leaf,c);require(proof['passed'],'fixture failed '+json.dumps(proof,sort_keys=True,separators=(',',':')))
    receipts={p:json.loads(gzip.decompress(raw[p])) for p in (SEAM,TAIL,RADIAL,SWIRL,PREVIOUS)};require(all(r['result']['status']=='returned' for r in receipts.values()),'receipt status');seam,tail,rad,sw,prev=(receipts[p]['result']['output'] for p in (SEAM,TAIL,RADIAL,SWIRL,PREVIOUS));selected=json.loads(raw[SELECT]);correction=json.loads(raw[CORRECTION]);family=seam['source_family'];atoms=seam['source_logs'];require(all(x['source_family']==family and x['source_logs']==atoms for x in (tail,rad,sw,prev)) and selected['source_family']==family,'common source family and logs');pins={v['path']:v['sha256'] for v in receipts[PREVIOUS]['source']['files']};require(all(pins[p]==hashes[p] for p in (RADIAL,SWIRL,CORRECTION,SELECT)),'leading and cutoff authenticated chain');tailpins={v['path']:v['sha256'] for v in receipts[TAIL]['source']['files']};require(tailpins[SEAM]==hashes[SEAM] and tail['old_errors_removed'] is False and tail['only_symmetric_remainder_bounds'] is True,'fixed inlet tail chain');require(correction['receipt_sha256']==hashes[SWIRL] and correction['numerical_result_changed'] is False and correction['input_id']=='F0_connection','amplitude path correction')
    v=seam['actual_E00_seam_and_collar'];contract=selected['actual_cutoff_complete_support_measure'];sc=iv(c,v['geometry']['sc']);bits=tuple(contract['source_sc_exact_mpf_tuple']);require(sc._mpi_==(bits,bits) and sc.a>0 and sw['sc']==v['geometry']['sc'] and contract['passed'] and contract['measure_is_Z_independent'] and '0<hb<1/2' in contract['inequality_proof'],'physical geometry and sc');index=dict(seam['source_index']);require(all(index[k]==record for k,record in tail['source_index'].items()),'tail source paths identical to fixed inlet');inputs={'K0':{'U':{},'F':{}},'J1':{'F':{}},'J2':{'F':{}}};input_counts={};expected_q=((0,'0.25'),('0.25','0.5'),('0.5',1))
    def source_term(t):
        ids=t['source_ids'];require(ids and all(i in index for i in ids),'indexed field sector');return {'coefficient':iv(c,t['coefficient']),'powers':t['powers'],'source_ids':ids,'source_paths':['source_index:'+i for i in ids],'error_classes':t.get('error_classes',[]),'uncertainty':t.get('uncertainty',False),'error_only':t.get('error_only',False),'symmetric_error':t.get('symmetric_error',False),'field_branch':t.get('branch')}
    for i,(lo,hi) in enumerate(expected_q):
        scell=v['q_cells'][i];tcell=tail['tail_cells'][i+1];a=iv(c,scell['q_interval']);b=iv(c,tcell['q_interval']);require(a._mpi_==b._mpi_ and a.a==c.mpf(lo) and a.b==c.mpf(hi),'matching whole q cells')
        for part,channel,stage,names in (('K0','U','k0_transfer',('U_from_Uc','U_from_Uxc')),('K0','F','k0_transfer',('F_from_Fc','F_from_Fxc')),('J1','F','first_potential',('F_from_Fc','F_from_Fxc')),('J2','F','tail',('F_from_Fc','F_from_Fxc'))):
            branches=tcell['branches'] if stage=='tail' else scell[stage]['branches'];rows={}
            for z in ('Z0','Z1'):
                items=[source_term(t) for name in names for t in branches[name][z]];require(items and (part!='J2' or all(t['error_only'] and t['symmetric_error'] for t in items)),'missing field or tail status');rows[z]=items;input_counts[f'{part}:{channel}:{i}:{z}']=len(items)
            inputs[part][channel][i]=rows
    markers=(('Cauchy_phase','symmetric_certified_error_not_field_value'),('Cauchy_phase','cauchy_phase_variation_error'),('mean_phase','U_phase_error_history_integral'),('cutoff_value','minus_sigma_factor'),('cutoff_value','hb_sigma_factor'),('cutoff_derivative','sigma_prime_factor'));lead={'U':{},'F':{}}
    for order in (0,1):
        z='Z'+str(order);u=[];f=[];key='(0, '+str(order)+')'
        for i,t in enumerate(rad['actual_radial_sector_enclosures']['U'][key]):
            sid=f'lead_U:{order}:{i}';origin=json.dumps(t['origin'],sort_keys=True,separators=(',',':'));classes=sorted({label for label,needle in markers if needle in origin});index[sid]={'receipt':RADIAL,'path':['result','output','actual_radial_sector_enclosures','U',key,i],'origin_sha256':sha(origin.encode())};u.append({'powers':t['powers'],'coefficient':iv(c,t['coefficient']),'source_ids':[sid],'source_paths':['source_index:'+sid],'uncertainty':bool(classes),'error_classes':classes})
        for i,t in enumerate(sw['F0_connection']['leading'][str(order)]):
            sid=f'lead_F:{order}:{i}';index[sid]={'receipt':SWIRL,'path':['result','output','F0_connection','leading',str(order),i],'parent_source_ids':t['origin']['source_ids'],'connection_path_correction':CORRECTION};f.append({'powers':t['powers'],'coefficient':iv(c,t['coefficient']),'source_ids':[sid],'source_paths':['source_index:'+sid],'uncertainty':t['origin']['uncertainty'],'error_classes':t['origin']['error_classes']})
        require(u and f,'whole-phase leading boxes');require(all(prev['source_index'][t['source_ids'][0]]==index[t['source_ids'][0]] for t in u+f),'same leading source paths as authenticated physical collar');
        for i in range(3):
            lead['U'].setdefault(i,{})[z]=u;lead['F'].setdefault(i,{})[z]=f
    cert={'hb_strictly_between_zero_and_half':True,'sc_phase_Z_independent':True,'hb_epsilon_Z_independent':True,'original_cutoff_mass_5_over_8':True,'original_cutoff_mass_by_cell':True};result=leaf.build(c,sc=sc,leading_U=lead['U'],leading_F=lead['F'],field_sectors=inputs,source_family=family,source_atoms=v['source_atoms'],certificate=cert);require(result['schema']=='n1_e00_physical_m15_cells_v2' and result['appended_to_ledger'] is False,'scope');parts={};counts={};seen=set()
    for part in ('K0','J1','J2'):
        channels=result['parts'][part]['channels'];parts[part]={}
        for channel in ('axial','swirl'):
            parts[part][channel]={}
            for order in (0,1):
                parts[part][channel]['Z'+str(order)]={}
                for cell in range(4):
                    rows=channels[channel][order][cell] if channel in channels else [];
                    for t in rows:
                        require(all(sid in index for sid in t['source_ids']),'unindexed projected source');require(part!='J2' or (t['error_only'] and t['symmetric_error']),'j2 is symmetric error');seen.update(t['source_ids'])
                    require(cell!=3 or not rows,'zero cutoff tail');parts[part][channel]['Z'+str(order)][str(cell)]=rows;counts[f'{part}:{channel}:Z{order}:{cell}']=len(rows)
    require(all(not parts[p]['axial'][z][str(i)] for p in ('J1','J2') for z in ('Z0','Z1') for i in range(4)),'no axial potential terms');return {'source_family':family,'source_logs':atoms,'source_index':index,'used_source_id_count':len(seen),'input_counts':input_counts,'fixture':proof,'sector_counts':counts,'parts':enc(parts),'q_cells':enc(result['q_cells']),'physical_external_powers':result['physical_external_powers'],'old_prefix_disjointness_proved':False,'old_errors_removed':False,'appended_to_ledger':False,'local_forcing_computed':False,'physical_pressure_computed':False,'complete_m15_computed':False,'terminal_closure_proved':False}
def main():
    h,s,b=preflight();entry=save(ENTRY,{'component':COMPONENT,'source_sha256':s,'brief_sha256':b,'input_hashes':h,'entered':True});out={'component':COMPONENT,'source_sha256':s,'entry_sha256':entry,'status':'failed'}
    try:out.update(status='returned',output=consume(h),error=None)
    except BaseException as exc:out['error']={'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()}
    digest=save(RESULT,out);print(json.dumps({'component':COMPONENT,'status':out['status'],'result_sha256':digest}));return 0 if out['status']=='returned' else 1
if __name__=='__main__':raise SystemExit(main())
