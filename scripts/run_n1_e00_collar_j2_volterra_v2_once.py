"""One-entry source-bound E00 j>=2 Volterra remainder consumer."""
from __future__ import annotations
import gzip,hashlib,importlib.util,json,sys,traceback
from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext
ROOT=Path(__file__).resolve().parents[1];BRIEF=ROOT.parent.parent/'RESUME_BRIEF_2026_10_16.json';COMPONENT='n1_e00_collar_j2_volterra_v2';DIR=ROOT/'tmp'/COMPONENT
SOURCE,ENTRY,RESULT=(DIR/n for n in ('source.json','entered.json','result.json'))
SEAM='experiments/root_st073/n1_e00_collar_fixed_seam_v4_receipt.json.gz';KCHI='experiments/root_st073/n1_core_kchi_first_v1_receipt.json.gz'
LEAF='src/openai_ns_reconstruction/n1_e00_collar_j2_volterra_v2.py';FIXTURE='src/openai_ns_reconstruction/n1_e00_collar_j2_volterra_fixture_v2.py';RUNNER='scripts/run_n1_e00_collar_j2_volterra_v2_once.py'
FAILURE='experiments/root_st073/n1_e00_collar_j2_volterra_v1_failure.json.gz'
REQUIRED=(SEAM,KCHI,FAILURE,LEAF,FIXTURE,RUNNER)
FIXED={FAILURE:'623cc3617e8c342b6fdef5e4a0cfc8855521e6202e2634ff0bb2fb5896174903',SEAM:'6f8cc4c4fb4baedb9d14274c2f8ddc92f4ae536766d3d911a807d6f626539766',KCHI:'b0d654203e8c4acabdd0b61ca7b79c1b7050314b0c003a0c2f8757bde94973a0'}
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
    c=MPIntervalContext();c.prec=1800;leaf=module(LEAF,'n1_e00_collar_j2_volterra_v2');fixture=module(FIXTURE,'n1_e00_collar_j2_volterra_fixture_v2');proof=fixture.run_fixture(leaf,c);require(proof['passed'],'fixture failed '+json.dumps(proof,sort_keys=True,separators=(',',':')))
    s=json.loads(gzip.decompress(raw[SEAM]));k=json.loads(gzip.decompress(raw[KCHI]));require(s['result']['status']==k['result']['status']=='returned','receipt status');o=s['result']['output'];v=o['actual_E00_seam_and_collar'];kc=k['result']['output'];require(o['input_receipt_sha256'][KCHI]==hashes[KCHI] and o['source_family']==kc['source_family'] and o['source_logs']==kc['source_logs'],'axis potential source binding');require(o['seam_identity_source_count']==1632 and o['saved_E00_inlet_preserved_not_reset'] and v['j2_remainder_included'] is False,'saved seam scope')
    index=o['source_index'];branches=v['seam']['k0_transfer']['branches'];inlet={};counts={}
    for field,branch in (('F1','F_from_Fc'),('Fx','Fx_from_Fxc')):
        inlet[field]={}
        for z,expected in (('Z0',88),('Z1',176)):
            rows=branches[branch][z];require(len(rows)==expected,'fixed source count');out=[]
            for t in rows:
                ids=t['source_ids'];require(len(ids)==1 and ids[0] in index and len(index[ids[0]]['parent_source_ids'])>=2,'inlet source index');p=t['powers'];require(len(p)==4 and p[0]==0 and p[2]==0,'formal input powers');require(t['factor_is_identity'] and t['identity_difference_exact_zero'],'fixed seam identity');out.append({'coefficient':iv(c,t['coefficient']),'powers':p,'source_ids':ids,'error_classes':t['error_classes'],'uncertainty':t['uncertainty'],'source_path':index[ids[0]]['path']})
            inlet[field][z]=out;counts[field+':'+z]=len(out)
    chi={}
    for i in (0,1):
        sid='chi:Z'+str(i);require(sid in index and index[sid]['receipt']==KCHI,'chi owner');chi['Z'+str(i)]=[{'coefficient':iv(c,kc['chi']['chi'][i]),'powers':[0,0,0,0],'source_ids':[sid],'error_classes':['whole_Z_axis_potential_interval'],'uncertainty':True}]
    sc=iv(c,v['geometry']['sc']);cert=v['certificate_passthrough'];require(cert['fixed_geometry'] and cert['axis_chi_shared'] and cert['rho_radius_Z_independent'],'saved geometry');result=leaf.build(c,inlet=inlet,chi=chi,hb=c.mpf([0,'0.5']),sc=sc,q_cells=[(0,0),(0,'0.25'),('0.25','0.5'),('0.5',1)],source_family=o['source_family'],source_atoms=v['source_atoms'],certificate=cert)
    require(result['schema']=='n1_e00_collar_j2_volterra_v2' and result['cell_count']==4 and result['j2_remainder_included'] is True and result['appended_to_ledger'] is False,'tail scope');allids={t['source_ids'][0] for rows in inlet.values() for ts in rows.values() for t in ts};seen=set();seam_count=0;term_count=0
    cells=[]
    for i,cell in enumerate(result['q_cells']):
        cb={}
        for name,branch in cell['branches'].items():
            cb[name]={key:val for key,val in branch.items() if key not in ('response',)}
            for z in ('Z0','Z1'):
                for t in branch[z]:
                    require(t['error_only'] is True and t['symmetric_error'] is True and t['tail_starts_at_j']==2 and all(sid in index for sid in t['source_ids']),'tail lineage');require(len(t['input_source_ids'])==1 and t['input_source_ids'][0] in allids,'input lineage');seen.add(t['input_source_ids'][0]);term_count+=1
                    if i==0:require(t['coefficient'].a==t['coefficient'].b==0,'seam tail nonzero');seam_count+=1
        cells.append({'q_interval':cell['q_interval'],'kernel':cell['kernel'],'seam':cell['seam'],'branches':cb})
    require(seen==allids,'missing input sector');return {'source_family':o['source_family'],'source_logs':o['source_logs'],'source_index':{sid:index[sid] for sid in allids|{'chi:Z0','chi:Z1'}},'input_counts':counts,'fixture':proof,'cell_count':len(cells),'tail_term_count':term_count,'zero_seam_term_count':seam_count,'tail_cells':enc(cells),'formal_hb_powers_added':result['formal_hb_powers_added'],'only_symmetric_remainder_bounds':True,'same_source_old_endpoint_identity_proved':False,'old_errors_removed':False,'appended_to_ledger':False,'local_forcing_computed':False,'complete_inlet':False,'m15_computed':False,'terminal_closure_proved':False}
def main():
    h,s,b=preflight();entry=save(ENTRY,{'component':COMPONENT,'source_sha256':s,'brief_sha256':b,'input_hashes':h,'entered':True});out={'component':COMPONENT,'source_sha256':s,'entry_sha256':entry,'status':'failed'}
    try:out.update(status='returned',output=consume(h),error=None)
    except BaseException as exc:out['error']={'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()}
    digest=save(RESULT,out);print(json.dumps({'component':COMPONENT,'status':out['status'],'result_sha256':digest}));return 0 if out['status']=='returned' else 1
if __name__=='__main__':raise SystemExit(main())
