"""One-entry same-source core first-word ordinary Z0/Z1 at the exit."""
from __future__ import annotations
import gzip,hashlib,importlib.util,json,sys,traceback
from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext
ROOT=Path(__file__).resolve().parents[1];BRIEF=ROOT.parent.parent/'RESUME_BRIEF_2026_10_16.json'
COMPONENT='n1_core_low_word_v3';DIR=ROOT/'tmp'/COMPONENT
SOURCE,ENTRY,RESULT=(DIR/n for n in ('source.json','entered.json','result.json'))
BASE='experiments/root_st073/lei_ren_part1_paper_compliant_'
CORE=BASE+'core_physical_field.json';CORE_CHECK=BASE+'core_physical_field_check.json';TRANSFER=BASE+'core_transfer.json';TRANSFER_CHECK=BASE+'core_transfer_check.json'
PREVIOUS='experiments/root_st073/n1_core_word1_z2_v2_receipt.json.gz';COLLAR='experiments/root_st073/n1_actual_collar_low_word_v1_receipt.json.gz'
MATRIX='src/openai_ns_reconstruction/n1_core_low_matrix_v3.py';GREEN='src/openai_ns_reconstruction/n1_core_low_word_green_v3.py'
MFIX='src/openai_ns_reconstruction/n1_core_low_matrix_fixture_v3.py';GFIX='src/openai_ns_reconstruction/n1_core_low_word_fixture_v3.py';RUNNER='scripts/run_n1_core_low_word_v3_once.py'
FAILURE='experiments/root_st073/n1_core_low_word_v1_failure.json.gz'
FAILURE2='experiments/root_st073/n1_core_low_word_v2_failure.json.gz'
REQUIRED=(CORE,CORE_CHECK,TRANSFER,TRANSFER_CHECK,PREVIOUS,COLLAR,FAILURE,FAILURE2,MATRIX,GREEN,MFIX,GFIX,RUNNER)
FIXED={FAILURE2:'49642c3b9ef562352d25d62e2817bef7648f54be1c411e43c2a1140f09c60888',FAILURE:'e316494fd0b1fee61a66591f25909c2f9e0520f274d283352ac3c1e9596c1b68',CORE:'ad000166f6fce1ea94a7beeccf3324536fbffca073487b873d5460133f380096',CORE_CHECK:'0814fdd1afc27cbd260ec713f834d05804efb7827d2d500cc3551aad8f94fb17',TRANSFER:'4d7ec9b08ee2d36970d22c2a3287b8619bb1fa3deebf5ea264ab0d44b870d497',TRANSFER_CHECK:'91883a7c3d92b1c01ca0d66c5852a5d65814aa364dc7825adb8dd9323c5cb848',PREVIOUS:'c53ff32ad3ce5288279953e84bc7417b71477c82b00d94466d0d2880adbc7ac9',COLLAR:'f4ed142fd3c8ec79b919af92ad12f3e160b06a6c8b13a932bb4fc17df27df234'}
def sha(b):return hashlib.sha256(b).hexdigest()
def require(v,msg):
    if not v:raise ValueError(msg)
def save(path,value):
    b=(json.dumps(value,sort_keys=True,separators=(',',':'))+'\n').encode()
    with path.open('xb') as f:f.write(b)
    return sha(b)
def preflight():
    require(not ENTRY.exists() and not RESULT.exists(),'already entered')
    b=BRIEF.read_bytes();d=json.loads(b);src=SOURCE.read_bytes();ck=d['source_checkpoint']
    require(d['component']==COMPONENT and d['phase']==COMPONENT+'_pending' and d['goal_active'] is True and d['full_goal_complete'] is False and d['one_entry_only'] is True and d['do_not_reenter'] is True,'brief guard')
    require(Path(ck['path']).resolve()==SOURCE.resolve() and ck['sha256']==sha(src),'source checkpoint')
    m=json.loads(src);require(m['component']==COMPONENT and m['one_entry_only'] is True,'manifest');hashes={}
    for row in m['files']:
        p=row['path'];require(p in REQUIRED and p not in hashes,'unexpected input');h=sha((ROOT/p).read_bytes())
        require(h==row['sha256'] and (p not in FIXED or h==FIXED[p]),'changed '+p);hashes[p]=h
    require(set(hashes)==set(REQUIRED),'missing input');return hashes,sha(src),sha(b)
def iv(c,r):
    lo=tuple(r['lower_exact_mpf_tuple']);hi=tuple(r['upper_exact_mpf_tuple']);require(len(lo)==len(hi)==4,'interval endpoints')
    x=c.make_mpf((lo,hi));require(x.a<=x.b,'reversed interval');return x
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
    core=json.loads(raw[CORE]);ck=json.loads(raw[CORE_CHECK]);transfer=json.loads(raw[TRANSFER]);tc=json.loads(raw[TRANSFER_CHECK]);old=json.loads(gzip.decompress(raw[PREVIOUS]));collar=json.loads(gzip.decompress(raw[COLLAR]))
    require(ck['all_passed'] and tc['all_passed'] and old['result']['status']=='returned' and collar['result']['status']=='returned','accepted input')
    require(ck['input_hashes'][Path(CORE).name]==hashes[CORE] and tc['input_hashes'][Path(TRANSFER).name]==hashes[TRANSFER] and core['input_hashes'][Path(TRANSFER).name]==hashes[TRANSFER],'input check binding')
    a=old['result']['output'];b=collar['result']['output'];family=a['source_family'];whole=core['whole_core']
    require(b['source_family']==family and all(core[k]==family[k] for k in ('actual_five_defect_family_sha256','implicit_source_sha256','datum_enclosure_sha256')),'source family')
    require(transfer['implicit_source_sha256']==core['implicit_source_sha256'] and transfer['datum_enclosure_sha256']==core['datum_enclosure_sha256'] and whole['nonlinear_field_enclosed_not_replaced_by_model'] and whole['actual_selected_Cstar_source_retained'],'same analytic source')
    c=MPIntervalContext();c.prec=1800;lam=iv(c,core['actual_Lambda']);delta=iv(c,core['actual_delta']);j=iv(c,transfer['required_j']);z=iv(c,whole['Z'])
    require(z.a==-1 and z.b==1 and lam.a>0 and j.a>0 and iv(c,transfer['Lambda'])._mpi_==lam._mpi_,'same core geometry')
    require(iv(c,a['source_logs']['F0_of_Z'])._mpi_==iv(c,whole['F0_exact_positive_log_enclosure'])._mpi_==iv(c,b['source_logs']['F0_anchor'])._mpi_ and iv(c,a['source_logs']['Lambda'])._mpi_==iv(c,transfer['logLambda'])._mpi_==iv(c,b['source_logs']['Lambda'])._mpi_,'same formal source logs')
    grids=[]
    for key in ('Phi_derivatives','Uz_derivatives','radial_average_Uz_derivatives'):
        grid=whole[key];needed=[(0,k) for k in range(3)]
        if key!='radial_average_Uz_derivatives':needed.extend((1,k) for k in range(2))
        grids.append({f'rho{i}_Z{k}':iv(c,grid[f'rho{i}_Z{k}']) for i,k in needed})
    saved=a['result']['outer_cells'];require(len(saved)==8 and a['result']['no_source_cells_collapsed'] and a['independent_word1_green_fixture']['passed'],'saved core W0 cells')
    outer=[]
    for cell in saved:
        words={order:{col:[{**v,'coefficient':iv(c,v['coefficient'])} for v in terms] for col,terms in cols.items()} for order,cols in cell['W0_ordinary_Z'].items() if order in ('0','1','2')}
        outer.append({'t':iv(c,cell['t']),'rho':iv(c,cell['rho']),'W0_ordinary_Z':words})
    matrix=module(MATRIX,'n1_core_low_matrix_v1');green=module(GREEN,'n1_core_low_word_green_v1');mf=module(MFIX,'n1_core_low_matrix_fixture');gf=module(GFIX,'n1_core_low_word_fixture')
    p1=mf.run(matrix,c);require(p1['passed'],'matrix fixture failed '+str([k for k,v in p1['checks'].items() if not v]))
    p2=gf.run(green,c);require(p2['passed'],'Green fixture failed '+str([k for k,v in p2['checks'].items() if not v]))
    matrices=[{'t':cell['t'],'rho':cell['rho'],**matrix.build(c,cell['rho'],z,delta,lam,j,*grids)} for cell in outer]
    result=green.build(c,outer,matrices)
    return {'source_family':family,'source_logs':b['source_logs'],'domain':'rho=4, real Z in [-1,1]',
      'input_receipts':{PREVIOUS:hashes[PREVIOUS],COLLAR:hashes[COLLAR]},'core_source_hash':hashes[CORE],'core_transfer_hash':hashes[TRANSFER],
      'matrix_fixture':p1,'green_fixture':p2,'matrix_cells':enc(matrices),'result':enc(result),
      'saved_W0_outer_cells_reused':True,'no_ancestor_program_rerun':True,'Z2_not_recomputed':True,'allword_tail_added':False,'m15_computed':False,'terminal_closure_proved':False}
def main():
    h,s,b=preflight();entry=save(ENTRY,{'component':COMPONENT,'source_sha256':s,'brief_sha256':b,'input_hashes':h,'entered':True})
    out={'component':COMPONENT,'source_sha256':s,'entry_sha256':entry,'status':'failed'}
    try:out.update(status='returned',output=consume(h),error=None)
    except BaseException as exc:out['error']={'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()}
    digest=save(RESULT,out);print(json.dumps({'component':COMPONENT,'status':out['status'],'result_sha256':digest}));return 0 if out['status']=='returned' else 1
if __name__=='__main__':raise SystemExit(main())
