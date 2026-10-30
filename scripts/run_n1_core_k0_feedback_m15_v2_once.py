"""One-entry saved-source core K0-prefix first-feedback profile and m15."""
from __future__ import annotations
import gzip, hashlib, importlib.util, json, sys, traceback
from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext
ROOT=Path(__file__).resolve().parents[1]
BRIEF=ROOT.parent.parent/'RESUME_BRIEF_2026_10_16.json'
COMPONENT='n1_core_k0_feedback_m15_v2'; DIR=ROOT/'tmp'/COMPONENT
SOURCE,ENTRY,RESULT=(DIR/n for n in ('source.json','entered.json','result.json'))
LOW='experiments/root_st073/n1_core_low_word_v3_receipt.json.gz'
CORE='experiments/root_st073/lei_ren_part1_paper_compliant_core_physical_field.json'
SWIRL='experiments/root_st073/n1_first_word_swirl_m15_v1_receipt.json.gz'
CORRECTION='experiments/root_st073/n1_first_word_swirl_m15_v1_provenance_correction.json'
PREVIOUS='experiments/root_st073/n1_core_w0_m15_v1_receipt.json.gz'
SERIES='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_regular_series_function.py'
LEAF='src/openai_ns_reconstruction/n1_core_w1_profile_v2.py'
FIXTURE='src/openai_ns_reconstruction/n1_core_w1_profile_fixture_v2.py'
ADAPTER='src/openai_ns_reconstruction/n1_core_w1_m15_adapter_v2.py'
AFIX='src/openai_ns_reconstruction/n1_core_w1_m15_fixture_v2.py'
RUNNER='scripts/run_n1_core_k0_feedback_m15_v2_once.py'
FAILURE='experiments/root_st073/n1_core_k0_feedback_m15_v1_failure.json.gz'
REQUIRED=(LOW,CORE,SWIRL,CORRECTION,PREVIOUS,SERIES,LEAF,FIXTURE,ADAPTER,AFIX,RUNNER,FAILURE)
FIXED={FAILURE:'974d9a88e5671113e42eb6192807fe9d83043f1279fe5751e2821041b44f34c5',LOW:'5a5fe7c57a2efe10e58089ced608e22937c54fcf54c5f224a42abe54e450c6d6',CORE:'ad000166f6fce1ea94a7beeccf3324536fbffca073487b873d5460133f380096',SWIRL:'a626e71fa2b0ddf40b79ce992016da2af8ac9955db3f17337029f3c86c2513b4',CORRECTION:'f2133bf5824b2f9636b430a33a2ba1d19fbf391a7f90a81de0c8269d1fe638c4',PREVIOUS:'157f5677742b71a469c67e0efd21de0d94d230b4d3c2968d82acd2972413bda7'}
def sha(b): return hashlib.sha256(b).hexdigest()
def require(v,msg):
    if not v: raise ValueError(msg)
def save(path,value):
    b=(json.dumps(value,sort_keys=True,separators=(',',':'))+'\n').encode()
    with path.open('xb') as f: f.write(b)
    return sha(b)
def preflight():
    require(not ENTRY.exists() and not RESULT.exists(),'already entered')
    brief=BRIEF.read_bytes(); d=json.loads(brief); src=SOURCE.read_bytes(); ck=d['source_checkpoint']
    require(d['component']==COMPONENT and d['phase']==COMPONENT+'_pending' and d['goal_active'] is True and d['full_goal_complete'] is False and d['one_entry_only'] is True and d['do_not_reenter'] is True,'brief guard')
    require(Path(ck['path']).resolve()==SOURCE.resolve() and ck['sha256']==sha(src),'source checkpoint')
    m=json.loads(src); require(m['component']==COMPONENT and m['one_entry_only'] is True,'manifest'); hashes={}
    for row in m['files']:
        p=row['path']; require(p in REQUIRED and p not in hashes,'unexpected input'); h=sha((ROOT/p).read_bytes())
        require(h==row['sha256'] and (p not in FIXED or h==FIXED[p]),'changed '+p); hashes[p]=h
    require(set(hashes)==set(REQUIRED),'missing input'); return hashes,sha(src),sha(brief)
def iv(c,r):
    lo=tuple(r['lower_exact_mpf_tuple']); hi=tuple(r['upper_exact_mpf_tuple']); require(len(lo)==len(hi)==4,'interval endpoints')
    x=c.make_mpf((lo,hi)); require(x.a<=x.b,'reversed interval'); return x
def enc(x):
    if hasattr(x,'_mpi_'): return {'lower_exact_mpf_tuple':list(x._mpi_[0]),'upper_exact_mpf_tuple':list(x._mpi_[1]),'directed':True}
    if isinstance(x,dict): return {str(k):enc(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)): return [enc(v) for v in x]
    if isinstance(x,(str,int,float,bool)) or x is None: return x
    raise TypeError(type(x).__name__)
def module(path,name):
    spec=importlib.util.spec_from_file_location(name,ROOT/path); m=importlib.util.module_from_spec(spec); sys.modules[name]=m; spec.loader.exec_module(m); return m
def consume(hashes):
    raw={}
    for p in REQUIRED:
        b=(ROOT/p).read_bytes(); require(sha(b)==hashes[p],'changed after entry '+p); raw[p]=b
    low=json.loads(gzip.decompress(raw[LOW])); core=json.loads(raw[CORE]); swirl=json.loads(gzip.decompress(raw[SWIRL])); correction=json.loads(raw[CORRECTION]); previous=json.loads(gzip.decompress(raw[PREVIOUS]))
    for receipt in (low,swirl,previous): require(receipt['result']['status']=='returned','input receipt status')
    a=low['result']['output']; f=swirl['result']['output']; prev=previous['result']['output']; family=a['source_family']; atoms=a['source_logs']
    require(f['source_family']==prev['source_family']==family and f['source_logs']==prev['source_logs']==atoms,'source identity')
    pins={v['path']:v['sha256'] for v in previous['source']['files']}; require(pins[LOW]==hashes[LOW] and pins[CORE]==hashes[CORE] and pins[SWIRL]==hashes[SWIRL] and pins[CORRECTION]==hashes[CORRECTION],'same previously bound inputs')
    require(correction['receipt_sha256']==hashes[SWIRL] and a['saved_W0_outer_cells_reused'] and a['no_ancestor_program_rerun'] and prev['appended_to_cumulative_m15'] is False,'previous finite source scope')
    series=raw[SERIES].decode('utf8'); require('Kchi=sum_j>=0(-2chi)^j' in series and 'exact_infinite_swirl_kernel' in series,'original Kchi series missing')
    c=MPIntervalContext(); c.prec=1800
    def convert(value):
        if isinstance(value,dict):
            if 'lower_exact_mpf_tuple' in value and 'upper_exact_mpf_tuple' in value: return iv(c,value)
            return {k:convert(v) for k,v in value.items()}
        if isinstance(value,list): return [convert(v) for v in value]
        return value
    outer=convert(a['result']['outer_cells']); matrices=convert(a['matrix_cells']); require(len(outer)==len(matrices)==8,'saved cells')
    leaf=module(LEAF,'n1_core_w1_profile_v2'); fixture=module(FIXTURE,'n1_core_w1_profile_fixture_v2'); adapter=module(ADAPTER,'n1_core_w1_m15_adapter_v2'); afix=module(AFIX,'n1_core_w1_m15_fixture_v2')
    proof=fixture.run_fixture(leaf.build,c); require(proof['passed'],'profile fixture failed '+str([k for k,v in proof['checks'].items() if not v]))
    p2=afix.run(adapter,c); require(p2['passed'],'m15 fixture failed '+str([k for k,v in p2['checks'].items() if not v]))
    profile=leaf.build(c,outer,matrices)
    comparisons={}
    def key(term,old):
        return (term['outer_cell'] if old else term['source_cell'],term['origin'],term['matrix_key'],str(term['matrix_order']),term['matrix_term_index'],str(term['W0_order']),term['W0_term_index'])
    for order in (0,1):
        for field in ('F1','Uz1'):
            old=convert(a['result']['W1_Z'+str(order)+'_core_exit'][field]); new=profile['W1_Z'+str(order)+'_core_exit'][field]
            oldmap={key(t,True):t for t in old}; newmap={key(t,False):t for t in new}
            keys_match=len(oldmap)==len(old) and len(newmap)==len(new) and set(oldmap)==set(newmap)
            contained=0
            if keys_match:
                for k,t in newmap.items():
                    o=oldmap[k]; x=t['coefficient']; y=o['coefficient']
                    if x.a>=y.a and x.b<=y.b and t['F0_power']==o['F0_power'] and t['Lambda_power']==o['Lambda_power']: contained+=1
            comparisons[f'Z{order}_{field}']={'old_count':len(old),'new_count':len(new),'term_keys_match':keys_match,'new_inside_old_count':contained,'point_or_function_equality_claimed':False}
    whole=core['whole_core']; require(iv(c,whole['Z']).a==-1 and iv(c,whole['Z']).b==1 and iv(c,atoms['F0_anchor'])._mpi_==iv(c,whole['F0_exact_positive_log_enclosure'])._mpi_,'whole core source')
    index={}; classes=['core_whole_source_interval']
    def leadterm(value,powers,sid,extra=None):
        return {'powers':powers,'coefficient':value,'origin':{'source_ids':[sid],'uncertainty':True,'error_classes':extra or classes}}
    U={}; F={}; phi=[]
    for order in (0,1):
        for field,out,key,powers in (('U',U,'Uz_derivatives',[0,0,0,0]),('Phi',F,'Phi_derivatives',[0,1,0,0])):
            sid=f'core_leading:{field}:{order}'; index[sid]={'receipt':CORE,'path':['whole_core',key,'rho0_Z'+str(order)]}
            v=iv(c,whole[key]['rho0_Z'+str(order)]); out[order]=[leadterm(v,powers,sid)]
            if field=='Phi': phi.append(v)
    ell=iv(c,f['F0_connection']['ell_over_lambda']); index['core_F0_connection']={'receipt':SWIRL,'path':['F0_connection','ell_over_lambda'],'path_correction':CORRECTION}
    F[1].append({'powers':[0,1,0,1],'coefficient':phi[0]*ell,'origin':{'source_ids':['core_leading:Phi:0','core_F0_connection'],'uncertainty':True,'error_classes':['core_whole_source_interval','F0_connection_interval']}})
    normalized=[]
    for cell in profile['target_cells']:
        w={}
        for order in (0,1):
            w[order]={}
            for field,col in (('F1','1'),('Uz1','2')):
                terms=[]
                for term in cell['W1_by_order']['Z'+str(order)][field]:
                    j=term['source_cell']; mn=term['matrix_name']; mo=str(term['matrix_order']); mk=term['matrix_key']; mi=term['matrix_term_index']; wo=str(term['W0_order']); wc=str(term['W0_column']); wi=term['W0_term_index']
                    mid=f'matrix:{j}:{mn}:{mo}:{mk}:{mi}'; wid=f'W0:{j}:{wo}:{wc}:{wi}'
                    index[mid]={'receipt':LOW,'path':['matrix_cells',j,mn,mo,mk,mi]}; index[wid]={'receipt':LOW,'path':['result','outer_cells',j,'W0_ordinary_Z',wo,wc,wi]}
                    terms.append({'powers':[0,term['F0_power'],0,term['Lambda_power']],'coefficient':term['coefficient'],'origin':{'source_ids':[mid,wid],'uncertainty':True,'error_classes':classes,'source_cell':j,'target_cell':cell['index']}})
                w[order][col]=terms
        normalized.append({'t':cell['t'],'rho':cell['rho'],'W1_ordinary_Z':w})
    contribution=adapter.build(c,profile=normalized,leading_U=U,leading_F=F,source_family=family,source_atoms=atoms,certificate={'geometry_epsilon_Z_independent':True})
    return {'source_family':family,'source_logs':atoms,'domain':'rho target in [0,4], real Z in [-1,1]','profile_fixture':proof,'m15_fixture':p2,'exit_comparison':comparisons,'profile':enc(profile),'m15':enc(contribution),'input_index':index,
      'source_kernel_scope':'K0-only nested finite contribution; not a full Kchi W1 enclosure','Kchi_minus_K0_included':False,'allword_tail_added':False,'g4_pressure_added':False,'appended_to_cumulative_m15':False,'complete_m15_computed':False,'terminal_closure_proved':False}
def main():
    h,s,b=preflight(); entry=save(ENTRY,{'component':COMPONENT,'source_sha256':s,'brief_sha256':b,'input_hashes':h,'entered':True}); out={'component':COMPONENT,'source_sha256':s,'entry_sha256':entry,'status':'failed'}
    try: out.update(status='returned',output=consume(h),error=None)
    except BaseException as exc: out['error']={'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()}
    digest=save(RESULT,out); print(json.dumps({'component':COMPONENT,'status':out['status'],'result_sha256':digest})); return 0 if out['status']=='returned' else 1
if __name__=='__main__': raise SystemExit(main())
