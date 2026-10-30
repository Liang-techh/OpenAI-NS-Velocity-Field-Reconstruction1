"""One-entry core W0 and leading pressure m15 components at rho=4."""
from __future__ import annotations
import gzip,hashlib,importlib.util,json,sys,traceback
from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext
ROOT=Path(__file__).resolve().parents[1];BRIEF=ROOT.parent.parent/'RESUME_BRIEF_2026_10_16.json';COMPONENT='n1_core_w0_m15_v1';DIR=ROOT/'tmp'/COMPONENT
SOURCE,ENTRY,RESULT=(DIR/n for n in ('source.json','entered.json','result.json'))
LOW='experiments/root_st073/n1_core_low_word_v3_receipt.json.gz';PREVIOUS='experiments/root_st073/n1_core_word1_z2_v2_receipt.json.gz';COLLAR='experiments/root_st073/n1_first_collar_g4_pressure_v1_receipt.json.gz';SWIRL='experiments/root_st073/n1_first_word_swirl_m15_v1_receipt.json.gz';CORRECTION='experiments/root_st073/n1_first_word_swirl_m15_v1_provenance_correction.json'
CORE='experiments/root_st073/lei_ren_part1_paper_compliant_core_physical_field.json';CHECK='experiments/root_st073/lei_ren_part1_paper_compliant_core_physical_field_check.json';TRANSFER='experiments/root_st073/lei_ren_part1_paper_compliant_core_transfer.json';TC='experiments/root_st073/lei_ren_part1_paper_compliant_core_transfer_check.json';CUT='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_inner_cutoff_functions.py';OP='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_inner_operator.py';TAIL='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_selected_words.json'
LEAF='src/openai_ns_reconstruction/n1_core_w0_m15_v1.py';FIXTURE='src/openai_ns_reconstruction/n1_core_w0_m15_fixture_v1.py';RUNNER='scripts/run_n1_core_w0_m15_v1_once.py'
REQUIRED=(LOW,PREVIOUS,COLLAR,SWIRL,CORRECTION,CORE,CHECK,TRANSFER,TC,CUT,OP,TAIL,LEAF,FIXTURE,RUNNER)
FIXED={LOW:'5a5fe7c57a2efe10e58089ced608e22937c54fcf54c5f224a42abe54e450c6d6',PREVIOUS:'c53ff32ad3ce5288279953e84bc7417b71477c82b00d94466d0d2880adbc7ac9',COLLAR:'c229cf188d306413f6edd0641a729321ad4567c2587ebf4701ada47720853660',SWIRL:'a626e71fa2b0ddf40b79ce992016da2af8ac9955db3f17337029f3c86c2513b4',CORRECTION:'f2133bf5824b2f9636b430a33a2ba1d19fbf391a7f90a81de0c8269d1fe638c4',CORE:'ad000166f6fce1ea94a7beeccf3324536fbffca073487b873d5460133f380096',CHECK:'0814fdd1afc27cbd260ec713f834d05804efb7827d2d500cc3551aad8f94fb17',TRANSFER:'4d7ec9b08ee2d36970d22c2a3287b8619bb1fa3deebf5ea264ab0d44b870d497',TC:'91883a7c3d92b1c01ca0d66c5852a5d65814aa364dc7825adb8dd9323c5cb848',CUT:'b9f8b854b51b9eb0eda338503c8d25e6a1b69307da54b75bf40c1bd212ff1b29',OP:'268217c6b33bc717472742fa64d76a853c6760a956abfad92ac45bcf385df014',TAIL:'0e04c87a011bfec74e0ceb0f822774a03495a08a839aaf4e378d21949736cf1e'}
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
    low=json.loads(gzip.decompress(raw[LOW]));previous=json.loads(gzip.decompress(raw[PREVIOUS]));collar=json.loads(gzip.decompress(raw[COLLAR]));swirl=json.loads(gzip.decompress(raw[SWIRL]));correction=json.loads(raw[CORRECTION]);core=json.loads(raw[CORE]);check=json.loads(raw[CHECK]);transfer=json.loads(raw[TRANSFER]);tc=json.loads(raw[TC]);selected=json.loads(raw[TAIL])
    for r in (low,previous,collar,swirl):require(r['result']['status']=='returned','accepted receipt')
    a=low['result']['output'];p=previous['result']['output'];q=collar['result']['output'];f=swirl['result']['output'];family=a['source_family'];atoms=a['source_logs']
    require(p['source_family']==q['source_family']==f['source_family']==selected['source_family']==family and q['source_logs']==f['source_logs']==atoms and a['input_receipts'][PREVIOUS]==hashes[PREVIOUS],'same source and prior receipts')
    swirl_pins={row['path']:row['sha256'] for row in swirl['source']['files']};require(swirl_pins[CORE]==hashes[CORE] and swirl_pins[TRANSFER]==hashes[TRANSFER],'F0 connection uses the same core and transfer')
    require(check['all_passed'] and tc['all_passed'] and check['input_hashes'][Path(CORE).name]==hashes[CORE] and tc['input_hashes'][Path(TRANSFER).name]==hashes[TRANSFER] and core['input_hashes'][Path(TRANSFER).name]==hashes[TRANSFER] and all(core[k]==family[k] for k in ('actual_five_defect_family_sha256','implicit_source_sha256','datum_enclosure_sha256')),'source/core checks')
    require(correction['receipt_sha256']==hashes[SWIRL] and correction['numerical_result_changed'] is False and a['saved_W0_outer_cells_reused'] and q['core_m15_inlet_added'] is False and a['allword_tail_added'] is False,'separate core and collar inputs')
    op=raw[OP].decode('utf8');cut=raw[CUT].decode('utf8');require("g={'4':eps(forcing['N1p']*(2*x))" in op and "pressure=a.add(a.scale(a.mul(f,F),2),system['N1p'])" in cut and 'a.neg(a.scale(pressure,g.mul(eps2,r)))' in cut and 'self.cutoff_rows=dict(core=[a.fixed(1),a.fixed(0),a.fixed(0)],first=first)' in cut,'physical g4 pressure, velocity density, and unit core cutoff')
    c=MPIntervalContext();c.prec=1800;whole=core['whole_core'];require(iv(c,whole['Z']).a==-1 and iv(c,whole['Z']).b==1 and iv(c,atoms['F0_anchor'])._mpi_==iv(c,whole['F0_exact_positive_log_enclosure'])._mpi_ and iv(c,atoms['Lambda'])._mpi_==iv(c,transfer['logLambda'])._mpi_ and iv(c,transfer['Lambda'])._mpi_==iv(c,core['actual_Lambda'])._mpi_,'exact source logs')
    index={};coreclass=['core_whole_source_interval']
    def term(value,powers,sid,classes=coreclass):
        return {'powers':powers,'coefficient':value,'origin':{'source_ids':[sid],'uncertainty':True,'error_classes':classes}}
    U={};F={};phi=[]
    for order in (0,1):
        for field,out,key,powers in (('U',U,'Uz_derivatives',[0,0,0,0]),('Phi',F,'Phi_derivatives',[0,1,0,0])):
            sid=f'core_leading:{field}:{order}';path=[key,'rho0_Z'+str(order)];index[sid]={'receipt':CORE,'path':['whole_core']+path}
            v=iv(c,whole[key]['rho0_Z'+str(order)]);out[order]=[term(v,powers,sid)]
            if field=='Phi':phi.append(v)
    ell=iv(c,f['F0_connection']['ell_over_lambda']);index['core_F0_connection']={'receipt':SWIRL,'path':['F0_connection','ell_over_lambda'],'path_correction':CORRECTION}
    F[1].append({'powers':[0,1,0,1],'coefficient':phi[0]*ell,'origin':{'source_ids':['core_leading:Phi:0','core_F0_connection'],'uncertainty':True,'error_classes':['core_whole_source_interval','F0_connection_interval']}})
    outer=[]
    for cell,item in enumerate(a['result']['outer_cells']):
        w={}
        for order in (0,1):
            w[order]={}
            for col in ('1','2'):
                w[order][col]=[]
                for i,t in enumerate(item['W0_ordinary_Z'][str(order)][col]):
                    sid=f'core_W0:{cell}:{order}:{col}:{i}';index[sid]={'receipt':LOW,'path':['result','outer_cells',cell,'W0_ordinary_Z',str(order),col,i]}
                    w[order][col].append(term(iv(c,t['coefficient']),[0,t['F0_power'],0,t['Lambda_power']],sid))
        outer.append({'t':iv(c,item['t']),'rho':iv(c,item['rho']),'W0_ordinary_Z':w})
    g4={}
    for order in (0,1):
        g4[order]=[]
        for i,t in enumerate(p['global_low_source'][str(order)]['4']):
            sid=f'core_g4:{order}:{i}';index[sid]={'receipt':PREVIOUS,'path':['global_low_source',str(order),'4',i]}
            g4[order].append(term(iv(c,t['coefficient']),[0,t['F0_power'],0,t['Lambda_power']],sid))
    leaf=module(LEAF,'n1_core_w0_m15_v1');fixture=module(FIXTURE,'n1_core_w0_m15_fixture_v1');proof=fixture.run(leaf,c);require(proof['passed'],'fixture failed '+str([k for k,v in proof['checks'].items() if not v]))
    cert={'g4_equals_two_epsilon_x_N1p':True,'geometry_epsilon_Z_independent':True,'actual_W0_outer_cells':True}
    result=leaf.build(c,outer_cells=outer,leading_U=U,leading_F=F,g4=g4,source_family=family,source_atoms=atoms,certificate=cert)
    seen=set();classes=set()
    for channel in ('axial','swirl'):
        for cells in result[channel].values():
            for terms in cells.values():
                for t in terms:seen.update(t['origin']['source_ids']);classes.update(t['origin']['error_classes'])
    for terms in result['pressure'].values():
        for t in terms:seen.update(t['origin']['source_ids']);classes.update(t['origin']['error_classes'])
    require(seen==set(index) and classes,'lost input provenance')
    return {'source_family':family,'source_logs':atoms,'domain':'rho=4, real Z in [-1,1]','fixture':proof,'input_index':index,'result':enc(result),'output_error_classes':sorted(classes),
      'F0_connection_reused_without_recomputation':True,'source_cells_not_point_values':True,'previous_collar_receipt':COLLAR,'allword_tail_receipt':TAIL,
      'core_W1_radial_history_missing':True,'pressure_direct_density_no_endpoint_added':True,'appended_to_cumulative_m15':False,'allword_tail_added':False,'complete_m15_computed':False,'terminal_closure_proved':False}
def main():
    h,s,b=preflight();entry=save(ENTRY,{'component':COMPONENT,'source_sha256':s,'brief_sha256':b,'input_hashes':h,'entered':True});out={'component':COMPONENT,'source_sha256':s,'entry_sha256':entry,'status':'failed'}
    try:out.update(status='returned',output=consume(h),error=None)
    except BaseException as exc:out['error']={'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()}
    digest=save(RESULT,out);print(json.dumps({'component':COMPONENT,'status':out['status'],'result_sha256':digest}));return 0 if out['status']=='returned' else 1
if __name__=='__main__':raise SystemExit(main())
