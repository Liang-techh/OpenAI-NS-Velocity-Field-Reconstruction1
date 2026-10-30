"""One-entry same-source first Kchi correction and separate kernel tail."""
from __future__ import annotations
import gzip,hashlib,importlib.util,json,sys,traceback,math
from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext
ROOT=Path(__file__).resolve().parents[1]; BRIEF=ROOT.parent.parent/'RESUME_BRIEF_2026_10_16.json'; COMPONENT='n1_core_kchi_first_v1'; DIR=ROOT/'tmp'/COMPONENT
SOURCE,ENTRY,RESULT=(DIR/n for n in ('source.json','entered.json','result.json'))
LOW='experiments/root_st073/n1_core_word1_z2_v2_receipt.json.gz';LOW3='experiments/root_st073/n1_core_low_word_v3_receipt.json.gz'; PREV='experiments/root_st073/n1_core_k0_feedback_m15_v2_receipt.json.gz'
CORE='experiments/root_st073/lei_ren_part1_paper_compliant_core_physical_field.json'; TRANSFER='experiments/root_st073/lei_ren_part1_paper_compliant_core_transfer.json'; CSRC='experiments/root_st073/lei_ren_part1_paper_compliant_core_physical_field.py'; SERIES='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_regular_series_function.py'; WEIGHT='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_weighted_propagator.py'
CHI='src/openai_ns_reconstruction/n1_axis_chi_jet_v1.py';CF='src/openai_ns_reconstruction/n1_axis_chi_jet_fixture_v1.py';LEAF='src/openai_ns_reconstruction/n1_kchi_first_correction_v1.py';FIX='src/openai_ns_reconstruction/n1_kchi_first_correction_fixture_v1.py';TAIL='src/openai_ns_reconstruction/n1_kchi_remainder_v1.py';TF='src/openai_ns_reconstruction/n1_kchi_remainder_fixture_v1.py';RUNNER='scripts/run_n1_core_kchi_first_v1_once.py'
REQUIRED=(LOW,LOW3,PREV,CORE,TRANSFER,CSRC,SERIES,WEIGHT,CHI,CF,LEAF,FIX,TAIL,TF,RUNNER)
FIXED={LOW:'c53ff32ad3ce5288279953e84bc7417b71477c82b00d94466d0d2880adbc7ac9',LOW3:'5a5fe7c57a2efe10e58089ced608e22937c54fcf54c5f224a42abe54e450c6d6',PREV:'ca947df9fdf7765d9211f428375d81956af783d89074f9e782ba56aeff91dad0',CORE:'ad000166f6fce1ea94a7beeccf3324536fbffca073487b873d5460133f380096',TRANSFER:'4d7ec9b08ee2d36970d22c2a3287b8619bb1fa3deebf5ea264ab0d44b870d497',CSRC:'19acdd5f7a506e7582427f004be8e2cd7ae22771eafc76b66ee44d053b9aec54',SERIES:'b64b7cae8529cfc5caa158fe06ac5660ec1d2217fc4d845d02002dac3a195555',WEIGHT:'5219a968018601d8c445815058f7c1c832398dc417fcb51c689688d03011ba62'}
def sha(b):return hashlib.sha256(b).hexdigest()
def require(ok,msg):
    if not ok:raise ValueError(msg)
def save(path,value):
    b=(json.dumps(value,sort_keys=True,separators=(',',':'))+'\n').encode()
    with path.open('xb') as f:f.write(b)
    return sha(b)
def preflight():
    require(not ENTRY.exists() and not RESULT.exists(),'already entered');b=BRIEF.read_bytes();d=json.loads(b);src=SOURCE.read_bytes();ck=d['source_checkpoint']
    require(d['component']==COMPONENT and d['phase']==COMPONENT+'_pending' and d['goal_active'] is True and d['full_goal_complete'] is False and d['one_entry_only'] is True and d['do_not_reenter'] is True,'brief guard')
    require(Path(ck['path']).resolve()==SOURCE.resolve() and ck['sha256']==sha(src),'checkpoint');manifest=json.loads(src);require(manifest['component']==COMPONENT and manifest['one_entry_only'] is True,'manifest');hashes={}
    for row in manifest['files']:
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
    low=json.loads(gzip.decompress(raw[LOW])); low3=json.loads(gzip.decompress(raw[LOW3])); prev=json.loads(gzip.decompress(raw[PREV])); core=json.loads(raw[CORE]); transfer=json.loads(raw[TRANSFER])
    for receipt in (low,low3,prev):require(receipt['result']['status']=='returned','input receipt status')
    a=low['result']['output']; b=low3['result']['output']; p=prev['result']['output']; family=a['source_family']; atoms=b['source_logs']
    require(b['source_family']==p['source_family']==family and p['source_logs']==atoms and b['input_receipts'][LOW]==hashes[LOW],'same source family and prior receipt')
    pins={v['path']:v['sha256'] for v in prev['source']['files']};require(pins[LOW3]==hashes[LOW3] and pins[CORE]==hashes[CORE],'previous source binding')
    require(core['input_hashes'][Path(TRANSFER).name]==hashes[TRANSFER] and all(core[k]==family[k] for k in ('actual_five_defect_family_sha256','implicit_source_sha256','datum_enclosure_sha256')) and transfer['implicit_source_sha256']==core['implicit_source_sha256'] and transfer['datum_enclosure_sha256']==core['datum_enclosure_sha256'],'core source owner')
    src=raw[CSRC].decode('utf8');series=raw[SERIES].decode('utf8');weight=raw[WEIGHT].decode('utf8')
    require('chi=h2/denominator' in src and 'u=4*z+self.j; d=one-z*z; L=one-(z*z)*self.delta; H=z*((1-self.delta)/2)+d*u' in src and 'Kchi=sum_j>=0(-2chi)^j' in series and 'K_j=T^j K0 is positive' in weight,'original chi and kernel definitions')
    c=MPIntervalContext();c.prec=1800;whole=core['whole_core'];z=iv(c,whole['Z']);delta=iv(c,core['actual_delta']);j=iv(c,transfer['required_j'])
    require(z.a==-1 and z.b==1 and j.a>0 and iv(c,core['actual_Lambda'])._mpi_==iv(c,transfer['Lambda'])._mpi_ and iv(c,atoms['F0_anchor'])._mpi_==iv(c,whole['F0_exact_positive_log_enclosure'])._mpi_,'domain and scale')
    chi_leaf=module(CHI,'n1_axis_chi_jet_v1');cf=module(CF,'n1_axis_chi_jet_fixture_v1');leaf=module(LEAF,'n1_kchi_first_correction_v1');fix=module(FIX,'n1_kchi_first_correction_fixture_v1');tail=module(TAIL,'n1_kchi_remainder_v1');tf=module(TF,'n1_kchi_remainder_fixture_v1')
    proofs={'chi':cf.run(chi_leaf,c),'kernel':fix.run_fixture(leaf.build,c),'tail':tf.run(tail,c)}
    for name,proof in proofs.items():require(proof['passed'],name+' fixture failed '+str([k for k,v in proof['checks'].items() if not v]))
    jets=chi_leaf.build(c,z=z,delta=delta,j=j);chi=jets['chi'];require(chi[0].a>=0 and chi[0].b<=1,'real chi range')
    index={};rows={}
    for order in (0,1,2):
        rows[order]=[]
        for i,t in enumerate(a['global_low_source'][str(order)]['5']):
            sid=f'g5:Z{order}:{i}';index[sid]={'receipt':LOW,'path':['global_low_source',str(order),'5',i]}
            rows[order].append({'powers':[0,t['F0_power'],0,t['Lambda_power']],'coefficient':iv(c,t['coefficient']),'origin':{'source_ids':[sid],'error_classes':['whole_core_source_interval'],'uncertainty':True}})
    chi_rows={}
    for order in (0,1,2):
        sid=f'chi:Z{order}';index[sid]={'inputs':[CORE,TRANSFER,CSRC],'definition':'chi=H^2/(H^2+(j/500)^2)','ordinary_Z_order':order}
        chi_rows[order]=[{'powers':[0,0,0,0],'coefficient':chi[order],'origin':{'source_ids':[sid],'error_classes':['whole_Z_chi_interval'],'uncertainty':True}}]
    result=leaf.build(c,source_rows=[rows for _ in range(8)],chi=chi_rows,source_family=family,source_atoms=atoms)
    kernel_tail=tail.build(c,a=c.mpf(2),chi_z1=chi[1],chi_z2=chi[2]); error={}
    for field in ('F1','d_x_F1'):
        error[field]={}
        for order in (0,1,2):
            terms=[]
            for source_order in range(order+1):
                for item in rows[source_order]:
                    factor=c.mpf(2*math.comb(order,source_order)); upper=factor*kernel_tail[field][order-source_order]*abs(item['coefficient']); sym=c.make_mpf(((-upper)._mpi_[0],upper._mpi_[1]))
                    terms.append({'powers':item['powers'],'coefficient':sym,'source_ids':item['origin']['source_ids']+['chi:Z0','chi:Z1','chi:Z2'],'error_classes':['Kchi_j_ge_2_kernel_tail','whole_core_source_interval'],'uncertainty':True,'source_order':source_order,'kernel_Z_order':order-source_order,'physical_source_width_bound':2})
            error[field][order]=terms
    return {'source_family':family,'source_logs':atoms,'domain':'x in [0,2], rho in [0,4], real Z in [-1,1]','fixtures':proofs,'chi':enc(jets),'first_correction':enc(result),'kernel_tail_bounds':enc(kernel_tail),'W0_separate_symmetric_kernel_tail':enc(error),'input_index':index,
      'source_cells_reuse_one_whole_core_enclosure':True,'source_boxes_not_point_values':True,'first_correction_signed':True,'Kchi_tail_is_uncertainty_not_signed_field':True,'tail_includes_feedback':False,'feedback_cross_terms_added':False,'appended_to_cumulative_m15':False,'complete_W0_computed':False,'complete_m15_computed':False,'terminal_closure_proved':False}
def main():
    h,s,b=preflight();entry=save(ENTRY,{'component':COMPONENT,'source_sha256':s,'brief_sha256':b,'input_hashes':h,'entered':True});out={'component':COMPONENT,'source_sha256':s,'entry_sha256':entry,'status':'failed'}
    try:out.update(status='returned',output=consume(h),error=None)
    except BaseException as exc:out['error']={'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()}
    digest=save(RESULT,out);print(json.dumps({'component':COMPONENT,'status':out['status'],'result_sha256':digest}));return 0 if out['status']=='returned' else 1
if __name__=='__main__':raise SystemExit(main())
