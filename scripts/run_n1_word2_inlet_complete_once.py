"""Join previously saved signed inlet source and new core W0_Z3 enclosures."""
from __future__ import annotations
import gzip,hashlib,importlib.util,json,sys,traceback
from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext

ROOT=Path(__file__).resolve().parents[1];BRIEF=ROOT.parent.parent/'RESUME_BRIEF_2026_10_16.json'
COMPONENT='n1_word2_inlet_complete_v1';DIR=ROOT/'tmp'/COMPONENT
SOURCE,ENTRY,RESULT=(DIR/n for n in ('source.json','entered.json','result.json'))
BASE='experiments/root_st073/lei_ren_part1_paper_compliant_'
CORE=BASE+'core_physical_field.json'; CORE_SOURCE=BASE+'core_physical_field.py'
NORM=BASE+'physical_norm_family.json'; AMP=BASE+'anchored_axis_amplitude.py'
OP='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_inner_operator.json.gz'
OP_SOURCE='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_inner_operator.py'
PREV='experiments/root_st073/n1_word2_source_inlet_v1_receipt.json.gz'
WORD='experiments/root_st073/n1_core_g3_word0_v1_receipt.json.gz'
LEAF='src/openai_ns_reconstruction/n1_word2_missing_source_join.py';RUNNER='scripts/run_n1_word2_inlet_complete_once.py'
REQUIRED=(CORE,CORE_SOURCE,NORM,AMP,OP,OP_SOURCE,PREV,WORD,LEAF,RUNNER)
FIXED={CORE:'ad000166f6fce1ea94a7beeccf3324536fbffca073487b873d5460133f380096',
 CORE_SOURCE:'19acdd5f7a506e7582427f004be8e2cd7ae22771eafc76b66ee44d053b9aec54',
 NORM:'e8b97192b8d026a8fa2b22848bb07bae38bc876261d909fb70419b4898efb7e5',
 AMP:'f3ad62f07b7631ccbfe67c787b4989c0d70db2ff49c96325ca3ab2be284ad87c',
 OP:'63b181010f1715d3b20589b06c1c2139ef63bcc73d7b0f3c5b7c43d9c5e03f52',
 OP_SOURCE:'268217c6b33bc717472742fa64d76a853c6760a956abfad92ac45bcf385df014',
 PREV:'fed8a184174e2922ea4097316ff503aba740f2423913d3e30b254eb69c2cc264',
 WORD:'a109cd2b88f3f6740ed4569cb9407d527e98a43287bc56609b41ff08c1052fc8'}


def sha(b):return hashlib.sha256(b).hexdigest()
def require(ok,message):
    if not ok:raise ValueError(message)
def exclusive(path,value):
    b=(json.dumps(value,sort_keys=True,separators=(',',':'))+'\n').encode()
    with path.open('xb') as f:f.write(b)
    return sha(b)


def preflight():
    require(not ENTRY.exists() and not RESULT.exists(),'already entered')
    b=BRIEF.read_bytes();d=json.loads(b)
    require(d.get('component')==COMPONENT and d.get('phase')==COMPONENT+'_pending'
        and d.get('goal_active') is True and d.get('full_goal_complete') is False
        and d.get('one_entry_only') is True and d.get('do_not_reenter') is True,'brief guard')
    ck=d['source_checkpoint'];source=SOURCE.read_bytes()
    require(Path(ck['path']).resolve()==SOURCE.resolve() and sha(source)==ck['sha256'],'checkpoint')
    manifest=json.loads(source);require(manifest['component']==COMPONENT and manifest['one_entry_only'] is True,'manifest')
    hashes={}
    for item in manifest['files']:
        p=item['path'];require(p in REQUIRED and p not in hashes,'unexpected input');h=sha((ROOT/p).read_bytes())
        require(h==item['sha256'] and (p not in FIXED or h==FIXED[p]),'changed '+p);hashes[p]=h
    require(set(hashes)==set(REQUIRED),'missing input');return hashes,sha(source),sha(b)


def iv(c,record):
    lo=tuple(record['lower_exact_mpf_tuple']);hi=tuple(record['upper_exact_mpf_tuple'])
    require(len(lo)==len(hi)==4,'exact endpoints');x=c.make_mpf((lo,hi));require(x.a<=x.b,'reversed interval');return x
def enc(x):
    if hasattr(x,'_mpi_'):return {'lower_exact_mpf_tuple':list(x._mpi_[0]),'upper_exact_mpf_tuple':list(x._mpi_[1]),'directed':True}
    if isinstance(x,dict):return {str(k):enc(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [enc(v) for v in x]
    if isinstance(x,(str,int,float,bool)) or x is None:return x
    raise TypeError('unsupported output '+type(x).__name__)


def consume(hashes):
    raw={}
    for p in REQUIRED:
        b=(ROOT/p).read_bytes();require(sha(b)==hashes[p],'changed after entry '+p);raw[p]=b
    core=json.loads(raw[CORE]);norm=json.loads(raw[NORM]);op=json.loads(gzip.decompress(raw[OP]))
    previous=json.loads(gzip.decompress(raw[PREV]));word=json.loads(gzip.decompress(raw[WORD]))
    require(previous['result']['status']=='returned' and word['result']['status']=='returned','accepted inputs')
    p=previous['result']['output'];w=word['result']['output']
    require(p['source_family']==w['source_family']==op['source_family'],'source family')
    for key in ('actual_five_defect_family_sha256','implicit_source_sha256','datum_enclosure_sha256'):
        require(core[key]==p['source_family'][key],'actual core family '+key)
    norm_name=Path(NORM).name
    require(core['input_hashes'][norm_name]==hashes[NORM] and op['input_hashes'][norm_name]==hashes[NORM],
            'same selected Cstar norm source must be pinned in both archives')
    transfer_name='lei_ren_part1_paper_compliant_core_transfer.json'
    require(core['input_hashes'][transfer_name]==op['input_hashes'][transfer_name], 'same Lambda transfer source')
    cs=raw[CORE_SOURCE].decode('utf-8');amps=raw[AMP].decode('utf-8');ops=raw[OP_SOURCE].decode('utf-8')
    require("self.logC=read(norm,'selected_logCstar')" in cs and 'logF0=-self.core.logC-self.core.Lambda*actual_G' in amps
        and "amplitude=self.source.amplitude.evaluate(z)" in ops and "amplitude['logF0']" in ops,
        'same selected amplitude source definitions')
    c=MPIntervalContext();c.prec=1800
    require(iv(c,core['selected_logCstar'])._mpi_==iv(c,norm['selected_logCstar'])._mpi_,'selected Cstar exact source mismatch')
    z=iv(c,p['selected_Z_interval']);require(z.a>=-1 and z.b<=1,'seam Z outside whole-core enclosure')
    require(w['source_domain']['Z']=='[-1,1]' and w['W0_Z3_green_enclosure']['coordinate']=='core endpoint rho=4, X=2'
        and p['coordinate']=='rho=4, s=0; core W0 and actual first-collar operator', 'core-to-first seam')
    require(w['independent_symbolic_fixture']['passed'] is True and w['W0_Z3_green_enclosure']['cell_count']==8,
        'whole-cell W0_Z3 input')
    require(w['same_actual_analytic_core'] is True and p['same_four_source_anchors'] is True,
        'source binding flags')
    spec=importlib.util.spec_from_file_location('n1_word2_missing_join_leaf',ROOT/LEAF)
    require(spec is not None and spec.loader is not None,'leaf import');module=importlib.util.module_from_spec(spec)
    sys.modules[spec.name]=module;spec.loader.exec_module(module)
    result=module.join(p['result']['rows'],w['W0_Z3_green_enclosure']['fields'],interval=lambda v:iv(c,v))
    return {'source_family':p['source_family'],'selected_Z_interval':p['selected_Z_interval'],
        'binding':{'same_source_family':True,'same_selected_Cstar_norm_hash':hashes[NORM],
          'same_selected_Cstar_exact_record':True,'same_native_amplitude_definition':True,'same_Lambda_transfer_source':True,
          'same_core_exit_and_first_inlet':True,'whole_core_Z_contains_selected_Z':True,
          'independent_positive_factors_not_collapsed':True},
        'W0_source_logs':w['source_logs'],'core_source_input_scope':w['source_domain'],
        'result':enc(result),'interval_overlap_used_as_identity':False,
        'W1_Z2_integrated':False,'m15_computed':False,'global_five_moment_closure':False}


def main():
    h,s,b=preflight();e=exclusive(ENTRY,{'component':COMPONENT,'source_sha256':s,'brief_sha256':b,'input_hashes':h,'entered':True})
    out={'component':COMPONENT,'source_sha256':s,'entry_sha256':e,'status':'failed'}
    try:out['output']=consume(h);out['status']='returned';out['error']=None
    except BaseException as exc:out['error']={'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()}
    digest=exclusive(RESULT,out);print(json.dumps({'component':COMPONENT,'status':out['status'],'result_sha256':digest}))
    return 0 if out['status']=='returned' else 1


if __name__=='__main__':raise SystemExit(main())
