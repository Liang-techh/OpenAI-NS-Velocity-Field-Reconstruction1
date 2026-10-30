"""One-entry signed core W1_Z2 enclosure using saved actual source intervals."""
from __future__ import annotations
import gzip,hashlib,importlib.util,json,sys,traceback
from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext

ROOT=Path(__file__).resolve().parents[1];BRIEF=ROOT.parent.parent/'RESUME_BRIEF_2026_10_16.json'
COMPONENT='n1_core_word1_z2_v2';DIR=ROOT/'tmp'/COMPONENT
SOURCE,ENTRY,RESULT=(DIR/n for n in ('source.json','entered.json','result.json'))
BASE='experiments/root_st073/lei_ren_part1_paper_compliant_'
CORE=BASE+'core_physical_field.json';CORE_CHECK=BASE+'core_physical_field_check.json';CORE_SOURCE=BASE+'core_physical_field.py'
TRANSFER=BASE+'core_transfer.json';TRANSFER_CHECK=BASE+'core_transfer_check.json';TRANSFER_SOURCE=BASE+'core_transfer.py'
OP_SOURCE='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_inner_operator.py'
INNER='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_inner_volterra_functions.py'
PRIOR='experiments/root_st073/n1_core_g3_word0_v1_receipt.json.gz'
FAILURE='experiments/root_st073/n1_core_word1_z2_v1_failure.json.gz'
LOW='src/openai_ns_reconstruction/n1_core_low_source_matrix_v2.py'
GREEN='src/openai_ns_reconstruction/n1_core_word1_z2_green_v2.py'
LOW_FIXTURE='src/openai_ns_reconstruction/n1_core_low_independent_fixture_v2.py'
GREEN_FIXTURE='src/openai_ns_reconstruction/n1_word1_green_fixture_v2.py'
RUNNER='scripts/run_n1_core_word1_z2_v2_once.py'
REQUIRED=(CORE,CORE_CHECK,CORE_SOURCE,TRANSFER,TRANSFER_CHECK,TRANSFER_SOURCE,OP_SOURCE,INNER,PRIOR,FAILURE,LOW,GREEN,LOW_FIXTURE,GREEN_FIXTURE,RUNNER)
FIXED={CORE:'ad000166f6fce1ea94a7beeccf3324536fbffca073487b873d5460133f380096',
 CORE_CHECK:'0814fdd1afc27cbd260ec713f834d05804efb7827d2d500cc3551aad8f94fb17',
 CORE_SOURCE:'19acdd5f7a506e7582427f004be8e2cd7ae22771eafc76b66ee44d053b9aec54',
 TRANSFER:'4d7ec9b08ee2d36970d22c2a3287b8619bb1fa3deebf5ea264ab0d44b870d497',
 TRANSFER_CHECK:'91883a7c3d92b1c01ca0d66c5852a5d65814aa364dc7825adb8dd9323c5cb848',
 TRANSFER_SOURCE:'2ad468af4815afe650eb3ebedb29aa392062c4b8a9c7b29e663540a1a2d199d6',
 OP_SOURCE:'268217c6b33bc717472742fa64d76a853c6760a956abfad92ac45bcf385df014',
 INNER:'ff82f4edb03922a551cba3e2df20217a55dc80ed42bd39ccb8b77b83f63b0208',
 PRIOR:'a109cd2b88f3f6740ed4569cb9407d527e98a43287bc56609b41ff08c1052fc8',
 FAILURE:'8becc35f9c3d86fd3ee4b187ff9eef8cec41b9edb415f1427260daeccb0793bd'}


def sha(b):return hashlib.sha256(b).hexdigest()
def require(ok,msg):
    if not ok:raise ValueError(msg)
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
    ck=d['source_checkpoint'];src=SOURCE.read_bytes()
    require(Path(ck['path']).resolve()==SOURCE.resolve() and sha(src)==ck['sha256'],'checkpoint')
    manifest=json.loads(src);require(manifest['component']==COMPONENT and manifest['one_entry_only'] is True,'manifest')
    hashes={}
    for row in manifest['files']:
        p=row['path'];require(p in REQUIRED and p not in hashes,'unexpected input');h=sha((ROOT/p).read_bytes())
        require(h==row['sha256'] and (p not in FIXED or h==FIXED[p]),'changed '+p);hashes[p]=h
    require(set(hashes)==set(REQUIRED),'missing input');return hashes,sha(src),sha(b)


def iv(c,record):
    lo=tuple(record['lower_exact_mpf_tuple']);hi=tuple(record['upper_exact_mpf_tuple'])
    require(len(lo)==len(hi)==4,'exact endpoints');x=c.make_mpf((lo,hi));require(x.a<=x.b,'reversed interval');return x
def enc(x):
    if hasattr(x,'_mpi_'):return {'lower_exact_mpf_tuple':list(x._mpi_[0]),'upper_exact_mpf_tuple':list(x._mpi_[1]),'directed':True}
    if isinstance(x,dict):return {str(k):enc(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [enc(v) for v in x]
    if isinstance(x,(str,int,float,bool)) or x is None:return x
    raise TypeError('unsupported output '+type(x).__name__)
def module(path,name):
    spec=importlib.util.spec_from_file_location(name,ROOT/path);require(spec is not None and spec.loader is not None,'leaf import')
    obj=importlib.util.module_from_spec(spec);sys.modules[spec.name]=obj;spec.loader.exec_module(obj);return obj


def consume(hashes):
    raw={}
    for p in REQUIRED:
        b=(ROOT/p).read_bytes();require(sha(b)==hashes[p],'changed after entry '+p);raw[p]=b
    core=json.loads(raw[CORE]);ck=json.loads(raw[CORE_CHECK]);transfer=json.loads(raw[TRANSFER]);tc=json.loads(raw[TRANSFER_CHECK]);prior=json.loads(gzip.decompress(raw[PRIOR]))
    require(ck['all_passed'] is True and tc['all_passed'] is True and prior['result']['status']=='returned','accepted inputs')
    require(ck['input_hashes'][Path(CORE).name]==hashes[CORE] and tc['input_hashes'][Path(TRANSFER).name]==hashes[TRANSFER]
        and core['input_hashes'][Path(TRANSFER).name]==hashes[TRANSFER],'source/check binding')
    prior_hashes={item['path']:item['sha256'] for item in prior['source']['files']}
    require(prior_hashes[CORE]==hashes[CORE] and prior_hashes[TRANSFER]==hashes[TRANSFER],'same saved whole-core source')
    previous=prior['result']['output'];family=previous['source_family']
    for key in ('actual_five_defect_family_sha256','implicit_source_sha256','datum_enclosure_sha256'):
        require(core[key]==family[key],'source family '+key)
    require(transfer['implicit_source_sha256']==core['implicit_source_sha256']
        and transfer['datum_enclosure_sha256']==core['datum_enclosure_sha256'],'transfer family')
    whole=core['whole_core'];require(core['whole_Z_analytic_core_profile_enclosures_through_order5_available'] is True
        and whole['nonlinear_field_enclosed_not_replaced_by_model'] is True
        and whole['actual_selected_Cstar_source_retained'] is True,'actual analytic core')
    cs=raw[CORE_SOURCE].decode('utf-8');ts=raw[TRANSFER_SOURCE].decode('utf-8');ops=raw[OP_SOURCE].decode('utf-8');ins=raw[INNER].decode('utf-8')
    require(all(s in cs for s in ("self.j=read(major,'required_j')",'self.sigma=self.j/500',
        'gradient=L*H/denominator',"whole=self.profiles([0,'4.1'],[-1,1])"))
        and 'lam=ctx.exp(logLambda);eps=ctx.exp(-logLambda)' in ts
        and "B0['5,5']=eps(B*x);B0['6,6']=eps(B*x)" in ops
        and 'endpoint_independent_of_Z=True' in ins and 'complete_history_not_reset_at_core_seam=True' in ins,
        'native source/operator/Green definitions')
    c=MPIntervalContext();c.prec=1800
    whole_rho=iv(c,whole['rho']);z=iv(c,whole['Z']);delta=iv(c,core['actual_delta']);lam=iv(c,core['actual_Lambda']);j=iv(c,transfer['required_j'])
    require(whole_rho.a==0 and whole_rho.b>=4 and z.a==-1 and z.b==1 and lam.a>0 and j.a>0,'core domain')
    require(iv(c,transfer['Lambda'])._mpi_==lam._mpi_,'same exact Lambda')
    require(iv(c,previous['source_logs']['F0_of_Z'])._mpi_==iv(c,whole['F0_exact_positive_log_enclosure'])._mpi_
        and iv(c,previous['source_logs']['Lambda'])._mpi_==iv(c,transfer['logLambda'])._mpi_,
        'same exact current-Z F0 and Lambda source-log enclosures')
    grids=[]
    for name in ('Phi_derivatives','Uz_derivatives','radial_average_Uz_derivatives'):
        grid=whole[name];grids.append({f'rho{i}_Z{k}':iv(c,grid[f'rho{i}_Z{k}']) for i in range(3) for k in range(6-i)})
    low=module(LOW,'n1_core_low_source_matrix_leaf');green=module(GREEN,'n1_core_word1_z2_green_leaf')
    low_fixture=module(LOW_FIXTURE,'n1_core_low_independent_fixture');green_fixture=module(GREEN_FIXTURE,'n1_word1_green_fixture')
    f1=low_fixture.run_fixture(low.build,c);require(f1['passed'] is True,'independent source/matrix fixture')
    f2=green_fixture.run_fixture(green.build,c);require(f2['passed'] is True,'independent nested Green fixture')
    old=[]
    require(len(previous['cells'])==8 and previous['source_domain']['Z']=='[-1,1]'
        and previous['same_actual_analytic_core'] is True
        and previous['independent_symbolic_fixture']['passed'] is True
        and prior['frozen_byte_sha256']['result']=='21dff2d36e4d0f1c43189be6ab2fe43eeade60245bbd1a25ece502c2d8b396c2', 'prior g3 cells')
    for i,cell in enumerate(previous['cells']):
        left=c.mpf(i)/8;right=c.mpf(i+1)/8;t=c.make_mpf((left._mpi_[0],right._mpi_[1]))
        require(cell['t']['directed'] is True and cell['rho']['directed'] is True
            and iv(c,cell['t'])._mpi_==t._mpi_ and iv(c,cell['rho'])._mpi_==(4*t*t)._mpi_,'prior g3 cell coordinate')
        rows={}
        for row in ('4','5','6'):
            expected={'4':[(0,0),(0,-1)],'5':[(1,-1)],'6':[(0,-1),(0,-2)]}[row]
            require([(v['F0_power'],v['Lambda_power']) for v in cell['g'][row]]==expected, 'prior g3 source factors')
            require(all(v['coefficient']['directed'] is True for v in cell['g'][row]),'prior g3 directed endpoints')
            rows[row]=[{'F0_power':v['F0_power'],'Lambda_power':v['Lambda_power'],'coefficient':iv(c,v['coefficient'])} for v in cell['g'][row]]
        old.append({'t':t,'rho':4*t*t,'g':rows})
    rho_all=c.mpf([0,4]);global_low=low.build(c,rho_all,z,delta,lam,j,*grids,source_only=True)['g']
    matrix=[]
    for i in range(8):
        left=c.mpf(i)/8;right=c.mpf(i+1)/8;t=c.make_mpf((left._mpi_[0],right._mpi_[1]));rho=4*t*t
        require(old[i]['t']._mpi_==t._mpi_ and old[i]['rho']._mpi_==rho._mpi_,'source/outer cell alignment')
        rows=low.build(c,rho,z,delta,lam,j,*grids,matrix_only=True)
        matrix.append({'t':t,'rho':rho,'C0':rows['C0'],'B1':rows['B1']})
    result=green.build(c,global_low,old,matrix)
    return {'source_family':family,'source_domain':{'archived_derivative_rho':'[0,4.1]',
          'integration_rho':'[0,4]','endpoint_rho':'4','Z':'[-1,1]',
          'cell_inputs':'same actual whole-core enclosures reused; not independently sampled cell values'},
        'source_logs':previous['source_logs'],'prior_g3_result_sha256':prior['frozen_byte_sha256']['result'],
        'input_provenance':'same saved whole-core source, source factors and eight prior g3 source cells',
        'independent_source_matrix_fixture':f1,'independent_word1_green_fixture':f2,
        'global_low_source':enc(global_low),'result':enc(result),
        'positive_source_factors_not_materialized':True,'whole_core_intervals_not_point_values':True,
        'prior_Cauchy_and_allword_tail_removed':False,'m15_computed':False,'global_five_moment_closure':False}


def main():
    h,s,b=preflight();e=exclusive(ENTRY,{'component':COMPONENT,'source_sha256':s,'brief_sha256':b,'input_hashes':h,'entered':True})
    out={'component':COMPONENT,'source_sha256':s,'entry_sha256':e,'status':'failed'}
    try:out['output']=consume(h);out['status']='returned';out['error']=None
    except BaseException as exc:out['error']={'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()}
    digest=exclusive(RESULT,out);print(json.dumps({'component':COMPONENT,'status':out['status'],'result_sha256':digest}))
    return 0 if out['status']=='returned' else 1


if __name__=='__main__':raise SystemExit(main())
