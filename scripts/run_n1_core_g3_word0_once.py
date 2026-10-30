"""One-entry real-core g_Z3 and W0_Z3 interval consumer; no ancestor replay."""
from __future__ import annotations
import gzip
import hashlib
import importlib.util
import json
import sys
import traceback
from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext

ROOT=Path(__file__).resolve().parents[1]
BRIEF=ROOT.parent.parent/'RESUME_BRIEF_2026_10_16.json'
COMPONENT='n1_core_g3_word0_v1'; DIR=ROOT/'tmp'/COMPONENT
SOURCE,ENTRY,RESULT=(DIR/name for name in ('source.json','entered.json','result.json'))
BASE='experiments/root_st073/lei_ren_part1_paper_compliant_'
CORE=BASE+'core_physical_field.json'; CORE_CHECK=BASE+'core_physical_field_check.json'; CORE_SOURCE=BASE+'core_physical_field.py'
TRANSFER=BASE+'core_transfer.json'; TRANSFER_CHECK=BASE+'core_transfer_check.json'; TRANSFER_SOURCE=BASE+'core_transfer.py'
INNER='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_inner_volterra_functions.py'
PRIOR='experiments/root_st073/n1_word2_source_inlet_v1_receipt.json.gz'
LEAF='src/openai_ns_reconstruction/n1_core_forcing_z3.py'
GREEN='src/openai_ns_reconstruction/n1_core_word0_z3_green.py'
FIXTURE='src/openai_ns_reconstruction/n1_core_z3_independent_fixture.py'
RUNNER='scripts/run_n1_core_g3_word0_once.py'
REQUIRED=(CORE,CORE_CHECK,CORE_SOURCE,TRANSFER,TRANSFER_CHECK,TRANSFER_SOURCE,INNER,PRIOR,LEAF,GREEN,FIXTURE,RUNNER)
FIXED={CORE:'ad000166f6fce1ea94a7beeccf3324536fbffca073487b873d5460133f380096',
 CORE_CHECK:'0814fdd1afc27cbd260ec713f834d05804efb7827d2d500cc3551aad8f94fb17',
 CORE_SOURCE:'19acdd5f7a506e7582427f004be8e2cd7ae22771eafc76b66ee44d053b9aec54',
 TRANSFER:'4d7ec9b08ee2d36970d22c2a3287b8619bb1fa3deebf5ea264ab0d44b870d497',
 TRANSFER_CHECK:'91883a7c3d92b1c01ca0d66c5852a5d65814aa364dc7825adb8dd9323c5cb848',
 TRANSFER_SOURCE:'2ad468af4815afe650eb3ebedb29aa392062c4b8a9c7b29e663540a1a2d199d6',
 INNER:'ff82f4edb03922a551cba3e2df20217a55dc80ed42bd39ccb8b77b83f63b0208',
 PRIOR:'fed8a184174e2922ea4097316ff503aba740f2423913d3e30b254eb69c2cc264'}


def sha(x):return hashlib.sha256(x).hexdigest()
def require(x,message):
    if not x:raise ValueError(message)
def exclusive(path,obj):
    b=(json.dumps(obj,sort_keys=True,separators=(',',':'))+'\n').encode()
    with path.open('xb') as f:f.write(b)
    return sha(b)


def preflight():
    require(not ENTRY.exists() and not RESULT.exists(),'already entered')
    b=BRIEF.read_bytes(); brief=json.loads(b)
    require(brief.get('component')==COMPONENT and brief.get('phase')==COMPONENT+'_pending'
        and brief.get('goal_active') is True and brief.get('full_goal_complete') is False
        and brief.get('one_entry_only') is True and brief.get('do_not_reenter') is True,'brief guard')
    ck=brief['source_checkpoint']; raw=SOURCE.read_bytes()
    require(Path(ck['path']).resolve()==SOURCE.resolve() and sha(raw)==ck['sha256'],'checkpoint')
    source=json.loads(raw);require(source['component']==COMPONENT and source['one_entry_only'] is True,'manifest')
    hashes={}
    for row in source['files']:
        p=row['path'];require(p in REQUIRED and p not in hashes,'unexpected input')
        digest=sha((ROOT/p).read_bytes())
        require(digest==row['sha256'] and (p not in FIXED or digest==FIXED[p]),'changed input '+p)
        hashes[p]=digest
    require(set(hashes)==set(REQUIRED),'missing input')
    return hashes,sha(raw),sha(b)


def iv(c,record):
    lo=tuple(record['lower_exact_mpf_tuple']); hi=tuple(record['upper_exact_mpf_tuple'])
    require(len(lo)==len(hi)==4,'four-entry directed endpoints required')
    x=c.make_mpf((lo,hi));require(x.a<=x.b,'reversed interval');return x


def encode(x):
    if hasattr(x,'_mpi_'):
        return {'lower_exact_mpf_tuple':list(x._mpi_[0]),'upper_exact_mpf_tuple':list(x._mpi_[1]),'directed':True}
    if isinstance(x,dict):return {str(k):encode(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [encode(v) for v in x]
    if isinstance(x,(str,int,float,bool)) or x is None:return x
    raise TypeError('unsupported output '+type(x).__name__)


def module(path,name):
    spec=importlib.util.spec_from_file_location(name,ROOT/path)
    require(spec is not None and spec.loader is not None,'leaf import '+path)
    obj=importlib.util.module_from_spec(spec);sys.modules[spec.name]=obj;spec.loader.exec_module(obj);return obj


def consume(hashes):
    raw={}
    for p in REQUIRED:
        data=(ROOT/p).read_bytes();require(sha(data)==hashes[p],'changed after entry '+p);raw[p]=data
    core=json.loads(raw[CORE]); check=json.loads(raw[CORE_CHECK]); transfer=json.loads(raw[TRANSFER]); tc=json.loads(raw[TRANSFER_CHECK]); prior=json.loads(gzip.decompress(raw[PRIOR]))
    require(check['all_passed'] is True and tc['all_passed'] is True and prior['result']['status']=='returned','accepted input checks')
    require(check['input_hashes'][Path(CORE).name]==hashes[CORE] and tc['input_hashes'][Path(TRANSFER).name]==hashes[TRANSFER]
        and core['input_hashes'][Path(TRANSFER).name]==hashes[TRANSFER], 'archive/check/transfer binding')
    family=prior['result']['output']['source_family']
    for key in ('actual_five_defect_family_sha256','implicit_source_sha256','datum_enclosure_sha256'):
        require(core[key]==family[key], 'source family mismatch '+key)
    require(transfer['implicit_source_sha256']==core['implicit_source_sha256'] and transfer['datum_enclosure_sha256']==core['datum_enclosure_sha256'], 'transfer family')
    whole=core['whole_core']
    require(core['whole_Z_analytic_core_profile_enclosures_through_order5_available'] is True
        and whole['nonlinear_field_enclosed_not_replaced_by_model'] is True
        and whole['actual_selected_Cstar_source_retained'] is True, 'actual analytic core enclosure')
    cs=raw[CORE_SOURCE].decode('utf-8'); ts=raw[TRANSFER_SOURCE].decode('utf-8'); inner=raw[INNER].decode('utf-8')
    require(all(s in cs for s in ('gradient=L*H/denominator', 'ell=[-self.Lambda*gradient[k] for k in range(5)]',
        "self.j=read(major,'required_j')", 'self.sigma=self.j/500', "whole=self.profiles([0,'4.1'],[-1,1])"))
        and 'lam=ctx.exp(logLambda);eps=ctx.exp(-logLambda)' in ts
        and 'endpoint_independent_of_Z=True' in inner and 'kernels=(K,log,g.neg(K),g.one,g.mul(t,t,t),t)' in inner,
        'source amplitude/Green/scale identities')
    c=MPIntervalContext(); c.prec=1800
    rho_all=iv(c,whole['rho']); z=iv(c,whole['Z']); delta=iv(c,core['actual_delta']); lam=iv(c,core['actual_Lambda']); j=iv(c,transfer['required_j'])
    require(rho_all.a==0 and rho_all.b>=4 and z.a==-1 and z.b==1 and lam.a>0 and j.a>0,'core coverage')
    require(iv(c,transfer['Lambda'])._mpi_==lam._mpi_, 'exact same Lambda source interval')
    Flog=iv(c,whole['F0_exact_positive_log_enclosure']); Llog=iv(c,transfer['logLambda'])
    grids=[]
    for name in ('Phi_derivatives','Uz_derivatives','radial_average_Uz_derivatives'):
        grid=whole[name]; grids.append({f'rho{i}_Z{k}':iv(c,grid[f'rho{i}_Z{k}']) for i in range(3) for k in range(6-i)})
    force=module(LEAF,'n1_core_forcing_z3_leaf'); green=module(GREEN,'n1_core_word0_z3_green_leaf'); fixture=module(FIXTURE,'n1_core_z3_independent_fixture')
    fixture_result=fixture.run_fixture(force.build_core_forcing_z3,green.integrate,c)
    require(fixture_result['passed'] is True,'independent symbolic fixture')
    cells=[]; proof=None
    for index in range(8):
        left=c.mpf(index)/8;right=c.mpf(index+1)/8;t=c.make_mpf((left._mpi_[0],right._mpi_[1]));rho=4*t*t
        rows=force.build_core_forcing_z3(c,rho,z,delta,lam,j,*grids,Flog)
        if proof is None:proof=rows['proof']
        cells.append({'t':t,'rho':rho,'g':{'4':rows['g4_Z3'],'5':rows['g5_Z3'],'6':rows['g6_Z3']}})
    word=green.integrate(c,cells)
    return {'source_family':family,'source_domain':{'rho':'[0,4.1]','Z':'[-1,1]',
          'cell_inputs':'same actual whole-core enclosures reused; not independently sampled cell values'},
        'source_logs':{'F0_of_Z':encode(Flog),'Lambda':encode(Llog)},
        'ordinary_Z_order':3,'independent_symbolic_fixture':fixture_result,'forcing_proof':encode(proof),
        'cells':encode(cells),'W0_Z3_green_enclosure':encode(word),
        'actual_function_values_not_midpoints':True,'same_actual_analytic_core':True,
        'W1_Z2_integrated':False,'first_collar_added':False,'all_word_tail_added':False,
        'm15_computed':False,'global_five_moment_closure':False}


def main():
    h,s,b=preflight();e=exclusive(ENTRY,{'component':COMPONENT,'source_sha256':s,'brief_sha256':b,'input_hashes':h,'entered':True})
    out={'component':COMPONENT,'source_sha256':s,'entry_sha256':e,'status':'failed'}
    try:out['output']=consume(h);out['status']='returned';out['error']=None
    except BaseException as exc:out['error']={'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()}
    digest=exclusive(RESULT,out);print(json.dumps({'component':COMPONENT,'status':out['status'],'result_sha256':digest}))
    return 0 if out['status']=='returned' else 1


if __name__=='__main__':raise SystemExit(main())
