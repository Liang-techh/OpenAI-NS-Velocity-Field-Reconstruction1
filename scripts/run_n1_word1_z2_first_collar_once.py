"""One-entry same-source no-reset first-collar W1_Z2 enclosure."""
from __future__ import annotations
import gzip,hashlib,importlib.util,json,sys,traceback
from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext

ROOT=Path(__file__).resolve().parents[1];BRIEF=ROOT.parent.parent/'RESUME_BRIEF_2026_10_16.json'
COMPONENT='n1_word1_z2_first_collar_v1';DIR=ROOT/'tmp'/COMPONENT
SOURCE,ENTRY,RESULT=(DIR/n for n in ('source.json','entered.json','result.json'))
BASE='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_'
SELECTED=BASE+'selected_words.json';COUPLED=BASE+'coupled_first_step.json.gz';COUPLED_CHECK=BASE+'coupled_first_step_check.json';COUPLED_SOURCE=BASE+'coupled_first_step.py'
INNER=BASE+'inner_volterra_functions.py';OP_SOURCE=BASE+'inner_operator.py'
CORE='experiments/root_st073/n1_core_word1_z2_v2_receipt.json.gz'
JOIN='experiments/root_st073/n1_word2_inlet_complete_v1_receipt.json.gz'
GAP='experiments/root_st073/n1_seam_gap_binding_v1_receipt.json'
NUMERIC='src/openai_ns_reconstruction/original_n1_selected_word_numeric.py'
LEAF='src/openai_ns_reconstruction/n1_word1_z2_first_collar.py'
FIXTURE='src/openai_ns_reconstruction/n1_word1_z2_first_collar_fixture.py'
RUNNER='scripts/run_n1_word1_z2_first_collar_once.py'
REQUIRED=(SELECTED,COUPLED,COUPLED_CHECK,COUPLED_SOURCE,INNER,OP_SOURCE,CORE,JOIN,GAP,NUMERIC,LEAF,FIXTURE,RUNNER)
FIXED={SELECTED:'0e04c87a011bfec74e0ceb0f822774a03495a08a839aaf4e378d21949736cf1e',
 COUPLED:'30b807af4147658ec0978ebd253cc1b2b31bed67bfed430a965e66461450add1',
 COUPLED_CHECK:'6a58caf1db41697620820830cbe7f76620874a8b429b351ad3d58f637e4e1416',
 COUPLED_SOURCE:'97913b82b8570b220921cd89918569a4019ff3495bdbea0548b3c53e9d8c93e7',
 INNER:'ff82f4edb03922a551cba3e2df20217a55dc80ed42bd39ccb8b77b83f63b0208',
 OP_SOURCE:'268217c6b33bc717472742fa64d76a853c6760a956abfad92ac45bcf385df014',
 CORE:'c53ff32ad3ce5288279953e84bc7417b71477c82b00d94466d0d2880adbc7ac9',
 JOIN:'bd1c753e41bcfced02e24090cd461299c1d6fd85dd70a52342f7ef750811f228',
 GAP:'a703914721b964dadca8d4ed3aa3d5df57f557c6974bdb54b9fd9bf3a0deaa11'}


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
    ck=d['source_checkpoint'];src=SOURCE.read_bytes()
    require(Path(ck['path']).resolve()==SOURCE.resolve() and sha(src)==ck['sha256'],'checkpoint')
    manifest=json.loads(src);require(manifest['component']==COMPONENT and manifest['one_entry_only'] is True,'manifest')
    hashes={}
    for item in manifest['files']:
        p=item['path'];require(p in REQUIRED and p not in hashes,'unexpected input');h=sha((ROOT/p).read_bytes())
        require(h==item['sha256'] and (p not in FIXED or h==FIXED[p]),'changed '+p);hashes[p]=h
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
    selected=json.loads(raw[SELECTED]);coupled=json.loads(gzip.decompress(raw[COUPLED]));check=json.loads(raw[COUPLED_CHECK])
    core=json.loads(gzip.decompress(raw[CORE]));join=json.loads(gzip.decompress(raw[JOIN]));gap=json.loads(raw[GAP])
    require(selected['checks']['passed'] is True and check['all_passed'] is True
        and core['result']['status']=='returned' and join['result']['status']=='returned' and gap['result']['status']=='returned','accepted inputs')
    require(check['input_hashes'][Path(COUPLED).name]==hashes[COUPLED],'coupled archive check binding')
    previous=core['result']['output'];binding=join['result']['output']['binding'];g=gap['result']['output']
    tail=selected['actual_selected_raw_word_tail'];words=selected['actual_selected_first_word_enclosure']
    measure=selected['actual_cutoff_complete_support_measure'];sample=coupled['samples']['actual_first_keep'];norm=sample['analytic_preconditioned_tail']
    family=previous['source_family']
    require(family==selected['source_family']==coupled['source_family']==sample['source_family']==g['source_family'], 'same source family')
    require(binding['same_selected_Cstar_exact_record'] is True and binding['same_Lambda_transfer_source'] is True
        and binding['same_core_exit_and_first_inlet'] is True,'same core amplitude and seam')
    require(previous['result']['coordinate']=='rho=4, core exit, real Z in [-1,1]'
        and previous['independent_source_matrix_fixture']['passed'] is True
        and previous['independent_word1_green_fixture']['passed'] is True,'core W1_Z2 result')
    require(measure['passed'] is True and measure['source_family']==family
        and measure['integration_q_interval']==['0','1'] and measure['measure_is_Z_independent'] is True
        and '0<hb<1/2' in measure['inequality_proof'], 'actual first collar geometry')
    require(words['passed'] is True and words['source_log_basis_names']==['hb','F0_anchor','Pstar_squared','Lambda']
        and words['source_family']==family and tail['passed'] is True,'hb source binding')
    require(norm['scope']=='uniform analytic series tail for canonical core/actual-first-collar problem; later completed-leading supports unbound'
        and norm['retained_preconditioned_indices']==[0,1]
        and norm['all_k_derivative_words_allowed'] is True,'actual analytic bound scope')
    require(all(g['checks'].values()) and g['prior_result_sha256']=='bf6ad485216d237166af8ea55261650925fbd3b64fee267af6c9bdd2910caaac', 'Cauchy gap binding')
    ntext=raw[NUMERIC].decode('utf-8');itext=raw[INNER].decode('utf-8');otext=raw[OP_SOURCE].decode('utf-8')
    require(all(x in ntext for x in ('H = Cp * M + C1 * M / loss','gap = loss / 2',
        'h0=C0 W0+B1 W0_Z is bounded on eta/8'))
        and 'complete_history_not_reset_at_core_seam=True' in itext
        and "logs=(self.source.logh,amplitude['logF0'],2*self.source.core.logP,self.source.core.logLambda)" in otext,
        'feedback bound, no-reset history and Z-independent hb')
    c=MPIntervalContext();c.prec=1800
    bits=tuple(measure['source_sc_exact_mpf_tuple']);require(len(bits)==4,'exact sc')
    sc=c.make_mpf((bits,bits));require(sc.a>0 and sc.b<c.mpf(1)/1024,'positive microscopic sc')
    M=iv(c,norm['initial_base_common_tube_bound']).b;Cp=iv(c,norm['off_kinematic_B0_norm_upper']).b
    C1=iv(c,norm['B1_norm_upper']).b;loss=iv(c,norm['total_available_Z_loss']).a
    radius=iv(c,norm['solution_Z_radius']);gap_record=sample['derivative_tail']['Cauchy_gap'];margin=iv(c,gap_record).a
    require(M>=0 and Cp>=0 and C1>=0 and loss>0 and margin>0,'analytic source positivity')
    require(iv(c,g['cauchy_gap'])._mpi_==iv(c,gap_record)._mpi_
        and iv(c,norm['total_available_Z_loss'])._mpi_==radius._mpi_, 'same eta/8 tube and eta/16 margin')
    H=(M*(Cp+C1/loss)).b;H2=(2*H/(margin*margin)).b
    require(H2>=0,'negative feedback bound')
    fields={}
    for field,terms in previous['result']['W1_Z2_core_exit'].items():
        fields[field]=[{**v,'coefficient':iv(c,v['coefficient'])} for v in terms]
    leaf=module(LEAF,'n1_word1_z2_collar_leaf');fixture=module(FIXTURE,'n1_word1_z2_collar_fixture')
    proof=fixture.run_fixture(leaf.enclose,c);require(proof['passed'] is True,'independent propagator and kernel fixture')
    result=leaf.enclose(c,fields,sc=sc,q=c.mpf([0,1]),feedback_Z2_upper=H2)
    return {'source_family':family,'source_logs':{'F0_of_Z':previous['source_logs']['F0_of_Z'],
          'Lambda':previous['source_logs']['Lambda'],'hb':words['source_log_bases'][0]},
        'same_source_binding':{'selected_Cstar':True,'Lambda_transfer':True,'core_first_seam':True,
          'actual_first_collar_uniform_norm':True,'Cauchy_gap_verified':True},
        'bound_derivation':{'H':'Mbase*(Cp+C1/loss)','H2':'2*H/(eta/16)^2',
          'kernel_change_of_variable':'u=2 exp(hb*v/2); du=hb*u*dv/2; K=u*(1-(u/X)^2)/2; log=u*log(X/u)',
          'Delta_x':'<=hb*s/(1-sc/16)', 'radius_not_a_signed_value':True},
        'directed_bounds':enc({'Mbase_upper':M,'Cp_upper':Cp,'C1_upper':C1,'loss_lower':loss,
          'gap_lower':margin,'feedback_upper':H,'feedback_Z2_upper':H2}),
        'independent_propagator_fixture':proof,'result':enc(result),
        'local_feedback_signed_value_computed':False,'allword_tail_added':False,
        'P1_is_not_m15':True,'m15_integral_computed':False,'global_five_moment_closure':False}


def main():
    h,s,b=preflight();e=exclusive(ENTRY,{'component':COMPONENT,'source_sha256':s,'brief_sha256':b,'input_hashes':h,'entered':True})
    out={'component':COMPONENT,'source_sha256':s,'entry_sha256':e,'status':'failed'}
    try:out['output']=consume(h);out['status']='returned';out['error']=None
    except BaseException as exc:out['error']={'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()}
    digest=exclusive(RESULT,out);print(json.dumps({'component':COMPONENT,'status':out['status'],'result_sha256':digest}))
    return 0 if out['status']=='returned' else 1


if __name__=='__main__':raise SystemExit(main())
