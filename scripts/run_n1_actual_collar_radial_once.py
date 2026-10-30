"""One-entry actual collar radial jets with source-specific flat cutoff."""
from __future__ import annotations
import gzip,hashlib,importlib.util,json,sys,traceback
from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext
ROOT=Path(__file__).resolve().parents[1];BRIEF=ROOT.parent.parent/'RESUME_BRIEF_2026_10_16.json'
COMPONENT='n1_actual_collar_radial_v1';DIR=ROOT/'tmp'/COMPONENT
SOURCE,ENTRY,RESULT=(DIR/n for n in ('source.json','entered.json','result.json'))
BASE='experiments/root_st073/lei_ren_part1_paper_compliant_'
CORE=BASE+'core_physical_field.json';CORE_CHECK=BASE+'core_physical_field_check.json';CORE_SOURCE=BASE+'core_physical_field.py'
FLAT_SOURCE=BASE+'flat_pulse_derivatives.py';COLLAR_SOURCE=BASE+'current_original_n1_lower_collar_functions.py'
PHASE='experiments/root_st073/n1_actual_collar_phase_v1_receipt.json.gz';MEAN='experiments/root_st073/n1_actual_collar_mean_v1_receipt.json.gz'
FLAT='src/openai_ns_reconstruction/n1_collar_flat_factor_certificate.py';FLAT_FIX='src/openai_ns_reconstruction/n1_collar_flat_factor_fixture.py'
RADIAL='src/openai_ns_reconstruction/n1_actual_collar_radial.py';RADIAL_FIX='src/openai_ns_reconstruction/n1_actual_collar_radial_fixture.py';RUNNER='scripts/run_n1_actual_collar_radial_once.py'
REQUIRED=(CORE,CORE_CHECK,CORE_SOURCE,FLAT_SOURCE,COLLAR_SOURCE,PHASE,MEAN,FLAT,FLAT_FIX,RADIAL,RADIAL_FIX,RUNNER)
FIXED={CORE:'ad000166f6fce1ea94a7beeccf3324536fbffca073487b873d5460133f380096',CORE_CHECK:'0814fdd1afc27cbd260ec713f834d05804efb7827d2d500cc3551aad8f94fb17',CORE_SOURCE:'19acdd5f7a506e7582427f004be8e2cd7ae22771eafc76b66ee44d053b9aec54',FLAT_SOURCE:'32f51123a47c2bafbb6c2f7056302bab996a019883615333875207d72a9bb7d7',COLLAR_SOURCE:'d00dfe285735844eabc985f99d94190caa1cfd187f668eb14a2686055cf175fa',PHASE:'47c2caaf3b8c7006984445681b39a1f01cf6f1bd4788384602bfa36ef3dd2b11',MEAN:'367a1e2c930907307bac66db309a3676807c1dff5dcf199c49024cb99a9789cf'}
def sha(b):return hashlib.sha256(b).hexdigest()
def require(ok,message):
    if not ok:raise ValueError(message)
def save(path,value):
    b=(json.dumps(value,sort_keys=True,separators=(',',':'))+'\n').encode()
    with path.open('xb') as f:f.write(b)
    return sha(b)
def preflight():
    require(not ENTRY.exists() and not RESULT.exists(),'already entered')
    b=BRIEF.read_bytes();d=json.loads(b);src=SOURCE.read_bytes();ck=d['source_checkpoint']
    require(d['component']==COMPONENT and d['phase']==COMPONENT+'_pending' and d['goal_active'] is True and d['full_goal_complete'] is False and d['one_entry_only'] is True and d['do_not_reenter'] is True,'brief guard')
    require(Path(ck['path']).resolve()==SOURCE.resolve() and ck['sha256']==sha(src),'source checkpoint')
    manifest=json.loads(src);require(manifest['component']==COMPONENT and manifest['one_entry_only'] is True,'manifest')
    hashes={}
    for item in manifest['files']:
        p=item['path'];require(p in REQUIRED and p not in hashes,'unexpected input')
        h=sha((ROOT/p).read_bytes());require(h==item['sha256'] and (p not in FIXED or h==FIXED[p]),'changed '+p);hashes[p]=h
    require(set(hashes)==set(REQUIRED),'source missing');return hashes,sha(src),sha(b)
def iv(c,record):
    lo=tuple(record['lower_exact_mpf_tuple']);hi=tuple(record['upper_exact_mpf_tuple'])
    require(len(lo)==len(hi)==4,'directed exact endpoints');x=c.make_mpf((lo,hi));require(x.a<=x.b,'reversed interval');return x
def enc(x):
    if hasattr(x,'_mpi_'):return {'lower_exact_mpf_tuple':list(x._mpi_[0]),'upper_exact_mpf_tuple':list(x._mpi_[1]),'directed':True}
    if isinstance(x,dict):return {str(k):enc(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [enc(v) for v in x]
    if isinstance(x,(str,int,float,bool)) or x is None:return x
    raise TypeError('unsupported output '+type(x).__name__)
def module(path,name):
    spec=importlib.util.spec_from_file_location(name,ROOT/path);m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m
def consume(hashes):
    raw={}
    for p in REQUIRED:
        b=(ROOT/p).read_bytes();require(sha(b)==hashes[p],'changed after entry '+p);raw[p]=b
    core=json.loads(raw[CORE]);cc=json.loads(raw[CORE_CHECK]);phase=json.loads(gzip.decompress(raw[PHASE]));mean=json.loads(gzip.decompress(raw[MEAN]))
    require(cc['all_passed'] is True and cc['input_hashes'][Path(CORE).name]==hashes[CORE]
      and phase['result']['status']=='returned' and mean['result']['status']=='returned','accepted actual inputs')
    p=phase['result']['output'];m=mean['result']['output'];family=p['source_family'];atoms=p['source_logs']
    require(m['source_family']==family and m['source_logs']==atoms and m['prior_result_sha256']==phase['result_sha256']
      and p['same_source_and_analytic_bound_receipts_verified'] is True and m['same_source_and_inlet_receipt_verified'] is True
      and p['actual_collar_F_and_U_enclosed_under_admitted_analytic_assumptions'] is True
      and m['actual_mean_enclosed_under_admitted_analytic_assumptions'] is True,'actual source binding')
    for key in ('actual_five_defect_family_sha256','implicit_source_sha256','datum_enclosure_sha256'):
        require(core[key]==family[key],'comparison core family '+key)
    whole=core['whole_core'];require(core['whole_Z_analytic_core_profile_enclosures_through_order5_available'] is True
      and whole['nonlinear_field_enclosed_not_replaced_by_model'] is True and whole['actual_selected_Cstar_source_retained'] is True,'admitted comparison grids')
    text=raw[FLAT_SOURCE].decode('utf8');lower=raw[COLLAR_SOURCE].decode('utf8');source=raw[CORE_SOURCE].decode('utf8')
    require('odds=1/(1-x)**2-1/x**2' in text and 'L1=2/(1-x)**3+2/x**3' in text
      and 'sigma=sigma or (lambda x:own.jets.complete.Gamma.sigma(g,x))' in lower
      and 'chi=lambda x:g.add(g.one,g.mul(g.sub(hb,g.one),sigma(x)))' in lower
      and 'Vgrid[gridkey(i,k)]=base+self.epsilon*psi[gridkey(i,k)]' in source,'actual flat cutoff and physical core U definitions')
    c=MPIntervalContext();c.prec=1800;sc=iv(c,p['sc']);loghb=iv(c,atoms['hb']);Z=iv(c,whole['Z'])
    require(iv(c,m['sc'])._mpi_==sc._mpi_ and Z.a==-1 and Z.b==1,'same exact sc and real Z')
    def terms(rows):return [{**v,'coefficient':iv(c,v['coefficient'])} for v in rows]
    actual_phi={(0,k):terms(p['enclosure']['fields']['Phi'][str(k)]) for k in range(6)}
    actual_U={(0,k):terms(p['enclosure']['fields']['U'][str(k)]) for k in range(6)}
    actual_mean={(0,k):terms(m['mean']['mean_Z'][str(k)]['all_terms']) for k in range(6)}
    cmp_phi={(i,k):iv(c,whole['Phi_derivatives']['rho'+str(i)+'_Z'+str(k)]) for i in range(3) for k in range(6-i)}
    lam=iv(c,core['actual_Lambda']);require(lam.a>0,'same positive Lambda')
    cmp_U={}
    for i in range(3):
        for k in range(6-i):
            v=iv(c,whole['Uz_derivatives']['rho'+str(i)+'_Z'+str(k)])
            cmp_U[(i,k)]=(v if i==0 else [{'powers':(0,0,0,-1),'coefficient':v*lam,
              'origin':('same_core_U_radial_interval_times_same_Lambda',i,k),
              'semantics':'saved_actual_core_interval_formally_refactored_not_point_value'}])
    flat=module(FLAT,'n1_flat_factor');flat_fix=module(FLAT_FIX,'n1_flat_factor_fixture')
    radial=module(RADIAL,'n1_actual_radial');radial_fix=module(RADIAL_FIX,'n1_actual_radial_fixture')
    f=flat_fix.run(flat,c);require(f['passed'] is True,'flat fixture failed '+str([k for k,v in f['checks'].items() if not v]))
    r=radial_fix.run(radial,c);require(r['passed'] is True,'radial fixture failed '+str([k for k,v in r['checks'].items() if not v]))
    source_contract={'exact_original_flat_sigma':True,'same_hb_log_and_sc':True,'hb_and_phase_Z_independent':True}
    factor=flat.build(c,sc=sc,log_hb=loghb,source_family=family,source_atoms=atoms,source_contract=source_contract)
    geom={k:True for k in ('hb_strictly_between_zero_and_half','sc_phase_Z_independent','hb_epsilon_Z_independent','sigma_Z_independent')}
    actual=radial.build_collar_laurent_radial(c,factor['rho'],Z,loghb,actual_Phi=actual_phi,actual_U=actual_U,own_mean=actual_mean,
      comparison_Phi=cmp_phi,comparison_U=cmp_U,chi=factor['chi'],chi_r=factor['chi_r'],source_family=family,source_atoms=atoms,geometry_certificate=geom,
      source_provenance={'phase_result_sha256':phase['result_sha256'],'mean_result_sha256':mean['result_sha256'],'comparison_core_sha256':hashes[CORE]})
    return {'source_family':family,'source_logs':atoms,'domain':'rho=4 exp(hb*s), 0<=s<=sc/4, real Z in [-1,1]',
      'source_receipts_and_conditional_common_core_checked':True,'flat_factor':enc(factor),'flat_fixture':f,'radial_fixture':r,
      'actual_radial_sector_enclosures':enc({'Phi':actual['Phi'],'U':actual['U'],'mean':actual['mean']}),
      'actual_radial_metadata':enc(actual['metadata']),
      'comparison_U_epsilon_already_in_saved_intervals':True,'same_Lambda_refactored_into_formal_inverse':True,
      'refactoring_preserves_interval_uncertainty_and_may_widen':True,
      'inherited_interval_errors_preserved_in_sector_origins':True,'actual_signed_point_values_computed':False,
      'm15_computed':False,'allword_tail_added':False,'terminal_closure_proved':False}
def main():
    hashes,source,b=preflight();entry=save(ENTRY,{'component':COMPONENT,'source_sha256':source,'brief_sha256':b,'input_hashes':hashes,'entered':True})
    out={'component':COMPONENT,'source_sha256':source,'entry_sha256':entry,'status':'failed'}
    try:out.update(status='returned',output=consume(hashes),error=None)
    except BaseException as exc:out['error']={'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()}
    digest=save(RESULT,out);print(json.dumps({'component':COMPONENT,'status':out['status'],'result_sha256':digest}))
    return 0 if out['status']=='returned' else 1
if __name__=='__main__':raise SystemExit(main())
