"""One-entry actual same-source collar Phi/U phase enclosures."""
from __future__ import annotations
import gzip,hashlib,importlib.util,json,sys,traceback
from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext
ROOT=Path(__file__).resolve().parents[1];BRIEF=ROOT.parent.parent/'RESUME_BRIEF_2026_10_16.json'
COMPONENT='n1_actual_collar_phase_v1';DIR=ROOT/'tmp'/COMPONENT
SOURCE,ENTRY,RESULT=(DIR/n for n in ('source.json','entered.json','result.json'))
BASE='experiments/root_st073/lei_ren_part1_paper_compliant_'
CORE=BASE+'core_physical_field.json';CORE_CHECK=BASE+'core_physical_field_check.json'
TRANSFER=BASE+'core_transfer.json';TRANSFER_CHECK=BASE+'core_transfer_check.json'
DOMAIN=BASE+'current_original_n1_inner_analytic_domain.json.gz';DOMAIN_CHECK=BASE+'current_original_n1_inner_analytic_domain_check.json';DOMAIN_SOURCE=BASE+'current_original_n1_inner_analytic_domain.py'
SELECTED=BASE+'current_original_n1_selected_words.json'
BINDING='experiments/root_st073/n1_collar_laurent_binding_v3_receipt.json.gz'
LEAF='src/openai_ns_reconstruction/n1_collar_actual_phase_enclosure.py';FIXTURE='src/openai_ns_reconstruction/n1_collar_actual_phase_fixture.py';RUNNER='scripts/run_n1_actual_collar_phase_once.py'
REQUIRED=(CORE,CORE_CHECK,TRANSFER,TRANSFER_CHECK,DOMAIN,DOMAIN_CHECK,DOMAIN_SOURCE,SELECTED,BINDING,LEAF,FIXTURE,RUNNER)
FIXED={CORE:'ad000166f6fce1ea94a7beeccf3324536fbffca073487b873d5460133f380096',CORE_CHECK:'0814fdd1afc27cbd260ec713f834d05804efb7827d2d500cc3551aad8f94fb17',TRANSFER:'4d7ec9b08ee2d36970d22c2a3287b8619bb1fa3deebf5ea264ab0d44b870d497',TRANSFER_CHECK:'91883a7c3d92b1c01ca0d66c5852a5d65814aa364dc7825adb8dd9323c5cb848',DOMAIN:'c9905f7a9d13d4f651d1f55036932ea6c69db4fc5fa890e67a9839e611af4db7',DOMAIN_CHECK:'01a4a5f497a224d3dd6ebcb7b44809df5b3702e8d32028170a764903238e0ae2',DOMAIN_SOURCE:'56507e29d1c069c6ca71c13722717dfaeef662265401652afb3b62232bfe73f3',SELECTED:'0e04c87a011bfec74e0ceb0f822774a03495a08a839aaf4e378d21949736cf1e',BINDING:'2545e3688129f4fd0733849f573f7f6684b843a1720300aea7fe898d90d1bc23'}
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
    core=json.loads(raw[CORE]);cc=json.loads(raw[CORE_CHECK]);transfer=json.loads(raw[TRANSFER]);tc=json.loads(raw[TRANSFER_CHECK])
    domain=json.loads(gzip.decompress(raw[DOMAIN]));dc=json.loads(raw[DOMAIN_CHECK]);selected=json.loads(raw[SELECTED]);binding=json.loads(gzip.decompress(raw[BINDING]))
    require(cc['all_passed'] is True and tc['all_passed'] is True and dc['all_passed'] is True
      and dc['current_original_n1_core_and_actual_first_collar_analytic_bounds_installed'] is True
      and selected['checks']['passed'] is True and binding['result']['status']=='returned','accepted source receipts')
    for check,path in ((cc,CORE),(tc,TRANSFER),(dc,DOMAIN)):
        require(check['input_hashes'][Path(path).name]==hashes[path],'source receipt input '+path)
    require(core['input_hashes'][Path(TRANSFER).name]==hashes[TRANSFER],'same core transfer')
    family=domain['source_family'];require(family==selected['source_family']==binding['result']['output']['saved_proof_binding']['source_family'],'actual source family')
    for key in ('actual_five_defect_family_sha256','implicit_source_sha256','datum_enclosure_sha256'):
        require(core[key]==family[key],'core family '+key)
    require(binding['result']['output']['saved_proof_binding']['same_core_family_as_n1_receipt'] is True
      and domain['source_inner_common_holomorphic_bounds_constructed'] is True
      and domain['all_fixed_radial_derivatives_holomorphic_on_same_domain_by_smooth_ODE'] is True
      and domain['actual_inner_geometry']['actual_bridge_not_replaced_by_comparison'] is True,'admitted collar holomorphic source')
    words=selected['actual_selected_first_word_enclosure'];measure=selected['actual_cutoff_complete_support_measure']
    require(words['passed'] is True and words['source_family']==family and measure['passed'] is True and measure['source_family']==family
      and measure['integration_q_interval']==['0','1'] and measure['measure_is_Z_independent'] is True
      and '0<hb<1/2' in measure['inequality_proof'],'same exact width and selected phase')
    text=raw[DOMAIN_SOURCE].decode('utf8')
    require('actual_collar_radial_bounds=export({key:collar[key] for key in (\'Phi\',\'Uz\',\'mean\')})' in text
      and 'ordinary_Z_derivative_rule=\'for fields on the Gg tube use k!/(eta/4)^k' in text
      and "exact_actual_collar_identities=dict(" in text,'actual analytic bound derivation')
    c=MPIntervalContext();c.prec=1800
    whole=core['whole_core'];rho=iv(c,whole['rho']);Z=iv(c,whole['Z'])
    require(core['whole_Z_analytic_core_profile_enclosures_through_order5_available'] is True
      and whole['nonlinear_field_enclosed_not_replaced_by_model'] is True
      and whole['actual_selected_Cstar_source_retained'] is True and rho.a<=0 and rho.b>=c.mpf('4.1').a and Z.a==-1 and Z.b==1,'same nonlinear whole-core enclosure')
    loghb=words['source_log_bases'][0];require(iv(c,loghb)._mpi_==iv(c,domain['actual_inner_geometry']['source_log_hb'])._mpi_,'same hb logarithm')
    require(iv(c,transfer['Lambda'])._mpi_==iv(c,core['actual_Lambda'])._mpi_,'same Lambda')
    bits=tuple(measure['source_sc_exact_mpf_tuple']);require(len(bits)==4,'exact sc bits');sc=c.make_mpf((bits,bits))
    require(sc.a>0 and sc.b<c.mpf(1)/1024,'same sc admissible')
    phase=c.make_mpf((c.mpf(0)._mpi_[0],(sc/4).b._mpi_[1]));radius=iv(c,domain['common_tubes']['Gg_Z_radius'])
    rho_max=iv(c,domain['actual_inner_geometry']['real_rho_upper_bound']);require(radius.a>0 and rho_max.a>=4,'Cauchy radius and geometry')
    corejets={name:{k:iv(c,whole[key]['rho0_Z'+str(k)]) for k in range(6)} for name,key in (('Phi','Phi_derivatives'),('U','Uz_derivatives'))}
    mean0={k:iv(c,whole['radial_average_Uz_derivatives']['rho0_Z'+str(k)]) for k in range(6)}
    bounds={}
    for name,key in (('Phi','Phi'),('U','Uz')):
        data=domain['actual_collar_radial_bounds'][key][1];bounds[name]={int(p):iv(c,v['coefficient']) for p,v in data.items()}
        require(bounds[name] and all(p>=0 and val.a>=0 for p,val in bounds[name].items()),'collar radial magnitude source')
    atoms={'hb':loghb,'F0_anchor':whole['F0_exact_positive_log_enclosure'],
      'Pstar_squared':words['source_log_bases'][2],'Lambda':transfer['logLambda']}
    certificate={k:True for k in ('source_family_checked','core_enclosure_contains_rho4','complex_radial_bound_covers_collar','analytic_Z_disks_cover_real_domain','hb_Z_independent','rho_upper_covers_collar')}
    leaf=module(LEAF,'n1_actual_phase_leaf');fixture=module(FIXTURE,'n1_actual_phase_fixture')
    proof=fixture.run(leaf,c);require(proof['passed'] is True,'synthetic fixture failed: '+str([k for k,v in proof['checks'].items() if not v]))
    result=leaf.build(c,sc=sc,phase=phase,core_inlet_Z=corejets,radial_tube_bounds=bounds,cauchy_radius=radius,
      rho_upper=rho_max,source_family=family,source_atoms=atoms,admitted_source_contract=certificate)
    return {'source_family':family,'same_source_and_analytic_bound_receipts_verified':True,
       'domain':'rho=4 exp(hb*s), 0<=s<=sc/4, real Z in [-1,1]',
       'source_logs':atoms,'sc':enc(sc),'Cauchy_radius':enc(radius),'rho_upper':enc(rho_max),
       'radial_tube_bounds':enc(bounds),'own_mean_core_inlet':enc(mean0),'fixture':proof,'enclosure':enc(result),
       'actual_signed_point_value_computed':False,'separate_cauchy_error_retained':True,
       'actual_collar_F_and_U_enclosed_under_admitted_analytic_assumptions':True,
       'allword_tail_added':False,'m15_computed':False,'terminal_closure_proved':False}
def main():
    hashes,source,b=preflight();entry=save(ENTRY,{'component':COMPONENT,'source_sha256':source,'brief_sha256':b,'input_hashes':hashes,'entered':True})
    out={'component':COMPONENT,'source_sha256':source,'entry_sha256':entry,'status':'failed'}
    try:out.update(status='returned',output=consume(hashes),error=None)
    except BaseException as exc:out['error']={'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()}
    digest=save(RESULT,out);print(json.dumps({'component':COMPONENT,'status':out['status'],'result_sha256':digest}))
    return 0 if out['status']=='returned' else 1
if __name__=='__main__':raise SystemExit(main())
