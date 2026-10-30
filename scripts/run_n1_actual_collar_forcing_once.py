"""One-entry actual source-factored collar forcing and dynamic matrices."""
from __future__ import annotations
import gzip,hashlib,importlib.util,json,sys,traceback
from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext
ROOT=Path(__file__).resolve().parents[1];BRIEF=ROOT.parent.parent/'RESUME_BRIEF_2026_10_16.json'
COMPONENT='n1_actual_collar_forcing_v1';DIR=ROOT/'tmp'/COMPONENT
SOURCE,ENTRY,RESULT=(DIR/n for n in ('source.json','entered.json','result.json'))
BASE='experiments/root_st073/lei_ren_part1_paper_compliant_'
CORE=BASE+'core_physical_field.json';CORE_CHECK=BASE+'core_physical_field_check.json';TRANSFER=BASE+'core_transfer.json';TRANSFER_CHECK=BASE+'core_transfer_check.json'
ORIGINAL=BASE+'current_original_n1_inner_volterra_functions.py'
RADIAL='experiments/root_st073/n1_actual_collar_radial_v1_receipt.json.gz'
LEAF='src/openai_ns_reconstruction/n1_actual_collar_forcing_compact.py';FIXTURE='src/openai_ns_reconstruction/n1_actual_collar_forcing_fixture.py';RUNNER='scripts/run_n1_actual_collar_forcing_once.py'
REQUIRED=(CORE,CORE_CHECK,TRANSFER,TRANSFER_CHECK,ORIGINAL,RADIAL,LEAF,FIXTURE,RUNNER)
FIXED={CORE:'ad000166f6fce1ea94a7beeccf3324536fbffca073487b873d5460133f380096',CORE_CHECK:'0814fdd1afc27cbd260ec713f834d05804efb7827d2d500cc3551aad8f94fb17',TRANSFER:'4d7ec9b08ee2d36970d22c2a3287b8619bb1fa3deebf5ea264ab0d44b870d497',TRANSFER_CHECK:'91883a7c3d92b1c01ca0d66c5852a5d65814aa364dc7825adb8dd9323c5cb848',ORIGINAL:'ff82f4edb03922a551cba3e2df20217a55dc80ed42bd39ccb8b77b83f63b0208',RADIAL:'27d0e2a97568c50ff62f373cbdf04dea262b4980392d7ac4b13b9c7db21752cf'}
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
    core=json.loads(raw[CORE]);cc=json.loads(raw[CORE_CHECK]);transfer=json.loads(raw[TRANSFER]);tc=json.loads(raw[TRANSFER_CHECK]);radial=json.loads(gzip.decompress(raw[RADIAL]))
    require(cc['all_passed'] is True and tc['all_passed'] is True and radial['result']['status']=='returned'
      and cc['input_hashes'][Path(CORE).name]==hashes[CORE] and tc['input_hashes'][Path(TRANSFER).name]==hashes[TRANSFER]
      and core['input_hashes'][Path(TRANSFER).name]==hashes[TRANSFER],'accepted source inputs')
    saved=radial['result']['output'];family=saved['source_family'];atoms=saved['source_logs']
    require(saved['source_receipts_and_conditional_common_core_checked'] is True
      and saved['flat_factor']['source_proof']['sigma_and_derivative_le_hb_squared'] is True
      and saved['inherited_interval_errors_preserved_in_sector_origins'] is True
      and saved['actual_signed_point_values_computed'] is False,'actual radial source and errors')
    for key in ('actual_five_defect_family_sha256','implicit_source_sha256','datum_enclosure_sha256'):
        require(core[key]==family[key],'same source family '+key)
    require(transfer['implicit_source_sha256']==core['implicit_source_sha256'] and transfer['datum_enclosure_sha256']==core['datum_enclosure_sha256'],'same transfer')
    original=raw[ORIGINAL].decode('utf8')
    require('def raw_matrix_expressions():' in original and "'6,1':-eps**2*8*z*r*f/L" in original
      and "'5,1':eps*2*H/L" in original and 'endpoint_independent_of_Z=True' in original,'original raw dynamic system')
    c=MPIntervalContext();c.prec=1800;rho=iv(c,saved['flat_factor']['rho']);x=iv(c,saved['flat_factor']['x'])
    z=c.mpf([-1,1]);delta=iv(c,core['actual_delta']);lam=iv(c,core['actual_Lambda']);j=iv(c,transfer['required_j']);logF0=iv(c,atoms['F0_anchor'])
    require(lam.a>0 and j.a>0 and iv(c,transfer['Lambda'])._mpi_==lam._mpi_,'same positive Lambda and j')
    source_index={};grids={};classification_counts={}
    markers=(('Cauchy_phase','symmetric_certified_error_not_field_value'),('Cauchy_phase','cauchy_phase_variation_error'),
      ('mean_phase','U_phase_error_history_integral'),('cutoff_value','minus_sigma_factor'),('cutoff_value','hb_sigma_factor'),
      ('cutoff_derivative','sigma_prime_factor'))
    for field in ('Phi','U','mean'):
        grid={};radial_limit=2 if field!='mean' else 1
        for i in range(radial_limit+1):
            for k in range(6-i):
                key='('+str(i)+', '+str(k)+')';terms=saved['actual_radial_sector_enclosures'][field][key];group=[]
                for index,t in enumerate(terms):
                    sid=field+':'+str(i)+':'+str(k)+':'+str(index);origin_json=json.dumps(t['origin'],sort_keys=True,separators=(',',':'))
                    classes=sorted({label for label,needle in markers if needle in origin_json});uncertain=bool(classes)
                    source_index[sid]={'field':field,'row':key,'term_index':index,'origin_sha256':sha(origin_json.encode()),'error_classes':classes}
                    for label in classes:classification_counts[label]=classification_counts.get(label,0)+1
                    powers=t['powers'];require(len(powers)==4 and all(isinstance(a,int) and not isinstance(a,bool) for a in powers),'explicit source powers')
                    group.append({'powers':powers,'coefficient':iv(c,t['coefficient']),
                      'origin':{'source_ids':[sid],'uncertainty':uncertain,'error_classes':classes}})
                grid[(i,k)]=group
        grids[field]=grid
    require(all(classification_counts.get(k,0)>0 for k in ('Cauchy_phase','cutoff_value','cutoff_derivative','mean_phase')),'independent uncertainty provenance classes')
    leaf=module(LEAF,'n1_actual_collar_forcing_leaf');fixture=module(FIXTURE,'n1_actual_collar_forcing_fixture')
    proof=fixture.run(leaf,c);require(proof['passed'] is True,'independent fixture failed '+str([k for k,v in proof['checks'].items() if not v]))
    result=leaf.build_first_collar_laurent_forcing_jets(c,rho,z,delta,lam,j,Phi_derivatives=grids['Phi'],Uz_derivatives=grids['U'],radial_average_Uz_derivatives=grids['mean'],log_F0=logF0,x=x)
    allowed=set(source_index);seen=set();uncertain=0
    collections=[terms for rows in result['sectors'].values() for terms in rows.values()]
    collections += [terms for group in result['matrix'].values() for rows in group.values() for terms in rows.values()]
    for terms in collections:
        for term in terms:
            origin=term['origin'];ids=origin.get('source_ids',[]);require(all(q in allowed for q in ids),'unmapped input provenance')
            seen.update(ids)
            if origin.get('uncertainty'):uncertain+=1;require(bool(origin.get('error_classes')),'lost error class')
    require(uncertain>0 and seen,'no source and uncertainty outputs')
    return {'source_family':family,'source_logs':atoms,'domain':saved['domain'],
      'same_source_radial_and_core_receipts_verified':True,'prior_radial_result_sha256':radial['result_sha256'],
      'source_index_receipt_sha256':hashes[RADIAL],'source_index':source_index,'input_uncertainty_class_counts':classification_counts,
      'fixture':proof,'g':enc(result['sectors']),'matrix':enc(result['matrix']),'calculation_metadata':enc(result['proof']),
      'output_uncertain_group_count':uncertain,'distinct_input_ids_referenced':len(seen),
      'same_power_groups_are_outward_interval_sums':True,'source_or_error_groups_not_interpreted_as_point_values':True,
      'no_endpoint_identity_or_error_cancellation_proved':True,'allword_tail_added':False,'m15_computed':False,'terminal_closure_proved':False}
def main():
    hashes,source,b=preflight();entry=save(ENTRY,{'component':COMPONENT,'source_sha256':source,'brief_sha256':b,'input_hashes':hashes,'entered':True})
    out={'component':COMPONENT,'source_sha256':source,'entry_sha256':entry,'status':'failed'}
    try:out.update(status='returned',output=consume(hashes),error=None)
    except BaseException as exc:out['error']={'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()}
    digest=save(RESULT,out);print(json.dumps({'component':COMPONENT,'status':out['status'],'result_sha256':digest}))
    return 0 if out['status']=='returned' else 1
if __name__=='__main__':raise SystemExit(main())
