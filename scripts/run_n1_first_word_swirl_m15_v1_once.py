"""One-entry physical first-collar swirl m15 contribution of finite W1."""
from __future__ import annotations
import gzip,hashlib,importlib.util,json,sys,traceback
from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext
ROOT=Path(__file__).resolve().parents[1];BRIEF=ROOT.parent.parent/'RESUME_BRIEF_2026_10_16.json';COMPONENT='n1_first_word_swirl_m15_v1';DIR=ROOT/'tmp'/COMPONENT
SOURCE,ENTRY,RESULT=(DIR/n for n in ('source.json','entered.json','result.json'))
JOIN='experiments/root_st073/n1_core_collar_low_join_v2_receipt.json.gz';RADIAL='experiments/root_st073/n1_actual_collar_radial_v1_receipt.json.gz';TAIL='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_selected_words.json'
INNER='src/openai_ns_reconstruction/original_n1_inner_m15_numeric.py';CUTOFF='src/openai_ns_reconstruction/original_n1_cutoff_numeric.py';GAMMA='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_Gamma_future_functions.py';VOLTERRA='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_inner_volterra_functions.py'
CORE='experiments/root_st073/lei_ren_part1_paper_compliant_core_physical_field.json';CC='experiments/root_st073/lei_ren_part1_paper_compliant_core_physical_field_check.json';TRANSFER='experiments/root_st073/lei_ren_part1_paper_compliant_core_transfer.json';TC='experiments/root_st073/lei_ren_part1_paper_compliant_core_transfer_check.json'
AXIAL='experiments/root_st073/n1_first_word_axial_m15_v1_receipt.json.gz';CUT_SOURCE='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_inner_cutoff_functions.py';CONNECTION_SOURCE='src/openai_ns_reconstruction/n1_first_collar_forcing_jets.py'
CONNECTION='src/openai_ns_reconstruction/n1_leading_F_connection_v1.py';CONNECTION_FIX='src/openai_ns_reconstruction/n1_leading_F_connection_fixture_v1.py'
LEAF='src/openai_ns_reconstruction/n1_first_word_swirl_m15_v1.py';FIXTURE='src/openai_ns_reconstruction/n1_first_word_swirl_m15_fixture_v1.py';RUNNER='scripts/run_n1_first_word_swirl_m15_v1_once.py'
REQUIRED=(JOIN,RADIAL,TAIL,INNER,CUTOFF,GAMMA,VOLTERRA,CORE,CC,TRANSFER,TC,AXIAL,CUT_SOURCE,CONNECTION_SOURCE,CONNECTION,CONNECTION_FIX,LEAF,FIXTURE,RUNNER)
FIXED={AXIAL:'02acdcb22d9b0f03b4201f502388814d54d49937124ebb2288e66b9706770ab5',CUT_SOURCE:'b9f8b854b51b9eb0eda338503c8d25e6a1b69307da54b75bf40c1bd212ff1b29',CONNECTION_SOURCE:'9f853464f5ba19291418981efedd850d4a9af2e933d339c1c832efd9f326fdbf',CORE:'ad000166f6fce1ea94a7beeccf3324536fbffca073487b873d5460133f380096',CC:'0814fdd1afc27cbd260ec713f834d05804efb7827d2d500cc3551aad8f94fb17',TRANSFER:'4d7ec9b08ee2d36970d22c2a3287b8619bb1fa3deebf5ea264ab0d44b870d497',TC:'91883a7c3d92b1c01ca0d66c5852a5d65814aa364dc7825adb8dd9323c5cb848',JOIN:'a2e25ccec22f97dc99d1eb80609f0d806d4276e7a75e763ca882d90532e7ce67',RADIAL:'27d0e2a97568c50ff62f373cbdf04dea262b4980392d7ac4b13b9c7db21752cf',TAIL:'0e04c87a011bfec74e0ceb0f822774a03495a08a839aaf4e378d21949736cf1e',INNER:'004341ab4eb1f4693995bc86a748eae4c20ff2b5406fb0d784c4416963f5988b',CUTOFF:'7a2004b378ee495c6c053d04a8394406937f3a6fca235a721138f1ce4661c731',GAMMA:'38d3a8e5abd964c42b135410485f6ad58dfe9f7656035b4b5c7b1ed48c2ebb01',VOLTERRA:'ff82f4edb03922a551cba3e2df20217a55dc80ed42bd39ccb8b77b83f63b0208'}
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
    join=json.loads(gzip.decompress(raw[JOIN]));radial=json.loads(gzip.decompress(raw[RADIAL]));axial=json.loads(gzip.decompress(raw[AXIAL]));selected=json.loads(raw[TAIL]);core=json.loads(raw[CORE]);cc=json.loads(raw[CC]);transfer=json.loads(raw[TRANSFER]);tc=json.loads(raw[TC])
    require(join['result']['status']=='returned' and radial['result']['status']=='returned' and axial['result']['status']=='returned','accepted receipts')
    a=join['result']['output'];r=radial['result']['output'];family=a['source_family'];atoms=a['source_logs']
    require(r['source_family']==selected['source_family']==family and r['source_logs']==atoms,'same family and four source logs')
    require(axial['result']['output']['source_family']==family and axial['result']['output']['source_logs']==atoms and axial['result']['output']['sc']==a['sc'] and axial['result']['output']['appended_to_cumulative_m15'] is False,'same previous axial component')
    require(cc['all_passed'] and tc['all_passed'] and cc['input_hashes'][Path(CORE).name]==hashes[CORE] and tc['input_hashes'][Path(TRANSFER).name]==hashes[TRANSFER] and core['input_hashes'][Path(TRANSFER).name]==hashes[TRANSFER] and all(core[k]==family[k] for k in ('actual_five_defect_family_sha256','implicit_source_sha256','datum_enclosure_sha256')) and transfer['implicit_source_sha256']==core['implicit_source_sha256'] and transfer['datum_enclosure_sha256']==core['datum_enclosure_sha256'],'same core and transfer source')
    require(r['source_receipts_and_conditional_common_core_checked'] and r['inherited_interval_errors_preserved_in_sector_origins'] and a['allword_tail_added'] is False and a['result']['local_added_once'] and a['result']['core_history_reset'] is False,'same actual field and finite word')
    measure=selected['actual_cutoff_complete_support_measure'];require(measure['passed'] and measure['source_family']==family and measure['measure_is_Z_independent'] and '0<hb<1/2' in measure['inequality_proof'],'geometry source')
    inner=raw[INNER].decode('utf8');cutoff=raw[CUTOFF].decode('utf8');gamma=raw[GAMMA].decode('utf8');volterra=raw[VOLTERRA].decode('utf8')
    require("density='2*rho*U_lead*u_cut'" in inner and "measure='hb*(sc/4)*dq'" in inner and "R='epsilon*rho'" in inner and "definition='1-sigma(4q-2)',only_f_and_u_multiplied=True" in cutoff and 'exp(-1/x^2)/(exp(-1/x^2)+exp(-1/(1-x)^2)) inside' in gamma and 'raw_physical_units_retained=True' in volterra,'original physical density, sigma and cutoff definitions')
    cut_source=raw[CUT_SOURCE].decode('utf8');connection_source=raw[CONNECTION_SOURCE].decode('utf8')
    require("pressure=a.add(a.scale(a.mul(f,F),2),system['N1p'])" in cut_source and 'a.neg(a.scale(pressure,g.mul(eps2,r)))' in cut_source and 'H = z * (one - delta) / 2 + d * (4 * z + j)' in connection_source and 'ell = -(Lambda * L * H) / denominator' in connection_source,'original swirl sign, scale and F0 connection')
    c=MPIntervalContext();c.prec=1800;sc=iv(c,a['sc']);bits=tuple(measure['source_sc_exact_mpf_tuple']);require(sc._mpi_==(bits,bits),'same exact sc')
    whole=core['whole_core'];lam=iv(c,core['actual_Lambda']);delta=iv(c,core['actual_delta']);j=iv(c,transfer['required_j']);z=iv(c,whole['Z'])
    require(z.a==-1 and z.b==1 and iv(c,transfer['Lambda'])._mpi_==lam._mpi_ and iv(c,atoms['Lambda'])._mpi_==iv(c,transfer['logLambda'])._mpi_ and iv(c,atoms['F0_anchor'])._mpi_==iv(c,whole['F0_exact_positive_log_enclosure'])._mpi_,'same exact F0 and Lambda sources')
    markers=(('Cauchy_phase','symmetric_certified_error_not_field_value'),('Cauchy_phase','cauchy_phase_variation_error'),('mean_phase','U_phase_error_history_integral'),('cutoff_value','minus_sigma_factor'),('cutoff_value','hb_sigma_factor'),('cutoff_derivative','sigma_prime_factor'))
    index={'F0_connection':{'receipt':CORE,'path':['actual_delta','whole_core','Z'],'transfer':TRANSFER,'formula':'F0_Z/F0=Lambda*(-L*H/(H^2+(j/500)^2))'}};phi={};word={}
    for order in (0,1):
        key='(0, '+str(order)+')';phi[order]=[]
        for i,t in enumerate(r['actual_radial_sector_enclosures']['Phi'][key]):
            sid=f'swirl_phi:{order}:{i}';origin_json=json.dumps(t['origin'],sort_keys=True,separators=(',',':'));classes=sorted({label for label,needle in markers if needle in origin_json})
            index[sid]={'receipt':RADIAL,'path':['actual_radial_sector_enclosures','Phi',key,i],'origin_sha256':sha(origin_json.encode()),'error_classes':classes}
            phi[order].append({'powers':t['powers'],'coefficient':iv(c,t['coefficient']),'origin':{'source_ids':[sid],'uncertainty':bool(classes),'error_classes':classes}})
        word[order]={}
        for branch,group in (('core','core_transport'),('local','local_contribution')):
            word[order][branch]=[]
            for i,t in enumerate(a['result'][group][str(order)]['F1']):
                sid=f'swirl_word:{order}:{branch}:{i}';origin=t['origin'];index[sid]={'receipt':JOIN,'path':['result',group,str(order),'F1',i],'parent_source_ids':origin['source_ids']}
                word[order][branch].append({'powers':t['powers'],'coefficient':iv(c,t['coefficient']),'origin':{'source_ids':[sid],'uncertainty':origin['uncertainty'],'error_classes':origin['error_classes']}})
    conn=module(CONNECTION,'n1_leading_F_connection');cf=module(CONNECTION_FIX,'n1_leading_F_connection_fixture');connection_proof=cf.run(conn,c);require(connection_proof['passed'],'connection fixture failed '+str([k for k,v in connection_proof['checks'].items() if not v]))
    connected=conn.build(c,z=z,delta=delta,j=j,phi=phi,connection_source_id='F0_connection');leading=connected['leading']
    leaf=module(LEAF,'n1_first_word_swirl_m15_v1');fixture=module(FIXTURE,'n1_first_word_swirl_m15_fixture_v1');proof=fixture.run(leaf,c);require(proof['passed'],'fixture failed '+str([k for k,v in proof['checks'].items() if not v]))
    cert={'hb_strictly_between_zero_and_half':True,'sc_phase_Z_independent':True,'hb_epsilon_Z_independent':True,'original_cutoff_mass_5_over_8':True}
    result=leaf.build(c,sc=sc,leading=leading,first_word=word,source_family=family,source_atoms=atoms,certificate=cert)
    seen=set();classes=set()
    for branches in result['swirl'].values():
        for terms in branches.values():
            for term in terms:
                o=term['origin'];require(all(k in index for k in o['source_ids']),'unmapped product source');seen.update(o['source_ids']);classes.update(o['error_classes'])
    require(seen==set(index) and classes,'missing provenance or error classes')
    return {'source_family':family,'source_logs':atoms,'sc':a['sc'],'domain':a['domain'],'fixture':proof,'F0_connection_fixture':connection_proof,'F0_connection':enc(connected),'input_index':index,'referenced_input_ids':len(seen),'output_error_classes':sorted(classes),'result':enc(result),
      'physical_epsilon_squared_times_raw_F1':True,'cutoff_integral_5_over_8_by_exact_symmetry':True,'source_z_derivatives_are_ordinary':True,
      'whole_phase_uniform_boxes_not_point_values':True,'first_word_W1_only':True,'appended_to_cumulative_m15':False,'allword_tail_added':False,'complete_m15_computed':False,'terminal_closure_proved':False,
      'prior_axial_receipt':AXIAL,'prior_Z2_receipt':a['prior_Z2_green_receipt'],'alternative_cauchy_receipt':a['alternative_cauchy_receipt'],'allword_tail_receipt':a['allword_tail_receipt']}
def main():
    h,s,b=preflight();entry=save(ENTRY,{'component':COMPONENT,'source_sha256':s,'brief_sha256':b,'input_hashes':h,'entered':True});out={'component':COMPONENT,'source_sha256':s,'entry_sha256':entry,'status':'failed'}
    try:out.update(status='returned',output=consume(h),error=None)
    except BaseException as exc:out['error']={'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()}
    digest=save(RESULT,out);print(json.dumps({'component':COMPONENT,'status':out['status'],'result_sha256':digest}));return 0 if out['status']=='returned' else 1
if __name__=='__main__':raise SystemExit(main())
