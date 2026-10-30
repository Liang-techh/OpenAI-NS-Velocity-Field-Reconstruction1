"""One-entry same-source receipt binding and new collar Laurent checks."""
from __future__ import annotations
import gzip,hashlib,importlib.util,json,sys,traceback
from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext
ROOT=Path(__file__).resolve().parents[1];BRIEF=ROOT.parent.parent/'RESUME_BRIEF_2026_10_16.json'
COMPONENT='n1_collar_laurent_binding_v3';DIR=ROOT/'tmp'/COMPONENT
SOURCE,ENTRY,RESULT=(DIR/n for n in ('source.json','entered.json','result.json'))
BASE='experiments/root_st073/lei_ren_part1_paper_compliant_'
JOIN=BASE+'current_original_lower_active_common_source_join_functions_check.json'
JOIN_DATA=BASE+'current_original_lower_active_common_source_join_functions.json.gz'
JOIN_CHECK_SOURCE=BASE+'current_original_lower_active_common_source_join_functions_check.py'
JOIN_SOURCE=BASE+'current_original_lower_active_common_source_join_functions.py'
INTERFACE=BASE+'current_core_first_interface_check.json'
CORE='experiments/root_st073/n1_core_word1_z2_v2_receipt.json.gz'
MEAN='src/openai_ns_reconstruction/n1_collar_own_mean_history_v3.py';MEAN_FIX='src/openai_ns_reconstruction/n1_collar_own_mean_fixture_v3.py'
RADIAL='src/openai_ns_reconstruction/n1_collar_laurent_radial_v3.py';RADIAL_FIX='src/openai_ns_reconstruction/n1_collar_laurent_radial_fixture_v3.py'
FORCE='src/openai_ns_reconstruction/n1_collar_laurent_forcing_v3.py';FORCE_FIX='src/openai_ns_reconstruction/n1_collar_laurent_forcing_fixture_v3.py'
RUNNER='scripts/run_n1_collar_laurent_binding_v3_once.py'
FAILURE='experiments/root_st073/n1_collar_laurent_binding_v1_failure.json.gz'
FAILURE2='experiments/root_st073/n1_collar_laurent_binding_v2_failure.json.gz'
REQUIRED=(JOIN,JOIN_DATA,JOIN_CHECK_SOURCE,JOIN_SOURCE,INTERFACE,CORE,MEAN,MEAN_FIX,RADIAL,RADIAL_FIX,FORCE,FORCE_FIX,RUNNER,FAILURE,FAILURE2)
FIXED={JOIN:'980475686b9d884b3b03d641023910e9213853d48c889456c3883fb50f703e4b',
 JOIN_DATA:'97528e90f3a0f439bef768f8284ed10696a4350b8905fa69260f83b987817cb6',
 JOIN_CHECK_SOURCE:'c7ac551c511227b20128f192f92a3a4892a8c174a503483fb8b8ea21a368336b',
 JOIN_SOURCE:'fd8f42e0a36581ea49fca395d88d6b278614a867fcdcca641226e9dcdd12f2dc',
 INTERFACE:'e3a631d85da4cbeb62b9f6bbd7345e1109349935185edab177d2252a9d010400',
 CORE:'c53ff32ad3ce5288279953e84bc7417b71477c82b00d94466d0d2880adbc7ac9',
 FAILURE:'4a1df0b0b664d7f0b3b9ca5fde8fa867db91cf1d63475b907effc8247610e63e',
 FAILURE2:'057197a162b9697c110d1308e21b164441fd931ef7a4545add54178a78794dd8'}
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
    require(d['component']==COMPONENT and d['phase']==COMPONENT+'_pending' and d['goal_active'] is True
      and d['full_goal_complete'] is False and d['one_entry_only'] is True and d['do_not_reenter'] is True,'brief guard')
    require(Path(ck['path']).resolve()==SOURCE.resolve() and ck['sha256']==sha(src),'source checkpoint')
    manifest=json.loads(src);require(manifest['component']==COMPONENT and manifest['one_entry_only'] is True,'manifest')
    hashes={}
    for item in manifest['files']:
        p=item['path'];require(p in REQUIRED and p not in hashes,'unexpected source')
        h=sha((ROOT/p).read_bytes());require(h==item['sha256'] and (p not in FIXED or h==FIXED[p]),'source changed: '+p);hashes[p]=h
    require(set(hashes)==set(REQUIRED),'source missing');return hashes,sha(src),sha(b)
def module(path,name):
    spec=importlib.util.spec_from_file_location(name,ROOT/path);obj=importlib.util.module_from_spec(spec)
    sys.modules[name]=obj;spec.loader.exec_module(obj);return obj
def bind_saved_proof(hashes):
    for p in REQUIRED:require(sha((ROOT/p).read_bytes())==hashes[p],'post-entry input changed: '+p)
    join=json.loads((ROOT/JOIN).read_bytes());interface=json.loads((ROOT/INTERFACE).read_bytes());core=json.loads(gzip.decompress((ROOT/CORE).read_bytes()))
    require(join['all_passed'] is True and interface['all_passed'] is True and core['result']['status']=='returned','accepted receipts')
    for p in (JOIN_DATA,JOIN_CHECK_SOURCE,JOIN_SOURCE,INTERFACE):
        require(join['input_hashes'][Path(p).name]==hashes[p],'joined source hash '+p)
    family=join['source_family'];prior=core['result']['output']['source_family']
    require(family==prior,'joined and n1 source families differ')
    for key in ('actual_five_defect_family_sha256','implicit_source_sha256','datum_enclosure_sha256'):
        require(interface[key]==family[key],'core interface family '+key)
    c=join['checks'];require(all(c[k]['passed'] is True for k in ('source','inlet_primitives','physical_direction','comparison_core','common_IVP','physical_recovery','ordinary_sources','trace_transfer')),'proof checks')
    require(c['source']['actual_core_bridge_lower_and_finite_N_owner_links'] is True
       and c['inlet_primitives']['M_at_rho4_equals_U0_plus_epsilon_mean_Psi_by_integral_linearity'] is True,'same owner and inlet')
    d=c['physical_direction'];require(d['actual_current_direction_and_three_sector_drive_value_identities']==2
       and d['admitted_integrated_actual_core_equations_and_zero_axis_constants_consumed'] is True
       and d['actual_core_substitutions_give_D_plus_2_Phi_rho_over_epsilon_Phi_zero'] is True
       and d['actual_core_substitutions_give_full_drive_plus_epsilon_rho_Psi_rho_zero'] is True,'conditional direction identities')
    v=c['common_IVP'];require(v['angular_positive_linear_IVP_uniqueness_with_nonconstant_chi'] is True
       and v['actual_left_and_right_ratio_ODE_and_unit_inlet_identities']==4
       and v['actual_lower_and_prescribed_axial_integral_function_identity'] is True
       and v['mean_and_absolute_pressure_follow_same_M_C_FTC_IVP'] is True,'conditional same IVP')
    rows=set(c['ordinary_sources']['actual_saved_field_ordinary_derivative_rows'])
    require(all(side+'_'+field+'_Z'+str(k) in rows for side in ('left','right') for field in ('U','m') for k in range(4)), 'ordinary Z0..3 source field rows')
    theorem=interface['original_integrated_stress_free_core_equations'];atoms=interface['scaled_atoms_are_the_original_physical_primitives']
    require(theorem['passed'] is True and bool(theorem['axis_integration_constant_proof']) and atoms['passed'] is True,'core integrated equations and physical primitives')
    return {'source_family':family,'same_core_family_as_n1_receipt':True,
      'conditional_on_admitted_core_equations_and_checked_source_owners':True,
      'conditional_angular_direction':'D_over_R=-2*Lambda*Phi_rho/Phi',
      'conditional_axial_direction':'Gtotal=-R*d_R(U_core)',
      'conditional_bridge_ratio_identity':'same positive ODE and unit inlet, after admitted core substitutions',
      'conditional_mean_identity':'M and C histories have same checked ODE and inlet',
      'ordinary_field_Z0_to_Z3_rows_present':True,'all_six_histories_independently_reproved':False,
      'actual_collar_numerical_values_computed':False,'terminal_endpoint_cancellation_proved':False}
def main():
    hashes,source,b=preflight();entry=save(ENTRY,{'component':COMPONENT,'source_sha256':source,'brief_sha256':b,'input_hashes':hashes,'entered':True})
    out={'component':COMPONENT,'source_sha256':source,'entry_sha256':entry,'status':'failed'}
    try:
        binding=bind_saved_proof(hashes)
        out['partial']={'saved_proof_binding':binding,'fixtures':{}}
        c=MPIntervalContext();c.prec=300
        mean=module(MEAN,'n1_own_mean');mf=module(MEAN_FIX,'n1_own_mean_fixture')
        radial=module(RADIAL,'n1_collar_radial');rf=module(RADIAL_FIX,'n1_collar_radial_fixture')
        force=module(FORCE,'n1_collar_laurent_force');ff=module(FORCE_FIX,'n1_collar_laurent_force_fixture')
        proofs={}
        for label,fixture,leaf in (('own_mean',mf,mean),('radial',rf,radial),('forcing',ff,force)):
            proof=fixture.run(leaf,c);require(proof['passed'] is True,label+' fixture failed: '+str([k for k,v in proof['checks'].items() if not v]));proofs[label]=proof;out['partial']['fixtures'][label]=proof
        out.update(status='returned',output={'saved_proof_binding':binding,'fixtures':proofs,
          'actual_collar_source_values_computed':False,'actual_source_jet_provider_constructed':False,
          'allword_tail_added':False,'m15_computed':False},error=None)
    except BaseException as exc:out['error']={'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()}
    digest=save(RESULT,out);print(json.dumps({'component':COMPONENT,'status':out['status'],'result_sha256':digest}))
    return 0 if out['status']=='returned' else 1
if __name__=='__main__':raise SystemExit(main())
