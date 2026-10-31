"""One-entry physical collar m15 from new Q/P homogeneous inlet sectors."""
from __future__ import annotations
import gzip,hashlib,importlib.util,json,sys,traceback
from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext
ROOT=Path(__file__).resolve().parents[1];BRIEF=ROOT.parent.parent/'RESUME_BRIEF_2026_10_16.json';COMPONENT='n1_collar_qp_j2_m15_cells_v6';DIR=ROOT/'tmp'/COMPONENT
SOURCE,ENTRY,RESULT=(DIR/n for n in ('source.json','entered.json','result.json'))
GAUGE='experiments/root_st073/n1_core_gauge_exit_v2_receipt.json.gz';QP='experiments/root_st073/n1_collar_qp_inlet_v2_receipt.json.gz';PRESS='experiments/root_st073/n1_first_collar_g4_pressure_v1_receipt.json.gz';RADIAL='experiments/root_st073/n1_actual_collar_radial_v1_receipt.json.gz';SWIRL='experiments/root_st073/n1_first_word_swirl_m15_v1_receipt.json.gz';CORRECTION='experiments/root_st073/n1_first_word_swirl_m15_v1_provenance_correction.json';SELECT='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_selected_words.json'
CUT_SOURCE='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_inner_cutoff_functions.py';GAMMA='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_Gamma_future_functions.py';CUTOFF='src/openai_ns_reconstruction/original_n1_cutoff_numeric.py'
LEAF='src/openai_ns_reconstruction/n1_collar_qp_j2_m15_cells_v6.py';FIXTURE='src/openai_ns_reconstruction/n1_collar_qp_j2_m15_cells_fixture_v6.py';RUNNER='scripts/run_n1_collar_qp_j2_m15_cells_v6_once.py'
KCHI='experiments/root_st073/n1_core_kchi_first_v1_receipt.json.gz';PREVIOUS='experiments/root_st073/n1_collar_qp_k0_m15_cells_v1_receipt.json.gz'
J1='experiments/root_st073/n1_collar_qp_j1_m15_cells_v1_receipt.json.gz'
V1_ARCHIVE='experiments/root_st073/n1_collar_qp_j2_m15_cells_v1_terminated.json.gz';VALIDATION='src/openai_ns_reconstruction/n1_collar_qp_j2_independent_validation_v1.py'
V2_ARCHIVE='experiments/root_st073/n1_collar_qp_j2_m15_cells_v2_failure.json.gz'
V3_ARCHIVE='experiments/root_st073/n1_collar_qp_j2_m15_cells_v3_failure.json.gz'
V4_ARCHIVE='experiments/root_st073/n1_collar_qp_j2_m15_cells_v4_failure.json.gz'
V5_ARCHIVE='experiments/root_st073/n1_collar_qp_j2_m15_cells_v5_failure.json.gz'
REQUIRED=(GAUGE,QP,PRESS,RADIAL,SWIRL,CORRECTION,SELECT,CUT_SOURCE,GAMMA,CUTOFF,KCHI,PREVIOUS,J1,V1_ARCHIVE,V2_ARCHIVE,V3_ARCHIVE,V4_ARCHIVE,V5_ARCHIVE,LEAF,FIXTURE,VALIDATION,RUNNER)
FIXED={V5_ARCHIVE:'b7b67d481bcfc0a478bb9819ee0ea3078b813e82a9d14cb74e5a18a34a488a38',V4_ARCHIVE:'66718963ac536ca9cb8a98301cb1021112fa6b8aa9b8554fce9cfb04555a7edd',V3_ARCHIVE:'9879e424df9cd9b6a4f0902cbacd7a5daa6e6875c4b33f36ec3ccdd4706bdc9f',V2_ARCHIVE:'7eafe44fd46d34d20c320359aa51d9f28d7d5304aa7fd5ef4bb00be353ba39c1',V1_ARCHIVE:'2592ea3a34e819d928d26551d467c9ea95a0ac43a038a37de0b92040c15fd5cf',J1:'d03bb6c2ff6b0d26a74637ba9fb3428f55a1192498972ad3b65c823b6fd18144',KCHI:'b0d654203e8c4acabdd0b61ca7b79c1b7050314b0c003a0c2f8757bde94973a0',PREVIOUS:'48409ed314e95fe0f3117f6c27ed74b5aa8b56822c810ba5050d2b000b08e5c2',GAUGE:'4c1e49a83dbc7c13e5d9f09f4601ce34b26b487405a2d14a2f066ee2dfe27921',QP:'ac35695b793c6839b6f8cd030b68c413beda76147ab2700652aaf94b005574a7',PRESS:'c229cf188d306413f6edd0641a729321ad4567c2587ebf4701ada47720853660',RADIAL:'27d0e2a97568c50ff62f373cbdf04dea262b4980392d7ac4b13b9c7db21752cf',SWIRL:'a626e71fa2b0ddf40b79ce992016da2af8ac9955db3f17337029f3c86c2513b4',CORRECTION:'f2133bf5824b2f9636b430a33a2ba1d19fbf391a7f90a81de0c8269d1fe638c4',SELECT:'0e04c87a011bfec74e0ceb0f822774a03495a08a839aaf4e378d21949736cf1e',CUT_SOURCE:'b9f8b854b51b9eb0eda338503c8d25e6a1b69307da54b75bf40c1bd212ff1b29',GAMMA:'38d3a8e5abd964c42b135410485f6ad58dfe9f7656035b4b5c7b1ed48c2ebb01',CUTOFF:'7a2004b378ee495c6c053d04a8394406937f3a6fca235a721138f1ce4661c731'}
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
    require(Path(ck['path']).resolve()==SOURCE.resolve() and ck['sha256']==sha(src),'checkpoint');m=json.loads(src);require(m['component']==COMPONENT and m['one_entry_only'] is True,'manifest');h={}
    for row in m['files']:
        p=row['path'];require(p in REQUIRED and p not in h,'unexpected input');v=sha((ROOT/p).read_bytes());require(v==row['sha256'] and (p not in FIXED or v==FIXED[p]),'changed '+p);h[p]=v
    require(set(h)==set(REQUIRED),'missing input');return h,sha(src),sha(b)
def iv(c,r):
    lo=tuple(r['lower_exact_mpf_tuple']);hi=tuple(r['upper_exact_mpf_tuple']);require(len(lo)==len(hi)==4,'interval endpoints');x=c.make_mpf((lo,hi));require(x.a<=x.b,'reversed interval');return x
def conv(c,x):
    if isinstance(x,dict):
        if 'lower_exact_mpf_tuple' in x and 'upper_exact_mpf_tuple' in x:return iv(c,x)
        return {k:conv(c,v) for k,v in x.items()}
    if isinstance(x,list):return [conv(c,v) for v in x]
    return x
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
    for path in REQUIRED:
        v=(ROOT/path).read_bytes();require(sha(v)==hashes[path],'changed after entry '+path);raw[path]=v
    c=MPIntervalContext();c.prec=1800;leaf=module(LEAF,'n1_collar_qp_j2_m15_cells_v6');fixture=module(FIXTURE,'n1_collar_qp_j2_m15_cells_fixture_v6');proof=fixture.run_fixture(leaf,c);require(proof['passed'],'fixture failed '+json.dumps(proof,sort_keys=True,separators=(',',':')));validation=module(VALIDATION,'n1_collar_qp_j2_independent_validation_v1');independent=validation.run_validation();require(independent['passed'],'independent validation failed '+json.dumps(independent,sort_keys=True,separators=(',',':')))
    receipts={path:json.loads(gzip.decompress(raw[path])) for path in (GAUGE,QP,PRESS,RADIAL,SWIRL,KCHI,PREVIOUS,J1)};selected=json.loads(raw[SELECT]);correction=json.loads(raw[CORRECTION]);require(all(x['result']['status']=='returned' for x in receipts.values()),'accepted receipt status')
    g,q,p,r,f,k,prev,j1=(receipts[path]['result']['output'] for path in (GAUGE,QP,PRESS,RADIAL,SWIRL,KCHI,PREVIOUS,J1));family=g['source_family'];atoms=g['source_logs'];require(all(x['source_family']==family and x['source_logs']==atoms for x in (q,p,r,f,k,prev,j1)) and selected['source_family']==family,'same family and atoms')
    qp={v['path']:v['sha256'] for v in receipts[QP]['source']['files']};require(qp[GAUGE]==hashes[GAUGE] and qp[PRESS]==hashes[PRESS] and qp[SELECT]==hashes[SELECT] and qp[KCHI]==hashes[KCHI] and qp['experiments/root_st073/n1_core_kchi_compensation_v2_receipt.json.gz']=='c862f0560e8cb4b053abbbfcc982e4b4c1b602a806628e2d709399878e1c3aaf' and qp['experiments/root_st073/n1_collar_inlet_j2_tail_v3_receipt.json.gz']=='dd28d5be7e090afc720f7b71d48ddbeb0e542e8e4d5097e86a86d78ab8c36891' and qp['experiments/root_st073/n1_collar_axis_potential_binding_v1_receipt.json.gz']=='5859d5cb4bb82d372d8ac2f3ea1629e07c9787ffbaca149c99d1aba5edb8c6e8','same QP source chain and chi owner');require(correction['receipt_sha256']==hashes[SWIRL] and correction['numerical_result_changed'] is False and correction['input_id']=='F0_connection','F0 source path correction')
    contract=selected['actual_cutoff_complete_support_measure'];require(contract['passed'] and contract['source_family']==family and contract['measure_is_Z_independent'] and '0<hb<1/2' in contract['inequality_proof'],'physical geometry');sc=iv(c,p['sc']);bits=tuple(contract['source_sc_exact_mpf_tuple']);require(sc._mpi_==(bits,bits) and sc.a>0 and f['sc']==p['sc'],'exact sc')
    src=raw[CUT_SOURCE].decode('utf8');gamma=raw[GAMMA].decode('utf8');cutoff=raw[CUTOFF].decode('utf8');require("pressure=a.add(a.scale(a.mul(f,F),2),system['N1p'])" in src and 'a.neg(a.scale(pressure,g.mul(eps2,r)))' in src and "definition='1-sigma(4q-2)',only_f_and_u_multiplied=True" in cutoff and 'exp(-1/x^2)/(exp(-1/x^2)+exp(-1/(1-x)^2)) inside' in gamma,'original density and cutoff')
    fields=('F1','Uz1','K1','P1','d_x_F1','d_x_Uz1');groups={};index={};input_counts={}
    def ids(value):
        if value is None:return []
        if isinstance(value,str):return [value]
        require(isinstance(value,(list,tuple)),'source identifiers');return list(value)
    def add_group(name,data,path,*,error=False,orders='Z',support=fields):
        require(name not in groups,'duplicate group');packet={field:{} for field in fields};input_counts[name]={}
        for order in (0,1):
            key=('Z' if orders=='Z' else '')+str(order);row=data[key];require(set(row).issubset(fields),'unexpected field')
            for field in fields:
                vals=row.get(field,[]);require(isinstance(vals,list) and (field in support or not vals),'source field');out=[]
                for i,item in enumerate(vals):
                    require(isinstance(item,dict) and 'coefficient' in item and 'powers' in item,'source term');source=conv(c,item);origin=source.get('origin');om=origin if isinstance(origin,dict) else {};sid=f'{name}:{key}:{field}:{i}';parents=ids(source.get('source_ids'))+ids(source.get('source_id'))+ids(om.get('source_ids'))+ids(om.get('source_id'));classes=ids(source.get('error_classes'))+ids(om.get('error_classes'));require(len(source['powers'])==4 and all(isinstance(v,int) and not isinstance(v,bool) for v in source['powers']),'formal powers')
                    if error:
                        interval=item['coefficient'];lo=interval['lower_exact_mpf_tuple'];hi=interval['upper_exact_mpf_tuple'];require(source.get('error_only') is True and len(lo)==len(hi)==4 and lo[1:]==hi[1:] and ((lo[0],hi[0])==(1,0) or (lo[1]==0 and hi[1]==0)),'symmetric error source')
                    record={'receipt':GAUGE,'path':path+[key,field,i],'parent_source_ids':parents,'parent_origin':item.get('origin'),'parent_origin_record':item.get('origin_record')};require(q['source_index'].get(sid)==record,'QP exact source path and parent correspondence');index[sid]=record;meta={'source_ids':[sid],'error_classes':classes,'uncertainty':True,'parent_group':name,'error_only':error,'symmetric_error':error};out.append({'powers':source['powers'],'coefficient':source['coefficient'],'origin':meta,'error_only':error,'symmetric_error':error})
                packet[field]['Z'+str(order)]=out;input_counts[name][field+':Z'+str(order)]=len(out)
        groups[name]=packet
    for name in ('E01','E10','E11'):
        data=g['raw_Q_fixed_exit'][name];allowed=fields if name=='E01' else ('F1','d_x_F1');require(set(data['structural_zero_fields'])==set(fields)-set(allowed),'Q zero support');add_group('Q_'+name,data['aggregate'],['result','output','raw_Q_fixed_exit',name,'aggregate'],support=allowed)
    for name in ('G0_P_W0','G0_P_DeltaW0','D_P_W0','D_P_DeltaW0'):
        add_group('P_'+name,g['saved_potential_compensation_fixed_exit'][name]['exit_aggregate'],['result','output','saved_potential_compensation_fixed_exit',name,'exit_aggregate'],support=('F1','d_x_F1'))
    for name in ('G0','D'):
        add_group('Q_inner_R_'+name,g['inner_W0_kernel_tail_fixed_exit'][name]['aggregate'],['result','output','inner_W0_kernel_tail_fixed_exit',name,'aggregate'],error=True,support=fields if name=='G0' else ('F1','d_x_F1'))
    for name in ('R_Q_W0','R_Q_Delta','R_Q_Rg'):
        add_group('Q_outer_'+name,g['outer_kernel_tail_fixed_exit'][name],['result','output','outer_kernel_tail_fixed_exit',name],error=True,orders='',support=('F1','d_x_F1'))
    for name in ('G0_P_Rg','D_P_Rg','R_P_W0','R_P_Delta','R_P_Rg'):
        add_group('P_tail_'+name,g['potential_compensation_kernel_tail_fixed_exit'][name],['result','output','potential_compensation_kernel_tail_fixed_exit',name],error=True,orders='',support=('F1','d_x_F1'))
    require(input_counts==q['input_counts'],'same 17 source group ownership as saved QP endpoint consumer')
    markers=(('Cauchy_phase','symmetric_certified_error_not_field_value'),('Cauchy_phase','cauchy_phase_variation_error'),('mean_phase','U_phase_error_history_integral'),('cutoff_value','minus_sigma_factor'),('cutoff_value','hb_sigma_factor'),('cutoff_derivative','sigma_prime_factor'))
    leading_F={}
    for order in (0,1):
        leading_F[order]=[]
        for i,t in enumerate(f['F0_connection']['leading'][str(order)]):
            sid=f'lead_F:{order}:{i}';index[sid]={'receipt':SWIRL,'path':['result','output','F0_connection','leading',str(order),i],'parent_source_ids':t['origin']['source_ids'],'connection_path_correction':CORRECTION};leading_F[order].append({'powers':t['powers'],'coefficient':iv(c,t['coefficient']),'origin':{'source_ids':[sid],'uncertainty':t['origin']['uncertainty'],'error_classes':t['origin']['error_classes']}})
    chi={}
    for order in (0,1):
        sid=f'chi:Z{order}';index[sid]={'receipt':KCHI,'path':['result','output','chi','chi',order]};chi['Z'+str(order)]=[{'powers':[0,0,0,0],'coefficient':iv(c,k['chi']['chi'][order]),'origin':{'source_ids':[sid],'error_classes':['whole_Z_axis_potential_interval'],'uncertainty':True}}]
    previous_pins={v['path']:v['sha256'] for v in receipts[PREVIOUS]['source']['files']};require(all(previous_pins[path]==hashes[path] for path in (GAUGE,QP,PRESS,RADIAL,SWIRL,CORRECTION,SELECT,CUT_SOURCE,GAMMA,CUTOFF)) and prev['input_counts']==input_counts,'same K0 physical input chain')
    j1_pins={v['path']:v['sha256'] for v in receipts[J1]['source']['files']};require(all(j1_pins[path]==hashes[path] for path in (GAUGE,QP,PRESS,RADIAL,SWIRL,CORRECTION,SELECT,CUT_SOURCE,GAMMA,CUTOFF,KCHI,PREVIOUS)) and j1['input_counts']==input_counts,'same first-potential source chain')
    cert={'hb_strictly_between_zero_and_half':True,'sc_phase_Z_independent':True,'hb_epsilon_Z_independent':True,'original_cutoff_mass_5_over_8':True,'original_cutoff_mass_by_cell':True,'strict_hb_admitted':True,'hb_sc_Z_independent':True,'rho_radius_Z_independent':True}
    reduced_groups={name:{'Fc':packet['F1'],'Fxc':packet['d_x_F1']} for name,packet in groups.items()}
    result=leaf.build(c,sc=sc,inlet_groups=reduced_groups,leading_F=leading_F,chi=chi,source_family=family,source_atoms=atoms,certificate=cert)
    # Keep the disjoint branch terms and group metadata once; the full source
    # index retains every inlet term, including fields unused by these channels.
    compact={k:v for k,v in result.items() if k not in ('groups','group_results')}
    compact['group_results']={name:{k:v for k,v in value.items() if k not in ('inlet','axial','swirl','axial_paths','swirl_paths','axial_geometry','swirl_geometry')} for name,value in result['group_results'].items()}
    seen=set()
    for group in compact['group_results'].values():
        for channel in ('swirl',):
            for cells in group['paths'][channel].values():
                for branches in cells.values():
                    for terms in branches.values():
                        for term in terms:
                            source_ids=term['origin']['source_ids'];require(all(sid in index for sid in source_ids),'unindexed source id');seen.update(source_ids)
    return {'source_family':family,'source_logs':atoms,'source_index':index,'input_counts':input_counts,'fixture':proof,'independent_synthetic_validation':independent,'actual_qp_j2_m15_cells':enc(compact),'used_source_id_count':len(seen),'previous_qp_endpoint_reference':QP,'previous_j1_cells_reference':J1,'previous_k0_cells_reference':PREVIOUS,'previous_finite_prefix_reference':PRESS,'QP_source_paths_and_parents_match':True,'gauge_inputs_are_fixed_exit_inlet_boxes':True,'old_prefix_disjointness_not_proved':True,'E00_not_recomputed':True,'prior_endpoint_first_potential_not_recomputed':True,'complete_collar_field':False,'feedback_word_tail_included':False,'complete_m15_computed':False,'appended_to_ledger':False,'terminal_closure_proved':False}
def main():
    h,s,b=preflight();entry=save(ENTRY,{'component':COMPONENT,'source_sha256':s,'brief_sha256':b,'input_hashes':h,'entered':True});out={'component':COMPONENT,'source_sha256':s,'entry_sha256':entry,'status':'failed'}
    try:out.update(status='returned',output=consume(h),error=None)
    except BaseException as exc:out['error']={'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()}
    digest=save(RESULT,out);print(json.dumps({'component':COMPONENT,'status':out['status'],'result_sha256':digest}));return 0 if out['status']=='returned' else 1
if __name__=='__main__':raise SystemExit(main())
