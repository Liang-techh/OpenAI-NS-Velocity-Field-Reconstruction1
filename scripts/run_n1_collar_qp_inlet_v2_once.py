"""One-entry no-reset collar transport of saved Q/P inlet sectors."""
from __future__ import annotations
import gzip,hashlib,importlib.util,json,sys,traceback
from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext
ROOT=Path(__file__).resolve().parents[1];BRIEF=ROOT.parent.parent/'RESUME_BRIEF_2026_10_16.json';COMPONENT='n1_collar_qp_inlet_v2';DIR=ROOT/'tmp'/COMPONENT
SOURCE,ENTRY,RESULT=(DIR/n for n in ('source.json','entered.json','result.json'))
GAUGE='experiments/root_st073/n1_core_gauge_exit_v2_receipt.json.gz';COMP='experiments/root_st073/n1_core_kchi_compensation_v2_receipt.json.gz';KCHI='experiments/root_st073/n1_core_kchi_first_v1_receipt.json.gz';TAIL='experiments/root_st073/n1_collar_inlet_j2_tail_v3_receipt.json.gz';BIND='experiments/root_st073/n1_collar_axis_potential_binding_v1_receipt.json.gz';PRESS='experiments/root_st073/n1_first_collar_g4_pressure_v1_receipt.json.gz';SELECT='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_selected_words.json'
LEAF='src/openai_ns_reconstruction/n1_collar_qp_transport_v2.py';FIXTURE='src/openai_ns_reconstruction/n1_collar_qp_transport_fixture_v2.py';RUNNER='scripts/run_n1_collar_qp_inlet_v2_once.py'
FAILURE='experiments/root_st073/n1_collar_qp_inlet_v1_failure.json.gz'
REQUIRED=(GAUGE,COMP,KCHI,TAIL,BIND,PRESS,SELECT,LEAF,FIXTURE,RUNNER,FAILURE)
FIXED={GAUGE:'4c1e49a83dbc7c13e5d9f09f4601ce34b26b487405a2d14a2f066ee2dfe27921',COMP:'c862f0560e8cb4b053abbbfcc982e4b4c1b602a806628e2d709399878e1c3aaf',KCHI:'b0d654203e8c4acabdd0b61ca7b79c1b7050314b0c003a0c2f8757bde94973a0',TAIL:'dd28d5be7e090afc720f7b71d48ddbeb0e542e8e4d5097e86a86d78ab8c36891',BIND:'5859d5cb4bb82d372d8ac2f3ea1629e07c9787ffbaca149c99d1aba5edb8c6e8',PRESS:'c229cf188d306413f6edd0641a729321ad4567c2587ebf4701ada47720853660',SELECT:'0e04c87a011bfec74e0ceb0f822774a03495a08a839aaf4e378d21949736cf1e',FAILURE:'eeaff792752337e8792fbfb4608daa82daa37ed8d7493c9dae1964b3527e745a'}
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
    c=MPIntervalContext();c.prec=1800;leaf=module(LEAF,'n1_collar_qp_transport_v2');fixture=module(FIXTURE,'n1_collar_qp_transport_fixture_v2');proof=fixture.run_fixture(leaf,c);require(proof['passed'],'fixture failed '+str([k for k,v in proof['checks'].items() if not v]))
    receipts={path:json.loads(gzip.decompress(raw[path])) for path in (GAUGE,COMP,KCHI,TAIL,BIND,PRESS)};selected=json.loads(raw[SELECT]);require(all(x['result']['status']=='returned' for x in receipts.values()),'accepted receipt status')
    g,p,k,t,bind,press=(receipts[path]['result']['output'] for path in (GAUGE,COMP,KCHI,TAIL,BIND,PRESS));family=g['source_family'];atoms=g['source_logs'];require(all(x['source_family']==family and x['source_logs']==atoms for x in (p,k,t,bind,press)) and selected['source_family']==family,'same source family and formal atoms')
    gp={v['path']:v['sha256'] for v in receipts[GAUGE]['source']['files']};tp={v['path']:v['sha256'] for v in receipts[TAIL]['source']['files']};require(gp[COMP]==hashes[COMP] and gp[KCHI]==hashes[KCHI] and all(tp[path]==hashes[path] for path in (COMP,KCHI,BIND,PRESS,SELECT)),'authenticated source chain')
    contract=selected['actual_cutoff_complete_support_measure'];require(contract['passed'] and contract['source_family']==family and contract['measure_is_Z_independent'] and '0<hb<1/2' in contract['inequality_proof'],'certified geometry');sc=iv(c,press['sc']);bits=tuple(contract['source_sc_exact_mpf_tuple']);require(sc._mpi_==(bits,bits) and sc.a>0,'same exact sc')
    require(bind['source_checks']['one_query_jet_shared_in_integration'] is True and t['partial_inlet_groups_only'] is True and g['feedback_word_tail_included'] is False,'prior source scope')
    fields=('F1','Uz1','K1','P1','d_x_F1','d_x_Uz1');nonspin=('Uz1','K1','P1','d_x_Uz1');groups={};index={};input_counts={}
    def ids(value):
        if value is None:return []
        if isinstance(value,str):return [value]
        require(isinstance(value,(list,tuple)),'source identifiers');return list(value)
    def add_group(name,data,path,*,error=False,orders='Z',support=fields):
        require(name not in groups,'duplicate group '+name);packet={field:{} for field in fields};input_counts[name]={}
        for order in (0,1):
            key=('Z' if orders=='Z' else '')+str(order);row=data[key];require(set(row).issubset(fields),'unexpected field '+name)
            for field in fields:
                vals=row.get(field,[]);require(isinstance(vals,list),'term array');require(field in support or not vals,'structural-zero field '+name+field);out=[]
                for i,item in enumerate(vals):
                    require(isinstance(item,dict) and 'coefficient' in item and 'powers' in item,'source term contract');source=conv(c,item);origin=source.get('origin');om=origin if isinstance(origin,dict) else {};sid=f'{name}:{key}:{field}:{i}';parents=ids(source.get('source_ids'))+ids(source.get('source_id'))+ids(om.get('source_ids'))+ids(om.get('source_id'));classes=ids(source.get('error_classes'))+ids(om.get('error_classes'));require(len(source['powers'])==4 and all(isinstance(v,int) and not isinstance(v,bool) for v in source['powers']),'formal powers')
                    if error:
                        interval=item['coefficient'];lo=interval['lower_exact_mpf_tuple'];hi=interval['upper_exact_mpf_tuple'];require(source.get('error_only') is True and len(lo)==len(hi)==4 and lo[1:]==hi[1:] and ((lo[0],hi[0])==(1,0) or (lo[1]==0 and hi[1]==0)),'symmetric error source')
                    index[sid]={'receipt':GAUGE,'path':path+[key,field,i],'parent_source_ids':parents,'parent_origin':item.get('origin'),'parent_origin_record':item.get('origin_record')};meta={'source_ids':parents+[sid],'error_classes':classes,'uncertainty':True,'parent_group':name,'error_only':error,'symmetric_error':error};term={'powers':source['powers'],'coefficient':source['coefficient'],'origin':meta,'error_only':error,'symmetric_error':error};out.append(term)
                packet[field]['Z'+str(order)]=out;input_counts[name][field+':Z'+str(order)]=len(out)
        groups[name]=packet
    for name in ('E01','E10','E11'):
        data=g['raw_Q_fixed_exit'][name];allowed=fields if name=='E01' else ('F1','d_x_F1');require(set(data['structural_zero_fields'])==set(fields)-set(allowed),'declared Q zero support');add_group('Q_'+name,data['aggregate'],['result','output','raw_Q_fixed_exit',name,'aggregate'],support=allowed)
    for name in ('G0_P_W0','G0_P_DeltaW0','D_P_W0','D_P_DeltaW0'):
        data=g['saved_potential_compensation_fixed_exit'][name];add_group('P_'+name,data['exit_aggregate'],['result','output','saved_potential_compensation_fixed_exit',name,'exit_aggregate'],support=('F1','d_x_F1'))
    for name in ('G0','D'):
        data=g['inner_W0_kernel_tail_fixed_exit'][name];allowed=fields if name=='G0' else ('F1','d_x_F1');add_group('Q_inner_R_'+name,data['aggregate'],['result','output','inner_W0_kernel_tail_fixed_exit',name,'aggregate'],error=True,support=allowed)
    for name in ('R_Q_W0','R_Q_Delta','R_Q_Rg'):
        add_group('Q_outer_'+name,g['outer_kernel_tail_fixed_exit'][name],['result','output','outer_kernel_tail_fixed_exit',name],error=True,orders='',support=('F1','d_x_F1'))
    for name in ('G0_P_Rg','D_P_Rg','R_P_W0','R_P_Delta','R_P_Rg'):
        add_group('P_tail_'+name,g['potential_compensation_kernel_tail_fixed_exit'][name],['result','output','potential_compensation_kernel_tail_fixed_exit',name],error=True,orders='',support=('F1','d_x_F1'))
    chi={}
    for order in (0,1):
        sid=f'chi:Z{order}';index[sid]={'receipt':KCHI,'path':['result','output','chi','chi',order]};chi['Z'+str(order)]=[{'powers':[0,0,0,0],'coefficient':iv(c,k['chi']['chi'][order]),'origin':{'source_ids':[sid],'error_classes':['whole_Z_axis_potential_interval'],'uncertainty':True}}]
    cert={'strict_hb_admitted':True,'hb_sc_Z_independent':True,'rho_radius_Z_independent':True,'source':'selected actual_cutoff_complete_support_measure','source_sha256':hashes[SELECT]};out=leaf.build(c,inlet_groups=groups,chi=chi,hb=c.mpf([0,'0.5']),sc=sc,source_family=family,source_atoms=['hb','F0_anchor','Pstar_squared','Lambda'],geometry_certificate=cert)
    require(set(out['group_results'])==set(groups),'transport group ownership')
    def compact_stage(stage):
        meta={k:v for k,v in stage.items() if k not in ('fields','response','branches')}
        meta['branches']={name:{k:v for k,v in packet.items() if k!='response'} for name,packet in stage['branches'].items()}
        return meta
    compact={'geometry':out['geometry'],'kernel_bounds':out['kernel_bounds'],'groups':{name:{**{k:v for k,v in value.items() if k not in ('inlet','k0','first_potential','j2_tail')},**{stage:compact_stage(value[stage]) for stage in ('k0','first_potential','j2_tail')}} for name,value in out['group_results'].items()}}
    return {'source_family':family,'source_logs':atoms,'source_index':index,'input_counts':input_counts,'fixture':proof,'transport':enc(compact),'previous_E00_reference':{'receipt':g['saved_baseline_E00_receipt'],'status':'previous_K0_collar_result_not_recomputed'},'previous_W0_groups_reference':{'receipt':TAIL,'path':['result','output','inlet_j2_tail']},'branch_owners_preserved':True,'formal_atoms_not_materialized':True,'E00_not_recomputed':True,'complete_collar_field':False,'feedback_word_tail_included':False,'m15_computed':False,'appended_to_ledger':False,'terminal_closure_proved':False}
def main():
    h,s,b=preflight();entry=save(ENTRY,{'component':COMPONENT,'source_sha256':s,'brief_sha256':b,'input_hashes':h,'entered':True});out={'component':COMPONENT,'source_sha256':s,'entry_sha256':entry,'status':'failed'}
    try:out.update(status='returned',output=consume(h),error=None)
    except BaseException as exc:out['error']={'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()}
    digest=save(RESULT,out);print(json.dumps({'component':COMPONENT,'status':out['status'],'result_sha256':digest}));return 0 if out['status']=='returned' else 1
if __name__=='__main__':raise SystemExit(main())
