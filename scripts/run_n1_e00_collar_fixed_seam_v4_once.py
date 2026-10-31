"""One-entry same-source fixed E00 inlet and no-reset collar transport."""
from __future__ import annotations
import gzip,hashlib,importlib.util,json,sys,traceback
from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext
ROOT=Path(__file__).resolve().parents[1];BRIEF=ROOT.parent.parent/'RESUME_BRIEF_2026_10_16.json';COMPONENT='n1_e00_collar_fixed_seam_v4';DIR=ROOT/'tmp'/COMPONENT
SOURCE,ENTRY,RESULT=(DIR/n for n in ('source.json','entered.json','result.json'))
E00='experiments/root_st073/n1_core_e00_fixed_missing_v1_receipt.json.gz';PROFILE='experiments/root_st073/n1_core_k0_feedback_m15_v2_receipt.json.gz';KCHI='experiments/root_st073/n1_core_kchi_first_v1_receipt.json.gz';BIND='experiments/root_st073/n1_collar_axis_potential_binding_v1_receipt.json.gz';PRESS='experiments/root_st073/n1_first_collar_g4_pressure_v1_receipt.json.gz';SELECT='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_selected_words.json';QP='experiments/root_st073/n1_collar_qp_j2_m15_cells_v6_receipt.json.gz'
LEAF='src/openai_ns_reconstruction/n1_e00_collar_fixed_seam_v4.py';FIXTURE='src/openai_ns_reconstruction/n1_e00_collar_fixed_seam_fixture_v4.py';RUNNER='scripts/run_n1_e00_collar_fixed_seam_v4_once.py'
FAILURE='experiments/root_st073/n1_e00_collar_fixed_seam_v1_failure.json.gz'
FAILURE2='experiments/root_st073/n1_e00_collar_fixed_seam_v2_failure.json.gz'
FAILURE3='experiments/root_st073/n1_e00_collar_fixed_seam_v3_failure.json.gz'
REQUIRED=(E00,PROFILE,KCHI,BIND,PRESS,SELECT,QP,FAILURE,FAILURE2,FAILURE3,LEAF,FIXTURE,RUNNER)
FIXED={FAILURE3:'e1880de571e505a9a4f92be016f5b84c605f1d0f499c6065d53730838a4f8af2',FAILURE2:'57cd598d34b227d368c8981d9791766ed0ffa19806b168f04cbef2bd549ec5d4',FAILURE:'44de3a7ce09e1aadd4f9b01c8cdbc589a7fd7ff351f04e0b87801c97cd6e19b0',E00:'39c8f92f0f2b9f3c284e163b053a4b428e1b09a1ae76677bae84fa4665f62dda',PROFILE:'ca947df9fdf7765d9211f428375d81956af783d89074f9e782ba56aeff91dad0',KCHI:'b0d654203e8c4acabdd0b61ca7b79c1b7050314b0c003a0c2f8757bde94973a0',BIND:'5859d5cb4bb82d372d8ac2f3ea1629e07c9787ffbaca149c99d1aba5edb8c6e8',PRESS:'c229cf188d306413f6edd0641a729321ad4567c2587ebf4701ada47720853660',SELECT:'0e04c87a011bfec74e0ceb0f822774a03495a08a839aaf4e378d21949736cf1e',QP:'4f894b94f8d2623f291d5ab3613c2d5fc9b98d93f95c0389a1ed22ce716a55c5'}
def sha(b):return hashlib.sha256(b).hexdigest()
def require(v,msg):
    if not v:raise ValueError(msg)
def save(path,value):
    b=(json.dumps(value,sort_keys=True,separators=(',',':'))+'\n').encode()
    with path.open('xb') as f:f.write(b)
    return sha(b)
def preflight():
    require(not ENTRY.exists() and not RESULT.exists(),'already entered');b=BRIEF.read_bytes();d=json.loads(b);src=SOURCE.read_bytes();ck=d['source_checkpoint'];require(d['component']==COMPONENT and d['phase']==COMPONENT+'_pending' and d['goal_active'] is True and d['full_goal_complete'] is False and d['one_entry_only'] is True and d['do_not_reenter'] is True,'brief guard');require(Path(ck['path']).resolve()==SOURCE.resolve() and ck['sha256']==sha(src),'checkpoint');m=json.loads(src);require(m['component']==COMPONENT and m['one_entry_only'] is True,'manifest');h={}
    for row in m['files']:
        p=row['path'];require(p in REQUIRED and p not in h,'unexpected input');v=sha((ROOT/p).read_bytes());require(v==row['sha256'] and (p not in FIXED or v==FIXED[p]),'changed '+p);h[p]=v
    require(set(h)==set(REQUIRED),'missing input');return h,sha(src),sha(b)
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
    for path in REQUIRED:
        v=(ROOT/path).read_bytes();require(sha(v)==hashes[path],'changed after entry '+path);raw[path]=v
    c=MPIntervalContext();c.prec=1800;leaf=module(LEAF,'n1_e00_collar_fixed_seam_v4');fixture=module(FIXTURE,'n1_e00_collar_fixed_seam_fixture_v4');proof=fixture.run_fixture(leaf,c);require(proof['passed'],'fixture failed '+json.dumps(proof,sort_keys=True,separators=(',',':')))
    receipts={path:json.loads(gzip.decompress(raw[path])) for path in (E00,PROFILE,KCHI,BIND,PRESS,QP)};selected=json.loads(raw[SELECT]);require(all(x['result']['status']=='returned' for x in receipts.values()),'accepted receipt status')
    e,p,k,b,press,q=(receipts[path]['result']['output'] for path in (E00,PROFILE,KCHI,BIND,PRESS,QP));family=e['source_family'];atoms=e['source_logs'];require(all(x['source_family']==family and x['source_logs']==atoms for x in (p,k,b,press,q)) and selected['source_family']==family,'common source family and formal atoms')
    pin=lambda r:{v['path']:v['sha256'] for v in r['source']['files']};require(pin(receipts[E00])[PROFILE]==hashes[PROFILE] and pin(receipts[E00])[QP]==hashes[QP] and pin(receipts[QP])[KCHI]==hashes[KCHI] and pin(receipts[BIND])[KCHI]==hashes[KCHI] and pin(receipts[QP])[PRESS]==hashes[PRESS] and pin(receipts[QP])[SELECT]==hashes[SELECT],'authenticated common source chain')
    require(b['source_checks']['one_query_jet_shared_in_integration'] is True and b['phase_shear_cutoff_identified_with_axis_potential'] is False,'axis chi owner');contract=selected['actual_cutoff_complete_support_measure'];require(contract['passed'] and contract['source_family']==family and contract['measure_is_Z_independent'] and '0<hb<1/2' in contract['inequality_proof'],'certified geometry');sc=iv(c,press['sc']);bits=tuple(contract['source_sc_exact_mpf_tuple']);require(sc._mpi_==(bits,bits) and sc.a>0,'same exact sc')
    fields=('F1','Uz1','K1','P1','d_x_F1','d_x_Uz1');inlet={field:{} for field in fields};index={};counts={}
    def ids(value):
        if value is None:return []
        if isinstance(value,str):return [value]
        require(isinstance(value,(list,tuple)),'source ids');return list(value)
    for order in (0,1):
        z='Z'+str(order)
        for field in fields:
            if field in ('F1','Uz1'):
                name='W1_Z'+str(order)+'_core_exit';items=p['profile'][name][field];path=['result','output','profile',name,field];receipt=PROFILE
            else:
                items=e['fixed_missing_fields']['fixed_fields'][z][field];path=['result','output','fixed_missing_fields','fixed_fields',z,field];receipt=E00
            expected=e['prior_fixed_F_U_reference'][z][field]['term_count'] if field in ('F1','Uz1') else None
            if expected is not None:require(expected==len(items),'same saved F/U field')
            require(len(items)>0,'nonzero inlet packet');counts[field+':'+z]=len(items);out=[]
            for i,item in enumerate(items):
                sid=f'E00:{z}:{field}:{i}';source=iv(c,item['coefficient']);f0=item['F0_power'];lam=item['Lambda_power'];require(isinstance(f0,int) and not isinstance(f0,bool) and isinstance(lam,int) and not isinstance(lam,bool),'formal powers');parent=ids(item.get('source_ids'))
                if field in ('F1','Uz1'):
                    cell=item['source_cell'];mn=item['matrix_name'];mo=str(item['matrix_order']);mk=item['matrix_key'];mi=item['matrix_term_index'];wo=str(item['W0_order']);wc=str(item['W0_column']);wi=item['W0_term_index'];parent=[f'matrix:{cell}:{mn}:{mo}:{mk}:{mi}',f'W0:{cell}:{wo}:{wc}:{wi}']
                require(len(parent)>=2 and all(pid in e['source_index'] for pid in parent),'indexed common E00 source paths')
                index[sid]={'receipt':receipt,'path':path+[i],'parent_source_ids':parent,'parent_origin':item.get('origin'),'source_cell':item.get('source_cell',item.get('outer_cell')),'matrix_name':item.get('matrix_name'),'matrix_key':item.get('matrix_key'),'matrix_order':item.get('matrix_order'),'matrix_term_index':item.get('matrix_term_index'),'W0_order':item.get('W0_order'),'W0_column':item.get('W0_column'),'W0_term_index':item.get('W0_term_index')}
                out.append({'source_ids':[sid],'powers':[0,f0,0,lam],'coefficient':source,'origin':{'source_ids':[sid],'uncertainty':True,'error_classes':['saved_core_E00_interval'],'parent_field':field,'parent_order':z}})
            inlet[field][z]=out
    chi={}
    for order in (0,1):
        sid=f'chi:Z{order}';index[sid]={'receipt':KCHI,'path':['result','output','chi','chi',order]};chi['Z'+str(order)]=[{'source_ids':[sid],'powers':[0,0,0,0],'coefficient':iv(c,k['chi']['chi'][order]),'origin':{'source_ids':[sid],'error_classes':['whole_Z_axis_potential_interval'],'uncertainty':True}}]
    cert={'strict_hb_admitted':True,'hb_sc_Z_independent':True,'rho_radius_Z_independent':True,'fixed_geometry':True,'axis_chi_shared':True,'source_family':family,'source':'authenticated selected geometry and axis potential receipts','owners':{'q':{'receipt':SELECT,'path':['actual_cutoff_complete_support_measure']},'x':{'receipt':E00,'path':['result','output','fixed_missing_fields','fixed_target']},'rho':{'receipt':SELECT,'path':['actual_cutoff_complete_support_measure']},'chi':{'receipt':BIND,'path':['result','output','source_checks']}}}
    result=leaf.build(c,sc=sc,hb=c.mpf([0,'0.5']),q_cells=[(0,'0.25'),('0.25','0.5'),('0.5',1)],inlet=inlet,chi=chi,source_family=family,source_atoms=['hb','F0_anchor','Pstar_squared','Lambda'],certificate=cert)
    source_by_id={term['origin']['source_ids'][0]:term for field in fields for z in ('Z0','Z1') for term in inlet[field][z]};seam=result['seam'];require(result['schema']=='n1_e00_collar_fixed_seam_v4' and seam['q'].a==seam['q'].b==0 and seam['x'].a==seam['x'].b==2 and seam['rho'].a==seam['rho'].b==4,'exact fixed seam geometry')
    matched=set();zero_routes=0
    for branch in seam['k0_transfer']['branches'].values():
        factor=branch['factor'];require(factor.a==factor.b and factor.a in (0,1),'seam factor is literal zero or one')
        for z in ('Z0','Z1'):
            for term in branch[z]:
                sids=term['source_ids'];require(len(sids)==1 and sids[0] in source_by_id,'seam source ownership');source=source_by_id[sids[0]]
                if factor.a==1:
                    require(term['coefficient']._mpi_==source['coefficient']._mpi_ and tuple(term['powers'])==tuple(source['powers']) and sids[0] not in matched,'identity route changes or duplicates inlet');matched.add(sids[0])
                else:require(term['coefficient'].a==term['coefficient'].b==0,'cross route not zero');zero_routes+=1
    require(matched==set(source_by_id),'missing seam identity source')
    for stage_name in ('first_potential','j2_tail'):
        for branch in seam[stage_name]['branches'].values():
            for z in ('Z0','Z1'):
                for term in branch[z]:require(term['coefficient'].a==term['coefficient'].b==0,'nonzero seam potential word')
    require(seam['local_forcing']['exact_zero_at_q0'] is True and result['j2_remainder_included'] is False,'seam source or tail scope')
    def stage(value):return {k:({name:{bk:bv for bk,bv in branch.items() if bk not in ('fields','response')} for name,branch in v.items()} if k=='branches' else v) for k,v in value.items() if k not in ('fields','response')}
    compact={k:v for k,v in result.items() if k not in ('inlet','chi_jets','source_term_ledger','q_cells','seam','q0_seam','k0_transfer','first_potential','j2_tail','local_forcing')}
    compact['seam']={k:(stage(v) if k in ('k0_transfer','first_potential','j2_tail','local_forcing') else v) for k,v in seam.items() if k not in ('source_term_ledger',)}
    compact['q_cells']=[{k:(stage(v) if k in ('k0_transfer','first_potential','j2_tail','local_forcing') else v) for k,v in cell.items()} for cell in result['q_cells']]
    return {'source_family':family,'source_logs':atoms,'source_index':index,'input_counts':counts,'fixture':proof,'seam_identity_source_count':len(matched),'seam_zero_cross_route_count':zero_routes,'actual_E00_seam_and_collar':enc(compact),'input_receipt_sha256':{path:hashes[path] for path in (E00,PROFILE,KCHI,BIND,PRESS,SELECT,QP)},'saved_E00_inlet_preserved_not_reset':True,'old_prefix_disjointness_not_proved':True,'complete_n1_inlet':False,'m15_computed':False,'appended_to_ledger':False,'terminal_closure_proved':False}
def main():
    h,s,b=preflight();entry=save(ENTRY,{'component':COMPONENT,'source_sha256':s,'brief_sha256':b,'input_hashes':h,'entered':True});out={'component':COMPONENT,'source_sha256':s,'entry_sha256':entry,'status':'failed'}
    try:out.update(status='returned',output=consume(h),error=None)
    except BaseException as exc:out['error']={'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()}
    digest=save(RESULT,out);print(json.dumps({'component':COMPONENT,'status':out['status'],'result_sha256':digest}));return 0 if out['status']=='returned' else 1
if __name__=='__main__':raise SystemExit(main())
