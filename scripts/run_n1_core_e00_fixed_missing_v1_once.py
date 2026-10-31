"""One-entry fixed-exit four-field supplement of the saved K0 first feedback."""
from __future__ import annotations
import gzip,hashlib,importlib.util,json,sys,traceback
from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext
ROOT=Path(__file__).resolve().parents[1];BRIEF=ROOT.parent.parent/'RESUME_BRIEF_2026_10_16.json';COMPONENT='n1_core_e00_fixed_missing_v1';DIR=ROOT/'tmp'/COMPONENT
SOURCE,ENTRY,RESULT=(DIR/n for n in ('source.json','entered.json','result.json'))
LOW='experiments/root_st073/n1_core_low_word_v3_receipt.json.gz';PROFILE='experiments/root_st073/n1_core_k0_feedback_m15_v2_receipt.json.gz';GAUGE='experiments/root_st073/n1_core_gauge_exit_v2_receipt.json.gz';QP='experiments/root_st073/n1_collar_qp_j2_m15_cells_v6_receipt.json.gz';OLD_LEAF='src/openai_ns_reconstruction/n1_core_low_word_green_v3.py';KERNEL_SRC='src/openai_ns_reconstruction/n1_g0_six_exit_v2.py'
LEAF='src/openai_ns_reconstruction/n1_core_e00_fixed_missing_v1.py';FIXTURE='src/openai_ns_reconstruction/n1_core_e00_fixed_missing_fixture_v1.py';RUNNER='scripts/run_n1_core_e00_fixed_missing_v1_once.py'
REQUIRED=(LOW,PROFILE,GAUGE,QP,OLD_LEAF,KERNEL_SRC,LEAF,FIXTURE,RUNNER)
FIXED={LOW:'5a5fe7c57a2efe10e58089ced608e22937c54fcf54c5f224a42abe54e450c6d6',PROFILE:'ca947df9fdf7765d9211f428375d81956af783d89074f9e782ba56aeff91dad0',GAUGE:'4c1e49a83dbc7c13e5d9f09f4601ce34b26b487405a2d14a2f066ee2dfe27921',QP:'4f894b94f8d2623f291d5ab3613c2d5fc9b98d93f95c0389a1ed22ce716a55c5',OLD_LEAF:'d34c36574f336eb3c6bcd9f3fd8f5f6a0253fef0fe30744cfc7d50b17afcd7fb',KERNEL_SRC:'0921023ed45b2bcd4471b2262a61e1019ff86ce9ec7265769de95a8d1f31e377'}
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
    c=MPIntervalContext();c.prec=1800;leaf=module(LEAF,'n1_core_e00_fixed_missing_v1');fixture=module(FIXTURE,'n1_core_e00_fixed_missing_fixture_v1');proof=fixture.run_fixture(leaf,c);require(proof['passed'],'fixture failed '+json.dumps(proof,sort_keys=True,separators=(',',':')))
    receipts={path:json.loads(gzip.decompress(raw[path])) for path in (LOW,PROFILE,GAUGE,QP)};require(all(v['result']['status']=='returned' for v in receipts.values()),'receipt status')
    low,profile,gauge,qp=(receipts[path]['result']['output'] for path in (LOW,PROFILE,GAUGE,QP));family=low['source_family'];atoms=low['source_logs'];require(all(x['source_family']==family and x['source_logs']==atoms for x in (profile,gauge,qp)),'common source family and formal atoms')
    pin=lambda r:{v['path']:v['sha256'] for v in r['source']['files']}
    require(pin(receipts[PROFILE])[LOW]==hashes[LOW] and pin(receipts[GAUGE])[LOW]==hashes[LOW] and pin(receipts[QP])[GAUGE]==hashes[GAUGE],'same saved source chain')
    src=raw[OLD_LEAF].decode('utf8');require('"K1": -K' in src and '"P1": c.mpf(1)' in src and '"d_x_F1": outer_t**3' in src and '"d_x_Uz1": outer_t' in src and 'kernel' in src,'saved physical kernel convention')
    outer=conv(c,low['result']['outer_cells']);matrices=conv(c,low['matrix_cells']);require(len(outer)==len(matrices)==8 and low['saved_W0_outer_cells_reused'],'eight saved cells')
    index={}
    for cell in range(8):
        for name in ('C0','B1'):
            for order,entries in matrices[cell][name].items():
                for key,terms in entries.items():
                    for i,t in enumerate(terms):
                        sid=f'matrix:{cell}:{name}:{order}:{key}:{i}';index[sid]={'receipt':LOW,'path':['result','output','matrix_cells',cell,name,order,key,i],'parent_origin':low['matrix_cells'][cell][name][order][key][i].get('origin')};t['source_id']=sid
        for order,cols in outer[cell]['W0_ordinary_Z'].items():
            if str(order) not in ('0','1','2'):continue
            for col,terms in cols.items():
                for i,t in enumerate(terms):
                    sid=f'W0:{cell}:{order}:{col}:{i}';index[sid]={'receipt':LOW,'path':['result','output','result','outer_cells',cell,'W0_ordinary_Z',order,col,i],'parent_origin':low['result']['outer_cells'][cell]['W0_ordinary_Z'][order][col][i].get('origin')};t['source_id']=sid
    result=leaf.build(c,outer_cells=outer,matrix_cells=matrices,source_family=family,source_logs=atoms)
    target=result['fixed_target'];require(result['schema']=='n1_core_e00_fixed_missing_v1' and result['source_schema']=='n1_core_low_word_v3' and result['source_cell_count']==8 and result['ordinary_Z_orders']==[0,1] and target['fixed'] is True and target['x'].a==target['x'].b==2 and target['t'].a==target['t'].b==1 and target['rho'].a==target['rho'].b==4,'fixed endpoint schema')
    prior_target=profile['profile']['exit_target'];require(prior_target['fixed'] is True and iv(c,prior_target['x']).a==iv(c,prior_target['x']).b==2 and iv(c,prior_target['t']).a==iv(c,prior_target['t']).b==1 and iv(c,prior_target['rho']).a==iv(c,prior_target['rho']).b==4,'prior F/U fixed endpoint')
    fields=('K1','P1','d_x_F1','d_x_Uz1');compact={k:v for k,v in result.items() if k not in ('E00_fixed','fixed_E00','E00_fixed_by_source_cell','E00_fixed_Z0','E00_fixed_Z1','E00_Z0_core_exit','E00_Z1_core_exit')};compact['fixed_fields']={};containment={}
    def term_key(t):return (t['outer_cell'],t['origin'],t['matrix_key'],str(t['matrix_order']),t['matrix_term_index'],str(t['W0_order']),t['W0_term_index'])
    for order in (0,1):
        row=result['E00_fixed_Z'+str(order)];require(set(row)==set(fields),'four field support');compact['fixed_fields']['Z'+str(order)]=row
        for field in fields:
            require(len(row[field])>0,'nonempty fixed field')
            for t in row[field]:require(t['matrix_source_id'] in index and t['W0_source_id'] in index and all(sid in index for sid in t['source_ids']),'unindexed provenance')
            old=conv(c,low['result']['W1_Z'+str(order)+'_core_exit'][field]);oldmap={term_key(t):t for t in old};newmap={term_key(t):t for t in row[field]};require(len(oldmap)==len(old) and len(newmap)==len(row[field]) and set(oldmap)==set(newmap),'saved unmerged source keys')
            inside=sum(1 for key,t in newmap.items() if t['F0_power']==oldmap[key]['F0_power'] and t['Lambda_power']==oldmap[key]['Lambda_power'] and t['coefficient'].a>=oldmap[key]['coefficient'].a and t['coefficient'].b<=oldmap[key]['coefficient'].b);require(inside==len(newmap),'fixed mass interval containment')
            containment['Z'+str(order)+'_'+field]={'count':len(newmap),'contained':inside,'proves_point_or_function_equality':False}
    compact['source_cell_masses']={cell:value['kernel_masses'] for cell,value in result['E00_fixed_by_source_cell'].items()}
    fixed_reference={};comparison={}
    for order in (0,1):
        name='W1_Z'+str(order)+'_core_exit';fixed_reference['Z'+str(order)]={}
        for field in ('F1','Uz1'):
            terms=profile['profile'][name][field];require(len(terms)>0,'saved fixed-target baseline');fixed_reference['Z'+str(order)][field]={'receipt':PROFILE,'path':['result','output','profile',name,field],'term_count':len(terms),'status':'saved_fixed_target_interval_sectors_not_recomputed'}
            comparison['Z'+str(order)+'_'+field]=profile['exit_comparison']['Z'+str(order)+'_'+field]
    return {'source_family':family,'source_logs':atoms,'source_index':index,'fixture':proof,'fixed_missing_fields':enc(compact),'fixed_vs_old_range_containment':containment,'prior_fixed_F_U_reference':fixed_reference,'prior_exit_containment_evidence':comparison,'input_receipt_sha256':{path:hashes[path] for path in (LOW,PROFILE,GAUGE,QP,OLD_LEAF,KERNEL_SRC)},'source_chain_authenticated':True,'old_prefix_disjointness_not_proved':True,'E00_is_only_K0_first_feedback':True,'F_U_not_recomputed':True,'exact_function_value_computed':False,'complete_n1_inlet':False,'m15_computed':False,'appended_to_ledger':False,'terminal_closure_proved':False}
def main():
    h,s,b=preflight();entry=save(ENTRY,{'component':COMPONENT,'source_sha256':s,'brief_sha256':b,'input_hashes':h,'entered':True});out={'component':COMPONENT,'source_sha256':s,'entry_sha256':entry,'status':'failed'}
    try:out.update(status='returned',output=consume(h),error=None)
    except BaseException as exc:out['error']={'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()}
    digest=save(RESULT,out);print(json.dumps({'component':COMPONENT,'status':out['status'],'result_sha256':digest}));return 0 if out['status']=='returned' else 1
if __name__=='__main__':raise SystemExit(main())
