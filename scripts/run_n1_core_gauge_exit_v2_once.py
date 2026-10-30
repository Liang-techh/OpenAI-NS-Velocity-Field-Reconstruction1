"""One-entry source-bound gauge evidence and fixed core raw-Q exits."""
from __future__ import annotations
import gzip,hashlib,importlib.util,json,sys,traceback
from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext
ROOT=Path(__file__).resolve().parents[1];BRIEF=ROOT.parent.parent/'RESUME_BRIEF_2026_10_16.json';COMPONENT='n1_core_gauge_exit_v2';DIR=ROOT/'tmp'/COMPONENT
SOURCE,ENTRY,RESULT=(DIR/n for n in ('source.json','entered.json','result.json'))
LOW='experiments/root_st073/n1_core_low_word_v3_receipt.json.gz';KCHI='experiments/root_st073/n1_core_kchi_first_v1_receipt.json.gz';COMP='experiments/root_st073/n1_core_kchi_compensation_v2_receipt.json.gz';BRIDGE='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_inner_numeric_bridge_check.json';BRIDGE_SRC='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_inner_numeric_bridge.py';BRIDGE_CHECK_SRC='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_inner_numeric_bridge_check.py'
WEIGHT='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_weighted_propagator.py';NATURAL='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_natural_component_enclosure.py';REGULAR='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_regular_series_function.py';RESIDUAL='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_weighted_feedback_bound.py';OP='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_inner_operator.py'
EXIT='src/openai_ns_reconstruction/n1_kchi_feedback_exit_v2.py';EF='src/openai_ns_reconstruction/n1_kchi_feedback_exit_fixture_v2.py';GAUGE='src/openai_ns_reconstruction/n1_gauge_conjugation_identity_v2.py';GF='src/openai_ns_reconstruction/n1_gauge_conjugation_identity_fixture_v2.py';RUNNER='scripts/run_n1_core_gauge_exit_v2_once.py'
FAILURE='experiments/root_st073/n1_core_gauge_exit_v1_failure.json.gz'
REQUIRED=(LOW,KCHI,COMP,BRIDGE,BRIDGE_SRC,BRIDGE_CHECK_SRC,WEIGHT,NATURAL,REGULAR,RESIDUAL,OP,EXIT,EF,GAUGE,GF,RUNNER,FAILURE)
FIXED={FAILURE:'127ce9804dbe7c26e0e848dfd912d61929da861d5cb808dd6246dd3233a2894c',LOW:'5a5fe7c57a2efe10e58089ced608e22937c54fcf54c5f224a42abe54e450c6d6',KCHI:'b0d654203e8c4acabdd0b61ca7b79c1b7050314b0c003a0c2f8757bde94973a0',COMP:'c862f0560e8cb4b053abbbfcc982e4b4c1b602a806628e2d709399878e1c3aaf',WEIGHT:'5219a968018601d8c445815058f7c1c832398dc417fcb51c689688d03011ba62',REGULAR:'b64b7cae8529cfc5caa158fe06ac5660ec1d2217fc4d845d02002dac3a195555',OP:'268217c6b33bc717472742fa64d76a853c6760a956abfad92ac45bcf385df014'}
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
        b=(ROOT/path).read_bytes();require(sha(b)==hashes[path],'changed after entry '+path);raw[path]=b
    c=MPIntervalContext();c.prec=1800;leaf=module(EXIT,'n1_kchi_feedback_exit_v2');ef=module(EF,'n1_kchi_feedback_exit_fixture_v2');gauge=module(GAUGE,'n1_gauge_conjugation_identity_v2');gf=module(GF,'n1_gauge_conjugation_identity_fixture_v2')
    proofs={'exit':ef.run_fixture(leaf.build_feedback_exit,c),'gauge':gf.run(gauge)}
    for name,v in proofs.items():require(v['passed'],name+' fixture failed '+str([q for q,x in v['checks'].items() if not x]))
    low=json.loads(gzip.decompress(raw[LOW]));kc=json.loads(gzip.decompress(raw[KCHI]));comp=json.loads(gzip.decompress(raw[COMP]));bridge=json.loads(raw[BRIDGE])
    for r in (low,kc,comp):require(r['result']['status']=='returned','source receipt status')
    a=low['result']['output'];k=kc['result']['output'];p=comp['result']['output'];family=a['source_family'];atoms=a['source_logs']
    require(k['source_family']==p['source_family']==family and k['source_logs']==p['source_logs']==atoms,'same source family and formal atoms')
    kpins={v['path']:v['sha256'] for v in kc['source']['files']};cpins={v['path']:v['sha256'] for v in comp['source']['files']};require(kpins[LOW]==cpins[LOW]==hashes[LOW] and cpins[KCHI]==hashes[KCHI],'authenticated previous source chain')
    src={path:raw[path].decode('utf8') for path in (WEIGHT,NATURAL,REGULAR,RESIDUAL,OP,BRIDGE_SRC,BRIDGE_CHECK_SRC)}
    require('S^-1 B0 S+S^-1 B1 S_Z' in src[WEIGHT] and 'potential=algebra.lift(-2*chi)' in src[WEIGHT] and "residual['5,1']=residual['5,1']-potential" in src[WEIGHT] and 'def conjugate_extra' in src[NATURAL] and "B0['5,1']=B0['5,1']+algebra.lift(-2*chi)" in src[REGULAR] and 'all 72 matrix entries' in src[BRIDGE_CHECK_SRC].lower(),'source identity formulas')
    saved=bridge.get('physical_conjugation_proof',{});owner=bridge.get('source_identity',{});pins=bridge.get('input_hashes',{}); pin_audit={}
    for path in (WEIGHT,NATURAL,REGULAR,RESIDUAL,OP,BRIDGE_SRC):
        name=Path(path).name;pin_audit[name]={'present':name in pins,'matches_current':pins.get(name)==hashes[path] if name in pins else False}
    historical={'check_sha256':hashes[BRIDGE],'all_passed':bridge.get('all_passed') is True,'source_family_matches':bridge.get('source_family')==family,'saved_source_identity_passed':owner.get('passed') is True,'saved_proof_passed':saved.get('passed') is True,'saved_matrix_entry_count':saved.get('all_exact_physical_conjugation_matrix_entries'),'saved_forcing_rows':saved.get('exact_physical_forcing_conjugation_rows'),'saved_relative_amplitude_rows':saved.get('original_relative_amplitude_ordinary_Z0_to_Z3_rows'),'source_pin_audit':pin_audit,'source_identity_reexecuted':False}
    index={};matrices=conv(c,a['matrix_cells']);w0=conv(c,a['result']['outer_cells']);dcells=k['first_correction']['target_cells'];tail=k['W0_separate_symmetric_kernel_tail'];require(len(matrices)==len(w0)==len(dcells)==8,'eight source cells')
    for cell in range(8):
        for name in ('C0','B1'):
            for order,entries in matrices[cell][name].items():
                for key,terms in entries.items():
                    for i,t in enumerate(terms):
                        sid=f'matrix:{cell}:{name}:{order}:{key}:{i}';index[sid]={'receipt':LOW,'path':['matrix_cells',cell,name,order,key,i]};t['source_id']=sid;t['origin']={'source_ids':[sid],'error_classes':['whole_core_matrix_interval'],'uncertainty':True}
        for order,cols in w0[cell]['W0_ordinary_Z'].items():
            if order not in ('0','1','2'):continue
            for col,terms in cols.items():
                for i,t in enumerate(terms):
                    sid=f'W0:{cell}:{order}:{col}:{i}';index[sid]={'receipt':LOW,'path':['result','outer_cells',cell,'W0_ordinary_Z',order,col,i]};t['source_id']=sid;t['origin']={'source_ids':[sid],'error_classes':['whole_core_W0_interval'],'uncertainty':True}
    delta=[];tailcells=[]
    for cell,rawcell in enumerate(dcells):
        t=w0[cell]['t'];require((2*t)._mpi_==iv(c,rawcell['x'])._mpi_,'same source cell');d={'t':t,'W0_ordinary_Z':{str(n):{str(col):[] for col in range(1,7)} for n in range(3)}};e={'t':t,'W0_ordinary_Z':{str(n):{str(col):[] for col in range(1,7)} for n in range(3)}}
        for field,col in (('F1','1'),('d_x_F1','5')):
            for order in (0,1,2):
                terms=conv(c,rawcell[field]['Z'+str(order)]);groups={}
                for pos,term in enumerate(terms):
                    key=(tuple(term['powers']),term['chi_order'],term['source_order']);g=groups.get(key)
                    if g is None:g={'coefficient':term['coefficient'],'ids':set(),'classes':set(),'parents':[]};groups[key]=g
                    else:g['coefficient']+=term['coefficient']
                    g['ids'].update(term['source_ids']);g['classes'].update(term['error_classes']);g['parents'].append(pos)
                vals=[]
                for i,(key,g) in enumerate(groups.items()):
                    sid=f'delta:{cell}:{field}:{order}:{i}';index[sid]={'receipt':KCHI,'path':['first_correction','target_cells',cell,field,'Z'+str(order)],'term_indices':g['parents']};vals.append({'powers':key[0],'coefficient':g['coefficient'],'source_id':sid,'origin':{'source_ids':sorted(g['ids'])+[sid],'error_classes':sorted(g['classes']),'uncertainty':True,'aggregation':'outward only'}})
                d['W0_ordinary_Z'][str(order)][col]=vals;ev=conv(c,tail[field][str(order)])
                for i,term in enumerate(ev):
                    sid=f'W0tail:{field}:{order}:{i}';index[sid]={'receipt':KCHI,'path':['W0_separate_symmetric_kernel_tail',field,str(order),i]};term['source_id']=sid;term['origin']={'source_ids':term['source_ids']+[sid],'error_classes':term['error_classes'],'uncertainty':True,'error_only':True};term['error_only']=True
                e['W0_ordinary_Z'][str(order)][col]=ev
        delta.append(d);tailcells.append(e)
    chi=conv(c,k['chi']['chi']);chirows={n:[{'powers':[0,0,0,0],'coefficient':chi[n],'source_id':f'chi:Z{n}','origin':{'source_ids':[f'chi:Z{n}'],'error_classes':['whole_Z_chi_interval'],'uncertainty':True}}] for n in (0,1)}
    for n in (0,1):index[f'chi:Z{n}']={'receipt':KCHI,'path':['chi','chi',n]}
    result=leaf.build_feedback_exit(c,w0,matrices,chirows,delta_w0_cells=delta,W0_kernel_tail_cells=tailcells)
    kt=conv(c,k['kernel_tail_bounds']);outer_error={}
    for kind,qmap in (('R_Q_W0',result['Q_W0_by_source_cell']),('R_Q_Delta',result['Q_DeltaW0_by_source_cell']),('R_Q_Rg',result['Q_W0_kernel_tail_by_source_cell'])):
        outer_error[kind]={}
        for order in (0,1):
            fields={}
            for field in ('F1','d_x_F1'):
                terms=[]
                for j in range(8):
                    for qorder in range(order+1):
                        for v in qmap[str(j)]['Z'+str(qorder)]['5']:
                            upper=c.mpf(1)/4*kt[field][order-qorder]*abs(v['coefficient']);sym=c.make_mpf(((-upper)._mpi_[0],upper._mpi_[1]));ids=list(v.get('source_ids',[]))+['chi:Z0','chi:Z1']
                            terms.append({'powers':v['powers'],'coefficient':sym,'source_ids':ids,'error_classes':list(v.get('error_classes',[]))+[kind,'outer_Kchi_kernel_tail'],'uncertainty':True,'error_only':True,'source_cell':j})
                fields[field]=terms
            outer_error[kind][order]=fields
    # P=2 chi E51 acts only on F1. The previously saved Pq for finite W0
    # and DeltaW0 is reused; P(Rg) is formed here as a symmetric error.
    pcomp=p['potential_compensation_exit'];pq_w0={};pq_delta={};pq_tail={}
    def symterm(value,powers,ids,classes,kind,cell):
        radius=abs(value).b; symmetric=c.make_mpf(((-radius)._mpi_[0],radius._mpi_[1]))
        return {'powers':powers,'coefficient':symmetric,'source_ids':ids,'error_classes':classes+[kind],'uncertainty':True,'error_only':True,'source_cell':cell}
    for j in range(8):
        pq_w0[str(j)]=conv(c,pcomp['G0_P_W0']['exit_by_source_cell'][str(j)]['P_q5'])
        pq_delta[str(j)]=conv(c,pcomp['G0_P_DeltaW0']['exit_by_source_cell'][str(j)]['P_q5'])
        pq_tail[str(j)]={'0':[],'1':[]}
        for order in (0,1):
            rules=((0,0),) if order==0 else ((1,0),(0,1))
            for chi_order,tail_order in rules:
                for i,v in enumerate(conv(c,tail['F1'][str(tail_order)])):
                    ids=list(v['source_ids'])+[f'W0tail:F1:{tail_order}:{i}',f'chi:Z{chi_order}']
                    pq_tail[str(j)][str(order)].append(symterm(2*chi[chi_order]*v['coefficient'],v['powers'],ids,list(v['error_classes']),'P_on_W0_kernel_tail',j))
    potential_error={name:{0:{'F1':[],'d_x_F1':[]},1:{'F1':[],'d_x_F1':[]}} for name in ('G0_P_Rg','D_P_Rg','R_P_W0','R_P_Delta','R_P_Rg')}
    for j in range(8):
        lo=c.mpf(j)/8;hi=c.mpf(j+1)/8;gm=leaf.cell_kernel_masses(c,c.mpf(1),lo,hi);dm=leaf.delta_kernel_masses(c,c.mpf(2),2*lo,2*hi)
        for order in (0,1):
            for field in ('F1','d_x_F1'):
                for v in pq_tail[str(j)][str(order)]:
                    potential_error['G0_P_Rg'][order][field].append(symterm(gm[field]*v['coefficient'],v['powers'],v['source_ids'],v['error_classes'],'G0_after_P_R',j))
                rules=((0,0),) if order==0 else ((1,0),(0,1))
                for chi_order,qorder in rules:
                    for v in pq_tail[str(j)][str(qorder)]:
                        potential_error['D_P_Rg'][order][field].append(symterm(-2*chi[chi_order]*dm[field]*v['coefficient'],v['powers'],v['source_ids']+[f'chi:Z{chi_order}'],v['error_classes'],'D_after_P_R',j))
                for name,packet in (('R_P_W0',pq_w0),('R_P_Delta',pq_delta),('R_P_Rg',pq_tail)):
                    for qorder in range(order+1):
                        for v in packet[str(j)][str(qorder)]:
                            bound=c.mpf(1)/4*kt[field][order-qorder]*abs(v['coefficient']);ids=list(v.get('source_ids',[]))+['chi:Z0','chi:Z1']
                            potential_error[name][order][field].append(symterm(bound,v['powers'],ids,list(v.get('error_classes',[])),'outer_Kchi_kernel_tail',j))
    compact={name:{'aggregate':branch['aggregate'],'exit_by_source_cell':branch['exit_by_source_cell'],'structural_zero_fields':branch['structural_zero_fields']} for name,branch in result['raw_Q_branches'].items()}
    errors={name:{'aggregate':branch['aggregate'],'exit_by_source_cell':branch['exit_by_source_cell']} for name,branch in result['symmetric_R_branches'].items()}
    return {'source_family':family,'source_logs':atoms,'historical_source_proof':historical,'fixtures':proofs,'input_index':index,'raw_Q_fixed_exit':enc(compact),'inner_W0_kernel_tail_fixed_exit':enc(errors),'outer_kernel_tail_fixed_exit':enc(outer_error),'potential_compensation_kernel_tail_fixed_exit':enc(potential_error),
      'saved_W0_first_correction_fixed_exit':p['signed_delta_W0_fixed_exit'],'saved_potential_compensation_fixed_exit':p['potential_compensation_exit'],
      'saved_baseline_E00_receipt':LOW,'baseline_E00_recomputed':False,'W0_baseline_last_cell_is_only_interval_enclosure':True,
      'P_compensation_R_terms_included':True,'feedback_word_tail_included':False,'full_six_component_inlet_computed':False,'appended_to_cumulative_m15':False,'complete_m15_computed':False,'terminal_closure_proved':False}
def main():
    h,s,b=preflight();entry=save(ENTRY,{'component':COMPONENT,'source_sha256':s,'brief_sha256':b,'input_hashes':h,'entered':True});out={'component':COMPONENT,'source_sha256':s,'entry_sha256':entry,'status':'failed'}
    try:out.update(status='returned',output=consume(h),error=None)
    except BaseException as exc:out['error']={'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()}
    digest=save(RESULT,out);print(json.dumps({'component':COMPONENT,'status':out['status'],'result_sha256':digest}));return 0 if out['status']=='returned' else 1
if __name__=='__main__':raise SystemExit(main())
