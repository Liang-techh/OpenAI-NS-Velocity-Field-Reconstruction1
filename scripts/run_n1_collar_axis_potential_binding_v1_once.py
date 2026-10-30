"""One-entry binding of the saved axis potential to the native collar kernel."""
from __future__ import annotations
import ast,gzip,hashlib,json,traceback
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];BRIEF=ROOT.parent.parent/'RESUME_BRIEF_2026_10_16.json';COMPONENT='n1_collar_axis_potential_binding_v1';DIR=ROOT/'tmp'/COMPONENT
SOURCE,ENTRY,RESULT=(DIR/n for n in ('source.json','entered.json','result.json'))
KCHI='experiments/root_st073/n1_core_kchi_first_v1_receipt.json.gz';BASE='experiments/root_st073/n1_core_six_inlet_v2_receipt.json.gz';GAUGE='experiments/root_st073/n1_core_gauge_exit_v2_receipt.json.gz';CV4='experiments/root_st073/n1_kchi_collar_volterra_v4_receipt.json.gz';WEIGHT='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_weighted_propagator.py';COREPY='experiments/root_st073/lei_ren_part1_paper_compliant_core_physical_field.py';SERIES='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_regular_series_function.py';COLLAR='src/openai_ns_reconstruction/n1_kchi_collar_volterra_v4.py';RUNNER='scripts/run_n1_collar_axis_potential_binding_v1_once.py'
REQUIRED=(KCHI,BASE,GAUGE,CV4,WEIGHT,COREPY,SERIES,COLLAR,RUNNER)
FIXED={KCHI:'b0d654203e8c4acabdd0b61ca7b79c1b7050314b0c003a0c2f8757bde94973a0',BASE:'59c0ff88652d83d7b92ba9c3a564ee90eac5ae794220b8b858c112936b76342d',GAUGE:'4c1e49a83dbc7c13e5d9f09f4601ce34b26b487405a2d14a2f066ee2dfe27921',CV4:'eb10c032856d4c44a43c5fe5c2b394663c85a164f4d35344a3c71f54fac6edb1',WEIGHT:'5219a968018601d8c445815058f7c1c832398dc417fcb51c689688d03011ba62',COREPY:'19acdd5f7a506e7582427f004be8e2cd7ae22771eafc76b66ee44d053b9aec54',SERIES:'b64b7cae8529cfc5caa158fe06ac5660ec1d2217fc4d845d02002dac3a195555'}
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
def methods(source):
    tree=ast.parse(source);cls=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=='CurrentOriginalN1WeightedPropagator');return {n.name:n for n in cls.body if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef))}
def consume(hashes):
    raw={}
    for p in REQUIRED:
        b=(ROOT/p).read_bytes();require(sha(b)==hashes[p],'changed after entry '+p);raw[p]=b
    kchi=json.loads(gzip.decompress(raw[KCHI]));base=json.loads(gzip.decompress(raw[BASE]));gauge=json.loads(gzip.decompress(raw[GAUGE]));cv4=json.loads(gzip.decompress(raw[CV4]));require(all(r['result']['status']=='returned' for r in (kchi,base,gauge,cv4)),'source receipt status')
    kp={v['path']:v['sha256'] for v in kchi['source']['files']};bp={v['path']:v['sha256'] for v in base['source']['files']};gp={v['path']:v['sha256'] for v in gauge['source']['files']};cp={v['path']:v['sha256'] for v in cv4['source']['files']}
    require(kp[WEIGHT]==hashes[WEIGHT] and kp[COREPY]==hashes[COREPY] and kp[SERIES]==hashes[SERIES] and bp[GAUGE]==hashes[GAUGE] and gp[KCHI]==hashes[KCHI] and cp[COLLAR]==hashes[COLLAR],'pinned source chain')
    k=kchi['result']['output'];b=base['result']['output'];g=gauge['result']['output'];require(k['source_family']==b['source_family']==g['source_family'] and k['source_logs']==b['source_logs']==g['source_logs'],'same source family and formal atom records')
    weight=raw[WEIGHT].decode('utf8');core=raw[COREPY].decode('utf8');series=raw[SERIES].decode('utf8');leaf=raw[COLLAR].decode('utf8');m=methods(weight);init=ast.unparse(m['__init__']);graph=ast.unparse(m['assert_graph']);operator=ast.unparse(m['_operator']);evaluate=ast.unparse(m['evaluate'])
    checks={'native_axis_from_same_physical_owner':'self.physical = self.inputs.n1.source.core' in init and 'self.physical is not self.inputs.n1.source.core' in graph,
      'native_axis_jet_uses_Z_only':"axis = self.physical.axis_inputs(z)" in operator and "chi = op.Jet(self.c, axis['chi'].coefficients[:3])" in operator,
      'principal_potential_is_minus_2chi':"potential = algebra.lift(-2 * chi)" in operator and "residual['5,1'] = residual['5,1'] - potential" in operator,
      'one_query_jet_shared_in_integration':"self._operator(z, q['actual_source_coordinate'], chart)" in evaluate and 'potential_kernel(chi, t, X, N, 1)' in evaluate,
      'both_core_and_first_collar_use_integral':"'core')" in evaluate and 'cell_integral(cell, xcell, core_end / core_cells, False)' in evaluate and 'cell_integral(cell, xcell, xcell * phase / (2 * collar_cells), True)' in evaluate,
      'axis_formula_pinned':'chi=h2/denominator' in core and 'u=4*z+self.j; d=one-z*z; L=one-(z*z)*self.delta; H=z*((1-self.delta)/2)+d*u' in core,
      'series_principal_sign_pinned':'Kchi=sum_j>=0(-2chi)^j' in series,
      'adapter_uses_same_kernel_equation':'t * (1 - (t / x) ** 2) / 2' in leaf and '(-2) * k * chi[j] * previous["F"][j]' in leaf,
      'previous_gauge_evidence_family_match':g['historical_source_proof']['source_family_matches'] is True and g['historical_source_proof']['saved_proof_passed'] is True}
    require(all(checks.values()),'source check failed '+str([key for key,value in checks.items() if not value]))
    jets=k['chi']['chi'];require(len(jets)>=2,'ordinary chi jet');chi={}
    for n in (0,1):
        v=jets[n];require(v['directed'] is True and len(v['lower_exact_mpf_tuple'])==len(v['upper_exact_mpf_tuple'])==4,'directed saved jet');chi[str(n)]={'coefficient':v,'source_id':f'chi:Z{n}','source_receipt':KCHI,'source_record_path':['result','output','chi','chi',n],'uncertainty_class':'whole_Z_axis_potential_interval'}
    return {'source_family':k['source_family'],'source_logs':k['source_logs'],'source_checks':checks,'saved_chi_Z0_Z1':chi,'domain':'real Z in [-1,1]','source_function_identity':'same axis_inputs(z) chi jet captured once by the native weighted propagator and used for core and first collar integration','ordinary_Z1_convention':'first Taylor-jet coefficient equals first ordinary derivative','seam_function_identity_scope':'for the native weighted axis gauge, same chi(Z) on both sides of the core/first collar seam','phase_shear_cutoff_identified_with_axis_potential':False,'actual_live_source_graph_reexecuted':False,'point_values_computed':False,'actual_collar_transport_run':False,'local_collar_forcing_bound':False,'feedback_word_tail_included':False,'terminal_closure_proved':False}
def main():
    h,s,b=preflight();entry=save(ENTRY,{'component':COMPONENT,'source_sha256':s,'brief_sha256':b,'input_hashes':h,'entered':True});out={'component':COMPONENT,'source_sha256':s,'entry_sha256':entry,'status':'failed'}
    try:out.update(status='returned',output=consume(h),error=None)
    except BaseException as exc:out['error']={'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()}
    digest=save(RESULT,out);print(json.dumps({'component':COMPONENT,'status':out['status'],'result_sha256':digest}));return 0 if out['status']=='returned' else 1
if __name__=='__main__':raise SystemExit(main())
