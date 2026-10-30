"""One-entry signed source-jet join at the core/first seam; no old replay."""
from __future__ import annotations
import gzip
import hashlib
import importlib.util
import json
import sys
import traceback
from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext

ROOT=Path(__file__).resolve().parents[1]
BRIEF=ROOT.parent.parent/'RESUME_BRIEF_2026_10_16.json'
COMPONENT='n1_word2_source_inlet_v1'
DIR=ROOT/'tmp'/COMPONENT
SOURCE,ENTRY,RESULT=(DIR/name for name in ('source.json','entered.json','result.json'))
BASE='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_'
OP=BASE+'inner_operator.json.gz'; OP_CHECK=BASE+'inner_operator_check.json'; OP_SOURCE=BASE+'inner_operator.py'
COUPLED=BASE+'coupled_first_step.json.gz'; COUPLED_CHECK=BASE+'coupled_first_step_check.json'; COUPLED_SOURCE=BASE+'coupled_first_step.py'
PRIOR='experiments/root_st073/n1_inlet_residual_join_v1_receipt.json.gz'
LEAF='src/openai_ns_reconstruction/n1_word2_source_integrand.py'
RUNNER='scripts/run_n1_word2_source_inlet_once.py'
REQUIRED=(OP,OP_CHECK,OP_SOURCE,COUPLED,COUPLED_CHECK,COUPLED_SOURCE,PRIOR,LEAF,RUNNER)
FIXED={OP:'63b181010f1715d3b20589b06c1c2139ef63bcc73d7b0f3c5b7c43d9c5e03f52',
    OP_CHECK:'846cfedcbd0355a4a8e6a5ed445558d367cd55c7459085dd9c5d74821d969bd4',
    OP_SOURCE:'268217c6b33bc717472742fa64d76a853c6760a956abfad92ac45bcf385df014',
    COUPLED:'30b807af4147658ec0978ebd253cc1b2b31bed67bfed430a965e66461450add1',
    COUPLED_CHECK:'6a58caf1db41697620820830cbe7f76620874a8b429b351ad3d58f637e4e1416',
    COUPLED_SOURCE:'97913b82b8570b220921cd89918569a4019ff3495bdbea0548b3c53e9d8c93e7',
    PRIOR:'6eb7552ea8253c6ee4ab8a0f0ed507fc8c24b5e083013b0a4d4dac865ce10c5c'}
FIELDS=('F1','Uz1','K1','P1','d_x_F1','d_x_Uz1')


def sha(x): return hashlib.sha256(x).hexdigest()
def require(x,message):
    if not x: raise ValueError(message)
def exclusive(path,obj):
    raw=(json.dumps(obj,sort_keys=True,separators=(',',':'))+'\n').encode()
    with path.open('xb') as f:f.write(raw)
    return sha(raw)


def preflight():
    require(not ENTRY.exists() and not RESULT.exists(),'already entered')
    b=BRIEF.read_bytes(); brief=json.loads(b)
    require(brief.get('component')==COMPONENT and brief.get('phase')==COMPONENT+'_pending'
        and brief.get('goal_active') is True and brief.get('full_goal_complete') is False
        and brief.get('one_entry_only') is True and brief.get('do_not_reenter') is True,'brief guard')
    ck=brief['source_checkpoint']; raw=SOURCE.read_bytes()
    require(Path(ck['path']).resolve()==SOURCE.resolve() and sha(raw)==ck['sha256'],'checkpoint')
    source=json.loads(raw); require(source['component']==COMPONENT and source['one_entry_only'] is True,'manifest')
    hashes={}
    for row in source['files']:
        p=row['path']; require(p in REQUIRED and p not in hashes,'unexpected input'); digest=sha((ROOT/p).read_bytes())
        require(digest==row['sha256'] and (p not in FIXED or digest==FIXED[p]),'changed input '+p); hashes[p]=digest
    require(set(hashes)==set(REQUIRED),'missing input')
    return hashes,sha(raw),sha(b)


def iv(c,record):
    lo=tuple(record['lower_exact_mpf_tuple']); hi=tuple(record['upper_exact_mpf_tuple'])
    require(len(lo)==len(hi)==4,'four-entry exact endpoints required')
    x=c.make_mpf((lo,hi)); require(x.a<=x.b,'reversed interval'); return x


def encode(x):
    if hasattr(x,'_mpi_'):
        return {'lower_exact_mpf_tuple':list(x._mpi_[0]),'upper_exact_mpf_tuple':list(x._mpi_[1]),'directed':True}
    if isinstance(x,dict):return {str(k):encode(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [encode(v) for v in x]
    if isinstance(x,(str,int,float,bool)) or x is None:return x
    raise TypeError('unsupported output '+type(x).__name__)


def consume(hashes):
    raw={}
    for p in REQUIRED:
        data=(ROOT/p).read_bytes(); require(sha(data)==hashes[p],'changed after entry '+p); raw[p]=data
    operator=json.loads(gzip.decompress(raw[OP])); coupled=json.loads(gzip.decompress(raw[COUPLED]));
    oc=json.loads(raw[OP_CHECK]); cc=json.loads(raw[COUPLED_CHECK]); prior=json.loads(gzip.decompress(raw[PRIOR]))
    require(oc['all_passed'] is True and cc['all_passed'] is True and prior['result']['status']=='returned','accepted inputs')
    require(oc['input_hashes'][Path(OP).name]==hashes[OP] and cc['input_hashes'][Path(COUPLED).name]==hashes[COUPLED], 'check/archive binding')
    a=operator['samples']['first_inlet']; w=coupled['samples']['core_exit']; earlier=prior['result']['output']
    require(a['source_family']==w['source_family']==operator['source_family']==coupled['source_family']==earlier['source_family'], 'source family')
    require(a['Z']==w['Z']==earlier['selected_Z_interval'], 'selected Z mismatch')
    require(a['source_log_basis_names']==w['source_log_basis_names']==['hb','F0_anchor','Pstar_squared','Lambda']
        and a['source_log_bases']==w['source_log_bases'], 'source log anchors mismatch')
    require(a['amplitude_anchor_fixed_during_jet_algebra'] is True
        and a['ordinary_Z_amplitude_derivatives_in_coefficients'] is True
        and a['outputs_are_enclosures_of_original_functions_not_selected_point_values'] is True,'physical jet semantics')
    require(w['computed_preconditioned_indices']==[0,1] and w['initial_Z_derivative_order']==2
        and w['actual_whole_cell_initial_function_used'] is True,'whole-cell W0 Z2 source')
    require(a['geometry']['chart']=='first' and w['endpoint']['chart']=='core','seam charts')
    c=MPIntervalContext(); c.prec=1800
    require(iv(c,a['geometry']['phase'])==0 and iv(c,a['geometry']['rho_enclosure_only'])==4
        and iv(c,w['endpoint']['coordinate'])==4,'physical seam coordinate')
    src=raw[COUPLED_SOURCE].decode('utf-8'); op_source=raw[OP_SOURCE].decode('utf-8')
    require('Wbase=(I+J)Gg; H=(I+J)G(B0-Akin+B1*d_Z)' in src
        and "initial_function={name:op.export(value,2)" in src
        and "ordinary_Z_derivatives=[math.factorial(k)*v" in op_source, 'actual source jet conventions')

    def terms(row,order):
        require(row['ordinary_Z_derivative_order']>=order,'insufficient ordinary Z derivative')
        out=[]
        for t in row['terms']:
            exponent=t['source_exponents']; require(len(exponent)==4 and all(type(n) is int for n in exponent),'source exponent')
            out.append({'source_exponents':exponent,'coefficient_interval':iv(c,t['ordinary_Z_derivatives'][order]),
                        'combined_source_log_interval':iv(c,t['combined_source_log'])})
        return out

    W0={order:{i:terms(w['initial_function'][field],order) for i,field in enumerate(FIELDS,1)} for order in range(3)}
    matrices={}
    for name in ('B0','B1'):
        matrices[name]={order:{} for order in range(3)}
        for key,row in a['regular_system'][name].items():
            pair=tuple(int(v) for v in key.split(',')); require(len(pair)==2,'matrix key')
            for order in range(3):matrices[name][order][pair]=terms(row,order)
    spec=importlib.util.spec_from_file_location('n1_word2_source_integrand_leaf',ROOT/LEAF)
    require(spec is not None and spec.loader is not None,'leaf import'); module=importlib.util.module_from_spec(spec)
    sys.modules[spec.name]=module; spec.loader.exec_module(module)
    result=module.assemble(matrices['B0'],matrices['B1'],W0)
    require(result['W0_ZZZ_computed'] is False and result['W1_ZZ_integrated'] is False,'scope overclaimed')
    return {'source_family':a['source_family'],'selected_Z_interval':a['Z'],
        'coordinate':'rho=4, s=0; core W0 and actual first-collar operator',
        'source_and_quadrature_intervals_retained':True,'same_four_source_anchors':True,
        'computed_first_word_integrand_not_integral':True,'result':encode(result)}


def main():
    hashes,s,b=preflight(); entry=exclusive(ENTRY,{'component':COMPONENT,'source_sha256':s,'brief_sha256':b,'input_hashes':hashes,'entered':True})
    out={'component':COMPONENT,'source_sha256':s,'entry_sha256':entry,'status':'failed'}
    try:out['output']=consume(hashes);out['status']='returned';out['error']=None
    except BaseException as exc:out['error']={'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()}
    digest=exclusive(RESULT,out);print(json.dumps({'component':COMPONENT,'status':out['status'],'result_sha256':digest}))
    return 0 if out['status']=='returned' else 1


if __name__=='__main__':raise SystemExit(main())
