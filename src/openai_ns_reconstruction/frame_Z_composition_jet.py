"""Inspect two saved Z substitution rows and derive the general C2 jet.

Only the two recorded substitution rows and their immediate variable/value
descriptors are accessed; expression subtrees and geometry closures are not
traversed. The symbolic jet theorem is conditional on source smoothness and
on the graph nodes denoting the corresponding compositions.
"""
import base64
import gzip
import hashlib
import json
import sympy as s

ARCHIVE = 'experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rp_native_C5_pulse_source.json.gz'
ARCHIVE_SHA = '442df70e5ab4051066b241e943c8d0531067d057e774a04dcc341dadbd703ebe'
REPORT = 'experiments/root_st073/lei_ren_part1_paper_compliant_current_frame_replacement_free_names.json'
REPORT_SHA = 'e26ce19fa5a566555e4e57324f3e0343f446cef5fe18bf1cbbed5c90393b7d97'
RESULT = 'tmp/frame_replacement_free_names_v3/result.json'
RESULT_SHA = '134ec429c5f0608a63894c20cd712d0e465c633e881a5d4e86fec1c532e74fff'
PROTECTED = frozenset((17,33572,25220,263344,263348,263349,262918,205147,263346))


def recover(root):
    data = (root / REPORT).read_bytes()
    if hashlib.sha256(data).hexdigest() != REPORT_SHA:
        raise ValueError('report hash changed')
    report = json.loads(data)
    raw_result = base64.b64decode(report['artifacts'][RESULT]['base64'], validate=True)
    if hashlib.sha256(raw_result).hexdigest() != RESULT_SHA:
        raise ValueError('saved result hash changed')
    prior = json.loads(raw_result)
    if prior.get('status') != 'returned' or prior['result'].get('saved_pair') != [31655,265772]:
        raise ValueError('unexpected saved replacement')
    contexts = prior['result']['saved_binder_contexts']
    binder_ids = {int(b[1]) for c in contexts for b in c['binder_stack'] if b[0] == 'substitution' and b[2] == 'Z'}
    if binder_ids != {31776,31944}:
        raise ValueError('saved Z binder ids changed')
    data = (root / ARCHIVE).read_bytes()
    if hashlib.sha256(data).hexdigest() != ARCHIVE_SHA:
        raise ValueError('raw gzip hash mismatch before decode')
    nodes = json.loads(gzip.decompress(data))['exact_graph_nodes']
    records = []
    for node in sorted(binder_ids):
        row = nodes[node]
        if row.get('operation') != 'function_substitution':
            raise ValueError('saved Z binder is not function_substitution')
        var, value, expression = (row.get(k) for k in ('variable','value','expression'))
        if any(type(i) is not int or not 0 <= i < len(nodes) for i in (var,value,expression)):
            raise ValueError('invalid saved substitution operand')
        if var in PROTECTED or value in PROTECTED:
            raise ValueError('protected immediate descriptor; cannot inspect')
        variable_desc, value_desc = nodes[var], nodes[value]
        if variable_desc != {'operation':'bound_variable','name':'Z'}:
            raise ValueError('saved variable is not the declared Z atom')
        records.append({'node':node,'variable_ref':var,'variable_descriptor':variable_desc,
                        'expression_ref':expression,'value_ref':value,'value_descriptor':value_desc,
                        'other_row_metadata':{k:v for k,v in row.items() if k not in ('operation','variable','expression','value')}})
    # Arbitrary quadratic Taylor jets suffice to derive the C2 chain rule;
    # z and q are independent before substituting q=a(z).
    z,q = s.symbols('z q')
    a0,a1,a2,f0,fz,fq,fzz,fzq,fqq = s.symbols('a0 a1 a2 f0 fz fq fzz fzq fqq')
    F = f0+fz*z+fq*(q-a0)+fzz*z**2/2+fzq*z*(q-a0)+fqq*(q-a0)**2/2
    a = a0+a1*z+a2*z**2/2
    H = F.subs(q,a)
    checks = {'value':s.expand(H.subs(z,0)-f0)==0,
              'first':s.expand(s.diff(H,z).subs(z,0)-(fz+fq*a1))==0,
              'second':s.expand(s.diff(H,z,2).subs(z,0)-(fzz+2*fzq*a1+fqq*a1**2+fq*a2))==0,
              'fixed_argument_first':s.expand((fz+fq*a1).subs(a1,0)-fz)==0,
              'fixed_argument_second':s.expand((fzz+2*fzq*a1+fqq*a1**2+fq*a2).subs({a1:0,a2:0})-fzz)==0}
    return {'schema':'frame_Z_composition_jet_v1','status':'descriptor_and_conditional_C2_jet',
            'source_archive_sha256':ARCHIVE_SHA,'previous_result_sha256':RESULT_SHA,
            'binder_records':records,'symbolic_checks':checks,'symbolic_checks_passed':all(checks.values()),
            'jet_formula':{'first':'F_z + F_q*a1','second':'F_zz + 2*F_zq*a1 + F_qq*a1^2 + F_q*a2'},
            'conditions':['same source functions and ordinary C2 jets','substitution denotes composition','argument jet is actual ordinary Z jet'],
            'expression_subtrees_traversed':False,'value_subtrees_traversed':False,'geometry_closure_traversed':False,
            'actual_argument_derivatives_proved':False,'actual_inlet_jet_transport_proved':False,
            'actual_capture_free_proved':False,'physical_endpoint_proved':False,'source_values_computed':False,
            'old_errors_removed':False,'runtime_model_metadata':None}
