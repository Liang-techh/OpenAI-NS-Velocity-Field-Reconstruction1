"""Bind the observed axis-context inlet port and distinguish its jet orders.

This leaf reads fixed saved native port lists and the saved mapping, without
walking expression edges or evaluating source functions. The axis/derivative
identities below are independent symbolic examples with arbitrary coefficients.
"""
import base64
import gzip
import hashlib
import json
import sympy as s

P='lei_ren_part1_paper_compliant_current_original_'
ARCHIVE='experiments/root_st073/'+P+'Rp_native_C5_pulse_source.json.gz'
ARCHIVE_SHA='442df70e5ab4051066b241e943c8d0531067d057e774a04dcc341dadbd703ebe'
REPORT='experiments/root_st073/lei_ren_part1_paper_compliant_current_frame_Z_composition_jet.json'
REPORT_SHA='1a075e0a10e9df6589e924f7ba25a58f2417a09c6d384ee319aece9b08054919'
RESULT='tmp/frame_Z_composition_jet_v1/result.json'
RESULT_SHA='9a45b8bc308e4ac81fffcc9887e2e1ed39c9897857cabb21032521ab0fe1e931'
INPUTS=(P+'n1_pulse_C3_functions.json.gz',P+'n1_pulse_pressure_C3_functions.json.gz')


def recover(root):
    report_raw=(root/REPORT).read_bytes()
    if hashlib.sha256(report_raw).hexdigest()!=REPORT_SHA:raise ValueError('axis report hash changed')
    report=json.loads(report_raw)
    result_raw=base64.b64decode(report['artifacts'][RESULT]['base64'],validate=True)
    if hashlib.sha256(result_raw).hexdigest()!=RESULT_SHA:raise ValueError('axis result hash changed')
    result=json.loads(result_raw)
    if result.get('status')!='returned':raise ValueError('axis result not returned')
    for record in result['result']['binder_records']:
        if record['variable_descriptor']!={'operation':'bound_variable','name':'Z'} or record['value_descriptor']!={'operation':'exact_rational','numerator':0,'denominator':1}:
            raise ValueError('actual axis descriptor differs from exact Z=0')
    if {r['node'] for r in result['result']['binder_records']}!={31776,31944}:
        raise ValueError('axis binder set changed')
    raw=(root/ARCHIVE).read_bytes()
    if hashlib.sha256(raw).hexdigest()!=ARCHIVE_SHA:raise ValueError('raw archive hash changed before decoding')
    archive=json.loads(gzip.decompress(raw)); frame=archive['Rp_native_C5_pulse_source_bindings']
    mapping=frame['input_substitution']; current=frame['current_frame']; found=[]
    if mapping.get('31655')!=265772:raise ValueError('target mapping changed')
    for name in INPUTS:
        packet=frame['inputs'].get(name,{})
        for key, old in packet.get('native',{}).items():
            if not isinstance(old,list):raise ValueError('native port is not a list')
            positions=[i for i,ref in enumerate(old) if ref==31655]
            if not positions:continue
            target=current.get(key)
            if not isinstance(target,list):raise ValueError('current native port missing')
            rows=[]
            for i,ref in enumerate(old):
                if i>=len(target):break
                mapped=ref if ref in (0,1) else mapping.get(str(ref))
                rows.append({'order':i,'original_ref':ref,'mapping_ref':mapped,
                             'current_port_ref':target[i],
                             'direct_saved_mapping_matches':mapped==target[i]})
            found.append({'input':name,'key':key,'target_positions':positions,'rows':rows})
    # At the fixed axis, E0 and differentiation do not commute in general.
    z,t=s.symbols('z t'); f0,f1,f2,a,b,c=s.symbols('f0 f1 f2 a b c')
    f=f0+f1*z+f2*z*z/2; axis=f.subs(z,0)
    tests={'axis_value':s.expand(axis-f0)==0,
           'differentiate_then_axis_order1':s.expand(s.diff(f,z).subs(z,0)-f1)==0,
           'differentiate_then_axis_order2':s.expand(s.diff(f,z,2).subs(z,0)-f2)==0,
           'axis_then_differentiate_order1':s.diff(axis,z)==0,
           'axis_then_differentiate_order2':s.diff(axis,z,2)==0,
           'outer_parameter_first':s.expand(s.diff((t*t).subs(t,a+b*z+c*z*z/2),z).subs(z,0)-2*a*b)==0,
           'outer_parameter_second':s.expand(s.diff((t*t).subs(t,a+b*z+c*z*z/2),z,2).subs(z,0)-(2*b*b+2*a*c))==0}
    return {'schema':'frame_axis_native_port_v1','status':'port_roles_and_conditional_axis_jet',
            'target_mapping':[31655,265772],'native_port_roles':found,'port_role_found':bool(found),
            'axis_binder_refs':[31776,31944],'symbolic_checks':tests,'symbolic_checks_passed':all(tests.values()),
            'source_values_computed':False,'expression_edges_traversed':False,'geometry_closure_traversed':False,
            'actual_port_is_derivative_jet_proved':False,'ordinary_Z_transport_proved':False,
            'actual_axis_composition_identity_proved':False,'actual_capture_free_proved':False,
            'physical_endpoint_proved':False,'old_errors_removed':False,'runtime_model_metadata':None}
