"""Bind saved axis-expression handles to saved Gamma history roles.

This reads authenticated literal function-binding records only. It does not
walk graph edges, evaluate functions, or rerun any ancestor constructor.
"""
import gzip
import hashlib
import json

P='lei_ren_part1_paper_compliant_current_original_'
RP='experiments/root_st073/'+P+'Rp_native_C5_pulse_source.json.gz'
GAMMA='experiments/root_st073/'+P+'n1_Gamma_future_functions.json.gz'
AXIS='tmp/frame_Z_composition_jet_v1/result.json'
PORT='tmp/frame_axis_native_port_v1/result.json'
PIN={RP:'442df70e5ab4051066b241e943c8d0531067d057e774a04dcc341dadbd703ebe',
     GAMMA:'dbc3dd702930b85599e7967b5e14d934f3c326a731cd1b9528990158ff74c612',
     AXIS:'9a45b8bc308e4ac81fffcc9887e2e1ed39c9897857cabb21032521ab0fe1e931'}


def read(root,name,gzip_encoded=False):
    raw=(root/name).read_bytes()
    if name in PIN and hashlib.sha256(raw).hexdigest()!=PIN[name]:raise ValueError('input hash changed: '+name)
    return json.loads(gzip.decompress(raw) if gzip_encoded else raw)


def recover(root):
    rp=read(root,RP,True); gamma=read(root,GAMMA,True); axis=read(root,AXIS); port=read(root,PORT)
    if axis.get('status')!='returned' or port.get('status')!='returned':raise ValueError('prior leaf did not return')
    histories=gamma['function_bindings']['histories']; native=gamma['function_bindings']['native_binding']
    records={r['node']:r for r in axis['result']['binder_records']}
    roles=[]
    for label,node in [('Xt',31776),('Xf',31944)]:
        jet=histories[label]
        value=jet[0] if isinstance(jet,list) and jet else None
        roles.append({'history':label,'history_jet':jet,'axis_node':node,
                      'saved_axis_expression_ref':records[node]['expression_ref'],
                      'same_saved_reference':value==records[node]['expression_ref']})
    gamma_file=P+'n1_Gamma_future_functions.json.gz'
    return {'schema':'frame_axis_gamma_roles_v1','status':'saved_role_binding_read',
            'source_family_matches':gamma.get('source_family')==rp.get('source_family'),
            'rp_input_gamma_pin':rp.get('input_hashes',{}).get(gamma_file),
            'gamma_raw_sha256':PIN[GAMMA],'axis_roles':roles,'gamma_native_binding':native,
            'axis_port_roles':port['result']['native_port_roles'],
            'expression_edges_traversed':False,'geometry_closure_traversed':False,
            'original_builders_called':False,'source_values_computed':False,
            'same_namespace_semantics_proved':False,'actual_derivative_transport_proved':False,
            'actual_axis_composition_identity_proved':False,'physical_endpoint_proved':False,
            'old_errors_removed':False,'runtime_model_metadata':None}
