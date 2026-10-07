"""Fresh six-edge physical bounds from checked primitive patch packets.

The old full33 view is not a numerical input. All signed normalization and
physical rows are regenerated from the current encoded primitive grids.
"""
import gzip
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_patch_physical_velocity_traces as transport
from lei_ren_part1_paper_compliant_global_physical_assembly import (
    CompliantGlobalPhysicalAssembly as BASE,cartesian_source_row,time_source_row,log_row,BASES)
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE,PREFIX,sha = transport.HERE,transport.PREFIX,transport.sha
NAME = PREFIX+'current_patch_physical_numeric_bounds.json'
RECEIPT = PREFIX+'current_patch_physical_numeric_bounds_check.json'
VIEWS_NAME = PREFIX+'current_patch_physical_numeric_bounds_views.json.gz'
PACKET_NAME = PREFIX+'actual_feedback_patch_mixed_C4.json'
NORM_NAME = PREFIX+'physical_norm_family.json'
PRESSURE_NAME = PREFIX+'pressure_source.json'


def inputs():
    accepted = transport._accepted(require_bound_definition=False)
    receipt = json.loads((HERE/transport.RECEIPT).read_bytes())
    if not receipt['all_passed'] or not receipt['current_six_patch_support_physical_spatial4_time1_source_joins_certified']:
        raise ValueError('Checked physical two-sided source theorem required')
    hashes = {**accepted['input_hashes'],**receipt['input_hashes'],transport.RECEIPT:sha(transport.RECEIPT)}
    for name,digest in hashes.items():
        if sha(name) != digest: raise ValueError('Changed physical/primitive source: '+name)
    data = {name:json.loads((HERE/name).read_bytes()) for name in (PACKET_NAME,NORM_NAME,PRESSURE_NAME)}
    for name in data:
        if hashes.get(name) != sha(name): raise ValueError('Primitive/scalar receipt is not bound by the current graph: '+name)
    primitive, norm, pressure = data[PACKET_NAME],data[NORM_NAME],data[PRESSURE_NAME]['compliant_source']
    for key,value in accepted['source_family'].items():
        if primitive.get(key) != value: raise ValueError('Foreign primitive patch source: '+key)
    if norm['implicit_source_sha256'] != accepted['source_family']['implicit_source_sha256'] or pressure['implicit_source_sha256'] != norm['implicit_source_sha256']:
        raise ValueError('Scalar source differs from the actual implicit patch family')
    parameters = pressure['implicit_source_definition']
    if parameters['Md'] != '40' or parameters['logPstar'] != 'exp(Md)+11':
        raise ValueError('Original Pstar defining scalar changed')
    if not norm['full_physical_C3_K_norms_certified_for_uniform_analytic_family']:
        raise ValueError('Admitted actual selected Cstar scalar required')
    c = MPIntervalContext();c.dps = 240
    logp = c.exp(c.mpf(parameters['Md']))+11
    logc = transport.read_interval(c,norm['selected_logCstar'])
    logrref = c.ln(110)+10*(logc+logp)
    delta = transport.read_interval(c,accepted['actual_delta'])
    asts = transport.SourceAST()
    for stem,method in (('global_physical_assembly','normalized_sources'),
                        ('global_physical_assembly','radius'),
                        ('core_physical_field','__init__')):
        asts.method(stem,method)
    asts.expression('global_physical_assembly','__init__','self.logRref',
        wanted='c.ln(110)+10*(self.logC+self.logP)')
    asts.expression('core_physical_field','__init__','self.logC',wanted="read(norm,'selected_logCstar')")
    field = SimpleNamespace(ctx=c,logP=logp)
    return dict(accepted=accepted,primitive=primitive,ctx=c,field=field,
        delta=delta,logP=logp,logRref=logrref,source_bindings=asts.bindings,
        hashes={**hashes,**asts.hashes,Path(__file__).name:sha(Path(__file__).name)})


def build_view(source,edge):
    c = source['ctx'];index = transport.EDGES.index(edge)
    packet = source['primitive']['source_bump_edge_packets'][index]
    x = c.mpf(edge)/40
    xl,xh = transport.endpoints(transport.read_interval(c,packet['x']))
    if not xl <= transport.endpoints(x)[0] <= transport.endpoints(x)[1] <= xh:
        raise ValueError('Primitive edge enclosure does not contain the exact rational edge')
    Z = transport.read_interval(c,packet['Z'])
    if transport.endpoints(Z) != (-1,1): raise ValueError('Whole axial source domain required')
    if not packet['same_actual_coefficient_family_and_P0_retained'] or not packet['ordinary_derivatives_not_Taylor_radial_coefficients']:
        raise ValueError('Actual source and ordinary derivative units required')
    raw = {label:{key:transport.read_interval(c,value) for key,value in rows.items()}
        for label,rows in packet['physical_velocity_pressure_y_Z_mixed4'].items()}
    if set(raw) != set(transport.LABELS) or any(len(rows) != 15 for rows in raw.values()):
        raise ValueError('Complete four-label mixed4 primitive rows required')
    logr = source['logRref']-6+c.ln(x)
    grids,logs,amplitudes = BASE.normalized_sources(source['field'],'actual_patch',Z,x,
        {'physical_velocity_pressure_y_Z_mixed4':raw},None,logr)
    lt = c.mpf([-3,-1]);angular = c.mpf([-1,1]);spatial = {};time = {}
    for i,j,b in transport.INDICES:
        key = 'x%d_y%d_z%d'%(i,j,b);spatial[key] = {}
        for component in transport.COMPONENTS:
            parts = cartesian_source_row(c,grids,component,i,j,b,Z,source['delta'],angular,angular,amplitudes)
            spatial[key][component] = {label:log_row(c,rows,logs,gamma,lt/2) for label,(rows,gamma) in parts.items()}
    basis = {'ux':{transport.UR:angular,transport.UT:-angular},
        'uy':{transport.UR:angular,transport.UT:angular},'uz':{transport.UZ:c.mpf(1)},'p':{transport.P:c.mpf(1)}}
    for component,labels in basis.items():
        time[component] = {}
        for label,factor in labels.items():
            terms,gamma = time_source_row(c,grids,label,Z,source['delta'],amplitudes)
            time[component][label] = log_row(c,[(powers,value*factor) for powers,value in terms],logs,gamma,lt/2)
    return dict(chart='actual_patch',exact_support_x=str(transport.s.Rational(edge,40)),
        Z=Z,requested_log_tau=lt,cos_theta=angular,sin_theta=angular,
        physical_spatial_cartesian_mixed4=spatial,first_fixed_x_physical_time_derivative=time,
        positive_source_log_bases=dict(zip(BASES,logs)),source_logR_enclosure=logr,
        physical_layout='four_label_Cartesian',
        **source['accepted']['source_family'],
        modified_velocity_pressure_dispatch_definition_sha256=source['accepted']['current_dispatch_definition'],
        primitive_source_packet=PACKET_NAME+'#/source_bump_edge_packets/'+str(index),
        common_physical_source_join_receipt=transport.RECEIPT,
        same_bound_applies_to_both_proved_source_traces=True,
        exact_zero_and_signed_source_sectors_retained=True,
        normalization_not_differentiated_twice=True,
        physical_log_bounds_use_lambda_ge_sqrt_tau=True,
        actual_lambda_not_identified_with_sqrt_tau=True,
        positive_source_exponentials_not_materialized=True,
        finite_physical_coordinate_domain='tau>0, |Z|<1, all angles; log_tau[-3,-1]; viscosity1',
        numeric_edge_field_values_resolved=False,all_global_velocity_interfaces_certified=False)


def definition(source):
    terms = dict(current_dispatch_definition=source['accepted']['current_dispatch_definition'],
        source_join_receipt=sha(transport.RECEIPT),primitive_packets=sha(PACKET_NAME),
        norm_scalars=sha(NORM_NAME),pressure_parameters=sha(PRESSURE_NAME),
        source=sha(Path(__file__).name))
    return hashlib.sha256(json.dumps(terms,sort_keys=True).encode()).hexdigest()


def run():
    source = inputs();views = {}
    for edge in transport.EDGES:
        views[str(edge)] = build_view(source,edge)
        print('Fresh primitive-source physical patch bound:',edge,'/40',flush=True)
    payload = (json.dumps(encode(pack(views)),separators=(',',':'))+'\n').encode()
    (HERE/VIEWS_NAME).write_bytes(gzip.compress(payload,compresslevel=6,mtime=0))
    hashes = {**source['hashes'],VIEWS_NAME:sha(VIEWS_NAME)}
    result = dict(source_family=source['accepted']['source_family'],
        fresh_patch_physical_bound_definition_sha256=definition(source),
        current_dispatch_definition=source['accepted']['current_dispatch_definition'],
        primitive_source=PACKET_NAME,exact_support_edges=list(transport.EDGES),complete_fresh_views=VIEWS_NAME,
        source_bindings=source['source_bindings'],
        current_six_patch_support_physical_numeric_trace_bounds_available=True,
        old_full33_numeric_view_used_as_numerical_input=False,
        current_modified_global_velocity_interfaces_certified=False,
        common_N_modified_cones_energy_recursion_full_NS_certified=False,input_hashes=hashes)
    (HERE/NAME).write_bytes((json.dumps(result,indent=2)+'\n').encode())
    return result


if __name__ == '__main__': run()
