"""Check fresh edge maps, domains, exact layouts and current definitions."""
import gzip
import json
from pathlib import Path
import mpmath as mp
import lei_ren_part1_paper_compliant_current_patch_physical_numeric_bounds as source


def validate_view(view,edge,definition,ctx=None):
    t = source.transport
    if view['modified_velocity_pressure_dispatch_definition_sha256'] != definition:
        raise ValueError('Foreign dispatcher definition in fresh patch bound')
    if view['exact_support_x'] != str(t.s.Rational(edge,40)):
        raise ValueError('Wrong actual support edge')
    if set(view['physical_spatial_cartesian_mixed4']) != {'x%d_y%d_z%d'%v for v in t.INDICES}:
        raise ValueError('Incomplete physical spatial index layout')
    spatial = view['physical_spatial_cartesian_mixed4'];time = view['first_fixed_x_physical_time_derivative']
    counts = 0
    for components in spatial.values():
        if set(components) != set(t.COMPONENTS): raise ValueError('Incomplete vector/pressure layout')
        for component,labels in components.items():
            expected = {t.UR,t.UT} if component in ('ux','uy') else {t.UZ} if component=='uz' else {t.P}
            if set(labels) != expected: raise ValueError('Lost original signed source sector')
            counts += len(labels)
    counts += sum(len(labels) for labels in time.values())
    if counts != 216: raise ValueError('Expected210 spatial plus6 time groups')
    for groups in ([labels for components in spatial.values() for labels in components.values()], [labels for labels in time.values()]):
        for labels in groups:
            for row in labels.values():
                if row['exact_zero']:
                    if row['terms'] or row['log_absolute_upper'] is not None: raise ValueError('Invalid exact-zero layout')
                elif not row['terms'] or row['log_absolute_upper'] is None:
                    raise ValueError('Missing signed source terms or finite log bound')
                elif ctx is not None and any(not mp.isfinite(value) for value in
                    t.endpoints(t.read_interval(ctx,row['log_absolute_upper']))):
                    raise ValueError('Nonfinite logarithmic derivative upper bound')
    return counts


def run():
    data = json.loads((source.HERE/source.NAME).read_bytes())
    for name,digest in data['input_hashes'].items():
        if source.sha(name) != digest: raise ValueError('Changed fresh patch bound input: '+name)
    current = source.inputs()
    if source.definition(current) != data['fresh_patch_physical_bound_definition_sha256']:
        raise ValueError('Changed primitive/scalar/map bound definition')
    views = json.loads(gzip.decompress((source.HERE/source.VIEWS_NAME).read_bytes()))
    if tuple(map(int,views)) != source.transport.EDGES: raise ValueError('Incomplete six-support export')
    total = 0
    for edge in source.transport.EDGES:
        view = views[str(edge)]
        total += validate_view(view,edge,current['accepted']['current_dispatch_definition'],current['ctx'])
        for key,value in current['accepted']['source_family'].items():
            if view[key] != value: raise ValueError('Foreign actual primitive source family')
        for key,expected in (('Z',(-1,1)),('requested_log_tau',(-3,-1)),('cos_theta',(-1,1)),('sin_theta',(-1,1))):
            if source.transport.endpoints(source.transport.read_interval(current['ctx'],view[key])) != expected:
                raise ValueError('Incorrect quantitative bound domain')
    # Independently regenerate both terminal support bounds from the raw
    # checked packets. In particular the last one retains its actual P0
    # and functional terminal closure rather than a fabricated zero field.
    for edge in (source.transport.EDGES[0],source.transport.EDGES[-1]):
        fresh = source.encode(source.pack(source.build_view(current,edge)))
        if fresh != views[str(edge)]: raise ValueError('Fresh edge map differs from encoded source/scalar transport')
    bad = dict(views['49'],modified_velocity_pressure_dispatch_definition_sha256=current['accepted']['saved_dispatch_definition'])
    try: validate_view(bad,49,current['accepted']['current_dispatch_definition'])
    except ValueError as error:
        if 'Foreign dispatcher' not in str(error): raise
    else: raise ArithmeticError('Old foreign embedded definition was accepted')
    if data['old_full33_numeric_view_used_as_numerical_input'] is not False or any(data[key] is not False for key in
        ('current_modified_global_velocity_interfaces_certified','common_N_modified_cones_energy_recursion_full_NS_certified')):
        raise ValueError('Fresh six-support bounds exceed completed scope')
    hashes = dict(data['input_hashes'])
    for name in (source.NAME,Path(__file__).name): hashes[name] = source.sha(name)
    receipt = dict(all_passed=True,source_family=current['accepted']['source_family'],
        fresh_patch_physical_bound_definition_sha256=source.definition(current),
        exact_support_edges=list(source.transport.EDGES),physical_contribution_groups=total,
        independently_regenerated_terminal_edge_maps=2,
        old_foreign_dispatch_definition_rejected=True,
        current_six_patch_support_physical_numeric_trace_bounds_available=True,
        old_full33_numeric_view_used_as_numerical_input=False,
        current_modified_global_velocity_interfaces_certified=False,
        common_N_modified_cones_energy_recursion_full_NS_certified=False,input_hashes=hashes)
    (source.HERE/source.RECEIPT).write_bytes((json.dumps(receipt,indent=2)+'\n').encode())
    print('Fresh six physical patch bounds PASS:',total,'signed contribution groups; current definition bound',flush=True)
    return receipt


if __name__ == '__main__': run()
