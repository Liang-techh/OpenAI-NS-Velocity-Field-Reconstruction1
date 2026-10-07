"""Transport six proved patch germs through the actual physical maps.

No defining graph is rebuilt. Function equality comes from the accepted
two-sided source theorem. Saved whole-patch bounds are admitted only with
the same dispatcher definition and on their explicit time sector.
"""
import gzip
import json
from pathlib import Path
from types import SimpleNamespace
import mpmath as mp
import sympy as s
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_current_patch_support_velocity_pressure import (
    HERE, PREFIX, EDGES, LABELS, sha, SourceAST)
from lei_ren_part1_paper_compliant_cartesian_field import (
    INDICES, COMPONENTS, CS, SN, cartesian_templates, angular_polynomial)
from lei_ren_part1_paper_compliant_global_physical_assembly import BASES
from lei_ren_part1_paper_compliant_pulse_physical_bounds import (
    UZ, UT, UR, P, physical_operators, interval_expression)
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

NAME = PREFIX+'current_patch_physical_velocity_traces.json'
RECEIPT = PREFIX+'current_patch_physical_velocity_traces_check.json'
SOURCE_NAME = PREFIX+'current_patch_support_velocity_pressure.json'
SOURCE_RECEIPT = PREFIX+'current_patch_support_velocity_pressure_check.json'
DISPATCH_RECEIPT = PREFIX+'current_modified_velocity_pressure_dispatch_check.json'
VIEWS_NAME = PREFIX+'current_modified_velocity_pressure_dispatch_views.json.gz'
PHYSICAL_RECEIPT = PREFIX+'current_modified_physical_velocity_check.json'
COVER_NAME = PREFIX+'current_global_tensor_cover.json'
RAW_NAMES = {UT:LABELS[0], UZ:LABELS[1], UR:LABELS[2], P:LABELS[3]}


def _accepted(require_bound_definition=True):
    records = {name:json.loads((HERE/name).read_bytes()) for name in
        (SOURCE_NAME, SOURCE_RECEIPT, DISPATCH_RECEIPT, PHYSICAL_RECEIPT, COVER_NAME)}
    source, dispatch, physical = (records[name] for name in
        (SOURCE_RECEIPT, DISPATCH_RECEIPT, PHYSICAL_RECEIPT))
    if not source['all_passed'] or not source['current_six_patch_support_velocity_pressure_mixed4_source_joins_certified']:
        raise ValueError('Checked independent actual patch source germs required')
    if not dispatch['all_passed'] or not physical['all_passed']:
        raise ValueError('Checked canonical native and physical velocity maps required')
    hashes = {}
    for name, record in records.items():
        for path, digest in record['input_hashes'].items():
            if path in hashes and hashes[path] != digest:
                raise ValueError('Conflicting current source definition: '+path)
            hashes[path] = digest
        hashes[name] = sha(name)
    for name, digest in hashes.items():
        if sha(name) != digest:
            raise ValueError('Changed accepted physical/source dependency: '+name)
    for key, value in source['source_family'].items():
        if dispatch.get(key) != value or physical.get(key) != value or records[COVER_NAME].get(key) != value:
            raise ValueError('Different source family, implicit controls or P0: '+key)
    if physical['input_hashes'].get(DISPATCH_RECEIPT) != sha(DISPATCH_RECEIPT):
        raise ValueError('Physical field does not consume the same checked native dispatcher')
    if dispatch['input_hashes'].get(VIEWS_NAME) != sha(VIEWS_NAME):
        raise ValueError('Saved native views are not bound by the accepted dispatcher')
    views = json.loads(gzip.decompress((HERE/VIEWS_NAME).read_bytes()))
    native = views['actual_patch']
    if len(native['source_pieces']) != 1 or native['source_pieces'][0]['source_branch'] != 'unchanged_same_checked_canonical_route':
        raise ValueError('One unchanged actual patch function required')
    view = native['source_pieces'][0]['complete_velocity_pressure_source']
    if view['chart'] != 'actual_patch' or view['physical_layout'] != 'four_label_Cartesian':
        raise ValueError('Saved view is not the actual patch Cartesian map')
    bound_definition_matches = native['modified_velocity_pressure_dispatch_definition_sha256'] == dispatch['modified_velocity_pressure_dispatch_definition_sha256']
    if require_bound_definition and not bound_definition_matches:
        raise ValueError('Saved view belongs to a different native dispatcher')
    for key, value in source['source_family'].items():
        if view.get(key) != value:
            raise ValueError('Saved patch bound has a foreign source family: '+key)
    c = MPIntervalContext(); c.dps = 240
    xlo, xhi = endpoints(read_interval(c, native['requested_coordinate']))
    if any(not xlo <= mp.mpf(edge)/40 <= xhi for edge in EDGES):
        raise ValueError('Saved whole-patch domain does not cover the six supports')
    if endpoints(read_interval(c, view['Z'])) != (-1,1) or endpoints(read_interval(c, view['requested_log_tau'])) != (-3,-1):
        raise ValueError('Saved patch source/time sector changed')
    if any(endpoints(read_interval(c, view[key])) != (-1,1) for key in ('cos_theta','sin_theta')):
        raise ValueError('Whole angular domain required')
    hypotheses = records[COVER_NAME]['actual_current_source_cover_hypotheses']
    dl, dh = endpoints(read_interval(c, hypotheses['source_delta']))
    if not hypotheses['passed'] or dl < 0 or dh >= 1:
        raise ValueError('Same actual nonsingular 0<=delta<1 hypothesis required')
    return dict(records=records, candidate_saved_view=view, input_hashes=hashes,
        saved_bound_definition_verified=bound_definition_matches,
        saved_dispatch_definition=native['modified_velocity_pressure_dispatch_definition_sha256'],
        current_dispatch_definition=dispatch['modified_velocity_pressure_dispatch_definition_sha256'],
        source_family=source['source_family'], actual_delta=hypotheses['source_delta'])


def exact_context():
    def exact(value):
        return s.Rational(str(value)) if isinstance(value,(str,int,float)) else value
    return SimpleNamespace(mpf=exact, sqrt=s.sqrt, ln=s.log, exp=s.exp)


def programs(asts):
    c = exact_context()
    env = dict(mp=c, UZ=UZ, UT=UT, UR=UR, P=P, BASES=BASES,
        MICRO=(), COMMON=(), PREPULSE=(), PULSE=(), POST=(),
        cartesian_templates=cartesian_templates, angular_polynomial=angular_polynomial,
        physical_operators=physical_operators, interval_expression=interval_expression)
    for name in ('zero_powers','shift','ordinary_terms','normalized_sources',
                 'cartesian_source_row','time_source_row'):
        asts.replay('global_physical_assembly',name,env)
    for stem, names in (
        ('cartesian_field',('transverse_step','cartesian_templates','angular_polynomial')),
        ('pulse_physical_bounds',('physical_operators','interval_expression'))):
        for name in names: asts.method(stem,name)
    # The live native mapper calls this same original base implementation.
    asts.method('current_modified_velocity_pressure_dispatch_operator','canonical_original_velocity')
    asts.expression('current_modified_velocity_pressure_dispatch_operator','canonical_original_velocity','physical',
        wanted='field.geometry.evaluate(chart,Z,value,log_tau=log_tau,theta=theta,axis=axis)')
    asts.method('global_physical_assembly','evaluate')
    return c, env


def physical_rows(asts, side, perturb=None, raw=None):
    """Replay fixed-unit normalization and exact Cartesian/time operators."""
    c, env = programs(asts)
    x = s.Symbol('same_positive_patch_x',positive=True)
    z, delta = s.symbols('same_Z same_delta',real=True)
    logp, logr = s.symbols('same_logPstar same_logR',real=True)
    if raw is None:
        raw = {name:{'y%d_Z%d'%(j,n):s.Symbol('%s_%s_y%d_Z%d'%(side,name,j,n),real=True)
                    for j in range(5) for n in range(5-j)} for name in LABELS}
    source_atoms = {name:dict(rows) for name,rows in raw.items()}
    if perturb == 'radial_unit':
        raw = {name:dict(rows) for name,rows in raw.items()}
        raw[LABELS[2]] = {key:value*s.sqrt(x) for key,value in raw[LABELS[2]].items()}
    if perturb == 'P0':
        raw = {name:dict(rows) for name,rows in raw.items()}
        raw[LABELS[3]]['y0_Z0'] += 1
    field = SimpleNamespace(ctx=c,logP=logp)
    grids, logs, amplitudes = env['normalized_sources'](field,'actual_patch',z,x,
        {'physical_velocity_pressure_y_Z_mixed4':raw},None,logr)
    spatial, time = {}, {}
    for i,j,b in INDICES:
        key = 'x%d_y%d_z%d'%(i,j,b)
        for component in COMPONENTS:
            for label, (terms,gamma) in env['cartesian_source_row'](c,grids,component,i,j,b,z,delta,CS,SN,amplitudes).items():
                spatial[key+'/'+component+'/'+label] = dict(terms=terms,gamma=gamma)
    basis = {'ux':{UR:CS,UT:-SN},'uy':{UR:SN,UT:CS},'uz':{UZ:1},'p':{P:1}}
    for component, labels in basis.items():
        for label, angular in labels.items():
            terms, gamma = env['time_source_row'](c,grids,label,z,delta,amplitudes)
            time[component+'/'+label] = dict(terms=[(powers,angular*value) for powers,value in terms],gamma=gamma)
    return dict(raw=raw,source_atoms=source_atoms,grids=grids,amplitudes=amplitudes,spatial=spatial,time=time,x=x,Z=z,delta=delta)


def merged_terms(row, substitution=None):
    out = {}
    for powers, value in row['terms']:
        if substitution: value = value.xreplace(substitution)
        out[powers] = out.get(powers,s.Integer(0))+value
    return {powers:s.cancel(value) for powers,value in out.items() if value != 0}


def compare(left,right):
    substitution = {right['source_atoms'][name][key]:value for name,rows in left['source_atoms'].items() for key,value in rows.items()}
    checks = {}
    for group in ('spatial','time'):
        if left[group].keys() != right[group].keys(): raise ValueError('Physical layout differs')
        for key,row in left[group].items():
            other = right[group][key]
            a,b = merged_terms(row),merged_terms(other,substitution)
            ok = row['gamma'] == other['gamma'] and all(s.cancel(a.get(power,0)-b.get(power,0)) == 0 for power in a.keys()|b.keys())
            checks[group+'/'+key] = ok
    return checks


def theorem():
    accepted = _accepted(require_bound_definition=False); asts = SourceAST()
    left,right = physical_rows(asts,'left'),physical_rows(asts,'right')
    checks = compare(left,right)
    if not all(checks.values()): raise ArithmeticError('Canonical physical pullbacks do not preserve source joins')
    # Retain the distinct lambda exponents and fixed basepoint units.
    for label, name in RAW_NAMES.items():
        for j in range(5):
            for n in range(5-j):
                value = left['grids'][label][j,n][0][1]
                expected = left['raw'][name]['y%d_Z%d'%(j,n)]/(s.sqrt(left['x']) if label==UR else 1)
                if s.cancel(value-expected) != 0: raise ArithmeticError('Canonical patch normalization differs')
    # Inventory comparison is diagnostic only if the saved definition is
    # stale. No numeric bound is admitted until exact lineage is supplied.
    spatial = accepted['candidate_saved_view']['physical_spatial_cartesian_mixed4']
    time = accepted['candidate_saved_view']['first_fixed_x_physical_time_derivative']
    saved_keys = {'spatial/'+index+'/'+component+'/'+label
        for index,components in spatial.items() for component,labels in components.items() for label in labels}
    saved_keys |= {'time/'+component+'/'+label for component,labels in time.items() for label in labels}
    if saved_keys != set(checks): raise ValueError('Saved physical bound layout differs from the actual operators')
    hashes = {**accepted['input_hashes'],**asts.hashes,Path(__file__).name:sha(Path(__file__).name)}
    return dict(source_family=accepted['source_family'], source_bindings=asts.bindings,
        exact_linear_physical_source_join_identities=checks,
        six_edges={str(edge):dict(exact_x=str(s.Rational(edge,40)),
            consumed_two_sided_source_theorem=SOURCE_NAME,
            spatial4_time1_source_identities_transport_by_same_linear_operator=True,
            candidate_common_bound_native_view=VIEWS_NAME+'#actual_patch/source_pieces/0/complete_velocity_pressure_source',
            numeric_common_bound_admitted=accepted['saved_bound_definition_verified']) for edge in EDGES},
        canonical_basepoint_radial_unit='raw Ur/sqrt(Rm/2) derivative rows divided by sqrt(x); unit is not differentiated again',
        same_actual_radius='R=Rm*x=Rref*exp(-6)*x',
        nonsingular_hypotheses='Rm>0; Pstar>0; lambda>0; 0<=delta<1; |Z|<1; L=1-delta*Z^2>0',
        actual_delta=accepted['actual_delta'],
        physical_coordinates='r=sqrt(nu)*lambda*sqrt(2*Rm*x); z=sqrt(nu)*Z*lambda^(1-delta); tau=lambda^2*(1-Z^2)',
        constant_viscosity_spatial_power='nu^((1-order)/2) for velocity; nu^((2-order)/2) for pressure',
        constant_viscosity_time_power='sqrt(nu) for velocity; nu for pressure',
        common_saved_bound_domain=dict(x='[1,e]',Z='[-1,1] source closure; finite physical points require |Z|<1',
            log_tau='[-3,-1]',theta='all angles',viscosity='1; positive constant nu supplied as a separate exact scale'),
        saved_numeric_bound_definition_verified=accepted['saved_bound_definition_verified'],
        saved_dispatch_definition=accepted['saved_dispatch_definition'],
        current_dispatch_definition=accepted['current_dispatch_definition'],
        function_equality_proved_before_common_upper_bounds=True,
        global_velocity_interfaces_energy_common_N_cones_recursion_full_NS_remain_open=True,
        passed=True,input_hashes=hashes)


def patch_physical_trace(edge):
    """Common factored bounds for both proved traces on the saved sector."""
    if edge not in EDGES: raise ValueError('One of the six actual support numerators required')
    accepted = _accepted()
    return dict(exact_x=str(s.Rational(edge,40)),same_two_sided_source_theorem=SOURCE_NAME,
        common_physical_velocity_pressure_bound=accepted['candidate_saved_view'],
        copied_bound_is_not_the_function_equality_proof=True,
        common_bound_is_conservative_whole_patch_not_a_resolved_edge_value=True,
        log_tau_domain='[-3,-1]',finite_physical_Z_domain='(-1,1)',
        all_global_velocity_interfaces_certified=False)


def run():
    proof = theorem()
    result = dict(exact_actual_patch_physical_velocity_trace_theorem=proof,
        current_six_patch_support_physical_spatial4_time1_source_joins_certified=True,
        current_six_patch_support_physical_numeric_trace_bounds_available=proof['saved_numeric_bound_definition_verified'],
        current_modified_global_velocity_interfaces_certified=False,
        common_N_modified_cones_energy_recursion_full_NS_certified=False,input_hashes=proof['input_hashes'])
    (HERE/NAME).write_bytes((json.dumps(result,indent=2)+'\n').encode())
    print('Six patch source joins transported to physical spatial4/time1; numeric bound admitted:',
        proof['saved_numeric_bound_definition_verified'],flush=True)
    return result


if __name__ == '__main__': run()
