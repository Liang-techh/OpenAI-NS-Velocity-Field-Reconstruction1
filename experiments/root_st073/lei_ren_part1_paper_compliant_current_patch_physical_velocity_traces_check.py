"""Independent physical derivatives, asymmetric units, and lineage guards."""
import json
from pathlib import Path
import sympy as s
import lei_ren_part1_paper_compliant_current_patch_physical_velocity_traces as source


def independent_derivatives():
    # An exactly differentiable physical field with all four nonzero
    # cylindrical sources. This tests moving basis, fixed radial units,
    # all35 spatial multiindices and the fixed-position time operator.
    # The fixture is an operator check, not an actual graph point sample.
    raw = {name:{'y%d_Z%d'%(j,n):
        (s.sqrt(s.Symbol('same_positive_patch_x',positive=True))*s.Rational(1,2)**j
         if name in source.LABELS[:1]+source.LABELS[2:3]
         else s.Symbol('same_positive_patch_x',positive=True)) if n==0 else s.Integer(0)
        for j in range(5) for n in range(5-j)} for name in source.LABELS}
    mapped = source.physical_rows(source.SourceAST(),'fixture',raw=raw)
    px,py,pz,t = s.symbols('physical_x physical_y physical_z physical_t',real=True)
    lam = s.sqrt(pz*pz+1-t)
    expected = {'ux':(px-3*py)/(2*lam**2),'uy':(py+3*px)/(2*lam**2),
        'uz':(px*px+py*py)/(4*lam**3),'p':9*(px*px+py*py)/(4*lam**4)}
    point = {px:s.Rational(3,4),py:s.Rational(1,3),pz:s.Rational(2,5),t:s.Rational(1,7)}
    l0 = lam.subs(point); radius = s.sqrt((px*px+py*py).subs(point))
    R = radius**2/(2*l0**2)
    substitution = {mapped['x']:R/2,mapped['Z']:point[pz]/l0,mapped['delta']:0,
        source.CS:point[px]/radius,source.SN:point[py]/radius}
    bases = (1,3,1,1,1,R,2)
    def evaluate(row):
        result = s.Integer(0)
        for powers, coefficient in source.merged_terms(row).items():
            amplitude = s.prod(s.sympify(base)**power for base,power in zip(bases,powers))
            result += coefficient.subs(substitution)*amplitude
        return result*l0**row['gamma'].subs(substitution)
    count = 0
    for i,j,b in source.INDICES:
        index = 'x%d_y%d_z%d'%(i,j,b)
        for component in source.COMPONENTS:
            prefix = index+'/'+component+'/'
            actual = sum((evaluate(row) for key,row in mapped['spatial'].items() if key.startswith(prefix)),s.Integer(0))
            independent = s.diff(expected[component],px,i,py,j,pz,b).subs(point)
            if s.simplify(actual-independent) != 0:
                raise ArithmeticError('Independent physical derivative failed: '+index+'/'+component)
            count += 1
    for component in source.COMPONENTS:
        actual = sum((evaluate(row) for key,row in mapped['time'].items() if key.startswith(component+'/')),s.Integer(0))
        if s.simplify(actual-s.diff(expected[component],t).subs(point)) != 0:
            raise ArithmeticError('Independent fixed-position time derivative failed: '+component)
        count += 1
    return count


def run():
    data = json.loads((source.HERE/source.NAME).read_bytes())
    for name,digest in data['input_hashes'].items():
        if source.sha(name) != digest: raise ValueError('Changed physical patch trace dependency: '+name)
    proof = source.theorem()
    if proof != data['exact_actual_patch_physical_velocity_trace_theorem']:
        raise ValueError('Current canonical physical transport theorem differs')
    left = source.physical_rows(source.SourceAST(),'left')
    negatives = {}
    for perturb in ('radial_unit','P0'):
        right = source.physical_rows(source.SourceAST(),'right',perturb)
        failed = [key for key,ok in source.compare(left,right).items() if not ok]
        if not failed: raise ArithmeticError('Asymmetric physical source/unit change was hidden: '+perturb)
        negatives[perturb] = dict(right_only_change_detected=True,physical_contributions_changed=len(failed))
    count = independent_derivatives()
    bound_available = proof['saved_numeric_bound_definition_verified']
    if not bound_available:
        try: source.patch_physical_trace(49)
        except ValueError as error:
            if 'different native dispatcher' not in str(error): raise
        else: raise ArithmeticError('Stale saved numeric bound definition was accepted')
    if any(data[key] is not False for key in ('current_modified_global_velocity_interfaces_certified',
        'common_N_modified_cones_energy_recursion_full_NS_certified')):
        raise ValueError('Local physical joins exceed completed scope')
    checks = proof['exact_linear_physical_source_join_identities']
    hashes = dict(data['input_hashes'])
    for name in (source.NAME,Path(__file__).name): hashes[name] = source.sha(name)
    receipt = dict(all_passed=True,source_family=proof['source_family'],
        exact_canonical_spatial_contribution_identities=sum(key.startswith('spatial/') for key in checks),
        exact_canonical_time_contribution_identities=sum(key.startswith('time/') for key in checks),
        six_actual_edges_bound_to_same_current_source_theorem=list(source.EDGES),
        independent_physical_operator_fixture_derivatives=count,
        asymmetric_right_physical_unit_and_pressure_negative_controls=negatives,
        stale_saved_numeric_bound_rejected=not bound_available,
        current_six_patch_support_physical_spatial4_time1_source_joins_certified=True,
        current_six_patch_support_physical_numeric_trace_bounds_available=bound_available,
        current_modified_global_velocity_interfaces_certified=False,
        common_N_modified_cones_energy_recursion_full_NS_certified=False,input_hashes=hashes)
    (source.HERE/source.RECEIPT).write_bytes((json.dumps(receipt,indent=2)+'\n').encode())
    print('Six physical patch joins PASS:',len(checks),'source contribution identities;',count,
        'independent physical derivatives; numeric bound admitted:',bound_available,flush=True)
    return receipt


if __name__ == '__main__': run()
