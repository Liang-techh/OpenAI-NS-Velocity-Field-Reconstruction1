"""Exact common-source radius algebra before any numerical enclosure."""
import functools
import sympy as s
from lei_ren_part1_paper_compliant_current_physical_tensor_locator_operator import (
    V,SYMBOLS,DOMAINS,original_native_radius_recipes,original_radial_cover_theorem)
from lei_ren_part1_paper_compliant_current_tensor_registry import ROUTES,ADJACENT
from lei_ren_part1_paper_compliant_current_O3_transition_background_tensor import SourceAST

OFFSET=s.Symbol('requested_radius_offset',real=True)
LOG_HB=s.Symbol('exact_source_log_hb',real=True)


@functools.lru_cache(maxsize=1)
def original_equalities():
    a=SYMBOLS
    ref=s.log(110)+10*(a['log_C']+a['log_P'])
    pulse=ref+a['log_P']+1+a['Tw']
    return {a['log_Ra']:s.log(4)+a['log_epsilon'],a['h_switch']:a['h_bridge'],
        a['log_Rref']:ref,a['log_gap']:10*(a['log_C']+a['log_P'])-a['T'],
        a['log_Rp']:pulse,a['log_Rv']:pulse+13/a['mu']}


@functools.lru_cache(maxsize=8192)
def canonical_radius(expression):
    """Only admitted exact source equalities; no numeric midpoint/cap."""
    a=SYMBOLS
    expression=expression.xreplace(original_equalities()).subs(a['log_P'],s.exp(a['Md'])+11)
    # The symbol denotes the exact source, never its numerical cap. Keep
    # that defining exp(logh) visible in every canonical radius expression.
    expression=expression.subs(a['h_bridge'],s.exp(LOG_HB))
    expression=s.expand_log(expression,force=False)
    # Sympy can leave log(110) in one endpoint and log(10)+log(11)
    # in the other. Normalize the small exact rational logarithms in the
    # original source to prime logs; interval subtraction must never decide
    # this identity. Larger user constants may remain unsimplified, safely
    # retaining candidate ambiguity rather than attempting costly factoring.
    replacements={}
    for node in s.preorder_traversal(expression):
        if node.func is not s.log:continue
        argument=node.args[0]
        if not argument.is_Rational or argument<=0:continue
        numerator,denominator=int(argument.p),int(argument.q)
        if max(numerator,denominator)>1000000:continue
        replacements[node]=sum((power*s.log(prime) for prime,power in s.factorint(numerator).items()),s.Integer(0))-sum(
            (power*s.log(prime) for prime,power in s.factorint(denominator).items()),s.Integer(0))
    return s.expand(expression.xreplace(replacements))


@functools.lru_cache(maxsize=1)
def source_geometry():
    recipes,proof=original_native_radius_recipes();geometry={}
    for region,expr in recipes.items():
        lo,hi=map(s.sympify,DOMAINS[region]);original=canonical_radius(expr)
        geometry[region]=dict(original=original,lower=(-s.oo if region=='core_positive_radius' else canonical_radius(expr.subs(V,lo))),
            upper=canonical_radius(expr.subs(V,hi)),domain=(lo,hi),
            base=canonical_radius(expr.subs(V,0)) if region not in ('core_positive_radius','actual_patch','O2_axial') else None,
            slope=s.simplify(s.diff(original,V)))
    return geometry,proof


@functools.lru_cache(maxsize=1)
def anchors():
    a=SYMBOLS;r=a['log_Rref'];p=a['log_Rp'];v=a['log_Rv'];rel=v+100+a['Lrel']
    return dict(epsilon=a['log_epsilon'],Ra=a['log_Ra'],R100=s.log(100),R110=s.log(110),
        Rsh=s.log(110)+a['T'],Rref=r,Rm=r-6,Rh=r-5,Rd=r+a['log_P'],Rw=r+a['log_P']+1,
        Rp=p,Rv=v,Rf=v+100,Rrel=rel,Rs=rel+1,Rq=rel+2+a['Ts'],
        Rt=rel+2+a['Ts']+a['wait'],Rtail=rel+5+a['Ts']+a['wait'])


@functools.lru_cache(maxsize=1)
def boundary_sources():
    geometry,_=source_geometry();regions=list(geometry);result={}
    for j,(_,name,_,_) in enumerate(ADJACENT):
        left=regions[j];right=regions[j+1]
        if s.simplify(geometry[left]['upper']-geometry[right]['lower'])!=0:
            raise ValueError('Original exact source radius endpoint changed: '+name)
        result[name]=dict(expression=geometry[left]['upper'],kind='adjacent',left_region=left,right_region=right)
    for edge in (49,51,59,61,69,71):
        result['patch_support_'+str(edge)]=dict(expression=canonical_radius(SYMBOLS['log_Rref']-6+s.log(s.Rational(edge,40))),
            kind='internal_support',region='actual_patch',exact_native_coordinate=s.Rational(edge,40))
    for center in (-3,-1):
        for side in (-1,1):
            edge=s.Rational(center)+side*s.Rational(3,20)
            result['pulse_end_'+str(edge)]=dict(expression=canonical_radius(SYMBOLS['log_Rv']+edge),
                kind='internal_support',region='pulse_end',exact_native_coordinate=edge)
            name='outer_angular_'+str(center)+'_'+str(side)
            result[name]=dict(expression=canonical_radius(anchors()['Rrel']+edge),
                kind='internal_support',region='outer_angular',exact_native_coordinate=edge)
    if len(result)!=46:raise ValueError('Exactly 32 adjacent plus 14 support source radii required')
    return result


@functools.lru_cache(maxsize=8192)
def inverse_source_expression(region,requested_logR):
    geometry,_=source_geometry();row=geometry[region];a=SYMBOLS
    if region in ('core_positive_radius','actual_patch'):
        base=canonical_radius(original_native_radius_recipes()[0][region]-s.log(V))
        log_native=s.simplify(requested_logR-base)
        return s.simplify(s.exp(log_native)),log_native
    if region=='O2_axial':
        delta=s.simplify(requested_logR-canonical_radius(a['log_Rref']))
        return s.simplify(s.log(delta)/a['Md']),delta
    # Cancellation precedes all interval arithmetic. In particular HB and
    # 1/mu cancel in microscopic/pulse inverses of source-correlated inputs.
    return s.cancel((requested_logR-row['base'])/row['slope']),None


@functools.lru_cache(maxsize=1)
def correlated_source_theorem():
    previous=original_radial_cover_theorem();equalities=original_equalities()
    if {str(k):str(v) for k,v in equalities.items()}!=previous['source_parameter_equalities_required']:
        raise ValueError('Correlated canonicalization requires the admitted original equalities')
    geometry,_=source_geometry();boundaries=boundary_sources();asts=SourceAST()
    for stem,method in (('current_actual_patch_background_tensor','current_patch_source_theorem'),
        ('current_pulse_end_background_tensor','support_interface'),
        ('current_angular_internal_background_tensor','interface')):
        asts.method(stem,method)
    # Independently normalize the reciprocal gap derivative to the positive
    # source width specified by the original d coordinate.
    gap=s.simplify(geometry['pulse_gap_end']['slope']-(1-4*SYMBOLS['mu'])/SYMBOLS['mu'])
    if gap!=0:raise ValueError('Original reciprocal gap slope differs')
    return dict(original_exact_equalities_before_enclosure={str(k):str(v) for k,v in equalities.items()},
        exact_positive_width_source='h_bridge=h_switch=exp(exact_source_log_hb)>0',
        source_recipe_and_cap_enclosure_are_distinct=True,
        all_32_adjacent_source_radius_equalities=True,all_14_original_support_source_radii_available=True,
        reciprocal_gap_positive_slope_identity=True,
        regional_original_derivative_sources={region:str(row['slope']) for region,row in geometry.items()},
        required_positive_sources=previous['positive_source_requirements'],
        lower_core_endpoint='R=0 is the admitted nonsingular analytic core extension',
        upper_exterior_endpoint='heat_exterior offset tends to infinity; original Gamma T/E identically zero',
        physical_log_radius_recipe='(log2+lognu+2*actual_log_lambda+exact_source_logR)/2',
        request_semantics='Source-dependent physical point family; not a fixed absolute coordinate or a resolved velocity point',
        boundary_source_inventory=list(boundaries),
        input_hashes={**previous['input_hashes'],**asts.hashes},passed=True)
