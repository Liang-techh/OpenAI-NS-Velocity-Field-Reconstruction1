"""Same-source core-inlet plus separate Cauchy phase-variation enclosures.

The source owner proves the input core jets and collar radial bounds. This
leaf uses FTC and Cauchy's estimate to enclose actual Phi and U on a first
collar phase interval, retaining hb as a formal factor. The variation is a
symmetric error radius, never an evaluated field or a signed point value.
"""
from __future__ import annotations
import math


def _radius(c,upper):
    upper=c.mpf(upper).b
    if upper<0:raise ValueError('nonnegative error upper required')
    return c.make_mpf(((-upper)._mpi_[0],upper._mpi_[1]))


def build(c, *, sc, phase, core_inlet_Z, radial_tube_bounds, cauchy_radius,
          rho_upper, source_family, source_atoms, admitted_source_contract):
    """Enclose actual Z0..Z5 with core data plus explicit hb variation errors.

    core_inlet_Z[field][k] must cover the unique source at rho=4 and all real
    Z in the chosen interval. radial_tube_bounds[field] maps formal hb powers
    to positive upper intervals for |d_rho(field)| on complex disks of radius
    at least cauchy_radius around every such Z. The caller must authenticate
    this analytic contract and the same source; this leaf cannot establish it.
    """
    sc,phase,rho_upper,radius=(c.mpf(x) for x in (sc,phase,rho_upper,cauchy_radius))
    if sc.a<=0 or sc.b>=c.mpf(1)/1024 or phase.a<0 or phase.b>(sc/4).a:
        raise ValueError('original microscopic collar phase required')
    if radius.a<=0 or rho_upper.a<4:raise ValueError('positive complex radius and rho upper required')
    if source_family is None or set(source_atoms)!=set(('hb','F0_anchor','Pstar_squared','Lambda')):
        raise ValueError('same-source atom records required')
    keys=('source_family_checked','core_enclosure_contains_rho4','complex_radial_bound_covers_collar',
          'analytic_Z_disks_cover_real_domain','hb_Z_independent','rho_upper_covers_collar')
    if any(admitted_source_contract.get(k) is not True for k in keys):raise ValueError('caller analytic source contract required')
    if set(core_inlet_Z)!=set(('Phi','U')) or set(radial_tube_bounds)!=set(('Phi','U')):
        raise ValueError('Phi and U source data required')
    out={}
    for field in ('Phi','U'):
        if set(core_inlet_Z[field])!=set(range(6)) or not radial_tube_bounds[field]:
            raise ValueError('core Z0..Z5 and radial source bounds required')
        bounds=radial_tube_bounds[field]
        if any(not isinstance(p,int) or isinstance(p,bool) or p<0 for p in bounds):
            raise ValueError('nonnegative collar radial-bound hb powers required')
        out[field]={}
        for k in range(6):
            base=c.mpf(core_inlet_Z[field][k])
            terms=[{'powers':(0,0,0,0),'coefficient':base,'origin':('core_inlet_source_enclosure',field,k),
                    'semantics':'actual_core_interval_not_selected_value'}]
            for power,bound in sorted(bounds.items()):
                bound=c.mpf(bound)
                if bound.a<0:raise ValueError('radial tube bound must be nonnegative')
                cap=(rho_upper.b*phase.b*bound.b*math.factorial(k)/(radius.a**k)).b
                terms.append({'powers':(power+1,0,0,0),'coefficient':_radius(c,cap),
                  'origin':('cauchy_phase_variation_error',field,k,power),
                  'semantics':'symmetric_certified_error_not_field_value'})
            out[field][k]=terms
    return {'source_family':source_family,'source_atoms':source_atoms,'fields':out,
      'ordinary_Z_orders':list(range(6)),'phase_interval':phase,
      'same_source_contract_independently_proved':False,
      'source_semantics':'core inlet interval plus separate formal-hb uncertainty',
      'analytic_estimate':'|d_Z^k (field(s)-field(0))| <= hb*s*rho_upper*k!*B/radius^k',
      'no_point_value_or_sign_claim':True,'allword_tail_added':False,'m15_computed':False}
