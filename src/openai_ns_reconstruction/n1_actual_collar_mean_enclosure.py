"""Actual collar own-mean enclosure from same-source whole-phase U ranges.

The original radial mean obeys m_s=hb*(U-m). Integrating its nonzero core
inlet once preserves the signed core source and keeps every Cauchy phase
uncertainty as a separate formal-width sector. No source value is selected.
"""
from __future__ import annotations


def _scale(terms,factor,hb_power,origin):
    out=[]
    for term in terms:
        if len(term['powers'])!=4 or not hasattr(term['coefficient'],'_mpi_'):
            raise ValueError('four formal powers and directed coefficient required')
        powers=list(term['powers']);powers[0]+=hb_power
        out.append({'powers':tuple(powers),'coefficient':term['coefficient']*factor,
          'origin':(origin,{k:v for k,v in term.items() if k not in ('powers','coefficient')}),
          'semantics':term['semantics']})
    return out


def build(c, *, sc, phase, core_mean_Z, U_Z, source_family, source_atoms, admitted_source_contract):
    """FTC enclosure for ordinary mean Z0..Z5, retaining inlet and errors.

    U_Z[k] must enclose the same actual U_Zk over the full phase interval.
    It is split into the core inlet sector and explicit Cauchy error sectors.
    Source identity, analyticity and geometry are established by the caller.
    """
    sc,phase=c.mpf(sc),c.mpf(phase)
    if sc.a<=0 or sc.b>=c.mpf(1)/1024 or phase.a<0 or phase.b>(sc/4).a:
        raise ValueError('original microscopic collar phase required')
    if source_family is None or set(source_atoms)!=set(('hb','F0_anchor','Pstar_squared','Lambda')):
        raise ValueError('same-source atom records required')
    keys=('same_source_checked','core_mean_inlet_actual','U_whole_phase_enclosure','hb_Z_independent','rho4_inlet_history')
    if any(admitted_source_contract.get(k) is not True for k in keys):raise ValueError('caller actual-source contract required')
    if set(core_mean_Z)!=set(range(6)) or set(U_Z)!=set(range(6)):
        raise ValueError('complete ordinary Z0..Z5 rows required')
    # exprel(-hb*s) = (1-exp(-hb*s))/(hb*s) in [1-s/4,1]
    # because 0<hb<1/2, and is defined as 1 at s=0.
    lower=(1-phase/4).a
    if lower<0:raise ValueError('phase outside positive exponential range')
    exprel=c.make_mpf((lower._mpi_[0],c.mpf(1)._mpi_[1]))
    weight=phase*exprel
    out={}
    for k in range(6):
        base=c.mpf(core_mean_Z[k]);inputs=U_Z[k]
        if not inputs or inputs[0].get('semantics')!='actual_core_interval_not_selected_value':
            raise ValueError('actual U core inlet sector required')
        for term in inputs[1:]:
            if term.get('semantics')!='symmetric_certified_error_not_field_value' or term['powers'][0]<1:
                raise ValueError('separate certified phase uncertainty required')
        inlet=[{'powers':(0,0,0,0),'coefficient':base,'origin':('actual_core_own_mean_inlet',k),
                'semantics':'actual_core_interval_not_selected_value'}]
        # m = m0 + hb*s*exprel*(U-m0); each occurrence is kept distinct.
        correction=_scale(inlet,-weight,1,('inlet_exponential_correction',k))
        drive=_scale(inputs[:1],weight,1,('actual_U_history_integral',k))
        errors=_scale(inputs[1:],weight,1,('U_phase_error_history_integral',k))
        out[k]={'inlet':inlet,'signed_history_terms':correction+drive,'phase_uncertainty_terms':errors,
                'all_terms':inlet+correction+drive+errors}
    return {'source_family':source_family,'source_atoms':source_atoms,'mean_Z':out,
       'ordinary_Z_orders':list(range(6)),'same_source_contract_independently_proved':False,
       'history_formula':'m0+hb*s*exprel(-hb*s)*(U-m0), interval U over the full phase',
       'input_phase_errors_retained_separately':True,'core_history_reset':False,
       'actual_point_value_computed':False,'m15_computed':False,'allword_tail_added':False}
