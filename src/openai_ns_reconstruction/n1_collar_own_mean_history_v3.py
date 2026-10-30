"""Own-mean history and radial identities with formal collar width.

For rho=4 exp(hb*s), an actual signed U and its actual core-exit own mean
satisfy m(s)=exp(-hb*s)m(0)+hb*integral_0^s exp(-hb*(s-v))*U(v)dv.
Input source provenance is retained, and no midpoint or bound is selected as
an actual value. Source identity and callback cell enclosures are caller duties.
"""
from __future__ import annotations


def _unit(c,phase):
    lower=1-phase/2
    if lower.a<0:raise ValueError('phase outside exponential bound')
    return c.make_mpf((lower.a._mpi_[0],c.mpf(1)._mpi_[1]))


def _check(t):
    if not all(k in t for k in ('powers','coefficient','origin')):raise ValueError('sector metadata required')
    if len(t['powers'])!=4 or any(not isinstance(v,int) or isinstance(v,bool) for v in t['powers']):raise ValueError('four integer powers required')
    if not hasattr(t['coefficient'],'_mpi_'):raise TypeError('directed coefficient required')


def _scale(terms,factor,hb_power,origin):
    out=[]
    for term in terms:
        _check(term);p=list(term['powers']);p[0]+=hb_power
        out.append({'powers':tuple(p),'coefficient':term['coefficient']*factor,
          'origin':(origin,{k:v for k,v in term.items() if k not in ('powers','coefficient')})})
    return out


def build(c, *, sc, phase, inlet_mean_Z, U_Z, U_r_Z, cells, source_family, source_atoms, geometry_certificate):
    """Return signed m, m_r, m_rr ordinary Z jets from actual own history.

    inlet_mean_Z[k] is the complete core-exit own mean jet, k=0..5.
    U_Z(v)[k] covers the entire phase interval v (k=0..5).
    U_r_Z[k] covers the endpoint phase interval (k=0..3). Each jet is a list
    of signed formal sectors in (hb,F0_anchor,Pstar_squared,Lambda) order.
    Callers establish same-source identity and source-cell interval validity.
    """
    sc,phase=c.mpf(sc),c.mpf(phase)
    if sc.a<=0 or sc.b>=c.mpf(1)/1024 or phase.a<0 or phase.b>(sc/4).a:
        raise ValueError('microscopic first-collar phase required')
    if not isinstance(cells,int) or isinstance(cells,bool) or cells<1:raise ValueError('positive cell count required')
    if source_family is None or set(source_atoms)!=set(('hb','F0_anchor','Pstar_squared','Lambda')):
        raise ValueError('four source atoms and family required')
    certs=('hb_strictly_between_zero_and_half','sc_phase_Z_independent','hb_epsilon_Z_independent')
    if any(geometry_certificate.get(k) is not True for k in certs):raise ValueError('geometry certificates required')
    if set(inlet_mean_Z)!=set(range(6)) or set(U_r_Z)!=set(range(4)):
        raise ValueError('complete ordinary inlet Z0..5 and U_r Z0..3 required')
    for jets in (inlet_mean_Z,U_r_Z):
        for terms in jets.values():
            for term in terms:_check(term)
    mean={k:_scale(inlet_mean_Z[k],_unit(c,phase),0,('own_mean_inlet',k)) for k in range(6)}
    for i in range(cells):
        left=c.mpf(i)/cells;right=c.mpf(i+1)/cells
        t=c.make_mpf((left.a._mpi_[0],right.b._mpi_[1]));v=phase*t
        data=U_Z(v)
        if set(data)!=set(range(6)):raise ValueError('actual U_Z0..Z5 on entire cell required')
        factor=phase*(right-left)*_unit(c,phase*(1-t))
        for k in range(6):mean[k].extend(_scale(data[k],factor,1,('own_mean_local',i,k)))
    endpoint=U_Z(phase)
    if set(endpoint)!=set(range(6)):raise ValueError('actual U endpoint Z0..Z5 required')
    inverse_rho=_unit(c,phase)/4
    mr={k:_scale(endpoint[k],inverse_rho,0,('mean_r_U',k))+
          _scale(mean[k],-inverse_rho,0,('mean_r_minus_mean',k)) for k in range(5)}
    mrr={k:_scale(U_r_Z[k],inverse_rho,0,('mean_rr_Ur',k))+
           _scale(mr[k],-2*inverse_rho,0,('mean_rr_minus_2mr',k)) for k in range(4)}
    return {'mean':mean,'mean_r':mr,'mean_rr':mrr,'source_family':source_family,
       'source_atoms':source_atoms,'caller_geometry_certificate':geometry_certificate,
       'same_source_identity_proved':False,'callback_enclosures_proved':False,
       'history_formula':'exp(-hb*s)*m(0)+hb*integral exp(-hb*(s-v))*U(v)dv',
       'rho_formula':'4*exp(hb*s)','core_history_reset':False,'signed_sectors_preserved':True,
       'ordinary_Z_fixed_phase_and_R':True,'m15_computed':False}
