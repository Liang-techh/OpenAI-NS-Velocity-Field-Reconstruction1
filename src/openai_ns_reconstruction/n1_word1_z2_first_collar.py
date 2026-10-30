"""Transport the core W1_Z2 and enclose its local first-collar increment.

The caller supplies actual same-source core sectors, an independently
certified local feedback-Z2 norm, and directed interval arithmetic.  hb is
positive and remains a formal factor.  The local radius is an error bound,
not a signed source value or an all-word tail.
"""
from __future__ import annotations


FIELDS=('F1','Uz1','K1','P1','d_x_F1','d_x_Uz1')


def enclose(c, core, *, sc, q, feedback_Z2_upper):
    """Cover 0<=q<=1 with s=sc*q/4, X=2*exp(hb*s/2).

    The inputs assert 0<hb<1/2 and 0<sc<1/1024.  Core terms have F0_power,
    Lambda_power and a directed ``coefficient`` interval.  This computes only
    the finite first Volterra word's ordinary Z2, at fixed phase s.
    """
    if not (sc.a>0 and sc.b<c.mpf(1)/1024 and q.a>=0 and q.b<=1
            and feedback_Z2_upper.a>=0):
        raise ValueError('positive microscopic width, q range and bound required')
    if set(core)!=set(FIELDS):
        raise ValueError('six original core W1_Z2 fields required')
    for field in FIELDS:
        if not core[field]:raise ValueError('empty core signed source field '+field)
    s=sc*q/4
    # 0 <= hb*s <= s/2.  Taylor inequalities bound only positive
    # multipliers; hb itself is never assigned a numerical representative.
    theta=s.b/2
    if not theta<c.mpf(1)/2:raise ValueError('invalid exponential domain')
    one=c.mpf(1)
    def positive_unit(lower):
        return c.make_mpf((lower._mpi_[0],one._mpi_[1]))
    decay=positive_unit((one-theta).a)
    decay_three_halves=positive_unit((one-3*theta/2).a)
    decay_half=positive_unit((one-theta/2).a)
    exprel_negative=positive_unit((one-theta/2).a)
    def emit(field, factor, hb_power, origin):
        return [{'F0_power':v['F0_power'],'Lambda_power':v['Lambda_power'],
            'hb_power':hb_power,'coefficient_interval':v['coefficient']*factor,
            'core_field':field,'origin':origin,
            'core_provenance':{key:value for key,value in v.items() if key!='coefficient'}} for v in core[field]]
    transported={
        'F1':emit('F1',one,0,'core F')+emit('d_x_F1',s*exprel_negative,1,'(1-exp(-hb*s))*core Fx'),
        'Uz1':emit('Uz1',one,0,'core U')+emit('d_x_Uz1',s,1,'hb*s*core Ux'),
        'K1':emit('K1',decay,0,'exp(-hb*s)*core K')+emit('d_x_Uz1',-s*exprel_negative,1,'-(1-exp(-hb*s))*core Ux'),
        'P1':emit('P1',one,0,'core P'),
        'd_x_F1':emit('d_x_F1',decay_three_halves,0,'exp(-3*hb*s/2)*core Fx'),
        'd_x_Uz1':emit('d_x_Uz1',decay_half,0,'exp(-hb*s/2)*core Ux'),
    }
    # X-2=hb*s*exprel(hb*s/2) <= hb*s/(1-theta/2).
    D=(one/(one-theta/2)).b
    dx_without_hb=(s.b*D).b
    local={}
    for field in FIELDS:
        power=2 if field in ('F1','Uz1','K1') else 1
        radius=(feedback_Z2_upper.b*(dx_without_hb**power)).b
        local[field]={'hb_power':power,'positive_coefficient_radius':radius,
            'meaning':'symmetric local first-collar Green error, not a field value'}
    return {'core_transport_signed_sectors':transported,'local_feedback_symmetric_radii':local,
        's_interval':s,'decay_enclosures':{'exp_minus_hb_s':decay,
            'exp_minus_3hb_s_over_2':decay_three_halves,'exp_minus_hb_s_over_2':decay_half,
            'exprel_minus_hb_s':exprel_negative},
        'dx_upper_without_hb':dx_without_hb,'local_feedback_Z2_upper':feedback_Z2_upper,
        'no_seam_reset':True,'core_history_transported_once':True,
        'local_collar_feedback_is_radius_not_signed_value':True,
        'allword_tail_added':False,'m15_integral_computed':False,
        'scope':'finite W1_Z2; first-collar q in supplied range; fixed phase ordinary Z'}
