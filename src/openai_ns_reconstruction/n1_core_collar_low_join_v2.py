"""No-reset low-word core transport and separate first-collar sectors."""
FIELDS=('F1','Uz1','K1','P1','d_x_F1','d_x_Uz1')
ATOM_NAMES=('hb','F0_anchor','Pstar_squared','Lambda')

def build(c,*,sc,phase,core,local,source_family,source_atoms,geometry_certificate):
    if source_family is None or set(source_atoms)!=set(ATOM_NAMES):raise ValueError('source family and four atoms required')
    for key in ('hb_strictly_between_zero_and_half','sc_phase_Z_independent','hb_epsilon_Z_independent'):
        if geometry_certificate.get(key) is not True:raise ValueError('missing geometry certificate '+key)
    sc=c.mpf(sc);phase=c.mpf(phase)
    if sc.a<=0 or sc.a!=sc.b or sc.b>=c.mpf(1)/1024 or phase.a<0 or phase.b>(sc/4).a:raise ValueError('invalid exact sc or phase')
    if set(core)!={0,1} or set(local)!={0,1}:raise ValueError('ordinary Z0 and Z1 required')
    def unit(lower):
        if lower.a<0:raise ValueError('exponential lower bound')
        return c.make_mpf((lower.a._mpi_[0],c.mpf(1)._mpi_[1]))
    # hb<1/2: exp(-a hb s) >= 1-a*s/2, and
    # (1-exp(-hb*s))/(hb*s) >= 1-s/4.
    e1=unit(1-phase/2);e3=unit(1-3*phase/4);eh=unit(1-phase/4);exprel=unit(1-phase/4)
    def sector(t,factor,hb,operation):
        sid=t['input_id'];outer=t.get('outer_cell');classes=['core_W1_outer_'+str(outer),'core_source_detail_not_separately_classified']
        return {'powers':[hb,t['F0_power'],0,t['Lambda_power']],'coefficient':t['coefficient']*factor,
          'origin':{'source_ids':[sid],'uncertainty':True,'error_classes':classes,'operation':operation}}
    def scale(terms,factor,hb,operation):return [sector(t,factor,hb,operation) for t in terms]
    def copy_local(t):
        origin=dict(t['origin']);ids=list(origin.get('source_ids',[]));sid=t.get('input_id')
        if sid is not None and sid not in ids:ids.append(sid)
        origin['source_ids']=ids
        return {'powers':list(t['powers']),'coefficient':t['coefficient'],'origin':origin}
    transported={};local_out={};finite={}
    for order in (0,1):
        if set(core[order])!=set(FIELDS) or set(local[order])!=set(FIELDS):raise ValueError('six component fields required')
        f,u,k,p,fx,ux=(core[order][name] for name in FIELDS)
        values=(scale(f,c.mpf(1),0,'F_diagonal')+scale(fx,phase*exprel,1,'F_cross'),
                scale(u,c.mpf(1),0,'U_diagonal')+scale(ux,phase,1,'U_cross'),
                scale(k,e1,0,'K_diagonal')+scale(ux,-phase*exprel,1,'K_cross'),
                scale(p,c.mpf(1),0,'P_diagonal'),scale(fx,e3,0,'Fx_diagonal'),scale(ux,eh,0,'Ux_diagonal'))
        transported[order]=dict(zip(FIELDS,values));local_out[order]={name:[copy_local(t) for t in local[order][name]] for name in FIELDS}
        finite[order]={name:transported[order][name]+local_out[order][name] for name in FIELDS}
    return {'source_family':source_family,'source_atoms':source_atoms,'source_power_basis':list(ATOM_NAMES),
      'core_transport':transported,'local_contribution':local_out,'finite_sum':finite,
      'ordinary_Z_orders':[0,1],'core_history_reset':False,'local_added_once':True,'allword_tail_added':False,
      'm15_computed':False,'hb_materialized':False,'same_source_cancellation_claimed':False,'finite_sum_is_unmerged_sector_list':True,
      'geometry_certificate_independently_verified':False,'source_callback_identity_proved':False}
