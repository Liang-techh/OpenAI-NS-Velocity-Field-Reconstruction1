"""Four-cell ordinary Z0/Z1 source enclosures of first-collar W0."""
FIELDS=('F1','Uz1','K1','P1','d_x_F1','d_x_Uz1')
ROWS=('5','6','6','4','5','6')
ATOMS=('hb','F0_anchor','Pstar_squared','Lambda')

def build(c,*,sc,core,g,source_family,source_atoms,certificate):
    if source_family is None or set(source_atoms)!=set(ATOMS):raise ValueError('four same-source atoms required')
    for k in ('hb_strictly_between_zero_and_half','sc_phase_Z_independent','hb_epsilon_Z_independent'):
        if certificate.get(k) is not True:raise ValueError('missing geometry certificate '+k)
    sc=c.mpf(sc)
    if sc.a<=0 or sc.a!=sc.b or sc.b>=c.mpf(1)/1024:raise ValueError('exact microscopic sc required')
    if set(core)!={0,1} or set(g)!={0,1}:raise ValueError('ordinary Z0/Z1 required')
    def box(lo,hi):return c.make_mpf((lo.a._mpi_[0],hi.b._mpi_[1]))
    def mapped(term):
        origin=dict(term['origin']);ids=list(origin.get('source_ids',[]));sid=term.get('input_id')
        if sid is not None and sid not in ids:ids.append(sid)
        if not ids:raise ValueError('source IDs required')
        origin['source_ids']=ids
        return {'powers':list(term['powers']),'coefficient':term['coefficient'],'origin':origin}
    core_out={};local={};cells=[]
    for order in (0,1):
        if set(core[order])!=set(FIELDS) or set(g[order])!={'4','5','6'}:raise ValueError('incomplete field/source rows')
        core_out[order]={field:[mapped(t) for t in core[order][field]] for field in FIELDS};local[order]={}
        for index in range(4):
            lo=sc*index/16;hi=sc*(index+1)/16;ratio=1/(1-hi/4)
            # u/2=exp(hb*v/2) is in [1,ratio], and
            # exprel(-hb*(s-v)) >= 1-hi/4.  The kernels are nonnegative
            # apart from K, with exact integral of (s-v) equal to s^2/2.
            f=box(lo**2/2*(1-hi/4),hi**2/2*ratio**2)
            u=box(lo**2/2,hi**2/2*ratio**2)
            p=box(lo,hi*ratio)
            fx=box(lo*(1-3*hi/4),hi*ratio)
            ux=box(lo*(1-hi/4),hi*ratio)
            factors=(f,u,-f,p,fx,ux);orders=(2,2,2,1,1,1)
            local[order][index]={}
            for field,row,factor,hbpower in zip(FIELDS,ROWS,factors,orders):
                terms=[]
                for t in g[order][row]:
                    q=mapped(t);q['powers'][0]+=hbpower;q['coefficient']*=factor
                    q['origin']={**q['origin'],'operation':'local_green_cell','phase_cell':index,'source_row':row}
                    terms.append(q)
                local[order][index][field]=terms
            if order==0:cells.append({'q_lower':c.mpf(index)/4,'q_upper':c.mpf(index+1)/4,'s_lower':lo,'s_upper':hi})
    return {'source_family':source_family,'source_atoms':source_atoms,'source_power_basis':list(ATOMS),
      'core_transport':core_out,'local':local,'q_cells':cells,'ordinary_Z_orders':[0,1],
      'finite_word_definition':'core_transport[field] + local[order][cell][field]',
      'g_uniform_on_entire_first_collar_required':True,'core_history_reset':False,'local_added_once_in_definition':True,
      'allword_tail_added':False,'m15_computed':False,'hb_materialized':False,'same_source_cancellation_claimed':False,
      'geometry_certificate_independently_verified':False}
