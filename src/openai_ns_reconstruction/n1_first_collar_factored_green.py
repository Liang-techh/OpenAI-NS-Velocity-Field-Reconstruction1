"""Source-factored, signed interval Green words on the first collar.

Each sector represents coefficient * hb**a * F0(Z)**b * Pstar2**c *
Lambda**d. The positive sources are not materialized. Source identity and the
validity of every caller-supplied interval extension remain caller duties.
Cellwise Green factors bound their positive exponential multipliers; no source
sector is collapsed, canceled, or replaced with an error radius.
"""
from __future__ import annotations

FIELDS=('F1','Uz1','K1','P1','d_x_F1','d_x_Uz1')
ROWS=('5','6','6','4','5','6')
C0_KEYS=frozenset(('4,1','5,1','5,2','5,3','5,5','6,1','6,2','6,3','6,4','6,6'))
B1_KEYS=frozenset(('5,1','5,2','5,3','6,2','6,3','6,4'))


def _check_sector(t):
    if set(('powers','coefficient','origin'))-set(t):raise ValueError('sector metadata missing')
    if len(t['powers'])!=4 or any(not isinstance(x,int) or isinstance(x,bool) for x in t['powers']):
        raise ValueError('four integer original source powers required')
    if not hasattr(t['coefficient'],'_mpi_'):raise TypeError('directed interval coefficient required')


def _scale(terms,factor,hb_power,origin):
    result=[]
    for t in terms:
        _check_sector(t);powers=list(t['powers']);powers[0]+=hb_power
        result.append({'powers':tuple(powers),'coefficient':t['coefficient']*factor,
          'origin':(origin,{k:v for k,v in t.items() if k not in ('powers','coefficient')})})
    return result


def _product(left,right,multiplier,origin):
    out=[]
    for a in left:
        _check_sector(a)
        for b in right:
            _check_sector(b)
            out.append({'powers':tuple(x+y for x,y in zip(a['powers'],b['powers'])),
              'coefficient':a['coefficient']*b['coefficient']*multiplier,
              'origin':(origin,{k:v for k,v in a.items() if k not in ('powers','coefficient')},
                        {k:v for k,v in b.items() if k not in ('powers','coefficient')})})
    return out


def _unit(c, lower):
    if lower.a<0:raise ValueError('microscopic exponential bound left positive range')
    return c.make_mpf((lower.a._mpi_[0],c.mpf(1)._mpi_[1]))


def _factors(c, phase):
    # 0 < hb < 1/2, so exp(-alpha hb phase)>=1-alpha phase/2.
    return (_unit(c,1-phase/2),_unit(c,1-3*phase/4),
            _unit(c,1-phase/4),_unit(c,1-phase/4))


def _transition(c,phase,initial):
    if len(initial)!=6:raise ValueError('six core components required')
    a,b,e,exprel=_factors(c,phase)
    f,u,k,p,fx,ux=initial
    return (list(f)+_scale(fx,phase*exprel,1,'F_transport_cross'),
       list(u)+_scale(ux,phase,1,'U_transport_cross'),
       _scale(k,a,0,'K_transport_diagonal')+_scale(ux,-phase*exprel,1,'K_transport_cross'),
       list(p),_scale(fx,b,0,'Fx_transport_diagonal'),_scale(ux,e,0,'Ux_transport_diagonal'))


def _partition(c,count):
    if not isinstance(count,int) or isinstance(count,bool) or count<1:raise ValueError('positive cell count required')
    return [(c.mpf(i)/count,c.mpf(i+1)/count) for i in range(count)]


def _green(c,phase,callback,count):
    """Return sector lists with actual cell width and formal hb measure."""
    fields={name:[] for name in FIELDS}
    for index,(left,right) in enumerate(_partition(c,count)):
        t=c.make_mpf((left.a._mpi_[0],right.b._mpi_[1]));v=phase*t;p=phase*(1-t)
        # u/2 = exp(hb*v/2) <= 1/(1-v/4), and exprel(-hb*p) is in [1-p/4,1].
        upper=(1/(1-v.b/4)).b
        uhalf=c.make_mpf((c.mpf(1)._mpi_[0],upper._mpi_[1]))
        a,b,e,exprel=_factors(c,p)
        # K=(u/2)*hb*p*exprel; log=(u/2)*hb*p.
        kernels=((uhalf*p*exprel,1),(uhalf*p,1),(-uhalf*p*exprel,1),
                 (c.mpf(1),0),(b,0),(e,0))
        source=callback(v)
        if set(source)!=set(('4','5','6')):raise ValueError('signed forcing rows 4,5,6 required')
        measure=phase*(right-left)*uhalf
        for field,row,(kernel,kpower) in zip(FIELDS,ROWS,kernels):
            fields[field].extend(_scale(source[row],measure*kernel,1+kpower,('local_green',index,row,field)))
    return fields


def enclose(c, *, sc, phase, core_W0_Z, core_W1_Z2, forcing_Z, matrix_Z,
            inner_cells, outer_cells, source_family, source_atoms, geometry_certificate):
    """Compute finite W1_Z2 sectors with a formal hb Green measure.

    Inputs are order-major source-sector lists using four powers [hb,F0,
    Pstar2,Lambda]. The callbacks must return signed interval extensions of
    the same source on each full phase/Z cell. This function establishes no
    source identity, leading-field existence, tail, or moment closure.
    """
    sc,phase=c.mpf(sc),c.mpf(phase)
    if sc.a<=0 or sc.b>=c.mpf(1)/1024 or phase.a<0 or phase.b>(sc/4).a:
        raise ValueError('first-collar phase and microscopic sc required')
    if source_family is None:raise ValueError('source family label required')
    if set(source_atoms)!=set(('hb','F0_anchor','Pstar_squared','Lambda')):
        raise ValueError('four same-source atom records required')
    required_geometry=('hb_strictly_between_zero_and_half','sc_phase_Z_independent','hb_epsilon_Z_independent')
    if any(geometry_certificate.get(k) is not True for k in required_geometry):
        raise ValueError('caller geometry certificates required')
    if set(core_W0_Z)!=set(range(4)) or any(len(core_W0_Z[k])!=6 for k in range(4)) or len(core_W1_Z2)!=6:
        raise ValueError('complete six-component core W0/W1 jets required')
    for values in (*core_W0_Z.values(),core_W1_Z2):
        for terms in values:
            for term in terms:_check_sector(term)
    _partition(c,inner_cells);_partition(c,outer_cells)
    def word0(v):
        out={}
        for order in range(4):
            initial=_transition(c,v,core_W0_Z[order])
            def body(w,order=order):
                data=forcing_Z(w)
                if set(data)!=set(range(4)):raise ValueError('g_Z0..Z3 required')
                return data[order]
            local=_green(c,v,body,inner_cells)
            out[order]=tuple(initial[i]+local[field] for i,field in enumerate(FIELDS))
        return out
    def feedback(v):
        word=word0(v);matrix=matrix_Z(v)
        if set(matrix)!=set(('C0','B1')):raise ValueError('C0/B1 required')
        for name,keys in (('C0',C0_KEYS),('B1',B1_KEYS)):
            if set(matrix[name])!=set(range(3)) or any(set(matrix[name][k])!=keys for k in range(3)):
                raise ValueError('raw matrix Z0..Z2 support required')
        result={str(k):[] for k in (4,5,6)}
        specs=(('C0',2,0,1),('C0',1,1,2),('C0',0,2,1),('B1',2,1,1),('B1',1,2,2),('B1',0,3,1))
        for name,mo,wo,factor in specs:
            for key,terms in matrix[name][mo].items():
                row,col=key.split(',')
                result[row].extend(_product(terms,word[wo][int(col)-1],factor,(name,mo,wo,key)))
        return result
    core=_transition(c,phase,core_W1_Z2)
    local=_green(c,phase,feedback,outer_cells)
    return {'source_family':source_family,'source_atoms':source_atoms,
       'caller_geometry_certificate':geometry_certificate,'geometry_certificate_independently_verified':False,
       'source_power_basis':['hb','F0_anchor','Pstar_squared','Lambda'],
       'core_transport':dict(zip(FIELDS,core)),'local_signed_sectors':local,
       'finite_W1_Z2':{name:core[i]+local[name] for i,name in enumerate(FIELDS)},
       'source_callback_identity_proved':False,'ordinary_Z_fixed_phase_and_R':True,
       'allword_tail_added':False,'m15_computed':False,'core_history_reset':False}
