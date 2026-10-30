"""Signed first-collar Green words from caller-bound *actual* source jets.

Ordinary Z derivatives are taken at fixed phase. Callers supply interval
extensions of the same-source forcing and dynamic matrices on each requested
phase interval, plus the complete core-exit jets. No midpoint, norm cap,
previous-error deletion, or all-word completion is performed here.

This numerical leaf materializes hb; it does not preserve formal source-log
sectors or bind exact source nodes. It is valid only after an independent
source-bound adapter proves its inputs enclose the same actual fields on every
requested interval. A synthetic run cannot be used as actual field evidence.
"""
from __future__ import annotations

FIELDS = ('F1', 'Uz1', 'K1', 'P1', 'd_x_F1', 'd_x_Uz1')
ROWS = ('5', '6', '6', '4', '5', '6')
C0_KEYS = frozenset(('4,1','5,1','5,2','5,3','5,5','6,1','6,2','6,3','6,4','6,6'))
B1_KEYS = frozenset(('5,1','5,2','5,3','6,2','6,3','6,4'))


def matrix_entries_by_order(raw):
    """Reorder the companion matrix leaf's ordinary jets without changing values."""
    if set(raw)!=set(('C0','B1')):
        raise ValueError('C0 and B1 matrices required')
    result={}
    for name,expected in (('C0',C0_KEYS),('B1',B1_KEYS)):
        if set(raw[name])!=expected or any(len(v)!=3 for v in raw[name].values()):
            raise ValueError('raw matrix supports and three-component jets required')
        result[name]={k:{key:value[k] for key,value in raw[name].items()} for k in range(3)}
    return result


def _exprel_negative(c, distance):
    # exp(-d) <= (1-exp(-d))/d <= 1 on d>=0, including d=0.
    if distance.a < 0:
        raise ValueError('negative phase distance')
    return c.make_mpf((c.exp(-distance).a._mpi_[0], c.mpf(1)._mpi_[1]))


def _transition(c, hb, phase, initial):
    d = hb*phase
    cross = d*_exprel_negative(c,d)
    a = c.exp(-d)
    b = c.exp(-3*d/2)
    e = c.exp(-d/2)
    f,u,k,p,fx,ux = initial
    return (f+cross*fx,u+d*ux,a*k-cross*ux,p,b*fx,e*ux)


def _partition(c, count):
    if not isinstance(count,int) or isinstance(count,bool) or count < 1:
        raise ValueError('positive integer partition count required')
    return [(c.mpf(i)/count,c.mpf(i+1)/count) for i in range(count)]


def _integral(c, hb, phase, source, count):
    """Enclose the six raw Green integrals by interval cell ranges.

    source(v) is a dict for signed rows 4,5,6 covering the entire v interval.
    The phase transform v=phase*t and positive measure hb*u*phase/2 are
    retained; ranges are integrated only against their actual cell widths.
    """
    terms={name:[] for name in FIELDS}
    for left,right in _partition(c,count):
        t=c.make_mpf((left.a._mpi_[0],right.b._mpi_[1]))
        v=phase*t
        u=2*c.exp(hb*v/2)
        distance=hb*phase*(1-t)
        ratio=c.exp(-distance/2)
        K=u*distance*_exprel_negative(c,distance)/2
        L=u*distance/2
        kernels=(K,L,-K,c.mpf(1),ratio**3,ratio)
        values=source(v)
        if set(values)!=set(('4','5','6')):
            raise ValueError('forcing rows 4, 5, 6 required')
        measure=hb*u*phase*(right-left)/2
        for name,row,kernel in zip(FIELDS,ROWS,kernels):
            terms[name].append(values[row]*kernel*measure)
    return terms


def _sum(c, terms):
    value=c.mpf(0)
    for term in terms:value=value+term
    return value


def enclose(c, *, hb, sc, phase, core_W0_Z, core_W1_Z2,
            forcing_Z, matrix_Z, inner_cells, outer_cells, source_family):
    """Enclose finite W0_Z0..3 and W1_Z2 using signed same-source jets.

    core_W0_Z[order][field] and core_W1_Z2[field] are six-vectors; forcing_Z
    returns {order: {row: interval}} for order 0..3, and matrix_Z returns
    {'C0': {order: {'row,col': interval}}, 'B1': ...} for orders 0..2.
    Each callback must return ranges valid for the full supplied phase interval
    and the same real Z range. All Z jets are ordinary derivatives, not Taylor
    coefficients. The caller is responsible for source identity and error bounds.
    The output contains per-cell signed terms as well as interval sums.
    """
    hb,sc,phase=(c.mpf(x) for x in (hb,sc,phase))
    if hb.a<=0 or hb.b>=c.mpf('0.5') or sc.a<=0 or sc.b>=c.mpf(1)/1024:
        raise ValueError('actual first-collar hb/sc admissibility required')
    if phase.a<0 or phase.b>(sc/4).a:
        raise ValueError('phase outside actual first collar')
    if source_family is None:
        raise ValueError('same source family must be identified')
    if set(core_W0_Z)!=set(range(4)) or any(len(core_W0_Z[k])!=6 for k in range(4)) or len(core_W1_Z2)!=6:
        raise ValueError('complete six-component core-exit W0 jets and W1_Z2 required')
    _partition(c,inner_cells);_partition(c,outer_cells)

    def word0(v):
        initial={k:_transition(c,hb,v,core_W0_Z[k]) for k in range(4)}
        out={}
        for order in range(4):
            def body(w, order=order):
                data=forcing_Z(w)
                if set(data)!=set(range(4)):
                    raise ValueError('forcing Z0..Z3 required')
                return data[order]
            local=_integral(c,hb,v,body,inner_cells)
            out[order]=tuple(initial[order][i]+_sum(c,local[name]) for i,name in enumerate(FIELDS))
        return out

    def feedback(v):
        word=word0(v);matrix=matrix_Z(v)
        if set(matrix)!=set(('C0','B1')):
            raise ValueError('actual C0 and B1 jets required')
        result={str(k):c.mpf(0) for k in (4,5,6)}
        for name,expected in (('C0',C0_KEYS),('B1',B1_KEYS)):
            if set(matrix[name])!=set(range(3)) or any(set(matrix[name][k])!=expected for k in range(3)):
                raise ValueError('matrix Z0..Z2 and original raw support required')
        specs=(('C0',2,0,1),('C0',1,1,2),('C0',0,2,1),
               ('B1',2,1,1),('B1',1,2,2),('B1',0,3,1))
        for name,mo,wo,factor in specs:
            for key,value in matrix[name][mo].items():
                row,column=key.split(',')
                result[row]+=factor*value*word[wo][int(column)-1]
        return result

    core=_transition(c,hb,phase,core_W1_Z2)
    local=_integral(c,hb,phase,feedback,outer_cells)
    total={field:core[i]+_sum(c,local[field]) for i,field in enumerate(FIELDS)}
    return {'source_family':source_family,'core_transport':dict(zip(FIELDS,core)),
        'local_signed_cell_terms':local,'finite_W1_Z2':total,
        'Z_derivative':'ordinary at fixed phase and fixed physical R',
        'formal_source_sectors_preserved':False,'exact_source_identity_proved':False,
        'local_terms_are_interval_enclosures':True,
        'caller_must_certify_source_callbacks':True,'allword_tail_added':False,
        'm15_computed':False,'core_history_reset':False}
