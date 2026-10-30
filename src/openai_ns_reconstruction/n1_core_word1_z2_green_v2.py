"""Whole-cell signed enclosures for the core W1 ordinary-Z2 Green word.

All source and matrix inputs are caller-bound enclosures of the same actual
core. Positive F0(Z) and Lambda factors stay formal and cell terms stay
separate. This module neither runs an ancestor nor chooses midpoint values.
"""
from __future__ import annotations


FIELDS=(('F1','5'),('Uz1','6'),('K1','6'),('P1','4'),('d_x_F1','5'),('d_x_Uz1','6'))
EXPECTED_C0={'4,1','5,1','5,2','5,3','5,5','6,1','6,2','6,3','6,4','6,6'}
EXPECTED_B1={'5,1','5,2','5,3','6,2','6,3','6,4'}


def _term(term, factor, **extra):
    return {'F0_power':term['F0_power'],'Lambda_power':term['Lambda_power'],
            'coefficient':term['coefficient']*factor,**extra}


def _product(a,b,multiplier,origin,column):
    return {'F0_power':a['F0_power']+b['F0_power'],
        'Lambda_power':a['Lambda_power']+b['Lambda_power'],
        'coefficient':multiplier*a['coefficient']*b['coefficient'],
        'origin':origin,'column':column,'W0_source_cell':b.get('source_cell')}


def build(c, low_source, saved_g3_cells, matrix_cells):
    """Enclose W0 jets on eight outer x cells, then W1_Z2 at rho=4.

    low_source[order][row] is a signed, global source enclosure for orders
    0,1,2. saved_g3_cells are the previously computed eight absolute source
    cells. matrix_cells contain native C0/B1 jets at each outer rho cell.
    Source/cell/tube identity and directed input are caller responsibilities.
    """
    if len(saved_g3_cells)!=8 or len(matrix_cells)!=8:
        raise ValueError('eight original source and eight outer cells required')
    output={name:[] for name,_ in FIELDS};outer=[]
    for index,cell in enumerate(matrix_cells):
        left=c.mpf(index)/8;right=c.mpf(index+1)/8
        t=cell['t'];rho=cell['rho'];x=2*t
        if t.a!=left.a or t.b!=right.b or rho._mpi_!=(4*t*t)._mpi_:
            raise ValueError('outer cell is not rho=4*t^2')
        for name,expected in (('C0',EXPECTED_C0),('B1',EXPECTED_B1)):
            for order in ('0','1','2'):
                if set(cell[name][order])!=expected:
                    raise ValueError('raw dynamic matrix support changed')
        x2=x**2
        masses={'F1':x2/8,'Uz1':x2/4,'K1':-x2/8,'P1':x,'d_x_F1':x/4,'d_x_Uz1':x/2}
        word={str(k):{} for k in range(4)}
        for col,(field,row) in enumerate(FIELDS,1):
            mass=masses[field]
            for order in ('0','1','2'):
                word[order][col]=[_term(v,mass,source_cell=None) for v in low_source[order][row]]
            signed_weight=c.mpf([0,abs(mass).b])
            if field=='K1':signed_weight=-signed_weight
            terms=[]
            for source_index,source in enumerate(saved_g3_cells):
                # Source x lies in [source_index/4,(source_index+1)/4].
                # A cell at or beyond every allowed endpoint has zero measure.
                if x.b<=c.mpf(source_index)/4:continue
                for v in source['g'][row]:
                    terms.append(_term(v,signed_weight,source_cell=source_index))
            word['3'][col]=terms
        specs=(('C0','2','0',1,'C0_ZZ W0'),('C0','1','1',2,'2 C0_Z W0_Z'),
               ('C0','0','2',1,'C0 W0_ZZ'),('B1','2','1',1,'B1_ZZ W0_Z'),
               ('B1','1','2',2,'2 B1_Z W0_ZZ'),('B1','0','3',1,'B1 W0_ZZZ'))
        source={'4':[],'5':[],'6':[]}
        for matrix,mo,wo,multiplier,origin in specs:
            for key,terms in cell[matrix][mo].items():
                row,col=(int(v) for v in key.split(','))
                if row not in (4,5,6) or col not in range(1,7):
                    raise ValueError('invalid raw matrix coordinate')
                for a in terms:
                    for b in word[wo][col]:
                        source[str(row)].append(_product(a,b,multiplier,origin,col))
        tlog=c.mpf([0,'0.5']) if index==0 else t*c.log(1/t)
        K=t*(1-t*t)
        kernels={'F1':K,'Uz1':2*tlog,'K1':-K,'P1':c.mpf(1),'d_x_F1':t**3,'d_x_Uz1':t}
        for field,row in FIELDS:
            for term in source[row]:
                output[field].append(_term(term,kernels[field]/4,outer_cell=index,
                    source_row=row,origin=term['origin'],column=term['column'],
                    W0_source_cell=term['W0_source_cell']))
        outer.append({'t':t,'rho':rho,'W0_ordinary_Z':word,'feedback_Z2':source})
    return {'outer_cells':outer,'W1_Z2_core_exit':output,'ordinary_Z_order':2,
        'coordinate':'rho=4, core exit, real Z in [-1,1]',
        'W0_low_order_source':'global signed enclosure times exact positive Green mass',
        'W0_Z3_source':'distinct previous source cells times individual nonnegative Green-weight bounds',
        'no_source_cells_collapsed':True,'axis_endpoint_continuation':True,
        'source_exponentials_materialized':False,'first_collar_added':False,
        'allword_tail_added':False,'m15_computed':False}
