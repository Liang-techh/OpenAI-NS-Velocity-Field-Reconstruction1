"""Bind saved core source enclosures to the exact first-collar inlet."""
FIELDS=('F1','Uz1','K1','P1','d_x_F1','d_x_Uz1')
ROWS=('5','6','6','4','5','6')

def build(c,low,g3fields,w1fields):
    masses=(c.mpf(1)/2,c.mpf(1),-c.mpf(1)/2,c.mpf(2),c.mpf(1)/2,c.mpf(1))
    if set(low)!={0,1,2} or set(g3fields)!=set(FIELDS) or set(w1fields)!=set(FIELDS):raise ValueError('incomplete inlet')
    index={};w0={};w1={}
    def sector(term,coefficient,sid,classes,ref):
        index[sid]=ref
        return {'powers':[0,term['F0_power'],0,term['Lambda_power']],'coefficient':coefficient,
          'origin':{'source_ids':[sid],'uncertainty':True,'error_classes':classes}}
    for k in range(3):
        w0[k]={}
        for field,row,mass in zip(FIELDS,ROWS,masses):
            result=[]
            for i,t in enumerate(low[k][row]):
                sid=f'core_low:{k}:{row}:{i}'
                result.append(sector(t,t['coefficient']*mass,sid,['core_whole_source_interval'],{'receipt':'core_word1_z2_v2','path':['global_low_source',str(k),row,i]}))
            w0[k][field]=result
    w0[3]={}
    for field,row in zip(FIELDS,ROWS):
        result=[];cells=set()
        for i,t in enumerate(g3fields[field]):
            cell=t['cell'];cells.add(cell);sid=f'core_g3:{field}:{i}'
            result.append(sector(t,t['coefficient_interval'],sid,[f'core_g3_cell_{cell}'],{'receipt':'core_g3_word0_v1','path':['W0_Z3_green_enclosure','fields',field,i],'source_cell':cell,'source_row':t['source_row']}))
        if len(cells)!=8:raise ValueError('missing g3 source cells '+field)
        w0[3][field]=result
    for field in FIELDS:
        result=[]
        for i,t in enumerate(w1fields[field]):
            sid=f'core_W1:{field}:{i}';classes=['core_W1_outer_'+str(t['outer_cell'])]
            if t.get('W0_source_cell') is not None:classes.append('core_g3_cell_'+str(t['W0_source_cell']))
            result.append(sector(t,t['coefficient'],sid,classes,{'receipt':'core_word1_z2_v2','path':['result','W1_Z2_core_exit',field,i],'source_row':t['source_row'],'outer_cell':t['outer_cell'],'W0_source_cell':t.get('W0_source_cell')}))
        w1[field]=result
    return {'core_W0_Z':w0,'core_W1_Z2':w1,'source_index':index,'Z_derivatives_are_ordinary':True,'g3_and_W1_already_integrated':True,'point_values_not_claimed':True}
