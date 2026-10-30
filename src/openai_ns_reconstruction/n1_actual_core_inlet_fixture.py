"""Synthetic independent checks of core inlet integration bookkeeping."""
def run(leaf,c):
    low={k:{r:[{'F0_power':1,'Lambda_power':-1,'coefficient':c.mpf(k+2)}] for r in ('4','5','6')} for k in range(3)}
    g3={f:[{'F0_power':1,'Lambda_power':-1,'coefficient_interval':c.mpf(i+1),'cell':i,'source_row':'5'} for i in range(8)] for f in leaf.FIELDS}
    w1={f:[{'F0_power':0,'Lambda_power':-2,'coefficient':c.mpf(7),'outer_cell':2,'source_row':'6','W0_source_cell':3}] for f in leaf.FIELDS}
    result=leaf.build(c,low,g3,w1);expected=(c.mpf(1)/2,c.mpf(1),-c.mpf(1)/2,c.mpf(2),c.mpf(1)/2,c.mpf(1));checks={}
    for i,f in enumerate(leaf.FIELDS):
        for k in range(3):checks[f'low_{f}_{k}']=result['core_W0_Z'][k][f][0]['coefficient']==(k+2)*expected[i]
        checks['g3_'+f]=all(t['coefficient']==c.mpf(j+1) and 'core_g3_cell_'+str(j) in t['origin']['error_classes'] for j,t in enumerate(result['core_W0_Z'][3][f]))
        checks['w1_'+f]=result['core_W1_Z2'][f][0]['coefficient']==c.mpf(7) and 'core_g3_cell_3' in result['core_W1_Z2'][f][0]['origin']['error_classes']
    checks['ordinary_Z']=result['Z_derivatives_are_ordinary'];checks['source_index']=all(t['origin']['source_ids'][0] in result['source_index'] for terms in result['core_W0_Z'][3].values() for t in terms)
    return {'passed':all(checks.values()),'check_count':len(checks),'checks':checks,'synthetic_only':True}
