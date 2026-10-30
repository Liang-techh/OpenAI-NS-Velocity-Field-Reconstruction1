"""Exact rational product-rule check for the source-factored F0 connection."""
def run(leaf,c):
    from fractions import Fraction as F
    def iv(q):return c.mpf(q.numerator)/q.denominator
    z=F(1,3);delta=F(1,5);j=F(7,3);lam=F(3);f0=F(2);p=F(7,4);pz=F(-5,6)
    source={0:[{'powers':[0,0,0,0],'coefficient':iv(p),'origin':{'source_ids':['p0'],'uncertainty':True,'error_classes':['source_box']}}],1:[{'powers':[0,0,0,0],'coefficient':iv(pz),'origin':{'source_ids':['p1'],'uncertainty':False,'error_classes':[]}}]}
    out=leaf.build(c,z=iv(z),delta=iv(delta),j=iv(j),phi=source,connection_source_id='ell')
    L=1-delta*z*z;H=z*(1-delta)/2+(1-z*z)*(4*z+j);ell_over=-L*H/(H*H+(j/500)**2)
    def value(terms):
        total=c.mpf(0)
        for t in terms:total+=t['coefficient']*iv(f0)**t['powers'][1]*iv(lam)**t['powers'][3]
        return total
    target0=iv(f0*p);target1=iv(f0*(pz+lam*ell_over*p));v0=value(out['leading'][0]);v1=value(out['leading'][1])
    checks={'value_contains':v0.a<=target0.a and v0.b>=target0.b,'ordinary_Z1_contains':v1.a<=target1.a and v1.b>=target1.b,
      'formal_connection_powers':out['leading'][1][1]['powers']==[0,1,0,1],
      'input_ids_retained':out['leading'][1][1]['origin']['source_ids']==['p0','ell'],
      'source_error_retained':'source_box' in out['leading'][1][1]['origin']['error_classes'],
      'connection_once':out['ordinary_Z1_uses_F0_connection_once'],'no_source_factor_materialized':out['positive_F0_and_Lambda_formal']}
    return {'passed':all(checks.values()),'checks':checks,'check_count':len(checks),'synthetic_only':True,'reference':'exact rational product rule and F0 logarithmic derivative'}
