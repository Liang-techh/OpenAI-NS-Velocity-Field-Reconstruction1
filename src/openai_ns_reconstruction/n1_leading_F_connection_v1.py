"""Physical F=F0*Phi ordinary Z0/Z1 sectors without materializing F0."""
def build(c,*,z,delta,j,phi,connection_source_id):
    z=c.mpf(z);delta=c.mpf(delta);j=c.mpf(j)
    if j.a<=0:raise ValueError('positive j required')
    L=1-delta*z**2;d=1-z**2;H=z*(1-delta)/2+d*(4*z+j)
    denominator=H**2+(j/500)**2
    if L.a<=0 or denominator.a<=0:raise ValueError('source denominator is not positive')
    ell_over_lambda=-L*H/denominator
    if set(phi)!={0,1}:raise ValueError('ordinary Phi Z0/Z1 required')
    def lift(term,extra_lambda,scale,connection):
        powers=list(term['powers']);powers[1]+=1;powers[3]+=extra_lambda
        old=term['origin'];ids=list(old['source_ids']);classes=list(old['error_classes'])
        if connection:
            ids.append(connection_source_id);classes.append('F0_connection_interval')
        return {'powers':powers,'coefficient':term['coefficient']*scale,
          'origin':{'source_ids':list(dict.fromkeys(ids)),'uncertainty':bool(old['uncertainty'] or connection),'error_classes':list(dict.fromkeys(classes))}}
    value=[lift(t,0,c.mpf(1),False) for t in phi[0]]
    derivative=[lift(t,0,c.mpf(1),False) for t in phi[1]]+[lift(t,1,ell_over_lambda,True) for t in phi[0]]
    return {'leading':{0:value,1:derivative},'ell_over_lambda':ell_over_lambda,'ordinary_Z1_uses_F0_connection_once':True,
      'positive_F0_and_Lambda_formal':True,'same_source_cancellation_claimed':False}
