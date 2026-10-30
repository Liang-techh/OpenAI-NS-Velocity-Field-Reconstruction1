"""Synthetic pointwise and interval checks for the original chi jet."""
def run(leaf,c):
    import mpmath as mp
    mp.mp.dps=80; delta=mp.mpf('0.2'); j=mp.mpf('0.7'); checks={}
    def f(z):
        h=z*(1-delta)/2+(1-z*z)*(4*z+j)
        return h*h/(h*h+(j/500)**2)
    for i,z in enumerate([mp.mpf('-1'),mp.mpf('-0.7'),mp.mpf('-0.2'),mp.mpf('0'),mp.mpf('0.2'),mp.mpf('0.7'),mp.mpf('1')]):
        padding=c.mpf(['-1e-60','1e-60']); zi=c.mpf(mp.nstr(z,80))+padding; di=c.mpf('0.2')+padding; ji=c.mpf('0.7')+padding; got=leaf.build(c,z=zi,delta=di,j=ji)['chi']
        for n in range(3):
            ref=mp.diff(f,z,n); lo=mp.make_mpf(got[n]._mpi_[0]); hi=mp.make_mpf(got[n]._mpi_[1]); checks[f'point_{i}_Z{n}']=bool(lo<=ref<=hi)
    box=leaf.build(c,z=c.mpf([-1,1]),delta=c.mpf('0.2')+c.mpf(['-1e-60','1e-60']),j=c.mpf('0.7')+c.mpf(['-1e-60','1e-60']))['chi']
    for i,z in enumerate([mp.mpf('-1'),mp.mpf('-0.7'),mp.mpf('-0.2'),mp.mpf('0'),mp.mpf('0.2'),mp.mpf('0.7'),mp.mpf('1')]):
        for n in range(3):
            ref=mp.diff(f,z,n);lo=mp.make_mpf(box[n]._mpi_[0]);hi=mp.make_mpf(box[n]._mpi_[1]);checks[f'box_{i}_Z{n}']=bool(lo<=ref<=hi)
    checks['chi_range']=bool(box[0].a>=0 and box[0].b<=1)
    return {'passed':all(checks.values()),'check_count':len(checks),'checks':checks,'synthetic_only':True}
