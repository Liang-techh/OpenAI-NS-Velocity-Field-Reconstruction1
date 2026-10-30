"""Synthetic first-word physical swirl m15 check with a nonconstant profile."""
def run(leaf,c):
    import mpmath as mp
    def term(value,sid,uncertain=False):return {'powers':[0,0,0,0],'coefficient':value,'origin':{'source_ids':[sid],'uncertainty':uncertain,'error_classes':['synthetic_error'] if uncertain else []}}
    leading={0:[term(c.mpf([2,3]),'F0')],1:[term(c.mpf(3),'F1',True)]}
    for group in leading.values():
        for t in group:t['powers'][1]=1
    first={0:{'core':[term(c.mpf([1,3]),'Wc0')],'local':[term(c.mpf([-1,0]),'Wl0',True)]},1:{'core':[term(c.mpf(2),'Wc1')],'local':[term(c.mpf(-1),'Wl1',True)]}}
    for branches in first.values():
        for group in branches.values():
            for t in group:t['powers'][1]=1;t['powers'][3]=-1
    sc=c.mpf(1)/2048;atoms={k:'synthetic_'+k for k in ('hb','F0_anchor','Pstar_squared','Lambda')}
    cert={'hb_strictly_between_zero_and_half':True,'sc_phase_Z_independent':True,'hb_epsilon_Z_independent':True,'original_cutoff_mass_5_over_8':True}
    result=leaf.build(c,sc=sc,leading=leading,first_word=first,source_family='synthetic',source_atoms=atoms,certificate=cert)
    def value(terms):
        total=c.mpf(0);bases=(c.mpf(1)/3,c.mpf(2),c.mpf(3),c.mpf(2))
        for t in terms:
            v=t['coefficient']
            for base,power in zip(bases,t['powers']):v*=base**power
            total+=v
        return total
    mp.mp.dps=75;hb=mp.mpf(1)/3;eps=mp.mpf(1)/2;width=mp.mpf(1)/8192
    def sigma(t):
        if t<=0:return mp.mpf(0)
        if t>=1:return mp.mpf(1)
        a=mp.exp(-1/t**2);b=mp.exp(-1/(1-t)**2);return a/(a+b)
    def integrand(q,order,branch):
        U=2+q;Uz=3;W=(1+2*q if branch=='core' else -1+q);Wz=(2 if branch=='core' else -1)
        product=U*W if order==0 else Uz*W+U*Wz
        return -2*eps**2*hb*width*16*mp.exp(2*hb*width*q)*(1-sigma(4*q-2))*product*2*(2*eps)
    checks={}
    mass=mp.quad(lambda q:1-sigma(4*q-2),[0,mp.mpf('.5'),mp.mpf('.75'),1]);checks['cutoff_mass']=abs(mass-mp.mpf(5)/8)<mp.mpf('1e-65')
    for order in (0,1):
        for branch in ('core','local'):
            terms=result['swirl'][order][branch];ref=c.mpf(mp.nstr(mp.quad(lambda q:integrand(q,order,branch),[0,mp.mpf('.5'),mp.mpf('.75'),1]),60));got=value(terms)
            checks[f'integral_{order}_{branch}']=bool(got.a<=ref.a and got.b>=ref.b)
            checks[f'physical_scales_{order}_{branch}']=all(t['powers'][0]==1 and t['powers'][1]==2 and t['powers'][3]==-3 for t in terms)
            checks[f'provenance_{order}_{branch}']=bool(terms) and all(t['origin']['source_ids'] for t in terms)
    checks['error_retained']=any('synthetic_error' in t['origin']['error_classes'] for fields in result['swirl'].values() for terms in fields.values() for t in terms)
    checks['tail_separate']=result['allword_tail_added'] is False
    checks['not_full_m15']=result['complete_m15_computed'] is False
    return {'passed':all(checks.values()),'checks':checks,'check_count':len(checks),'synthetic_only':True,'reference':'independent high precision direct cutoff integral of nonconstant q-dependent fields; numerical quadrature not certified'}
