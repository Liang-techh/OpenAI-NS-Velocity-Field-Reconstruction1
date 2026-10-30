"""Independent nonconstant synthetic W0 m15 integrals by phase cell."""
def run(leaf,c):
    import mpmath as mp
    def term(value,sid,powers=(0,0,0,0),uncertain=False):
        return {'powers':list(powers),'coefficient':value,'origin':{'source_ids':[sid],'uncertainty':uncertain,'error_classes':['synthetic_error'] if uncertain else []}}
    U={0:[term(c.mpf([2,3]),'U')],1:[term(c.mpf(3),'UZ',uncertain=True)]}
    F={0:[term(c.mpf([1,2]),'F',(0,1,0,0))],1:[term(c.mpf(2),'FZ',(0,1,0,0),True)]}
    w={0:{'Uz1':{'core':[term(c.mpf([1,2]),'Uc')],'local':{}},'F1':{'core':[term(c.mpf([1,2]),'Fc',(0,1,0,-1))],'local':{}}},
       1:{'Uz1':{'core':[term(c.mpf(2),'Ucz')],'local':{}},'F1':{'core':[term(c.mpf(1),'Fcz',(0,1,0,-1))],'local':{}}}}
    for i in range(4):
        q=c.mpf([i/4,(i+1)/4]);w[0]['Uz1']['local'][i]=[term(-1+q,f'Ul{i}',uncertain=True)]
        w[1]['Uz1']['local'][i]=[term(c.mpf(-1),f'Ulz{i}',uncertain=True)]
        w[0]['F1']['local'][i]=[term(q,f'Fl{i}',(0,1,0,-1),True)]
        w[1]['F1']['local'][i]=[term(c.mpf(1),f'Flz{i}',(0,1,0,-1),True)]
    sc=c.mpf(1)/2048;atoms={k:'synthetic_'+k for k in ('hb','F0_anchor','Pstar_squared','Lambda')}
    cert={'hb_strictly_between_zero_and_half':True,'sc_phase_Z_independent':True,'hb_epsilon_Z_independent':True,'original_cutoff_mass_5_over_8':True,'original_cutoff_mass_by_cell':True}
    out=leaf.build(c,sc=sc,leading_U=U,leading_F=F,w0=w,source_family='synthetic',source_atoms=atoms,certificate=cert)
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
    def integrand(q,order,channel,branch):
        eta=1-sigma(4*q-2)
        if channel=='axial':
            lead=2+q;leadZ=3;word=(1+q if branch=='core' else -1+q);wordZ=(2 if branch=='core' else -1)
            product=lead*word if order==0 else leadZ*word+lead*wordZ
            return eps*hb*width*8*mp.exp(hb*width*q)*eta*product
        lead=2*(1+q);leadZ=4;word=(2-q if branch=='core' else q);wordZ=1
        product=lead*word if order==0 else leadZ*word+lead*wordZ
        return -eps**2*hb*width*32*mp.exp(2*hb*width*q)*eta*product
    checks={}
    for channel in ('axial','swirl'):
        for order in (0,1):
            for branch in ('core',0,1,2):
                terms=out[channel][order]['core'] if branch=='core' else out[channel][order]['local'][branch]
                limits=([0,mp.mpf('.5'),mp.mpf('.75'),1] if branch=='core' else [mp.mpf(branch)/4,mp.mpf(branch+1)/4])
                ref=c.mpf(mp.nstr(mp.quad(lambda q:integrand(q,order,channel,'core' if branch=='core' else 'local'),limits),60));got=value(terms)
                checks[f'integral_{channel}_{order}_{branch}']=bool(got.a<=ref.a and got.b>=ref.b)
                checks[f'physical_{channel}_{order}_{branch}']=all(t['powers'][0]==1 and t['powers'][3]==(-1 if channel=='axial' else -3) for t in terms)
                checks[f'ids_{channel}_{order}_{branch}']=bool(terms) and all(t['origin']['source_ids'] for t in terms)
            checks[f'cutoff_{channel}_{order}']=set(out[channel][order]['local'])=={0,1,2}
    checks['no_tail']=out['allword_tail_added'] is False
    checks['no_full_m15']=out['complete_m15_computed'] is False
    return {'passed':all(checks.values()),'checks':checks,'check_count':len(checks),'synthetic_only':True,'reference':'independent high precision physical integrals of nonconstant phase fields; quadrature not certified'}
