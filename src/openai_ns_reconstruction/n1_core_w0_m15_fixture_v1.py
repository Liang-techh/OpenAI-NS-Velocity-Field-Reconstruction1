"""Synthetic variable-radius tests of physical core W0 m15 components."""
def run(leaf,c):
    import mpmath as mp
    def term(v,sid,powers=(0,0,0,0),error=False):return {'powers':list(powers),'coefficient':v,'origin':{'source_ids':[sid],'uncertainty':error,'error_classes':['synthetic_error'] if error else []}}
    U={0:[term(c.mpf([2,3]),'U')],1:[term(c.mpf(3),'UZ',error=True)]};F={0:[term(c.mpf([1,2]),'F',(0,1,0,0))],1:[term(c.mpf(2),'FZ',(0,1,0,0),True)]}
    outer=[]
    for i in range(8):
        t=c.mpf([i/8,(i+1)/8]);w={0:{'1':[term(2-t,f'F0:{i}',(0,1,0,-1),True)],'2':[term(1+t,f'U0:{i}')]},1:{'1':[term(c.mpf(1),f'F1:{i}',(0,1,0,-1),True)],'2':[term(c.mpf(2),f'U1:{i}')]}}
        outer.append({'t':t,'rho':4*t*t,'W0_ordinary_Z':w})
    g4={0:[term(c.mpf([1,2]),'g4:0')],1:[term(c.mpf([-2,-1]),'g4:1',error=True)]}
    atoms={k:'synthetic_'+k for k in ('hb','F0_anchor','Pstar_squared','Lambda')};cert={'g4_equals_two_epsilon_x_N1p':True,'geometry_epsilon_Z_independent':True,'actual_W0_outer_cells':True}
    out=leaf.build(c,outer_cells=outer,leading_U=U,leading_F=F,g4=g4,source_family='synthetic',source_atoms=atoms,certificate=cert)
    def value(terms):
        total=c.mpf(0);bases=(c.mpf(1)/3,c.mpf(2),c.mpf(3),c.mpf(2))
        for t in terms:
            v=t['coefficient']
            for base,p in zip(bases,t['powers']):v*=base**p
            total+=v
        return total
    mp.mp.dps=75;eps=mp.mpf(1)/2;checks={}
    for order in (0,1):
        for i in range(8):
            lo=mp.mpf(i)/8;hi=mp.mpf(i+1)/8
            def source(t,channel):
                if channel=='axial':
                    U0=2+t;Uz=3;W0=1+t;Wz=2;product=U0*W0 if order==0 else Uz*W0+U0*Wz
                    return eps*16*t*product
                F0=2*(1+t);Fz=4;W0=2-t;Wz=1;product=F0*W0 if order==0 else Fz*W0+F0*Wz
                return -eps**2*64*t**3*product
            for channel in ('axial','swirl'):
                terms=out[channel][order][i];ref=c.mpf(mp.nstr(mp.quad(lambda t:source(t,channel),[lo,hi]),60));got=value(terms)
                checks[f'{channel}_{order}_{i}_integral']=bool(got.a<=ref.a and got.b>=ref.b)
                checks[f'{channel}_{order}_{i}_ids']=bool(terms) and all(x['origin']['source_ids'] for x in terms)
                checks[f'{channel}_{order}_{i}_physical_power']=all(x['powers'][0]==0 and x['powers'][3]==(-1 if channel=='axial' else -3) for x in terms)
        ref=c.mpf(mp.nstr(mp.quad(lambda t:-8*eps*t*t*(1+t if order==0 else -2+t),[0,1]),60));got=value(out['pressure'][order]);checks[f'pressure_{order}']=bool(got.a<=ref.a and got.b>=ref.b)
    checks['no_tail']=out['allword_tail_added'] is False;checks['no_W1']=out['W1_included'] is False;checks['no_full_m15']=out['complete_m15_computed'] is False
    return {'passed':all(checks.values()),'check_count':len(checks),'checks':checks,'synthetic_only':True,'reference':'independent high-precision variable-radius integrals; no physical-source proof'}
