"""Independent synthetic physical pressure-density integration check."""
def run(leaf,c):
    import mpmath as mp
    def term(box,sid):return {'powers':[0,0,0,0],'coefficient':c.mpf(box),'origin':{'source_ids':[sid],'uncertainty':True,'error_classes':['synthetic_source']}}
    g4={0:[term([1,2],'g40')],1:[term([-2,-1],'g41')]};sc=c.mpf(1)/2048
    atoms={k:'synthetic_'+k for k in ('hb','F0_anchor','Pstar_squared','Lambda')}
    cert={'hb_strictly_between_zero_and_half':True,'sc_phase_Z_independent':True,'hb_epsilon_Z_independent':True,'g4_equals_two_epsilon_x_N1p':True}
    result=leaf.build(c,sc=sc,g4=g4,source_family='synthetic',source_atoms=atoms,certificate=cert)
    def value(terms):
        total=c.mpf(0);bases=(c.mpf(1)/3,c.mpf(2),c.mpf(3),c.mpf(2))
        for t in terms:
            v=t['coefficient']
            for base,power in zip(bases,t['powers']):v*=base**power
            total+=v
        return total
    mp.mp.dps=75;hb=mp.mpf(1)/3;eps=mp.mpf(1)/2;S=mp.mpf(1)/8192
    checks={}
    for order in (0,1):
        def integrand(q):
            x=2*mp.exp(hb*S*q/2);g=(1+q if order==0 else -2+q)
            return -eps*hb*S*x**3*g/2
        ref=c.mpf(mp.nstr(mp.quad(integrand,[0,1]),60));got=value(result['pressure'][order])
        checks[f'direct_integral_{order}']=bool(got.a<=ref.a and got.b>=ref.b)
        checks[f'physical_factor_{order}']=all(t['powers']==(1,0,0,-1) or t['powers']==[1,0,0,-1] for t in result['pressure'][order])
        checks[f'provenance_{order}']=all('g4'+str(order) in t['origin']['source_ids'] for t in result['pressure'][order])
    checks['uncut_last_cell']=result['velocity_cutoff_applied'] is False
    checks['no_endpoint_duplication']=result['explicit_pressure_endpoint_added'] is False
    checks['no_allword']=result['allword_tail_added'] is False
    checks['not_complete']=result['complete_m15_computed'] is False
    return {'passed':all(checks.values()),'checks':checks,'check_count':len(checks),'synthetic_only':True,'reference':'independent high-precision physical g4 pressure-density integral; quadrature is not certified'}
