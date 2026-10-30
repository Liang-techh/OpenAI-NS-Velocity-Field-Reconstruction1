"""Synthetic checks for no-reset core transport plus local collar terms."""
def run(leaf,c):
    import mpmath as mp
    from fractions import Fraction
    fields=('F1','Uz1','K1','P1','d_x_F1','d_x_Uz1')
    core={order:{field:[{'F0_power':0,'Lambda_power':0,'coefficient':c.mpf((order+2)*(i+1))/13,'input_id':f'core:{order}:{field}:0','outer_cell':i}] for i,field in enumerate(fields)} for order in (0,1)}
    local={order:{field:[{'powers':[1,0,0,0],'coefficient':c.mpf((-1 if i%2 else 1)*(i+order+1))/19,'origin':{'source_ids':[f'local:{order}:{field}'],'uncertainty':True,'error_classes':['synthetic_local']}}] for i,field in enumerate(fields)} for order in (0,1)}
    sc=c.mpf(1)/2048;phase=sc/4;atoms={k:'synthetic_'+k for k in ('hb','F0_anchor','Pstar_squared','Lambda')};geometry={'hb_strictly_between_zero_and_half':True,'sc_phase_Z_independent':True,'hb_epsilon_Z_independent':True}
    out=leaf.build(c,sc=sc,phase=phase,core=core,local=local,source_family='synthetic',source_atoms=atoms,geometry_certificate=geometry)
    fullphase=c.make_mpf((c.mpf(0)._mpi_[0],phase.b._mpi_[1]));whole=leaf.build(c,sc=sc,phase=fullphase,core=core,local=local,source_family='synthetic',source_atoms=atoms,geometry_certificate=geometry)
    def val(terms):
        total=c.mpf(0)
        for term in terms:total+=term['coefficient']*(c.mpf(1)/3)**term['powers'][0]
        return total
    mp.mp.dps=80;hb=mp.mpf(1)/3;s=mp.mpf(1)/8192
    def transport(order):
        f,u,k,p,fx,ux=[mp.mpf((order+2)*(i+1))/13 for i in range(6)];e=mp.exp(-hb*s)
        return (f+(1-e)*fx,u+hb*s*ux,e*k-(1-e)*ux,p,mp.exp(-3*hb*s/2)*fx,mp.exp(-hb*s/2)*ux)
    checks={}
    for order in (0,1):
        reference=transport(order)
        for i,field in enumerate(fields):
            local_value=hb*mp.mpf((-1 if i%2 else 1)*(i+order+1))/19
            if i in (1,3):
                rational=Fraction((order+2)*(i+1),13)+Fraction((-1 if i%2 else 1)*(i+order+1),57)
                if i==1:rational+=Fraction(1,3*8192)*Fraction((order+2)*6,13)
                ref=c.mpf(rational.numerator)/rational.denominator
            else:ref=c.mpf(mp.nstr(reference[i]+local_value,65))
            value=val(out['finite_sum'][order][field]);wide=val(whole['finite_sum'][order][field])
            checks[f'transport_{order}_{field}']=bool(value.a<=ref.a and value.b>=ref.b)
            checks[f'full_phase_{order}_{field}']=bool(wide.a<=ref.a and wide.b>=ref.b)
            checks[f'local_once_{order}_{field}']=len(out['local_contribution'][order][field])==1 and out['local_contribution'][order][field][0]['origin']['source_ids']==[f'local:{order}:{field}']
            checks[f'core_id_{order}_{field}']=any(f'core:{order}:{field}:0' in t['origin']['source_ids'] for t in out['core_transport'][order][field])
    checks['tail_separate']=out['allword_tail_added'] is False
    checks['no_m15']=out['m15_computed'] is False
    return {'passed':all(checks.values()),'checks':checks,'check_count':len(checks),'synthetic_only':True,'reference':'independent high precision exact homogeneous transport and explicit local term; no physical-source identity assumed'}
