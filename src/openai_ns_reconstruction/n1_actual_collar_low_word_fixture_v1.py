"""Independent synthetic direct-integral checks of local ordinary Z0/Z1."""
def run(leaf,c):
    import mpmath as mp
    fields=('F1','Uz1','K1','P1','d_x_F1','d_x_Uz1');rows=('5','6','6','4','5','6')
    ckeys={'4,1','5,1','5,2','5,3','5,5','6,1','6,2','6,3','6,4','6,6'}
    bkeys={'5,1','5,2','5,3','6,2','6,3','6,4'}
    def value(j,i):return (j+2)*(i+1)*(1 if (j+i)%2==0 else -1)
    def sec(j,i):return {'powers':[0,0,0,0],'coefficient':c.mpf(value(j,i))/17,'origin':{'source_ids':[f'synthetic:{j}:{i}'],'uncertainty':j%2==0,'error_classes':['synthetic_interval'] if j%2==0 else []}}
    core={j:{field:[sec(j,i)] for i,field in enumerate(fields)} for j in range(4)}
    g={j:{r:[sec(j+4,int(r))] for r in ('4','5','6')} for j in range(4)}
    present={'C0':{'4,1','5,1','6,2'},'B1':{'5,1','6,2'}}
    matrix={name:{j:{key:([sec(j+(10 if name=='C0' else 13),i)] if key in present[name] else []) for i,key in enumerate(sorted(keys))} for j in range(3)} for name,keys in (('C0',ckeys),('B1',bkeys))}
    sc=c.mpf(1)/2048;end=sc/4;atoms={x:'synthetic_'+x for x in ('hb','F0_anchor','Pstar_squared','Lambda')};geometry={'hb_strictly_between_zero_and_half':True,'sc_phase_Z_independent':True,'hb_epsilon_Z_independent':True}
    out=leaf.enclose_low(c,sc=sc,phase=end,core_W0_Z=core,forcing_Z=g,matrix_Z=matrix,inner_cells=2,outer_cells=2,source_family='synthetic',source_atoms=atoms,geometry_certificate=geometry)
    phase=c.make_mpf((c.mpf(0)._mpi_[0],end.b._mpi_[1]))
    entire=leaf.enclose_low(c,sc=sc,phase=phase,core_W0_Z=core,forcing_Z=g,matrix_Z=matrix,inner_cells=1,outer_cells=1,source_family='synthetic',source_atoms=atoms,geometry_certificate=geometry)
    def evalterms(terms):
        total=c.mpf(0)
        for term in terms:total+=term['coefficient']*(c.mpf(1)/3)**term['powers'][0]
        return total
    mp.mp.dps=75;hb=mp.mpf(1)/3;top=mp.mpf(1)/8192
    def val(j,i):return mp.mpf(value(j,i))/17
    def word(v,j):
        d=hb*v;e=mp.exp(-d);f,u,k,p,fx,ux=[val(j,i) for i in range(6)]
        t=(f+(1-e)*fx,u+d*ux,e*k-(1-e)*ux,p,mp.exp(-3*d/2)*fx,mp.exp(-d/2)*ux)
        x=2*mp.exp(d/2);a=(x*x-4)/4-(x**4-16)/(8*x*x);b=x*x/4-2*mp.log(x/2)-1
        masses=(a,b,-a,x-2,(x**4-16)/(4*x**3),(x*x-4)/(2*x))
        return [t[i]+masses[i]*val(j+4,int(row)) for i,row in enumerate(rows)]
    def source(v,order):
        w=[word(v,j) for j in range(3)];out={r:mp.mpf(0) for r in ('4','5','6')}
        for name,keys in (('C0',ckeys),('B1',bkeys)):
            for i,key in enumerate(sorted(keys)):
                if key not in present[name]:continue
                row,col=key.split(',');col=int(col)-1;shift=0 if name=='C0' else 1;offset=10 if name=='C0' else 13
                # Coefficient of z^order in a product of Taylor polynomials.
                out[row]+=sum(val(offset+j,i)*w[order-j+shift][col] for j in range(order+1))
        return out
    def integrand(v,i,order):
        u=2*mp.exp(hb*v/2);d=hb*(top-v);ratio=mp.exp(-d/2);k=u*(1-mp.exp(-d))/2
        kernels=(k,u*d/2,-k,mp.mpf(1),ratio**3,ratio)
        return source(v,order)[rows[i]]*kernels[i]*hb*u/2
    checks={}
    for order in (0,1):
        for i,field in enumerate(fields):
            ref=c.mpf(mp.nstr(mp.quad(lambda v:integrand(v,i,order),[0,top]),60));got=evalterms(out['local_signed_sectors'][order][field]);whole=evalterms(entire['local_signed_sectors'][order][field])
            checks[f'direct_{order}_{field}']=bool(got.a<=ref.a and got.b>=ref.b)
            checks[f'full_{order}_{field}']=bool(whole.a<=ref.a and whole.b>=ref.b)
            checks[f'provenance_{order}_{field}']=all(t['origin']['source_ids'] for t in out['local_signed_sectors'][order][field])
    checks['fields_and_rows']=leaf.FIELDS==fields and leaf.ROWS==rows
    checks['separate_missing_core_inlet']=out['core_W1_inlet_included'] is False
    checks['tail_not_added']=out['allword_tail_added'] is False
    return {'passed':all(checks.values()),'checks':checks,'check_count':len(checks),'synthetic_only':True,'reference':'independent high-precision direct integration and polynomial product coefficients; quadrature is not certified'}
