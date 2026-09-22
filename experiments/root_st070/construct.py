"""Freeze a smooth eta-dependent annular control before independent diagnostics."""
from monotone_annulus import *
from datetime import datetime,timezone
import time

def moments(fam,co,e):
    f=fam.fixed(fam.x,e)-fam.K@co[:4]
    u=fam.axial_base(fam.x,e)+fam.B@co[4:]
    _,_,inner,target,scale=fam.base.setup_eta(e);w=fam.w;x=fam.x
    vals=inner+np.array([w@u,w@(2*x*f),w@(2*x*f*u),w@(u*u-x*f*f),w@(f*f)])
    mass=fam.fixed(4.,e)-sum(co[:4])-float(profile(4.,e,fam.c,fam.h,n=96)['F'])
    return (vals-target)/scale,mass

def run():
    t=time.monotonic();fam=MonotoneFamily(n=40)
    save('evidence/construction_registration.json',dict(utc=datetime.now(timezone.utc).isoformat(),c=fam.c,centers=fam.centers,widths=fam.widths,ucenters=fam.ucenters,uwidth=fam.uwidth,positive_density_floor=fam.floor,eta_levels=[17,33,65,129],n_per_segment=40,cdf_order=64,calibration_tolerance=1e-8,eta_calibration=np.linspace(-.5,.5,37),independent_seed=9227091,full_fd_seed=9227092,scope='Algebraic positive-derivative mixture and constrained axial norm construction; no PDE fit or wave construction'))
    records=[]
    for N in [17,33,65,129]:
        e=.5*np.cos(np.pi*np.arange(N)/(N-1));co=[];rows=[]
        for et in e:
            aa,rr=fam.solve_eta(float(et));co.append(aa);rows.append(rr)
        np.savez_compressed(ROOT/f'data/solved_nodes_N{N}.npz',eta=e,coefficients=np.array(co))
        C=ch.chebfit(2*e,np.array(co),N-1)
        test=[];seams=[]
        for et in np.linspace(-.5,.5,37):
            cc=ch.chebval(2*et,C);err,mass=moments(fam,cc,et);test.append(err);seams.append(mass)
        maxerr=float(np.max(abs(np.array(test))));rec=dict(nodes=N,max_calibration_error=maxerr,mass_seam_error=float(np.max(abs(np.asarray(seams)))),elapsed=time.monotonic()-t,rows=rows)
        records.append(rec);print('N',N,maxerr,rec['mass_seam_error'],rec['elapsed'],flush=True)
        meta=dict(id='ST070-S',schema='monotone_five_moment_strip_v1',c=fam.c,h=fam.h,nu=.01,Xc=.25,Xb=4.,eta_max=.5,pde_validated=False,global_field_ready=False,stress_realizable=False,source_correspondence_verified=False,blowup_proved=False,construction='Smooth nonnegative radial derivative mixture; five integral conditions, not dynamics')
        np.savez_compressed(ROOT/f'data/ST070_N{N}.npz',cheb=C,c=np.array(fam.c),h=np.array(fam.h),metadata=np.array(json.dumps(meta)))
        if maxerr<1e-8:
            import shutil
            shutil.copyfile(ROOT/f'data/ST070_N{N}.npz',ROOT/'data/ST070-S.npz');break
    else:raise RuntimeError('No passing interpolation level')
    save('evidence/construction.json',dict(records=records,elapsed=time.monotonic()-t,selected_nodes=N))
    save('evidence/freeze.json',dict(utc=datetime.now(timezone.utc).isoformat(),candidate_sha256=sha(ROOT/'data/ST070-S.npz'),source_sha256={name:sha(ROOT/name) for name in ['construct.py','monotone_annulus.py','moment_bound.py']},independent_eta_seed=9227091,independent_fd_seed=9227092,pde_validated=False))
if __name__=='__main__':run()
