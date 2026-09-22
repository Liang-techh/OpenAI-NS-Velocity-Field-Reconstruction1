"""ST064-S: scale-native kinematic scaffold, NOT a validated NS solution.
Physical evolution is evaluated, not camera-scaled. Its exponents/profile are
prescribed; this does NOT demonstrate dynamically generated scale recurrence.
"""
from pathlib import Path
import hashlib,json
import numpy as np
from profile import Profile,bump
ROOT=Path(__file__).resolve().parent

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
class ScaleField:
    def __init__(self,config=None):
        self.config=json.loads((ROOT/'data/model.json').read_text()) if config is None else dict(config)
        c=self.config
        if c['pde_validated'] or c['dynamic_recursion_validated']:raise ValueError('Unsupported scientific promotion')
        pp=ROOT/'data/profile.npz'
        if sha(pp)!=c['profile_sha256']:raise ValueError('Profile hash mismatch')
        a=np.load(pp,allow_pickle=False);self.profile=Profile(a['F'],a['G'],a['P'])
        self.h=c['h'];self.nu=c['nu'];self.tau0=c['tau0'];self.T=c['T'];self.kmax=c['kmax'];self.force_params=np.array(c['force'])
    def scales(self,k):
        k=float(k)
        if not np.isfinite(k) or not 0<=k<=self.kmax:raise ValueError('Scale outside frozen window')
        q=2.**(-k);ar=q**.5;az=q**(.5-self.h);tau=self.tau0*q
        return dict(k=k,q=q,tau=tau,t=self.T-tau,ar=ar,az=az,ur_amp=q**(-.5),utheta_amp=q**(-.5-self.h),uz_amp=q**(-.5-self.h))
    def k_from_time(self,t):
        tau=self.T-float(t)
        if not np.isfinite(tau) or tau<=0:raise ValueError('Endpoint not in domain')
        k=np.log2(self.tau0/tau)
        if k<-1e-12 or k>self.kmax+1e-12:raise ValueError('Time outside NEW window')
        return float(np.clip(k,0,self.kmax))
    def cylindrical(self,s,z,k):
        sc=self.scales(k);q,ar,az,tau=[sc[n] for n in ('q','ar','az','tau')]
        s,z=np.broadcast_arrays(np.asarray(s,float),np.asarray(z,float));s=s.ravel();z=z.ravel();ss=s/ar**2;zz=z/az
        d=self.profile.jets(ss,zz);fa=q**-1;fb=q**(-1-self.h);fc=az/q;fp=q**(-1-2*self.h)
        D={}
        for tag,fac in (('A',fa),('B',fb),('C',fc),('P',fp)):
            D[tag]=fac*d[tag]
            for suf,den in (('s',ar**2),('z',az),('ss',ar**4),('zz',az**2)):
                if tag+suf in d:D[tag+suf]=fac/den*d[tag+suf]
        for tag,fac,power in (('A',fa,1),('B',fb,1+self.h),('C',fc,.5+self.h)):
            D[tag+'t']=fac/tau*(power*d[tag]+ss*d[tag+'s']+(.5-self.h)*zz*d[tag+'z'])
        br=bump(s/4);bz=bump(z*z/4);b=br[:,0]*bz[:,0];bs=br[:,1]*bz[:,0]/4;b_z=br[:,0]*bz[:,1]*z/2
        a,c=self.force_params;g=float(bump(np.array([(2*sc['t']-1)**2]))[0,0])
        D['fa']=-g*a*(b+z*b_z);D['fb']=g*c*(b+s*bs);D['fc']=2*g*a*z*(b+s*bs)
        A,B,C=D['A'],D['B'],D['C']
        D['AL']=8*D['As']+4*s*D['Ass']+D['Azz'];D['BL']=8*D['Bs']+4*s*D['Bss']+D['Bzz'];D['CL']=4*D['Cs']+4*s*D['Css']+D['Czz']
        D['Ra']=D['At']+A*A+2*s*A*D['As']+C*D['Az']-B*B+2*D['Ps']-self.nu*D['AL']-D['fa']
        D['Rb']=D['Bt']+2*A*B+2*s*A*D['Bs']+C*D['Bz']-self.nu*D['BL']-D['fb']
        D['Rc']=D['Ct']+2*s*A*D['Cs']+C*D['Cz']+D['Pz']-self.nu*D['CL']-D['fc']
        D['div']=2*A+2*s*D['As']+D['Cz'];return D
    def evaluate(self,points,k):
        x=np.asarray(points,float)
        if x.ndim!=2 or x.shape[1]!=3 or not np.isfinite(x).all():raise ValueError('Finite Nx3 points required')
        if len(x)==0:return {n:np.empty((0,3)) if n!='p' else np.empty(0) for n in ['u','p','omega','residual','force','ut','gradp','lap']}
        outputs=[]
        for start in range(0,len(x),2048):
            X,Y,Z=x[start:start+2048].T;s=X*X+Y*Y;D=self.cylindrical(s,Z,k)
            def vec(a,b,c):return np.c_[X*D[a]-Y*D[b],Y*D[a]+X*D[b],D[c]]
            outputs.append(dict(u=vec('A','B','C'),p=D['P'],omega=np.c_[-X*D['Bz']-Y*(D['Az']-2*D['Cs']),-Y*D['Bz']+X*(D['Az']-2*D['Cs']),2*(D['B']+s*D['Bs'])],
                residual=vec('Ra','Rb','Rc'),force=vec('fa','fb','fc'),ut=vec('At','Bt','Ct'),lap=vec('AL','BL','CL'),gradp=np.c_[2*X*D['Ps'],2*Y*D['Ps'],D['Pz']],div=D['div']))
        return {n:np.concatenate([o[n] for o in outputs]) for n in outputs[0]}
    def velocity(self,points,k):
        pts=np.asarray(points,float);s0=self.scales(k);q=s0['q'];ar=s0['ar'];az=s0['az'];X,Y,Z=pts.T
        a,b,c=self.profile.abc((X*X+Y*Y)/ar**2,Z/az)
        A=a/q;B=b*q**(-1-self.h)
        return np.c_[X*A-Y*B,Y*A+X*B,az/q*c]
    def fields_time(self,x,t):
        D=self.evaluate(x,self.k_from_time(t));return D['u'],D['p']
    def forcing(self,x,t):
        X,Y,Z=np.asarray(x).T;s=X*X+Y*Y;br=bump(s/4);bz=bump(Z*Z/4);b=br[:,0]*bz[:,0];bs=br[:,1]*bz[:,0]/4;bzv=br[:,0]*bz[:,1]*Z/2
        a,c=self.force_params;g=bump(np.array([(2*t-1)**2]))[0,0];A=-a*(b+Z*bzv);B=c*(b+s*bs);C=2*a*Z*(b+s*bs)
        return g*np.c_[X*A-Y*B,Y*A+X*B,C]
