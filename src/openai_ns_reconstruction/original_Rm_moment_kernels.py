"""Positive native patch kernels for beta and survival-tail pressure modes."""
from fractions import Fraction
import hashlib

from .original_core_numeric import ROOT,_bounds


def powers(code):
    result=[]
    for i in range(4):
        result.extend([i]*(code%3));code//=3
    if code or len(result)>2:raise ValueError('Only original linear/quadratic four-shape modes')
    return result


class OriginalRmMomentKernels:
    def __init__(self,rows,*,cells=32,partial_cells=128):
        if type(cells) is not int or not 8<=cells<=128:raise ValueError('Use8..128 positive native cells')
        if type(partial_cells) is not int or not 64<=partial_cells<=512:raise ValueError('Use64..512 partial bump cells')
        self.rows=rows;self.patch=rows.patch;self.ctx,self.iv,self.audit=rows.ctx,rows.iv,rows.audit
        self.cells,self.partial_cells=cells,partial_cells
        self.centers=(Fraction(5,4),Fraction(7,4));self._cdfs={};self._moments={};self._packets={}
        name='src/openai_ns_reconstruction/original_Rm_moment_kernels.py'
        self.code_hashes={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()}

    @staticmethod
    def number(ctx,q):return ctx.mpf(q.numerator)/q.denominator

    def cdf(self,i,q):
        key=i,q
        if key not in self._cdfs:
            c,iv,p=self.ctx,self.iv,self.patch;center=self.centers[i];a,b=center-Fraction(1,40),center+Fraction(1,40)
            if q<=a:out=c.mpf(0),iv.mpf(0)
            elif q>=b:out=c.mpf(1),iv.mpf(1)
            else:
                out=p._point_weight(self.number(c,q),2*i,'0',1),p.bump.partial_weight(
                    self.number(iv,q),center,'0',1,cells=self.partial_cells)
            self._cdfs[key]=out
        return self._cdfs[key]

    def beta(self,ctx,i,x,directed):
        center=self.number(ctx,self.centers[i])
        if directed:return self.patch.bump.beta(x-center)
        u=(x-center)/ctx.mpf('.025')
        if abs(u)>=1:return ctx.mpf(0),ctx.mpf(0)
        value=ctx.exp(-1/(1-u*u))/(ctx.mpf('.025')*self.patch.normp)
        return value,-2*u*value/(ctx.mpf('.025')*(1-u*u)**2)

    def tail(self,i,q):
        p,t=self.cdf(i,q);lo,hi=_bounds(self.audit,1-t)
        return 1-p,self.iv.mpf([max(0,lo),min(1,hi)])

    def moments(self,i,q):
        """A=integral x beta; B=integral x beta²; J=integral T; L/S=integral T/T² dx/x."""
        key=i,q
        if key in self._moments:return self._moments[key]
        c,iv=self.ctx,self.iv;center=self.centers[i];a,b=center-Fraction(1,40),center+Fraction(1,40)
        end=b if q is None or q>=b else q
        if end<=a:
            xp,xt=self.number(c,end),self.number(iv,end)
            out=dict(point=dict(A=c.mpf(0),B=c.mpf(0),J=xp-1,L=c.log(xp),S=c.log(xp)),
                     truth=dict(A=iv.mpf(0),B=iv.mpf(0),J=xt-1,L=iv.log(xt),S=iv.log(xt)),cells=0)
        else:
            ap,at=self.number(c,a),self.number(iv,a)
            point=dict(A=c.mpf(0),B=c.mpf(0),J=ap-1,L=c.log(ap),S=c.log(ap))
            truth=dict(A=iv.mpf(0),B=iv.mpf(0),J=at-1,L=iv.log(at),S=iv.log(at));used=0
            for k in range(self.cells):
                left=a+(b-a)*Fraction(k,self.cells);right=min(end,a+(b-a)*Fraction(k+1,self.cells))
                if right<=left:break
                midpoint=(left+right)/2;xp=self.number(c,midpoint)
                xl,xr=self.number(iv,left),self.number(iv,right)
                box=iv.mpf([_bounds(self.audit,xl)[0],_bounds(self.audit,xr)[1]])
                width=self.number(iv,right-left);wp=self.number(c,right-left)
                bp,_=self.beta(c,i,xp,False);bt,_=self.beta(iv,i,box,True)
                tp,_=self.tail(i,midpoint);_,tl=self.tail(i,left);_,tr=self.tail(i,right)
                tailbox=iv.mpf([_bounds(self.audit,tr)[0],_bounds(self.audit,tl)[1]])
                logwidth=iv.log(xr/xl);logpoint=c.log(self.number(c,right)/self.number(c,left))
                point['A']+=xp*bp*wp;point['B']+=xp*bp*bp*wp;point['J']+=tp*wp
                point['L']+=tp*logpoint;point['S']+=tp*tp*logpoint
                truth['A']+=box*bt*width;truth['B']+=box*bt*bt*width;truth['J']+=tailbox*width
                truth['L']+=tailbox*logwidth;truth['S']+=tailbox*tailbox*logwidth;used+=1
            if end==b:
                # Exact normalization and symmetry, including beta² symmetry.
                cp,ct=self.number(c,center),self.number(iv,center)
                point['A'],truth['A']=cp,ct;point['J'],truth['J']=cp-1,ct-1
                point['B']=cp*self.patch._point_weight(c.mpf(2),2*i,'0',2)
                truth['B']=ct*self.patch.bump.full_weights[center,Fraction(0),2]
            out=dict(point=point,truth=truth,cells=used)
        self._moments[key]=out;return out

    def evaluate(self,x):
        if x in self._packets:return self._packets[x]
        geo=self.rows.geometry(x);q=None if x=='Rh' else Fraction(x);c,iv=self.ctx,self.iv
        leaf=[self.moments(i,q) for i in range(2)];endpoint=[]
        for i in range(2):
            if q is None:Ip,It,Tp,Tt=c.mpf(1),iv.mpf(1),c.mpf(0),iv.mpf(0)
            else:(Ip,It),(Tp,Tt)=self.cdf(i,q),self.tail(i,q)
            endpoint.append(dict(point=dict(I=Ip,T=Tp),truth=dict(I=It,T=Tt)))
        result={}
        for label,ctx,xx in (('point',c,geo['x_point']),('truth',iv,geo['x'])):
            values={0:(xx*xx-1)/2};L=[p[label] for p in leaf];E=[p[label] for p in endpoint]
            for i in range(2):values[3**i]=L[i]['A'];values[3**(i+2)]=L[i]['J']
            for i in range(4):
                for j in range(i,4):
                    if i==j and i<2:value=L[i]['B']
                    elif (i,j)==(0,1) or (i,j)==(1,2):value=ctx.mpf(0)
                    elif (i,j)==(0,2):value=(1-E[0]['T']*E[0]['T'])/2
                    elif (i,j)==(1,3):value=(1-E[1]['T']*E[1]['T'])/2
                    elif (i,j)==(0,3):value=E[0]['I']
                    elif i==j:value=L[i-2]['S']
                    elif (i,j)==(2,3):value=L[0]['L']
                    else:raise AssertionError('Unmapped native shape product')
                    values[3**i+3**j]=value
            result[label]=values
        result.update(geometry=geo,endpoint=endpoint,leaf=leaf,
            contract=dict(physical_measure='R*dR=Rm²*x*dx',four_shapes=['beta0','beta2','T0/x','T2/x'],
                complete_positive_prefix_from_one=True,original_normalized_bump_and_survival_tail=True,
                disjoint_support_cross_identities=True,full_symmetry_moments_exact=True,
                tail_log_moments_direct_positive_cells=True,point_quadrature_only_proxy=True,
                original_fixed_control_scale_not_epsilon=True))
        self._packets[x]=result;return result
