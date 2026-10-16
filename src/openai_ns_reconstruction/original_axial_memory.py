"""Positive original O2-axial convolution memories on exponential radius scales.

Peak normalization resolves the nonzero endpoint memory even when the local
axial cutoff vanishes. No physical-radius exponential or float field is used.
"""
from .original_core_numeric import _bounds


class OriginalAxialMemory:
    def __init__(self,ctx,iv,J,*,relative_target='1e-8'):
        from lei_ren_part1_paper_compliant_flat_pulse_derivatives import sigma_jets
        self.ctx,self.iv,self.J,self.sigma_jets=ctx,iv,J,sigma_jets
        self.target=ctx.mpf(relative_target);self.cache={};self.roots={}
        if not 0<self.target<ctx.mpf('.01'):raise ValueError('Positive small convolution target required')
        a=iv.mpf('1e-4');margin=6*a**4/(1-a)**4+80*a**4/(1-a)**3+80*a-6
        if _bounds(ctx,margin)[1]>=0:raise ArithmeticError('Small-r log concavity required')
        bounds=[]
        for lo,hi in (('1e-4','.001'),('.001','.01'),('.01','.1'),('.1','.5')):
            l,r=iv.mpf(lo),iv.mpf(hi)
            bounds.append((2/(1-r)**3+2/l**3)*iv.exp(-40*(1-r))/40)
        for lo,hi in (('.5','.75'),('.75','1')):
            derivative=self.sigma_jets(iv,iv.mpf([ctx.mpf(lo),ctx.mpf(hi)]))[1]
            upper=max(abs(v) for v in _bounds(ctx,derivative))
            bounds.append(2*iv.mpf(upper)*iv.exp(-40*(1-iv.mpf(hi)))/40)
        self.outside_t_upper=max(_bounds(ctx,v)[1] for v in bounds)
        if self.outside_t_upper>=ctx.mpf('.25'):raise ArithmeticError('Global decreasing log-tail guard required')
        self.global_proof=dict(small_r_concavity_margin=margin,small_r_boundary='1e-4',
            outside_log_beta_q_derivative_upper=self.outside_t_upper,
            all_powers_1_2_log_integrand_slope_le_minus_half_outside_small_r=True)

    def _t(self,r):
        iv=self.iv;odds=1/(1-r)**2-1/r**2
        return (2/(1-r)**3+2/r**3)*iv.exp(-40*(1-r))/(40*(1+iv.exp(odds)))

    def _root(self,power):
        if power in self.roots:return self.roots[power]
        c,iv=self.ctx,self.iv;lo=c.mpf(0);hi=c.mpf('1e-4')
        for _ in range(220):
            mid=(lo+hi)/2;v=self._t(iv.mpf(mid))*power
            a,b=_bounds(c,v)
            if a>1:lo=mid
            elif b<1:hi=mid
            else:break
        out=iv.mpf([lo,hi]);self.roots[power]=out;return out

    def _ln(self,j):
        ctx=j.ctx;n=j.order
        if not n:return self.J(ctx,[ctx.ln(j[0])])
        derivative=self.J(ctx,[(k+1)*j[k+1] for k in range(n)])
        quotient=derivative/j.truncate(n-1)
        return self.J(ctx,[ctx.ln(j[0])]+[quotient[k-1]/k for k in range(1,n+1)])

    def _log_beta(self,ctx,r):
        # The checked sigmoid jet treats flat endpoints without singular
        # reciprocals. On the left, log-sigmoid keeps huge negative odds.
        lo=_bounds(self.ctx,r[0])[0] if ctx is self.iv else r[0]
        if ctx is self.iv and lo>=self.ctx.mpf('.25'):
            sj=self.sigma_jets(ctx,r[0]);zero=self.J(ctx,[ctx.mpf(0)]+list(r.coefficients[1:]))
            result=self.J.constant(ctx,0,r.order)
            for k in range(r.order+1):result+=zero**k*sj[k]
            return self._ln(result)
        odds=(1-r)**-2-r**-2
        if lo>=self.ctx.mpf('.5'):return -self._ln(1+(-odds).exp())
        return odds-self._ln(1+odds.exp())

    def evaluate(self,s,power):
        if power not in (1,2):raise ValueError('Original beta power1/2 required')
        key=(s,power)
        if key in self.cache:return self.cache[key]
        c,iv,J=self.ctx,self.iv,self.J;sp=c.mpf(s);si=iv.mpf(s)
        if not c.isfinite(sp) or not 0<=sp<=1:raise ValueError('Original axial native0<=s<=1 required')
        dp=c.expm1(40*sp);dt=iv.expm1(40*si);U=c.exp(40*sp);Ui=iv.exp(40*si)
        if sp==0:
            out=dict(point=c.mpf(1),directed=iv.mpf(1),proof=dict(zero_width_actual_history=True,power=power))
            self.cache[key]=out;return out
        r0=1-sp;r0i=1-si;root=self._root(power);rl,rh=_bounds(c,root);rs=(rl+rh)/2
        ordinary=_bounds(c,r0i)[0]>=c.mpf('1e-4') or r0>0 and _bounds(c,self._t(r0i)*power)[1]<=c.mpf('.5')
        def logs(ctx,z,qorigin,scale,order=4):
            u=(ctx.mpf(U) if ctx is c else Ui)-qorigin-scale*z
            initial=ctx.mpf(r0) if ctx is c else r0i
            rvalue=initial-ctx.ln(1-(qorigin+scale*z)/(ctx.mpf(U) if ctx is c else Ui))/40
            ratio=scale/u
            r=J(ctx,[rvalue]+[ratio**k/(40*k) for k in range(1,order+1)])
            return self._log_beta(ctx,r)*power
        if ordinary:
            origin=c.mpf(0);scale=c.mpf(1)
            peak=logs(iv,iv.mpf(0),iv.mpf(0),iv.mpf(1))[0]
            anchor=_bounds(c,peak)[1]
            mass=-iv.expm1(-dt);masslo=_bounds(c,mass)[0]
            lower=iv.exp(iv.mpf(_bounds(c,peak)[0])-iv.mpf(anchor))*mass
            right=min(dp,max(c.mpf(80),-2*c.ln(self.target*_bounds(c,lower)[0]/100)))
            left=c.mpf(0);mode='exact exponential kernel; local beta Taylor3/fourth remainder'
            tail=iv.mpf(0) if right==dp else 2*iv.exp(-iv.mpf(right)/2)
            endpoint_error=iv.mpf(0)
            if right==dp:
                dlo,dhi=_bounds(c,dt);radius=max(abs(dlo-dp),abs(dhi-dp))
                endpoint_error=iv.mpf(radius)*iv.exp(-iv.mpf(min(dlo,dp))-iv.mpf(anchor))
            budget=self.target*_bounds(c,lower)[0]/4
        else:
            # The saddle is fixed in r for each power; its q position changes
            # with the query. It may lie just left of the finite domain.
            origin=max(c.mpf(0),U-c.exp(40*(1-rs)))
            Ustar=iv.exp(40*(1-root));qstar=Ui-Ustar
            logstar=self._log_beta(iv,J.constant(iv,root,4))[0]*power-qstar
            anchor=_bounds(c,logstar)[1]
            derivative=logs(iv,iv.mpf(0),iv.mpf(origin),iv.mpf(1));curvature=2*derivative[2]
            if _bounds(c,curvature)[1]>=0:raise ArithmeticError('Strict peak curvature required')
            scale=1/c.sqrt(-sum(_bounds(c,curvature),c.mpf(0))/2)
            slope0=derivative[1]-1
            if origin==0 and _bounds(c,slope0)[1]<0:
                # The unconstrained saddle lies left of q=0. The maximum on
                # the actual domain is its left endpoint; use its decay scale
                # instead of a much larger Gaussian saddle width.
                anchor=_bounds(c,derivative[0])[1]
                scale=min(scale,1/(-sum(_bounds(c,slope0),c.mpf(0))/2))
            window=max(c.mpf(32),-c.ln(self.target/1000))
            left=max(-window,-origin/scale);right=min(window,(dp-origin)/scale)
            def log_integrand(z):return logs(iv,iv.mpf(z),iv.mpf(origin),iv.mpf(scale))-J.variable(iv,iv.mpf(z),4)*scale-origin-anchor
            tail=iv.mpf(0)
            if origin+scale*left>0:
                L=log_integrand(left)
                if _bounds(c,L[1])[0]<=0:raise ArithmeticError('Positive left log-tail slope required')
                tail+=iv.exp(L[0])/L[1]
            if origin+scale*right<dp:
                L=log_integrand(right);slope=min(-_bounds(c,L[1])[1],scale/2)
                if slope<=0:raise ArithmeticError('Negative right log-tail slope required')
                tail+=iv.exp(L[0])/iv.mpf(slope)
            L0=log_integrand(c.mpf(0));span=c.mpf('.1')
            derivativebox=logs(iv,iv.mpf([c.mpf(0),span]),iv.mpf(origin),iv.mpf(scale))
            derivmax=max(abs(v) for v in _bounds(c,derivativebox[1]-scale))
            lower=iv.mpf(span)*iv.exp(iv.mpf(_bounds(c,L0[0])[0])-iv.mpf(span)*iv.mpf(derivmax))
            budget=self.target*_bounds(c,lower)[0]/4
            endpoint_error=iv.mpf(0)
            # If a scaled boundary is clipped at a finite endpoint, the
            # binary point endpoints need not equal the exact source limits.
            # Global exp(L-anchor)<=1 bounds this missing/excess width too.
            if left==-origin/scale:
                delta=iv.mpf(origin)+iv.mpf(scale)*iv.mpf(left)
                endpoint_error+=iv.mpf(max(abs(v) for v in _bounds(c,delta)))/iv.mpf(scale)
            if right==(dp-origin)/scale:
                delta=iv.mpf(origin)+iv.mpf(scale)*iv.mpf(right)-dt
                endpoint_error+=iv.mpf(max(abs(v) for v in _bounds(c,delta)))/iv.mpf(scale)
            mode='peak-normalized Taylor2/fourth remainder with tangent tails'
        if budget<=0:raise ArithmeticError('Positive normalized convolution budget required')
        def cell(a,b):
            ai,bi=iv.mpf(a),iv.mpf(b);mid=(ai+bi)/2;h=(bi-ai)/2;box=iv.mpf([a,b])
            mp=(a+b)/2;hp=(b-a)/2
            center=logs(iv,mid,iv.mpf(origin),iv.mpf(scale))-anchor
            whole=logs(iv,box,iv.mpf(origin),iv.mpf(scale))-anchor
            pcenter=logs(c,mp,origin,scale)-anchor
            if ordinary:
                gj=center.exp();pjet=pcenter.exp();bound=whole.exp()[4]
                cap=max(abs(v) for v in _bounds(c,bound));mass=iv.exp(-ai)*(-iv.expm1(-(bi-ai)))
                moments=[iv.exp(h)-iv.exp(-h)]
                pmom=[c.exp(hp)-c.exp(-hp)]
                for n in range(1,4):
                    moments.append((-h)**n*iv.exp(h)-h**n*iv.exp(-h)+n*moments[-1])
                    pmom.append((-hp)**n*c.exp(hp)-hp**n*c.exp(-hp)+n*pmom[-1])
                truth=sum((gj[n]*moments[n] for n in range(4)),iv.mpf(0))*iv.exp(-mid)
                point=sum((pjet[n]*pmom[n] for n in range(4)),c.mpf(0))*c.exp(-mp)
                error=iv.mpf(cap)*h**4*mass
            else:
                z=J.variable(iv,mid,4);zz=J.variable(iv,box,4);zp=J.variable(c,mp,4)
                center-=z*scale+origin;whole-=zz*scale+origin;pcenter-=zp*scale+origin
                gj=center.exp();pjet=pcenter.exp()
                # exp(L-anchor)<=1 globally by the source maximum proof.
                normalized=J(iv,[iv.mpf(0)]+list(whole.coefficients[1:])).exp()
                cap=max(abs(v) for v in _bounds(c,normalized[4]))
                # A local mean-value log bound keeps the Gaussian tails from
                # inheriting the peak value1 on every quadrature cell.
                slope=max(abs(v) for v in _bounds(c,whole[1]))
                local_upper=_bounds(c,iv.mpf(_bounds(c,center[0])[1])+h*iv.mpf(slope))[1]
                cap=_bounds(c,iv.mpf(cap)*iv.exp(iv.mpf(min(c.mpf(0),local_upper))))[1]
                truth=2*h*gj[0]+2*h**3*gj[2]/3
                point=2*hp*pjet[0]+2*hp**3*pjet[2]/3
                error=2*h**5*iv.mpf(cap)/5
            return point,truth,error
        stack=[(left,right,0)];integral=iv.mpf(0);point=c.mpf(0);accepted=0;maxdepth=0
        while stack:
            a,b,depth=stack.pop();p,t,error=cell(a,b);err=_bounds(c,error)[1]
            if err>budget*(b-a)/(right-left):
                if depth>=24:raise ArithmeticError('Convolution adaptive remainder budget not reached')
                mid=(a+b)/2;stack.extend(((mid,b,depth+1),(a,mid,depth+1)));continue
            lo,hi=_bounds(c,t+iv.mpf([-err,err]))
            integral+=iv.mpf([max(c.mpf(0),lo),max(c.mpf(0),hi)])
            point+=p;accepted+=1;maxdepth=max(maxdepth,depth)
        th=_bounds(c,tail)[1];ee=_bounds(c,endpoint_error)[1]
        integral+=iv.mpf([0,th])+iv.mpf([-ee,ee])
        result=iv.exp(-dt)+iv.exp(iv.mpf(anchor))*iv.mpf(scale)*integral
        proxy=c.exp(-dp)+c.exp(anchor)*scale*point
        lo,hi=_bounds(c,result)
        if lo<=0:raise ArithmeticError('Positive original convolution memory must be retained')
        out=dict(point=proxy,directed=result,proof=dict(power=power,method=mode,accepted_cells=accepted,max_depth=maxdepth,
            source_global_tail_proof=self.global_proof,saddle_r_interval=root,normalized_log_anchor=anchor,
            q_origin=origin,q_scale=scale,integration_left=left,integration_right=right,
            normalized_tail_upper=tail,finite_endpoint_rounding_bound=endpoint_error,
            positive_actual_memory_lower=lo,local_axial_zero_does_not_zero_memory=True,
            physical_measure='dq=d(d-v); original exponential Jacobian already absorbed',
            actual_absolute_memory_exponential_evaluated_with_unbounded_exponent_mpmath=True))
        self.cache[key]=out;return out
