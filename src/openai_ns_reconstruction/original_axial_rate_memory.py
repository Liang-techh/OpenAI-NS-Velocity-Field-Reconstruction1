"""Positive O2 axial beta kernels with original exponential rates one and two.

This extends the accepted rate-one source integration without modifying its
pinned implementation. Cell masses retain a positive source distribution for
integration of the squared own mean; homogeneous and forced parts stay separate.
"""
from .original_core_numeric import _bounds
from .original_axial_memory import OriginalAxialMemory


class OriginalAxialRateMemory(OriginalAxialMemory):
    def __init__(self,ctx,iv,J,*,kernel_rate=1,relative_target='1e-8'):
        if kernel_rate not in (1,2):raise ValueError('Original pressure kernel rate one or two required')
        super().__init__(ctx,iv,J,relative_target=relative_target)
        self.rate,self.rate_iv=ctx.mpf(kernel_rate),iv.mpf(kernel_rate)
        self._squares={}

    def _root(self,power):
        if power in self.roots:return self.roots[power]
        c,iv=self.ctx,self.iv;lo=c.mpf(0);hi=c.mpf('1e-4')
        for _ in range(220):
            mid=(lo+hi)/2;v=self._t(iv.mpf(mid))*power
            a,b=_bounds(c,v)
            if a>self.rate:lo=mid
            elif b<self.rate:hi=mid
            else:break
        out=iv.mpf([lo,hi]);self.roots[power]=out;return out

    def evaluate(self,s,power):
        if power not in (1,2):raise ValueError('Original beta power1/2 required')
        key=(s,power)
        if key in self.cache:return self.cache[key]
        c,iv,J=self.ctx,self.iv,self.J;sp=c.mpf(s);si=iv.mpf(s);ratep,ratei=self.rate,self.rate_iv
        if not c.isfinite(sp) or not 0<=sp<=1:raise ValueError('Original axial native0<=s<=1 required')
        dp=c.expm1(40*sp);dt=iv.expm1(40*si);U=c.exp(40*sp);Ui=iv.exp(40*si)
        if sp==0:
            out=dict(point=c.mpf(1),directed=iv.mpf(1),proof=dict(zero_width_actual_history=True,power=power))
            self.cache[key]=out;return out
        r0=1-sp;r0i=1-si;root=self._root(power);rl,rh=_bounds(c,root);rs=(rl+rh)/2
        ordinary=_bounds(c,r0i)[0]>=c.mpf('1e-4') or r0>0 and _bounds(c,self._t(r0i)*power)[1]<=ratep/2
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
            mass=-iv.expm1(-ratei*dt)/ratei;masslo=_bounds(c,mass)[0]
            lower=iv.exp(iv.mpf(_bounds(c,peak)[0])-iv.mpf(anchor))*mass
            right=min(dp,max(80/ratep,-2/ratep*c.ln(self.target*_bounds(c,lower)[0]/100)))
            left=c.mpf(0);mode='exact exponential kernel; local beta Taylor3/fourth remainder'
            left_tail=iv.mpf(0);right_tail=iv.mpf(0) if right==dp else 2/ratei*iv.exp(-ratei*iv.mpf(right)/2)
            tail=left_tail+right_tail
            endpoint_error=iv.mpf(0)
            if right==dp:
                dlo,dhi=_bounds(c,dt);radius=max(abs(dlo-dp),abs(dhi-dp))
                endpoint_error=iv.mpf(radius)*iv.exp(-ratei*iv.mpf(min(dlo,dp))-iv.mpf(anchor))
            budget=self.target*_bounds(c,lower)[0]/4
        else:
            # The saddle is fixed in r for each power; its q position changes
            # with the query. It may lie just left of the finite domain.
            origin=max(c.mpf(0),U-c.exp(40*(1-rs)))
            Ustar=iv.exp(40*(1-root));qstar=Ui-Ustar
            logstar=self._log_beta(iv,J.constant(iv,root,4))[0]*power-ratei*qstar
            anchor=_bounds(c,logstar)[1]
            derivative=logs(iv,iv.mpf(0),iv.mpf(origin),iv.mpf(1));curvature=2*derivative[2]
            if _bounds(c,curvature)[1]>=0:raise ArithmeticError('Strict peak curvature required')
            scale=1/c.sqrt(-sum(_bounds(c,curvature),c.mpf(0))/2)
            slope0=derivative[1]-ratei
            if origin==0 and _bounds(c,slope0)[1]<0:
                # The unconstrained saddle lies left of q=0. The maximum on
                # the actual domain is its left endpoint; use its decay scale
                # instead of a much larger Gaussian saddle width.
                anchor=_bounds(c,derivative[0])[1]
                scale=min(scale,1/(-sum(_bounds(c,slope0),c.mpf(0))/2))
            window=max(c.mpf(32),-c.ln(self.target/1000))
            left=max(-window,-origin/scale);right=min(window,(dp-origin)/scale)
            def log_integrand(z):return logs(iv,iv.mpf(z),iv.mpf(origin),iv.mpf(scale))-J.variable(iv,iv.mpf(z),4)*(scale*ratep)-origin*ratep-anchor
            left_tail=iv.mpf(0);right_tail=iv.mpf(0)
            if origin+scale*left>0:
                L=log_integrand(left)
                if _bounds(c,L[1])[0]<=0:raise ArithmeticError('Positive left log-tail slope required')
                left_tail=iv.exp(L[0])/L[1]
            if origin+scale*right<dp:
                L=log_integrand(right);slope=min(-_bounds(c,L[1])[1],scale*ratep/2)
                if slope<=0:raise ArithmeticError('Negative right log-tail slope required')
                right_tail=iv.exp(L[0])/iv.mpf(slope)
            tail=left_tail+right_tail
            L0=log_integrand(c.mpf(0));span=c.mpf('.1')
            derivativebox=logs(iv,iv.mpf([c.mpf(0),span]),iv.mpf(origin),iv.mpf(scale))
            derivmax=max(abs(v) for v in _bounds(c,derivativebox[1]-scale*ratep))
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
                cap=max(abs(v) for v in _bounds(c,bound));mass=iv.exp(-ratei*ai)*(-iv.expm1(-ratei*(bi-ai)))/ratei
                moments=[iv.exp(ratei*h)*(-iv.expm1(-2*ratei*h))/ratei]
                pmom=[c.exp(ratep*hp)*(-c.expm1(-2*ratep*hp))/ratep]
                for n in range(1,4):
                    moments.append(((-h)**n*iv.exp(ratei*h)-h**n*iv.exp(-ratei*h)+n*moments[-1])/ratei)
                    pmom.append(((-hp)**n*c.exp(ratep*hp)-hp**n*c.exp(-ratep*hp)+n*pmom[-1])/ratep)
                truth=sum((gj[n]*moments[n] for n in range(4)),iv.mpf(0))*iv.exp(-ratei*mid)
                point=sum((pjet[n]*pmom[n] for n in range(4)),c.mpf(0))*c.exp(-ratep*mp)
                error=iv.mpf(cap)*h**4*mass
            else:
                z=J.variable(iv,mid,4);zz=J.variable(iv,box,4);zp=J.variable(c,mp,4)
                center-=(z*scale+origin)*ratep;whole-=(zz*scale+origin)*ratep;pcenter-=(zp*scale+origin)*ratep
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
            metadata=dict(left=a,right=b)
            if ordinary:metadata.update(coefficients=list(gj.coefficients[:4]),point_coefficients=list(pjet.coefficients[:4]),fourth_absolute_cap=iv.mpf(cap))
            return point,truth,error,metadata
        stack=[(left,right,0)];integral=iv.mpf(0);point=c.mpf(0);accepted=0;maxdepth=0;cell_records=[]
        while stack:
            a,b,depth=stack.pop();p,t,error,record=cell(a,b);err=_bounds(c,error)[1]
            if err>budget*(b-a)/(right-left):
                if depth>=24:raise ArithmeticError('Convolution adaptive remainder budget not reached')
                mid=(a+b)/2;stack.extend(((mid,b,depth+1),(a,mid,depth+1)));continue
            lo,hi=_bounds(c,t+iv.mpf([-err,err]))
            piece=iv.mpf([max(c.mpf(0),lo),max(c.mpf(0),hi)])
            integral+=piece;record.update(point=p,truth=piece);cell_records.append(record)
            point+=p;accepted+=1;maxdepth=max(maxdepth,depth)
        th=_bounds(c,tail)[1];ee=_bounds(c,endpoint_error)[1]
        integral+=iv.mpf([0,th])+iv.mpf([-ee,ee])
        result=iv.exp(-ratei*dt)+iv.exp(iv.mpf(anchor))*iv.mpf(scale)*integral
        proxy=c.exp(-ratep*dp)+c.exp(anchor)*scale*point
        lo,hi=_bounds(c,result)
        if lo<=0:raise ArithmeticError('Positive original convolution memory must be retained')
        out=dict(point=proxy,directed=result,proof=dict(power=power,kernel_rate=int(ratep),method=mode,accepted_cells=accepted,max_depth=maxdepth,
            source_global_tail_proof=self.global_proof,saddle_r_interval=root,normalized_log_anchor=anchor,
            q_origin=origin,q_scale=scale,integration_left=left,integration_right=right,
            normalized_tail_upper=tail,finite_endpoint_rounding_bound=endpoint_error,
            positive_actual_memory_lower=lo,local_axial_zero_does_not_zero_memory=True,
            physical_measure='dq=d(d-v); original exponential Jacobian already absorbed',
            actual_absolute_memory_exponential_evaluated_with_unbounded_exponent_mpmath=True))
        out['distribution']=dict(ordinary=ordinary,anchor=anchor,scale=scale,q_origin=origin,
            left=left,right=right,actual_distance=dt,point_distance=dp,
            source_cells=sorted(cell_records,key=lambda row:row['left']),
            normalized_source_integral=integral,normalized_source_point=point,
            normalized_left_tail=left_tail,normalized_right_tail=right_tail,
            normalized_endpoint_mass_error=endpoint_error,
            normalized_homogeneous=iv.exp(-ratei*dt-iv.mpf(anchor))/iv.mpf(scale),
            normalized_point_homogeneous=c.exp(-ratep*dp-anchor)/scale,
            lower_endpoint_clipped=left==-origin/scale,
            upper_endpoint_clipped=right==(dp-origin)/scale)
        self.cache[key]=out;return out

    def squared_mean_integral(self,s):
        """Return a factored enclosure of integral exp(2d)*m(d)^2 dd.

        With q=D-d, exp(-q)m(D-q) is the future CDF of the positive
        exp(-q)beta(D-q) source, plus exp(-D). This avoids evaluating a
        second nested convolution at every pressure quadrature node.
        """
        if self.rate!=1:raise ValueError('The original own mean has rate one')
        if s in self._squares:return self._squares[s]
        c,iv=self.ctx,self.iv
        if c.mpf(s)==0:
            return dict(point=c.mpf(0),truth=iv.mpf(0),log_scale_point=c.mpf(0),log_scale=iv.mpf(0),
                proof=dict(exact_empty_integral=True))
        source=self.evaluate(s,1)['distribution'];cells=source['source_cells']
        origin,scale,anchor=(source[n] for n in ('q_origin','scale','anchor'))
        D,Dp=source['actual_distance'],source['point_distance']
        qleft=origin+scale*source['left'];qright=origin+scale*source['right']
        def positive(v):
            lo,hi=_bounds(c,v);return iv.mpf([max(c.mpf(0),lo),max(c.mpf(0),hi)])
        def upper(v):return iv.mpf(max(c.mpf(0),_bounds(c,v)[1]))
        ee=upper(source['normalized_endpoint_mass_error'])
        perturb=iv.mpf([-_bounds(c,ee)[1],_bounds(c,ee)[1]])
        total=positive(source['normalized_homogeneous']+source['normalized_source_integral'])
        height=upper(total)
        future=source['normalized_homogeneous']+source['normalized_right_tail']
        futurep=source['normalized_point_homogeneous']
        result=iv.mpf(0);proxy=c.mpf(0);cell_error=iv.mpf(0)

        def moments(ctx,h,rate,degree):
            out=[ctx.exp(rate*h)*(-ctx.expm1(-2*rate*h))/rate]
            for k in range(1,degree+1):
                out.append(((-h)**k*ctx.exp(rate*h)-h**k*ctx.exp(-rate*h)+k*out[-1])/rate)
            return out

        def square_primitive(ctx,row,F):
            a,b=ctx.mpf(row['left']),ctx.mpf(row['right']);mid=(a+b)/2;h=(b-a)/2
            co=row['coefficients' if ctx is iv else 'point_coefficients']
            # A=B+B'+B''+B'''; (exp(-q)A)'=-exp(-q)B.
            import math
            ac=[sum((co[k]*math.factorial(k)/math.factorial(j) for k in range(j,4)),ctx.mpf(0)) for j in range(4)]
            C=F-ctx.exp(-b)*sum((ac[j]*h**j for j in range(4)),ctx.mpf(0))
            mu1,mu2=moments(ctx,h,ctx.mpf(1),3),moments(ctx,h,ctx.mpf(2),6)
            first=ctx.exp(-mid)*sum((ac[j]*mu1[j] for j in range(4)),ctx.mpf(0))
            second=ctx.exp(-2*mid)*sum((ac[j]*ac[k]*mu2[j+k] for j in range(4) for k in range(4)),ctx.mpf(0))
            return C*C*(b-a)+2*C*first+second

        for row in reversed(cells):
            F=positive(future+perturb);width=iv.mpf(scale)*(iv.mpf(row['right'])-iv.mpf(row['left']))
            if source['ordinary']:
                val=square_primitive(iv,row,F);point=square_primitive(c,row,futurep)
                a,b=iv.mpf(row['left']),iv.mpf(row['right']);h=(b-a)/2
                eta=row['fourth_absolute_cap']*h**4
                mass1=iv.exp(-a)*(-iv.expm1(-(b-a)))
                mass2=iv.exp(-2*a)*(-iv.expm1(-2*(b-a)))/2
                error=upper(2*height*eta*mass1+3*eta*eta*mass2)
                err=_bounds(c,error)[1];result+=positive(val+iv.mpf([-err,err]));cell_error+=error
                proxy+=max(c.mpf(0),point)
            else:
                # The CDF is monotone inside every accepted positive cell.
                # This is a source-adaptive enclosure, not a uniform s grid.
                low=_bounds(c,F)[0];high=_bounds(c,positive(F+row['truth']))[1]
                result+=width*iv.mpf([low*low,high*high])
                proxy+=scale*(row['right']-row['left'])*(futurep+row['point']/2)**2
            future+=row['truth'];futurep+=row['point']
        atleft=positive(future+perturb)
        leftwidth=max(c.mpf(0),qleft)
        result+=iv.mpf(leftwidth)*iv.mpf([_bounds(c,atleft)[0]**2,_bounds(c,height)[1]**2])
        proxy+=leftwidth*futurep**2
        rightwidth=positive(D-iv.mpf(qright))
        rightheight=upper(source['normalized_homogeneous']+source['normalized_right_tail']+ee)
        result+=rightwidth*iv.mpf([0,_bounds(c,rightheight)[1]**2])
        # Clipped binary endpoints may differ from exact physical endpoints.
        length_error=iv.mpf(0)
        if source['lower_endpoint_clipped']:
            length_error+=upper(abs(iv.mpf(origin)+iv.mpf(scale)*iv.mpf(source['left'])))
        if source['upper_endpoint_clipped']:
            length_error+=upper(abs(iv.mpf(origin)+iv.mpf(scale)*iv.mpf(source['right'])-D))
        rounding=upper(length_error*height*height);err=_bounds(c,rounding)[1]
        result=positive(result+iv.mpf([-err,err]))
        out=dict(point=proxy,truth=result,log_scale_point=2*Dp+2*anchor+2*c.log(scale),
            log_scale=2*D+2*iv.mpf(anchor)+2*iv.log(iv.mpf(scale)),
            proof=dict(positive_future_source_CDF=True,source_cells=len(cells),
                ordinary_exponential_polynomial_primitive=source['ordinary'],
                peak_monotone_cell_enclosure=not source['ordinary'],
                normalized_own_mean_square_integral=result,normalized_Taylor_square_error=cell_error,
                endpoint_length_rounding_error=rounding,both_positive_source_tails_retained=True,
                nested_own_mean_quadrature_avoided=True))
        self._squares[s]=out;return out
