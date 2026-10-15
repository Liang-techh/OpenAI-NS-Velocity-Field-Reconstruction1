"""Numerical leading+n1 core velocity as source-bounded exponential sums.

The absolute F0 exponential is never materialized. Original tiny meridional
corrections, delta*time terms and temporal orders remain separate. The
Cartesian entry point embeds a directed physical inverse in every numeric
backend's input boxes before reusing its original-source error bounds.
"""
import copy
import hashlib
import json
import math

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from .original_core_log_amplitude import OriginalCoreLogAmplitude
from .original_core_numeric import ROOT,PREFIX,_box,_bounds
from .original_n1_core_volterra_numeric import OriginalN1CoreVolterraNumeric


class OriginalCoreVelocity:
    """Physical core u,v,w functions in native or Cartesian coordinates.

    Output is a signed sum of coefficient*exp(separate logarithmic terms),
    with directed errors on every coefficient/log. It is not a float vector
    or an admitted infinite temporal sum. Core-only spatial support applies.
    """
    def __init__(self,*,precision=70,amplitude_precision=110,directory=None):
        self.n1 = OriginalN1CoreVolterraNumeric(precision=precision,quadrature_nodes=192,directory=directory)
        self.amplitude = OriginalCoreLogAmplitude(precision=amplitude_precision,directory=directory)
        if self.n1.family != self.amplitude.family:raise ValueError("One original source family required")
        self.family = self.n1.family
        self.ctx,self.iv = self.amplitude.ctx,self.amplitude.iv
        self.hashes = dict(self.n1.hashes)
        for name,digest in self.amplitude.hashes.items():
            if self.hashes.get(name,digest) != digest:raise ValueError("Inconsistent original source pin")
            self.hashes[name] = digest
        self.code_hashes = dict(self.n1.code_hashes)
        for name in ('original_core_log_amplitude.py','original_n1_core_volterra_numeric.py'):
            path = ROOT/'src/openai_ns_reconstruction'/name
            self.code_hashes[path.relative_to(ROOT).as_posix()] = hashlib.sha256(path.read_bytes()).hexdigest()
        report_path = self.n1.core.directory/(PREFIX+'compliant_current_original_n1_core_volterra_numeric.json')
        report = json.loads(report_path.read_text(encoding='utf-8'))
        if not report['checks']['passed'] or report['source_family'] != self.family:
            raise ValueError("Checked full-core n1 numerical implementation required")
        for relative,digest in report['implementation_hashes'].items():
            if hashlib.sha256((ROOT/relative).read_bytes()).hexdigest() != digest:
                raise ValueError("Changed n1 implementation: "+relative)
            self.code_hashes[relative] = digest
        self.hashes[report_path.name] = hashlib.sha256(report_path.read_bytes()).hexdigest()
        self._leading_cache,self._n1_cache = {},{}

    def _restore(self,row,key='source_enclosure_exact_mpf_tuples'):
        return self.iv.mpf([mp.make_mpf(tuple(v)) for v in row[key]])

    def _leading(self,Z,rho):
        leading = self.n1.leading
        if Z not in self._leading_cache:self._leading_cache[Z] = leading.coefficients(Z,degree=18)
        packet = self._leading_cache[Z]
        sample = leading.evaluate(packet,rho)
        c,iv = self.ctx,self.iv
        def profile(key):
            row = sample['profiles'][key]
            return c.mpf(row['approximation']),self._restore(row)
        out = {key:profile(key) for key in ('Phi','Psi','mean_Psi')}
        z,zi = c.mpf(Z),iv.mpf(Z)
        r,ri = c.mpf(rho),iv.mpf(rho)
        de,dei = self.amplitude.delta,self.amplitude.delta_iv
        U0,U0i = 4*z+self.amplitude.j,4*zi+self.amplitude.j_iv
        out['U0'] = U0,U0i
        # This derivative is the same original mean, not finite differences
        # of a rounded Uz value in which epsilon*Psi disappeared.
        mz = sum((c.mpf(packet['point']['Uz'][n][1])*r**n/(c.mpf(leading.epsilon)*(n+1))
            for n in range(1,19)),c.mpf(0))
        mzi = sum((iv.mpf(packet['directed']['Uz'][n][1])*ri**n/(iv.mpf(leading.epsilon_iv)*(n+1))
            for n in range(1,19)),iv.mpf(0))
        if r:
            correction = _box(iv,leading.major['scaled_map_size_upper'])/(1-_box(iv,leading.major['scaled_map_Lipschitz_upper']))
            tail = correction*leading.tail_factor(iv,degree=18,radial_order=0,axial_order=1,
                radius=ri,h=_box(iv,leading.tube['Xh_parameter']))['tail_per_Xh_norm']
            high = _bounds(c,tail)[1]
            mzi += iv.mpf([-high,high])
        P,Pi = out['Psi'];M,Mi = out['mean_Psi']
        out['Qaxis'] = ((1+de)*z*U0-4*(1-z*z))/(1-de*z*z),((1+dei)*zi*U0i-4*(1-zi*zi))/(1-dei*zi*zi)
        out['Qpsi'] = (2*z*P-(1-de)*z*M-(1-z*z)*mz)/(1-de*z*z),(
            2*zi*Pi-(1-dei)*zi*Mi-(1-zi*zi)*mzi)/(1-dei*zi*zi)
        return out

    def _n1(self,Z,rho):
        key = Z,rho
        if key not in self._n1_cache:
            self._n1_cache[key] = self.n1.evaluate(Z,rho,degree=18,kernel_degree=24)
        packet = self._n1_cache[key]
        out = {}
        for key in ('f','u','mean_u','q'):
            row = packet['normalized_original_n1'][key]
            row = row[0] if isinstance(row,list) else row
            out[key] = self.ctx.mpf(row['approximation']),self._restore(row,'original_source_enclosure_exact_mpf_tuples')
        return out

    def _summary(self,pair):return self.amplitude.summarize(*pair)

    def _time(self,Z,log_tau,scaled_time):
        c,iv = self.ctx,self.iv
        z,zi = c.mpf(Z),iv.mpf(Z)
        half_log_d = -c.log(1-z*z)/2,-iv.ln(1-zi*zi)/2
        if scaled_time is None:
            if type(log_tau) is not str or not c.isfinite(c.mpf(log_tau)) or c.mpf(log_tau)>0:
                raise ValueError("Finite decimal log(tau)<=0 required")
            a = c.mpf(log_tau)/2,iv.mpf(log_tau)/2
            da = self.amplitude.delta*a[0],self.amplitude.delta_iv*a[1]
            definition = dict(kind='log_tau',exact_log_tau=log_tau)
        else:
            if (log_tau is not None or type(scaled_time) is not str
                    or not c.isfinite(c.mpf(scaled_time)) or c.mpf(scaled_time)<0):
                raise ValueError("Use log_tau=None and a finite decimal S=-delta*log(tau)/2>=0")
            S,Si = c.mpf(scaled_time),iv.mpf(scaled_time)
            a = -S/self.amplitude.delta,-Si/self.amplitude.delta_iv
            # Cancel delta exactly before boxing the time product.
            da = -S,-Si
            definition = dict(kind='source_scaled_time',exact_S=scaled_time,
                log_tau_definition='-2*S/original_delta',delta_cancelled_before_time_comparison=True)
        db = self.amplitude.delta*half_log_d[0],self.amplitude.delta_iv*half_log_d[1]
        return [a,half_log_d],[da,db],definition

    def _assemble(self,Z,rho,log_tau,scaled_time,rotation):
        c,iv = self.ctx,self.iv
        z,r = c.mpf(Z),c.mpf(rho)
        if not (-1<z<1 and 0<=r<=4):raise ValueError("Original core requires |Z|<1 and 0<=rho<=4")
        profiles,n1 = self._leading(Z,rho),self._n1(Z,rho)
        ell,de,time = self._time(Z,log_tau,scaled_time)
        ll = self.amplitude.ll,self.amplitude.ll_iv
        amplitude = self.amplitude.global_amplitude(Z)
        F0 = [dict(negative_of_exact_mpf_tuple=self.amplitude.selected_logC_tuple),
              dict(amplitude['log_F0_relative_to_anchor'])]
        def scale(pair,value):return pair[0]*value,pair[1]*value
        def logs(a,b,L=0,geometry=None,swirl=False):
            out = [self._summary(scale(q,a)) for q in ell if a]
            out += [self._summary(scale(q,b)) for q in de if b]
            if L:out.append(self._summary(scale(ll,L)))
            if geometry is not None:out.append(self._summary(geometry))
            return copy.deepcopy(F0)+out if swirl else out
        def term(name,order,pair,exponents):
            return dict(name=name,temporal_order=order,coefficient=self._summary(pair),log_scale_terms=exponents)
        def zero(name,order):return term(name,order,(c.mpf(0),iv.mpf(0)),[])
        axial = [term('leading_axis',0,profiles['U0'],logs(-1,-1)),
            term('leading_axial_correction',0,profiles['Psi'],logs(-1,-1,-1)),
            term('n1_axial',1,n1['u'],logs(-1,1,-1))]
        if r:
            ri = iv.mpf(rho)
            radial_geometry = (c.log(r/2)/2,iv.ln(ri/2)/2)
            swirl_geometry = (c.log(2*r)/2,iv.ln(2*ri)/2)
            radial = [term('leading_radial_axis',0,profiles['Qaxis'],logs(-1,0,c.mpf('-0.5'),radial_geometry)),
                term('leading_radial_correction',0,profiles['Qpsi'],logs(-1,0,c.mpf('-1.5'),radial_geometry)),
                term('n1_radial',1,n1['q'],logs(-1,2,c.mpf('-1.5'),radial_geometry))]
            swirl = [term('leading_swirl',0,profiles['Phi'],logs(-1,-1,c.mpf('-0.5'),swirl_geometry,True)),
                term('n1_swirl',1,n1['f'],logs(-1,1,c.mpf('0.5'),swirl_geometry,True))]
            mass_geometry = c.log(r),iv.ln(ri)
            psi = [term('leading_mass_axis',0,profiles['U0'],logs(1,-1,-1,mass_geometry)),
                term('leading_mass_correction',0,profiles['mean_Psi'],logs(1,-1,-2,mass_geometry)),
                term('n1_mass',1,n1['mean_u'],logs(1,1,-2,mass_geometry))]
            radius = term('physical_radius',None,(c.mpf(1),iv.mpf(1)),logs(1,0,c.mpf('-0.5'),swirl_geometry))
        else:
            radial = [zero('leading_radial_axis',0),zero('leading_radial_correction',0),zero('n1_radial',1)]
            swirl = [zero('leading_swirl',0),zero('n1_swirl',1)]
            psi = [zero('leading_mass_axis',0),zero('leading_mass_correction',0),zero('n1_mass',1)]
            radius = zero('physical_radius',None)
        cos,sin = rotation
        def rotate(packet,pair):
            out = copy.deepcopy(packet)
            coeff = packet['coefficient']
            out['coefficient'] = self._summary((c.mpf(coeff['approximation'])*pair[0],self._restore(coeff)*pair[1]))
            return out
        negsin = -sin[0],-sin[1]
        Cartesian = dict(u=[rotate(q,cos) for q in radial]+[rotate(q,negsin) for q in swirl],
            v=[rotate(q,sin) for q in radial]+[rotate(q,cos) for q in swirl],w=copy.deepcopy(axial))
        comparison = dict(available=False,reason='both physical swirl terms vanish on the axis')
        if r:
            f,fi = n1['f'];P,Pi = profiles['Phi']
            flo,fhi = _bounds(c,fi)
            if f and (flo>0 or fhi<0) and _bounds(c,Pi)[0]>0:
                logpoint = self.amplitude.ll+2*sum((q[0] for q in de),c.mpf(0))+c.log(abs(f)/P)
                logbox = self.amplitude.ll_iv+2*sum((q[1] for q in de),iv.mpf(0))+iv.ln(abs(fi)/Pi)
                lo,hi = _bounds(c,logbox)
                comparison = dict(available=True,log_absolute_ratio=self.amplitude.summarize(logpoint,logbox),
                    classification='n1_smaller' if hi<0 else 'n1_larger' if lo>0 else 'unresolved',
                    formula='logLambda+2*delta*log(lambda)+log(abs(f/Phi))',
                    F0_cancelled_exactly=True,proves_higher_temporal_recursion=False)
        return dict(native=dict(Z=Z,rho=rho),time=time,source_family=self.family,
            coordinates=dict(r=radius,x=rotate(radius,cos),y=rotate(radius,sin),
                z=term('physical_z',None,(z,iv.mpf(Z)),logs(1,-1)),
                t=dict(kind='one_minus_exp_of_log_tau',definition=time)),
            cylindrical=dict(ur=radial,utheta=swirl,uz=axial),Cartesian=Cartesian,
            streamfunction=psi,swirl_order_comparison=comparison,
            value_semantics='sum coefficient*exp(sum separate log_scale_terms); all coefficient/log intervals retained',
            coefficient_orders_included=[0,1],shared_streamfunction_meridional_recovery=True,
            tiny_leading_corrections_and_delta_time_terms_kept_separate=True,
            original_core_coefficient_and_log_errors_available=True,
            useful_combined_physical_relative_error_certified=False,
            absolute_F0_or_physical_velocity_exponentials_materialized=False,
            whole_global_temporal_sum_admitted=False,global_field_available=False,
            genuine_n_ge_2_temporal_recursion=False,full_corrected_NS_solution=False)

    def evaluate(self,Z,rho,*,log_tau='-1',scaled_time=None,azimuth='0'):
        c,iv = self.ctx,self.iv
        if type(Z) is not str or type(rho) is not str or type(azimuth) is not str or not c.isfinite(c.mpf(azimuth)):
            raise ValueError("Use finite decimal native coordinates/angle")
        phi,phii = c.mpf(azimuth),iv.mpf(azimuth)
        return self._assemble(Z,rho,log_tau,scaled_time,((c.cos(phi),iv.cos(phii)),(c.sin(phi),iv.sin(phii))))

    def cartesian(self,x,y,z,*,log_tau='-1'):
        """Physical x,y,z -> core velocity, after directed inverse containment."""
        return self._cartesian(x,y,z,log_tau,None)

    def at_time(self,x,y,z,t):
        """Physical x,y,z,t -> u,v,w; use log_tau for unresolved near-T times."""
        return self._cartesian(x,y,z,None,t)

    def _cartesian(self,x,y,z,log_tau,t):
        if any(type(v) is not str for v in (x,y,z,t if t is not None else log_tau)):
            raise ValueError("Exact decimal strings required")
        c = mp.mp.clone();c.dps = 320
        iv = MPIntervalContext();iv.dps = 320
        xp,yp,zp = (c.mpf(v) for v in (x,y,z))
        xi,yi,zi = (iv.mpf(v) for v in (x,y,z))
        if t is None:
            lt,lti = c.mpf(log_tau),iv.mpf(log_tau)
        else:
            tp,ti = c.mpf(t),iv.mpf(t)
            if not c.isfinite(tp) or not 0<=tp<1:
                raise ValueError("Physical time requires 0<=t<1; use log_tau near T")
            lt,lti = c.log1p(-tp),iv.ln(1-ti)
            if not all(c.isfinite(v) for v in _bounds(c,lti)):
                raise ValueError("Physical time unresolved; use log_tau or native scaled time")
            log_tau = c.nstr(lt,300)
        if not all(c.isfinite(v) for v in (xp,yp,zp,lt)) or lt>0:raise ValueError("Finite coordinates/log(tau)<=0 required")
        lp = iv.exp(40)+11
        delta = iv.exp(-4*lp-30)
        tolerance = c.mpf('1e-280')
        if zp == 0:
            Zbox = iv.mpf(0);Ztext = '0';iterations = 0
        else:
            target = iv.ln(abs(zi))
            def H(value):
                q = iv.mpf(value)
                return iv.ln(q)+(1-delta)*(lti-iv.ln(1-q*q))/2-target
            lo,hi = c.mpf(0),c.mpf(1)
            iterations = 0
            while hi-lo>tolerance:
                mid = (lo+hi)/2
                a,b = _bounds(c,H(mid))
                if b<0:lo=mid
                elif a>0:hi=mid
                else:
                    left,right = mid-tolerance/4,mid+tolerance/4
                    if not (0<left<right<1 and _bounds(c,H(left))[1]<0 and _bounds(c,H(right))[0]>0):
                        raise ValueError("Physical inverse unresolved at current precision")
                    lo,hi=left,right
                iterations += 1
                if iterations>1200:raise ValueError("Physical inverse iteration limit")
            if zp<0:lo,hi=-hi,-lo
            Zbox = iv.mpf([lo,hi]);Ztext=c.nstr((lo+hi)/2,300)
        r2 = xi*xi+yi*yi
        if xp == yp == 0:
            rhobox=iv.mpf(0);rhotext='0'
            rotation=((self.ctx.mpf(1),self.iv.mpf(1)),(self.ctx.mpf(0),self.iv.mpf(0)))
        else:
            logrho = iv.ln(r2)-iv.ln(2)-lti+iv.ln(1-Zbox*Zbox)+4*lp+1000
            if _bounds(c,logrho)[1]>c.log(4):raise ValueError("Physical point is outside the original core")
            # Avoid constructing a Python exponent integer with astronomical
            # bit length; such inputs must use the native factored interface.
            if max(v._mpf_[2] for v in _bounds(c,abs(logrho)))>4096:
                raise ValueError("Use native factored coordinates for this unmaterializable radius")
            rhobox = iv.exp(logrho)
            l,h = _bounds(c,rhobox);rhotext=c.nstr((l+h)/2,300)
            rr = iv.sqrt(r2)
            rotation=tuple((self.ctx.mpf(xp/c.sqrt(xp*xp+yp*yp) if name=='x' else yp/c.sqrt(xp*xp+yp*yp)),
                self.iv.mpf([mp.make_mpf(v) for v in (value/rr)._mpi_])) for name,value in (('x',xi),('y',yi)))
        for context in (self.n1.iv,self.iv):
            for text,box in ((Ztext,Zbox),(rhotext,rhobox),(log_tau,lti)):
                a,b=_bounds(c,context.mpf(text));lo,hi=_bounds(c,box)
                if not a<=lo<=hi<=b:
                    raise ValueError("Directed inverse must be contained in both source input boxes")
        out=self._assemble(Ztext,rhotext,log_tau,None,rotation)
        out['physical_cartesian_request']=dict(x=x,y=y,z=z,
            **(dict(exact_log_tau=log_tau) if t is None else dict(exact_t=t)))
        if t is not None:
            out['time']=dict(kind='exact_physical_t',exact_t=t,log_tau_numeric_proxy=log_tau,
                source_time_boxes_contain_true_log_tau=True)
            out['coordinates']['t']=dict(kind='exact_decimal_physical_time',value=t)
        out['physical_inverse']=dict(directed=True,bisection_iterations=iterations,
            source_input_boxes_contain_entire_inverse=True,
            source_time_boxes_contain_true_log_tau=True,
            Z_exact_mpf_tuples=[list(v) for v in Zbox._mpi_],rho_exact_mpf_tuples=[list(v) for v in rhobox._mpi_],
            monotonicity='d(z_physical)/dZ=lambda^(1-delta)*(1-delta*Z^2)/(1-Z^2)>0',
            no_uncontrolled_interpolation_from_native_point=True)
        return out
