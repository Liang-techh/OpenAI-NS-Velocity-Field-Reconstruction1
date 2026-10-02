"""Whole-Z analytic core jets and nonsingular physical Cartesian core map.

Evaluates enclosures of the SAME selected-Cstar analytic fixed point, using
its explicit leading pair plus the admitted Xh correction norm. A factorial
tail bounds the infinite Bessel model. No local old coefficient state is
extrapolated, and no model is substituted for the nonlinear solution.
"""
import functools
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_cartesian_field import CompliantCartesianField,INDICES
from lei_ren_part1_paper_compliant_core_uniform_bounds import embedding
from lei_ren_part1_paper_compliant_pressure_source import CompliantPressureDatum
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_'
X,Y,Z,D,B=s.symbols('X Y Z delta beta',real=True)
Q,F,V,PD,PI='radial_recovery_Q','Fcore_over_F0','axial_velocity','axis_pressure_over_Pstar_squared','pressure_increment_over_epsilon_F0_squared'
BASES={'ux':{Q:X/2,F:-Y},'uy':{Q:Y/2,F:X},'uz':{V:s.Integer(1)},'p':{PD:s.Integer(1),PI:s.Integer(1)}}


def gridkey(i,k):return 'rho'+str(i)+'_Z'+str(k)


def symmetric(c,bound):
    hi=endpoints(bound)[1]
    if hi<0:raise ValueError('Positive correction bound required')
    return c.mpf([-hi,hi])


def intersection(c,a,b):
    lo=max(endpoints(a)[0],endpoints(b)[0]); hi=min(endpoints(a)[1],endpoints(b)[1])
    if lo>hi:raise ArithmeticError('Independent analytic core enclosures do not intersect')
    return c.mpf([lo,hi])


def model_phi_jets(c,chi,rho,degree=80):
    """All rho/Z mixed derivatives<=5 of B(chi*rho/2), including infinite tail."""
    if degree<6 or endpoints(rho)[0]<0 or endpoints(rho)[1]>endpoints(c.mpf('4.1'))[1]:
        raise ValueError('Core rho in[0,4.1] and factorial-tail degree>=6 required')
    M=sum((c.mpf(max(abs(v) for v in endpoints(chi[k]))) for k in range(1,6)),c.mpf(0))
    powers=[IntervalTaylor.constant(c,1,5)]
    for _ in range(degree):powers.append(powers[-1]*chi)
    rows={}; tails={}; rmax=c.mpf(endpoints(rho)[1])
    for i in range(6):
        value=IntervalTaylor.constant(c,0,5)
        for n in range(i,degree+1):
            value+=powers[n]*((-1)**n*rho**(n-i)/(2**n*math.factorial(n-i)*math.factorial(n+1)))
        for k in range(6-i):
            n=degree+1
            first=(rmax**(n-i)/(2**n*math.factorial(n-i)*math.factorial(n+1))
                   *(k+1)*(n+1)**k*(1+M)**k)
            ratio=rmax/2*(c.mpf(n+2)/(n+1))**k/((n+1-i)*(n+2))
            if endpoints(ratio)[1]>=1:raise ArithmeticError('Factorial model tail ratio failed')
            tail=first/(1-ratio)*math.factorial(k)
            rows[gridkey(i,k)]=value[k]*math.factorial(k)+symmetric(c,tail)
            tails[gridkey(i,k)]=tail
    # The alternating cubic/quadratic proof is independent of interval summation.
    q=rmax/2; floor=1-q/2+q*q/12-q**3/144
    rows[gridkey(0,0)]=intersection(c,rows[gridkey(0,0)],c.mpf([endpoints(floor)[0],1]))
    return rows,tails


def relative_amplitude_jet(c,ell,multiplier=1):
    rows=[c.mpf(1)]
    for n in range(1,6):rows.append(sum((multiplier*ell[j]*rows[n-1-j] for j in range(n)),c.mpf(0))/n)
    return [value*math.factorial(k) for k,value in enumerate(rows)]


def square_grid(c,grid,order=4):
    return {gridkey(i,k):sum((math.comb(i,a)*math.comb(k,b)*grid[gridkey(a,b)]*grid[gridkey(i-a,k-b)]
        for a in range(i+1) for b in range(k+1)),c.mpf(0)) for i in range(order+1) for k in range(order+1-i)}


def core_step(row,direction,gamma=None):
    """No inverse-r factor: normalized X,Y satisfy rho=(X^2+Y^2)/2."""
    out={}
    def add(key,value):out[key]=out.get(key,s.Integer(0))+value
    for (i,k),coefficient in row.items():
        if direction in ('x','y'):
            variable=X if direction=='x' else Y
            add((i,k),s.diff(coefficient,variable)); add((i+1,k),variable*coefficient)
        elif direction=='z':
            L=1-D*Z**2; d=1-Z**2
            add((i,k),(gamma*Z*coefficient+d*s.diff(coefficient,Z)-Z*(X*s.diff(coefficient,X)+Y*s.diff(coefficient,Y)))/L)
            add((i,k+1),d*coefficient/L); add((i+1,k),-Z*(X**2+Y**2)*coefficient/L)
        else:raise ValueError('Physical core derivative direction required')
    return {key:s.cancel(value) for key,value in out.items() if value!=0}


@functools.lru_cache(maxsize=1)
def core_templates():
    out={}
    for component,bases in BASES.items():
        for label,seed in bases.items():
            for nx,ny,nz in INDICES:
                row={(0,0):seed}
                for _ in range(nx):row=core_step(row,'x')
                for _ in range(ny):row=core_step(row,'y')
                for b in range(nz):row=core_step(row,'z',B-nx-ny+b*(D-1))
                out[(component,label,nx,ny,nz)]=row
    return out


def core_time_template(seed):
    L=1-D*Z**2
    return {(0,0):(-B*seed/2+(1-D)*Z*s.diff(seed,Z)/2+(X*s.diff(seed,X)+Y*s.diff(seed,Y))/2)/L,
            (0,1):(1-D)*Z*seed/(2*L),
            (1,0):(X**2+Y**2)*seed/(2*L)}


def coefficient_value(c,expression,x,y,z,delta,beta):
    values={X:x,Y:y,Z:z,D:delta,B:beta}
    if expression in values:return values[expression]
    if expression.is_Rational:return c.mpf(int(expression.p))/int(expression.q)
    if expression.is_Add:return sum((coefficient_value(c,v,x,y,z,delta,beta) for v in expression.args),c.mpf(0))
    if expression.is_Mul:
        value=c.mpf(1)
        for v in expression.args:value*=coefficient_value(c,v,x,y,z,delta,beta)
        return value
    if expression.is_Pow and expression.args[1].is_Integer:return coefficient_value(c,expression.args[0],x,y,z,delta,beta)**int(expression.args[1])
    raise ValueError('Unsupported exact core operator expression')


class CompliantCorePhysicalField:
    def __init__(self):
        self.outer=CompliantCartesianField(); self.ctx=c=self.outer.ctx
        self.delta=self.outer.delta; self.family=self.outer.family; self.source=self.outer.source
        self.hashes=dict(self.outer.hashes); self.records={}
        for part in ('cartesian_field_check','physical_energy_check'):
            name=PREFIX+'compliant_'+part+'.json'; check=json.loads((HERE/name).read_bytes())
            if not check['all_passed'] or check['implicit_source_sha256']!=self.source or check['actual_five_defect_family_sha256']!=self.family:
                raise ValueError('Accepted outer physical map/energy prerequisites required')
            for path,digest in check['input_hashes'].items():
                if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:raise ValueError('Accepted physical source changed: '+path)
            self.hashes.update(check['input_hashes']); self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        for part in ('core_transfer','core_transfer_check','core_uniform_bounds','physical_norm_family','shared_analytic_tube','shared_linear_resolvent'):
            name=PREFIX+(part if part.startswith('shared_') else 'compliant_'+part)+'.json'
            record=json.loads((HERE/name).read_bytes()); self.records[part]=record
            for path,digest in record.get('input_hashes',{}).items():
                if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:raise ValueError('Core physical source changed: '+path)
            self.hashes.update(record.get('input_hashes',{})); self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        major=self.records['core_transfer']; norm=self.records['physical_norm_family']; uniform=self.records['core_uniform_bounds']
        analytic_family=major['analytic_core_family_sha256']
        if (norm['base_analytic_core_family_sha256']!=analytic_family
                or uniform['analytic_core_family_sha256']!=analytic_family
                or self.records['shared_analytic_tube']['analytic_core_family_sha256']!=analytic_family
                or major['implicit_source_sha256']!=self.source):
            raise ValueError('Core, tube and selected-Cstar analytic family differ')
        if not self.records['core_transfer_check']['all_passed'] or not major['contraction_proved'] or not norm['uniform_analytic_fixed_point_admission_extended']:
            raise ValueError('Selected-Cstar analytic fixed point admission required')
        if not uniform['full_real_axis_normalized_swirl_positive'] or norm['implicit_source_sha256']!=self.source:
            raise ValueError('Positive same-source normalized core required')
        read=lambda r,key:read_interval(c,r[key])
        self.logLambda=read(major,'logLambda'); self.Lambda=read(major,'Lambda'); self.epsilon=read(major,'epsilon')
        self.logP=read(major,'logPstar'); self.logC=read(norm,'selected_logCstar')
        self.j=read(major,'required_j'); self.h=read(self.records['shared_analytic_tube'],'Xh_parameter')
        self.sigma=self.j/500
        q=read(major,'scaled_map_Lipschitz_upper')
        if endpoints(q)[1]>=1:raise ValueError('Admitted analytic contraction factor required')
        # The size majorant is already uniform over the ball. The extra
        # factor is conservative and also dominates an a-posteriori bound.
        self.correction=read(major,'scaled_map_size_upper')/(1-q)
        self.phi_floor=read(uniform,'actual_Phi_lower'); self.phi_ceiling=read(uniform,'actual_Phi_upper')
        self.Gbar=read(major,'Gupper_in_logC_definition')
        self.datum=CompliantPressureDatum('40',160)
        if self.datum.source_sha!=self.source or self.datum.datum_sha!=major['datum_enclosure_sha256']:
            raise ValueError('Core physical datum differs from accepted analytic source')
        self.hashes.update(self.datum.input_hashes); self.hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        self.cache={}

    def axis_inputs(self,Zvalue):
        c=self.ctx; z=IntervalTaylor.variable(c,Zvalue,6); one=IntervalTaylor.constant(c,1,6)
        u=4*z+self.j; d=one-z*z; L=one-(z*z)*self.delta; H=z*((1-self.delta)/2)+d*u
        h2=H*H; square=list(h2.coefficients); square[0]=H[0]**2; h2=IntervalTaylor(c,square)
        denominator=h2+self.sigma**2; chi=h2/denominator
        lo,hi=endpoints(h2[0]); sl,sh=endpoints(self.sigma**2)
        chi0=c.mpf([endpoints(c.mpf(lo)/(c.mpf(lo)+sh))[0],endpoints(c.mpf(hi)/(c.mpf(hi)+sl))[1]])
        chi=IntervalTaylor(c,[chi0]+list(chi.coefficients[1:])).truncate(5)
        datum=self.datum.normalized_jets(endpoints(c.mpf(Zvalue)),6)
        pressure=IntervalTaylor(c,[c.mpf(endpoints(v)) for v in datum['normalized_pressure_coefficients']])
        physicalP=pressure*c.exp(2*self.logP)
        pz=IntervalTaylor(c,[(k+1)*physicalP[k+1] for k in range(6)])
        g=-((1-2*z*u)*u*((1+self.delta)/2)+4*H+d*pz-z*physicalP*(2*(1+self.delta)))
        slope=-g/(2*L)
        gradient=L*H/denominator
        ell=[-self.Lambda*gradient[k] for k in range(5)]
        return dict(z=z,u=u,d=d,L=L,H=H,chi=chi,slope=slope.truncate(5),
                    beta=((z*self.j+3-self.delta/2)/L).truncate(5),pressure=pressure,
                    relative_F0_derivatives=relative_amplitude_jet(c,ell),
                    relative_F0_squared_derivatives=relative_amplitude_jet(c,ell,2))

    def normalized_jets(self,rho,Zvalue):
        c=self.ctx; rho=c.mpf(rho); z=c.mpf(Zvalue)
        if endpoints(rho)[0]<0 or endpoints(rho)[1]>endpoints(c.mpf('4.1'))[1] or endpoints(z)[0]<-1 or endpoints(z)[1]>1:
            raise ValueError('Original analytic core domain rho in[0,4.1], Z in[-1,1] required')
        key=(rho._mpi_,z._mpi_)
        if key in self.cache:return self.cache[key]
        source=self.axis_inputs(z); phi,tails=model_phi_jets(c,source['chi'],rho)
        psi={}; average={}; errors={}; rmax=c.mpf(endpoints(rho)[1]); axis=endpoints(rho)[1]==0
        for i in range(6):
            for k in range(6-i):
                index=gridkey(i,k)
                if i==0:
                    err=self.correction*rmax*embedding(c,self.h,rmax,1,k)
                    avgerr=err/2
                else:
                    err=self.correction*embedding(c,self.h,rmax,i,k); avgerr=err/(i+1)
                errors[index]=err
                phi[index]+=symmetric(c,err)
                model=(source['slope'][k]*math.factorial(k))*(rho if i==0 else 1 if i==1 else 0)
                psi[index]=model+symmetric(c,err); average[index]=model/2+symmetric(c,avgerr)
                if axis and i==0:
                    phi[index]=c.mpf(1 if k==0 else 0); psi[index]=average[index]=c.mpf(0)
                elif axis and i==1:
                    phi[index]=-(source['chi'][k]+self.epsilon*source['beta'][k])*math.factorial(k)/4
                    psi[index]=model; average[index]=model/2
        phi[gridkey(0,0)]=intersection(c,phi[gridkey(0,0)],c.mpf([endpoints(self.phi_floor)[0],endpoints(self.phi_ceiling)[1]]))
        Vgrid={}; Mgrid={}
        for i in range(6):
            for k in range(6-i):
                base=source['u'][k]*math.factorial(k) if i==0 else c.mpf(0)
                Vgrid[gridkey(i,k)]=base+self.epsilon*psi[gridkey(i,k)]
                Mgrid[gridkey(i,k)]=base+self.epsilon*average[gridkey(i,k)]
        result=dict(rho=rho,Z=z,Phi=phi,Psi=psi,Uz=Vgrid,Mz_over_R=Mgrid,
                    explicit_model_factorial_tail=tails,nonlinear_correction_bounds=errors,source=source)
        self.cache[key]=result; return result

    def profiles(self,rho,Zvalue):
        c=self.ctx; packet=self.normalized_jets(rho,Zvalue); rho=packet['rho']; z=packet['Z']; source=packet['source']
        phi=packet['Phi']; vel=packet['Uz']; mean=packet['Mz_over_R']; square=square_grid(c,phi)
        covering=self.normalized_jets(c.mpf([0,endpoints(rho)[1]]),z)
        square_cover=square_grid(c,covering['Phi']); K={}
        for i in range(5):
            for k in range(5-i):K[gridkey(i,k)]=rho*square_cover[gridkey(0,k)] if i==0 else square[gridkey(i-1,k)]
        coefficients=(2*source['z']/source['L'],source['z']*(-(1-self.delta))/source['L'],-source['d']/source['L'])
        grids={Q:{},F:{},V:{},PD:{},PI:{}}
        for i in range(5):
            for k in range(5-i):
                index=gridkey(i,k)
                grids[Q][index]=sum((math.comb(k,j)*math.factorial(j)*(coefficients[0][j]*vel[gridkey(i,k-j)]
                    +coefficients[1][j]*mean[gridkey(i,k-j)]+coefficients[2][j]*mean[gridkey(i,k-j+1)]) for j in range(k+1)),c.mpf(0))
                grids[F][index]=sum((math.comb(k,j)*source['relative_F0_derivatives'][j]*phi[gridkey(i,k-j)] for j in range(k+1)),c.mpf(0))
                grids[V][index]=vel[index]
                grids[PD][index]=source['pressure'][k]*math.factorial(k) if i==0 else c.mpf(0)
                grids[PI][index]=sum((math.comb(k,j)*source['relative_F0_squared_derivatives'][j]*K[gridkey(i,k-j)] for j in range(k+1)),c.mpf(0))
        return dict(rho=rho,Z=z,ordinary_mixed_profile_grids=grids,
            Phi_derivatives=phi,Uz_derivatives=vel,radial_average_Uz_derivatives=mean,
            explicit_model_factorial_tail=packet['explicit_model_factorial_tail'],nonlinear_correction_bounds=packet['nonlinear_correction_bounds'],
            F0_exact_positive_log_enclosure=c.mpf([endpoints(-self.logC-self.Lambda*self.Gbar)[0],endpoints(-self.logC)[1]]),
            actual_selected_Cstar_source_retained=True,nonlinear_field_enclosed_not_replaced_by_model=True,
            original_P0_and_centrifugal_pressure_increment_retained=True)

    def physical_map(self,packet,axis=False):
        c=self.ctx; grids=packet['ordinary_mixed_profile_grids']; z=packet['Z']
        x=y=c.mpf(0) if axis else c.mpf([-endpoints(c.sqrt(c.mpf('8.2')))[1],endpoints(c.sqrt(c.mpf('8.2')))[1]])
        beta={Q:c.mpf(-1),F:-1-self.delta,V:-1-self.delta,PD:-2-2*self.delta,PI:-2-2*self.delta}
        rows={}; scales={}
        for nx,ny,nz in INDICES:
            index='x'+str(nx)+'_y'+str(ny)+'_z'+str(nz); rows[index]={}; N=nx+ny
            for component,bases in BASES.items():
                rows[index][component]={}
                for label in bases:
                    value=sum((coefficient_value(c,coefficient,x,y,z,self.delta,beta[label])*grids[label][gridkey(i,k)]
                               for (i,k),coefficient in core_templates()[(component,label,nx,ny,nz)].items()),c.mpf(0))
                    norm=c.mpf(max(abs(v) for v in endpoints(value))); zero=endpoints(norm)[1]==0
                    scale_key=label+'_transverse'+str(N)+'_z'+str(nz)
                    if scale_key not in scales:
                        lp=(c.mpf(N-1)/2 if label in (Q,F) else c.mpf(N)/2-(1 if label==PI else 0))*self.logLambda
                        amp=-self.logC if label==F else -2*self.logC if label==PI else 2*self.logP if label==PD else c.mpf(0)
                        gamma=beta[label]-N+nz*(self.delta-1)
                        scales[scale_key]=dict(logLambda_term=lp,amplitude_log_upper=amp,physical_lambda_exponent=gamma,
                            log_tau_sector_terms={lt:gamma*c.mpf(lt)/2 for lt in ('-1','-10','-100')},
                            formal_F0_retained_for_swirl_and_pressure=True)
                    rows[index][component][label]=dict(bracket=value,absolute_upper=norm,exactly_zero=zero,scale_key=scale_key)
        times={}
        for component,bases in BASES.items():
            times[component]={}
            for label,seed in bases.items():
                value=sum((coefficient_value(c,coef,x,y,z,self.delta,beta[label])*grids[label][gridkey(i,k)]
                           for (i,k),coef in core_time_template(seed).items()),c.mpf(0))
                norm=c.mpf(max(abs(v) for v in endpoints(value)))
                source_scale=scales[label+'_transverse0_z0']
                gamma=beta[label]-2
                times[component][label]=dict(bracket=value,absolute_upper=norm,exactly_zero=endpoints(norm)[1]==0,
                    logLambda_term=source_scale['logLambda_term'],amplitude_log_upper=source_scale['amplitude_log_upper'],
                    physical_lambda_exponent=gamma,log_tau_sector_terms={lt:gamma*c.mpf(lt)/2 for lt in ('-1','-10','-100')})
        return dict(cartesian_spatial_multiindices=rows,shared_physical_prefactor_bounds=scales,
                    first_fixed_x_physical_time_derivative=times,
                    includes_axis_without_inverse_radius=True,axis_only=axis)

    def local_energy(self,whole):
        c=self.ctx; grids=whole['ordinary_mixed_profile_grids']
        norm=lambda value:c.mpf(max(abs(v) for v in endpoints(value)))
        q=norm(grids[Q][gridkey(0,0)]); phi=norm(whole['Phi_derivatives'][gridkey(0,0)]); v=norm(grids[V][gridkey(0,0)])
        rho=c.mpf('4.1'); alpha_r=(3-self.delta)/2; alpha_p=(3-3*self.delta)/2
        masses=dict(radial=dict(logLambda_term=-2*self.logLambda,amplitude_term=c.mpf(0),
                bounded_profile_term=2*c.ln(rho*q)-c.ln(4)),
            swirl=dict(logLambda_term=-2*self.logLambda,amplitude_term=-2*self.logC,
                bounded_profile_term=2*c.ln(rho*phi)),
            axial=dict(logLambda_term=-self.logLambda,amplitude_term=c.mpf(0),
                bounded_profile_term=c.ln(rho)+2*c.ln(v)))
        rows={}
        for name in masses:
            alpha=alpha_r if name=='radial' else alpha_p; power=alpha-1
            rows[name]=dict(radial_integral_log_upper_parts=masses[name],
                fixed_physical_strip_energy_log_terms={lt:c.ln(2*c.pi)+alpha*c.ln(3)-c.mpf(lt) for lt in ('-1','-10','-100')},
                similarity_sector_energy_log_terms={z:{lt:c.ln(2*c.pi*c.mpf(z))+power*c.mpf(lt)-alpha*c.ln(1-c.mpf(z)**2)
                    for lt in ('-1','-10','-100')} for z in ('.5','.9','.99')},
                fixed_physical_strip_spacetime_energy_log_term=c.ln(2*c.pi)+alpha*c.ln(3)+c.ln(99))
        return dict(component_contributions=rows,
            domain='0<=rho<=4.1 (core including short analytic extension), |z_phys|<=1, fixed positive tau; spacetime exp(-100)<=tau<=exp(-1)',
            source_inequalities='Ir<=epsilon^2*rho_max^2*supQ^2/4; Itheta<=epsilon^2*rho_max^2*supPhi^2*F0_upper^2; Iz<=epsilon*rho_max*supUz^2',
            core_fixed_positive_time_local_kinetic_energy_bounded=True,
            core_local_spacetime_energy_away_from_terminal_time_bounded=True,
            uniform_terminal_time_energy_certified=False,full_background_physical_energy_integral_certified=False)

    def report(self):
        c=self.ctx; whole=self.profiles([0,'4.1'],[-1,1]); axis=self.profiles(0,[-1,1])
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            datum_enclosure_sha256=self.datum.datum_sha,uniform_Cstar_family_sha256=self.records['physical_norm_family']['uniform_Cstar_family_sha256'],
            selected_logCstar=self.logC,actual_Lambda=self.Lambda,actual_epsilon_core=self.epsilon,actual_delta=self.delta,
            normalized_core_definition='rho=Lambda*R; F=F0*Phi, Uz=4Z+j+epsilon*Psi, F0=exp(-selected_logCstar-Lambda*G); unique admitted analytic fixed point',
            model_definition='Phi0=B(chi*rho/2), Psi0=-rho*g/(2L); Phi/Psi = model + actual admitted Xh correction',
            radial_recovery='Ur=sqrt(R/2)*Q, Q=[2Z*Uz-(1-delta)Z*(Mz/R)-d*dZ(Mz/R)]/L',
            pressure_recovery='P=Pstar^2*P0_normalized+epsilon*F0^2*V(Phi^2)',
            nonsingular_cartesian_definition='X=sqrt(Lambda)*x/lambda,Y=sqrt(Lambda)*y/lambda; ux=sqrt(epsilon)*[lambda^-1*X*Q/2-lambda^(-1-delta)*Y*F0*Phi], uy analogous',
            core_domain='rho in[0,4.1], Z in[-1,1]; physical map tau>0, |Z|<1; Ra=4/Lambda',
            whole_core=whole,whole_axis=axis,whole_core_physical_map=self.physical_map(whole),whole_axis_physical_map=self.physical_map(axis,True),
            core_local_physical_energy=self.local_energy(whole),
            samples=[self.profiles(r,z) for z in ('-1','0','.5','1') for r in ('0','2','4')],
            original_axis_values_and_first_radial_slopes_restored=True,
            whole_Z_analytic_core_profile_enclosures_through_order5_available=True,
            whole_core_and_axis_cartesian_spatial4_enclosures_available=True,
            whole_core_and_axis_first_physical_time_derivative_enclosures_available=True,
            core_local_physical_energy_bounds_available=True,
            local_old_finite_coefficient_state_used=False,selected_point_coefficients_recomputed=False,
            old_finite_coefficients_extrapolated=False,actual_F0_not_materialized=True,
            core_inner_annulus_interfaces_certified=False,full_cartesian_vector_derivatives_certified=False,
            physical_energy_integral_certified=False,whole_outer_cone_certified=False,temporal_recursion=False,
            next_dependency='Core physical local energy and inner/annulus source dispatch/interfaces; whole-field physical divergence/residual; admissible stress and independent flat remainder',input_hashes=self.hashes)


def run():
    with mp.workdps(280):result=CompliantCorePhysicalField().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(result),indent=2)+'\n',encoding='utf8')
    print('Same selected-Cstar whole-Z analytic core/axis jets and nonsingular Cartesian spatial4 enclosures generated',flush=True)
    return result


if __name__=='__main__':run()
