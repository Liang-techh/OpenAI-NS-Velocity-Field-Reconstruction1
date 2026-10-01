"""Same-family O.1/O.2 outer field in local normalized moment units.

The inlet is the actual functional repair's exact terminal identities.
Rref stays formal; all velocities, moments, pressure and stress quantities
use its selected family. Huge axial turnoff masses use endpoint-weighted
positive kernels with a retained nonzero tail, never exp(exp(logCstar)).
This module does not install the later outer moment repairs/heat exterior.
"""
import hashlib
import json
from fractions import Fraction
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_shared_five_moment_repair import SharedFiveMomentRepair, pack
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_interval_outer_slope_field import transition_integrals, fraction_box
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_'


def stable_sigma(c,t):
    """The existing flat cutoff, evaluating both sides without subtraction."""
    lo,hi=endpoints(c.mpf(t))
    def point(v):
        if v<=0:return c.mpf(0),c.mpf(0)
        if v>=1:return c.mpf(1),c.mpf(0)
        x=c.mpf(v);a=c.exp(-1/x**2);b=c.exp(-1/(1-x)**2)
        value=a/(a+b);derivative=value*(b/(a+b))*(2/x**3+2/(1-x)**3)
        return value,derivative
    vl,dl=point(lo);vh,dh=point(hi)
    value=c.mpf([max(mp.mpf(0),endpoints(vl)[0]),min(mp.mpf(1),endpoints(vh)[1])])
    if lo==hi:
        derivative=c.mpf([max(mp.mpf(0),endpoints(dl)[0]),min(mp.mpf(8),endpoints(dh)[1])])
    else:derivative=c.mpf([0,8])
    return value,derivative


def turnoff_kernels(c,y,Md,cells=256,window=800):
    """Enclose int_1^y exp(s-y) B(log(s)/Md)^j ds, j=1,2.

    Substitute r=y-s. Each cell uses its exact exponential mass and
    monotone endpoint bounds for B. The omitted far interval is bounded
    by its positive exponential mass; it is never set to zero.
    """
    if not isinstance(cells,int) or cells<1 or not isinstance(window,int) or window<1:
        raise ValueError('Positive integer cells and window required')
    y=c.mpf(y);md=c.mpf(Md)
    if endpoints(y)[0]<1:raise ValueError('y>=1 required')
    def B(s):
        sl,sh=endpoints(s)
        def val(x):
            if x<=1:return c.mpf(1)
            phase=1-c.ln(c.mpf(x))/md
            return stable_sigma(c,phase)[0]
        return c.mpf([endpoints(val(sh))[0],endpoints(val(sl))[1]])
    # Integrate to the lower bound of the upper endpoint. Any endpoint
    # uncertainty is covered separately, so interval y is supported.
    length=y-1;ll,lh=endpoints(length);limit=min(ll,mp.mpf(window))
    mass=[c.mpf(0),c.mpf(0)]
    for i in range(cells):
        a=c.mpf(limit)*i/cells;b=c.mpf(limit)*(i+1)/cells
        left=B(y-b);right=B(y-a)
        box=c.mpf([max(mp.mpf(0),endpoints(right)[0]),min(mp.mpf(1),endpoints(left)[1])])
        weight=c.exp(-a)-c.exp(-b)
        mass[0]+=weight*box;mass[1]+=weight*box**2
    tail_upper=c.exp(-c.mpf(limit))-c.exp(-c.mpf(lh))
    tail=c.mpf([0,max(mp.mpf(0),endpoints(tail_upper)[1])])
    return dict(B_mass=mass[0]+tail,B_squared_mass=mass[1]+tail,
                retained_far_tail=tail,window=window,cells=cells,
                exponential_cell_weights_integrated_exactly=True,
                integrals_normalized_at_current_radius=True)


class SharedOuterInitial:
    def __init__(self):
        self.repair=SharedFiveMomentRepair();self.ctx=c=self.repair.ctx
        names=('shared_five_moment_repair','shared_five_moment_repair_check')
        self.hashes=dict(self.repair.hashes)
        self.records={n:json.loads((HERE/(PREFIX+n+'.json')).read_bytes()) for n in names}
        for n,record in self.records.items():
            for source,digest in record['input_hashes'].items():
                if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:
                    raise ValueError('Selected outer inlet dependency changed: '+source)
                self.hashes[source]=digest
            self.hashes[PREFIX+n+'.json']=hashlib.sha256((HERE/(PREFIX+n+'.json')).read_bytes()).hexdigest()
        repair,check=[self.records[n] for n in names]
        if not (repair['actual_implicit_functional_five_moment_identities_analytically_certified']
                and check['actual_implicit_functional_five_moment_closure_independently_checked']
                and repair['actual_five_defect_family_sha256']==check['actual_five_defect_family_sha256']):
            raise ValueError('Actual same-family repaired terminal identities required')
        self.family=repair['actual_five_defect_family_sha256']
        self.datum=self.repair.datum;self.params=self.datum.parameters
        self.delta=self.repair.delta;self.invP2=self.repair.invP2
        self.logC=read_interval(c,self.repair.records['shared_physical_norm_family']['selected_logCstar'])
        for name in (Path(__file__).name,PREFIX+'interval_outer_slope_field.py',
                     PREFIX+'interval_long_reshape_field.py',PREFIX+'interval_taylor.py'):
            self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()

    def coordinate(self,Z):
        c=self.ctx;Z=c.mpf(Z)
        if endpoints(Z)[0]<-1 or endpoints(Z)[1]>1:raise ValueError('Z in [-1,1] required')
        z=IntervalTaylor(c,[Z,c.mpf(1)])
        q=IntervalTaylor(c,[1+Z**2,2*Z]);return z,q.reciprocal()

    def packet(self,Z,y,u,uy,V,Vy,m,h,k,e,p,stage,extra=None):
        c=self.ctx;z,qi=self.coordinate(Z);Z=z[0];dz=1-Z**2;L=1-self.delta*Z**2
        jets=(u,uy,V,Vy,m,h,k,e,p)
        if any(not isinstance(v,IntervalTaylor) or v.order!=1 for v in jets):
            raise TypeError('Normalized outer axial jet lost')
        if endpoints(L)[0]<=0 or endpoints(u[0])[0]<=0:
            raise ArithmeticError('Outer L/u positivity unresolved; subdivide the input box')
        p0rows=self.datum.normalized_jets(endpoints(Z),1)['normalized_pressure_coefficients']
        p0=IntervalTaylor(c,[c.mpf(endpoints(v)) for v in p0rows]);P=p+p0
        W=1-(1-self.delta)*Z*m[0]-dz*m[1]
        angular=(1-self.delta/2)*h[0]-(1-self.delta)*Z*h[1]/2-dz*k[1]+(2*self.delta-1)*Z*k[0]
        Q=-W+angular/u[0]
        N=(-W*V[0]+(1-self.delta)*(m[0]-Z*m[1])/2)*self.invP2
        N+=2*self.delta*Z*e[0]-dz*e[1]+2*(1+self.delta)*Z*P[0]-dz*P[1]
        Ur=(2*Z*V[0]-(1-self.delta)*Z*m[0]-dz*m[1])/L
        a=1-2*uy[0]/u[0];bP=2*Vy[0]/u[0]
        if endpoints(a)[0]<=0:
            raise ArithmeticError('Outer angular-shear denominator unresolved')
        # Q can have a broad interval on wide axial boxes; do not turn a
        # diagnostic denominator into an unsupported cone certificate.
        nonzeroQ=not endpoints(Q)[0]<=0<=endpoints(Q)[1]
        bw=2*N*Vy[0]/(u[0]**2*Q) if nonzeroQ else None
        kappa=a+bP**2*self.invP2/a
        result=dict(stage=stage,y=y,Z=Z,Utheta_over_Pstar=list(u.coefficients),
            Utheta_y_over_Pstar=list(uy.coefficients),Uz=list(V.coefficients),Uz_y=list(Vy.coefficients),
            Mz_over_R=list(m.coefficients),Mtheta_over_sqrt2_R_3half_Pstar=list(h.coefficients),
            Mtheta_z_over_sqrt2_R_3half_Pstar=list(k.coefficients),Mztheta_over_R_Pstar_squared=list(e.coefficients),
            Mp_over_Pstar_squared=list(p.coefficients),P0_over_Pstar_squared=list(p0.coefficients),
            P_over_Pstar_squared=list(P.coefficients),Ur_over_sqrt_R_over_2=Ur,
            Q=Q,N_over_Pstar_squared=N,angular_shear_a=a,axial_shear_b_times_Pstar=bP,
            axial_stress_product_bw=bw,kappa=kappa,Ur_Z_available=False,
            normalized_local_radial_moment_rhs=dict(m=list((V-m).coefficients),h=list((u-h*c.mpf('1.5')).coefficients),
                k=list((u*V-k*c.mpf('1.5')).coefficients),
                e=list((V*V*self.invP2-u*u/2-e).coefficients),p=list((u*u/2).coefficients)),
            terminal_identity_source='actual F30 functional five-moment closure, not midpoint controls',
            pressure_schedule_Md=self.params.Md,implicit_source_sha256=self.datum.source_sha,
            source_pressure_not_changed=True,physical_Rref_not_materialized=True,
            all_five_moment_histories_retained=True,retained_axial_order=1,
            outer_moment_repairs_installed=False,whole_outer_cone_certified=False,
            exact_heat_exterior_installed=False,temporal_recursion=False)
        if extra:result.update(extra)
        return result

    def reference(self,Z,offset):
        c=self.ctx;y=fraction_box(c,offset)
        if endpoints(y)[0]<endpoints(c.ln(2)-6)[0] or endpoints(y)[1]>0:
            raise ValueError('ln(2)-6 <= relative reference offset <=0 required')
        z,qi=self.coordinate(Z);u=qi*c.exp(y/10);V=z*4
        m=V;h=u*c.mpf('.625');k=h*z*4
        e=z*z*(16*self.invP2)-u*u*c.mpf(5)/12;p=u*u*c.mpf('2.5')
        return self.packet(Z,y,u,u/10,V,z*0,m,h,k,e,p,'O.1 selected-family reference')

    def slope(self,Z,y,cells=256):
        c=self.ctx;yy=fraction_box(c,y);z,qi=self.coordinate(Z)
        J,integrals=transition_integrals(c,y,cells)
        factor=c.exp(yy/10-c.mpf('.6')*J);slope=c.mpf('.1')-c.mpf('.6')*stable_sigma(c,yy)[0]
        u=qi*factor;V=z*4;m=V
        h=qi*((c.mpf('.625')+integrals[0])*c.exp(-c.mpf('1.5')*yy));k=h*z*4
        e=z*z*(16*self.invP2)-qi*qi*((c.mpf(5)/12+integrals[2]/2)*c.exp(-yy))
        p=qi*qi*(c.mpf('2.5')+integrals[1]/2)
        return self.packet(Z,yy,u,u*slope,V,z*0,m,h,k,e,p,'O.2 selected-family slope transition',
                           dict(J=J,dimensionless_increment_integrals=integrals,directed_integral_cells=cells))

    def axial(self,Z,phase=None,buffer_offset=None,cells=256,window=800):
        c=self.ctx;z,qi=self.coordinate(Z);md=c.mpf(self.params.Md)
        if (phase is None)==(buffer_offset is None):raise ValueError('Choose phase or buffer_offset')
        if phase is not None:
            qq=fraction_box(c,phase)
            if endpoints(qq)[0]<0 or endpoints(qq)[1]>1:raise ValueError('phase in [0,1] required')
            y=c.exp(md*qq);B,dB=stable_sigma(c,1-qq);By=-dB/(md*y)
        else:
            offset=fraction_box(c,buffer_offset)
            if endpoints(offset)[0]<0 or endpoints(offset)[1]>11:raise ValueError('buffer offset in [0,11] required')
            y=c.exp(md)+offset;B=c.mpf(0);By=c.mpf(0)
        inlet=self.slope(Z,1,cells);t=y-1;decay=c.exp(-t);root_decay=c.exp(-t/2);decay3=c.exp(-c.mpf('1.5')*t)
        get=lambda key:IntervalTaylor(c,inlet[key])
        u1=get('Utheta_over_Pstar');u=u1*root_decay;V=z*(4*B);Vy=z*(4*By)
        K=turnoff_kernels(c,y,self.params.Md,cells,window)
        m=get('Mz_over_R')*decay+z*(4*K['B_mass'])
        h=get('Mtheta_over_sqrt2_R_3half_Pstar')*decay3+u1*(root_decay-decay3)
        k=get('Mtheta_z_over_sqrt2_R_3half_Pstar')*decay3+u1*z*(4*root_decay*K['B_mass'])
        e=get('Mztheta_over_R_Pstar_squared')*decay+z*z*(16*self.invP2*K['B_squared_mass'])-u1*u1*(t*decay/2)
        p=get('Mp_over_Pstar_squared')+u1*u1*((1-decay)/2)
        return self.packet(Z,y,u,-u/2,V,Vy,m,h,k,e,p,'O.2 selected-family axial turnoff / buffer',
            dict(phase=str(phase) if phase is not None else None,buffer_offset=str(buffer_offset) if buffer_offset is not None else None,
                 cutoff=B,cutoff_y=By,normalized_cutoff_integrals=K,
                 accumulated_axial_moment_not_reset_after_Uz_zero=True))

    def report(self):
        with mp.workdps(210):
            samples=[self.reference('.5',v) for v in ('-5','-1','0')]
            samples += [self.slope('.5',v) for v in ('0','.5','1')]
            samples += [self.axial('.5',phase=v) for v in ('0','.5','1')]
            samples += [self.axial('.5',buffer_offset='11')]
            return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.datum.source_sha,
                datum_enclosure_sha256=self.datum.datum_sha,actual_delta=self.delta,selected_logCstar=self.logC,
                selected_radius_definition='Rref=110*(Cstar*Pstar)^10; symbolic exact relative offsets',
                normalized_moment_units=['Mz/R','Mtheta/(sqrt(2)*R^1.5*Pstar)',
                    'Mtheta_z/(sqrt(2)*R^1.5*Pstar)','Mztheta/(R*Pstar^2)','Mp/Pstar^2'],
                real_axial_domain=['-1','1'],samples=samples,
                same_family_reference_and_initial_outer_callable=True,
                same_analytic_preheat_pressure_and_actual_delta_retained=True,
                independent_outer_repairs_and_heat_still_required=True,
                complete_corrected_outer_built=False,heat_exterior_matched=False,
                global_admissible_stress_lift_constructed=False,temporal_recursion=False,
                input_hashes=self.hashes)


def run():
    field=SharedOuterInitial();result=field.report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf-8')
    print('Same-family reference, outer slope and axial turnoff callable: 10 samples; positive kernel tails retained',flush=True)
    return result


if __name__=='__main__':run()
