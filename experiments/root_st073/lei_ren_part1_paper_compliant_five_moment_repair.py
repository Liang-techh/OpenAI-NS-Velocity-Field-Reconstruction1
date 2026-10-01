"""Actual implicit functional five-bump repair and callable profile enclosures.

Signed defects enclose the true source primitives on the entire axial domain.
The numerical inverse encloses their unique coefficient family; it does not
replace unresolved source functions by midpoints or zero tails. Fields,
partial moments, pressure and normalized stresses retain the same data.
"""
# Recomputed for the distinct compliant pressure source.
# Formula origin: lei_ren_part1_paper_shared_five_moment_repair.py; legacy source/receipts remain unchanged.

import hashlib
import json
from fractions import Fraction
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_five_defect_admission import signed_jets
from lei_ren_part1_paper_interval_five_bump_inverse import certify,weights
from lei_ren_part1_paper_interval_exit_continuation_enclosure import restore_value
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_compliant_pressure_source import CompliantPressureDatum as LogarithmicPressureDatum
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_'


class SharedFiveMomentRepair:
    def __init__(self):
        self.ctx=c=MPIntervalContext();c.dps=160
        names=['compliant_five_defect_admission','shared_bump_constants','compliant_reference_join_bounds',
               'compliant_physical_norm_family','shared_analytic_tube']
        self.records={n:json.loads((HERE/(PREFIX+n+'.json')).read_bytes()) for n in names}
        self.hashes={}
        for record in self.records.values():
            for n,digest in record['input_hashes'].items():
                if hashlib.sha256((HERE/n).read_bytes()).hexdigest()!=digest:
                    raise ValueError('Actual repair dependency changed: '+n)
                self.hashes[n]=digest
        self.admit=self.records['compliant_five_defect_admission']
        fixed=self.records['shared_bump_constants'];join=self.records['compliant_reference_join_bounds']
        norms=self.records['compliant_physical_norm_family']
        if (not norms['radius_9_17_compatibility_certified_for_selected_analytic_family']
                or endpoints(read_interval(c,norms['radius_log_margins']['Rm16_log_margin']))[0]<=0
                or norms['uniform_Cstar_family_sha256']!=self.admit['uniform_Cstar_family_sha256']):
            raise ValueError('Selected same-family Rm>=16 proof required')
        if (not self.admit['complete_10_19_actual_functional_defect_test_certified']
                or self.admit['reference_join_family_sha256']!=join['reference_join_family_sha256']):
            raise ValueError('Actual functional defects not admitted')
        wn=PREFIX+'bump_integral_enclosures_check.json'
        self.wraw=json.loads((HERE/wn).read_bytes())
        if not self.wraw['directed_integrals_certified']:
            raise ValueError('True fixed bump integrals not admitted')
        for n,key in ((PREFIX+'bump_integral_enclosures.py','module_sha256'),
                      (PREFIX+'bump_integral_enclosures_check.py','check_source_sha256')):
            if hashlib.sha256((HERE/n).read_bytes()).hexdigest()!=self.wraw[key]:
                raise ValueError('Fixed bump quadrature changed')
            self.hashes[n]=self.wraw[key]
        self.datum=LogarithmicPressureDatum('40',precision=160)
        if (self.datum.source_sha!=self.admit['implicit_source_sha256']
                or self.datum.datum_sha!=self.admit['datum_enclosure_sha256']):
            raise ValueError('Same axis pressure datum changed')
        self.hashes.update(self.datum.input_hashes)
        with mp.workdps(200):
            get=lambda r,n:read_interval(c,r[n])
            hi=lambda v:c.mpf(endpoints(v)[1])
            self.j=get(self.admit,'required_j')
            self.rho=get(self.admit,'actual_v1_minus_4Z_minus_j_C1_upper')
            self.logP=get(self.admit,'logPstar');self.invP2=get(self.admit,'inverse_Pstar_squared')
            self.kernels={n:read_interval(c,v) for n,v in self.admit['exact_signed_kernel_enclosures'].items()}
            self.tail=get(self.admit,'actual_tail_C1_positive_cap')
            self.e=get(self.admit,'complete_actual_functional_e_upper')
            self.CA=fixed['CA'];self.CQ=fixed['CQ'];self.CS=fixed['CS'];self.KN=join['KN']
            self.t=hi(2*self.CA*self.e)
            self.tstar=get(fixed,'t_star')
            if endpoints(self.t)[1]>=endpoints(self.tstar)[0]:
                raise ArithmeticError('Corrected coefficient C1 t_star gate failed')
            # The tube stores an upper bound, not the physical parameter.
            # Recover delta from the unchanged source datum, as the core did.
            self.delta=c.mpf(endpoints(self.datum.parameters.delta))
            tube_delta=get(self.records['shared_analytic_tube'],'delta')
            if not (0<endpoints(self.delta)[0]<=endpoints(self.delta)[1]
                    <=endpoints(tube_delta)[1]<=mp.mpf('.001')):
                raise ValueError('Fixed CS elementary delta hypothesis failed')
            self.W=weights(c,self.wraw)
            if any(endpoints(v)!=endpoints(read_interval(c,fixed['fixed_matrix'][i][j]))
                   for i,row in enumerate(self.W['L']) for j,v in enumerate(row)):
                raise ValueError('Actual correction matrix differs from fixed CA matrix')
            # Prove anisotropic box support, separately from the scalar e gate.
            self.axial_scale=self.t
            self.angular_scale_coefficient=10**10*self.axial_scale**2
            self.angular_scale=self.angular_scale_coefficient*self.invP2
            eta=get(self.admit,'actual_axial_eta_C1_upper')
            source_coeff=12*c.exp(c.mpf('1.2'))*self.kernels['d4']*eta**2+3*c.mpf('1e-800')
            self.angular_source_ratio=hi(self.CA*source_coeff/self.angular_scale_coefficient)
            if endpoints(self.angular_source_ratio)[1]>=mp.mpf('.1'):
                raise ArithmeticError('Angular forcing not covered by anisotropic coefficient box')
            self.scales=[self.axial_scale]*2+[self.angular_scale]*3
            self.normalization=restore_value(c,self.wraw['normalization'])
            self.beta_sup=get(fixed,'bump_sup_bound');self.beta_deriv_sup=get(fixed,'bump_first_derivative_sup_bound')
            self.full_weights={(Fraction(r['center']),Fraction(r['power']),r['multiplicity']):restore_value(c,r['weight_interval'])
                               for r in self.wraw['weight_records'].values()}
        self.hashes.update({PREFIX+n+'.json':hashlib.sha256((HERE/(PREFIX+n+'.json')).read_bytes()).hexdigest() for n in names})
        self.hashes.update({n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in
            (Path(__file__).name,wn,PREFIX+'interval_five_bump_inverse.py',PREFIX+'compliant_five_defect_admission.py')})

    def coefficients(self,Z,tightening=12):
        c=self.ctx
        with mp.workdps(200):
            Z=c.mpf(Z)
            d,invAm2=signed_jets(c,Z,self.j,self.rho,self.logP,self.kernels,self.tail)
            result=certify(c,self.W,d,invAm2,scales=self.scales,tightening=tightening)
            if not result['certified']:
                raise ArithmeticError('True actual five-bump map failed strict inclusion/contraction')
            return result,d,invAm2

    def beta(self,offset):
        c=self.ctx;offset=c.mpf(offset);radius=c.mpf('.025')
        left,right=endpoints(offset)
        if right<=-endpoints(radius)[1] or left>=endpoints(radius)[1]:
            return c.mpf(0),c.mpf(0)
        if left<=-endpoints(radius)[0] or right>=endpoints(radius)[0]:
            return c.mpf([0,endpoints(self.beta_sup)[1]]),c.mpf([-endpoints(self.beta_deriv_sup)[1],endpoints(self.beta_deriv_sup)[1]])
        den=1-(offset/radius)**2
        value=c.exp(-1/den)/(radius*self.normalization)
        derivative=-2*offset*value/(radius**2*den**2)
        return value,derivative

    def partial_weight(self,x,center,power,multiplicity=1,cells=64):
        """Exact full integral or directed closed-cell partial integral."""
        c=self.ctx;x=c.mpf(x);center=Fraction(center);power=Fraction(power)
        cc=c.mpf(center.numerator)/center.denominator;radius=c.mpf('.025')
        def endpoint(upper):
            upper=c.mpf(upper);left=cc-radius;right=cc+radius
            if endpoints(upper)[1]<=endpoints(left)[0]:
                return c.mpf(0)
            if endpoints(upper)[0]>=endpoints(right)[1]:
                if power==0 and multiplicity==1:return c.mpf(1)
                return self.full_weights[center,power,multiplicity]
            # Include the tiny directed uncertainty of the support endpoint.
            length=upper-left
            if endpoints(length)[0]<0:
                length=c.mpf([0,max(mp.mpf(0),endpoints(length)[1])])
            total=c.mpf(0)
            pp=c.mpf(power.numerator)/power.denominator
            for i in range(cells):
                a=left+length*i/cells;b=left+length*(i+1)/cells
                box=c.mpf([endpoints(a)[0],endpoints(b)[1]])
                val,_=self.beta(box-cc)
                total+=box**pp*val**multiplicity*length/cells
            upper_bound=(c.mpf(1) if power==0 and multiplicity==1
                         else self.full_weights[center,power,multiplicity])
            return c.mpf([max(mp.mpf(0),endpoints(total)[0]),
                          min(endpoints(total)[1],endpoints(upper_bound)[1])])
        a,b=endpoints(x);low=endpoint(a);high=endpoint(b)
        return c.mpf([endpoints(low)[0],endpoints(high)[1]])

    def evaluate_patch(self,x,Z,partial_cells=64):
        """Enclose corrected fields and actual partial moments in scaled units.

        R=Rm*x is left implicit. Ur is divided bysqrt(R/2); u and pressure
        byPstar andPstar^2. Only C1 Z data is retained, so Ur_Z is unavailable.
        """
        c=self.ctx
        with mp.workdps(200):
            x=c.mpf(x);Z=c.mpf(Z)
            if endpoints(x)[0]<1 or endpoints(x)[1]>endpoints(c.exp(1))[1]:
                raise ValueError('Reference repair interval1<=x<=e required')
            inverse,d,invAm2=self.coefficients(Z)
            h=inverse['controls'];z=IntervalTaylor(c,[Z,c.mpf(1)])
            q=IntervalTaylor(c,[1+Z**2,2*Z])
            am=q.reciprocal()*c.exp(c.mpf('-.6'))
            f=z*0;fx=z*0;g=z*0;gx=z*0
            centers=('5/4','3/2','7/4')
            vals=[self.beta(x-c.mpf(Fraction(k).numerator)/Fraction(k).denominator) for k in centers]
            for k in range(3):
                f+=h[k+2]*vals[k][0];fx+=h[k+2]*vals[k][1]
            for i,k in enumerate((0,2)):
                g+=h[i]*vals[k][0];gx+=h[i]*vals[k][1]
            H=f+x**c.mpf('.1')
            u=am*H;uy=am*(fx*x+x**c.mpf('.1')/10)
            V=z*4+g;Vy=gx*x
            cache={}
            def w(k,p,m=1):
                key=(k,p,m)
                if key not in cache:
                    cache[key]=self.partial_weight(x,centers[k],p,m,partial_cells)
                return cache[key]
            dg=sum((h[i]*w(k,'0') for i,k in enumerate((0,2))),z*0)
            dt=sum((h[k+2]*w(k,'.5') for k in range(3)),z*0)
            dm=sum((h[i]*w(k,'.6')+h[i]*h[k+2]*w(k,'.5',2)
                    for i,k in enumerate((0,2))),z*0)
            de=sum((invAm2*h[i]*h[i]*w(k,'0',2) for i,k in enumerate((0,2))),z*0)
            de-=sum((h[k+2]*w(k,'.1')+h[k+2]*h[k+2]*w(k,'0',2)/2 for k in range(3)),z*0)
            dp=sum((h[k+2]*w(k,'-.9')+h[k+2]*h[k+2]*w(k,'-1',2)/2 for k in range(3)),z*0)
            mass=z*4+(d[0]+dg)/x
            theta=d[2]+dt+x**c.mpf('1.6')*c.mpf('.625')
            mixed=theta*z*4+d[1]+dm
            energy=d[3]+de-x**c.mpf('1.2')*c.mpf(5)/12
            pressure_moment=d[4]+dp+x**c.mpf('.2')*c.mpf('2.5')
            if not all(isinstance(v,IntervalTaylor) and v.order==1
                       for v in (H,u,uy,V,Vy,mass,theta,mixed,energy,pressure_moment)):
                raise TypeError('Axial jet structure lost in reference moment recovery')
            pressure=self.datum.normalized_jets(endpoints(Z),1)
            p0=IntervalTaylor(c,[c.mpf(endpoints(v)) for v in pressure['normalized_pressure_coefficients']])
            P=p0+am*am*pressure_moment
            dz=1-Z**2;L=1-self.delta*Z**2
            if any(endpoints(v)[0]<=0 for v in (L,H[0],u[0])):
                raise ArithmeticError('Patch denominator unresolved; subdivide the input box')
            Ur=(2*Z*V[0]-(1-self.delta)*Z*mass[0]-dz*mass[1])/L
            W=1-(1-self.delta)*Z*mass[0]-dz*mass[1]
            zeta_am=-2*Z/(1+Z**2)
            angular=((1-self.delta/2)*theta[0]-(1-self.delta)*Z*(theta[1]+zeta_am*theta[0])/2
                    -dz*(mixed[1]+zeta_am*mixed[0])+(2*self.delta-1)*Z*mixed[0])
            Q=-W+angular/(x**c.mpf('1.5')*H[0])
            energy_scaled=(z*8*mass-z*z*16)*self.invP2+am*am*energy/x
            N_scaled=(-W*V[0]+(1-self.delta)*(mass[0]-Z*mass[1])/2)*self.invP2
            N_scaled+=2*self.delta*Z*energy_scaled[0]-dz*energy_scaled[1]
            N_scaled+=2*(1+self.delta)*Z*P[0]-dz*P[1]
            a=1-2*uy[0]/u[0]
            if endpoints(Q)[0]<=0 or endpoints(a)[0]<=0:
                raise ArithmeticError('Patch Q/a positivity unresolved; subdivide the input box')
            b_times_P=2*Vy[0]/u[0]
            bw=2*N_scaled*Vy[0]/(u[0]**2*Q)
            kappa=a+(b_times_P**2*self.invP2)/a
            # Rm>=16 is independently admitted; use it as a conservative
            # physical D floor, not an artificially selected actual radius.
            cone_margin=16*x*(Q/L)*(a-bw)-2*a
            return dict(x=x,Z=Z,Utheta_over_Pstar=list(u.coefficients),Uz=list(V.coefficients),
                Utheta_y_over_Pstar=list(uy.coefficients),Uz_y=list(Vy.coefficients),
                Ur_over_sqrt_R_over_2=Ur,Ur_Z_available=False,
                Mz_over_R=list(mass.coefficients),
                normalized_Mtheta=list(theta.coefficients),normalized_Mtheta_z=list(mixed.coefficients),
                normalized_centered_energy=list(energy.coefficients),Mp_over_Am_squared=list(pressure_moment.coefficients),
                P_over_Pstar_squared=list(P.coefficients),Q=Q,N_over_Pstar_squared=N_scaled,
                angular_shear_a=a,axial_shear_b_times_Pstar=b_times_P,axial_stress_product_bw=bw,kappa=kappa,
                relaxed_Da_minus_bw_minus_2a_lower_enclosure=cone_margin,
                partial_weights={str(k):v for k,v in cache.items()},
                actual_coefficients=inverse['controls'],actual_defects=[list(v.coefficients) for v in d],
                moment_inheritance='reference primitives plus true initial defects plus actual partial bump integrals',
                normalized_units=True,physical_Rm_or_Cstar_not_materialized=True,
                retained_axial_order=1,terminal_moments_not_reset=True,
                source_pressure_not_changed=True,temporal_recursion=False)

    def report(self):
        c=self.ctx
        with mp.workdps(200):
            whole,d,invAm2=self.coefficients([-1,1])
            CS_t=self.CS*self.t
            banach_ratio=4*self.CA**2*self.CQ*self.e
            bw=256*self.KN*CS_t
            kappa=c.mpf('.9')+(16*CS_t)**2/c.mpf('.7')
            cone_margin=8*(c.mpf('.7')-bw)-c.mpf('1.8')
            if (endpoints(CS_t)[1]>=mp.mpf('.01') or endpoints(banach_ratio)[1]>=mp.mpf('.5')
                    or endpoints(bw)[1]>=mp.mpf('.1') or endpoints(kappa)[1]>=1
                    or endpoints(cone_margin)[0]<=0):
                raise ArithmeticError('Corrected functional cone/velocity bounds failed')
            samples=[self.evaluate_patch(x,'.5') for x in ('1','1.25','1.5','1.75','2')]
            self.hashes[PREFIX+'logarithmic_pressure_datum.py']=hashlib.sha256((HERE/(PREFIX+'logarithmic_pressure_datum.py')).read_bytes()).hexdigest()
            result=dict(actual_five_defect_family_sha256=self.admit['actual_five_defect_family_sha256'],
                reference_join_family_sha256=self.admit['reference_join_family_sha256'],
                uniform_Cstar_family_sha256=self.admit['uniform_Cstar_family_sha256'],
                implicit_source_sha256=self.admit['implicit_source_sha256'],datum_enclosure_sha256=self.admit['datum_enclosure_sha256'],
                actual_source_delta_enclosure=self.delta,tube_delta_upper_not_used_as_parameter=True,
                real_axial_domain=['-1','1'],correction_support_x=['49/40','71/40'],available_interval_x=['1','e'],
                coefficient_order=['c1','c2','xi1','xi2','xi3'],actual_whole_axis_coefficient_inverse=whole,
                complete_actual_defect_C1_e_upper=self.e,coefficient_C1_t_upper=self.t,
                global_C1_banach_contraction_upper=banach_ratio,anisotropic_angular_source_ratio=self.angular_source_ratio,
                axial_box_scale=self.axial_scale,angular_box_scale=self.angular_scale,
                CS_times_coefficient_C1_bound=CS_t,corrected_bw_upper=bw,corrected_kappa_upper=kappa,
                corrected_relaxed_cone_lower_margin=cone_margin,normalized_patch_samples=samples,
                actual_implicit_functional_five_moment_repair_specified=True,
                actual_implicit_functional_five_moment_identities_analytically_certified=True,
                numerical_coefficient_and_partial_moment_enclosures_callable=True,
                corrected_inner_Ra_Rh_relaxed_cone_analytically_certified=True,
                original_core_and_inner_admissible_collar_preserved=True,axis_pressure_P0_preserved=True,
                signed_tails_or_source_data_projected_to_midpoints=False,
                terminal_moments_reset_or_pressure_tail_fitted=False,
                physical_full_field_evaluator_complete=False,exact_coefficient_point_functions_evaluated=False,
                independent_five_primitive_to_bump_map_check_pending=True,
                corrected_outer_at_selected_radius_built=False,heat_exterior_matched=False,
                admissible_stress_lift_constructed=False,full_Section9_parameter_admission=False,temporal_recursion=False,
                proof=dict(
                    actual_target='source-bound actual five axis defects in10.3 units; true scalar restoration integrals and positive retained tails, not fitted data',
                    global_C1='fixed inverse CA and quadratic CQ on the summedC1 Banach ball giveunique h with||h||C1<=2CAe; the actual smooth source data gives smooth h by invertible pointwise derivative',
                    source_norm_vs_box='the C1 e bound belongs to the correlated true source functions. Independent Taylor coefficient boxes enlarge that class; their arbitrary box norm is not used as the analytic Banach source norm. The separate directed inverse is certified on the enlarged boxes.',
                    interval_inverse='strict directed self-map and uniform pointwise contraction enclose this same actual h family andits derivative; fixed preconditioner midpoints never replace mapweights, defects oramplitude',
                    angular_scale='the explicit source coefficient covers d3,d5 andall nonaxiald4 tails plusfull Am^-2 E1^2 C1; CA times this coefficient divided bytheangularbox coefficient is<.1',
                    moments='the five actual normalized changes are10.8; h solvesA h+Q(h,h)=-d, so their centered defects vanish identically. Invertible centering restores each original moment. Numericalpartialintegrals are enclosure diagnostics; total identities follow fromthemap.',
                    fields='u/Pstar=(exp(-.6)/(1+Z^2))*(x^.1+f_h),V=4Z+g_h; all cumulativeprimitives start fromreference moments+actualdefects atRm; P=P0+Mp; Ur follows the same Mz primitive.',
                    cone='e<e_star gives t<=t_star andCS t<.01; Section10 scalarQ/N barriers retain their inheriteddata, D>=8, bw<.1 andkappa<1 imply the stored relaxed margin',
                    outer_scope='identities restore the prescribed reference target moments at2Rm..Rh; actual completed outer/heat construction atselectedRref isstillpending; no globaladmissible or temporalclaim'),
                input_hashes=self.hashes)
            return result


def pack(value):
    if isinstance(value,IntervalTaylor):return dict(coefficients=list(value.coefficients))
    if isinstance(value,dict):return {k:pack(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [pack(v) for v in value]
    return value


def run():
    repair=SharedFiveMomentRepair();result=repair.report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf-8')
    print('Actual implicit whole-axis five-moment repair and partial field enclosures PASS',flush=True)
    print('Coefficient C1 t<=',mp.nstr(endpoints(result['coefficient_C1_t_upper'])[1],14),
          'corrected |bw|<=',mp.nstr(endpoints(result['corrected_bw_upper'])[1],14),flush=True)
    return result


if __name__=='__main__':
    run()
