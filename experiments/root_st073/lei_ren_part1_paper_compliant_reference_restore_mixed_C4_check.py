"""Independent physical primitive/mixed4 checks for actual axial restore."""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
import sympy as s
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_reference_restore_mixed_C4 import (
    CompliantReferenceRestoreMixedC4,source_mixed,IntervalTaylor)
from lei_ren_part1_paper_compliant_actual_patch_mixed_C4 import CompliantActualPatchMixedC4
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import sigma_jets,positive_exp
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
NAME='lei_ren_part1_paper_compliant_reference_restore_mixed_C4.json'


def structural_checks():
    y=s.symbols('y',real=True);z,E,alpha=s.symbols('z E alpha',real=True)
    R=s.Function('R')(y);u=s.Function('u')(y)
    em,eh,ek,aa,eb,ep=[s.Function(name)(y) for name in ('em','eh','ek','aa','eb','ep')]
    V=4*z+E*alpha;m=4*z+em;H=s.Rational(5,8)+eh;K=4*z*H+ek
    A=16*z*z+8*z*em+aa;b=s.Rational(5,6)+eb;p=5+ep
    rules={s.diff(R,y):R,s.diff(u,y):u/10,s.diff(em,y):E*alpha-em,
        s.diff(eh,y):-s.Rational(8,5)*eh,s.diff(ek,y):E*alpha-s.Rational(8,5)*ek,
        s.diff(aa,y):(E*alpha)**2-aa,s.diff(eb,y):-s.Rational(6,5)*eb,s.diff(ep,y):-ep/5}
    moments=[R*m,s.sqrt(2)*R**s.Rational(3,2)*u*H,s.sqrt(2)*R**s.Rational(3,2)*u*K,
        R*A-R*u*u*b/2,u*u*p/2]
    expected=[R*V,s.sqrt(2)*R**s.Rational(3,2)*u,s.sqrt(2)*R**s.Rational(3,2)*u*V,
        R*(V**2-u*u/2),u*u/2]
    if any(s.simplify(s.diff(row,y).xreplace(rules)-rhs)!=0 for row,rhs in zip(moments,expected)):
        raise ArithmeticError('Centered source equations change physical primitive RHSs')
    Q=s.Function('Q')(y)
    for k in range(5):
        expected=s.exp(y/2)*sum(s.binomial(k,j)*s.Rational(1,2)**(k-j)*s.diff(Q,y,j) for j in range(k+1))
        if s.simplify(s.diff(s.exp(y/2)*Q,y,k)-expected)!=0:raise ArithmeticError('True radial velocity prefactor derivative failed')
    # These are identities of the SAME source functions, not tests of
    # intersecting interval boxes. Transport from Rm to Rh cancels the
    # original defect normalization at Rh exactly, including the swirl row.
    length=s.symbols('length',real=True);Am=s.symbols('Am',positive=True)
    names=(em,ek,eh,aa,eb,ep)
    rates=(1,s.Rational(8,5),s.Rational(8,5),1,s.Rational(6,5),s.Rational(1,5))
    terminal=[row*s.exp(-rate*length) for row,rate in zip(names,rates)]
    defects=[terminal[0]*s.exp(length),terminal[1]*s.exp(s.Rational(8,5)*length),
        terminal[2]*s.exp(s.Rational(8,5)*length),
        terminal[3]*s.exp(length)/Am**2-terminal[4]*s.exp(s.Rational(6,5)*length)/2,
        terminal[5]*s.exp(length/5)/2]
    inlet=[em,ek,eh,aa/Am**2-eb/2,ep/2]
    if any(s.simplify(a-b)!=0 for a,b in zip(defects,inlet)):
        raise ArithmeticError('Rh defect transport does not recover the exact Rm source')
    theta=s.Rational(5,8)+eh
    actual=[4*z+em,theta,4*z*theta+ek,aa/Am**2-(s.Rational(5,6)+eb)/2,(5+ep)/2]
    patch=[4*z+defects[0],s.Rational(5,8)+defects[2],
        4*z*(s.Rational(5,8)+defects[2])+defects[1],defects[3]-s.Rational(5,12),defects[4]+s.Rational(5,2)]
    if any(s.simplify(a-b)!=0 for a,b in zip(actual,patch)):
        raise ArithmeticError('Actual five physical Rm inlet histories changed')
    return dict(exact_physical_primitive_y_RHS_identities=5,exact_radial_velocity_prefactor_identities=5,
        exact_Rh_defect_to_Rm_transport_identities=5,exact_actual_Rm_primitive_inlet_identities=5,
        Rm_mixed_join_proof='identical source histories, original P0 and swirl; open zero-bump neighborhood gives identical primitive RHSs, hence all mixed derivatives through4',passed=True)


def independent_physical_fixture():
    with mp.workdps(85):
        c=MPIntervalContext();c.dps=110;tol=mp.mpf('1e-65');y=mp.mpf('.3');z=mp.mpf('.2')
        Pstar=mp.mpf('2.4');delta=mp.mpf('.0005');coeff=[mp.mpf(v) for v in (1,-1,2,-1)]+[mp.mpf(1)/3,-mp.mpf(1)/5]
        squared=[sum((coeff[i]*coeff[k-i] for i in range(len(coeff)) if 0<=k-i<len(coeff)),mp.mpf(0)) for k in range(11)]
        alpha=lambda yy:sum(a*yy**i for i,a in enumerate(coeff))
        def integrate(k,polynomial,yy):
            # Exact independent antiderivatives of exp(k*s)*polynomials.
            return sum(a*(mp.exp(k*yy)*sum((-1)**j*mp.factorial(i)/mp.factorial(i-j)*yy**(i-j)/k**(j+1)
                for j in range(i+1))-(-1)**i*mp.factorial(i)/k**(i+1)) for i,a in enumerate(polynomial))
        E=lambda zz:mp.mpf('.01')+zz/500+zz**3/3000
        initial=lambda zz:[mp.mpf(i+1)/100+zz/1000+zz**2/2000+zz**5/10000 for i in range(6)]
        u0=lambda zz:Pstar*mp.exp(mp.mpf('-.8'))/(1+zz*zz)
        datum=lambda zz:mp.mpf('.03')+zz/100+zz**4/1000
        def physical(yy,zz):
            R=mp.exp(yy);u=u0(zz)*mp.exp(yy/10);ee=E(zz);v=4*zz+ee*alpha(yy)
            em,eh,ek,aa,eb,ep=initial(zz)
            I1=integrate(1,coeff,yy);I16=integrate(mp.mpf('1.6'),coeff,yy);I2=integrate(1,squared,yy)
            M=4*zz*R+em+ee*I1
            T=mp.sqrt(2)*u0(zz)*(mp.mpf('.625')+eh+(mp.exp(mp.mpf('1.6')*yy)-1)/mp.mpf('1.6'))
            J=4*zz*T+mp.sqrt(2)*u0(zz)*(ek+ee*I16)
            A=16*zz*zz*R+8*zz*(em+ee*I1)+aa+ee*ee*I2
            B=u0(zz)**2*(mp.mpf(5)/6+eb+(mp.exp(mp.mpf('1.2')*yy)-1)/mp.mpf('1.2'))
            Mp=u0(zz)**2*(5+ep+(mp.exp(yy/5)-1)*5)/2
            m=M/R
            mZ=mp.diff(lambda q:4*q+(initial(q)[0]+E(q)*I1)/R,zz)
            Q=(2*zz*v-(1-delta)*zz*m-(1-zz*zz)*mZ)/(1-delta*zz*zz)
            return dict(Utheta=u,Uz=v,Ur=mp.sqrt(R/2)*Q,P_over_Pstar2=datum(zz)+Mp/Pstar**2,
                Mtheta=T,Mtheta_z=J,Mz=M,Mztheta_over_Pstar2=(A-B/2)/Pstar**2,Mp_over_Pstar2=Mp/Pstar**2)
        def jet(fn):return IntervalTaylor(c,[c.mpf([v-tol,v+tol])/math.factorial(n)
            for n in range(6) for v in (mp.diff(fn,z,n),)])
        I1=integrate(1,coeff,y);I16=integrate(mp.mpf('1.6'),coeff,y);I2=integrate(1,squared,y)
        centered=dict(mean_error=jet(lambda zz:(initial(zz)[0]+E(zz)*I1)*mp.exp(-y)),
            angular_error=jet(lambda zz:initial(zz)[1]*mp.exp(mp.mpf('-1.6')*y)),
            mixed_error=jet(lambda zz:(initial(zz)[2]+E(zz)*I16)*mp.exp(mp.mpf('-1.6')*y)),
            axial_square=jet(lambda zz:(initial(zz)[3]+E(zz)**2*I2)*mp.exp(-y)),
            swirl_error=jet(lambda zz:initial(zz)[4]*mp.exp(mp.mpf('-1.2')*y)),
            pressure_error=jet(lambda zz:initial(zz)[5]*mp.exp(-y/5)))
        logu=jet(lambda zz:mp.mpf('-.8')+y/10-mp.log(1+zz*zz))
        derivatives=[c.mpf(mp.diff(alpha,y,k)) for k in range(5)]
        packet=source_mixed(c,c.mpf(z),c.mpf(delta),jet(E),derivatives,centered,logu,jet(datum),c.mpf(1/Pstar**2))
        R0=mp.exp(y);baseu=u0(z)*mp.exp(y/10)
        mapping={'Utheta_over_current_Utheta':('Utheta',baseu),'Uz':('Uz',1),
            'Ur_over_current_sqrt_R_over_2':('Ur',mp.sqrt(R0/2)), 'P_over_Pstar2':('P_over_Pstar2',1),
            'Mtheta_over_current_sqrt2_R_1p5_Utheta':('Mtheta',mp.sqrt(2)*R0**mp.mpf('1.5')*baseu),
            'Mtheta_z_over_current_sqrt2_R_1p5_Utheta':('Mtheta_z',mp.sqrt(2)*R0**mp.mpf('1.5')*baseu),
            'Mz_over_current_R':('Mz',R0),'Mztheta_over_current_R_Pstar2':('Mztheta_over_Pstar2',R0),'Mp_over_Pstar2':('Mp_over_Pstar2',1)}
        count=0
        for group in ('physical_velocity_pressure_y_Z_mixed4','physical_five_primitive_y_Z_mixed4'):
            for name,grid in packet[group].items():
                label,scale=mapping[name]
                for key,value in grid.items():
                    k,n=[int(v[1:]) for v in key.split('_')]
                    expected=mp.diff(lambda yy,zz:physical(yy,zz)[label],(y,z),(k,n))/scale
                    lo,hi=endpoints(value)
                    if not lo-tol*10000<=expected<=hi+tol*10000:raise ArithmeticError('Independent physical restore derivative failed: '+name+' '+key)
                    count+=1
        return dict(independent_physical_closed_primitive_mixed_derivatives=count,
            nonconstant_cutoff_derivatives_and_nonzero_histories=True,finite_fixture_only=True,actual_source_admission=False,passed=True)


def run():
    with mp.workdps(280):
        raw=json.loads((HERE/NAME).read_bytes());provider=CompliantReferenceRestoreMixedC4();c=provider.ctx
        for name,digest in raw['input_hashes'].items():
            if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:raise ValueError('Actual reference mixed source changed: '+name)
        packets=[raw[name] for name in ('whole_reference','whole_restoration','whole_postrestore','actual_Rsh_exit',
            'actual_Rz_reference_side','actual_restore_inlet','actual_restore_exit','actual_postrestore_inlet','actual_Rm_exit')]+raw['interior_packets']
        counts={name:0 for name in ('physical_velocity_pressure_y_Z_mixed4','physical_five_primitive_y_Z_mixed4')}
        for packet in packets:
            for group in counts:
                for grid in packet[group].values():
                    if set(grid)!={'y'+str(k)+'_Z'+str(n) for k in range(5) for n in range(5-k)}:raise ValueError('Reference mixed grid incomplete')
                    for value in grid.values():
                        lo,hi=endpoints(read_interval(c,value))
                        if not lo<=hi or not mp.isfinite(lo) or not mp.isfinite(hi):raise ArithmeticError('Nonfinite reference mixed derivative')
                        counts[group]+=1
            if not packet['positive_exponential_cap_is_enclosure_only'] or not packet['actual_reference_moment_histories_and_P0_retained']:
                raise ValueError('Exact positive source or actual histories lost')
            if packet['Rsh_reshape_join_certified'] or packet['full_inner_interfaces_certified'] or packet['full_cartesian_vector_derivatives_certified'] or packet['temporal_recursion']:
                raise ValueError('Unbuilt full spatial/temporal scope promoted')
        for a,b in (('actual_Rz_reference_side','actual_restore_inlet'),('actual_restore_exit','actual_postrestore_inlet')):
            for group in counts:
                if raw[a][group]!=raw[b][group]:raise ArithmeticError('Exact source-bound mixed endpoint histories changed')
        for t in (0,1):
            jets=sigma_jets(c,c.mpf(t))
            if endpoints(jets[0])!=(mp.mpf(t),mp.mpf(t)) or any(endpoints(jets[k])!=(mp.mpf(0),mp.mpf(0)) for k in range(1,5)):
                raise ArithmeticError('Original flat cutoff endpoint derivatives changed')
        patch_provider=CompliantActualPatchMixedC4()
        prior=json.loads((HERE/'lei_ren_part1_paper_compliant_actual_patch_mixed_C4_check.json').read_bytes())
        if not prior['exact_coordinate_and_join_checks']['actual_Rm_join_uses_open_zero_correction_neighborhood']:
            raise ValueError('Original open zero-bump Rm source join missing')
        if patch_provider.family!=provider.family or patch_provider.source!=provider.source:
            raise ValueError('Rm join uses different actual source families')
        patch=patch_provider.evaluate(1,[-1,1])
        ref=raw['actual_Rm_exit']['physical_velocity_pressure_y_Z_mixed4']
        names={'Uz':'Uz','Ur_over_current_sqrt_R_over_2':'Ur_over_sqrt_Rm_over_2','P_over_Pstar2':'P_over_Pstar2'}
        for left,right in names.items():
            for key,value in ref[left].items():
                lo,hi=endpoints(read_interval(c,value));a,b=endpoints(patch['physical_velocity_pressure_y_Z_mixed4'][right][key])
                if max(lo,a)>min(hi,b):raise ArithmeticError('Actual Rm reference/patch mixed data inconsistent')
        if not raw['postrestore_unpatched_interval_stops_at_Rm']:raise ValueError('Unpatched source incorrectly crosses actual patch')
        result=dict(actual_five_defect_family_sha256=provider.family,implicit_source_sha256=provider.source,
            structural_checks=structural_checks(),independent_physical_fixture=independent_physical_fixture(),
            actual_mixed_bounds_checked=counts,exact_Rz_and_restore_exit_mixed_packets_retained=True,
            original_flat_cutoff_endpoint_derivatives_checked=10,
            actual_reference_restore_mixed4_available=True,Rz_restore_end_Rm_functional_joins_certified=True,
            Rm_patch_data_consistency_complements_functional_source_join=True,
            exact_positive_source_amplitudes_retained=True,postrestore_unpatched_interval_stops_at_Rm=True,
            Rsh_reshape_join_certified=False,full_inner_interfaces_certified=False,
            full_cartesian_vector_derivatives_certified=False,admissible_stress_lift_constructed=False,
            temporal_recursion=False,all_passed=True,
            input_hashes={**raw['input_hashes'],NAME:hashlib.sha256((HERE/NAME).read_bytes()).hexdigest(),
                Path(__file__).name:hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('Actual reference/restore mixed4, exact flat interfaces and original source units PASS',flush=True)
    return result


if __name__=='__main__':run()
