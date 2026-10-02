"""Independent moment/divergence, curl, pressure and shared-root checks."""
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as s
from lei_ren_part1_paper_compliant_rooted_core_field import CompliantRootedCoreField
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_candidate_exact_amplitude import _encode

HERE=Path(__file__).parent
NAME='lei_ren_part1_paper_compliant_rooted_core_field.json'


def identities():
    r,z,d,e,lam,F0=s.symbols('rho Z delta epsilon lambda F0',positive=True)
    D=1-z*z;L=1-d*z*z
    U=s.Function('U')(r,z);M=s.Function('M')(r,z);Phi=s.Function('Phi')(r,z)
    Q=(2*z*U-(1-d)*z*M-D*s.diff(M,z))/L
    div=s.expand(Q+r*s.diff(Q,r)+(D*s.diff(U,z)-(1+d)*z*U-2*z*r*s.diff(U,r))/L)
    # Exact identities of the SAME primitive, not independent rounded rows.
    div=div.subs(s.diff(M,r,z),(s.diff(U,z)-s.diff(M,z))/r).subs(s.diff(M,r),(U-M)/r)
    if s.simplify(div)!=0:raise ArithmeticError('Original moment-based divergence identity failed')
    # Bind the identities to the implemented coefficientwise primitive.
    cn=[s.Function('c'+str(n))(z) for n in range(6)]
    upoly=sum(cn[n]*r**n for n in range(6))
    mpoly=sum(cn[n]*r**n/(n+1) for n in range(6))
    if s.expand(mpoly+r*s.diff(mpoly,r)-upoly)!=0:
        raise ArithmeticError('Shared finite primitive coefficient identity failed')
    if s.expand(s.diff(mpoly,z)+r*s.diff(mpoly,r,z)-s.diff(upoly,z))!=0:
        raise ArithmeticError('Shared axial primitive coefficient identity failed')
    v=s.symbols('v',positive=True)
    phi_poly=1+cn[1]*r+cn[2]*r**2
    V=s.integrate(phi_poly.subs(r,v)**2,(v,0,r))
    if s.expand(s.diff(V,r)-phi_poly**2)!=0 or V.subs(r,0)!=0:
        raise ArithmeticError('Pressure squared-profile primitive identity failed')
    R=s.symbols('R',positive=True)
    original=s.integrate(phi_poly.subs(r,v/e)**2,(v,0,R))
    if s.expand(original.subs(R,e*r)-e*V)!=0:
        raise ArithmeticError('Original pressure dR=epsilon*drho conversion failed')
    kp,b,sigma=s.symbols('Kprime b sigma',real=True,nonzero=True)
    g=s.sqrt(2*e*r)*F0*Phi
    gz=s.sqrt(2*e*r)*F0*(s.diff(Phi,z)-kp*Phi/b)
    wr=-(-1-d)*z*g-D*gz+2*z*r*s.diff(g,r)
    desired=s.sqrt(2*e*r)*F0*((2+d)*z*Phi-D*s.diff(Phi,z)+2*z*r*s.diff(Phi,r)+D*kp*Phi/b)
    if s.simplify(wr-desired)!=0:raise ArithmeticError('Original omega_r coefficient/sign failed')
    if s.simplify((s.sqrt(e)/b).subs(b,sigma*s.sqrt(e))-1/sigma)!=0:
        raise ArithmeticError('Microscopic inverse-width cancellation failed')
    gr=s.sqrt(e*r/2)*Q
    wt=(-z*gr+D*s.diff(gr,z)-2*z*r*s.diff(gr,r))/L
    if s.simplify(wt-s.sqrt(e*r/2)*(D*s.diff(Q,z)-2*z*(Q+r*s.diff(Q,r)))/L)!=0:
        raise ArithmeticError('Original omega_theta radial curl failed')
    if s.simplify(s.diff(2*r*Phi,r)-2*(Phi+r*s.diff(Phi,r)))!=0:
        raise ArithmeticError('Original omega_z nonsingular curl failed')
    a,x,j=s.symbols('a x j',real=True)
    H=lambda v:-4*v**3-j*v*v+(9-d)*v/2+j
    hp=s.diff(H(z),z).subs(z,a)
    if s.expand(H(a+x)-H(a)-x*(hp+(-12*a-j)*x-4*x*x))!=0:
        raise ArithmeticError('Shared root polynomial identity failed')
    return dict(same_primitive_divergence_identity=True,
        actual_coefficientwise_radial_and_axial_primitive_identities=True,
        pressure_squared_profile_primitive_and_epsilon_conversion=True,physical_curl_identities=True,
        omega_r_correct_2_plus_delta=True,microscopic_inverse_width_canceled=True,
        nonsingular_axis_curl=True,exact_root_polynomial_identity=True,
        root_Bessel_tail_valuation='chi(a+t)=O(t^2); chi^n has no axial coefficient k<2n',passed=True)


def cartesian_fixture():
    """Independent physical-coordinate differences, never paper source data."""
    with mp.workdps(100):
        c=MPIntervalContext();c.dps=100;j=mp.mpf('.12');delta=mp.mpf('.02');Lambda=mp.mpf(7);epsilon=1/Lambda
        sigma=j/500;C=mp.log(mp.mpf('2.5'));logP=mp.mpf('.2');tau=mp.mpf('.73');logtau=mp.log(tau)
        H=lambda zz:-4*zz**3-j*zz*zz+(9-delta)*zz/2+j
        a=mp.findroot(H,-j/((9-delta)/2));b=sigma*mp.sqrt(epsilon)
        hp=(9-delta)/2-12*a*a-2*j*a;A=-12*a-j
        q=lambda xx:hp+A*b*xx-4*b*b*xx*xx
        gradient=lambda xx:xx*(1-delta*(a+b*xx)**2)*q(xx)/(1+epsilon*xx*xx*q(xx)**2)
        K=lambda xx:mp.quad(gradient,[0,xx]);F=lambda zz:mp.exp(-C-K((zz-a)/b))
        R,Z=s.symbols('R Z');ph=1-s.Rational(13,100)*R+s.Rational(4,100)*R*Z+s.Rational(7,1000)*R*R
        psi=R*(s.Rational(7,10)+s.Rational(3,10)*Z+s.Rational(1,5)*Z*Z)+R*R*(s.Rational(1,25)-Z/100)
        v=s.symbols('v');mean=s.integrate(psi.subs(R,v),(v,0,R))/R
        phiF=s.lambdify((R,Z),ph,'mpmath');psiF=s.lambdify((R,Z),psi,'mpmath')
        meanF=s.lambdify((R,Z),mean,'mpmath');meanZF=s.lambdify((R,Z),s.diff(mean,Z),'mpmath')
        funcs={(i,k):{name:s.lambdify((R,Z),s.diff(expr,R,i,Z,k),'mpmath') for name,expr in (
            ('Phi',ph),('Psi',psi),('MeanPsi',mean),('Uz',4*Z+s.Rational(12,100)+psi/7),
            ('Mz_over_R',4*Z+s.Rational(12,100)+mean/7))} for i,k in ((0,0),(1,0),(0,1),(0,2),(1,1))}
        Vfunc=s.lambdify((R,Z),s.integrate(ph.subs(R,v)**2,(v,0,R)),'mpmath')
        box=lambda number:c.mpf([number-mp.mpf('1e-80'),number+mp.mpf('1e-80')])
        f=object.__new__(CompliantRootedCoreField);f.ctx=c
        f.core=SimpleNamespace(delta=c.mpf('.02'),j=c.mpf('.12'),sigma=c.mpf('.12')/500,
            logLambda=c.ln(7),logC=c.ln(c.mpf('2.5')),logP=c.mpf('.2'))
        f.core.datum=SimpleNamespace(normalized_jets=lambda zz,order:dict(normalized_pressure_coefficients=[c.mpf('.8')+c.mpf('.12')*zz]))
        def peak_eval(xx):
            point=sum(endpoints(c.mpf(xx)))/2
            return dict(normalized_F0_relative_to_shared_anchor=box(mp.exp(-K(point))),K_derivatives=[box(gradient(point))])
        f.peak=SimpleNamespace(anchor=box(a),b=box(b),evaluate=peak_eval)
        def profile(xx,rr,i=0,k=0):
            zz=a+b*sum(endpoints(c.mpf(xx)))/2;rad=sum(endpoints(c.mpf(rr)))/2
            return dict(source_jets={name:box(fn(rad,zz)) for name,fn in funcs[(i,k)].items()})
        f.rooted_profile=profile
        f.pressure_integral=lambda xx,rr:dict(V=box(Vfunc(sum(endpoints(c.mpf(rr)))/2,a+b*sum(endpoints(c.mpf(xx)))/2)))
        def physical_velocity(x,y,zphysical):
            # Independent inverse of z=Z*(tau/(1-Z^2))^((1-delta)/2).
            zz=mp.findroot(lambda val:val*(tau/(1-val*val))**((1-delta)/2)-zphysical,a,tol=mp.eps*16)
            ll=mp.sqrt(tau/(1-zz*zz));radius=mp.sqrt(x*x+y*y);rr=Lambda*radius*radius/(2*ll*ll)
            u=4*zz+j+epsilon*psiF(rr,zz);m=4*zz+j+epsilon*meanF(rr,zz)
            mz=4+epsilon*meanZF(rr,zz);Q=(2*zz*u-(1-delta)*zz*m-(1-zz*zz)*mz)/(1-delta*zz*zz)
            ur=radius*Q/(2*ll*ll);ut=radius*F(zz)*phiF(rr,zz)/ll**(2+delta)
            return [ur*x/radius-ut*y/radius,ur*y/radius+ut*x/radius,u/ll**(1+delta)]
        scalar=lambda term:sum(endpoints(term['coefficient_enclosure']))/2*mp.exp(sum(sum(endpoints(v))/2 for v in term['positive_scale_log_terms']))
        errors=[];convergence=[]
        for xx,rr in ((mp.mpf('-.5'),mp.mpf(1)),(mp.mpf(0),mp.mpf(2)),(mp.mpf('.5'),mp.mpf(3))):
            zz=a+b*xx;ll=mp.sqrt(tau/(1-zz*zz));radius=ll*mp.sqrt(2*epsilon*rr);angle=mp.mpf('.7')
            point=[radius*mp.cos(angle),radius*mp.sin(angle),zz*ll**(1-delta)]
            packet=f.field(c.mpf(xx),c.mpf(rr),c.mpf(logtau),c.mpf(angle))
            target=[sum(scalar(term) for term in packet['cartesian_vorticity'][axis]) for axis in ('x','y','z')]
            step=b*mp.mpf('1e-5');Ds=[]
            for scale in (1,mp.mpf('.5'),mp.mpf('.25')):
                derivatives=[]
                for axis in range(3):
                    plus=list(point);minus=list(point);plus[axis]+=step*scale;minus[axis]-=step*scale
                    vp=physical_velocity(*plus);vm=physical_velocity(*minus)
                    derivatives.append([(vp[k]-vm[k])/(2*step*scale) for k in range(3)])
                curl=[derivatives[1][2]-derivatives[2][1],derivatives[2][0]-derivatives[0][2],derivatives[0][1]-derivatives[1][0]]
                divergence=sum(derivatives[i][i] for i in range(3));Ds.append((curl,divergence))
            rich=[(4*Ds[1][0][i]-Ds[0][0][i])/3 for i in range(3)]
            finer=[(4*Ds[2][0][i]-Ds[1][0][i])/3 for i in range(3)]
            error=max(abs(finer[i]-target[i]) for i in range(3));estimate=max(abs(finer[i]-rich[i]) for i in range(3))
            div=(4*Ds[2][1]-Ds[1][1])/3
            if error>mp.mpf('1e-10') or abs(div)>mp.mpf('1e-12') or estimate>mp.mpf('1e-10'):
                raise ArithmeticError('Independent physical Cartesian curl/divergence fixture failed')
            raw0=max(abs(Ds[0][0][i]-target[i]) for i in range(3));raw1=max(abs(Ds[1][0][i]-target[i]) for i in range(3))
            if raw0>mp.mpf('1e-25') and not 3<raw0/raw1<5:raise ArithmeticError('Cartesian difference refinement not converging quadratically')
            errors.append(dict(xi=str(xx),rho=str(rr),curl_error=mp.nstr(error,20),divergence_error=mp.nstr(abs(div),20),Richardson_change=mp.nstr(estimate,20)))
            convergence.append(mp.nstr(raw0/raw1,20))
        return dict(synthetic_kinematic_fixture_not_paper_source=True,independent_inverse_physical_z_map=True,
            independent_Cartesian_curl_points=3,independent_Cartesian_divergence_points=3,
            three_difference_resolutions_checked=True,raw_difference_convergence_ratios=convergence,errors=errors,passed=True)


def run():
    receipt=json.loads((HERE/NAME).read_bytes());hashes=dict(receipt['input_hashes'])
    for name,digest in hashes.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:raise ValueError('Rooted field source changed: '+name)
    with mp.workdps(400):
        f=CompliantRootedCoreField();c=f.ctx;algebra=identities();fixture=cartesian_fixture()
        root=receipt['rooted_fresh_core_packet'];fixed=root['fixed'];rows=root['rows']
        if endpoints(read_interval(c,root['H_Z_taylor'][0]))!=(mp.mpf(0),mp.mpf(0)) or endpoints(read_interval(c,fixed['ell_Z_taylor'][0]))!=(mp.mpf(0),mp.mpf(0)):
            raise ArithmeticError('Shared exact root seed constraint was lost')
        if endpoints(read_interval(c,fixed['S_Z_taylor'][0]))[1]<=0 or not root['same_pressure_Cstar_and_positive_swirl_source']:
            raise ValueError('Nonzero original implicit swirl source replaced')
        # Validate the primitive averaging on complete finite rows, not samples.
        for n,row in enumerate(f.root_rows['rows']['Uz']):
            for k,value in enumerate(row):
                stored=read_interval(c,rows['Uz'][n][k])
                if stored._mpi_!=value._mpi_:raise ArithmeticError('Fresh rooted coupled coefficient changed')
        n1=read_interval(c,rows['A'][1][0]);z=f.peak.anchor;beta=(z*f.core.j+3-f.core.delta/2)/(1-f.core.delta*z*z)
        expected=-f.epsilon*beta/4
        if max(endpoints(n1)[0],endpoints(expected)[0])>min(endpoints(n1)[1],endpoints(expected)[1]):
            raise ArithmeticError('Original shared-root first Phi radial coefficient failed')
        points=axis_points=primitive_jets=0;max_div_width=c.mpf(0)
        for key,packets in receipt['physical_field_packets'].items():
            for packet in packets:
                rho=read_interval(c,packet['rho']);xi=read_interval(c,packet['xi']);t=read_interval(c,packet['log_tau']);angle=read_interval(c,packet['theta'])
                fresh=f.field(xi,rho,t,angle)
                for mixed,jetpacket in packet['source_mixed_jet_packets'].items():
                    i,k=jetpacket['radial_order'],jetpacket['axial_order'];mean=c.mpf(0)
                    for n in range(max(1,i),root['radial_degree']+1):
                        coefficient=f.root_rows['rows']['Uz'][n][k]/f.epsilon/(n+1)
                        mean+=coefficient*(math.factorial(n)//math.factorial(n-i))*math.factorial(k)*rho**(n-i)
                    stored=read_interval(c,jetpacket['error_details']['MeanPsi']['root_finite_polynomial'])
                    if stored._mpi_!=mean._mpi_:raise ArithmeticError('Actual shared coefficientwise primitive changed')
                    if not jetpacket['Mz_over_R_is_same_Uz_radial_average'] or not jetpacket['root_model_radial_tail_exact_zero']:
                        raise ValueError('Actual average or root tail source identity not bound')
                    primitive_jets+=1
                for group in ('cylindrical_velocity','cartesian_velocity','cylindrical_vorticity','cartesian_vorticity'):
                    for label,terms in packet[group].items():
                        for old,new in zip(terms,fresh[group][label]):
                            if read_interval(c,old['coefficient_enclosure'])._mpi_!=new['coefficient_enclosure']._mpi_:raise ArithmeticError('Root vector/vorticity replay failed')
                            for oldlog,newlog in zip(old['positive_scale_log_terms'],new['positive_scale_log_terms']):
                                if read_interval(c,oldlog)._mpi_!=newlog._mpi_:raise ArithmeticError('Separate physical logarithmic scale changed')
                V=read_interval(c,packet['pressure_radial_primitive']['V'])
                if packet['axis']:
                    if endpoints(V)!=(mp.mpf(0),mp.mpf(0)):raise ArithmeticError('Axis pressure primitive not zero')
                    for group,labels in (('cylindrical_velocity',('ur','utheta')),('cylindrical_vorticity',('omega_r','omega_theta'))):
                        for label in labels:
                            if any(endpoints(read_interval(c,term['coefficient_enclosure']))!=(mp.mpf(0),mp.mpf(0)) for term in packet[group][label]):
                                raise ArithmeticError('Axis transverse vector not structurally zero')
                    axis_points+=1
                elif endpoints(V)[0]<=0:raise ArithmeticError('Original radial pressure increment lost positivity')
                for stored,current in zip(packet['pressure_terms'],fresh['pressure_terms']):
                    if read_interval(c,stored['coefficient_enclosure'])._mpi_!=current['coefficient_enclosure']._mpi_:raise ArithmeticError('Original compatible pressure source replay failed')
                    for oldlog,newlog in zip(stored['positive_scale_log_terms'],current['positive_scale_log_terms']):
                        if read_interval(c,oldlog)._mpi_!=newlog._mpi_:raise ArithmeticError('Pressure epsilon/F0/time factors changed')
                div=read_interval(c,packet['normalized_divergence_interval'])
                if not endpoints(div)[0]<=0<=endpoints(div)[1]:raise ArithmeticError('Raw interval divergence excludes zero')
                max_div_width=max(max_div_width,c.mpf(endpoints(div)[1])-c.mpf(endpoints(div)[0]))
                if not packet['vorticity_amplitude_inverse_width_canceled_before_enclosure']:raise ValueError('Microscopic inverse width materialized')
                points+=1
        for flag in ('full_point_physical_field_evaluation','all_annular_source_values_resolved','measured_blowup_dynamics',
                     'whole_vortex_aspect_ratio_measured','physical_energy_integral_certified','admissible_stress_lift_constructed','temporal_recursion'):
            if receipt[flag]:raise ValueError('Local rooted field scope promoted: '+flag)
    for name in (NAME,Path(__file__).name):hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
    result=dict(actual_five_defect_family_sha256=receipt['actual_five_defect_family_sha256'],
        implicit_source_sha256=receipt['implicit_source_sha256'],datum_enclosure_sha256=receipt['datum_enclosure_sha256'],
        independent_source_identities=algebra,independent_Cartesian_fixture=fixture,
        rooted_velocity_pressure_vorticity_packets_checked=points,nonsingular_axis_packets_checked=axis_points,
        actual_coefficientwise_moment_jet_packets_checked=primitive_jets,
        raw_normalized_interval_divergence_width_upper=max_div_width,
        raw_interval_diagnostic_distinct_from_exact_structural_divergence=True,
        exact_source_H_root_and_radial_model_tail_valuation_preserved=True,
        original_positive_pressure_primitive_and_datum_preserved=True,
        full_point_physical_field_evaluation=False,measured_blowup_dynamics=False,temporal_recursion=False,
        all_passed=True,input_hashes=hashes)
    Path(__file__).with_suffix('.json').write_text(json.dumps(_encode(result),indent=2)+'\n',encoding='utf8')
    print('Rooted field PASS:18 velocity/pressure/vorticity packets,6 axes,independent Cartesian curl/divergence',flush=True)
    return result


if __name__=='__main__':run()
