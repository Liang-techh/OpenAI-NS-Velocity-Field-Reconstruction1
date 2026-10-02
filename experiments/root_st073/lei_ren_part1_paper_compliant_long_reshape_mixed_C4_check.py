"""Independent actual long-reshape physical mixed4 and interface checks."""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
import sympy as s
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_long_reshape_mixed_C4 import (
    CompliantLongReshapeMixedC4,reshape_mixed,exponential_derivatives,IntervalTaylor)
from lei_ren_part1_paper_compliant_reference_restore_mixed_C4 import CompliantReferenceRestoreMixedC4
from lei_ren_part1_paper_compliant_reference_restore_profiles_check import source_binding
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import sigma_jets
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
NAME='lei_ren_part1_paper_compliant_long_reshape_mixed_C4.json'


def symbolic_checks():
    y=s.symbols('y',real=True);L=s.Function('L')(y);B,T=s.symbols('B T',nonzero=True)
    q=[s.diff(L,y,k) for k in range(1,5)];rows=exponential_derivatives(q)
    for k,row in enumerate(rows):
        if s.simplify(s.diff(s.exp(L),y,k)/s.exp(L)-row)!=0:
            raise ArithmeticError('Ordinary exponential Bell derivative changed')
    m,eh,ek,A,b,p=[s.symbols(n) for n in ('m','eh','ek','A','b','p')]
    z,H=s.symbols('z H')
    # The centered provider is an algebraic coordinate change of the same
    # actual Rsh primitive histories; the raw V source is shared.
    centered=[m-4*z,H-s.Rational(5,8),ek-4*z*H,A-8*z*m+16*z*z,b-s.Rational(5,6),p-5]
    recovered=[4*z+centered[0],s.Rational(5,8)+centered[1],4*z*H+centered[2],
        16*z*z+8*z*centered[0]+centered[3],s.Rational(5,6)+centered[4],5+centered[5]]
    if any(s.simplify(a-bb)!=0 for a,bb in zip(recovered,(m,H,ek,A,b,p))):
        raise ArithmeticError('Rsh actual moment coordinate change failed')
    sig=s.Function('sig');ell=y/10+B*(1-sig(y/T))
    for point in (0,T):
        for k in range(1,5):
            value=s.diff(ell,y,k).subs(y,point)
            # Substitute the original flat endpoint jet in the exact
            # chain-rule row. No cap width is differentiated.
            target=s.Rational(1,10) if k==1 else 0
            flat={node:0 for node in value.atoms(s.Subs) if node.expr.has(s.Derivative)}
            if s.simplify(value.xreplace(flat)-target)!=0:raise ArithmeticError('Flat reshape endpoint derivative mismatch')
    return dict(exact_exponential_Bell_derivative_identities=5,
        exact_Rsh_actual_centered_history_identities=6,exact_flat_endpoint_log_velocity_y_identities=8,
        functional_join_proof='same actual histories, P0 and raw V; sigma flatness gives identical physical primitive RHS derivatives; R110 original power is supplied on R109..R110',passed=True)


def independent_physical_fixture():
    with mp.workdps(65):
        c=MPIntervalContext();c.dps=95;tol=mp.mpf('1e-45');y=mp.mpf('.3');z=mp.mpf('.2')
        T=mp.mpf('1.7');Pstar=mp.mpf('2.4');delta=mp.mpf('.0005')
        B=lambda zz:mp.mpf('.015')+zz/500+zz**2/300+zz**5/1000
        sigma=lambda t:t-t*t+2*t**3-t**4+t**5/3
        ell=lambda yy,zz:mp.mpf('-.8')-mp.log(1+zz*zz)+yy/10+B(zz)*(1-sigma(yy/T))
        u=lambda yy,zz:Pstar*mp.exp(ell(yy,zz))
        V=lambda zz:4*zz+mp.mpf('.01')+zz**3/2000
        initial=lambda i,zz:mp.mpf(i+1)/20+zz/1000+zz**2/500+zz**5/10000
        p0=lambda zz:mp.mpf('.03')+zz/100+zz**4/1000
        integral=lambda fn,yy,zz:mp.quad(lambda t:fn(t,zz),[0,yy])
        rhs_theta=lambda yy,zz:mp.sqrt(2)*mp.exp(mp.mpf('1.5')*yy)*u(yy,zz)
        rhs_swirl=lambda yy,zz:mp.exp(yy)*u(yy,zz)**2
        rhs_pressure=lambda yy,zz:u(yy,zz)**2/2
        def mean(yy,zz):return (initial(2,zz)+V(zz)*(mp.exp(yy)-1))/mp.exp(yy)
        def ur(yy,zz):
            m=mean(yy,zz);mZ=mp.diff(lambda q:mean(yy,q),zz)
            Q=(2*zz*V(zz)-(1-delta)*zz*m-(1-zz*zz)*mZ)/(1-delta*zz*zz)
            return mp.exp(yy/2)/mp.sqrt(2)*Q
        def jet(fn):
            return IntervalTaylor(c,[c.mpf([v-tol,v+tol])/math.factorial(n)
                for n in range(6) for v in (mp.diff(fn,z,n),)])
        logu=jet(lambda zz:ell(y,zz))
        log_y=[jet(lambda zz,k=k:mp.diff(lambda yy:ell(yy,zz),y,k)) for k in range(1,5)]
        # Differentiate under full finite integrals; do not numerically
        # differentiate quadrature output at exponentially tiny Z steps.
        # This is an independent integral fixture, not the source Bell bound.
        coefficients=lambda fn:[mp.diff(fn,z,n)/math.factorial(n) for n in range(6)]
        integrated=lambda fn:[mp.quad(lambda yy:mp.diff(lambda zz:fn(yy,zz),z,n),[0,y])/math.factorial(n) for n in range(6)]
        add=lambda a,b:[x+q for x,q in zip(a,b)]
        product=lambda a,b:[sum(a[j]*b[n-j] for j in range(n+1)) for n in range(6)]
        scale=lambda a,f:[q*f for q in a]
        itheta=integrated(rhs_theta);iswirl=integrated(rhs_swirl);ipressure=integrated(rhs_pressure)
        vp=coefficients(V);R=mp.exp(y)
        primitive_values=dict(Mtheta=add(coefficients(lambda zz:initial(0,zz)),itheta),
            Mtheta_z=add(coefficients(lambda zz:initial(1,zz)),product(vp,itheta)),
            Mz=add(coefficients(lambda zz:initial(2,zz)),scale(vp,R-1)),
            axial=add(coefficients(lambda zz:initial(3,zz)),scale(product(vp,vp),R-1)),
            swirl=add(coefficients(lambda zz:initial(4,zz)),iswirl),
            Mp=add(coefficients(lambda zz:initial(5,zz)),ipressure))
        cover=lambda values:IntervalTaylor(c,[c.mpf([q-tol,q+tol]) for q in values])
        true={name:cover(values) for name,values in primitive_values.items()};ujet=jet(lambda zz:u(y,zz))
        shapes=dict(theta=true['Mtheta']/(ujet*(c.sqrt(2)*c.mpf(R)**c.mpf('1.5'))),
            theta_z=true['Mtheta_z']/(ujet*(c.sqrt(2)*c.mpf(R)**c.mpf('1.5'))),
            mean=true['Mz']/c.mpf(R),axial=true['axial']/c.mpf(R),
            swirl=true['swirl']/(ujet*ujet*c.mpf(R)),pressure=(true['Mp']*2)/(ujet*ujet))
        packet=reshape_mixed(c,c.mpf(z),c.mpf(delta),logu,log_y,jet(V),shapes,jet(p0),c.mpf(1/Pstar**2))
        R=mp.exp(y);uu=u(y,z)
        scales={'Utheta_over_current_Utheta':uu,'Uz':1,'Ur_over_current_sqrt_R_over_2':mp.sqrt(R/2),
            'P_over_Pstar2':1,'Mtheta_over_current_sqrt2_R_1p5_Utheta':mp.sqrt(2)*R**mp.mpf('1.5')*uu,
            'Mtheta_z_over_current_sqrt2_R_1p5_Utheta':mp.sqrt(2)*R**mp.mpf('1.5')*uu,
            'Mz_over_current_R':R,'Mztheta_over_current_R_Pstar2':R*Pstar**2,'Mp_over_Pstar2':Pstar**2}
        def original(name,yy,zz):
            if name=='Utheta_over_current_Utheta':return u(yy,zz)
            if name=='Uz':return V(zz)
            if name=='Ur_over_current_sqrt_R_over_2':return ur(yy,zz)
            raise ValueError('Full primitives use independent differentiated integrals')
        point_values={
            'Mtheta_over_current_sqrt2_R_1p5_Utheta':primitive_values['Mtheta'],
            'Mtheta_z_over_current_sqrt2_R_1p5_Utheta':primitive_values['Mtheta_z'],
            'Mz_over_current_R':primitive_values['Mz'],
            'Mztheta_over_current_R_Pstar2':add(primitive_values['axial'],scale(primitive_values['swirl'],-mp.mpf('.5'))),
            'Mp_over_Pstar2':primitive_values['Mp'],
            'P_over_Pstar2':add(coefficients(p0),scale(primitive_values['Mp'],1/Pstar**2))}
        def primitive_rhs(name,yy,zz):
            if name.startswith('Mtheta_z'):return rhs_theta(yy,zz)*V(zz)
            if name.startswith('Mtheta_'):return rhs_theta(yy,zz)
            if name.startswith('Mztheta_'):return mp.exp(yy)*V(zz)**2-rhs_swirl(yy,zz)/2
            if name.startswith('Mz_'):return mp.exp(yy)*V(zz)
            return rhs_pressure(yy,zz)/(Pstar**2 if name=='P_over_Pstar2' else 1)
        count=0
        for group in ('physical_velocity_pressure_y_Z_mixed4','physical_five_primitive_y_Z_mixed4'):
            for name,grid in packet[group].items():
                for key,value in grid.items():
                    k,n=[int(v[1:]) for v in key.split('_')]
                    if k and (group=='physical_five_primitive_y_Z_mixed4' or name=='P_over_Pstar2'):
                        expected=mp.diff(lambda yy,zz:primitive_rhs(name,yy,zz),(y,z),(k-1,n))
                    elif name in point_values:expected=point_values[name][n]*math.factorial(n)
                    else:expected=mp.diff(lambda yy,zz:original(name,yy,zz),(y,z),(k,n))
                    expected/=scales[name];lo,hi=endpoints(value)
                    if not lo-tol*10000<=expected<=hi+tol*10000:
                        raise ArithmeticError('Independent physical reshape derivative failed: '+name+' '+key)
                    count+=1
        return dict(independent_full_physical_primitive_mixed_derivatives=count,
            finite_full_integrals=True,nonconstant_log_velocity_and_nonzero_histories=True,
            finite_fixture_only=True,actual_source_admission=False,passed=True)


def run():
    with mp.workdps(280):
        raw=json.loads((HERE/NAME).read_bytes());provider=CompliantLongReshapeMixedC4();c=provider.ctx
        for path,digest in raw['input_hashes'].items():
            if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:raise ValueError('Actual reshape mixed source changed: '+path)
        packets=[raw[key] for key in ('whole_reshape','actual_R110_inlet','actual_Rsh_exit','whole_local_R109_R110_power','actual_R110_power_side')]+raw['interior_packets']
        counts={name:0 for name in ('physical_velocity_pressure_y_Z_mixed4','physical_five_primitive_y_Z_mixed4')}
        for packet in packets:
            for group in counts:
                for grid in packet[group].values():
                    if set(grid)!={'y'+str(k)+'_Z'+str(n) for k in range(5) for n in range(5-k)}:raise ValueError('Actual reshape mixed grid incomplete')
                    for value in grid.values():
                        lo,hi=endpoints(read_interval(c,value))
                        if not lo<=hi or not mp.isfinite(lo) or not mp.isfinite(hi):raise ArithmeticError('Nonfinite reshape derivative')
                        counts[group]+=1
            if not packet['positive_source_cap_is_enclosure_only'] or not packet['large_derivative_factors_combined_before_positive_source_cap']:
                raise ValueError('Positive source was replaced or capped before large factors')
            if any(packet[name] for name in ('microscopic_switch_mixed4_certified','full_inner_interfaces_certified','full_cartesian_vector_derivatives_certified','admissible_stress_lift_constructed','temporal_recursion')):
                raise ValueError('Unbuilt full/microscopic scope promoted')
        for endpoint in (0,1):
            jets=sigma_jets(c,c.mpf(endpoint))
            if endpoints(jets[0])!=(mp.mpf(endpoint),mp.mpf(endpoint)) or any(endpoints(jets[k])!=(mp.mpf(0),mp.mpf(0)) for k in range(1,5)):
                raise ArithmeticError('Original flat reshape endpoint source changed')
        reference=CompliantReferenceRestoreMixedC4()
        if reference.family!=provider.family or reference.source!=provider.source:raise ValueError('Rsh actual family/source differs')
        # Accepted reference inputs explicitly inherit this actual Rsh parent.
        if canonical_source(encode_parent(raw['actual_Rsh_exit']['actual_inherited_axial5_packet']))!=canonical_source(encode_parent(reference.reference.inputs([-1,1])['parent'])):
            raise ValueError('Reference has not inherited the SAME actual Rsh histories')
        source_proof=actual_axial_binding(provider,reference,raw)
        power=raw['actual_R110_power_side']['actual_inherited_axial5_packet']
        original=provider.reshape.inputs([-1,1])['inlet']
        for name in ('F_actual_over_F0_axial5_coefficients','Uz_actual_axial5_coefficients',
                     'actual_moment_shape_axial5_coefficients','pressure_axis_axial5_coefficients'):
            if canonical_source(power[name])!=canonical_source(encode_parent(original[name])):
                raise ValueError('R110 power/reshape raw source history is not identical')
        for group in counts:
            if canonical_source(raw['actual_R110_inlet'][group])!=canonical_source(raw['actual_R110_power_side'][group]):
                raise ValueError('R110 two-sided physical mixed grids do not have exact directed endpoints')
        capcount=0
        for proof in raw['factored_positive_source_cap_proofs']:
            if endpoints(read_interval(c,proof['log_magnitude_upper']))[1]>endpoints(read_interval(c,proof['log_cap']))[0] or not proof['exact_source_not_replaced']:
                raise ArithmeticError('Actual positive source cap before large derivative factors is invalid')
            capcount+=1
        from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
        from lei_ren_part1_paper_compliant_five_moment_repair import pack
        right=encode(pack(reference.reference_branch([-1,1],0)))
        for left,right_packet in ((raw['actual_R110_inlet'],raw['actual_R110_power_side']),(raw['actual_Rsh_exit'],right)):
            for group in counts:
                for name,grid in left[group].items():
                    for key,value in grid.items():
                        lo,hi=endpoints(read_interval(c,value));a,b=endpoints(read_interval(c,right_packet[group][name][key]))
                        if max(lo,a)>min(hi,b):raise ArithmeticError('Actual functional join enclosure inconsistent: '+name+' '+key)
        if endpoints(provider.reshape.T)!=endpoints(400*provider.reshape.A):raise ValueError('Original selected T changed')
        result=dict(actual_five_defect_family_sha256=provider.family,implicit_source_sha256=provider.source,
            actual_mixed_bounds_checked=counts,symbolic_checks=symbolic_checks(),independent_physical_fixture=independent_physical_fixture(),
            exact_original_cutoff_endpoint_checks=10,same_actual_Rsh_parent_retained=True,
            exact_shared_axial_source_binding=source_proof,actual_factored_positive_source_caps_checked=capcount,
            same_actual_R110_power_parent_retained=True,join_overlap_is_diagnostic_only=True,
            exact_R110_two_sided_physical_mixed_grids_retained=True,
            original_selected_T_and_inverse_T_factors_retained=True,
            actual_long_reshape_all_mixed_derivatives_total_order_le4_available=True,
            Rsh_reference_mixed4_join_certified=True,R110_postswitch_power_reshape_local_mixed4_join_certified=True,
            microscopic_switch_mixed4_certified=False,full_inner_interfaces_certified=False,
            full_cartesian_vector_derivatives_certified=False,admissible_stress_lift_constructed=False,temporal_recursion=False,all_passed=True,
            input_hashes={**raw['input_hashes'],NAME:hashlib.sha256((HERE/NAME).read_bytes()).hexdigest(),Path(__file__).name:hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('Actual long reshape mixed4, full physical integrals and local functional joins PASS',flush=True)
    return result


def encode_parent(packet):
    from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
    from lei_ren_part1_paper_compliant_five_moment_repair import pack
    return encode(pack(packet))


def canonical_source(value):
    # Exact directed endpoints are authoritative; width/display strings may
    # have been printed outside the producer's high precision context.
    if isinstance(value,dict):
        if 'lower_exact_mpf_tuple' in value:return (value['lower_exact_mpf_tuple'],value['upper_exact_mpf_tuple'])
        return {name:canonical_source(row) for name,row in value.items()}
    if isinstance(value,list):return [canonical_source(row) for row in value]
    return value


def actual_axial_binding(provider,reference,raw):
    graph=raw['shared_exact_axial_source']
    if graph!=provider.shared_axial_source or not graph['caps_are_not_source_integrals']:
        raise ValueError('Actual signed integral source graph changed')
    terms=graph['formal_signed_integrals']
    actual_bridge=provider.reshape.switch.bridge.actual([-1,1],provider.reshape.switch.bridge.r/100)
    bridge_sources=actual_bridge['source_integral_definitions']
    switch_sources=provider.reshape.switch.phase([-1,1],1)
    # Canonicalize the exact ORIGINAL source definitions independently of
    # the producer graph. Each sign, weight, quotient, radius, bound and
    # physical drive scale must match, not merely its source-file hash.
    if bridge_sources['V']!='V=v-integral_0^y chi*(phi_actual/phi_bar)*(R*hydro+R*Pstar^2*pressure+R^2*F0^2*swirl) dt' or bridge_sources['chi']!='1-(1-hb)*sigma(y/hb); hb=cstar*K^-100':
        raise ValueError('Original bridge source cannot be canonicalized')
    if switch_sources['exact_Uz_source']!='V100 - hb^2*integral_0^min(phase,1)(1-sigma(t))*(phi_actual/barphi)*drive(100exp(hb*t),Z)dt':
        raise ValueError('Original first-switch source cannot be canonicalized')
    drive=['R*hydro','R*Pstar^2*pressure','R^2*F0^2*swirl']
    expected=dict(I_bridge=dict(variable='s=log(R/Ra)',bounds=['0','log(100/Ra)'],sign=-1,
        chi='1-(1-hb)*sigma(s/hb)',quotient='phi_actual(R,Z)/phi_bar(R,Z)',
        radius='R=Ra*exp(s)',drive_terms=drive,exact_original_V_source=bridge_sources['V'],
        exact_original_chi_source=bridge_sources['chi']),
        I_first_switch=dict(variable='t',bounds=['0','1'],sign=-1,factor='hb^2',weight='1-sigma(t)',
        quotient='phi_actual(R,Z)/phi_bar(R,Z)',radius='R=100*exp(hb*t)',
        drive='same current-radius axial direction sum; not the angular Dbar',drive_terms=drive,
        exact_original_V_source=switch_sources['exact_Uz_source']))
    if terms!=expected:raise ValueError('Signed formal source AST differs from original source integrals')
    if graph['axial_Taylor_orders']!=list(range(6)) or graph['ordinary_derivative_factorials']!=[math.factorial(k) for k in range(6)]:
        raise ValueError('Ordinary axial source derivative units changed')
    if terms['I_bridge']['sign']!=-1 or terms['I_first_switch']['sign']!=-1:
        raise ValueError('Actual axial signed integral source lost its minus sign')
    for name in ('I_bridge','I_first_switch'):
        if terms[name]['drive_terms']!=['R*hydro','R*Pstar^2*pressure','R^2*F0^2*swirl']:
            raise ValueError('Actual axial hydro/pressure/swirl source scale omitted')
    if terms['I_first_switch']['weight']!='1-sigma(t)' or terms['I_first_switch']['factor']!='hb^2':
        raise ValueError('Original first-switch signed source changed')
    for name,digest in graph['source_bindings'].items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:raise ValueError('Actual formal integral source changed')
    # This independently binds the numerical covers to the original core,
    # actual bridge and first-switch functions. The following exact source
    # graph retains the actual signed integrals themselves, not their caps.
    receipt=source_binding(reference.reference,reference.reference.inputs([-1,1]))
    if not receipt['exact_actual_core_bridge_switch_Rsh_axial_source_chain']:raise ValueError('Actual axial source chain missing')
    Z=s.symbols('Z',real=True);j,eps=s.symbols('j eps',real=True)
    psi=s.Function('Psi')(Z);Ibridge=s.Function('I_bridge')(Z);Iswitch=s.Function('I_first_switch')(Z)
    def expression(row):
        if 'formal_integral' in row:
            if row['shared_source']!=graph['shared_source_namespace']:raise ValueError('Independent copy of an actual source integral')
            return Ibridge if row['formal_integral']=='I_bridge' else Iswitch
        if row.get('source')=='same selected j':return j
        if row['op']=='sum':return sum(expression(arg) for arg in row['args'])
        return 4*Z if row['args']==[4,'Z'] else eps*psi
    v100=expression(graph['V100']);v110=expression(graph['V110']);E=expression(graph['E_V110_minus_4Z'])
    count=0
    for k in range(6):
        for left,right in ((v100,4*Z+j+eps*psi+Ibridge),(v110,v100+Iswitch),(v110,4*Z+E)):
            if s.simplify(s.diff(left-right,Z,k)/math.factorial(k))!=0:raise ArithmeticError('Actual shared signed-integral axial source identity failed')
            count+=1
    switch=provider.reshape.switch.phase([-1,1],2)
    if not switch['Uz_constant_from_phase1_onward']:raise ValueError('Second switch introduced a new axial source')
    return dict(exact_common_signed_integral_Taylor_identities=count,axial_orders=list(range(6)),
        canonical_signed_integral_ASTs_match_original_sources=True,
        original_core_bridge_first_switch_and_sigma_source_hashes_bound=True,
        same_actual_v_cover_and_Rsh_history_bound=True,caps_are_enclosures_only=True,
        formal_integrals_numerically_reconstructed=False,passed=True)


if __name__=='__main__':run()
