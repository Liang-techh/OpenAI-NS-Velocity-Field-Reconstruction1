"""Common physical Cartesian source bounds, core/axis through heat exterior.

Every chart's already differentiated physical profile is divided by a FIXED
basepoint unit. The coordinate operators consume those derivatives directly;
normalization factors are never differentiated again. Microscopic rows use
their uncapped source terms and exact hb^-k conversion before log bounds.
The result is a signed factored enclosure, not a point coefficient choice.
"""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
from lei_ren_part1_paper_compliant_source_dispatcher import CompliantSourceDispatcher,ROUTES
from lei_ren_part1_paper_compliant_cartesian_field import (
    cartesian_templates,angular_polynomial,CS,SN,INDICES)
from lei_ren_part1_paper_compliant_pulse_physical_bounds import (
    physical_operators,interval_expression,UZ,UT,UR,P)
from lei_ren_part1_paper_compliant_inner_bridge_profiles import square,IntervalTaylor
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'
COMPONENTS=('ux','uy','uz','p')
MICRO=('bridge_first','bridge_second','switch_first','switch_second')
COMMON=('bridge_macro','switch_power','reshape','inner_reference','axial_restore','restore_buffer')
PREPULSE=('Rh_reference','O2_slope','O2_axial','O2_buffer','O3_slope_mu','O3_power')
PULSE=('pulse_entrance','pulse_main','pulse_exit','pulse_gap','pulse_gap_end','pulse_end')
POST=('flatten','outer_power','outer_angular','steep_entry','steep_power','steep_exit','waiting','heat_collar','heat_exterior')
BASES=('loghb','logPstar','logF0base','logphi','log_velocity_unit','logR','log2')


def zero_powers():return (mp.mpf(0),)*len(BASES)


def shift(powers,changes):
    result=list(powers)
    for index,value in changes.items():result[index]+=mp.mpf(value)
    return tuple(result)


def log_row(c,terms,logs,lambda_exponent,loglambda_bound,uncapped_source_rows=True):
    """Triangle bound AFTER every exact source factor/coordinate operator.

    Terms with identical factor exponents share one coefficient. Signed
    coefficients are retained. No giant positive/negative exponential is
    materialized and no numerical cap is substituted for a source scale.
    """
    merged={}
    for powers,coefficient in terms:
        merged[powers]=merged.get(powers,c.mpf(0))+coefficient
    rows=[];upper=[]
    for powers,coefficient in sorted(merged.items()):
        magnitude=max(abs(value) for value in endpoints(coefficient))
        if magnitude==0:continue
        combined=sum((logs[index]*power for index,power in enumerate(powers) if power),c.mpf(0))
        logupper=combined+lambda_exponent*loglambda_bound+c.ln(c.mpf(magnitude))
        upper.append(endpoints(logupper)[1])
        rows.append(dict(source_log_exponents=list(powers),signed_coefficient=coefficient,
            log_absolute_upper=logupper))
    bound=c.mpf(max(upper))+c.ln(len(upper)) if upper else None
    return dict(terms=rows,physical_lambda_exponent=lambda_exponent,
        log_absolute_upper=bound,exact_zero=not upper,
        original_source_factors_combined_before_enclosure=uncapped_source_rows,
        global_source_factors_combined_before_final_physical_bound=True,
        source_row_mode='uncapped_factored_rows' if uncapped_source_rows else 'provider_prebounded_mixed_rows',
        positive_source_exponentials_not_materialized=True)


def ordinary_terms(grid):
    return {(int(key.split('_')[0][1:]),int(key.split('_')[1][1:])):[(zero_powers(),value)] for key,value in grid.items()}


def micro_terms(packet,component,phase_to_y=True):
    """Expand the TRUE uncapped phase row; cancel hb before physical map.

    Original four bases are loghb,2logPstar,2logF0,2log(u/Pstar).
    Substitute 2log(u/Pstar)=log(2R)+2logF0+2logphi-2logPstar
    algebraically. Do not subtract those enormous numeric log boxes.
    """
    rows={}
    prefix=component+'/'
    for row in packet['final_factored_physical_row_ledgers']:
        if not row['physical_row'].startswith(prefix):continue
        key=row['physical_row'][len(prefix):];k,n=(int(value[1:]) for value in key.split('_'))
        terms=[]
        for term in row['terms']:
            h,p,f,u=term['source_exponents']
            powers=(mp.mpf(h-k if phase_to_y else h),mp.mpf(2*p-2*u),mp.mpf(2*f+2*u),mp.mpf(2*u),mp.mpf(0),mp.mpf(u),mp.mpf(u))
            terms.append((powers,term['final_ordinary_coefficient']))
        rows[k,n]=terms
    if len(rows)!=15:raise ValueError('Full uncapped phase mixed4 source ledger required: '+component)
    return rows


def cartesian_source_row(c,grids,component,i,j,b,Z,delta,cosine,sine,amplitudes):
    """Exact linear pullback of source terms including the moving basis."""
    N=i+j
    beta={UR:c.mpf(-1),UT:-1-delta,UZ:-1-delta,P:-2-2*delta}
    # Contributions can have different lambda exponents (radial vs swirl).
    by_label={}
    for (label,a,q),angular in cartesian_templates()[(component,i,j,b)].items():
        angular_factor=angular_polynomial(c,angular,cosine,sine)*c.sqrt(2)**(a-q)
        for (k,n),expression in physical_operators()[a,b].items():
            coefficient=angular_factor*interval_expression(c,expression,Z,delta,beta[label])
            powers_to_add=dict(amplitudes[label]);powers_to_add[5]=powers_to_add.get(5,0)-mp.mpf(N)/2
            for powers,source_coefficient in grids[label][k,n]:
                by_label.setdefault(label,[]).append((shift(powers,powers_to_add),source_coefficient*coefficient))
    return {label:(rows,beta[label]-N+b*(delta-1)) for label,rows in by_label.items()}


def time_source_row(c,grids,label,Z,delta,amplitudes):
    """Fixed physical x: lambda^-2 times original Lemma2.1 time operator."""
    beta=c.mpf(-1) if label==UR else (-2-2*delta if label==P else -1-delta)
    L=1-delta*Z**2;terms=[]
    for index,coefficient in (((0,0),-beta/(2*L)),((0,1),(1-delta)*Z/(2*L)),((1,0),1/L)):
        for powers,value in grids[label][index]:terms.append((shift(powers,amplitudes[label]),value*coefficient))
    return terms,beta-2


class CompliantGlobalPhysicalAssembly:
    def __init__(self):
        self.dispatch=CompliantSourceDispatcher();self.hashes={}
        name=PREFIX+'source_dispatcher.json';receipt=json.loads((HERE/name).read_bytes())
        if not (receipt['all_passed'] and receipt['all_routes_exercised']):raise ValueError('Accepted same-source complete chart dispatcher required')
        for path,digest in receipt['input_hashes'].items():
            if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:raise ValueError('Global physical source changed: '+path)
        self.hashes.update(receipt['input_hashes']);self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        self.family=receipt['actual_five_defect_family_sha256'];self.source=receipt['implicit_source_sha256']
        self.core=self.dispatch.provider('core');self.ctx=c=self.core.ctx;self.delta=self.core.delta
        self.pre=self.dispatch.provider('Rh_reference');self.params=self.pre.params
        self.logP=c.mpf(endpoints(self.params.logPstar));self.logC=self.core.logC
        self.logRref=c.ln(110)+10*(self.logC+self.logP)
        self.logRp=self.logRref+self.logP+1+self.params.Tw
        self.hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    def radius(self,chart,value,packet,provider):
        c=self.ctx;v=c.mpf(value)
        if chart.startswith('bridge_'):
            base=c.ln(provider.bridge.r);h=packet['width_enclosure_is_not_source']
            logR=(base+h*v if chart in MICRO else base+(c.ln(100)-base)*v+2*h*(1-v))
            if chart=='bridge_macro' and endpoints(v)==(mp.mpf(1),mp.mpf(1)):logR=c.ln(100)
            logR=c.mpf([max(endpoints(base)[0],endpoints(logR)[0]),min(endpoints(c.ln(100))[1],endpoints(logR)[1])])
            source=dict(original_radius_tree=packet['source_radius_tree'],
                exact_log_source='logRa+hb*s' if chart in MICRO else 'logRa+fraction*(log100-logRa)+2hb*(1-fraction)',
                exact_positive_hb_log=packet['source_width_log'],hb_numerical_enclosure_only=h,
                microscopic_radius_variation_retained_formally=True,absolute_radius_not_rounded_as_source=True,
                numeric_logR_bound_uses_positive_hb_enclosure=True,
                formal_hb_radius_correlation_evaluated=False,
                numerical_radius_bound_scope='Conservative cap-based enclosure; formal positive hb source retained separately, its correlations not evaluated')
            return logR,source
        if chart.startswith('switch_'):
            base=c.ln(100);h=packet['width_enclosure_is_not_source']
            logR=(base+h*v if chart in MICRO else base+c.ln(c.mpf(110)/100)*v+2*h*(1-v))
            if chart=='switch_power' and endpoints(v)==(mp.mpf(1),mp.mpf(1)):logR=c.ln(110)
            logR=c.mpf([max(endpoints(base)[0],endpoints(logR)[0]),min(endpoints(c.ln(110))[1],endpoints(logR)[1])])
            source=dict(original_radius_tree=packet.get('formal_log_radius_tree',packet.get('formal_radius_tree')),
                exact_log_source='log100+hb*s' if chart in MICRO else 'log100+fraction*log(110/100)+2hb*(1-fraction)',
                exact_positive_hb_log=packet['exact_positive_width_log'],hb_numerical_enclosure_only=h,
                microscopic_radius_variation_retained_formally=True,absolute_radius_not_rounded_as_source=True,
                numeric_logR_bound_uses_positive_hb_enclosure=True,
                formal_hb_radius_correlation_evaluated=False,
                numerical_radius_bound_scope='Conservative cap-based enclosure; formal positive hb source retained separately, its correlations not evaluated')
            return logR,source
        if chart=='reshape':return c.ln(110)+provider.reshape.T*v,'R=110*exp(T*phase), T=400*Abar'
        if chart=='inner_reference':
            ref=provider.reference
            return c.ln(110)+ref.reshape.T+(ref.loggap-8)*v,'R=110*exp(T+phase*(10*(logCstar+logPstar)-T-8))'
        if chart=='axial_restore':return self.logRref-8+v,'R=Rref*exp(-8+t)'
        if chart=='restore_buffer':return self.logRref+v,'R=Rref*exp(offset)'
        if chart=='actual_patch':return self.logRref-6+c.ln(v),'R=Rm*x=Rref*exp(-6)*x'
        if chart in ('Rh_reference','O2_slope'):return self.logRref+v,'R=Rref*exp(offset or y)'
        if chart in ('O2_axial','O2_buffer'):return self.logRref+packet['actual_y'],packet['exact_radius_source']
        if chart=='O3_slope_mu':return self.logRref+self.logP+v,'R=Rd*exp(t), log(Rd/Rref)=logPstar'
        if chart=='O3_power':return self.logRref+self.logP+1+self.params.Tw*v,'R=Rw*exp(Tw*phase), Rw=e*Rd'
        if chart in PULSE:
            if chart=='pulse_entrance':offset=v
            elif chart in ('pulse_main','pulse_exit','pulse_gap'):offset=v/self.params.mu
            else:offset=13/self.params.mu+v
            return self.logRp+offset,'R=Rp*exp(original pulse coordinate); Rv=Rp*exp(13/mu)'
        if chart in POST:
            if chart=='flatten':offset=v
            elif chart in ('outer_power','outer_angular'):
                outer=provider;offset=100+(outer.Lrel-4)*v if chart=='outer_power' else 100+outer.Lrel+v
            elif chart.startswith('steep_') or chart=='waiting':
                steep=provider;origin=100+steep.outer.Lrel
                offset=(origin+v if chart=='steep_entry' else origin+1+steep.Ts*v if chart=='steep_power'
                    else origin+1+steep.Ts+v if chart=='steep_exit' else origin+2+steep.Ts+steep.wait*v)
            else:
                # Pressure/stress companions retain the original heat
                # source object; radius offsets still come from that source.
                heat=getattr(provider,'heat',provider);steep=heat.steep
                offset=100+steep.outer.Lrel+2+steep.Ts+steep.wait+v
            return self.logRp+13/self.params.mu+offset,'R=Rv*exp(exact original flatten/angular/steep/waiting/collar offset)'
        raise ValueError('No original radius source for chart: '+chart)

    def normalized_sources(self,chart,Z,value,packet,provider,logR):
        c=self.ctx;logs=[c.mpf(0),self.logP,c.mpf(0),c.mpf(0),c.mpf(0),logR,c.ln(2)]
        amplitudes={UZ:{},UT:{},UR:{5:mp.mpf('.5'),6:mp.mpf('-.5')},P:{1:mp.mpf(2)}}
        if chart in MICRO or chart=='bridge_macro':
            logs[0]=packet.get('source_width_log',packet.get('exact_positive_width_log'))
            logs[2]=packet['factored_source_log_bases'][2]/2
            parent=packet.get('actual_parent_axial5_packet')
            logs[3]=c.ln(parent['F_actual_over_F0_axial5_coefficients'][0])
            names={UZ:'Uz',UT:'Utheta_over_current_Utheta',UR:'Ur_over_current_sqrt_R_over_2',P:'P_over_Pstar2'}
            grids={label:micro_terms(packet,name,phase_to_y=chart in MICRO) for label,name in names.items()}
            amplitudes[UT]={2:mp.mpf(1),3:mp.mpf(1),5:mp.mpf('.5'),6:mp.mpf('.5')}
        elif chart in COMMON:
            raw=packet['physical_velocity_pressure_y_Z_mixed4'];parent=packet.get('actual_inherited_axial5_packet')
            if chart=='switch_power':
                phi=packet['actual_postswitch_phi_axial5'][0]
                logs[4]=c.ln(2)/2+logR/2+provider.logF0+c.ln(phi)
            else:logs[4]=self.logP+parent['log_Utheta_over_Pstar_axial5_coefficients'][0]
            names={UZ:'Uz',UT:'Utheta_over_current_Utheta',UR:'Ur_over_current_sqrt_R_over_2',P:'P_over_Pstar2'}
            grids={label:ordinary_terms(raw[name]) for label,name in names.items()};amplitudes[UT]={4:mp.mpf(1)}
        elif chart=='actual_patch' or chart in PREPULSE:
            raw=packet['physical_velocity_pressure_y_Z_mixed4'];ur='Ur_over_sqrt_Rm_over_2' if chart=='actual_patch' else 'Ur_over_current_sqrt_R_over_2'
            names={UZ:'Uz',UT:'Utheta_over_Pstar',UR:ur,P:'P_over_Pstar2'}
            grids={label:ordinary_terms(raw[name]) for label,name in names.items()};amplitudes[UT]={1:mp.mpf(1)}
            if chart=='actual_patch':
                fixed_to_current=c.sqrt(c.mpf(value))
                grids[UR]={index:[(powers,coefficient/fixed_to_current) for powers,coefficient in terms] for index,terms in grids[UR].items()}
        else:
            raw=packet['physical_mixed_derivatives_total_order_le4'];grids={label:ordinary_terms(raw[label]) for label in (UZ,UT,UR,P)}
            if chart in PULSE:
                if chart=='pulse_entrance':offset=c.mpf(value)
                elif chart in ('pulse_main','pulse_exit','pulse_gap'):offset=c.mpf(value)/self.params.mu
                else:offset=13/self.params.mu+c.mpf(value)
                logs[4]=self.logP-(c.mpf('.5')+self.params.mu)*offset
                amplitudes[UT]={4:mp.mpf(1)};amplitudes[UZ]={4:mp.mpf(1)}
                amplitudes[UR].update({4:mp.mpf(1)})
            else:
                flatten=self.dispatch.provider('flatten')
                logs[4]=self.logP+sum(flatten.logEv2_parts.values(),c.mpf(0))/2
                amplitudes.update({label:{4:mp.mpf(1)} for label in (UZ,UT,UR)})
        return grids,tuple(logs),amplitudes

    def evaluate(self,chart,Z,coordinate,log_tau='-1',theta='0',axis=False):
        c=self.ctx;Z=c.mpf(Z);logtau=c.mpf(log_tau)
        if any(not mp.isfinite(v) for v in endpoints(logtau)):raise ValueError('Finite log(tau), tau>0 required')
        if axis and (chart!='core' or endpoints(c.mpf(coordinate))!=(mp.mpf(0),mp.mpf(0))):raise ValueError('Nonsingular axis requires core rho=0')
        provider=self.dispatch.provider(chart);dispatched=self.dispatch.evaluate(chart,Z,coordinate);packet=dispatched['source_packet']
        if chart=='core':
            native=self.core.physical_map(packet,axis=axis)
            # The native core mapper retains its general lambda exponents.
            # Evaluate the requested time sector here, rather than returning
            # only its three historically materialized sample sectors.
            bounds={}
            for index,components in native['cartesian_spatial_multiindices'].items():
                bounds[index]={}
                for component,parts in components.items():
                    bounds[index][component]={}
                    for label,row in parts.items():
                        scale=native['shared_physical_prefactor_bounds'][row['scale_key']]
                        norm=endpoints(row['absolute_upper'])[1]
                        bounds[index][component][label]=(None if not norm else
                            c.ln(c.mpf(norm))+scale['logLambda_term']+scale['amplitude_log_upper']+scale['physical_lambda_exponent']*logtau/2)
            timebounds={}
            for component,parts in native['first_fixed_x_physical_time_derivative'].items():
                timebounds[component]={label:(None if endpoints(row['absolute_upper'])[1]==0 else
                    c.ln(row['absolute_upper'])+row['logLambda_term']+row['amplitude_log_upper']+row['physical_lambda_exponent']*logtau/2)
                    for label,row in parts.items()}
            return dict(chart=chart,actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
                core_axis_nonsingular_native_map=native,requested_core_spatial_log_bounds=bounds,
                requested_core_time_log_bounds=timebounds,
                physical_lambda_relation='lambda^2*(1-Z^2)=tau; original nonsingular core Cartesian map',
                output_kind='source bounds in original core factorization; no point coefficient selection',
                requested_log_tau=logtau,whole_core_includes_axis_without_inverse_radius=True,
                full_point_physical_field_evaluation=False,admissible_stress_lift_constructed=False,temporal_recursion=False)
        logR,radius_source=self.radius(chart,coordinate,packet,provider)
        grids,logs,amplitudes=self.normalized_sources(chart,Z,coordinate,packet,provider,logR)
        uncapped_source_rows=chart in MICRO or chart=='bridge_macro'
        if theta is None:cosine=sine=c.mpf([-1,1])
        else:cosine=c.cos(c.mpf(theta));sine=c.sin(c.mpf(theta))
        # Negative physical lambda powers admit a uniform upper bound even
        # when the source Z box touches +/-1 (physical infinity, not axis).
        loglambda_bound=logtau/2
        spatial={}
        for i,j,b in INDICES:
            key='x'+str(i)+'_y'+str(j)+'_z'+str(b);spatial[key]={}
            for component in COMPONENTS:
                parts=cartesian_source_row(c,grids,component,i,j,b,Z,self.delta,cosine,sine,amplitudes)
                spatial[key][component]={label:log_row(c,terms,logs,gamma,loglambda_bound,uncapped_source_rows) for label,(terms,gamma) in parts.items()}
        # Time rotation differentiates no theta: fixed physical x leaves
        # the cylindrical basis fixed, while the implicit lambda varies.
        time={}
        basis={'ux':{UR:cosine,UT:-sine},'uy':{UR:sine,UT:cosine},'uz':{UZ:c.mpf(1)},'p':{P:c.mpf(1)}}
        for component,seeds in basis.items():
            time[component]={}
            for label,angular in seeds.items():
                terms,gamma=time_source_row(c,grids,label,Z,self.delta,amplitudes)
                time[component][label]=log_row(c,[(powers,value*angular) for powers,value in terms],logs,gamma,loglambda_bound,uncapped_source_rows)
        physical_log_coordinates=dict(log_r=logR/2+c.ln(2)/2+loglambda_bound,
            log_r_uses_lambda_lower_bound_only=True,exact_log_r='log(lambda)+log(2R)/2',
            exact_z='Z*lambda^(1-delta)',lambda_relation='lambda^2*(1-Z^2)=tau')
        lo,hi=endpoints(Z)
        if -1<lo and hi<1:
            loglambda=(logtau-c.ln(1-Z**2))/2
            physical_log_coordinates.update(log_lambda_enclosure=loglambda,log_r=logR/2+c.ln(2)/2+loglambda,
                log_r_uses_lambda_lower_bound_only=False)
        return dict(chart=chart,actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            Z=Z,requested_log_tau=logtau,cos_theta=cosine,sin_theta=sine,
            original_radius_source=radius_source,source_logR_enclosure=logR,positive_source_log_bases=dict(zip(BASES,logs)),
            physical_spatial_cartesian_mixed4=spatial,first_fixed_x_physical_time_derivative=time,
            physical_log_coordinates=physical_log_coordinates,
            mapping=dict(lambda_relation='lambda^2-lambda^(2delta)*z^2=tau; Z=z/lambda^(1-delta)',
                radial='R=r^2/(2lambda^2)',velocity='ur=lambda^-1 Ur; utheta/uz=lambda^(-1-delta) Utheta/Uz',
                pressure='p=lambda^(-2-2delta)P'),
            microscope_phase_to_logR_applied_before_source_bound=chart in MICRO,
            global_source_row_mode='uncapped_factored_rows' if uncapped_source_rows else 'provider_prebounded_mixed_rows',
            moving_cylindrical_basis_differentiated=True,normalization_not_differentiated_twice=True,
            physical_log_bounds_use_lambda_ge_sqrt_tau=True,
            chart_spatial4_and_first_time1_source_map_available=True,
            output_kind='signed factored bounds on the SAME physical source; not point coefficients or measured dynamics',
            full_point_physical_field_evaluation=False,physical_energy_integral_certified=False,
            admissible_stress_lift_constructed=False,temporal_recursion=False)

    def report(self):
        c=self.ctx;domains={'core':[0,4],'bridge_first':[0,1],'bridge_second':[1,2],'bridge_macro':[0,1],
            'switch_first':[0,1],'switch_second':[1,2],'switch_power':[0,1],'reshape':[0,1],
            'inner_reference':[0,1],'axial_restore':[0,1],'restore_buffer':[-7,-6],'actual_patch':[1,endpoints(c.exp(1))[1]],
            'Rh_reference':[-5,0],'O2_slope':[0,1],'O2_axial':[0,1],'O2_buffer':[0,11],'O3_slope_mu':[0,1],'O3_power':[0,1],
            'pulse_entrance':[0,endpoints(c.mpf('.02')/self.params.mu)[0]],'pulse_main':['.02','10'],
            'pulse_exit':[10,11],'pulse_gap':[11,'12.0001'],'pulse_gap_end':[-endpoints(1/self.params.mu)[0],-4],
            'pulse_end':[-4,0],'flatten':[0,100],'outer_power':[0,1],'outer_angular':[-4,0],
            'steep_entry':[0,1],'steep_power':[0,1],'steep_exit':[0,1],'waiting':[0,1],
            'heat_collar':[0,3],'heat_exterior':[3,mp.inf]}
        charts={}
        for chart in ROUTES:
            charts[chart]=self.evaluate(chart,[-1,1],domains[chart],theta=None)
            print('Physical source chart mapped: '+chart,flush=True)
        axis=self.evaluate('core',[-1,1],0,axis=True,theta=None)
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            whole_source_chart_physical_maps=charts,nonsingular_axis_map=axis,
            all_33_original_source_charts_physical_spatial4_time1_mapped=True,
            exact_original_chart_coordinates_and_amplitudes_retained=True,
            full_point_physical_field_evaluation=False,full_background_NS_validation=False,
            physical_energy_integral_certified=False,admissible_stress_lift_constructed=False,
            independently_bounded_flat_remainder=False,temporal_recursion=False,input_hashes=self.hashes)


def run():
    with mp.workdps(280):result=CompliantGlobalPhysicalAssembly().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Common physical source assembly generated; genuine recursion remains open',flush=True)
    return result


if __name__=='__main__':run()
