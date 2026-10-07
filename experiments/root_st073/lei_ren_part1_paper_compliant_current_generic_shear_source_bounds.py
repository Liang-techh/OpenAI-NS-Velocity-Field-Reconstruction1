"""Current whole-cover source norms and admitted noncore quotient derivative bounds.

Bounds are logarithms of absolute ordinary derivatives. No physical radius,
Pstar, positive microscopic width, or reciprocal small source is evaluated.
Checked whole covers bound functions; they never define point values. A
quotient requires a separate current analytic positive-function certificate.
"""
import json
import ast
import math
from pathlib import Path
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_patch_relaxed_inputs as patch

packets=patch.packets;HERE,PREFIX,sha=patch.HERE,patch.PREFIX,patch.sha
NAME=PREFIX+'current_generic_shear_source_bounds.json'
RECEIPT=PREFIX+'current_generic_shear_source_bounds_check.json'
GATE='current_original_noncore_generic_quotient_derivative_log_bounds_certified'
OPEN=patch.OPEN
MIDDLE=('reshape','inner_reference','axial_restore','restore_buffer','actual_patch')
O2=('Rh_reference','O2_slope','O2_axial','O2_buffer')
INNER=('bridge_first','bridge_second','bridge_macro','switch_first','switch_second','switch_power')
ORDERS=tuple((j,k) for j in range(3) for k in range(2))


class LogUpper:
    """|f|<=exp(log_upper), or exact zero. Sum uses max+log(term count)."""
    def __init__(self,c,log=None):
        self.ctx=c;self.log=log
        if log is not None:
            lo,hi=packets.recovery.endpoints(log)
            if not all(mp.isfinite(v) for v in (lo,hi)):
                raise ValueError('Finite logarithmic source bound required')
            self.log=c.mpf(hi)

    @classmethod
    def constant(cls,c,value):
        v=c.mpf(value);lo,hi=packets.recovery.endpoints(v)
        upper=max(abs(lo),abs(hi))
        return cls(c,None if upper==0 else c.ln(c.mpf(upper)))

    def __mul__(self,other):
        if self.ctx is not other.ctx:raise ValueError('Same directed source context required')
        return LogUpper(self.ctx,None if self.log is None or other.log is None else self.log+other.log)

    @classmethod
    def add(cls,c,values):
        values=[v for v in values if v.log is not None]
        if any(v.ctx is not c for v in values):raise ValueError('Same directed source context required')
        if not values:return cls(c)
        largest=max(packets.recovery.endpoints(v.log)[1] for v in values)
        return cls(c,c.mpf(largest)+c.ln(len(values)))

    def divide_positive(self,log_lower):
        return LogUpper(self.ctx,None if self.log is None else self.log-log_lower)

    def record(self):
        return dict(exact_zero=self.log is None,log_absolute_upper=self.log,
            source_value_or_exponential_not_materialized=True)


def modal_partial_bound(packet,row,k,radius_power=0):
    """Bound D_Z^k of a factored row, preserving its fixed source log bases."""
    if row.algebra is not packet.algebra or type(k) is not int or k<0 or k>row.order:
        raise ValueError('Same source factor algebra and available axial derivative required')
    c=packet.algebra.ctx;values=[]
    for exponents,jet in row.terms.items():
        coefficient=jet.coefficients[k]*math.factorial(k)
        value=LogUpper.constant(c,coefficient)
        if value.log is None:continue
        shift=sum((c.mpf(str(power))*base for power,base in zip(exponents,packet.algebra.logs)),c.mpf(0))
        shift+=radius_power*packet.provenance['logR_cover']
        values.append(LogUpper(c,value.log+shift))
    return LogUpper.add(c,values)


def row_table(packet,rows,orders=ORDERS,radius_power=0):
    return {(j,k):modal_partial_bound(packet,rows[j],k,radius_power) for j,k in orders}


def quotient_table(c,numerator,denominator,log_denominator_lower,orders=ORDERS):
    """Exact differentiated n=d*g, bounded using true d>=exp(log_lower)>0."""
    lo,hi=packets.recovery.endpoints(log_denominator_lower)
    if not all(mp.isfinite(v) for v in (lo,hi)):
        raise ValueError('Finite actual positive-source log lower bound required')
    result={}
    for j,k in sorted(orders,key=lambda key:(sum(key),key)):
        terms=[numerator[(j,k)]]
        for i in range(j+1):
            for ell in range(k+1):
                if i==ell==0:continue
                weight=LogUpper.constant(c,math.comb(j,i)*math.comb(k,ell))
                terms.append(weight*denominator[(i,ell)]*result[(j-i,k-ell)])
        result[(j,k)]=LogUpper.add(c,terms).divide_positive(log_denominator_lower)
    return result


def exact_quotient_theorem():
    y,Z=s.symbols('ordinary_logR Z');n=s.Function('full_signed_numerator')(y,Z)
    d=s.Function('actual_positive_denominator')(y,Z);g=s.Function('actual_quotient')(y,Z)
    checks={}
    for j,k in ORDERS:
        left=s.diff(d*g,y,j,Z,k)
        right=sum(math.comb(j,i)*math.comb(k,ell)*s.diff(d,y,i,Z,ell)*s.diff(g,y,j-i,Z,k-ell)
            for i in range(j+1) for ell in range(k+1))
        if s.expand(left-right)!=0:raise ArithmeticError('Mixed quotient Leibniz identity failed')
        checks['y%d_Z%d'%(j,k)]=True
    R=s.exp(y);f=s.Function('same_inertial_shape')(y)
    physical=[s.diff(s.exp(y/2)*f,y,j)/s.exp(y/2) for j in range(3)]
    for j in range(3):
        converted=sum(math.comb(j,i)*s.Rational(1,2)**(j-i)*physical[i] for i in range(j+1))
        if s.simplify(converted-s.diff(R*f,y,j)/R)!=0:
            raise ArithmeticError('Original inertial radius half-shift reapplied incorrectly')
        checks['physical_half_shift_to_R_numerator_y%d'%j]=True
    return dict(passed=True,exact_mixed_derivative_identities=checks,
        derivatives_are_ordinary_logR_and_ordinary_Z=True,
        numerator_radius_D_y_shift_applied_exactly_once=True)


def outer_positive_source_theorem():
    asts=packets.recovery.numeric.transport.SourceAST()
    wanted={('coordinates','qi'):'(1+square(z)).reciprocal()',
        ('reference','u'):'qi*c.exp(y/10)',
        ('slope','factor'):"c.exp(y/10-c.mpf('.6')*J)",
        ('slope','u'):'qi*factor',('axial','u'):'u1*root',
        ('axial','logu'):"parent['log_Utheta_over_Pstar_base_source']-t/2",
        ('axial','t'):'y-1',('axial','root'):'c.exp(-t/2)'}
    for (method,target),value in wanted.items():asts.expression('pre_pulse_mixed_C4',method,target,wanted=value)
    asts.method('pre_pulse_mixed_C4','slope_masses')
    asts.method('pre_pulse_mixed_C4','turnoff_derivatives')
    return dict(passed=True,original_actual_O2_profile_AST_bindings=asts.bindings,
        axial_factor_qi_lower='1/(1+Z^2)>=1/2, same original qi source',
        slope_integral_bound='0<=J(y)=integral_0^y sigma(s)ds<=1 for 0<=y<=1',
        axial_extent_source='exp(Md)=logPstar-11 from the same checked implicit source definition',
        original_slope_a_lower='.8+1.2*sigma(y)>=.8',original_axial_buffer_a_exact='2',
        input_hashes=asts.hashes)


def inner_positive_source_theorem():
    asts=packets.recovery.numeric.transport.SourceAST()
    wanted={('inner_bridge_profiles','inputs','at'):'self.core.normalized_jets(4,Z)',
        ('inner_bridge_profiles','actual','inputs'):'self.inputs(Z)',
        ('inner_bridge_profiles','actual','floor'):'endpoints(self.r/100)[0]',
        ('inner_bridge_profiles','actual','relative'):'IntervalTaylor(c,ell).exp()',
        ('inner_bridge_profiles','actual','actualphi'):"inputs['phi'].truncate(5)*relative",
        ('inner_switch_profiles','phase','phi'):"inp['phi']*ell.exp()",
        ('inner_switch_profiles','post','phi'):"(inp['phi']*correction.exp())*c.exp(-y*c.mpf(2)/5)",
        ('frozen_comparison_field','__init__','self.r'):'4*self.core.epsilon',
        ('inner_bridge_profiles','__init__','self.r'):'self.frozen.r',
        ('core_transfer','run','logLambda'):'4*logP+1000',
        ('core_transfer','run','eps'):'ctx.exp(-logLambda)'}
    for (stem,method,target),value in wanted.items():asts.expression(stem,method,target,wanted=value)
    def keyword(stem,method,name):
        node=asts.method(stem,method)
        found=[kw.value for item in ast.walk(node) if isinstance(item,ast.Return) and isinstance(item.value,ast.Call)
            for kw in item.value.keywords if kw.arg==name]
        if len(found)!=1:raise ValueError('Unique original defining source keyword required: '+name)
        return found[0]
    defining=keyword('inner_bridge_profiles','actual','source_integral_definitions')
    actual={kw.arg:ast.literal_eval(kw.value) for kw in defining.keywords}
    if actual['F']!='F=f*exp(-.5*integral_0^y chi*Dbar dt)' or actual['chi']!='1-(1-hb)*sigma(y/hb); hb=cstar*K^-100':
        raise ValueError('Actual cumulative log(F/f), comparison direction or shared width source changed')
    angular=ast.literal_eval(keyword('inner_switch_profiles','phase','angular_source'))
    expected='-.5*hb^2*[integral_0^min(phase,1)Dbar(100exp(hb*t),Z)dt + integral_0^max(phase-1,0)(1-sigma(t))*Dbar(100exp(hb*(1+t)),Z)dt] -.4*hb*integral_0^max(phase-1,0)sigma(t)dt'
    if angular!=expected:raise ValueError('Original actual two-switch cumulative angular source changed')
    asts.bindings['inner_bridge_profiles.actual.same_exact_cumulative_logF_definition']=True
    asts.bindings['inner_switch_profiles.phase.same_exact_cumulative_logF_definition']=True
    return dict(passed=True,original_actual_core_bridge_switch_AST_bindings=asts.bindings,
        bridge_inlet='rho4 source is Ra=4*epsilon; epsilon=exp(-logLambda), logLambda=4logPstar+1000',
        normalized_inverse_core_bound='I=log upper norm of 1/(Cstar*Fcore), hence log f(Ra)>=-logCstar-I',
        actual_bridge_ratio='actual K1 cumulative log(F/f) bound applied to the exact same integral of chi*Dbar, not comparison log(Fbar/f)',
        actual_two_switch_ratio='original same-width signed integral in phase() consumes the checked cumulative two-switch log bound',
        final_power_ratio='post() retains exact negative .4log(R/100) plus inherited switch correction; no moment or amplitude reset',
        input_hashes=asts.hashes)


def encode_table(table):
    return {'y%d_Z%d'%key:value.record() for key,value in sorted(table.items())}


class CurrentGenericSourceBounds:
    def __init__(self,service=None):
        self.patch=patch.CurrentPatchRelaxedInputs(service);self.service=self.patch.service
        self.ctx=self.patch.ctx;self.family=self.patch.family
        record=json.loads((HERE/patch.RECEIPT).read_bytes())
        if not record['all_passed'] or not record[patch.GATE] or record['source_family']!=self.family:
            raise ValueError('Same checked current whole middle-source relaxed input required')
        self.service.bind_hashes(record['input_hashes']);self.service.bind_hashes({patch.RECEIPT:sha(patch.RECEIPT)})
        self.theorem=exact_quotient_theorem();self.outer=outer_positive_source_theorem();self.inner=inner_positive_source_theorem()
        self.service.bind_hashes(self.outer['input_hashes']);self.service.bind_hashes(self.inner['input_hashes'])
        self.service.bind_hashes({Path(__file__).name:sha(Path(__file__).name)})
        self.conditions=['current whole original relaxed input through Rh; same pressure/source/five-history functions',
            'source bounds use checked whole covers, never point-source functions or selected cap values',
            'ordinary y rows already normalized; radial velocity half-shift not used as a shear derivative']
        self.lower=self.positive_denominators()

    def positive_denominators(self):
        c=self.ctx;read=lambda v:packets.interval(c,v)
        lower=lambda v:c.mpf(packets.recovery.endpoints(v)[0])
        source=self.patch.restore.reshape;join=source.rows['reference_join_bounds']
        norm=source.inner.rows['physical_norm_family']
        logC=read(norm['selected_logCstar']);logP=self.service.data['logP']
        A=read(join['Abar']);a=c.mpf('.7')
        # Same B/T expression as the checked reshape theorem: B>=-2A,
        # 0<=1-sigma<=1, y=log(R/110)>=0 and 1+Z^2<=2.
        shape=lower(-2*A-logC-logP-c.ln(2))
        # Original exact reference continuation starts at y=T=400A;
        # its normalized log swirl is y/10-logC-logP-log(1+Z^2).
        reference=lower(read(join['T'])/10-logC-logP-c.ln(2))
        restore_min=self.patch.restore.proof['original_u_over_Pstar_lower']
        patch_min=self.patch.proof['original_normalized_swirl_lower']
        ell={'reshape':shape,'inner_reference':reference,
            'axial_restore':lower(c.ln(restore_min)),
            'restore_buffer':lower(c.ln(restore_min)),
            'actual_patch':lower(c.ln(patch_min))}
        extent=logP-11
        ell.update(Rh_reference=lower(-c.mpf('.5')-c.ln(2)),
            O2_slope=lower(-c.mpf('.6')-c.ln(2)),
            O2_axial=lower(-extent/2-c.ln(2)),
            O2_buffer=lower(-extent/2-c.mpf('5.5')-c.ln(2)))
        original=source.inner.rows['global_exit_certificate'];ledger=source.inner.rows['K1_ledger']
        if original['field_prescription']['actual']!='chi=1-(1-epsilon)*sigma(y/h); logF=logf-.5 integralchi*Dbar dy; V=v-integralchi sqrt(R/2)F Ebar dy, Ra<R<=100':
            raise ValueError('Same original exact cumulative inner function required')
        if 'log(F/f)C2<=CqK^8(h+epsilon)/2' not in ledger['proof']['later_weight']:
            raise ValueError('K1 bound must be for actual cumulative log(F/f), not comparison-only data')
        if norm['physical_amplitude_units']['normalized_swirl']!='Fhat=Cstar*Fcore=exp(-Lambda*G)*Phi':
            raise ValueError('Exact original normalized inverse-core function units required')
        inverse_core_log=read(norm['normalized_core_log_weighted_norms']['3']['inverse_Cstar_Fcore'])
        ratio=read(ledger['decreasing_K_smallness_bounds']['actual_log_C2'])
        switch=read(join['short_switch_relative_log_C2_upper'])
        logRa=lower(c.ln(4)-(4*logP+1000))
        inner_E=lower((c.ln(2)+logRa)/2-logC-logP-inverse_core_log-ratio-switch-c.mpf('.4')*c.ln(c.mpf('1.1')))
        inner_a=lower(source.inner.proof['source_log_a_lower_bound'])
        ell.update({chart:inner_E for chart in INNER})
        records={}
        for chart,value in ell.items():
            amin=a if chart in MIDDLE else c.mpf('.8')
            logamin=inner_a if chart in INNER else lower(c.ln(amin))
            records[chart]=dict(log_E_positive_lower=value,log_C_positive_lower=lower(value+logamin),
                log_actual_a_positive_lower=logamin,source_function_positivity_not_inferred_from_saved_box=True)
        return dict(source_charts=records,
            reshape_proof='same checked actual B/T source, ||B||C2<=2A, y>=0; log E>=-2A-logCstar-logPstar-log2',
            reference_proof='same original positive power continuation, y>=T=400A; log E>=T/10-logCstar-logPstar-log2',
            restoration_and_patch_proof='same checked actual original normalized swirl lower bounds on whole source domains',
            actual_original_O2_source_theorem=self.outer,
            actual_original_inner_source_theorem=self.inner,
            inner_source_norm_components=dict(logRa=logRa,inverse_normalized_core_log_norm=inverse_core_log,
                whole_actual_cumulative_logF_ratio_C2_upper=ratio,two_switch_cumulative_logF_C2_upper=switch),
            inner_proof='same exact core rho4 inlet f>=exp(-logCstar-I), actual log(F/f)>=-L through100; original two-switch loss<=Lsw and post-power loss<=.4log1.1; R>=Ra. Positive a lower is the already checked exact hb/(2Kmax) bound',
            O2_buffer_full_source_domain='offset in[0,11], including the last two units beyond the smaller relaxed-buffer interval',
            minimum_a_from_current_original_whole_middle_theorems='.7',
            unresolved_source_charts=[chart for chart in packets.CHARTS if chart not in records],
            log_bounds_not_point_profile_definitions=True)

    def chart(self,chart):
        if chart not in packets.CHARTS:raise ValueError('Known original current source chart required')
        packet=self.service.saved(chart);c=self.ctx
        coverage=packet.provenance
        if packets.recovery.endpoints(coverage['Z_box'])!=packets.recovery.endpoints(c.mpf([-1,1])):
            raise ValueError('Whole axial source cover required for norm admission')
        E,V=packet.velocity['theta'],packet.velocity['axial']
        C=tuple(E[j]-2*E[j+1] for j in range(3));B=tuple(2*V[j+1] for j in range(3))
        field=packet.recover_original(self.service.data['delta'])
        stress=field['full_signed_stress_ordinary_y_rows'];shift=(0,.5,0,0)
        def inertial_rows(component):
            physical=[left+packet.algebra.shift(right,shift) for left,right in
                zip(stress['inertial_'+component+'_linear'],stress['inertial_'+component+'_quadratic'])]
            # Existing rows carry +1/2 from sqrt(R/2). R*shape carries +1:
            # apply only the remaining +1/2, then bind the formal R factor.
            return tuple(sum((physical[i]*(math.comb(j,i)*c.mpf('.5')**(j-i))
                for i in range(j+1)),packet.algebra.lift(0)) for j in range(3))
        nt,nz=inertial_rows('theta'),inertial_rows('axial')
        tables=dict(E=row_table(packet,E),V=row_table(packet,V),C=row_table(packet,C),B=row_table(packet,B),
            inertial_theta_numerator=row_table(packet,nt,radius_power=1),
            inertial_axial_numerator=row_table(packet,nz,radius_power=1))
        histories={key:row_table(packet,rows) for key,rows in packet.histories.items()}
        absolute=row_table(packet,packet.absolute_pressure)
        P0={(0,k):modal_partial_bound(packet,packet.P0,k) for k in range(2)}
        quotient={};positivity=None
        if chart in self.lower['source_charts']:
            positivity=self.lower['source_charts'][chart]
            for key,num,den,ell in (
                ('a','C','E','log_E_positive_lower'),('b','B','E','log_E_positive_lower'),
                ('p1','inertial_theta_numerator','E','log_E_positive_lower'),
                ('p2','inertial_axial_numerator','E','log_E_positive_lower'),
                ('t0','B','C','log_C_positive_lower')):
                quotient[key]=quotient_table(c,tables[num],tables[den],positivity[ell])
        return dict(chart=chart,source_family=self.family,original_source_provenance=coverage,
            actual_denominator_source_certificate=positivity,
            ordinary_mixed_source_log_norms={key:encode_table(table) for key,table in tables.items()},
            own_original_five_history_log_norms={key:encode_table(table) for key,table in histories.items()},
            absolute_pressure_log_norms=encode_table(absolute),separate_original_P0_log_norms=encode_table(P0),
            admitted_original_quotient_log_norms={key:encode_table(table) for key,table in quotient.items()},
            quotient_derivative_orders=[list(key) for key in ORDERS],
            quotient_derivatives_bound_on_checked_whole_chart=bool(quotient),
            numerator_all_signed_inertial_pressure_energy_meridional_terms_retained=True,
            original_radial_prefactor_inertial_shift_applied_once=True,
            original_width_inverse_and_source_factor_correlation_retained=True,
            saved_core_positive_sector_does_not_include_axis=chart=='core',
            source_point_values_or_inverse_width_radius_amplitude_not_materialized=True,
            phase_held_loop_jets_installed=False,whole_generic_scales_instantiated=False,
            **dict.fromkeys(OPEN,False))

    def run(self):
        records={chart:self.chart(chart) for chart in packets.CHARTS}
        result=dict(source_family=self.family,source_function_attachment_conditions=self.conditions,
            exact_quotient_and_radius_derivative_theorem=self.theorem,
            positive_noncore_source_denominator_theorem=self.lower,current_original_source_log_bound_charts=records,
            admitted_quotient_source_charts=list(self.lower['source_charts']),raw_whole_cover_norm_charts=list(records),
            **{GATE:True},**dict.fromkeys(OPEN,False),
            whole_upstream_source_derivative_norms_certified=False,whole_generic_scales_instantiated=False,
            phase_held_loop_primitive_derivative_bounds_certified=False,source_graph_ancestor_constructors_called=False,
            scope='Current whole-cover raw source norms on16 saved charts and current analytic positive-denominator quotient bounds through y2/Z1 on all15 noncore charts; core positive saved sector remains raw-only. No whole-loop domain/scale/N/modified field/recursion admission.',
            input_hashes=self.service.hashes)
        (HERE/NAME).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
        print('Current whole-cover source log norms and admitted15-chart quotient derivatives PASS',flush=True)
        return result


def run():return CurrentGenericSourceBounds().run()


if __name__=='__main__':run()
