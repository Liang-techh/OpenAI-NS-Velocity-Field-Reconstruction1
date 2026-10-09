"""Actual second-switch functions and post-power transport through R110.

Actual first-exit inlets and anchored pressure/amplitude are reused. The
original shifted radius, true sigma weights and microscopic errors remain.
This is conditional two-frame source closure, not whole-axis/global closure.
"""
from fractions import Fraction
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_first_switch_functions as first
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_compliant_inner_bridge_profiles import logarithm

fields=first.fields;endpoint=first.endpoint;moments=first.moments
HERE,PREFIX,sha=first.HERE,first.PREFIX,first.sha
prior,base,ep=first.prior,first.base,first.ep
NAME=PREFIX+'current_original_second_switch_R110.json'
RECEIPT=PREFIX+'current_original_second_switch_R110_check.json'
GATE='original_actual_second_switch_and_two_frame_R110_function_transport_installed'


def source_bindings():
    from lei_ren_part1_paper_compliant_actual_switch_mixed_C4 import assignment_source_bindings
    controls=assignment_source_bindings('microswitch_mixed_C4','switch_controls',{
        'Vphase':'[quotient*0]*4'})
    transport=assignment_source_bindings('inner_switch_profiles','power_transport',{
        'angular':'difference(c.mpf(2)/5,2)*c.mpf(5)/4',
        'swirl':'difference(c.mpf(4)/5,2)*c.mpf(5)/6',
        'pressure':'difference(c.mpf(4)/5,1)*5'})
    logshape=assignment_source_bindings('inner_switch_profiles','inlet',{
        'B':"G*(-self.core.Lambda)+logarithm(phi)+logarithm(1+square(z))+c.ln(c.mpf(220))/2"})
    return dict(passed=True,original_second_V_zero_assignment=controls,
        original_positive_post_power_assignments=transport,original_log_shape_assignment=logshape,
        second_comparison_radius='100*exp(hb*(1+t)); fixed first-exit offset is retained',
        anchored_G_recipe_replaces_the_enclosure_only_Gbar_range_not_the_defining_formula=True)


class SecondSwitchFunctions:
    def __init__(self,first_owner,actual_fields,actual_histories):
        if not isinstance(first_owner,first.FirstSwitchFunctions):raise ValueError('Original first-switch mode owner required')
        self.first=first_owner;self.flow=f=first_owner.flow;self.c=f.c
        if set(actual_fields)!=set(('phi','V')) or set(actual_histories)!=set(moments.RATES):
            raise ValueError('Actual first-exit fields and all six incoming histories required')
        for row in (*actual_fields.values(),*actual_histories.values()):
            if len(row)!=6 or any(v.scale.bases is not f.logs or v.ledger is not f.ledger for v in row):
                raise ValueError('Same first-exit source basis/ledger and ordinary Z0..5 required')
        self.inlet_fields=actual_fields;self.inlet=actual_histories;self.cache={};self.posts={}
        # Current second-switch radius starts at100*exp(h), not100.
        self.D=[f.scale(row,first_owner.scalar_series(1-j,self.c.mpf(1))[0])
                for j,row in enumerate(first_owner.D)]
        self.D0=f.add(*self.D)

    def Ddelta(self,t):
        f=self.flow
        return f.add(*[f.scale(row,self.first.scalar_series(1-j,t,minus_one=True)[0])
                       for j,row in enumerate(self.D)])

    def angular(self,t,complement_prefix,weighted_D):
        f=self.flow;c=self.c
        sig=t-complement_prefix
        sig=c.mpf([max(mp.mpf(0),ep(sig)[0]),min(ep(t)[1],ep(sig)[1])])
        deterministic=[f.h*(-c.mpf(2)/5)*sig]+[f.scalar(0)]*5
        ell=f.add(deterministic,f.scale(weighted_D,-f.h*f.h*c.mpf('.5')))
        norm=f.norm(ell);upper=f.small_upper(norm)
        if ep(upper)[1]>mp.mpf('.5'):raise ValueError('Complete second-switch angular norm must be <=1/2')
        eta=norm*norm*norm*(c.exp(upper)/6)
        em=f.add(ell,f.scale(f.multiply(ell,ell),c.mpf('.5')),f.error_rows(eta))
        delta=f.multiply(self.inlet_fields['phi'],em)
        return dict(ell=ell,phi=f.add(self.inlet_fields['phi'],delta),phi_delta=delta,
                    sigma_prefix_mass=sig,nonzero_complete_exponential_error_norm=eta)

    def evaluate(self,phase):
        q=first.fraction(phase)
        if q in self.cache:return self.cache[q]
        op=self.first;f=self.flow;c=self.c;s=op.cv(q);h=f.h;V=self.inlet_fields['V']
        cuts=[q*Fraction(n,op.cells) for n in range(op.cells+1)] if q else [Fraction(0)]
        prefix=c.mpf(0);weighted_delta=[f.scalar(0)]*6;cells=[]
        for left,right in zip(cuts[:-1],cuts[1:]):
            t=op.box(left,right);mass=op.mass(left,right);next_prefix=prefix+mass
            if right==1:next_prefix=c.mpf('.5')
            Pcover=c.mpf([max(mp.mpf(0),ep(prefix)[0]),min(ep(s)[1],ep(next_prefix)[1])])
            delta=self.Ddelta(t)
            WD=f.add(f.scale(self.D0,Pcover),weighted_delta,f.scale(delta,c.mpf([0,ep(mass)[1]])))
            ang=self.angular(t,Pcover,WD)
            cells.append(dict(left=left,right=right,phi=ang['phi'],cutoff_mass=mass,weighted_D=WD))
            weighted_delta=f.add(weighted_delta,f.scale(delta,mass));prefix=next_prefix
        if q==1:prefix=c.mpf('.5')
        weighted_D=f.add(f.scale(self.D0,prefix),weighted_delta)
        at=self.angular(s,prefix,weighted_D)
        if q==0:at['phi']=self.inlet_fields['phi']
        phi=at['phi'];source_end=dict(H=f.scale(phi,2),M=V,K=f.scale(f.multiply(phi,V),2),
            A=f.multiply(V,V),B=f.multiply(phi,phi),C=f.multiply(phi,phi))
        histories={};increments={};ODE={};evidence={}
        for name,rate in moments.RATES.items():
            dm,decay_error=op.scalar_series(-rate,s,minus_one=True)
            if name in ('M','A'):
                # Exact constant-V equilibrium difference retains memory and
                # avoids subtracting independent positive kernel enclosures.
                change=f.scale(f.add(self.inlet[name],f.scale(source_end[name],-c.mpf(1)/rate)),dm)
                integral=None
            else:
                pieces=[]
                for cell in cells:
                    pc=cell['phi'];sources=dict(H=f.scale(pc,2),K=f.scale(f.multiply(pc,V),2),
                                               B=f.multiply(pc,pc),C=f.multiply(pc,pc))
                    width=op.cv(cell['right']-cell['left'])
                    backward,_=op.scalar_series(-rate,s-op.cv(cell['right']))
                    avg,_=op.scalar_series(-rate,width,average=True)
                    pieces.append(f.scale(sources[name],h*width*backward*avg))
                integral=f.add(*pieces) if pieces else [f.scalar(0)]*6
                change=f.add(f.scale(self.inlet[name],dm),integral)
            actual=f.add(self.inlet[name],change)
            if q==0:actual=self.inlet[name]
            histories[name]=actual;increments[name]=change
            ODE[name]=f.scale(f.add(source_end[name],f.scale(actual,-rate)),h)
            evidence[name]=dict(rate=rate,actual_incoming_decay_minus_one=dm,
                nonzero_complete_decay_error=decay_error,complete_signed_kernel_integral=integral,
                constant_V_equilibrium_identity_used=name in ('M','A'),microscopic_jacobian_hb_applied_once=True)
        sig=first.sigma_jets(c,s)[0];Dend=f.add(self.D0,self.Ddelta(s))
        angular_rate=f.add(f.scale(Dend,h*h*(-c.mpf('.5'))*(1-sig)),
                          [h*(-c.mpf(2)/5)*sig]+[f.scalar(0)]*5)
        phi_phase=f.multiply(phi,angular_rate)
        displacement,radial_error=op.scalar_series(1,1+s,minus_one=True)
        result=dict(actual_fields=dict(phi=phi,V=V),actual_six_histories=histories,
            signed_phi_increment=at['phi_delta'],signed_six_history_increments=increments,
            actual_phase_ODE_rows=dict(phi=phi_phase,V=[f.scalar(0)]*6,**ODE),
            log_phi_over_first_exit=at['ell'],sigma_prefix_mass=at['sigma_prefix_mass'],
            complement_prefix_mass=prefix,complete_weighted_original_D_integral=weighted_D,
            complete_angular_exponential_error_norm=at['nonzero_complete_exponential_error_norm'],
            history_integral_evidence=evidence,
            weighted_source_cells=[dict(left=[x['left'].numerator,x['left'].denominator],
                right=[x['right'].numerator,x['right'].denominator],cutoff_mass=x['cutoff_mass'],
                complete_weighted_D_prefix=x['weighted_D']) for x in cells],
            geometry=dict(phase=[q.numerator,q.denominator],radius=f.scalar(100)+displacement*100,
                signed_nonzero_radius_minus100=displacement*100,log_radius_over100=h*(1+s),
                complete_radial_exponential_error=radial_error,
                exact_source='R=100*exp(hb*(1+t)); dy=hb*dt'),
            actual_first_exit_V_preserved_exactly=True,ordinary_Z_orders=list(range(6)),
            shifted_original_comparison_radius_retained=True,actual_moments_not_used_to_redefine_direction=True)
        self.cache[q]=result;return result

    def post(self,radius):
        if type(radius) is int:Rq=Fraction(radius)
        elif isinstance(radius,tuple) and len(radius)==2 and all(type(v) is int for v in radius) and radius[1]>0:
            Rq=Fraction(*radius)
        else:raise ValueError('Exact integer or rational fixed post-switch radius required')
        if not 100<Rq<=110:raise ValueError('Fixed post-switch radius in(100,110] required')
        if Rq in self.posts:return self.posts[Rq]
        op=self.first;f=self.flow;c=self.c;R=op.cv(Rq);y=c.ln(R/100)
        if ep(y)[0]<=ep(2*f.hupper)[1]:raise ValueError('Fixed radius must exceed the original R2=100*exp(2hb)')
        incoming=self.evaluate((1,1));phi2=incoming['actual_fields']['phi'];V=incoming['actual_fields']['V']
        initial=incoming['actual_six_histories'];powers={};corrections={};errors={}
        for key,p in (('one',c.mpf(1)),('two',c.mpf(2)),('angular',c.mpf(2)/5),('squared',c.mpf(4)/5)):
            em,error=op.scalar_series(2*p,c.mpf(1),minus_one=True)
            baseline=f.scalar(c.exp(-p*y));powers[key]=baseline*(1+em)
            corrections[key]=baseline*em;errors[key]=error
        def positive(row):
            # The original common theta has0<theta<=1 by the R>=R2 guard.
            # Intersect only its exact positive weight, never select a value.
            value=f.ordinary_cover(row);lo,hi=ep(value)
            if hi<0:raise ArithmeticError('Post-power source contradicts positive kernel identity')
            return f.scalar(c.mpf([max(mp.mpf(0),lo),hi]))
        angular=positive((powers['angular']-powers['two'])*(c.mpf(5)/4))
        swirl=positive((powers['squared']-powers['two'])*(c.mpf(5)/6))
        pressure=positive((powers['squared']-powers['one'])*5)
        axial=positive(1-powers['one']);phi=f.scale(phi2,powers['angular'])
        histories=dict(H=f.add(f.scale(initial['H'],powers['two']),f.scale(phi2,angular)),
            M=f.add(f.scale(initial['M'],powers['one']),f.scale(V,axial)),
            K=f.add(f.scale(initial['K'],powers['two']),f.scale(f.multiply(phi2,V),angular)),
            A=f.add(f.scale(initial['A'],powers['one']),f.scale(f.multiply(V,V),axial)),
            B=f.add(f.scale(initial['B'],powers['two']),f.scale(f.multiply(phi2,phi2),swirl)),
            C=f.add(f.scale(initial['C'],powers['one']),f.scale(f.multiply(phi2,phi2),pressure)))
        sources=dict(H=f.scale(phi,2),M=V,K=f.scale(f.multiply(phi,V),2),A=f.multiply(V,V),
                     B=f.multiply(phi,phi),C=f.multiply(phi,phi))
        result=dict(actual_fields=dict(phi=phi,V=V),actual_six_histories=histories,
            actual_log_radius_ODE_rows=dict(phi=f.scale(phi,-c.mpf(2)/5),V=[f.scalar(0)]*6,
                **{name:f.add(sources[name],f.scale(row,-moments.RATES[name])) for name,row in histories.items()}),
            geometry=dict(radius=f.scalar(R),exact_fixed_radius=[Rq.numerator,Rq.denominator],
                original_source_R2='100*exp(2hb)',log_radius_over_R2=f.scalar(y)-f.h*2,
                R_ge_original_R2_checked=True),original_theta_powers=powers,
            signed_nonzero_micro_theta_corrections=corrections,complete_theta_exponential_errors=errors,
            original_positive_kernel_weights=dict(angular=angular,swirl=swirl,pressure=pressure,axial=axial),
            actual_second_exit_incoming_histories=initial,actual_first_exit_V_preserved_exactly=True,
            original_post_power_a='4/5',original_post_power_b=0,ordinary_Z_orders=list(range(6)))
        self.posts[Rq]=result;return result


def physical_at_radius(flow,R,Z,delta,p0,ratio,ratio2,actual_fields,histories):
    value=endpoint.recover_endpoint(flow,Z,delta,p0,ratio,ratio2,actual_fields,histories)
    f=flow;c=f.c;R=c.mpf(R);q=R/100;F0=f.factor((0,0,.5,0,0));F02=f.factor((0,0,1,0,0))
    value['radius']=f.scalar(R)
    value['physical_velocity_axial_coefficients']['Ur']=[v*c.sqrt(q) for v in value['physical_velocity_axial_coefficients']['Ur']]
    value['physical_velocity_axial_coefficients']['Utheta']=f.scale(value['physical_velocity_axial_coefficients']['Utheta'],c.sqrt(q))
    increment=f.scale(value['pressure_increment_true_axial5_divided_by_R_F0_squared'],f.scalar(R)*F02)
    value['physical_pressure_radial_increment_axial5']=increment
    value['physical_total_pressure_axial5']=f.add(value['physical_pressure_axis_axial5'],increment)
    r1=f.jet(ratio);r2=f.jet(ratio2)
    value['physical_cumulative_moment_axial5']=dict(
        Mtheta=f.scale(f.multiply(histories['H'],r1),f.scalar(R*R)*F0),
        Mz=f.scale(histories['M'],R),Mtheta_z=f.scale(f.multiply(histories['K'],r1),f.scalar(R*R)*F0),
        Mztheta=f.add(f.scale(histories['A'],R),f.scale(f.multiply(histories['B'],r2),f.scalar(-R*R)*F02)),Mp=increment)
    value['same_original_endpoint_first_switch_schema_ready']=False
    value['same_original_fixed_R_physical_recovery_ready']=True
    return value


class OriginalSecondSwitchR110:
    mode='genuine_original_second_switch_and_R110_actual_source_functions'
    def __init__(self,dps=500):
        self.upstream=first.OriginalFirstSwitchFunctions(dps);self.c=self.upstream.c
        self.family=self.upstream.family;self.hashes=dict(self.upstream.hashes);self.owners={};self.saved=None
        for name in (first.NAME,first.RECEIPT):
            row=json.loads((HERE/name).read_bytes())
            if not row.get(first.GATE) or name==first.RECEIPT and not row.get('all_passed'):
                raise ValueError('Accepted actual first-switch functions required')
            if row['source_family']!=self.family:raise ValueError('Same actual first-switch source family required')
            for path,digest in row['input_hashes'].items():fields.previous.bind(self.hashes,path,digest)
            fields.previous.bind(self.hashes,name,sha(name))
            if name==first.NAME:self.saved=row
        fields.previous.bind(self.hashes,Path(__file__).name,sha(Path(__file__).name))
        self.bindings=source_bindings()

    def owner(self,label):
        if label in self.owners:return self.owners[label]
        prior_owner=self.upstream.owner(label);flow=prior_owner.flow
        source=self.saved['packets'][label][-1]
        if (source['source_frame']!=label or source['source_family']!=self.family
                or source['implicit_source_sha256']!=self.upstream.upstream.upstream.fields.source
                or source['datum_enclosure_sha256']!=self.upstream.upstream.upstream.fields.datum
                or source['function_evaluation']['geometry']['phase']!=[1,1]):
            raise ValueError('Same original complete first-switch endpoint required')
        value=source['function_evaluation'];restore=lambda rows:[endpoint.restore_row(flow,v) for v in rows]
        actual_fields={name:restore(row) for name,row in value['actual_fields'].items()}
        histories={name:restore(row) for name,row in value['actual_six_histories'].items()}
        self.owners[label]=SecondSwitchFunctions(prior_owner,actual_fields,histories);return self.owners[label]

    def evaluate(self,label,phase):
        with mp.workdps(self.c.dps+40):value=self.owner(label).evaluate(phase)
        return dict(mode=self.mode,source_frame=label,source_family=self.family,
            function_evaluation=fields.serialized(value),whole_axis_functions_installed=False,
            original_upstream_micro_function_providers_complete=False)

    def inlet(self,label,radius=110):
        with mp.workdps(self.c.dps+40):
            op=self.owner(label);value=op.post(radius);f=op.flow;c=self.c
            _,proof,_=self.upstream.upstream.owner(label)
            records=self.upstream.upstream.upstream.fields.records
            Lambda=fields.previous.read_interval(c,records['core_transfer']['Lambda'])
            gradient=proof['original_gradient_coefficients']
            def ratios(multiplier):
                rows=[c.mpf(1)]
                for n in range(1,6):rows.append(sum((-multiplier*Lambda*gradient[j]*rows[n-1-j]
                                                    for j in range(n)),c.mpf(0))/n)
                return IntervalTaylor(c,rows)
            p0=IntervalTaylor(c,proof['original_P0_coefficients'][:6])
            delta=c.exp(fields.previous.read_interval(c,records['pressure_source']['compliant_source']['parameter_bounds']['log_delta']))
            R=value['geometry']['radius'].coefficient
            physical=physical_at_radius(f,R,proof['Z'],delta,p0,ratios(1),ratios(2),
                value['actual_fields'],value['actual_six_histories'])
            phi=IntervalTaylor(c,[f.ordinary_cover(row) for row in value['actual_fields']['phi']])
            if ep(phi[0])[0]<=0:raise ValueError('Positive actual fixed-radius phi required for original anchored log shape')
            G=IntervalTaylor(c,[fields.previous.read_interval(c,proof['actual_anchored_G'])]
                +[gradient[n-1]/n for n in range(1,6)])
            z=IntervalTaylor.variable(c,proof['Z'],5)
            B=G*(-Lambda)+logarithm(phi)+logarithm(1+z*z)+c.ln(2*R)/2
        return dict(mode=self.mode,source_frame=label,source_family=self.family,
            post_power_function_evaluation=fields.serialized(value),physical_function_evaluation=fields.serialized(physical),
            original_anchored_log_shape_B_axial5=list(B.coefficients),
            log_shape_identity='B=-Lambda*G_anchored+log(phi_R)+.5log(2R)+log(1+Z^2)',
            actual_anchored_G_not_Gbar_used=True,logCstar_cancelled_before_enclosure=True,
            original_P0_retained_separately=True,whole_axis_functions_installed=False,
            original_upstream_micro_function_providers_complete=False)


def run():
    began=time.monotonic();owner=OriginalSecondSwitchR110()
    result=dict(**{GATE:True},source_family=owner.family,original_source_bindings=owner.bindings,
        second_switch_packets={label:[owner.evaluate(label,q) for q in ((0,1),(1,2),(1,1))] for label in ('0','.5')},
        post_power_packets={label:[owner.inlet(label,R) for R in (105,110)] for label in ('0','.5')},
        **dict.fromkeys(fields.previous.OPEN,False),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
        scope='Original actual second-switch functions and true inherited power transport to105/110, physical Q/velocity/pressure and anchored log B at native0,.5. Conditional admitted upstream inlets; micro providers/whole Z/global reconstruction stay open.')
    (HERE/NAME).write_text(json.dumps(base.encoded(result),indent=2)+'\n',encoding='utf8')
    print('Actual original second switch and source-owned R110 functions recovered',flush=True);return result


if __name__=='__main__':run()
