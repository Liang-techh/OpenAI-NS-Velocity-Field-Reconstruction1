"""Actual original Rh source radial jets and endpoint-retaining C0 averaging.

Linear leading densities are odd in the original inverse phase. Their
primitive and slow radial derivative, arbitrary phase endpoints and the
actual N-dependent nonlinear remainder are all bounded before integration.
Ordinary Z contributions retain the predecessor's independent enclosure;
the mixed radial/Z averaging evaluator remains required for C1 tightening.
"""
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_reference_whole_cell_integrals as whole
import lei_ren_part1_paper_compliant_current_native_signed_averaging as averaging

points=whole.points;base=whole.base;ep=whole.ep
HERE,PREFIX,sha=whole.HERE,whole.PREFIX,whole.sha
NAME=PREFIX+'current_original_reference_phase_averaged_integrals.json.gz'
RECEIPT=PREFIX+'current_original_reference_phase_averaged_integrals_check.json'
GATE='original_reference_actual_source_radial_jets_and_endpoint_retaining_C0_Nminus2_integrals_enclosed'


def second_exponential_remainder(c,argument):
    """R2(x)=int_0^1(1-t)exp(t*x)dt with a directed integrated tail."""
    x=c.mpf(argument)
    if max(abs(v) for v in ep(x))>1:raise ValueError('Original bounded A/N required')
    if ep(x)==(0,0):return c.mpf(1)/2
    result=c.mpf(1)/2;power=c.mpf(1);M=64
    for k in range(1,M+1):power*=x;result+=power/c.factorial(k+2)
    radius=c.mpf(max(abs(v) for v in ep(x)))
    tail=ep(c.exp(1)*radius**(M+1)/c.factorial(M+3))[1]
    return result+c.mpf((-tail,tail))


def magnitude_bound(atlas,value):
    """A source-factor bound, never a selected signed function value.

    The finite offset contains the radial cell; its upper endpoint is used
    only in a supremum bound. Original constant factor powers remain exact.
    """
    if value.zero:return atlas.scalar(0)
    magnitude=max(abs(v) for v in ep(value.coefficient))
    scale=whole.prior.FormalScale(atlas.bases,value.scale.powers,atlas.ctx.mpf(ep(value.scale.offset)[1]))
    return whole.prior.ScaledEnclosure(scale,magnitude,atlas.ledger)


def symmetric_bound(atlas,bound):
    if bound.zero:return bound
    magnitude=max(abs(v) for v in ep(bound.coefficient))
    scale=whole.prior.FormalScale(atlas.bases,bound.scale.powers,atlas.ctx.mpf(ep(bound.scale.offset)[1]))
    return whole.prior.ScaledEnclosure(scale,atlas.ctx.mpf((-magnitude,magnitude)),atlas.ledger)


def tighter_source_bound(atlas,centered,direct):
    if centered.zero:return centered,'exact_centered_zero'
    if direct.zero:return direct,'exact_direct_zero'
    c=atlas.ctx
    logratio=(centered.scale-direct.scale).evaluate()+c.ln(c.mpf(ep(centered.coefficient)[1]))-c.ln(c.mpf(ep(direct.coefficient)[1]))
    if ep(logratio)[0]>0:return direct,'direct_smaller_by_collected_source_factor_ratio'
    return centered,'centered_valid_bound_retained'


def leading_pairs(g,E,V,A,B):
    EA=g.c1mul(E,A);EEA=g.c1mul(E,EA)
    return dict(m=B,h=EA,k=g.c1add(g.c1mul(V,EA),g.c1mul(E,B)),
        e=g.c1add(g.c1scale(g.constant(2),g.c1mul(V,B)),g.c1scale(g.constant(-1),EEA)),p=EEA)


class OriginalReferencePhaseAveraging:
    def __init__(self,*,Z):
        self.owner=whole.OriginalReferenceWholeCells(Z=Z);self.atlas=a=self.owner.atlas
        self.ctx=c=a.ctx;self.family=self.owner.family;self.Z=self.owner.Z
        receipt=json.loads((HERE/whole.RECEIPT).read_bytes())
        parity=json.loads((HERE/averaging.RECEIPT).read_bytes())
        for checked,gate in ((receipt,whole.GATE),(parity,averaging.GATE)):
            if not checked.get('all_passed') or not checked.get(gate) or checked['source_family']!=self.family:
                raise ValueError('Accepted same-source whole-cell and exact centering evidence required')
        self.hashes=dict(self.owner.hashes)
        for checked,name in ((receipt,whole.RECEIPT),(parity,averaging.RECEIPT)):
            for filename,digest in {**checked['input_hashes'],name:sha(name)}.items():
                if sha(filename)!=digest:raise ValueError('Changed original averaging dependency: '+filename)
                if filename in self.hashes and self.hashes[filename]!=digest:raise ValueError('Original averaging families disagree')
                self.hashes[filename]=digest
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.radial_templates={};self.radial_compiled={};self.radial_pressure={}
        inputs=self.owner.templates['inputs'];radial=inputs[1:5]
        for name in ('E','V','b','p1','p2'):
            rows=[];compiled=[];pressure=[]
            for powers,expression in self.owner.templates['rows'][(name,0)]:
                r=powers[0]
                derivative=r*expression+sum(s.diff(expression,x)*rate*x for x,rate in
                    zip(radial,(s.Rational(1,10),s.Rational(1,10),s.Rational(1,5),s.Rational(1,5)),strict=True))
                rows.append((powers,derivative))
                compiled.append(s.lambdify(inputs,derivative,modules=[{'mpf':c.mpf},'mpmath']))
                pressure.append(tuple(s.lambdify(inputs,s.diff(derivative,P),modules=[{'mpf':c.mpf},'mpmath'])
                                      for P in self.owner.templates['pressure_symbols']))
            self.radial_templates[name]=tuple(rows);self.radial_compiled[name]=tuple(compiled)
            self.radial_pressure[name]=tuple(pressure)
        self.cache={}

    def radial_query(self,left,right):
        query=self.owner.roots(left,right);a=self.atlas;c=self.ctx;y=query['coordinate']
        z=a.rational(self.Z);q=1+z*z;f=c.exp(y/10);alpha=a.copy_interval(self.owner.inputs.alpha_enclosure)
        values=(z,f,c.mpf(5)/8*f,c.mpf(5)/12*f*f,c.mpf(5)/2*f*f,
                -alpha/q**2,4*alpha*z/q**3,alpha*(4-20*z*z)/q**4)
        roots={};records=[]
        for name,rows in self.radial_templates.items():
            derivative=a.scalar(0);terms=[]
            for i,(powers,expr) in enumerate(rows):
                coefficient=c.mpf(self.radial_compiled[name][i](*values))
                derivative=a.add(derivative,a.term(powers,coefficient,coordinate=y))
                late=c.mpf(0)
                for sensitivity,bound in zip(self.radial_pressure[name][i],(5,10,44),strict=True):
                    late+=abs(c.mpf(sensitivity(*values)))*c.exp(c.mpf(3)/5)*bound/(2*q*q)
                if ep(late)[1]:
                    r,p,d,ell=powers;upper=ep(late)[1]
                    derivative=a.add(derivative,a.term((r,p-1,d,ell),c.mpf((-upper,upper)),coordinate=y))
                terms.append(dict(original_factor_powers=powers,ordinary_y_coefficient=coefficient,
                    late_pressure_error_Pstar_power=-1 if ep(late)[1] else None,positive_late_error_budget=late))
            roots[name]={(0,0):query['roots'][name][(0,0)],(0,1):derivative}
            records.append(dict(original_input=name,ordinary_native_y_terms=terms))
        for name in ('a','t0'):roots[name]={(0,0):query['roots'][name][(0,0)],(0,1):a.scalar(0)}
        query.update(radial_roots=roots,radial_source_records=records)
        return query

    def cell(self,left,right,*,N):
        N=points.candidate_N(N);left=points.reference.reference_coordinate(left);right=points.reference.reference_coordinate(right)
        key=(left,right,N)
        if key in self.cache:return self.cache[key]
        a=self.atlas;c=self.ctx;g=whole.DirectedCoefficientAlgebra(a)
        with mp.workdps(c.dps+40):
            query=self.radial_query(left,right);kernel=query['kernel']
            if kernel.geometry!='signed_Mobius':raise ValueError('Original strict signed radial source required')
            coordinate=c.mpf((0,1));chart='E';primitive=kernel.primitives(coordinate,chart)
            # The accepted first-direction implicit chain rule is variable
            # agnostic. Its derivative input slot receives actual slow y
            # roots here, not ordinary Z rows or a total N*phase derivative.
            slow,proof=points.slow.slow_values(kernel,query['radial_roots'],coordinate,chart)
            BY,bounded=whole.bounded_signed_implicit_B_Z(a,kernel,query['radial_roots'],coordinate,chart)
            pair=lambda name:points.exact.source.C1Function(query['radial_roots'][name][(0,0)],query['radial_roots'][name][(0,1)])
            E,V=pair('E'),pair('V')
            A=points.exact.source.C1Function(primitive['A'],slow['A_Z_slow'])
            B=points.exact.source.C1Function(primitive['B_over_Pstar'],BY)
            leading=leading_pairs(g,E,V,A,B)
            exact,F=points.exact.coefficient_pairs(g,E,V,A,B,g.constant(N))
            x=whole.conditioned.bounded(A.value*(c.mpf(1)/N))
            Q=E.value*base.current.square(A.value)*a.scalar(second_exponential_remainder(c,x))
            FF=base.current.square(F.value);BB=base.current.square(B.value)
            remainder=dict(m=a.scalar(0),h=Q,k=a.add(V.value*Q,F.value*B.value),
                e=a.sum((-E.value*Q,BB,-FF*(c.mpf(1)/2))),p=a.add(E.value*Q,FF*(c.mpf(1)/2)))
            direct={name:a.add(exact[-1][name].value*(c.mpf(1)/N),exact[-2][name].value*(c.mpf(1)/N**2))
                    for name in points.exact.RATES}
            caps={name:dict(leading=magnitude_bound(a,row.value),slow_y=magnitude_bound(a,row.Z),
                nonlinear_remainder=magnitude_bound(a,remainder[name]),direct_density=magnitude_bound(a,direct[name]))
                for name,row in leading.items()}
            endpoints={which:self.owner.dispatcher.reference.owner.radius.evaluate(y=y,N=N)
                for which,y in (('left',left),('right',right))}
        record=dict(source_family=self.family,exact_radial_cell=[str(left),str(right)],candidate_N=N,
            original_closed_slow_y_source_terms=query['radial_source_records'],
            radial_derivative_contract=dict(source_coordinate='native reference log radius y',
                fixed_original_phi=True,true_phase_total_y_includes_N_f_phi_not_used_here=True,
                exact_profile_rates=['1/10','1/10','1/5','1/5'],radius_factor_derivative='r times original term',
                original_P0_y_exact_zero=True,original_L_y_exact_zero=True,
                reused_first_direction_slot='ordinary y roots in generic first-direction chain-rule slots',
                slow_implicit_proof=proof,bounded_implicit_product=bounded),
            exact_leading_signed_coefficients={name:dict(C0=v.value.record(),slow_y=v.Z.record()) for name,v in leading.items()},
            actual_N_dependent_nonlinear_remainder={name:v.record() for name,v in remainder.items()},
            direct_original_finite_N_density_enclosures={name:v.record() for name,v in direct.items()},
            original_full_phase_bound_inputs={name:{kind:v.record() for kind,v in row.items()} for name,row in caps.items()},
            actual_original_endpoint_phases=endpoints,original_phase_linear_native_y_slope=N,
            reflection_leading_mean_exact_zero=True,nonlinear_finite_N_mean_not_zeroed=True,
            native_radius_Jacobian=1,ordinary_y_function_derivatives_not_derivatives_of_caps=True)
        result=dict(caps=caps,record=record,query=query,leading=leading,remainder=remainder,F=F)
        self.cache[key]=result;return result

    def integrate(self,*,count,N):
        N=points.candidate_N(N)
        if type(count) is not int or not 1<=count<=4096:raise ValueError('Exact radial partition required')
        a=self.atlas;c=self.ctx;total={name:a.scalar(0) for name in points.exact.RATES}
        direct={name:a.scalar(0) for name in points.exact.RATES};cells=[]
        with mp.workdps(c.dps+40):
            for i in range(count):
                left=-5+s.Rational(5*i,count);right=-5+s.Rational(5*(i+1),count)
                query=self.cell(left,right,N=N);contributions={}
                for name,rate in points.exact.RATES.items():
                    rate=c.mpf(rate.numerator)/rate.denominator
                    wl,wr=c.exp(rate*a.rational(left)),c.exp(rate*a.rational(right))
                    mass=a.rational(right-left) if ep(rate)==(0,0) else (wr-wl)/rate
                    if ep(mass)[0]<=0:raise ArithmeticError('Original positive own-rate mass required')
                    cap=query['caps'][name];G=cap['leading']*(c.mpf(1)/2);Gy=cap['slow_y']*(c.mpf(1)/2)
                    each_cell_endpoint=G*(wl+wr)
                    # G is one original closed reference-source function on
                    # this entire window, with G(0)=G(1)=0 in fast phase.
                    # Actual internal traces cancel exactly, even across a
                    # modulo wrap. No cross-chart seam is presumed here.
                    endpoint=G*((wl if i==0 else c.mpf(0))+(wr if i==count-1 else c.mpf(0)))
                    slow=a.add(Gy,G*rate)*mass;rest=cap['nonlinear_remainder']*mass
                    bound=a.sum((endpoint,slow,rest))*(c.mpf(1)/N**2)
                    total[name]=a.add(total[name],magnitude_bound(a,bound))
                    direct[name]=a.add(direct[name],magnitude_bound(a,cap['direct_density']*mass))
                    contributions[name]=dict(each_cell_endpoint_absolute_bound_before_telescoping=each_cell_endpoint.record(),
                        retained_global_endpoint_term=endpoint.record(),slow_and_kernel_term=slow.record(),
                        actual_nonlinear_term=rest.record(),positive_own_rate_mass=mass,N_power=-2,
                        source_IBP_C0_enclosure=symmetric_bound(a,bound).record())
                cells.append(dict(source=query['record'],contributions=contributions))
        selected={name:tighter_source_bound(a,total[name],direct[name]) for name in total}
        return dict(source_family=self.family,original_Z_exact=str(self.Z),candidate_N=N,exact_window=['-5','0'],
            exact_radial_cells=count,C0_phase_averaged_contribution_enclosures={name:symmetric_bound(a,v).record() for name,v in total.items()},
            C0_direct_source_range_contribution_enclosures={name:symmetric_bound(a,v).record() for name,v in direct.items()},
            C0_effective_contribution_enclosures={name:symmetric_bound(a,v[0]).record() for name,v in selected.items()},
            effective_C0_bound_selection={name:v[1] for name,v in selected.items()},
            full_cell_source_and_IBP_records=cells,actual_global_endpoint_terms_retained=True,
            exact_internal_same_closed_source_endpoint_cancellation=True,
            no_cross_chart_seam_cancellation_assumed=True,
            actual_N_dependent_nonlinear_remainder_retained=True,integer_period_count_not_enumerated=True,
            direct_or_incoming_histories_not_assumed_zero=True,original_P0_not_reset=True,
            C1_Z_tightening_not_claimed_without_mixed_yZ_oracle=True,
            actual_terminal_defect_or_control_closure=False)


def run():
    begin=time.monotonic();owner=OriginalReferencePhaseAveraging(Z='.37');levels=[]
    for count,N in ((4,160),(16,160),(16,320)):
        levels.append(owner.integrate(count=count,N=N))
        print('Original reflection/IBP radial contribution:',count,'cells, N',N,flush=True)
    raw=gzip.decompress((HERE/whole.NAME).read_bytes());previous=json.loads(raw)
    retained=[dict(source_family=v['source_family'],original_Z_exact=v['original_Z_exact'],candidate_N=v['candidate_N'],
        cells=v['exact_original_cells'],ordinary_Z_contribution_enclosures={name:row['Z'] for name,row in
            v['actual_finite_N_five_density_contribution_enclosures'].items()})
        for v in previous['original_reference_whole_window_refinements']]
    report=dict(**{GATE:True},source_family=owner.family,actual_original_reference_phase_averaged_levels=levels,
        predecessor_original_Z_integrals_retained=retained,
        exact_reflection_and_own_rate_IBP_source=averaging.NAME,
        fixed_nonzero_Z_source_window_only=True,actual_original_radial_source_jets_installed=True,
        mixed_yZ_oracle_and_Z_averaging_installed=False,
        all_17_chart_or_24_cell_oracle_installed=False,actual_five_controls_installed=False,
        current_whole_N_selected=False,**dict.fromkeys(whole.point.source.inertial.profiles.loop.OPEN,False),
        input_hashes=owner.hashes,execution_seconds=time.monotonic()-begin,
        scope='Actual original Rh_reference slow radial jets and reflection-centered endpoint/own-rate C0 integration at candidate N, retaining the exact N-dependent nonlinear remainder. Predecessor ordinary Z contribution enclosures retained separately. Not mixed-yZ averaging, whole-Z terminal closure, all-route oracle, numerical controls/global N or recursive velocity.')
    (HERE/NAME).write_bytes(gzip.compress(json.dumps(base.encoded(report),indent=2).encode()+b'\n',compresslevel=9,mtime=0))
    return report


if __name__=='__main__':run()
