"""True-phase means of the original finite-N O2 nonlinear five densities.

Reflection is an identity of the same original inverse function. Paired
integrands are integrated over true phase, then slow radius with original
own-rate masses. The oscillatory spatial contribution is still separate.
"""
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_O2_correlated_qy_density_jets as mixed

integrals=mixed.integrals;phase=mixed.phase;source=mixed.source
base,prior,ep=mixed.base,mixed.prior,mixed.ep
MixedJet=mixed.MixedJet;C0,Y,Z,YZ=mixed.ORDERS
HERE,PREFIX,sha=mixed.HERE,mixed.PREFIX,mixed.sha
NAME=PREFIX+'current_original_O2_reflected_phase_means.json.gz'
RECEIPT=PREFIX+'current_original_O2_reflected_phase_means_check.json'
GATE='original_O2_true_phase_five_finite_N_mean_C0_Z_and_Duhamel_bias_enclosed'
PHASE_PARTS=4


def stable_hyperbolic_changes(a,x):
    """Cosh(x)-1 and sinh(x), preserving tiny native x before conversion."""
    c=a.ctx;finite=integrals.positive.bounded(x[C0])
    if ep(finite)[0]<-1 or ep(finite)[1]>1:raise ValueError('Actual |x|<=1 required before stable hyperbolic conversion')
    square=finite**2;power=c.mpf(1);cm=c.mpf(0);sm=c.mpf(0)
    for k in range(32):
        cm+=power/c.factorial(2*k+2);sm+=power/c.factorial(2*k+1);power*=square
    ct=ep(c.exp(1)/c.factorial(66))[1];st=ep(c.exp(1)/c.factorial(65))[1]
    cm=base.conditioned.clipped(c,cm+c.mpf((0,ct)),c.mpf(1)/2,(c.exp(1)+c.exp(-1))/2-1)
    sm=base.conditioned.clipped(c,sm+c.mpf((0,st)),1,(c.exp(1)-c.exp(-1))/2)
    coshm1=base.current.square(x[C0])*cm;sinh=x[C0]*sm
    cosh=a.scalar(1+integrals.positive.bounded(coshm1))
    coshjet=MixedJet(a,{C0:coshm1,Y:x[Y]*sinh,Z:x[Z]*sinh,
        YZ:a.add(x[YZ]*sinh,x[Y]*x[Z]*cosh)})
    sinhjet=MixedJet(a,{C0:sinh,Y:x[Y]*cosh,Z:x[Z]*cosh,
        YZ:a.add(x[YZ]*cosh,x[Y]*x[Z]*sinh)})
    return coshjet,sinhjet


def reflected_pair_densities(a,E,V,A,B,*,N):
    """Exact half-sum at phi and1-phi; these are still phase functions."""
    integrals.candidate_N(N);inverse=a.ctx.mpf(1)/N;x=A*inverse;z=B*inverse
    cm,sh=stable_hyperbolic_changes(a,x);cm2,unused=stable_hyperbolic_changes(a,x*2)
    pp=E*E*cm2*(a.ctx.mpf(1)/2)
    return dict(m=MixedJet.constant(a,0),h=E*cm,k=V*E*cm+E*z*sh,e=z*z-pp,p=pp)


def reflection_theorem():
    phi,psi,T1,T2,a,E,V,q,nu,N=s.symbols('phi psi T1 T2 a E V q nu N',real=True)
    # Accepted original t0=0 moments: integral t=0, integral t^2=4*pi*q^2.
    A=a*(phi-psi/(2*s.pi))/2;B=-a*E*T1/(4*s.pi)
    reflected={phi:1-phi,psi:2*s.pi-psi,T1:-T1,T2:4*s.pi*q*q-T2}
    assert s.expand(A.xreplace(reflected)+A)==0 and s.expand(B.xreplace(reflected)+B)==0
    Phi=(psi+T2)/(2*s.pi*(1+2*q*q))
    assert s.cancel(Phi.xreplace(reflected)-(1-Phi))==0
    x,z=s.symbols('x z',real=True)
    def densities(A,B):
        de=E*(s.exp(A)-1);dv=B
        return dict(m=dv,h=de,k=V*de+E*dv+de*dv,
            e=2*V*dv+dv*dv-E*de-de*de/2,p=E*de+de*de/2)
    left,right=densities(x,z),densities(-x,-z)
    pair=dict(m=0,h=E*(s.cosh(x)-1),k=V*E*(s.cosh(x)-1)+E*z*s.sinh(x),
        e=z*z-E*E*(s.cosh(2*x)-1)/2,p=E*E*(s.cosh(2*x)-1)/2)
    for key in integrals.KEYS:assert s.simplify(((left[key]+right[key])/2-pair[key]).rewrite(s.exp))==0
    return dict(passed=True,exact_original_phase_A_B_reflection_identities=3,exact_five_nonlinear_pair_identities=5,
        original_t0_zero_full_period_moments=['integral_0^2pi t dpsi=0','integral_0^2pi t^2 dpsi=4*pi*q^2'],
        exact_true_phase_mean='M_j=2*integral_0^1/2 pair_j(y,Z,phi)dphi',
        ordinary_slow_y_Z_yZ_differentiate_same_function_reflection=True,
        original_q_y_not_set_to_zero=True,phase_held_fixed_not_total_rapid_derivative=True,
        exact_mean_m_zero_and_nonnegative_h_p_and_e_plus_p=True,
        nonlinear_bias_not_removed_by_linear_A_B_zero_means=True)


class OriginalO2ReflectedMeans:
    def __init__(self):
        self.phase=mixed.CorrelatedO2MixedDensity();self.source=self.phase.source
        self.ctx=self.phase.ctx;self.family=self.phase.family;self.hashes=dict(self.phase.hashes)
        receipt=json.loads((HERE/mixed.RECEIPT).read_bytes())
        if not receipt.get('all_passed') or not receipt.get(mixed.GATE) or receipt['source_family']!=self.family:
            raise ValueError('Accepted same-source correlated O2 mixed density interface required')
        for name,digest in {**receipt['input_hashes'],mixed.RECEIPT:sha(mixed.RECEIPT)}.items():
            if sha(name)!=digest:raise ValueError('Original O2 mean dependency changed: '+name)
            if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Original O2 mean closures disagree')
            self.hashes[name]=digest
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.count=min(self.source.source.parent.parent.levels)
        self.out=source.O2MixedAtlas(self.source.owner.inputs.frame,lower=-1,upper=1,logq=self.ctx.mpf(0))
        self.Z_unit=prior.ScaledEnclosure(prior.FormalScale(self.out.bases,(11,10,-2,0,0)),1,self.out.ledger)
        self.issued_means={};self.theorem=reflection_theorem()

    def pair(self,frame,left,right,*,N):
        integrals.candidate_N(N);self.source.describe(frame)
        left,right=s.Rational(left),s.Rational(right)
        if not 0<=left<right<=s.Rational(1,2):raise ValueError('One ordered exact true-half-phase cell required')
        a=frame.roots['q'].atlas;c=a.ctx
        target=c.mpf((ep(mixed.rational(c,left))[0],ep(mixed.rational(c,right))[1]))
        inverse=frame.kernel.evaluate(target,bits=12)
        if inverse['status']!='enclosed':raise ArithmeticError('Original full phase interval inverse not enclosed')
        selected=inverse['selected_inverse'];coordinate=c.mpf(selected['coordinate_interval'])
        if frame.branch=='regular':
            if selected['chart']!='psi':raise ValueError('Original regular psi inverse required')
        elif selected['chart']!='E':coordinate=frame.kernel.angles(coordinate,selected['chart'])[1]
        primitive=self.phase.primitive(frame,coordinate,phase=target)
        # Kernel C0 uses the equivalent canceled original primitive formula.
        # These stronger value ranges do not define or differentiate a cap.
        A,B=primitive['A'],primitive['B']
        A=MixedJet(a,{**A.rows,C0:integrals.restore_value(a,inverse['primitives']['A'])})
        B=MixedJet(a,{**B.rows,C0:integrals.restore_value(a,inverse['primitives']['B_over_Pstar'])})
        pair=reflected_pair_densities(a,frame.roots['E'],frame.roots['V'],A,B,N=N)
        record=dict(source_family=self.family,source_level=frame.count,source_index=frame.index,branch=frame.branch,
            exact_y_cell=frame.record['exact_y_cell'],exact_outer_Z_bounds=frame.record['exact_outer_Z_bounds'],
            actual_original_predicate=frame.record['actual_original_source_predicate'],
            exact_true_phase_cell=[str(left),str(right)],exact_reflected_true_phase_cell=[str(1-right),str(1-left)],
            selected_original_interval_inverse=selected,original_phase_interval=target,
            actual_backend_coordinate_fraction=coordinate,actual_backend_chart='psi' if frame.branch=='regular' else 'E',
            source_defined_same_function_A_B_C0_y_Z_yZ=dict(A=A.record(),B=B.record()),
            actual_reflected_pair_density_C0_y_Z_yZ={key:jet.record() for key,jet in pair.items()},
            exact_pair_to_full_mean_weight=str(2*(right-left)),
            independent_angle_phase_rectangle_is_outer_true_inverse_graph_cover=True,
            stronger_C0_from_canceled_original_primitive_not_selected_function=True,
            slow_jets_are_actual_derivatives_not_derivatives_of_C0_cover=True,
            paired_same_original_function_not_two_independent_caps=True,
            nonzero_original_V_term_in_k_retained=True,actual_positive_q_eta_and_source_family_retained=True,
            pair_integrand_not_yet_full_phase_mean=True)
        return dict(pair=pair,record=record)

    def phase_mean(self,frame,*,N):
        self.source.describe(frame);a=frame.roots['q'].atlas;c=a.ctx
        mean={key:MixedJet.constant(a,0) for key in integrals.KEYS};records=[]
        for index in range(PHASE_PARTS):
            left,right=s.Rational(index,2*PHASE_PARTS),s.Rational(index+1,2*PHASE_PARTS)
            result=self.pair(frame,left,right,N=N);weight=c.mpf(1)/PHASE_PARTS
            for key in integrals.KEYS:mean[key]=mean[key]+result['pair'][key]*weight
            records.append(result['record'])
        exports=[{key:integrals.export_range(a,jet[order],derivative_unit=order==Z) for key,jet in mean.items()} for order in (C0,Z)]
        if ep(exports[0]['h'])[0]<0 or ep(exports[0]['p'])[0]<0:
            raise ArithmeticError('Original positive E pair drift sign lost')
        return dict(mean=mean,export=exports,record=dict(branch=frame.branch,source_family=self.family,
            exact_y_cell=frame.record['exact_y_cell'],exact_outer_Z_bounds=frame.record['exact_outer_Z_bounds'],
            actual_original_predicate=frame.record['actual_original_source_predicate'],
            complete_original_half_phase_partition=records,exact_full_phase_measure='1',
            actual_true_phase_mean_C0_y_Z_yZ={key:jet.record() for key,jet in mean.items()},
            native_mean_values_not_selected_functions=True,normalized_mean_C0_Z_over_source_unit=exports,
            mean_m_and_all_its_slow_jets_exact_zero=True,mean_h_and_p_nonnegative=True,
            nonlinear_pressure_energy_bias_and_nonzero_V_mixed_moment_retained=True,
            true_phase_integral_range_not_point_quadrature=True))

    def integrate(self,*,N=160):
        integrals.candidate_N(N);c=self.ctx;began=time.monotonic()
        totals=[{key:c.mpf(0) for key in integrals.KEYS} for unused in range(2)];rows=[]
        with mp.workdps(c.dps+40):
            for index in range(self.count):
                branches=[];exports=[{key:[] for key in integrals.KEYS} for unused in range(2)]
                for branch in integrals.BRANCHES:
                    frame=self.source.frame(self.count,index,branch=branch);mean=self.phase_mean(frame,N=N)
                    branches.append(mean['record'])
                    for j in range(2):
                        for key in integrals.KEYS:exports[j][key].append(mean['export'][j][key])
                hull=lambda values:c.mpf((min(ep(v)[0] for v in values),max(ep(v)[1] for v in values)))
                union=[{key:hull(values) for key,values in group.items()} for group in exports]
                left,right=c.mpf(index)/self.count,c.mpf(index+1)/self.count
                unused,masses,decay=integrals.five.masses(c,c.mpf(1)/self.count,left,right)
                increments=[{key:group[key]*masses[key] for key in integrals.KEYS} for group in union]
                for total,add in zip(totals,increments,strict=True):
                    for key in integrals.KEYS:total[key]+=add[key]
                rows.append(dict(exact_y_cell=[str(s.Rational(index,self.count)),str(s.Rational(index+1,self.count))],
                    conditional_source_true_phase_mean_records=branches,normalized_mean_C0_Z_unions=union,
                    positive_own_rate_final_endpoint_masses=masses,normalized_mean_Duhamel_contributions=increments,
                    branch_overlap_hulled_not_summed=True,phase_reflection_weight_applied_once=True,
                    no_extra_R_Jacobian_or_period_count=True))
                if (index+1)%8==0:print('True original O2 reflected phase means:',index+1,'/',self.count,flush=True)
            values={key:self.out.scalar(totals[0][key]) for key in integrals.KEYS}
            derivatives={key:self.Z_unit*totals[1][key] for key in integrals.KEYS}
            report=dict(source_family=self.family,explicit_candidate_N=N,exact_y_window=['0','1'],exact_Z_range=['-1','1'],
                source_cells=self.count,source_predicate_frames=self.count*3,
                actual_true_phase_interval_inverse_queries=self.count*3*PHASE_PARTS,
                exact_half_phase_parts=PHASE_PARTS,own_rates=integrals.five.RATES,normalized_own_units=integrals.five.UNITS,
                true_phase_mean_Duhamel_C0={key:v.record() for key,v in values.items()},
                true_phase_mean_Duhamel_ordinary_Z={key:v.record() for key,v in derivatives.items()},
                mean_Z_unit=self.Z_unit.record(),all_normalized_C0_Z_mean_bias_contributions=totals,
                original_source_true_phase_mean_records=rows,original_phase_reflection_theorem=self.theorem,
                original_P0_datum_sha256=self.family['datum_enclosure_sha256'],
                exact_mean_m_zero_and_h_p_nonnegative=True,full_original_true_phase_means_enclosed=True,
                finite_N_nonlinear_bias_retained_and_integrated=True,
                oscillatory_spatial_density_integral_remainder_enclosed=False,
                mean_Duhamel_contribution_not_complete_oscillatory_contribution=True,
                same_source_incoming_history_or_P0_not_reset=True,functional_terminal_identity_solved=False,
                current_whole_N_selected=False,execution_seconds=time.monotonic()-began)
        result=dict(report=report,values=values,Z=derivatives);self.issued_means[id(result)]=result;return result


def run():
    began=time.monotonic();owner=OriginalO2ReflectedMeans();result=owner.integrate(N=160)
    report=dict(**{GATE:True},source_family=owner.family,actual_original_O2_true_phase_means_and_finite_N_bias=result['report'],
        true_phase_nonlinear_mean_C0_Z_layer_installed=True,actual_finite_N_mean_Duhamel_contributions_installed=True,
        sharp_O2_phase_averaging_installed=False,oscillatory_spatial_integral_remainder_enclosed=False,
        actual_all_route_incoming_histories_installed=False,functional_terminal_identity_solved=False,
        current_whole_N_selected=False,all_17_chart_or_24_cell_oracle_installed=False,actual_five_controls_installed=False,
        **dict.fromkeys(source.ordered.base.point.source.inertial.profiles.loop.OPEN,False),
        input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
        scope='Actual original O2 same-function reflection, true phase interval inverse partition, nonlinear five phase mean C0/Z and own-rate mean Duhamel contributions at candidate N160 over full original y/Z predicate union. Retains finite-N bias; oscillatory spatial contribution, sharp averaging, actual all-route incoming, terminal/global N, recursion and full reconstruction remain open.')
    (HERE/NAME).write_bytes(gzip.compress(json.dumps(base.encoded(report),indent=2).encode()+b'\n',compresslevel=9,mtime=0))
    print('Original O2 true nonlinear phase means and finite-N Duhamel bias generated',flush=True);return report


if __name__=='__main__':run()
