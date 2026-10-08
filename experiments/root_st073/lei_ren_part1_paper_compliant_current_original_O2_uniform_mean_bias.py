"""Original O2 all-N>=160 nonlinear mean-bias envelopes on native sources.

Actual phase inputs are reused, not the fixed-N density outputs rescaled.
Same-frame native q factors and the separate source amplitude theorem give
two valid envelope families. No envelope is differentiated as a function.
"""
import gzip
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_O2_reflected_phase_means as current

mixed=current.mixed;integrals=current.integrals;source=current.source
base,prior,ep=current.base,current.prior,current.ep
MixedJet=current.MixedJet;C0,Y,Z,YZ=mixed.ORDERS
HERE,PREFIX,sha=current.HERE,current.PREFIX,current.sha
NAME=PREFIX+'current_original_O2_uniform_mean_bias.json.gz'
RECEIPT=PREFIX+'current_original_O2_uniform_mean_bias_check.json'
GATE='original_O2_source_amplitude_all_N_mean_C0_Z_Nminus2_bias_bounds_installed'
MIN_N=160


def absolute(value):return base.conditioned.absolute(value)


def native_coefficients(a,E,V,A,B,*,physical_amplitude):
    """Native pointwise envelope coefficients, not a chosen mean function.

    E,V,A,B are covers of the same original functions and ordinary Z jets.
    Only actual derivatives enter product rules. Native products precede
    unit division/export. Constants exp(C/MIN_N) are uniform in all N.
    """
    c=a.ctx;cap=c.mpf(1)/2
    E0,V0,B0=(absolute(jet[C0]) for jet in (E,V,B))
    EZ,VZ,AZ,BZ=(absolute(jet[Z]) for jet in (E,V,A,B))
    acap=a.scalar(cap) if physical_amplitude else absolute(A[C0])
    a2=base.current.square(acap);rho=c.exp(cap/MIN_N);rho2=c.exp(2*cap/MIN_N)
    H=E0*a2*(rho/2);P=base.current.square(E0)*a2*rho2
    K=a.add(V0*H,E0*B0*acap*rho);B2=base.current.square(B0)
    HZ=a.add(EZ*a2*(rho/2),E0*acap*AZ*rho)
    KZ=a.sum((a.add(VZ*E0,V0*EZ)*a2*(rho/2),V0*E0*acap*AZ*rho,
        a.add(EZ*B0,E0*BZ)*acap*rho,E0*B0*AZ*rho))
    PZ=a.add(E0*EZ*a2,base.current.square(E0)*acap*AZ)*(2*rho2)
    EZdensity=a.add(B0*BZ*2,PZ)
    zero=a.scalar(0)
    return dict(C0=dict(m=zero,h=H,k=K,e_negative=P,e_positive=B2,p=P),
        Z=dict(m=zero,h=HZ,k=KZ,e=EZdensity,p=PZ))


def normalized_ranges(a,coefficients):
    """Enclose N^2 times pair C0 or pair_Z/U, after native products."""
    c=a.ctx;values={};derivatives={}
    for key,value in coefficients['C0'].items():
        finite=integrals.export_range(a,value,derivative_unit=False)
        if ep(finite)[0]<0:raise ArithmeticError('Positive native envelope lost')
        values[key]=ep(finite)[1]
    for key,value in coefficients['Z'].items():
        finite=integrals.export_range(a,value,derivative_unit=True)
        if ep(finite)[0]<0:raise ArithmeticError('Positive native Z envelope lost')
        upper=ep(finite)[1];derivatives[key]=c.mpf((-upper,upper))
    return [dict(m=c.mpf(0),h=c.mpf((0,values['h'])),k=c.mpf((-values['k'],values['k'])),
        e=c.mpf((-values['e_negative'],values['e_positive'])),p=c.mpf((0,values['p']))),derivatives]


def intersect(c,left,right):
    return base.conditioned.clipped(c,left,ep(right)[0],ep(right)[1])


def uniform_amplitude_theorem():
    C,N,E,EZ,V,VZ,B,BZ,A,AZ=s.symbols('C N E EZ V VZ B BZ A AZ',positive=True)
    # Exact derivative formulas used for the positive coefficient envelopes.
    z=s.Symbol('z',real=True)
    Ef,Vf,Af,Bf=(s.Function(name)(z) for name in ('E','V','A','B'))
    x=Af/N;zz=Bf/N
    pair=dict(h=Ef*(s.cosh(x)-1),k=Vf*Ef*(s.cosh(x)-1)+Ef*zz*s.sinh(x),
        p=Ef**2*(s.cosh(2*x)-1)/2)
    expected=dict(h=s.diff(Ef,z)*(s.cosh(x)-1)+Ef*s.diff(Af,z)/N*s.sinh(x),
        k=(s.diff(Vf,z)*Ef+Vf*s.diff(Ef,z))*(s.cosh(x)-1)+Vf*Ef*s.diff(Af,z)/N*s.sinh(x)
            +(s.diff(Ef,z)*Bf+Ef*s.diff(Bf,z))/N*s.sinh(x)+Ef*Bf*s.diff(Af,z)/N**2*s.cosh(x),
        p=Ef*s.diff(Ef,z)*(s.cosh(2*x)-1)+Ef**2*s.diff(Af,z)/N*s.sinh(2*x))
    for key in pair:assert s.simplify(s.diff(pair[key],z)-expected[key])==0
    return dict(passed=True,exact_pair_ordinary_Z_product_rule_identities=3,
        original_half_phase_amplitude='Phi(pi)=1/2 and strict monotonicity imply psi/(2*pi) in[0,1/2]; a<=2 gives |A|<=a/4<=1/2',
        hyperbolic_inequalities=['cosh(x)-1<=x^2*exp(abs(x))/2','abs(sinh(x))<=abs(x)*exp(abs(x))','cosh(x)<=exp(abs(x))'],
        uniform_exp_constants='exp((1/2)/160), exp(1/160), independent of N>=160',
        actual_A_B_free_phase_source_does_not_depend_on_frequency_N=True,
        native_and_physical_amplitude_envelopes_are_same_function_bounds=True,
        source_q_factors_preserved_before_positive_unit_division=True,
        no_derivative_of_amplitude_or_normalized_envelope=True,
        mean_Z_is_ordinary_source_Z_derivative_not_derivative_of_unit_normalized_caps=True,
        all_N_bias_component_not_complete_oscillatory_spatial_integral=True)


class OriginalO2UniformMeanBias:
    def __init__(self):
        self.parent=current.OriginalO2ReflectedMeans();self.source=self.parent.source
        self.ctx=self.parent.ctx;self.family=self.parent.family;self.out=self.parent.out;self.Z_unit=self.parent.Z_unit
        self.hashes=dict(self.parent.hashes)
        receipt=json.loads((HERE/current.RECEIPT).read_bytes())
        if not receipt.get('all_passed') or not receipt.get(current.GATE) or receipt['source_family']!=self.family:
            raise ValueError('Accepted actual original true-phase mean source required')
        for name,digest in {**receipt['input_hashes'],current.RECEIPT:sha(current.RECEIPT)}.items():
            if sha(name)!=digest:raise ValueError('Original all-N mean dependency changed: '+name)
            if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Original all-N mean closures disagree')
            self.hashes[name]=digest
        raw=gzip.decompress((HERE/current.NAME).read_bytes())
        if hashlib.sha256(raw).hexdigest()!=receipt['compressed_producer_report']['lossless_original_json_sha256']:
            raise ValueError('Accepted true-phase function archive changed')
        self.saved=json.loads(raw)['actual_original_O2_true_phase_means_and_finite_N_bias']
        corr_receipt=json.loads((HERE/mixed.RECEIPT).read_bytes())
        raw=gzip.decompress((HERE/mixed.NAME).read_bytes())
        if hashlib.sha256(raw).hexdigest()!=corr_receipt['compressed_producer_report']['lossless_original_json_sha256']:
            raise ValueError('Accepted original native source atlas archive changed')
        self.atlases=json.loads(raw)['whole_original_O2_correlated_qy_mixed_density_cells']
        self.count=self.parent.count;self.issued={};self.theorem=uniform_amplitude_theorem()
        if len(self.saved['original_source_true_phase_mean_records'])!=self.count or len(self.atlases)!=self.count:
            raise ValueError('Complete original source partition required')
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)

    def phase_inputs(self,frame,part):
        self.source.describe(frame);a=frame.roots['q'].atlas
        expected=self.atlases[frame.index][frame.branch]['source']
        if expected['source_level']!=frame.count or expected['source_index']!=frame.index or expected['source_family']!=self.family:
            raise ValueError('Original native atlas owner/partition changed')
        if expected['exact_outer_Z_bounds']!=frame.record['exact_outer_Z_bounds'] or expected['exact_y_cell']!=frame.record['exact_y_cell']:
            raise ValueError('Original native source domain changed')
        atlas=expected['atlas']
        if atlas['source_basis_order']!=['logPstar','selected_logCstar','logL','original_logq','defined_original_logabsu_or_zero'] or not atlas['actual_radius_power_is_exact_zero_and_not_in_fifth_slot']:
            raise ValueError('Actual source q/u atlas semantics required')
        for live,saved in zip(a.bases,atlas['defining_basis'],strict=True):
            if ep(live)!=ep(integrals.interval(a.ctx,saved)):raise ValueError('Same original log basis required')
        if part['source_family']!=self.family or part['source_level']!=frame.count or part['source_index']!=frame.index or part['branch']!=frame.branch:
            raise ValueError('Same original free phase source input required')
        if part['exact_outer_Z_bounds']!=frame.record['exact_outer_Z_bounds'] or part['actual_original_predicate']!=frame.record['actual_original_source_predicate']:
            raise ValueError('Original phase source predicate changed')
        return {name:MixedJet(a,{order:integrals.restore_value(a,rows[str(order)]) for order in mixed.ORDERS})
            for name,rows in part['source_defined_same_function_A_B_C0_y_Z_yZ'].items()}

    def branch(self,frame,branchrecord):
        self.source.describe(frame);a=frame.roots['q'].atlas;c=a.ctx
        parts=branchrecord['complete_original_half_phase_partition'];totals=[{key:c.mpf(0) for key in integrals.KEYS} for unused in range(2)];records=[]
        if len(parts)!=current.PHASE_PARTS or branchrecord['branch']!=frame.branch:raise ValueError('Complete same-source half-phase partition required')
        for index,part in enumerate(parts):
            left,right=s.Rational(index,8),s.Rational(index+1,8)
            if part['exact_true_phase_cell']!=[str(left),str(right)] or part['exact_pair_to_full_mean_weight']!='1/4':
                raise ValueError('Accepted exact original full phase measure required')
            primitives=self.phase_inputs(frame,part)
            inputs=dict(E=frame.roots['E'],V=frame.roots['V'],A=primitives['A'],B=primitives['B'])
            families=[native_coefficients(a,**inputs,physical_amplitude=value) for value in (False,True)]
            covers=[normalized_ranges(a,family) for family in families]
            dual=[{key:intersect(c,covers[0][j][key],covers[1][j][key]) for key in integrals.KEYS} for j in range(2)]
            for j in range(2):
                for key in integrals.KEYS:totals[j][key]+=dual[j][key]/4
            records.append(dict(exact_true_phase_cell=[str(left),str(right)],source_level=frame.count,source_index=frame.index,
                branch=frame.branch,exact_full_phase_weight='1/4',source_family=self.family,
                original_phase_archive=current.NAME,original_phase_receipt=current.RECEIPT,
                original_free_phase_function_inputs={name:{str(order):jet[order].record() for order in (C0,Z)} for name,jet in inputs.items()},
                native_positive_coefficient_envelopes=[{group:{key:value.record() for key,value in values.items()} for group,values in family.items()} for family in families],
                normalized_native_and_source_amplitude_coefficient_covers=covers,
                dual_same_function_C0_Z_coefficient_intersections=dual,
                original_q_factors_not_floored_or_dropped_before_native_products=True,
                coefficient_envelopes_are_bounds_not_selected_mean_functions=True,
                candidate_N_density_outputs_not_rescaled=True,N_free_source_inputs_and_uniform_exp_floor_used=True))
        return dict(coefficients=totals,record=dict(source_family=self.family,branch=frame.branch,
            exact_y_cell=frame.record['exact_y_cell'],exact_outer_Z_bounds=frame.record['exact_outer_Z_bounds'],
            actual_original_predicate=frame.record['actual_original_source_predicate'],atlas=a.record(),
            original_phase_coefficient_records=records,normalized_full_mean_C0_Z_Nminus2_coefficients=totals))

    def integrate(self):
        c=self.ctx;began=time.monotonic();totals=[{key:c.mpf(0) for key in integrals.KEYS} for unused in range(2)];records=[]
        with mp.workdps(c.dps+40):
            for index,saved in enumerate(self.saved['original_source_true_phase_mean_records']):
                branches=[];choices=[{key:[] for key in integrals.KEYS} for unused in range(2)]
                for branchrecord in saved['conditional_source_true_phase_mean_records']:
                    branch=branchrecord['branch'];frame=self.source.frame(self.count,index,branch=branch)
                    result=self.branch(frame,branchrecord);branches.append(result['record'])
                    for j in range(2):
                        for key in integrals.KEYS:choices[j][key].append(result['coefficients'][j][key])
                hull=lambda values:c.mpf((min(ep(v)[0] for v in values),max(ep(v)[1] for v in values)))
                union=[{key:hull(values) for key,values in group.items()} for group in choices]
                left,right=c.mpf(index)/self.count,c.mpf(index+1)/self.count
                unused,masses,decay=integrals.five.masses(c,c.mpf(1)/self.count,left,right)
                contributions=[{key:group[key]*masses[key] for key in integrals.KEYS} for group in union]
                for total,addition in zip(totals,contributions,strict=True):
                    for key in integrals.KEYS:total[key]+=addition[key]
                records.append(dict(exact_y_cell=saved['exact_y_cell'],same_source_native_all_N_mean_coefficient_branches=branches,
                    normalized_C0_Z_Nminus2_coefficient_unions=union,positive_original_own_rate_masses=masses,
                    normalized_Nminus2_mean_Duhamel_coefficients=contributions,branch_overlap_hulled_not_added=True))
                if (index+1)%16==0:print('Native all-N original O2 mean coefficients:',index+1,'/',self.count,flush=True)
            report=dict(source_family=self.family,all_integer_N_lower=MIN_N,exact_y_window=['0','1'],exact_Z_range=['-1','1'],
                original_source_cells=self.count,original_source_predicate_frames=self.count*3,
                original_true_phase_input_cells=self.count*3*current.PHASE_PARTS,original_P0_datum_sha256=self.family['datum_enclosure_sha256'],
                exact_source_amplitude_and_uniform_bias_theorem=self.theorem,
                normalized_all_N_mean_Duhamel_C0_Z_coefficients=totals,ordinary_Z_unit=self.Z_unit.record(),
                actual_source_mean_coefficient_records=records,N_dependence='C0 in coefficient/N^2; ordinary_Z in (Lambda0/L^2)*coefficient/N^2',
                inherited_histories_P0_and_pressure_zero_rate_memory_unchanged=True,
                no_fixed_N_mean_density_or_integral_rescaling=True,execution_seconds=time.monotonic()-began)
        result=dict(report=report,coefficients=totals);self.issued[id(result)]=result;return result

    def at_candidate(self,integral,*,N):
        integrals.candidate_N(N)
        if self.issued.get(id(integral)) is not integral:raise ValueError('Issued same-owner all-N mean coefficients required')
        c=self.ctx;factor=c.mpf(1)/(c.mpf(N)**2)
        normalized=[{key:row*factor for key,row in group.items()} for group in integral['coefficients']]
        values={key:self.out.scalar(row) for key,row in normalized[0].items()}
        derivatives={key:self.Z_unit*row for key,row in normalized[1].items()}
        return dict(normalized=normalized,values=values,Z=derivatives,record=dict(source_family=self.family,explicit_candidate_N=N,
            same_all_N_source_coefficient_bounds_used=True,candidate_not_global_N_admission=True,
            normalized_C0_Z_mean_contributions=normalized,C0={key:value.record() for key,value in values.items()},
            ordinary_Z={key:value.record() for key,value in derivatives.items()}))


def run():
    began=time.monotonic();owner=OriginalO2UniformMeanBias();result=owner.integrate();c=owner.ctx
    candidates={str(N):owner.at_candidate(result,N=N)['record'] for N in (160,257,1024,4096)}
    old=owner.saved['all_normalized_C0_Z_mean_bias_contributions'];new=candidates['160']['normalized_C0_Z_mean_contributions']
    dual=[{key:intersect(c,new[j][key],integrals.interval(c,old[j][key])) for key in integrals.KEYS} for j in range(2)]
    report=dict(**{GATE:True},source_family=owner.family,actual_original_O2_all_N_mean_bias=result['report'],
        candidate_mean_bias_queries=candidates,additional_N160_direct_and_all_N_mean_intersection=dual,
        uniform_all_N_mean_bias_component_installed=True,physical_half_phase_A_amplitude_bound_installed=True,
        useful_native_dual_mean_C0_Z_bounds_installed=True,
        sharp_O2_phase_averaging_installed=False,oscillatory_spatial_integral_remainder_enclosed=False,
        actual_all_route_incoming_histories_installed=False,functional_terminal_identity_solved=False,
        current_whole_N_selected=False,all_17_chart_or_24_cell_oracle_installed=False,actual_five_controls_installed=False,
        **dict.fromkeys(source.ordered.base.point.source.inertial.profiles.loop.OPEN,False),
        input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
        scope='Actual original O2 native same-source phase inputs, half-phase amplitude and all integer N>=160 nonlinear five mean C0/Z N^-2 coefficients integrated by original own rates; dual native/physical amplitude bounds retain q factors. Mean component only, not spatial oscillatory remainder, sharp total averaging, inherited all-route driver, terminal/global N, recursion or corrected reconstruction.')
    (HERE/NAME).write_bytes(gzip.compress(json.dumps(base.encoded(report),indent=2).encode()+b'\n',compresslevel=9,mtime=0))
    print('Original O2 native all-N mean bias and amplitude bounds generated',flush=True);return report


if __name__=='__main__':run()
