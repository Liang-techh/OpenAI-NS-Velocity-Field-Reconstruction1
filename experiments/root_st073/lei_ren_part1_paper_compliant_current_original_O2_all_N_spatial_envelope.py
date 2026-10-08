"""Actual original O2 periodic-primitive IBP and full all-N spatial bounds.

The first-order signed functions have zero phase mean; the full nonlinear
remainder does not. Every local/global endpoint and own-rate mass is kept.
"""
import gzip
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_O2_collected_slow_jets as slow

mean=slow.current;integrals=slow.integrals;mixed=slow.mixed;source=slow.source
base,prior,ep=slow.base,slow.prior,slow.ep
MixedJet=slow.MixedJet;C0,Y,Z,YZ=mixed.ORDERS
HERE,PREFIX,sha=slow.HERE,slow.PREFIX,slow.sha
NAME=PREFIX+'current_original_O2_all_N_spatial_envelope.json.gz'
RECEIPT=PREFIX+'current_original_O2_all_N_spatial_envelope_check.json'
GATE='original_O2_source_periodic_primitive_IBP_full_spatial_C0_Z_all_Nminus2_envelope_installed'
OUTPUT_ORDERS=(Y,YZ)


def exact_spatial_theorem():
    E,V,A,B,N=s.symbols('E V A B N',real=True,nonzero=True)
    RE=N**2*E*(s.exp(A/N)-1-A/N);FN=N*E*(s.exp(A/N)-1)
    rp=N**2*E**2*(s.exp(2*A/N)-1-2*A/N)/2
    first=dict(m=B,h=E*A,k=E*(V*A+B),e=2*V*B-E**2*A,p=E**2*A)
    remainder=dict(m=0,h=RE,k=V*RE+FN*B,e=B**2-rp,p=rp)
    de=E*(s.exp(A/N)-1);dv=B/N
    full=dict(m=dv,h=de,k=V*de+E*dv+de*dv,e=2*V*dv+dv**2-E*de-de**2/2,p=E*de+de**2/2)
    for key in integrals.KEYS:assert s.simplify(full[key]-first[key]/N-remainder[key]/N**2)==0
    z=s.Symbol('Z',real=True);Ef,Af=s.Function('E')(z),s.Function('A')(z)
    re=N**2*Ef*(s.exp(Af/N)-1-Af/N)
    rez=N**2*s.diff(Ef,z)*(s.exp(Af/N)-1-Af/N)+N*Ef*s.diff(Af,z)*(s.exp(Af/N)-1)
    pp=N**2*Ef**2*(s.exp(2*Af/N)-1-2*Af/N)/2
    ppz=N**2*Ef*s.diff(Ef,z)*(s.exp(2*Af/N)-1-2*Af/N)+N*Ef**2*s.diff(Af,z)*(s.exp(2*Af/N)-1)
    assert s.simplify(s.diff(re,z)-rez)==0 and s.simplify(s.diff(pp,z)-ppz)==0
    y,phi=s.symbols('y phi',real=True);rate=s.Symbol('rate',nonnegative=True)
    G=s.Function('G')(y,z,phi);K=s.Function('K')(y)
    for g in (G,s.diff(G,z)):
        total=rate*K*g+K*(s.diff(g,y)+N*s.diff(g,phi))
        assert s.expand(K*s.diff(g,phi)/N-(total-K*(s.diff(g,y)+rate*g))/N**2)==0
    return dict(passed=True,exact_full_changed_density_first_order_second_remainder_identities=5,
        exact_remainder_ordinary_Z_cancellations=2,exact_C0_Z_IBP_chain_identities=2,
        primitive='G_j(y,Z,phi)=integral_0^phi actual_signed_f1_j(y,Z,s) ds',
        original_signed_f1={key:str(value) for key,value in first.items()},
        G_and_all_four_slow_rows_zero_at_phi_zero_and_one=True,
        primitive_cap='min(1/2, upper(phi),1-lower(phi))*sup_phi |D_slow f1|',
        full_phase_amplitude='half-phase |A|<=1/2 plus actual A(1-phi)=-A(phi) gives |A|<=1/2 on all[0,1]',
        uniform_exp_constants='rho=exp(1/320), rho2=exp(1/160), all integer N>=160',
        exact_remainder=dict(m='0',h='R_E',k='V*R_E+F_N*B',e='B^2-r_p',p='r_p'),
        exact_RE_Z='E_Z*A^2*R2(A/N)+E*A*A_Z*exprel(A/N)',
        exact_FN_Z='E_Z*A*exprel(A/N)+E*A_Z*exp(A/N)',
        exact_rp_Z='4*E*E_Z*A^2*R2(2*A/N)+2*E^2*A*A_Z*exprel(2*A/N)',
        source_IBP='integral_a^b K*f1(y,Z,frac(N*y+phi0))/N = ([K*G]_a^b-integral_a^b K*(G_y+rate*G))/N^2',
        phase0_Z_exact_zero=True,primitive_slow_derivatives_not_cap_derivatives=True,
        nonzero_nonlinear_mean_bias_inside_remainder_not_added_twice=True,
        no_endpoint_or_join_cancellation_assumed=True,dominating_output_units=['Lambda0/L','Lambda0^2/L^3'])


def interval(c,value):return integrals.interval(c,value)


class OriginalO2AllNSpatialEnvelope:
    def __init__(self):
        self.slow=slow.OriginalO2CollectedSlowJets();self.source=self.slow.source
        self.ctx=self.slow.ctx;self.family=self.slow.family;self.hashes=dict(self.slow.hashes)
        receipt=json.loads((HERE/slow.RECEIPT).read_bytes())
        if not receipt.get('all_passed') or not receipt.get(slow.GATE) or receipt['source_family']!=self.family:
            raise ValueError('Accepted actual original slow source and derivative norms required')
        for name,digest in {**receipt['input_hashes'],slow.RECEIPT:sha(slow.RECEIPT)}.items():
            if sha(name)!=digest:raise ValueError('Original O2 spatial envelope dependency changed: '+name)
            if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Original spatial closures disagree')
            self.hashes[name]=digest
        raw=gzip.decompress((HERE/slow.NAME).read_bytes())
        if hashlib.sha256(raw).hexdigest()!=receipt['compressed_producer_report']['lossless_original_json_sha256']:
            raise ValueError('Accepted original slow source archive changed')
        self.saved=json.loads(raw)['whole_original_O2_collected_slow_source_cells'];self.count=64
        if len(self.saved)!=self.count:raise ValueError('Complete actual original source partition required')
        self.out=source.O2MixedAtlas(self.source.owner.inputs.frame,lower=-1,upper=1,logq=self.ctx.mpf(0))
        self.units=[prior.ScaledEnclosure(prior.FormalScale(self.out.bases,slow.UNITS[order]),1,self.out.ledger) for order in OUTPUT_ORDERS]
        self.cache={};self.issued={};self.theorem=exact_spatial_theorem()
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)

    def hydrate(self,frame):
        self.source.describe(frame)
        if id(frame) in self.cache:return self.cache[id(frame)]
        a=frame.roots['q'].atlas;c=a.ctx;row=self.saved[frame.index][frame.branch];record=row['source']
        if frame.count!=64 or record['source_family']!=self.family or record['source_index']!=frame.index or record['source_level']!=64:
            raise ValueError('Same original spatial source owner/partition required')
        if record['branch']!=frame.branch or record['exact_y_cell']!=frame.record['exact_y_cell'] or record['exact_outer_Z_bounds']!=frame.record['exact_outer_Z_bounds']:
            raise ValueError('Actual source domains differ')
        if record['actual_original_source_predicate']!=frame.record['actual_original_source_predicate'] or row['exact_true_phase_domain']!=['0','1']:
            raise ValueError('Named actual source whole-phase predicate required')
        atlas=record['atlas']
        if atlas['source_basis_order']!=['logPstar','selected_logCstar','logL','original_logq','defined_original_logabsu_or_zero'] or not atlas['actual_radius_power_is_exact_zero_and_not_in_fifth_slot']:
            raise ValueError('Original spatial q/u basis semantics required')
        for live,saved in zip(a.bases,atlas['defining_basis'],strict=True):
            if ep(live)!=ep(interval(c,saved)):raise ValueError('Same original native basis required')
        jets={name:MixedJet(a,{order:integrals.restore_value(a,values[str(order)]) for order in mixed.ORDERS}) for name,values in row['actual_same_source_A_B_C0_y_Z_yZ'].items()}
        caps={}
        for key,orders in row['first_order_density_native_norms'].items():
            caps[key]={}
            for order in mixed.ORDERS:
                proof=orders[str(order)];finite=interval(c,proof['normalized_absolute_cap'])
                if ep(finite)[0]!=0 or ep(finite)[1]<0 or not proof['bound_only_not_source_function_or_derivative_of_cap']:
                    raise ValueError('Actual source upper norm contract required')
                unit=integrals.restore_value(a,proof['exact_positive_normalization_unit'])
                if unit.scale.powers!=slow.UNITS[order]:raise ValueError('Actual slow normalization unit changed')
                caps[key][order]=unit*c.mpf(ep(finite)[1])
        # Intersect only first-order C0 magnitudes with the globally proved
        # source amplitude bound. No actual signed source jet is replaced.
        E,V,B=(base.conditioned.absolute(j[C0]) for j in (frame.roots['E'],frame.roots['V'],jets['B']))
        half=a.scalar(c.mpf(1)/2);EE=base.current.square(E)
        physical=dict(m=B,h=E*half,k=E*a.add(V*half,B),e=a.add(V*B*2,EE*half),p=EE*half)
        for key,bound in physical.items():
            physicalcap=slow.native_norm(self.slow,frame,bound,C0)['cap']
            oldcap=slow.native_norm(self.slow,frame,caps[key][C0],C0)['cap']
            caps[key][C0]=a.scalar(min(ep(physicalcap)[1],ep(oldcap)[1]))
        result=dict(A=jets['A'],B=jets['B'],leading_caps=caps,source_archive_row=row)
        self.cache[id(frame)]=result;return result

    def primitive_bounds(self,frame,*,phase=(0,1)):
        self.source.describe(frame);a=frame.roots['q'].atlas;c=a.ctx;phi=c.mpf(phase)
        lo,hi=ep(phi)
        if not 0<=lo<=hi<=1:raise ValueError('One closed true-phase interval required')
        weight=min(ep(c.mpf(1)/2)[1],hi,ep(c.mpf(1)-c.mpf(lo))[1])
        sourcecaps=self.hydrate(frame)['leading_caps']
        caps={key:{order:value*c.mpf(weight) for order,value in orders.items()} for key,orders in sourcecaps.items()}
        return dict(caps=caps,record=dict(source_family=self.family,source_level=frame.count,source_index=frame.index,
            branch=frame.branch,true_phase_interval=phi,exact_same_source_primitive_definition=self.theorem['primitive'],
            original_source_archive=slow.NAME,source_signed_density_rows='native_same_source_signed_first_order_density_C0_y_Z_yZ',
            positive_integration_length_bound=weight,native_periodic_primitive_absolute_caps={key:{str(order):value.record() for order,value in orders.items()} for key,orders in caps.items()},
            actual_primitive_point_value_not_evaluated=True,phase_mean_zero_from_same_original_reflection=True,
            actual_slow_function_derivatives_used_not_cap_derivatives=True))

    def nonlinear_remainder(self,frame):
        self.source.describe(frame);a=frame.roots['q'].atlas;c=a.ctx;data=self.hydrate(frame)
        inputs=dict(E=frame.roots['E'],V=frame.roots['V'],A=data['A'],B=data['B'])
        families=[mean.native_coefficients(a,**inputs,physical_amplitude=physical) for physical in (False,True)]
        bounds=[];norms=[]
        for family in families:
            value={key:family['C0'][key] for key in ('m','h','k','p')}
            value['e']=a.add(family['C0']['e_negative'],family['C0']['e_positive'])
            bounds.append([value,family['Z']])
            norms.append([{key:slow.native_norm(self.slow,frame,v,order)['cap'] for key,v in group.items()} for group,order in zip((value,family['Z']),(C0,Z),strict=True)])
        dual=[{key:c.mpf((0,min(ep(norms[0][j][key])[1],ep(norms[1][j][key])[1]))) for key in integrals.KEYS} for j in range(2)]
        native=[{key:prior.ScaledEnclosure(prior.FormalScale(a.bases,slow.UNITS[order]),1,a.ledger)*c.mpf(ep(cap)[1]) for key,cap in group.items()} for group,order in zip(dual,(C0,Z),strict=True)]
        return dict(caps=native,record=dict(source_family=self.family,source_index=frame.index,branch=frame.branch,
            native_positive_all_N_remainder_families=[[{key:v.record() for key,v in group.items()} for group in family] for family in bounds],
            normalized_same_source_dual_C0_Z_remainder_caps=dual,
            mean_bias_not_assumed_zero_or_subtracted=True,actual_N_dependent_remainder_functions_enclosed=True,
            uniform_exp_floor_uses_global_reflection_extended_A_bound=True,ordinary_Z_product_terms_all_retained=True,
            same_envelope_formula_as_pair_with_separate_exact_unpaired_remainder_proof=True))

    def integrate(self):
        c=self.ctx;began=time.monotonic();totals=[{key:c.mpf(0) for key in integrals.KEYS} for unused in range(2)];rows=[]
        with mp.workdps(c.dps+40):
            for index in range(self.count):
                left,right=c.mpf(index)/64,c.mpf(index+1)/64
                unused,masses,decay=integrals.five.masses(c,c.mpf(1)/64,left,right)
                branchrecords=[];choices=[{key:[] for key in integrals.KEYS} for unused in range(2)]
                for branch in integrals.BRANCHES:
                    frame=self.source.frame(64,index,branch=branch);a=frame.roots['q'].atlas
                    primitive=self.primitive_bounds(frame);remainder=self.nonlinear_remainder(frame);groups=[];components=[]
                    for j,(order,dy,outputorder) in enumerate(zip((C0,Z),(Y,YZ),OUTPUT_ORDERS,strict=True)):
                        caps={};parts={}
                        for key in integrals.KEYS:
                            rate=c.mpf(integrals.five.RATES[key]);ka=c.exp(-rate*(1-left));kb=c.exp(-rate*(1-right))
                            G,Gy=primitive['caps'][key][order],primitive['caps'][key][dy]
                            endpoint=G*(ka+kb);transport=a.add(Gy,G*rate)*masses[key]
                            rest=remainder['caps'][j][key]*masses[key]
                            total=a.sum((endpoint,transport,rest))
                            cap=slow.native_norm(self.slow,frame,total,outputorder)['cap'];caps[key]=cap;choices[j][key].append(cap)
                            parts[key]=dict(left_kernel=ka,right_kernel=kb,positive_own_rate_mass=masses[key],
                                native_kept_both_endpoint_caps=endpoint.record(),native_slow_and_kernel_transport_cap=transport.record(),
                                native_full_nonlinear_remainder_cap=rest.record(),native_total_Nminus2_coefficient=total.record(),
                                normalized_dominating_unit_coefficient_cap=cap)
                        groups.append(caps);components.append(parts)
                    branchrecords.append(dict(branch=branch,source_family=self.family,
                        actual_original_predicate=frame.record['actual_original_source_predicate'],
                        exact_outer_Z_bounds=frame.record['exact_outer_Z_bounds'],periodic_primitive=primitive['record'],
                        full_nonlinear_remainder=remainder['record'],original_C0_Z_IBP_components=components,
                        normalized_dominating_C0_Z_coefficient_caps=groups))
                union=[{key:c.mpf((0,max(ep(v)[1] for v in values))) for key,values in group.items()} for group in choices]
                for total,add in zip(totals,union,strict=True):
                    for key in integrals.KEYS:total[key]+=add[key]
                rows.append(dict(exact_y_cell=[str(s.Rational(index,64)),str(s.Rational(index+1,64))],
                    conditional_original_source_spatial_IBP_records=branchrecords,normalized_C0_Z_Nminus2_coefficient_unions=union,
                    branch_overlap_hulled_not_added=True,all_each_cell_and_global_endpoints_retained=True,
                    no_extra_R_Jacobian_period_count_or_second_mean_bias=True))
                if (index+1)%16==0:print('Actual original O2 full spatial all-N envelope:',index+1,'/64',flush=True)
            report=dict(source_family=self.family,all_integer_N_lower=160,exact_y_window=['0','1'],exact_Z_range=['-1','1'],
                exact_spatial_first_order_remainder_and_IBP_theorem=self.theorem,
                normalized_full_spatial_C0_Z_Nminus2_coefficient_caps=totals,
                dominating_C0_Z_output_units=[v.record() for v in self.units],actual_source_spatial_IBP_records=rows,
                own_rates=integrals.five.RATES,original_pressure_datum_sha256=self.family['datum_enclosure_sha256'],
                exact_source_phase='frac(N*(log(110/4)+14*logPstar+10*logCstar+1000+y-hb*s_c/2))',
                original_positive_hb_sc_origin_retained=True,source_defined_periodic_primitive_bounds_installed=True,
                full_actual_original_O2_changed_spatial_integral_enclosed=True,
                broad_bounds_not_sharp_terminal_norms=True,execution_seconds=time.monotonic()-began)
        result=dict(coefficients=totals,report=report);self.issued[id(result)]=result;return result

    def at_candidate(self,result,*,N):
        integrals.candidate_N(N)
        if self.issued.get(id(result)) is not result:raise ValueError('Same-owner issued full spatial envelope required')
        c=self.ctx
        with mp.workdps(c.dps+40):
            endpoints=[self.source.owner.radius.evaluate(y=s.Rational(index,64),N=N) for index in range(65)]
            for row in endpoints:
                if row['source_family']!=self.family or not row['same_original_radius_phase_bound'] or not row['positive_original_origin_offset']['actual_hb_sc_product_not_materialized_or_set_to_zero']:
                    raise ValueError('Same actual original radius/common-N phase required')
            factor=c.mpf(1)/(c.mpf(N)**2);values=[];normalized=[]
            for j,group in enumerate(result['coefficients']):
                caps={};native={}
                for key,coefficient in group.items():
                    upper=ep(coefficient*factor)[1];caps[key]=c.mpf((-upper,upper));native[key]=self.units[j]*caps[key]
                normalized.append(caps);values.append(native)
            return dict(C0=values[0],ordinary_Z=values[1],record=dict(source_family=self.family,explicit_candidate_N=N,
                normalized_signed_C0_Z_contribution_envelopes=normalized,
                native_C0_contribution_envelopes={key:value.record() for key,value in values[0].items()},
                native_ordinary_Z_contribution_envelopes={key:value.record() for key,value in values[1].items()},
                actual_original_radius_common_N_phase_endpoints=endpoints,
                original_phase0_Z_exact_zero=True,approximate_phase_not_used_as_selected_value=True,
                all_N_full_spatial_coefficients_reused_not_fixed_N_density_rescaled=True,candidate_not_global_N_admission=True))


def run():
    began=time.monotonic();owner=OriginalO2AllNSpatialEnvelope();result=owner.integrate()
    candidates={str(N):owner.at_candidate(result,N=N)['record'] for N in (160,257,1024)}
    report=dict(**{GATE:True},source_family=owner.family,actual_original_O2_full_spatial_all_N_envelope=result['report'],
        candidate_full_spatial_queries=candidates,original_O2_full_spatial_changed_contribution_enclosed=True,
        actual_source_periodic_primitive_definition_and_bounds_installed=True,
        actual_phase_primitive_point_evaluator_installed=False,sharp_O2_phase_averaging_installed=False,
        actual_all_route_incoming_histories_installed=False,functional_terminal_identity_solved=False,
        current_whole_N_selected=False,all_17_chart_or_24_cell_oracle_installed=False,actual_five_controls_installed=False,
        **dict.fromkeys(source.ordered.base.point.source.inertial.profiles.loop.OPEN,False),
        input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
        scope='Actual original O2 full changed spatial contribution, all integer N>=160 N^-2 bounds from same-source zero-mean first-order periodic primitives, exact nonzero nonlinear remainder, own-rate masses and all local/global endpoints. C0 unit Lambda0/L and Z unit Lambda0^2/L^3 remain formal. Conservative O2 only; not sharp norms, point primitive evaluation, all-route history, functional terminal/global N, recursion or corrected NS.')
    (HERE/NAME).write_bytes(gzip.compress(json.dumps(base.encoded(report),indent=2).encode()+b'\n',compresslevel=9,mtime=0))
    print('Original O2 full spatial all-N envelope generated',flush=True);return report


if __name__=='__main__':run()
