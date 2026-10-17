"""Same-generator swirl and all-sector source jets over saved original data.

The exact positive absolute F0 squared amplitude remains a source atom.
Its logarithmic enclosure certifies a numerical approximation error, and
never selects the amplitude or substitutes an underflow zero for the field.
"""
from __future__ import annotations

import math
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from mpmath.libmp import fzero, mpf_cmp

from .original_hydro_source_point import (
    ROOT, DIRECTORY, FIRST, BRIDGE, DOMAIN, MASS, sha, digest, load,
    exact_iv, upper, ordinary_norm, width_errors, four_product,
)
from .original_leading_viscous_kernel import decode_pair, encode_pair
from .original_modal_source_point import _function_ast_guard
from .original_signed_C1_source import macro_nonlinear_enclosure, main_coefficient
from .original_source_jet_kernel import ordinary_product, ordinary_quotient

PARENT='lei_ren_part1_paper_compliant_current_original_hydro_source_C1_numeric.json'
PARENT_SHA256='660f377e75998bc631298d431fb2388fc89101f7b9cf58e5dd566bf6b2262033'
PRESSURE='lei_ren_part1_paper_compliant_current_original_signed_C1_source_numeric.json'
ANCHORED='lei_ren_part1_paper_compliant_anchored_axis_amplitude.json'


def relative_amplitude_from_dressed_moment(dressed,raw,*,normalized_base_proved):
    """Recover ordinary relative amplitude Z0..3 from one same-source moment.

    r0=1 follows from F0(Zbase)^2/F0(Zbase)^2, before interval division.
    Higher rows solve the ordinary product rule.  Both source moments stay
    enclosures; no endpoint, bound or midpoint supplies a signed value.
    """
    if normalized_base_proved is not True:raise ValueError('Exact normalized amplitude source identity required')
    if len(dressed)!=4 or len(raw)!=4 or any(v is None for v in list(dressed)+list(raw)):
        raise ValueError('Same-source dressed/raw ordinary moment rows through Z3 required')
    lower=raw[0]._mpi_[0] if hasattr(raw[0],'_mpi_') else raw[0]._mpf_
    if mpf_cmp(lower,fzero)<=0:raise ValueError('Positive raw moment denominator enclosure required')
    result=[raw[0]*0+1]
    for k in range(1,4):
        result.append((dressed[k]-sum((math.comb(k,j)*result[j]*raw[k-j] for j in range(k)),raw[0]*0))/raw[0])
    return result


def relative_amplitude_from_gradient(gradient,Lambda,*,source_proved):
    """Relative F0 squared jet from ordinary G'/G''/G''' and exact Lambda."""
    if source_proved is not True:raise ValueError('The same anchored amplitude gradient source is required')
    if len(gradient)!=3 or any(v is None for v in gradient):raise ValueError('Ordinary G gradient rows through Z2 required')
    a,b,d=[-2*Lambda*v for v in gradient]
    return [a*0+1,a,b+a*a,d+3*a*b+a*a*a]


def finite_width_swirl_comparison(ctx,iv,core,domain,weight,hcap,logh,*,source_proved,error_cap='1e-1000'):
    """Native2 comparison Phi/B/C in the ordinary order-three jet ring."""
    if source_proved is not True:raise ValueError('Same-core comparison and own-moment sources required')
    if hcap.a<=0 or (4*iv.exp(2*hcap)).b>=iv.mpf('4.1').a:
        raise ValueError('Comparison must stay in the admitted core continuation tube')
    d=domain['joint_core_domain'];rz=exact_iv(ctx,iv,d['Cauchy_Z_radius']);rr=iv.mpf(d['Cauchy_rho_radius'])
    major=exact_iv(ctx,iv,d['Xh_norm_majorants']['Phi'])*d['core_value_factor']
    floor=exact_iv(ctx,iv,domain['denominator_bounds']['Phi_core_modulus_lower'])
    tail=major*sum((weight**k/rz**k for k in range(1,4)),iv.mpf(0))
    gap=floor-tail
    if gap.a<=0:raise ValueError('Positive inverse-jet norm margin required')
    cover=major+tail
    log_radial_norm=4*iv.exp(2*hcap)*cover/rr/gap
    phi0=ordinary_norm(ctx,iv,core['Phi'],weight)
    phi_cover=phi0*iv.exp(2*hcap*log_radial_norm)
    changes={'Phi':2*log_radial_norm*phi_cover,
             'B':2*(2*ordinary_norm(ctx,iv,core['B'],weight)+phi_cover**2),
             'C':2*(ordinary_norm(ctx,iv,core['C'],weight)+phi_cover**2)}
    truth={};proofs={}
    for name,bound in changes.items():
        errors,proof=width_errors(ctx,iv,bound,weight,logh,orders=3,
                                  source_proved=source_proved,error_cap=error_cap)
        truth[name]=[v+e for v,e in zip(core[name],errors)];proofs[name]=proof
    return truth,{'passed':True,'phase':'2','core_seed_phase':'0 inlet',
                  'inverse_jet_margin_exact_mpf_tuples':[list(v) for v in gap._mpi_],
                  'log_core_radial_norm_upper_exact_mpf_tuples':[list(v) for v in log_radial_norm._mpi_],
                  'normalized_own_ODEs':['B_y+2B=Phi_bar^2','C_y+C=Phi_bar^2'],
                  'log_Phi_positive_IBP_source':True,'ordinary_error_proofs':proofs,
                  'source_width_materialized_or_selected':False}


def scaled_atom_enclosure(ctx,iv,coefficient,log_amplitude,*,source_proved,error_cap='1e-1000'):
    """Enclose A*coefficient without exponentiating unrepresentable log(A)."""
    if source_proved is not True:raise ValueError('Exact positive amplitude source and log enclosure required')
    if len(coefficient)!=3:raise ValueError('Complete dressed ordinary source rows through Z2 required')
    cap=iv.mpf(error_cap)
    if cap.a<=0 or cap.b>=1:raise ValueError('Positive approximation error cap smaller than one required')
    proofs=[];truth=[]
    for k,row in enumerate(coefficient):
        bound=upper(ctx,row)
        combined=None if bound==0 else log_amplitude+iv.log(iv.mpf(bound))
        if combined is not None and mpf_cmp(combined._mpi_[1],iv.log(cap)._mpi_[0])>=0:
            raise ValueError('The original positive amplitude does not admit this numerical approximation error')
        truth.append(iv.mpf(0) if bound==0 else iv.mpf([-1,1])*cap)
        proofs.append({'ordinary_Z_order':k,'coefficient_absolute_upper_exact_mpf_tuple':list(bound._mpf_),
                       'combined_amplitude_log_bound_exact_mpf_tuples':None if combined is None else [list(v) for v in combined._mpi_],
                       'numerical_absolute_error_cap':error_cap,'amplitude_exponentiated_or_selected':False,'passed':True})
    return [ctx.mpf(0)]*3,truth,{'passed':True,'orders':proofs,
        'zero_is_controlled_numerical_approximation_only':True,'actual_positive_amplitude_is_zero':False,
        'underflow_used_to_define_field':False}


class OriginalSwirlSourcePoint:
    def __init__(self,*,directory=DIRECTORY):
        self.directory=Path(directory)
        if sha(self.directory/PARENT)!=PARENT_SHA256:raise ValueError('Pinned accepted hydro parent required')
        self.parent=load(self.directory/PARENT)
        if self.parent['checks']['passed'] is not True:raise ValueError('Passed parent required')
        self.input_hashes=dict(self.parent['input_hashes']);self.input_hashes[PARENT]=PARENT_SHA256
        self.implementation_hashes=dict(self.parent['implementation_hashes'])
        if any(name not in self.input_hashes for name in (FIRST,BRIDGE,DOMAIN,MASS,PRESSURE,ANCHORED)):
            raise ValueError('Every hydrated source must have an accepted input pin')
        for base,group in ((self.directory,self.input_hashes),(ROOT,self.implementation_hashes)):
            for name,expected in group.items():
                if sha(base/name)!=expected:raise ValueError('Changed accepted pin: '+name)
        self.first=load(self.directory/FIRST);self.bridge=load(self.directory/BRIDGE)
        self.domain=load(self.directory/DOMAIN);self.mass=load(self.directory/MASS)
        self.pressure=load(self.directory/PRESSURE)
        self.anchored=load(self.directory/ANCHORED)
        if digest(self.mass['interior_mass'])!=self.parent['inherited_mass_packet_digest']:
            raise ValueError('Saved original mass changed')
        if digest(self.pressure['source_jet'])!=self.parent['inherited_pressure_source_digest']:
            raise ValueError('Saved pressure source changed')
        self.ctx=mp.mp.clone();self.ctx.dps=180
        self.iv=MPIntervalContext();self.iv.dps=220
        self.weight=exact_iv(self.ctx,self.iv,self.bridge['axial_jet_weight'])
        self.hcap=exact_iv(self.ctx,self.iv,self.bridge['width_enclosure_cap'])
        self.logh=exact_iv(self.ctx,self.iv,self.bridge['source_log_hb_enclosure'])
        self.logA=exact_iv(self.ctx,self.iv,self.bridge['swirl_scale_log_enclosure'])

    def source_guards(self):
        base=ROOT/'experiments/root_st073';src=ROOT/'src/openai_ns_reconstruction'
        specs=(
            ('swirl_dressing_direction',base/'lei_ren_part1_paper_compliant_inner_bridge_profiles.py','direction',None,
             ["b=dress(moments[MZT]['swirl'],ratios2)","p=dress(moments[MP],ratios2)",
              "(-(z*b)*(2*delta)+d*derivative(b)+(z*p)*(2*(1+delta))-d*derivative(p))/(2*L)"]),
            ('swirl_mode_sources',base/'lei_ren_part1_paper_compliant_macro_signed_integrals.py','_mode_rows',None,
             ["[base['drive_swirl'],one['drive_swirl']-base['drive_swirl'],two['drive_swirl']-base['drive_swirl']]",
              "delta_one={MTH:_zero(phi),MZ:delta[MZ],MTHZ:_zero(phi),MZT:{'axial':delta[MZT]['axial'],'swirl':_zero(phi)},MP:delta[MP]}",
              "delta_two={MTH:delta[MTH],MZ:_zero(phi),MTHZ:delta[MTHZ],MZT:{'axial':_zero(phi),'swirl':delta[MZT]['swirl']},MP:_zero(phi)}"]),
            ('relative_amplitude_seed',src/'original_first_collar_jets.py','inlet_jets','OriginalFirstCollarJets',
             ["ratios=[c.mpf(1)]","logarithmic_derivative=[2*v/eps for v in fixed['ell_Z_taylor'][:order]]",
              "ratios.append(sum((logarithmic_derivative[i]*ratios[k-1-i] for i in range(k)),c.mpf(0))/k)"]),
            ('saved_dressed_C',src/'original_first_collar_jets.py','evaluate','OriginalFirstCollarJets',
             ["values['C']*inlet[label]['F0_squared_ratio']"]),
            ('swirl_own_histories',base/'lei_ren_part1_paper_compliant_bridge_mixed_C4.py','formal_comparison_source','CompliantBridgeMixedC4',
             ["('B',2,'phi_bar(t,Z)^2')","('C',1,'phi_bar(t,Z)^2')",
              "dict(initial='log(Phi_core(4,Z))',integrand='alpha(t)*D_t log(Phi_core(4exp(t),Z))',bounds=['0','y'])"]),
            ('absolute_swirl_source',base/'lei_ren_part1_paper_compliant_actual_bridge_integrals.py','__init__','CompliantActualBridgeIntegrals',
             ["self.swirl_log=c.mpf([endpoints(-2*self.core.logC-2*self.core.Lambda*self.core.Gbar)[0],endpoints(-2*self.core.logC)[1]])"]),
            ('anchored_amplitude',base/'lei_ren_part1_paper_compliant_anchored_axis_amplitude.py','evaluate','CompliantAnchoredAxisAmplitude',
             ["logF0=-self.core.logC-self.core.Lambda*actual_G","self.gradient_jets(z)",
              "actual_G=c.mpf(0) if anchored else c.mpf([max(mp.mpf(0),endpoints(G.real)[0]),endpoints(G.real)[1]])"]),
            ('anchored_gradient',base/'lei_ren_part1_paper_compliant_anchored_axis_amplitude.py','gradient_jets','CompliantAnchoredAxisAmplitude',
             ["H=-4*z**3-z*z*self.params['j']+z*self.params['a']+self.params['j']",
              "((one-z*z*self.params['delta'])*H/IntervalTaylor(c,coefficients)).coefficients"]),
            ('anchored_parameters',base/'lei_ren_part1_paper_compliant_anchored_axis_amplitude.py','__init__','CompliantAnchoredAxisAmplitude',
             ["self.params=dict(j=self.core.j,delta=self.core.delta,sigma=self.core.sigma,a=(9-self.core.delta)/2)"]),
            ('same_core_gradient',src/'original_core_numeric.py','_seed','OriginalCoreNumeric',
             ["H=variable*(1-delta)/2+(one-variable*variable)*U","gradient=L*H/(H*H+(j/500)**2)",
              "[-v for v in gradient.coefficients]"]),
            ('same_core_sigma',base/'lei_ren_part1_paper_compliant_core_physical_field.py','__init__','CompliantCorePhysicalField',
             ["self.sigma=self.j/500"]),
        )
        guards={name:_function_ast_guard(path,fn,fragments,class_name=cls) for name,path,fn,cls,fragments in specs}
        for name,path,_,_,_ in specs:
            expected=self.input_hashes.get(path.name) if path.parent==base else self.implementation_hashes.get(path.relative_to(ROOT).as_posix())
            if not expected or sha(path)!=expected:raise ValueError('Unbound source AST: '+name)
        checks={'inherited_hydro_generator_and_Cauchy_guards':self.parent['source_guards']['passed'],
                'same_pressure_and_mass':self.pressure['source_family']==self.parent['source_family']==self.first['source_family'],
                'actual_positive_swirl_source':self.bridge['exact_swirl_scale_source']==
                    'F0(Z)^2=Cstar^-2*exp(-2Lambda*G(Z)); Gbar is a bound only',
                'same_formal_axis_source':self.domain['original_axis_amplitude']['source_formula']==
                    'exp(-selected_logCstar-Lambda*G(Z))',
                'saved_relative_C_base_amplitude_once':all(s['pressure']['original_amplitude_applied_once']
                                                         for s in self.first['samples']),
                'anchored_same_family':all(self.anchored[k]==v for k,v in self.parent['source_family'].items()),
                'anchored_same_selected_Cstar':self.anchored['same_selected_Cstar_and_current_source'] is True,
                'anchored_G_directed_source':self.anchored['actual_anchored_G_resolved_with_directed_error'] is True}
        return {'passed':all(v['passed'] for v in guards.values()) and all(checks.values()),'AST':guards,'checks':checks}

    def source_row(self,ztext,guards):
        from .original_swirl_source_kernel import dress_ordinary_jet,swirl_moment_mode_numerator_delta_sectors_ordinary
        c,iv=self.ctx,self.iv
        sample=next(s for s in self.first['samples'] if s['Z']==ztext and s['phase']=='0')
        def decode4(raw):
            if len(raw)<4:raise ValueError('Missing saved ordinary Z3 seed')
            pairs=[decode_pair(c,iv,row) for row in raw[:4]]
            return [a*math.factorial(k) for k,(a,b) in enumerate(pairs)],[b*math.factorial(k) for k,(a,b) in enumerate(pairs)]
        point={};core={}
        for name in ('Phi','B','C'):
            point[name],core[name]=decode4(sample['normalized_fields_and_histories'][name]['inlet'])
        cdp,cdt=decode4(sample['pressure']['C_inlet_dressed_by_relative_F0_squared'])
        rp=relative_amplitude_from_dressed_moment(cdp,point['C'],normalized_base_proved=guards['passed'])
        moment_rt=relative_amplitude_from_dressed_moment(cdt,core['C'],normalized_base_proved=guards['passed'])
        truth,transport=finite_width_swirl_comparison(c,iv,core,self.domain,self.weight,self.hcap,self.logh,
                                                     source_proved=guards['passed'])
        logeps,lei=decode_pair(c,iv,sample['source_scales']['log_epsilon'])
        anchor=self.anchored['anchored_amplitude_packets'][ztext]
        if not all(anchor[k] is True for k in ('exact_source_anchor_retained','root_parameter_uncertainty_propagated',
                                                'log_branches_certified','exp_logF0_not_materialized','positive_F0_retained_formally')):
            raise ValueError('Certified same-source anchored primitive required')
        gradient=[exact_iv(c,iv,v)*math.factorial(k) for k,v in enumerate(anchor['G_gradient_taylor_coefficients'][:3])]
        rt=relative_amplitude_from_gradient(gradient,iv.exp(-lei),source_proved=guards['passed'])
        logA=2*exact_iv(c,iv,anchor['logF0'])
        if mpf_cmp(logA._mpi_[0],self.logA._mpi_[0])<0 or mpf_cmp(logA._mpi_[1],self.logA._mpi_[1])>0:
            raise ValueError('Saved anchored amplitude lies outside the original global scale enclosure')
        index=('-.3','0','.3').index(ztext)
        qproof=self.pressure['source_jet']['comparison_q']['correlated_shared_phi_proof'][index]
        if not qproof['passed']:raise ValueError('Saved same-Phi quotient proof required')
        qp=[c.mpf(1),c.mpf(0),c.mpf(0)]
        qt=[iv.mpf(1 if k==0 else 0)+iv.mpf([-1,1])*iv.mpf(qproof['orders'][k]['representable_absolute_error_cap']) for k in range(3)]
        rhoproof=self.parent['hydro_source_rows'][index]['finite_width_comparison']['rho_error_proof'][0]
        if not rhoproof['passed']:raise ValueError('Saved native2 radius proof required')
        rho=iv.mpf(4)+iv.mpf([0,1])*iv.mpf(rhoproof['representable_absolute_error_cap'])
        results=[];packets=[]
        for ctx,values,ratio,logepsilon,radius,q in ((c,point,rp,logeps,c.mpf(4),qp),(iv,truth,rt,lei,rho,qt)):
            phi2=four_product(values['Phi'],values['Phi'])
            db=[a-b/2 for a,b in zip(values['B'],phi2)];dc=[a-b for a,b in zip(values['C'],phi2)]
            target=dress_ordinary_jet(ratio,phi2)
            bd=dress_ordinary_jet(ratio,db);cd=dress_ordinary_jet(ratio,dc)
            z=ctx.mpf(ztext);delta=ctx.exp(logepsilon+970)
            denominator=[2*(1-delta*z*z),-4*delta*z,-4*delta]
            modes=swirl_moment_mode_numerator_delta_sectors_ordinary(z,target,bd,cd)
            pieces={};numerators={}
            for j,sectors in enumerate(modes):
                for degree,numerator in sectors.items():
                    relative=ordinary_product(q,ordinary_quotient([-v for v in numerator],denominator))
                    scale=100**(2-j)*radius**j*ctx.exp((j-1+degree)*logepsilon+970*degree)
                    pieces[(j,degree)]=[v*scale for v in relative];numerators[(j,degree)]=numerator
            total=[sum((v[k] for v in pieces.values()),ctx.mpf(0)) for k in range(3)]
            results.append(total);packets.append((pieces,numerators))
        relative_point,relative_truth=results
        relative_truth,nonlinear=macro_nonlinear_enclosure(c,iv,relative_truth,self.weight,self.hcap,
                                                         macro_norm_bound_proved=guards['passed'])
        ap,at,scaleproof=scaled_atom_enclosure(c,iv,relative_truth,logA,source_proved=guards['passed'])
        pieces=[]
        for key in sorted(packets[0][0]):
            j,degree=key
            pieces.append({'mode_j':j,'delta_degree':degree,'epsilon_power':j-1+degree,
                           'source_numerator_ordinary_rows':[encode_pair(c,iv,c.mpf(a),iv.mpf(b)) for a,b in
                               zip(packets[0][1][key],packets[1][1][key])],
                           'relative_coefficient_ordinary_rows':[encode_pair(c,iv,a,b) for a,b in
                               zip(packets[0][0][key],packets[1][0][key])]})
        return {'Z':ztext,'ordinary_Z_orders':[0,1,2],
                'relative_amplitude_ordinary_Z0_Z3':[encode_pair(c,iv,a,b) for a,b in zip(rp,rt)],
                'relative_amplitude_approximation_source':'same-source dressed C/raw C with exact relative base one',
                'relative_amplitude_truth_source':'saved anchored G gradient Bell jet, Lambda=1/epsilon',
                'independent_moment_quotient_ordinary_Z0_Z3':[encode_pair(c,iv,a,b) for a,b in zip(rp,moment_rt)],
                'anchored_gradient_source_rows':anchor['G_gradient_taylor_coefficients'][:3],
                'relative_source_rows':[encode_pair(c,iv,a,b) for a,b in zip(relative_point,relative_truth)],
                'absolute_source_rows':[encode_pair(c,iv,a,b) for a,b in zip(ap,at)],
                'factored_actual_source':'exact positive F0(Zbase)^2 times relative_source_ordinary_jet',
                'absolute_source_atom':{'definition':self.bridge['exact_swirl_scale_source'],
                    'log_enclosure_exact_mpf_tuples':[list(v) for v in logA._mpi_],
                    'global_log_enclosure_exact_mpf_tuples':[list(v) for v in self.logA._mpi_],
                    'anchored_G_source_enclosure':anchor['G'],'selected_logCstar_source_enclosure':anchor['selected_logCstar'],
                    'anchored_primitive_definition':anchor['original_anchored_primitive'],
                    'anchored_packet_digest':digest(anchor),'pointwise_log_inside_global_enclosure':True,
                    'selected_value':None,'actual_amplitude_positive':True,'local_base_frozen':True,
                    'relative_derivatives_already_dressed':True,'Gbar_selected_as_G':False},
                'finite_width_comparison':transport,'amplitude_log_product_proof':scaleproof,
                'mode_delta_pieces':pieces,'macro_nonlinear_remainder':nonlinear,
                'source_quotient_proof_digest':digest(qproof),'native2_radius_proof_digest':digest(rhoproof),
                'all_three_radial_modes_retained':True,'external_amplitude_differentiated_twice':False,
                'approximation_is_exact_source_value':False,'actual_swirl_set_to_zero':False}

    def all_sector_rows(self,swirl):
        c,iv=self.ctx,self.iv
        pressure={r['Z']:r for r in self.pressure['source_jet']['rows']}
        hydro={r['Z']:r for r in self.parent['hydro_source_rows']}
        output=[]
        for row in swirl:
            groups=[[decode_pair(c,iv,v) for v in group] for group in
                    (pressure[row['Z']]['source_rows'],hydro[row['Z']]['source_rows'],row['absolute_source_rows'])]
            values=[encode_pair(c,iv,sum((g[k][0] for g in groups),c.mpf(0)),
                               sum((g[k][1] for g in groups),iv.mpf(0))) for k in range(3)]
            output.append({'Z':row['Z'],'ordinary_Z_orders':[0,1,2],'source_rows':values,
                           'components':['pressure','hydro','swirl'],'all_sectors_enclosed_at_query':True,
                           'swirl_positive_source_atom_retained':True,'weighted_primitive_complete':False})
        return output

    def integral_rows(self,swirl):
        c,iv=self.ctx,self.iv
        masses={c.mpf(v['phase']):v for v in self.mass['interior_mass']['rows']}
        sw={v['Z']:v for v in swirl};output=[]
        for old in self.parent['signed_source_integral_rows']:
            z,phase=old['Z'],old['phase'];sample=next(s for s in self.first['samples'] if s['Z']==z and s['phase']=='0')
            le,lei=decode_pair(c,iv,sample['source_scales']['log_epsilon'])
            delta,di=c.exp(le+970),iv.exp(lei+970);mass,mi=decode_pair(c,iv,masses[c.mpf(phase)])
            def coefficient(key):
                pairs=[decode_pair(c,iv,v) for v in sw[z][key]]
                p=main_coefficient(c.mpf(z),delta,[v[0] for v in pairs],mass,mass_Z_independent=True)
                t=main_coefficient(iv.mpf(z),di,[v[1] for v in pairs],mi,mass_Z_independent=True)
                return [encode_pair(c,iv,a,b) for a,b in zip(p,t)]
            factored=coefficient('relative_source_rows');absolute=coefficient('absolute_source_rows')
            combined=[]
            for known,extra in zip(old['pressure_plus_hydro_main_component_C0_C1'],absolute):
                kp,kt=decode_pair(c,iv,known);sp,st=decode_pair(c,iv,extra)
                combined.append(encode_pair(c,iv,kp+sp,kt+st))
            output.append({'Z':z,'phase':phase,'ordinary_Z_orders':[0,1],'raw_number':46,
                           'swirl_main_component_relative_C0_C1':factored,
                           'swirl_main_component_absolute_C0_C1':absolute,
                           'all_sector_main_component_C0_C1':combined,'physical_scale':old['physical_scale'],
                           'absolute_swirl_positive_atom':sw[z]['absolute_source_atom'],
                           'mass_Z_independent':True,'all_three_sectors_enclosed_at_query':True,
                           'raw47_weighted_remainder_retained':True,'complete_weighted_primitive':False})
        return output

    def build(self):
        guards=self.source_guards()
        if guards['passed'] is not True:raise ValueError('Swirl source guards failed: '+str(guards))
        swirl=[self.source_row(z,guards) for z in ('-.3','0','.3')]
        return {'source_family':self.parent['source_family'],'parent_report':PARENT,'parent_report_sha256':PARENT_SHA256,
                'source_guards':guards,'swirl_source_rows':swirl,'all_sector_source_rows':self.all_sector_rows(swirl),
                'signed_source_integral_rows':self.integral_rows(swirl),
                'inherited_m15_packet':self.parent['inherited_m15_packet'],
                'inherited_m15_packet_digest':self.parent['inherited_m15_packet_digest'],
                'inherited_mass_packet_digest':self.parent['inherited_mass_packet_digest'],
                'inherited_hydro_source_digest':digest(self.parent['hydro_source_rows']),
                'inherited_hydro_integrals_digest':digest(self.parent['signed_source_integral_rows']),
                'inherited_pressure_source_digest':self.parent['inherited_pressure_source_digest'],
                'scope':{'all_source_sectors_enclosed_at_three_Z_queries':True,
                         'all_sector_main_C0_C1_enclosed_at_19_queries':True,
                         'absolute_swirl_amplitude_selected':False,'complete_weighted_source_primitive':False},
                'contract':{'accepted_ancestor_replayed':False,'mass_recomputed':False,'actual_swirl_set_to_zero':False,
                            'positive_cap_or_Gbar_as_signed_value':False,'source_width_selected':False,
                            'R100_modes_reused_as_R0_modes':False,'raw213_reduced183_and_all_errors_retained':True},
                'limits':dict(self.parent['limits'],full_goal_complete=False)}
