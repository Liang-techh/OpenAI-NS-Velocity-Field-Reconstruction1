"""Whole original near-midplane source in a factored signed axial coordinate.

Physical Z=zeta/(Pstar^11*Cstar^10) is retained formally. Every Z power is
collected with the original radius/source powers before any bounded conversion.
Actual pressure jets enclose ONE original remainder; odd first-jet errors keep
their Z carrier. Original coefficient, inverse, Fourier and rate programs are
reused. Ordinary Z derivatives retain their native large factor.
"""
from dataclasses import dataclass
import gzip
import hashlib
import json
from pathlib import Path
from types import MappingProxyType
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_reference_midplane_integrals as midplane

whole=midplane.whole;points=whole.points;point=whole.point;base=whole.base;prior=whole.prior;ep=whole.ep
HERE,PREFIX,sha=whole.HERE,whole.PREFIX,whole.sha
NAME=PREFIX+'current_original_reference_near_midplane_integrals.json.gz'
RECEIPT=PREFIX+'current_original_reference_near_midplane_integrals_check.json'
GATE='original_whole_near_midplane_zeta_source_C0_Z_and_reference_integrals_enclosed'


class NearMidplaneAtlas(whole.factors.OriginalSourceFactorAtlas):
    """One original-family variable-L atlas for an exact finite zeta window."""
    def __init__(self,frame,*,zeta_lower,zeta_upper,dps=260):
        lo,hi=(point.source.exact_rational(v) for v in (zeta_lower,zeta_upper))
        if not -1<=lo<hi<=1:raise ValueError('Strict exact zeta window inside[-1,1] required')
        super().__init__(frame,Z=0,dps=dps)
        c=self.ctx
        with mp.workdps(c.dps+40):
            self.zeta_bounds=(lo,hi);self.Z=None
            zeta=c.mpf((ep(self.rational(lo))[0],ep(self.rational(hi))[1]))
            logLambda=11*self.bases[0]+10*self.bases[1]
            if ep(logLambda)[0]<=0:raise ArithmeticError('Original positive large axial carrier required')
            # Only bounded Q/L budgets are converted. Z itself, and all its
            # amplified occurrences in source rows, remain exact formal powers.
            z2=prior.ScaledEnclosure(prior.FormalScale(self.bases,(-22,-20,0,0,0)),zeta**2,self.ledger)
            dz2=self.parameter('delta')*z2
            Q=1+whole.conditioned.bounded(z2);L=1-whole.conditioned.bounded(dz2)
            if ep(L)[0]<=0 or ep(Q)[0]<1:raise ArithmeticError('Original Q/L positivity lost')
            self.bases=self.bases[:2]+(c.ln(L),)+self.bases[3:]
            self.Q=Q;self.L=L;self.logLambda=logLambda
            self.budget=dict(zeta_closed_interval=zeta,Q=Q,L=L,log_original_L=c.ln(L),
                physical_Z_squared_formal_powers=[-22,-20,0,0,0],
                delta_Z_squared_formal_powers=[-26,-20,0,0,0],
                original_Q_definition='1+(zeta/Lambda0)^2',
                original_L_definition='1-delta*(zeta/Lambda0)^2',
                L_and_Q_whole_zeta_hulls_not_exact_function_correlation=True,
                positive_small_tail_budgets_not_set_to_zero=True)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)

    def zterm(self,powers,coefficient,*,coordinate,zeta,Z_power):
        if type(Z_power) is not int or Z_power<0:raise ValueError('Nonnegative original Z polynomial power required')
        value=self.term(powers,coefficient*self.ctx.mpf(zeta)**Z_power,coordinate=coordinate)
        shift=prior.FormalScale(self.bases,(-11*Z_power,-10*Z_power,0,0,0))
        return prior.ScaledEnclosure(value.scale+shift,value.coefficient,self.ledger)

    def physical_Z(self,zeta):
        return prior.ScaledEnclosure(prior.FormalScale(self.bases,(-11,-10,0,0,0)),
                                     self.ctx.mpf(zeta),self.ledger)

    def rebase_piece_value(self,*args,**kwargs):
        raise ValueError('Fixed-Z point frames require an explicit near-midplane containment adapter')

    def record(self):
        return dict(source_family=self.family,coordinate_type='original_factored_axial_zeta_interval',
            exact_zeta_bounds=[str(v) for v in self.zeta_bounds],canonical_basis_order=self.order,
            defining_basis=self.bases,physical_Z='zeta/Lambda0',
            Lambda0='Pstar^11*Cstar^10',original_source_Z_power_shift='(-11*k,-10*k,0,0,0)',
            original_delta_identity='logdelta=-4logPstar-30',
            original_radius_identity='logR=log110+10logCstar+10logPstar+y',
            variable_L_positive_budget=self.budget,
            source_family_and_original_parameter_definitions_unchanged=True,
            amplified_Z_factors_collected_before_bounded_conversion=True,
            physical_Z_not_replaced_by_zero_or_native_float=True)


def pressure_odd_remainder_theorem(original_owner):
    """Original all-stage derivative factors sharpen only the odd jet cover."""
    pressure=point.pressure;operator=original_owner.inputs.frame.owner.pressure
    partition=operator.partition;z,Q=partition['z'],partition['q'];identities={}
    beta2=set(pressure.source.inertial.pressure.BETA2);beta0=set(pressure.source.inertial.pressure.BETA0)
    for name,density in operator.original_densities.items():
        if name in beta2:factor=-4*z/Q
        elif name in beta0:factor=s.Integer(0)
        else:
            if name!='z_flatten':raise ValueError('Unsupported original pressure parity sector')
            factor=4*z*(partition['sigma'](partition['t']/100)-1)/Q
        if s.simplify(s.diff(density,z)-factor*density)!=0:
            raise ValueError('Original pressure density derivative factor changed: '+name)
        identities[name]=dict(original_derivative_factor=s.srepr(factor),
            bound_by_four_abs_Z_over_Q=True)
    return dict(passed=True,original_stage_derivative_factor_identities=identities,
        bounds_apply_to_one_defined_original_even_remainder_and_its_jets=True,
        separate_jet_boxes_do_not_preserve_joint_arithmetic_correlation=True,
        normalized_pressure_identity='P0/Pstar^2=-alpha/Q^2+R0; (P0/Pstar^2)_Z=Z*(4*alpha/Q^3+H_R)',
        odd_carrier_definition='H_R=R0_Z/Z, with original smooth even extension at zero',
        source_R0_absolute_bound='5*exp(3/5)/(2*Pstar*Q^2)',
        source_H_R_absolute_bound='10*exp(3/5)/(Pstar*Q^3)',
        source_R0_ZZ_absolute_bound='44*exp(3/5)/(2*Pstar*Q^2)',
        H_R_not_identified_with_R0_ZZ=True,
        finite_alpha_quadrature_error_in_same_defining_coefficient_enclosure=True,
        derivation='Each original density derivative factor has magnitude<=4*abs(Z)/Q, including removed finite axial tail; apply accepted order0 factor5 to that same source integral. No derivative of an independently selected bound.',
        original_midplane_parity=original_owner.parity,
        original_all_late_error_theorem=original_owner.inputs.pressure_proof)


@dataclass(frozen=True)
class OriginalNearMidplaneFrame:
    owner:object
    family:object
    left:object
    right:object
    zeta_bounds:tuple
    roots:object
    kernel:object
    record:dict


class OriginalReferenceNearMidplane(whole.OriginalReferenceWholeCells):
    def __init__(self,*,zeta_lower='-1/100000000',zeta_upper='1/100000000'):
        self.midplane=midplane.OriginalReferenceMidplaneWholeCells()
        self.dispatcher=self.midplane.dispatcher;self.family=self.midplane.family
        self.inputs=self.midplane.inputs;self.scales=self.midplane.scales;self.templates=self.inputs.templates
        self.atlas=a=NearMidplaneAtlas(self.inputs.frame,zeta_lower=zeta_lower,zeta_upper=zeta_upper)
        self.ctx=c=a.ctx;self.Z=None;self.cache={};self.frames={};self.source_records={}
        self.pressure_theorem=pressure_odd_remainder_theorem(self.midplane)
        receipt=json.loads((HERE/midplane.RECEIPT).read_bytes())
        if not receipt.get('all_passed') or not receipt.get(midplane.GATE) or receipt['source_family']!=self.family:
            raise ValueError('Accepted original full-pressure midplane source required')
        self.hashes={**self.midplane.hashes,**a.hashes}
        for name,digest in {**receipt['input_hashes'],midplane.RECEIPT:sha(midplane.RECEIPT)}.items():
            if sha(name)!=digest:raise ValueError('Original near-midplane dependency changed: '+name)
            if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Original near-midplane source families differ')
            self.hashes[name]=digest
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.z=self.templates['inputs'][0];self.Qsymbol=s.Symbol('original_Q',positive=True)
        self.alpha=s.Symbol('original_alpha',real=True);self.terms={};self.compiler_proof={}
        self.arguments=(*self.templates['inputs'][1:5],self.alpha,self.Qsymbol)
        p0,p1,p2=self.templates['pressure_symbols'];z=self.z;Q=self.Qsymbol
        baseline={p0:-self.alpha/Q**2,p1:4*self.alpha*z/Q**3,p2:self.alpha*(4-20*z*z)/Q**4}
        error_factors=(s.Rational(5,2)/Q**2,10*z/Q**3,22/Q**2)
        for key,rows in self.templates['rows'].items():
            compiled=[];proofs=[]
            for powers,expr in rows:
                for P in self.templates['pressure_symbols']:
                    for S in self.templates['pressure_symbols']:
                        if s.diff(expr,P,S)!=0:raise ValueError('Original pressure-affine template required')
                pieces=[(None,expr.xreplace(baseline))]
                pieces.extend((order,s.diff(expr,P)*factor)
                    for order,(P,factor) in enumerate(zip(self.templates['pressure_symbols'],error_factors,strict=True)))
                for error_order,value in pieces:
                    if value==0:continue
                    terms,proof=self.Z_polynomial(value)
                    actual_powers=powers if error_order is None else (powers[0],powers[1]-1,powers[2],powers[3])
                    for k,coefficient in terms:
                        fn=s.lambdify(self.arguments,coefficient,modules=[{'mpf':c.mpf},'mpmath'])
                        compiled.append(dict(powers=actual_powers,k=k,fn=fn,expr=coefficient,error_order=error_order))
                    proofs.append(dict(original_factor_powers=powers,pressure_error_order=error_order,**proof))
            self.terms[key]=tuple(compiled);self.compiler_proof[str(key)]=proofs

    def Z_polynomial(self,value):
        # Denominators become bounded Q functions; signed Z powers remain
        # symbolic until they have cancelled the native source carrier.
        numerator,denominator=s.fraction(s.factor(s.together(value)))
        denominator=s.factor(denominator).subs(1+self.z**2,self.Qsymbol)
        if denominator.has(self.z):raise ValueError('Original source denominator needs a non-polynomial Z route')
        polynomial=s.Poly(s.expand(numerator),self.z)
        terms=tuple((int(k[0]),coefficient/denominator) for k,coefficient in polynomial.terms())
        rebuilt=sum(coefficient*self.z**k for k,coefficient in terms)
        if s.cancel((rebuilt-value).subs(self.Qsymbol,1+self.z**2))!=0:
            raise ValueError('Original Z polynomial source reconstruction failed')
        return terms,dict(reconstruction_after_Q_equals_one_plus_Z_squared_exact=True,
            original_Z_powers=[k for k,coefficient in terms],
            denominator_Z_independent_after_original_Q_binding=True)

    def source_frame(self,left,right):
        with mp.workdps(self.ctx.dps+40):
            left=points.reference.reference_coordinate(left);right=points.reference.reference_coordinate(right)
            if not left<right:raise ValueError('Strict original radial source cell required')
            a=self.atlas;c=self.ctx;y=c.mpf((ep(a.rational(left))[0],ep(a.rational(right))[1]))
            zl,zh=a.zeta_bounds;zeta=c.mpf((ep(a.rational(zl))[0],ep(a.rational(zh))[1]))
            f=c.exp(y/10);values=(f,c.mpf(5)/8*f,c.mpf(5)/12*f*f,c.mpf(5)/2*f*f,
                a.copy_interval(self.inputs.alpha_enclosure),a.Q)
            rows={};records=[]
            for key,terms in self.terms.items():
                result=a.scalar(0);source_terms=[]
                for term in terms:
                    coefficient=c.mpf(term['fn'](*values))
                    if term['error_order'] is not None:coefficient=coefficient*c.exp(c.mpf(3)/5)*c.mpf((-1,1))
                    value=a.zterm(term['powers'],coefficient,coordinate=y,zeta=zeta,Z_power=term['k'])
                    result=a.add(result,value)
                    source_terms.append(dict(original_factor_powers=term['powers'],original_Z_power=term['k'],
                        finite_coefficient_enclosure=coefficient,pressure_remainder_order=term['error_order'],
                        exact_collected_source_value=value.record()))
                rows[key]=result;records.append(dict(input=key[0],ordinary_Z_order=key[1],terms=source_terms))
            roots={name:{(0,k):rows[(name,k)] for k in (0,1)} for name in ('E','V','b','p1','p2')}
            roots['a']={(0,0):a.scalar(c.mpf(4)/5),(0,1):a.scalar(0)}
            roots['t0']={(0,0):a.scalar(0),(0,1):a.scalar(0)}
            eta=a.copy_interval(self.scales.logs['eta']);logamin=a.copy_interval(self.scales.logs['a_min'])
            q=base.current.q_enclosure(roots['a'][(0,0)],roots['a'][(0,0)]-2,eta,logamin)['q']
            dstar=a.copy_interval(self.scales.logs['d_star'])
            d=prior.ScaledEnclosure(prior.FormalScale(a.bases,offset=dstar),1,a.ledger)
            u=(roots['p2'][(0,0)]*q).positive_divide(d,dstar)
            ubox=whole.conditioned.bounded(u)
            if max(abs(v) for v in ep(ubox))>c.mpf('.25'):
                raise ArithmeticError('Actual whole zeta source exceeds original regular-u guard')
            kernel=whole.conditioned.PositiveLogQPhase(dict(q=q,roots=roots,original_u_source=u),dstar)
            if kernel.geometry!='small_r_series':raise ArithmeticError('Original whole regular branch required')
            record=dict(source_family=self.family,exact_reference_cell=[str(left),str(right)],
                exact_zeta_interval=[str(v) for v in a.zeta_bounds],physical_Z_map=a.physical_Z(zeta).record(),
                original_whole_cell_terms=records,actual_original_u_range=ubox,
                actual_original_roots={name:{str(k):v.record() for k,v in row.items()} for name,row in roots.items()},
                original_pressure_remainder_carrier=self.pressure_theorem,
                original_source_compiler_identities=self.compiler_proof,
                bounded_Q_and_positive_variable_L=a.budget,
                whole_y_and_zeta_source_enclosures_not_samples=True,
                ordinary_Z_jets_from_original_templates_not_derivatives_of_caps=True,
                physical_Z_not_zeroed_before_source_amplification=True)
            frame=OriginalNearMidplaneFrame(self,self.family,left,right,a.zeta_bounds,
                MappingProxyType({name:MappingProxyType(row) for name,row in roots.items()}),kernel,record)
            self.frames[id(frame)]=frame;return frame

    def primitives(self,frame,coordinate):
        if type(frame) is not OriginalNearMidplaneFrame or self.frames.get(id(frame)) is not frame or frame.owner is not self:
            raise ValueError('Issued original near-midplane frame required')
        with mp.workdps(self.ctx.dps+40):
            primitive=frame.kernel.primitives(coordinate,'psi')
            slow,proof=points.slow.slow_values(frame.kernel,frame.roots,coordinate,'psi')
            return dict(C0=primitive,Z=slow,proof=proof)

    def roots(self,left,right):
        frame=self.source_frame(left,right)
        self.source_records[(frame.left,frame.right)]=frame.record
        return dict(roots=frame.roots,kernel=frame.kernel,source_terms=frame.record['original_whole_cell_terms'],
            coordinate=self.ctx.mpf((ep(self.atlas.rational(frame.left))[0],ep(self.atlas.rational(frame.right))[1])),
            left=frame.left,right=frame.right,frame=frame)

    def cell(self,left,right,*,N,bits=32):
        with mp.workdps(self.ctx.dps+40):
            result=super().cell(left,right,N=N,bits=bits)
            source=self.source_records[tuple(map(points.reference.reference_coordinate,result['record']['exact_reference_cell']))]
            result['record'].update(exact_zeta_interval=[str(v) for v in self.atlas.zeta_bounds],
                original_physical_Z_map='Z=zeta/(Pstar^11*Cstar^10)',
                whole_near_midplane_source_and_C0_Z_primitives=True,
                original_pressure_odd_carrier_and_even_errors_retained=True,
                variable_L_and_Q_positive_budgets_retained=True,
                native_large_Z_derivative_factor_retained=True,
                whole_zeta_actual_original_u_range=source['actual_original_u_range'],
                whole_zeta_physical_Z_map=source['physical_Z_map'],
                whole_zeta_positive_Q_variable_L=source['bounded_Q_and_positive_variable_L'])
            return result

    def integrate(self,*,count,N):
        with mp.workdps(self.ctx.dps+40):
            result=super().integrate(count=count,N=N);report=result['report']
            report.pop('original_Z_exact')
            report.update(exact_zeta_interval=[str(v) for v in self.atlas.zeta_bounds],
                original_physical_Z_map='Z=zeta/(Pstar^11*Cstar^10)',
                whole_near_midplane_reference_C0_Z_contributions_installed=True,
                original_regular_source_and_Fourier_Z_functions_used=True,
                signed_implicit_Z_products_bounded_before_multiplication=False,
                native_large_Z_derivative_factor_not_capped=True,
                actual_source_functions_uniform_over_this_zeta_window=True,
                high_precision_source_and_integral_exports_retained=True,
                whole_Z_terminal_or_all_route_closure=False)
            return result


def run():
    begin=time.monotonic();owner=OriginalReferenceNearMidplane();levels=[]
    for count,N in ((4,160),(16,160),(16,16384)):
        levels.append(owner.integrate(count=count,N=N)['report'])
        print('Original whole near-midplane source/integrals:',count,'cells, N',N,flush=True)
    wider=OriginalReferenceNearMidplane(zeta_lower='-1/1000000',zeta_upper='1/1000000')
    if wider.family!=owner.family or wider.hashes!=owner.hashes:
        raise ValueError('Expanded neighborhood must use the identical original family and inputs')
    for count,N in ((16,160),(16,16384)):
        levels.append(wider.integrate(count=count,N=N)['report'])
        print('Original hundredfold expanded near-midplane source/integrals:',count,'cells, N',N,flush=True)
    report=dict(**{GATE:True},source_family=owner.family,
        actual_original_near_midplane_reference_levels=levels,
        actual_variable_L_source_atlas=owner.atlas.record(),
        original_pressure_odd_remainder_theorem=owner.pressure_theorem,
        original_Z_polynomial_compiler_identities=owner.compiler_proof,
        full_closed_zeta_neighborhood_not_midplane_samples=True,
        accepted_signed_zeta_windows=[['-1/100000000','1/100000000'],['-1/1000000','1/1000000']],
        wider_neighborhood_half_width_ratio=100,
        actual_C0_Z_source_and_regular_primitives_installed=True,
        original_large_Z_amplification_and_even_pressure_errors_retained=True,
        all_17_chart_or_24_cell_oracle_installed=False,actual_five_controls_installed=False,current_whole_N_selected=False,
        **dict.fromkeys(point.source.inertial.profiles.loop.OPEN,False),
        input_hashes=owner.hashes,execution_seconds=time.monotonic()-begin,
        scope='Whole original Rh_reference[-5,0] C0/Z source and regular primitive/integral enclosures on signed zeta[-1e-8,1e-8] and hundredfold expanded[-1e-6,1e-6], physical Z=zeta/(Pstar^11*Cstar^10). Full-pressure odd derivative carrier, even jet errors, finite alpha enclosure, variable positive L/Q hulls, exact source/Z power collection and actual N coefficients. Ordinary Z native factor remains. Separate pressure jet boxes retain provenance but not joint arithmetic correlation. Not full[-1,1] Z, terminal/all-route control/global N/recursive corrected NS.')
    (HERE/NAME).write_bytes(gzip.compress(json.dumps(base.encoded(report),indent=2).encode()+b'\n',compresslevel=9,mtime=0))
    return report


if __name__=='__main__':run()
