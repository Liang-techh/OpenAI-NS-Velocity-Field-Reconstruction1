"""Full original O2 p2/Z carrier with pressure parity and varying profiles.

The original pressure remainder has an odd first jet, not an arbitrary
constant error divided by Z. No reference profile ratio or q lower cap is
borrowed. This certificate supplies source geometry, not terminal controls.
"""
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_O2_mixed_source_jets as source

base,prior,ep=source.base,source.prior,source.ep
HERE,PREFIX,sha=source.HERE,source.PREFIX,source.sha
C0,Y,Z,YZ=source.C0,source.Y,source.Z,source.YZ
NAME=PREFIX+'current_original_O2_axial_carrier.json.gz'
RECEIPT=PREFIX+'current_original_O2_axial_carrier_check.json'
GATE='original_O2_whole_axial_p2_carrier_and_pressure_parity_bound'


def original_pressure_odd_carrier(owner):
    inputs=owner.owner.inputs;operator=inputs.frame.owner.pressure;p=operator.partition
    z,Q=p['z'],p['q'];pressure=source.ordered.base.point.pressure
    beta2=set(pressure.source.inertial.pressure.BETA2);beta0=set(pressure.source.inertial.pressure.BETA0)
    identities={}
    for name,density in operator.original_densities.items():
        if name in beta2:factor=-4*z/Q
        elif name in beta0:factor=s.Integer(0)
        else:
            if name!='z_flatten':raise ValueError('Unsupported original pressure parity sector')
            factor=4*z*(p['sigma'](p['t']/100)-1)/Q
        if s.simplify(s.diff(density,z)-factor*density)!=0:
            raise ValueError('Original pressure first-jet source changed: '+name)
        identities[name]=dict(exact_derivative_factor=s.srepr(factor),bounded_by_four_abs_Z_over_Q=True)
    if not inputs.pressure_proof['passed'] or inputs.pressure_proof['total_remainder_factors']!=[5,10,44]:
        raise ValueError('Accepted full original pressure remainder theorem required')
    return dict(passed=True,source_family=owner.family,original_density_first_jet_identities=identities,
        original_datum_independent_of_source_y=True,original_remainder_even_and_first_jet_odd=True,
        smooth_odd_jet_carrier_H_R='R0_Z/Z, extended smoothly at Z=0; H_R is not R0_ZZ',
        error_factors_after_exp_three_fifths_over_Pstar=['(5/2)/Q^2','10*Z/Q^3','22/Q^2'],
        same_original_removed_axial_tail_and_all_late_stages_retained=True,
        original_late_pressure_source_proof=inputs.pressure_proof,
        derivation='All original density first derivatives are bounded by4*abs(Z)/Q times density, including removed finite axial tail. Apply the accepted order0 factor5 to the same defining integrals, not to a selected error cap.')


def even_rational_horner(expr,z,x):
    expr=s.cancel(expr);numerator,denominator=s.fraction(expr)
    def polynomial(value):
        poly=s.Poly(value,z);terms=[]
        for (power,),coefficient in poly.terms():
            if power%2:raise ValueError('Full original carrier must be even in Z')
            terms.append(coefficient*x**(power//2))
        return sum(terms,s.Integer(0))
    return s.horner(polynomial(numerator),x)/s.factor(polynomial(denominator))


class OriginalO2AxialCarrier:
    def __init__(self):
        self.source=source.OriginalO2MixedSources();self.owner=self.source.owner;self.family=self.source.family;self.ctx=self.source.ctx
        receipt=json.loads((HERE/source.RECEIPT).read_bytes())
        if not receipt.get('all_passed') or not receipt.get(source.GATE) or receipt['source_family']!=self.family:
            raise ValueError('Accepted full original O2 mixed source required')
        self.hashes=dict(self.source.hashes)
        for name,digest in {**receipt['input_hashes'],source.RECEIPT:sha(source.RECEIPT)}.items():
            if sha(name)!=digest:raise ValueError('Original O2 carrier dependency changed: '+name)
            if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Original O2 carrier closures disagree')
            self.hashes[name]=digest
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.pressure_proof=original_pressure_odd_carrier(self.source)
        templates=self.source.original;z,f,H,D,P,p0,p0Z,p0ZZ=templates['inputs']
        W=self.source.symbols[Y][-2];a=self.source.a_symbol;Q=1+z*z
        baseline={p0:-(P+W)/Q**2,p0Z:4*(P+W)*z/Q**3,p0ZZ:(4-20*z*z)*(P+W)/Q**4}
        errors=(s.Rational(5,2)/Q**2,10*z/Q**3,22/Q**2)
        x=s.Symbol('original_Z_squared',nonnegative=True);self.symbols=(f,H,D,P,W,a,x)
        self.compiled={};self.template_proof=[];self.expressions={}
        for direction,rows in ((C0,templates['rows']),(Y,self.source.derivative_templates)):
            compiled=[];expressions=[]
            for powers,expr in rows[('p2',0)]:
                pieces=[(None,expr.subs(baseline))]
                for order,(symbol,error) in enumerate(zip((p0,p0Z,p0ZZ),errors,strict=True)):
                    pieces.append((order,s.diff(expr,symbol)*error))
                for order,value in pieces:
                    if value==0:continue
                    quotient=s.cancel(value/z)
                    if s.simplify(quotient.subs(z,-z)-quotient)!=0:
                        raise ValueError('Original full p2 pressure carrier parity changed')
                    finite=even_rational_horner(quotient,z,x)
                    if s.cancel(finite.subs(x,z*z)*z-value)!=0:
                        raise ValueError('Original odd carrier factorization failed')
                    actual=powers if order is None else (powers[0],powers[1]-1,powers[2],powers[3])
                    compiled.append((actual,s.lambdify(self.symbols,finite,modules=[{'mpf':self.ctx.mpf},'mpmath']),order))
                    expressions.append((actual,finite,order))
                    self.template_proof.append(dict(ordinary_y_order=direction[0],original_factor_powers=powers,
                        pressure_error_order=order,actual_source_factor_powers=actual,
                        exact_odd_Z_factor_and_even_quotient=True,
                        carrier_has_finite_smooth_extension_at_axis=True))
            self.compiled[direction]=tuple(compiled);self.expressions[direction]=tuple(expressions)

    def cell(self,count,index,lower,upper):
        if count not in self.source.parent.parent.levels or type(index) is not int or not 0<=index<count:
            raise ValueError('Accepted whole original O2 mass cell required')
        level,_=self.source.parent.parent.levels[count];profiles=level['whole_source_cells'][index]['whole_original_radial_profile_covers']
        left,right=s.Rational(index,count),s.Rational(index+1,count);c=self.ctx
        eta=source.ordered.interval(c,self.owner.scales.logs['eta'])
        logq=source.ordered.hull(c,source.positive.endpoint_logq(c,right,eta),source.positive.endpoint_logq(c,left,eta))
        atlas=source.O2MixedAtlas(self.owner.inputs.frame,lower=lower,upper=upper,logq=logq);c=atlas.ctx
        with mp.workdps(c.dps+40):
            y=c.mpf((ep(atlas.rational(left))[0],ep(atlas.rational(right))[1]))
            z=c.mpf((ep(atlas.rational(s.Rational(lower)))[0],ep(atlas.rational(s.Rational(upper)))[1]))
            values=tuple(atlas.copy_interval(profiles[name]) for name in ('f','H','D','P','remaining_pressure_mass','a'))+(z*z,)
            result={}
            for direction,terms in self.compiled.items():
                contributions=[]
                for powers,fn,error_order in terms:
                    coefficient=c.mpf(fn(*values))
                    if error_order is not None:coefficient*=c.exp(c.mpf(3)/5)*c.mpf((-1,1))
                    contributions.append(source.axial.regular.finite_offset_anchor(atlas,atlas.term(powers,coefficient,coordinate=y)))
                value=atlas.sum(contributions)
                result[direction]=source.axial.normalize(atlas,value,(11,10,0,0,0))
            return dict(exact_y_cell=[str(left),str(right)],exact_abs_Z_cell=[str(lower),str(upper)],
                p2_over_Z_divided_by_Lambda0=result[C0],p2_y_over_Z_divided_by_Lambda0=result[Y],
                original_pressure_odd_first_jet_not_constant_divided_by_Z=True,
                actual_original_f_H_D_P_a_and_remaining_pressure_mass=True,
                full_original_late_pressure_errors_retained=True)

    def certificate(self,*,axial_cells=16):
        if type(axial_cells) is not int or axial_cells<4:raise ValueError('At least four whole axial cells required')
        count=min(self.source.parent.parent.levels);records=[]
        for index in range(count):
            for j in range(axial_cells):
                row=self.cell(count,index,s.Rational(j,axial_cells),s.Rational(j+1,axial_cells))
                if ep(row['p2_over_Z_divided_by_Lambda0'])[1]>=0:
                    raise ArithmeticError('Actual O2 full carrier sign needs refinement: '+str((index,j,ep(row['p2_over_Z_divided_by_Lambda0']))))
                records.append(row)
        c=self.ctx
        with mp.workdps(c.dps+40):
            g=source.axial.hull(c,[row['p2_over_Z_divided_by_Lambda0'] for row in records])
            gy=source.axial.hull(c,[row['p2_y_over_Z_divided_by_Lambda0'] for row in records])
            return dict(passed=True,source_family=self.family,original_y_window=['0','1'],original_Z_window=['-1','1'],
                sign=-1,Lambda0='Pstar^11*Cstar^10',original_mass_source_level=count,axial_cells_per_radial_cell=axial_cells,
                g_normalized=g,g_y_normalized=gy,g_y_over_g=gy/g,
                all_negative_Z_from_exact_even_carrier_and_y_derivative_parity=True,
                all_Z_zero_from_same_smooth_carrier_extension=True,whole_original_source_partition=records,
                source_identities=['p2=Z*g','p2_y=Z*g_y','u=p2*q/dstar','u_y/u=g_y/g+q_y/q (where u!=0)'],
                actual_q_y_not_set_to_zero=True,reference_profile_ratios_not_used=True,
                signed_p2_sign_opposite_physical_Z=True,
                regular_signed_predicate_frames_or_mixed_primitives_installed=False)


def run():
    began=time.monotonic();owner=OriginalO2AxialCarrier();certificate=owner.certificate()
    report=dict(**{GATE:True},source_family=owner.family,original_pressure_odd_carrier_theorem=owner.pressure_proof,
        exact_original_coefficient_carrier_templates=owner.template_proof,actual_O2_whole_axial_carrier_certificate=certificate,
        actual_nonconstant_original_O2_source_carrier_enclosed=True,
        regular_signed_predicate_frames_or_mixed_primitives_installed=False,
        actual_changed_five_integrals_installed=False,all_17_chart_or_24_cell_oracle_installed=False,
        actual_five_controls_installed=False,current_whole_N_selected=False,
        **dict.fromkeys(source.ordered.base.point.source.inertial.profiles.loop.OPEN,False),
        input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
        scope='Full original O2 p2/Z and p2_y/Z carrier over continuous y[0,1]/Z[-1,1], original pressure odd first-jet errors, actual varying profiles. Source sign/ratio certificate only, not mixed signed phase, five integrals, terminal/global N or full reconstruction.')
    (HERE/NAME).write_bytes(gzip.compress(json.dumps(base.encoded(report),indent=2).encode()+b'\n',compresslevel=9,mtime=0))
    print('Actual whole original O2 carrier:',len(certificate['whole_original_source_partition']),'source rectangles; g sign negative',flush=True)
    return report


if __name__=='__main__':run()
