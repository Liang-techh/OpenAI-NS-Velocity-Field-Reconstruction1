"""Actual generic-loop shear/inertial inputs without resolving source factors.

F=Utheta/sqrt(2R)>0 is the paper shear unit, not its scalar velocity
amplitude. I/F is recovered from all signed linear/quadratic inertial
sectors. Quotients keep their true denominators; broad coefficient covers
are never inverted or chosen as defining point values. Radius factors and
the original radius tree stay separate from the four fixed source bases.
"""
import gzip
import json
from dataclasses import dataclass
from pathlib import Path
import sympy as s
import lei_ren_part1_paper_compliant_current_generic_shear_source_packets as packets

HERE,PREFIX,sha=packets.HERE,packets.PREFIX,packets.sha
NAME=PREFIX+'current_generic_shear_inputs.json'
RECEIPT=PREFIX+'current_generic_shear_inputs_check.json'
VIEWS=PREFIX+'current_generic_shear_inputs_views.json.gz'
GATE='current_original_generic_shear_input_expressions_available'
OPEN=packets.OPEN


class RadiusPolynomial:
    """Finite signed R-power sum with unchanged four-factor axial jet coefficients."""
    def __init__(self,algebra,terms):
        if any(type(k) is not int or k<0 for k in terms):raise ValueError('Nonnegative integer original radius powers required')
        self.algebra=algebra;self.ctx=algebra.ctx
        self.terms={k:algebra.lift(v) for k,v in terms.items() if algebra.lift(v).terms}

    def coerce(self,value):
        if isinstance(value,RadiusPolynomial):
            if value.algebra is not self.algebra:raise ValueError('Original source factor bases differ')
            return value
        return RadiusPolynomial(self.algebra,{0:self.algebra.lift(value)})

    def __add__(self,value):
        other=self.coerce(value);terms=dict(self.terms)
        for power,row in other.terms.items():terms[power]=terms[power]+row if power in terms else row
        return RadiusPolynomial(self.algebra,terms)
    __radd__=__add__
    def __neg__(self):return RadiusPolynomial(self.algebra,{k:-v for k,v in self.terms.items()})
    def __sub__(self,value):return self+-self.coerce(value)
    def __rsub__(self,value):return self.coerce(value)+-self
    def __mul__(self,value):
        other=self.coerce(value);terms={}
        for a,left in self.terms.items():
            for b,right in other.terms.items():
                key=a+b;row=left*right
                terms[key]=terms[key]+row if key in terms else row
        return RadiusPolynomial(self.algebra,terms)
    __rmul__=__mul__
    def __pow__(self,n):
        if type(n) is not int or n<0:raise ValueError('Nonnegative integer source polynomial power required')
        out=self.coerce(1)
        for _ in range(n):out=out*self
        return out
    def record(self):return dict(original_radius_power_terms=[dict(radius_power=k,coefficient=v) for k,v in sorted(self.terms.items())],
        no_radius_or_source_factor_materialization=True,no_signed_sector_pruning=True)


@dataclass(frozen=True)
class SourceQuotient:
    numerator: RadiusPolynomial
    denominator: RadiusPolynomial
    denominator_condition: str
    def __post_init__(self):
        if self.numerator.algebra is not self.denominator.algebra:raise ValueError('Quotient source bases differ')
        if not self.denominator.terms:raise ValueError('Identically zero source denominator')
    def record(self):return dict(numerator=self.numerator.record(),denominator=self.denominator.record(),
        denominator_condition=self.denominator_condition,interval_denominator_inversion_performed=False,
        ratio_is_source_expression_not_point_value=True)


@dataclass(frozen=True)
class GenericSourceInputs:
    packet: packets.CurrentSourcePacket
    E: RadiusPolynomial
    C: RadiusPolynomial
    B: RadiusPolynomial
    inertial_theta: RadiusPolynomial
    inertial_axial: RadiusPolynomial
    denominator: RadiusPolynomial
    kappa_numerator: RadiusPolynomial
    kappa_excess_numerator: RadiusPolynomial
    stronger_numerator: RadiusPolynomial
    direction_numerator: RadiusPolynomial
    transverse_numerator: RadiusPolynomial

    def quotients(self):
        e='E=Utheta/Pstar>0 as an original source function; a>0 not inferred from coefficient caps'
        ce='C*E=a*E^2>0 requires the actual source a>0 proof'
        return dict(a=SourceQuotient(self.C,self.E,e),b=SourceQuotient(self.B,self.E,e),
            p1=SourceQuotient(self.inertial_theta,self.E,e),p2=SourceQuotient(self.inertial_axial,self.E,e),
            t0=SourceQuotient(-self.B,self.C,'C=a*E>0 requires the actual source a>0 proof'),
            kappa=SourceQuotient(self.kappa_numerator,self.denominator,ce),
            kappa_minus2=SourceQuotient(self.kappa_excess_numerator,self.denominator,ce),
            H0_minus2=SourceQuotient(self.stronger_numerator,self.denominator,ce),
            D=SourceQuotient(self.direction_numerator,self.denominator,ce),
            J=SourceQuotient(self.transverse_numerator,self.denominator,ce))

    def quadratic_numerator(self):
        """Deferred exact full signed quadratic, including radius/nonlinear terms."""
        return 2*self.direction_numerator**2*self.denominator-self.kappa_excess_numerator*self.transverse_numerator**2

    def record(self):
        return dict(chart=self.packet.chart,source_family=self.packet.source_family,
            original_source_provenance=self.packet.provenance,
            fixed_source_log_basis_names=packets.LOG_NAMES,fixed_source_log_bases=self.packet.algebra.logs,
            original_physical_radius_tree=self.packet.provenance['original_radius_source'],
            paper_shear_unit='F=Utheta/sqrt(2R)=S*E/sqrt(2R)',
            shear_over_F='(-a,b)',inertial_over_F='(p1,p2), before adding S/F',
            source_numerators=dict(E=self.E.record(),C=self.C.record(),B=self.B.record(),
                inertial_theta=self.inertial_theta.record(),inertial_axial=self.inertial_axial.record(),
                positive_denominator=self.denominator.record(),kappa=self.kappa_numerator.record(),
                kappa_minus2=self.kappa_excess_numerator.record(),H0_minus2=self.stronger_numerator.record(),
                D=self.direction_numerator.record(),J=self.transverse_numerator.record()),
            quotient_definitions={k:v.record() for k,v in self.quotients().items()},
            full_quadratic_recipe=dict(numerator='2*D_numerator^2*positive_denominator-kappa_minus2_numerator*J_numerator^2',
                denominator='positive_denominator^3',evaluated_without_point_factor_values=True,
                expansion_deferred_not_terms_discarded=True),
            full_relaxed_conditions=dict(positive_profile='E>0',positive_shear='C>0',
                kappa_le2='H0_minus2_numerator>0',
                kappa_gt2='D_numerator>0 and full_quadratic_numerator>0',
                generic_stronger_uniform_input='H0_minus2_numerator>0 also for kappa>2'),
            complete_signed_energy_mixed_pressure_and_meridional_terms_retained=True,
            P0_is_same_original_axis_datum=True,
            whole_box_quotient_bounds_certified=False,original_relaxed_cone_certified_by_this_packet=False,
            **{GATE:True},**dict.fromkeys(OPEN,False))


def from_packet(packet,delta):
    """Original actual source algebra; no phase or normalization is reapplied."""
    field=packet.recover_original(delta);algebra=packet.algebra
    constant=lambda row:RadiusPolynomial(algebra,{0:row})
    with_R=lambda row:RadiusPolynomial(algebra,{1:row})
    E=constant(packet.velocity['theta'][0])
    C=constant(packet.velocity['theta'][0]-2*packet.velocity['theta'][1])
    B=constant(2*packet.velocity['axial'][1])
    # sqrt(R/2)*S and sqrt(R/2)*S^2 divide by S*E/sqrt(2R).
    # The quotient is R*(linear+S*quadratic)/E in both components.
    It=with_R(field['inertial_theta_linear']+algebra.shift(field['inertial_theta_quadratic'],(0,.5,0,0)))
    Iz=with_R(field['inertial_axial_linear']+algebra.shift(field['inertial_axial_quadratic'],(0,.5,0,0)))
    den=C*E;kap=C*C+B*B;excess=kap-2*den
    Hnum=C*It-B*Iz;stronger=Hnum-2*den;D=Hnum-kap;J=C*Iz+B*It
    return GenericSourceInputs(packet,E,C,B,It,Iz,den,kap,excess,stronger,D,J)


def exact_theorem():
    E,C,B,R,S=s.symbols('E C B R S',positive=True)
    itl,itq,izl,izq=s.symbols('It_linear It_quadratic Iz_linear Iz_quadratic',real=True)
    # B is allowed signed; positivity is used only for E,C,R,S.
    B=s.Symbol('signed_B',real=True)
    a,b=C/E,B/E;F=S*E/s.sqrt(2*R)
    nt=R*(itl+S*itq);nz=R*(izl+S*izq);p1,p2=nt/E,nz/E
    den=C*E;kap=C*C+B*B;Hnum=C*nt-B*nz;Dnum=Hnum-kap;Jnum=C*nz+B*nt
    t0=-b/a;kappa=a+b*b/a;H0=p1+p2*t0;D=H0-kappa;J=p2-p1*t0
    checks={}
    def zero(name,left,right):
        if s.cancel(left-right)!=0:raise ArithmeticError('Actual generic input identity: '+name)
        checks[name]=True
    zero('same_full_theta_inertial_over_F',s.sqrt(R/2)*(S*itl+S*S*itq)/F,p1)
    zero('same_full_axial_inertial_over_F',s.sqrt(R/2)*(S*izl+S*S*izq)/F,p2)
    Ey,Vy=s.symbols('actual_E_y actual_V_y',real=True)
    zero('same_signed_angular_shear',((2*Ey-E)*S/s.sqrt(2*R))/F,-a.subs(C,E-2*Ey))
    zero('same_signed_axial_shear',(2*Vy*S/s.sqrt(2*R))/F,b.subs(B,2*Vy))
    zero('same_t0_signed_source',t0,-B/C)
    zero('same_kappa_division_free',kappa,kap/den)
    zero('same_H0_minus2_division_free',H0-2,(Hnum-2*den)/den)
    zero('same_signed_direction_division_free',D,Dnum/den)
    zero('same_transverse_division_free',J,Jnum/den)
    zero('same_complete_quadratic_division_free',2*D*D-(kappa-2)*J*J,
        (2*Dnum*Dnum*den-(kap-2*den)*Jnum*Jnum)/den**3)
    theta,axial=p1-a,p2+b
    zero('same_T_to_I_theta',theta+a,p1);zero('same_T_to_I_axial',axial-b,p2)
    zero('same_original_signed_D_projection',theta-b*axial/a,D)
    zero('same_original_signed_J_projection',axial+b*theta/a,J)
    zero('same_stronger_relaxed_H0_margin',theta-b*axial/a+kappa-2,H0-2)
    return dict(passed=True,identities=checks,
        paper_input_is_I_over_F_not_T_over_F=True,
        original_five_history_full_recovery_used=True,
        no_shear_only_or_angular_only_substitution=True,
        all_source_factors_and_exact_radius_retained=True,
        source_positivity_and_box_bounds_still_required=True)


class CurrentGenericShearInputs:
    def __init__(self,service=None):
        self.service=packets.CurrentSourcePackets() if service is None else service
        admitted=json.loads((HERE/packets.RECEIPT).read_bytes())
        if not admitted.get('all_passed') or not admitted.get(packets.GATE) or admitted.get('source_family')!=self.service.family:
            raise ValueError('Same checked actual source packet interface required')
        self.service.bind_hashes(admitted['input_hashes']);self.service.bind_hashes({packets.RECEIPT:sha(packets.RECEIPT)})
        self.service.bind_hashes({Path(__file__).name:sha(Path(__file__).name)})
    def saved(self,chart,view=None):return from_packet(self.service.saved(chart,view),self.service.data['delta'])
    def query(self,chart,Z,coordinate,log_tau='-1',theta=None,viscosity='1'):
        packet=self.service.query(chart,Z,coordinate,log_tau,theta,viscosity)
        return from_packet(packet,packet.algebra.ctx.mpf(self.service.data['delta']))


def run():
    owner=CurrentGenericShearInputs();records={};inventory={}
    for chart in packets.CHARTS:
        inputs=owner.saved(chart);records[chart]=inputs.record()
        inventory[chart]=dict(saved_view=inputs.packet.provenance.get('view'),
            quotient_inputs=list(inputs.quotients()),full_quadratic_deferred=True,
            source_positivity_and_whole_box_bounds_not_inferred=True)
        if inputs.packet.algebra.proofs or inputs.packet.algebra.final_rows:raise ValueError('Production source factors were resolved')
        print('Actual full generic shear inputs: '+chart,flush=True)
    (HERE/VIEWS).write_bytes(gzip.compress((json.dumps(packets.encode(records),separators=(',',':'))+'\n').encode(),compresslevel=6,mtime=0))
    owner.service.bind_hashes({VIEWS:sha(VIEWS)})
    result=dict(source_family=owner.service.family,current_original_chart_count=len(records),
        actual_source_input_inventory=inventory,full_input_cover_views=VIEWS,exact_full_source_input_theorem=exact_theorem(),
        **{GATE:True},**dict.fromkeys(OPEN,False),
        ancestor_constructors_called=False,source_factors_resolved=False,
        scope='Actual full inertial/shear source expressions and division-free original relaxed invariants for saved boxes. Not quotient norm bounds, complete relaxed admission, a loop installation or a new N.',
        input_hashes=owner.service.hashes)
    (HERE/NAME).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()
