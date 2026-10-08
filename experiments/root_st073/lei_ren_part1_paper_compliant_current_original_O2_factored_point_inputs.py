"""Complete original O2 finite point coefficients with source factors/errors.

R, Pstar, delta and L=1-delta*Z² remain correlated exact source factors.
True defining quadrature supplies finite radial point coefficients; directed
original integrals bound their error. The checked actual pressure coefficient
is reused as a Z-independent constant, not as a saved field-cover value.
This does not numerically expand native factors or install a phase oracle.
"""
from dataclasses import dataclass
from fractions import Fraction
import json
from pathlib import Path
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as s
import lei_ren_part1_paper_compliant_current_original_pressure_point_jets as pressure
from lei_ren_part1_paper_interval_outer_slope_field import transition_integrals
from lei_ren_part1_paper_interval_long_reshape_field import sigma_value_derivative
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

source=pressure.source;HERE,PREFIX,sha=source.HERE,source.PREFIX,source.sha
NAME=PREFIX+'current_original_O2_factored_point_inputs.json'
RECEIPT=PREFIX+'current_original_O2_factored_point_inputs_check.json'
GATE='original_O2_full_factored_point_inputs_with_directed_finite_coefficient_and_pressure_tail_errors'


@dataclass(frozen=True)
class FactoredPointTerm:
    # Powers of same original (R,Pstar,delta,L); L=1-delta*Z².
    factor_powers: tuple
    coefficient: object
    coefficient_error_upper: object
    late_pressure_error_log_upper: object=None

    def record(self):
        return dict(original_R_Pstar_delta_L_powers=list(self.factor_powers),
            approximate_source_point_coefficient=self.coefficient,
            directed_finite_coefficient_absolute_error_upper=self.coefficient_error_upper,
            late_pressure_error=dict(zero=self.late_pressure_error_log_upper is None,
                log_upper=self.late_pressure_error_log_upper))


@dataclass(frozen=True)
class FactoredPointRow:
    terms: tuple
    original_y: object
    original_Z: object

    def unscaled_scalar(self):
        if any(any(term.factor_powers) for term in self.terms):
            raise ValueError('Original correlated source factors cannot be cast to a native scalar')
        return sum((term.coefficient for term in self.terms),0)

    def record(self):
        return dict(terms=[term.record() for term in self.terms],source_factor_order=['R','Pstar','delta','L'],
            original_y=str(self.original_y),original_Z=str(self.original_Z),
            exact_L_definition='1-original_delta*original_Z²',
            error_contract='sum |R^r Pstar^p delta^d L^ell|*(coefficient_error+exp(late_pressure_error_log)); exact-zero late budget adds0',
            original_factors_not_materialized=True)


def finite_coefficient_templates(frame):
    t=frame.owner.template;z=t['z'];ps=t['Pstar'];delta=t['delta'];R=t['R'];L=1-delta*z*z
    f,H,D,P=s.symbols('f H D P',real=True)
    p0,p0Z,p0ZZ=s.symbols('P0 P0_Z P0_ZZ',real=True)
    replace={t['f']:f,t['H']:H,t['D']:D,t['P']:P,
        t['P0']:p0,s.diff(t['P0'],z):p0Z,s.diff(t['P0'],z,2):p0ZZ}
    inputs=(z,f,H,D,P,p0,p0Z,p0ZZ);rows={};identities={}
    for key in ('p1','p2'):
        for order in (0,1):
            original=t[key if not order else key+'_Z'].xreplace(replace)
            extra=ps if key=='p2' else 1
            numerator=s.cancel(original*extra*L**(order+1)/R)
            polynomial=s.Poly(numerator,ps,delta)
            terms=[]
            for (ps_power,delta_power),coefficient in polynomial.terms():
                assert ps_power in ((0,2) if key=='p2' else (0,))
                assert delta_power<=order+1
                powers=(1,ps_power-(1 if key=='p2' else 0),delta_power,-order-1)
                terms.append((powers,coefficient))
            rebuilt=sum((coefficient*R**powers[0]*ps**powers[1]*delta**powers[2]*L**powers[3]
                for powers,coefficient in terms),s.Integer(0))
            assert s.cancel(rebuilt-original)==0,(key,order)
            rows[(key,order)]=tuple(terms);identities[key+'_Z'+str(order)]=True
    for key,expression in (('E',f/(1+z*z)),('V',4*z),('a',s.Symbol('a',real=True)),('b',s.Integer(0))):
        if key=='a':continue
        for order in (0,1):
            value=s.diff(expression,z,order)
            rows[(key,order)]=() if value==0 else (((0,-1 if key=='V' else 0,0,0),value),)
    return dict(inputs=inputs,pressure_symbols=(p0,p0Z,p0ZZ),rows=rows,
        source_template_reconstruction_identities=identities,
        source_delta_and_L_factors_preserved_without_approximation=True)


class OriginalO2FactoredPointInputs:
    mode='original_O2_factored_point_coefficients_with_directed_error_budgets'
    def __init__(self,dps=50,cells=4096):
        if type(cells) is not int or cells<16:raise ValueError('At least 16 directed radial cells required')
        self.frame=source.OriginalO2SourceParameterFrame(dps);self.family=self.frame.family
        self.ctx=c=self.frame.owner.profiles.ctx;self.interval=iv=MPIntervalContext();iv.dps=dps+40
        accepted=json.loads((HERE/pressure.RECEIPT).read_bytes())
        if not accepted['all_passed'] or not accepted[pressure.GATE] or accepted['source_family']!=self.family:
            raise ValueError('Checked actual original normalized pressure point/error service required')
        self.hashes={**accepted['input_hashes'],pressure.RECEIPT:sha(pressure.RECEIPT)}
        for filename,digest in self.hashes.items():
            if sha(filename)!=digest:raise ValueError('Original factored-input dependency changed: '+filename)
        record=json.loads((HERE/pressure.NAME).read_bytes())
        self.pressure_proof=accepted['original_late_source_error_theorem']
        if not self.pressure_proof['passed'] or self.pressure_proof!=record['original_late_source_error_theorem']:
            raise ValueError('Same checked all-late-stage pressure error theorem required')
        # The accepted actual defining quadrature constant is independent of
        # Z/y/Md. Its exact saved point bits and directed enclosure are reused.
        self.alpha=c.make_mpf(tuple(record['actual_defining_quadrature_alpha']['exact_mpf_tuple']))
        bounds=record['directed_alpha_enclosure']
        lo=c.make_mpf(tuple(bounds['lower']['exact_mpf_tuple']));hi=c.make_mpf(tuple(bounds['upper']['exact_mpf_tuple']))
        self.alpha_enclosure=iv.mpf([lo,hi]);assert lo<=self.alpha<=hi
        self.cells=cells;self.radial_cache={};self.point_cache={}
        self.templates=finite_coefficient_templates(self.frame)
        self.compiled={};self.pressure_linear={}
        for key,terms in self.templates['rows'].items():
            self.compiled[key]=tuple((powers,
                s.lambdify(self.templates['inputs'],expr,modules=[{'mpf':c.mpf},'mpmath']),
                s.lambdify(self.templates['inputs'],expr,modules=[{'mpf':iv.mpf},'mpmath'])) for powers,expr in terms)
            self.pressure_linear[key]=tuple(tuple(s.lambdify(self.templates['inputs'],s.diff(expr,p0),
                modules=[{'mpf':iv.mpf},'mpmath']) for p0 in self.templates['pressure_symbols']) for powers,expr in terms)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)

    def radial(self,y):
        y=source.exact_coordinate(y)
        if y not in self.radial_cache:
            c=self.ctx;iv=self.interval
            point=self.frame.owner.profiles.radial(c.mpf(int(y.p))/int(y.q))
            with mp.workdps(iv.dps+40):
                J,mass=transition_integrals(iv,Fraction(int(y.p),int(y.q)),self.cells)
                iy=iv.mpf(int(y.p))/int(y.q)
                boxes=dict(f=iv.exp(iy/10-iv.mpf(3)/5*J),
                    H=(iv.mpf(5)/8+mass[0])*iv.exp(-3*iy/2),
                    D=(iv.mpf(5)/12+mass[2]/2)*iv.exp(-iy),P=iv.mpf(5)/2+mass[1]/2,
                    a=iv.mpf(4)/5+iv.mpf(6)/5*sigma_value_derivative(iv,iy)[0])
                assert endpoints(boxes['f'])[0]>0
            self.radial_cache[y]=(point,boxes)
        return self.radial_cache[y]

    def evaluate(self,*,y,Z):
        y=source.exact_coordinate(y);Z=pressure.exact_Z(Z);cache_key=(y,Z)
        if cache_key in self.point_cache:return self.point_cache[cache_key]
        c=self.ctx;iv=self.interval;radial,boxes=self.radial(y)
        with mp.workdps(iv.dps+40):
            z=c.mpf(int(Z.p))/int(Z.q);iz=iv.mpf(int(Z.p))/int(Z.q);q=1+z*z;iq=1+iz*iz
            factors=(q**-2,-4*z/q**3,(-4+20*z*z)/q**4)
            ifactors=(iq**-2,-4*iz/iq**3,(-4+20*iz*iz)/iq**4)
            pvalues=tuple(-self.alpha*factor for factor in factors)
            pboxes=tuple(-self.alpha_enclosure*factor for factor in ifactors)
            tail_logs=tuple(None if order==1 and Z==0 else
                endpoints(iv.mpf(3)/5-(iv.exp(40)+11)+iv.ln(bound)-iv.ln(2)-2*iv.ln(iq))[1]
                for order,bound in enumerate((5,10,44)))
            values=(z,*(radial[key] for key in ('f','H','D','P')),*pvalues)
            intervals=(iz,*(boxes[key] for key in ('f','H','D','P')),*pboxes)
            rows={}
            for key,terms in self.compiled.items():
                evaluated=[]
                for index,(powers,point_fn,interval_fn) in enumerate(terms):
                    coefficient=c.mpf(point_fn(*values));box=iv.mpf(interval_fn(*intervals))
                    coefficient_error=endpoints(abs(box-iv.mpf(coefficient)))[1]
                    logs=[]
                    for sensitivity,tail in zip(self.pressure_linear[key][index],tail_logs):
                        multiplier=endpoints(abs(iv.mpf(sensitivity(*intervals))))[1]
                        if multiplier and tail is not None:logs.append(endpoints(iv.ln(iv.mpf(multiplier))+iv.mpf(tail))[1])
                    late_log=None if not logs else endpoints(iv.mpf(max(logs))+iv.ln(len(logs)))[1]
                    if coefficient or coefficient_error or late_log is not None:
                        evaluated.append(FactoredPointTerm(powers,coefficient,coefficient_error,late_log))
                rows[key]=FactoredPointRow(tuple(evaluated),y,Z)
            a=c.mpf(radial['a']);aerror=endpoints(abs(boxes['a']-iv.mpf(a)))[1]
            rows[('a',0)]=FactoredPointRow((FactoredPointTerm((0,0,0,0),a,aerror),),y,Z)
            rows[('a',1)]=FactoredPointRow((),y,Z)
        result=dict(source_family=self.family,mode=self.mode,original_y_exact=str(y),original_Z_exact=str(Z),
            inputs={key:(rows[(key,0)],rows[(key,1)]) for key in ('E','V','a','b','p1','p2')},
            exact_source_factor_contract=dict(original_logR=s.srepr(self.frame.radius_log(y)),
                original_logPstar=s.srepr(self.frame.definitions['logPstar']),
                original_logdelta=s.srepr(self.frame.definitions['log_delta']),
                original_L='1-original_delta*Z²',source_delta_and_radius_Z_independent=True),
            directed_radial_and_pressure_coefficient_errors_installed=True,
            full_original_O2_I_over_F_point_coefficients_installed=True,
            conditioned_native_factor_and_phase_evaluation_installed=False,
            numerical_original_source_point_or_integral_oracle_installed=False)
        self.point_cache[cache_key]=result
        return result


def record_query(query):
    return {**query,'inputs':{name:[row.record() for row in pair] for name,pair in query['inputs'].items()}}


def run():
    began=time.monotonic();owner=OriginalO2FactoredPointInputs()
    samples=[record_query(owner.evaluate(y=y,Z=z)) for y,z in (('0','0'),('.53','-.6'),('.53','0'),('.53','.37'),('1','.7'))]
    report=dict(**{GATE:True},source_family=owner.family,mode=owner.mode,
        source_template_reconstruction_identities=owner.templates['source_template_reconstruction_identities'],
        original_full_factored_O2_point_queries=samples,directed_radial_cells=owner.cells,
        original_pressure_defining_coefficient_cache_reused_not_field_covers=True,
        original_delta_L_R_Pstar_correlations_and_positive_axial_sector_retained=True,
        finite_coefficient_and_late_pressure_errors_separately_factored=True,
        full_original_O2_I_over_F_point_coefficients_installed=True,
        conditioned_native_factor_and_phase_evaluation_installed=False,
        original_p1_p2_scalar_point_values_installed=False,
        numerical_original_source_point_or_integral_oracle_installed=False,
        actual_five_controls_installed=False,current_whole_N_selected=False,
        **dict.fromkeys(source.inertial.profiles.loop.OPEN,False),input_hashes=owner.hashes,
        execution_seconds=time.monotonic()-began,
        scope='Complete O2 E/V/a/b/p1/p2 finite point coefficients and ordinary Z rows with original correlated R/Pstar/delta/L factors, directed radial/pressure errors and nonzero late pressure budgets. Native factors/phase and all-chart numerical oracle are not evaluated; no installed controls, global N or corrected field is claimed.')
    (HERE/NAME).write_text(json.dumps(source.inertial.profiles.loop.encoded(report),indent=2)+'\n',encoding='utf8')
    print('Full original O2 factored input coefficients and directed error budgets connected',flush=True)
    return report


if __name__=='__main__':run()
