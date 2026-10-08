"""One exact source parameter frame for O2 profiles, pressure and inertial inputs.

The selected logCstar is the builder's explicit dyadic parameter choice.
It is represented compactly rather than expanding its enormous integer.
Pressure/radius enclosure endpoints are never selected as field values.
This installs source definitions, not a conditioned numerical phase oracle.
"""
import ast
import hashlib
import json
from pathlib import Path
import time
import sympy as s
import lei_ren_part1_paper_compliant_current_original_O2_inertial_point_functions as inertial

HERE,PREFIX,sha=inertial.HERE,inertial.PREFIX,inertial.sha
NAME=PREFIX+'current_original_O2_source_parameter_frame.json'
RECEIPT=PREFIX+'current_original_O2_source_parameter_frame_check.json'
GATE='original_O2_correlated_source_parameter_frame_and_selected_radius_definition_connected'
PRESSURE_SOURCE=PREFIX+'pressure_source.json'
PHYSICAL_FAMILY=PREFIX+'physical_norm_family.json'


class ExactPositiveDyadic(s.Function):
    """Exact mantissa*2**exponent, with no integer or float expansion.

    The two integer arguments are its complete defining data. SymPy keeps
    this source constant opaque, including during substitution/evalf.
    """
    nargs=2
    is_positive=True
    is_real=True
    is_finite=True

    @classmethod
    def eval(cls,mantissa,exponent):
        if not isinstance(mantissa,s.Integer) or not isinstance(exponent,s.Integer) or mantissa<=0:
            raise ValueError('Positive exact integer dyadic source data required')
        return None

    def _eval_evalf(self,prec):
        return None


class ExactSourceExponential(s.Function):
    """The exact positive exp(argument), protected from automatic evalf.

    Native nested exponentials can exceed addressable integer storage.
    Keeping this explicit source operation prevents SymPy assumptions
    from attempting an implicit numerical expansion during substitution.
    """
    nargs=1
    is_positive=True
    is_real=True
    is_finite=True

    @classmethod
    def eval(cls,argument):
        if argument.is_real is False:raise ValueError('Real source exponent required')
        return None

    def _eval_evalf(self,prec):
        return None


def source_assignment(filename,method,target,wanted,cls=None):
    tree=ast.parse((HERE/filename).read_text(encoding='utf8'))
    body=next(n.body for n in tree.body if isinstance(n,ast.ClassDef) and n.name==cls) if cls else tree.body
    fn=next(n for n in body if isinstance(n,ast.FunctionDef) and n.name==method)
    expected=ast.dump(ast.parse(wanted,mode='eval').body)
    matches=[n for n in ast.walk(fn) if isinstance(n,ast.Assign)
        and any(ast.unparse(t)==target for t in n.targets) and ast.dump(n.value)==expected]
    if len(matches)!=1:raise ValueError('Original parameter definition changed: '+filename+':'+target)
    return dict(source=filename,method=method,target=target,expression=wanted,passed=True)


def exact_rational(value):
    if hasattr(value,'_mpf_'):
        sign,mantissa,exponent,bits=value._mpf_
        if bits<0:raise ValueError('Finite source coordinate required')
        if abs(exponent)>10000:raise ValueError('Finite O2 coordinate with bounded representation required')
        value=s.Integer((-1 if sign else 1)*mantissa)*s.Integer(2)**exponent
    else:
        try:value=s.Rational(value)
        except (ValueError,TypeError):raise ValueError('Finite rational source coordinate required') from None
    if value.is_real is not True or value.is_finite is not True:
        raise ValueError('Finite real source coordinate required')
    return value


def exact_coordinate(value):
    value=exact_rational(value)
    if not 0<=value<=1:
        raise ValueError('Exact original O2 coordinate y in[0,1] required')
    return value


class OriginalO2SourceParameterFrame:
    def __init__(self,dps=50):
        self.owner=inertial.OriginalO2InertialPointFunctions(dps);self.family=self.owner.family
        admitted=json.loads((HERE/inertial.RECEIPT).read_bytes())
        if not admitted['all_passed'] or not admitted[inertial.GATE] or admitted['source_family']!=self.family:
            raise ValueError('Checked same-family full original O2 inertial functions required')
        self.hashes=dict(admitted['input_hashes'])
        for name,digest in self.hashes.items():
            if sha(name)!=digest:raise ValueError('Accepted original parameter dependency differs: '+name)
        self.hashes[inertial.RECEIPT]=sha(inertial.RECEIPT)
        sources=(PRESSURE_SOURCE,PHYSICAL_FAMILY,PREFIX+'pressure_source.py',
            'lei_ren_part1_paper_logarithmic_outer_parameters.py',PREFIX+'physical_norm_family.py',
            PREFIX+'current_heat_physical_assembly.py')
        if any(self.hashes.get(name)!=sha(name) for name in sources):
            raise ValueError('Original definition and selected-parameter sources must match accepted hashes')
        source=json.loads((HERE/PRESSURE_SOURCE).read_bytes())['compliant_source']
        definition=source['implicit_source_definition']
        canonical=lambda obj:hashlib.sha256(json.dumps(obj,sort_keys=True).encode()).hexdigest()
        if canonical(definition)!=self.family['implicit_source_sha256'] or source['implicit_source_sha256']!=self.family['implicit_source_sha256']:
            raise ValueError('Exact pressure-source definition belongs to another family')
        if source['datum_enclosure_sha256']!=self.family['datum_enclosure_sha256']:
            raise ValueError('Original pressure datum family differs')
        expected=dict(Md='40',logPstar='exp(Md)+11',c_mu='.001',c_delta='.001',c_epsilon='.001',
            delta='min(1e-200,exp(-4logPstar-30))',Tw='-60log(mu)',Ts='4log(2/delta)',Tf=100,
            waiting='unique positive root of the continuous raw preheat waiting equation',
            angular_profile='Section 6.1 reference-plus-outer ansatz, H replaced by 1',
            cutoff_and_schedule_python_sha256=sha('lei_ren_part1_paper_outer.py'))
        if definition!=expected:raise ValueError('The current pinned original parameter definition differs')
        norm=json.loads((HERE/PHYSICAL_FAMILY).read_bytes());nd=norm['definition']
        if nd['pressure_source_sha256']!=self.family['implicit_source_sha256'] or nd['pressure_datum_sha256']!=self.family['datum_enclosure_sha256']:
            raise ValueError('Selected Cstar/radius uses another pressure family')
        if canonical(nd)!=norm['uniform_Cstar_family_sha256']:
            raise ValueError('Selected physical parameter definition hash differs')
        if nd['selection']!='logCstar=2*max(explicit logarithmic lower bounds)':
            raise ValueError('Explicit selected Cstar prescription required')
        chosen=norm['selected_logCstar']
        if chosen!=nd['selected_logC'] or chosen['lower_exact_mpf_tuple']!=chosen['upper_exact_mpf_tuple']:
            raise ValueError('An explicitly chosen singleton dyadic logCstar is required, not an interval point selection')
        sign,mantissa,exponent,bits=self.selected_logCstar_mpf_tuple=tuple(chosen['lower_exact_mpf_tuple'])
        if sign!=0 or type(mantissa) is not int or mantissa<=0 or type(exponent) is not int or mantissa.bit_length()!=bits:
            raise ValueError('Invalid selected logCstar dyadic definition')
        self.selected_physical_family_sha256=norm['uniform_Cstar_family_sha256']
        self.source_bindings=self.bind_definitions()
        Md=s.Integer(40);logP=s.exp(Md)+11;logmu=-s.log(1000)-4*logP
        # exp(40)>=1+40+40²/2=841, whereas log(10)<3 because
        # exp(3)>1+3+9/2+27/6=13. Hence 4logP+30>200log(10).
        self.delta_branch_proof=dict(passed=True,exp_Md_polynomial_lower=841,
            four_logP_plus30_lower=3438,two_hundred_log10_strict_upper=600,
            selected_branch='exp(-4logPstar-30)',no_rounded_branch_comparison=True)
        logdelta=-4*logP-30;delta=ExactSourceExponential(logdelta)
        epsilon=ExactSourceExponential(logdelta-s.log(1000))
        logC=ExactPositiveDyadic(s.Integer(mantissa),s.Integer(exponent))
        self.definitions=dict(Md=Md,logCstar=logC,logPstar=logP,
            Pstar=ExactSourceExponential(logP),log_mu=logmu,mu=ExactSourceExponential(logmu),log_delta=logdelta,delta=delta,
            log_epsilon=logdelta-s.log(1000),epsilon=epsilon,Td=s.exp(Md)+10,
            yd=logP,Tw=-60*logmu,Ts=4*(s.log(2)-logdelta),Tf=s.Integer(100),L=-30*logmu,
            logRref=s.log(110)+10*(logC+logP))
        self.definitions['Rref']=ExactSourceExponential(self.definitions['logRref'])
        p=self.owner.pressure.partition;t=self.owner.template
        self.substitutions={p[key]:self.definitions[key] for key in ('mu','delta','epsilon','yd','Tw','L','Ts')}
        self.substitutions.update({t['delta']:delta,t['Pstar']:self.definitions['Pstar']})
        self.definitions['W']=self.owner.pressure.raw_waiting_root.xreplace(self.substitutions)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)

    def bind_definitions(self):
        rows=[]
        filename='lei_ren_part1_paper_logarithmic_outer_parameters.py'
        for target,expression in (
            ('self.logPstar','c.exp(self.md)+11'),('self.log_mu',"c.ln(c.mpf('.001'))-4*self.logPstar"),
            ('choice','-4*self.logPstar-30'),('self.log_delta','choice'),
            ('self.Tw','-60*self.log_mu'),('self.Ts','4*(c.ln(2)-self.log_delta)'),('self.yd','c.exp(self.md)+11')):
            rows.append(source_assignment(filename,'__init__',target,expression,'LogarithmicOuterParameters'))
        rows.append(source_assignment(PREFIX+'pressure_source.py','__init__','self.log_epsilon',
            "c.ln(c.mpf('.001'))+self.log_delta",'CompliantOuterParameters'))
        rows.append(source_assignment(PREFIX+'physical_norm_family.py','run','selected_logC',
            'b.hi(2*max(endpoints(v)[1] for v in restrictions.values()))'))
        rows.append(source_assignment(PREFIX+'current_heat_physical_assembly.py','__init__','self.logRref',
            'c.ln(110)+10*(self.logC+self.logP)','CurrentHeatPhysicalAssembly'))
        outer='lei_ren_part1_paper_outer.py'
        for target,expression in (('self.y_rel','self.y_f-Decimal(30)*self.log_mu'),
            ('self.y_d','Decimal(1)+self.Td')):
            rows.append(source_assignment(outer,'__init__',target,expression,'PaperOuterSchedule'))
        return dict(passed=True,original_parameter_AST_assignments=rows,
            pressure_definition_canonical_hash_identified=True,selected_logCstar_definition_hash_identified=True,
            pressure_radius_and_waiting_enclosure_endpoints_not_used_as_values=True)

    def radius_log(self,y):
        return self.definitions['logRref']+exact_coordinate(y)

    def pressure_at(self,Z,order=0):
        return self.owner.pressure_at(Z,order).xreplace(self.substitutions)

    def functions(self,y):
        y=exact_coordinate(y)
        # Rational coordinates enter the approximate quadrature with all
        # available digits; radius uses the same exact requested coordinate.
        c=self.owner.profiles.ctx
        coefficient_y=c.mpf(int(y.p))/int(y.q)
        result=self.owner.functions(coefficient_y)
        mapping={**self.substitutions,self.owner.template['R']:ExactSourceExponential(self.radius_log(y))}
        for key in ('p1','p2'):result[key]=tuple(row.xreplace(mapping) for row in result[key])
        for key in ('E','V','a','b'):result[key]=result[key].xreplace(mapping)
        result.update(original_y_exact=y,original_R=mapping[self.owner.template['R']],
            original_Pstar=self.definitions['Pstar'],original_delta=self.definitions['delta'],
            original_logR=self.radius_log(y),selected_physical_family_sha256=self.selected_physical_family_sha256,
            exact_correlated_original_parameter_frame_installed=True,
            original_pressure_integrals_evaluated_numerically=False,
            conditioned_native_phase_or_scalar_oracle_installed=False,
            original_requested_coordinate_shared_by_radius_and_radial_quadrature=True)
        return result


def run():
    began=time.monotonic();frame=OriginalO2SourceParameterFrame();got=frame.functions('.53')
    record=dict(**{GATE:True},source_family=frame.family,
        selected_physical_family_sha256=frame.selected_physical_family_sha256,
        source_definition_binding=frame.source_bindings,delta_branch_proof=frame.delta_branch_proof,
        selected_logCstar_exact_mpf_tuple=list(frame.selected_logCstar_mpf_tuple),
        exact_source_parameter_definitions={key:s.srepr(value) for key,value in frame.definitions.items()},
        source_bound_O2_point_expressions=dict(y=str(got['original_y_exact']),
            p1=[s.srepr(row) for row in got['p1']],p2=[s.srepr(row) for row in got['p2']],
            original_logR=s.srepr(got['original_logR'])),
        native_radius_and_scale_source_definitions_selected=True,
        original_raw_waiting_root_definition_bound=True,
        astronomical_dyadic_integer_or_radius_expansion_performed=False,
        original_pressure_integrals_evaluated_numerically=False,
        conditioned_native_phase_or_scalar_oracle_installed=False,
        quadrature_and_roundoff_certified=False,
        original_p1_p2_scalar_point_values_installed=False,
        numerical_original_source_point_or_integral_oracle_installed=False,
        actual_five_controls_installed=False,current_whole_N_selected=False,
        **dict.fromkeys(inertial.profiles.loop.OPEN,False),input_hashes=frame.hashes,
        execution_seconds=time.monotonic()-began,
        scope='One original selected parameter/radius definition frame, bound to the existing source and selected dyadic logCstar family. Full O2 p1/p2 expression functions share these exact parameters and raw pressure integrals. Approximate radial quadrature, unevaluated pressure integrals and unconditioned phase remain; no native numerical oracle or completed field is claimed.')
    (HERE/NAME).write_text(json.dumps(record,indent=2)+'\n',encoding='utf8')
    print('Original O2 correlated parameter/radius source frame connected',flush=True)
    return record


if __name__=='__main__':run()
