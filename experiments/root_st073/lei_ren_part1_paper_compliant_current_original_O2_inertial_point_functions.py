"""Full original O2 inertial input functions with realized pressure jets.

The O2 radial coefficients are approximate defining-integral evaluations.
R, Pstar and delta remain original source parameters. The pressure symbols
are bound to the exact fourteen-stage original preheat operator and its
prescribed raw waiting root. No interval pressure midpoint is selected.
"""
import json
import math
from pathlib import Path
from types import SimpleNamespace
import time
import sympy as s
import lei_ren_part1_paper_compliant_current_original_O2_slope_point_profiles as profiles
import lei_ren_part1_paper_compliant_current_pre_pulse_stress_operator as stress
import lei_ren_part1_paper_compliant_current_pressure_terminal_closure as pressure

HERE,PREFIX,sha=profiles.HERE,profiles.PREFIX,profiles.sha
NAME=PREFIX+'current_original_O2_inertial_point_functions.json'
RECEIPT=PREFIX+'current_original_O2_inertial_point_functions_check.json'
GATE='original_O2_full_inertial_point_function_templates_and_exact_P0_Z_jets_connected'


def full_inertial_template():
    """Replay the original full stress program on the original O2 functions."""
    y,z=s.symbols('original_y Z',real=True);R,Pstar=s.symbols('original_R original_Pstar',positive=True)
    delta=s.Symbol('original_delta',positive=True)
    f,H,D,P=(s.Function('original_O2_'+name)(y) for name in ('f','H','D','P'))
    P0=s.Function('original_normalized_P0')(z);C=1/(1+z*z)
    E=C*f;Vraw=4*z
    h=C*H;k=Vraw*h;energy=Vraw**2/Pstar**2-C*C*D;cumulative_p=C*C*P
    histories=dict(m=Vraw,h=h,k=k,e=energy,p=cumulative_p)
    rows=lambda value:[s.diff(value,y,j) for j in range(5)]
    product=lambda a,b:[sum(math.comb(j,i)*a[i]*b[j-i] for i in range(j+1)) for j in range(5)]
    shifted=lambda values,rate:[sum(math.comb(j,i)*rate**(j-i)*values[i] for i in range(j+1)) for j in range(4)]
    asts=stress.SourceAST()
    raw=asts.replay('current_pre_pulse_stress_operator','raw_pre_stress_rows',
        dict(axial_derivative=lambda value:s.diff(value,z),product_rows=product,shifted_rows=shifted))
    ctx=SimpleNamespace(mpf=lambda value:s.Rational(str(value)))
    output=raw(ctx,delta,z,rows(E),rows(Vraw),{name:rows(value) for name,value in histories.items()},rows(P0+cumulative_p))
    theta=sum(part['shape'][0] for name,part in output['theta'].items() if name!='variable_radial_shear')
    axial=sum(Pstar**part['mode'][1]*part['shape'][0] for name,part in output['axial'].items() if name!='axial_radial_shear')
    p1=s.factor(R*theta/E);p2=s.factor(R*axial/(Pstar*E))
    p1_Z=s.diff(p1,z);p2_Z=s.diff(p2,z)
    # Preserve the original inertial/shear distinction and common radius.
    a=s.Symbol('original_a',positive=True);b=s.Integer(0)
    return dict(y=y,z=z,R=R,Pstar=Pstar,delta=delta,f=f,H=H,D=D,P=P,P0=P0,
        E=E,V=Vraw/Pstar,histories=histories,p1=p1,p2=p2,p1_Z=p1_Z,p2_Z=p2_Z,
        a=a,b=b,source_inertial_theta=theta,source_inertial_axial=axial,
        original_program_AST_hashes=asts.hashes,
        pressure_derivative_orders_required=(0,1,2),
        source_conditions='original R>0,Pstar>0,0<delta<1,Z in[-1,1]; radius and scales Z-independent')


class OriginalO2InertialPointFunctions:
    def __init__(self,dps=50):
        self.profiles=profiles.OriginalO2SlopePointProfiles(dps);self.family=self.profiles.family
        checked=json.loads((HERE/profiles.RECEIPT).read_bytes())
        if not checked['all_passed'] or not checked[profiles.GATE] or checked['source_family']!=self.family:
            raise ValueError('Accepted original O2 point coefficient service required')
        self.hashes={**checked['input_hashes'],profiles.RECEIPT:sha(profiles.RECEIPT)}
        admitted=json.loads((HERE/pressure.RECEIPT).read_bytes())
        if not admitted['all_passed'] or not admitted[pressure.GATES[0]] or {key:admitted[key] for key in self.family}!=self.family:
            raise ValueError('Same-family original fourteen-stage pressure function identity required')
        proof=admitted['original_pressure_function_identification']
        if not proof['passed']:raise ValueError('Original pressure function identification must pass')
        for name,digest in admitted['input_hashes'].items():
            if sha(name)!=digest:raise ValueError('Original pressure-function source changed: '+name)
            self.hashes[name]=digest
        for name in (Path(pressure.__file__).name,):
            if admitted['input_hashes'].get(name)!=sha(name):raise ValueError('Original pressure/stress program differs from accepted source')
            self.hashes[name]=sha(name)
        self.hashes[pressure.RECEIPT]=sha(pressure.RECEIPT)
        self.pressure=pressure.ExactOriginalPreheatPressureOperator()
        self.template=full_inertial_template()
        allN=json.loads((HERE/profiles.ALLN_CHECK).read_bytes())
        for name,digest in self.template['original_program_AST_hashes'].items():
            if allN['input_hashes'].get(name)!=digest:raise ValueError('Original inertial replay source differs: '+name)
            self.hashes[name]=digest
        self.pressure_jet_cache={};self.hashes[Path(__file__).name]=sha(Path(__file__).name)

    def pressure_Z_jet(self,order):
        """Exact normalized P0 derivative; fixed original integration bounds."""
        if type(order) is not int or not 0<=order<=2:raise ValueError('P0 Z orders 0..2 are required here')
        if order not in self.pressure_jet_cache:
            p=self.pressure.partition;z=p['z'];t=p['t']
            expression=-sum((s.Integral(s.diff(density,z,order),(t,*p['domains'][name]))
                for name,density in self.pressure.original_densities.items()),s.Integer(0))
            self.pressure_jet_cache[order]=expression.subs(p['W'],self.pressure.raw_waiting_root).xreplace(
                {p['delta']:self.template['delta']})
        return self.pressure_jet_cache[order]

    def pressure_at(self,Z,order=0):
        Z=s.sympify(Z)
        if Z.is_number and (Z.is_real is not True or Z.is_finite is not True or not -1<=Z<=1):
            raise ValueError('Original pressure Z domain is [-1,1]')
        return self.pressure_Z_jet(order).subs(self.pressure.partition['z'],Z)

    def functions(self,y,*,resolve_pressure=True):
        """C1(Z) full inertial expressions at one original radial coordinate.

        Radial quadrature coefficients remain approximate. Exact original
        pressure integrals and scale parameters remain unmaterialized.
        """
        radial=self.profiles.radial(y);c=self.profiles.ctx;t=self.template;z=t['z']
        representative=lambda value:s.Rational(value._mpf_[1]*(-1 if value._mpf_[0] else 1))*s.Integer(2)**value._mpf_[2]
        subs={t[key]:representative(radial[key]) for key in ('f','H','D','P')}
        pairs={name:(t[name].subs(subs),t[name+'_Z'].subs(subs)) for name in ('p1','p2')}
        if resolve_pressure:
            jets={t['P0']:self.pressure_at(z),s.diff(t['P0'],z):self.pressure_at(z,1),
                  s.diff(t['P0'],z,2):self.pressure_at(z,2)}
            pairs={name:tuple(row.xreplace(jets) for row in pair) for name,pair in pairs.items()}
        return dict(source_family=self.family,original_y=radial['y'],Z_variable=z,
            p1=pairs['p1'],p2=pairs['p2'],E=t['E'].subs(subs),V=t['V'],
            a=s.Rational(4,5)+s.Rational(6,5)*representative(profiles.loop.flat_step(c,radial['y'])),b=s.Integer(0),
            original_R=t['R'],original_Pstar=t['Pstar'],original_delta=t['delta'],
            pressure_jet_derivatives_from_original_fourteen_stage_integral=resolve_pressure,
            original_raw_waiting_root_constraint_retained=True,
            pressure_and_inertial_delta_are_the_same_source_symbol=True,
            approximate_radial_coefficients=True,scalar_p1_p2_values_or_native_parameters_selected=False,
            numerical_error_certified=False)


def run():
    began=time.monotonic();owner=OriginalO2InertialPointFunctions();got=owner.functions('.53',resolve_pressure=False)
    encode=lambda expression:s.srepr(expression)
    record=dict(**{GATE:True},source_family=owner.family,
        exact_original_pressure_stage_count=len(owner.pressure.original_densities),
        exact_normalized_P0_Z_jet_expressions={str(k):encode(owner.pressure_Z_jet(k)) for k in range(3)},
        original_full_O2_inertial_template={key:encode(owner.template[key]) for key in ('p1','p1_Z','p2','p2_Z')},
        original_O2_point_inertial_functions=dict(y=profiles.loop.encoded(got['original_y']),
            p1=[encode(row) for row in got['p1']],p2=[encode(row) for row in got['p2']],
            pressure_symbols_bound_to_exact_integral_jet_expressions=True,radial_coefficients_approximate=True),
        original_full_inertial_program_replayed=True,shear_terms_excluded_from_I_over_F=True,
        original_same_R_factor_and_Pstar_powers_retained=True,
        pressure_and_inertial_delta_are_the_same_source_symbol=True,
        full_energy_pressure_and_nonlinear_meridional_terms_retained=True,
        pressure_P0_P0_Z_P0_ZZ_function_symbols_realized_by_original_integrals=True,
        original_datum_interval_jets_or_mass_caps_selected_as_values=False,
        original_p1_p2_expression_functions_connected=True,original_p1_p2_scalar_point_values_installed=False,
        native_source_scale_and_phase_point_oracle_installed=False,
        numerical_original_source_point_or_integral_oracle_installed=False,
        actual_five_controls_installed=False,certified_actual_fixed_point_tail_installed=False,
        actual_terminal_Z_function_closure_installed=False,current_whole_N_selected=False,
        **dict.fromkeys(profiles.loop.OPEN,False),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
        scope='Full original O2 I/F functions and Z rows from original stress program with exact fourteen-stage P0/P0_Z/P0_ZZ integral realization. Radial coefficients approximate; native scale/pressure/phase scalar values and certified numerical errors remain open.')
    (HERE/NAME).write_text(json.dumps(record,indent=2)+'\n',encoding='utf8')
    print('Original O2 full inertial functions and exact preheat pressure Z jets connected',flush=True)
    return record


if __name__=='__main__':run()
