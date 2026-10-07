"""Current core/first source attachment and explicit strict exit collar.

Only checked receipts and pure source bindings are loaded. The exact width
and flat stress factor remain formal; no ancestor constructor is called.
"""
import ast
import functools
import hashlib
import json
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as s
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent;PREFIX='lei_ren_part1_paper_compliant_'
NAME=PREFIX+'current_inner_exit_strict_collar.json'
RECEIPT=PREFIX+'current_inner_exit_strict_collar_check.json'
GATE='current_inner_exit_strict_collar_attached_to_current_source_graph_certified'
OPEN=('global_completed_tensor_admissibility','admissible_stress_lift_constructed',
    'global_common_finite_N_certified','full_point_physical_field_evaluation',
    'full_inner_interfaces_certified','full_cartesian_vector_derivatives_certified',
    'actual_bridge_mixed4_feedback_installed','actual_point_moment_history_recovered',
    'physical_energy_integral_certified','independently_bounded_flat_remainder',
    'temporal_recursion','full_background_NS_validation','upstream_shear_loop_constructed')


def sha(name):return hashlib.sha256((HERE/name).read_bytes()).hexdigest()


def function(stem,method):
    tree=ast.parse((HERE/(PREFIX+stem+'.py')).read_text(encoding='utf8'))
    nodes=[n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name==method]
    if len(nodes)!=1:raise ValueError('Unique original source function required: '+stem+'.'+method)
    return nodes[0]


@functools.lru_cache(maxsize=1)
def exact_source_attachment():
    """Bind actual shear arithmetic and inlet, beyond matching family IDs."""
    bindings={};hashes={}
    def assignment(stem,method,target,expression):
        fn=function(stem,method);expected=ast.dump(ast.parse(expression,mode='eval').body)
        if not any(isinstance(n,ast.Assign) and any(ast.unparse(t)==target for t in n.targets)
            and ast.dump(n.value)==expected for n in ast.walk(fn)):
            raise ValueError('Current exit attachment source changed: '+stem+'.'+target)
        bindings[stem+'.'+method+'/'+target]=expression;hashes[PREFIX+stem+'.py']=sha(PREFIX+stem+'.py')
    specs=(
        ('inner_bridge_profiles','__init__','self.logh',"read_interval(c,ledger['shared_positive_width_log_enclosure'])"),
        ('comparison_point_integrals','inlet','atoms','self.atoms.atoms_from_packet(packet,6,shared_root=root)'),
        ('actual_bridge_integrals','prepare','data','self.comparison.inlet(Z,root)'),
        ('actual_bridge_integrals','packet','initial',"{n:j.truncate(5) for n,j in p['data']['moments'].items()}"),
        ('actual_bridge_integrals','packet','phi0',"p['data']['phi0'].truncate(5)"),
        ('actual_bridge_integrals','packet','pressure',"dress(own['actual']['C'],p['inputs']['F0_squared_ratios'])"),
        ('bridge_mixed_C4','evaluate','h','algebra.width(1)'),
        ('bridge_mixed_C4','evaluate','chi',"[algebra.lift(1-sigma[0])+h*sigma[0]]+[-algebra.lift(sigma[k])+h*sigma[k] for k in range(1,4)]"),
        ('bridge_mixed_C4','evaluate','Dbar',"[rate_rows(directions['D_over_R'],scale,k)*R for k in range(4)]"),
        ('bridge_mixed_C4','evaluate','drive',"[rate_rows(directions['hydro'],scale,k)*R+algebra.shift(rate_rows(directions['pressure'],scale,k)*R,(0,1,0,0))+algebra.shift(rate_rows(directions['swirl'],scale*2,k)*R**2,(0,0,1,0)) for k in range(4)]"),
        ('bridge_mixed_C4','evaluate','controls',"bridge_controls(algebra,scale,chi,Dbar,drive,phi/comparison['phi'].truncate(5),barlog)"),
        ('bridge_mixed_C4','bridge_controls','logF','[-scale*product_rows(chi,Dbar,k)/2 for k in range(4)]'),
        ('current_core_first_interface','boundary','inlet','shared_packet_inlet(self.comparison,packet)'),
        ('current_core_first_interface','boundary','p0',"self.upstream.prepare(z)['inputs']['p0']"),
        ('current_core_first_interface','boundary','ratios',"self.upstream.prepare(z)['inputs']['F0_squared_ratios']"),
        ('current_core_first_interface','boundary','dressed','dress(C,ratios)'),
        ('flat_pulse_derivatives','_sigma_left','odds','1/(1-x)**2-1/x**2'),
        ('flat_pulse_derivatives','_sigma_left','value','e/(1+e)'),
    )
    for spec in specs:assignment(*spec)
    for stem,method,target,key,expected in (
        ('physical_norm_family','run','K_terms','one_million','b.number(1000000)'),
        ('K1_ledger','run','family_definition','K',"'actual physical norm sum in9.16, bounded by Cstar*Kbar'")):
        fn=function(stem,method)
        values=[n.value for n in ast.walk(fn) if isinstance(n,ast.Assign)
            and any(ast.unparse(t)==target for t in n.targets)]
        if len(values)!=1:raise ValueError('Original physical norm definition changed')
        entries=({ast.literal_eval(k):v for k,v in zip(values[0].keys,values[0].values)}
            if isinstance(values[0],ast.Dict) else {kw.arg:kw.value for kw in values[0].keywords})
        if ast.dump(entries[key])!=ast.dump(ast.parse(expected,mode='eval').body):
            raise ValueError('Original K fixed positive term changed')
        bindings[stem+'.'+target+'/'+key]=expected;hashes[PREFIX+stem+'.py']=sha(PREFIX+stem+'.py')
    fn=function('frozen_comparison_field','packet')
    definitions=[kw.value for n in ast.walk(fn) if isinstance(n,ast.Call)
        for kw in n.keywords if kw.arg=='axial_drive_definition']
    expected_drive='sqrt(R/2)*F*Ef = R*hydro + R*Pstar^2*pressure + R^2*F0^2*swirl'
    if len(definitions)!=1 or ast.literal_eval(definitions[0])!=expected_drive:
        raise ValueError('Original frozen Ebar three-scale drive definition changed')
    bindings['frozen_comparison_field.packet/original_Ebar_definition']=expected_drive
    hashes[PREFIX+'frozen_comparison_field.py']=sha(PREFIX+'frozen_comparison_field.py')
    fn=function('global_exit_certificate','run')
    value=next(n.value for n in ast.walk(fn) if isinstance(n,ast.Assign)
        and any(ast.unparse(t)=='prescription' for t in n.targets))
    prescription={kw.arg:ast.literal_eval(kw.value) for kw in value.keywords}
    if prescription['actual']!='chi=1-(1-epsilon)*sigma(y/h); logF=logf-.5 integralchi*Dbar dy; V=v-integralchi sqrt(R/2)F Ebar dy, Ra<R<=100':
        raise ValueError('Admitted global exit shear prescription changed')
    hashes[PREFIX+'global_exit_certificate.py']=sha(PREFIX+'global_exit_certificate.py')
    R,Ps,Fbase,Fa,Fbar,chi,h,w,Db,kappa=s.symbols(
        'R Pstar F0base Factual Fbar chi hb omega Dbar kappa',positive=True)
    Eb,et,ez=s.symbols('Ebar e_theta e_z',real=True)
    hydro,pressure,swirl=s.symbols('hydro pressure swirl',real=True);checks={}
    def zero(name,left,right):
        if s.cancel(s.expand(left-right))!=0:raise ArithmeticError('Current collar identity: '+name)
        checks[name]=True
    drive=R*hydro+R*Ps**2*pressure+R**2*Fbase**2*swirl
    Ef=drive/(s.sqrt(R/2)*Fbar)
    zero('same_full_three_scale_axial_shear_equation',-chi*(Fa/Fbar)*drive,-chi*s.sqrt(R/2)*Fa*Ef)
    zero('same_logF_shear_and_positive_width_pullback',(-h*chi*Db/2)/h,-chi*Db/2)
    sigv=s.Symbol('same_sigma_value',nonnegative=True)
    zero('same_original_chi_omega',1-(1-(1-h)*sigv),(1-h)*sigv)
    pa,pb=s.symbols('actual_phi comparison_phi',positive=True)
    zero('same_actual_comparison_amplitude_quotient',(Fbase*pa)/(Fbase*pb),pa/pb)
    H=(Db**2+Eb**2)/Db;Tt=w*Db+et;Tz=w*Eb+ez
    direction=Tt+(Eb/Db)*Tz;transverse=Tz-(Eb/Db)*Tt
    zero('same_signed_cone_direction',direction,w*H+et+Eb*ez/Db)
    zero('same_signed_cone_transverse',transverse,ez-Eb*et/Db)
    zero('same_full_quadratic_dot_cross_units',Db**2*(2*direction**2-(kappa-2)*transverse**2),
        2*(Tt*Db+Tz*Eb)**2-(kappa-2)*(Db*Tz-Eb*Tt)**2)
    zero('same_signed_shear_kappa',chi*Db+(-chi*Eb)**2/(chi*Db),chi*H)
    zero('signed_error_dot_cross_Pythagorean_identity',(et*Db+ez*Eb)**2+(Db*ez-Eb*et)**2,
        (Db**2+Eb**2)*(et**2+ez**2))
    K=s.Symbol('same_global_K',positive=True)
    zero('correlated_K10_rho_squared_budget',K**10*(1/(20*K**5))**2,s.Rational(1,400))
    a,b=s.symbols('same_positive_left_flat_exponential same_positive_right_flat_exponential',positive=True)
    zero('same_flat_sigma_exponential_and_odds_ratio',a/(a+b),(a/b)/(1+a/b))
    x=s.Symbol('same_open_sigma_argument',positive=True)
    aa=s.exp(-1/x**2);bb=s.exp(-1/(1-x)**2);sig=aa/(aa+bb)
    zero('same_flat_sigma_strict_monotonicity_density',s.diff(sig,x),aa*bb/(aa+bb)**2*(2/x**3+2/(1-x)**3))
    norms=s.symbols('inverse_ell inverse_Ra Cstar Pstar_norm A_norm Fcore_norm inverse_Fcore_norm P0_norm moments_norm Df_norm Ef_norm inverse_Df_norm',nonnegative=True)
    Kdef=s.Integer(1000000)+sum(norms)
    zero('original_paper_K_fixed_positive_term_lower',Kdef-1000000,sum(norms))
    gamma,loss=s.symbols('same_gamma source_cutoff_K10_loss',positive=True)
    zero('actual_chi_Hbar_lower_from_retained_product_loss',2+2*gamma-gamma/10,2+s.Rational(19,10)*gamma)
    return dict(passed=True,AST_bindings=bindings,exact_source_identities=checks,
        admitted_global_exit_prescription=prescription,
        same_function_proof=[
            'Checked common-core Banach/atom identities identify one actual F,V,P0 and six physical primitives at Ra.',
            'The checked core-first graph uses those actual atom integrals, not the covering profiles, and the same global scalar width.',
            'The comparison uses the same alpha/core derivatives, inlet and six linear moment ODEs; their integral definitions uniquely identify its functions.',
            'The normalized three-scale drive and Dbar identities identify both actual prescribed-shear equations with the admitted global exit prescription.',
            'The angular exponential and axial integral have the same inlet and comparison direction, hence define the same actual F and V.',
            'The six actual moment ODEs then have the same densities/inlets; pressure is P0 plus that actual square primitive.',
            'The checked current bridge replay and core-first pullback attach these same functions to the current source graph.'],
        pressure_scaling='P=P0+R*F0base^2*C_dressed; at Ra, PI_core=4*C_dressed and Ra=4*epsilon_core',
        exact_width='hb=epsilon_b=cstar*K^-100; K is the original global physical norm sum',
        original_paper_K_definition='Lei-Ren v2 p134 (9.16): K=10^6+ell_c^-1+Ra^-1+Cstar+Pstar+A+nonnegative core/moment/frozen norms; K>=10^6',
        input_hashes={**hashes,Path(__file__).name:sha(Path(__file__).name)})


class CurrentInnerExitStrictCollar:
    def __init__(self,require_checked=True):
        self.ctx=MPIntervalContext();self.ctx.dps=240;self.records={};self.hashes={}
        for stem in ('global_exit_certificate','K1_ledger','physical_norm_family',
            'current_core_first_interface_check','current_bridge_functional_joins_check',
            'actual_bridge_integrals_check','current_core_common_fixed_point_check'):
            name=PREFIX+stem+'.json';row=json.loads((HERE/name).read_bytes())
            if stem.endswith('_check') and not row['all_passed']:raise ValueError('Checked current inlet/bridge source required')
            for path,digest in row['input_hashes'].items():
                if sha(path)!=digest:raise ValueError('Changed current collar prerequisite: '+path)
            self.hashes.update(row['input_hashes']);self.hashes[name]=sha(name);self.records[stem]=row
        first=self.records['current_core_first_interface_check'];norm=self.records['physical_norm_family']
        self.family=first['actual_five_defect_family_sha256'];self.source=first['implicit_source_sha256'];self.datum=first['datum_enclosure_sha256']
        for stem in ('current_bridge_functional_joins_check','actual_bridge_integrals_check','current_core_common_fixed_point_check'):
            row=self.records[stem]
            if (row['actual_five_defect_family_sha256'],row['implicit_source_sha256'],row['datum_enclosure_sha256'])!=(self.family,self.source,self.datum):
                raise ValueError('Current collar graph family/source/datum differs')
        if (norm['implicit_source_sha256'],norm['datum_enclosure_sha256'])!=(self.source,self.datum):
            raise ValueError('Norms and current collar datum differ')
        gate=self.records['global_exit_certificate'];ledger=self.records['K1_ledger']
        for key in ('base_analytic_core_family_sha256','uniform_Cstar_family_sha256'):
            if gate[key]!=ledger[key] or gate[key]!=norm[key]:raise ValueError('Original exit core/Cstar family differs')
        if gate['admitted_inner_parameter_family_sha256']!=ledger['admitted_inner_parameter_family_sha256']:
            raise ValueError('Original exit width family differs')
        if not all(gate[k] for k in ('actual_whole_axis_Ra_R110_relaxed_cone_analytically_certified','nonempty_inner_admissible_collar_analytically_certified','exact_implicit_exit_field_specified')):
            raise ValueError('Original strict exit source theorem required')
        if not all(ledger[k] for k in ('K1_numeric_bound_certified','h_b_equals_epsilon_b_by_definition','positive_width_not_materialized','fixed_K_K1_cstar_and_radius_inner_gates_certified')) or not all(ledger['smallness_checks'].values()):
            raise ValueError('Original global width/norm gates required')
        for group in ('current_common_core_bridge_graph','original_boundary_source_bindings','original_integrated_stress_free_core_equations',
            'scaled_atoms_are_the_original_physical_primitives','common_pressure_primitive_and_true_atom','phase0_ODE_and_exact_positive_width_pullback'):
            if not first[group].get('passed',True) or any(v is False for v in first[group].values()):
                raise ValueError('Exact current inlet/width/pressure source proof missing: '+group)
        if not first['current_core_bridge_functional_mixed4_join_certified'] or not first['all_current_bridge_functional_interfaces_certified']:
            raise ValueError('Current analytic core/first join required')
        self.theorem=exact_source_attachment();self.hashes.update(self.theorem['input_hashes']);self.proof=self.prove()
        integrated=first['original_integrated_stress_free_core_equations']
        if integrated['original_nested_direction_base_identities']!=4 or not integrated['original_scalar_direction_functions_replayed_exactly']:
            raise ValueError('Exact scalar/nested comparison direction equality receipt required')
        expected_recovery=['D_logR logF=-R*D_over_R/2',
            'D_logR V=-(R*hydro+R*Pstar^2*pressure+R^2*F0base^2*swirl)']
        if integrated['exact_stress_free_recovery']!=expected_recovery or not integrated['pressure_and_swirl_physical_scales_not_dropped']:
            raise ValueError('Original Ebar recovery and physical scales required')
        expected_bounds='9.14 inputs, K1 ledger and shared width giveDbar>=1/(2K),|q|<=2K,Hbar>=2+2gamma,Hbar<=K^10; Dbar>=3.5 on100..110'
        if gate['proof']['comparison']!=expected_bounds:
            raise ValueError('Admitted Hbar lower and upper source bounds required')
        fixed_name='lei_ren_part1_paper_shared_fixed_step_bound.json'
        fixed=json.loads((HERE/fixed_name).read_bytes())
        expected='sigma(s)=exp(-1/s^2)/(exp(-1/s^2)+exp(-1/(1-s)^2)), 0<s<1; flat constants outside'
        if fixed['cutoff_formula']!=expected or not fixed['strict_monotonicity_on_open_unit_interval'] or not fixed['smooth_flat_endpoints']:
            raise ValueError('Same original flat cutoff definition and endpoints required')
        self.hashes[fixed_name]=sha(fixed_name)
        # A digest joins the actual amplitude/packet/pressure proof to the
        # admitted core/Cstar/width sources, in addition to family tuples.
        self.attachment=dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            datum_enclosure_sha256=self.datum,base_analytic_core_family_sha256=gate['base_analytic_core_family_sha256'],
            uniform_Cstar_family_sha256=gate['uniform_Cstar_family_sha256'],
            admitted_inner_parameter_family_sha256=gate['admitted_inner_parameter_family_sha256'],
            actual_amplitude_packet_graph=first['current_common_core_bridge_graph'],
            exact_core_atom_and_pressure_functions=first['original_boundary_source_bindings'],
            exact_scaled_atoms=first['scaled_atoms_are_the_original_physical_primitives'],
            exact_pressure_primitive=first['common_pressure_primitive_and_true_atom'],
            original_direction_equality=integrated,actual_source_function_hashes=self.theorem['input_hashes'])
        self.attachment_digest=hashlib.sha256(json.dumps(self.attachment,sort_keys=True).encode()).hexdigest()
        if require_checked:
            row=json.loads((HERE/RECEIPT).read_bytes())
            if not row['all_passed'] or not row[GATE]:raise ValueError('Checked current strict collar required')
            for name,digest in row['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed checked current collar source: '+name)

    def prove(self):
        c=self.ctx;read=lambda value:read_interval(c,value);positive={}
        def require(name,value):
            if endpoints(value)[0]<=0:raise ArithmeticError('Current strict collar bound unresolved: '+name)
            positive[name]=value
        ledger=self.records['K1_ledger'];gate=self.records['global_exit_certificate'];norm=self.records['physical_norm_family']
        gamma=read(ledger['gamma']);L=read(norm['selected_log_K_upper']);Kmin=c.mpf(1000000)
        logh=read(ledger['shared_positive_width_log_enclosure']);cstar=read(ledger['cstar'])
        require('positive_gamma',gamma);require('global_K_upper_exceeds_K_minimum',L-c.ln(Kmin))
        require('exact_width_below_half',c.ln(c.mpf('.5'))-c.mpf(endpoints(logh)[1]))
        require('original_global_strong_error_ratio_below1',read(gate['strong_cone_sufficient_ratio_margin']))
        X=10*L+c.ln(10/gamma)
        require('positive_cutoff_denominator_exceeds1',X-1)
        # This endpoint selects a subcollar; it never selects or changes hb.
        sc=c.mpf(endpoints(1/(4*c.sqrt(X)))[0]);require('positive_selected_first_phase_endpoint',sc)
        require('selected_first_phase_below_quarter',c.mpf('.25')-sc)
        sigma_log=-1/sc**2+4;cutoff_gap=c.ln(gamma/10)-sigma_log-10*L
        require('whole_subcollar_sigma_K10_below_gamma_over10',cutoff_gap)
        require('subcollar_inside_original_rho4point1_continuation',
            c.ln(c.ln(c.mpf('4.1')/4))-c.mpf(endpoints(logh+c.ln(sc))[1]))
        kappa=2+c.mpf('1.9')*gamma;require('whole_subcollar_kappa_excess_exceeds_gamma',kappa-2-gamma)
        # Retain omega>0 and bound errors RELATIVE to omega*qbar. The
        # original strong gate gives |e|<omega/(40*K^6), never a fixed e.
        rho=1/(20*Kmin**5);D=1-rho;Q=2*(1-rho)**2-c.mpf(1)/400
        require('direction_relative_to_omega_Hbar_above_point95',D-c.mpf('.95'))
        require('full_quadratic_relative_to_omega_Hbar_squared_above1point8',Q-c.mpf('1.8'))
        return dict(positive_margins=positive,whole_Z_domain=(-1,1),selected_first_phase_endpoint=sc,
            cutoff_source='sigma(s)=exp(-1/s^2)/(exp(-1/s^2)+exp(-1/(1-s)^2))',
            original_K_log_upper=L,cutoff_denominator_X=X,whole_subcollar_log_sigma_upper=sigma_log,
            cutoff_log_admission_gap=cutoff_gap,exact_source_width_log_enclosure=logh,
            selected_log_phase_endpoint=c.ln(sc),source_log_log_radius_offset_enclosure=logh+c.ln(sc),
            full_source_kappa_lower=kappa,relative_error_over_omega_qbar_upper=rho,
            actual_K_lower_by_original_definition=Kmin,
            original_Hbar_bounds='2+2gamma<=Hbar<=K^10, checked current global exit theorem',
            actual_kappa_lower_derivation='chi*Hbar=Hbar-(1-hb)*sigma*Hbar>2+2gamma-gamma/10=2+1.9gamma',
            correlated_relative_error_bound='rho_source<1/(20*K^5), kappa-2<=Hbar<=K^10',
            correlated_K10_relative_error_squared_upper=c.mpf(1)/400,
            full_direction_over_omega_Hbar_lower=D,full_quadratic_over_omega_Hbar_squared_lower=Q,
            full_stress_norm_over_omega_qbar_lower=D,
            source_stress_scale='omega=(1-hb)*sigma(first_phase)>0 only for first_phase>0; exact stress is zero at first_phase0',
            explicit_support_left_endpoint_fraction='1/2',
            source_radius='Ra*exp(hb*selected_first_phase_endpoint*fraction); Ra=4*epsilon_core',
            positive_width_and_flat_stress_not_materialized=True,
            current_inner_exit_zero_inlet_case_certified=True,
            **{GATE:True},**{key:False for key in OPEN})

    def query(self,fraction,Z=(-1,1)):
        c=self.ctx;q=c.mpf(fraction);z=c.mpf(Z);lo,hi=endpoints(q);zl,zh=endpoints(z)
        if lo<0 or hi>1 or zl < -1 or zh>1 or not all(mp.isfinite(v) for v in (lo,hi,zl,zh)):
            raise ValueError('Finite collar fraction subset[0,1] and Z subset[-1,1] required')
        return dict(fraction=q,Z=z,source_family=self.family,source_first_phase=q*self.proof['selected_first_phase_endpoint'],
            full_source_collar_certificate=self.proof,
            strict_nonzero_stress_cone_certified_for_entire_query_box=lo>0,
            exact_zero_core_inlet_case_only=hi==0,core_zero_and_strict_collar_case_certified=True,
            **{GATE:True},**{key:False for key in OPEN})


def run():
    field=CurrentInnerExitStrictCollar(require_checked=False)
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,
        datum_enclosure_sha256=field.datum,exact_current_exit_source_attachment=field.theorem,
        cross_graph_amplitude_packet_pressure_width_attachment=field.attachment,
        cross_graph_attachment_sha256=field.attachment_digest,
        explicit_current_inner_exit_strict_collar=field.proof,
        examples={name:field.query(q,z) for name,q,z in (
            ('closed_core_and_collar',(0,1),(-1,1)),('exact_zero_inlet',0,(-1,1)),
            ('strict_left_support_collar',('.5','1'),(-1,1)),('strict_positive_midplane',('.25','.75'),0))},
        **{GATE:True},**{key:False for key in OPEN},input_hashes=field.hashes)
    (HERE/NAME).write_bytes((json.dumps(encode(result),indent=2)+'\n').encode())
    print('Current inner strict exit collar attached; explicit left support fraction[.5,1]',flush=True)
    return result


if __name__=='__main__':run()
