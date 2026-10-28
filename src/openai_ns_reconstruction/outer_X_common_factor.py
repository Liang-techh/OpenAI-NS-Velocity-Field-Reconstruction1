"""Six-chart leading E/h common-factor theorem with explicit source gates.

No original chart integrals, source graphs or ancestor builders are run.
The scalar Duhamel terms stand for exact unevaluated defining integrals.
"""
import ast
import base64
import gzip
import hashlib
import json
import sympy as s

P='lei_ren_part1_paper_compliant_current_original_'
OUTER_SRC='experiments/root_st073/'+P+'true_outer_C6_C5_functions.py'
OUTER_SHA='4464a934d303a95a0da157e9e289f6ed8a04d5ff0adcb0d663857f0d9c57dcd3'
RP='experiments/root_st073/'+P+'Rp_native_C5_pulse_source.json.gz'
RP_SHA='442df70e5ab4051066b241e943c8d0531067d057e774a04dcc341dadbd703ebe'
OUTER=P+'true_outer_C6_C5_functions.json.gz'
REPORT='experiments/root_st073/lei_ren_part1_paper_compliant_current_compact_X_cleared_residual.json'
REPORT_SHA='fff6e2b5d9848de525207964546fc84f4f8c9024b1569a938b44735848b4af76'
RESULT='tmp/compact_X_cleared_residual_v1/result.json'
RESULT_SHA='9cd2f4400545eb02c966da06012e790b87ae86dba201bd3f638b25b568ecf14e'
CHARTS=('Rh_reference','O2_slope','O2_axial','O2_buffer','O3_transition','O3_power')


def source_bytes(root,name,sha):
    raw=(root/name).read_bytes()
    if hashlib.sha256(raw).hexdigest()!=sha:raise ValueError('changed input: '+name)
    return raw


def scalar_stage(previous_E,previous_h,scale,decay,forcing_integral):
    """One exact Duhamel step on scalar coefficients, preserving memory."""
    return previous_E*scale,previous_h*decay+forcing_integral


def recover(root):
    raw=source_bytes(root,OUTER_SRC,OUTER_SHA); code=raw.decode('utf-8').replace('\r\n','\n')
    tree=ast.parse(code)
    cls=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=='CurrentOriginalTrueOuterC6C5Functions')
    init=next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name=='__init__')
    assignments=[n.value for n in ast.walk(init) if isinstance(n,ast.Assign)
                 and any(isinstance(t,ast.Name) and t.id=='E' for t in n.targets)]
    expected={"self.leading_Rh['fields']['E']",
              "self.background_endpoints['Rh_reference']['fields']['E']",
              "self.background_endpoints['O2_slope']['fields']['E']",
              "self.background_endpoints['O2_axial']['fields']['E']",
              "self.background_endpoints['O2_buffer']['fields']['E']",
              "self.background_endpoints['O3_transition']['fields']['E']"}
    if len(assignments)!=6:raise ValueError('unexpected number of outer E assignments')
    predecessors=[]; factor_syntax=[]
    for expr in assignments:
        if not isinstance(expr,ast.Call) or not isinstance(expr.func,ast.Attribute) or not isinstance(expr.func.value,ast.Name) or expr.func.value.id!='a' or expr.func.attr!='scale' or len(expr.args)!=2:
            raise ValueError('outer E is not a two-argument scale')
        predecessors.append(ast.unparse(expr.args[0]))
        if any((isinstance(n,ast.Name) and n.id=='Z') or (isinstance(n,ast.Constant) and n.value=='Z') for n in ast.walk(expr.args[1])):
            raise ValueError('explicit Z found in scalar factor syntax')
        factor_syntax.append(ast.unparse(expr.args[1]))
    if set(predecessors)!=expected:raise ValueError('outer E predecessor chain differs')
    anchors=["return dict(m=V,h=E,k=a.mul(E,V)",
             "out[key]=a.add(a.scale(incoming[key],decay),local[key])",
             "H,parts,memory=transport_geometry(g,background_density(g,E,V),incoming,x,domain[0],offset,jac,name,6)",
             "incoming=self.background_endpoints[name]['history']",
             "canonical_normalized_Rh='E=exp(-1/2)/(1+Z^2); V=m=4Z/Pstar; h=5E/8;"]
    if any(x not in code for x in anchors):raise ValueError('outer source definition anchor missing')
    report=json.loads(source_bytes(root,REPORT,REPORT_SHA));result_raw=base64.b64decode(report['artifacts'][RESULT]['base64'],validate=True)
    if hashlib.sha256(result_raw).hexdigest()!=RESULT_SHA:raise ValueError('prior compact result changed')
    prior=json.loads(result_raw)
    if prior.get('status')!='returned':raise ValueError('prior compact result did not return')
    prior=prior['result']
    rp=json.loads(gzip.decompress(source_bytes(root,RP,RP_SHA)))
    outer=rp['Rp_native_C5_pulse_source_bindings']['inputs'][OUTER]
    roles=[]
    for chart in CHARTS:
        bg=outer['background'][chart]
        roles.append({'chart':chart,'E_rows':bg['fields']['E'],'h_rows':bg['fields']['h'],
                      'phase_ref':bg['phase'],'domain_refs':bg['domain'],
                      'log_offset_ref':bg['log_offset'],'jacobian_ref':bg['physical_measure']})
    seed=outer['seed_identity']
    if 'h=5E/8' not in seed['canonical_normalized_Rh']:
        raise ValueError('saved Rh seed differs')
    # Mathematical induction with opaque, Z-independent scalar Duhamel
    # integrals J_i; no source integral is computed or numerically replaced.
    z=s.symbols('z',real=True);a=s.exp(-s.Rational(1,2))/(1+z*z)
    factors=s.symbols('f0:6',nonzero=True);decays=s.symbols('d0:6');J=s.symbols('J0:6')
    F=s.Integer(1);C=s.Rational(5,8);checks={};steps=[]
    for i,chart in enumerate(CHARTS):
        F,C=scalar_stage(F,C,factors[i],decays[i],J[i])
        e=a*F;h=a*C
        checks['common_factor_'+chart]=s.cancel(h/e-C/F)==0
        steps.append({'chart':chart,'scalar_E_factor':str(F),'scalar_h_history':str(C)})
    checks.update({'final_ratio_Z'+str(n):s.simplify(s.diff((a*C)/(a*F),z,n))==0 for n in range(1,6)})
    Fc,Cc=s.symbols('F_compact C_compact',nonzero=True)
    Eb,Hb=a*Fc,a*Cc;dE,dh=s.symbols('delta_E delta_h')
    checks['terminal_ratio_if_corrections_zero']=s.cancel(((Hb+dh)/(Eb+dE)).subs({dh:0,dE:0})-Cc/Fc)==0
    return {'schema':'outer_X_common_factor_v1','status':'source_model_conditional_factorization',
            'source_ast_scale_count':len(assignments),'source_predecessors':predecessors,
            'source_factor_syntax':factor_syntax,'saved_chart_roles':roles,'saved_Rh_seed':seed['canonical_normalized_Rh'],
            'source_family':rp['source_family'],'inductive_scalar_steps':steps,'checks':checks,'checks_passed':all(checks.values()),
            'prior_terminal_evidence':prior['inherited_source_evidence'],
            'conditions':['Rh seed is the same all-Z source function','sigma, mu and radial geometry are independent of Z',
                'each source scale and Duhamel integral has its declared mathematical semantics',
                'all chart endpoint substitutions preserve the same source functions',
                'compact terminal delta_E and delta_h vanish as functions on the full Z interval',
                'the source jet rows are ordinary derivatives and denominators are nonzero'],
            'conditional_conclusion':'background h/E is Z-independent; modified terminal ratio agrees if both terminal corrections vanish',
            'actual_outer_graph_factorization_proved':False,'actual_terminal_function_identity_proved':False,
            'actual_compact_P4_P5_zero_proved':False,'actual_current_X4_X5_zero_proved':False,
            'source_values_computed':False,'integrals_evaluated':False,'old_errors_removed':False,
            'runtime_model_metadata':None}
