"""Stateless original O2 slope point profiles in the common Pstar basis.

Actual defining integrals, not saved-cover values, produce finite point
coefficients. Pstar^{-1} and Pstar^{-2} remain exact formal source factors;
no astronomical original exponential is materialized. The all-N E/V roles
are bound to the accepted original graph. Quadrature remains approximate.
This partial point provider does not supply p1/p2, P0 or the full oracle.
"""
import ast
from dataclasses import dataclass
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_generic_shear_loop as loop

HERE,PREFIX,sha=loop.HERE,loop.PREFIX,loop.sha
NAME=PREFIX+'current_original_O2_slope_point_profiles.json'
RECEIPT=PREFIX+'current_original_O2_slope_point_profiles_check.json'
GATE='original_O2_slope_point_coefficients_and_common_basis_E_V_dispatch_implemented'
VIEWS=PREFIX+'current_generic_loop_function_sources_views.json.gz'
ALLN=PREFIX+'current_native_Rc_all_N_function_controls.json'
ALLN_CHECK=PREFIX+'current_native_Rc_all_N_function_controls_check.json'


@dataclass(frozen=True)
class PstarPoint:
    """Sum of true point coefficients times exact original Pstar powers.

    All coefficients are approximate point evaluations. They are never
    enclosures or fixed cover midpoints. Factor exponents are exact integers.
    """
    terms: tuple
    approximate_coefficients: bool=True
    original_factor_basis: str='same original positive Pstar'

    def unscaled_scalar(self):
        if any(power!=0 for power,unused in self.terms):
            raise ValueError('Original Pstar factor is unresolved; scalar cast forbidden')
        return sum((value for unused,value in self.terms),0)

    def record(self):
        return dict(terms=[dict(original_Pstar_power=power,point_coefficient=value) for power,value in self.terms],
            approximate_point_coefficients=True,source_caps_or_midpoints_used=False,
            original_positive_Pstar_factor_not_materialized=True)


@dataclass(frozen=True)
class PointC1:
    value: PstarPoint
    Z: PstarPoint

    def record(self):return dict(C0=self.value.record(),Z=self.Z.record())


def original_recipe_binding():
    name=PREFIX+'pre_pulse_mixed_C4.py';tree=ast.parse((HERE/name).read_text(encoding='utf8'))
    cls=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=='CompliantPrePulseMixedC4')
    fn=next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name=='slope')
    assignments={ast.unparse(t):n.value for n in ast.walk(fn) if isinstance(n,ast.Assign) for t in n.targets}
    expected={
        '(J, mass)':'slope_masses(c,y,self.cells)',
        'factor':"c.exp(y/10-c.mpf('.6')*J)",
        'u':'qi*factor','V':'z*4',
        'h':"qi*(c.mpf('.625')+mass[0])*c.exp(-c.mpf('1.5')*y)",
        'hist':"dict(m=V,h=h,k=h*V,e=square(V)*self.invP2-square(qi)*(c.mpf(5)/12+mass[2]/2)*c.exp(-y),p=square(qi)*(c.mpf('2.5')+mass[1]/2))"}
    for key,value in expected.items():
        if key not in assignments or ast.dump(assignments[key])!=ast.dump(ast.parse(value,mode='eval').body):
            raise ValueError('Original O2 slope source recipe changed: '+key)
    massfile='lei_ren_part1_paper_interval_outer_slope_field.py'
    mt=ast.parse((HERE/massfile).read_text(encoding='utf8'))
    mf=next(n for n in mt.body if isinstance(n,ast.FunctionDef) and n.name=='transition_integrals')
    ma={ast.unparse(t):n.value for n in ast.walk(mf) if isinstance(n,ast.Assign) for t in n.targets}
    for key,value in (('rates',"(c.mpf('1.6'),c.mpf('.2'),c.mpf('1.2'))"),('powers','(1,2,2)')):
        if ast.dump(ma[key])!=ast.dump(ast.parse(value,mode='eval').body):raise ValueError('Original slope mass definition changed')
    want=ast.parse("dy*c.exp(rate*scell-c.mpf('.6')*power*jcell)",mode='eval').body
    if not any(isinstance(n,ast.AugAssign) and ast.unparse(n.target)=='masses[k]' and ast.dump(n.value)==ast.dump(want) for n in ast.walk(mf)):
        raise ValueError('Original slope integral density changed')
    return dict(passed=True,original_source_assignments_bound=len(expected),original_mass_densities=3,
        input_hashes={name:sha(name),massfile:sha(massfile),Path(loop.__file__).name:sha(Path(loop.__file__).name)})


def original_background_leaf_identity():
    """Bind the actual source chain from pre-pulse profiles to E/V leaves."""
    hashes={};identities=[]
    def method(stem,name,cls=None):
        filename=PREFIX+stem+'.py';hashes[filename]=sha(filename)
        tree=ast.parse((HERE/filename).read_text(encoding='utf8'))
        body=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name==cls).body if cls else tree.body
        return next(n for n in body if isinstance(n,ast.FunctionDef) and n.name==name)
    def assignment(fn,key,value):
        rows={ast.unparse(t):n.value for n in ast.walk(fn) if isinstance(n,ast.Assign) for t in n.targets}
        if key not in rows or ast.dump(rows[key])!=ast.dump(ast.parse(value,mode='eval').body):
            raise ValueError('Original background source chain changed: '+key)
        identities.append(key)
    fn=method('pre_pulse_mixed_C4','packet','CompliantPrePulseMixedC4')
    for key,value in (('Utheta_over_Pstar_axial5_coefficients','list(u.coefficients)'),('Uz_ordinary_y_derivative_axial5','V')):
        if not any(isinstance(n,ast.keyword) and n.arg==key and ast.dump(n.value)==ast.dump(ast.parse(value,mode='eval').body) for n in ast.walk(fn)):
            raise ValueError('Original point source packet changed: '+key)
        identities.append(key)
    fn=method('current_pre_pulse_stress_operator','raw_pre_velocity_rows')
    assignment(fn,'u0',"IntervalTaylor(c,[c.mpf(endpoints(value)) for value in pre['Utheta_over_Pstar_axial5_coefficients']])")
    assignment(fn,'u','[u0*value for value in exponential_derivatives(logjet)]')
    assignment(fn,'V',"[copy_jet(c,value) for value in pre['Uz_ordinary_y_derivative_axial5']]")
    fn=method('long_reshape_mixed_C4','exponential_derivatives')
    assignment(fn,'one','log_derivatives[0]*0+1');assignment(fn,'out','[one]')
    fn=method('current_native_generic_source_packets','query','NativeGenericSourcePackets')
    assignment(fn,'velocity',"{key:tuple(algebra.shift(row,packets.INVERSE_S) if key in ('axial','radial') else row for row in values) for key,values in native_v.items()}")
    fn=method('current_native_generic_source_packets','original_raw','NativeGenericSourcePackets')
    assignment(fn,'logs',"raw['algebra'].logs if 'algebra' in raw else (c.mpf(0),2*self.seed.logP,c.mpf(0),c.mpf(0))")
    fn=method('current_generic_shear_signed_jets','source_leaves')
    assignment(fn,'E',"packet.velocity['theta']");assignment(fn,'V',"packet.velocity['axial']")
    filename=PREFIX+'current_generic_shear_source_packets.py';hashes[filename]=sha(filename)
    tree=ast.parse((HERE/filename).read_text(encoding='utf8'))
    shift=next(n.value for n in tree.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='INVERSE_S' for t in n.targets))
    if ast.literal_eval(shift)!=(0,-.5,0,0):raise ValueError('Original common Pstar shift changed')
    return dict(passed=True,source_chain_AST_identities=identities,
        original_E_source_leaf_equals_original_u_and_its_ordinary_Z_row=True,
        original_V_source_leaf_equals_4Z_times_original_Pstar_inverse_and_its_Z_row=True,
        coefficient_to_ordinary_Z_order0_1_factorials_are_one=True,
        common_inverse_factor_identity='exp((-1/2)*(2*logP))=Pstar^-1',
        source_function_identity_not_equality_of_enclosure_endpoints=True,input_hashes=hashes)


class OriginalO2SlopePointProfiles:
    mode='original_point_coefficients_with_formal_Pstar'
    def __init__(self,dps=50):
        if type(dps) is not int or dps<40:raise ValueError('At least 40 decimal digits required')
        self.ctx=c=mp.mp.clone();c.dps=dps;self.J_cache={};self.radial_cache={}
        self.recipe=original_recipe_binding();self.leaf_identity=original_background_leaf_identity()
        self.hashes={**self.recipe['input_hashes'],**self.leaf_identity['input_hashes']}
        checked=json.loads((HERE/ALLN_CHECK).read_bytes());manifest=json.loads((HERE/ALLN).read_bytes())
        if not checked['all_passed'] or checked['source_family']!=manifest['source_family']:
            raise ValueError('Accepted original all-N source family required')
        if checked['input_hashes'][ALLN]!=sha(ALLN):raise ValueError('All-N original source graph changed')
        for name,digest in self.hashes.items():
            if checked['input_hashes'].get(name)!=digest:
                raise ValueError('Source recipe/leaf path differs from the accepted original receipt: '+name)
        self.source_graph_sha=manifest['input_hashes'][VIEWS]
        if self.source_graph_sha!=sha(VIEWS):raise ValueError('Original loop source graph changed')
        self.family=manifest['source_family'];view=json.loads(gzip.decompress((HERE/VIEWS).read_bytes()))['O2_slope']
        if view['source_family']!=self.family:raise ValueError('O2 source belongs to another family')
        dag=view['original_signed_input_graph']['jet_expression_dag'];nodes=view['function_graph_nodes']
        self.roles={}
        for role,order in (('all_N_original_E_C0','y0_Z0'),('all_N_original_E_Z','y0_Z1')):
            root=dag['nodes'][dag['roots']['E'][order]]
            if root!={'operation':'source_derivative','name':'E_'+order}:
                raise ValueError('Original E graph root is not its bound derivative leaf')
            self.roles[role]=('original_signed_input_graph.jet_expression_dag',dag['roots']['E'][order])
        for role,label in (('all_N_original_V_C0','V_y0_Z0'),('all_N_original_V_Z','V_y0_Z1')):
            ids=[i for i,row in enumerate(nodes) if row.get('operation')=='source_derivative' and row.get('name')==label]
            if len(ids)!=1:raise ValueError('Unique original axial source leaf required')
            self.roles[role]=('function_graph_nodes',ids[0])
        self.hashes.update({name:sha(name) for name in (ALLN,ALLN_CHECK,VIEWS,Path(__file__).name)})

    def coordinate(self,y):
        c=self.ctx;y=c.mpf(y)
        if not c.isfinite(y) or not 0<=y<=1:raise ValueError('Original O2 slope coordinate y in[0,1] required')
        return y

    def J(self,y):
        y=self.coordinate(y);c=self.ctx;key=(c.prec,y._mpf_)
        if key not in self.J_cache:
            if y==0:value=c.mpf(0)
            elif y==1:value=c.mpf('.5')
            else:value=c.quad(lambda s:loop.flat_step(c,s),[0,y/2,y])
            self.J_cache[key]=value
        return self.J_cache[key]

    def radial(self,y):
        y=self.coordinate(y);c=self.ctx;key=(c.prec,y._mpf_)
        if key not in self.radial_cache:
            J=self.J(y);f=c.exp(y/10-c.mpf('.6')*J);sigma=loop.flat_step(c,y)
            masses=[c.quad(lambda s:c.exp(rate*s-c.mpf('.6')*power*self.J(s)),[0,y/2,y]) if y else c.mpf(0)
                for rate,power in ((c.mpf('1.6'),1),(c.mpf('.2'),2),(c.mpf('1.2'),2))]
            self.radial_cache[key]=dict(y=y,J=J,masses=masses,f=f,log_E_y=c.mpf('.1')-c.mpf('.6')*sigma,
                a=c.mpf('.8')+c.mpf('1.2')*sigma,
                H=(c.mpf('.625')+masses[0])*c.exp(-c.mpf('1.5')*y),
                D=(c.mpf(5)/12+masses[2]/2)*c.exp(-y),P=c.mpf('2.5')+masses[1]/2)
        return self.radial_cache[key]

    def evaluate(self,*,Z,y):
        c=self.ctx;z=c.mpf(Z)
        if not c.isfinite(z) or not -1<=z<=1:raise ValueError('Original Z in[-1,1] required')
        row=self.radial(y);C=1/(1+z*z);CZ=-2*z*C*C
        pp=lambda terms:PstarPoint(tuple(sorted((p,v) for p,v in terms.items() if v!=0)))
        pair=lambda terms,jets:PointC1(pp(terms),pp(jets))
        E=pair({0:C*row['f']},{0:CZ*row['f']});V=pair({-1:4*z},{-1:c.mpf(4)})
        h=pair({0:C*row['H']},{0:CZ*row['H']})
        histories=dict(m=V,h=h,
            k=pair({-1:4*z*C*row['H']},{-1:4*(C+z*CZ)*row['H']}),
            e=pair({-2:16*z*z,0:-C*C*row['D']},{-2:32*z,0:-2*C*CZ*row['D']}),
            p=pair({0:C*C*row['P']},{0:2*C*CZ*row['P']}))
        return dict(chart='O2_slope',Z=z,y=row['y'],source_family=self.family,E=E,V=V,histories=histories,
            a=pair({0:row['a']},{}),b=pair({},{}),
            log_E_y=row['log_E_y'],original_defining_integrals=dict(J=row['J'],masses=row['masses']),
            absolute_pressure=dict(original_analytic_P0_reference=dict(operation='same_original_analytic_pressure_datum',
                source_family=self.family,datum_enclosure_sha256=self.family['datum_enclosure_sha256'],
                point_value_and_Z_derivative_not_evaluated=True),cumulative_p=histories['p'],pressure_datum_reset=False),
            original_p1_p2_and_loop_A_B_not_evaluated=True,
            mathematical_original_Pstar_basis_preserved=True,
            approximate_coefficients=True,certified_numeric_error_bound=False)

    def background(self,*,Z,y):
        """Cheap E/V/a/b source query: cumulative masses are not needed."""
        c=self.ctx;z=c.mpf(Z);y=self.coordinate(y)
        if not c.isfinite(z) or not -1<=z<=1:raise ValueError('Original Z in[-1,1] required')
        C=1/(1+z*z);CZ=-2*z*C*C;f=c.exp(y/10-c.mpf('.6')*self.J(y))
        pp=lambda terms:PstarPoint(tuple(sorted((p,v) for p,v in terms.items() if v!=0)))
        pair=lambda terms,jets:PointC1(pp(terms),pp(jets))
        return dict(E=pair({0:C*f},{0:CZ*f}),V=pair({-1:4*z},{-1:c.mpf(4)}),
            a=pair({0:c.mpf('.8')+c.mpf('1.2')*loop.flat_step(c,y)},{}),b=pair({},{}))

    def dispatch_original_background(self,row,*,coordinate,Z):
        role=row.get('function_role')
        if row.get('operation')!='original_function_graph' or row.get('chart')!='O2_slope' or role not in self.roles:
            raise ValueError('Only original O2 slope background E/V roles are implemented')
        if row.get('graph_sha256')!=self.source_graph_sha or (row.get('source_graph_namespace'),row.get('source_node'))!=self.roles[role]:
            raise ValueError('Original source role/namespace/node binding differs')
        got=self.background(Z=Z,y=coordinate);key='E' if '_E_' in role else 'V'
        return getattr(got[key],'Z' if role.endswith('_Z') else 'value')


def run():
    began=time.monotonic();owner=OriginalO2SlopePointProfiles();frame=owner.evaluate(Z='.37',y='.53')
    record=dict(**{GATE:True},source_family=owner.family,original_recipe_binding=owner.recipe,
        original_E_V_graph_leaf_identity=owner.leaf_identity,
        original_source_recipe_and_leaf_chain_hashes_compared_to_accepted_receipt=True,
        exact_original_background_role_bindings={k:dict(namespace=v[0],source_node=v[1]) for k,v in owner.roles.items()},
        original_O2_slope_point=dict(Z=frame['Z'],y=frame['y'],E=frame['E'].record(),V=frame['V'].record(),
            histories={k:v.record() for k,v in frame['histories'].items()},a=frame['a'].record(),b=frame['b'].record(),
            original_defining_integrals=frame['original_defining_integrals']),
        original_finite_coefficients_from_defining_integrals_not_cover_endpoints=True,
        background_E_V_queries_skip_unneeded_cumulative_mass_quadrature=True,
        original_nonzero_axial_source_kept_with_exact_Pstar_inverse_factor=True,
        original_five_history_functions_and_ordinary_Z_rows_evaluated=True,
        original_absolute_pressure_datum_kept_separate=True,
        original_analytic_P0_reference=frame['absolute_pressure']['original_analytic_P0_reference'],
        original_p1_p2_and_phase_provider_installed=False,
        full_original_O2_source_point_provider_installed=False,
        quadrature_and_roundoff_certified=False,numerical_original_source_point_or_integral_oracle_installed=False,
        actual_five_controls_installed=False,certified_actual_fixed_point_tail_installed=False,
        actual_terminal_Z_function_closure_installed=False,current_whole_N_selected=False,
        **dict.fromkeys(loop.OPEN,False),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
        scope='Actual original O2 slope defining profile integrals with approximate point coefficients and exact common Pstar powers. Four background E/V source roles bound. Full p1/p2, analytic P0 point service, native phase, certified integrals, controls and global field remain open.')
    (HERE/NAME).write_text(json.dumps(loop.encoded(record),indent=2)+'\n',encoding='utf8')
    print('Original O2 slope point E/V and five history coefficients evaluated without source ancestors',flush=True)
    return record


if __name__=='__main__':run()
