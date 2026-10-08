"""True Rh_reference point coefficients and actual-radius C0/Z primitives.

The original reference source has closed radial profiles, so no defining
slope quadrature or ancestor constructor is rerun. Original inertial and
pressure programs, source factors, selected scales and inverse kernels are
reused. This is a second chart of the actual source oracle, not all 17.
"""
import ast
import copy
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_O2_conditioned_slow_Z as slow

base=slow.base;point=base.point;prior=base.prior;ep=point.endpoints
HERE,PREFIX,sha=base.HERE,base.PREFIX,base.sha
NAME=PREFIX+'current_original_reference_point_oracle.json'
RECEIPT=PREFIX+'current_original_reference_point_oracle_check.json'
GATE='original_Rh_reference_true_point_C0_Z_inputs_and_radius_phase_primitives_connected'


def reference_coordinate(value):
    value=point.source.exact_rational(value)
    if not -5<=value<=0:raise ValueError('Original Rh_reference offset in[-5,0] required')
    return value


def reference_recipe_binding():
    name=PREFIX+'pre_pulse_mixed_C4.py'
    tree=ast.parse((HERE/name).read_text(encoding='utf8'))
    cls=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=='CompliantPrePulseMixedC4')
    fn=next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name=='reference')
    rows={ast.unparse(t):n.value for n in ast.walk(fn) if isinstance(n,ast.Assign) for t in n.targets}
    wanted={'u':'qi*c.exp(y/10)','V':'z*4','h':"u*c.mpf('.625')",
        'hist':"dict(m=V,h=h,k=h*V,e=square(V)*self.invP2-square(u)*c.mpf(5)/12,p=square(u)*c.mpf('2.5'))"}
    for key,value in wanted.items():
        if ast.dump(rows[key])!=ast.dump(ast.parse(value,mode='eval').body):
            raise ValueError('Original reference source changed: '+key)
    return dict(passed=True,source=name,sha256=sha(name),closed_reference_assignments= wanted,
        exact_logarithmic_shear='a=1-2*(1/10)=4/5',b_exact_zero=True,
        original_full_five_histories_and_separate_P0_retained=True)


def project_method(module,clsname,method,expected_coordinates,replace_radius=False):
    """Change chart domain and equivalent radius metadata, never source math."""
    tree=ast.parse(Path(module.__file__).read_text(encoding='utf8'))
    cls=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name==clsname)
    original=next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name==method)
    fn=copy.deepcopy(original);changes=[];radius=[]
    class Project(ast.NodeTransformer):
        def visit_Call(self,node):
            self.generic_visit(node)
            name=ast.unparse(node.func)
            if name in ('source.exact_coordinate','point.source.exact_coordinate'):
                changes.append(name);node.func=ast.Name(id='reference_coordinate',ctx=ast.Load())
            elif replace_radius and name=='self.frame.radius_log':
                radius.append(ast.unparse(node))
                node=ast.parse("self.frame.definitions['logRref']+y",mode='eval').body
            return node
    Project().visit(fn)
    if len(changes)!=expected_coordinates or len(radius)!=int(replace_radius):
        raise ValueError('Original point method domain/metadata sites changed')
    # The reverse mapping certifies that all original mathematical operations
    # outside chart validation and the equivalent radius metadata are unchanged.
    restored=copy.deepcopy(fn)
    class Restore(ast.NodeTransformer):
        def visit_Call(self,node):
            self.generic_visit(node)
            if isinstance(node.func,ast.Name) and node.func.id=='reference_coordinate':
                node.func=ast.parse(changes.pop(0),mode='eval').body
            return node
        def visit_BinOp(self,node):
            if replace_radius and ast.dump(node)==ast.dump(ast.parse("self.frame.definitions['logRref']+y",mode='eval').body):
                return ast.parse(radius[0],mode='eval').body
            return self.generic_visit(node)
    Restore().visit(restored)
    if ast.dump(restored)!=ast.dump(original):raise ValueError('Original point math changed')
    env=dict(vars(module),reference_coordinate=reference_coordinate)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),
        '<original point method; reference domain and equivalent radius metadata>','exec'),env)
    return env[method]


class ReferenceFactoredInputs(point.OriginalO2FactoredPointInputs):
    mode='original_Rh_reference_closed_point_coefficients_with_directed_errors'
    def radial(self,y):
        y=reference_coordinate(y)
        if y not in self.radial_cache:
            c=self.ctx;iv=self.interval;cy=c.mpf(int(y.p))/int(y.q);iy=iv.mpf(int(y.p))/int(y.q)
            f=c.exp(cy/10);F=iv.exp(iy/10)
            values=dict(f=f,H=c.mpf(5)/8*f,D=c.mpf(5)/12*f*f,P=c.mpf(5)/2*f*f,a=c.mpf(4)/5)
            boxes=dict(f=F,H=iv.mpf(5)/8*F,D=iv.mpf(5)/12*F**2,P=iv.mpf(5)/2*F**2,a=iv.mpf(4)/5)
            self.radial_cache[y]=(values,boxes)
        return self.radial_cache[y]


ReferenceFactoredInputs.evaluate=project_method(point,'OriginalO2FactoredPointInputs','evaluate',1,True)


class ReferenceRadiusPhase:
    def __init__(self,original):self.original=original;self.family=original.family
    def evaluate(self,*,y,N):
        y=reference_coordinate(y);original=self.original.evaluate(y=0,N=N)
        c=base.MPIntervalContext();c.dps=max(260,self.original.dps+N.bit_length())
        rows=[]
        with mp.workdps(c.dps+40):
            shift=c.mpf(N)*c.mpf(int(y.p))/int(y.q)
            for box in original['true_original_phase_directed_boxes']:
                projection=base.radius.phase.ordinary_mod_one(c,c.mpf([box['lower'],box['upper']])+shift)
                rows.extend(base.radius.interval_record(v) for v in projection['boxes'])
            p=mp.mp.clone();p.dps=c.dps
            approximate=p.mpf(original['approximate_original_fractional_phase'])+p.mpf(N)*int(y.p)/int(y.q)
            approximate-=p.floor(approximate)
        return dict(original,original_y_exact=str(y),approximate_original_fractional_phase=approximate,
            true_original_phase_directed_boxes=rows,reference_phase_from_exact_N_offset=True,
            exact_source_phase='frac(N*(log(110/4)+14*logPstar+10*logCstar+1000+offset-hb*s_c/2))',
            reference_R_equals_Rref_exp_offset=True,positive_origin_budget_not_zeroed=True,
            phase_arithmetic_error_contract='Directed base phase plus exact rational N*offset; every new add/divide is directed',
            periodic_projection_full_period=any(ep(c.mpf([v['lower'],v['upper']]))==(0,1) for v in rows))


reference_query=project_method(base,'OriginalO2ConditionedPrimitives','query',1)


class OriginalReferencePointOracle:
    mode='original_Rh_reference_factored_true_point_and_radius_phase_C0_Z_enclosures'
    def __init__(self):
        self.owner=base.OriginalO2ConditionedPrimitives()
        old=self.owner.inputs
        inputs=ReferenceFactoredInputs.__new__(ReferenceFactoredInputs)
        inputs.__dict__=dict(old.__dict__,radial_cache={},point_cache={})
        self.owner.inputs=inputs;self.owner.radius=ReferenceRadiusPhase(self.owner.radius)
        self.family=self.owner.family;self.ctx=self.owner.ctx;self.recipe=reference_recipe_binding()
        self.hashes=dict(self.owner.hashes)
        receipt=json.loads((HERE/slow.RECEIPT).read_bytes())
        if not receipt.get('all_passed') or not receipt.get(slow.GATE) or receipt['source_family']!=self.family:
            raise ValueError('Accepted same-family original phase-held Z kernel required')
        for name,digest in {**receipt['input_hashes'],slow.RECEIPT:sha(slow.RECEIPT)}.items():
            if sha(name)!=digest:raise ValueError('Reference oracle ancestor changed: '+name)
            if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Original family closures disagree')
            self.hashes[name]=digest
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.cache={}
    def evaluate(self,*,offset,Z,N,bits=80):
        offset=reference_coordinate(offset);Z=point.pressure.exact_Z(Z)
        if type(bits) is not int or not 4<=bits<=256:raise ValueError('Inverse bits in[4,256] required')
        key=(offset,Z,N,bits)
        if key in self.cache:return self.cache[key]
        c=self.ctx;phase=self.owner.radius.evaluate(y=offset,N=N)
        with mp.workdps(c.dps+40):
            query=reference_query(self.owner,y=offset,Z=Z);kernel=query['kernel'];rows=[]
            for box in phase['true_original_phase_directed_boxes']:
                target=c.mpf([box['lower'],box['upper']]);inverse=kernel.evaluate(target,bits=bits)
                if inverse['status']!='enclosed':raise ArithmeticError('Refine original reference signed inputs')
                selected=inverse['selected_inverse']
                values,proof=slow.slow_values(kernel,query['roots'],selected['coordinate_interval'],selected['chart'])
                inverse.update(free_phase_parameter_not_spatial_phase=False,
                    original_common_N_and_radius_phase_bound=True)
                rows.append(dict(C0=inverse,slow_Z={k:v.record() for k,v in values.items()},
                    derivative_contract=proof))
        record=dict(source_family=self.family,mode=self.mode,chart='Rh_reference',
            original_offset_exact=str(offset),original_Z_exact=str(Z),explicit_candidate_N=N,
            original_factored_point_inputs=point.record_query(query['point']),
            actual_original_radius_phase=phase,actual_phase_held_Z_enclosures=rows,
            source_factor_basis=query['basis_contract'],numerical_arithmetic_ledger=query['ledger'],
            original_closed_reference_recipe=self.recipe,
            true_reference_C0_Z_source_and_primitives_installed=True,
            true_radius_phase_Z_exact_zero=True,closed_profiles_require_no_defining_quadrature=True,
            source_caps_or_midpoints_selected_as_field_values=False,
            full_17_chart_numeric_oracle_installed=False,actual_five_controls_installed=False,
            current_whole_N_selected=False,**dict.fromkeys(point.source.inertial.profiles.loop.OPEN,False))
        self.cache[key]=record;return record


def run():
    begin=time.monotonic();owner=OriginalReferencePointOracle();samples=[]
    for offset,Z,N in (('-5','-.37',7),('-2.337','0',7),('-2.337','.37',7),('0','.37',7)):
        row=owner.evaluate(offset=offset,Z=Z,N=N);samples.append(row)
        print('True reference C0/Z point:',offset,Z,N,row['actual_phase_held_Z_enclosures'][0]['derivative_contract']['branch'],flush=True)
    result=dict(**{GATE:True},source_family=owner.family,mode=owner.mode,
        actual_original_reference_point_queries=samples,reference_recipe=owner.recipe,
        original_inertial_pressure_inverse_and_Z_programs_reused=True,
        original_reference_to_O2_slope_source_join='offset=0 and original slope y=0 share f=1,H=5/8,D=5/12,P=5/2,a=4/5,R=Rref',
        no_original_slope_quadrature_or_ancestor_producer_reexecuted=True,
        true_reference_C0_Z_source_and_primitives_installed=True,
        full_17_chart_numeric_oracle_installed=False,actual_five_controls_installed=False,
        current_whole_N_selected=False,**dict.fromkeys(point.source.inertial.profiles.loop.OPEN,False),
        input_hashes=owner.hashes,execution_seconds=time.monotonic()-begin,
        scope='A second original source chart: closed Rh_reference E/V/a/b/p1/p2 C0/Z point coefficients and error budgets feed same selected scales, true candidate-N radius phase, original inverse/A-B and phase-held Z. All-chart integral/control/global-N and recursive corrected field remain open.')
    (HERE/NAME).write_text(json.dumps(base.encoded(result),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()
