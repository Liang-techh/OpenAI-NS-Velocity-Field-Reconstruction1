"""Original full background tensor and NS remainder on live current-Rp rows.

Replay the paper stress AST with all actual moments and absolute pressure.
Expand shared source products and apply exact FTC/theta laws before bounds.
Point decomposition does not certify a stress cone or a flat remainder.
"""
import ast
import copy
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
import gzip
import json
from pathlib import Path
import time

import sympy as s
import lei_ren_part1_paper_compliant_current_original_Rp_axial_vorticity as axial
import lei_ren_part1_paper_compliant_pulse_high_jets as pulse_high
from lei_ren_part1_paper_compliant_pulse_high_jets_check import source_identities as radial_source_identities
from lei_ren_part1_paper_compliant_pulse_end_physical_C2 import original_full_stress, full_physical_identities
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

HERE, PREFIX, sha, ends = axial.HERE, axial.PREFIX, axial.sha, axial.ends
signed, correlated, box, centered = axial.signed, axial.correlated, axial.box, axial.centered
NAME = PREFIX+'current_original_Rp_background_stress.json.gz'
RECEIPT = PREFIX+'current_original_Rp_background_stress_check.json'
GATES = ('current_original_Rp_full_paper_stress_actual_history_rows_installed',
    'current_original_Rp_completed_tensor_and_separate_remainder_point_installed',
    'current_original_Rp_source_bound_signed_background_momentum_decomposition_installed')
OPEN = axial.OPEN
Y, Z, DELTA, MU = s.symbols('y Z delta mu', real=True)
R = s.Symbol('R', positive=True)
CS, SN = s.symbols('cos_theta sin_theta', real=True)
LABELS = ('Utheta', 'Uz', 'Ur', 'pressure', 'Mtheta', 'Mz', 'Mtheta_z', 'Mztheta')
F = {name:s.Function(name)(Y, Z) for name in LABELS}


def dy(value):
    return s.diff(value, Y)+R*s.diff(value, R)


def dz(value, beta):
    return (beta*Z*value+(1-Z**2)*s.diff(value, Z)-2*Z*dy(value))/(1-DELTA*Z**2)


def dt(value, beta):
    return (-beta*value/2+(1-DELTA)*Z*s.diff(value, Z)/2+dy(value))/(1-DELTA*Z**2)


def reduce_source(value, theta_law=True):
    """Original cumulative FTC equations; independent P0 survives at y order 0."""
    density = dict(Mtheta=s.sqrt(2)*R**s.Rational(3,2)*F['Utheta'], Mz=R*F['Uz'],
        Mtheta_z=s.sqrt(2)*R**s.Rational(3,2)*F['Utheta']*F['Uz'],
        Mztheta=R*(F['Uz']**2-F['Utheta']**2/2), pressure=F['Utheta']**2/2)
    rules = {}
    for derivative in value.atoms(s.Derivative):
        if derivative.expr not in F.values():
            raise ValueError('Original source derivatives required')
        label = str(derivative.expr.func)
        powers = dict(derivative.variable_count); k, n = powers.get(Y, 0), powers.get(Z, 0)
        if k and label in density:
            result = density[label]
            for _ in range(k-1):result=dy(result)
            rules[derivative]=s.diff(result,Z,n)
    value=value.xreplace(rules)
    if theta_law:
        rules={}
        for derivative in value.atoms(s.Derivative):
            if derivative.expr==F['Utheta']:
                powers=dict(derivative.variable_count); k,n=powers.get(Y,0),powers.get(Z,0)
                if k:rules[derivative]=(-(s.Rational(1,2)+MU))**k*s.diff(F['Utheta'],Z,n)
        value=value.xreplace(rules)
    return value


@lru_cache(maxsize=1)
def source_operators():
    Ut,Uz,Ur,P=(F[name] for name in ('Utheta','Uz','Ur','pressure'))
    moments={key:F[name] for key,name in dict(theta='Mtheta',z='Mz',theta_z='Mtheta_z',z_theta='Mztheta',p='pressure').items()}
    actual,hashes=original_full_stress(dict(a=Ut,b=Uz,ay=dy(Ut),az=s.diff(Ut,Z),by=dy(Uz),bz=s.diff(Uz,Z),
        dt=DELTA,z=Z,R=R,root=s.sqrt(2*R),d=1-Z**2,L=1-DELTA*Z**2,m=moments,
        mz={key:s.diff(value,Z) for key,value in moments.items()},p=P,pz=s.diff(P,Z)))
    Tt,Tz=actual['Itheta']+actual['Stheta'],actual['Iz']+actual['Sz']
    bt=-2-DELTA; bv=-1-DELTA
    divt=s.sqrt(2/R)*(dy(Tt)+Tt); divz=s.sqrt(2/R)*(dy(Tz)+Tz/2)
    diagonal=s.sqrt(2*R)*dz(Tz,bt)
    er0=dt(Ur,-1)+s.sqrt(2/R)*Ur*dy(Ur)+Uz*dz(Ur,-1)-2/R*(dy(dy(Ur))-Ur/4)
    er1=-dz(dz(Ur,-1),DELTA-2)
    et=-dz(dz(Ut,bv),bv+DELTA-1); ez=-dz(dz(Uz,bv),bv+DELTA-1)
    rt=dt(Ut,bv)+s.sqrt(2/R)*Ur*(dy(Ut)+Ut/2)+Uz*dz(Ut,bv)-2/R*(dy(dy(Ut))-Ut/4)
    rz=dt(Uz,bv)+s.sqrt(2/R)*Ur*dy(Uz)+Uz*dz(Uz,bv)+dz(P,-2-2*DELTA)-2/R*dy(dy(Uz))
    # The actual recovery/FTC identities establish R_B = -div(T) + E.
    identities={}
    for name,expr in (('theta',rt+divt),('axial',rz+divz)):
        expr=expr.subs(Ur,actual['Ur']).doit()
        if s.cancel(s.expand(reduce_source(expr,theta_law=False)))!=0:
            raise ArithmeticError('Original full background decomposition differs: '+name)
        identities[name+'_momentum_decomposition']=True
    pressure_centrifugal=s.sqrt(2/R)*dy(P)-Ut**2/s.sqrt(2*R)
    if s.cancel(s.expand(reduce_source(pressure_centrifugal)))!=0:
        raise ArithmeticError('Original radial pressure/centrifugal cancellation differs')
    identities['radial_pressure_centrifugal_source_cancellation']=True
    zero=[]
    result=dict(cylindrical_stress=dict(r_theta=[(bt,Tt)],r_z=[(bt,Tz)],theta_theta=[(-2,diagonal)],rr=zero,theta_z=zero,zz=zero),
        cylindrical_divergence=dict(radial=zero,theta=[(-3-DELTA,divt)],axial=[(-3-DELTA,divz)]),
        cylindrical_remainder=dict(radial=[(-3,er0),(-3+2*DELTA,er1)],theta=[(-3+DELTA,et)],axial=[(-3+DELTA,ez)]),
        cylindrical_momentum_residual=dict(radial=[(-3,er0),(-3+2*DELTA,er1)],
            theta=[(-3-DELTA,rt),(-3+DELTA,et)],axial=[(-3-DELTA,rz),(-3+DELTA,ez)]),
        Cartesian_stress=dict(xx=[(bt,-2*CS*SN*Tt),(-2,SN**2*diagonal)],
            xy=[(bt,(CS**2-SN**2)*Tt),(-2,-CS*SN*diagonal)],xz=[(bt,CS*Tz)],
            yy=[(bt,2*CS*SN*Tt),(-2,CS**2*diagonal)],yz=[(bt,SN*Tz)],zz=zero))
    for output,source in (('Cartesian_divergence','cylindrical_divergence'),('Cartesian_remainder','cylindrical_remainder'),
            ('Cartesian_momentum_residual','cylindrical_momentum_residual')):
        v=result[source]
        result[output]=dict(x=[(b,CS*e) for b,e in v['radial']]+[(b,-SN*e) for b,e in v['theta']],
            y=[(b,SN*e) for b,e in v['radial']]+[(b,CS*e) for b,e in v['theta']],z=list(v['axial']))
    definitions={section:{name:[dict(lambda_exponent=s.srepr(s.sympify(b)),source_expression=s.srepr(e)) for b,e in parts]
        for name,parts in rows.items()} for section,rows in result.items()}
    return result,definitions,identities,hashes


def source_symbols(expression):
    replacements={}; rows={}
    for atom in expression.atoms(s.Derivative):
        if atom.expr not in F.values():raise ValueError('Actual mixed source primitive required')
        powers=dict(atom.variable_count); k,n=int(powers.get(Y,0)),int(powers.get(Z,0))
        if set(powers)-{Y,Z} or k+n>4:raise ValueError('Available ordinary source4 derivative required')
        label=str(atom.expr.func); symbol=s.Symbol('S_'+label+'_'+str(k)+'_'+str(n))
        replacements[atom]=symbol; rows[symbol]=(label,k,n)
    expression=expression.xreplace(replacements)
    replacements={}
    for name,atom in F.items():
        if expression.has(atom):
            symbol=s.Symbol('S_'+name+'_0_0'); replacements[atom]=symbol; rows[symbol]=(name,0,0)
    return expression.xreplace(replacements),rows


def operator_definitions(operators):
    return {section:{name:[dict(lambda_exponent=s.srepr(s.sympify(beta)),source_expression=s.srepr(expr))
        for beta,expr in parts] for name,parts in rows.items()} for section,rows in operators.items()}


def radial_recovery_binding():
    """Bind the live normalized recovery to the original paper equation."""
    name=PREFIX+'pulse_high_jets.py'
    tree=ast.parse((HERE/name).read_bytes())
    cls=next(node for node in tree.body if isinstance(node,ast.ClassDef) and node.name=='CompliantPulseHighJets')
    fn=next(node for node in cls.body if isinstance(node,ast.FunctionDef) and node.name=='radial')
    assignments={target.id:node.value for node in ast.walk(fn) if isinstance(node,ast.Assign)
        for target in node.targets if isinstance(target,ast.Name)}
    expected=dict(z='IntervalTaylor(c, [Z, 1]+[0]*max(0, order-1))',
        q='IntervalTaylor(c, [1+Z**2, 2*Z, 1]+[0]*max(0, order-2)).truncate(order)',
        d='IntervalTaylor(c, [1-Z**2, -2*Z, -1]+[0]*max(0, order-2)).truncate(order)',
        L='IntervalTaylor(c, [1-self.delta*Z**2, -2*self.delta*Z, -self.delta]+[0]*max(0, order-2)).truncate(order)',
        b='B.truncate(order)',m='m1.truncate(order)',mz='derivative(m1)',
        numerator='z*b*2-z*m*(1-self.delta)-d*(mz-z*m/q*2)')
    for key,value in expected.items():
        if ast.dump(assignments[key])!=ast.dump(ast.parse(value,mode='eval').body):
            raise ValueError('Actual original radial recovery differs: '+key)
    returned=[node.value for node in ast.walk(fn) if isinstance(node,ast.Return)]
    if len(returned)!=1 or ast.dump(returned[0])!=ast.dump(ast.parse('numerator/L',mode='eval').body):
        raise ValueError('Original radial recovery denominator differs')
    mixed_name=PREFIX+'pulse_mixed_C4.py'
    tree=ast.parse((HERE/mixed_name).read_bytes())
    fn=next(node for node in tree.body if isinstance(node,ast.FunctionDef) and node.name=='transport_mixed')
    assignments={target.id:node.value for node in ast.walk(fn) if isinstance(node,ast.Assign)
        for target in node.targets if isinstance(target,ast.Name)}
    recovery=assignments['A']
    if ast.dump(recovery)!=ast.dump(ast.parse('[field.radial(Z,b,m) for b,m in zip(Brows,m1)]',mode='eval').body):
        raise ValueError('Actual radial mixed rows must consume the original recovery')
    physical=assignments['physical']
    radial=next(value for key,value in zip(physical.keys,physical.values)
        if isinstance(key,ast.Constant) and key.value==correlated.physical.mixed.UR)
    if ast.dump(radial)!=ast.dump(ast.parse('[u.truncate(4)*binomial_product(A,-mu,k) for k in range(5)]',mode='eval').body):
        raise ValueError('Original radial product derivatives differ')
    proof=radial_source_identities()
    keys=('production_paper_3_9_numerator','radial_log_derivative_from_paper_3_8','physical_radial_Z_product_rule')
    if not all(proof[key] for key in keys):raise ValueError('Original radial source identities required')
    return dict(source_file=name,source_sha256=sha(name),mixed_source_file=mixed_name,mixed_source_sha256=sha(mixed_name),
        proof_source_file=PREFIX+'pulse_high_jets_check.py',proof_source_sha256=sha(PREFIX+'pulse_high_jets_check.py'),
        normalized_recovery_AST={key:ast.dump(assignments_value) for key,assignments_value in
            ((key,ast.parse(value,mode='eval').body) for key,value in expected.items())},
        recovery_call_AST=ast.dump(recovery),radial_product_AST=ast.dump(radial),
        original_radial_source_identities={key:proof[key] for key in keys})


@dataclass(frozen=True,eq=False)
class OriginalBackgroundStress:
    chart:str


class CurrentOriginalRpBackgroundStress:
    @source_precision
    def __init__(self,before,require_checked=True):
        if type(before) is not axial.CurrentOriginalRpAxialVorticity or not before.acceptance_loaded:
            raise ValueError('Accepted current original source/axial curl owner required')
        self.before,self.differential,self.product,self.arithmetic=before,before.before,before.product,before.arithmetic
        self.graph,self.ctx,self.family_record=before.graph,before.ctx,before.family_record
        self.operators,self.definitions,self.identities,hashes=source_operators()
        self._definitions=copy.deepcopy(self.definitions)
        self.radial_law=radial_recovery_binding()
        self._radial_callback=pulse_high.CompliantPulseHighJets.radial
        self.mathematics=full_physical_identities()
        if not all(self.mathematics['identities'].values()):raise ValueError('Original full physical source theorem required')
        self.hashes=dict(before.hashes)
        for name,digest in {**hashes,**self.mathematics['input_hashes']}.items():
            if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Original stress source conflict '+name)
            self.hashes[name]=digest
        for name in (axial.NAME,axial.RECEIPT,Path(__file__).name,PREFIX+'pulse_end_physical_C2.py'):
            self.hashes[name]=sha(name)
        self.acceptance_loaded=False;self._fields={}
        self.assert_graph()
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not all(receipt[key] for key in GATES) or any(receipt[key] for key in OPEN) or \
                    receipt['source_family']!=self.family_record or receipt['original_operator_definitions']!=self.definitions or \
                    receipt['actual_radial_recovery_source_binding']!=self.radial_law:
                raise ValueError('Original background stress receipt/source/scope differs')
            for name in (Path(__file__).name,Path(__file__).stem+'_check.py',NAME):
                if receipt['input_hashes'].get(name)!=sha(name):raise ValueError('Unbound background source '+name)
            box.pulse.radius.post.selected.inlet.add_hashes(self.hashes,receipt['input_hashes'])
            self.acceptance_loaded=True

    @source_precision
    def assert_graph(self):
        self.before.assert_graph()
        checks=dict(accepted_original_differential_and_theta_source=self.before.acceptance_loaded,
            one_live_original_graph_and_context=self.graph is self.product.graph and self.ctx is self.product.ctx,
            full_stress_source_definitions_unchanged=self.definitions==self._definitions,
            actual_source_operator_expressions_unchanged=operator_definitions(self.operators)==self._definitions,
            full_original_decomposition_identities=all(self.identities.values()),
            full_original_physical_operator_identities=all(self.mathematics['identities'].values()),
            original_radial_recovery_AST_and_identities_unchanged=radial_recovery_binding()==self.radial_law,
            actual_live_radial_recovery_callback=self.before.refined.pulse.radial.__func__ is self._radial_callback
                is pulse_high.CompliantPulseHighJets.radial,
            actual_radial_derivative_callable=self._radial_callback.__globals__['derivative'] is pulse_high.derivative,
            actual_radial_Taylor_type=self._radial_callback.__globals__['IntervalTaylor'] is pulse_high.IntervalTaylor)
        if not all(checks.values()):raise ValueError('Original background stress source differs: '+str(checks))
        return checks

    def _validate_source(self,source):
        evidence=source.get('velocity_source_evidence',{})
        if 'original_forward_source_packet' not in source or \
                evidence.get('source')!='actual_current_pulse_physical_radial_product_rows' or \
                evidence.get('physical_radial_prefactors_already_differentiated') is not True:
            raise ValueError('Actual original radial recovery source packet required')

    def _exact(self,expr,coordinates):
        g=self.graph
        if expr==Z:return coordinates['Z']
        if expr==DELTA:return self.before.delta
        if expr==MU:return self.before.mu
        if expr==CS:return coordinates['cosine']
        if expr==SN:return coordinates['sine']
        if expr.is_Rational:return g.constant(Fraction(int(expr.p),int(expr.q)))
        if expr.is_Add:return g.add(*(self._exact(v,coordinates) for v in expr.args))
        if expr.is_Mul:return g.mul(*(self._exact(v,coordinates) for v in expr.args))
        if expr.is_Pow and expr.args[1].is_Integer:
            base=self._exact(expr.args[0],coordinates); exponent=int(expr.args[1])
            return g.mul(*([base]*exponent)) if exponent>=0 else g.node('exact_operator_integer_power',
                base=base.node,exponent=exponent,nonzero_certificate='original 1-delta*Z^2 positive')
        if expr.is_Pow and expr.args[0].is_positive and expr.args[0].is_Rational and expr.args[1].is_Rational:
            return g.unary('exp',g.mul(self._exact(expr.args[1],coordinates),g.unary('log',self._exact(expr.args[0],coordinates))))
        raise ValueError('Closed exact source operator coefficient required: '+str(expr))

    def _numeric(self,expr,bindings):
        if expr in bindings:return bindings[expr]
        c=self.ctx
        if expr.is_Rational:return c.mpf(int(expr.p))/int(expr.q)
        if expr.is_Add:return sum((self._numeric(v,bindings) for v in expr.args),c.mpf(0))
        if expr.is_Mul:
            result=c.mpf(1)
            for value in expr.args:result*=self._numeric(value,bindings)
            return result
        if expr.is_Pow and expr.args[1].is_Integer:return self._numeric(expr.args[0],bindings)**int(expr.args[1])
        if expr.is_Pow and expr.args[0].is_positive and expr.args[0].is_Rational and expr.args[1].is_Rational:
            return c.exp(self._numeric(expr.args[1],bindings)*c.ln(self._numeric(expr.args[0],bindings)))
        raise ValueError('Directed exact coefficient expression required: '+str(expr))

    def _enclose(self,parts,source,reader,coordinates,target):
        combined={}
        for beta,expr in parts:
            beta=s.sympify(beta)
            combined[beta]=combined.get(beta,s.Integer(0))+expr
        primitive_terms={}; total_symbolic_terms=0
        for beta,expr in combined.items():
            expression,variables=source_symbols(reduce_source(expr))
            for term in s.Add.make_args(s.expand(expression)):
                if term==0:continue
                powers=term.as_powers_dict(); primitive=[]; coefficient=term
                for symbol,key in sorted(variables.items(),key=lambda value:str(value[0])):
                    exponent=powers.get(symbol,s.Integer(0))
                    if not exponent.is_Integer or exponent<0:raise ValueError('Shared source monomial required')
                    if exponent:
                        primitive.append((*key,int(exponent))); coefficient/=symbol**exponent
                coefficient=s.cancel(coefficient)
                # cancel()/expand() can leave R inside a rational denominator.
                # Extract its exact homogeneous degree before creating scales.
                rp=s.cancel(R*s.diff(coefficient,R)/coefficient)
                if not rp.is_Rational:raise ValueError('Exact homogeneous original radial scale power required')
                coefficient=s.simplify(coefficient/R**rp)
                if coefficient.free_symbols-{Z,DELTA,MU,CS,SN}:raise ValueError('Unbound primitive coefficient '+str(coefficient))
                key=(beta,Fraction(int(rp.p),int(rp.q)),tuple(primitive))
                primitive_terms[key]=primitive_terms.get(key,s.Integer(0))+coefficient
                total_symbolic_terms+=1
        angle=reader.at(coordinates['theta'])
        numeric={Z:reader.at(coordinates['Z']),DELTA:reader.at(self.before.delta),MU:reader.at(self.before.mu),
            CS:self.ctx.cos(angle),SN:self.ctx.sin(angle)}
        groups={};ledger=[];cancelled=0
        for (beta,rp,primitive),coefficient in primitive_terms.items():
            coefficient=s.cancel(coefficient)
            if coefficient==0:cancelled+=1;continue
            log=self.graph.add(self.graph.mul(self.graph.constant(rp),coordinates['logR']),
                self.graph.mul(self._exact(beta,coordinates),coordinates['loglambda']))
            value=self.ctx.mpf(1); typed_units={'R':rp};dependencies=[]
            for label,k,n,exponent in primitive:
                row=source['log_radius_mixed_rows'][label]['y%d_Z%d'%(k,n)]
                if type(row) is not correlated.physical.mixed.FactorizedMixedSourceRow or row.name!=label or \
                        row.derivative!=(k,n) or row.powers[-1]!=0 or row.coefficients.order!=0 or row.coefficients.ctx is not self.product.amplitude.ctx:
                    raise ValueError('Actual ordinary shared source primitive required')
                value*=self.ctx.mpf(row.coefficients[0])**exponent
                log=self.graph.add(log,self.graph.mul(self.graph.constant(exponent),self.graph.add(*(ref for _,ref in row.log_scale_parts))))
                for name,power in zip(row.source_units,row.powers):typed_units[name]=typed_units.get(name,Fraction(0))+exponent*power
                dependencies.append(dict(source=label,ordinary_derivative=(k,n),power=exponent,
                    actual_source_units=row.source_units,
                    actual_source_powers=[signed.rational_record(q) for q in row.powers]))
            directed_operator=self._numeric(coefficient,numeric);exact_operator=self._exact(coefficient,coordinates)
            value*=directed_operator
            polynomial=reader.polynomial(log)
            key=(tuple(sorted(polynomial.items())),tuple(sorted(typed_units.items())),s.srepr(beta))
            if key not in groups:
                canonical=signed.polynomial_function(self.graph,polynomial)
                groups[key]=dict(log_function=canonical,log_bound=reader.at(canonical),coefficient=self.ctx.mpf(0),proofs=[])
            groups[key]['coefficient']+=value
            ledger.append(dict(actual_shared_source_product=dependencies,exact_operator=exact_operator.node,
                exact_operator_expression=s.srepr(coefficient),directed_operator_enclosure=directed_operator,
                radial_power=signed.rational_record(rp),lambda_exponent=s.srepr(beta),
                exact_product_log_scale=log.node,signed_product_coefficient=value,
                shared_primitive_monomial_collected_before_enclosure=True))
        result=signed.enclose_factored_sum(self.graph,self.ctx,reader,list(groups.values()),target)
        if not result['exact_zero_enclosure'] and result.get('common_scale_coefficient_enclosure') is not None:
            center=self.product.center(reader,box.pulse.radius.FunctionRef(self.graph,result['exact_reference_log_scale_function']))
            accuracy=centered.centered_relative_budget(self.ctx,result,target,center)
        else:
            accuracy=dict(ordinary_numeric_relative_width_satisfied=result['exact_zero_enclosure'],
                ordinary_numeric_delivery_target_satisfied=result['exact_zero_enclosure'])
        result.update(physical_accuracy=accuracy,actual_source_product_ledger=ledger,
            expanded_source_monomial_count=total_symbolic_terms,exact_zero_source_monomials_cancelled=cancelled,
            nonlinear_shared_source_products_collected_before_interval_bounds=True,
            background_cone_flat_global_and_recursion_gates_unpromoted=True)
        return result

    def _source_dag(self,views):
        pending=[]; snapshots={}
        for rows in views.values():
            for result in rows.values():
                pending.extend(item[key] for item in result['actual_source_product_ledger']
                    for key in ('exact_operator','exact_product_log_scale'))
                pending.extend(item[key] for item in result.get('ratio_terms',[])
                    for key in ('exact_ratio_function','exact_log_ratio_function'))
        while pending:
            node=pending.pop()
            if node in snapshots:continue
            data=self.graph.nodes[node];snapshots[node]=copy.deepcopy(data);op=data['operation']
            if op in ('sum','product'):pending.extend(data['arguments'])
            elif op in ('negative','analytic_unary'):pending.append(data['argument'])
            elif op=='positive_quotient':pending.extend((data['numerator'],data['denominator']))
            elif op=='exact_operator_integer_power':pending.append(data['base'])
        return snapshots

    @source_precision
    def evaluate(self,field,relative_width_target='1/1000'):
        self.assert_graph();record=self.differential._require(field)
        if record['chart']!='pulse_exit':raise ValueError('Current source-bound stress adapter admits pulse_exit')
        target=Fraction(relative_width_target)
        if not 0<target<1:raise ValueError('Exact relative-width target required')
        delivery=record['delivery'];source=self.product.before.source(delivery)
        self._validate_source(source)
        reader=self.product.reader(delivery)
        _,_,request=self.product.amplitude.before._validate_delivery(delivery)
        views={section:{name:self._enclose(parts,source,reader,request.forward_coordinates,target)
            for name,parts in rows.items()} for section,rows in self.operators.items()}
        value=OriginalBackgroundStress(record['chart'])
        snapshot=self.product.amplitude.before._fingerprint(box.report(source))
        self._fields[id(value)]=(value,field,source,snapshot,views,self._source_dag(views))
        return value

    @source_precision
    def report(self,value):
        self.assert_graph();entry=self._fields.get(id(value))
        if type(value) is not OriginalBackgroundStress or entry is None or entry[0] is not value or value.chart!=entry[1].chart:
            raise ValueError('Live background stress field issued by this owner required')
        field=self.differential._require(entry[1]);source=self.product.before.source(field['delivery'])
        self._validate_source(source)
        if source is not entry[2] or self.product.amplitude.before._fingerprint(box.report(source))!=entry[3]:
            raise ValueError('Original complete moment/pressure source changed')
        if any(self.graph.nodes[node]!=data for node,data in entry[5].items()):
            raise ValueError('Original background source product/operator DAG changed')
        return dict(chart=value.chart,source_family=copy.deepcopy(self.family_record),**copy.deepcopy(entry[4]),
            original_momentum_decomposition_identities=copy.deepcopy(self.identities),original_full_physical_theorem=copy.deepcopy(self.mathematics),
            actual_independent_P0_and_remaining_absolute_pressure_retained=True,
            incoming_full_moments_and_nonlinear_cross_products_retained=True,
            completed_diagonal_r_partial_z_Tz_retained=True,
            actual_radial_recovery_source_binding=copy.deepcopy(self.radial_law),
            convention='R_B=-div(T_B)+E_B',original_viscosity=1,
            regional_remainder_not_claimed_flat=True,unrestricted_physical_point_API=False,
            full_certified_physical_accuracy=False,**dict.fromkeys(GATES,self.acceptance_loaded),**dict.fromkeys(OPEN,False))


@source_precision
def run(before,fields):
    began=time.monotonic();owner=CurrentOriginalRpBackgroundStress(before,require_checked=False)
    values={name:owner.evaluate(field) for name,field in fields.items()}
    views={name:owner.report(value) for name,value in values.items()}
    result=dict(source_family=owner.family_record,actual_background_stress=views,source_assertions=owner.assert_graph(),
        actual_radial_recovery_source_binding=owner.radial_law,
        original_operator_definitions=owner.definitions,input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
        **dict.fromkeys(GATES+OPEN,False))
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(correlated.report(result),separators=(',',':'))+'\n').encode(),mtime=0))
    return owner,values,views
