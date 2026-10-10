"""Direct current raw similarity rows -> physical Cartesian/time rows.

The incoming rows already differentiate the full similarity components.
Keep their exact source factors and signed ordinary derivative enclosures;
append only the original coordinate-operator factors. No old physical
provider, rounded absolute radius, or point coefficient is constructed.
This off-axis chart-parametric interface is not an arbitrary x/y/z/t oracle.
"""
from dataclasses import dataclass
from fractions import Fraction
import gzip
import json
from pathlib import Path
import time
import sympy as s

import lei_ren_part1_paper_compliant_current_original_Rp_postpulse_mixed_seams as chain
from lei_ren_part1_paper_compliant_cartesian_field import (
    cartesian_templates, angular_polynomial, CS, SN, INDICES)
from lei_ren_part1_paper_compliant_pulse_physical_bounds import (
    physical_operators, interval_expression, ZSYM, DSYM, BSYM, UR, UT, UZ, P)
from lei_ren_part1_paper_compliant_global_physical_assembly_check import source_identities
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

mixed=chain.mixed
HERE,PREFIX,sha=chain.HERE,chain.PREFIX,chain.sha
NAME=PREFIX+'current_original_Rp_physical_source_map.json.gz'
RECEIPT=PREFIX+'current_original_Rp_physical_source_map_check.json'
GATES=('current_original_Rp_same_graph_physical_coordinate_functions_installed',
       'current_original_Rp_scaled_Cartesian_spatial4_source_rows_installed',
       'current_original_Rp_scaled_fixed_x_time1_source_rows_installed')
OPEN=tuple(dict.fromkeys(chain.OPEN+('full_point_physical_field_evaluation',
    'arbitrary_physical_coordinate_inverse_installed','physical_axis_limit_installed',
    'physical_energy_integral_certified','full_background_NS_validation',
    'admissible_stress_lift_constructed','independently_bounded_flat_remainder','temporal_recursion')))
COMPONENTS=('ux','uy','uz','p')
SOURCE={UR:'Ur',UT:'Utheta',UZ:'Uz',P:'pressure'}


def operator_definitions():
    """Exact pure operator definitions; no numerical provider construction."""
    return dict(cylindrical={str(a)+','+str(b):{str(k)+','+str(n):s.srepr(v)
        for (k,n),v in row.items()} for (a,b),row in physical_operators().items()},
        Cartesian={','.join(map(str,key)):{'|'.join(map(str,index)):s.srepr(v)
        for index,v in row.items()} for key,row in cartesian_templates().items()},
        time_coefficients=['-beta/(2*(1-delta*Z**2))',
            '(1-delta)*Z/(2*(1-delta*Z**2))','1/(1-delta*Z**2)'],
        input_is_full_ordinary_similarity_derivative=True,
        no_second_factorial_or_amplitude_derivative=True,
        current_Ur_already_has_original_sqrt_R_over_2_factor=True)


def exact_expression(g,expression,bindings):
    if expression in bindings:return bindings[expression]
    if expression.is_Rational:return g.constant(Fraction(int(expression.p),int(expression.q)))
    if expression.is_Add:return g.add(*(exact_expression(g,v,bindings) for v in expression.args))
    if expression.is_Mul:return g.mul(*(exact_expression(g,v,bindings) for v in expression.args))
    if expression.is_Pow and expression.args[1].is_Integer:
        base=exact_expression(g,expression.args[0],bindings);power=int(expression.args[1])
        if power>=0:return g.mul(*([base]*power))
        # cancel() may write a denominator as +/- L**b. The sign is part
        # of the expression, so do not label that base a positive quotient.
        return g.node('exact_operator_integer_power',base=base.node,exponent=power,
            nonzero_certificate='original operator denominator; 1-delta*Z^2 > 0')
    raise ValueError('Unsupported exact physical operator: '+str(expression))


@dataclass(frozen=True)
class PhysicalSourceTerm:
    source_label: str
    source_row: mixed.FactorizedMixedSourceRow
    operator_function: object
    operator_coefficient: mixed.IntervalTaylor
    radial_power: Fraction
    lambda_exponent: object
    log_scale_parts: tuple
    signed_coefficient: mixed.IntervalTaylor


@dataclass(frozen=True)
class PhysicalSourceRow:
    component: str
    derivative: tuple
    terms: tuple

    @source_precision
    def groups(self):
        """Add signed terms only when all exact scale/unit nodes coincide."""
        grouped={}
        for term in self.terms:
            key=(term.source_row.source_units,term.source_row.powers,
                 tuple((name,value.node) for name,value in term.log_scale_parts))
            if key not in grouped:
                grouped[key]=dict(scale_units=term.source_row.source_units,
                    exact_source_scale_powers=term.source_row.powers,
                    log_scale_parts=term.log_scale_parts,signed_coefficient=term.signed_coefficient,
                    contributors=[term])
            else:
                grouped[key]['signed_coefficient']+=term.signed_coefficient
                grouped[key]['contributors'].append(term)
        return tuple(grouped.values())

    @source_precision
    def report(self):
        rows=[]
        for group in self.groups():
            value=group['signed_coefficient'];lo,hi=mixed.pulse.radius.post.selected.inlet.endpoints(value[0])
            rows.append(dict(scale_units=group['scale_units'],
                exact_source_scale_powers=[dict(numerator=v.numerator,denominator=v.denominator)
                    for v in group['exact_source_scale_powers']],
                positive_scale_log_parts={name:ref.node for name,ref in group['log_scale_parts']},
                signed_scaled_sum_enclosure=value,
                directed_scaled_sum_width_bound=value.ctx.mpf(hi)-value.ctx.mpf(lo),
                original_source_terms=[dict(source_label=t.source_label,
                    ordinary_similarity_derivative=t.source_row.derivative,
                    exact_operator_coefficient_function=t.operator_function.node,
                    exact_physical_lambda_exponent_function=t.lambda_exponent.node)
                    for t in group['contributors']]))
        return dict(component=self.component,physical_derivative=self.derivative,
            signed_exact_scale_groups=rows,source_term_count=len(self.terms),
            full_source_derivatives_not_redifferentiated=True,
            coefficient_enclosures_not_selected_as_point_values=True,
            numerical_absolute_point_value_installed=False)


class CurrentOriginalRpPhysicalSourceMap:
    @source_precision
    def __init__(self,before=None,require_checked=True):
        self.before=before if before is not None else chain.CurrentOriginalRpPostpulseMixedSeams()
        if type(self.before) is not chain.CurrentOriginalRpPostpulseMixedSeams or not self.before.acceptance_loaded:
            raise ValueError('Accepted typed current fourteen-interface owner required')
        self.transport=self.before.transport;self.owner=self.transport.owner;self.radius=self.before.radius
        self.graph=self.before.graph;self.ctx=self.before.ctx;self.family_record=self.before.family_record
        self.delta_function=self.radius.functions['delta'];self.delta=self.owner.selected.seed.future.delta
        self.original_N_definition=self.radius.frame.bridge.leading.bridge.source_integer_and_Rc_binding['original_N_definition']
        self.hashes=dict(self.before.hashes);self.definitions=operator_definitions()
        canonical_name=PREFIX+'global_physical_assembly_check.json'
        canonical=json.loads((HERE/canonical_name).read_bytes())
        self.canonical_source_identities=source_identities()
        if not canonical['all_passed'] or canonical['source_and_divergence_identities']!=self.canonical_source_identities:
            raise ValueError('Accepted unchanged original linear physical identities required')
        if (canonical['actual_five_defect_family_sha256'],canonical['implicit_source_sha256'])!=(self.before.post.family,self.before.post.source):
            raise ValueError('Canonical physical operator family differs')
        self.canonical_fixture=canonical['independent_implicit_physical_fixture']
        if not self.canonical_fixture['passed'] or self.canonical_fixture['independent_implicit_root_cartesian_derivatives']!=140 \
                or self.canonical_fixture['independent_fixed_x_time_derivatives']!=4:
            raise ValueError('Accepted independent original implicit-coordinate fixture required')
        mixed.pulse.radius.post.selected.inlet.add_hashes(self.hashes,canonical['input_hashes'])
        for name in (chain.NAME,chain.RECEIPT,Path(__file__).name,canonical_name,
                PREFIX+'cartesian_field.py',PREFIX+'pulse_physical_bounds.py',
                PREFIX+'global_physical_assembly.py',PREFIX+'global_physical_assembly_check.py'):
            self.hashes[name]=sha(name)
        self.acceptance_loaded=False;self.call_trace=[];self.assert_graph()
        if require_checked:
            record=json.loads((HERE/RECEIPT).read_bytes())
            if not record['all_passed'] or not all(record[k] for k in GATES) or any(record[k] for k in OPEN):
                raise ValueError('Current scaled physical map receipt/scope differs')
            if record['source_family']!=self.family_record or record['original_linear_operator_definitions']!=self.definitions:
                raise ValueError('Current physical family or exact operator definitions differ')
            mixed.pulse.radius.post.selected.inlet.add_hashes(self.hashes,record['input_hashes'])
            self.acceptance_loaded=True

    def assert_graph(self):
        g=self.graph;f=self.radius.functions
        expected_delta=g.unary('exp',g.add(g.mul(g.constant(-4),f['logP']),g.constant(-30)))
        ends=mixed.pulse.radius.post.selected.inlet.endpoints
        lo,hi=ends(self.delta)
        result=dict(accepted_same_current_fourteen_interface_owner=self.before.acceptance_loaded
                and all(self.before.assert_graph().values()),
            same_actual_mixed_transport_and_radius_graph=self.transport is self.before.transport
                and self.graph is self.transport.graph is self.radius.graph,
            original_exact_delta_function_not_an_enclosure=self.delta_function==expected_delta,
            same_original_delta_enclosure_context=self.delta.ctx is self.ctx
                and self.delta._mpi_==self.radius.before.before.pulse.delta._mpi_ and 0<lo<=hi<1,
            exact_original_finite_correction_integer_retained=self.original_N_definition==
                self.radius.frame.bridge.leading.bridge.source_integer_and_Rc_binding['original_N_definition']
                and self.original_N_definition['definition']=='N=2^(2^J)'
                and self.original_N_definition['exact_positive_integer'],
            original_cartesian_moving_basis_and_cylindrical_operators=self.definitions==operator_definitions(),
            inherited_analytic_P0_and_complete_Gamma_future=all(self.before.assert_graph().values()))
        if not all(result.values()):raise ValueError('Current physical source graph differs: '+str(result))
        return result

    def source_view(self,view,Z):
        """Validate a live typed observation; never decode receipt rows.

        evaluate() obtains this observation from the accepted owner. The
        optional observed-view run path reuses caller-owned live objects
        from that same current runtime, with scale/schema/context guards.
        It is not a serialized point-value import facility.
        """
        self.assert_graph();chart=view['chart']
        if chart not in mixed.CHARTS:raise ValueError('Actual current Rp chart required')
        z=mixed.pulse.radius.exact_coordinate(Z)
        if not -1<z<1:raise ValueError('Strict interior Z required for finite physical coordinates')
        c=self.ctx;zc=c.mpf(z.numerator)/z.denominator
        if type(view['Z']) is not type(zc) or view['Z'].ctx is not c \
                or (view['Z']-zc).a>0 or (view['Z']-zc).b<0:
            raise ValueError('Same exact axial argument and current context required')
        geometry=view['geometry'];coordinate=geometry.get('exact_native_coordinate')
        if coordinate is None:raise ValueError('This chart point interface requires an exact rational native coordinate')
        value=Fraction(coordinate['numerator'],coordinate['denominator'])
        if geometry!=self.radius.geometry(chart,value):raise ValueError('Exact current native radius geometry required')
        if not view['velocity_source_evidence'] or not view['radial_velocity_source_order4_retained']:
            raise ValueError('Actual full radial/axial similarity source rows required')
        for name in SOURCE.values():
            rows=view['log_radius_mixed_rows'][name]
            base=view['original_factorized_values'][name]
            if len(rows)!=15:raise ValueError('Complete ordinary mixed4 grid required: '+name)
            for k in range(5):
                for n in range(5-k):
                    row=rows['y%d_Z%d'%(k,n)]
                    if type(row) is not mixed.FactorizedMixedSourceRow or row.name!=name \
                            or row.derivative!=(k,n) or row.powers[-1]!=0 \
                            or type(row.coefficients) is not mixed.IntervalTaylor \
                            or row.coefficients.order!=0 or row.coefficients.ctx is not c \
                            or any(ref.graph is not self.graph for _,ref in row.log_scale_parts):
                        raise ValueError('Typed ordinary full-source row required: '+name)
                    powers=(tuple(map(Fraction,(0,2,2 if chart in mixed.pulse.PULSE else 0)))
                        if name=='pressure' and k else base.powers)
                    expected=self.transport.factor(name,chart,value,geometry,powers,row.coefficients,k,n,False)
                    if row.powers!=expected.powers or row.source_units!=expected.source_units or row.log_scale_parts!=expected.log_scale_parts:
                        raise ValueError('Original source units and exact lazy scales required')
            if rows['y0_Z0'].powers[:-1]!=base.powers:
                raise ValueError('Original absolute source-value units required')
        return value,z,zc

    @source_precision
    def coordinates(self,geometry,Z,log_tau,theta):
        g=self.graph;ref=mixed.pulse.radius.FunctionRef
        z=g.constant(Z);angle=mixed.pulse.radius.exact_coordinate(theta)
        if type(log_tau) is ref:
            if log_tau.graph is not g:raise ValueError('Same current graph log-time function required')
            lt=log_tau
        else:lt=g.constant(mixed.pulse.radius.exact_coordinate(log_tau))
        d=g.sub(g.one,g.mul(z,z));loglambda=g.mul(g.constant('1/2'),g.sub(lt,g.unary('log',d)))
        lr=ref(g,geometry['logR']);half=g.constant('1/2');delta=self.delta_function
        logr=g.add(loglambda,g.mul(half,g.add(g.unary('log',g.constant(2)),lr)))
        r=g.unary('exp',logr);theta_ref=g.constant(angle)
        cs=g.unary('cos',theta_ref);sn=g.unary('sin',theta_ref)
        axial=g.mul(z,g.unary('exp',g.mul(g.sub(g.one,delta),loglambda)))
        return dict(log_tau=lt,loglambda=loglambda,lambda_value=g.unary('exp',loglambda),
            logR=lr,R=ref(g,geometry['R']),log_r=logr,r=r,theta=theta_ref,cosine=cs,sine=sn,
            x=g.mul(r,cs),y=g.mul(r,sn),z=axial,t=g.sub(g.one,g.unary('exp',lt)),
            tau=g.unary('exp',lt),Z=z,delta=delta,
            L=g.sub(g.one,g.mul(delta,z,z)))

    def beta(self,label):
        g=self.graph;d=self.delta_function
        if label==UR:return g.constant(-1),self.ctx.mpf(-1)
        if label==P:return g.sub(g.constant(-2),g.mul(g.constant(2),d)),-2-2*self.delta
        return g.sub(g.constant(-1),d),-1-self.delta

    def term(self,view,label,index,operator,coefficient,radial_power,gamma,coordinates):
        g=self.graph;row=view['log_radius_mixed_rows'][SOURCE[label]]['y%d_Z%d'%index]
        scale=row.log_scale_parts+(('physical_radial',g.mul(g.constant(radial_power),coordinates['logR'])),
            ('physical_lambda',g.mul(gamma,coordinates['loglambda'])))
        multiplier=mixed.IntervalTaylor.constant(self.ctx,coefficient,0)
        return PhysicalSourceTerm(label,row,operator,multiplier,Fraction(radial_power),gamma,scale,
            row.coefficients*multiplier)

    @source_precision
    def map_source_view(self,view,Z,log_tau='-10',theta='7/10'):
        coordinate,z,zc=self.source_view(view,Z);g=self.graph;c=self.ctx
        coords=self.coordinates(view['geometry'],z,log_tau,theta)
        angle=mixed.pulse.radius.exact_coordinate(theta);tc=c.mpf(angle.numerator)/angle.denominator
        cosine,sine=c.cos(tc),c.sin(tc);L=1-self.delta*zc*zc
        if mixed.pulse.radius.post.selected.inlet.endpoints(L)[0]<=0:raise ValueError('Positive physical L required')
        log2=g.unary('log',g.constant(2));spatial={};times={}
        for component in COMPONENTS:
            spatial[component]={}
            for i,j,b in INDICES:
                terms=[];degree=i+j
                for (label,a,q),angular in cartesian_templates()[component,i,j,b].items():
                    beta,beta_box=self.beta(label)
                    gamma=g.add(beta,g.constant(-degree),g.mul(g.constant(b),g.sub(self.delta_function,g.one)))
                    gamma_box=beta_box-degree+b*(self.delta-1)
                    if mixed.pulse.radius.post.selected.inlet.endpoints(gamma_box)[1]>=0:
                        raise ArithmeticError('Original physical exponent must be negative')
                    angular_exact=exact_expression(g,angular,{CS:coords['cosine'],SN:coords['sine']})
                    two_exact=g.unary('exp',g.mul(g.constant(Fraction(a-q,2)),log2))
                    angular_box=angular_polynomial(c,angular,cosine,sine)*c.sqrt(2)**(a-q)
                    for index,expression in physical_operators()[a,b].items():
                        operator=g.mul(angular_exact,two_exact,exact_expression(g,expression,
                            {ZSYM:coords['Z'],DSYM:self.delta_function,BSYM:beta}))
                        coefficient=angular_box*interval_expression(c,expression,zc,self.delta,beta_box)
                        terms.append(self.term(view,label,index,operator,coefficient,Fraction(-degree,2),gamma,coords))
                spatial[component]['x%d_y%d_z%d'%(i,j,b)]=PhysicalSourceRow(component,(i,j,b),tuple(terms))
            terms=[]
            # At fixed physical x the moving angular basis is fixed in time.
            for (label,a,q),angular in cartesian_templates()[component,0,0,0].items():
                beta,beta_box=self.beta(label);gamma=g.add(beta,g.constant(-2))
                angular_exact=exact_expression(g,angular,{CS:coords['cosine'],SN:coords['sine']})
                angular_box=angular_polynomial(c,angular,cosine,sine)
                factors=(((0,0),g.quotient(g.neg(beta),g.mul(g.constant(2),coords['L']),'positive 2L'),-beta_box/(2*L)),
                    ((0,1),g.quotient(g.mul(g.sub(g.one,self.delta_function),coords['Z']),
                        g.mul(g.constant(2),coords['L']),'positive 2L'),(1-self.delta)*zc/(2*L)),
                    ((1,0),g.quotient(g.one,coords['L'],'positive L'),1/L))
                for index,exact,box in factors:
                    terms.append(self.term(view,label,index,g.mul(angular_exact,exact),angular_box*box,Fraction(0),gamma,coords))
            times[component]=PhysicalSourceRow(component,('t',),tuple(terms))
        self.call_trace.append(dict(chart=view['chart'],exact_native_coordinate=dict(numerator=coordinate.numerator,denominator=coordinate.denominator),
            spatial_rows=140,fixed_physical_x_time_rows=4,full_current_raw_source_consumed=True,
            ordinary_rows_not_native_Jacobian_rows=True,no_extra_Ur_sqrt2=True))
        return dict(chart=view['chart'],exact_native_coordinate=dict(numerator=coordinate.numerator,denominator=coordinate.denominator),
            exact_Z=dict(numerator=z.numerator,denominator=z.denominator),
            coordinates=coords,coordinate_domain='r>0, -1<Z<1, tau=exp(log_tau)>0, t=1-tau',
            Cartesian_spatial_rows=spatial,fixed_x_time_rows=times,
            physical_velocity_pressure={key:spatial[key]['x0_y0_z0'] for key in COMPONENTS},
            independent_P0_retained_in_absolute_pressure=True,
            original_N_definition=self.original_N_definition,
            scale_scope='exact positive source factors times signed directed derivative enclosures; no finite radius or point coefficient chosen',
            **dict.fromkeys(GATES,self.acceptance_loaded),**dict.fromkeys(OPEN,False))

    @source_precision
    def evaluate(self,chart,Z,coordinate,log_tau='-10',theta='7/10'):
        view=self.transport.evaluate(chart,Z,coordinate)
        return self.map_source_view(view,Z,log_tau,theta)


def report(view):
    return {**{key:value for key,value in view.items() if key not in
        ('coordinates','Cartesian_spatial_rows','fixed_x_time_rows','physical_velocity_pressure')},
        'exact_physical_coordinate_function_nodes':{key:value.node for key,value in view['coordinates'].items()},
        'Cartesian_spatial_rows':{key:{index:row.report() for index,row in rows.items()}
            for key,rows in view['Cartesian_spatial_rows'].items()},
        'fixed_x_time_rows':{key:row.report() for key,row in view['fixed_x_time_rows'].items()},
        'physical_velocity_pressure_rows':'x0_y0_z0 of ux,uy,uz,p; same typed rows, no numeric midpoint'}


@source_precision
def run(before=None,observed_source_views=None):
    began=time.monotonic();owner=CurrentOriginalRpPhysicalSourceMap(before,require_checked=False)
    sources=observed_source_views
    if sources is None:sources={chart:owner.transport.evaluate(chart,'.521',coordinate)
        for _,chart,coordinate in mixed.pulse.VIEWS}
    if set(sources)!=set(mixed.CHARTS):raise ValueError('One live typed observation of every current chart required')
    views={chart:owner.map_source_view(view,'.521','-10','7/10') for chart,view in sources.items()}
    result=dict(source_family=owner.family_record,candidate_current_physical_map_constructed=True,
        actual_source_graph=owner.assert_graph(),original_linear_operator_definitions=owner.definitions,
        accepted_independent_original_physical_fixture=owner.canonical_fixture,
        original_source_and_divergence_identities=owner.canonical_source_identities,
        actual_current_fifteen_chart_physical_views={key:report(view) for key,view in views.items()},
        actual_source_call_trace=owner.call_trace,original_N_definition=owner.original_N_definition,
        exact_physical_expression_graph=owner.graph.nodes,
        source_scope='Current fifteen Rp-to-heat charts; 140 Cartesian spatial rows and 4 fixed-x time rows per chart. Signed scaled enclosures, interior off-axis chart parameters. Arbitrary physical point inverse, axis, global norms, stress and temporal recursion remain open',
        **dict.fromkeys(GATES+OPEN,False),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began)
    data=json.dumps(mixed.pulse.raw.packed(result),separators=(',',':'))+'\n'
    (HERE/NAME).write_bytes(gzip.compress(data.encode(),mtime=0))
    print('CURRENT_ORIGINAL_RP_PHYSICAL_SOURCE_MAP fifteen live charts mapped',flush=True)
    return owner,views


if __name__=='__main__':run()
