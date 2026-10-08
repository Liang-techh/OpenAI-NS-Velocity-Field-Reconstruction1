"""Execute the accepted original O2 signed five-rate graph on source points.

Native C0/phase-held Z primitives bind graph nodes. Stable expm1 retains
tiny factored changes. These are normalized own-rate densities per logR,
not an integral or a globally admitted common frequency/field.
"""
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_O2_conditioned_slow_Z as slow
import lei_ren_part1_paper_compliant_current_generic_shear_moment_recovery as recovery

base=slow.base;HERE,PREFIX,sha=slow.HERE,slow.PREFIX,slow.sha;ep=slow.ep
NAME=PREFIX+'current_original_O2_signed_densities.json'
RECEIPT=PREFIX+'current_original_O2_signed_densities_check.json'
GATE='original_O2_true_radius_signed_five_own_rate_densities_and_Z_enclosures_connected'


class BoundDensityGraph:
    """Interpret the accepted graph, stopping at actual source/kernel bindings."""
    def __init__(self,graph,query,primitives,jets,N):
        if type(N) is not int or N<1:raise ValueError('One explicit positive integer candidate N required')
        self.graph=graph;self.nodes=graph['function_graph_nodes'];self.kernel=query['kernel']
        self.c=self.kernel.c;self.scalar=self.kernel.scalar;self.bindings={};self.cache={};self.expcache={}
        self.operations={};self.ledger=query['ledger'];groot=graph['roots'];self.has_Z=jets is not None
        source=graph['original_signed_input_graph']['jet_expression_dag']['roots']
        for order in ((0,1) if self.has_Z else (0,)):
            self.bindings[source['E']['y0_Z%d'%order]]=query['roots']['E'][(0,order)]
            label='V_y0_Z%d'%order
            ids=[i for i,node in enumerate(self.nodes) if node.get('operation')=='source_derivative' and node.get('name')==label]
            if len(ids)!=1:raise ValueError('Unique original axial source leaf required: '+label)
            self.bindings[ids[0]]=query['roots']['V'][(0,order)]
        self.bindings.update({groot['A']:primitives['A'],groot['B_over_Pstar']:primitives['B_over_Pstar'],groot['N']:self.scalar(N)})
        if self.has_Z:self.bindings.update({groot['A_Z_slow']:jets['A_Z_slow'],groot['B_Z_slow']:jets['B_Z_slow']})
        if self.nodes[groot['N']]!={'operation':'shared_positive_integer_parameter','name':'common_N',
            'domain':'N>=1','same_for_all_charts':True}:raise ValueError('Original common-N node contract differs')
        contract=graph['common_phase_binding']
        if not contract['total_Z_phase_derivative_exact_zero'] or not contract['slow_derivatives_hold_phi_fixed']:
            raise ValueError('Original phase-held Z convention required')
        self.transport=graph['own_defect_transport_contract']
        if self.transport['rates']!=recovery.RATES:raise ValueError('Original own-rate normalization differs')

    def modulation_exp(self,index):
        if index in self.expcache:return self.expcache[index]
        x=self.evaluate(index);c=self.c
        # Original O2 a=.8+1.2*sigma, 0<=sigma<=1, and
        # A=a*(phi-psi/(2*pi))/2 imply |A/N|<=1 for every N>=1.
        # Intersection is a function range proof, never a selected value.
        finite=base.conditioned.clipped(c,base.conditioned.bounded_value(x),-1,1)
        M=64;mean=c.mpf(1);power=c.mpf(1)
        for k in range(1,M+1):
            power*=finite;mean+=power/c.factorial(k+1)
        tail=ep(c.exp(1)/c.factorial(M+2))[1]
        mean+=c.mpf((-tail,tail))
        mean=base.conditioned.clipped(c,mean,c.exp(-1),c.exp(1))
        changed=x*mean
        result=(self.scalar(1)+changed,changed)
        self.ledger['bounded_modulation_exp_taylor_tail_operations']=self.ledger.get('bounded_modulation_exp_taylor_tail_operations',0)+1
        self.expcache[index]=result;return result

    def evaluate(self,index):
        if index in self.bindings:return self.bindings[index]
        if index in self.cache:return self.cache[index]
        node=self.nodes[index];op=node['operation'];self.operations[op]=self.operations.get(op,0)+1
        if op=='constant':value=self.scalar(node['value'])
        elif op=='mathematical_pi':value=self.scalar(self.c.pi)
        elif op=='sum':value=sum((self.evaluate(i) for i in node['arguments']),self.scalar(0))
        elif op=='product':
            value=self.scalar(1)
            for i in node['arguments']:value*=self.evaluate(i)
        elif op=='negative':value=-self.evaluate(node['argument'])
        elif op=='positive_function_quotient':
            if node['positive_certificate'] not in ('common_N_positive_integer','exact_positive_integer2','exact_positive_2pi'):
                raise ValueError('Unexpected unbound original density denominator')
            numerator=self.evaluate(node['numerator']);denominator=self.evaluate(node['denominator'])
            lower=ep(denominator.coefficient)[0]
            if lower<=0:raise ArithmeticError('Exact positive graph denominator lost')
            value=numerator.positive_divide(denominator,denominator.scale.evaluate()+self.c.ln(self.c.mpf(lower)))
        elif op=='analytic_unary' and node['name'] in ('exp','expm1'):
            positive,changed=self.modulation_exp(node['argument']);value=positive if node['name']=='exp' else changed
        else:raise ValueError('Density subgraph reached an unbound source/inverse operation: '+op)
        self.cache[index]=value;return value

    def values(self):
        graph=self.graph;roots=graph['roots']
        velocities={key:self.evaluate(roots[key]) for key in ('E_N','V_N','delta_E','delta_V')}
        return dict(velocities=velocities,densities={k:self.evaluate(i) for k,i in graph['five_signed_increment_rate_roots'].items()})

    def outputs(self):
        if not self.has_Z:raise ValueError('Actual Z primitive/source bindings required for Z outputs')
        graph=self.graph;Z=graph['five_signed_increment_rate_first_derivatives']['Z'];result=self.values()
        velocities=result['velocities']
        velocities.update(E_N_Z=self.evaluate(graph['changed_velocity_first_ordinary_derivatives']['theta']['Z']),
            V_N_Z=self.evaluate(graph['changed_velocity_first_ordinary_derivatives']['axial']['Z']),
            delta_E_Z=self.evaluate(Z['h']),delta_V_Z=self.evaluate(Z['m']))
        result['density_Z']={k:self.evaluate(i) for k,i in Z.items()};return result

    def record(self):
        return dict(accepted_function_graph_executed_not_density_formula_rewritten=True,
            exact_node_bindings={str(i):value.record() for i,value in self.bindings.items()},
            evaluated_operation_counts=self.operations,source_E_V_and_Z_bound=True,
            actual_C0_bound=True,actual_phase_held_Z_bound=self.has_Z,common_N_not_replaced_by_pointwise_frequency=True,
            stable_expm1='expm1(x)=x*integral_0^1 exp(t*x)dt; mean exponential uses 64 terms and directed exp(1)/66! tail',
            tiny_original_factored_x_not_zeroed=True,original_O2_modulation_argument_absolute_bound=1,
            bound_proof='a=.8+1.2*sigma in[.8,2]; A=a*(phi-psi/(2*pi))/2; N>=1',
            own_rate_units_with_common_S_Pstar=recovery.UNITS,own_rates=recovery.RATES,
            integration_coordinate='y=log(R/Rref); no additional R/Jacobian in normalized density',
            actual_inlet_histories_and_P0_required_for_own_transport=True)


class OriginalO2SignedDensities:
    mode='original_O2_true_radius_signed_density_and_Z_enclosures'
    def __init__(self):
        admitted=json.loads((HERE/slow.RECEIPT).read_bytes())
        if not admitted.get('all_passed') or not admitted.get(slow.GATE):raise ValueError('Accepted actual original O2 slow-Z receipt required')
        self.source=slow.OriginalO2ConditionedSlowZ();self.family=self.source.family;self.owner=self.source.owner
        if admitted['source_family']!=self.family:raise ValueError('O2 density and derivative families differ')
        self.hashes=dict(self.source.hashes)
        for name,digest in {**admitted['input_hashes'],slow.RECEIPT:sha(slow.RECEIPT)}.items():
            if sha(name)!=digest:raise ValueError('Original density source dependency changed: '+name)
            if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Density/derivative dependencies disagree')
            self.hashes[name]=digest
        self.graph=self.owner.scales.graph
        if self.graph['source_family']!=self.family:raise ValueError('Accepted O2 graph family differs')
        for name in (Path(__file__).name,Path(recovery.__file__).name):self.hashes[name]=sha(name)

    def evaluate(self,*,y,Z,N,bits=80):
        if type(bits) is not int or not 4<=bits<=256:raise ValueError('Explicit inverse bits in[4,256] required')
        phase=self.owner.radius.evaluate(y=y,N=N);c=self.owner.ctx
        key=(base.point.source.exact_coordinate(y),base.point.pressure.exact_Z(Z))
        reused=key in self.owner.inputs.point_cache
        with mp.workdps(c.dps+40):
            query=self.owner.query(y=y,Z=Z);kernel=query['kernel'];records=[]
            for box in phase['true_original_phase_directed_boxes']:
                C0=kernel.evaluate(c.mpf([box['lower'],box['upper']]),bits=bits)
                if C0['status']!='enclosed':raise ArithmeticError('Actual density source inverse requires refinement')
                C0.update(free_phase_parameter_not_spatial_phase=False,original_common_N_and_radius_phase_bound=True)
                selected=C0['selected_inverse'];coordinate=selected['coordinate_interval'];chart=selected['chart']
                primitives=kernel.primitives(coordinate,chart)
                jets,derivative=slow.slow_values(kernel,query['roots'],coordinate,chart)
                interpreter=BoundDensityGraph(self.graph,query,primitives,jets,N);got=interpreter.outputs()
                records.append(dict(C0=C0,phase_held_Z_contract=derivative,
                    modulated_normalized_velocity_and_Z_enclosures={k:v.record() for k,v in got['velocities'].items()},
                    five_signed_own_rate_density_enclosures={k:v.record() for k,v in got['densities'].items()},
                    five_signed_own_rate_density_Z_enclosures={k:v.record() for k,v in got['density_Z'].items()},
                    actual_graph_execution_contract=interpreter.record()))
        return dict(source_family=self.family,mode=self.mode,original_y_exact=str(key[0]),original_Z_exact=str(key[1]),
            explicit_candidate_N=N,actual_original_radius_phase=phase,source_factor_basis=query['basis_contract'],
            actual_original_O2_factored_inputs=base.point.record_query(query['point']),
            actual_signed_density_point_queries=records,numerical_arithmetic_ledger=query['ledger'],
            full_original_five_signed_density_values_and_Z_installed_on_this_O2_point=True,
            source_caps_or_midpoints_used_as_field_values=False,cached_actual_defining_point_coefficients_and_errors_reused=reused,
            actual_changed_five_moment_integral_evaluated=False,numerical_original_source_point_or_integral_oracle_installed=False,
            actual_five_controls_installed=False,current_whole_N_selected=False)


def run():
    began=time.monotonic();owner=OriginalO2SignedDensities()
    precursor=json.loads((HERE/slow.NAME).read_bytes());records=[]
    for saved in precursor['actual_original_O2_slow_Z_point_queries']:
        y,Z,N=(saved[key] for key in ('original_y_exact','original_Z_exact','explicit_candidate_N'))
        row=owner.evaluate(y=y,Z=Z,N=N);records.append(row)
        print('True-radius original O2 signed densities:',y,Z,N,flush=True)
    report=dict(**{GATE:True},source_family=owner.family,mode=owner.mode,
        actual_original_O2_signed_density_queries=records,
        full_original_five_signed_density_values_and_Z_installed_on_actual_O2_points=True,
        accepted_exact_density_function_graph_executed=True,
        original_normalized_own_rates_and_common_velocity_units_preserved=True,
        no_defining_quadratures_or_ancestor_constructors_reexecuted=True,
        actual_changed_five_moment_integral_evaluated=False,numerical_original_source_point_or_integral_oracle_installed=False,
        actual_five_controls_installed=False,current_whole_N_selected=False,
        **dict.fromkeys(base.point.source.inertial.profiles.loop.OPEN,False),
        input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
        scope='Actual original O2 signed five own-rate density and Z point enclosures, via accepted exact graph and true radius phase. No continuous-cell integral, own transported histories, all-chart oracle, controls, common N or recursion.')
    (HERE/NAME).write_text(json.dumps(base.encoded(report),indent=2)+'\n',encoding='utf8')
    return report


if __name__=='__main__':run()
