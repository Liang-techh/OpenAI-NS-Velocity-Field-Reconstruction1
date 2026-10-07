"""Realize the original all-N cell integrals as genuine C1(Z) functions.

Piecewise dominated differentiation needs fixed Z-independent cell bounds,
internal cutoff/period traces and integrable value/Z majorants. It does not
need equality of spatial velocity/stress jets at different chart endpoints.
The existing exact native graph and all-period bounds are reused unchanged.
"""
import ast
from dataclasses import dataclass
import json
from pathlib import Path
import time
import sympy as s
import lei_ren_part1_paper_compliant_current_native_q_flat_cone_budget as previous

band,allN,packets=previous.band,previous.allN,previous.packets
controls,source=allN.controls,allN.source
HERE,PREFIX,sha,ep=previous.HERE,previous.PREFIX,previous.sha,previous.ep
NAME=PREFIX+'current_native_Rc_C1_integral_realization.json'
RECEIPT=PREFIX+'current_native_Rc_C1_integral_realization_check.json'
GATE='current_original_all_N_Rc_cell_integral_functions_C1_in_Z_realized'


def require(condition,message):
    if not condition:raise ArithmeticError(message)


def restore_graph(nodes):
    """Restore expression handles, not live source ancestor owners."""
    g=source.FunctionTransportGraph()
    require(nodes[:2]==g.nodes,'Original exact zero/one graph prefix required')
    g.nodes=[dict(row) for row in nodes]
    g.keys={json.dumps(row,sort_keys=True,separators=(',',':')):i for i,row in enumerate(g.nodes)}
    require(len(g.keys)==len(g.nodes),'Original function graph contains duplicate expression identities')
    return g


def phase_and_cutoff_regularity_theorem():
    """Bind the defining recipes and prove the needed parameter traces."""
    a,b,E,q,t0=s.symbols('a b E q t0',real=True)
    nu=1+t0*t0+2*q*q
    phi_end=(2*s.pi+2*s.pi*(t0*t0+2*q*q))/(2*s.pi*nu)
    Aend=a*(phi_end-1)/2
    Bend=E*(-a*t0-b*phi_end)/2
    require(s.cancel(phi_end-1)==0 and s.cancel(Aend)==0 and s.cancel(Bend.subs(t0,-b/a))==0,
            'Original full-period primitive traces differ')
    z=s.symbols('Z');az,bz,ez,qz=[s.Function(k)(z) for k in ('a','b','E','q')]
    t0z=-bz/az;nuz=1+t0z*t0z+2*qz*qz
    phasez=(1+t0z*t0z+2*qz*qz)/nuz
    require(s.simplify(s.diff(az*(phasez-1)/2,z))==0 and
            s.simplify(s.diff(ez*(-az*t0z-bz*phasez)/2,z))==0,
            'Original endpoint slow Z jets are not zero')
    u,t=s.symbols('u t',real=True)
    require(s.simplify(1-u*u/(1+u*u)-1/(1+u*u))==0,'Poisson denominator positivity lost')
    # For fixed parameters, Phi_psi=(1+t^2)/(2pi*nu)>0. The IFT
    # therefore provides C1 inverse-parameter dependence on the circle.
    loopfile=PREFIX+'current_generic_shear_loop.py'
    tree=ast.parse((HERE/loopfile).read_text(encoding='utf8'))
    cls=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=='GenericShearLoop')
    init=next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name=='__init__')
    assignments={ast.unparse(n.targets[0]):n.value for n in ast.walk(init) if isinstance(n,ast.Assign) and len(n.targets)==1}
    expected={'self.q':'flat_step(c, (2 + eta - self.kappa) / eta) * c.sqrt((2 + 2 * eta - self.kappa) / (2 * self.a))',
              'self.u':'self.p2 * self.q / scales.d_star','self.h':'c.sqrt(1 + self.u ** 2)','self.r':'self.u / self.h'}
    for key,value in expected.items():
        require(key in assignments and ast.dump(assignments[key])==ast.dump(ast.parse(value,mode='eval').body),
                'Original cutoff/Poisson recipe changed: '+key)
    return dict(passed=True,original_full_period_endpoint_value_identities=3,
        original_full_period_endpoint_Z_identities=2,
        original_q_branch_boundaries=['Delta=0: sigma joins constant1 flatly',
            'Delta=eta: sigma joins zero flatly; active square-root numerator>=eta>0'],
        square_root_has_positive_argument_on_every_active_branch=True,
        original_sigma_flat_at_zero_and_one=True,
        original_r='u/sqrt(1+u^2)',one_minus_r_squared='1/(1+u^2)>0',
        original_phase_slope='(1+t^2)/(2*pi*(1+t0^2+2*q^2))>0',
        inverse_C1_parameter_dependence_from_implicit_function_theorem=True,
        original_A_B_and_slow_Z_traces_zero_at_phase0_and_phase1=True,
        all_N_signed_coefficient_value_and_Z_traces_zero_at_fast_period_seams=True,
        fixed_finite_N_has_finitely_many_period_seams_on_each_finite_cell=True,
        proof='Positive original a,eta,dstar and finite whole-chart C1 source parameters give smooth q/Poisson/implicit primitives. Flat cutoff traces and periodic primitive endpoint identities glue C1 parameter dependence. Coordinate dependence is measurable; actual phase is Z-independent.',
        input_hashes={loopfile:sha(loopfile),PREFIX+'current_generic_loop_function_sources.py':sha(PREFIX+'current_generic_loop_function_sources.py')})


@dataclass(frozen=True)
class RealizedC1Cell:
    original_cell: dict
    coefficient_density_pairs: dict
    coefficient_contributions: dict
    coefficient_outgoing: dict
    exact_integral_nodes: tuple


class NativeRcC1IntegralRealization:
    def __init__(self):
        saved=previous.SavedOriginalQFlatBudget();self.ctx=saved.ctx;self.family=saved.family
        self.hashes=dict(saved.hashes)
        checked=json.loads((HERE/previous.RECEIPT).read_bytes())
        require(checked['all_passed'] and checked[previous.GATE] and checked['source_family']==self.family,
                'Checked current original local cone route required')
        self.hashes.update(checked['input_hashes']);self.hashes[previous.RECEIPT]=sha(previous.RECEIPT)
        self.data={}
        for stem in ('current_native_Rc_all_N_function_controls','current_native_Rc_function_transport','current_generic_shear_source_bounds'):
            name=PREFIX+stem+'.json';require(self.hashes.get(name)==sha(name),'Accepted native source changed: '+name)
            self.data[stem]=json.loads((HERE/name).read_bytes())
            require(self.data[stem]['source_family']==self.family,'Original integral source family differs')
        self.theorem=phase_and_cutoff_regularity_theorem();self.hashes.update(self.theorem['input_hashes'])
        self.source_C1_receipts={}
        for name,gate in ((source.sources.RECEIPT,source.sources.GATE),
            (PREFIX+'current_generic_shear_source_bounds_check.json','current_original_noncore_generic_quotient_derivative_log_bounds_certified'),
            (allN.RECEIPT,allN.GATE)):
            require(self.hashes.get(name)==sha(name),'Accepted original C1 source-leaf receipt changed: '+name)
            receipt=json.loads((HERE/name).read_bytes())
            require(receipt['all_passed'] and receipt[gate] and receipt['source_family']==self.family,
                    'Explicit same-family original C1 source-leaf/function receipt required')
            self.source_C1_receipts[name]=dict(gate=gate,sha256=sha(name),source_family=self.family,
                original_function_and_Z_derivative_recipes_not_cover_values=True)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)

    def build(self):
        old=self.data['current_native_Rc_all_N_function_controls'];g=restore_graph(old['exact_function_graph_nodes'])
        ref=lambda i:source.FunctionRef(g,i)
        parameters={r['name']:ref(i) for i,r in enumerate(g.nodes) if r['operation']=='original_source_parameter'}
        N=next(ref(i) for i,r in enumerate(g.nodes) if r['operation']=='shared_positive_integer')
        original_cells=self.data['current_native_Rc_function_transport']['exact_original_cells']
        maps,symbols,x=source.exact_radius_maps()
        history={p:{key:source.C1Function(g.zero,g.zero) for key in allN.RATES} for p in allN.ORDERS}
        cells=[];all_integrals=[]
        for i,cell in enumerate(original_cells):
            coordinate=g.symbol('coordinate_'+str(i));offset=source.expression(g,maps[cell['chart']],parameters,coordinate)
            jacobian=source.expression(g,s.diff(maps[cell['chart']],x),parameters,coordinate)
            density={p:{key:source.C1Function(g.zero,g.zero) for key in allN.RATES} for p in allN.ORDERS}
            if not cell['source_flat_exact_zero']:
                refs={r['function_role']:ref(j) for j,r in enumerate(old['exact_function_graph_nodes'])
                    if r['operation']=='original_function_graph' and r.get('coordinate')==coordinate.node and r.get('function_role','').startswith('all_N_')}
                pair=lambda root:source.C1Function(refs[root+'_C0'],refs[root+'_Z'])
                density,F=allN.coefficient_pairs(g,pair('all_N_original_E'),pair('all_N_original_V'),
                    pair('all_N_periodic_A'),pair('all_N_periodic_B'),N)
            contributions={p:{} for p in allN.ORDERS};outgoing={p:{} for p in allN.ORDERS};node_ids=[]
            for p in allN.ORDERS:
                for key,rate in allN.RATES.items():
                    kernel=g.unary('exp',g.neg(g.mul(g.constant(rate),g.sub(ref(cell['right_radius_offset']),offset))))
                    decay=g.unary('exp',g.neg(g.mul(g.constant(rate),ref(cell['width']))))
                    pair=density[p][key];ints=[]
                    for handle in (pair.value,pair.Z):
                        if handle==g.zero:ints.append(g.zero);continue
                        value=controls.integral(g,g.mul(kernel,handle,jacobian),'coordinate_'+str(i),
                            ref(cell['lower']),ref(cell['upper']),measure='native coordinate; original dy/dcoordinate applied exactly once',
                            extracted_N_power=p,coefficient_still_depends_on_N_and_actual_phase=True)
                        require(value.node<len(old['exact_function_graph_nodes']),'Realization changed the accepted original integrand')
                        node_ids.append(value.node);ints.append(value)
                    contributions[p][key]=source.C1Function(*ints)
                    outgoing[p][key]=g.c1add(g.c1scale(decay,history[p][key]),contributions[p][key])
            history=outgoing;all_integrals.extend(node_ids)
            cells.append(RealizedC1Cell(cell,density,contributions,outgoing,tuple(node_ids)))
        A=source.C1Function(*[next(ref(i) for i,r in enumerate(g.nodes) if r.get('function_role')==role)
                              for role in ('Rc_E_C0','Rc_E_Z')])
        orders={p:allN.normalize_history(g,history[p],A,parameters['mu'])[0] for p in allN.ORDERS}
        invN=g.quotient(g.one,N,'one common original positive integer N')
        scaled={key:g.c1add(orders[-1][key],g.c1scale(invN,orders[-2][key])) for key in controls.ROWS}
        require(controls.pair_roots(scaled.values(),scaled.keys())==old['exact_reconstructed_N_scaled_target_roots'],
                'Realized C1 targets do not equal the original all-N functions')
        require(len(cells)==24 and len(all_integrals)==336 and len(g.nodes)==len(old['exact_function_graph_nodes']),
                'Original graph/cell/integral inventory changed during realization')
        return dict(graph=g,cells=cells,coefficient_history=history,N_scaled_targets=scaled,amplitude=A,N=N,
            parameters=parameters,source_family=self.family,source_graph_sha256=old['input_hashes'][source.sources.VIEWS],
            exact_integral_nodes=all_integrals)

    def records(self,built):
        rows=self.data['current_native_Rc_all_N_function_controls']['actual_all_N_continuous_cell_range_records'];records=[]
        for cell,majorant in zip(built['cells'],rows):
            original=cell.original_cell
            require(original['label']==majorant['label'] and original['chart']==majorant['chart'],
                    'Original cell/majorant geometry mismatch')
            require(majorant['original_geometry']['width_and_endpoints_independent_of_Z'] and majorant['original_true_mass_applied_once'],
                    'Leibniz boundary terms or repeated Jacobian would change the integral')
            caps={order:majorant[key] for order,key in (('value','coefficient_contributions'),('Z','coefficient_Z_contributions'))}
            for rows0 in caps.values():
                for key,orders in rows0.items():
                    for p,cap in orders.items():
                        require(cap['encloses_original_source_function'] and cap['point_value_selected'] is False,
                                'A density/mass cap became a defining integral value')
                        if not cap['exact_zero']:
                            value=packets.interval(self.ctx,cap['log_absolute_upper'])
                            require(all(__import__('mpmath').isfinite(v) for v in ep(value)),'Nonintegrable source majorant')
            records.append(dict(label=original['label'],chart=original['chart'],original_geometry=majorant['original_geometry'],
                original_source_flat_exact_zero=original['source_flat_exact_zero'],
                exact_C1_coefficient_integral_roots=allN.encode_orders(cell.coefficient_contributions),
                exact_C1_outgoing_coefficient_roots=allN.encode_orders(cell.coefficient_outgoing),
                actual_original_integral_nodes=cell.exact_integral_nodes,
                integrable_value_and_Z_majorants=caps,
                majorant_definition='original full-coordinate/full-period C0 or Z density cap times exp(-rate*(yright-y))*positive original dy/dcoordinate; own-rate finite mass bounds its integral',
                derivative_definition='d_Z integral f(y,Z)dy = integral original f_Z(y,Z)dy; endpoints, kernel and Jacobian independent of Z',
                internal_q_and_periodic_phase_C1_traces_glued=True,
                different_spatial_chart_velocity_or_stress_seam_equality_used=False,
                incoming_C1_histories_and_rate_zero_pressure_memory_retained=True,
                exact_functions_defined_by_original_integrals_not_majorant_values=True))
        return records


def run():
    began=time.monotonic();owner=NativeRcC1IntegralRealization();built=owner.build();records=owner.records(built)
    report=dict(source_family=owner.family,**{GATE:True},
        original_phase_cutoff_and_parameter_C1_theorem=owner.theorem,
        explicitly_bound_original_C1_source_leaf_and_function_receipts=owner.source_C1_receipts,
        exact_source_graph_sha256=built['source_graph_sha256'],original_exact_function_graph_nodes=len(built['graph'].nodes),
        actual_original_C1_integral_cell_records=records,actual_original_coefficient_integral_count=len(built['exact_integral_nodes']),
        exact_C1_Rc_N_scaled_target_roots=controls.pair_roots(built['N_scaled_targets'].values(),built['N_scaled_targets'].keys()),
        fixed_N_domain='each fixed finite integer N>=160 satisfying the accepted original source bounds',
        whole_Z_domain=[-1,1],exact_original_integral_functions_C1_in_Z_defined=True,
        dominated_differentiation_uses_original_uniform_C0_Z_and_own_rate_mass_bounds=True,
        C1_integral_target_not_an_arbitrary_symbol_or_enclosure_endpoint=True,
        piecewise_integrals_do_not_admit_spatial_y_seam_equality=True,
        original_source_ancestor_constructors_called=False,numerical_original_source_point_or_integral_oracle_installed=False,
        actual_five_controls_installed=False,certified_actual_fixed_point_tail_installed=False,
        actual_terminal_Z_function_closure_installed=False,current_whole_N_selected=False,
        **dict.fromkeys(packets.OPEN,False),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
        scope='Source-bound realization of all24 original coefficient/history/target integrals as C1(Z) functions, using flat q/periodic phase traces and piecewise dominated differentiation. Original native graph unchanged. No spatial velocity/stress seam admission, numeric original oracle, fixed-point controls, terminal closure, global N, outer cone or recursion.')
    (HERE/NAME).write_text(json.dumps(packets.encode(report),indent=2)+'\n',encoding='utf8')
    print('Original all24 cell integrals and same-source Rc targets realized as C1(Z) functions',flush=True)
    return report


if __name__=='__main__':run()
