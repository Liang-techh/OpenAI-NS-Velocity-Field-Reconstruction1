"""Issued Rh-reference integral callbacks using genuine whole-cell sources.

The accepted defining coefficient functions, ordinary-Z rows and canonical
factor algebra supply all five own-rate integrals. The current precise phase
is applied to the same original radius. This is a partial fixed-nonzero-Z
callback, not an all-chart source oracle or a solved control system.
"""
from dataclasses import dataclass
from fractions import Fraction
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_point_source_leaves as leaves
import lei_ren_part1_paper_compliant_current_original_reference_whole_cell_integrals as whole
import lei_ren_part1_paper_compliant_current_native_Rc_all_N_function_controls_check as identities

HERE,PREFIX,sha=leaves.HERE,leaves.PREFIX,leaves.sha
precise=leaves.precise;transport=leaves.transport;ep=leaves.ep
NAME=PREFIX+'current_original_reference_integral_callback.json'
RECEIPT=PREFIX+'current_original_reference_integral_callback_check.json'
GATE='actual_original_Rh_reference_graph_bound_C0_Z_integral_callbacks_installed'


def symbolic_node(provider,index,source=None):
    """Collect exact radius differences before any numerical evaluation."""
    g=provider.built['graph'];parameters={ref.node:precise.sy.Symbol(name,real=True)
        for name,ref in provider.built['parameters'].items()};memo={}
    density=sy.Symbol('issued_original_signed_density',real=True)
    def visit(i):
        if i in memo:return memo[i]
        row=g.nodes[i];op=row['operation']
        if source is not None and row is source:value=density
        elif i in parameters:value=parameters[i]
        elif op=='exact_rational':value=sy.Rational(row['numerator'],row['denominator'])
        elif op=='bound_variable':value=sy.Symbol(row['name'],real=True)
        elif op=='sum':value=sy.Add(*(visit(q) for q in row['arguments']))
        elif op=='product':value=sy.Mul(*(visit(q) for q in row['arguments']))
        elif op=='negative':value=-visit(row['argument'])
        elif op=='positive_quotient':value=visit(row['numerator'])/visit(row['denominator'])
        elif op=='analytic_unary' and row['name'] in ('exp','log'):
            value=(sy.exp if row['name']=='exp' else sy.log)(visit(row['argument']))
        else:raise ValueError('Unbound original integral expression: '+op)
        memo[i]=value;return value
    return visit(index),density


def exact_density_binding(provider):
    """Prove the reused coefficient program equals these issued density roots."""
    exact=whole.points.exact;g=transport.FunctionTransportGraph()
    N=g.node('shared_positive_integer',name='N',lower=160)
    pairs={};outside={};inside={};bindings=[]
    view=provider.views['Rh_reference']
    E,EZ,V,VZ,A,AZ,B,BZ=sy.symbols('E E_Z V V_Z mathcal_A mathcal_A_Z mathcal_B mathcal_B_Z',real=True)
    source=view['original_signed_input_graph']['jet_expression_dag']['roots']
    original_ids={'E':(source['E']['y0_Z0'],source['E']['y0_Z1']),
        'A':(view['roots']['A'],view['roots']['A_Z_slow']),
        'B':(view['roots']['B_over_Pstar'],view['roots']['B_Z_slow'])}
    axial=[]
    for label in ('V_y0_Z0','V_y0_Z1'):
        matches=[i for i,row in enumerate(view['function_graph_nodes']) if row.get('operation')=='source_derivative' and row.get('name')==label]
        if len(matches)!=1:raise ValueError('One defining original axial source binding required')
        axial.append(matches[0])
    original_ids['V']=tuple(axial)
    for name,symbols in (('E',(E,EZ)),('V',(V,VZ)),('A',(A,AZ)),('B',(B,BZ))):
        handles=[g.symbol(str(symbol)) for symbol in symbols]
        pairs[name]=transport.C1Function(*handles)
        for handle,symbol,original in zip(handles,symbols,original_ids[name]):
            outside[handle.node]=symbol;inside[original]=symbol
            bindings.append(dict(input=name,symbol=str(symbol),actual_original_source_node=original))
    coefficients,F=exact.coefficient_pairs(g,pairs['E'],pairs['V'],pairs['A'],pairs['B'],N)
    old=identities.symbolic_walk(view['function_graph_nodes'],inside)
    new=identities.symbolic_walk(g.nodes,outside);Ns=sy.Symbol('N',positive=True,integer=True);roots=[]
    for key in transport.RATES:
        for jet,attribute in (('C0','value'),('Z','Z')):
            issued=provider.rows['Rh_reference','density_'+key+'_'+jet];provider.require_row(issued)
            reconstructed=sum(new(getattr(coefficients[order][key],attribute).node)*Ns**order for order in exact.ORDERS)
            if sy.simplify(old(issued['source_node'])-reconstructed)!=0:
                raise ValueError('Whole-cell coefficient program differs from issued signed density')
            roots.append(dict(role=issued['function_role'],source_node=issued['source_node'],exact_identity_passed=True))
    return dict(passed=True,actual_issued_C0_Z_full_signed_density_identities=roots,
        original_E_V_A_B_C0_Z_atom_bindings=bindings,
        reused_coefficient_pairs_source=Path(exact.__file__).name,reused_whole_cell_source=Path(whole.__file__).name,
        accepted_whole_cell_enclosure_receipt=whole.RECEIPT,
        retained_exact_N_dependent_exprel_and_exp_functions=True,
        exprel_identity_continuously_extended_at_A_zero=True,
        original_B_over_Pstar_normalization_no_additional_division=True)


class ExactPhaseReferenceCells(whole.OriginalReferenceWholeCells):
    def __init__(self,provider,*,Z):
        self.provider=provider
        super().__init__(Z=Z)
        if self.family!=provider.family:raise ValueError('Same genuine reference source family required')
        original=self.dispatcher.reference.owner.inputs.frame
        current=provider.reference.owner.inputs.frame
        if original.selected_logCstar_mpf_tuple!=current.selected_logCstar_mpf_tuple or any(
                original.definitions[k]!=current.definitions[k] for k in ('logPstar','log_delta','logRref','L')):
            raise ValueError('Same defining reference parameter/factor recipes required')
        old_phase=self.dispatcher.reference.owner.radius.original
        binder=provider.phase.binder
        if tuple(old_phase.sc_tuple)!=tuple(ep(binder.sc)[0]._mpf_) or ep(binder.sc)[0]!=ep(binder.sc)[1]:
            raise ValueError('Same exact selected microscopic origin endpoint required')
        if old_phase.domain['formal_radii'].get('Ra')!='4*epsilon_core' or old_phase.domain['formal_radii'].get('r_minus')!='Ra*exp(hb*s_c/2)':
            raise ValueError('Unchanged original r_minus recipe required')
        P,C,B,sc,x=[provider.phase.symbols[k] for k in ('logP','logC','hbB','sc')]+[provider.phase.x]
        original_logR=sy.log(110)+10*(C+P)+x
        original_log_r_minus=sy.log(4)-4*P-1000+B*sc/2
        offset=provider.phase.maps['Rh_reference']
        if sy.simplify(offset-(original_logR-original_log_r_minus))!=0 or sy.diff(offset,x)!=1:
            raise ValueError('Original reference phase argument/coordinate Jacobian differs')
        if 'Rh_reference' not in binder.identity['exact_source_radius_minus_same_left_radius_identities']:
            raise ValueError('Accepted original radius identity required')
        self.phase_argument_identity=dict(passed=True,old_and_current_exact_argument_equal_not_endpoint_overlap=True,
            original_logR=sy.srepr(original_logR),original_log_r_minus=sy.srepr(original_log_r_minus),
            same_selected_sc_mpf_tuple=list(old_phase.sc_tuple),
            original_hb_width_source_attachment=old_phase.collar['exact_current_exit_source_attachment'],
            same_origin_recipe_source_hash=sha(leaves.base.radius.phase.COLLAR),
            accepted_current_native_radius_identity=True,reference_Jacobian_exact_one=True,
            phase_cell_increment_exact='N*(right-left)',original_microscopic_origin_not_zeroed=True)

    def phase_boxes(self,left,right,N):
        row=self.provider.rows['Rh_reference','density_m_C0'];c=self.ctx
        phase=self.provider.phase.phase_for_source(row,self.provider.built,coordinate=str(left),N=N)
        boxes=[c.mpf(ep(value)) for value in phase['actual_phase_directed_boxes']]
        shift=c.mpf((0,ep(c.mpf(N)*self.atlas.rational(right-left))[1]))
        projection=precise.periodic_add(c,boxes,shift)
        return projection['boxes'],dict(phase,exact_source_cell_phase_increment='N*(right-left)',
            reference_native_radius_Jacobian_exact_one=True,
            actual_cell_phase_is_original_point_origin_plus_exact_coordinate_width=True)


@dataclass(frozen=True)
class OriginalReferenceIntegralFrame:
    Z: object
    N: int
    count: int
    bits: int
    owner: object
    values: dict
    coefficients: dict
    cells: tuple
    record: dict


class OriginalReferenceIntegralCallback:
    mode='original_reference_integral_partial'
    def __init__(self,provider=None):
        self.provider=leaves.OriginalPointSourceLeaves() if provider is None else provider
        if type(self.provider) is not leaves.OriginalPointSourceLeaves:
            raise TypeError('Genuine issued point-source provider required')
        self.family=self.source_family=self.provider.family
        self.source_graph_sha256=self.provider.source_graph_sha256
        self.hashes=dict(self.provider.hashes);self.owners={};self.cache={};self.issued={}
        for module in (leaves,whole):
            accepted=json.loads((HERE/module.RECEIPT).read_bytes())
            if not accepted.get('all_passed') or not accepted.get(module.GATE) or accepted['source_family']!=self.family:
                raise ValueError('Accepted same-source point/cell prerequisite required')
            self.bind_hashes({**accepted['input_hashes'],module.RECEIPT:sha(module.RECEIPT)})
        self.density_identity=exact_density_binding(self.provider)
        self.bind_hashes({Path(identities.__file__).name:sha(Path(identities.__file__).name)})
        saved=json.loads((HERE/transport.NAME).read_bytes())
        cells=[cell for cell in saved['exact_original_cells'] if cell['chart']=='Rh_reference']
        if len(cells)!=1:raise ValueError('One exact original reference route cell required')
        self.original_cell=cells[0];g=self.provider.built['graph'];self.roles={};self.rows={};self.contracts=[]
        for key,pair in self.original_cell['contributions'].items():
            for jet,label in (('C0','value'),('Z','Z')):
                index=pair[label];row=g.nodes[index];self.roles[index]=(key,jet)
                self.rows[key,jet]=row
                self.contracts.append(self.require_integral(row))
        self.bind_hashes({Path(__file__).name:sha(Path(__file__).name)})

    def bind_hashes(self,closure):
        for name,digest in closure.items():
            if sha(name)!=digest:raise ValueError('Original integral dependency changed: '+name)
            if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Original integral closures disagree: '+name)
            self.hashes[name]=digest

    def require_integral(self,row):
        provider=self.provider;g=provider.built['graph']
        if leaves.digest_rows(g.nodes)!=provider.graph_digest:
            raise ValueError('Original integral graph must remain unchanged')
        matches=[i for i,item in enumerate(g.nodes) if row is item]
        if len(matches)!=1 or matches[0] not in self.roles:
            raise ValueError('Issued Rh_reference definite-integral row required')
        index=matches[0];key,jet=self.roles[index]
        if row.get('operation')!='definite_integral' or not row.get('exact_function_integral'):
            raise ValueError('Exact defining integral required')
        if row['lower']!=self.original_cell['lower'] or row['upper']!=self.original_cell['upper']:
            raise ValueError('Unchanged original graph endpoints required')
        lower,_=symbolic_node(provider,row['lower']);upper,_=symbolic_node(provider,row['upper'])
        if lower!=-5 or upper!=0:raise ValueError('Exact reference window[-5,0] required')
        source=provider.rows['Rh_reference','density_'+key+'_'+jet];provider.require_row(source)
        coordinate=g.nodes[source['coordinate']]
        if coordinate!={'operation':'bound_variable','name':row['variable']}:
            raise ValueError('Same original bound coordinate required')
        x=sy.Symbol(row['variable'],real=True);rate=transport.RATES[key]
        expression,density=symbolic_node(provider,row['integrand'],source)
        if sy.simplify(expression-sy.exp(sy.Rational(rate.numerator,rate.denominator)*x)*density)!=0:
            raise ValueError('Original own-rate kernel/source/Jacobian binding differs')
        if row['measure']!='native coordinate; original dy/dcoordinate applied exactly once':
            raise ValueError('Original coordinate measure required')
        return dict(integral_node=index,source_node=source['source_node'],source_role=source['function_role'],
            chart='Rh_reference',ordinary_order=jet,key=key,exact_endpoints=['-5','0'],
            exact_reduced_integrand='exp(rate*coordinate)*issued_original_signed_density',
            rate=str(rate),original_coordinate_Jacobian=1,
            huge_common_radius_terms_canceled_symbolically=True)

    def frame(self,*,Z,N,count=4,bits=32):
        z=leaves.base.point.pressure.exact_Z(Z)
        if z==0:raise ValueError('Nonzero reference Z required by the accepted whole-cell signed route')
        if type(N) is not int or N<160 or N.bit_length()>4096:
            raise ValueError('Explicit original graph integer N>=160, at most4096 bits required')
        if type(count) is not int or not 1<=count<=256:raise ValueError('Bounded exact partition count in[1,256] required')
        if type(bits) is not int or not 4<=bits<=256:raise ValueError('Inverse bits in[4,256] required')
        self.require_integral(self.rows['m','C0']);cache_key=(z,N,count,bits)
        if cache_key in self.cache:return self.cache[cache_key]
        if z not in self.owners:
            owner=ExactPhaseReferenceCells(self.provider,Z=str(z));self.bind_hashes(owner.hashes);self.owners[z]=owner
        owner=self.owners[z];c=owner.ctx;a=owner.atlas
        total={order:{key:{jet:a.scalar(0) for jet in ('C0','Z')} for key in transport.RATES}
            for order in whole.points.exact.ORDERS};cells=[];began=time.monotonic()
        with mp.workdps(c.dps+40):
            for i in range(count):
                left=sy.Rational(-5)+sy.Rational(5*i,count);right=sy.Rational(-5)+sy.Rational(5*(i+1),count)
                query=owner.cell(left,right,N=N,bits=bits);masses={}
                for key,rate in transport.RATES.items():
                    rr=c.mpf(rate.numerator)/rate.denominator
                    mass=a.rational(right-left) if not rate else (c.exp(rr*a.rational(right))-c.exp(rr*a.rational(left)))/rr
                    if ep(mass)[0]<=0:raise ArithmeticError('Positive exact own-rate mass lost')
                    masses[key]=mass
                    for order in total:
                        for jet in ('C0','Z'):
                            total[order][key][jet]=a.add(total[order][key][jet],query['coefficients'][order][key][jet]*mass)
                cells.append(dict(source=query['record'],exact_positive_own_rate_masses=masses,
                    original_log_radius_Jacobian_applied_once=1,nonlinear_exact_N_coefficients_enclosed_before_phase_union=True))
            full={key:{jet:a.add(total[-1][key][jet]*(c.mpf(1)/N),total[-2][key][jet]*(c.mpf(1)/N**2))
                for jet in ('C0','Z')} for key in transport.RATES}
            ordinary={};unmaterialized=[]
            for key in full:
                ordinary[key]={}
                for jet,value in full[key].items():
                    try:ordinary[key][jet]=leaves.base.conditioned.bounded_value(value)
                    except ArithmeticError:unmaterialized.append([key,jet])
        record=dict(source_family=self.family,source_graph_sha256=self.source_graph_sha256,
            original_Z_exact=str(z),explicit_candidate_N=N,exact_original_window=['-5','0'],
            exact_source_cells=count,inverse_bits=bits,canonical_atlas=a.record(),
            exact_original_integral_bindings=self.contracts,
            local_coefficient_program_equals_issued_original_density=self.density_identity,
            old_and_current_original_phase_argument_identity=owner.phase_argument_identity,
            exact_N_coefficient_own_rate_integrals={str(order):{key:{jet:v.record() for jet,v in pair.items()}
                for key,pair in rows.items()} for order,rows in total.items()},
            actual_native_C0_Z_integral_values={key:{jet:v.record() for jet,v in pair.items()} for key,pair in full.items()},
            ordinary_directed_integral_intervals=ordinary,unmaterializable_values_retained_factored=unmaterialized,
            genuine_whole_cell_source_phase_and_error_records=cells,
            defining_source_not_sampled_quadrature=True,point_caps_or_midpoints_not_source_values=True,
            original_order_minus1_and_minus2_N_dependent_functions_retained=True,
            source_pressure_late_errors_and_positive_auxiliary_scales_retained=True,
            incoming_corrections_and_original_P0_not_reset=True,
            integral_enclosure_error_contract='Defining closed source cell intervals, directed original pressure errors, exact radius-phase covers, certified inverse and original ordinary-Z ranges, then positive kernel masses; no unproved point quadrature remainder.',
            complete_target_averaging_not_claimed_from_local_integral=True,
            full_scalar_control_evaluator_compatibility_installed=False,
            full_17_chart_source_or_24_cell_integral_oracle_installed=False,
            actual_five_controls_installed=False,current_whole_N_selected=False,
            execution_seconds=time.monotonic()-began)
        frame=OriginalReferenceIntegralFrame(z,N,count,bits,owner,full,total,tuple(cells),record)
        self.cache[cache_key]=frame;self.issued[id(frame)]=frame;return frame

    def dispatch(self,row,frame):
        contract=self.require_integral(row)
        if type(frame) is not OriginalReferenceIntegralFrame or self.issued.get(id(frame)) is not frame:
            raise ValueError('Issued genuine reference integral frame required')
        return frame.values[contract['key']][contract['ordinary_order']]

    def integral_factored(self,row,*,Z,N,count=4,bits=32):
        self.require_integral(row)
        return self.dispatch(row,self.frame(Z=Z,N=N,count=count,bits=bits))

    def integrate(self,row,*,Z,N,count=4,bits=32):
        """Standalone issued graph callback; full evaluator still unavailable."""
        self.require_integral(row);frame=self.frame(Z=Z,N=N,count=count,bits=bits)
        with mp.workdps(frame.owner.ctx.dps+40):result=leaves.base.conditioned.bounded_value(self.dispatch(row,frame))
        if not hasattr(result,'_mpi_') or hasattr(result,'scale'):raise TypeError('Ordinary directed integral interval required')
        return result

    def parameter(self,name):raise NotImplementedError('Complete original parameter evaluator remains open')
    def source(self,*args,**kwargs):raise NotImplementedError('All-chart source callback remains open')


@precise.phase.native.inlet.source_precision
def run(*,return_live=False):
    began=time.monotonic();owner=OriginalReferenceIntegralCallback();frames=[]
    for count,N in ((4,160),(16,160),(4,257)):
        frame=owner.frame(Z='.37',N=N,count=count);frames.append(frame)
        for row in owner.rows.values():owner.integrate(row,Z='.37',N=N,count=count)
        print('Genuine issued original reference integrals:',count,'cells, N',N,flush=True)
    result=dict(**{GATE:True},source_family=owner.family,source_graph_sha256=owner.source_graph_sha256,
        mode=owner.mode,exact_original_integral_bindings=owner.contracts,
        actual_original_reference_integral_frames=[frame.record for frame in frames],
        actual_original_graph_integral_callbacks=10*len(frames),
        original_source_integral_graph_unchanged=leaves.digest_rows(owner.provider.built['graph'].nodes)==owner.provider.graph_digest,
        full_17_chart_source_or_24_cell_integral_oracle_installed=False,
        full_scalar_control_evaluator_compatibility_installed=False,
        actual_five_controls_installed=False,actual_terminal_Z_function_closure_installed=False,
        current_whole_N_selected=False,**dict.fromkeys(precise.phase.packets.OPEN,False),
        input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
        scope='Actual issued Rh_reference[-5,0] five C0/Z definite-integral callbacks at fixed nonzero Z using genuine whole-cell source intervals, current precise phase and retained exact N-dependent coefficients. Partial provider only; full source/evaluator, all24 integrals, controls/global N, recursion and corrected NS remain open.')
    (HERE/NAME).write_bytes(json.dumps(precise.encode(result),indent=2).encode()+b'\n')
    return (result,owner,frames) if return_live else result


if __name__=='__main__':run()
