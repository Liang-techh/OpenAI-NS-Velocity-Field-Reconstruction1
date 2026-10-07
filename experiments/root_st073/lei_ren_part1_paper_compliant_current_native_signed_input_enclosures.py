"""Live signed Section 11 root enclosures without enormous exponentials.

The returned intervals enclose analytic source functions. No interval midpoint
or upper bound defines a function value. A positive-function theorem can cut
an otherwise sign-indefinite denominator enclosure. Conservative radius-log
covers and the four source logs remain separate factors; exponentials below a cutoff retain a
directed nonzero tail enclosure. This is a numerical root backend, before the
conditioned phase inverse, changed-history integrals and common N selection.
"""
import ast
import copy
import dataclasses
import importlib
import json
import math
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_native_generic_source_packets as native
import lei_ren_part1_paper_compliant_current_generic_shear_signed_jets as signed
from lei_ren_part1_paper_compliant_global_physical_assembly import CompliantGlobalPhysicalAssembly
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import sigma_jets

packets=native.packets; HERE,PREFIX,sha=native.HERE,native.PREFIX,native.sha
NAME=PREFIX+'current_native_signed_input_enclosures.json'
RECEIPT=PREFIX+'current_native_signed_input_enclosures_check.json'
GATE='current_live_signed_generic_root_interval_backend_executed'
ZERO=(0,0,0,0,0)
ep=packets.recovery.endpoints


class FormalScale:
    """A collected linear form; matching factors cancel before evaluation."""
    def __init__(self,bases,powers=ZERO,offset=0):
        self.bases=bases; self.ctx=bases[0].ctx
        self.powers=tuple(powers); self.offset=self.ctx.mpf(offset)
        if len(self.powers)!=5 or any(type(x) not in (int,float) or 2*x!=int(2*x) for x in self.powers):
            raise ValueError('Five exact integer/half-integer factor powers required')
    def pair(self,other):
        if other.bases is not self.bases:raise ValueError('One fixed source/radius basis required')
    def __add__(self,other):
        self.pair(other)
        return FormalScale(self.bases,tuple(a+b for a,b in zip(self.powers,other.powers)),self.offset+other.offset)
    def __neg__(self):return FormalScale(self.bases,tuple(-a for a in self.powers),-self.offset)
    def __sub__(self,other):return self+-other
    def evaluate(self):
        # Zero powers are removed first, including exact radius cancellation.
        return sum((b*p for b,p in zip(self.bases,self.powers) if p),self.offset)
    def record(self):return dict(source_exponents=list(self.powers[:4]),radius_power=self.powers[4],additional_log_interval=self.offset)


class ScaledEnclosure:
    """Signed coefficient interval times exp(a formal log), never a point."""
    def __init__(self,scale,coefficient,ledger):
        self.scale=scale; self.ctx=scale.ctx; self.coefficient=self.ctx.mpf(coefficient); self.ledger=ledger
        if any(not mp.isfinite(x) for x in ep(self.coefficient)+ep(scale.offset)):
            raise ValueError('Finite coefficient and formal offset required')
        magnitude=max(abs(x) for x in ep(self.coefficient))
        if magnitude:
            logmag=self.ctx.ln(self.ctx.mpf(magnitude))
            threshold=self.ctx.dps*mp.log(10)
            if ep(logmag)[0]>threshold or ep(logmag)[1]<-threshold:
                # This coordinate normalization retains the entire signed box.
                # It prevents a tiny source mu from being swamped by a tail cap
                # expressed in unrelated unit coordinates.
                self.coefficient=self.coefficient/self.ctx.mpf(magnitude)
                self.scale=FormalScale(scale.bases,scale.powers,scale.offset+logmag)
                ledger['coefficient_log_coordinate_rescalings']=ledger.get('coefficient_log_coordinate_rescalings',0)+1
    @property
    def zero(self):return ep(self.coefficient)==(0,0)
    def scalar(self,value):return ScaledEnclosure(FormalScale(self.scale.bases),value,self.ledger)
    def coerce(self,other):
        if not isinstance(other,ScaledEnclosure):other=self.scalar(other)
        self.scale.pair(other.scale)
        if other.ledger is not self.ledger:raise ValueError('Same directed arithmetic ledger required')
        return other
    def __neg__(self):return ScaledEnclosure(self.scale,-self.coefficient,self.ledger)
    def bounded_exp(self,log):
        lo,hi=ep(log); c=self.ctx
        cutoff=-2*c.dps*c.ln(10)
        if hi<ep(cutoff)[0]:
            # This is a tail enclosure, not a substitution of zero.
            self.ledger['directed_small_exponential_tails']+=1
            return c.mpf([0,ep(c.exp(cutoff))[1]])
        if hi>2*c.dps*mp.log(10):
            raise ArithmeticError('Positive relative exponential too large; choose a dominating source scale')
        if lo<ep(cutoff)[0]:return c.mpf([0,ep(c.exp(c.mpf(hi)))[1]])
        return c.exp(log)
    def __add__(self,other):
        other=self.coerce(other)
        if self.zero:return other
        if other.zero:return self
        # Compare only to choose coordinates; no endpoint becomes a field value.
        left,right=self,other
        if ep(left.scale.evaluate())[1]<ep(right.scale.evaluate())[1]:left,right=right,left
        relative=(right.scale-left.scale).evaluate()
        if ep(relative)[1]>2*self.ctx.dps*mp.log(10):
            # A theorem-intersected scale can be an independent wide log box.
            # Its endpoints cannot share the defining source correlations.
            # A common upper log is merely an arithmetic coordinate, with
            # directed tails; it is never selected as a source value.
            reference=FormalScale(self.scale.bases,offset=self.ctx.mpf(max(
                ep(left.scale.evaluate())[1],ep(right.scale.evaluate())[1])))
            a=left.coefficient*left.bounded_exp(left.scale.evaluate()-reference.offset)
            b=right.coefficient*right.bounded_exp(right.scale.evaluate()-reference.offset)
            self.ledger['directed_independent_log_rescalings']+=1
            return ScaledEnclosure(reference,a+b,self.ledger)
        factor=left.bounded_exp(relative)
        return ScaledEnclosure(left.scale,left.coefficient+right.coefficient*factor,self.ledger)
    __radd__=__add__
    def __sub__(self,other):return self+-self.coerce(other)
    def __rsub__(self,other):return self.coerce(other)+-self
    def __mul__(self,other):
        other=self.coerce(other)
        return ScaledEnclosure(self.scale+other.scale,self.coefficient*other.coefficient,self.ledger)
    __rmul__=__mul__
    def positive_divide(self,other,log_positive_lower):
        other=self.coerce(other); lo,hi=ep(other.coefficient); c=self.ctx
        if hi<=0:raise ArithmeticError('Source cover contradicts positive denominator theorem')
        if lo>0:
            return ScaledEnclosure(self.scale-other.scale,self.coefficient/other.coefficient,self.ledger)
        lower=ep(c.mpf(log_positive_lower))[0]
        upper=ep(other.scale.evaluate()+c.ln(c.mpf(hi)))[1]
        if lower>upper:raise ArithmeticError('Positive source lower bound exceeds enclosing upper bound')
        # The theorem applies to the original function, not to each coefficient.
        # Keep the positive logarithmic interval without exp(huge) or 1/zero.
        self.ledger['positive_function_denominator_intersections']+=1
        denominator_log=FormalScale(self.scale.bases,offset=c.mpf([lower,upper]))
        return ScaledEnclosure(self.scale-denominator_log,self.coefficient,self.ledger)
    def positive_intersection(self,lower):
        lo,hi=ep(self.coefficient)
        if lo>0:return self
        if hi<=0:raise ArithmeticError('Original positive root has no positive enclosing values')
        loglo=ep(self.ctx.mpf(lower))[0]
        loghi=ep(self.scale.evaluate()+self.ctx.ln(self.ctx.mpf(hi)))[1]
        if loglo>loghi:raise ArithmeticError('Positive root theorem and source cover disagree')
        self.ledger['positive_function_root_intersections']+=1
        return ScaledEnclosure(FormalScale(self.scale.bases,offset=self.ctx.mpf([loglo,loghi])),1,self.ledger)
    def finite_interval(self,max_log=1000):
        """Only bounded exponentials may be materialized for downstream kernels."""
        if self.zero:return self.ctx.mpf(0)
        lo,hi=ep(self.scale.evaluate())
        if lo < -max_log or hi > max_log:
            raise ArithmeticError('Root remains logarithmically factored; scalar materialization forbidden')
        return self.coefficient*self.ctx.exp(self.scale.evaluate())
    def record(self):
        lo,hi=ep(self.coefficient); magnitude=max(abs(lo),abs(hi))
        return dict(coefficient_interval=self.coefficient,formal_positive_scale=self.scale.record(),
            exact_zero=self.zero,sign='zero' if self.zero else 'positive' if lo>0 else 'negative' if hi<0 else 'undetermined',
            log_absolute_upper=None if not magnitude else self.scale.evaluate()+self.ctx.ln(self.ctx.mpf(magnitude)),
            point_value_selected=False,encloses_original_source_function=True)


def compile_factored_amplitude_adapter(stem,name,algebra,logu):
    """Original row formula with exp(logu) stored as a formal positive unit.

Only arithmetic order is changed at named scalar-unit/P0 sites. The
accepted arbitrary frozen-unit normalization theorem already treats this unit
symbolically. Its source amplitude and its square are never capped here.
"""
    module=importlib.import_module(PREFIX+stem)
    tree=ast.parse(Path(module.__file__).read_text(encoding='utf8'))
    fn=copy.deepcopy(next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name==name))
    changes=[]
    class Orient(ast.NodeTransformer):
        def visit_BinOp(self,node):
            self.generic_visit(node)
            if isinstance(node.op,ast.Mult) and isinstance(node.right,ast.Name) and node.right.id in ('invP2','u0','u20'):
                changes.append(ast.unparse(node)); node.left,node.right=node.right,node.left
            elif isinstance(node.op,ast.Add) and isinstance(node.left,ast.Name) and node.left.id=='p0':
                changes.append(ast.unparse(node)); node.left,node.right=node.right,node.left
            return node
    Orient().visit(fn)
    expected=['A[j] * invP2',"p0 + histories['p'][0]"] if stem=='current_reshape_stress_operator' else ['amp * u0','amp2 * u20','A[j] * invP2',"p0 + histories['p'][0]"]
    if changes!=expected:raise ValueError('Original amplitude adapter multiplication sites changed')
    calls=[]
    def unit(ctx,value):
        if ctx is not algebra.ctx:raise ValueError('Original source amplitude context required')
        matches=[p for p in (1,2) if ctx.mpf(value)._mpi_==(logu*p)._mpi_]
        if len(matches)!=1:raise ValueError('Only the same actual frozen amplitude and square may be deferred')
        power=matches[0];calls.append(power)
        return algebra.shift(1,(0,0,0,power/2))
    env=dict(vars(module))
    env['positive_exp']=unit
    env['scaled_positive_source']=lambda ctx,value,row,proofs:unit(ctx,value)*row
    code=compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),'<original source amplitude deferred>','exec')
    exec(code,env)
    return env[name],dict(source_module=stem,source_function=name,
        commutative_scalar_factor_sites_reoriented=changes,
        positive_exponential_calls=calls,original_row_formula_unchanged=True,
        frozen_positive_source_log=logu,positive_source_cap_not_called=True)


class NativeSignedInputEnclosures:
    def __init__(self,backend):
        if type(backend) is not native.NativeGenericSourcePackets:raise ValueError('Existing genuine17-chart native backend required')
        self.native=backend; self.service=backend.service; self.ctx=backend.ctx; self.family=backend.family
        if backend.seed.radius.__func__ is not CompliantGlobalPhysicalAssembly.radius:
            raise ValueError('Unchanged original radius binding required')
        receipt=json.loads((HERE/signed.current.RECEIPT).read_bytes())
        signed.validate_source_receipt(receipt,self.family)
        self.service.bind_hashes(receipt['input_hashes'])
        self.service.bind_hashes({signed.current.RECEIPT:sha(signed.current.RECEIPT)})
        self.inventory=json.loads((HERE/signed.current.NAME).read_bytes())['current_actual_loop_jet_log_bounds_by_chart']
        self.scales=json.loads((HERE/(PREFIX+'current_generic_shear_uniform_inputs.json')).read_bytes())['current_actual_logarithmic_loop_scales']
        self.service.bind_hashes({name:sha(name) for name in (signed.current.NAME,
            PREFIX+'current_generic_shear_uniform_inputs.json',Path(__file__).name,
            PREFIX+'global_physical_assembly.py',PREFIX+'current_generic_shear_signed_jets.py')})
        self.queries=[]
    def decode(self,value):
        if isinstance(value,dict):
            if 'lower' in value and 'upper' in value:return packets.interval(self.ctx,value)
            return {k:self.decode(v) for k,v in value.items()}
        if isinstance(value,list):return [self.decode(v) for v in value]
        return value
    @native.inlet.source_precision
    def packet(self,chart,Z,coordinate):
        base=self.native.query(chart,Z,coordinate)
        source=self.native.original_raw(chart,Z,coordinate)
        logR,radius=self.native.seed.radius(chart,source['coordinate'],source['source_packet'],source['native_provider'])
        provenance=dict(base.provenance,coordinate_contract=packets.CHARTS[chart][3] if chart in packets.CHARTS
            else 'original O3 selector; ordinary logR rows already converted',
            original_radius_source=radius,logR_cover=logR,
            radius_cover_is_original_conservative_enclosure_not_a_point=True,
            formal_microscopic_offset_retained_as_source_metadata_only=True,
            arithmetic_radius_basis='original conservative logR cover; formal radius tree not evaluated',
            formal_radius_correlations_used_in_arithmetic=False)
        adapter=None
        if chart in ('switch_power','reshape','inner_reference','axial_restore','restore_buffer'):
            parent=source['source_packet']['actual_inherited_axial5_packet']
            logu=packets.jet(self.ctx,parent['log_Utheta_over_Pstar_axial5_coefficients'])[0]
            algebra=packets.FactoredAlgebra(self.ctx,(self.ctx.mpf(0),2*self.native.seed.logP,self.ctx.mpf(0),2*logu),[])
            stem='current_reshape_stress_operator' if chart in ('switch_power','reshape') else 'current_restore_stress_operator'
            name='raw_reshape_rows' if stem=='current_reshape_stress_operator' else 'raw_restore_rows'
            fn,adapter=compile_factored_amplitude_adapter(stem,name,algebra,logu)
            raw=fn(self.ctx,source['source_packet'],algebra.shift(1,(0,-1,0,0)))
            def rows(values):return tuple(packets.decode_row(algebra,packets.factored_rows_record(v)) for v in values)
            nv={key:rows(v) for key,v in raw['velocity'].items()}; nm={key:rows(v) for key,v in raw['histories'].items()}
            vel={key:tuple(algebra.shift(v,packets.INVERSE_S) if key in ('axial','radial') else v for v in values) for key,values in nv.items()}
            hist={key:tuple(algebra.shift(v,packets.INVERSE_S) if key in ('m','k') else v for v in values) for key,values in nm.items()}
            base=packets.CurrentSourcePacket(chart,base.source_family,algebra,base.Z,vel,hist,rows(raw['absolute_pressure']),
                algebra.lift(source['P0']),nv,nm,provenance)
            provenance.update(amplitude_adapter=adapter,source_factor_resolution_performed=False,
                original_capped_amplitude_covers_replaced_by_exact_formal_unit=True)
        else:
            provenance.update(plain_native_rows_are_prebounded_covers=source['native_provider_group'] in ('patch','pre','O3'),
                original_hidden_factors_not_invented=True)
        return dataclasses.replace(base,provenance=provenance)
    def leaf(self,polynomial,bases,ledger):
        result=ScaledEnclosure(FormalScale(bases),0,ledger)
        for radius,row in polynomial.terms.items():
            for mode,jet in row.terms.items():
                value=ScaledEnclosure(FormalScale(bases,mode+(radius,)),jet[0],ledger)
                result=result+value
        return result
    @native.inlet.source_precision
    def query(self,chart,Z,coordinate):
        packet=self.packet(chart,Z,coordinate)
        positive=self.decode(self.inventory[chart]['actual_positive_denominator_theorem'])
        definition=signed.from_packet(packet,self.service.data['delta'],positive)
        leaves=signed.source_leaves(packet,self.service.data['delta'])
        bases=tuple(packet.algebra.logs)+(packet.provenance['logR_cover'],)
        ledger=dict(directed_small_exponential_tails=0,positive_function_denominator_intersections=0,
            positive_function_root_intersections=0,directed_independent_log_rescalings=0)
        leaf_values={name:self.leaf(value,bases,ledger) for name,value in leaves.items()}
        graph=definition['jet_expression_dag'];values=[]
        certificate={key:self.decode(v['log_positive_lower']) for key,v in graph['positive_denominator_certificates'].items()}
        scalar=lambda v:ScaledEnclosure(FormalScale(bases),v,ledger)
        for node in graph['nodes']:
            op=node['operation']
            if op=='constant':value=scalar(node['value'])
            elif op=='source_derivative':value=leaf_values[node['name']]
            elif op=='sum':value=sum((values[k] for k in node['arguments']),scalar(0))
            elif op=='negative':value=-values[node['argument']]
            elif op=='product':
                value=scalar(1)
                for k in node['arguments']:value=value*values[k]
            elif op=='positive_function_quotient':
                value=values[node['numerator']].positive_divide(values[node['denominator']],certificate[node['positive_certificate']])
            else:raise ValueError('Known signed source arithmetic instruction required')
            values.append(value)
        roots={name:{order:values[k] for order,k in row.items()} for name,row in graph['roots'].items()}
        roots['E']['y0_Z0']=roots['E']['y0_Z0'].positive_intersection(positive['log_E_positive_lower'])
        roots['a']['y0_Z0']=roots['a']['y0_Z0'].positive_intersection(positive['log_actual_a_positive_lower'])
        exact_excess=None
        if chart.startswith('O3_'):
            mu=packets.interval(self.ctx,positive['actual_positive_mu']);v=self.ctx.mpf(coordinate)
            sigma=sigma_jets(self.ctx,v) if chart=='O3_slope_mu' else None
            for j,k in signed.ORDERS:
                delta=0 if k else (2*mu*sigma[j]*math.factorial(j) if sigma is not None else 2*mu if j==0 else 0)
                roots['kappa_minus2']['y%d_Z%d'%(j,k)]=scalar(delta)
            exact_excess=dict(original_source_definition=positive['exact_a_minus2_source'],
                original_positive_mu=mu,rounded_two_plus_mu_not_subtracted=True)
        eta=ScaledEnclosure(FormalScale(bases,offset=packets.interval(self.ctx,self.scales['selected_positive_eta_log'])),1,ledger)
        branch_difference=roots['kappa_minus2']['y0_Z0']-eta
        lo,hi=ep(branch_difference.coefficient)
        branch='flat' if lo>=0 else 'active' if hi<0 else 'requires_source_box_refinement'
        result=dict(chart=chart,source_family=self.family,source_provenance=packet.provenance,
            signed_ordinary_root_enclosures={name:{order:v.record() for order,v in row.items()} for name,row in roots.items()},
            root_orders=[list(k) for k in signed.ORDERS],fixed_source_log_bases=packet.algebra.logs,
            original_radius_log_cover=packet.provenance['logR_cover'],
            numerical_arithmetic_ledger=ledger,original_correlated_O3_excess=exact_excess,
            branch_against_actual_eta=branch,branch_difference_enclosure=branch_difference.record(),
            source_caps_or_midpoints_used_as_defining_values=False,
            original_radius_or_source_amplitude_exponential_materialized=False,
            positive_denominator_theorems_used_as_function_constraints=True,
            inverse_phase_or_signed_density_backend_installed=False,
            current_whole_N_selected=False,actual_integral_functions_installed=False,
            **dict.fromkeys(packets.OPEN,False))
        self.queries.append(dict(chart=chart,Z_box=packet.provenance['Z_box'],coordinate_box=packet.provenance['coordinate_box'],
            branch=branch,root_count=sum(len(v) for v in roots.values()),ledger=dict(ledger)))
        return dict(packet=packet,roots=roots,record=result,definition=definition,leaf_values=leaf_values)


def run(bridge=None):
    began=time.monotonic()
    if bridge is None:bridge,_=native.inlet.native_bridge_owner()
    with native.inlet.CheckedSourceRuntime():
        owner=NativeSignedInputEnclosures(native.NativeGenericSourcePackets(bridge))
        points=dict(bridge_first='.1337',bridge_second='1.831',bridge_macro='.537',switch_first='.537',
            switch_second='1.337',switch_power='.537',reshape='.537',inner_reference='.537',axial_restore='.537',
            restore_buffer='-6.337',actual_patch='1.337',Rh_reference='-2.337',O2_slope='.537',O2_axial='.1337',
            O2_buffer='5.337',O3_slope_mu='.537',O3_power=owner.ctx.mpf('.537')/owner.native.pre.params.Tw)
        records={}
        for chart in native.DOMAINS:
            records[chart]=owner.query(chart,(-1,1),points[chart])['record']
            print('Live signed source roots:',chart,records[chart]['branch_against_actual_eta'],flush=True)
    result=dict(source_family=owner.family,**{GATE:True},live_native_chart_count=len(records),
        live_signed_ordinary_root_enclosure_count=sum(v['root_count'] for v in owner.queries),
        current_native_signed_source_root_records=records,query_inventory=owner.queries,
        numerical_caps_used_as_field_values=False,ancestor_source_graph_rebuilt=False,
        selected_positive_eta_is_same_original_log_source=True,
        inverse_phase_and_changed_density_functions_installed=False,
        actual_defect_integral_functions_installed=False,actual_repair_control_functions_installed=False,
        **dict.fromkeys(packets.OPEN,False),execution_seconds=time.monotonic()-began,
        scope='Live signed y2/Z1 source-root interval enclosures from one native17-chart seed, collected log factors and positive-function denominator constraints; exact O3 excess retained separately. Not phase inversion, common N, changed signed-density integration, repair, global cone or coefficient recursion.',
        input_hashes=dict(owner.service.hashes))
    (HERE/NAME).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()
