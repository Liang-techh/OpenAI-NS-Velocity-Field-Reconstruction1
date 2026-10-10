"""Current fifteen-chart affine radius inverse, with exact cancellation.

Source-owned parameter bounds enclose expressions; they do not define them.
All possibly overlapping charts are returned. No native coordinate midpoint,
velocity packet, seam value or huge absolute radius is materialized.
"""
from fractions import Fraction
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_Rp_physical_inverse as inverse
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

physical=inverse.physical
radius_source=physical.mixed.pulse.radius
HERE,PREFIX,sha=inverse.HERE,inverse.PREFIX,inverse.sha
NAME=PREFIX+'current_original_Rp_native_chart_locator.json.gz'
RECEIPT=PREFIX+'current_original_Rp_native_chart_locator_check.json'
GATES=('current_original_Rp_native_affine_inverse_functions_installed',
       'current_original_Rp_source_bound_logRp_enclosure_installed',
       'current_original_Rp_native_chart_location_installed')
OPEN=tuple(key for key in inverse.OPEN if key not in GATES)
ends=inverse.ends
CHARTS=tuple(physical.mixed.CHARTS)
ANCHORS=dict(pulse_entrance='1',pulse_main='5',pulse_exit='21/2',pulse_gap='23/2',
    pulse_gap_end='-5',pulse_end='-2',flatten='50',outer_power='1/2',outer_angular='-2',
    steep_entry='1/2',steep_power='1/2',steep_exit='1/2',waiting='1/2',heat_collar='3/2',heat_exterior='4')


class CancelledGraphBounds:
    """Expand only exact rational arithmetic before directed evaluation.

    Integrals and analytic operations stay exact atomic source functions.
    In particular hbB*sc cancels before either unknown parameter is needed.
    """
    def __init__(self,graph,ctx,bindings,allowed_exponentials):
        self.graph=graph;self.ctx=ctx;self.bindings=bindings;self.polys={};self.memo={}
        self.allowed_exponentials=frozenset(allowed_exponentials)

    def polynomial(self,i):
        i=i.node if hasattr(i,'node') else i
        if i in self.polys:return self.polys[i]
        n=self.graph.nodes[i];op=n['operation'];out={}
        def add(p,sign=1):
            for key,value in p.items():out[key]=out.get(key,Fraction(0))+sign*value
        if op=='exact_rational':out[()]=Fraction(n['numerator'],n['denominator'])
        elif op=='sum':
            for j in n['arguments']:add(self.polynomial(j))
        elif op=='negative':add(self.polynomial(n['argument']),-1)
        elif op=='product':
            out={():Fraction(1)}
            for j in n['arguments']:
                product={}
                for a,av in out.items():
                    for b,bv in self.polynomial(j).items():
                        key=tuple(sorted(a+b));product[key]=product.get(key,Fraction(0))+av*bv
                out=product
        else:out[(i,)]=Fraction(1)
        out={key:value for key,value in out.items() if value}
        self.polys[i]=out;return out

    def at(self,i):
        i=i.node if hasattr(i,'node') else i
        if i in self.bindings:return self.ctx.mpf(self.bindings[i])
        if i in self.memo:return self.memo[i]
        c=self.ctx;n=self.graph.nodes[i];op=n['operation']
        if op in ('exact_rational','sum','product','negative'):
            value=c.mpf(0)
            for key,coefficient in self.polynomial(i).items():
                term=inverse.box(c,coefficient)
                for atom in key:term*=self.at(atom)
                value+=term
        elif op=='positive_quotient':
            denominator=self.at(n['denominator'])
            if ends(denominator)[0]<=0:raise ValueError('Positive current graph denominator required')
            value=self.at(n['numerator'])/denominator
        elif op=='analytic_unary' and n['name']=='log':
            argument=self.at(n['argument'])
            if ends(argument)[0]<=0:raise ValueError('Positive current source logarithm required')
            value=c.ln(argument)
        elif op=='analytic_unary' and n['name']=='exp':
            # Only original finite parameter exponentials are evaluated here.
            # Absolute R/lambda exponentials and root nodes require bindings.
            if i not in self.allowed_exponentials:
                raise ValueError('Only whitelisted original finite parameter exponentials are admitted')
            argument=self.at(n['argument'])
            value=c.exp(argument)
        else:raise ValueError('No admitted current source bound for graph node '+str(i)+'/'+op)
        if not all(mp.isfinite(v) for v in ends(value)):raise ValueError('Finite directed log geometry required')
        self.memo[i]=value;return value


class CurrentOriginalRpNativeChartLocator:
    @source_precision
    def __init__(self,before=None,require_checked=True):
        self.before=before if before is not None else inverse.CurrentOriginalRpPhysicalInverse()
        if type(self.before) is not inverse.CurrentOriginalRpPhysicalInverse or not self.before.acceptance_loaded:
            raise ValueError('Accepted current physical inverse required')
        self.graph=self.before.graph;self.ctx=self.before.ctx;self.radius=self.before.radius
        self.family_record=self.before.family_record;self.hashes=dict(self.before.hashes)
        self.closed=self.before.before.owner.closed
        self.acceptance_loaded=False;self.call_trace=[];self._source_logR_nodes=set()
        self.bindings=self.source_bindings();g=self.graph
        self.source_exponentials={g.unary('exp',g.constant(40)).node,self.radius.functions['mu'].node}
        self.maps={}
        for chart in CHARTS:
            actual=self.radius._map(chart,g.zero)
            lower,upper=self.domain(chart)
            self.maps[chart]=dict(origin=self.radius.logRp,pulse_term=actual['pulse_term'],
                offset_zero=actual['offset'],base=actual['logR'],
                jacobian=actual['native_to_log_radius_jacobian'],lower=lower,upper=upper)
        self.assert_graph()
        for name in (inverse.NAME,inverse.RECEIPT,Path(__file__).name):self.hashes[name]=sha(name)
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not all(receipt[k] for k in GATES) or any(receipt[k] for k in OPEN):
                raise ValueError('Current locator receipt or scope differs')
            if receipt['source_family']!=self.family_record:raise ValueError('Current locator family differs')
            physical.mixed.pulse.radius.post.selected.inlet.add_hashes(self.hashes,receipt['input_hashes'])
            self.acceptance_loaded=True

    def source_bindings(self):
        self.closed.assert_graph();r=self.closed.exact.repair;f=self.radius.functions
        if not self.closed.absolute_binding['same_exact_selected_logC_source'] or not \
                self.closed.absolute_binding['actual_absolute_Rp_matches_existing_repair_reference']:
            raise ValueError('Accepted current absolute source radius identity required')
        selected=r.heat.logC
        if selected._mpi_[0]!=selected._mpi_[1]:raise ValueError('Exact selected logC source datum required')
        leaf=self.graph.nodes[self.radius.parameters['logC'].node]
        if leaf['operation']!='current_original_source_parameter' or leaf['name']!='logC' or \
                leaf['source_family']!=self.family_record or not leaf['quantity_is_exact_original_parameter_not_range_endpoint']:
            raise ValueError('Current exact logC parameter graph required')
        self.hashes[leaf['current_source_receipt']]=sha(leaf['current_source_receipt'])
        if self.hashes[leaf['current_source_receipt']]!=leaf['current_source_receipt_sha256']:
            raise ValueError('Exact selected parameter source receipt changed')
        post=self.radius.before
        return {self.radius.parameters['logC'].node:selected,
            f['mu'].node:post.before.pulse.mu,f['Lrel'].node:post.outer.Lrel,
            f['Ts'].node:post.steep.Ts,f['waiting'].node:post.steep.wait}

    def assert_graph(self):
        same=self.source_bindings()
        result=dict(accepted_current_inverse=self.before.acceptance_loaded and all(self.before.assert_graph().values()),
            same_current_source_graph=self.graph is self.before.graph is self.radius.graph,
            accepted_current_absolute_radius_binding=self.closed.acceptance_loaded and all(self.closed.assert_graph().values()),
            exact_current_parameter_bounds_unchanged=set(same)==set(self.bindings)
                and all(same[key]._mpi_==self.bindings[key]._mpi_ for key in same),
            only_original_finite_parameter_exponentials=self.source_exponentials=={
                self.graph.unary('exp',self.graph.constant(40)).node,self.radius.functions['mu'].node},
            all_fifteen_actual_source_maps=set(self.maps)==set(self.radius.maps)==set(CHARTS))
        for chart,row in self.maps.items():
            actual=self.radius._map(chart,self.graph.zero)
            if (row['base'],row['jacobian'],row['origin'],row['pulse_term'],row['offset_zero'],row['lower'],row['upper'])!=(
                    actual['logR'],actual['native_to_log_radius_jacobian'],self.radius.logRp,
                    actual['pulse_term'],actual['offset'],*self.domain(chart)):
                result['all_fifteen_actual_source_maps']=False
        if not all(result.values()):raise ValueError('Current native locator graph differs: '+str(result))
        return result

    def domain(self,chart):
        g=self.graph;mu=self.radius.functions['mu']
        ranges={'pulse_main':('1/50',10),'pulse_exit':(10,11),'pulse_gap':(11,12),
            'pulse_end':(-4,0),'flatten':(0,100),'outer_angular':(-4,0),
            'heat_collar':(0,3),'heat_exterior':(3,None)}
        if chart=='pulse_entrance':return g.zero,g.quotient(g.constant('1/50'),mu,'accepted current positive mu')
        if chart=='pulse_gap_end':return g.neg(g.quotient(g.one,mu,'accepted current positive mu')),g.constant(-4)
        left,right=ranges.get(chart,(0,1))
        return g.constant(left),None if right is None else g.constant(right)

    def reader(self,extra=None):
        return CancelledGraphBounds(self.graph,self.ctx,{**self.bindings,**(extra or {})},self.source_exponentials)

    def source_log_radius(self,chart,coordinate):
        """Register a radius at an admitted exact rational source coordinate."""
        self.assert_graph();geometry=self.radius.geometry(chart,coordinate)
        value=radius_source.FunctionRef(self.graph,geometry['logR'])
        self._source_logR_nodes.add(value.node);return value

    def source_endpoint_log_radius(self,chart,side):
        """Register the original exact source endpoint, including +/-1/mu."""
        self.assert_graph()
        if chart not in self.maps or side not in ('left','right'):raise ValueError('Exact current chart endpoint required')
        endpoint=self.maps[chart]['lower' if side=='left' else 'upper']
        if endpoint is None:raise ValueError('Infinite tail has no finite right endpoint')
        value=self.radius._map(chart,endpoint)['logR']
        self._source_logR_nodes.add(value.node);return value

    @source_precision
    def locate_function(self,logR):
        """Trusted exact graph evaluation; no caller-provided numeric bounds."""
        if type(logR) is not radius_source.FunctionRef or logR.graph is not self.graph:
            raise TypeError('Registered current source radius function required')
        if logR.node not in self._source_logR_nodes:
            raise ValueError('Use this locator source_log_radius or source_endpoint_log_radius')
        return self._locate(logR,self.reader(),None)

    def _locate(self,logR,reader,physical_Z_resolved):
        self.assert_graph();g=self.graph;c=self.ctx
        if type(logR) is not radius_source.FunctionRef or logR.graph is not g:
            raise TypeError('Same current graph log radius required')
        if reader.graph is not g or reader.ctx is not c:raise ValueError('Same current directed graph reader required')
        absolute=reader.at(logR);rows={};possible=[];contained=[]
        for chart,route in self.maps.items():
            difference=g.sub(logR,route['base'])
            vref=g.quotient(difference,route['jacobian'],'accepted current strictly positive native radius Jacobian')
            delta=reader.at(difference);jac=reader.at(route['jacobian'])
            if ends(jac)[0]<=0:raise ValueError('Positive native radius Jacobian required')
            native=delta/jac;lo,hi=ends(native)
            left=reader.at(route['lower']);ll,lh=ends(left)
            right=None if route['upper'] is None else reader.at(route['upper'])
            rl,rh=(mp.inf,mp.inf) if right is None else ends(right)
            excluded=hi<ll or lo>rh
            inside=not excluded and lo>=lh and hi<=rl
            if not excluded:possible.append(chart)
            if inside:contained.append(chart)
            clipped=None if excluded else c.mpf([max(lo,ll),min(hi,rh)])
            rows[chart]=dict(exact_native_coordinate_function=vref,
                exact_native_domain_lower=route['lower'],exact_native_domain_upper=route['upper'],
                exact_source_log_radius_origin=route['origin'],exact_source_pulse_term=route['pulse_term'],
                exact_source_offset_zero=route['offset_zero'],exact_native_Jacobian=route['jacobian'],
                directed_native_coordinate=native,conditional_chart_native_intersection=clipped,
                native_difference_enclosure=delta,native_Jacobian_enclosure=jac,
                domain_lower_enclosure=left,domain_upper_enclosure=right,
                excluded_by_directed_bounds=excluded,entire_coordinate_enclosure_in_chart=inside,
                interior_resolved=inside and lo>lh and hi<rl,
                unbounded_heat_tail=chart=='heat_exterior',
                exact_common_factors_cancelled_before_numeric_arithmetic=True,
                clipped_box_is_conditional_chart_overlap_not_selected_coordinate=True)
        before=reader.at(g.sub(logR,self.radius.logRp))
        strictly_before=ends(before)[1]<0
        result=dict(source_family=self.family_record,exact_input_logR=logR,
            source_logR_enclosure=absolute,current_logRp_enclosure=reader.at(self.radius.logRp),
            logR_minus_current_logRp=before,chart_rows=rows,possible_charts=possible,contained_charts=contained,
            definitely_before_current_Rp=strictly_before,
            classification='before_current_Rp' if strictly_before else 'contained' if contained else 'unresolved_overlap',
            physical_Z_numeric_interior_resolved=physical_Z_resolved,
            native_source_interval_callback_available=False,source_point_value_installed=False,
            exact_seam_source_call_performed=False,overlap_charts_not_rounded_or_discarded=True,
            original_N_definition=self.before.original_N_definition,
            **dict.fromkeys(GATES,self.acceptance_loaded),**dict.fromkeys(OPEN,False))
        if not possible and not strictly_before:raise ArithmeticError('Current full affine cover lost a finite radius')
        self.call_trace.append(dict(classification=result['classification'],possible_charts=possible,
            no_source_packet_or_velocity_value_called=True))
        return result

    @source_precision
    def locate_inverse(self,view):
        """Revalidate an unchanged live inverse; reject imported scalar rows."""
        refs=view['coordinate_functions'];mapping=view['directed_inverse_mapping'];before=self.before
        proofs=[proof for proof in before._solve_witnesses.values() if proof[1] is mapping and proof[2] is refs]
        if len(proofs)!=1:raise ValueError('Live physical inverse from this accepted owner required')
        current=before.assemble(view['input_kind'],view['exact_input_record'],refs,mapping,
            view['physical_log_radius_enclosure'],view['theta_enclosure'],witness=proofs[0][0])
        if physical.mixed.pulse.raw.packed(inverse.report(current))!=physical.mixed.pulse.raw.packed(inverse.report(view)):
            raise ValueError('Unchanged current inverse view required')
        reader=self.reader({refs['log_lambda'].node:mapping['actual_log_lambda'],
            refs['log_r'].node:view['physical_log_radius_enclosure']})
        return self._locate(refs['logR'],reader,current['strict_numerical_Z_interior_resolved'])


def report(result):
    def nodes(value):
        if type(value) is radius_source.FunctionRef:return dict(exact_current_function_node=value.node)
        if isinstance(value,dict):return {key:nodes(item) for key,item in value.items()}
        if isinstance(value,(tuple,list)):return [nodes(item) for item in value]
        return value
    return physical.mixed.pulse.raw.packed(nodes(result))


@source_precision
def run(before=None):
    began=time.monotonic();owner=CurrentOriginalRpNativeChartLocator(before,require_checked=False)
    anchors={chart:owner.locate_function(owner.source_log_radius(chart,value))
        for chart,value in ANCHORS.items()}
    points={name:owner.before.cartesian(*args) for name,args in inverse.CASES.items()}
    points.update({name:owner.before.log_cylindrical(*args) for name,args in inverse.LOG_CASES.items()})
    located={name:owner.locate_inverse(view) for name,view in points.items()}
    result=dict(source_family=owner.family_record,actual_current_graph=owner.assert_graph(),
        source_bound_logRp_enclosure=owner.reader().at(owner.radius.logRp),
        actual_fifteen_source_affine_anchor_locations={key:report(value) for key,value in anchors.items()},
        actual_nine_physical_inverse_locations={key:report(value) for key,value in located.items()},
        exact_current_native_locator_graph=owner.graph.nodes,actual_source_call_trace=owner.call_trace,
        **dict.fromkeys(GATES+OPEN,False),input_hashes=owner.hashes,
        execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(report(result),separators=(',',':'))+'\n').encode(),mtime=0))
    print('CURRENT_ORIGINAL_RP_NATIVE_LOCATOR fifteen affine maps and physical inputs constructed',flush=True)
    return owner,anchors,points,located
