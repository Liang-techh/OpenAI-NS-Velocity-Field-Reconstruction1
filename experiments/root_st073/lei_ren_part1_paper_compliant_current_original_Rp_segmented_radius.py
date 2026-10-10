"""Absolute source radii for the current selected pulse and postpulse chain.

The accepted logRp, inverse-mu pulse length and finite offset remain separate
function expressions. Waiting is its true integral-defined function, never
a root-box endpoint. Numerical outputs are directed local geometry bounds;
no enormous absolute radius or Cartesian point field is materialized.
"""
import copy
from fractions import Fraction
import gzip
import json
from pathlib import Path
import time
import sympy as s
import lei_ren_part1_paper_compliant_current_original_Rp_postpulse_source as post
import lei_ren_part1_paper_compliant_current_original_C3_Rp_native_frame as frame_source
import lei_ren_part1_paper_compliant_current_original_outer_leading_C3_identity as leading_source
from lei_ren_part1_paper_compliant_current_original_outer_leading_C3_identity_check import Interpreter
from lei_ren_part1_paper_compliant_current_native_Rc_function_transport import FunctionRef
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

HERE,PREFIX,sha=post.HERE,post.PREFIX,post.sha
NAME=PREFIX+'current_original_Rp_segmented_radius.json.gz'
RECEIPT=PREFIX+'current_original_Rp_segmented_radius_check.json'
GATES=('current_original_Rp_absolute_segmented_radius_functions_installed',
       'current_original_Rp_actual_postpulse_source_radius_caller_installed')
OPEN=post.OPEN
VIEWS=(('flatten',100),('outer_power',0),('outer_power',1),('outer_angular',-4),
    ('outer_angular',0),('steep_entry',0),('steep_entry',1),('steep_power',0),
    ('steep_power',1),('steep_exit',0),('steep_exit',1),('waiting',0),('waiting',1),
    ('heat_collar',0),('heat_collar',3),('heat_exterior',3),('heat_exterior',4))


def exact_coordinate(value):
    if type(value) is int or type(value) is str or isinstance(value,Fraction):return Fraction(value)
    raise TypeError('An explicit exact rational native coordinate is required')


def ref(g,value):return FunctionRef(g,value.node)


def waiting_functions(owner):
    """Literal native affine terminal and uncorrected waiting construction.

    Angular repair coefficients act later and do not redefine this waiting
    datum. All kernels below are mathematical integrals of original shapes.
    """
    g=owner.graph;leading=leading_source
    mu=owner.parameters['mu'];logP=owner.parameters['logP']
    rate=g.sub(g.one,mu);logdelta=g.add(g.mul(g.constant(-4),logP),g.constant(-30))
    delta=g.unary('exp',logdelta);epsilon=g.quotient(delta,g.constant(1000),'positive constant 1000')
    k=g.sub(g.one,g.mul(g.constant('1/2'),delta));eq=g.quotient(g.one,k,'admitted positive restore rate')
    aeq=g.quotient(g.one,rate,'admitted positive pulse rate')
    logmu=g.unary('log',mu);L=g.mul(g.constant(-30),logmu)
    Ts=g.mul(g.constant(4),g.sub(g.unary('log',g.constant(2)),logdelta))
    frame=owner.frame.functions['actual_identified_Rp_native_frame_C3']
    Xp=g.node('function_substitution',expression=frame_source.target.rows(frame['X'])[0].node,
        variable=g.symbol('Z').node,value=g.zero.node,Z_independent_substitution=True)
    decay=g.unary('exp',g.neg(g.quotient(g.mul(g.constant(13),rate),mu,'admitted positive mu')))
    Xv=g.add(aeq,g.mul(g.sub(Xp,aeq),decay))
    u=g.symbol('current_Rp_wait_kernel_u');v=g.symbol('current_Rp_wait_kernel_v')
    def sigma(x):return leading.sigma(g,x)
    def integral(body,name,upper,lower=None):return leading.integral(g,body,name,upper,lower)
    J=integral(sigma(v),'current_Rp_wait_kernel_v',u)
    F=g.unary('exp',g.neg(g.mul(g.unary('log',g.constant(2)),
        sigma(g.quotient(u,g.constant(100),'positive flatten length 100')))))
    flatten_mass=integral(g.mul(F,g.unary('exp',g.neg(g.mul(rate,g.sub(g.constant(100),u))))),
        'current_Rp_wait_kernel_u',g.constant(100))
    Xf=g.mul(g.constant(2),g.add(g.mul(Xv,g.unary('exp',g.mul(g.constant(-100),rate))),flatten_mass))
    Xrel=g.add(aeq,g.mul(g.sub(Xf,aeq),g.unary('exp',g.neg(g.mul(rate,L)))))
    Iin=integral(g.unary('exp',g.mul(rate,g.sub(u,J))),'current_Rp_wait_kernel_u',g.one)
    Iout=integral(g.unary('exp',g.mul(k,J)),'current_Rp_wait_kernel_u',g.one)
    Xt=g.mul(g.add(g.mul(g.add(Xrel,Iin),g.unary('exp',g.mul(g.constant('-1/2'),rate))),Ts,Iout),
        g.unary('exp',g.mul(g.constant('-1/2'),k)))
    edge=g.node('exact_original_positive_flat_edge',argument=g.quotient(g.sub(g.constant(3),u),g.constant(2),
        'positive constant 2').node,definition='exp(-1/x^2) for x>0, zero for x<=0',
        source_module=PREFIX+'outer_angular_candidate.py',source_function='collar_preheat_integral.f')
    sig=sigma(u)
    collarJ=integral(g.mul(g.unary('exp',g.mul(k,u)),g.add(g.sub(g.one,sig),g.mul(sig,edge))),
        'current_Rp_wait_kernel_u',g.constant(3))
    logone=g.unary('log',g.sub(g.one,epsilon))
    numerator=g.add(g.unary('log',g.sub(Xt,eq)),logone,g.neg(g.unary('log',epsilon)),
        g.neg(g.unary('log',g.add(eq,collarJ))))
    waiting=g.quotient(numerator,k,'accepted current positive waiting-source restore rate')
    return dict(mu=mu,logP=logP,logmu=logmu,logdelta=logdelta,delta=delta,epsilon=epsilon,
        rate=rate,k=k,Lrel=L,Ts=Ts,Xp=Xp,Xv=Xv,flatten_mass=flatten_mass,Iin=Iin,Iout=Iout,
        collarJ=collarJ,Xt=Xt,waiting_logone=logone,waiting=waiting)


class RadiusInterpreter(Interpreter):
    def at(self,i):
        i=i.node if hasattr(i,'node') else i
        if i in self.bindings:return self.bindings[i]
        if self.g.nodes[i]['operation']=='exact_original_positive_flat_edge':
            n=self.g.nodes[i];assert n['definition']=='exp(-1/x^2) for x>0, zero for x<=0'
            x=self.at(n['argument'])
            return s.Piecewise((s.exp(-1/x**2),x>0),(0,True))
        return super().at(i)


class CurrentOriginalRpSegmentedRadius:
    @source_precision
    def __init__(self,before=None,require_checked=True):
        self.before=before if before is not None else post.CurrentOriginalRpPostpulseSource()
        if type(self.before) is not post.CurrentOriginalRpPostpulseSource or not self.before.acceptance_loaded:
            raise ValueError('Accepted typed current postpulse owner required')
        self.before.assert_graph();self.ctx=self.before.ctx
        self.family_record=self.before.family_record;self.hashes=dict(self.before.hashes)
        self.frame=self.before.before.inlet.exact.frame
        self.graph=copy.deepcopy(self.frame.graph);self.prefix=copy.deepcopy(self.graph.nodes)
        self.parameters={key:ref(self.graph,value) for key,value in self.frame.bridge.leading.parameters.items()}
        self.functions=waiting_functions(self);g=self.graph;f=self.functions
        self.logRp=ref(g,self.before.before.current_Rp_physical_log_radius)
        self.pulse_length=g.quotient(g.constant(13),f['mu'],'accepted current positive mu')
        self.native=g.symbol('current_Rp_radius_coordinate')
        self.maps={chart:self._map(chart,self.native) for chart in self.before.registry}
        self.call_trace=[];self.acceptance_loaded=False
        for name in (post.NAME,post.RECEIPT,Path(__file__).name):self.hashes[name]=sha(name)
        self.assert_graph()
        if require_checked:
            record=json.loads((HERE/RECEIPT).read_bytes())
            if not record['all_passed'] or not all(record[key] for key in GATES) or any(record[key] for key in OPEN):
                raise ValueError('Accepted current segmented-radius receipt required')
            if record['source_family']!=self.family_record:raise ValueError('Segmented radius family differs')
            post.selected.inlet.add_hashes(self.hashes,record['input_hashes'])
            self.acceptance_loaded=True

    def assert_graph(self):
        b=self.before;f=b.before;g=self.graph
        graph=dict(current_postpulse_owner=b.acceptance_loaded and all(b.assert_graph().values()),
            same_accepted_original_logRp=self.logRp.node==f.current_Rp_physical_log_radius.node,
            unchanged_accepted_graph_prefix=g.nodes[:len(self.prefix)]==self.prefix==self.frame.graph.nodes,
            separate_expression_graph=g is not self.frame.graph,
            current_waiting_owner=b.heat.future.angular is b.steep.future.angular is f.future.angular,
            current_unique_repair=b.heat.repair is f.future.repair is f.angular4.repair,
            current_native_geometry_context=b.outer.ctx is b.steep.ctx is b.heat.ctx is self.ctx,
            all_fifteen_radius_routes=set(self.maps)==set(b.registry),
            no_original_N_replaced=True)
        if not all(graph.values()):raise ValueError('Current radius source graph differs: '+str(graph))
        return graph

    def _map(self,chart,v):
        g=self.graph;f=self.functions;one=g.one;L=f['Lrel'];T=f['Ts'];W=f['waiting']
        p=g.zero;jac=one
        if chart=='pulse_entrance':offset=v
        elif chart in ('pulse_main','pulse_exit','pulse_gap'):
            offset=g.quotient(v,f['mu'],'accepted current positive mu');jac=g.quotient(one,f['mu'],'accepted current positive mu')
        elif chart in ('pulse_gap_end','pulse_end'):p=self.pulse_length;offset=v
        else:
            p=self.pulse_length
            if chart=='flatten':offset=v
            elif chart=='outer_power':
                jac=g.sub(L,g.constant(4));offset=g.add(g.constant(100),g.mul(jac,v))
            elif chart=='outer_angular':offset=g.add(g.constant(100),L,v)
            elif chart=='steep_entry':offset=g.add(g.constant(100),L,v)
            elif chart=='steep_power':jac=T;offset=g.add(g.constant(101),L,g.mul(T,v))
            elif chart=='steep_exit':offset=g.add(g.constant(101),L,T,v)
            elif chart=='waiting':jac=W;offset=g.add(g.constant(102),L,T,g.mul(W,v))
            elif chart in ('heat_collar','heat_exterior'):offset=g.add(g.constant(102),L,T,W,v)
            else:raise ValueError('Unknown current radius chart: '+chart)
        logR=g.add(self.logRp,p,offset)
        return dict(origin=self.logRp,pulse_term=p,offset=offset,logR=logR,R=g.unary('exp',logR),
            native_to_log_radius_jacobian=jac)

    def geometry(self,chart,coordinate):
        self.assert_graph()
        if chart not in self.maps:raise ValueError('Unknown current radius chart')
        value=exact_coordinate(coordinate)
        domains={'pulse_main':(Fraction(1,50),10),'pulse_exit':(10,11),'pulse_gap':(11,12),
            'pulse_end':(-4,0),'flatten':(0,100),'outer_angular':(-4,0),
            'heat_collar':(0,3),'heat_exterior':(3,None)}
        if chart=='pulse_entrance':
            upper=post.selected.inlet.endpoints(self.ctx.mpf('.02')/self.before.before.pulse.mu)[0]
            if value<0 or self.ctx.mpf(value.numerator)/value.denominator>upper:
                raise ValueError('Current pulse entrance source domain required')
        elif chart=='pulse_gap_end':
            lower=post.selected.inlet.endpoints(-1/self.before.before.pulse.mu)[1]
            if value>-4 or self.ctx.mpf(value.numerator)/value.denominator<lower:
                raise ValueError('Current pulse gap-end source domain required')
        else:
            left,right=domains.get(chart,(0,1))
            if value<left or right is not None and value>right:raise ValueError('Current original radius domain required')
        actual=self._map(chart,self.graph.constant(value))
        return dict(chart=chart,exact_native_coordinate=dict(numerator=value.numerator,denominator=value.denominator),
            **{key:ref.node for key,ref in actual.items()},
            absolute_radius_is_a_lazy_function=True,short_offset_not_added_to_huge_absolute_value=True,
            original_waiting_integrals_not_replaced_by_root_box=True,
            numerical_absolute_radius_point_value_installed=False)

    @source_precision
    def local_step(self,chart,left,right):
        """A local log-radius increment with its true directed Jacobian bound.

        Absolute logRp and the shared pulse term cancel before arithmetic.
        Enclosures bound functions and are never used as their definitions.
        """
        a,b=exact_coordinate(left),exact_coordinate(right);self.geometry(chart,a);self.geometry(chart,b)
        c=self.ctx;before=self.before
        bound=(1/before.before.pulse.mu if chart in ('pulse_main','pulse_exit','pulse_gap') else
            before.outer.Lrel-4 if chart=='outer_power' else before.steep.Ts if chart=='steep_power' else
            before.steep.wait if chart=='waiting' else c.mpf(1))
        if post.selected.inlet.endpoints(bound)[0]<=0:raise ArithmeticError('Positive source coordinate Jacobian required')
        step=c.mpf((b-a).numerator)/(b-a).denominator
        exact=self.graph.mul(self.maps[chart]['native_to_log_radius_jacobian'],self.graph.constant(b-a))
        return dict(exact_log_radius_increment=exact.node,native_to_log_radius_jacobian_bound=bound,
            log_radius_increment_bound=bound*step,constant_origin_cancelled_before_arithmetic=True,
            enclosure_not_used_as_defining_function=True)

    @source_precision
    def evaluate(self,chart,Z,coordinate):
        result=dict(self.before.evaluate(chart,Z,coordinate))
        result['current_absolute_source_geometry']=self.geometry(chart,coordinate)
        self.call_trace.append(dict(chart=chart,actual_source_provider_unchanged=True,
            actual_absolute_source_radius_mapper_called=True))
        result.update(**dict.fromkeys(GATES,self.acceptance_loaded),**dict.fromkeys(OPEN,False))
        return result


@source_precision
def run(before=None):
    began=time.monotonic();owner=CurrentOriginalRpSegmentedRadius(before,require_checked=False)
    views=[]
    for chart,coordinate in VIEWS:
        view=owner.evaluate(chart,'.427',coordinate)
        views.append(dict(chart=chart,coordinate=coordinate,geometry=view['current_absolute_source_geometry']))
        print('Actual current segmented radius:',chart,coordinate,flush=True)
    local={chart:owner.local_step(chart,3,'7/2') if chart=='heat_exterior' else owner.local_step(chart,0,'1/2') for chart in
        ('flatten','outer_power','steep_entry','steep_power','steep_exit','waiting','heat_collar','heat_exterior')}
    result=dict(source_family=owner.family_record,candidate_current_segmented_radius_constructed=True,
        actual_same_object_graph=owner.assert_graph(),accepted_original_graph_prefix_length=len(owner.prefix),
        exact_function_graph_nodes=owner.graph.nodes,exact_waiting_functions={key:value.node for key,value in owner.functions.items()},
        current_absolute_radius_maps={chart:{key:value.node for key,value in row.items()} for chart,row in owner.maps.items()},
        actual_seventeen_current_radius_calls=views,actual_source_call_trace=owner.call_trace,
        actual_local_radius_steps=local,original_logRp_node=owner.logRp.node,
        no_numeric_absolute_radius_or_new_N_materialized=True,
        **dict.fromkeys(GATES+OPEN,False),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began)
    data=json.dumps(post.selected.inlet.encode(post.selected.inlet.pack(result)),separators=(',',':'))+'\n'
    (HERE/NAME).write_bytes(gzip.compress(data.encode(),mtime=0))
    return owner


if __name__=='__main__':run()
