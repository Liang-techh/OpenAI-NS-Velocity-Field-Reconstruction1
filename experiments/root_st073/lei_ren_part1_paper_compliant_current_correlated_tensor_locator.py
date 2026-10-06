"""Source-correlated physical radius requests on the checked complete graph.

Absolute input boxes remain supported by the previous locator. This API adds
physical point families expressed relative to the actual construction's
common radius anchors; it preserves their correlations before inversion.
"""
from dataclasses import dataclass
import functools
import json
import re
from pathlib import Path
import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_current_physical_tensor_locator import (
    CurrentPhysicalTensorLocator,HERE,PREFIX,sha,pack,encode,endpoints,source_precision,
    accepted,_verify_hashes,OPEN)
from lei_ren_part1_paper_compliant_current_physical_tensor_locator_operator import (
    V,SYMBOLS,DOMAINS,implicit_log_coordinate_map,intersect,nonpositive_exp)
from lei_ren_part1_paper_compliant_current_correlated_radius_operator import (
    OFFSET,LOG_HB,canonical_radius,source_geometry,anchors,boundary_sources,
    inverse_source_expression,correlated_source_theorem)

NAME=PREFIX+'current_correlated_tensor_locator.json'
RECEIPT=PREFIX+'current_correlated_tensor_locator_check.json'
GATES=('current_correlated_source_radius_inverse_locator_available',
       'current_46_exact_source_radius_boundary_trace_selectors_available')


@dataclass(frozen=True)
class SourceRadiusRequest:
    expression:s.Expr
    offset_box:object
    family:str
    source:str
    datum:str
    original_radius_program_sha256:str


def exact_expression(value):
    if isinstance(value,s.Basic):return value
    if isinstance(value,int):return s.Integer(value)
    if isinstance(value,str) and re.fullmatch(r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?',value):
        return s.Rational(value)
    raise ValueError('Exact Sympy radius expression, integer or exact decimal string required')


class CurrentCorrelatedTensorLocator:
    @source_precision
    def __init__(self,locator=None,require_checked=True):
        self.locator=locator if locator is not None else CurrentPhysicalTensorLocator()
        if type(self.locator) is not CurrentPhysicalTensorLocator or not self.locator.acceptance_loaded:
            raise ValueError('Checked actual physical coordinate/radius locator required')
        self.registry=self.locator.registry;self.physical=self.locator.physical;self.ctx=self.locator.ctx
        self.family=self.locator.family;self.source=self.locator.source;self.datum_sha=self.locator.datum_sha
        self.geometry,self.radius_program_proof=source_geometry();self.boundaries=boundary_sources()
        if set(self.boundaries)!=set(self.registry.adjacent)|set(self.registry.internal):
            raise ValueError('Every source boundary must have its admitted complete trace route')
        self.proof=correlated_source_theorem()
        self.radius_program_sha=self.locator.radius_program_proof['input_hashes'][PREFIX+'global_physical_assembly.py']
        self.hashes={**self.locator.hashes,**self.proof['input_hashes']}
        for stem in ('current_correlated_radius_operator','current_correlated_tensor_locator'):
            self.hashes[PREFIX+stem+'.py']=sha(PREFIX+stem+'.py')
        self.acceptance_loaded=False;self.assert_graph()
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source,GATES[0]);_verify_hashes(receipt)
            if receipt['datum_enclosure_sha256']!=self.datum_sha or not all(receipt[k] for k in GATES) or any(receipt[k] for k in OPEN):
                raise ValueError('Correlated source radius admission exceeds scope')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.acceptance_loaded=True

    def assert_graph(self):
        self.locator.assert_graph()
        if not self.locator.acceptance_loaded or self.registry is not self.locator.registry or self.geometry is not source_geometry()[0] or self.boundaries is not boundary_sources():
            raise ValueError('Foreign correlated radius/operator/physical graph')

    @source_precision
    def radius_expression(self,expression,offset_box=None):
        self.assert_graph();expr=exact_expression(expression)
        allowed=set(SYMBOLS.values())|{LOG_HB,OFFSET}
        if expr.free_symbols-allowed or expr.has(s.Float,s.zoo,s.nan) or expr==s.oo:
            raise ValueError('Only the original exact source radius parameters and requested offset are allowed')
        if (OFFSET in expr.free_symbols)!=(offset_box is not None):raise ValueError('Exactly one enclosure for the requested offset variable required')
        box=None if offset_box is None else self.ctx.mpf(offset_box)
        if box is not None and any(not mp.isfinite(v) for v in endpoints(box)):raise ValueError('Finite requested radius offset box required')
        request=SourceRadiusRequest(canonical_radius(expr),box,self.family,self.source,self.datum_sha,self.radius_program_sha)
        # Domain validation for original logarithms/reciprocals. HB is
        # positive by its exact finite log source, even when its cap touches0.
        self._validate_request(request)
        return request

    def _validate_request(self,request):
        values=self._values(request)
        for node in s.preorder_traversal(request.expression):
            if node.func is s.log:
                # An arbitrary box touching zero cannot define a finite log
                # source. Exact exp(logh)>0 still passes even when its numeric
                # cap touches zero; no rounded lower bound substitutes for it.
                if self._sign(node.args[0],values)!=1:raise ValueError('Strictly positive exact logarithm argument required for source radius')
            if isinstance(node,s.Pow) and node.exp.is_negative and self._sign(node.base,values) not in (-1,1):
                raise ValueError('Nonzero exact source denominator required')
            if isinstance(node,s.Pow) and not node.exp.is_Integer and endpoints(self._eval(node.base,values))[0]<0:
                raise ValueError('Nonnegative real fractional-power source base required')
        self._eval(request.expression,values)

    @source_precision
    def radius(self,anchor,offset=0):
        if anchor not in anchors():raise ValueError('Unknown original source radius anchor')
        if isinstance(offset,(tuple,list)):
            return self.radius_expression(anchors()[anchor]+OFFSET,offset)
        return self.radius_expression(anchors()[anchor]+exact_expression(offset))

    def _values(self,request):
        if type(request) is not SourceRadiusRequest or (request.family,request.source,request.datum,request.original_radius_program_sha256)!=(self.family,self.source,self.datum_sha,self.radius_program_sha):
            raise ValueError('Foreign source radius request')
        allowed=set(SYMBOLS.values())|{LOG_HB,OFFSET}
        if not isinstance(request.expression,s.Basic) or request.expression.free_symbols-allowed or request.expression.has(s.Float,s.zoo,s.nan) or request.expression==s.oo or request.expression!=canonical_radius(request.expression):
            raise ValueError('Original canonical exact source radius expression required')
        if (OFFSET in request.expression.free_symbols)!=(request.offset_box is not None):raise ValueError('Source radius offset binding differs')
        if request.offset_box is not None and any(not mp.isfinite(v) for v in endpoints(request.offset_box)):
            raise ValueError('Finite source radius offset binding required')
        return {**self.locator.values,LOG_HB:self.locator.value_bindings['exact_positive_hb_log'],
            **({OFFSET:request.offset_box} if request.offset_box is not None else {})}

    def _eval(self,expression,values):
        c=self.ctx
        if expression in values:return values[expression]
        if expression in (s.oo,-s.oo):return c.mpf('inf' if expression==s.oo else '-inf')
        if expression is s.E:return c.exp(c.mpf(1))
        if expression.is_Rational:return c.mpf(str(expression.p))/int(expression.q)
        if isinstance(expression,s.Add):return sum((self._eval(arg,values) for arg in expression.args),c.mpf(0))
        if isinstance(expression,s.Mul):
            result=c.mpf(1)
            for arg in expression.args:result*=self._eval(arg,values)
            return result
        if isinstance(expression,s.Pow):return self._eval(expression.base,values)**self._eval(expression.exp,values)
        if expression.func is s.log:return c.ln(self._eval(expression.args[0],values))
        if expression.func is s.exp:
            arg=self._eval(expression.args[0],values);lo,hi=endpoints(arg)
            if hi<=0:return nonpositive_exp(c,arg)
            if hi>1000000:
                # Original finite source expression retained. This bound is
                # only an enclosure; no gigantic exponent integer is formed.
                lower=0 if lo<0 else endpoints(c.exp(c.mpf(min(lo,mp.mpf(1000000)))))[0]
                return c.mpf([lower,mp.inf])
            if lo<-1024:return c.mpf([0,endpoints(c.exp(c.mpf(hi)))[1]])
            return c.exp(arg)
        raise ValueError('Unsupported original correlated radius expression: '+str(expression))

    def _sign(self,expression,values):
        expression=s.sympify(expression)
        if expression==0:return 0
        if expression in (s.oo,-s.oo):return 1 if expression==s.oo else -1
        if expression in (SYMBOLS['h_bridge'],SYMBOLS['h_switch']):return 1
        if expression.func is s.exp:return 1
        if isinstance(expression,s.Mul):
            signs=[self._sign(arg,values) for arg in expression.args]
            if 0 in signs:return 0
            if None not in signs:
                result=1
                for sign in signs:result*=sign
                return result
        if isinstance(expression,s.Pow) and expression.exp.is_Integer:
            sign=self._sign(expression.base,values)
            if sign in (-1,1):return 1 if int(expression.exp)%2==0 else sign
        box=self._eval(expression,values);lo,hi=endpoints(box)
        if lo>0:return 1
        if hi<0:return -1
        factored=s.factor(expression)
        if factored!=expression and isinstance(factored,s.Mul):return self._sign(factored,values)
        return None

    def _candidate(self,region,request,values):
        c=self.ctx;row=self.geometry[region];expr=request.expression
        left=s.expand(expr-row['lower']);right=s.expand(expr-row['upper'])
        ls=self._sign(left,values);rs=self._sign(right,values)
        if ls==-1 or rs==1:return None
        lo,hi=row['domain'];domain=c.mpf([endpoints(self._eval(lo,values))[0],endpoints(self._eval(hi,values))[1]])
        native,extra=inverse_source_expression(region,expr)
        status='original_exact_correlated_source_inverse'
        if region in ('core_positive_radius','actual_patch'):
            log_domain=c.ln(domain);log_box=intersect(c,self._eval(extra,values),log_domain)
            if log_box is None:return None
            shift=endpoints(log_box)[1];coordinate=c.exp(c.mpf(shift))*nonpositive_exp(c,log_box-shift)
        elif region=='O2_axial':
            argument=intersect(c,self._eval(extra,values),c.mpf([1,endpoints(c.exp(values[SYMBOLS['Md']]))[1]]))
            if argument is None:return None
            # Symbolic log(exp(Md*t))/Md can cancel exactly before the box.
            coordinate=self._eval(native,values) if native.func is not s.log and not native.has(s.log) else c.ln(argument)/values[SYMBOLS['Md']]
        else:
            coordinate=self._eval(native,values)
        coordinate=intersect(c,coordinate,domain)
        if coordinate is None:
            # An unresolved sign/denominator enclosure must not erase a
            # possible source region. Exact source domains remain the bound.
            if ls is None or rs is None:coordinate=domain;status='source_inverse_precision_limited_full_legal_domain'
            else:raise ValueError('Correlated source inverse contradicts exact radius distances')
        return dict(region=region,native_coordinate_enclosure=coordinate,exact_native_coordinate_source=str(native),
            exact_lower_radius_distance=str(left),exact_upper_radius_distance=str(right),
            lower_distance_sign=ls,upper_distance_sign=rs,inverse_status=status,
            cancellation_precedes_parameter_enclosure=True,source_cap_not_promoted_to_field=True)

    @source_precision
    def locate(self,request,z_phys,log_tau,theta=0,viscosity=1):
        self.assert_graph();self._validate_request(request);values=self._values(request);c=self.ctx;angle=c.mpf(theta)
        if any(not mp.isfinite(v) for v in endpoints(angle)):raise ValueError('Finite physical angle required')
        mapping=implicit_log_coordinate_map(c,z_phys,log_tau,viscosity,self.physical.delta)
        if request.expression==-s.oo:
            candidates=[dict(region='core_positive_radius',native_coordinate_enclosure=c.mpf(0),exact_native_coordinate_source='0',
                inverse_status='exact_nonsingular_axis',cancellation_precedes_parameter_enclosure=True,source_cap_not_promoted_to_field=True)]
        else:candidates=[item for region in self.geometry if (item:=self._candidate(region,request,values)) is not None]
        if not candidates:raise ValueError('No original source region encloses the correlated physical radius')
        matches=[name for name,row in self.boundaries.items() if s.expand(request.expression-row['expression'])==0]
        logR=self._eval(request.expression,values)
        lr=(c.ln(2)+c.ln(mapping['physical_viscosity'])+2*mapping['actual_log_lambda']+logR)/2
        return dict(**mapping,theta=angle,source_logR_enclosure=logR,physical_log_radius_enclosure=lr,
            exact_source_logR=str(request.expression),requested_offset_enclosure=request.offset_box,
            exact_physical_log_radius_recipe='(log2+lognu+2*actual_log_lambda+exact_source_logR)/2',
            request_semantics='Source-dependent physical point family; not a fixed absolute coordinate',
            candidates=candidates,candidate_region_count=len(candidates),region_ambiguity_retained=len(candidates)>1,
            exact_source_boundary_matches=matches,
            source_family=self.family,implicit_source_sha256=self.source,datum_enclosure_sha256=self.datum_sha,
            **dict.fromkeys(OPEN,False))

    def _replay(self,request,location):
        expected=self.locate(request,location['physical_z'],location['requested_log_tau'],location['theta'],location['physical_viscosity'])
        if encode(pack(expected))!=encode(pack(location)):raise ValueError('Foreign or changed correlated source locator result')

    @source_precision
    def tensor(self,request,z_phys,log_tau,theta=0,viscosity=1):
        location=self.locate(request,z_phys,log_tau,theta,viscosity);c=self.ctx;views=[]
        for candidate in location['candidates']:
            region=candidate['region'];v=candidate['native_coordinate_enclosure']
            if region=='core_positive_radius':
                root=c.sqrt(2*v);X=root*c.cos(location['theta']);Y=root*c.sin(location['theta'])
                source=self.locator.core_program(self.registry.owners['axis'],location['Z'],X,Y,
                    location['requested_log_tau'],location['physical_viscosity'],rho_source=v)
                view=self.registry._view(region,source,'Cartesian_core')
                view['joint_polar_Cartesian_source_constraint']=self.locator.core_program_proof
            elif region=='heat_exterior':view=self.registry.exterior(location['Z'],location['requested_log_tau'],location['theta'],location['physical_viscosity'])
            else:view=self.registry.native(region,location['Z'],v,location['requested_log_tau'],location['theta'],location['physical_viscosity'])
            views.append(dict(candidate=candidate,original_full_native_view=view,
                actual_lambda_signed_component_groups=self.locator._refined_groups(view,location['actual_log_lambda'])))
        return dict(location=location,candidate_full_tensor_views=views,
            source_radius_correlation_retained_before_native_enclosure=True,
            return_kind='Complete original physical T/E source-function enclosure for the correlated point family',
            **dict.fromkeys(OPEN,False))

    @source_precision
    def traces(self,request,location):
        self._replay(request,location);result={}
        for name in location['exact_source_boundary_matches']:
            result[name]=self.registry.trace(name,location['Z'],location['requested_log_tau'],location['theta'],location['physical_viscosity'])
        return dict(exact_original_source_radius_trace_routes=result,actual_log_lambda=location['actual_log_lambda'],
            actual_log_tau_retained=True,common_trace_bound_retained_as_conservative_native_sqrt_tau=True,
            source_function_equality_not_radius_interval_overlap=True,**dict.fromkeys(OPEN,False))

    def manifest(self):
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,datum_enclosure_sha256=self.datum_sha,
            current_exact_correlated_radius_source_theorem=self.proof,
            original_source_radius_anchors={name:str(value) for name,value in anchors().items()},
            current_exact_source_radius_boundary_registry={name:{key:(str(value) if isinstance(value,s.Basic) else value)
                for key,value in row.items()} for name,row in self.boundaries.items()},
            actual_current_tensor_regions_available=list(self.geometry),
            scope='Exact original correlated radius expressions before enclosure, 33 inverse candidate sources and 46 exact radius trace selectors on the checked complete tensor graph. Requests denote source-dependent physical point families. Arbitrary absolute-coordinate seam/global physical cover certification, resolved velocity/pressure points, cone/lift, energy, corrected NS and completed flatness/recursion remain open.',
            input_hashes=self.hashes,**dict.fromkeys(GATES+OPEN,False))


@source_precision
def run(field=None):
    field=field if field is not None else CurrentCorrelatedTensorLocator(require_checked=False)
    field.assert_graph();result=field.manifest()
    (HERE/NAME).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current exact correlated radius tensor locator: common source cancellation, 33 inverses and 46 exact trace selectors',flush=True)
    return result


if __name__=='__main__':run()
