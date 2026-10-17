"""Small exact scalar expression graph over pinned original C5 functions.

References remain functions from their original source graph. No enclosure,
cap or decimal approximation is registered as an original function value.
"""
from __future__ import annotations

from fractions import Fraction
import json


class ScalarExpression:
    __slots__=('graph','node')
    def __init__(self,graph,node):self.graph,self.node=graph,node
    def __add__(self,other):return self.graph.binary('sum',self,other)
    def __radd__(self,other):return self.graph.binary('sum',other,self)
    def __sub__(self,other):return self.graph.binary('difference',self,other)
    def __rsub__(self,other):return self.graph.binary('difference',other,self)
    def __mul__(self,other):return self.graph.binary('product',self,other)
    def __rmul__(self,other):return self.graph.binary('product',other,self)
    def __truediv__(self,other):return self.graph.binary('quotient',self,other)
    def __rtruediv__(self,other):return self.graph.binary('quotient',other,self)
    def __neg__(self):return self.graph.unary('negative',self)
    def __pow__(self,power):
        if type(power) is not int:raise TypeError('Only exact integer powers are used here')
        if power==0:return self.graph.constant(1)
        if power<0:return self.graph.constant(1)/(self**(-power))
        if power==1:return self
        return self.graph.insert({'operation':'integer_power','argument':self.node,'power':power})


class SignedDensityExpressionGraph:
    def __init__(self):
        self.nodes=[];self._index={}
        self.zero=self.constant(0);self.one=self.constant(1)

    def insert(self,node):
        key=json.dumps(node,sort_keys=True,separators=(',',':'))
        if key not in self._index:self._index[key]=len(self.nodes);self.nodes.append(node)
        return ScalarExpression(self,self._index[key])

    def constant(self,value):
        q=Fraction(value)
        return self.insert({'operation':'exact_rational','numerator':q.numerator,'denominator':q.denominator})

    def coerce(self,value):
        if isinstance(value,ScalarExpression):
            if value.graph is not self:raise ValueError('Original source namespaces must not be mixed')
            return value
        return self.constant(value)

    def reference(self,report,sha256,node,*,source_nodes):
        if type(node) is not int or not 0<=node<len(source_nodes):raise ValueError('Original function node missing')
        original=source_nodes[node]
        if not isinstance(original,dict) or 'operation' not in original:raise ValueError('Original source function node required')
        if any(t in original['operation'].lower() for t in ('cap','enclosure','majorant','bound_only',
            'magnitude','upper_bound','log_bound','norm_bound','certificate','proved_')):
            raise ValueError('A magnitude bound cannot be a source function')
        return self.insert({'operation':'original_function_reference','report':report,'report_sha256':sha256,
            'original_node':node,'original_operation':original['operation']})

    def unary(self,operation,value):
        x=self.coerce(value)
        if operation=='negative':
            node=self.nodes[x.node]
            if node['operation']=='exact_rational':return self.constant(-Fraction(node['numerator'],node['denominator']))
            if node['operation']=='negative':return ScalarExpression(self,node['argument'])
        return self.insert({'operation':operation,'argument':x.node})

    def binary(self,operation,left,right):
        a,b=self.coerce(left),self.coerce(right)
        A,B=self.nodes[a.node],self.nodes[b.node]
        rational=lambda v:Fraction(v['numerator'],v['denominator'])
        if A['operation']==B['operation']=='exact_rational':
            x,y=rational(A),rational(B)
            return self.constant({'sum':lambda:x+y,'difference':lambda:x-y,'product':lambda:x*y,'quotient':lambda:x/y}[operation]())
        zero=lambda v:v['operation']=='exact_rational' and v['numerator']==0
        one=lambda v:v['operation']=='exact_rational' and v['numerator']==v['denominator']
        if operation=='sum':
            if zero(A):return b
            if zero(B):return a
        if operation=='difference':
            if a.node==b.node:return self.zero
            if zero(B):return a
            if zero(A):return -b
        if operation=='product':
            if zero(A) or zero(B):return self.zero
            if one(A):return b
            if one(B):return a
        if operation=='quotient':
            if one(B):return a
            # Do not cancel x/x without a separate nonzero source theorem.
        args=[a.node,b.node]
        if operation in ('sum','product'):args.sort()
        return self.insert({'operation':operation,'arguments':args})

    def evaluate(self,bindings,*,through=None):
        """Only for NEW bounded tests, with independently supplied values."""
        values=[]
        for node in self.nodes if through is None else self.nodes[:through+1]:
            op=node['operation']
            if op=='exact_rational':v=Fraction(node['numerator'],node['denominator'])
            elif op=='original_function_reference':v=bindings[(node['report'],node['original_node'])]
            elif op=='negative':v=-values[node['argument']]
            elif op=='integer_power':v=values[node['argument']]**node['power']
            elif op=='exp':raise ValueError('Original physical exponentials are never evaluated here')
            else:
                a,b=(values[i] for i in node['arguments'])
                v={'sum':lambda:a+b,'difference':lambda:a-b,'product':lambda:a*b,'quotient':lambda:a/b}[op]()
            values.append(v)
        return values

    def packet(self):
        return {'exact_local_expression_nodes':self.nodes,'local_nodes':len(self.nodes),
            'original_graph_nodes_added_or_replaced':False,'bounds_or_midpoints_registered_as_values':False,
            'original_sources_numerically_executed':False,'physical_scale_materialized':False}
