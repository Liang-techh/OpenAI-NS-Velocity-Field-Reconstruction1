"""Free-variable cache for nested exact repair-band function quadrature."""
import lei_ren_part1_paper_compliant_current_native_Rc_functional_controls as controls


class BandFunctionEvaluator(controls.FunctionEvaluator):
    """Respect integral binding and retain only relevant variables in keys.

    Fixed raw-bump integrals and Picard coefficients are reused at outer
    partial-history quadrature points. The explicit original/synthetic
    oracle contract and rejection of range-as-value data are unchanged.
    This is an evaluation cache, not a quadrature or fixed-point certificate.
    """
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        self.free_variables=[]
        for row in self.built['graph'].nodes:
            op=row['operation']
            union=lambda ids:frozenset().union(*(self.free_variables[i] for i in ids))
            if op=='bound_variable':free=frozenset((row['name'],))
            elif op in ('exact_rational','original_source_parameter','shared_positive_integer'):free=frozenset()
            elif op in ('sum','product'):free=union(row['arguments'])
            elif op in ('negative','analytic_unary','compact_raw_beta'):free=self.free_variables[row['argument']]
            elif op=='positive_quotient':free=union((row['numerator'],row['denominator']))
            elif op=='original_function_graph':free=union((row['coordinate'],row['phase']))
            elif op=='definite_integral':
                free=(self.free_variables[row['integrand']]-{row['variable']})|union((row['lower'],row['upper']))
            else:raise ValueError('Unknown function operation for free-variable cache: '+op)
            self.free_variables.append(free)

    def walk(self,index,variables):
        return super().walk(index,{key:value for key,value in variables.items() if key in self.free_variables[index]})
