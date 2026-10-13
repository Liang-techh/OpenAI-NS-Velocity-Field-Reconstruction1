"""Independent physical-measure and inverse-action checks for new n1 bumps."""
import gzip
import json
import math
from pathlib import Path

import mpmath as mp
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_n1_reserved_linear_operator as current

HERE,PREFIX,sha,ep=current.HERE,current.PREFIX,current.sha,current.ep
NAME=current.RECEIPT


def contains(box,value):
    a,b=ep(box)
    if not a<=value<=b:raise ArithmeticError('Independent value left directed enclosure')


def determinant_identity():
    h0,hm,a,b,mu,d=sy.symbols('H0 Hmu a b mu d',positive=True)
    hd=sy.Symbol('Hd')
    A=[h0*sy.exp(q/2) for q in (a,b)]
    D=[sy.exp(q/2)*(hd+hm*(sy.exp(-mu*q)-1)/mu) for q in (a,b)]
    expected=h0*hm*sy.exp((a+b)/2)*sy.exp(-mu*a)*(1-sy.exp(-mu*(b-a)))/mu
    if sy.simplify(-(A[0]*D[1]-A[1]*D[0])-expected)!=0:
        raise ArithmeticError('New translated axial divided determinant identity failed')
    x,y,z=sy.symbols('x y z')
    matrix=sy.Matrix([[1,x,x*x],[1,y,y*y],[1,z,z*z]])
    if sy.expand(matrix.det()-(y-x)*(z-x)*(z-y))!=0:
        raise ArithmeticError('New angular translated Vandermonde identity failed')
    # k=(1,-1-mu,-mu): the first two differences are negative and the third positive.
    if sy.simplify((-2-mu)*(-1-mu))!=(mu+1)*(mu+2):
        raise ArithmeticError('Angular determinant sign failed')
    return dict(axial_exact_difference_identity=True,angular_exact_Vandermonde_identity=True,
        angular_determinant_positive_for_mu_nonnegative=True,
        original_axial_determinant_negative_for_mu_positive=True)


def physical_fixture():
    """Direct paper-weight integrals with nonunit Rc, Astar and positive mu."""
    with mp.workdps(75):
        Rc,A,mu=mp.mpf('3.7'),mp.mpf('2.3'),mp.mpf('.08')
        ell=mp.log(2)/400;centers=[j*mp.log(2)/40 for j in current.CENTERS]
        beta=lambda t:mp.exp(-1/(1-t*t)) if abs(t)<1 else mp.mpf(0)
        integrate=lambda f:mp.quad(f,[-1,-.75,0,.75,1])
        Jhalf=integrate(lambda t:beta(t)*mp.exp(ell*t/2))
        def entry(j,k):
            ci=centers[j]
            def value(t):
                s=ci+ell*t;R=Rc*mp.exp(s)
                b=mp.sqrt(2/Rc)*beta(t)/(ell*mp.exp(ci/2)*Jhalf)
                return b*(2*R)**k*R*ell
            return integrate(value)
        Bz=mp.matrix([[entry(j,k) for j in current.AXIAL] for k in (0,-mu)])
        Btheta=mp.matrix([[entry(j,k) for j in current.ANGULAR] for k in (.5,-1.5-mu,-.5-mu)])
        for j in range(5):
            if abs(entry(j,-mp.mpf('.5'))-1)>mp.mpf('1e-65'):
                raise ArithmeticError('Physical dR/sqrt(2R) normalization failed')
        c=current.MPIntervalContext();c.dps=85
        numeric=current.numeric_weights(c,c.mpf('.08'))
        sr=mp.sqrt(2*Rc)
        for j in range(2):
            contains(numeric['axial_divided_rows'][0][j],Bz[0,j]/sr)
            contains(numeric['axial_divided_rows'][1][j],((2*Rc)**mu*Bz[1,j]-Bz[0,j])/(mu*sr))
        factors=[2*Rc,(2*Rc)**(-1-mu),(2*Rc)**(-mu)]
        for i in range(3):
            for j in range(3):contains(numeric['angular_rows'][i][j],Btheta[i,j]/factors[i])
        contains(numeric['positive_negative_axial_divided_determinant'],
            -mp.det(Bz)*(2*Rc)**mu/(mu*sr*sr))
        contains(numeric['positive_angular_determinant'],mp.det(Btheta)/(math.prod(factors)))
        moments=[mp.mpf(s) for s in ('.37','-.29','.41','-.23','.19')]
        rhsz=mp.matrix([-moments[0],-moments[3]/A])
        rhst=mp.matrix([-moments[1],-moments[2]/(2*A),moments[4]/A])
        alpha=mp.lu_solve(Bz,rhsz);theta=mp.lu_solve(Btheta,rhst)
        # Physical changes are the paper identities, independent of the producer transform.
        changes=[(Bz*alpha)[0],(Btheta*theta)[0],2*A*(Btheta*theta)[1],
            A*(Bz*alpha)[1],-A*(Btheta*theta)[2]]
        if max(abs(a+b) for a,b in zip(changes,moments))>mp.mpf('1e-65'):
            raise ArithmeticError('Paper five-moment signs or factor two failed')
        return dict(nonunit_Rc='3.7',nonunit_Astar='2.3',finite_mu='.08',
            physical_normalizations=5,physical_matrix_integrals=13,
            fresh_determinant_formulas=True,paper_m13_factor_two=True,paper_m15_minus_sign=True)


class Interpreter:
    """Numerical evaluation of defining graph functions, without range selection."""
    def __init__(self,g,bindings):self.g=g;self.bindings=bindings;self.cache={}
    def at(self,node,local=None):
        if local is not None and node in local:return local[node]
        if node in self.bindings:return self.bindings[node]
        if local is None and node in self.cache:return self.cache[node]
        row=self.g.nodes[node];op=row['operation'];at=lambda i:self.at(i,local)
        if op=='exact_rational':v=mp.mpf(row['numerator'])/row['denominator']
        elif op=='sum':v=sum(map(at,row['arguments']),mp.mpf(0))
        elif op=='product':v=math.prod(map(at,row['arguments']))
        elif op=='negative':v=-at(row['argument'])
        elif op=='positive_quotient':v=at(row['numerator'])/at(row['denominator'])
        elif op=='analytic_unary':
            a=at(row['argument']);name=row['name']
            v=mp.expm1(a)/a if name=='exprel' and a else (mp.mpf(1) if name=='exprel' else getattr(mp,name)(a))
        elif op=='compact_raw_beta':
            t=at(row['argument']);v=mp.exp(-1/(1-t*t)) if abs(t)<1 else mp.mpf(0)
        elif op=='definite_integral':
            variable=next(i for i,r in enumerate(self.g.nodes) if r==dict(operation='bound_variable',name=row['variable']))
            f=lambda x:self.at(row['integrand'],{**(local or {}),variable:x})
            v=mp.quad(f,[at(row['lower']),-.75,0,.75,at(row['upper'])])
        else:raise ValueError('Unbound defining node '+str(row))
        if local is None:self.cache[node]=v
        return v


def function_fixture(owner):
    """Exercise the production physical C1 transform with genuine fixture functions."""
    with mp.workdps(75):
        g=owner.graph;W=owner.functions;b=owner.background;ref=lambda n:current.source.FunctionRef(g,n)
        original=mp.mpf('.08');z=mp.mpf('.371');Rc=mp.mpf('3.7');A=mp.mpf('2.3')
        bindings={W['mu'].node:original,W['log2Rc'].node:mp.log(2*Rc),
            b['Astar_ordinary_Z0_Z3'][0]:A*(1+z*z),b['Astar_ordinary_Z0_Z3'][1]:2*A*z}
        rows={}
        definitions=[lambda x:.37+x*x,lambda x:-.29+2*x,lambda x:.41-x,
            lambda x:-.23+3*x*x,lambda x:.19+x*x*x]
        for key,fn in zip(('m11','m12','m13','m14','m15'),definitions):
            value=g.symbol('independent_fixture_'+key);deriv=g.symbol('independent_fixture_'+key+'_Z')
            rows[key]=current.source.C1Function(value,deriv);bindings[value.node]=fn(z);bindings[deriv.node]=mp.diff(fn,z)
        coefficients=current.coefficient_action(g,W,b,rows)
        interp=Interpreter(g,bindings);numerical=[interp.at(p.value.node) for p in coefficients]
        def direct(at):
            M=mp.matrix([fn(at) for fn in definitions]);a=A*(1+at*at)
            interval_interp=Interpreter(g,{W['mu'].node:original})
            Bz=mp.matrix([[interval_interp.at(i.node) for i in r] for r in W['axial_rows']])
            Bt=mp.matrix([[interval_interp.at(i.node) for i in r] for r in W['angular_rows']])
            sr=mp.sqrt(2*Rc);y0=-M[0]/sr
            yd=(M[0]-(2*Rc)**original*M[3]/a)/(original*sr)
            rhs=[-M[1]/(2*Rc),-M[2]*(2*Rc)**(1+original)/(2*a),M[4]*(2*Rc)**original/a]
            return list(mp.lu_solve(Bz,mp.matrix([y0,yd])))+list(mp.lu_solve(Bt,mp.matrix(rhs)))
        exact=direct(z)
        if max(abs(a-b) for a,b in zip(numerical,exact))>mp.mpf('1e-60'):
            raise ArithmeticError('Production physical inverse action failed')
        for j,p in enumerate(coefficients):
            if abs(interp.at(p.Z.node)-mp.diff(lambda zz:direct(zz)[j],z))>mp.mpf('1e-57'):
                raise ArithmeticError('Production C1 quotient/product derivative failed')
        try:current.coefficient_action(g,W,b,{k:mp.mpf(0) for k in rows})
        except TypeError:pass
        else:raise ArithmeticError('Coefficient action accepted interval/scalar values as functions')
        return dict(genuine_fixture_C1_action=True,amplitude_derivative_retained=True,
            caps_rejected_as_functions=True,actual_defect_functions_not_selected_by_fixture=True)


@current.precision
def run(owner,raw):
    prefix_len=raw['accepted_background_graph_prefix_length']
    if raw['exact_graph_nodes'][:prefix_len]!=owner.reports[current.quiet.band.NAME]['exact_graph_nodes']:
        raise ValueError('Immutable accepted background graph changed')
    if raw['support_certificate']!=current.exact_supports():raise ValueError('New support geometry changed')
    for name,digest in raw['input_hashes'].items():
        if sha(name)!=digest:raise ValueError('Changed checked input '+name)
    if raw[current.GATE] or any(raw[k] for k in current.OPEN):raise ValueError('Raw scope widened')
    if raw['selected_coefficient_functions_installed'] or raw['complete_pressure_C1_accuracy_installed']:
        raise ValueError('Unfinished function/pressure layer marked complete')
    c=owner.c;vals=owner.values
    if ep(vals['mu'])[0]<=0 or ep(vals['base_normalizer'])[0]<=0:
        raise ArithmeticError('Actual original positive mu/normalizer lost')
    for k in ('positive_negative_axial_divided_determinant','positive_angular_determinant'):
        if ep(vals[k])[0]<=0:raise ArithmeticError('Actual new linear matrix invertibility failed')
    if raw['actual_normalized_matrix_enclosures']!=current.record(vals):raise ValueError('Matrix archive differs')
    # Restore the same radius independently from the accepted native-frame expression.
    g=owner.graph;bg=owner.background;params=bg['parameters']
    bindings={i:sy.Symbol(k,real=True) for k,i in params.items()}
    expected=current.arithmetic(g.nodes,bg['logRc'],bindings)
    f=owner.reports[current.frame.NAME]
    native=current.arithmetic(f['exact_graph_nodes'],f['actual_Rp_native_frame_functions']['original_logRc'],bindings)
    if sy.simplify(expected-native)!=0:raise ValueError('Actual original microscopic Rc identity changed')
    identities=determinant_identity();physical=physical_fixture()
    # Fixture nodes belong only to the checker. Never archive or install them.
    graph_nodes=list(g.nodes);graph_keys=dict(g.keys)
    functions=function_fixture(owner)
    g.nodes=graph_nodes;g.keys=graph_keys
    hashes={**owner.hashes,NAME.replace('_check.json','_check.py'):sha(NAME.replace('_check.json','_check.py')),current.NAME:sha(current.NAME)}
    receipt=dict(all_passed=True,source_family=owner.family,**{current.GATE:True},
        **dict.fromkeys(current.OPEN,False),determinant_identities=identities,
        nonunit_physical_fixture=physical,C1_function_fixture=functions,
        immutable_background_prefix=True,actual_Rc_identity=True,original_amplitude_positive_on_all_four_Z_cells=True,
        new_supports_strictly_disjoint=True,actual_mu_positive_and_retained=True,
        source_functions_and_bounds_distinct=True,full_defect_C1_adapter_still_open=True,input_hashes=hashes)
    (HERE/NAME).write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8',newline='\n')
    print('Original n1 normalized physical linear repair operator checks passed',flush=True)
    return receipt


if __name__=='__main__':
    owner,raw=current.run();run(owner,raw)
