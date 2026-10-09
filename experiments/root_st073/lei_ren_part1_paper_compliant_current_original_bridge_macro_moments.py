"""Actual frozen-macro six moments from the same original coupled F/V.

Finite signed exponential polynomials are integrated before enclosure.
The complete exponential has a separate nonzero remainder. Native actual
micro histories and positive incoming decay remain, with no moment reset.
"""
import json
import math
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_bridge_macro_functions as fields

HERE,PREFIX,sha=fields.HERE,fields.PREFIX,fields.sha
base,prior,ep=fields.base,fields.prior,fields.ep
NAME=PREFIX+'current_original_bridge_macro_moments.json'
RECEIPT=PREFIX+'current_original_bridge_macro_moments_check.json'
GATE='original_frozen_macro_actual_six_moment_function_integrals_installed'
RATES=dict(H=2,M=1,K=2,A=1,B=2,C=1)


class ExponentialPolynomial:
    """sum_(a,k) coefficient_Z_jet * t^k * exp(a*t), integer a,k."""
    def __init__(self,flow,terms=None):
        self.flow=flow;self.terms={}
        for key,row in (terms or {}).items():
            if not isinstance(key,tuple) or len(key)!=2 or any(type(x) is not int for x in key) or key[1]<0:
                raise ValueError('Integer exponential degree and nonnegative polynomial degree required')
            if len(row)!=6 or any(v.scale.bases is not flow.logs or v.ledger is not flow.ledger for v in row):
                raise ValueError('Same source basis and ordinary Z0..5 coefficient rows required')
            if any(not v.zero for v in row):self.terms[key]=row

    @classmethod
    def constant(cls,flow,row):return cls(flow,{(0,0):row})

    def pair(self,other):
        if not isinstance(other,ExponentialPolynomial) or other.flow is not self.flow:
            raise ValueError('One defining flow owner required')

    def __add__(self,other):
        self.pair(other);terms=dict(self.terms)
        for key,row in other.terms.items():terms[key]=self.flow.add(terms[key],row) if key in terms else row
        return type(self)(self.flow,terms)

    def __neg__(self):return self.scale(-1)
    def __sub__(self,other):return self+-other

    def scale(self,value):return type(self)(self.flow,{key:self.flow.scale(row,value) for key,row in self.terms.items()})

    def __mul__(self,other):
        self.pair(other);terms={}
        for (a,k),row in self.terms.items():
            for (b,l),second in other.terms.items():
                key=(a+b,k+l);product=self.flow.multiply(row,second)
                terms[key]=self.flow.add(terms[key],product) if key in terms else product
        return type(self)(self.flow,terms)

    def primitive(self):
        """Exact source primitive with value zero at t=0, all resonances."""
        f=self.flow;terms={}
        def add(key,row):terms[key]=f.add(terms[key],row) if key in terms else row
        for (a,k),row in self.terms.items():
            if a==0:add((0,k+1),f.scale(row,f.c.mpf(1)/(k+1)))
            else:
                for j in range(k+1):
                    coefficient=f.c.mpf((-1)**j*math.factorial(k)//math.factorial(k-j))/a**(j+1)
                    add((a,k-j),f.scale(row,coefficient))
                add((0,0),f.scale(row,-f.c.mpf((-1)**k*math.factorial(k))/a**(k+1)))
        return type(self)(f,terms)

    def factor_at(self,a,R0,R1):
        # Collect source radius factors before any logarithmic evaluation.
        return self.flow.radial_power(R1,a)*self.flow.radial_power(R0,-a)

    def evaluate(self,S,R0,R1):
        f=self.flow
        return f.add(*[f.scale(row,self.factor_at(a,R0,R1)*S**k) for (a,k),row in self.terms.items()])

    def mass(self,rate,a,k,S,R0,R1,empty=False):
        f=self.flow;c=f.c
        if type(rate) is not int or rate not in (1,2):raise ValueError('Original history rate1 or2 required')
        if empty:return f.scalar(0)
        b=a+rate;incoming=self.factor_at(-rate,R0,R1)
        if b==0:return incoming*(S**(k+1)/(k+1))
        endpoint=sum((c.mpf((-1)**j*math.factorial(k)//math.factorial(k-j))*S**(k-j)/b**(j+1)
                      for j in range(k+1)),c.mpf(0))
        return self.factor_at(a,R0,R1)*endpoint-incoming*(c.mpf((-1)**k*math.factorial(k))/b**(k+1))

    def volterra(self,rate,S,R0,R1,empty=False):
        return self.flow.add(*[self.flow.scale(row,self.mass(rate,a,k,S,R0,R1,empty))
                               for (a,k),row in self.terms.items()])


class ActualMacroMoments:
    def __init__(self,flow,inlet):
        if not flow.sources_set or set(inlet)!=set(RATES):raise ValueError('Defining F/V and all six actual inlet histories required')
        self.flow=flow;self.inlet={}
        for key,row in inlet.items():
            if len(row)!=6 or any(v.scale.bases is not flow.logs or v.ledger is not flow.ledger for v in row):
                raise ValueError('Actual inlet rows must share the source basis/ledger')
            self.inlet[key]=row

    def evaluate(self,fraction):
        f=self.flow;c=f.c;S,R0,R1,empty=f.geometry(fraction)
        endpoint=f.evaluate(fraction)
        B=endpoint['full_prefix_log_jet_norm'];upper=f.small_upper(B)
        if ep(upper)[1]>.5:raise ValueError('Full weighted exponential prefix must be <=1/2')
        exp_error=B*B*B*(c.exp(upper)/6)
        const=lambda row:ExponentialPolynomial.constant(f,row)
        one=const([f.scalar(1)]+[f.scalar(0)]*5)
        ell=const(f.ell_in)
        for l,row in enumerate(f.d):
            primitive=ExponentialPolynomial(f,{(1-l,0):f.scale(row,R0)}).primitive()
            ell=ell-primitive.scale(f.h*c.mpf('.5'))
        exp_poly=one+ell+(ell*ell).scale(c.mpf('.5'))
        phi=const(f.phi0)*exp_poly
        Vin=f.add(f.V0,f.deltaV_in);V=const(Vin)
        source_mass=f.scalar(0)
        for part,p in (('hydro',1),('pressure',1),('swirl',2)):
            scale=f.factor((1,int(part=='pressure'),int(part=='swirl'),0,0))
            for j,row in enumerate(f.drive[part]):
                product=f.multiply(f.q,row)
                drive=ExponentialPolynomial(f,{(p-j,0):f.scale(product,f.radial_power(R0,p))})
                V=V-(drive*exp_poly).primitive().scale(scale)
                I=f.I(p,j,S,R0,R1,empty)
                mass=prior.ScaledEnclosure(I.scale,fields.magnitude(c,I.coefficient),f.ledger)
                source_mass+=f.norm(product)*mass*scale
        exp_poly_norm=f.scalar(1)+B+B*B*c.mpf('.5')
        phi_norm=f.norm(f.phi0)*exp_poly_norm
        V_norm=f.norm(Vin)+source_mass*exp_poly_norm
        phi_error=f.norm(f.phi0)*exp_error;V_error=source_mass*exp_error
        source_polys=dict(H=phi.scale(2),M=V,K=(phi*V).scale(2),A=V*V,B=phi*phi,C=phi*phi)
        source_errors=dict(H=phi_error*2,M=V_error,
            K=(phi_norm*V_error+V_norm*phi_error+phi_error*V_error)*2,
            A=V_norm*V_error*2+V_error*V_error,
            B=phi_norm*phi_error*2+phi_error*phi_error,
            C=phi_norm*phi_error*2+phi_error*phi_error)
        field_phi=f.add(phi.evaluate(S,R0,R1),f.error_rows(phi_error))
        field_V=f.add(V.evaluate(S,R0,R1),f.error_rows(V_error))
        source_end=dict(H=f.scale(field_phi,2),M=field_V,K=f.scale(f.multiply(field_phi,field_V),2),
            A=f.multiply(field_V,field_V),B=f.multiply(field_phi,field_phi),C=f.multiply(field_phi,field_phi))
        histories={};derivatives={};evidence={}
        for name,rate in RATES.items():
            poly=source_polys[name];decay=poly.factor_at(-rate,R0,R1)
            signed=poly.volterra(rate,S,R0,R1,empty)
            # The positive original Volterra mass is a bound on the complete
            # source remainder, with no short-body or infinite-tail reset.
            W=f.scalar(0) if empty else (f.scalar(1)-decay)*(c.mpf(1)/rate)
            mass=prior.ScaledEnclosure(W.scale,fields.magnitude(c,W.coefficient),f.ledger)
            integral_error=source_errors[name]*mass
            inherited=f.scale(self.inlet[name],decay)
            actual=f.add(inherited,signed,f.error_rows(integral_error))
            if empty:actual=self.inlet[name]
            histories[name]=actual
            derivatives[name]=f.add(source_end[name],f.scale(actual,-rate))
            evidence[name]=dict(rate=rate,source={'H':'2phi','M':'V','K':'2phiV','A':'V^2','B':'phi^2','C':'phi^2'}[name],
                original_incoming_decay=decay,original_positive_kernel_mass=W,
                actual_inlet_rows=self.inlet[name],inherited_rows=inherited,
                signed_complete_polynomial_integral=signed,complete_source_error_norm=source_errors[name],
                complete_integral_error_norm=integral_error,
                polynomial_term_count=len(poly.terms),max_polynomial_degree=max((k for a,k in poly.terms),default=0),
                source_products_formed_before_enclosure=True,original_finite_interval_not_shortened=True)
        return dict(actual_six_moment_functions=histories,actual_radial_ODE_derivative_functions=derivatives,
            same_original_macro_field_functions=dict(phi=field_phi,V=field_V),
            complete_integral_evidence=evidence,
            complete_exponential_error_norm=exp_error,complete_field_error_norms=dict(phi=phi_error,V=V_error),
            geometry=dict(fraction=list(fraction),S=S,R0=R0,R1=R1),
            original_volterra_ODEs_preserved=True,comparison_moments_not_substituted=True,
            analytic_P0_remains_separate=True,ordinary_Z_orders=list(range(6)),
            kernel_coordinate='t=log(R/R0); dR/R=dt; original positive kernel applied once',
            source_basis_and_micro_inlet_errors_retained=True)


class OriginalBridgeMacroMoments:
    mode='genuine_actual_frozen_macro_six_history_functions_with_micro_inlet_errors'
    def __init__(self,dps=500):
        self.fields=fields.OriginalBridgeMacroFunctions(dps);self.c=self.fields.c
        self.family=self.fields.family;self.hashes=dict(self.fields.hashes);self.owners={}
        for name in (fields.NAME,fields.RECEIPT):
            row=json.loads((HERE/name).read_bytes())
            if not row.get(fields.GATE) or name==fields.RECEIPT and not row.get('all_passed'):
                raise ValueError('Accepted complete original macro field functions required')
            for path,digest in row['input_hashes'].items():fields.previous.bind(self.hashes,path,digest)
            fields.previous.bind(self.hashes,name,sha(name))
        fields.previous.bind(self.hashes,Path(__file__).name,sha(Path(__file__).name))

    def owner(self,label):
        if label in self.owners:return self.owners[label]
        flow,proof=self.fields.owner(label)
        packet=self.fields.records['actual_bridge_integrals']['packets'][label]['second_exit']
        raw=packet['actual_own_six_moments_axial5']
        inlet={name:[flow.scalar(fields.previous.read_interval(self.c,v)) for v in rows] for name,rows in raw.items()}
        owner=ActualMacroMoments(flow,inlet)
        self.owners[label]=(owner,proof);return owner,proof

    def evaluate(self,label,fraction):
        owner,proof=self.owner(label)
        with mp.workdps(self.c.dps+40):result=owner.evaluate(fraction)
        return dict(mode=self.mode,source_frame=label,evaluation=fields.serialized(result),
            source_family=self.family,implicit_source_sha256=self.fields.source,datum_enclosure_sha256=self.fields.datum,
            source_function_proof_reference='same OriginalBridgeMacroFunctions source owner, phi/V and actual second_exit histories',
            actual_micro_function_provider_installed=False,whole_axis_functions_installed=False,
            actual_R100_R110_switch_function_installed=False)


def run():
    began=time.monotonic();owner=OriginalBridgeMacroMoments()
    packets={label:[owner.evaluate(label,q) for q in ((0,1),(1,2),(1,1))] for label in ('0','.5')}
    result=dict(**{GATE:True},source_family=owner.family,packets=packets,
        **dict.fromkeys(fields.previous.OPEN,False),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
        scope='Actual frozen macro H/M/K/A/B/C defining finite Volterra integrals and radial ODE rows at native0,.5. Complete nonlinear F/V and actual micro inlet errors retained. Micro functions, whole Z and switches remain open.')
    (HERE/NAME).write_text(json.dumps(base.encoded(result),indent=2)+'\n',encoding='utf8')
    print('Actual original frozen macro six moment functions integrated',flush=True);return result


if __name__=='__main__':run()
