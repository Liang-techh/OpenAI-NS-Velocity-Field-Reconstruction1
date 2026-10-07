"""Own five-history transport, pressure and radial/full stress recovery.

One common velocity unit S is kept formal. Radius-normalized histories
have a stable Duhamel transport and do not reset outside modulation. The
operator accepts full source jet covers; it does not manufacture covers
for the generic loop or infer source functions from cached enclosures.
Current original O2 cache conversion is executable without constructors.
"""
import ast
import gzip
import json
import math
from pathlib import Path
from types import SimpleNamespace
import mpmath as mp
import sympy as s
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_generic_shear_loop as loop
import lei_ren_part1_paper_compliant_current_patch_physical_numeric_bounds as numeric
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_compliant_collar_stress_C3 import axial_derivative
from lei_ren_part1_paper_compliant_long_reshape_mixed_C4 import exponential_derivatives
from lei_ren_part1_paper_compliant_current_O3_independent_repair_operator import exp_average

HERE,PREFIX,sha=loop.HERE,loop.PREFIX,loop.sha
NAME=PREFIX+'current_generic_shear_moment_recovery.json'
RECEIPT=PREFIX+'current_generic_shear_moment_recovery_check.json'
VIEWS=PREFIX+'current_generic_shear_moment_recovery_views.json.gz'
GATE='generic_shear_own_moment_transport_and_recovery_implemented'
CACHE_GATE='current_original_O2_common_velocity_unit_recovery_attached'
OPEN=loop.OPEN
RATES=dict(m=1,h='1.5',k='1.5',e=1,p=0)
UNITS=dict(m='Mz/(R*S)',h='Mtheta/(sqrt(2)*R^1.5*S)',
           k='Mtheta_z/(sqrt(2)*R^1.5*S^2)',e='Mztheta/(R*S^2)',p='Mp/S^2')
endpoints=numeric.transport.endpoints
read=numeric.transport.read_interval


def encode(value):
    if isinstance(value,dict):return {k:encode(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)):return [encode(v) for v in value]
    if isinstance(value,IntervalTaylor):return dict(coefficients=[encode(v) for v in value.coefficients])
    if hasattr(value,'_mpi_'):
        lo,hi=endpoints(value)
        return dict(lower=mp.nstr(lo,85),upper=mp.nstr(hi,85),
                    lower_exact_mpf_tuple=list(lo._mpf_),upper_exact_mpf_tuple=list(hi._mpf_))
    return value


def check_jets(c,values):
    if not values or any(type(v) is not IntervalTaylor or v.ctx is not c for v in values):
        raise ValueError('Source jets must use the same interval Taylor context')
    if min(v.order for v in values)<1:
        raise ValueError('At least one actual axial derivative required')
    if any(not mp.isfinite(x) for v in values for q in v.coefficients for x in endpoints(q)):
        raise ValueError('Finite source jet coefficients required')


def history_densities(E,V):
    return dict(m=V,h=E,k=E*V,e=V*V-E*E/2,p=E*E/2)


def increment_densities(E,V,delta_E,delta_V):
    """Full signed source differences, including nonzero original V."""
    e,u=delta_E,delta_V
    return dict(m=u,h=e,k=V*e+E*u+e*u,
                e=2*V*u+u*u-E*e-e*e/2,p=E*e+e*e/2)


def zero_histories(jet):
    return {key:jet*0 for key in RATES}


def shifted_rows(c,rows,rate):
    return [sum((rows[k]*math.comb(j,k)*c.mpf(rate)**(j-k) for k in range(j+1)),rows[0]*0)
            for j in range(len(rows))]


def ordinary_product_row(a,b,j):
    return sum((a[k]*b[j-k]*math.comb(j,k) for k in range(j+1)),a[0]*0)


class GenericMomentRecovery:
    def __init__(self,c,*,source_family,P0,original,defect=None):
        if set(source_family)!=set(('actual_five_defect_family_sha256','implicit_source_sha256','datum_enclosure_sha256')):
            raise ValueError('Explicit current family/source/axis datum required')
        if set(original)!=set(RATES) or (defect is not None and set(defect)!=set(RATES)):
            raise ValueError('All five inlet histories required')
        defect=zero_histories(P0) if defect is None else defect
        check_jets(c,[P0,*original.values(),*defect.values()])
        self.ctx=c;self.family=dict(source_family);self.P0=P0
        self.original=dict(original);self.defect=dict(defect)

    def own(self):
        return {key:self.original[key]+self.defect[key] for key in RATES}

    def advance(self,*,width,E,V,delta_E,delta_V,source_family):
        """Enclose a continuous Duhamel step using entire-cell source covers.

        E,V and increments must cover the same actual functions and axial
        derivatives throughout this logR cell. They are not samples. The
        width is a Z-independent exact/source-enclosed log-radius length.
        No phase sampling, frequency clipping or pressure resetting occurs.
        """
        c=self.ctx;w=c.mpf(width);wl,wh=endpoints(w)
        if source_family!=self.family:raise ValueError('Changed moment-source family or pressure datum')
        if wl<0 or not all(mp.isfinite(v) for v in (wl,wh)):
            raise ValueError('Finite nonnegative log-radius cell width required')
        check_jets(c,[E,V,delta_E,delta_V])
        old=history_densities(E,V);inc=increment_densities(E,V,delta_E,delta_V)
        original={};defect={}
        for key,r in RATES.items():
            rate=c.mpf(r);decay=c.exp(-rate*w)
            # w*integral_0^1 exp(-rate*w*t)dt, including w=0 exactly.
            # exp_average keeps microscopic positive widths; 1-exp is
            # never formed by subtraction from a rounded unit.
            mass=w if r==0 else w*exp_average(c,-rate*w)
            original[key]=self.original[key]*decay+old[key]*mass
            defect[key]=self.defect[key]*decay+inc[key]*mass
        return GenericMomentRecovery(c,source_family=self.family,P0=self.P0,
                                     original=original,defect=defect)

    def field(self,*,Z,delta,E,V,E_y,V_y):
        """Recover from own histories; E,V and y derivatives are changed profiles."""
        c=self.ctx;check_jets(c,[Z,E,V,E_y,V_y]);de=c.mpf(delta)
        if endpoints(de)[0]<0 or endpoints(de)[1]>=1:
            raise ValueError('Similarity parameter must satisfy 0<=delta<1')
        zl,zh=endpoints(Z[0])
        if zl < -1 or zh > 1:raise ValueError('Whole profile axial domain [-1,1] required')
        M=self.own();m,h,k,e,p=(M[key] for key in ('m','h','k','e','p'))
        L=1-Z*Z*de;d=1-Z*Z
        pressure=self.P0+p
        transport=m*(Z*(1-de))+axial_derivative(m)*d
        Q=(V*(2*Z)-transport)/L
        theta_linear=(-E+h*(1-de/2)-axial_derivative(h)*(Z*((1-de)/2)))/L
        theta_quadratic=(k*(Z*(2*de-1))-axial_derivative(k)*d+E*transport)/L
        axial_linear=(-V+(m-axial_derivative(m)*Z)*((1-de)/2))/L
        axial_quadratic=(V*transport+e*(Z*(2*de))-axial_derivative(e)*d
                         +pressure*(Z*(2*(1+de)))-axial_derivative(pressure)*d)/L
        return dict(source_family=self.family,own_normalized_five_histories=M,
                    original_histories=self.original,transported_defects=self.defect,
                    original_axis_pressure_over_S_squared=self.P0,
                    absolute_pressure_over_S_squared=pressure,
                    Utheta_over_S=E,Uz_over_S=V,Ur_over_S_sqrt_R_over_2=Q,
                    inertial_theta_linear=theta_linear,inertial_theta_quadratic=theta_quadratic,
                    inertial_axial_linear=axial_linear,inertial_axial_quadratic=axial_quadratic,
                    shear_theta=2*E_y-E,shear_axial=2*V_y,
                    formal_stress_modes=dict(inertial_linear='sqrt(R/2)*S',
                        inertial_quadratic='sqrt(R/2)*S^2',shear='S/sqrt(2R)'),
                    pressure_radial_balance='d_y(P/S^2)=E^2/2',
                    radial_incompressibility='L*d_y(R*S*Q)=R*S*((1+delta)*Z*V-d*V_Z+2*Z*V_y)',
                    cell_ranges_are_covers_not_defining_values=True,
                    **dict.fromkeys(OPEN,False))

    def field_rows(self,*,Z,delta,E_rows,V_rows):
        """Own history/velocity mixed4 and signed stress mixed3 source rows.

        Inputs are actual ordinary logR derivatives, with axial Taylor jets.
        The +1/2 radial/inertial and -1/2 shear prefactors are differentiated
        exactly once. Cell/source covers still do not define a point field.
        """
        if len(E_rows)<5 or len(V_rows)<5:
            raise ValueError('Actual ordinary logR profile rows zero through four required')
        c=self.ctx;check_jets(c,[Z,*E_rows[:5],*V_rows[:5]])
        result=self.field(Z=Z,delta=delta,E=E_rows[0],V=V_rows[0],E_y=E_rows[1],V_y=V_rows[1])
        de=c.mpf(delta);L=1-Z*Z*de;d=1-Z*Z
        M={key:[value] for key,value in self.own().items()}
        for j in range(4):
            M['m'].append(V_rows[j]-M['m'][j])
            M['h'].append(E_rows[j]-M['h'][j]*c.mpf('1.5'))
            M['k'].append(ordinary_product_row(E_rows,V_rows,j)-M['k'][j]*c.mpf('1.5'))
            M['e'].append(ordinary_product_row(V_rows,V_rows,j)-ordinary_product_row(E_rows,E_rows,j)/2-M['e'][j])
            M['p'].append(ordinary_product_row(E_rows,E_rows,j)/2)
        m,h,k,e,p=(M[key] for key in ('m','h','k','e','p'))
        transport=[m[j]*(Z*(1-de))+axial_derivative(m[j])*d for j in range(5)]
        Q=[(V_rows[j]*(2*Z)-transport[j])/L for j in range(5)]
        P=[self.P0+p[0]]+p[1:]
        theta_linear=[(-E_rows[j]+h[j]*(1-de/2)-axial_derivative(h[j])*(Z*((1-de)/2)))/L for j in range(4)]
        theta_quadratic=[(k[j]*(Z*(2*de-1))-axial_derivative(k[j])*d+ordinary_product_row(E_rows,transport,j))/L for j in range(4)]
        axial_linear=[(-V_rows[j]+(m[j]-axial_derivative(m[j])*Z)*((1-de)/2))/L for j in range(4)]
        axial_quadratic=[(ordinary_product_row(V_rows,transport,j)+e[j]*(Z*(2*de))-axial_derivative(e[j])*d
                           +P[j]*(Z*(2*(1+de)))-axial_derivative(P[j])*d)/L for j in range(4)]
        shear_theta=[2*E_rows[j+1]-E_rows[j] for j in range(4)]
        shear_axial=[2*V_rows[j+1] for j in range(4)]
        result.update(own_normalized_history_ordinary_y_rows=M,
            physical_velocity_pressure_ordinary_y_rows=dict(theta=E_rows[:5],axial=V_rows[:5],
                radial=shifted_rows(c,Q,'.5'),pressure=P),
            physical_primitive_ordinary_y_rows={key:shifted_rows(c,rows,rate) for (key,rows),rate in zip(M.items(),RATES.values())},
            full_signed_stress_ordinary_y_rows=dict(
                inertial_theta_linear=shifted_rows(c,theta_linear,'.5'),
                inertial_theta_quadratic=shifted_rows(c,theta_quadratic,'.5'),
                inertial_axial_linear=shifted_rows(c,axial_linear,'.5'),
                inertial_axial_quadratic=shifted_rows(c,axial_quadratic,'.5'),
                shear_theta=shifted_rows(c,shear_theta,'-.5'),shear_axial=shifted_rows(c,shear_axial,'-.5')),
            actual_profile_y_rows_supplied_not_chart_derivatives=True,
            original_axis_pressure_only_in_absolute_row_zero=True,
            physical_radial_and_stress_prefactor_shifts_applied_once=True)
        return result


def exact_theorem():
    checks={};asts=numeric.transport.SourceAST()
    def zero(name,a,b):
        if s.cancel(s.expand(a-b))!=0:raise ArithmeticError('Generic moment recovery identity: '+name)
        checks[name]=True
    y,z,delta,S=s.symbols('y Z delta S',real=True)
    E,V,u,v=s.symbols('E V delta_E delta_V',real=True)
    m,h,k,e,p=s.symbols('m h k e p',real=True)
    invS=s.Symbol('inverse_S')
    original=asts.method('pre_pulse_mixed_C4','physical_mixed')
    wanted=("rows['m'].append(V[j]-rows['m'][j])",
            "rows['h'].append(U[j]-rows['h'][j]*c.mpf('1.5'))",
            "rows['k'].append(product_rows(U,V,j)-rows['k'][j]*c.mpf('1.5'))",
            "rows['e'].append(product_rows(V,V,j)*invP2-product_rows(U,U,j)/2-rows['e'][j])",
            "rows['p'].append(product_rows(U,U,j)/2)")
    calls={ast.dump(node) for node in ast.walk(original) if isinstance(node,ast.Call)}
    for text in wanted:
        if ast.dump(ast.parse(text,mode='eval').body) not in calls:
            raise ValueError('Original five-history source program changed')
    checks['current_source_five_ODEs_bound']=True
    old=dict(m=V,h=E,k=E*V,e=V*V-E*E/2,p=E*E/2)
    new=dict(m=V+v,h=E+u,k=(E+u)*(V+v),e=(V+v)**2-(E+u)**2/2,p=(E+u)**2/2)
    inc=dict(m=v,h=u,k=V*u+E*v+u*v,e=2*V*v+v*v-E*u-u*u/2,p=E*u+u*u/2)
    for key in RATES:zero('own_full_nonzero_source_density_'+key,new[key]-old[key],inc[key])
    for key,rate in RATES.items():
        r=s.Rational(str(rate));D=s.Function('delta_'+key)(y)
        A=s.Function('integrated_'+key)(y);f=s.Function('source_'+key)(y)
        zero('Duhamel_FTC_'+key,s.diff(s.exp(-r*y)*A,y).subs(s.diff(A,y),s.exp(r*y)*f),f-r*s.exp(-r*y)*A)
    # Ordinary current raw pre source uses physical Uz/m/k. Express it in
    # common S units, and replay its full signed inertial/shear operator.
    E,V,m,h,k,e,P=(s.Function(name)(z) for name in ('E','V','m','h','k','e','absolute_P'))
    Ey,Vy=s.symbols('ordinary_E_y ordinary_V_y',real=True)
    ctx=SimpleNamespace(mpf=lambda x:s.Rational(str(x)))
    def products(a,b,count=None):
        if count is None:return [products(a,b,j) for j in range(min(len(a),len(b)))]
        return sum(a[j]*b[count-j]*s.binomial(count,j) for j in range(count+1))
    raw=asts.replay('current_pre_pulse_stress_operator','raw_pre_stress_rows',dict(
        axial_derivative=lambda value:s.diff(value,z),product_rows=products,
        shifted_rows=lambda rows,rate:rows))
    zeros=lambda value:[value,s.Integer(0),s.Integer(0),s.Integer(0),s.Integer(0)]
    erows=[E,Ey,0,0,0];vrows=[S*V,S*Vy,0,0,0]
    parts=raw(ctx,delta,z,erows,vrows,dict(m=zeros(S*m),h=zeros(h),k=zeros(S*k),e=zeros(e)),zeros(P))
    L=1-delta*z*z;d=1-z*z;transport=(1-delta)*z*m+d*s.diff(m,z)
    Itl=(-E+(1-delta/2)*h-z*(1-delta)/2*s.diff(h,z))/L
    Itq=((2*delta-1)*z*k-d*s.diff(k,z)+E*transport)/L
    Izl=(-V+(1-delta)/2*(m-z*s.diff(m,z)))/L
    Izq=(V*transport+2*delta*z*e-d*s.diff(e,z)+2*(1+delta)*z*P-d*s.diff(P,z))/L
    theta=sum(row['shape'][0] for name,row in parts['theta'].items() if name!='variable_radial_shear')
    axial=sum(row['shape'][0]*(S*S if row['mode'][1]==2 else 1) for name,row in parts['axial'].items() if name!='axial_radial_shear')
    zero('same_current_full_inertial_theta_units',S*theta,S*Itl+S*S*Itq)
    zero('same_current_full_inertial_axial_units',axial,S*Izl+S*S*Izq)
    zero('same_current_full_angular_shear',parts['theta']['variable_radial_shear']['shape'][0],2*Ey-E)
    zero('same_current_full_axial_shear',parts['axial']['axial_radial_shear']['shape'][0],2*S*Vy)
    V=s.Function('V')(y,z);m=s.Function('m')(y,z)
    Q=(2*z*V-(1-delta)*z*m-d*s.diff(m,z))/L
    zero('own_moment_radial_divergence_exact',
         (L*(Q+s.diff(Q,y))).subs({s.diff(m,y,z):s.diff(V,z)-s.diff(m,z),
                                  s.diff(m,y):V-m},simultaneous=True).doit(),
         (1+delta)*z*V-d*s.diff(V,z)+2*z*s.diff(V,y))
    return dict(passed=True,identities=checks,common_units=UNITS,
                current_primitive_conversion='m=raw_m/S, h=raw_h, k=raw_k/S, e=raw_e, p=raw_p',
                pressure='same P0/S^2 + own_p; never exterior-renormalize before repair',
                full_stress_not_local_increment=True,input_hashes=asts.hashes)


class CurrentO2RecoveryCache:
    """Convert saved entire-chart covers, never arbitrary point evaluation."""
    def __init__(self):
        self.data=numeric.inputs();self.ctx=self.data['ctx'];self.hashes=dict(self.data['hashes'])
        name=PREFIX+'current_O2_background_tensor.json.gz';check=PREFIX+'current_O2_background_tensor_check.json'
        receipt=json.loads((HERE/check).read_bytes())
        if not receipt['all_passed'] or not receipt['current_actual_four_O2_completed_tensor_joins_certified']:
            raise ValueError('Checked current O2 source/tensor joins required')
        for file,digest in receipt['input_hashes'].items():
            if sha(file)!=digest:raise ValueError('Changed current source: '+file)
        self.family=self.data['accepted']['source_family']
        if any(receipt.get(key)!=value for key,value in self.family.items()):
            raise ValueError('Current O2 source/pressure family differs')
        self.hashes.update(receipt['input_hashes']);self.hashes[check]=sha(check);self.hashes[name]=sha(name)
        self.records=json.loads(gzip.decompress((HERE/name).read_bytes()))['current_actual_O2_tensor_views']
        self.invS=self.ctx.exp(-self.data['logP'])

    def jet(self,value):
        coeffs=value['coefficients'] if isinstance(value,dict) else value
        return IntervalTaylor(self.ctx,[read(self.ctx,v) for v in coeffs])

    def recover(self,chart):
        if chart not in ('Rh_reference','O2_slope','O2_axial','O2_buffer'):
            raise ValueError('Saved original O2 whole chart required')
        c=self.ctx;pre=self.records[chart+'_whole']['actual_upstream_original_pre_O2_source']
        if not pre['actual_five_histories_and_analytic_pressure_retained']:
            raise ValueError('Actual original histories/P0 required')
        expected=dict(Utheta='Pstar',Uz='1',Ur='sqrt(current R/2)',pressure='Pstar^2',
            Mz='current R',Mtheta='sqrt2*current R^1.5*Pstar',
            Mtheta_z='sqrt2*current R^1.5*Pstar',Mztheta='current R*Pstar^2',Mp='Pstar^2')
        if pre['exact_formal_prefactors']!=expected:raise ValueError('Current O2 raw source units changed')
        velocities=pre['physical_velocity_pressure_y_derivative_Taylor']
        E0=self.jet(pre['Utheta_over_Pstar_axial5_coefficients'])
        logrows=[self.jet(value) for value in pre['log_Utheta_ordinary_y_derivatives']]
        Erows=[E0*row for row in exponential_derivatives(logrows)]
        Vrows=[self.jet(value)*self.invS for value in pre['Uz_ordinary_y_derivative_axial5']]
        original={key:self.jet(pre['actual_normalized_primitive_y_derivative_axial5'][key][0])*
                  (self.invS if key in ('m','k') else 1) for key in RATES}
        P0=self.jet(pre['original_P0_axial5_coefficients'])
        Z=IntervalTaylor.variable(c,read(c,pre['Z']),5)
        recovery=GenericMomentRecovery(c,source_family=self.family,P0=P0,original=original)
        result=recovery.field_rows(Z=Z,delta=self.data['delta'],E_rows=Erows,V_rows=Vrows)
        result.update(chart=chart,coverage_coordinate=pre['coverage_coordinate'],
            covers_saved_original_functions_only=True,arbitrary_point_field_not_evaluated=True,
            modified_generic_loop_not_substituted_for_original=True)
        return result


def run():
    theorem=exact_theorem();cache=CurrentO2RecoveryCache()
    parent=json.loads((HERE/loop.RECEIPT).read_bytes())
    if not parent['all_passed']:raise ValueError('Checked generic loop kernel required')
    hashes={**cache.hashes,**theorem['input_hashes'],**parent['input_hashes'],
            loop.RECEIPT:sha(loop.RECEIPT),Path(__file__).name:sha(Path(__file__).name)}
    for name,digest in hashes.items():
        if sha(name)!=digest:raise ValueError('Generic recovery dependency changed: '+name)
    covers={key:cache.recover(key) for key in ('Rh_reference','O2_slope','O2_axial','O2_buffer')}
    (HERE/VIEWS).write_bytes(gzip.compress((json.dumps(encode(covers),indent=2)+'\n').encode(),mtime=0))
    hashes[VIEWS]=sha(VIEWS)
    result=dict(source_family=cache.family,exact_own_five_history_and_full_recovery_theorem=theorem,
        current_original_O2_common_unit_cover_views=VIEWS,current_original_O2_saved_domains=list(covers),
        **{GATE:True,CACHE_GATE:True},**dict.fromkeys(OPEN,False),
        scope='Whole-source cover/own-history transport operator and full pressure/radial/signed stress recovery; actual saved original O2 common-unit conversion. Generic loop current full packets and changed whole-family cumulative values remain uninstalled.',
        input_hashes=hashes)
    (HERE/NAME).write_text(json.dumps(encode(result),indent=2)+'\n',encoding='utf8')
    print('Generic own-history transport and full recovery; current original O2 cache attached',flush=True)
    return result


if __name__=='__main__':run()
