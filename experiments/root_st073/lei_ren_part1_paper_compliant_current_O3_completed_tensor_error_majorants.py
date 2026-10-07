"""Variable-N errors of the original completed tensor and leading remainder.

Replay the actual signed source programs, subtract the unchanged original
source, and bound every differential polynomial. Original nonzero radial
velocity and moments remain in cross terms. No graph constructors or saved
fixed-N modified values are used. These bounds do not admit a signed cone.
"""
import ast
import copy
import functools
import gzip
import json
import math
from pathlib import Path
from types import SimpleNamespace
import mpmath as mp
import sympy as s
from sympy.core.function import AppliedUndef
import lei_ren_part1_paper_compliant_current_O3_recovered_error_majorants as recovered
import lei_ren_part1_paper_compliant_pulse_end_physical_C2 as physical
from lei_ren_part1_paper_compliant_pulse_physical_bounds import (
    physical_operators,interval_expression,ZSYM,DSYM,BSYM)
from lei_ren_part1_paper_interval_taylor import IntervalTaylor

HERE,PREFIX,sha=recovered.HERE,recovered.PREFIX,recovered.sha
NAME=PREFIX+'current_O3_completed_tensor_error_majorants.json'
RECEIPT=PREFIX+'current_O3_completed_tensor_error_majorants_check.json'
VIEWS_NAME=PREFIX+'current_O3_completed_tensor_error_majorants_views.json.gz'
BASE_RECEIPT=PREFIX+'current_modified_pre_stress_check.json'
BASE_VIEWS=PREFIX+'current_modified_pre_stress_views.json.gz'
parameters,endpoints,upper=recovered.parameters,recovered.endpoints,recovered.upper
OPEN=tuple(k for k in recovered.OPEN if k!='completed_signed_tensor_error_bounds_available')+(
    'whole_signed_tensor_cones_certified','global_flat_remainder_certified')
Y=s.Symbol('logR',real=True);Z=ZSYM;D=DSYM;EPS=s.Symbol('formal_Ad',positive=True)
FUNCTIONS={name:s.Function(name)(Y,Z) for name in (
    'U0','M0','R0','H0','K0','E0','P0','du','V','dm','dh','dk','de','dp','dr')}


class Formula:
    """Exact source function; no interval sampling or rational simplification."""
    def __init__(self,value):self.expr=s.sympify(value)
    @staticmethod
    def val(v):return v.expr if isinstance(v,Formula) else s.sympify(v)
    def __add__(self,v):return Formula(self.expr+self.val(v))
    __radd__=__add__
    def __sub__(self,v):return Formula(self.expr-self.val(v))
    def __rsub__(self,v):return Formula(self.val(v)-self.expr)
    def __mul__(self,v):return Formula(self.expr*self.val(v))
    __rmul__=__mul__
    def __truediv__(self,v):return Formula(self.expr/self.val(v))
    def __neg__(self):return Formula(-self.expr)
    def __pow__(self,v):return Formula(self.expr**v)


def source_view(row,epsilon):
    """Formal Ad powers in exactly the actual Pstar source units."""
    part=lambda p,v:dict(Pstar_power=p,ordinary_logR_rows=v)
    add=lambda a,b:[x+y for x,y in zip(a,b)]
    mul=lambda a,p:[x*epsilon**p for x in a]
    return dict(modified_cylindrical_velocity_source_log_sectors=dict(
        theta=[part(1,add(row['U0'],mul(row['du'],1)))],
        axial=[part(1,mul(row['V'],1))],
        radial=[part(0,row['R0']),part(1,mul(row['dr'],1))]),
        modified_five_histories_in_original_normalized_units_Pstar_sectors=dict(
            m=[part(0,row['M0']),part(1,mul(row['dm'],1))],
            h=[part(0,add(row['H0'],mul(row['dh'],1)))],
            k=[part(0,row['K0']),part(1,mul(row['dk'],2))],
            e=[part(0,add(row['E0'],mul(row['de'],2)))],
            p=[part(0,add(row['P0'],mul(row['dp'],2)))]),
        modified_absolute_pressure_over_Pstar2_ordinary_logR_rows=add(row['P0'],mul(row['dp'],2)))


@functools.lru_cache(maxsize=1)
def source_model():
    """Exact difference of original source operators, including all products."""
    asts=recovered.histories.SourceAST();c=SimpleNamespace(mpf=lambda v:s.Rational(str(v)))
    derivative=lambda value,n=1:Formula(s.diff(value.expr,Z,n))
    env=dict(math=math,axial_derivative=derivative,axial_n=derivative,
        physical_operators=physical_operators,product_rows=physical.product_rows,
        shifted_rows=physical.shifted_rows)
    asts.replay('pulse_end_physical_C2','axial_operator_rows',env)
    asts.replay('current_modified_pre_stress_operator','sectors',env)
    stress=asts.replay('current_modified_pre_stress_operator','modified_pre_stress_rows',env)
    remainder=asts.replay('current_modified_pre_stress_operator','modified_pre_remainder_sectors',env)
    rows={name:[Formula(sum(math.comb(j,i)*s.Rational(1,2)**(j-i)*s.diff(f,Y,i)
        for i in range(j+1))) if name in ('R0','dr') else Formula(s.diff(f,Y,j))
        for j in range(5)] for name,f in FUNCTIONS.items()}
    changed=source_view(rows,EPS);baseline=source_view(rows,s.Integer(0))
    results={};checks={};cancellations={}
    for kind,fn,order in (('stress',stress,3),('remainder',remainder,2)):
        new=fn(c,D,Formula(Z),changed);old=fn(c,D,Formula(Z),baseline)
        results[kind]={};cancellations[kind]=[]
        for label,parts in new.items():
            results[kind][label]={}
            for name,part in parts.items():
                field='full_derivative_rows' if kind=='stress' else 'rows'
                expressions=[s.expand(a.expr-b.expr) for a,b in zip(part[field],old[label][name][field])]
                polys=[s.Poly(e,EPS) for e in expressions]
                if any(p.degree()>2 for p in polys if not p.is_zero):raise ValueError('Unexpected Ad degree')
                if any(p.nth(0)!=0 for p in polys):raise ArithmeticError('Unchanged baseline did not cancel')
                powers=[a for a in (1,2) if any(p.nth(a)!=0 for p in polys)]
                if not powers:cancellations[kind].append(label+'/'+name)
                for a in powers:
                    expr=polys[0].nth(a)
                    # The stress program includes the radius shift once;
                    # the remainder already includes its own internal shifts.
                    for j in range(order+1):
                        rate=s.Rational(str(part['mode'][0]))
                        want=sum(math.comb(j,i)*rate**(j-i)*s.diff(expr,Y,i) for i in range(j+1))
                        difference=s.expand(want-polys[j].nth(a))
                        if difference!=0 and s.cancel(difference)!=0:
                            raise ArithmeticError('Source rows differ: %s/%s/%s/Ad%d/y%d'%(kind,label,name,a,j))
                        checks[kind+'/'+label+'/'+name+'/Ad'+str(a)+'/y'+str(j)]=True
                    out=dict(mode=part['mode'],error_Ad_power=a,expression=expr,original_sector=name)
                    if kind=='remainder':out.update(beta=part['beta'],normalization_half=part['normalization_half'])
                    results[kind][label][name+'_Ad'+str(a)]=out
    # Baseline functions that cancel may not leak into the difference caps.
    atoms=set().union(*(p['expression'].atoms(AppliedUndef) for k in results.values() for q in k.values() for p in q.values()))
    if atoms.intersection(FUNCTIONS[n] for n in ('H0','K0','E0','P0')):
        raise ArithmeticError('Unchanged history/axis pressure leaked into difference')
    checks['nonzero_original_swirl_radial_velocity_and_meridional_history_retained']=all(FUNCTIONS[n] in atoms for n in ('U0','M0','R0'))
    q,tau=s.symbols('loglambda tau',real=True)
    if s.simplify(s.exp(2*q)*(1-Z**2)-tau).subs(q,s.log(tau/(1-Z**2))/2).simplify()!=0:
        raise ArithmeticError('Actual lambda-time correlation failed')
    checks['lambda_squared_times_one_minus_Z_squared_equals_tau']=True
    return results,dict(passed=all(checks.values()),identities=checks,exact_zero_baseline_difference_sectors=cancellations,
        source_error_sector_counts={k:{l:len(v) for l,v in r.items()} for k,r in results.items()},
        formal_Ad_is_factor_bookkeeping_not_a_fixed_N_perturbation=True,input_hashes=asts.hashes)


def indices(order):return ((j,k) for j in range(order+1) for k in range(order+1-j))


@functools.lru_cache(maxsize=None)
def polynomial_terms(expr,j,k,rate=0):
    """Differentiate the exact source, then collect its jet monomials."""
    rate=s.Rational(str(rate))
    value=sum(math.comb(j,i)*rate**(j-i)*s.diff(expr,Y,i,Z,k) for i in range(j+1))
    derivatives=value.atoms(s.Derivative)
    # Replace outer derivatives before the underlying applied functions.
    atoms=sorted(derivatives,key=str)+sorted(value.atoms(AppliedUndef),key=str)
    if not atoms:return [(s.cancel(value),[])]
    symbols=s.symbols('jet0:'+str(len(atoms)))
    lookup={}
    for atom,var in zip(atoms,symbols):
        fn=atom.expr if isinstance(atom,s.Derivative) else atom
        if fn not in FUNCTIONS.values():raise ValueError('Foreign source function')
        orders=dict(atom.variable_count) if isinstance(atom,s.Derivative) else {}
        if set(orders)-{Y,Z}:raise ValueError('Foreign source derivative')
        lookup[var]=(fn.func.__name__,int(orders.get(Y,0)),int(orders.get(Z,0)))
    poly=s.Poly(value.xreplace(dict(zip(atoms,symbols))),*symbols)
    terms=[]
    for powers,coefficient in poly.terms():
        terms.append((s.cancel(coefficient),[(lookup[v],p) for v,p in zip(symbols,powers) if p]))
    return terms


def bound_expression(c,expr,caps,z,delta,order,rate=0):
    out={}
    for j,k in indices(order):
        result=c.mpf(0)
        for coefficient,monomial in polynomial_terms(expr,j,k,rate):
            value=upper(c,abs(interval_expression(c,coefficient,z,delta,c.mpf(0))))
            for (name,a,b),power in monomial:value*=caps[name][a,b]**power
            result+=value
        out[j,k]=upper(c,result)
    return out


def compile_error_lift():
    """Keep original physical/completion algebra; change only error inputs."""
    asts=recovered.histories.SourceAST();fn=copy.deepcopy(asts.method('pulse_end_physical_C2','lift_physical_packet'))
    fn.decorator_list=[];changes=[]
    for node in ast.walk(fn):
        if isinstance(node,ast.Assign) and ast.unparse(node)=='source_errors = pulse_remainder_sectors(c, delta, None, z, velocity)':
            node.value=ast.parse("packet['actual_error_remainder_sectors']",mode='eval').body
            changes.append('remainder_difference_input')
        if isinstance(node,ast.Assign) and ast.unparse(node)=='parts = factor_logs(c, packet, sector[\'mode\'], sector[\'normalization_half\'])':
            node.value=ast.parse('error_factor_logs(c,packet,sector)',mode='eval').body
            changes.append('remainder_original_Ad_log_factor')
    if sorted(changes)!=['remainder_difference_input','remainder_original_Ad_log_factor']:
        raise ValueError('Unreviewed physical lift adaptation')
    def error_factor_logs(c,packet,sector):
        return dict(physical.factor_logs(c,packet,sector['mode'],sector['normalization_half']),
            original_Rd_amplitude=sector['error_Ad_power']*packet['exact_logAd'])
    env={**vars(physical),'error_factor_logs':error_factor_logs}
    exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),
        '<original completed tensor lift; error inputs only>','exec'),env)
    for stem in ('pulse_physical_bounds','collar_physical_C2','paper_interval_taylor'):
        name=('lei_ren_part1_'+stem+'.py') if stem=='paper_interval_taylor' else PREFIX+stem+'.py'
        # IntervalTaylor has a paper_ prefix, without compliant_.
        if stem=='paper_interval_taylor':name='lei_ren_part1_paper_interval_taylor.py'
        asts.hashes[name]=sha(name)
    return env[fn.name],dict(exact_AST_adaptations=changes,input_hashes=asts.hashes,
        original_physical_derivatives_diagonal_divergence_and_Cartesian_formulas_unchanged=True,
        mapper_receives_twice_actual_log_lambda=True)


class CurrentCompletedTensorErrorMajorants:
    def __init__(self,require_checked=True):
        self.recovered=recovered.CurrentRecoveredErrorMajorants();self.ctx=self.recovered.ctx
        self.data=self.recovered.data;self.mu=self.recovered.mu;self.delta=self.recovered.delta
        self.model,self.theorem=source_model();self.lift,self.lift_binding=compile_error_lift()
        receipt=json.loads((HERE/BASE_RECEIPT).read_bytes())
        if not receipt['all_passed'] or not receipt['full_actual_source_and_original_paper_programs_bound']:
            raise ValueError('Original arbitrary-source full stress binding required')
        for name,digest in receipt['input_hashes'].items():
            if sha(name)!=digest:raise ValueError('Changed original baseline source: '+name)
        if receipt['input_hashes'].get(BASE_VIEWS)!=sha(BASE_VIEWS):raise ValueError('Baseline packet is not checked')
        self.hashes={**self.recovered.hashes,recovered.RECEIPT:sha(recovered.RECEIPT),
            **receipt['input_hashes'],BASE_RECEIPT:sha(BASE_RECEIPT),**self.theorem['input_hashes'],
            **self.lift_binding['input_hashes'],Path(__file__).name:sha(Path(__file__).name)}
        for name,digest in self.hashes.items():
            if sha(name)!=digest:raise ValueError('Changed tensor error source: '+name)
        views=json.loads(gzip.decompress((HERE/BASE_VIEWS).read_bytes()));self.baselines={}
        for chart,key in (('O2','O2_full_buffer'),('O3','O3_full_transition'),('quiet','quiet_full_repair')):
            actual=views[key]['actual_source'];pre=actual['current_original_pre_source']
            for name,value in self.data['source']['accepted']['source_family'].items():
                if actual.get(name)!=value:raise ValueError('Foreign original baseline family')
            if not pre['radial_prefactors_differentiated_before_mixed_grid']:
                raise ValueError('Original native radial derivatives not included')
            read=parameters.numeric.transport.read_interval;c=self.ctx
            expected_chart={'O2':'O2_11_unit_buffer','O3':'O3_slope_mu','quiet':'O3_power_to_Rp'}[chart]
            coverage=read(c,actual['coordinate']);expected=c.mpf({'O2':(0,11),'O3':(0,1),'quiet':(1,2)}[chart])
            if pre['chart']!=expected_chart or coverage._mpi_!=expected._mpi_ or read(c,pre['Z'])._mpi_!=c.mpf([-1,1])._mpi_:
                raise ValueError('Original baseline does not cover the current whole source chart')
            source=dict(U0=pre['physical_velocity_pressure_y_Z_mixed4']['Utheta_over_Pstar'],
                R0=pre['physical_velocity_pressure_y_Z_mixed4']['Ur_over_current_sqrt_R_over_2'],
                M0=pre['physical_five_primitive_y_Z_mixed4']['Mz_over_current_R'])
            uz=pre['physical_velocity_pressure_y_Z_mixed4']['Uz']
            if any(read(c,v)._mpi_!=c.mpf(0)._mpi_ for v in uz.values()):raise ValueError('Nonzero original axial source needs extra cross sectors')
            # Saved primitive/radial grids differentiate their physical R
            # prefactors. Undo those shifts before using normalized functions.
            self.baselines[chart]={}
            for n,g in source.items():
                rate=c.mpf('.5') if n=='R0' else c.mpf(1) if n=='M0' else c.mpf(0)
                self.baselines[chart][n]={(j,k):upper(c,abs(sum(
                    math.comb(j,i)*(-rate)**(j-i)*read(c,g['y%d_Z%d'%(i,k)])
                    for i in range(j+1)))) for j,k in indices(4)}
        if require_checked:
            checked=json.loads((HERE/RECEIPT).read_bytes())
            if not checked['all_passed']:raise ValueError('Checked completed error bounds required')
            for name,digest in checked['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed checked tensor error dependency: '+name)

    def inputs(self,chart,coordinate,N):
        c=self.ctx
        if chart=='quiet':
            q=self.recovered.quiet_repair(N,coordinate);y=q['coordinate']
            f2=c.exp(-1-3*self.mu/2);G=q['local_disjoint_bump_logR_derivative_caps']
            theta=[upper(c,f2*self.mu*self.recovered.repair.R*v/N) for v in G]
            axial=[upper(c,f2*c.sqrt(self.mu)*self.recovered.repair.R*v/N) for v in G]
        elif chart in ('O2','O3'):
            q=self.recovered.modulation(coordinate,N);lo,hi=endpoints(q['coordinate'])
            if chart=='O2' and hi>0 or chart=='O3' and lo<0:raise ValueError('Split at actual O2/O3 source seam')
            t=q['coordinate'];norm=self.recovered.norms
            old=[upper(c,parameters.maximum(c,norm['T'][j],c.exp(-t/2)/2**j)) for j in range(5)]
            F=[[old[j]*math.factorial(k) for k in range(6)] for j in range(5)]
            profiles=parameters.normalized_envelopes(c,self.mu,N,q['local_cutoff_derivative_caps'],F,norm['L1'])
            theta=[upper(c,row[0]) for row in profiles['theta_increment_over_Pstar_over_Ad_majorants']]
            axial=[upper(c,row[0]) for row in profiles['modified_axial_over_Pstar_over_Ad_majorants']]
        else:raise ValueError('Expected O2, O3 or quiet current source chart')
        caps=dict(self.baselines[chart])
        for name,rows in (('du',theta),('V',axial)):
            caps[name]={(j,k):upper(c,rows[j]*math.factorial(k)) for j,k in indices(4)}
        for name,key in (('dm','m'),('dh','h'),('dk','k'),('de','e'),('dp','p')):
            rows=q['five_normalized_history_error_ordinary_logR4_axial5'][key]
            caps[name]={(j,k):rows[j][k] for j,k in indices(4)}
        rows=q['radial_error_over_Pstar_Ad_sqrtRover2_ordinary_logR4_axial5']
        caps['dr']={(j,k):upper(c,sum(math.comb(j,i)*c.mpf('.5')**(j-i)*rows[i][k]
            for i in range(j+1))) for j,k in indices(4)}
        return q,caps

    def query(self,chart,coordinate,N,Z=(-1,1),log_lambda=(-3,-1),theta=None,viscosity='1'):
        c=self.ctx;z=c.mpf(Z);q,caps=self.inputs(chart,coordinate,N);ll=c.mpf(log_lambda)
        lo,hi=endpoints(z)
        if lo< -1 or hi>1 or not all(mp.isfinite(v) for v in endpoints(z)+endpoints(ll)):
            raise ValueError('Finite Z subset[-1,1] and log(lambda) required')
        signed=lambda v:c.mpf([-1,1])*upper(c,v)
        columns={};errors={};native={}
        for kind,labels in self.model.items():
            native[kind]={}
            for label,parts in labels.items():
                native[kind][label]={}
                for name,part in parts.items():
                    order=3 if kind=='stress' else 2
                    bounded=bound_expression(c,part['expression'],caps,z,self.delta,order,part['mode'][0])
                    native[kind][label][name]={'y%d_Z%d'%jk:v for jk,v in bounded.items()}
                    if kind=='stress':
                        rp,bp,_,_=part['mode']
                        columns.setdefault(label,{})[name]=dict(mode=part['mode'],error_Ad_power=part['error_Ad_power'],
                            full_stress_mixed3_coefficient_enclosures={'s%d_Z%d'%jk:signed(v) for jk,v in bounded.items()},
                            exact_source_log_parts=dict(source_logR=c.mpf(str(rp))*q['native_logR'],
                                logPstar=bp*self.data['logP'],original_Rd_amplitude=part['error_Ad_power']*self.data['logAd'],
                                normalization=-c.ln(2)/2))
                    else:
                        rows=[IntervalTaylor(c,[signed(bounded[j,k])/math.factorial(k) for k in range(3-j)]) for j in range(3)]
                        errors.setdefault(label,{})[name]=dict(mode=part['mode'],error_Ad_power=part['error_Ad_power'],
                            rows=rows,beta=interval_expression(c,s.sympify(part['beta']),z,self.delta,c.mpf(0)),
                            normalization_half=part['normalization_half'])
        packet=dict(Z=z,s=q['coordinate'],exact_logR=q['native_logR'],exact_logAd=self.data['logAd'],
            exact_pulse_reference_logB_parts=dict(logPstar=self.data['logP']),exact_logD=c.mpf(0),exact_logH=c.mpf(0),
            full_meridional_stress_log_sectors=columns,actual_error_remainder_sectors=errors)
        # Original scalar abs() uses global mp precision. Keep it above the
        # interval context so none of its positive upper endpoints round down.
        with mp.workdps(c.dps+30):physical_bounds=self.lift(c,packet,self.delta,None,2*ll,theta,viscosity)
        physical_bounds.pop('requested_log_tau')
        physical_bounds['actual_log_lambda']=ll
        return dict(chart=chart,coordinate=q['coordinate'],finite_integer_N=N,Z=z,
            source_family=q['source_family'],native_error_absolute_mixed_caps=native,physical_error_bounds=physical_bounds,
            bound_kind='modified source minus the same unchanged original source',
            axial_closure_is_formal_source_bound=True,finite_physical_time_requires_strict_axial_domain=True,
            cached_modified_fixed_N_values_used=False,nonzero_original_radial_and_meridional_cross_terms_retained=True,
            same_original_pressure_datum_retained=True,completed_signed_tensor_error_bounds_available=True,
            **{k:False for k in OPEN})

    def at_physical_time(self,chart,coordinate,N,Z='0',log_tau='-1',theta=None,viscosity='1'):
        c=self.ctx;z=c.mpf(Z);lt=c.mpf(log_tau);lo,hi=endpoints(z)
        if lo<=-1 or hi>=1 or not all(mp.isfinite(v) for v in endpoints(z)+endpoints(lt)):
            raise ValueError('Finite physical time requires strict |Z|<1')
        ll=(lt-c.ln(1-z*z))/2
        result=self.query(chart,coordinate,N,z,ll,theta,viscosity)
        result['physical_error_bounds']['requested_log_tau']=lt
        result['physical_error_bounds']['exact_lambda_time_correlation_used']=True
        return result


def run():
    field=CurrentCompletedTensorErrorMajorants(require_checked=False)
    examples={name:field.query(chart,coordinate,N) for name,chart,coordinate,N in (
        ('before_support','O2','-3',1),('left_taper','O2','-1.99',37),
        ('transition','O3',('0','.5'),10**12),('after_cutoff','O3','.75',37),
        ('quiet_whole','quiet',(0,1),field.recovered.repair.quiet_threshold),
        ('terminal','quiet',('.9','1'),field.recovered.repair.quiet_threshold))}
    examples['actual_time']=field.at_physical_time('O3','.2',10**12,Z=('.3','.6'),log_tau=('-3','-1'))
    (HERE/VIEWS_NAME).write_bytes(gzip.compress((json.dumps(parameters.encoded(examples),indent=2)+'\n').encode(),mtime=0))
    result=dict(source_family=field.data['source']['accepted']['source_family'],
        exact_original_source_difference_theorem=field.theorem,original_completed_physical_lift_binding=field.lift_binding,
        example_names=list(examples),complete_error_views=VIEWS_NAME,completed_signed_tensor_error_bounds_available=True,
        variable_N_original_cross_terms_diagonal_divergence_and_remainder_bounds_available=True,
        **{k:False for k in OPEN},input_hashes={**field.hashes,VIEWS_NAME:sha(VIEWS_NAME)})
    (HERE/NAME).write_bytes((json.dumps(parameters.encoded(result),indent=2)+'\n').encode())
    print('Completed tensor error bounds generated; original cross terms, mixed3 stress/mixed2 remainder',flush=True)
    return result


if __name__=='__main__':run()
