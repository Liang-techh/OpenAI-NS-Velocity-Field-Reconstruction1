"""Focused original entrance physical checks and independent Cartesian derivatives."""
import ast
import copy
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
import sympy as s
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_pulse_entrance_physical_C2 import (
    CompliantPulseEntrancePhysicalC2,entrance_physical_binding,compiled_main_exit_lift,
    FALSE_FLAGS,PREFIX,DOMAIN,source_precision)
from lei_ren_part1_paper_compliant_pulse_entrance_similarity_C4 import (
    CompliantPulseEntranceSimilarityC4,compiled_entrance_exporter)
from lei_ren_part1_paper_compliant_pulse_main_exit_fixture_integrals import run as fixture_integrals
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval

HERE=Path(__file__).parent


def independent_entrance_fixture_setup():
    """Exact local derivatives of original startup primitives, not a global polynomial field."""
    c=MPIntervalContext();c.dps=130
    v,z=s.symbols('y Z',real=True)
    mu=s.Rational(1,5);delta=s.Rational(1,2);bp=s.Rational(1,2)+mu
    B0,D0,Rp,Xp=map(s.Rational,('.45','.07','2.8','.3'))
    C=1/(1+z*z);xx=mu*v;p=1+2*mu
    ap=s.Rational(105,100)+3*z/100-z*z/50
    incoming=(s.Rational(1,10)*(z+z**3),s.Rational(7,100)*(z+z**3))
    ein=-s.Rational(15,1000)+(z*z+2*z**4+z**6)/500
    J0=s.Rational(1,5)+z*z/20;P0=-7*C*C/20+z/50
    xi0,z0=mp.mpf('.01'),mp.mpf('.31');v0=xi0/mp.mpf('.2')
    constant=lambda value:s.Float(mp.nstr(value,115),115)
    sigma=lambda x:mp.mpf(0) if x<=0 else mp.mpf(1) if x>=1 else 1/(1+mp.exp(1/x**2-1/(1-x)**2))
    gp_cache={}
    def gp(x):
        key=mp.nstr(x,110)
        if key not in gp_cache:
            gp_cache[key]=mp.quad(lambda b:sigma(50*b),[0,x/2,x],method='gauss-legendre')
        return gp_cache[key]
    g0=gp(xi0)
    gder=[g0]+[mp.mpf('.2')**n*50**(n-1)*mp.diff(sigma,mp.mpf('.5'),n-1) for n in range(1,9)]
    local=v-s.Rational(1,20)
    gpoly=sum(constant(value)*local**n/math.factorial(n) for n,value in enumerate(gder))
    kernels=[];kernel_values=[]
    for i in (1,2):
        rate=mp.mpf('.5')-i*mp.mpf('.2')
        head=mp.quad(lambda b:mp.exp(rate*b/mp.mpf('.2'))*sigma(50*b),[0,xi0/2,xi0])
        k0=(g0-mp.exp(-rate*v0)*head)/rate;kernel_values.append(k0)
        rows=[k0]
        for n in range(8):rows.append(gder[n]-rate*rows[n])
        kernels.append(sum(constant(value)*local**n/math.factorial(n) for n,value in enumerate(rows)))
    partial_energy=mp.quad(lambda a:mp.exp(-2*a)*gp(a)**2,[0,xi0/2,xi0],method='gauss-legendre')
    vals=fixture_integrals()['values'];Kpulse=s.Float(vals['Kpulse'],115)
    energy0=s.exp(s.Rational(1,50))*(ein+ap*ap*constant(partial_energy)/mu-(1-s.exp(-s.Rational(1,50)))/(4*mu))
    energy_rows=[energy0]
    for n in range(8):
        square=sum(math.comb(n,k)*constant(gder[k])*constant(gder[n-k]) for k in range(n+1))
        energy_rows.append(2*mu*energy_rows[n]+ap*ap*square-(s.Rational(1,2) if n==0 else 0))
    energy=sum(value*local**n/math.factorial(n) for n,value in enumerate(energy_rows))
    future=s.exp(26)/mu*(ap*ap*Kpulse-(1-s.exp(-26))/4+mu*ein)+D0*D0*J0
    m,n=[ap*kernels[i-1]+s.exp(-(s.Rational(1,2)-i*mu)*v)*incoming[i-1] for i in (1,2)]
    Bh=ap*gpoly
    B=B0*s.exp(-bp*v);R=Rp*s.exp(v);H=s.exp(-(1-mu)*v)
    Ut=B*C;Uz=Ut*Bh
    Mt=s.sqrt(2)*R**s.Rational(3,2)*Ut*(1/(1-mu)+(Xp-1/(1-mu))*H)
    Mz=R*Ut*m;Mtz=s.sqrt(2)*R**s.Rational(3,2)*Ut*Ut*n;Mzt=R*Ut*Ut*energy
    P=B*B*(-C*C/(2*p)+s.exp(-p*(13-xx)/mu)*(P0+C*C/(2*p)))
    d,L=1-z*z,1-delta*z*z
    transport=-R+(1-delta)*z*Mz+d*s.diff(Mz,z)
    It=Ut*transport/(L*s.sqrt(2*R))+((1-delta/2)*Mt
        -(1-delta)*z*s.diff(Mt,z)/2-d*s.diff(Mtz,z)+(2*delta-1)*z*Mtz)/(2*L*R)
    Iz=(transport*Uz+(1-delta)*(Mz-z*s.diff(Mz,z))/2+2*delta*z*Mzt
        -d*s.diff(Mzt,z)+R*(2*(1+delta)*z*P-d*s.diff(P,z)))/(L*s.sqrt(2*R))
    Tt=It+(2*s.diff(Ut,v)-Ut)/s.sqrt(2*R)
    Tz=Iz+s.sqrt(2*R)*s.diff(Uz,v)/R
    Ur=(2*z*R*Uz-(1-delta)*z*Mz-d*s.diff(Mz,z))/(L*s.sqrt(2*R))
    direct={name:s.lambdify((v,z),expr,'mpmath',cse=True) for name,expr in
        (('Ur',Ur),('Ut',Ut),('Uz',Uz),('P',P),('Tt',Tt),('Tz',Tz))}
    def jet(expr):
        return IntervalTaylor(c,[c.mpf(s.lambdify(z,s.diff(expr,z,k),'mpmath')(z0)/math.factorial(k)) for k in range(6)])
    sim=object.__new__(CompliantPulseEntranceSimilarityC4)
    sim.ctx=c;sim.mu=c.mpf('.2');sim.delta=c.mpf('.5');sim.Xp=c.mpf('.3');sim.cells=128
    sim.ap=jet(ap);sim.incoming=[jet(expr) for expr in incoming];sim.incoming_energy=jet(ein)
    sim.future=jet(future);sim.J0=jet(J0);sim.P0=jet(P0)
    sim.Pin=c.mpf('.12');sim.logP=c.ln(c.mpf('.45'));sim.logU=c.mpf(0)
    sim.logRp=c.ln(c.mpf('2.8'));sim.finite=c.ln(c.mpf('.07'))+1/sim.mu
    sim.exporter,_=compiled_entrance_exporter()
    def exact_kernel(ctx,mm,x,rate,cells):
        index=0 if endpoints(rate)[0]>mp.mpf('.2') else 1
        return dict(enclosure=c.mpf(kernel_values[index]))
    sim.exporter.__globals__['partial_linear_kernel']=exact_kernel
    sim.exporter.__globals__['partial_future_energy']=lambda *args:c.mpf(mp.mpf(vals['Kpulse'])-partial_energy)
    fixture=object.__new__(CompliantPulseEntrancePhysicalC2)
    fixture.ctx=c;fixture.similarity=sim;fixture.lift,_=compiled_main_exit_lift()
    print('Independent original entrance integral and local derivative fixture constructed',flush=True)
    return c,fixture,direct,xi0,z0,v0,Rp


def independent_entrance_Cartesian_oracle():
    """Reuse unchanged independent Cartesian operators with the original entrance fixture."""
    path=HERE/(PREFIX+'pulse_main_exit_physical_C2_check.py')
    tree=ast.parse(path.read_text(encoding='utf8'))
    original=next(node for node in tree.body if isinstance(node,ast.FunctionDef) and node.name=='independent_main_Cartesian_oracle')
    withnode=original.body[-1]
    start=next(i for i,node in enumerate(withnode.body) if isinstance(node,ast.Assign)
        and any(isinstance(target,ast.Name) and target.id=='checks' for target in node.targets))
    tail=copy.deepcopy(withnode.body[start:])
    class EntranceScope(ast.NodeTransformer):
        def visit_Attribute(self,node):
            node=self.generic_visit(node)
            if isinstance(node.value,ast.Name) and node.value.id=='fixture' and node.attr=='main_exit':node.attr='entrance'
            return node
        def visit_Constant(self,node):
            if isinstance(node.value,str):
                node.value=node.value.replace('Independent main Cartesian oracle','Independent entrance Cartesian oracle').replace('nonzero main error','nonzero entrance error')
            return node
        def visit_keyword(self,node):
            node=self.generic_visit(node)
            if node.arg=='original_startup_integrals_and_exact_plateau_primitives_used':
                node.arg='original_startup_integrals_and_exact_local_derivatives_used'
            return node
    tail=[EntranceScope().visit(node) for node in tail]
    setup=ast.parse('c,fixture,direct,xi0,z0,v0,Rp=independent_entrance_fixture_setup()').body
    fn=ast.FunctionDef(name='replayed_entrance_oracle',args=ast.arguments(posonlyargs=[],args=[],kwonlyargs=[],kw_defaults=[],defaults=[]),
        body=[ast.With(items=[ast.withitem(context_expr=ast.parse('mp.workdps(100)',mode='eval').body)],body=setup+tail)],
        decorator_list=[])
    module=ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[]))
    env=dict(globals());exec(compile(module,str(path)+' [entrance fixture]','exec'),env)
    result=env['replayed_entrance_oracle']()
    result.update(Cartesian_operator_formulas_reused_without_arithmetic_changes=True,
        local_polynomial_used_only_to_represent_original_derivatives_at_one_point=True,
        original_entrance_xi='.01',derivative_jet_order=8,
        Cartesian_operator_source_sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    return result


@source_precision
def run():
    companion=CompliantPulseEntrancePhysicalC2()
    name=PREFIX+'pulse_entrance_physical_C2.json';raw=(HERE/name).read_bytes();record=json.loads(raw)
    actual=companion.report()
    if encode(pack(actual))!=record:raise ValueError('Current entrance physical report differs')
    for path,digest in record['input_hashes'].items():
        if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:raise ValueError('Entrance physical source changed: '+path)
    if record['domain']!=DOMAIN:raise ValueError('Original entrance domain required')
    if record['source_and_physical_join_binding']!=entrance_physical_binding(companion.records):raise ValueError('Physical source join differs')
    c=MPIntervalContext();c.dps=300;finite=zeros=0
    def rowcheck(row):
        nonlocal finite,zeros
        lo,hi=endpoints(read_interval(c,row['signed_coefficient']))
        if not all(mp.isfinite(v) for v in (lo,hi)):raise ArithmeticError('Nonfinite physical row')
        if row['exact_zero']:
            if lo or hi or row['log_absolute_upper'] is not None:raise ArithmeticError('Incorrect structural zero')
            zeros+=1
        else:
            if not all(mp.isfinite(v) for v in endpoints(read_interval(c,row['log_absolute_upper']))):raise ArithmeticError('Nonfinite physical log bound')
            finite+=1
    groups=('physical_cylindrical_stress_mixed3','physical_cylindrical_stress_divergence_mixed2',
        'physical_three_component_remainder_mixed2','completed_theta_theta_stress_mixed2')
    for pointkey in ('whole_original_entrance','original_inlet','early_ordinary_y_chart','common_entrance_main'):
        point=record[pointkey]
        for flag in ('actual_pulse_entrance_physical_decomposition_constructed',
            'original_inlet_and_entrance_main_completed_physical_interfaces_verified',
            'all_local_incoming_nonlinear_history_cross_products_retained',
            'frozen_source_logs_not_differentiated_again','equivalent_forward_energy_not_double_counted'):
            if not point[flag]:raise ValueError('Missing entrance physical scope '+flag)
        for group in groups:
            order=3 if group=='physical_cylindrical_stress_mixed3' else 2
            required={'r'+str(i)+'_z'+str(j) for i in range(order+1) for j in range(order+1-i)}
            sectors=point[group]
            grids=sectors.values() if group=='completed_theta_theta_stress_mixed2' else (
                grid for parts in sectors.values() for grid in parts.values())
            for grid in grids:
                if set(grid)!=required:raise ValueError('Incomplete physical mixed derivative grid')
                for row in grid.values():rowcheck(row)
        for label in ('exact_completed_tensor_radial_divergence','exact_physical_divergence'):
            if endpoints(read_interval(c,point[label]))!=(0,0):raise ArithmeticError('Lost physical structural zero')
        for flag in FALSE_FLAGS:
            if point[flag] or record[flag]:raise ValueError('Physical scope overclaimed '+flag)
    inlet=record['original_inlet']['physical_three_component_remainder_mixed2']
    if any(not row['exact_zero'] for grid in inlet['axial'].values() for row in grid.values()):raise ArithmeticError('Flat inlet axial remainder is not zero')
    for label in ('radial','theta'):
        if not any(not row['exact_zero'] for grid in inlet[label].values() for row in grid.values()):raise ArithmeticError('Nonzero inlet history removed: '+label)
    main=companion.main_exit([-1,1],'.02',theta=None,viscosity='.7')
    for group in groups:
        if encode(pack(main[group]))!=record['common_entrance_main'][group]:raise ValueError('Same physical xi=.02 source differs: '+group)
    print('Entrance completed physical rows and original inlet/main joins admitted',flush=True)
    oracle=independent_entrance_Cartesian_oracle()
    hashes=dict(record['input_hashes']);hashes[name]=hashlib.sha256(raw).hexdigest()
    hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    for stem in ('pulse_main_exit_fixture_integrals.py','pulse_main_exit_fixture_integrals.json','pulse_main_exit_physical_C2_check.py'):
        path=PREFIX+stem;hashes[path]=hashlib.sha256((HERE/path).read_bytes()).hexdigest()
    result=dict(all_passed=True,actual_five_defect_family_sha256=companion.family,implicit_source_sha256=companion.source,
        domain=DOMAIN,input_hashes=hashes,actual_pulse_entrance_physical_decomposition_constructed=True,
        original_inlet_and_entrance_main_completed_physical_interfaces_verified=True,
        actual_finite_signed_physical_rows=finite,exact_physical_zero_rows=zeros,
        current_whole_entrance_physical_report_recomputed=True,
        actual_source_and_join_identities=len(record['source_and_physical_join_binding']['identities']),
        admitted_full_physical_operator_identities=record['source_and_physical_join_binding']['admitted_full_physical_operator_identities'],
        independent_entrance_Cartesian_oracle=oracle,source_caps_used_as_defining_field_values=False,
        **{flag:False for flag in FALSE_FLAGS})
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('PASS original entrance physical tensor and full remainder; entrance cone/global/recursion pending',flush=True)
    return result


if __name__=='__main__':run()
