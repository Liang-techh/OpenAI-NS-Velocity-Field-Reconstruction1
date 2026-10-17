"""Exact conditional weighted Green kernel; no actual complex n1 norm yet.

The rejected approximate log-bound prototype is not used. Actual source
sup norms, outward arithmetic and the Xh embedding remain required inputs.
"""
from __future__ import annotations

import ast
from fractions import Fraction
import math
from pathlib import Path
from .original_finite_N_signed_source import load,sha,digest,require

PARENT='lei_ren_part1_paper_compliant_current_original_radial_IBP.json'
PARENT_SHA256='67945b4e9a8f2b1cdaa9a255b42a3bd9e5377f514cb12409ba41e9f17f957af0'
RAW='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_inner_volterra_functions.py'
POWERS=(2,2,2,2,1,1)
SOURCE_POWERS={4:1,5:0,6:0}


def word_majorant(Msrc,C0,C1,R,sigma,order):
    """Conditional word bound; the caller must supply genuine sup bounds."""
    require(type(order) is int and order>=1,'Positive integer word order required')
    values=[]
    for value in (Msrc,C0,C1,R,sigma):
        if type(value) is int:value=Fraction(value)
        require(isinstance(value,Fraction) or hasattr(value,'_mpi_'),
            'Exact rational or outward interval arithmetic required')
        values.append(value)
    Msrc,C0,C1,R,sigma=values
    return Msrc*(R/2)**order*(C0+C1*order/sigma)**order/math.factorial(order)


def cauchy_tail_factor(degree,*,target=Fraction(4),radius=Fraction(9,2)):
    require(type(degree) is int and degree>=0 and 0<=target<radius,'Invalid Taylor domain')
    q=target/radius
    return q**(degree+1)/(1-q)


def _degree(n):
    """Minimum x power and parity; leading fields are conditional even series."""
    if isinstance(n,ast.Name):return (1,1) if n.id=='x' else ((2,0) if n.id=='r' else (0,0))
    if isinstance(n,ast.Constant) and isinstance(n.value,(int,float)):return 0,0
    if isinstance(n,ast.UnaryOp):return _degree(n.operand)
    require(isinstance(n,ast.BinOp),'Unsupported matrix expression: '+ast.dump(n))
    a,pa=_degree(n.left);b,pb=_degree(n.right)
    if isinstance(n.op,(ast.Add,ast.Sub)):
        require(pa==pb,'Mixed radial parity');return min(a,b),pa
    if isinstance(n.op,ast.Mult):return a+b,(pa+pb)%2
    if isinstance(n.op,ast.Div):
        require(b==pb==0,'Radial denominator requires an axis proof');return a,pa
    if isinstance(n.op,ast.Pow):
        require(isinstance(n.right,ast.Constant) and type(n.right.value) is int and n.right.value>=0,
            'Nonnegative integer power required')
        return a*n.right.value,pa*n.right.value%2
    raise ValueError('Unsupported matrix algebra')


def _entries(root,parent):
    path=root/RAW;pin=parent['input_hashes'][path.name]
    require(sha(path)==pin,'Actual raw operator changed')
    tree=ast.parse(path.read_text(encoding='utf-8'))
    f=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='raw_matrix_expressions']
    require(len(f)==1,'Unique actual raw matrix source required')
    dictionaries={target.id:n.value for n in f[0].body if isinstance(n,ast.Assign)
        and isinstance(n.value,ast.Dict) for target in n.targets
        if isinstance(target,ast.Name) and target.id in ('B0','B1')}
    require(set(dictionaries)=={'B0','B1'},'Actual matrix literals required')
    rows=[]
    for label,d in dictionaries.items():
        for key,value in zip(d.keys,d.values):
            require(isinstance(key,ast.Constant) and isinstance(key.value,str),'Literal entry required')
            i,j=map(int,key.value.split(','))
            if i not in SOURCE_POWERS:
                require(label=='B0' and key.value in ('1,5','2,6','3,6'),'Unexpected kinematic row');continue
            degree,parity=_degree(value);shift=POWERS[j-1]-SOURCE_POWERS[i]
            require(degree+shift>=2 and (parity+shift)%2==0,'Weighted entry is not rho times a regular coefficient')
            rows.append(dict(matrix=label,entry=key.value,expression=ast.unparse(value),
                raw_min_x_degree=degree,axis_weight_shift=shift,
                after_division_by_rho_min_x_degree=degree+shift-2,even_in_x_after_weighting=True))
    require(len(rows)==16,'All sixteen off-kinematic entries required')
    return dict(passed=True,path=RAW,source_sha256=pin,entries=rows,source_executed=False,
        identity='x^(p_j-q_i)*B_ij=rho*A_ij; division means analytic axis continuation')


def exact_controls():
    import sympy as s
    t,u,k=s.symbols('t u k',positive=True)
    kernels=[t*(1-t*t)/2,t*s.log(1/t),-t*(1-t*t)/2,t,t**3,t]
    constants=[s.integrate((-v if i==2 else v),(t,0,1)) for i,v in enumerate(kernels)]
    expected=[s.Rational(1,8),s.Rational(1,4),s.Rational(1,8),s.Rational(1,2),s.Rational(1,4),s.Rational(1,2)]
    # bar_h(rho*t^2)=rho*t^2*A*v, ds=2*rho*t*dt.
    transformed=[s.simplify(s.expand_log(v.subs(t,s.sqrt(u))*s.sqrt(u)/2,force=True)) for v in kernels]
    wanted=[u*(1-u)/4,u*s.log(1/u)/4,-u*(1-u)/4,u/2,u*u/2,u/2]
    f=u*s.log(1/u)/4
    checks=dict(six_source_integral_constants=constants==expected,
        six_exact_rho_segment_kernels=all(s.simplify(a-b)==0 for a,b in zip(transformed,wanted)),
        log_kernel_stationary_point=s.simplify(s.diff(f,u).subs(u,1/s.E))==0,
        log_kernel_maximum=s.simplify(f.subs(u,1/s.E)-1/(4*s.E))==0,
        factorial_log_integral=s.simplify(s.integrate(s.log(t),(t,1,k))-(k*s.log(k)-k+1))==0,
        raw_degree16_gives_weighted_degree15=cauchy_tail_factor(15)==9*Fraction(8,9)**16)
    require(all(checks.values()),'New exact weighted Green control failed')
    return dict(passed=True,checks=checks,source_constants=[str(v) for v in constants],
        rho_segment_kernels=[str(v) for v in transformed],uniform_operator_kernel_upper='1/2',
        domination='0<=u<=1; u(1-u)<=1/4; u*log(1/u)<=1/e',source_values_evaluated=False)


def build(root,parent):
    root=Path(root);path=root/'experiments/root_st073'/PARENT
    require(sha(path)==PARENT_SHA256 and digest(parent)==digest(load(path)),'Actual accepted parent required')
    entries=_entries(root,parent);controls=exact_controls()
    return dict(passed=True,proof_scope='Actual operator parity and exact conditional weighted Green/word kernel only',
        parent_report=PARENT,parent_report_sha256=PARENT_SHA256,source_only=True,
        axis_weights=list(POWERS),normalized_source=['g4/x','g5','g6'],
        normalized_unknown='v=(F1/x^2,Uz1/x^2,K1/x^2,P1/x^2,F1x/x,Uz1x/x)',
        zero_axis_extension_by_parity=True,actual_weighted_matrix_entries=entries,exact_controls=controls,
        source_norm='Msrc=max(||g4/x||/2,||g5||/4,||g6||/2)',
        operator_norms='C_a=max_i sum_j sup|A^a_ij| on common complex tube',
        word_bound='Msrc*(R/2)^k*(C0+k*C1/sigma)^k/k!',
        Z_loss_allocation='total gap sigma, sigma/k for each of k applications',
        factorial_proof='sum(log(j),j=1..k)>=integral_1^k log(t)dt>=k*log(k)-k',
        geometric_beta='e*R/2*(C0+C1/sigma)',conditional_solution_norm='||v||<=Msrc/(1-beta)<=2*Msrc if beta<=1/2',
        kinematic_preconditioning_already_in_six_kernels=True,extra_I_plus_J_factor_applied=False,
        complex_path='s=rho*u; absolute nested paths give |rho|^k/k!',
        radial_tail=dict(weighted_degree_N='M_v*(r/R)^(N+1)/(1-r/R)',
            raw_base_degree_N='r*M_v*(r/R)^N/(1-r/R), N>=1',
            prefix_map='W1..W4=rho*v; raw degree16 gives weighted degree15',
            odd_rows='W5,W6=x*v5,x*v6; divide by x before rho Cauchy',
            odd_raw_weighted_degree_N='sqrt(r)*M_v*(r/R)^(N+1)/(1-r/R)',
            odd_raw_from_base_degree_N='sqrt(r)*M_v*(r/R)^N/(1-r/R), N>=1',
            target='4',candidate_radius='9/2',ratio='8/9',prefactor='9'),
        actual_complex_source_sup_bounds_installed=False,actual_directed_solution_norm_installed=False,
        actual_radial_Taylor_remainder_certified=False,actual_collar_used_as_analytic_core_extension=False,
        accepted_owner_fixture_recurrence_or_quadrature_replayed=False,physical_scale_or_N_materialized=False,
        missing_inputs=['Holomorphy of actual coefficient fields and nonvanishing of L and other denominators on the common complex tube',
            'Actual Xh weighted-norm embedding into a complex-rho disk sup norm',
            'Outward actual n0 source and A0/A1 norms on nested complex Z tubes',
            'Outward Msrc and beta<1 at the chosen R>4',
            'Actual coefficient scaling and ordinary-Z Cauchy remainder binding'])
