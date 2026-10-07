"""Original paired Poisson parameter derivatives before whole-period bounds.

P=h^-1 W1 and H=h^-2 W2 retain u=(p2/dstar)*q. Pairing exact
derivatives removes unneeded q*u_Z products, with genuine linear q_Z kept.
"""
import ast
import inspect
import json
from pathlib import Path
import textwrap
import time
import sympy as sy
import lei_ren_part1_paper_compliant_current_native_periodic_C1_support_transport as accepted
import lei_ren_part1_paper_compliant_current_native_Rc_range_loss as loss

cutoff=accepted.cutoff;cover=accepted.cover;density=accepted.density;native=accepted.native;packets=accepted.packets
HERE,PREFIX,sha=accepted.HERE,accepted.PREFIX,accepted.sha;ep=accepted.ep;ZERO=accepted.ZERO;DZ=accepted.DZ
NAME=PREFIX+'current_native_paired_C1_transport.json'
RECEIPT=PREFIX+'current_native_paired_C1_transport_check.json'
GATE='current_original_paired_Poisson_Z_support_and_fixed_N_transport_executed'


def paired_parameter_theorem():
    r,z,w,E,psi=sy.symbols('r z sinpsi E psi',real=True);s=1-r*r;D=1-2*r*z+r*r
    cosE=((1+r*r)*z-2*r)/D;sinE=s*w/D;Er=2*w/D
    # W1=(E-psi)/(2r), P=sqrt(s)*W1, u*r_u=r*s.
    P=sy.sqrt(s)*(E-psi)/(2*r)
    Pr=sy.diff(P,r)+sy.diff(P,E)*Er
    if sy.simplify(P+r*s*Pr-s**sy.Rational(3,2)*w/D)!=0:
        raise ArithmeticError('Original paired P derivative identity failed')
    H=((2-3*s)*E+s*psi+2*r*sy.Symbol('sinE'))/(4*r*r)
    Hr=sy.diff(H,r)+sy.diff(H,E)*Er+sy.diff(H,sy.Symbol('sinE'))*cosE*Er
    paired=(2*H+r*s*Hr).subs(sy.Symbol('sinE'),sinE)
    target=E+r*s*w/D+s*s*w*(z-r)/(D*D)
    if sy.factor(paired-target)!=0:raise ArithmeticError('Original paired H derivative identity failed')
    if sy.expand(D-((z-r)**2+w*w)).subs(w*w,1-z*z)!=0:
        raise ArithmeticError('Original complex Poisson denominator identity failed')
    return dict(passed=True,original_exact_paired_derivative_identities=3,
        P_plus_u_Pu='s^(3/2)*sin(psi)/D; absolute upper<=sqrt(s)<=1',
        twoH_plus_u_Hu='E+r*s*sin(psi)/D+s²*sin(psi)*(cos(psi)-r)/D²',
        twoH_plus_u_Hu_absolute_upper='2pi+3',
        proof_P='D²-s²*sin²(psi)=((1+r²)*cos(psi)-2r)²>=0',
        proof_H='|E|<=2pi; |r*s*sin(psi)/D|<=1; 2|sin(psi)*(cos(psi)-r)|<=D; s²/D<=(1+|r|)²<=4',
        original_signed_r_and_zero_r_limits_retained=True,
        P_absolute_upper='4pi; inherited original W1 bound and h^-1<=1',
        Pu_absolute_upper='20pi; inherited W1_u<16pi and |(h^-1)_u|<=1',
        Hu_absolute_upper='1000pi; accepted original whole-period parameter theorem',
        qP_Z='q_Z*(P+u*P_u)+q²*c_Z*P_u',
        q2H_Z='(q²_Z/2)*(2H+u*H_u)+q³*c_Z*H_u',
        c_Z='original p2_Z/dstar; dstar has no slow derivatives',
        genuine_linear_q_Z_remains=True,outer_bounds_not_defining_field_values=True)


def paired_parameter_bounds(Q,QZ,R,RZ,roots,dstar,dstar_log,pi):
    cZ=accepted.serial.absolute_upper(roots['p2'][DZ].positive_divide(dstar,dstar_log))
    return QZ+R*cZ*(20*pi),RZ*((2*pi+3)/2)+R*Q*cZ*(1000*pi),cZ


def compile_paired_support():
    tree=ast.parse(inspect.getsource(accepted.periodic_C1_support))
    changes=[('(TZ*Q+T*QZ+T*Q*UZ)*(4*pi)+T*Q*UZ*(16*pi)','TZ*Q*(4*pi)+T*paired_qP_Z'),
        ('T*TZ*(4*pi)+firstZ*4+RZ*(400*pi)+R*UZ*(4000*pi)','T*TZ*(4*pi)+firstZ*4+paired_q2H_Z*4'),
        ('EZ*B_cap+EU*(TZ*A_cap+T*AZ_cap+AZ*Q*2+AU*QZ*2+AU*Q*UZ*10+AU*(T+1)*L*(1/(4*pi)))',
         'EZ*B_cap+EU*(TZ*A_cap+T*AZ_cap+AZ*Q*2+AU*paired_qP_Z*(1/(2*pi))+AU*(T+1)*L*(1/(4*pi)))')]
    tree,counts=density.replace_expressions(tree,changes);function=tree.body[0];inserts=0
    for i,node in enumerate(function.body):
        if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='firstZ' for t in node.targets):
            function.body[i:i]=ast.parse('paired_qP_Z,paired_q2H_Z,paired_c_Z=paired_parameter_bounds(Q,QZ,R,RZ,roots,dstar,dstar_log,pi)').body
            inserts+=1;break
    if inserts!=1:raise ValueError('Original whole-period derivative insertion site changed')
    ast.fix_missing_locations(tree);scope=dict(vars(accepted));scope['paired_parameter_bounds']=paired_parameter_bounds
    exec(compile(tree,'<original-whole-period-Z-support-with-paired-Poisson-factors>','exec'),scope)
    return scope['periodic_C1_support'],counts


PAIRED_SUPPORT,SUPPORT_ADAPTER_COUNTS=compile_paired_support()


def paired_C1_support(roots,qrows,q2rows,log_a_lower,dstar_log):
    got=PAIRED_SUPPORT(roots,qrows,q2rows,log_a_lower,dstar_log)
    got['record'].update(original_paired_Poisson_parameter_theorem=paired_parameter_theorem(),
        original_whole_period_Z_support_AST_replacements=SUPPORT_ADAPTER_COUNTS,
        original_unpaired_q_times_u_Z_bounds_replaced_by_exact_paired_derivatives=True,
        original_linear_q_Z_and_nonzero_p2_Z_retained=True)
    return got


def compile_paired_query(source_query,kernel,log_a_lower,support_rows):
    tree=ast.parse(inspect.getsource(accepted.compile_periodic_query))
    tree,counts=density.replace_expressions(tree,[
        ("periodic_C1_support(source['roots'],qrows,q2rows,log_a_lower,dstar)",
         "paired_C1_support(source['roots'],qrows,q2rows,log_a_lower,dstar)")])
    scope=dict(vars(accepted));scope['paired_C1_support']=paired_C1_support
    exec(compile(tree,'<original-periodic-query-with-exact-paired-parameter-support>','exec'),scope)
    return scope['compile_periodic_query'](source_query,kernel,log_a_lower,support_rows)


class NativePairedC1Density(accepted.NativePeriodicC1Density):
    def __init__(self,owner):
        super().__init__(owner)
        self.hashes={**self.hashes,Path(__file__).name:sha(Path(__file__).name),Path(loss.__file__).name:sha(Path(loss.__file__).name)}
        self.service.bind_hashes(self.hashes)

    @native.inlet.source_precision
    def spatial_query(self,chart,Z,coordinate,N):
        tree=ast.parse(textwrap.dedent(inspect.getsource(accepted.NativePeriodicC1Density.spatial_query)))
        scope=dict(vars(accepted));scope['compile_periodic_query']=compile_paired_query
        exec(compile(tree,'<original-cutoff-query-with-paired-Poisson-Z-support>','exec'),scope)
        got=scope['spatial_query'](self,chart,Z,coordinate,N)
        got['record'].update(original_paired_Poisson_support_before_nonlinear_density=True,
            original_paired_Poisson_parameter_theorem=paired_parameter_theorem())
        return got


class NativePairedC1Oracle(accepted.NativePeriodicC1Oracle):
    def __init__(self,role_owner,built=None):
        super().__init__(role_owner,built);self.owner=NativePairedC1Density(self.original_density_owner)
        self.hashes={**self.hashes,**self.owner.hashes};self.service.bind_hashes(self.hashes)

    @native.inlet.source_precision
    def density_frame(self,*,chart,Z,coordinate,N):
        frame=super().density_frame(chart=chart,Z=Z,coordinate=coordinate,N=N)
        frame.record['phase_solver_backend']='original branch-local q squared and exact paired Poisson whole-period Z supports'
        return frame


class NativePairedC1Transport(accepted.NativePeriodicC1Transport):
    def __init__(self,role_owner):
        super().__init__(role_owner)
        checked=json.loads((HERE/accepted.RECEIPT).read_bytes())
        if not checked.get('all_passed') or not checked.get(accepted.GATE) or checked['source_family']!=self.family:
            raise ValueError('Checked original whole-period Z support baseline required')
        self.oracle=NativePairedC1Oracle(role_owner,self.built)
        self.hashes={**self.hashes,**self.oracle.hashes,**checked['input_hashes'],accepted.RECEIPT:sha(accepted.RECEIPT),
            Path(__file__).name:sha(Path(__file__).name)};self.service.bind_hashes(self.hashes)


def range_frontier(owner,live):
    """Current weighted contribution ranks, using applied density supports."""
    if len(live['cells'])!=24 or live['target'] is None:raise ValueError('Actual original full24 route required')
    c=owner.ctx;coords=owner.coordinates;target=live['target'];suffix={key:coords.scalar(1) for key in loss.current.RATES}
    rows=[];dominant={};primitives={};count=0
    for cell in reversed(live['cells']):
        values={key:cell['values'][key]*suffix[key] for key in loss.current.RATES}
        jets={key:cell['Z_derivatives'][key]*suffix[key] for key in loss.current.RATES}
        got=loss.current.fixed_N_target_rows(values,jets,target['A'],target['AZ'],target['logA'],coords.scalar(target['mu']),target['logmu'])
        label=cell['record']['label'];per={}
        for group in ('values','Z_derivatives'):
            for key,value in got[group].items():
                upper=loss.upper(value);name=key+('' if group=='values' else '_Z')
                per[name]=dict(exact_zero=value.zero,log_absolute_upper=None if upper is None else c.mpf(upper))
                if upper is not None and (name not in dominant or upper>dominant[name][0]):dominant[name]=(upper,label)
        if cell['frame'] is not None:
            source=cell['frame'].record['actual_original_spatial_source']
            for support in source['original_C0_A_B_source_support_ranges_before_nonlinear_density']:
                actual={key:(dict(exact_zero=True) if support[key]['exact_original_zero'] else
                    support[key]['original_C0_range'] if support[key]['original_tighter_formal_range_retained']
                    else support[key]['used_C0_range']) for key in ('A','B_over_Pstar')}
                actual.update(A_Z=support['original_A_Z'],B_Z_over_Pstar=support['original_B_Z_over_Pstar'])
                for key,value in actual.items():
                    count+=1
                    if value['exact_zero']:continue
                    upper=ep(packets.interval(c,value['log_absolute_upper']))[1]
                    if key not in primitives or upper>primitives[key][0]:primitives[key]=(upper,label)
        rows.append(dict(label=label,chart=cell['record']['chart'],normalized_final_target_contribution_ranges=per,
            final_suffix_decay={key:value.record() for key,value in suffix.items()}))
        suffix={key:cell['factors'][key]['decay']*suffix[key] for key in loss.current.RATES}
    return dict(original_continuous_cells=24,weighted_original_cell_ranges=list(reversed(rows)),
        dominant_normalized_target_contribution_cells={key:dict(label=label,log_absolute_upper=c.mpf(upper)) for key,(upper,label) in dominant.items()},
        dominant_actual_density_primitive_range_cells={key:dict(label=label,log_absolute_upper=c.mpf(upper)) for key,(upper,label) in primitives.items()},
        actual_density_primitive_ranges_inspected=count,applied_C0_and_Z_support_ranges_used=True,
        per_cell_target_hulls_do_not_prove_joint_cancellation=True,source_range_bounds_not_actual_NS_residual=True)


@native.inlet.source_precision
def run(role_owner,*,return_live=False,baseline_live=None):
    began=time.monotonic();owner=NativePairedC1Transport(role_owner);got=owner.route()
    if got['history'] is None or len(got['cells'])!=24:raise ArithmeticError('All24 original continuous cells required')
    baseline=json.loads((HERE/accepted.NAME).read_bytes());comparison=accepted.accepted.accepted.compare_targets(owner.ctx,got['record'],baseline)
    frontier=range_frontier(owner,got)
    old_frontier=None if baseline_live is None else range_frontier(*baseline_live)
    result=dict(**got['record'],**{GATE:True},original_paired_parameter_theorem=paired_parameter_theorem(),
        comparison_with_checked_unpaired_periodic_Z_target_ranges=comparison,
        strict_target_absolute_upper_reductions=sum(row['strict_absolute_upper_reduction'] for row in comparison.values()),
        current_weighted_target_range_frontier=frontier,checked_unpaired_weighted_target_range_frontier=old_frontier,
        original_linear_q_Z_sensitivity_retained=True,useful_repair_contraction_or_actual_controls_established=False,
        execution_seconds=time.monotonic()-began,input_hashes=owner.hashes,
        scope='Exact original paired qP/q²H Z identities and whole-period support before density, complete24-cell original fixed-N transport and refreshed actual weighted source contributions. Genuine linear q_Z remains; no actual controls/global N/closure/recursion admission.')
    (HERE/NAME).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    print('Original paired Poisson supports transported through24 cells;strict target upper reductions:',result['strict_target_absolute_upper_reductions'],flush=True)
    return (result,owner,got) if return_live else result
