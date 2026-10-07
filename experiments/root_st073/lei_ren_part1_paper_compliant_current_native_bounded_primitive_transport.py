"""Original C0 primitive support ranges before nonlinear density transport.

Theorems bound source ranges, never define field values. Tighter original
tiny formal factors survive; all original A_Z/B_Z and E/V rows stay intact.
"""
import json
from pathlib import Path
import time
import sympy as sy
import lei_ren_part1_paper_compliant_current_native_cutoff_range_transport as accepted

density=accepted.density;cutoff=density.cutoff;packets=accepted.packets;native=accepted.native
HERE,PREFIX,sha=accepted.HERE,accepted.PREFIX,accepted.sha;ep=accepted.ep
NAME=PREFIX+'current_native_bounded_primitive_transport.json'
RECEIPT=PREFIX+'current_native_bounded_primitive_transport_check.json'
GATE='current_original_supported_C0_primitives_before_fixed_N_C1_transport_executed'


def primitive_support_theorem():
    a,b,q,eta,Delta,sigma=sy.symbols('a b q eta Delta sigma',real=True)
    kappa=a+b*b/a;v=a*(1+b*b/(a*a)+2*q*q)
    if sy.simplify(v-kappa-2*a*q*q)!=0:raise ArithmeticError('Original v/nu normalization identity differs')
    active_v=2+Delta+sigma*sigma*(2*eta-Delta)
    if sy.expand(3-active_v-(1-2*eta+(1-sigma*sigma)*(2*eta-Delta)))!=0:
        raise ArithmeticError('Original active v<=3 decomposition differs')
    return dict(passed=True,original_A_support=cutoff.primitive_identity_checks(),
        original_B_identity='B/E=-(a*T1/(2*pi)+b*phi)/2; T1=integral_0^psi t',
        original_direction='t=t0+2*q*w/h; w=sum_{n>=1} r^(n-1)*cos(n*psi), |r|<1',
        exact_period_mean_w=0,exact_period_mean_w_squared='1/(2*(1-r^2))=h^2/2',
        exact_period_mean_t_squared='t0^2+2*q^2=nu-1',
        Cauchy_Schwarz_argument='|a*T1/(2*pi)|<=sqrt(a*(v-a)); b^2<=a*(v-a); sqrt(a*(v-a))<=v/2',
        active_v_argument='q!=0 => Delta<eta; v=kappa+sigma^2*(2*eta-Delta)<=2+2*eta<=3',
        bounds='|A|<=5/4; B=E*beta with |beta|<=3/2, uniformly over original p2/sign/phase',
        exact_flat_case='q=0 => original A=B=0',
        derivative_rows_not_inferred_from_C0_bounds=True,bounds_not_defining_field_values=True)


def signed_C0_support_range(original,bound):
    """Keep a tighter original coordinate, else a signed source-support range.

    Both inputs enclose the same source. No midpoint or bound endpoint is
    selected as a source value. The support product B=E*beta stays factored.
    """
    original.coerce(bound);c=original.ctx
    if original.zero:return original,dict(original_tighter_formal_range_retained=True,exact_original_zero=True)
    old=original.record();new=bound.record();lo,hi=ep(original.coefficient)
    old_upper=ep(old['log_absolute_upper'])[1];support_upper=ep(new['log_absolute_upper'])[1]
    if lo>0 or hi<0:
        minimum=min(abs(lo),abs(hi))
        lower=ep(original.scale.evaluate()+c.ln(c.mpf(minimum)))[0]
        if lower>support_upper:raise ArithmeticError('Original C0 range contradicts the proved source-support theorem')
    if old_upper<=support_upper:
        return original,dict(original_tighter_formal_range_retained=True,exact_original_zero=False,
            original_C0_range=old,source_support_outer_range=new)
    lower,upper=ep(bound.coefficient)
    if lo>0:lower=max(lower,0)
    if hi<0:upper=min(upper,0)
    got=cutoff.prior.ScaledEnclosure(bound.scale,c.mpf((lower,upper)),original.ledger)
    return got,dict(original_tighter_formal_range_retained=False,exact_original_zero=False,
        original_C0_range=old,source_support_outer_range=new,used_C0_range=got.record(),
        original_source_sign_retained=True,outer_range_coordinate_not_new_source_value=True)


class NativeBoundedPrimitiveDensity(density.NativeCutoffDensityCover):
    def __init__(self,owner):
        super().__init__(owner);self.original_kernel=self.kernel;self.support_rows=[]
        def supported_kernel(E,E_Z,V,V_Z,primitives,N):
            A,ar=signed_C0_support_range(primitives['A'],primitives['A'].scalar(self.ctx.mpf((-1.25,1.25))))
            B,br=signed_C0_support_range(primitives['B_over_Pstar'],E*self.ctx.mpf((-1.5,1.5)))
            changed=dict(primitives,A=A,B_over_Pstar=B)
            if changed['A_Z'] is not primitives['A_Z'] or changed['B_Z_over_Pstar'] is not primitives['B_Z_over_Pstar']:
                raise ValueError('Original primitive ordinary Z rows must remain identical objects')
            self.support_rows.append(dict(A=ar,B_over_Pstar=br,original_A_Z=primitives['A_Z'].record(),
                original_B_Z_over_Pstar=primitives['B_Z_over_Pstar'].record(),
                original_A_Z_and_B_Z_same_objects=True,B_support_retains_original_E_factor=True,
                original_E_E_Z_V_V_Z_and_all_cross_terms_retained=True))
            return self.original_kernel(E,E_Z,V,V_Z,changed,N)
        self.kernel=supported_kernel
        self.hashes={**self.hashes,Path(__file__).name:sha(Path(__file__).name)};self.service.bind_hashes(self.hashes)

    @native.inlet.source_precision
    def spatial_query(self,chart,Z,coordinate,N):
        start=len(self.support_rows);got=super().spatial_query(chart,Z,coordinate,N)
        rows=self.support_rows[start:];record=got['record']
        record['original_A_over_N_formal_factor_and_A_Z_preserved']=all(
            row['A']['original_tighter_formal_range_retained'] for row in rows)
        record.update(original_C0_A_B_source_support_ranges_before_nonlinear_density=rows,
            original_A_C0_arithmetic_coordinates_may_be_source_theorem_restricted=True,
            original_A_Z_B_Z_and_nonzero_E_V_rows_unchanged=True,
            tighter_original_tiny_formal_factors_not_replaced_by_scalar_caps=True,
            exact_original_C0_primitive_support_theorem=primitive_support_theorem())
        return got


class NativeBoundedPrimitiveOracle(density.NativeCutoffFactoredOracle):
    def __init__(self,role_owner,built=None):
        super().__init__(role_owner,built);self.owner=NativeBoundedPrimitiveDensity(self.original_density_owner)
        self.hashes={**self.hashes,**self.owner.hashes};self.service.bind_hashes(self.hashes)

    @native.inlet.source_precision
    def density_frame(self,*,chart,Z,coordinate,N):
        frame=super().density_frame(chart=chart,Z=Z,coordinate=coordinate,N=N)
        frame.record['phase_solver_backend']='original cutoff/signed-u first jets with proved C0 A/B support ranges before density'
        return frame


class NativeBoundedPrimitiveTransport(accepted.NativeCutoffRangeTransport):
    def __init__(self,role_owner):
        super().__init__(role_owner)
        checked=json.loads((HERE/accepted.RECEIPT).read_bytes())
        if not checked.get('all_passed') or not checked.get(accepted.GATE) or checked['source_family']!=self.family:
            raise ValueError('Checked same-original complete fixed-N range baseline required')
        self.oracle=NativeBoundedPrimitiveOracle(role_owner,self.built)
        self.hashes={**self.hashes,**self.oracle.hashes,**checked['input_hashes'],accepted.RECEIPT:sha(accepted.RECEIPT),
            Path(__file__).name:sha(Path(__file__).name)}
        self.service.bind_hashes(self.hashes)


def compare_targets(c,new,baseline):
    rows={}
    for group in ('original_Rc_target_C0_ranges','original_Rc_target_Z_ranges'):
        for key,value in new[group].items():
            old=baseline[group][key];label=key+('' if group.endswith('C0_ranges') else '_Z')
            before=None if old['exact_zero'] else packets.interval(c,old['log_absolute_upper'])
            after=None if value['exact_zero'] else value['log_absolute_upper']
            rows[label]=dict(original_log_absolute_upper=before,bounded_primitive_log_absolute_upper=after,
                strict_absolute_upper_reduction=before is not None and (after is None or ep(after)[1]<ep(before)[0]),
                log_upper_reduction=None if before is None or after is None else before-after)
    return rows


@native.inlet.source_precision
def run(role_owner,*,return_live=False):
    began=time.monotonic();owner=NativeBoundedPrimitiveTransport(role_owner);got=owner.route()
    baseline=json.loads((HERE/accepted.NAME).read_bytes())
    if got['history'] is None or len(got['cells'])!=24:raise ArithmeticError('All original24 bounded-primitive cells required')
    comparison=compare_targets(owner.ctx,got['record'],baseline)
    result=dict(**got['record'],**{GATE:True},original_primitive_support_theorem=primitive_support_theorem(),
        comparison_with_checked_original_fixed_N_target_ranges=comparison,
        strict_target_absolute_upper_reductions=sum(row['strict_absolute_upper_reduction'] for row in comparison.values()),
        original_derivative_rows_retained_not_capped=True,
        useful_repair_contraction_or_actual_controls_established=False,
        execution_seconds=time.monotonic()-began,input_hashes=owner.hashes,
        scope='Proved original C0 A/B support ranges before nonlinear density and complete fixed-N24-cell C1 transport. Tighter original tiny factors and all original derivative/cross rows retained. Target range improvement is not actual controls/global N/closure or recursion.')
    (HERE/NAME).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    print('Original bounded primitives transported through24 cells;strict target upper reductions:',result['strict_target_absolute_upper_reductions'],flush=True)
    return (result,owner,got) if return_live else result
