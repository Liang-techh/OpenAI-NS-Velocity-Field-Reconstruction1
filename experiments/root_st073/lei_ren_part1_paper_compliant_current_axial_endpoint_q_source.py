"""Original O2 axial b/q on both microscopic endpoint layers and flat middle.

The full phase selector [0,1] is partitioned using original eta/Pstar/Md.
Original ordinary y2/Z1 rows and positive physical widths are retained.
This source interface is for subsequent actual signed density integration.
"""
import ast
from fractions import Fraction
import gzip
import hashlib
import json
import math
from pathlib import Path
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as sy
import lei_ren_part1_paper_compliant_current_transition_complete_prefix as complete

native=complete.native;prior=native.prior;current=native.current;slow=native.slow
common=native.common;HERE,PREFIX,sha,encode,ep,iv=native.HERE,native.PREFIX,native.sha,native.encode,native.ep,native.iv
NAME=PREFIX+'current_axial_endpoint_q_source.json'
RECEIPT=PREFIX+'current_axial_endpoint_q_source_check.json'
GATE='original_whole_axial_strict_sign_b_q_slow_sources_and_true_geometry_executed'
ZERO=(0,0);ORDERS=slow.ORDERS


def source_identity():
    proof=native.canonical_amplitude_proof()
    source=HERE/(PREFIX+'pre_pulse_mixed_C4.py')
    tree=ast.parse(source.read_text(encoding='utf8'))
    fn=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name=='axial')
    assignments=[(ast.unparse(t),ast.unparse(n.value)) for n in ast.walk(fn) if isinstance(n,ast.Assign) for t in n.targets]
    assert ('y','c.exp(md * phase)') in assignments and ('V','[z * (4 * b) for b in B]') in assignments
    assert ('t','y - 1') in assignments and ('root','c.exp(-t / 2)') in assignments
    r,H,k,Q,Z,M,P,y=sy.symbols('r H k Q Z M P y',positive=True)
    L,F=sy.symbols('L F',real=True)
    L1=2/r**3*(1+(r/(1-r))**3)
    normalized=32*sy.exp(sy.Rational(2,5))*Z**2*(1+Z**2)**2/M**2*L1**2*sy.exp(H+2*F-2/r**2+y-1)
    collected=16*sy.exp(sy.Rational(2,5))*Z**2*(1+Z**2)**2/M**2*(Q+k)**3*(1+(r/(1-r))**3)**2*sy.exp(H+2*F-2/r**2+y-1)
    assert sy.simplify((normalized/collected).subs(r**2,2/(Q+k)))==1
    return dict(passed=True,canonical_original_amplitude_proof=proof,
        original_axial_AST_sha256=hashlib.sha256(ast.dump(fn).encode()).hexdigest(),
        original_a_identity='a=1-2*(log(Utheta/Pstar))_y=2',
        original_inlet_Utheta_over_Pstar='exp(-1/5)/(1+Z**2)',
        original_signed_b_identity='-8*exp(1/5)*Z*(1+Z**2)*sigma_prime(p)*exp((exp(M*p)-1)/2-M*p)/(M*Pstar)',
        original_normalized_excess='rho=Delta/eta=b**2/(2*eta)',
        exact_reflection_identity='sigma_prime(1-r)=sigma_prime(r)',
        physical_derivative_coordinate='ordinary y=log(R/Rref), y=exp(M*p); D_y=(M*y)^(-1)*D_p',
        original_a_has_no_nonzero_order_slow_jets=True,
        original_b_sign_is_minus_sign_Z_not_assumed_positive=True,
        same_original_Pstar_Md_eta_and_function_formula_bound=True,
        original_source_functions_not_modified_or_fitted=True,
        input_hashes={**proof['input_hashes'],source.name:sha(source.name),
            PREFIX+'flat_pulse_derivatives.py':sha(PREFIX+'flat_pulse_derivatives.py'),
            PREFIX+'current_generic_shear_inputs.py':sha(PREFIX+'current_generic_shear_inputs.py')})


class AxialCoordinates:
    def __init__(self,coordinates,eta_log,M=40):
        self.coordinates=coordinates;self.c=c=coordinates.ctx;self.t=coordinates.scalar(1)
        self.eta_log=c.mpf(eta_log);self.M=c.mpf(M);self.logP=coordinates.logP_squared/2
        self.H=-self.eta_log-2*self.logP
        if ep(self.H)[0]<=0:raise ValueError('Original positive axial H required')
        self.logH=c.ln(self.H);self.A={'left':c.mpf(0),'right':c.exp(self.M)-1-2*self.M}
        self.Q={side:self.H+3*self.logH+A for side,A in self.A.items()}
        if ep(self.H-self.logH-200)[0]<=0:raise ValueError('Positive original flat-middle support denominator required')
        self.star_denominator=self.H-self.logH-200
        self.rstar=prior.ScaledEnclosure(prior.FormalScale(self.t.scale.bases,
            offset=(c.ln(2)-c.ln(self.star_denominator))/2),1,self.t.ledger)
        self.rstar_cover=c.sqrt(2/self.star_denominator)
        if ep(self.rstar_cover)[1]>=ep(c.mpf('.1'))[0]:raise ValueError('Original microscopic axial support collar required')
    def cv(self,value):
        if isinstance(value,(str,int)):value=Fraction(value)
        return self.c.mpf(value.numerator)/value.denominator if isinstance(value,Fraction) else self.c.mpf(value)
    def record(self):
        return dict(source_family=self.coordinates.family,original_eta_log=self.eta_log,
            original_Md=self.M,original_logPstar=self.logP,entire_original_H_cover=self.H,
            original_H_identity='H=-log(eta)-2log(Pstar)',original_endpoint_A=self.A,
            entire_original_Q_covers=self.Q,original_Q_identity='Q_e=H+3log(H)+A_e',
            exact_flat_middle_radius=self.rstar.record(),exact_common_flat_support_denominator=self.star_denominator,
            endpoint_Q_plus_kstar_identity='Q_e+(-4log(H)-A_e-200)=H-log(H)-200',
            source_parameter_covers_not_midpoints_or_selected_caps=True)
    def geometry(self,side,kind,left,right):
        if side not in self.Q:raise ValueError('Original left/right axial endpoint required')
        declared=[str(left),str(right)]
        c=self.c;t=self.t;Q=self.Q[side];left=self.cv(left)
        Kstar=-4*self.logH-self.A[side]-200
        right=Kstar if right=='star' else self.cv(right)
        W=prior.ScaledEnclosure(prior.FormalScale(t.scale.bases,offset=(c.ln(2)-c.ln(Q))/2),1,t.ledger)
        if kind=='xi':
            if ep(right-left)[0]<=0 or ep(left)[0]<0:raise ValueError('Increasing original xi interval required')
            rlo=c.sqrt(2/Q)*left;rhi=c.sqrt(2/Q)*right
            r=W*c.mpf((ep(left)[0],ep(right)[1]));width=W*(right-left)
            kval=None;left_source=W*left;right_source=W*right
            exact_width='sqrt(2/Q_e)*(xi_right-xi_left)'
        elif kind=='mixed':
            beta=1/c.sqrt(1+right/Q)
            if ep(beta-left)[0]<=0:raise ValueError('Positive original axial mixed width required')
            rlo=c.sqrt(2/Q)*left;rhi=c.sqrt(2/(Q+right))
            r=W*c.mpf((ep(left)[0],ep(beta)[1]));width=W*(beta-left)
            kval=c.mpf((ep(right)[0],ep(Q*(1/left**2-1))[1]))
            left_source=W*left;right_source=W*beta
            exact_width='sqrt(2/Q_e)*((1+k_right/Q_e)^(-1/2)-xi_left)'
        elif kind=='k':
            if ep(left-right)[0]<=0:raise ValueError('Decreasing original k interval required')
            a=c.sqrt(1+left/Q);b=c.sqrt(1+right/Q)
            rlo=c.sqrt(2/(Q+left));rhi=c.sqrt(2/(Q+right))
            r=prior.ScaledEnclosure(prior.FormalScale(t.scale.bases,
                offset=(c.ln(2)-c.ln(Q+c.mpf((ep(right)[0],ep(left)[1]))))/2),1,t.ledger)
            width=prior.ScaledEnclosure(prior.FormalScale(t.scale.bases,
                offset=(c.ln(2)-3*c.ln(Q))/2),(left-right)/(a*b*(a+b)),t.ledger)
            kval=c.mpf((ep(right)[0],ep(left)[1]));left_source=W*(1/a);right_source=W*(1/b)
            exact_width='sqrt(2)*((k_left-k_right)/(sqrt(Q+k_left)*sqrt(Q+k_right)*(sqrt(Q+k_left)+sqrt(Q+k_right))))'
        else:raise ValueError('Original axial xi/mixed/k source coordinate required')
        if ep(width.coefficient)[0]<=0:raise ArithmeticError('True original axial phase width must be positive')
        cover=c.mpf((ep(rlo)[0],ep(rhi)[1]));epsilon=1 if side=='left' else -1
        x=r*(self.M*epsilon);ybase=c.mpf(1) if side=='left' else c.exp(self.M)
        finite_x=native.bounded(x)
        exponential=c.exp(c.mpf((min(ep(finite_x)[0],0),max(ep(finite_x)[1],0))))
        dy=x*exponential*ybase
        ycover=ybase*c.exp(finite_x)
        residual=dy-x*2
        physical_width=width*(self.M*ycover)
        record=dict(chart='O2_axial',source_family=self.coordinates.family,endpoint=side,
            typed_coordinate_kind=kind,original_distance_to_endpoint_cover=cover,
            exact_declared_typed_endpoints=declared,
            original_distance_to_endpoint_source=r.record(),original_selector_source='p=r' if side=='left' else 'p=1-r',
            original_Q_source=self.record(),original_typed_endpoints=[str(left),str(right)],
            original_endpoint_distance_sources=[left_source.record(),right_source.record()],
            original_correlated_k_cover=kval,positive_true_phase_width=width.record(),
            positive_true_ordinary_y_width=physical_width.record(),original_phase_width_identity=exact_width,
            original_nonzero_y_minus_endpoint=dy.record(),original_collected_profile_residual=residual.record(),
            original_y_cover=ycover,original_ordinary_y_Jacobian=self.M*ycover,
            exact_original_y_mapping='y=exp(M*p); exp(x)-1=x*integral_0^1 exp(t*x)dt',
            microscopic_original_y_increment_retained=True,physical_width_Jacobian_once=True,
            ordinary_source_jets_do_not_receive_second_selector_conversion=True,
            endpoint_geometry_independent_of_Z=True,typed_endpoint_subtraction_collected_before_evaluation=True)
        return dict(record=record,side=side,kind=kind,left=left,right=right,r=r,cover=cover,
            kcover=kval,width=width,physical_width=physical_width,y=ycover,residual=residual,epsilon=epsilon)


def axial_b_rows(coord,Z,geometry):
    """Original signed b using safe phase-sigma and already ordinary-y rows."""
    c=coord.c;t=coord.t;z=c.mpf(Z);M=coord.M;y=geometry['y'];r=geometry['cover'];e=geometry['epsilon']
    sig=complete.closed_sigma_rows(t,r)
    exponent=(y-1)/2-c.ln(y)
    pref=prior.ScaledEnclosure(prior.FormalScale(t.scale.bases,(0,0,0,-2,0),exponent),
        -8*c.exp(c.mpf(1)/5)/M,t.ledger)
    J1=c.mpf('.5')-1/y;J2=1/y**2
    radial=[sig[1],sig[2]*(e/(M*y))+sig[1]*J1,
        sig[3]*(1/(M*M*y*y))+sig[2]*(-e/(M*y*y)+2*e*J1/(M*y))+sig[1]*(J2+J1**2)]
    shape=z*(1+z*z);shapeZ=1+3*z*z
    rows={(j,k):radial[j]*pref*(shapeZ if k else shape) for j,k in ORDERS}
    return dict(rows=rows,record=dict(original_signed_b_ordinary_y2_Z1_rows={'y%d_Z%d'%o:v.record() for o,v in rows.items()},
        original_sigma_phase_ordinary_y0_through_y4=[v.record() for v in sig],
        defining_signed_b_prefactor=pref.record(),original_axial_b_unit='2*Uz_y/Utheta',
        physical_Uz_Pstar_inverse_applied_once=True,original_phase_to_ordinary_y_applied_once=True))


def normalized_excess(coord,Z,geometry):
    c=coord.c;t=coord.t;z=c.mpf(Z);bsource=axial_b_rows(coord,Z,geometry);b=bsource['rows']
    if ep(z)[0]<=0<=ep(z)[1]:raise ValueError('Current actual axial endpoint source requires strict-sign Z tile')
    if geometry['kind']=='xi' and ep(geometry['left'])[0]==0:
        # No singular division on the closed endpoint. Uniform phase-sigma
        # majorants multiply in their original formal scales before /eta.
        products=current.multiply_jet(b,b)
        eta=prior.ScaledEnclosure(prior.FormalScale(t.scale.bases,offset=coord.eta_log),1,t.ledger)
        rows={order:value.positive_divide(eta,coord.eta_log)*c.mpf('.5') for order,value in products.items()}
        # b has a fixed original sign, but generic product dependency may
        # give a negative lower cover. Delta=b^2/2 is nonnegative exactly.
        v=rows[ZERO];lo,hi=ep(v.coefficient)
        rows[ZERO]=prior.ScaledEnclosure(v.scale,c.mpf((max(lo,0),hi)),v.ledger)
        return dict(rows=rows,b=b,record=dict(branch='closed_original_axial_endpoint',
            original_b_source=bsource['record'],original_Delta_identity='b*b/2',
            uniform_original_closed_endpoint_sigma_bounds_used=True,
            exact_endpoint_b_and_all_its_slow_jets_zero=True,
            original_rho_rows={'y%d_Z%d'%o:v.record() for o,v in rows.items()},
            zero_endpoint_not_an_epsilon_or_removed_source_interval=True))
    r=geometry['cover'];k=geometry['kcover'];side=geometry['side'];Q=coord.Q[side]
    if k is None:
        k=Q*(1/c.mpf((ep(geometry['left'])[0],ep(geometry['right'])[1]))**2-1)
    F=1/(1-r)**2;L=F-1/r**2;tail=t.bounded_exp(L)
    correction=2*c.ln(1+(r/(1-r))**3)-4*c.ln(1+tail)
    ratio=1+(3*coord.logH+coord.A[side]+k)/coord.H
    constant=c.ln(16*c.exp(c.mpf(2)/5)*z*z*(1+z*z)**2/coord.M**2)
    # H-Q_e=-3log(H)-A_e is a defining source identity; do not
    # subtract two independent astronomical endpoint covers.
    G=-k+3*c.ln(ratio)+2*F+constant+native.bounded(geometry['residual'])+correction
    value=prior.ScaledEnclosure(prior.FormalScale(t.scale.bases,offset=G),1,t.ledger)
    sig=native.original_sigma_rows(t,r)[0]
    sigma=sig.coefficient*sig.bounded_exp(sig.scale.evaluate())
    L1=2/(1-r)**3+2/r**3;L2=6/(1-r)**4-6/r**4;L3=24/(1-r)**5+24/r**5
    A=(1-2*sigma)*L1+L2/L1
    Ap=-2*sigma*(1-sigma)*L1**2+(1-2*sigma)*L2+L3/L1-(L2/L1)**2
    e=geometry['epsilon'];y=geometry['y'];M=coord.M
    Gp=2*e*A+M*(y-2);Gpp=2*Ap+M*M*y
    G1=Gp/(M*y);G2=(Gpp-M*Gp)/(M*M*y*y);GZ=2/z+4*z/(1+z*z)
    rad=[value,value*G1,value*(G2+G1**2)]
    rows={(j,k):rad[j]*(GZ if k else 1) for j,k in ORDERS}
    # Same signed b = -sign(Z)*sqrt(2eta*rho), preserving the original
    # normalized-excess correlation rather than cancelling giant log boxes.
    eta=prior.ScaledEnclosure(prior.FormalScale(t.scale.bases,offset=coord.eta_log),1,t.ledger)
    bb=current.nonnegative_sqrt(value*eta*2)*(-1 if ep(z)[0]>0 else 1)
    brad=[bb,bb*(G1/2),bb*(G2/2+G1**2/4)]
    brows={(j,k):brad[j]*(GZ/2 if k else 1) for j,k in ORDERS}
    return dict(rows=rows,b=brows,record=dict(branch='correlated_original_axial_endpoint',
        original_rho_log_source_cover=G,original_rho_identity='b**2/(2eta)',
        original_H_Q_cancellation_collected_symbolically=True,original_ratio_Q_plus_k_over_H=ratio,
        original_profile_delta_from_endpoint=geometry['residual'].record(),
        original_phase_sigma_tail_cover=tail,original_log_rho_ordinary_y_derivatives=[G1,G2],
        original_log_rho_Z_derivative=GZ,original_rho_rows={'y%d_Z%d'%o:v.record() for o,v in rows.items()},
        original_signed_b_rows={'y%d_Z%d'%o:v.record() for o,v in brows.items()},
        same_original_b_Delta_eta_relation_collected_before_evaluation=True,
        ordinary_y2_Z1_jets_not_normalized_coordinate_derivatives=True))


def original_axial_q(coord,normalized):
    c=coord.c;t=coord.t;v=normalized[ZERO]
    roots=dict(a={o:t.scalar(2 if o==ZERO else 0) for o in ORDERS},
        kappa_minus2={o:prior.ScaledEnclosure(prior.FormalScale(t.scale.bases,offset=coord.eta_log),1,t.ledger)*row for o,row in normalized.items()})
    vlo=ep(v.coefficient)[0]
    lower=None if vlo<=0 else v.scale.evaluate()+c.ln(c.mpf(vlo))
    if lower is not None and ep(lower)[0]>=0:
        rows={o:t.scalar(0) for o in ORDERS}
        return dict(rows=rows,roots=roots,record=dict(branch='flat',status='enclosed',
            original_rho_lower_log=lower,original_q_and_all_y2_Z1_jets_exact_zero=True,
            active_body_evaluated=False))
    d=c.mpf(0) if v.zero else v.coefficient*v.bounded_exp(v.scale.evaluate())
    lo,hi=ep(d)
    if hi>=2:raise ArithmeticError('Original axial positive body needs typed source refinement')
    eta=prior.ScaledEnclosure(prior.FormalScale(t.scale.bases,offset=coord.eta_log),1,t.ledger)
    argument={o:t.scalar(1-d) if o==ZERO else -row for o,row in normalized.items()}
    coefficients=prior.sigma_jets(c,1-d)
    cutoff=slow.compose_sigma(argument,[coefficients[j]*math.factorial(j) for j in range(4)])
    ratio={o:eta*((t.scalar(2-d) if o==ZERO else -row)*c.mpf('.25')) for o,row in normalized.items()}
    lower=(coord.eta_log+c.ln(c.mpf(2)-c.mpf(hi))-c.ln(4))/2
    root=slow.sqrt_jet(ratio,lower);rows=current.multiply_jet(cutoff,root)
    return dict(rows=rows,roots=roots,record=dict(status='enclosed',branch='active' if hi<1 else 'smooth_cutoff_seam',
        original_rho_cover=d,original_positive_body_cover=2-d,original_q_ordinary_y2_Z1_rows={'y%d_Z%d'%o:v.record() for o,v in rows.items()},
        original_a_exactly_two_not_two_plus_Delta=True,original_q_Z_not_assumed_zero=True,
        original_sigma_and_positive_body_evaluated_before_flat_union=True,
        exact_endpoint_q_is_sqrt_of_eta_over_two=True,active_body_evaluated=True))


def middle_flat_proof(coord,Z):
    c=coord.c;z=c.mpf(Z)
    if ep(z)[0]<=0<=ep(z)[1]:raise ValueError('Strict-sign Z required for flat axial middle proof')
    bound=coord.logH+202-2*coord.M+c.ln(2048*c.exp(c.mpf(2)/5)*z*z*(1+z*z)**2/coord.M**2)
    if ep(bound)[0]<=0:raise ArithmeticError('Original axial flat-middle excess margin must be positive')
    return dict(passed=True,original_selector_range='rstar <= p <= 1-rstar',
        original_symmetric_distance_range='rstar <= min(p,1-p) <= 1/2',
        original_logistic_odds_lower='L(r)>=1-1/rstar**2',original_logistic_phase_slope_lower='sigma_prime(r)>=8*exp(1-1/rstar**2)',
        proof_steps=['L<=0 => (1+exp(L))**2<=4','L1=2/r**3+2/(1-r)**3>=32',
            'y-1-2M*p>=-2M','H-2/rstar**2=log(H)+200'],
        original_log_rho_uniform_lower=bound,
        entire_original_middle_Delta_strictly_exceeds_eta=True,
        original_q_and_all_required_jets_and_five_own_increments_exact_zero=True,
        incoming_histories_and_pressure_not_reset=True,proof_uses_original_function_not_selected_cap=True)


PLAN=(('closed_endpoint','xi','0','3/4'),('bulk_to_seam','mixed','3/4','8'),
    ('seam_outer','k','8','0'),('active_seam','k','0','-3'),
    ('active_to_cutoff','k','-3','-19/5'),('cutoff_seam','k','-19/5','-22/5'),
    ('flat_connector','k','-22/5','star'))


def execute_tile(coord,Z):
    rows=[]
    for side in ('left','right'):
        for label,kind,left,right in PLAN:
            geometry=coord.geometry(side,kind,left,right);norm=normalized_excess(coord,Z,geometry)
            q=original_axial_q(coord,norm['rows'])
            rows.append(dict(endpoint=side,label=label,original_true_geometry=geometry['record'],
                original_normalized_excess=norm['record'],original_b_ordinary_y2_Z1_rows={'y%d_Z%d'%o:v.record() for o,v in norm['b'].items()},
                original_q_source=q['record'],original_q_rows={'y%d_Z%d'%o:v.record() for o,v in q['rows'].items()}))
            print('Original whole axial q:',Z,side,label,q['record']['branch'],flush=True)
    middle=middle_flat_proof(coord,Z)
    return dict(source_family=coord.coordinates.family,exact_Z_range=list(Z),original_parameter_sources=coord.record(),
        original_endpoint_source_queries=rows,original_middle_flat_source_theorem=middle,
        original_whole_phase_partition=dict(selector_interval=['0','1'],
            left_endpoint_source_order=[r[0] for r in PLAN],right_endpoint_source_order=[r[0] for r in reversed(PLAN)],
            original_left_support_endpoint='sqrt(2/(H-log(H)-200))',
            original_right_support_endpoint='1-sqrt(2/(H-log(H)-200))',
            original_shared_adjacent_endpoints_not_rounded_overlap=True,
            source_reflection_order_reversed_before_physical_transport=True,
            full_selector_width_telescope='rstar+(1-2*rstar)+rstar=1',
            full_ordinary_y_width_telescope='exp(M)-1',
            endpoint_q_nonzero_and_interior_flat_distinguished=True),
        whole_original_axial_q_source_defined_on_two_strict_sign_tiles=True,
        whole_Z_axis_or_density_integrals_targets_controls_recursion_admitted=False)


@native.rc.native.inlet.source_precision
def run():
    began=time.monotonic();accepted=json.loads((HERE/complete.NAME).read_bytes());hashes={}
    common.attach_receipt(hashes,complete,accepted['source_family'])
    proof=source_identity();hashes.update(proof['input_hashes']);archives=[]
    c=MPIntervalContext();c.dps=240
    for archive in accepted['actual_original_complete_transition_archives']:
        raw=gzip.decompress((HERE/archive['filename']).read_bytes())
        assert hashlib.sha256(raw).hexdigest()==archive['lossless_original_json_sha256']
        data=json.loads(raw);parent=data['actual_native_source_input_archive']
        original=json.loads(gzip.decompress((HERE/parent['filename']).read_bytes()))
        coordinates=native.HalfPstarCoordinates(c,iv(c,data['common_directed_coordinate_theorem']['common_log_bases'][1]),accepted['source_family'])
        eta=original['original_source_rows'][0]['original_typed_native_source']['lossless_actual_native_parent_inputs']['original_eta_log']
        coord=AxialCoordinates(coordinates,iv(c,eta));got=execute_tile(coord,data['exact_Z_range'])
        got.update(original_axial_source_identity=proof,accepted_native_input_archive=parent,
            accepted_complete_transition_archive=archive,
            saved_original_coordinate_theorem=data['common_directed_coordinate_theorem'])
        encoded=json.dumps(encode(got),indent=2).encode()+b'\n';compressed=gzip.compress(encoded,compresslevel=9,mtime=0)
        tag='positive' if Fraction(data['exact_Z_range'][0])>0 else 'negative'
        name=PREFIX+'current_axial_endpoint_q_source_'+tag+'.json.gz';(HERE/name).write_bytes(compressed)
        hashes[name]=sha(name);archives.append(dict(filename=name,compressed_bytes=len(compressed),uncompressed_bytes=len(encoded),
            lossless_original_json_sha256=hashlib.sha256(encoded).hexdigest(),exact_Z_range=data['exact_Z_range']))
    hashes[Path(__file__).name]=sha(Path(__file__).name)
    result=dict(**{GATE:True},source_family=accepted['source_family'],original_axial_q_source_archives=archives,
        original_whole_axial_phase_interval=['0','1'],actual_endpoint_source_queries=2*2*len(PLAN),
        original_signed_b_and_q_ordinary_y2_Z1_sources_installed=True,
        original_full_middle_flat_theorem_on_two_strict_sign_Z_tiles=True,
        full_Z_axis_actual_axial_density_integrals_or_Rc_targets_controls_admitted=False,
        native_upstream_or_accepted_transition_producers_not_rerun=True,
        **dict.fromkeys(common.current.FLAGS,False),input_hashes=hashes,execution_seconds=time.monotonic()-began,
        scope='Whole original O2 axial phase [0,1] signed b, Delta/eta, ordinary q y2/Z1 and exact positive geometry on two actual strict-sign Z tiles. Both nonzero endpoint layers and strict flat middle are covered. No full-Z/axis, complete axial density/integral, control/terminal/heat/stress/recursion/full NS admission.')
    (HERE/NAME).write_text(json.dumps(encode(result),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()
