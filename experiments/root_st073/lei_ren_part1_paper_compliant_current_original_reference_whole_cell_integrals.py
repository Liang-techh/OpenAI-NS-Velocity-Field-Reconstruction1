"""Whole Rh_reference original C0/Z five-density contribution enclosures.

Closed source formulas are evaluated on full exact radial cells. Original
finite-N phase unions and kernels are retained; nonlinear density coefficients
are evaluated before phase union. Canonical factors permit actual cross-cell
Duhamel sums. Incoming corrections and P0 are never reset. This is a coarse
fixed-Z source-window contribution, not functional terminal closure.
"""
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_source_factor_atlas as factors
import lei_ren_part1_paper_compliant_current_original_O2_positive_logq_cells as conditioned

points=factors.points;base=points.base;point=points.point;prior=points.prior
HERE,PREFIX,sha=factors.HERE,factors.PREFIX,factors.sha;ep=points.ep
NAME=PREFIX+'current_original_reference_whole_cell_integrals.json.gz'
RECEIPT=PREFIX+'current_original_reference_whole_cell_integrals_check.json'
GATE='original_Rh_reference_whole_window_actual_N_C0_Z_five_density_contributions_enclosed'


class DirectedCoefficientAlgebra:
    """Interpret the unchanged original coefficient-pair program."""
    def __init__(self,atlas):
        self.atlas=atlas;self.ctx=atlas.ctx;self.zero=atlas.scalar(0);self.one=atlas.scalar(1)
    def constant(self,value):
        q=Fraction(value);return self.atlas.scalar(self.ctx.mpf(q.numerator)/q.denominator)
    def add(self,*values):return self.atlas.sum(values)
    def neg(self,value):return -value
    def mul(self,*values):
        if len(values)==2 and values[0] is values[1]:return base.current.square(values[0])
        result=self.one
        for value in values:result=result*value
        return result
    def quotient(self,numerator,denominator,certificate):
        if ep(denominator.coefficient)[0]<=0:raise ArithmeticError('Original positive denominator required')
        lower=ep(denominator.scale.evaluate()+self.ctx.ln(self.ctx.mpf(ep(denominator.coefficient)[0])))[0]
        return numerator.positive_divide(denominator,lower)
    def unary(self,name,value):
        argument=conditioned.bounded(value)
        if max(abs(x) for x in ep(argument))>1:raise ArithmeticError('Original bounded modulation required')
        if name=='exprel':return self.atlas.scalar(points.directed_exprel(self.ctx,argument))
        if name=='exp':return self.atlas.scalar(self.ctx.exp(argument))
        raise ValueError('Unsupported original coefficient unary')
    c1add=points.exact.source.FunctionTransportGraph.c1add
    c1mul=points.exact.source.FunctionTransportGraph.c1mul
    c1scale=points.exact.source.FunctionTransportGraph.c1scale


def bounded_signed_implicit_B_Z(atlas,kernel,roots,coordinate,chart):
    """Bound the original t*psi_Z product before interval multiplication.

    At fixed true phase psi_Z=-T2_Z/(1+t^2). The rational factor
    t/(1+t^2) belongs to[-1/2,1/2] for every real original t, including
    its native narrow peak. This bounds a function, not a selected cap value.
    The original T1, fixed-angle T1_Z and T2_Z retain their source factors.
    """
    if kernel.geometry!='signed_Mobius' or chart!='E':
        raise ValueError('Strict signed original E chart required')
    if any(not roots[key][(0,1)].zero for key in ('a','b')) or not kernel.t0.zero:
        raise ValueError('Exact original a_Z=b_Z=t0=0 required')
    c=atlas.ctx;scalar=atlas.scalar;add=atlas.add
    q,a,r,rho,s_source,hinv=(kernel.q,kernel.a,kernel.r,kernel.rho,kernel.s_source,kernel.hinv)
    psi_fraction,chi_fraction=kernel.angles(c.mpf(coordinate),chart)
    psi,chi=2*c.pi*psi_fraction,2*c.pi*chi_fraction;difference=chi-psi
    uZ=(roots['p2'][(0,1)]*q).positive_divide(kernel.dstar,kernel.dstar.scale.evaluate())
    K=uZ*hinv;rZ=K*s_source;SZ=-rZ*(2*r)
    sinchi=scalar(c.sin(chi));chiZ=K*(2*c.sin(chi))
    leading=scalar(4*c.cos(chi/2)**2) if kernel.sign>0 else scalar(4*c.sin(chi/2)**2)
    lam=add(add(leading,rho*((-2 if kernel.sign>0 else 2)*c.cos(chi))),-s_source*3)
    H=atlas.sum((scalar(2*chi),s_source*(psi-3*chi),sinchi*(2*r)))
    T1=q*hinv*(difference/r)
    T1Z=q*hinv*(1/r)*add(chiZ,-K*(difference/r))
    T2Z=base.current.square(q)*(1/r**2)*atlas.sum((
        SZ*(psi-3*chi),rZ*sinchi*2,chiZ*lam,-rZ*H*(2/r)))
    implicit_product=-T2Z*c.mpf(['-.5','.5'])
    BZ=-a*add(roots['E'][(0,1)]*T1,kernel.E*add(T1Z,implicit_product))*(1/(4*c.pi))
    return BZ,dict(original_identity='t*psi_Z=-T2_Z*t/(1+t^2)',
        signed_rational_factor_interval=c.mpf(['-.5','.5']),
        exact_nonnegative_certificates=['(t-1)^2>=0','(t+1)^2>=0'],
        T1_fixed_angle=T1.record(),T1_Z_fixed_angle=T1Z.record(),T2_Z_fixed_angle=T2Z.record(),
        source_factors_retained=True,positive_implicit_denominator_lower=1,
        selected_endpoint_midpoint_or_source_cap_value=False)


class OriginalReferenceWholeCells:
    def __init__(self,*,Z):
        self.dispatcher=points.TwoChartOriginalPointSources();self.family=self.dispatcher.source_family
        self.atlas=factors.OriginalSourceFactorAtlas(self.dispatcher.reference.owner.inputs.frame,Z=Z)
        self.ctx=c=self.atlas.ctx;self.Z=self.atlas.Z
        # This implementation currently requires a strict nonmidplane source.
        if self.Z==0:raise ValueError('Whole-window midplane derivative requires a separate signed transition integrator')
        self.inputs=self.dispatcher.reference.owner.inputs;self.scales=self.dispatcher.reference.owner.scales
        self.templates=self.inputs.templates;self.compiled={};self.sensitivity={}
        for key,terms in self.templates['rows'].items():
            self.compiled[key]=tuple((powers,s.lambdify(self.templates['inputs'],expression,modules=[{'mpf':c.mpf},'mpmath']))
                for powers,expression in terms)
            self.sensitivity[key]=tuple(tuple(s.lambdify(self.templates['inputs'],s.diff(expression,P0),
                modules=[{'mpf':c.mpf},'mpmath']) for P0 in self.templates['pressure_symbols']) for powers,expression in terms)
        receipt=json.loads((HERE/conditioned.RECEIPT).read_bytes())
        if not receipt.get('all_passed') or not receipt.get(conditioned.GATE) or receipt['source_family']!=self.family:
            raise ValueError('Accepted original positive-factor phase/primitives kernel required')
        self.hashes={**self.dispatcher.hashes,**self.atlas.hashes}
        for name,digest in {**receipt['input_hashes'],conditioned.RECEIPT:sha(conditioned.RECEIPT)}.items():
            if sha(name)!=digest:raise ValueError('Whole-cell original kernel dependency changed: '+name)
            if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Whole-cell source families disagree')
            self.hashes[name]=digest
        self.hashes[Path(__file__).name]=sha(Path(__file__).name);self.cache={}
    def roots(self,left,right):
        left=points.reference.reference_coordinate(left);right=points.reference.reference_coordinate(right)
        if not left<right:raise ValueError('Strict original radial cell required')
        c=self.ctx;a=self.atlas;lo=a.rational(left);hi=a.rational(right);y=c.mpf((ep(lo)[0],ep(hi)[1]))
        z=a.rational(self.Z);q=1+z*z;f=c.exp(y/10)
        radial=(f,c.mpf(5)/8*f,c.mpf(5)/12*f**2,c.mpf(5)/2*f**2)
        alpha=a.copy_interval(self.inputs.alpha_enclosure)
        pressure=(-alpha/q**2,4*alpha*z/q**3,alpha*(4-20*z*z)/q**4)
        intervals=(z,*radial,*pressure);rows={};source_terms=[]
        for key,terms in self.compiled.items():
            result=a.scalar(0);records=[]
            for index,(powers,fn) in enumerate(terms):
                coefficient=c.mpf(fn(*intervals))
                result=a.add(result,a.term(powers,coefficient,coordinate=y))
                # Original all-late normalized P0 jet budgets have an exact
                # Pstar^-1 factor. Keep it before physical amplification.
                budget=c.mpf(0)
                for order,(sensitivity,bound) in enumerate(zip(self.sensitivity[key][index],(5,10,44),strict=True)):
                    multiplier=abs(c.mpf(sensitivity(*intervals)))
                    budget+=multiplier*c.exp(c.mpf(3)/5)*bound/(2*q**2)
                late=ep(budget)[1]
                if late:
                    r,p,d,ell=powers
                    result=a.add(result,a.term((r,p-1,d,ell),c.mpf((-late,late)),coordinate=y))
                records.append(dict(original_R_Pstar_delta_L_powers=powers,
                    exact_closed_source_coefficient_enclosure=coefficient,
                    normalized_pressure_late_error_Pstar_power=-1 if late else None,
                    positive_late_error_finite_budget_upper=late))
            rows[key]=result;source_terms.append(dict(input=key[0],ordinary_Z_order=key[1],terms=records))
        roots={name:{(0,k):rows[(name,k)] for k in (0,1)} for name in ('E','V','b','p1','p2')}
        roots['a']={(0,0):a.scalar(c.mpf(4)/5),(0,1):a.scalar(0)}
        roots['t0']={(0,0):a.scalar(0),(0,1):a.scalar(0)}
        eta=a.copy_interval(self.scales.logs['eta']);logamin=a.copy_interval(self.scales.logs['a_min'])
        qsource=base.current.q_enclosure(roots['a'][(0,0)],roots['a'][(0,0)]-2,eta,logamin)
        dstar=a.copy_interval(self.scales.logs['d_star'])
        d=a.scalar(1);d=prior.ScaledEnclosure(prior.FormalScale(a.bases,offset=dstar),1,a.ledger)
        u=(roots['p2'][(0,0)]*qsource['q']).positive_divide(d,dstar)
        kernel=conditioned.PositiveLogQPhase(dict(q=qsource['q'],roots=roots,original_u_source=u),dstar)
        if kernel.geometry=='requires_signed_source_refinement':raise ArithmeticError('Whole reference cell needs source sign subdivision')
        return dict(roots=roots,kernel=kernel,source_terms=source_terms,coordinate=y,left=left,right=right)
    def phase_boxes(self,left,right,N):
        phase=self.dispatcher.reference.owner.radius.evaluate(y=left,N=N);c=self.ctx
        width=self.atlas.rational(right-left);shift=c.mpf((0,ep(c.mpf(N)*width)[1]));boxes=[]
        for old in phase['true_original_phase_directed_boxes']:
            cover=base.radius.phase.ordinary_mod_one(c,c.mpf((old['lower'],old['upper']))+shift)
            boxes.extend(cover['boxes'])
        return boxes,phase
    def cell(self,left,right,*,N,bits=32):
        N=points.candidate_N(N);left=points.reference.reference_coordinate(left);right=points.reference.reference_coordinate(right)
        key=(left,right,N,bits)
        if key in self.cache:return self.cache[key]
        c=self.ctx;a=self.atlas;g=DirectedCoefficientAlgebra(a)
        with mp.workdps(c.dps+40):
            query=self.roots(left,right);kernel=query['kernel'];boxes,phase=self.phase_boxes(left,right,N);evaluated=[]
            for phi in boxes:
                if ep(phi)==(0,1):
                    # Original strict monotonicity and periodic endpoint
                    # identities imply the inverse whole-period range.
                    chart='E' if kernel.geometry=='signed_Mobius' else 'psi';coordinate=c.mpf((0,1))
                    inverse=dict(status='enclosed',chart=chart,coordinate_interval=coordinate,
                        actual_source_phase=phi,proof='original strict inverse maps full period[0,1] onto[0,1]')
                else:
                    inverse0=kernel.evaluate(phi,bits=bits)
                    if inverse0['status']!='enclosed':raise ArithmeticError('Refine whole-cell actual phase/source')
                    selected=inverse0['selected_inverse'];chart=selected['chart'];coordinate=selected['coordinate_interval']
                    inverse=dict(inverse0,chart=chart,coordinate_interval=coordinate)
                primitive=kernel.primitives(coordinate,chart)
                jets,derivative=points.slow.slow_values(kernel,query['roots'],coordinate,chart)
                if kernel.geometry=='signed_Mobius' and chart=='E':
                    jets['B_Z_slow'],proof=bounded_signed_implicit_B_Z(a,kernel,query['roots'],coordinate,chart)
                    derivative=dict(derivative,bounded_original_implicit_product=proof)
                pair=lambda name:points.exact.source.C1Function(query['roots'][name][(0,0)],query['roots'][name][(0,1)])
                A=points.exact.source.C1Function(primitive['A'],jets['A_Z_slow'])
                B=points.exact.source.C1Function(primitive['B_over_Pstar'],jets['B_Z_slow'])
                coefficients,F=points.exact.coefficient_pairs(g,pair('E'),pair('V'),A,B,g.constant(N))
                evaluated.append(dict(coefficients=coefficients,inverse=inverse,ordinary_Z=derivative))
            coefficients={order:{key:{jet:points.union([getattr(p['coefficients'][order][key],attr) for p in evaluated])
                for jet,attr in (('C0','value'),('Z','Z'))} for key in points.exact.RATES} for order in points.exact.ORDERS}
        record=dict(source_family=self.family,exact_reference_cell=[str(left),str(right)],candidate_N=N,
            whole_cell_closed_source_terms=query['source_terms'],whole_cell_original_roots={key:{str(k):v.record() for k,v in row.items()} for key,row in query['roots'].items()},
            actual_phase_boxes=boxes,actual_phase_left_origin=phase,
            original_inverse_and_Z_piece_records=[dict(inverse=p['inverse'],ordinary_Z=p['ordinary_Z']) for p in evaluated],
            original_q_positive=not kernel.q.zero,source_geometry=kernel.geometry,
            whole_cell_N_dependent_coefficients={str(o):{k:{j:v.record() for j,v in pair.items()} for k,pair in rows.items()} for o,rows in coefficients.items()},
            closed_source_not_point_sample_extrapolation=True,original_pressure_late_Pstar_factor_collected=True,
            nonlinear_original_coefficient_precedes_phase_union=True,true_radius_phase_Z_independent=True)
        result=dict(coefficients=coefficients,record=record);self.cache[key]=result;return result
    def integrate(self,*,count,N):
        N=points.candidate_N(N)
        if type(count) is not int or count<1 or count>4096:raise ValueError('Exact radial partition count in[1,4096] required')
        c=self.ctx;a=self.atlas;total={o:{key:{j:a.scalar(0) for j in ('C0','Z')} for key in points.exact.RATES} for o in points.exact.ORDERS}
        cells=[]
        with mp.workdps(c.dps+40):
            for i in range(count):
                left=s.Rational(-5)+s.Rational(5*i,count);right=s.Rational(-5)+s.Rational(5*(i+1),count)
                query=self.cell(left,right,N=N);masses={}
                for key,rate in points.exact.RATES.items():
                    rate=c.mpf(rate.numerator)/rate.denominator if isinstance(rate,Fraction) else c.mpf(str(rate))
                    mass=a.rational(right-left) if ep(rate)==(0,0) else (c.exp(rate*a.rational(right))-c.exp(rate*a.rational(left)))/rate
                    if ep(mass)[0]<=0:raise ArithmeticError('Strictly positive original Duhamel cell mass required')
                    masses[key]=mass
                    for order in points.exact.ORDERS:
                        for jet in ('C0','Z'):total[order][key][jet]=a.add(total[order][key][jet],query['coefficients'][order][key][jet]*mass)
                cells.append(dict(source=query['record'],original_positive_own_rate_masses=masses,
                    native_radius_Jacobian=1,physical_kernel='exp(-rate*(0-y)); original reference log-radius coordinate y'))
            full={key:{jet:a.add(total[-1][key][jet]*(c.mpf(1)/N),total[-2][key][jet]*(c.mpf(1)/N**2))
                for jet in ('C0','Z')} for key in points.exact.RATES}
        report=dict(source_family=self.family,original_Z_exact=str(self.Z),candidate_N=N,exact_source_window=['-5','0'],
            exact_original_cells=count,atlas=self.atlas.record(),
            coefficient_own_rate_integral_enclosures={str(o):{k:{j:v.record() for j,v in pair.items()} for k,pair in rows.items()} for o,rows in total.items()},
            actual_finite_N_five_density_contribution_enclosures={k:{j:v.record() for j,v in pair.items()} for k,pair in full.items()},
            full_whole_cell_source_phase_records=cells,ordinary_Z_integral_rows_are_true_source_derivative_enclosures=True,
            original_pressure_rate_zero_memory_retained=True,
            signed_implicit_Z_products_bounded_before_multiplication=True,
            incoming_global_correction_histories_not_assumed_zero=True,original_P0_kept_separate_and_not_reset=True,
            contribution_only_not_terminal_defects=True,all_17_chart_or_24_cell_oracle_installed=False,
            actual_five_controls_installed=False,current_whole_N_selected=False,
            **dict.fromkeys(point.source.inertial.profiles.loop.OPEN,False))
        return dict(report=report,values=full,coefficients=total)


def run():
    begin=time.monotonic();owner=OriginalReferenceWholeCells(Z='.37');integrals=[]
    for count in (4,16):
        got=owner.integrate(count=count,N=160);integrals.append(got['report'])
        print('Whole original reference C0/Z contribution:',count,'cells, N160',flush=True)
    report=dict(**{GATE:True},source_family=owner.family,original_reference_whole_window_refinements=integrals,
        canonical_source_atlas=owner.atlas.record(),actual_whole_reference_C0_Z_density_integrals_enclosed=True,
        exact_original_delta_radius_correlations_collected=True,
        signed_implicit_Z_products_bounded_before_multiplication=True,
        canonical_arithmetic_ledger=dict(owner.atlas.ledger),
        fixed_Z_source_window_only=True,incoming_global_histories_or_P0_not_reset=True,
        all_17_chart_or_24_cell_oracle_installed=False,actual_five_controls_installed=False,
        current_whole_N_selected=False,**dict.fromkeys(point.source.inertial.profiles.loop.OPEN,False),
        input_hashes=owner.hashes,execution_seconds=time.monotonic()-begin,
        scope='Canonical original factors and complete Rh_reference[-5,0] actual candidate-N five-density own-rate contribution enclosures with ordinary Z rows at fixed nonzero Z. Whole cells and actual phase unions, not point quadrature. Coarse local contribution only; no unknown incoming reset, full-Z terminal closure, all-route controls/global-N/recursive corrected field.')
    raw=json.dumps(base.encoded(report),indent=2).encode()+b'\n'
    (HERE/NAME).write_bytes(gzip.compress(raw,compresslevel=9,mtime=0))
    return report


if __name__=='__main__':run()
