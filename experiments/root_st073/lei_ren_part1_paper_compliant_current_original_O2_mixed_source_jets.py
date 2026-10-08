"""Original whole O2 slope mixed source and varying q/nu rows.

Continuous original mass cells supply C0/y/Z/yZ functions, not chosen fields.
q is a positive source log even at the flat profile endpoint y=1.
No mixed phase primitive or terminal controller is claimed by this module.
"""
from dataclasses import dataclass
import gzip
import json
from pathlib import Path
from types import MappingProxyType
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_O2_positive_logq_cells as positive
import lei_ren_part1_paper_compliant_current_original_reference_axial_predicate_C1 as axial

ordered=positive.ordered;base=positive.base;prior=base.prior;ep=positive.ep
HERE,PREFIX,sha=positive.HERE,positive.PREFIX,positive.sha
C0,Y,Z,YZ=axial.mixed.ORDERS;MixedJet=axial.MixedJet
NAME=PREFIX+'current_original_O2_mixed_source_jets.json.gz'
RECEIPT=PREFIX+'current_original_O2_mixed_source_jets_check.json'
GATE='original_whole_O2_slope_four_source_jets_and_nonconstant_q_nu_bound'


class O2MixedAtlas(axial.AxialSourceAtlas):
    """Collected original radius; reserved slot3 is defined source logq."""
    def __init__(self,frame,*,lower,upper,logq):
        lo,hi=(ordered.base.point.pressure.exact_Z(v) for v in (lower,upper))
        if not -1<=lo<=hi<=1:raise ValueError('Closed original O2 axial interval required')
        axial.whole.factors.OriginalSourceFactorAtlas.__init__(self,frame,Z=0)
        c=self.ctx
        with mp.workdps(c.dps+40):
            self.Z=None;self.bounds=(lo,hi);self.carrier=(0,0)
            z=c.mpf((ep(self.rational(lo))[0],ep(self.rational(hi))[1]))
            self.Q=1+z**2
            self.L=1-axial.whole.conditioned.bounded(self.parameter('delta'))*z**2
            self.bases=self.bases[:2]+(c.ln(self.L),c.mpf(ep(logq)),c.mpf(0))
            if ep(self.L)[0]<=0 or ep(self.Q)[0]<1 or ep(self.bases[3])[1]>0:
                raise ValueError('Positive original Q/L and q<=1 source ranges required')
            if ep(self.bases[0])[0]<=0 or ep(self.bases[1])[0]<=0 or ep(self.bases[2])[1]>0:
                raise ValueError('Original positive P/C and nonpositive logL basis required')

    def add(self,left,right):
        left.scale.pair(right.scale)
        if left.ledger is not self.ledger or right.ledger is not self.ledger:
            raise ValueError('One O2 source ledger required')
        if left.zero:return right
        if right.zero:return left
        lp,rp=left.scale.powers,right.scale.powers
        if lp[4]!=0 or rp[4]!=0:raise ValueError('Original O2 radius must already be collected')
        powers=(max(lp[0],rp[0]),max(lp[1],rp[1]),min(lp[2],rp[2]),min(lp[3],rp[3]),0)
        anchor=prior.FormalScale(self.bases,powers)
        values=[v.coefficient*v.bounded_exp((v.scale-anchor).evaluate()) for v in (left,right)]
        self.ledger['O2_collected_P_C_L_q_additions']=self.ledger.get('O2_collected_P_C_L_q_additions',0)+1
        return prior.ScaledEnclosure(anchor,values[0]+values[1],self.ledger)

    def record(self):
        return dict(source_family=self.family,exact_Z_bounds=[str(v) for v in self.bounds],
            source_basis_order=['logPstar','selected_logCstar','logL','original_logq','zero'],
            defining_basis=self.bases,original_Q=self.Q,original_L=self.L,
            exact_source_radius_collected_before_native_log_evaluation=True,
            original_q_is_positive_log_function_not_selected_field_value=True)


@dataclass(frozen=True)
class OriginalO2MixedFrame:
    owner:object
    family:object
    count:int
    index:int
    left:object
    right:object
    roots:object
    record:dict


def radial_source_derivatives(templates):
    z,f,H,D,P,p0,p0Z,p0ZZ=templates['inputs'];a=s.Symbol('original_a',real=True)
    rates={f:(1-a)*f/2,H:f-s.Rational(3,2)*H,D:f*f/2-D,P:f*f/2}
    rows={}
    for key,terms in templates['rows'].items():
        rows[key]=tuple((powers,powers[0]*expr+sum(s.diff(expr,x)*rate for x,rate in rates.items()))
            for powers,expr in terms)
    return dict(templates,rows=rows),a,dict(passed=True,
        exact_original_radial_ODEs=['f_y=(1-a)*f/2','H_y=f-3H/2','D_y=f^2/2-D','P_y=f^2/2'],
        original_pressure_datum_y_and_yZ_exact_zero=True,
        native_radius_factor_y_derivative_retained=True,
        original_Z_rows_differentiated_radially_no_second_L_Z_applied=True)


class OriginalO2MixedSources:
    def __init__(self):
        self.parent=positive.OriginalO2PositiveLogQCells()
        self.owner=self.parent.owner;self.family=self.parent.family;self.ctx=self.owner.ctx
        accepted=json.loads((HERE/positive.RECEIPT).read_bytes())
        if not accepted.get('all_passed') or not accepted.get(positive.GATE):
            raise ValueError('Accepted original ordered positive-q cell source required')
        if accepted['source_family']!=self.family:raise ValueError('O2 mixed source family differs')
        self.hashes=dict(self.parent.hashes)
        for name,digest in {**accepted['input_hashes'],positive.RECEIPT:sha(positive.RECEIPT)}.items():
            if sha(name)!=digest:raise ValueError('Original O2 source dependency changed: '+name)
            if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Original O2 hash closures disagree')
            self.hashes[name]=digest
        for module in (axial,):
            receipt=json.loads((HERE/module.RECEIPT).read_bytes())
            if not receipt.get('all_passed') or not receipt.get(module.GATE) or receipt['source_family']!=self.family:
                raise ValueError('Same-family accepted canonical mixed arithmetic required')
            for name,digest in {**receipt['input_hashes'],module.RECEIPT:sha(module.RECEIPT)}.items():
                if sha(name)!=digest:raise ValueError('Accepted mixed arithmetic dependency changed: '+name)
                if name in self.hashes and self.hashes[name]!=digest:
                    raise ValueError('Mixed arithmetic and O2 source closures disagree')
                self.hashes[name]=digest
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.original=self.owner.inputs.templates
        derivative,self.a_symbol,self.derivative_proof=radial_source_derivatives(self.original)
        self.derivative_templates=derivative['rows']
        drows,symbols,self.pressure_correlation=ordered.pressure_correlation(derivative)
        self.rows={C0:self.parent.parent.rows,Y:drows}
        self.symbols={C0:self.parent.parent.symbols,Y:(*symbols,self.a_symbol)}
        self.compiled={}
        c=self.ctx
        for direction,rows in self.rows.items():
            self.compiled[direction]={key:tuple((powers,
                s.lambdify(self.symbols[direction],expr,modules=[{'mpf':c.mpf},'mpmath']),
                tuple(s.lambdify(self.symbols[direction],sensitivity,modules=[{'mpf':c.mpf},'mpmath'])
                    for sensitivity in sensitivities)) for powers,expr,sensitivities in terms)
                for key,terms in rows.items()}
        self.frames={};self.cache={}

    def source_frame(self,count,index,*,Z_lower,Z_upper):
        if count not in self.parent.parent.levels or type(index) is not int or not 0<=index<count:
            raise ValueError('Accepted original source level/cell required')
        lower=ordered.base.point.pressure.exact_Z(Z_lower);upper=ordered.base.point.pressure.exact_Z(Z_upper)
        key=(count,index,lower,upper)
        if key in self.cache:return self.cache[key]
        c=self.ctx
        with mp.workdps(c.dps+40):
            level,_=self.parent.parent.levels[count];cell=level['whole_source_cells'][index]
            left,right=s.Rational(index,count),s.Rational(index+1,count)
            eta=ordered.interval(c,self.owner.scales.logs['eta'])
            loglo=positive.endpoint_logq(c,right,eta);loghi=positive.endpoint_logq(c,left,eta)
            logq=ordered.hull(c,loglo,loghi)
            a=O2MixedAtlas(self.owner.inputs.frame,lower=lower,upper=upper,logq=logq);ctx=a.ctx
            y=ctx.mpf((ep(a.rational(left))[0],ep(a.rational(right))[1]))
            z=ctx.mpf((ep(a.rational(lower))[0],ep(a.rational(upper))[1]))
            profiles={name:a.copy_interval(value) for name,value in cell['whole_original_radial_profile_covers'].items()}
            values=(z,*(profiles[name] for name in ('f','H','D','P','remaining_pressure_mass')))
            rows={};records=[]
            for direction,compiled in self.compiled.items():
                arguments=values if direction==C0 else (*values,profiles['a'])
                for (name,order),terms in compiled.items():
                    value=a.scalar(0);source_terms=[]
                    for powers,fn,sensitivities in terms:
                        coefficient=ctx.mpf(fn(*arguments))
                        native=axial.regular.finite_offset_anchor(a,a.term(powers,coefficient,coordinate=y))
                        value=a.add(value,native);errors=[]
                        for pressure_order,(sensitivity,budget) in enumerate(zip(sensitivities,(5,10,44),strict=True)):
                            sensitivity=ctx.mpf(sensitivity(*arguments))
                            if lower==upper==0 and pressure_order==1:continue
                            error=abs(sensitivity)*ctx.exp(ctx.mpf(3)/5)*budget/(2*a.Q*a.Q)
                            if ep(error)[1]:
                                r,p,d,ell=powers;bound=ep(error)[1]
                                tail=axial.regular.finite_offset_anchor(a,a.term((r,p-1,d,ell),
                                    ctx.mpf((-bound,bound)),coordinate=y))
                                value=a.add(value,tail);errors.append(dict(pressure_order=pressure_order,
                                    finite_sensitivity=sensitivity,full_original_error_source=tail.record()))
                        source_terms.append(dict(original_factor_powers=powers,finite_original_coefficient=coefficient,
                            original_term=native.record(),full_pressure_error_terms=errors))
                    rows.setdefault(name,{})[(direction[0],order)]=value
                    records.append(dict(input=name,ordinary_y_order=direction[0],ordinary_Z_order=order,terms=source_terms))
            roots={name:MixedJet(a,row) for name,row in rows.items()}
            # Same original flat cutoff, actual ordinary derivative, not a
            # derivative of the saved profile hull. Its monotonicity is exact.
            sigma=prior.sigma_jets(ctx,y)
            ay_box=ctx.mpf(6)/5*sigma[1]
            ay_box=axial.whole.conditioned.base.conditioned.clipped(ctx,ay_box,0,ep(ay_box)[1])
            ay=a.scalar(ay_box);av=a.scalar(profiles['a'])
            roots['a']=MixedJet(a,{C0:av,Y:ay,Z:a.scalar(0),YZ:a.scalar(0)})
            roots['b']=MixedJet.constant(a,a.scalar(0))
            roots['t0']=MixedJet.constant(a,a.scalar(0))
            q=prior.ScaledEnclosure(prior.FormalScale(a.bases,(0,0,0,1,0)),1,a.ledger)
            eta_source=prior.ScaledEnclosure(prior.FormalScale(a.bases,offset=a.copy_interval(eta)),1,a.ledger)
            one_plus_eta=a.add(a.scalar(1),eta_source)
            a2=base.current.square(av);lower_a=a.ctx.ln(ctx.mpf(4)/5)
            qy=-(one_plus_eta*ay).positive_divide(a2*q*2,lower_a*2+logq+ctx.ln(2))
            # This exact cancellation avoids differentiating a rounded q^2.
            nu0=(one_plus_eta*2).positive_divide(av,lower_a)
            nuy=-(one_plus_eta*ay*2).positive_divide(a2,lower_a*2)
            roots['q']=MixedJet(a,{C0:q,Y:qy,Z:a.scalar(0),YZ:a.scalar(0)})
            roots['nu']=MixedJet(a,{C0:nu0,Y:nuy,Z:a.scalar(0),YZ:a.scalar(0)})
            record=dict(source_family=self.family,source_level=count,source_index=index,
                exact_y_cell=[str(left),str(right)],exact_Z_range=[str(lower),str(upper)],
                actual_original_mixed_source_terms=records,actual_root_jets={name:row.record() for name,row in roots.items()},
                original_radial_derivative_recipe=self.derivative_proof,
                pressure_prefix_suffix_function_correlation=self.pressure_correlation,
                original_a_derivative_source='a_y=(6/5)*sigma_y; sigma_jets original flat source, ordinary first derivative',
                original_q_definition='q^2=(2+2*eta-a)/(2*a); a<=2 implies original active cutoff sigma=1',
                actual_qy_definition='q_y=-(1+eta)*a_y/(2*a^2*q)',
                actual_nu_definitions=['nu=1+2*q^2=2*(1+eta)/a','nu_y=-2*(1+eta)*a_y/a^2'],
                q_and_nu_Z_yZ_exact_zero_only_on_original_O2_slope=True,
                original_b_and_t0_exact_zero_only_on_original_O2_slope=True,
                original_positive_eta_source=eta_source.record(),native_positive_logq=logq,
                q_at_endpoint_one_positive_source_not_exact_flat=True,
                actual_parameter_derivatives_not_reset_to_reference_constants=True,
                original_P0_datum_sha256=self.family['datum_enclosure_sha256'],
                ordinary_Z_yZ_not_scaled_coordinate_derivatives=True,atlas=a.record(),
                separate_source_hulls_not_exact_selected_compatible_field=True,
                mixed_phase_inverse_or_primitives_installed=False,actual_five_controls_installed=False)
            result=OriginalO2MixedFrame(self,self.family,count,index,left,right,MappingProxyType(roots),record)
            self.frames[id(result)]=result;self.cache[key]=result
            return result

    def describe(self,frame):
        if type(frame) is not OriginalO2MixedFrame or frame.owner is not self or self.frames.get(id(frame)) is not frame:
            raise ValueError('Issued same original O2 mixed source frame required')
        return frame.record


def run():
    begin=time.monotonic();owner=OriginalO2MixedSources()
    count=min(owner.parent.parent.levels);records=[]
    for index in range(count):
        records.append(owner.describe(owner.source_frame(count,index,Z_lower=-1,Z_upper=1)))
    extras=[]
    for index,zl,zh in ((0,0,0),(count//2,0,0),(count-1,0,0),
        (count//2,s.Rational(36,100),s.Rational(38,100)),
        (count-1,s.Rational(-38,100),s.Rational(-36,100))):
        extras.append(owner.describe(owner.source_frame(count,index,Z_lower=zl,Z_upper=zh)))
    report=dict(**{GATE:True},source_family=owner.family,ordered_source_cells=count,
        exact_original_O2_slope_domain=dict(y=['0','1'],Z=['-1','1']),
        actual_whole_O2_mixed_source_cells=records,actual_axis_and_signed_source_cells=extras,
        actual_original_a_q_nu_nonconstant_y_rows_installed=True,
        original_pressure_and_incoming_mass_source_unchanged=True,
        whole_continuous_source_domains_not_native_samples=True,
        mixed_phase_inverse_or_primitives_installed=False,actual_changed_five_integrals_installed=False,
        all_17_chart_or_24_cell_oracle_installed=False,actual_five_controls_installed=False,current_whole_N_selected=False,
        **dict.fromkeys(ordered.base.point.source.inertial.profiles.loop.OPEN,False),
        input_hashes=owner.hashes,execution_seconds=time.monotonic()-begin,
        scope='Complete original O2 slope continuous y[0,1]/Z[-1,1] C0/y/Z/yZ roots and actual varying a,q,nu source rows, full pressure prefix/suffix and late errors, positive q at y1. Source interface only; next mixed phase/Fourier/Mobius/implicit products and density transport remain open. Not terminal/all-route/controls/global N/recursion/corrected NS.')
    (HERE/NAME).write_bytes(gzip.compress(json.dumps(positive.ordered.base.encoded(report),indent=2).encode()+b'\n',
        compresslevel=9,mtime=0))
    print('Whole original O2 mixed source cells:',count,'plus5 axis/signed cells',flush=True)
    return report


if __name__=='__main__':run()
