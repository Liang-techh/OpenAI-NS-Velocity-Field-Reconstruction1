"""Original frozen-macro F/V integrals with signed nonlinear transport.

The known comparison coefficients are recovered from defining width-polynomial
source jets. Complete exponential integrals use exact I/J masses and a
nonzero third-width-order remainder. Native micro inlet uncertainty remains.
"""
import json
import math
from pathlib import Path
from types import SimpleNamespace
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_original_reshape_terminal_kernels as previous
import lei_ren_part1_paper_compliant_actual_bridge_integrals as bridge
import lei_ren_part1_paper_compliant_macro_signed_integrals as modes
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_logarithmic_pressure_datum import LogarithmicPressureDatum
from lei_ren_part1_paper_compliant_inner_bridge_profiles import logarithm

HERE,PREFIX,sha=previous.HERE,previous.PREFIX,previous.sha
prior,base,ep=previous.prior,previous.base,previous.ep
NAME=PREFIX+'current_original_bridge_macro_functions.json'
RECEIPT=PREFIX+'current_original_bridge_macro_functions_check.json'
GATE='original_frozen_macro_complete_nonlinear_FV_integral_functions_installed'
PARTS=('hydro','pressure','swirl')


def magnitude(c,value):
    return c.mpf(max(abs(x) for x in ep(value)))


class MacroFlow:
    """Truncated ordinary Z algebra; all radial endpoints are Z-independent."""
    def __init__(self,c,logh,logP2,logF02,logRa,Y,weight):
        self.c=c;self.weight=c.mpf(weight);self.Y=c.mpf(Y)
        self.logs=tuple(c.mpf(v) for v in (logh,logP2,logF02,logRa,0))
        exact_Y=c.ln(100)-self.logs[3]
        if self.Y._mpi_!=exact_Y._mpi_:raise ValueError('Original Y must be the same ln(100)-logRa source expression')
        if ep(self.weight)[0]<=0 or ep(self.Y)[0]<=0:
            raise ValueError('Positive fixed axial weight and macro logarithmic length required')
        if any(not mp.isfinite(x) for row in self.logs+(self.Y,self.weight) for x in ep(row)):
            raise ValueError('Finite fixed source log records required')
        self.ledger=dict(directed_small_exponential_tails=0,positive_function_denominator_intersections=0,
            positive_function_root_intersections=0,directed_independent_log_rescalings=0)
        self.scalar=lambda v:prior.ScaledEnclosure(prior.FormalScale(self.logs),v,self.ledger)
        self.h=self.factor((1,0,0,0,0))
        self.hupper=self.small_upper(self.h)
        if ep(self.hupper)[1]>=ep(self.Y/4)[0]:raise ValueError('Original2hb smoothing precedes macro')
        self.sources_set=False

    def factor(self,powers,offset=0):
        return prior.ScaledEnclosure(prior.FormalScale(self.logs,powers,offset),1,self.ledger)

    def small_upper(self,row):
        if row.zero:return self.c.mpf(0)
        log=row.scale.evaluate()+self.c.ln(magnitude(self.c,row.coefficient))
        hi=ep(log)[1]
        if hi>mp.log(2):raise ValueError('Exponential perturbation is not admitted small')
        cut=-2*self.c.dps*self.c.ln(10)
        # An outward tail bound only; never a source value or selected width.
        return self.c.exp(self.c.mpf(max(hi,ep(cut)[1])))

    def ordinary_cover(self,row):
        if row.zero:return self.c.mpf(0)
        log=row.scale.evaluate();lo,hi=ep(log);cut=-2*self.c.dps*self.c.ln(10)
        if hi<ep(cut)[0]:
            self.ledger['directed_small_exponential_tails']+=1
            return row.coefficient*self.c.mpf([0,ep(self.c.exp(cut))[1]])
        return row.coefficient*self.c.exp(log)

    def jet(self,value):
        if not isinstance(value,IntervalTaylor) or value.ctx is not self.c or value.order!=5:
            raise ValueError('Same-context ordinary source Taylor jets0..5 required')
        if any(not mp.isfinite(x) for row in value.coefficients for x in ep(row)):
            raise ValueError('Finite source jet coefficients required')
        return [self.scalar(v) for v in value.coefficients]

    def add(self,*rows):
        return [sum((row[n] for row in rows),self.scalar(0)) for n in range(6)]

    def multiply(self,a,b):
        return [sum((a[j]*b[n-j] for j in range(n+1)),self.scalar(0)) for n in range(6)]

    def scale(self,rows,value):return [v*value for v in rows]

    def norm(self,rows):
        return sum((prior.ScaledEnclosure(v.scale,magnitude(self.c,v.coefficient),self.ledger)*self.weight**n
                    for n,v in enumerate(rows)),self.scalar(0))

    def error_rows(self,bound):
        return [prior.ScaledEnclosure(bound.scale,bound.coefficient*self.c.mpf([-1,1])/self.weight**n,self.ledger)
                for n in range(6)]

    def term(self,rows,width_power,part=None):
        powers=[width_power,0,0,0,0]
        if part=='pressure':powers[1]=1
        elif part=='swirl':powers[2]=1
        elif part not in (None,'hydro'):raise ValueError('Original source part required')
        return self.scale(self.jet(rows),self.factor(tuple(powers)))

    def set_sources(self,directions,drives,quotient,ell_in,phi0,V0,deltaV_in):
        if len(directions)!=3 or set(drives)!=set(PARTS) or any(len(v)!=3 for v in drives.values()):
            raise ValueError('Original three fixed exponential comparison modes required')
        self.d=[self.jet(v) for v in directions]
        self.drive={key:[self.jet(v) for v in value] for key,value in drives.items()}
        self.q=self.jet(quotient);self.phi0=self.jet(phi0);self.V0=self.jet(V0)
        for row in (ell_in,deltaV_in):
            if len(row)!=6 or any(v.scale.bases is not self.logs or v.ledger is not self.ledger for v in row):
                raise ValueError('One source basis/ledger for inlet Taylor rows required')
        self.ell_in=ell_in;self.deltaV_in=deltaV_in;self.sources_set=True

    def geometry(self,fraction):
        if not isinstance(fraction,tuple) or len(fraction)!=2 or any(type(v) is not int for v in fraction):
            raise ValueError('An exact rational macro fraction pair is required')
        n,d=fraction
        if d<=0 or not 0<=n<=d:raise ValueError('Macro fraction must lie in[0,1]')
        q=self.c.mpf(n)/d
        S=q*(self.Y-2*self.c.mpf([0,ep(self.hupper)[1]]))
        R0=self.factor((0,0,0,1,0),self.c.mpf([0,ep(2*self.hupper)[1]]))
        if n==d:R1=self.scalar(100)
        elif n==0:R1=R0
        else:R1=self.factor((0,0,0,1,0),q*self.Y+self.c.mpf([0,ep(2*self.hupper*(1-q))[1]]))
        return S,R0,R1,n==0

    def radial_power(self,row,power):
        return prior.ScaledEnclosure(prior.FormalScale(self.logs,
            tuple(v*power for v in row.scale.powers),row.scale.offset*power),row.coefficient**power,self.ledger)

    def I(self,p,j,S,R0,R1,empty=False):
        if empty:return self.scalar(0)
        a=p-j
        if a==0:return self.radial_power(R0,p)*S
        return (self.radial_power(R0,j)*self.radial_power(R1,a)-self.radial_power(R0,p))*(self.c.mpf(1)/a)

    def J(self,p,j,l,S,R0,R1,empty=False):
        if empty:return self.scalar(0)
        if l!=1:
            return (self.I(p+1,j+l,S,R0,R1)-R0*self.I(p,j,S,R0,R1))*(self.c.mpf(1)/(1-l))
        a=p-j
        if a==0:return self.radial_power(R0,p+1)*(S*S/2)
        return (self.radial_power(R0,j+1)*self.radial_power(R1,a)*(a*S-1)+self.radial_power(R0,p+1))*(self.c.mpf(1)/(a*a))

    def evaluate(self,fraction):
        if not self.sources_set:raise ValueError('Defining source jets must be attached')
        c=self.c;S,R0,R1,empty=self.geometry(fraction)
        # Positive original radial masses also bound every shorter prefix.
        I={(p,j):self.I(p,j,S,R0,R1,empty) for p in (1,2) for j in range(3)}
        ell=self.add(self.ell_in,*[self.scale(row,self.h*I[1,l]*(-c.mpf('.5')))
                                  for l,row in enumerate(self.d)])
        prefix=self.norm(self.ell_in)+sum((self.norm(row)*self.h*prior.ScaledEnclosure(
            I[1,l].scale,magnitude(c,I[1,l].coefficient),self.ledger)*c.mpf('.5')
            for l,row in enumerate(self.d)),self.scalar(0))
        upper=self.small_upper(prefix)
        if ep(upper)[1]>c.mpf('.5'):raise ValueError('Full prefix weighted log jet norm must be <=1/2')
        remainder=prefix*prefix*(c.exp(upper)/2)
        delta_phi=self.multiply(self.phi0,self.add(ell,self.error_rows(remainder)))
        kernels={};parts={};errors={}
        for part,p in (('hydro',1),('pressure',1),('swirl',2)):
            sums=[];err=self.scalar(0);kernel_rows=[]
            for j in range(3):
                unit=[I[p,j]]+[self.scalar(0)]*5
                linear=self.add(unit,self.scale(self.ell_in,I[p,j]),
                    *[self.scale(row,self.h*self.J(p,j,l,S,R0,R1,empty)*(-c.mpf('.5')))
                      for l,row in enumerate(self.d)])
                mass=prior.ScaledEnclosure(I[p,j].scale,magnitude(c,I[p,j].coefficient),self.ledger)
                kernel_error=remainder*mass
                kernel=self.add(linear,self.error_rows(kernel_error))
                kernel_rows.append(dict(coefficients=kernel,signed_first_exponential_terms=linear,
                                        full_exponential_absolute_error=kernel_error))
                product=self.multiply(self.q,self.drive[part][j])
                sums.append(self.multiply(product,kernel))
                err+=self.norm(product)*kernel_error
            scale=self.factor((1,int(part=='pressure'),int(part=='swirl'),0,0))
            parts[part]=self.scale(self.add(*sums),-scale)
            errors[part]=err*scale
            kernels[part]=kernel_rows
        deltaV=self.add(self.deltaV_in,*parts.values())
        return dict(ell=ell,delta_phi=delta_phi,phi=self.add(self.phi0,delta_phi),
            deltaV=deltaV,V=self.add(self.V0,deltaV),kernels=kernels,parts=parts,
            macro_absolute_error_norms=errors,full_prefix_log_jet_norm=prefix,
            exponential_remainder_norm=remainder,geometry=dict(fraction=list(fraction),S=S,
                R0=R0,R1=R1,fixed_radial_endpoints_Z_independent=True),
            source_contract='q=phi_core_exit/barphi2; ell includes micro angular prefix exactly once',
            ordinary_Z_orders=list(range(6)),full_original_macro_exponential_integral_enclosed=True)


def serialized(value):
    if isinstance(value,prior.ScaledEnclosure):return value.record()
    if isinstance(value,dict):return {k:serialized(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [serialized(v) for v in value]
    return value


class OriginalBridgeMacroFunctions:
    mode='genuine_original_frozen_macro_FV_functions_with_source_micro_inlet_errors'
    def __init__(self,dps=500):
        if type(dps) is not int or dps<500:raise ValueError('At least500 source digits required')
        self.c=c=MPIntervalContext();c.dps=dps;self.hashes={};self.records={}
        for stem in ('actual_bridge_integrals','comparison_point_integrals','anchored_axis_amplitude'):
            name=PREFIX+stem+'_check.json';checked=json.loads((HERE/name).read_bytes())
            if not checked.get('all_passed'):raise ValueError('Accepted defining source required: '+stem)
            for path,digest in checked['input_hashes'].items():previous.bind(self.hashes,path,digest)
            previous.bind(self.hashes,name,sha(name))
            name=PREFIX+stem+'.json';self.records[stem]=json.loads((HERE/name).read_bytes())
            previous.bind(self.hashes,name,sha(name))
        original=self.records['actual_bridge_integrals']
        self.family=original['actual_five_defect_family_sha256'];self.source=original['implicit_source_sha256'];self.datum=original['datum_enclosure_sha256']
        for row in self.records.values():
            for key,want in (('actual_five_defect_family_sha256',self.family),('implicit_source_sha256',self.source),('datum_enclosure_sha256',self.datum)):
                if row[key]!=want:raise ValueError('One original family/source/P0 required')
        for stem in ('pressure_source','core_transfer'):
            name=PREFIX+stem+'.json';self.records[stem]=json.loads((HERE/name).read_bytes());previous.bind(self.hashes,name,sha(name))
        self.bindings=bridge.bridge_control_bindings()
        # Replay the original pressure jet method from its admitted fourteen-atom
        # cache, without constructing or recomputing outer waiting/pressure data.
        ps=self.records['pressure_source']['compliant_source'];datum=LogarithmicPressureDatum.__new__(LogarithmicPressureDatum)
        datum.ctx=c;datum.m2=previous.read_interval(c,ps['fixed_beta2_mass_normalized']);datum.m0=previous.read_interval(c,ps['fixed_beta0_mass_normalized'])
        datum.rho=previous.read_interval(c,ps['complex_strip_half_width']);datum.flatten_complex_upper=previous.read_interval(c,ps['flatten_complex_mass_upper'])
        datum.stages={'z_flatten':{'mass':previous.read_interval(c,ps['stages']['z_flatten']['mass'])}}
        datum.parameters=SimpleNamespace(logPstar=previous.read_interval(c,ps['parameter_bounds']['logPstar']))
        datum.source_sha=self.source;datum.datum_sha=self.datum;datum.input_hashes=ps['input_hashes'];self.pressure=datum
        if ps['implicit_source_sha256']!=self.source or ps['datum_enclosure_sha256']!=self.datum:
            raise ValueError('Original pressure cache source identity required')
        self.owners={}
        for module in (previous,bridge,modes):previous.bind(self.hashes,Path(module.__file__).name,sha(Path(module.__file__).name))
        previous.bind(self.hashes,Path(__file__).name,sha(Path(__file__).name))

    def owner(self,label):
        if label not in ('0','.5'):raise ValueError('Admitted comparison/amplitude source frames0,.5 only')
        if label in self.owners:return self.owners[label]
        c=self.c;original=self.records['actual_bridge_integrals'];comp=self.records['comparison_point_integrals']
        frame=original['packets'][label];entry=comp['comparison_point_packets'][label]['macro'][0]
        amplitude=self.records['anchored_axis_amplitude']['anchored_amplitude_packets'][label]
        Z=previous.read_interval(c,entry['Z'])
        if Z._mpi_!=previous.read_interval(c,frame['core']['Z'])._mpi_ or Z._mpi_!=previous.read_interval(c,amplitude['Z'])._mpi_:
            raise ValueError('Exact same original source axial coordinate required')
        core=self.records['core_transfer'];logLambda=previous.read_interval(c,core['logLambda']);logRa=c.ln(4)-logLambda
        logh=previous.read_interval(c,original['source_log_hb_enclosure']);logF02=2*previous.read_interval(c,amplitude['logF0'])
        logP2=2*previous.read_interval(c,core['logPstar']);weight=previous.read_interval(c,original['axial_jet_weight'])
        if ep(previous.read_interval(c,entry['coordinate']))!=(0,0):raise ValueError('Original frozen comparison phase0 required')
        if previous.read_interval(c,entry['axial_jet_weight'])._mpi_!=weight._mpi_:
            raise ValueError('One original source weighted Taylor norm required')
        flow=MacroFlow(c,logh,logP2,logF02,logRa,c.ln(100)-logRa,weight)
        def jet(rows):return IntervalTaylor(c,[previous.read_interval(c,row) for row in rows])
        def poly(packet):
            coefficients=[jet(row) for row in packet['signed_hb_power_axial_coefficients']]
            if any(v.order!=6 for v in coefficients):raise ValueError('Original comparison axial6 source required')
            row=[sum((flow.scalar(v[k])*flow.factor((n,0,0,0,0)) for n,v in enumerate(coefficients)),flow.scalar(0)) for k in range(7)]
            remainder=previous.read_interval(c,packet['omitted_hb_cubed_weighted_axial_jet_norm_upper'])
            bound=flow.h*flow.h*flow.h*remainder
            row=[v+prior.ScaledEnclosure(bound.scale,bound.coefficient*c.mpf([-1,1])/flow.weight**k,flow.ledger) for k,v in enumerate(row)]
            values=[flow.ordinary_cover(v) for v in row]
            return IntervalTaylor(c,values),dict(signed_width_polynomial=packet,
                directed_nonzero_remainder_norm=bound.record(),materialized_as_source_enclosure_only=True)
        phi,phi_proof=poly(entry['comparison_phi']);V,V_proof=poly(entry['comparison_raw_V'])
        moment_rows={};moment_proofs={}
        for name,packet in entry['comparison_own_six_moments'].items():moment_rows[name],moment_proofs[name]=poly(packet)
        Lambda=previous.read_interval(c,core['Lambda']);gradient=[previous.read_interval(c,v) for v in amplitude['G_gradient_taylor_coefficients']]
        def ratios(multiplier):
            rows=[c.mpf(1)]
            for n in range(1,7):rows.append(sum((-multiplier*Lambda*gradient[j]*rows[n-1-j] for j in range(n)),c.mpf(0))/n)
            return [v*math.factorial(n) for n,v in enumerate(rows)]
        p0=IntervalTaylor(c,[c.mpf(ep(v)) for v in self.pressure.normalized_jets(ep(Z),6)['normalized_pressure_coefficients']])
        delta=c.exp(previous.read_interval(c,self.records['pressure_source']['compliant_source']['parameter_bounds']['log_delta']))
        inputs=dict(p0=p0,F0_ratios=ratios(1),F0_squared_ratios=ratios(2))
        known=modes._mode_rows(SimpleNamespace(ctx=c,delta=delta),Z,inputs,phi,V,bridge.named_moments(moment_rows))
        phi0=jet(frame['core']['actual_phi_axial5']);V0=jet(frame['core']['actual_raw_V_axial5'])
        def inlet_terms(terms,part=None):
            rows=[]
            for term in terms:
                rows.append(flow.term(jet(term['normalized_axial_coefficients']),term['width_power'],part))
                norm=previous.read_interval(c,term['nonlinear_prefix_error_weighted_norm_per_next_hb'])
                powers=[term['width_power']+1,int(part=='pressure'),int(part=='swirl'),0,0]
                rows.append(flow.error_rows(flow.factor(tuple(powers))*norm))
            return flow.add(*rows)
        inlet=frame['second_exit'];ell=inlet_terms(inlet['signed_angular_integral_terms'])
        dV=flow.add(*[inlet_terms(inlet['signed_axial_integral_terms'][part],part) for part in PARTS])
        flow.set_sources([v.truncate(5) for v in known['D_over_R']],
            {part:[v.truncate(5) for v in known['drive_'+part]] for part in PARTS},
            phi0/phi.truncate(5),ell,phi0,V0,dV)
        proof=dict(Z=Z,source_family=self.family,implicit_source_sha256=self.source,datum_enclosure_sha256=self.datum,
            original_comparison_source_namespace=comp['original_comparison_source']['source_namespace'],
            original_source_bindings=self.bindings,comparison_phi=phi_proof,comparison_V=V_proof,
            comparison_moments=moment_proofs,actual_anchored_G=amplitude['G'],
            log_F0_squared=logF02,original_gradient_coefficients=gradient,
            actual_G_not_Gbar_used=True,original_P0_coefficients=list(p0.coefficients),
            original_mode_rows_recomputed_from_defining_finite_width_comparison=True,
            normalized_swirl_rows_include_true_F0_squared_derivatives_once=True,
            quotient='Phi_core_exit/barphi2; full micro log prefix enters exponential separately',
            original_micro_inlet_errors_retained=True,source_width_not_selected=True,
            no_original_ancestor_constructors_or_producers_executed=True)
        self.owners[label]=(flow,proof);return flow,proof

    def evaluate(self,label,fraction):
        flow,proof=self.owner(label)
        with mp.workdps(self.c.dps+40):value=flow.evaluate(fraction)
        return dict(mode=self.mode,source_proof=proof,function_evaluation=serialized(value),
            actual_micro_function_provider_installed=False,actual_bridge_own_moment_functions_closed=False,
            actual_R100_R110_switch_function_installed=False,whole_axis_functions_installed=False)


def run():
    began=time.monotonic();owner=OriginalBridgeMacroFunctions()
    packets={label:[owner.evaluate(label,q) for q in ((0,1),(1,2),(1,1))] for label in ('0','.5')}
    result=dict(**{GATE:True},source_family=owner.family,packets=packets,
        **dict.fromkeys(previous.OPEN,False),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
        scope='Original complete nonlinear frozen macro F/V integral functions at two genuine source frames; exact first exponential correction and nonzero remainder. Micro inlet errors remain, own moments and first switch not closed.')
    (HERE/NAME).write_text(json.dumps(base.encoded(result),indent=2)+'\n',encoding='utf8')
    print('Original nonlinear macro F/V source functions evaluated',flush=True);return result


if __name__=='__main__':run()
