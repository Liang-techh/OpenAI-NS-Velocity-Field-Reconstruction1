"""Original reference mixed source jets and correlated C1 phase integration.

Four actual root rows feed fixed-angle and implicit mixed chain rules. A
source-linked weighted-curvature inequality bounds the complete product,
preserving the microscopic factors rather than multiplying independent
Poisson peak covers. Fixed nonzero Z only; the full-Z route remains open.
"""
from dataclasses import dataclass
import gzip
import hashlib
import json
from pathlib import Path
from types import MappingProxyType
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_reference_phase_averaged_integrals as phase

whole=phase.whole;points=phase.points;base=phase.base;ep=phase.ep
HERE,PREFIX,sha=phase.HERE,phase.PREFIX,phase.sha
NAME=PREFIX+'current_original_reference_mixed_C1_integrals.json.gz'
RECEIPT=PREFIX+'current_original_reference_mixed_C1_integrals_check.json'
GATE='original_reference_actual_yZ_source_and_correlated_phase_C1_own_rate_integrals_enclosed'
C0,Y,Z,YZ=(0,0),(1,0),(0,1),(1,1);ORDERS=(C0,Y,Z,YZ)


class MixedJet:
    """Directed enclosures of four ordinary function derivatives."""
    def __init__(self,atlas,rows):
        if set(rows)!=set(ORDERS):raise ValueError('C0,y,Z,yZ source rows required')
        for value in rows.values():
            if type(value) is not whole.prior.ScaledEnclosure or value.scale.bases is not atlas.bases or value.ledger is not atlas.ledger:
                raise ValueError('One original atlas context/basis/ledger required')
        self.atlas=atlas;self.rows=MappingProxyType(dict(rows))
    def __getitem__(self,key):return self.rows[key]
    @classmethod
    def constant(cls,atlas,value):
        value=value if isinstance(value,whole.prior.ScaledEnclosure) else atlas.scalar(value)
        return cls(atlas,{key:value if key==C0 else atlas.scalar(0) for key in ORDERS})
    def coerce(self,other):
        if not isinstance(other,MixedJet):other=self.constant(self.atlas,other)
        if other.atlas is not self.atlas:raise ValueError('Same original mixed-jet owner required')
        return other
    def __neg__(self):return MixedJet(self.atlas,{key:-v for key,v in self.rows.items()})
    def __add__(self,other):
        other=self.coerce(other);a=self.atlas
        return MixedJet(a,{key:a.add(self[key],other[key]) for key in ORDERS})
    __radd__=__add__
    def __sub__(self,other):return self+-self.coerce(other)
    def __mul__(self,other):
        other=self.coerce(other);a=self.atlas
        return MixedJet(a,{(j,k):a.sum(self[(i,l)]*other[(j-i,k-l)]
            for i in range(j+1) for l in range(k+1)) for j,k in ORDERS})
    __rmul__=__mul__
    def reciprocal(self):
        a=self.atlas;c=a.ctx;value=self[C0];lo,hi=ep(value.coefficient)
        if not lo*hi>0:raise ValueError('Strict original signed nonzero function required')
        sign=1 if lo>0 else -1;positive=value*sign
        lower=positive.scale.evaluate()+c.ln(c.mpf(ep(positive.coefficient)[0]))
        inv=a.scalar(1).positive_divide(positive,lower)*sign
        return MixedJet(a,{C0:inv,Y:-self[Y]*inv*inv,Z:-self[Z]*inv*inv,
            YZ:a.add(self[Y]*self[Z]*inv*inv*inv*2,-self[YZ]*inv*inv)})
    def sincos(self):
        a=self.atlas;c=a.ctx;x=self[C0].finite_interval();sn,cs=a.scalar(c.sin(x)),a.scalar(c.cos(x))
        sine=MixedJet(a,{C0:sn,Y:cs*self[Y],Z:cs*self[Z],YZ:a.add(cs*self[YZ],-sn*self[Y]*self[Z])})
        cosine=MixedJet(a,{C0:cs,Y:-sn*self[Y],Z:-sn*self[Z],YZ:a.add(-sn*self[YZ],-cs*self[Y]*self[Z])})
        return sine,cosine
    def record(self):return {str(key):value.record() for key,value in self.rows.items()}


@dataclass(frozen=True)
class OriginalMixedFrame:
    owner:object
    family:dict
    left:object
    right:object
    Z:object
    roots:object
    query:object
    record:dict


def weighted_curvature_certificate(atlas,kernel):
    """Enclose complete source-linked mixed curvature products.

    n=r+cos(chi), h=1/sqrt(1+u^2), s=h^2, Dchi=s+2*r*n.
    T2_i=K_i*J, K_i=u_i*h,
    J=2*q^2/r^2*(sin(chi)*(Dchi-3*s)+s*(chi-psi)/r).
    The bound uses these identities BEFORE absolute values.
    """
    c=atlas.ctx;q=kernel.q.finite_interval();rmin=min(abs(v) for v in ep(kernel.r))
    if ep(q)[0]<c.mpf('.5') or not rmin>0:raise ValueError('Original q>=1/2 and strict signed r required')
    qu=c.mpf(ep(q)[1]);rl=c.mpf(rmin);D0=2*c.pi/rl
    C25=2*(20*c.power(5,c.mpf(5)/4)+c.power(5,c.mpf(3)/4)*D0*D0)
    C2=2*(20*c.power(5,c.mpf(3)/2)+c.power(5,c.mpf(3)/4)*D0*D0)
    inverse=16*qu**5/rl**4*C25;primitive=8*qu**5/rl**4*C2
    sym=lambda value:atlas.scalar(c.mpf((-ep(value)[1],ep(value)[1])))
    return dict(inverse=sym(inverse),primitive=sym(primitive),record=dict(
        original_K_i='u_i/h_inverse; equivalently u_i*hinv',
        original_T2_i='K_i*2*q^2/r^2*(sinchi*(Dchi-3*s)+s*(chi-psi)/r)',
        original_t_psi='-2*q*sinchi*Dchi/hinv^3',
        original_direction_first='t_i=K_i*(2*coschi-r)*t-2*q*K_i*hinv',
        source_identities=['s=hinv^2=1-r^2','Dchi=s+2*r*n','n=r+coschi','t=2*q*n/hinv'],
        inequalities=['sinchi^2<=s+2*abs(n)','Dchi<=s+2*abs(n)',
            'abs(Dchi-3*s)<=4*s+2*abs(n)','abs(chi-psi)<=2*pi',
            'k=abs(n)/hinv; 1+4*q^2*k^2>=1+k^2',
            '1+2*k<=sqrt(5)*sqrt(1+k^2)','4+2*k<=sqrt(20)*sqrt(1+k^2)',
            'hinv^(3/2)*(1+k^2)^(1/4)<=5^(1/4)',
            'abs(2*t)/(1+t^2)^3<=2/(1+t^2)^(5/2)',
            'abs(1-t^2)/(1+t^2)^3<=1/(1+t^2)^2'],
        bound_definition=dict(inverse='|2*t*t_psi*J^2/(1+t^2)^3|<=16*q_upper^5/r_lower^4*C25',
            primitive='|(1-t^2)*t_psi*J^2/(1+t^2)^3|<=8*q_upper^5/r_lower^4*C2'),
        original_q_directed_range=q,strict_r_magnitude_lower=rmin,C25=C25,C2=C2,
        inverse_weighted_curvature_absolute_bound=inverse,primitive_weighted_curvature_absolute_bound=primitive,
        microscopic_source_factors_collected_before_bounding=True,
        cap_not_selected_as_field_point_value=True))


def fixed_and_implicit_mixed(atlas,kernel,roots,coordinate):
    """Actual original signed E-chart mixed function enclosures."""
    a=atlas;c=a.ctx;const=lambda value:MixedJet.constant(a,value)
    if kernel.geometry!='signed_Mobius':raise ValueError('Original strict signed source required')
    for name in ('a','b','t0'):
        if any(not roots[name][key].zero for key in (Y,Z,YZ)):raise ValueError('Original constant reference parameter identities required')
    q=kernel.q;hinv=kernel.hinv;s0=kernel.s_source;r0=a.scalar(kernel.r)
    ui={key:(roots['p2'][key]*q).positive_divide(kernel.dstar,kernel.dstar.scale.evaluate()) for key in (Y,Z,YZ)}
    K={key:ui[key]*hinv for key in (Y,Z,YZ)};KY,KZ,KYZ=(K[key] for key in (Y,Z,YZ))
    rr=MixedJet(a,{C0:r0,Y:KY*s0,Z:KZ*s0,YZ:s0*a.add(KYZ,-r0*KY*KZ*3)})
    ss=MixedJet(a,{C0:s0,Y:-r0*KY*s0*2,Z:-r0*KZ*s0*2,
        YZ:-s0*a.add(r0*KYZ,a.add(s0,-base.current.square(r0)*3)*KY*KZ)*2})
    hh=MixedJet(a,{C0:hinv,Y:-r0*KY*hinv,Z:-r0*KZ*hinv,
        YZ:hinv*a.add(a.add(base.current.square(r0)*3,-a.scalar(1))*KY*KZ,-r0*KYZ)})
    psif,chif=kernel.angles(c.mpf(coordinate),'E');psi=const(2*c.pi*psif);chi0=a.scalar(2*c.pi*chif)
    sn=a.scalar(c.sin(2*c.pi*chif));cs=a.scalar(c.cos(2*c.pi*chif))
    chi=MixedJet(a,{C0:chi0,Y:KY*sn*2,Z:KZ*sn*2,
        YZ:sn*a.sum((KYZ,-r0*KY*KZ,KY*KZ*cs*2))*2})
    sine,cosine=chi.sincos();inv_r=rr.reciprocal();difference=chi-psi
    T1=const(q)*hh*inv_r*difference
    H=(const(2)-ss*3)*chi+ss*psi+rr*sine*2
    T2=const(base.current.square(q))*inv_r*inv_r*H
    beta={key:K[key]*a.add(cosine[C0]*2,-r0) for key in (Y,Z)}
    gamma={key:-q*K[key]*hinv*2 for key in (Y,Z)}
    R0=a.scalar(c.mpf((0,1)));Rt=a.scalar(c.mpf(['-.5','.5']))
    Rtt=a.scalar(c.mpf((0,c.mpf('.5'))));R2t=a.scalar(c.mpf((-1,1)))
    G2=a.scalar(c.mpf((-1,c.mpf(1)/8)));G2t=Rt
    curve=weighted_curvature_certificate(a,kernel)
    psiYZ=a.sum((-R0*T2[YZ],a.add(beta[Y]*Rtt,gamma[Y]*R2t)*T2[Z],
        a.add(beta[Z]*Rtt,gamma[Z]*R2t)*T2[Y],-curve['inverse']*KY*KZ))
    totalY=a.add(T1[Y],-Rt*T2[Y]);totalZ=a.add(T1[Z],-Rt*T2[Z])
    totalYZ=a.sum((T1[YZ],-Rt*T2[YZ],a.add(beta[Y]*G2t,gamma[Y]*G2)*T2[Z],
        a.add(beta[Z]*G2t,gamma[Z]*G2)*T2[Y],curve['primitive']*KY*KZ))
    total=MixedJet(a,{C0:T1[C0],Y:totalY,Z:totalZ,YZ:totalYZ})
    Bmixed=-(roots['a']*roots['E']*total)*(c.mpf(1)/(4*c.pi))
    AYZ=-roots['a'][C0]*psiYZ*(c.mpf(1)/(4*c.pi))
    primitive=kernel.primitives(coordinate,'E');first={}
    for order in (Y,Z):
        directional={name:{(0,0):row[C0],(0,1):row[order]} for name,row in roots.items()}
        values,proof=points.slow.slow_values(kernel,directional,coordinate,'E')
        B1,Bproof=whole.bounded_signed_implicit_B_Z(a,kernel,directional,coordinate,'E')
        first[order]=dict(A=values['A_Z_slow'],B=B1,proof=proof,Bproof=Bproof)
    A=MixedJet(a,{C0:primitive['A'],Y:first[Y]['A'],Z:first[Z]['A'],YZ:AYZ})
    B=MixedJet(a,{C0:primitive['B_over_Pstar'],Y:first[Y]['B'],Z:first[Z]['B'],YZ:Bmixed[YZ]})
    return dict(A=A,B=B,record=dict(original_fixed_angle_T1=T1.record(),original_fixed_angle_T2=T2.record(),
        actual_source_K_rows={str(key):v.record() for key,v in K.items()},
        original_r_s_hinv_mixed_jets=dict(r=rr.record(),s=ss.record(),hinv=hh.record()),
        fixed_angle_chi_mixed=chi.record(),original_total_T1_mixed=total.record(),
        original_inverse_psi_yZ=psiYZ.record(),correlated_weighted_curvature=curve['record'],
        exact_rational_function_ranges=dict(inv_one_plus_t_squared=[0,1],t_over_one_plus_t_squared=['-1/2','1/2'],
            two_t_squared_over_denominator_squared=[0,'1/2'],two_t_over_denominator_squared=[-1,1],
            t_squared_minus_one_over_denominator_squared=[-1,'1/8'],t_times_t_squared_minus_one_over_denominator_squared=['-1/2','1/2']),
        original_A_and_B_mixed=dict(A=A.record(),B=B.record()),
        both_fixed_angle_first_cross_terms_and_curvature_retained=True,
        phase_held_fixed_and_original_phi_Z_exact_zero=True,
        source_enclosures_not_saved_cap_or_midpoint_values=True))


class OriginalReferenceMixedC1:
    def __init__(self,*,Z):
        self.parent=phase.OriginalReferencePhaseAveraging(Z=Z);self.atlas=a=self.parent.atlas
        self.ctx=c=a.ctx;self.family=self.parent.family;self.Z=self.parent.Z;self.frames={};self.cache={}
        receipt=json.loads((HERE/phase.RECEIPT).read_bytes())
        if not receipt.get('all_passed') or not receipt.get(phase.GATE) or receipt['source_family']!=self.family:
            raise ValueError('Accepted original radial/phase-averaging source required')
        self.hashes=dict(self.parent.hashes)
        for name,digest in {**receipt['input_hashes'],phase.RECEIPT:sha(phase.RECEIPT)}.items():
            if sha(name)!=digest:raise ValueError('Original mixed source dependency changed: '+name)
            if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Mixed source families disagree')
            self.hashes[name]=digest
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        template=self.parent.owner.templates;inputs=template['inputs'];radial=inputs[1:5]
        self.templates={};self.compiled={};self.pressure={}
        for name in ('E','V','b','p1','p2'):
            rows=[];compiled=[];pressure=[]
            for powers,expr in template['rows'][(name,1)]:
                derivative=powers[0]*expr+sum(s.diff(expr,x)*rate*x for x,rate in
                    zip(radial,(s.Rational(1,10),s.Rational(1,10),s.Rational(1,5),s.Rational(1,5)),strict=True))
                rows.append((powers,derivative));compiled.append(s.lambdify(inputs,derivative,modules=[{'mpf':c.mpf},'mpmath']))
                pressure.append(tuple(s.lambdify(inputs,s.diff(derivative,P),modules=[{'mpf':c.mpf},'mpmath'])
                    for P in template['pressure_symbols']))
            self.templates[name]=tuple(rows);self.compiled[name]=tuple(compiled);self.pressure[name]=tuple(pressure)

    def source_frame(self,left,right):
        with mp.workdps(self.ctx.dps+40):
            return self._source_frame(left,right)

    def _source_frame(self,left,right):
        query=self.parent.radial_query(left,right);a=self.atlas;c=self.ctx;y=query['coordinate']
        z=a.rational(self.Z);q=1+z*z;f=c.exp(y/10);alpha=a.copy_interval(self.parent.owner.inputs.alpha_enclosure)
        values=(z,f,c.mpf(5)/8*f,c.mpf(5)/12*f*f,c.mpf(5)/2*f*f,
            -alpha/q**2,4*alpha*z/q**3,alpha*(4-20*z*z)/q**4)
        roots={};records=[]
        for name,rows in self.templates.items():
            mixed=a.scalar(0);terms=[]
            for i,(powers,expr) in enumerate(rows):
                coefficient=c.mpf(self.compiled[name][i](*values));mixed=a.add(mixed,a.term(powers,coefficient,coordinate=y))
                late=c.mpf(0)
                for sensitivity,bound in zip(self.pressure[name][i],(5,10,44),strict=True):
                    late+=abs(c.mpf(sensitivity(*values)))*c.exp(c.mpf(3)/5)*bound/(2*q*q)
                if ep(late)[1]:
                    r,p,d,ell=powers;upper=ep(late)[1]
                    mixed=a.add(mixed,a.term((r,p-1,d,ell),c.mpf((-upper,upper)),coordinate=y))
                terms.append(dict(original_factor_powers=powers,ordinary_yZ_coefficient=coefficient,
                    pressure_late_Pstar_factor=-1 if ep(late)[1] else None,positive_late_budget=late))
            roots[name]=MixedJet(a,{C0:query['roots'][name][(0,0)],Y:query['radial_roots'][name][(0,1)],
                Z:query['roots'][name][(0,1)],YZ:mixed})
            records.append(dict(original_input=name,ordinary_mixed_source_terms=terms))
        for name in ('a','t0'):roots[name]=MixedJet.constant(a,query['roots'][name][(0,0)])
        frame=OriginalMixedFrame(self,self.family,query['left'],query['right'],self.Z,MappingProxyType(roots),query,
            dict(source_family=self.family,original_Z_exact=str(self.Z),exact_reference_cell=[str(query['left']),str(query['right'])],
                original_mixed_root_terms=records,actual_root_jets={name:row.record() for name,row in roots.items()},
                ordinary_Z_rows_differentiated_radially=True,original_L_Z_not_applied_twice=True,
                shared_P0_y_and_yZ_exact_zero=True,source_pressure_late_factors_retained=True))
        self.frames[id(frame)]=frame;return frame

    def primitive(self,frame,coordinate):
        if type(frame) is not OriginalMixedFrame or self.frames.get(id(frame)) is not frame or frame.owner is not self:
            raise ValueError('Issued original mixed source frame required')
        with mp.workdps(self.ctx.dps+40):
            return fixed_and_implicit_mixed(self.atlas,frame.query['kernel'],frame.roots,coordinate)

    def cell(self,left,right,*,N):
        N=points.candidate_N(N);left=points.reference.reference_coordinate(left);right=points.reference.reference_coordinate(right)
        key=(left,right,N)
        if key in self.cache:return self.cache[key]
        a=self.atlas;c=self.ctx
        with mp.workdps(c.dps+40):
            frame=self.source_frame(left,right);primitive=self.primitive(frame,c.mpf((0,1)));A,B=primitive['A'],primitive['B']
            E,V=frame.roots['E'],frame.roots['V'];EA=E*A;EEA=E*EA
            leading=dict(m=B,h=EA,k=V*EA+E*B,e=V*B*2-EEA,p=EEA)
            g=whole.DirectedCoefficientAlgebra(a);pair=lambda row:points.exact.source.C1Function(row[C0],row[Z])
            coefficients,F=points.exact.coefficient_pairs(g,pair(E),pair(V),pair(A),pair(B),g.constant(N))
            x=whole.conditioned.bounded(A[C0]*(c.mpf(1)/N));R2=a.scalar(phase.second_exponential_remainder(c,x))
            M=a.scalar(points.directed_exprel(c,x));Q=E[C0]*base.current.square(A[C0])*R2
            QZ=a.add(E[Z]*base.current.square(A[C0])*R2,E[C0]*A[C0]*A[Z]*M)
            FF=base.current.square(F.value);BB=base.current.square(B[C0])
            rest={C0:dict(m=a.scalar(0),h=Q,k=a.add(V[C0]*Q,F.value*B[C0]),
                e=a.sum((-E[C0]*Q,BB,-FF*(c.mpf(1)/2))),p=a.add(E[C0]*Q,FF*(c.mpf(1)/2))),
                Z:dict(m=a.scalar(0),h=QZ,k=a.sum((V[Z]*Q,V[C0]*QZ,F.Z*B[C0],F.value*B[Z])),
                    e=a.sum((-E[Z]*Q,-E[C0]*QZ,B[C0]*B[Z]*2,-F.value*F.Z)),
                    p=a.sum((E[Z]*Q,E[C0]*QZ,F.value*F.Z)))}
            direct={jet:{name:a.add(getattr(coefficients[-1][name],attr)*(c.mpf(1)/N),
                getattr(coefficients[-2][name],attr)*(c.mpf(1)/N**2)) for name in points.exact.RATES}
                for jet,attr in ((C0,'value'),(Z,'Z'))}
            caps={jet:{name:dict(leading=phase.magnitude_bound(a,row[jet]),
                slow_y=phase.magnitude_bound(a,row[Y if jet==C0 else YZ]),
                nonlinear_remainder=phase.magnitude_bound(a,rest[jet][name]),
                direct_density=phase.magnitude_bound(a,direct[jet][name])) for name,row in leading.items()}
                for jet in (C0,Z)}
            endpoints={which:self.parent.owner.dispatcher.reference.owner.radius.evaluate(y=y,N=N)
                for which,y in (('left',left),('right',right))}
        with mp.workdps(c.dps+40):
            record=dict(source=frame.record,actual_mixed_primitive=primitive['record'],candidate_N=N,
                leading_density_mixed_jets={name:row.record() for name,row in leading.items()},
                actual_N_nonlinear_C0_Z_remainder={str(jet):{name:row.record() for name,row in rows.items()} for jet,rows in rest.items()},
                actual_original_endpoint_phases=endpoints,
                C0_Z_averaging_caps={str(jet):{name:{kind:row.record() for kind,row in pair.items()} for name,pair in rows.items()} for jet,rows in caps.items()},
                complete_Z_product_rules_and_exact_N_remainder_retained=True,
                phase_Z_exact_zero=True,native_log_radius_Jacobian=1)
            result=dict(record=record,caps=caps);self.cache[key]=result;return result

    def integrate(self,*,count,N):
        N=points.candidate_N(N)
        if type(count) is not int or not 1<=count<=4096:raise ValueError('Exact reference partition required')
        a=self.atlas;c=self.ctx;total={jet:{name:a.scalar(0) for name in points.exact.RATES} for jet in (C0,Z)}
        direct={jet:{name:a.scalar(0) for name in points.exact.RATES} for jet in (C0,Z)};cells=[]
        with mp.workdps(c.dps+40):
            for i in range(count):
                left=-5+s.Rational(5*i,count);right=-5+s.Rational(5*(i+1),count);got=self.cell(left,right,N=N)
                contributions={}
                for name,rate in points.exact.RATES.items():
                    rate=c.mpf(rate.numerator)/rate.denominator;wl,wr=c.exp(rate*a.rational(left)),c.exp(rate*a.rational(right))
                    mass=a.rational(right-left) if ep(rate)==(0,0) else (wr-wl)/rate
                    if ep(mass)[0]<=0:raise ArithmeticError('Original positive own-rate mass required')
                    rows={}
                    for jet in (C0,Z):
                        caps=got['caps'][jet][name];G=caps['leading']*(c.mpf(1)/2);Gy=caps['slow_y']*(c.mpf(1)/2)
                        endpoint=G*((wl if i==0 else c.mpf(0))+(wr if i==count-1 else c.mpf(0)))
                        slow=a.add(Gy,G*rate)*mass;rest=caps['nonlinear_remainder']*mass
                        bound=a.sum((endpoint,slow,rest))*(c.mpf(1)/N**2)
                        total[jet][name]=a.add(total[jet][name],phase.magnitude_bound(a,bound))
                        direct[jet][name]=a.add(direct[jet][name],phase.magnitude_bound(a,caps['direct_density']*mass))
                        rows[str(jet)]=dict(global_endpoint_term=endpoint.record(),slow_and_kernel_term=slow.record(),
                            exact_N_remainder_term=rest.record(),source_IBP_enclosure=phase.symmetric_bound(a,bound).record())
                    contributions[name]=dict(positive_own_rate_mass=mass,C0_Z_rows=rows,N_power=-2)
                cells.append(dict(source=got['record'],contributions=contributions))
        # The legacy magnitude helpers use real endpoint abs/negation.
        # Keep export inside the higher-precision context too: every directed
        # endpoint is an exact dyadic with fewer bits than this work precision.
        with mp.workdps(c.dps+40):
            selected={jet:{name:phase.tighter_source_bound(a,total[jet][name],direct[jet][name]) for name in total[jet]} for jet in total}
            encode=lambda rows:{str(jet):{name:phase.symmetric_bound(a,value).record() for name,value in pair.items()} for jet,pair in rows.items()}
            return dict(source_family=self.family,original_Z_exact=str(self.Z),candidate_N=N,exact_radial_cells=count,
                source_window=['-5','0'],C0_Z_phase_averaged_contribution_enclosures=encode(total),
                C0_Z_direct_contribution_enclosures=encode(direct),
                C0_Z_effective_contribution_enclosures=encode({jet:{name:value[0] for name,value in pair.items()} for jet,pair in selected.items()}),
                effective_bound_selection={str(jet):{name:value[1] for name,value in pair.items()} for jet,pair in selected.items()},
                whole_cell_mixed_source_and_IBP_records=cells,actual_global_C0_Z_endpoint_terms_retained=True,
                same_original_source_internal_C0_Z_traces_cancel=True,no_cross_chart_seam_assumed=True,
                incoming_histories_and_P0_not_reset=True,pressure_zero_rate_memory_retained=True,
                original_interval_export_work_precision_retained=True,
                actual_mixed_source_and_C1_averaging_installed_on_this_fixed_Z_reference_window=True,
                whole_Z_terminal_or_all_route_closure=False)


def run():
    begin=time.monotonic();owner=OriginalReferenceMixedC1(Z='.37');levels=[]
    for count,N in ((4,160),(16,160),(16,320),(16,16384)):
        levels.append(owner.integrate(count=count,N=N))
        print('Original mixed C1 reference contribution:',count,'cells, N',N,flush=True)
    report=dict(**{GATE:True},source_family=owner.family,actual_original_mixed_C1_reference_levels=levels,
        actual_mixed_yZ_source_and_reference_C1_averaging_installed=True,
        actual_correlated_weighted_curvature_function_enclosures_installed=True,fixed_nonzero_Z_only=True,
        all_17_chart_or_24_cell_oracle_installed=False,actual_five_controls_installed=False,current_whole_N_selected=False,
        **dict.fromkeys(whole.point.source.inertial.profiles.loop.OPEN,False),input_hashes=owner.hashes,
        execution_seconds=time.monotonic()-begin,
        scope='Actual original reference C0,y,Z,yZ roots and phase-held mixed primitives, preserving both inverse cross terms and source-correlated curvature. C0/Z own-rate endpoint-retaining N-dependent integration at fixed nonzero Z; dual direct/averaged bounds. Not whole-Z/midplane, terminal/all-route controls, global frequency or recursive corrected field.')
    (HERE/NAME).write_bytes(gzip.compress(json.dumps(base.encoded(report),indent=2).encode()+b'\n',compresslevel=9,mtime=0))
    return report


if __name__=='__main__':run()
