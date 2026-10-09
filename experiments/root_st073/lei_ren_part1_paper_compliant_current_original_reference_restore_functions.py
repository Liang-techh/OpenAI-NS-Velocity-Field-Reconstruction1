"""Actual Rsh histories through reference and original axial restoration to Rm.

Same-family accepted terminal functions replace broad legacy inlets. Signed
reference tails and full original scalar restoration kernels are retained.
Normalized and physical mixed y/Z rows are provided, before finite-N integrals
or the implicit active patch. Native frames0,.5 remain conditional upstream.
"""
import ast
from fractions import Fraction
import json
import math
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_long_reshape_endpoint as previous
import lei_ren_part1_paper_compliant_reference_restore_profiles as original
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import sigma_jets
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_compliant_inner_bridge_profiles import logarithm

fields,endpoint=previous.fields,previous.endpoint
prior,base,ep=previous.prior,previous.base,previous.ep
HERE,PREFIX,sha=previous.HERE,previous.PREFIX,previous.sha
NAME=PREFIX+'current_original_reference_restore_functions.json'
RECEIPT=PREFIX+'current_original_reference_restore_functions_check.json'
GATE='original_actual_two_frame_Rsh_reference_restoration_Rm_functions_installed'
ADMISSION=PREFIX+'five_defect_admission.json'
ADMISSION_CHECK=PREFIX+'actual_moment_patch_check.json'
RATES=dict(mean_error=1,angular_error='1.6',mixed_error='1.6',axial_square=1,swirl_error='1.2',pressure_error='.2')


def source_bindings():
    from lei_ren_part1_paper_compliant_actual_switch_mixed_C4 import assignment_source_bindings
    bindings=assignment_source_bindings('reference_restore_profiles','inputs',{
        'centered':"dict(mean_error=moments['mean']-z*4,angular_error=moments['theta']-c.mpf(5)/8,mixed_error=moments['theta_z']-(z*moments['theta'])*4,axial_square=moments['axial']-(z*moments['mean'])*8+square(z)*16,swirl_error=moments['swirl']-c.mpf(5)/6,pressure_error=moments['pressure']-5)"})
    bindings.update(assignment_source_bindings('reference_restore_profiles','restoration',{
        'V':"inp['z']*4+inp['E']*alpha",'logu':"-logarithm(1+square(inp['z']))+(-8+t)/10"}))
    bindings.update(assignment_source_bindings('reference_restore_profiles','terminal',{
        'rates':"dict(mean_error=1,angular_error='1.6',mixed_error='1.6',axial_square=1,swirl_error='1.2',pressure_error='.2')",
        'logu':"-logarithm(1+square(inp['z']))+offset/10"}))
    bindings.update(assignment_source_bindings('reference_restore_mixed_C4','source_mixed',{
        'Q':"[(2*z*V[k]-(z*m[k])*(1-delta)-d*derivative(m[k]))/L for k in range(5)]",
        'physical':"dict(Utheta_over_current_Utheta=[amp*c.mpf('.1')**k for k in range(5)],Uz=V,Ur_over_current_sqrt_R_over_2=[binomial_rate(Q,c.mpf('.5'),k) for k in range(5)],P_over_Pstar2=[p0+ratio*p[0]/2]+[ratio*binomial_rate(p,c.mpf('.2'),k)/2 for k in range(1,5)])",
        'primitives':"dict(Mtheta_over_current_sqrt2_R_1p5_Utheta=[amp*binomial_rate(H,c.mpf('1.6'),k) for k in range(5)],Mtheta_z_over_current_sqrt2_R_1p5_Utheta=[amp*binomial_rate(K,c.mpf('1.6'),k) for k in range(5)],Mz_over_current_R=[binomial_rate(m,c.mpf(1),k) for k in range(5)],Mztheta_over_current_R_Pstar2=[binomial_rate(A,c.mpf(1),k)*invP2-ratio*binomial_rate(b,c.mpf('1.2'),k)/2 for k in range(5)],Mp_over_Pstar2=[ratio*binomial_rate(p,c.mpf('.2'),k)/2 for k in range(5)])"}))
    expected={
        'reference_centered':"dict(mean_error=E+scale(initial['mean_error']-E,1),angular_error=scale(initial['angular_error'],'1.6'),mixed_error=scale(initial['mixed_error'],'1.6')+E*((1-d)/c.mpf('1.6')),axial_square=square(E)+scale(initial['axial_square']-square(E),1),swirl_error=scale(initial['swirl_error'],'1.2'),pressure_error=scale(initial['pressure_error'],'.2'))",
        'restore_centered':"dict(mean_error=initial['mean_error']*c.exp(-t)+E*kernels['mean'],angular_error=initial['angular_error']*c.exp(-c.mpf('1.6')*t),mixed_error=initial['mixed_error']*c.exp(-c.mpf('1.6')*t)+E*kernels['mixed'],axial_square=initial['axial_square']*c.exp(-t)+square(E)*kernels['square'],swirl_error=initial['swirl_error']*c.exp(-c.mpf('1.2')*t),pressure_error=initial['pressure_error']*c.exp(-c.mpf('.2')*t))"}
    tree=ast.parse(Path(original.__file__).read_text(encoding='utf8'))
    for name,text in expected.items():
        fn=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name==name)
        want=ast.dump(ast.parse(text,mode='eval').body)
        if sum(ast.dump(n.value)==want for n in ast.walk(fn) if isinstance(n,ast.Return))!=1:
            raise ValueError('Original centered source transport changed: '+name)
    return dict(passed=True,original_assignment_bindings=bindings,original_centered_return_bindings=list(expected),
        scalar_endpoint_kernel='J(k,j;1)=exp(-k)*admitted_integral_0^1 exp(k*s)*(1-sigma(s))^j ds',
        actual_E_source='same actual V_Rsh minus4Z, not a fitted centered cap',
        exact_radius_identity='Rref=110exp(10(logCstar+logPstar)); Rz/Rrestore/Rm offsets=-8/-7/-6',
        exact_Rz_amplitude_identity='T/10-logCstar-logPstar+(10(logCstar+logPstar)-T-8)/10=-.8',
        same_P0_and_fixed_T_function_identity_not_interval_overlap=True,
        old_patch_receipt_used_only_for_scalar_cutoff_integrals_not_patch_closure=True)


def fraction(value):
    if type(value) is not tuple or len(value)!=2 or any(type(v) is not int for v in value) or value[1]<=0:
        raise ValueError('Exact rational coordinate pair required')
    return Fraction(*value)


class RestorationKernelProvider:
    def __init__(self,c,endpoint_integrals):
        if set(endpoint_integrals)!=set(('linear_1','linear_8_over_5','square_1')):
            raise ValueError('All three original admitted cutoff integrals required')
        self.c=c;self.integrals={name:c.mpf(row) for name,row in endpoint_integrals.items()};self.cache={}
        if any(not 0<ep(row)[0]<=ep(row)[1] or not all(mp.isfinite(v) for v in ep(row)) for row in self.integrals.values()):
            raise ValueError('Finite positive restoration source intervals required')

    def evaluate(self,phase):
        q=fraction(phase)
        if not 0<=q<=1:raise ValueError('Original restoration phase in[0,1] required')
        if q in self.cache:return self.cache[q]
        c=self.c;t=c.mpf(q.numerator)/q.denominator
        if q==0:row=dict.fromkeys(('mean','mixed','square'),c.mpf(0))
        elif q==1:row={name:self.integrals[part]*c.exp(-c.mpf(rate))
            for name,part,rate in (('mean','linear_1',1),('mixed','linear_8_over_5','1.6'),('square','square_1',1))}
        else:row=original.restoration_kernels(c,t,cells=128)
        self.cache[q]=row;return row


class ActualReferenceRestoreFunctions:
    def __init__(self,flow,Z,delta,actual_shapes,V,P0,T,logC,kernel_provider):
        f=self.flow=flow;self.c=c=f.c;self.Z=c.mpf(Z);self.delta=c.mpf(delta)
        if ep(self.Z)[0]<-1 or ep(self.Z)[1]>1 or ep(self.delta)[0]<0 or ep(self.delta)[1]>=1:
            raise ValueError('Original |Z|<=1 and0<=delta<1 required')
        if set(actual_shapes)!=set(previous.NAMES):raise ValueError('All six actual Rsh shapes required')
        for row in (*actual_shapes.values(),V,P0):
            if len(row)!=6 or any(v.scale.bases is not f.logs or v.ledger is not f.ledger for v in row):
                raise ValueError('Same actual Rsh basis/ledger and ordinary Z0..5 required')
        if not isinstance(kernel_provider,RestorationKernelProvider) or kernel_provider.c is not c:
            raise ValueError('Same-context original scalar cutoff provider required')
        self.provider=kernel_provider;self.initial=actual_shapes;self.V=V;self.P0=P0
        self.T=c.mpf(T);self.logC=c.mpf(logC);self.logP=f.logs[1]/2
        self.logref=10*(self.logC+self.logP);self.gap=self.logref-self.T-8
        if ep(self.T)[0]<=0 or ep(self.gap)[0]<=0 or any(not mp.isfinite(v) for row in (self.T,self.logC,self.gap) for v in ep(row)):
            raise ValueError('Original Rsh<Rz and finite positive fixed T required')
        self.z=IntervalTaylor.variable(c,self.Z,5);self.zrows=f.jet(self.z);self.Vref=f.scale(self.zrows,4)
        self.E=f.add(V,f.scale(self.Vref,-1));self.E2=f.multiply(self.E,self.E)
        zero=[f.scalar(0)]*6;unit=lambda value:[f.scalar(value)]+zero[1:]
        s=actual_shapes
        self.centered=dict(mean_error=f.add(s['mean'],f.scale(self.Vref,-1)),
            angular_error=f.add(s['theta'],unit(-c.mpf(5)/8)),
            mixed_error=f.add(s['theta_z'],f.scale(f.multiply(self.zrows,s['theta']),-4)),
            axial_square=f.add(s['axial'],f.scale(f.multiply(self.zrows,s['mean']),-8),f.scale(f.multiply(self.zrows,self.zrows),16)),
            swirl_error=f.add(s['swirl'],unit(-c.mpf(5)/6)),pressure_error=f.add(s['pressure'],unit(-5)))
        self.references={};self.restores={};self.posts={}

    def scalar_exp(self,log):return self.flow.factor((0,0,0,0,0),log)
    def cv(self,q):return self.c.mpf(q.numerator)/q.denominator

    def reference(self,phase):
        q=fraction(phase)
        if not 0<=q<=1:raise ValueError('Original reference phase in[0,1] required')
        if q in self.references:return self.references[q]
        f=self.flow;c=self.c;t=self.cv(q);g=self.gap*t;decays={name:self.scalar_exp(-c.mpf(rate)*g) for name,rate in RATES.items()}
        initial=self.centered
        if q==0:cent=initial
        else:
            cent=dict(mean_error=f.add(self.E,f.scale(f.add(initial['mean_error'],f.scale(self.E,-1)),decays['mean_error'])),
                angular_error=f.scale(initial['angular_error'],decays['angular_error']),
                mixed_error=f.add(f.scale(initial['mixed_error'],decays['mixed_error']),f.scale(self.E,(1-decays['mixed_error'])*(c.mpf(5)/8))),
                axial_square=f.add(self.E2,f.scale(f.add(initial['axial_square'],f.scale(self.E2,-1)),decays['axial_square'])),
                swirl_error=f.scale(initial['swirl_error'],decays['swirl_error']),pressure_error=f.scale(initial['pressure_error'],decays['pressure_error']))
        logq=logarithm(1+self.z*self.z)
        logu=-logq+((self.T/10-self.logC-self.logP)*(1-t)-c.mpf('.8')*t)
        logR=c.ln(110)+self.T*(1-t)+(self.logref-8)*t
        value=self.packet(cent,self.V,[c.mpf(1)]+[c.mpf(0)]*4,logu,logR,
            'Rsh_to_Rz',offset=-8 if q==1 else None,shapes_override=self.initial if q==0 else None)
        value.update(geometry=dict(phase=[q.numerator,q.denominator],log_gap_from_Rsh=g,
            original_reference_length=self.gap,log_radius=logR,log_radius_Jacobian=1,
            exact_source='Rsh*exp(phase*(10(logCstar+logPstar)-T-8))'),
            signed_nonzero_incoming_reference_decays=decays,actual_Rsh_centered_inlet=initial)
        self.references[q]=value;return value

    def restoration(self,phase):
        q=fraction(phase)
        if not 0<=q<=1:raise ValueError('Original restoration phase in[0,1] required')
        if q in self.restores:return self.restores[q]
        f=self.flow;c=self.c;t=self.cv(q);initial=self.reference((1,1))['actual_centered_histories'];K=self.provider.evaluate(phase)
        if q==0:cent=initial
        else:
            cent={name:f.scale(row,c.exp(-c.mpf(RATES[name])*t)) for name,row in initial.items()}
            cent['mean_error']=f.add(cent['mean_error'],f.scale(self.E,K['mean']))
            cent['mixed_error']=f.add(cent['mixed_error'],f.scale(self.E,K['mixed']))
            cent['axial_square']=f.add(cent['axial_square'],f.scale(self.E2,K['square']))
        cutoff=sigma_jets(c,t);alpha=[1-cutoff[0]]+[-cutoff[k]*math.factorial(k) for k in range(1,5)]
        V=self.V if q==0 else self.Vref if q==1 else f.add(self.Vref,f.scale(self.E,alpha[0]))
        offset=-8+t;logu=-logarithm(1+self.z*self.z)+offset/10;logR=c.ln(110)+self.logref+offset
        value=self.packet(cent,V,alpha,logu,logR,'original_axial_restoration',offset=offset)
        value.update(geometry=dict(phase=[q.numerator,q.denominator],exact_reference_offset=offset,
            log_radius=logR,log_radius_Jacobian=1,exact_source='Rz*exp(t); Rz=Rref*exp(-8)'),
            full_original_restoration_kernels=K,actual_Rz_centered_inlet=initial,
            endpoint_integrals_reused_from_original_admission=q==1)
        self.restores[q]=value;return value

    def postrestore(self,offset):
        q=fraction(offset)
        if not -7<=q<=-6:raise ValueError('Unpatched postrestore stops at Rm, offsets[-7,-6]')
        if q in self.posts:return self.posts[q]
        f=self.flow;c=self.c;x=self.cv(q);s=x+7;initial=self.restoration((1,1))['actual_centered_histories']
        cent=initial if q==-7 else {name:f.scale(row,c.exp(-c.mpf(RATES[name])*s)) for name,row in initial.items()}
        logu=-logarithm(1+self.z*self.z)+x/10;logR=c.ln(110)+self.logref+x
        value=self.packet(cent,self.Vref,[c.mpf(0)]*5,logu,logR,'restore_exit_to_Rm',offset=x)
        value.update(geometry=dict(exact_reference_offset=[q.numerator,q.denominator],log_radius=logR,
            log_radius_Jacobian=1,exact_source='Rref*exp(offset), offsets[-7,-6] before active patch'),
            actual_restoration_exit_centered_inlet=initial)
        self.posts[q]=value;return value

    def packet(self,centered,V,alpha,logu,logR,chart,*,offset=None,shapes_override=None):
        f=self.flow;c=self.c;zero=[f.scalar(0)]*6;unit=lambda value:[f.scalar(value)]+zero[1:]
        H=f.add(centered['angular_error'],unit(c.mpf(5)/8))
        shapes=dict(theta=H,theta_z=f.add(f.scale(f.multiply(self.zrows,H),4),centered['mixed_error']),
            mean=f.add(self.Vref,centered['mean_error']),
            axial=f.add(f.scale(f.multiply(self.zrows,self.zrows),16),f.scale(f.multiply(self.zrows,centered['mean_error']),8),centered['axial_square']),
            swirl=f.add(centered['swirl_error'],unit(c.mpf(5)/6)),pressure=f.add(centered['pressure_error'],unit(5)))
        if shapes_override is not None:shapes=shapes_override
        mismatch=[f.scale(self.E,v) for v in alpha];Vy=[V]+mismatch[1:]
        yr={name:[row] for name,row in centered.items()}
        for k in range(4):
            yr['mean_error'].append(f.add(mismatch[k],f.scale(yr['mean_error'][k],-1)))
            yr['angular_error'].append(f.scale(yr['angular_error'][k],-c.mpf('1.6')))
            yr['mixed_error'].append(f.add(mismatch[k],f.scale(yr['mixed_error'][k],-c.mpf('1.6'))))
            square=f.add(*[f.scale(f.multiply(mismatch[j],mismatch[k-j]),math.comb(k,j)) for j in range(k+1)])
            yr['axial_square'].append(f.add(square,f.scale(yr['axial_square'][k],-1)))
            yr['swirl_error'].append(f.scale(yr['swirl_error'][k],-c.mpf('1.2')))
            yr['pressure_error'].append(f.scale(yr['pressure_error'][k],-c.mpf('.2')))
        sy={name:[shapes[name]] for name in previous.NAMES}
        for k in range(1,5):
            sy['theta'].append(yr['angular_error'][k])
            sy['theta_z'].append(f.add(f.scale(f.multiply(self.zrows,yr['angular_error'][k]),4),yr['mixed_error'][k]))
            sy['mean'].append(yr['mean_error'][k])
            sy['axial'].append(f.add(f.scale(f.multiply(self.zrows,yr['mean_error'][k]),8),yr['axial_square'][k]))
            sy['swirl'].append(yr['swirl_error'][k]);sy['pressure'].append(yr['pressure_error'][k])
        Q=[previous.radial_Q(f,self.Z,self.delta,Vy[k],sy['mean'][k]) for k in range(5)]
        Qfull=[row+[f.scalar(0)] for row in Q]
        ratio=IntervalTaylor(c,[c.mpf(0)]+list(logu.coefficients[1:])).exp();ratio2=IntervalTaylor(c,[c.mpf(0)]+[2*v for v in logu.coefficients[1:]]).exp()
        if offset is None:
            R=self.scalar_exp(logR);rootR=self.scalar_exp(logR/2);U=self.scalar_exp(logu[0]+self.logP)
        else:
            # At exact reference offsets, C/P source factors collect before
            # multiplication. logR/10 and logU cannot erase microscopic data.
            R=f.factor((0,5,0,0,0),10*self.logC+c.ln(110)+offset)
            rootR=f.factor((0,2.5,0,0,0),5*self.logC+c.ln(110)/2+offset/2)
            U=f.factor((0,.5,0,0,0),logu[0])
        U2=U*U;P2=f.factor((0,1,0,0,0));axis=f.scale(self.P0,P2)
        def binrate(rows,rate,k):return f.add(*[f.scale(rows[j],math.comb(k,j)*c.mpf(rate)**(k-j)) for j in range(k+1)])
        dress=lambda row,power:f.multiply(row,f.jet(ratio if power==1 else ratio2))
        velocity=dict(Ur=[f.scale(binrate(Qfull,'.5',k),rootR*(1/c.sqrt(2)))[:5] for k in range(5)],
            Utheta=[f.scale(f.jet(ratio),U*c.mpf('.1')**k) for k in range(5)],Uz=Vy)
        increments=[f.scale(dress(binrate(sy['pressure'],'.2',k),2),U2*c.mpf('.5')) for k in range(5)]
        pressure=[f.add(axis,increments[0])]+increments[1:]
        primitive=dict(Mtheta=[f.scale(dress(binrate(sy['theta'],'1.6',k),1),R*rootR*U*c.sqrt(2)) for k in range(5)],
            Mtheta_z=[f.scale(dress(binrate(sy['theta_z'],'1.6',k),1),R*rootR*U*c.sqrt(2)) for k in range(5)],
            Mz=[f.scale(binrate(sy['mean'],1,k),R) for k in range(5)],
            Mztheta=[f.add(f.scale(binrate(sy['axial'],1,k),R),f.scale(dress(binrate(sy['swirl'],'1.2',k),2),-R*U2*c.mpf('.5'))) for k in range(5)],Mp=increments)
        grid=lambda rows:{'y'+str(k)+'_Z'+str(n):row[n]*math.factorial(n) for k,row in enumerate(rows) for n in range(5-k)}
        return dict(chart=chart,actual_centered_histories=centered,actual_normalized_six_history_shapes=shapes,
            actual_centered_y_derivative_axial5=yr,actual_normalized_shape_y_derivative_axial5=sy,
            actual_Q_y_derivative_axial4=Q,actual_E_V110_minus4Z=self.E,actual_velocity_V=V,
            original_P0_normalized_axial5=self.P0,log_Utheta_over_Pstar_axial5=list(logu.coefficients),
            actual_alpha_ordinary_y_derivatives=alpha,
            physical_velocity_y_Z_mixed4={name:grid(rows) for name,rows in velocity.items()},
            physical_pressure_y_Z_mixed4=grid(pressure),physical_five_primitive_y_Z_mixed4={name:grid(rows) for name,rows in primitive.items()},
            physical_velocity_axial_coefficients={name:rows[0] for name,rows in velocity.items()},
            physical_pressure_axis_axial5=axis,physical_pressure_radial_increment_axial5=increments[0],physical_total_pressure_axial5=pressure[0],
            physical_cumulative_moment_axial5={name:rows[0] for name,rows in primitive.items()},
            formal_positive_geometry=dict(R=R,Utheta=U,sqrtR=rootR),
            actual_incoming_reference_restoration_memory_retained=True,separate_original_P0_retained=True,
            same_fixed_T_Z_exact_zero=True,radial_prefactors_differentiated_before_mixed_grid=True,
            actual_source_coefficients_not_selected_from_caps=True,positive_amplitudes_not_materialized=True,
            actual_finite_N_five_density_integrals_installed=False,actual_active_patch_connected=False)


class OriginalReferenceRestoreFunctions:
    mode='genuine_original_actual_Rsh_reference_and_restoration_to_Rm_functions'
    def __init__(self,dps=500):
        self.upstream=previous.OriginalLongReshapeEndpoint(dps);self.c=c=self.upstream.c
        self.family=self.upstream.family;self.hashes=dict(self.upstream.hashes);self.saved=None;self.owners={}
        for name in (previous.NAME,previous.RECEIPT):
            row=json.loads((HERE/name).read_bytes())
            if not row.get(previous.GATE) or name==previous.RECEIPT and not row.get('all_passed'):
                raise ValueError('Accepted actual finite long-reshape terminal functions required')
            if row['source_family']!=self.family:raise ValueError('Same actual Rsh family required')
            for path,digest in row['input_hashes'].items():fields.previous.bind(self.hashes,path,digest)
            fields.previous.bind(self.hashes,name,sha(name))
            if name==previous.NAME:self.saved=row
        admitted=json.loads((HERE/ADMISSION).read_bytes());checked=json.loads((HERE/ADMISSION_CHECK).read_bytes())
        if not checked.get('all_passed') or not checked.get('original_restoration_integrals_directly_source_bound'):
            raise ValueError('Original full cutoff integral receipt required')
        src=self.upstream.upstream.upstream.upstream.upstream.fields
        for key,want in (('actual_five_defect_family_sha256',self.family),('implicit_source_sha256',src.source),('datum_enclosure_sha256',src.datum)):
            if admitted[key]!=want:raise ValueError('Same restoration family/source/pressure datum required')
        if checked['actual_five_defect_family_sha256']!=self.family or checked['implicit_source_sha256']!=src.source:
            raise ValueError('Original scalar restoration receipt source differs')
        for path,digest in checked['input_hashes'].items():fields.previous.bind(self.hashes,path,digest)
        for name in (ADMISSION,ADMISSION_CHECK,Path(original.__file__).name,PREFIX+'reference_restore_mixed_C4.py',Path(__file__).name):fields.previous.bind(self.hashes,name,sha(name))
        self.provider=RestorationKernelProvider(c,{name:fields.previous.read_interval(c,row) for name,row in admitted['directed_restoration_integrals'].items()})
        self.bindings=source_bindings()

    def owner(self,label):
        if label in self.owners:return self.owners[label]
        if label not in ('0','.5'):raise ValueError('Admitted actual Rsh source frames0,.5 only')
        packet=self.saved['terminal_packets'][label];row=packet['function_evaluation']
        if packet['source_family']!=self.family or packet['source_frame']!=label:
            raise ValueError('Same actual long-reshape terminal inlet required')
        upstream=self.upstream.owner(label);f=upstream.flow;restore=lambda values:[endpoint.restore_row(f,v) for v in values]
        shapes={name:restore(values) for name,values in row['actual_terminal_normalized_six_history_shapes'].items()}
        V=restore(row['physical_velocity_axial_coefficients']['Uz']);P0=restore(row['original_P0_normalized_axial5'])
        for n,v in enumerate(P0):
            if v.record()!=f.scalar(upstream.p0[n]).record():raise ValueError('Canonical analytic P0 tuple identity required')
        if fields.previous.read_interval(self.c,row['source_geometry']['original_fixed_T'])._mpi_!=self.upstream.T._mpi_:
            raise ValueError('Same original fixed T source tuple required')
        self.owners[label]=ActualReferenceRestoreFunctions(f,upstream.Z,upstream.delta,shapes,V,P0,self.upstream.T,self.upstream.logC,self.provider)
        return self.owners[label]

    def evaluate(self,label,chart,coordinate):
        with mp.workdps(self.c.dps+40):
            op=self.owner(label)
            if chart=='reference':value=op.reference(coordinate)
            elif chart=='restoration':value=op.restoration(coordinate)
            elif chart=='postrestore':value=op.postrestore(coordinate)
            else:raise ValueError('Original reference/restoration/postrestore chart required')
        return dict(mode=self.mode,source_family=self.family,source_frame=label,function_evaluation=fields.serialized(value),
            whole_axis_functions_installed=False,original_upstream_micro_function_providers_complete=False,
            no_original_ancestor_producers_or_full_checkers_executed=True)


def run():
    began=time.monotonic();owner=OriginalReferenceRestoreFunctions()
    choices=(('reference',(1,1)),('restoration',(1,2)),('restoration',(1,1)),('postrestore',(-6,1)))
    report=dict(**{GATE:True},source_family=owner.family,original_source_bindings=owner.bindings,
        packets={label:[owner.evaluate(label,chart,q) for chart,q in choices] for label in ('0','.5')},
        **dict.fromkeys(fields.previous.OPEN,False),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
        scope='Actual two-frame Rsh histories through original reference, full-cutoff axial restoration and postrestore to Rm. Signed incoming tails, analytic P0, exact source geometry and factored physical mixed4 rows remain. Finite-N density integrals, active patch, whole Z/global N and real n-recursion stay open.')
    (HERE/NAME).write_text(json.dumps(base.encoded(report),indent=2)+'\n',encoding='utf8')
    print('Actual original Rsh reference/restoration functions and physical mixed4 rows through Rm generated',flush=True);return report


if __name__=='__main__':run()
