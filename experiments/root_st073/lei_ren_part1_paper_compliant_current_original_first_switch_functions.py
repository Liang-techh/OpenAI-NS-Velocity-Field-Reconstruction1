"""Actual first-switch F/V and six histories from the factored R100 inlet.

Known comparison modes continue at their true current radius. Directed
weighted-sigma cell integrals enclose the complete nonlinear source, with
explicit exponential remainders and the original microscopic Jacobian.
"""
from fractions import Fraction
import ast
import json
import math
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_R100_endpoint as endpoint
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import sigma_jets

fields=endpoint.fields;moments=endpoint.moments
HERE,PREFIX,sha=endpoint.HERE,endpoint.PREFIX,endpoint.sha
prior,base,ep=endpoint.prior,endpoint.base,endpoint.ep
NAME=PREFIX+'current_original_first_switch_functions.json'
RECEIPT=PREFIX+'current_original_first_switch_functions_check.json'
GATE='original_actual_first_switch_FV_and_six_history_function_integrals_installed'


def fixed_comparison_denominator_binding():
    """Bind q/core_phi=1/barphi to the accepted defining source, not caps."""
    path=Path(fields.__file__);tree=ast.parse(path.read_text(encoding='utf8'))
    owner=next(n for n in ast.walk(tree) if isinstance(n,ast.ClassDef) and n.name=='OriginalBridgeMacroFunctions')
    method=next(n for n in owner.body if isinstance(n,ast.FunctionDef) and n.name=='owner')
    calls=[n for n in ast.walk(method) if isinstance(n,ast.Call) and ast.unparse(n.func)=='flow.set_sources']
    if len(calls)!=1 or len(calls[0].args)!=7:
        raise ValueError('Original macro defining source attachment changed')
    want=lambda s:ast.dump(ast.parse(s,mode='eval').body)
    if ast.dump(calls[0].args[2])!=want('phi0/phi.truncate(5)') or ast.dump(calls[0].args[4])!=want('phi0'):
        raise ValueError('Same original core phi/comparison quotient identity required')
    from lei_ren_part1_paper_compliant_actual_switch_mixed_C4 import assignment_source_bindings
    original=assignment_source_bindings('switch_signed_integrals','packet',{
        'comparison':'self.bridge.comparison(z,self.bridge.r/100)',
        'phi':"comparison['phi'].truncate(6)",
        'quotient':'actualphi/phi.truncate(5)'})
    frozen=assignment_source_bindings('comparison_point_integrals','macro',{
        '(phi, V)':"(end['phi'],end['V'])"})
    return dict(passed=True,macro_source_q_equals_phi0_over_same_fixed_barphi=True,
        original_R100_fixed_comparison_denominator_assignments=original,
        frozen_macro_constant_phi_and_V_assignments=frozen,
        identity='actual_phi*(q/core_phi)=actual_phi/barphi; q=core_phi/barphi',
        comparison_phi_constant_through_original_frozen_macro=True,
        same_owner_source_family_basis_ledger_and_dependency_hashes_required=True)


def fraction(value):
    if not isinstance(value,tuple) or len(value)!=2 or any(type(v) is not int for v in value):
        raise ValueError('Exact rational first-switch phase pair required')
    a,b=value
    if b<=0 or not 0<=a<=b:raise ValueError('First-switch phase in[0,1] required')
    return Fraction(a,b)


def cutoff_mass(c,left,right,subcells=8):
    """Complete original complement integral, with local Simpson error."""
    a=c.mpf(left);b=c.mpf(right)
    if ep(a)[0]<0 or ep(b)[1]>1 or ep(a)[1]>ep(b)[0]:
        raise ValueError('Ordered first-switch cutoff endpoints required')
    if ep(a)==ep(b):return c.mpf(0)
    if ep(a)==(0,0) and ep(b)==(1,1):return c.mpf('.5')
    step=(b-a)/subcells;total=c.mpf(0);error=c.mpf(0)
    for n in range(subcells):
        l=a+step*n;r=b if n==subcells-1 else a+step*(n+1);mid=(l+r)/2
        box=c.mpf([max(mp.mpf(0),ep(l)[0]),min(mp.mpf(1),ep(r)[1])])
        values=[1-sigma_jets(c,x)[0] for x in (l,mid,r)]
        total+=(r-l)*(values[0]+4*values[1]+values[2])/6
        error+=(r-l)**5*fields.magnitude(c,sigma_jets(c,box)[4])*24/2880
    lo=max(mp.mpf(0),ep(total-error)[0]);hi=min(ep(b-a)[1],ep(total+error)[1])
    if lo>hi:raise ArithmeticError('Original cutoff mass and positive range disagree')
    return c.mpf([lo,hi])


class FirstSwitchFunctions:
    def __init__(self,flow,actual_fields,actual_histories,cells=8):
        if type(cells) is not int or not 2<=cells<=128:
            raise ValueError('Between2 and128 complete integration cells required')
        if not flow.sources_set or ep(flow.hupper)[1]>mp.mpf('.01'):
            raise ValueError('Original source modes and small positive width required')
        if set(actual_fields)!=set(('phi','V')) or set(actual_histories)!=set(moments.RATES):
            raise ValueError('Actual R100 fields and six own histories required')
        for row in (*actual_fields.values(),*actual_histories.values()):
            if len(row)!=6 or any(v.scale.bases is not flow.logs or v.ledger is not flow.ledger for v in row):
                raise ValueError('One R100 source basis/ledger and ordinary Z0..5 required')
        self.flow=flow;self.c=flow.c;self.inlet_fields=actual_fields;self.inlet=actual_histories
        self.cells=cells;self.mass_cache={};self.phase_cache={}
        _,R0,_,_=flow.geometry((1,1));R100=flow.scalar(100)
        # D and drive rows are unweighted modes at the second micro exit.
        # R0^j*100^(p-j) converts them ONCE to current R100 comparison modes.
        def modes(rows,p):return [flow.scale(row,flow.radial_power(R0,j)*flow.radial_power(R100,p-j))
                                  for j,row in enumerate(rows)]
        self.D=modes(flow.d,1)
        self.drives={part:modes(flow.drive[part],2 if part=='swirl' else 1) for part in fields.PARTS}
        core_phi=fields.IntervalTaylor(self.c,[flow.ordinary_cover(v) for v in flow.phi0])
        if ep(core_phi[0])[0]<=0:raise ValueError('Positive original core phi source required')
        self.inv_barphi=flow.multiply(flow.q,flow.jet(core_phi.reciprocal()))
        self.qin=flow.multiply(actual_fields['phi'],self.inv_barphi)
        self.G0={part:flow.add(*rows) for part,rows in self.drives.items()}
        self.force0={part:flow.multiply(self.qin,row) for part,row in self.G0.items()}
        self.scales={part:flow.factor((0,int(part=='pressure'),int(part=='swirl'),0,0)) for part in fields.PARTS}

    def cv(self,q):return self.c.mpf(q.numerator)/q.denominator

    def box(self,left,right):return self.c.mpf([max(mp.mpf(0),ep(self.cv(left))[0]),min(mp.mpf(1),ep(self.cv(right))[1])])

    def scalar_series(self,a,t,average=False,minus_one=False):
        """exp(a*h*t) or its average; h stays in the formal source basis."""
        f=self.flow;c=self.c;u=f.h*(c.mpf(a)*t);power=f.scalar(1)
        result=f.scalar(0) if minus_one else power
        for k in range(1,4):
            power=power*u
            result+=power*(c.mpf(1)/math.factorial(k+1 if average else k))
        upper=abs(c.mpf(a))*f.hupper*fields.magnitude(c,t)
        error=f.h*f.h*f.h*f.h*(abs(c.mpf(a))**4*fields.magnitude(c,t)**4*c.exp(upper)/math.factorial(5 if average else 4))
        return result+prior.ScaledEnclosure(error.scale,error.coefficient*c.mpf([-1,1]),f.ledger),error

    def angular(self,t):
        f=self.flow;c=self.c;terms=[];errors=[]
        for j,row in enumerate(self.D):
            avg,error=self.scalar_series(1-j,t,average=True)
            terms.append(f.scale(row,f.h*f.h*(-t/2)*avg))
            errors.append(f.norm(row)*f.h*f.h*(fields.magnitude(c,t)/2)*error)
        ell=f.add(*terms);norm=f.norm(ell);upper=f.small_upper(norm)
        if ep(upper)[1]>mp.mpf('.5'):raise ValueError('Complete first-switch log jet norm must be <=1/2')
        eta=norm*norm*norm*(c.exp(upper)/6)
        expminus=f.add(ell,f.scale(f.multiply(ell,ell),c.mpf('.5')),f.error_rows(eta))
        phi_delta=f.multiply(self.inlet_fields['phi'],expminus)
        phi=f.add(self.inlet_fields['phi'],phi_delta)
        return dict(ell=ell,expminus=expminus,phi_delta=phi_delta,phi=phi,
                    nonzero_exponential_error_norm=eta,nonzero_radial_mode_integral_error_norms=errors)

    def force_delta(self,t,angular):
        f=self.flow;deltas={};G={};Gdelta={}
        for part,p in (('hydro',1),('pressure',1),('swirl',2)):
            gd=f.add(*[f.scale(row,self.scalar_series(p-j,t,minus_one=True)[0])
                       for j,row in enumerate(self.drives[part])])
            Gdelta[part]=gd;G[part]=f.add(self.G0[part],gd)
            # Keep the shared R100 baseline before enclosure. No subtraction
            # of independent phi/G boxes fabricates a width-independent error.
            change=f.add(gd,f.multiply(self.G0[part],angular['expminus']),
                         f.multiply(gd,angular['expminus']))
            deltas[part]=f.multiply(self.qin,change)
        return deltas,G,Gdelta

    def mass(self,left,right):
        key=(left,right)
        if key not in self.mass_cache:
            self.mass_cache[key]=cutoff_mass(self.c,self.cv(left),self.cv(right))
        return self.mass_cache[key]

    def evaluate(self,phase):
        q=fraction(phase)
        if q in self.phase_cache:return self.phase_cache[q]
        f=self.flow;c=self.c;s=self.cv(q);h=f.h
        at=self.angular(s);end_delta,end_G,_=self.force_delta(s,at)
        cuts=[q*Fraction(n,self.cells) for n in range(self.cells+1)] if q else [Fraction(0)]
        cells=[];prefix=c.mpf(0);delta_prefix={part:[f.scalar(0)]*6 for part in fields.PARTS}
        for left,right in zip(cuts[:-1],cuts[1:]):
            t=self.box(left,right);ang=self.angular(t);dforce,G,_=self.force_delta(t,ang)
            mass=self.mass(left,right);next_prefix=prefix+mass
            if right==1:next_prefix=c.mpf('.5')
            Pcover=c.mpf([max(mp.mpf(0),ep(prefix)[0]),min(ep(s)[1],ep(next_prefix)[1])])
            # The complete actual V at every point in this cell consists of
            # constant source times its true sigma prefix plus a signed
            # integrated source difference and a positive partial-cell mass.
            Vparts={part:f.scale(f.add(f.scale(self.force0[part],Pcover),delta_prefix[part],
                f.scale(dforce[part],c.mpf([0,ep(mass)[1]]))),-h*h*self.scales[part])
                for part in fields.PARTS}
            V=f.add(self.inlet_fields['V'],*Vparts.values())
            cells.append(dict(left=left,right=right,t=t,phi=ang['phi'],V=V,
                weighted_cutoff_mass=mass,source_difference=dforce,original_current_drives=G))
            for part in fields.PARTS:delta_prefix[part]=f.add(delta_prefix[part],f.scale(dforce[part],mass))
            prefix=next_prefix
        if q==0:prefix=c.mpf(0)
        if q==1:prefix=c.mpf('.5')
        Vparts={part:f.scale(f.add(f.scale(self.force0[part],prefix),delta_prefix[part]),-h*h*self.scales[part])
                for part in fields.PARTS}
        deltaV=f.add(*Vparts.values());V=f.add(self.inlet_fields['V'],deltaV)
        if q==0:V=self.inlet_fields['V'];at['phi']=self.inlet_fields['phi']
        source_end=dict(H=f.scale(at['phi'],2),M=V,K=f.scale(f.multiply(at['phi'],V),2),
            A=f.multiply(V,V),B=f.multiply(at['phi'],at['phi']),C=f.multiply(at['phi'],at['phi']))
        histories={};increments={};derivatives={};evidence={}
        for name,rate in moments.RATES.items():
            dm,decay_error=self.scalar_series(-rate,s,minus_one=True)
            incoming_delta=f.scale(self.inlet[name],dm);pieces=[];kernel_errors=[]
            for cell in cells:
                phi,v=cell['phi'],cell['V']
                sources=dict(H=f.scale(phi,2),M=v,K=f.scale(f.multiply(phi,v),2),
                             A=f.multiply(v,v),B=f.multiply(phi,phi),C=f.multiply(phi,phi))
                width=self.cv(cell['right']-cell['left'])
                backward,err1=self.scalar_series(-rate,s-self.cv(cell['right']))
                avg,err2=self.scalar_series(-rate,width,average=True)
                mass=h*width*backward*avg
                pieces.append(f.scale(sources[name],mass));kernel_errors.append((err1,err2))
            integral=f.add(*pieces) if pieces else [f.scalar(0)]*6
            change=f.add(incoming_delta,integral);actual=f.add(self.inlet[name],change)
            if q==0:actual=self.inlet[name]
            histories[name]=actual;increments[name]=change
            derivatives[name]=f.scale(f.add(source_end[name],f.scale(actual,-rate)),h)
            evidence[name]=dict(rate=rate,retained_actual_incoming_decay_minus_one=dm,
                incoming_delta=incoming_delta,signed_complete_positive_kernel_integral=integral,
                complete_incoming_decay_error=decay_error,complete_kernel_errors=kernel_errors,
                microscopic_jacobian_hb_applied_once=True,signed_source_products_before_enclosure=True)
        Dend=f.add(*[f.scale(row,self.scalar_series(1-j,s)[0]) for j,row in enumerate(self.D)])
        phi_phase=f.scale(f.multiply(Dend,at['phi']),-h*h*c.mpf('.5'))
        complement=1-sigma_jets(c,s)[0]
        V_phase=f.add(*[f.scale(f.add(self.force0[part],end_delta[part]),-h*h*self.scales[part]*complement)
                        for part in fields.PARTS])
        displacement,radial_error=self.scalar_series(1,s,minus_one=True)
        result=dict(actual_fields=dict(phi=at['phi'],V=V),actual_six_histories=histories,
            signed_field_increments=dict(phi=at['phi_delta'],V=deltaV,parts=Vparts),
            signed_six_history_increments=increments,
            actual_phase_ODE_rows=dict(phi=phi_phase,V=V_phase,**derivatives),
            actual_first_switch_log_phi_over_R100=at['ell'],
            angular_full_exponential_error_norm=at['nonzero_exponential_error_norm'],
            angular_full_radial_mode_integral_error_norms=at['nonzero_radial_mode_integral_error_norms'],
            history_integral_evidence=evidence,weighted_sigma_prefix_mass=prefix,
            complete_weighted_sigma_cells=[dict(exact_left=[x['left'].numerator,x['left'].denominator],
                exact_right=[x['right'].numerator,x['right'].denominator],
                mass=x['weighted_cutoff_mass'],signed_source_difference=x['source_difference']) for x in cells],
            geometry=dict(phase=[q.numerator,q.denominator],radius=f.scalar(100)+displacement*100,
                signed_nonzero_radius_minus100=displacement*100,log_radius_over100=h*s,
                complete_radial_exponential_error=radial_error,exact_source='R=100*exp(hb*s); dy=hb*ds'),
            ordinary_Z_orders=list(range(6)),true_sigma_function_used=True,
            comparison_direction_not_redefined_by_actual_moments=True,
            complete_finite_integrals_not_leading_half_mass_reset=True,
            source_cell_range_quadrature_not_point_fits=True)
        self.phase_cache[q]=result;return result


class OriginalFirstSwitchFunctions:
    mode='genuine_original_actual_first_switch_function_integrals_with_factored_inlets'
    def __init__(self,dps=500,cells=8):
        self.upstream=endpoint.OriginalR100Endpoint(dps);self.c=self.upstream.c
        self.family=self.upstream.family;self.hashes=dict(self.upstream.hashes);self.cells=cells;self.owners={}
        for name in (endpoint.NAME,endpoint.RECEIPT):
            row=json.loads((HERE/name).read_bytes())
            if not row.get(endpoint.GATE) or name==endpoint.RECEIPT and not row.get('all_passed'):
                raise ValueError('Accepted factored actual R100 source required')
            for path,digest in row['input_hashes'].items():fields.previous.bind(self.hashes,path,digest)
            fields.previous.bind(self.hashes,name,sha(name))
        from lei_ren_part1_paper_compliant_first_switch_leading import switch_control_source_bridge
        self.bindings=switch_control_source_bridge()
        self.denominator_binding=fixed_comparison_denominator_binding()
        for path,digest in self.bindings['input_hashes'].items():fields.previous.bind(self.hashes,path,digest)
        for name in ('first_switch_leading','flat_pulse_derivatives','switch_signed_integrals'):
            path=PREFIX+name+'.py';fields.previous.bind(self.hashes,path,sha(path))
        fields.previous.bind(self.hashes,Path(__file__).name,sha(Path(__file__).name))

    def owner(self,label):
        if label not in self.owners:
            flow,proof,value=self.upstream.owner(label)
            self.owners[label]=FirstSwitchFunctions(flow,value['normalized_actual_fields'],
                value['normalized_actual_own_six_moments'],self.cells)
        return self.owners[label]

    def evaluate(self,label,phase):
        with mp.workdps(self.c.dps+40):value=self.owner(label).evaluate(phase)
        return dict(mode=self.mode,source_frame=label,source_family=self.family,
            implicit_source_sha256=self.upstream.upstream.fields.source,
            datum_enclosure_sha256=self.upstream.upstream.fields.datum,
            function_evaluation=fields.serialized(value),actual_second_switch_function_installed=False,
            actual_micro_function_provider_installed=False,whole_axis_functions_installed=False)


def run():
    began=time.monotonic();owner=OriginalFirstSwitchFunctions()
    result=dict(**{GATE:True},source_family=owner.family,source_control_bindings=owner.bindings,
        fixed_comparison_denominator_source_binding=owner.denominator_binding,
        packets={label:[owner.evaluate(label,q) for q in ((0,1),(1,2),(1,1))] for label in ('0','.5')},
        **dict.fromkeys(fields.previous.OPEN,False),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
        scope='Actual original first-switch phi/V and six history functions with complete weighted-sigma cell integrals, ordinary Z0..5, defining phase ODEs and factored R100 inlets. Second switch, whole bridge/Z, R110 and global reconstruction remain open.')
    (HERE/NAME).write_text(json.dumps(base.encoded(result),indent=2)+'\n',encoding='utf8')
    print('Original actual first-switch F/V and six histories integrated',flush=True);return result


if __name__=='__main__':run()
