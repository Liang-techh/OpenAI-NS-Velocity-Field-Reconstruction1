"""Current function attachment and independent active/transition/flat kernels."""
import gzip
import json
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_generic_loop_function_sources as source
from lei_ren_part1_paper_compliant_current_generic_shear_loop import GenericLoopScales,GenericShearLoop,flat_step

packets=source.packets


class FixtureEvaluator:
    """Explicit modest original-source oracle; saved coefficient covers unused."""
    def __init__(self,record,values,Vvalues,loop,phase,N):
        self.record=record;self.nodes=record['function_graph_nodes'];self.values=values;self.Vvalues=Vvalues
        self.loop=loop;self.c=loop.ctx;self.phase=phase;self.N=N;self.memo={}
        self.t1=record['roots']['T1_at_free_angle'];self.t2=record['roots']['T2_at_free_angle']
    def eval(self,index,angle=None):
        key=(index,None if angle is None else angle._mpf_)
        if key in self.memo:return self.memo[key]
        n=self.nodes[index];op=n['operation'];c=self.c
        ev=lambda q:self.eval(q,angle)
        if op=='constant':v=c.mpf(n['value'])
        elif op=='mathematical_pi':v=c.pi
        elif op=='source_derivative':
            name=n['name'];v=self.Vvalues[name] if name.startswith('V_') else self.values[name]
        elif op=='sum':v=sum((ev(i) for i in n['arguments']),c.mpf(0))
        elif op=='negative':v=-ev(n['argument'])
        elif op=='product':
            v=c.mpf(1)
            for i in n['arguments']:v*=ev(i)
        elif op=='positive_function_quotient':
            den=ev(n['denominator'])
            if den<=0:raise ArithmeticError('Active fixture functional denominator not positive')
            v=ev(n['numerator'])/den
        elif op=='function_variable':v=self.phase if n['name']=='phi' else angle
        elif op=='shared_positive_integer_parameter':v=c.mpf(self.N)
        elif op=='original_positive_log_scale':v=self.loop.scales.eta if n['name']=='eta' else self.loop.scales.d_star
        elif op=='flat_zero_branch':
            v=c.mpf(0) if ev(n['Delta'])>=ev(n['eta']) else ev(n['active_body'])
        elif op=='analytic_unary':
            x=ev(n['argument']);name=n['name']
            if name=='original_flat_sigma':v=flat_step(c,x)
            elif name=='original_flat_sigma_prime':v=c.diff(lambda z:flat_step(c,z),x)
            elif name=='positive_sqrt':
                if x<=0:raise ArithmeticError('Inactive square root was evaluated')
                v=c.sqrt(x)
            elif name=='absolute':v=abs(x)
            else:v=getattr(c,name)(x)
        elif op=='Poisson_half_angle_squared':
            r=ev(n['signed_r']);a=ev(n['angle']);v=(c.sin(a/2) if r>=0 else c.cos(a/2))**2
        elif op=='correlated_cos_minus_r':
            r=ev(n['signed_r']);a=ev(n['angle']);rho=ev(n['one_minus_abs_r'])
            v=rho-2*c.sin(a/2)**2 if r>=0 else 2*c.cos(a/2)**2-rho
        elif op=='angle_integral':
            upper=ev(n['upper_angle'])
            # Reuse the independent exact antiderivatives for the base t,t².
            # Slow derivative integrands are independently integrated here.
            if index==self.t1:v=self.loop.integrals(upper)[0]
            elif index==self.t2:v=self.loop.integrals(upper)[1]
            elif upper==0:v=c.mpf(0)
            else:v=c.quad(lambda psi:self.eval(n['integrand'],psi),[0,upper/2,upper])
        elif op=='monotone_phase_inverse':
            if ev(n['flat_Delta'])>=ev(n['flat_eta']):v=2*c.pi*self.phase
            else:v=self.loop.angle_at_phase(self.phase)
            residual=self.eval(n['phase_function'],v)-self.phase
            if abs(residual)>c.mpf('1e-60'):raise ArithmeticError('Graph phase binding differs from independently inverted phase')
        elif op=='substitute_inverse_angle':v=self.eval(n['body'],ev(n['inverse_angle']))
        else:raise ValueError('Unsupported fixture functional operation: '+op)
        self.memo[key]=v;return v


def independent_fixture():
    c=mp.mp.clone();c.dps=80;iv=MPIntervalContext();iv.dps=100
    comparisons=0;branches=[];fd_allowance=c.mpf('2e-7')
    for kind in ('active','cutoff_transition','flat'):
        scales=GenericLoopScales(a_min='.7',margin_min='1',boundary_kappa_excess_min='.1',t0_abs_max='1',p1_abs_max='10',p2_abs_max='3',dps=80)
        c=scales.ctx
        a0=c.mpf('.8') if kind=='active' else 2+scales.eta/2 if kind=='cutoff_transition' else c.mpf(3)
        ay=c.mpf('.02') if kind!='cutoff_transition' else scales.eta/50
        az=c.mpf('.015') if kind!='cutoff_transition' else scales.eta/80
        def fields(y,Z):
            E=(1+c.mpf('.04')*Z)*c.exp((1-a0)*y/2-ay*y*y/4-az*Z*y/2)
            V=c.mpf('.3')+c.mpf('.02')*Z
            V+=c.mpf('.0005')*y*y+c.mpf('.001')*y*Z if kind=='cutoff_transition' else c.mpf('.05')*y+c.mpf('.01')*y*Z
            return E,V
        def values(y,Z):
            E,V=fields(y,Z);aa=a0+ay*y+az*Z
            Vy=c.diff(lambda t:fields(t,Z)[1],y);bb=2*Vy/E
            return dict(E=E,a=aa,b=bb,p1=c.mpf(8),p2=c.mpf('.5')+c.mpf('.01')*y+c.mpf('.02')*Z,
                t0=-bb/aa,kappa_minus2=aa+bb*bb/aa-2)
        nodes=[];roots={};oracle={}
        for name in ('a','b','p1','p2','E','t0','kappa_minus2'):
            roots[name]={}
            for order in (source.ZERO,*source.FIRST):
                label=name+'_y%d_Z%d'%order;roots[name]['y%d_Z%d'%order]=len(nodes)
                nodes.append(dict(operation='source_derivative',name=label))
                oracle[label]=values(0,0)[name] if order==source.ZERO else c.diff(lambda t:values(t,0)[name],0) if order==(1,0) else c.diff(lambda t:values(0,t)[name],0)
        algebra=packets.FactoredAlgebra(iv,(iv.mpf(0),iv.mpf(0),iv.mpf(0),iv.mpf(0)),[])
        Vvalues={};Vrows={}
        for order in (source.ZERO,*source.FIRST):
            value=fields(0,0)[1] if order==source.ZERO else c.diff(lambda t:fields(t,0)[1],0) if order==(1,0) else c.diff(lambda t:fields(0,t)[1],0)
            Vvalues['V_y%d_Z%d'%order]=value;Vrows[order]=algebra.lift(iv.mpf(str(value)),0)
        inputs=dict(chart='fixture',source_family={'fixture':True},
            jet_expression_dag=dict(nodes=nodes,roots=roots,source_derivative_leaves={}),
            original_O3_small_excess_definition=None,original_O3_positive_mu=None,source_provenance={},original_radius_source={'fixture':True})
        record=source.build(inputs,Vrows,dict(eta={},d_star={}),dict(positive_conditioning=dict(log_lambda_positive_lower='fixture')))
        base=values(0,0);loop=GenericShearLoop(scales,a=base['a'],b=base['b'],p1=base['p1'],p2=base['p2'],Utheta=base['E'])
        h=c.mpf('1e-5');N=11
        def reference(y,Z,phi):
            par=values(y,Z);ee,vv=fields(y,Z)
            owner=GenericShearLoop(scales,a=par['a'],b=par['b'],p1=par['p1'],p2=par['p2'],Utheta=ee)
            point=owner.evaluate(phi);r=ee*c.expm1(point['A']/N);u=point['B']/N
            flux=dict(m=u,h=r,k=vv*r+ee*u+r*u,e=2*vv*u+u*u-ee*r-r*r/2,p=ee*r+r*r/2)
            return dict(E_N=ee*c.exp(point['A']/N),V_N=vv+u,delta_E=r,delta_V=u,A=point['A'],B_over_Pstar=point['B']),flux
        for phi in (c.mpf('.31'),):
            evaluator=FixtureEvaluator(record,oracle,Vvalues,loop,phi,N);ref,flux=reference(0,0,phi)
            for name in ref:
                actual=evaluator.eval(record['roots'][name])
                if abs(actual-ref[name])>c.mpf('1e-60')*(1+abs(ref[name])):raise ArithmeticError('Independent loop function graph value differs: '+kind+' '+name)
                comparisons+=1
            for name,node in record['five_signed_increment_rate_roots'].items():
                if abs(evaluator.eval(node)-flux[name])>c.mpf('1e-60')*(1+abs(flux[name])):raise ArithmeticError('Changed signed source rate differs')
                comparisons+=1
            yp,fp=reference(h,0,phi+N*h);ym,fm=reference(-h,0,phi-N*h)
            ypp,fpp=reference(2*h,0,phi+2*N*h);ymm,fmm=reference(-2*h,0,phi-2*N*h)
            zp,gp=reference(0,h,phi);zm,gm=reference(0,-h,phi)
            zpp,gpp=reference(0,2*h,phi);zmm,gmm=reference(0,-2*h,phi)
            # Fourth-order centered stencil resolves the fast phase without
            # relaxing the analytic first-derivative comparison allowance.
            derivative=lambda pp,p,m,mm,name:(-pp[name]+8*p[name]-8*m[name]+mm[name])/(12*h)
            for key,p,m,pp,mm in (('y',yp,ym,ypp,ymm),('Z',zp,zm,zpp,zmm)):
                for component,name in (('theta','E_N'),('axial','V_N')):
                    actual=evaluator.eval(record['changed_velocity_first_ordinary_derivatives'][component][key])
                    expected=derivative(pp,p,m,mm,name)
                    if abs(actual-expected)>fd_allowance*(1+abs(expected)):raise ArithmeticError('Actual total-y/Z changed profile derivative differs: '+kind+' '+key+' '+name+' actual='+c.nstr(actual,22)+' expected='+c.nstr(expected,22)+' diff='+c.nstr(actual-expected,15))
                    comparisons+=1
            for key,p,m,pp,mm in (('y',fp,fm,fpp,fmm),('Z',gp,gm,gpp,gmm)):
                for name,node in record['five_signed_increment_rate_first_derivatives'][key].items():
                    actual=evaluator.eval(node);expected=derivative(pp,p,m,mm,name)
                    if abs(actual-expected)>fd_allowance*(1+abs(expected)):raise ArithmeticError('Actual first changed signed source derivative differs: '+kind+' '+name+' '+key)
                    comparisons+=1
            for name,expected in (('a_N',1-2*evaluator.eval(record['changed_velocity_first_ordinary_derivatives']['theta']['y'])/ref['E_N']),
                ('b_N',2*evaluator.eval(record['changed_velocity_first_ordinary_derivatives']['axial']['y'])/ref['E_N'])):
                if abs(evaluator.eval(record['roots'][name])-expected)>c.mpf('1e-55')*(1+abs(expected)):raise ArithmeticError('Modulated a_N/b_N differs from actual derivative')
                comparisons+=1
            branches.append(kind)
    return dict(comparisons=comparisons,branches=branches,
        direct_point_allowance='1e-60 relative/absolute, modest fixtures only',
        first_finite_difference_step='fourth-order centered stencil,1e-5 step,2e-7 allowance; phi changes byN*step for total y',
        actual_source_values_or_eta_not_materialized=True,
        independent_existing_exact_antiderivative_and_inverse_oracle_used=True,
        first_slow_integrands_independently_quadrature_evaluated=True)


def run():
    manifest=json.loads((source.HERE/source.NAME).read_bytes())
    for name,digest in manifest['input_hashes'].items():
        if source.sha(name)!=digest:raise ValueError('Changed actual function source prerequisite: '+name)
    records=json.loads(gzip.decompress((source.HERE/source.VIEWS).read_bytes()))
    provider=source.CurrentLoopFunctionSources();count=0;terms=0
    if set(records)!=set(provider.bounds) or len(records)!=17:raise ArithmeticError('Actual function-source cover incomplete')
    for chart,record in records.items():
        actual=provider.chart(chart)
        if packets.encode(actual)!=record:raise ValueError('Actual analytic function graph differs: '+chart)
        if set(actual['five_signed_increment_rate_roots'])!=set(('m','h','k','e','p')):raise ArithmeticError('Full changed five source rates missing')
        if actual['own_defect_transport_contract']['quiet_gap_rule'].endswith('preserves memory') is not True:raise ArithmeticError('Quiet pressure memory reset')
        operations={q['operation'] for q in actual['function_graph_nodes']}
        for op in ('flat_zero_branch','angle_integral','monotone_phase_inverse','substitute_inverse_angle','shared_positive_integer_parameter'):
            if op not in operations:raise ArithmeticError('Defining analytic loop operation missing')
        for node in actual['function_graph_nodes']:
            if node['operation']=='flat_zero_branch' and not node['active_body_not_evaluated_on_flat_branch']:raise ArithmeticError('Flat branch eager evaluation would divide by inactive gamma')
            if node['operation']=='shared_positive_integer_parameter' and (node['name']!='common_N' or not node['same_for_all_charts']):raise ArithmeticError('Chartwise frequency reset')
        freq=actual['actual_source_frequency_majorants'];count+=len(actual['function_graph_nodes'])
        for key,row in freq['five_increment_rate_log_polynomials'].items():
            if any(term['N_power']>=0 for term in row['terms']):raise ArithmeticError('Value source lost1/N majorant')
            terms+=len(row['terms'])
        for row in freq['five_increment_rate_first_log_polynomials']['Z'].values():
            if any(term['N_power']>=0 for term in row['terms']):raise ArithmeticError('Z source lost slow derivative order')
        if chart=='O3_power' and not all(row['exact_zero'] for row in freq['five_increment_rate_log_polynomials'].values()):raise ArithmeticError('Original power flat loop source must be zero')
        if chart!='O3_power' and not any(term['N_power']==0 for row in freq['five_increment_rate_first_log_polynomials']['y'].values() for term in row['terms']):raise ArithmeticError('Fast y source derivatives incorrectly forced to1/N')
    if any(manifest[k] for k in source.OPEN) or manifest['actual_changed_five_moment_transport_integrated']:raise ArithmeticError('Function definitions promoted missing transport/physical stages')
    if packets.encode(provider.theorem)!=manifest['exact_loop_modulation_and_signed_increment_theorem']:raise ValueError('Current function-source theorem differs')
    fixture=independent_fixture()
    result=dict(all_passed=True,source_family=provider.family,**{source.GATE:True},**dict.fromkeys(source.OPEN,False),
        actual_source_charts=17,function_graph_nodes=count,five_value_source_frequency_terms=terms,
        exact_loop_modulation_density_and_Poisson_bindings=len(provider.theorem['exact_identities']),
        independent_active_transition_flat_fixture=fixture,
        flat_active_bodies_are_lazy=True,nonzero_original_V_and_full_signed_cross_terms_retained=True,
        single_N_and_global_phase_preserved=True,quiet_gap_defect_and_pressure_memory_preserved=True,
        original_power_zero_increment_source_admitted=True,
        signed_current_point_loop_or_inverse_jets_installed=False,actual_changed_five_moment_transport_integrated=False,
        current_whole_N_selected=False,full_velocity_mixed4_or_stress_mixed3_admitted=False,
        source_graph_ancestor_constructors_called=False,
        input_hashes={**manifest['input_hashes'],source.NAME:source.sha(source.NAME),Path(__file__).name:source.sha(Path(__file__).name)})
    (source.HERE/source.RECEIPT).write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('Actual loop and five signed source functions PASS:17 charts; total-y/Z independent kernels PASS',flush=True)
    return result


if __name__=='__main__':run()
