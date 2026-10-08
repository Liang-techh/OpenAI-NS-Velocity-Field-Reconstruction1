"""Adaptive, source-bound point phase for all17 original radius charts.

Exact selected dyadic periods are reduced before analytic arithmetic.
Actual microscopic widths stay in the source expression and are bounded
by signed numerical errors, not rounded to zero. No field value is selected.
"""
from fractions import Fraction
import json
from pathlib import Path
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as sy
import lei_ren_part1_paper_compliant_current_native_spatial_phase as phase
import lei_ren_part1_paper_compliant_current_native_Rc_function_transport as functions
import lei_ren_part1_paper_compliant_current_original_O2_source_parameter_frame as parameters

HERE,PREFIX,sha=phase.HERE,phase.PREFIX,phase.sha;ep=phase.ep
NAME=PREFIX+'current_original_all_chart_point_phase.json'
RECEIPT=PREFIX+'current_original_all_chart_point_phase_check.json'
GATE='actual_original17_chart_point_phase_with_adaptive_analytic_and_signed_microscopic_errors_installed'


def encode(v):
    if hasattr(v,'_mpi_'):
        lo,hi=ep(v)
        return dict(lower=mp.nstr(lo,85),upper=mp.nstr(hi,85),lower_exact_mpf_tuple=list(lo._mpf_),upper_exact_mpf_tuple=list(hi._mpf_))
    if hasattr(v,'_mpf_'):return dict(decimal=mp.nstr(v,85),exact_mpf_tuple=list(v._mpf_))
    if isinstance(v,dict):return {k:encode(x) for k,x in v.items()}
    if isinstance(v,(tuple,list)):return [encode(x) for x in v]
    return v


def rational(value):
    q=parameters.exact_rational(value)
    if max(abs(int(q.p)).bit_length(),int(q.q).bit_length())>8192:
        raise ValueError('Exact source coordinate representation limited to8192 bits')
    return q


def periodic_add(c,boxes,value):
    """Directed circular Minkowski sum; never hull across the 0/1 seam."""
    pieces=[]
    for box in boxes:
        projection=phase.ordinary_mod_one(c,box+value)
        if projection['full_period']:return dict(full_period=True,boxes=[c.mpf((0,1))])
        pieces.extend(projection['boxes'])
    ordered=sorted((ep(box) for box in pieces),key=lambda pair:pair[0]);merged=[]
    for lo,hi in ordered:
        if merged and lo<=merged[-1][1]:merged[-1]=(merged[-1][0],max(hi,merged[-1][1]))
        else:merged.append((lo,hi))
    return dict(full_period=len(merged)==1 and merged[0][0]<=0 and merged[0][1]>=1,
        boxes=[c.mpf(pair) for pair in merged])


class OriginalAllChartPointPhase:
    def __init__(self,binder):
        if type(binder) is not phase.NativeSpatialPhase:raise TypeError('Genuine accepted original radius binder required')
        self.binder=binder;self.family=binder.family;self.hashes=dict(binder.service.hashes);self.cache={}
        for module in (phase,functions):
            accepted=json.loads((HERE/module.RECEIPT).read_bytes())
            if not accepted['all_passed'] or not accepted[module.GATE] or accepted['source_family']!=self.family:
                raise ValueError('Same-family accepted original radius/function theorem required')
            for name,digest in {**accepted['input_hashes'],module.RECEIPT:sha(module.RECEIPT)}.items():
                if sha(name)!=digest:raise ValueError('Original phase source bytes changed: '+name)
                if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Source closures disagree: '+name)
                self.hashes[name]=digest
        if rational(binder.seed.params.Md)!=40:raise ValueError('Original Md=40 required')
        filename='lei_ren_part1_paper_logarithmic_outer_parameters.py'
        self.analytic_recipe_proof=[parameters.source_assignment(filename,'__init__',name,value,'LogarithmicOuterParameters')
            for name,value in (('self.logPstar','c.exp(self.md)+11'),('self.log_mu',"c.ln(c.mpf('.001'))-4*self.logPstar"),('self.Tw','-60*self.log_mu'))]
        self.hashes[filename]=sha(filename)
        self.maps,self.symbols,self.x=functions.exact_radius_maps()
        if set(self.maps)!=set(phase.native.DOMAINS):raise ValueError('All17 original radius charts required')
        check=MPIntervalContext();check.dps=max(300,binder.ctx.dps+60,binder.seed.params.ctx.dps+60)
        with mp.workdps(check.dps+40):
            self.live_analytic_parameter_binding={}
            for name,value in (('logP',binder.seed.logP),('Tw',binder.seed.params.Tw)):
                computed=self.analytic(self.symbols[name],check);lo,hi=ep(check.mpf(value));clo,chi=ep(computed)
                if not lo<=clo<=chi<=hi:raise ValueError('Live original analytic parameter disagrees: '+name)
                self.live_analytic_parameter_binding[name]=dict(original_native_enclosure=check.mpf(value),
                    defining_recipe_refinement=computed,same_source_refinement_contained=True)
        self.source_graph_sha256=sha(functions.sources.VIEWS)
        for name in (functions.sources.VIEWS,Path(functions.__file__).name,Path(parameters.__file__).name,Path(__file__).name):self.hashes[name]=sha(name)

    def coordinate(self,chart,value,c):
        if chart not in self.maps:raise ValueError('Original radius chart required')
        special=None
        if isinstance(value,dict):
            if set(value)=={'selected_sc_multiple'} and chart=='bridge_first':
                q=rational(value['selected_sc_multiple']);expr=self.symbols['sc']*q;special='selected_sc_multiple'
            elif set(value)=={'original_power_offset'} and chart=='O3_power':
                q=rational(value['original_power_offset']);expr=q/self.symbols['Tw'];special='original_power_offset'
            elif set(value)=={'original_patch_log_offset'} and chart=='actual_patch':
                q=rational(value['original_patch_log_offset'])
                if not 0<=q<=1:raise ValueError('Original patch logarithmic offset in[0,1] required')
                expr=sy.exp(q);special='original_patch_log_offset'
            else:raise ValueError('Unsupported original source coordinate expression')
        else:q=rational(value);expr=q
        coord=self.analytic(expr,c)
        lo,hi=phase.native.DOMAINS[chart]
        upper=c.exp(1) if hi=='e' else c.mpf(hi)
        # The exp(offset) patch convention proves this exact endpoint domain
        # algebraically; all other requests need their directed domain proof.
        if special!='original_patch_log_offset' and (ep(coord)[0]<ep(c.mpf(lo))[1] or ep(coord)[1]>ep(upper)[0]):
            raise ValueError('Original native chart coordinate outside its certified domain')
        return expr,coord,dict(kind=special or 'exact_rational',exact_coefficient=str(q),source_expression=sy.srepr(expr))

    def analytic(self,expr,c):
        if expr.is_Rational:return c.mpf(int(expr.p))/int(expr.q)
        if expr==self.symbols['Md']:return c.mpf(40)
        if expr==self.symbols['sc']:
            if ep(self.binder.sc)[0]!=ep(self.binder.sc)[1]:raise ValueError('Source-selected singleton sc required')
            return c.mpf(self.binder.sc)
        if expr==self.symbols['logP']:return c.exp(40)+11
        if expr==self.symbols['Tw']:return 60*c.ln(1000)+240*(c.exp(40)+11)
        if expr==sy.E:return c.exp(1)
        if expr.func is sy.Add:return sum((self.analytic(t,c) for t in expr.args),c.mpf(0))
        if expr.func is sy.Mul:
            out=c.mpf(1)
            for t in expr.args:out*=self.analytic(t,c)
            return out
        if expr.func is sy.Pow and expr.args[1].is_Integer:return self.analytic(expr.args[0],c)**int(expr.args[1])
        if expr.func is sy.exp:return c.exp(self.analytic(expr.args[0],c))
        if expr.func is sy.log:
            arg=self.analytic(expr.args[0],c)
            if ep(arg)[0]<=0:raise ArithmeticError('Original analytic logarithm lost positivity')
            return c.ln(arg)
        raise ValueError('Unbound original analytic radius operation: '+str(expr))

    def decomposition(self,chart,expr):
        original=sy.expand(self.maps[chart].subs(self.x,expr));terms={}
        for name in ('logP','logC','T','Tw','hbB','hbS'):terms[name]=sy.expand(original).coeff(self.symbols[name])
        constant=original.subs({self.symbols[name]:0 for name in terms})
        if sy.expand(original-constant-sum(self.symbols[name]*q for name,q in terms.items()))!=0:
            raise ArithmeticError('Original affine radius decomposition failed')
        if self.symbols['sc'] in constant.free_symbols:raise ArithmeticError('Microscopic origin escaped its width factor')
        return original,constant,terms

    def evaluate(self,*,chart,coordinate,N,decimal_digits=80):
        if type(N) is not int or N<1 or N.bit_length()>4096:raise ValueError('Explicit positive candidate integer N up to4096 bits required')
        if type(decimal_digits) is not int or not 40<=decimal_digits<=1000:raise ValueError('Phase digits in[40,1000] required')
        key=(chart,json.dumps(coordinate,sort_keys=True) if isinstance(coordinate,dict) else str(coordinate),N,decimal_digits)
        if key in self.cache:return self.cache[key]
        digitsN=(N.bit_length()*30103+99999)//100000
        c=MPIntervalContext();c.dps=max(260,decimal_digits+digitsN+60)
        with mp.workdps(c.dps+40):
            expr,coord,coord_record=self.coordinate(chart,coordinate,c)
            original,constant,terms=self.decomposition(chart,expr)
            projection=dict(full_period=False,boxes=[c.mpf(0)]);proofs=[]
            # Reduce each ordinary analytic component separately before adding
            # it to the selected enormous dyadic constants.
            analytic=[('constant',self.analytic(constant,c))]
            for name in ('logP','Tw'):
                if terms[name]!=0:analytic.append((name,self.analytic(terms[name]*self.symbols[name],c)))
            for name,value in analytic:
                wrapped=phase.ordinary_mod_one(c,value*N)
                proofs.append(dict(component=name,method='adaptive_directed_analytic_modulus',original_component=value,projection=wrapped))
                if wrapped['full_period']:raise ArithmeticError('Analytic phase enclosure needs more precision')
                pieces=[]
                for component in wrapped['boxes']:
                    added=periodic_add(c,projection['boxes'],component)
                    if added['full_period']:raise ArithmeticError('Analytic phase sum needs more precision')
                    pieces.extend(added['boxes'])
                # Adding zero merges the directed circular cover without taking
                # a hull through a periodic seam.
                projection=periodic_add(c,pieces,c.mpf(0))
            for name in ('logC','T'):
                if terms[name]==0:continue
                coefficient=terms[name]
                if coefficient.is_Rational is not True:raise ValueError('Exact rational multiplier required for selected dyadic modulus')
                multiplier=Fraction(int(coefficient.p),int(coefficient.q))*N
                selected=c.mpf(self.binder.fixed[name]);value,proof=phase.binary_mod_one(c,selected,multiplier)
                projection=periodic_add(c,projection['boxes'],value)
                proofs.append(dict(component=name,method='exact_selected_dyadic_modulus',**proof))
            micro=[];budget=c.mpf(0);epsilon=c.mpf(10)**-(decimal_digits+20)
            for name,kind,logwidth in (('hbB','bridge',self.binder.loghB),('hbS','switch',self.binder.loghS)):
                if terms[name]==0:continue
                coefficient=self.analytic(terms[name],c);lo,hi=ep(coefficient);maxabs=max(abs(lo),abs(hi))
                if maxabs==0:
                    enclosure=c.mpf(0);logupper=None;zero=True
                else:
                    logupper=c.ln(N)+c.mpf(logwidth)+c.ln(c.mpf(maxabs));zero=False
                    if ep(logupper)[1]>=ep(c.ln(epsilon))[0]:
                        raise ArithmeticError('Original microscopic phase exceeds the chosen error budget; refine explicitly')
                    enclosure=c.mpf((0,ep(epsilon)[1])) if lo>=0 else c.mpf((-ep(epsilon)[1],0)) if hi<=0 else c.mpf((-ep(epsilon)[1],ep(epsilon)[1]))
                budget+=enclosure
                micro.append(dict(original_width_kind=kind,original_positive_width_log=c.mpf(logwidth),
                    exact_coefficient_source_expression=sy.srepr(terms[name]),coefficient_enclosure=coefficient,
                    original_N_width_coefficient_phase_log_upper=logupper,numerical_phase_error_enclosure=enclosure,
                    exact_zero_coefficient=zero,actual_width_not_zeroed=True,actual_width_not_materialized=True,
                    directed_absolute_error_comparison_passed=True))
            projection=periodic_add(c,projection['boxes'],budget)
            if projection['full_period']:raise ArithmeticError('Analytic phase enclosure needs more precision; no selected phase fallback')
            result=dict(source_family=self.family,chart=chart,explicit_candidate_N=N,coordinate=coord_record,
                directed_native_coordinate=coord,source_graph_sha256=self.source_graph_sha256,
                exact_original_radius_minus_inlet_expression=sy.srepr(original),
                exact_source_phase='fractional_part(N*original_log_radius_minus_inlet)',
                actual_phase_directed_boxes=projection['boxes'],periodic_projection=projection,
                adaptive_interval_digits=c.dps,requested_phase_digits=decimal_digits,
                original_regular_component_proofs=proofs,original_signed_microscopic_phase_components=micro,
                total_microscopic_phase_error_enclosure=budget,source_selected_sc_same_as_original=True,
                selected_logC_T_integer_cycles_not_materialized=True,analytic_logP_Tw_not_selected_as_values=True,
                original_radius_Z_derivative_exact_zero=True,original_width_Jacobians_unchanged=True,
                source_phase_not_an_independent_free_angle=True,source_field_values_not_selected_from_caps=True,
                candidate_N_not_global_frequency_admission=True,current_whole_N_selected=False,
                full_original_source_point_or_integral_oracle_installed=False)
        self.cache[key]=result;return result

    def phase_for_source(self,row,built,*,coordinate,N,decimal_digits=80):
        if type(N) is not int or N<160:raise ValueError('Original graph frequency requires integer N>=160')
        if built['source_family']!=self.family or built['source_graph_sha256']!=self.source_graph_sha256:
            raise ValueError('Same original source family and graph required')
        if row.get('operation')!='original_function_graph' or row.get('graph_sha256')!=self.source_graph_sha256:
            raise ValueError('Bound original function-source reference required')
        g=built['graph']
        if type(g) is not functions.FunctionTransportGraph or not any(row is node for node in g.nodes):
            raise ValueError('Issued original function-source row required')
        g.ids([built['N'],*built['parameters'].values()]);chart=row['chart']
        if row.get('graph_file')!=functions.sources.VIEWS or row['shared_N']!=built['N'].node:
            raise ValueError('Same original shared-N radius source required')
        phase_node=g.nodes[row['phase']]
        if phase_node.get('operation')!='analytic_unary' or phase_node.get('name')!='fractional_part':
            raise ValueError('Original fractional radius phase required')
        Nsymbol=sy.Symbol('shared_N',positive=True,integer=True);memo={}
        bound={q.node:self.symbols[name] for name,q in built['parameters'].items() if name in self.symbols}
        def symbolic(index):
            if index in memo:return memo[index]
            item=g.nodes[index];op=item['operation']
            if index in bound:out=bound[index]
            elif index==built['N'].node:
                if op!='shared_positive_integer' or item.get('lower')!=160:raise ValueError('Original N>=160 binding required')
                out=Nsymbol
            elif op=='exact_rational':out=sy.Rational(item['numerator'],item['denominator'])
            elif op=='bound_variable':out=sy.Symbol(item['name'],real=True)
            elif op=='sum':out=sy.Add(*(symbolic(q) for q in item['arguments']))
            elif op=='product':out=sy.Mul(*(symbolic(q) for q in item['arguments']))
            elif op=='negative':out=-symbolic(item['argument'])
            elif op=='positive_quotient':out=symbolic(item['numerator'])/symbolic(item['denominator'])
            elif op=='analytic_unary' and item['name'] in ('exp','log'):
                out=(sy.exp if item['name']=='exp' else sy.log)(symbolic(item['argument']))
            else:raise ValueError('Unbound source radius graph operation')
            memo[index]=out;return out
        expected=Nsymbol*self.maps[chart].subs(self.x,symbolic(row['coordinate']))
        if sy.simplify(symbolic(phase_node['argument'])-expected)!=0:
            raise ValueError('Same original shared-N radius phase required')
        if g.nodes[row['coordinate']]['operation']!='bound_variable':
            requested,_,_=self.coordinate(chart,coordinate,self.binder.ctx)
            if sy.simplify(requested-symbolic(row['coordinate']))!=0:
                raise ValueError('Fixed original source coordinate must be preserved')
        return self.evaluate(chart=chart,coordinate=coordinate,N=N,decimal_digits=decimal_digits)


@phase.native.inlet.source_precision
def run(*,return_live=False):
    began=time.monotonic();bridge,construction=phase.native.inlet.native_bridge_owner()
    with phase.native.inlet.CheckedSourceRuntime() as runtime:
        binder=phase.NativeSpatialPhase(phase.native.NativeGenericSourcePackets(bridge));owner=OriginalAllChartPointPhase(binder)
        points=dict(bridge_first='.1337',bridge_second='1.831',bridge_macro='.537',switch_first='.537',switch_second='1.337',
            switch_power='.537',reshape='.537',inner_reference='.1337',axial_restore='.537',restore_buffer='-6.337',
            actual_patch='1.337',Rh_reference='-2.337',O2_slope='.537',O2_axial='.1337',O2_buffer='5.337',
            O3_slope_mu='.537',O3_power={'original_power_offset':'2'})
        records={chart:owner.evaluate(chart=chart,coordinate=x,N=257) for chart,x in points.items()}
        large={chart:owner.evaluate(chart=chart,coordinate=points[chart],N=(1<<1200)+37) for chart in ('inner_reference','O2_axial','O3_power')}
        endpoints=dict(zero_inlet=owner.evaluate(chart='bridge_first',coordinate={'selected_sc_multiple':'1/2'},N=257),
            flat_collar_end=owner.evaluate(chart='bridge_first',coordinate={'selected_sc_multiple':'3/4'},N=257),
            patch_exit=owner.evaluate(chart='actual_patch',coordinate={'original_patch_log_offset':'1'},N=257),
            reference_inlet=owner.evaluate(chart='Rh_reference',coordinate='-5',N=257))
    result=dict(**{GATE:True},source_family=owner.family,source_graph_sha256=owner.source_graph_sha256,
        original_all17_true_point_phase_records=records,large_explicit_N_point_phase_records=large,
        original_source_expression_endpoint_records=endpoints,original_analytic_parameter_recipe=owner.analytic_recipe_proof,
        live_original_analytic_parameter_binding=owner.live_analytic_parameter_binding,
        original_bridge_construction=construction,original_source_runtime=runtime.record(),
        full_original_source_point_or_integral_oracle_installed=False,actual_five_controls_installed=False,
        actual_terminal_Z_function_closure_installed=False,current_whole_N_selected=False,
        **dict.fromkeys(phase.packets.OPEN,False),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
        scope='Fresh original17-chart arbitrary exact-point radius phase with selected dyadic periods reduced exactly, adaptive defining analytic parameter arithmetic and signed actual microscopic-source error budgets. Source leaf values/integrals/controls, globally compatible N, recursion and corrected NS remain open.')
    (HERE/NAME).write_bytes(json.dumps(encode(result),indent=2).encode()+b'\n')
    print('Original17-chart adaptive point phase with retained source-width errors generated',flush=True)
    return (result,owner) if return_live else result


if __name__=='__main__':run()
