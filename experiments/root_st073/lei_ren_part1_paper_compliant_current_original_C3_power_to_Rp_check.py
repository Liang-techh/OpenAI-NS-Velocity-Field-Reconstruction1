"""Independent compact-repaired exit, quiet y4/Z3 flow and native frame checks."""
import json
import math
from pathlib import Path
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_C3_power_to_Rp as current
import lei_ren_part1_paper_compliant_current_original_Rh_C3_continuation_check as rh_check


def encoded(value):return json.loads(json.dumps(current.target.encoded(current.ranges.record(value))))


def symbolic(g,node,bindings,memo=None):
    memo={} if memo is None else memo
    if node in bindings:return bindings[node]
    if node in memo:return memo[node]
    row=g.nodes[node]
    if row['operation']=='analytic_unary' and row['name']=='exprel':
        q=symbolic(g,row['argument'],bindings,memo);result=(s.exp(q)-1)/q
    else:
        # Rh's exact graph interpreter recursively uses its own function.
        # Resolve exprel children to explicit analytic expressions first.
        nested=dict(bindings)
        def prepare(i):
            n=g.nodes[i];op=n['operation']
            if i in nested:return
            if op=='analytic_unary' and n['name']=='exprel':
                nested[i]=symbolic(g,i,bindings,memo);return
            keys={'sum':('arguments',),'product':('arguments',),'negative':('argument',),
                'positive_quotient':('numerator','denominator'),'analytic_unary':('argument',),
                'function_substitution':('expression','variable','value')}.get(op,())
            for key in keys:
                value=n[key]
                for child in value if isinstance(value,list) else [value]:prepare(child)
        prepare(node);result=rh_check.symbolic(g,node,nested)
    memo[node]=result;return result


def compact_exit(field):
    g=field.graph;bf=field.band_functions;f=field.functions;z=s.Symbol('Z');x=s.Symbol('repair_x');bindings={f['original_band_x'].node:x}
    count=0
    for key in current.current.RATES:
        L,D=(s.Function(prefix+key)(x,z) for prefix in ('leading_','correction_'))
        for rows,expr in ((bf['original_leading_power']['leading_histories'][key],L),(bf['partial_correction_histories'][key],D)):
            for j,node in enumerate(rows):bindings[node]=s.diff(expr,z,j)
        certificate=g.nodes[bf['relative_terminal_zero_certificates'][key]]
        assert certificate['operation']=='proved_C3_relative_terminal_zero_identity'
        assert certificate['expression_rows']==bf['relative_terminal_actual'][key]
        assert certificate['zero_rows']==[g.zero.node]*4 and certificate['relative_not_absolute_exterior']
        assert certificate['source_family']==field.identity
        for q,l,d in zip(*(current.target.rows(f[name][key]) for name in
            ('actual_repaired_exit_complete_C3','original_leading_exit_C3','actual_repaired_exit_relative_C3'))):
            assert s.expand(symbolic(g,q.node,bindings)-symbolic(g,l.node,bindings)-symbolic(g,d.node,bindings))==0
            count+=1
    assert count==20
    # Verify the actual compact bump arguments at x=2. This is a genuine
    # support calculation, not an assumption that the finite-N inlet vanishes.
    compact=set();visited=set()
    def walk(i):
        if i in visited:return
        visited.add(i);n=g.nodes[i];op=n['operation']
        if op in ('compact_raw_beta','compact_raw_beta_derivative','compact_raw_beta_second_derivative'):
            compact.add(i);return
        for key in {'sum':('arguments',),'product':('arguments',),'negative':('argument',),
            'positive_quotient':('numerator','denominator'),'analytic_unary':('argument',)}.get(op,()):
            value=n[key]
            for child in value if isinstance(value,list) else [value]:walk(child)
    for profile in bf['profiles_at_band_x']['radial_y_rows']:
        for key in ('F','G','E','V','original_E'):
            for node in profile[key]:walk(node)
    beta_bindings={f['original_band_x'].node:s.Integer(2)};arguments=[]
    for node in compact:
        argument=s.simplify(rh_check.symbolic(g,g.nodes[node]['argument'],beta_bindings))
        assert argument.is_Rational and abs(argument)>1,argument
        arguments.append(str(argument));beta_bindings[node]=s.Integer(0)
    assert len(compact)==9 and set(arguments)=={'8','20','32'}
    known=('exact_rational','bound_variable','sum','product','negative','positive_quotient','analytic_unary')
    for i,n in enumerate(g.nodes):
        if n['operation'] not in known and i not in beta_bindings:beta_bindings[i]=s.Symbol('source_leaf_'+str(i))
    profile_rows=0
    for profile in bf['profiles_at_band_x']['radial_y_rows']:
        for E,E0,V in zip(profile['E'],profile['original_E'],profile['V']):
            assert s.expand(symbolic(g,E,beta_bindings)-symbolic(g,E0,beta_bindings))==0
            assert symbolic(g,V,beta_bindings)==0;profile_rows+=2
    assert profile_rows==24
    return dict(actual_complete_minus_leading_equals_certified_relative_exit_rows=count,
        actual_compact_support_nodes=len(compact),actual_beta_arguments_at_2Rc=sorted(set(arguments)),
        actual_profile_y0_y1_y2_Z0_Z1_Z2_Z3_flat_exit_rows=profile_rows,
        same_source_compact_beta_all_higher_radial_derivatives_flat=True,
        genuine_compact_repair_consumed_instead_of_homogeneous_Rc_shortcut=True,
        complete_absolute_histories_and_independent_P0_not_zeroed=True)


def quiet_calculus(field):
    g=field.graph;f=field.functions;z=s.Symbol('Z');v=s.Symbol('current_C3_post_2Rc_power_s');mu=s.Symbol('mu',positive=True)
    A=s.Function('A2')(z);own={key:s.Function('H2_'+key)(z) for key in current.current.RATES}
    bindings={f['parameters']['mu'].node:mu}
    for j,q in enumerate(current.target.rows(f['actual_repaired_exit_E_C3'])):bindings[q.node]=s.diff(A,z,j)
    for key,row in f['actual_repaired_exit_complete_C3'].items():
        for j,q in enumerate(current.target.rows(row)):bindings[q.node]=s.diff(own[key],z,j)
    alpha=s.Rational(1,2)+mu;E=A*s.exp(-alpha*v);mass=lambda k:(1-s.exp(-k*v))/k
    expected=dict(m=own['m']*s.exp(-v),h=own['h']*s.exp(-3*v/2)+A*(s.exp(-alpha*v)-s.exp(-3*v/2))/(1-mu),
        k=own['k']*s.exp(-3*v/2),e=s.exp(-v)*(own['e']-A*A*mass(2*mu)/2),
        p=own['p']+A*A*mass(1+2*mu)/2)
    flow=f['actual_quiet_flow_C3'];history_rows=profile_rows=0
    for n,row in enumerate(flow['complete_histories_y0_y1_y2_y3_y4_C3']):
        for key,q in row.items():
            for j,ref in enumerate(current.target.rows(q)):
                assert s.simplify(symbolic(g,ref.node,bindings)-s.diff(expected[key],v,n,z,j))==0,(key,n,j)
                history_rows+=1
    for n,row in enumerate(flow['profiles_y0_y1_y2_y3_y4_C3']):
        for j,ref in enumerate(current.target.rows(row['E'])):
            assert s.simplify(symbolic(g,ref.node,bindings)-s.diff(E,v,n,z,j))==0;profile_rows+=1
        assert all(q==g.zero for q in current.target.rows(row['V']))
    drivers=dict(m=0,h=E,k=0,e=-E*E/2,p=E*E/2)
    for key,rate in current.current.RATES.items():
        assert s.simplify(s.diff(expected[key],v)+s.Rational(str(rate))*expected[key]-drivers[key])==0
        assert s.limit(expected[key],v,0)==own[key]
    # Original native complete-state semigroup: splitting the quiet interval
    # preserves absolute memory, including pressure rate zero.
    t=s.Symbol('t',positive=True);replace={A:E,**{own[key]:expected[key] for key in own}}
    for key,expr in expected.items():
        second=expr.subs(v,t).subs(replace,simultaneous=True)
        assert s.simplify(s.expand_power_exp(second-expr.subs(v,v+t)))==0,key
    assert history_rows==100 and profile_rows==20
    return dict(independent_quiet_ordinary_y0_to_y4_Z0_to_Z3_history_rows=history_rows,
        independent_quiet_y0_to_y4_Z0_to_Z3_swirl_rows=profile_rows,
        all_five_actual_ODE_and_inlet_identities=True,all_five_complete_state_semigroup_identities=True,
        relative_zero_propagates_without_zeroing_absolute_pressure=True)


def frame_units_and_geometry(field):
    g=field.graph;f=field.functions;z=s.Symbol('Z');E=s.Function('Ep')(z);H={key:s.Function('Hp_'+key)(z) for key in current.current.RATES};P0=s.Function('P0')(z)
    bindings={}
    for rows,expr in [(f['actual_terminal_Rp_swirl_C3'],E),(f['original_independent_P0_C3'],P0)]+[(f['actual_terminal_Rp_histories_C3'][key],H[key]) for key in H]:
        for j,ref in enumerate(current.target.rows(rows)):bindings[ref.node]=s.diff(expr,z,j)
    expected=dict(u=E,m1=H['m']/E,m2=H['k']/E**2,X=H['h']/E,energy=H['e']/E**2,Mp=H['p'],P0=P0,pressure=P0+H['p'])
    count=0
    for key,wanted in expected.items():
        for j,q in enumerate(current.target.rows(f['actual_C3_pulse_input_frame'][key])):
            assert s.cancel(symbolic(g,q.node,bindings)-s.diff(wanted,z,j))==0,(key,j);count+=1
    # Independently identify normalized common histories with the native paper units.
    R,S,U,M,J,I,En,Cp=s.symbols('R S U M J I En Cp',positive=True)
    norm={E:U/S,H['m']:M/(R*S),H['k']:J/(s.sqrt(2)*R**s.Rational(3,2)*S*S),
        H['h']:I/(s.sqrt(2)*R**s.Rational(3,2)*S),H['e']:En/(R*S*S),H['p']:Cp/(S*S)}
    raw=dict(m1=M/(R*U),m2=J/(s.sqrt(2)*R**s.Rational(3,2)*U*U),X=I/(s.sqrt(2)*R**s.Rational(3,2)*U),energy=En/(R*U*U),Mp=Cp/(S*S))
    for key,wanted in raw.items():assert s.cancel(expected[key].subs(norm,simultaneous=True)-wanted)==0
    mu,offset=s.symbols('mu offset',positive=True);v=s.Symbol('current_C3_post_2Rc_power_s')
    geo_bind={f['parameters']['mu'].node:mu,f['original_Rc_offset'].node:offset}
    length=-60*s.log(mu)-2-s.log(2)
    assert s.simplify(symbolic(g,f['quiet_length'].node,geo_bind)-length)==0
    assert s.simplify(s.expand_power_exp(symbolic(g,f['Rp'].node,geo_bind)-s.exp(offset-60*s.log(mu)-2)))==0
    assert s.simplify(symbolic(g,f['R'].node,geo_bind)-2*s.exp(offset+v))==0
    assert count==32
    return dict(independent_native_pulse_C3_quotient_and_pressure_rows=count,
        original_common_to_native_paper_unit_identities=True,extra_Pstar_division_excluded=True,
        exact_Rc_2Rc_Rp_and_quiet_length_identity=True,selected_native_pulse_inlet_not_yet_identified=True)


def ranges_and_admission(field,raw):
    comparisons=0;c=field.c
    read=lambda rec:current.ranges.LogUpper(c,None if rec['exact_zero'] else current.current.packets.interval(c,rec['log_absolute_upper']))
    zero=current.ranges.LogUpper.constant(c,0);constant=lambda q:current.ranges.LogUpper.constant(c,q)
    add=lambda *q:current.ranges.LogUpper.add(c,q)
    def product(a,b):return [add(*(constant(math.comb(j,k))*a[k]*b[j-k] for k in range(j+1))) for j in range(4)]
    def dominates(saved,expected):
        nonlocal comparisons
        for got,wanted in zip(saved,expected):
            if wanted.log is None:assert got['exact_zero']
            else:
                assert not got['exact_zero']
                assert current.ep(current.current.packets.interval(c,got['log_absolute_upper']))[1]>=current.ep(wanted.log)[1]
            comparisons+=1
    # These positive semigroup masses are analytic inequality witnesses:
    # 0<=1-exp(-r*s)<=1 for r,s>=0. Hence h mass<=2/3, e mass<=1,
    # and the exact E^2 future integral mass<=1/(1+2mu)<=1.
    r,v,mu=s.symbols('r v mu',positive=True)
    assert s.simplify(s.integrate(s.exp(-r*(v-s.Symbol('t'))),(s.Symbol('t'),0,v))-(1-s.exp(-r*v))/r)==0
    assert s.simplify(s.integrate(s.exp(-(1+2*mu)*s.Symbol('t')),(s.Symbol('t'),0,s.oo))-1/(1+2*mu))==0
    for row in raw['actual_four_Z_C3_power_ranges']:
        assert encoded(field.quantitative_range(row['exact_Z_cell']))==row
        ends=tuple(row['exact_Z_cell']);prior=field.rows[ends];length=current.current.packets.interval(c,row['actual_positive_quiet_length']);mu_iv=current.current.packets.interval(c,row['actual_positive_mu'])
        assert current.ep(length)[0]>0 and current.ep(mu_iv)[0]>0 and current.ep(mu_iv)[1]<current.ep(c.mpf(1)/6)[0]
        alpha=c.mpf('.5')+mu_iv
        A=[read(rec) for rec in prior['actual_endpoint_source']['ordinary_terminal_amplitude_C3']];AA=product(A,A)
        inlet={key:[read(rec) for rec in records] for key,records in prior['actual_complete_history_C3_bounds'].items()}
        hist={key:list(q) for key,q in inlet.items()}
        hist['h']=[add(h,constant(c.mpf(2)/3)*a) for h,a in zip(hist['h'],A)]
        for key in ('e','p'):hist[key]=[add(h,constant(c.mpf('.5'))*aa) for h,aa in zip(hist[key],AA)]
        for n in range(5):
            profiles=[constant(alpha**n)*a for a in A]
            dominates(row['actual_whole_quiet_profile_y0_y1_y2_y3_y4_C3_bounds'][n]['E'],profiles)
            dominates(row['actual_whole_quiet_profile_y0_y1_y2_y3_y4_C3_bounds'][n]['V'],[zero]*4)
            for key in hist:dominates(row['actual_whole_quiet_history_y0_y1_y2_y3_y4_C3_bounds'][n][key],hist[key])
            if n==0:base={key:list(q) for key,q in hist.items()}
            if n<4:
                density=dict(m=[zero]*4,h=profiles,k=[zero]*4,
                    e=[constant((2*alpha)**n/2)*aa for aa in AA],p=[constant((2*alpha)**n/2)*aa for aa in AA])
                hist={key:[add(d,constant(c.mpf(rate.numerator)/rate.denominator)*h) for d,h in zip(density[key],hist[key])]
                    for key,rate in current.current.RATES.items()}
        # Independent ordinary-j quotient recurrence from d*q=n, with the
        # actual positive denominator lower and all amplitude derivatives.
        lower=current.current.packets.interval(c,row['actual_terminal_swirl_log_lower'])
        def quotient(n,d,loglower):
            out=[]
            for j in range(4):
                numerator=add(n[j],*(constant(math.comb(j,k))*d[k]*out[j-k] for k in range(1,j+1)))
                out.append(numerator.divide_positive(loglower))
            return out
        P0=[read(rec) for rec in prior['actual_endpoint_source']['ordinary_independent_P0_C3']]
        frame=dict(u=A,m1=quotient(base['m'],A,lower),m2=quotient(base['k'],AA,2*lower),
            X=quotient(base['h'],A,lower),energy=quotient(base['e'],AA,2*lower),Mp=base['p'],P0=P0,
            pressure=[add(a,b) for a,b in zip(P0,base['p'])])
        for key,wanted in frame.items():dominates(row['actual_C3_pulse_frame_bounds'][key],wanted)
        assert row['relative_exit_zero_does_not_zero_complete_histories'] and row['range_caps_not_function_values']
        assert len(row['actual_whole_quiet_history_y0_y1_y2_y3_y4_C3_bounds'])==5
        assert row['positive_swirl_denominators_are_actual_source_functions']
    assert len(raw['actual_four_Z_C3_power_ranges'])==4
    return dict(actual_four_original_axial_cells=4,whole_quiet_profile_history_y4_Z3_and_frame_ranges=True,
        independent_semigroup_mass_and_product_quotient_majorant_comparisons=comparisons,
        independent_positive_own_rate_and_full_square_mass_inequality_witnesses=True,
        exact_source_amplitude_positive_lower_and_full_square_integral_mass_used=True,
        logarithmic_caps_remain_distinct_from_function_values=True)


def run(field=None):
    began=time.monotonic();raw=current.outer.rh.read(current.NAME)
    for name,digest in raw['input_hashes'].items():assert current.sha(name)==digest,name
    assert raw['candidate_actual_C3_power_to_Rp_constructed'] and not any(raw[k] for k in current.GATES+current.OPEN)
    with mp.workdps(540):
        field=field if field is not None else current.CurrentC3PowerToRp(require_checked=False)
        assert not field.acceptance_loaded and field.graph.nodes==raw['exact_graph_nodes']
        assert field.graph.nodes[:len(field.prefix)]==field.prefix
        assert encoded(field.functions)==raw['actual_C3_repaired_exit_power_and_pulse_functions']
        assert raw['source_bindings']==current.source_bindings() and raw['source_family']==field.identity
        assert encoded(field.range_source_binding)==raw['actual_range_source_node_binding']
        assert field.range_source_binding['same_source_recipe_ordinary_endpoint_amplitude_and_normalization']
        assert field.range_source_binding['actual_and_canonical_amplitude_recipe_row_identity']
        assert field.range_source_binding['actual_mu_graph_expression_enclosed_by_every_original_source_cell']
        exitproof=compact_exit(field);calculus=quiet_calculus(field);frame=frame_units_and_geometry(field);bounds=ranges_and_admission(field,raw)
        assert not raw['selected_native_pulse_constructor_consumes_current_C3_frame']
        assert not raw['current_outer_leading_endpoint_band_seed_function_identity_installed']
        result=dict(all_passed=True,**dict.fromkeys(current.GATES,True),source_family=field.identity,
            independent_actual_compact_repaired_exit=exitproof,independent_quiet_mixed_calculus=calculus,
            independent_native_frame_and_geometry=frame,actual_directed_ranges=bounds,**dict.fromkeys(current.OPEN,False),
            selected_native_pulse_constructor_consumes_current_C3_frame=False,
            current_outer_leading_endpoint_band_seed_function_identity_installed=False,
            current_numeric_point_field_oracle_installed=False,global_physical_time_Cartesian_heat_cone_and_temporal_recursion_installed=False,
            input_hashes={**field.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
            execution_seconds=time.monotonic()-began)
        (current.HERE/current.RECEIPT).write_text(json.dumps(current.target.encoded(result),indent=2)+'\n',encoding='utf8')
    print('Actual compact-repaired C3 exit, quiet y4/Z3 semigroup and Rp pulse frame checks passed',flush=True)
    return result


if __name__=='__main__':run()
