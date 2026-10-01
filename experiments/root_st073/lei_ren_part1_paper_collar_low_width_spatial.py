"""Cellwise enclosures for both collar width coefficients on all s in [0,2].

This is a continuous spatial coefficient enclosure at Z=.3, conditional on
the finite input and retained pressure powers. Higher width orders are open.
"""
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_collar_interval_receipt_reader import read_inlet
from lei_ren_part1_paper_coherent_pressure_error_transfer import accepted_profile
from lei_ren_part1_paper_schedule_endpoint_enclosures import ScheduleEndpointEnclosures,endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_high_order_preheat_integrals import HighOrderPreheatIntegrals,switch_taylor
from lei_ren_part1_paper_interval_taylor import IntervalTaylor,integrate_symmetric
from lei_ren_part1_paper_collar_first_width_interval_endpoint import first_width_coefficients
from lei_ren_part1_paper_collar_second_width_moments import second_width_spatial_coefficients
from lei_ren_part1_paper_collar_width_physical_endpoint import physical_moment_coefficients


def weighted_switch_grid(e,panels=128,order=12):
    iv=e.iv;h=iv.mpf(1)/panels;radius=h/2;edge=e.sigma_interval(h)
    total=iv.mpf([0,endpoints(h*h*edge/2)[1]]);grid=[iv.mpf(0),total]
    for i in range(1,panels-1):
        center=(iv.mpf(i)+iv.mpf('.5'))/panels
        cell=iv.mpf([endpoints(iv.mpf(i)/panels)[0],endpoints(iv.mpf(i+1)/panels)[1]])
        def density(x):return IntervalTaylor.variable(iv,x,order)*switch_taylor(iv,x,order)
        total+=integrate_symmetric(density(center),density(cell),radius);grid.append(total)
    total+=h-h*h/2-iv.mpf([0,endpoints(h*edge)[1]]);grid.append(total)
    return grid


def primitive_nodes(e,cells_per_unit=32,panels=128):
    if not isinstance(cells_per_unit,int) or cells_per_unit<4 or cells_per_unit%4:
        raise ValueError('Require at least four cells per unit, divisible by four')
    if panels%cells_per_unit:raise ValueError('Primitive grid must align with spatial cell endpoints')
    iv=e.iv;J_sigma=HighOrderPreheatIntegrals(e).primitive_grid(panels,12)
    Q_sigma=weighted_switch_grid(e,panels,12);nodes=[]
    for i in range(2*cells_per_unit+1):
        s=iv.mpf(i)/cells_per_unit
        if i<=cells_per_unit:
            n=i*(panels//cells_per_unit);J=s-J_sigma[n];K=s*s/2-s*J_sigma[n]+Q_sigma[n]
        else:J=iv.mpf('.5');K=Q_sigma[-1]+(s-1)/2
        jl,jh=endpoints(J);kl,kh=endpoints(K)
        J=iv.mpf([max(mp.mpf(0),jl),min(mp.mpf('.5'),jh)])
        K=iv.mpf([max(mp.mpf(0),kl),min(endpoints(s*s/2)[1],endpoints(s/2)[1],kh)])
        nodes.append((s,J,K))
    return nodes


def aggregate(jet):
    return {slot:dict(nominal=getattr(jet,slot).evaluate(1,0).nominal,
                      pressure_error=getattr(jet,slot).evaluate(1,0).difference)
            for slot in ('value','tangent','second')}


def run(cells_per_unit=32):
    with mp.workdps(520):
        inlet,provenance=read_inlet();iv=inlet['ctx'];profile,_=accepted_profile()
        e=ScheduleEndpointEnclosures(profile.schedule);nodes=primitive_nodes(e,cells_per_unit)
        rows=[];maxima={};coverage_end=mp.mpf(0)
        for index,(left,right) in enumerate(zip(nodes,nodes[1:])):
            sl=endpoints(left[0])[0];sr=endpoints(right[0])[1]
            if sl>coverage_end:raise AssertionError('Spatial coverage gap')
            coverage_end=sr;s=iv.mpf([sl,sr])
            J=iv.mpf([endpoints(left[1])[0],endpoints(right[1])[1]])
            K=iv.mpf([endpoints(left[2])[0],endpoints(right[2])[1]])
            first=first_width_coefficients(inlet,s=s,switch_integral=J)
            second=second_width_spatial_coefficients(inlet,s=s,switch_primitive=J,integrated_primitive=K)
            data={str(order):{name:aggregate(jet) for name,jet in values.items()}
                  for order,values in ((1,first),(2,second))}
            moments={str(order):{name:aggregate(jet) for name,jet in physical_moment_coefficients(inlet,values).items()}
                     for order,values in ((1,first),(2,second))}
            for order,values in data.items():
                for name,slots in values.items():
                    for slot,value in slots.items():
                        lo,hi=endpoints(value['pressure_error']);key=f'{order}:{name}:{slot}'
                        maxima[key]=max(maxima.get(key,mp.mpf(0)),abs(lo),abs(hi))
            rows.append(dict(index=index,s=s,J=J,K=K,coefficient_bounds=data,moment_coefficient_bounds=moments))
        if coverage_end!=2:raise AssertionError('Spatial coverage incomplete')
        # Independent switch quadrature checks directed primitives at interior
        # nodes; monotonicity, not sampled interpolation, covers whole cells.
        from lei_ren_part1_paper_axial_primitive import _sigma_mp
        checks=[]
        for text in ('.25','.5','.75','1','2'):
            s=mp.mpf(text);edge=min(s,mp.mpf(1))
            with mp.workdps(90):
                J=mp.quad(lambda x:_sigma_mp(1-x),[0,edge/2,edge])
                K=mp.quad(lambda x:(s-x)*_sigma_mp(1-x),[0,edge/2,edge])
            node=nodes[int(s*cells_per_unit)]
            for val,box in ((J,node[1]),(K,node[2])):
                lo,hi=endpoints(box)
                if not lo<=val<=hi:raise AssertionError('Independent primitive outside interval')
            checks.append(dict(s=text,J=J,K=K,contained=True))
        report=dict(source=provenance,Z='.3',s_domain=['0','2'],cells=len(rows),
            whole_s_cells_enclosed=True,analytic_width_orders=[1,2],
            all_eight_states_enclosed=True,axial_derivative_slots=[0,1,2],
            all_five_unnormalized_moment_coefficients_enclosed=True,
            collar_scope='Section9.25 inner exit collar; not exterior heat collar',
            retained_pressure_parameter_order=9,aggregation_pressure_value=1,
            pressure_error_maxima=maxima,independent_primitive_checks=checks,spatial_cells=rows,
            axis_generation_or_source_errors_enclosed=False,omitted_pressure_orders_enclosed=False,
            omitted_width_orders_enclosed=False,full_collar_ODE_enclosed=False,
            global_axial_domain_enclosed=False,full_spatial_divergence_certified=False,temporal_recursion=False)
        Path(__file__).with_suffix('.json').write_text(json.dumps(encode(report),indent=2)+'\n')
        print('whole collar coefficient cells',len(rows),'s=[0,2]',flush=True)
        for key in ('1:u:value','2:u:value','1:g:value','2:g:value'):
            print(key,'conditional pressure error',mp.nstr(maxima[key],15),flush=True)


if __name__=='__main__':run()
