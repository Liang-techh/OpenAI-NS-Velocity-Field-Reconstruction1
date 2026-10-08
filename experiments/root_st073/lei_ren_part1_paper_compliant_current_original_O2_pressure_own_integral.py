"""One genuine original O2 pressure own-rate contribution at exact Z=0.

Ordered directed defining-J rectangles cover whole source cells. Actual
radius phase covers are split at periods and passed to the original inverse.
Signed pressure density hulls integrate the full y in[0,1] window at rate0.
This contribution does not reset incoming pressure or install five controls.
"""
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_O2_signed_densities as density

base=density.base;HERE,PREFIX,sha=density.HERE,density.PREFIX,density.sha;ep=density.ep
NAME=PREFIX+'current_original_O2_pressure_own_integral.json'
RECEIPT=PREFIX+'current_original_O2_pressure_own_integral_check.json'
GATE='one_actual_original_O2_midplane_pressure_own_rate_integral_enclosed'


def intersection(c,left,right):
    lo=max(ep(left)[0],ep(right)[0]);hi=min(ep(left)[1],ep(right)[1])
    if lo>hi:raise ArithmeticError('Source identity and directed boxes disagree')
    return c.mpf((lo,hi))


def midplane_identity(templates):
    z,f,*_=templates['inputs'];pZ=templates['pressure_symbols'][1]
    subs={z:0,pZ:0}
    for key in ('p2','V','b'):
        assert all(sy.simplify(expr.subs(subs))==0 for powers,expr in templates['rows'][(key,0)])
    E=templates['rows'][('E',0)]
    assert len(E)==1 and E[0][0]==(0,0,0,0) and sy.simplify(E[0][1].subs(subs)-f)==0
    a,q,psi=sy.symbols('a q psi',real=True);nu=1+2*q*q
    T1=2*q*sy.sin(psi);T2=2*q*q*psi+q*q*sy.sin(2*psi)
    Phi=(psi+T2)/(2*sy.pi*nu)
    assert sy.simplify(a*(Phi-psi/(2*sy.pi))/2-a*q*q*sy.sin(2*psi)/(4*sy.pi*nu))==0
    return dict(passed=True,full_original_inertial_template_p2_V_b_zero_under_exact_pressure_odd_symmetry=True,
        full_original_E_equals_defining_f=True,exact_midplane_T1_T2_A_B_equivalence=True,
        original_pressure_symmetry='P0_Z(0)=0 from accepted original fourteen-stage pressure function',
        B_not_zeroed_with_original_V=True)


def phase_boxes(c,origin,N,left,right):
    result=[]
    for box in origin['true_original_phase_directed_boxes']:
        raw=c.mpf([box['lower'],box['upper']])+c.mpf([ep(left)[0],ep(right)[1]])*N
        lo,hi=ep(raw)
        if hi-lo>=1:return [c.mpf([0,1])]
        for k in range(int(mp.floor(lo)),int(mp.floor(hi))+1):
            lower=max(lo-k,mp.mpf(0));upper=min(hi-k,mp.mpf(1))
            if lower<=upper:result.append(c.mpf((lower,upper)))
    if not result:raise ArithmeticError('True source phase cover missing')
    return result


class OriginalO2PressureOwnIntegral:
    def __init__(self):
        admitted=json.loads((HERE/density.RECEIPT).read_bytes())
        if not admitted.get('all_passed') or not admitted.get(density.GATE):raise ValueError('Accepted actual signed density receipt required')
        self.density=density.OriginalO2SignedDensities();self.owner=self.density.owner;self.family=self.density.family
        if admitted['source_family']!=self.family:raise ValueError('Original pressure contribution family differs')
        self.hashes=dict(self.density.hashes)
        for name,digest in {**admitted['input_hashes'],density.RECEIPT:sha(density.RECEIPT)}.items():
            if sha(name)!=digest:raise ValueError('Original integral source dependency changed: '+name)
            if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Density/integral source dependencies disagree')
            self.hashes[name]=digest
        pressure=base.point.pressure;data=json.loads((HERE/pressure.NAME).read_bytes())
        saved=next(row for row in data['actual_original_pressure_point_queries'] if row['original_Z_exact']=='0')
        odd=saved['normalized_pressure_ordinary_Z_jets'][1]
        if not odd['exact_symmetry_zero'] or not odd['finite_end_and_all_late_source_error']['zero']:
            raise ValueError('Full original P0_Z(0)=0 symmetry required')
        self.identity=midplane_identity(self.owner.inputs.templates)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)

    def integrate(self,*,N=7,cells=128,bits=24):
        if type(cells) is not int or cells<16:raise ValueError('At least16 ordered whole-source cells required')
        if type(bits) is not int or not 8<=bits<=256:raise ValueError('Explicit inverse bits in[8,256] required')
        c=self.owner.ctx;origin=self.owner.radius.evaluate(y=0,N=N);began=time.monotonic()
        with mp.workdps(c.dps+40):
            grid=[c.mpf(i)/cells for i in range(cells+1)];width=c.mpf(1)/cells
            sigma=[base.point.sigma_value_derivative(c,x)[0] for x in grid]
            sigma_cells=[c.mpf((ep(sigma[i])[0],ep(sigma[i+1])[1])) for i in range(cells)]
            prefix=[c.mpf(0)]
            for row in sigma_cells:prefix.append(prefix[-1]+width*row)
            suffix=[c.mpf(0)]*(cells+1)
            for i in reversed(range(cells)):suffix[i]=suffix[i+1]+width*sigma_cells[i]
            J=[intersection(c,prefix[i],c.mpf('.5')-suffix[i]) for i in range(cells+1)]
            J[0]=c.mpf(0);J[-1]=c.mpf('.5')
            definitions=self.owner.inputs.frame.definitions
            logP=c.exp(40)+11;logdelta=-4*logP-30
            logC=c.mpf(mp.mp.make_mpf(self.owner.inputs.frame.selected_logCstar_mpf_tuple))
            logRref=c.ln(110)+10*(logC+logP)
            records=[];integral=c.mpf(0);wraps=0
            for i in range(cells):
                left,right=grid[i],grid[i+1];ys=c.mpf((ep(left)[0],ep(right)[1]))
                Jcell=c.mpf((ep(J[i])[0],ep(J[i+1])[1]))
                f=c.exp(ys/10-c.mpf(3)/5*Jcell)
                a=base.conditioned.clipped(c,c.mpf(4)/5+c.mpf(6)/5*sigma_cells[i],c.mpf(4)/5,2)
                bases=(logP,logdelta,c.mpf(0),c.mpf(0),logRref+ys)
                ledger=dict(directed_small_exponential_tails=0,positive_function_denominator_intersections=0,
                    positive_function_root_intersections=0,directed_independent_log_rescalings=0)
                scalar=lambda value:base.prior.ScaledEnclosure(base.prior.FormalScale(bases),value,ledger)
                eta=base.prior.ScaledEnclosure(base.prior.FormalScale(bases,offset=self.owner.scales.logs['eta']),1,ledger)
                eta_cover=base.conditioned.bounded_value(eta)
                q2=(2+2*eta_cover-a)/(2*a)
                q2=base.conditioned.clipped(c,q2,0,1+2*eta_cover)
                q=scalar(c.sqrt(q2))
                if q.zero:raise ArithmeticError('Positive original eta must not become an exact flat branch')
                roots={key:{(0,0):scalar(value)} for key,value in dict(a=a,E=f,V=0,p2=0,t0=0).items()}
                kernel=base.conditioned.ConditionedPhase(dict(q=q,roots=roots),self.owner.scales.logs['d_star'])
                phases=phase_boxes(c,origin,N,left,right);parts=[]
                for phi in phases:
                    inverse=kernel.evaluate(phi,bits=bits)
                    if inverse['status']!='enclosed':raise ArithmeticError('Original whole-cell inverse needs refinement')
                    psi=2*c.pi*inverse['selected_inverse']['coordinate_interval']
                    # Cancel identical linear angles by the exact source identity
                    # before interval arithmetic; no defining formula changes.
                    A=(roots['a'][(0,0)]*base.current.square(q)).positive_divide(kernel.nu,0)*(c.sin(2*psi)/(4*c.pi))
                    B=-roots['a'][(0,0)]*roots['E'][(0,0)]*q*(c.sin(psi)/(2*c.pi))
                    query=dict(kernel=kernel,roots=roots,ledger=ledger)
                    interpreter=density.BoundDensityGraph(self.density.graph,query,dict(A=A,B_over_Pstar=B),None,N)
                    pressure_density=interpreter.evaluate(self.density.graph['five_signed_increment_rate_roots']['p'])
                    parts.append(pressure_density.finite_interval())
                hull=c.mpf((min(ep(row)[0] for row in parts),max(ep(row)[1] for row in parts)))
                contribution=width*hull;integral+=contribution;wraps+=len(phases)>1
                records.append(dict(exact_y_cell=[str(i)+'/'+str(cells),str(i+1)+'/'+str(cells)],
                    J_source_cell=Jcell,f_source_cell=f,a_source_cell=a,q_source_cell=q.record(),
                    positive_eta_source=eta.record(),true_common_N_phase_boxes=phases,
                    signed_pressure_density_whole_cell=hull,signed_rate0_integral_contribution=contribution,
                    source_caps_or_midpoints_selected_as_field_values=False,
                    entire_cell_source_and_inverse_enclosure_not_point_quadrature=True))
                if (i+1)%max(64,cells//8)==0:print('Actual O2 pressure integral:',cells,i+1,flush=True)
        return dict(source_family=self.family,original_Z_exact='0',explicit_candidate_N=N,
            exact_y_window=['0','1'],ordered_source_cells=cells,inverse_bits=bits,
            original_true_radius_phase_at_y0=origin,
            directed_signed_pressure_own_rate0_integral=integral,
            normalized_units='delta_Mp_over_Pstar_squared; original Pstar² remains a formal physical factor',
            original_own_rate=0,original_inlet_pressure_and_P0_not_reset=True,
            contribution_only_requires_adding_incoming_pressure_history=True,
            whole_source_cells=records,phase_wrap_cell_count=wraps,
            exact_J_endpoint_values=[0,'1/2'],source_J_ordered_rectangles_and_endpoint_symmetry_used=True,
            exact_midplane_p2_symmetry_bound=self.identity,
            positive_original_eta_source_retained_q_cover_not_an_exact_q_point=True,
            no_additional_R_or_Jacobian_applied_to_normalized_density=True,
            full_five_controls_or_Z_functional_terminal_identity_installed=False,
            execution_seconds=time.monotonic()-began)


def run():
    began=time.monotonic();owner=OriginalO2PressureOwnIntegral();records=[]
    for cells in (64,256,2048,8192):records.append(owner.integrate(cells=cells))
    report=dict(**{GATE:True},source_family=owner.family,
        mode='actual_original_O2_midplane_normalized_pressure_own_rate0_integral_enclosures',
        original_full_midplane_template_and_primitive_identity=owner.identity,
        actual_original_pressure_integral_refinements=records,
        one_actual_original_O2_pressure_rate0_integral_enclosed=True,
        point_sampling_not_used_as_whole_cell_integral_proof=True,
        actual_changed_five_moment_integral_evaluated=False,numerical_original_source_point_or_integral_oracle_installed=False,
        actual_five_controls_installed=False,current_whole_N_selected=False,
        **dict.fromkeys(base.point.source.inertial.profiles.loop.OPEN,False),
        input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
        scope='One actual original O2 pressure own-rate0 contribution over y in[0,1] at exact Z0 and diagnostic N7. Whole-cell defining-source/phase/inverse/error enclosures; incoming pressure memory/P0, full five transport, controls, global N and recursion remain open.')
    (HERE/NAME).write_text(json.dumps(base.encoded(report),indent=2)+'\n',encoding='utf8')
    return report


if __name__=='__main__':run()
