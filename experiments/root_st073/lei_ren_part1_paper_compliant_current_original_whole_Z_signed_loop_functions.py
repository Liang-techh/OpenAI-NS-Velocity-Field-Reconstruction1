"""Live original signed loop C1 functions with current whole-Z radius phase.

Function evaluations are directed enclosures of the defining inverse, not
selected values from support bounds. Only Z/phase derivatives are exported.
The local integral interface retains its supplied correction-only incoming;
the full inner-to-Rc source oracle and terminal controls remain open.
"""
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
from types import FunctionType, SimpleNamespace
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_original_whole_Z_O3_Rc_finite_N as current
import lei_ren_part1_paper_compliant_current_native_signed_u_phase_cover as signed

phase=current.reference.phase
radius=phase.radius_phase
HERE,PREFIX,sha,bind,ep=current.HERE,current.PREFIX,current.sha,current.bind,current.ep
CELLS,RATES,C0,Z,OPEN=current.CELLS,current.RATES,current.C0,current.Z,current.OPEN
NAME=PREFIX+'current_original_whole_Z_signed_loop_functions.json.gz'
RECEIPT=PREFIX+'current_original_whole_Z_signed_loop_functions_check.json'
GATE='current_original_whole_Z_live_signed_inverse_Z_functions_and_local_integrals_installed'
serialized,encode=current.serialized,current.encode
OUTPUTS=phase.OUTPUTS
CHARTS=('Rh_reference','O2_slope','O2_axial','O2_buffer','O3_transition','O3_power')


def exact(value):
    if not isinstance(value,tuple) or len(value)!=2 or any(type(v) is not int for v in value) or value[1]<=0:
        raise ValueError('Exact rational source coordinate required')
    return Fraction(*value)


def union(rows):
    if not rows:raise ValueError('Nonempty source function alternatives required')
    return signed.density.local.same_source_union(rows)


def branch_Z_functions(source,qrows,dstar_log,phi,branch):
    """Reuse the admitted Z-only body with a conditional signed-u geometry."""
    namespace=dict(phase.Z_FIRST.__globals__)
    original_phase=namespace['phase']
    namespace['phase']=SimpleNamespace(**{**vars(original_phase),
        'ConditionedPhase':lambda src,log:signed.BranchConditionedPhase(src,log,branch)})
    fn=FunctionType(phase.Z_FIRST.__code__,namespace,'current_original_signed_Z_function',phase.Z_FIRST.__defaults__)
    got=fn(source,qrows,dstar_log,phi)
    if got['values'] is not None and set(got['values'])!=set(OUTPUTS):
        raise ValueError('Only genuine Z/phase primitive function slots may be exported')
    got['record']['original_Z_only_body_and_derivative_algebra_unchanged']=True
    got['record']['conditional_signed_u_geometry']=branch['name']
    got['record']['source_roots_and_q_Z_not_conditioned_or_selected']=True
    return got


class CurrentOuterRadiusPhase:
    def __init__(self,op,identity,N):
        self.op,self.identity,self.N=op,identity,N
        self.base=phase.RmRadiusPhase(op,identity['actual_five_defect_family_sha256'],identity)
        self.cache={}

    def offset(self,chart,left,right=None):
        right=left if right is None else right
        l,r=exact(left),exact(right)
        domains={'Rh_reference':(-5,0),'O2_slope':(0,1),'O2_axial':(0,1),'O2_buffer':(0,11),
            'O3_transition':(0,1),'O3_power':(0,2)}
        if chart not in domains or not domains[chart][0]<=l<=r<=domains[chart][1]:
            raise ValueError('Ordered actual outer source coordinate/domain required')
        key=(chart,l,r)
        if key in self.cache:return self.cache[key]
        c=MPIntervalContext();c.dps=max(260,80+(self.N.bit_length()*30103+99999)//100000+60)
        with mp.workdps(c.dps+40):
            cv=lambda q:c.mpf(q.numerator)/q.denominator
            t=c.mpf((ep(cv(l))[0],ep(cv(r))[1]))
            offset={'Rh_reference':lambda:6+t,'O2_slope':lambda:6+t,
                'O2_axial':lambda:6+c.exp(40*t),'O2_buffer':lambda:6+c.exp(40)+t,
                'O3_transition':lambda:17+c.exp(40)+t,'O3_power':lambda:18+c.exp(40)+t}[chart]()
            delta=signed.first.spatial.ordinary_mod_one(c,self.N*offset)
            base=self.base.point((1,1),self.N,decimal_digits=80)
            pieces=[];full=delta['full_period']
            if full:pieces=[c.mpf((0,1))]
            else:
                for before in base['phase_boxes']:
                    for part in delta['boxes']:
                        got=radius.periodic_add(c,[c.mpf(before)],part)
                        full=full or got['full_period'];pieces.extend(got['boxes'])
                if full:pieces=[c.mpf((0,1))]
                else:pieces=radius.periodic_add(c,pieces,c.mpf(0))['boxes']
            result=dict(source_identity=self.identity,candidate_N=self.N,actual_chart=chart,
                exact_left=list(left),exact_right=list(right),actual_Rm_phase_binding=base,
                actual_original_log_offset=offset,offset_periodic_projection=delta,
                actual_phase_boxes=[self.op.c.mpf(part) for part in pieces],full_period=full,
                actual_phase_Z_exact_zero=True,actual_same_source_Rm_radius=self.op.Rm_factor,
                global_phase_not_restarted=True,phase_not_supplied_as_free_angle=True,
                adaptive_analytic_interval_digits=c.dps,
                exact_original_global_phase='frac(N*(logRm+original_outer_offset-logRa-hb*s_c/2))')
        self.cache[key]=result;return result


class WholeZSignedLoopFunctions:
    def __init__(self,dps=500):
        self.source=current.WholeZO3RcFiniteN(dps);self.c=self.source.c
        self.identity,self.N=self.source.identity,self.source.N;self.hashes=dict(self.source.hashes)
        checked=json.loads((HERE/current.RECEIPT).read_bytes())
        if not checked.get('all_passed') or not checked.get(current.GATE) or checked['source_family']!=self.identity \
                or checked['candidate_N']!=self.N:
            raise ValueError('Checked current whole-Z original Rc source required')
        for name,digest in checked['input_hashes'].items():bind(self.hashes,name,digest)
        bind(self.hashes,current.RECEIPT,sha(current.RECEIPT))
        for module in (current,phase,phase.first,phase.first.phase,signed,signed.density.local,radius):
            name=Path(module.__file__).name;bind(self.hashes,name,sha(name))
        bind(self.hashes,Path(__file__).name,sha(Path(__file__).name))
        self.radius={};self.cache={}

    def owner(self,ends):return self.source.owner(ends)

    def source_query(self,ends,chart,left,right=None):
        right=left if right is None else right
        if chart not in CHARTS:raise ValueError('Issued original outer source chart required')
        if chart.startswith('O3_'):return self.source.query(ends,chart.removeprefix('O3_'),left,right)
        axial=self.source.source
        if chart in ('O2_axial','O2_buffer'):return axial.query(ends,chart.removeprefix('O2_'),left,right)
        slope=axial.source
        if chart=='O2_slope':return slope.query(ends,left,right)
        return slope.source.query(ends,left,right)

    def radius_query(self,ends,chart,left,right=None):
        key=tuple(ends)
        if key not in self.radius:self.radius[key]=CurrentOuterRadiusPhase(self.owner(ends),self.identity,self.N)
        return self.radius[key].offset(chart,left,right)

    def query(self,ends,chart,left,right=None):
        right=left if right is None else right;key=(tuple(ends),chart,left,right)
        if key in self.cache:return self.cache[key]
        geometry=self.radius_query(ends,chart,left,right);packet=self.source_query(ends,chart,left,right)
        op=self.owner(ends);f=op.flow
        roots=packet['original_roots'];qrows=packet['original_q_C0_Z']
        source=dict(q=qrows[C0],roots=roots)
        current.reference.parameters.same_source(f,[row for rows in roots.values() for row in rows.values()]+list(qrows.values()))
        if packet['exact_common_P0_axial5'] is not op.P0 or not packet['actual_phase_Z_exact_zero']:
            raise ValueError('Same independent P0 and source-radius Z identity required')
        original_u,branches,empty=signed.signed_u_branches(source,self.source.dstar_log)
        alternatives=[]
        for phi in geometry['actual_phase_boxes']:
            for branch in branches:
                got=branch_Z_functions(source,qrows,self.source.dstar_log,phi,branch)
                if got['values'] is None:raise ArithmeticError('Refine actual source/phase for signed inverse: '+got['record']['status'])
                E=packet['original_generic_source']['common_velocity_E_axial5']
                V=packet['original_generic_source']['common_velocity_V_axial5']
                A=got['values']['A'];logA=None if A.zero else ep(A.record()['log_absolute_upper'])[1]
                if logA is not None and logA>=ep(self.c.ln(self.N))[0]:raise ArithmeticError('Actual signed function A/N budget needs refinement')
                nonlinear=current.reference.phase.densities.density_Z_kernels(E[0],E[1],V[0],V[1],got['values'],self.N)
                alternatives.append(dict(actual_phase_box=phi,conditional_signed_u_branch=branch['name'],
                    original_inverse_and_Z_function_proof=got['record'],actual_signed_primitive_C0_Z_phi=got['values'],
                    actual_nonlinear_density_C0_Z=nonlinear,actual_A_N_admitted=True))
        primitives={name:union([row['actual_signed_primitive_C0_Z_phi'][name] for row in alternatives]) for name in OUTPUTS}
        density={part:{name:union([row['actual_nonlinear_density_C0_Z'][part][name] for row in alternatives]) for name in RATES}
            for part in ('kernels','Z_derivatives')}
        E,V=packet['original_generic_source']['common_velocity_E_axial5'],packet['original_generic_source']['common_velocity_V_axial5']
        result=dict(source_identity=self.identity,exact_Z_cell=list(ends),candidate_N=self.N,actual_chart=chart,
            exact_left=list(left),exact_right=list(right),exact_common_P0_axial5=op.P0,
            actual_source_geometry=geometry,actual_original_source_packet=packet,
            actual_signed_E_V_C0_Z=dict(E=list(E[:2]),V=list(V[:2])),
            actual_signed_inverse_function_alternatives=alternatives,original_u=original_u,
            conditional_u_branches_proved_empty=empty,actual_signed_primitive_C0_Z_phi=primitives,
            actual_signed_nonlinear_density_C0_Z=density,
            original_inverse_functions_evaluated_not_all_u_support_caps=True,
            conditional_branches_hulled_after_nonlinear_evaluation_not_added=True,
            actual_spatial_Z_is_phase_held_Z_by_radius_Z_zero=True,no_y_correction_jet_exported=True,
            source_intervals_not_selected_field_values=True,**dict.fromkeys(OPEN,False))
        self.cache[key]=result;return result

    def local_integral(self,ends,chart,left,right):
        """A signed C1 function enclosure under fixed Z-independent limits."""
        source=self.query(ends,chart,left,right);op=self.owner(ends);f,c=op.flow,op.c
        l,r=exact(left),exact(right)
        if r<=l:raise ValueError('Positive actual integral interval required')
        cv=lambda q:c.mpf(q.numerator)/q.denominator
        if chart=='O2_axial':width=c.exp(40*cv(l))*c.expm1(40*cv(r-l))
        else:width=cv(r-l)
        density=source['actual_signed_nonlinear_density_C0_Z'];integrals={};masses={}
        for name,rate in RATES.items():
            mass=current.original.incoming.weighted.terminal.local.positive_kernel_mass(c,width,rate)
            integrals[name]=[density[part][name]*mass for part in ('kernels','Z_derivatives')];masses[name]=mass
        return dict(source_identity=self.identity,exact_Z_cell=list(ends),candidate_N=self.N,actual_chart=chart,
            exact_left=list(left),exact_right=list(right),exact_common_P0_axial5=op.P0,
            actual_source_function=source,actual_physical_log_width=width,actual_positive_own_rate_masses=masses,
            actual_signed_local_integral_C0_Z=integrals,
            integral_recipe='int_left^right exp(-rate*(y_right-y))*signed_density(y,Z,N) dy',
            live_inverse_function_extension_under_Z_independent_integral=True,
            directed_rectangle_function_enclosure_not_selected_integral_value=True,
            physical_Jacobian_applied_once=True,incoming_correction_not_supplied_or_reset=True,
            actual_full_prefix_C1_defect_functions_installed=False,actual_terminal_controls_installed=False,
            **dict.fromkeys(OPEN,False))


def run():
    began=time.monotonic()
    with mp.workdps(540):
        owner=WholeZSignedLoopFunctions();queries=[];integrals=[]
        for ends in CELLS:
            for chart,point in (('Rh_reference',(-3,1)),('O2_slope',(1,4)),('O2_axial',(1,2)),
                    ('O2_buffer',(5,1)),('O3_transition',(1,2)),('O3_power',(1,1))):
                queries.append(serialized(owner.query(ends,chart,point)))
                print('Actual whole-Z signed inverse function: '+str(ends)+' '+chart,flush=True)
            for chart,left,right in (('O2_slope',(1,4),(1,2)),('O3_power',(0,1),(1,1))):
                integrals.append(serialized(owner.local_integral(ends,chart,left,right)))
                print('Actual whole-Z signed function local integral: '+str(ends)+' '+chart,flush=True)
        report=dict(**{GATE:True},source_family=owner.identity,candidate_N=owner.N,
            exact_Z_partition=CELLS,current_outer_charts=CHARTS,
            actual_signed_inverse_function_queries=queries,actual_signed_local_integral_function_queries=integrals,
            genuine_Z_only_source_binding=phase.BINDINGS,
            actual_full_prefix_C1_defect_functions_installed=False,actual_terminal_controls_installed=False,
            actual_global_frequency_admitted=False,**dict.fromkeys(OPEN,False),input_hashes=owner.hashes,
            execution_seconds=time.monotonic()-began,
            scope='Live original signed implicit inverse C0/Z/phase functions and local C1 integral '
                  'enclosures on current whole-Z outer source with actual N radius phase. Full prefix '
                  'source/integral oracle, sharp normalized Rc defects/controls/closure, higher correction '
                  'jets/global N/heat/cone and temporal recursion remain open.')
        (HERE/NAME).write_bytes(gzip.compress((json.dumps(encode(report),separators=(',',':'))+'\n').encode(),compresslevel=6,mtime=0))
    print('Current whole-Z actual signed inverse functions and local C1 integrals installed',flush=True);return report


if __name__=='__main__':run()
