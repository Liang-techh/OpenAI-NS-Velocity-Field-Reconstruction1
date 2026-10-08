"""Whole axial rectangles for original O2 five-history value/Z transport.

Every source cell is evaluated over the entire declared Z interval. The
accepted Z-independent true radial phase cover and positive own-rate masses
are retained. Original histories and incoming defects remain separate.
"""
import ast
import gzip
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_O2_density_Z_integrals as prior

current=prior.current;base=prior.base;ep=prior.ep
HERE,PREFIX,sha=prior.HERE,prior.PREFIX,prior.sha
NAME=PREFIX+'current_original_O2_continuous_Z_transport.json'
RECEIPT=PREFIX+'current_original_O2_continuous_Z_transport_check.json'
GATE='original_O2_continuous_strict_sign_Z_five_value_and_derivative_transport_enclosed'
KEYS=tuple(current.five.RATES)
FLAGS=('numerical_original_source_point_or_integral_oracle_installed',
       'actual_five_controls_installed','current_whole_N_selected',
       *base.point.source.inertial.profiles.loop.OPEN)


def original_inlet_contract():
    """Bind the original slope history assignment, then derive its inlet jet."""
    filename=PREFIX+'pre_pulse_mixed_C4.py'
    tree=ast.parse((HERE/filename).read_text(encoding='utf8'))
    method=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name=='slope')
    assignment=next(n for n in ast.walk(method) if isinstance(n,ast.Assign)
        and any(isinstance(t,ast.Name) and t.id=='hist' for t in n.targets))
    expected="dict(m=V,h=h,k=h*V,e=square(V)*self.invP2-square(qi)*(c.mpf(5)/12+mass[2]/2)*c.exp(-y),p=square(qi)*(c.mpf('2.5')+mass[1]/2))"
    if ast.dump(assignment.value)!=ast.dump(ast.parse(expected,mode='eval').body):
        raise ValueError('Original O2 cumulative history source changed')
    def require_statement(scope,statement):
        wanted=ast.dump(ast.parse(statement).body[0])
        if sum(ast.dump(n)==wanted for n in ast.walk(scope))!=1:
            raise ValueError('Original inlet source dataflow changed: '+statement)
    coordinates=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name=='coordinates')
    for statement in ('z=IntervalTaylor.variable(c,Z,5)','qi=(1+square(z)).reciprocal()','return z,qi'):
        require_statement(coordinates,statement)
    for statement in ('z,qi=self.coordinates(Z)','J,mass=slope_masses(c,y,self.cells)',
            "factor=c.exp(y/10-c.mpf('.6')*J)",'u=qi*factor','V=z*4',
            "h=qi*(c.mpf('.625')+mass[0])*c.exp(-c.mpf('1.5')*y)"):
        require_statement(method,statement)
    bindings={PREFIX+'pre_pulse_mixed_C4.py':('self.invP2=self.initial.invP2',),
        PREFIX+'outer_initial.py':('self.invP2=self.repair.invP2',),
        PREFIX+'five_moment_repair.py':("self.invP2=get(self.admit,'inverse_Pstar_squared')",),
        PREFIX+'five_defect_admission.py':('invP2=c.exp(-2*logP)',),
        'lei_ren_part1_paper_logarithmic_outer_parameters.py':('self.logPstar=c.exp(self.md)+11',),
        PREFIX+'current_original_O2_source_parameter_frame.py':('Md=s.Integer(40)','logP=s.exp(Md)+11')}
    for source,statements in bindings.items():
        source_tree=ast.parse((HERE/source).read_text(encoding='utf8'))
        for statement in statements:require_statement(source_tree,statement)
    defining=ast.parse((HERE/'lei_ren_part1_paper_interval_outer_slope_field.py').read_text(encoding='utf8'))
    transition=next(n for n in ast.walk(defining) if isinstance(n,ast.FunctionDef) and n.name=='transition_integrals')
    for statement in ('y=Fraction(y)','j=c.mpf(0)','masses=[c.mpf(0) for _ in range(3)]',
            'dy=fraction_box(c,y/cells)','next_j=j+dy*sj',
            "masses[k]+=dy*c.exp(rate*scell-c.mpf('.6')*power*jcell)",'j=next_j',
            "if y==1:j=c.mpf('.5')",'return j,masses'):
        require_statement(transition,statement)
    # At y0, dy=0: the bound dataflow preserves zero J and every mass.
    # There is no invented explicit y0 return branch in the original source.
    Z,P=sy.symbols('Z Pstar',real=True,nonzero=True);C=1/(1+Z*Z)
    values=dict(m=4*Z/P,h=sy.Rational(5,8)*C,k=sy.Rational(5,2)*Z*C/P,
        e=16*Z*Z/P**2-sy.Rational(5,12)*C*C,p=sy.Rational(5,2)*C*C)
    derivatives=dict(m=4/P,h=-sy.Rational(5,4)*Z*C*C,
        k=sy.Rational(5,2)*(1-Z*Z)*C*C/P,
        e=32*Z/P**2+sy.Rational(5,3)*Z*C**3,p=-10*Z*C**3)
    for key in KEYS:assert sy.cancel(sy.diff(values[key],Z)-derivatives[key])==0
    return dict(passed=True,original_slope_history_AST_bound=True,
        source_filename=filename,source_sha256=sha(filename),
        original_defining_J_and_three_masses_exact_zero_at_y0=True,
        normalized_common_velocity_unit='S=Pstar',
        ordinary_inlet_Z_symbolic_identities=True,
        original_inlet_functions={k:str(v) for k,v in values.items()},
        original_inlet_Z_functions={k:str(v) for k,v in derivatives.items()},
        pressure_datum_separate_from_original_cumulative_p=True)


def original_inlet_source(kernel,z):
    """Native source-defined nonmidplane inlet; positive Pstar^-1 is formal."""
    c=kernel.c;bases=kernel.q.scale.bases;ledger=kernel.q.ledger
    scalar=kernel.scalar;C=1/(1+z*z)
    inverse=base.prior.ScaledEnclosure(base.prior.FormalScale(bases,(-1,0,0,0,0)),c.mpf(1),ledger)
    inv2=prior.slow.square(inverse)
    old=dict(m=inverse*(4*z),h=scalar(c.mpf(5)/8*C),k=inverse*(c.mpf(5)/2*z*C),
        e=inv2*(16*z*z)-scalar(c.mpf(5)/12*C*C),p=scalar(c.mpf(5)/2*C*C))
    derivative=dict(m=inverse*4,h=scalar(-c.mpf(5)/4*z*C*C),
        k=inverse*(c.mpf(5)/2*(1-z*z)*C*C),
        e=inv2*(32*z)+scalar(c.mpf(5)/3*z*C**3),p=scalar(-10*z*C**3))
    return dict(original_inlet={k:current.bounded(v) for k,v in old.items()},
        original_inlet_Z={k:current.bounded(v) for k,v in derivative.items()},
        native_original_inlet_sources={k:v.record() for k,v in old.items()},
        native_original_inlet_Z_sources={k:v.record() for k,v in derivative.items()},
        microscopic_positive_inverse_Pstar_retained=True,
        whole_Z_rational_source_functions_not_axial_samples=True)


def source_history_densities(part):
    """The original unmodulated density functions and ordinary derivatives."""
    roots=part['roots'];E,V=(roots[k][(0,0)] for k in ('E','V'))
    EZ,VZ=(roots[k][(0,1)] for k in ('E','V'))
    square=prior.slow.square
    values=dict(m=V,h=E,k=E*V,e=square(V)-square(E)*.5,p=square(E)*.5)
    derivatives=dict(m=VZ,h=EZ,k=EZ*V+E*VZ,e=V*VZ*2-E*EZ,p=E*EZ)
    return values,derivatives


def require_identity(record,source_family,original_P0_datum_sha256):
    if source_family!=record['source_family'] or original_P0_datum_sha256!=source_family['datum_enclosure_sha256']:
        raise ValueError('Same original source family and separate P0 datum required')
    if record['original_P0_datum_sha256']!=original_P0_datum_sha256:
        raise ValueError('Recorded original pressure datum differs')
    if record['exact_y_window']!=['0','1'] or record['own_rates']!=current.five.RATES:
        raise ValueError('Original unit y window and five own rates required')
    if record['normalized_own_units']!=current.five.UNITS:
        raise ValueError('Original common-unit five histories required')


def explicit_ranges(c,values,label):
    if set(values)!=set(KEYS):raise ValueError('All five explicit '+label+' function covers required')
    result={k:(current.interval(c,v) if isinstance(v,dict) else c.mpf(v)) for k,v in values.items()}
    if any(not mp.isfinite(x) for v in result.values() for x in ep(v)):
        raise ValueError('Finite '+label+' whole-Z function covers required')
    return result


def apply_history_transport(c,record,*,original_inlet,original_inlet_Z,
        incoming,incoming_Z,source_family,original_P0_datum_sha256):
    """Affine ordinary-C1 transport of explicitly supplied actual histories.

    Inputs must enclose their source functions on record's full Z interval.
    This arithmetic interface cannot certify an arbitrary caller's functions.
    In particular, no missing original or incoming history defaults to zero.
    """
    require_identity(record,source_family,original_P0_datum_sha256)
    old=explicit_ranges(c,original_inlet,'original inlet')
    oldZ=explicit_ranges(c,original_inlet_Z,'original inlet ordinary Z derivative')
    inc=explicit_ranges(c,incoming,'incoming defect')
    incZ=explicit_ranges(c,incoming_Z,'incoming defect ordinary Z derivative')
    original={};originalZ={};defect={};defectZ={};own={};ownZ={}
    for key,rate in current.five.RATES.items():
        decay=c.exp(-c.mpf(rate))
        original[key]=decay*old[key]+current.interval(c,record['five_unmodulated_original_integral_contributions'][key])
        originalZ[key]=decay*oldZ[key]+current.interval(c,record['five_unmodulated_original_Z_integral_contributions'][key])
        defect[key]=decay*inc[key]+current.interval(c,record['five_original_C0_integral_contributions'][key])
        defectZ[key]=decay*incZ[key]+current.interval(c,record['five_genuine_ordinary_Z_integral_contributions'][key])
        own[key]=original[key]+defect[key];ownZ[key]=originalZ[key]+defectZ[key]
    return dict(source_family=source_family,exact_Z_range=record['exact_Z_range'],
        exact_y_window=record['exact_y_window'],original_P0_datum_sha256=original_P0_datum_sha256,
        original_histories_at_y1=original,original_history_Z_at_y1=originalZ,
        transported_incoming_and_modulation_defects=defect,transported_defect_Z=defectZ,
        own_five_histories_at_y1=own,own_five_ordinary_Z_at_y1=ownZ,
        input_function_covers_are_explicit_caller_obligations=True,
        actual_incoming_histories_not_set_to_zero=True,P0_not_reset_or_added_to_cumulative_pressure=True,
        functional_terminal_identity_solved=False,**dict.fromkeys(FLAGS,False))


def apply_actual_incoming(c,record,*,incoming,incoming_Z,source_family,original_P0_datum_sha256):
    """Use the source-bound original inlet and preserve explicit incoming C1 data."""
    inlet=record['source_defined_original_inlet']
    return apply_history_transport(c,record,original_inlet=inlet['original_inlet'],
        original_inlet_Z=inlet['original_inlet_Z'],incoming=incoming,incoming_Z=incoming_Z,
        source_family=source_family,original_P0_datum_sha256=original_P0_datum_sha256)


class OriginalO2ContinuousZTransport:
    def __init__(self):
        accepted=json.loads((HERE/prior.RECEIPT).read_bytes())
        if not accepted.get('all_passed') or not accepted.get(prior.GATE):
            raise ValueError('Accepted genuine original density-Z integrals required')
        self.parent=prior.OriginalO2DensityZIntegrals();self.c=self.parent.c;self.family=self.parent.family
        self.hashes=dict(self.parent.hashes)
        for name,digest in {**accepted['input_hashes'],prior.RECEIPT:sha(prior.RECEIPT)}.items():
            if sha(name)!=digest:raise ValueError('Continuous-Z dependency changed: '+name)
            if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Source closure differs')
            self.hashes[name]=digest
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.inlet_contract=original_inlet_contract()

    def integrate(self,count,*,Z_lower,Z_upper,N=7,bits=24):
        c=self.c;zl,zh=(base.point.pressure.exact_Z(v) for v in (Z_lower,Z_upper))
        if not (-1<=zl<zh<=1) or zl<=0<=zh:
            raise ValueError('Nondegenerate strict-sign original Z interval in [-1,1] required')
        level,phase=self.parent.parent.parent.levels[count]
        if N!=phase['explicit_candidate_N']:raise ValueError('Same actual common candidate phase cover required')
        begin=time.monotonic();totals=[{k:c.mpf(0) for k in KEYS} for unused in range(4)];rows=[];inlet=None
        with mp.workdps(c.dps+40):
            z=c.mpf((ep(c.mpf(int(zl.p))/int(zl.q))[0],ep(c.mpf(int(zh.p))/int(zh.q))[1]))
            for i,phase_cell in enumerate(phase['whole_source_cells']):
                query=self.parent.parent.query(count,i,Z_lower=str(zl),Z_upper=str(zh))
                if not query['pieces']:raise ArithmeticError('Whole-Z source needs signed refinement at '+str(i)+'/'+str(count))
                covers=[{k:[] for k in KEYS} for unused in range(4)];proofs=[]
                for source_index,original_part in enumerate(query['pieces']):
                    part=prior.correlated_p2_Z_source(original_part,query['record']);kernel=part['kernel']
                    if inlet is None:inlet=original_inlet_source(kernel,z)
                    old,oldZ=source_history_densities(part)
                    for k in KEYS:
                        covers[2][k].append(current.bounded(old[k]));covers[3][k].append(current.bounded(oldZ[k]))
                    for phase_cover in phase_cell['true_common_N_phase_boxes']:
                        phi=current.interval(c,phase_cover);selected=kernel.evaluate(phi,bits=bits)['selected_inverse']
                        primitive=kernel.primitives(selected['coordinate_interval'],selected['chart'])
                        derivative=prior.derivative_inverse(kernel,phi,bits)
                        jets,proof=prior.slow_values(part,derivative['coordinate_interval'],derivative['chart'])
                        graph=current.CellDensityGraph(self.parent.parent.owner.scales.graph,part,primitive,jets,N)
                        got=graph.outputs()
                        for key in KEYS:
                            covers[0][key].append(current.bounded(got['densities'][key]))
                            covers[1][key].append(current.bounded(got['density_Z'][key]))
                        proofs.append(dict(source_piece_index=source_index,actual_true_phase_cover=phase_cover,
                            correlated_original_p2_Z_carrier=part['p2_Z_common_carrier_binding'],
                            original_inverse=selected,original_derivative_inverse=derivative,
                            original_slow_Z_contract=proof,
                            native_original_primitive_Z_enclosures={k:v.record() for k,v in jets.items()},
                            density_Z_is_genuine_ordinary_source_derivative=True))
                hull=lambda vs:c.mpf((min(ep(v)[0] for v in vs),max(ep(v)[1] for v in vs)))
                hulls=[{k:hull(v) for k,v in cover.items()} for cover in covers]
                unused,masses,decay=current.five.masses(c,c.mpf(1)/count,c.mpf(i)/count,c.mpf(i+1)/count)
                adds=[{k:cover[k]*masses[k] for k in KEYS} for cover in hulls]
                for total,add in zip(totals,adds,strict=True):
                    for key in KEYS:total[key]+=add[key]
                rows.append(dict(source=query['record'],actual_source_phase_held_Z_records=proofs,
                    five_signed_C0_density_hulls=hulls[0],five_genuine_signed_ordinary_Z_density_hulls=hulls[1],
                    five_original_unmodulated_density_hulls=hulls[2],five_original_unmodulated_density_Z_hulls=hulls[3],
                    positive_own_rate_final_endpoint_masses=masses,five_C0_contributions=adds[0],
                    five_ordinary_Z_contributions=adds[1],five_original_unmodulated_contributions=adds[2],
                    five_original_unmodulated_Z_contributions=adds[3],
                    entire_source_y_Z_rectangle_enclosed=True,phase_and_masses_Z_independent=True))
                if (i+1)%max(8,count//8)==0:
                    print('Continuous original Z transport:',count,str(zl),str(zh),i+1,flush=True)
            original={k:c.exp(-c.mpf(current.five.RATES[k]))*inlet['original_inlet'][k]+totals[2][k] for k in KEYS}
            originalZ={k:c.exp(-c.mpf(current.five.RATES[k]))*inlet['original_inlet_Z'][k]+totals[3][k] for k in KEYS}
        return dict(source_family=self.family,ordered_source_cells=count,exact_Z_range=[str(zl),str(zh)],
            exact_y_window=['0','1'],explicit_candidate_N=N,inverse_bits=bits,own_rates=current.five.RATES,
            normalized_own_units=current.five.UNITS,five_original_C0_integral_contributions=totals[0],
            five_genuine_ordinary_Z_integral_contributions=totals[1],
            five_unmodulated_original_integral_contributions=totals[2],
            five_unmodulated_original_Z_integral_contributions=totals[3],
            source_defined_original_inlet=inlet,original_nonzero_Z_inlet_contract=self.inlet_contract,
            source_defined_original_histories_at_y1=original,source_defined_original_history_Z_at_y1=originalZ,
            whole_original_density_and_density_Z_source_records=rows,
            original_P0_datum_sha256=self.family['datum_enclosure_sha256'],
            entire_continuous_Z_interval_not_samples=True,actual_incoming_histories_not_set_to_zero=True,
            complete_original_y_window_all_source_and_phase_pieces_integrated=True,
            Z_interval_terminal_identities_or_global_controls_installed=False,current_whole_N_selected=False,
            execution_seconds=time.monotonic()-begin)


def run():
    begin=time.monotonic();owner=OriginalO2ContinuousZTransport();archives=[]
    for tag,zl,zh in (('positive','.36','.38'),('negative','-.38','-.36')):
        for count in (256,2048):
            report=owner.integrate(count,Z_lower=zl,Z_upper=zh)
            name=PREFIX+'current_original_O2_continuous_Z_transport_'+tag+'_'+str(count)+'.json.gz'
            raw=json.dumps(base.encoded(report),indent=2).encode('utf8')+b'\n'
            compressed=gzip.compress(raw,compresslevel=9,mtime=0)
            if len(compressed)>=100*1024*1024:raise ArithmeticError('Evidence exceeds GitHub single-file limit; split losslessly')
            (HERE/name).write_bytes(compressed)
            archives.append(dict(filename=name,compressed_bytes=len(compressed),uncompressed_bytes=len(raw),
                lossless_original_json_sha256=hashlib.sha256(raw).hexdigest(),exact_Z_range=report['exact_Z_range'],
                ordered_source_cells=count,five_original_C0_integral_contributions=report['five_original_C0_integral_contributions'],
                five_genuine_ordinary_Z_integral_contributions=report['five_genuine_ordinary_Z_integral_contributions'],
                five_unmodulated_original_integral_contributions=report['five_unmodulated_original_integral_contributions'],
                five_unmodulated_original_Z_integral_contributions=report['five_unmodulated_original_Z_integral_contributions']))
            print('Archived continuous Z rectangle',tag,count,len(compressed),flush=True)
    result=dict(**{GATE:True},source_family=owner.family,continuous_Z_integral_archives=archives,
        exact_original_source_derivative_identities=owner.parent.identities,
        entire_source_y_Z_rectangles_and_original_phase_unions=True,
        actual_original_nonzero_Z_inlet_functions_installed=True,
        explicit_actual_incoming_value_and_derivative_interface_installed=True,
        functional_terminal_identity_solved=False,**dict.fromkeys(FLAGS,False),
        input_hashes=owner.hashes,execution_seconds=time.monotonic()-begin,
        scope='Original O2 continuous strict-sign Z rectangles, candidate N7, five original and modulation densities and ordinary-Z integrals. Explicit incoming C1 affine transport. No crossing-Z atlas, functional terminal repair, global N, matching/stress/recursion or full NS.')
    (HERE/NAME).write_text(json.dumps(base.encoded(result),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()
