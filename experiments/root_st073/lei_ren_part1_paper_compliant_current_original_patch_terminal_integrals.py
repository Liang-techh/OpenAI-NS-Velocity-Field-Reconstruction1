"""Genuine terminal-patch whole-source integrals and affine history extension.

The original actual_patch integrand is restricted to its closed terminal
region. Whole coefficient/phase/error covers and positive dlog(x) masses
enclose its integral. Unknown incoming correction functions remain explicit.
"""
import ast
import copy
from dataclasses import dataclass
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_patch_terminal_points as points
import lei_ren_part1_paper_compliant_current_original_reference_integral_callback as reference
import lei_ren_part1_paper_compliant_current_original_partial_history_functions as partial

whole=reference.whole;leaves=points.leaves;precise=points.precise;base=points.base
HERE,PREFIX,sha,ep=points.HERE,points.PREFIX,points.sha,points.ep
NAME=PREFIX+'current_original_patch_terminal_integrals.json.gz'
RECEIPT=PREFIX+'current_original_patch_terminal_integrals_check.json'
GATE='actual_original_patch_terminal_whole_source_integrals_and_affine_history_extension_installed'


def exact_log_width(c,left,right):
    """Retain exact microscopic rational differences before log1p."""
    if left==right:return c.mpf(0)
    if left.is_Rational and right.is_Rational:
        q=(right-left)/left
        return c.log1p(c.mpf(int(q.p))/int(q.q))
    if not left.is_Rational and not right.is_Rational:
        q=points.log_coordinate(right)-points.log_coordinate(left)
        return c.mpf(int(q.p))/int(q.q)
    return points.patch_offset(c,right)-points.patch_offset(c,left)


def mass(c,left,right,endpoint,rate):
    if not left<right or right>endpoint:raise ValueError('Strict source cell below requested endpoint required')
    width=exact_log_width(c,left,right);distance=exact_log_width(c,right,endpoint)
    r=c.mpf(rate.numerator)/rate.denominator
    result=width if not rate else c.exp(-r*distance)*(-c.expm1(-r*width))/r
    if ep(result)[0]<=0:raise ArithmeticError('Original positive log-radius cell mass lost; refine endpoints')
    return result


def project_whole_method(method):
    tree=ast.parse(Path(whole.__file__).read_text(encoding='utf8'))
    cls=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=='OriginalReferenceWholeCells')
    original=next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name==method)
    fn=copy.deepcopy(original);changed=[]
    targets={"points.reference.reference_coordinate(left)":"terminal_coordinate(left)",
        "points.reference.reference_coordinate(right)":"terminal_coordinate(right)"}
    if method=='roots':targets.update({'a.rational(left)':'patch_offset(c,left)','a.rational(right)':'patch_offset(c,right)'})
    maps={ast.dump(ast.parse(a,mode='eval').body):ast.parse(b,mode='eval').body for a,b in targets.items()}
    class Project(ast.NodeTransformer):
        def visit(self,node):
            key=ast.dump(node)
            if key in maps:
                replacement=copy.deepcopy(maps[key]);changed.append((copy.deepcopy(node),ast.dump(replacement)));return replacement
            return super().visit(node)
    Project().visit(fn)
    if len(changed)!={'roots':4,'cell':2}[method]:raise ValueError('Original whole-cell coordinate sites changed')
    reverse={b:a for a,b in changed};restored=copy.deepcopy(fn)
    class Restore(ast.NodeTransformer):
        def visit(self,node):
            key=ast.dump(node)
            if key in reverse:return copy.deepcopy(reverse[key])
            return super().visit(node)
    Restore().visit(restored)
    if ast.dump(restored)!=ast.dump(original):raise ValueError('Original whole-source enclosure math changed')
    env=dict(vars(whole),terminal_coordinate=points.terminal_coordinate,patch_offset=points.patch_offset)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),
        '<original whole-source method; terminal patch coordinate>','exec'),env)
    return env[method]


class OriginalPatchTerminalWholeCells:
    def __init__(self,provider,*,Z):
        if type(provider) is not points.OriginalPatchTerminalPointLeaves:raise TypeError('Genuine accepted terminal source provider required')
        self.provider=provider;self.family=provider.family;self.inputs=provider.owner.inputs;self.scales=provider.owner.scales
        self.atlas=whole.factors.OriginalSourceFactorAtlas(self.inputs.frame,Z=Z)
        self.ctx=c=self.atlas.ctx;self.Z=self.atlas.Z
        if self.Z==0:raise ValueError('Strict nonzero-Z whole-cell signed route required')
        self.templates=self.inputs.templates;self.compiled={};self.sensitivity={}
        for key,terms in self.templates['rows'].items():
            self.compiled[key]=tuple((powers,s.lambdify(self.templates['inputs'],expr,modules=[{'mpf':c.mpf},'mpmath'])) for powers,expr in terms)
            self.sensitivity[key]=tuple(tuple(s.lambdify(self.templates['inputs'],s.diff(expr,P0),modules=[{'mpf':c.mpf},'mpmath'])
                for P0 in self.templates['pressure_symbols']) for powers,expr in terms)
        self.hashes={**provider.hashes,**self.atlas.hashes};self.cache={}

    roots=project_whole_method('roots')
    original_cell=project_whole_method('cell')

    def phase_boxes(self,left,right,N):
        p=self.provider;c=self.ctx;row=p.rows['actual_patch','density_m_C0']
        phase=p.phase.phase_for_source(row,p.built,coordinate=points.phase_coordinate(left),N=N)
        width=exact_log_width(c,left,right)
        if ep(width)[0]<=0:raise ArithmeticError('Strict original log-radius width required')
        shift=c.mpf((0,ep(c.mpf(N)*width)[1]))
        boxes=[c.mpf(ep(v)) for v in phase['actual_phase_directed_boxes']]
        projection=precise.periodic_add(c,boxes,shift)
        return projection['boxes'],dict(phase,actual_cell_phase_increment=N*width,
            exact_source_phase_increment='N*log(right/left)',original_patch_log_radius_measure='dx/x',
            true_source_phase_at_left_plus_complete_log_radius_width=True)

    def cell(self,left,right,*,N,bits=32):
        self.provider.require_row(self.provider.rows['actual_patch','density_m_C0'])
        if type(N) is not int or N<160 or N.bit_length()>4096:raise ValueError('Explicit original N>=160 required')
        if type(bits) is not int or not 4<=bits<=256:raise ValueError('Inverse bits in[4,256] required')
        left=points.terminal_coordinate(left);right=points.terminal_coordinate(right)
        result=self.original_cell(left,right,N=N,bits=bits);record=result['record']
        record['exact_patch_x_cell']=record.pop('exact_reference_cell',record.get('exact_patch_x_cell'))
        record.update(chart='actual_patch',source_graph_sha256=self.provider.source_graph_sha256,
            exact_current_R_offset_bounds=[str(points.log_coordinate(left)-6),str(points.log_coordinate(right)-6)],
            original_terminal_recipe=self.provider.recipe,
            original_current_R_normalization_and_inertial_compiler_retained=True)
        return result


def density_binding(provider):
    tree=ast.parse(Path(reference.__file__).read_text(encoding='utf8'))
    original=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='exact_density_binding')
    fn=copy.deepcopy(original);count=0
    class Project(ast.NodeTransformer):
        def visit_Constant(self,node):
            nonlocal count
            if node.value=='Rh_reference':count+=1;return ast.Constant(value='actual_patch')
            return node
    Project().visit(fn)
    if count!=2:raise ValueError('Original density source sites changed')
    env=dict(vars(reference));exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),
        '<same original coefficient identity; actual_patch source graph>','exec'),env)
    return env['exact_density_binding'](provider)


@dataclass(frozen=True)
class OriginalPatchTerminalIntegralFrame:
    left: object
    right: object
    Z: object
    N: int
    owner: object
    forcing: dict
    coefficients: dict
    decay: dict
    cells: tuple
    record: dict


class OriginalPatchTerminalIntegrals:
    mode='original_patch_terminal_integral_partial'
    def __init__(self,provider=None):
        self.provider=points.OriginalPatchTerminalPointLeaves() if provider is None else provider
        if type(self.provider) is not points.OriginalPatchTerminalPointLeaves:raise TypeError('Genuine terminal source required')
        self.family=self.source_family=self.provider.family;self.source_graph_sha256=self.provider.source_graph_sha256
        self.hashes=dict(self.provider.hashes);self.owners={};self.cache={};self.issued={}
        for module in (points,whole):
            record=json.loads((HERE/module.RECEIPT).read_bytes())
            if not record.get('all_passed') or not record.get(module.GATE) or record['source_family']!=self.family:raise ValueError('Accepted same-source cell prerequisite required')
            self.bind_hashes({**record['input_hashes'],module.RECEIPT:sha(module.RECEIPT)})
        self.density_identity=density_binding(self.provider);self.contracts=self.kernel_bindings()
        self.bind_hashes({Path(reference.__file__).name:sha(Path(reference.__file__).name),Path(__file__).name:sha(Path(__file__).name)})

    def bind_hashes(self,closure):
        for name,digest in closure.items():
            if sha(name)!=digest:raise ValueError('Terminal integral source changed: '+name)
            if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Terminal integral closures disagree: '+name)
            self.hashes[name]=digest

    def kernel_bindings(self):
        p=self.provider;g=p.built['graph'];saved=json.loads((HERE/leaves.transport.NAME).read_bytes())
        cells=[row for row in saved['exact_original_cells'] if row['chart']=='actual_patch']
        if len(cells)!=1:raise ValueError('One original actual_patch contribution required')
        self.original_cell=cells[0];result=[]
        for key,pair in self.original_cell['contributions'].items():
            for order,label in (('C0','value'),('Z','Z')):
                row=g.nodes[pair[label]];source=p.rows['actual_patch','density_'+key+'_'+order];p.require_row(source)
                if row['operation']!='definite_integral' or row['lower']!=self.original_cell['lower'] or row['upper']!=self.original_cell['upper']:raise ValueError('Original defining endpoints required')
                lower,_=reference.symbolic_node(p,row['lower']);upper,_=reference.symbolic_node(p,row['upper'])
                if lower!=1 or upper!=s.E:raise ValueError('Original full patch window[1,e] required')
                if g.nodes[source['coordinate']]!={'operation':'bound_variable','name':row['variable']}:raise ValueError('Same bound source coordinate required')
                expr,F=reference.symbolic_node(p,row['integrand'],source);old=s.Symbol(row['variable'],real=True);x=s.Symbol(row['variable'],positive=True)
                expr=expr.subs(old,x);rate=leaves.transport.RATES[key];r=s.Rational(rate.numerator,rate.denominator)
                B=s.Symbol('requested_patch_right',positive=True)
                expected=s.exp(-r*(1-s.log(x)))*F/x
                if s.simplify(expr-expected)!=0 or s.simplify(expr*s.exp(r*(1-s.log(B)))-s.exp(-r*(s.log(B)-s.log(x)))*F/x)!=0:raise ValueError('Original patch density/kernel/dx/x differs')
                if row['measure']!='native coordinate; original dy/dcoordinate applied exactly once':raise ValueError('Original source measure required')
                result.append(dict(original_integral_node=pair[label],source_node=source['source_node'],role=source['function_role'],
                    key=key,ordinary_order=order,rate=str(rate),original_full_endpoints=['1','e'],
                    original_full_integrand='exp(-rate*(1-log(x)))*issued_density/x',
                    localized_integrand='exp(-rate*(log(b)-log(x)))*issued_density/x',
                    original_dlogx_Jacobian_applied_exactly_once=True,source_kernel_localization_identity=True,
                    tail_subinterval_not_original_full_integral=True))
        return result

    def request(self,Z,N):
        z=point_Z=base.point.pressure.exact_Z(Z)
        if not z:raise ValueError('Nonzero fixed Z required by signed whole-cell route')
        if type(N) is not int or N<160 or N.bit_length()>4096:raise ValueError('Original candidate integer N>=160 required')
        self.provider.require_row(self.provider.rows['actual_patch','density_m_C0']);return point_Z

    def partition(self,left,right,count):
        if type(count) is not int or not 1<=count<=256:raise ValueError('Exact partition level in[1,256] required')
        if left>right:raise ValueError('Ordered terminal interval required')
        # Partition choices are exact rational geometry, never field samples.
        grid=[points.EDGE+(s.Rational(8,3)-points.EDGE)*s.Rational(i,count) for i in range(count+1)]
        inside=[x for x in grid if left<x<right]
        return [left,*inside,right] if left<right else [left]

    def subinterval(self,left,right,*,Z,N,count=4,bits=32):
        a=points.terminal_coordinate(left);b=points.terminal_coordinate(right);z=self.request(Z,N)
        if type(bits) is not int or not 4<=bits<=256:raise ValueError('Inverse bits in[4,256] required')
        knots=self.partition(a,b,count);key=a,b,z,N,count,bits
        if key in self.cache:return self.cache[key]
        if z not in self.owners:
            self.owners[z]=OriginalPatchTerminalWholeCells(self.provider,Z=str(z));self.bind_hashes(self.owners[z].hashes)
        owner=self.owners[z];c=owner.ctx;atlas=owner.atlas
        total={o:{k:{j:atlas.scalar(0) for j in ('C0','Z')} for k in leaves.transport.RATES} for o in whole.points.exact.ORDERS}
        records=[]
        with mp.workdps(c.dps+40):
            decay={k:c.exp(-(c.mpf(rate.numerator)/rate.denominator)*exact_log_width(c,a,b)) for k,rate in leaves.transport.RATES.items()}
            for l,r in zip(knots,knots[1:]):
                query=owner.cell(l,r,N=N,bits=bits);masses={}
                for k,rate in leaves.transport.RATES.items():
                    masses[k]=mass(c,l,r,b,rate)
                    for o in total:
                        for j in ('C0','Z'):
                            value=query['coefficients'][o][k][j]
                            if value.ctx is not c or value.ledger is not atlas.ledger or value.scale.bases is not atlas.bases:raise ValueError('One canonical source atlas required')
                            total[o][k][j]=atlas.add(total[o][k][j],value*masses[k])
                records.append(dict(source=query['record'],exact_positive_requested_endpoint_masses=masses))
            forcing={};full={}
            for k in leaves.transport.RATES:
                full[k]={}
                for j in ('C0','Z'):
                    value=atlas.add(total[-1][k][j]*(c.mpf(1)/N),total[-2][k][j]*(c.mpf(1)/N**2));full[k][j]=value
                    forcing[k,j]=base.conditioned.bounded_value(value)
        record=dict(source_family=self.family,source_graph_sha256=self.source_graph_sha256,
            exact_patch_x_interval=[str(a),str(b)],original_Z_exact=str(z),explicit_candidate_N=N,count=count,bits=bits,
            exact_partition_knots=[str(x) for x in knots],source_kernel_bindings=self.contracts,
            local_coefficient_program_equals_issued_actual_patch_density=self.density_identity,
            canonical_atlas=atlas.record(),actual_whole_source_phase_cells=records,
            actual_ordinary_C0_Z_forcing={k:{j:forcing[k,j] for j in ('C0','Z')} for k in leaves.transport.RATES},
            directed_homogeneous_multipliers=decay,exact_N_dependent_coefficient_integrals={str(o):{k:{j:v.record() for j,v in pair.items()} for k,pair in row.items()} for o,row in total.items()},
            unknown_incoming_correction_at_exact_left_endpoint_retained=True,separate_original_P0_P0_Z_unchanged=True,
            true_whole_cell_sources_not_point_quadrature=True,empty_interval_has_zero_forcing_and_unit_decay=a==b,
            full_original_patch_integral_installed=False,whole_Z_contract_installed=False)
        frame=OriginalPatchTerminalIntegralFrame(a,b,z,N,owner,forcing,total,decay,tuple(records),record)
        self.cache[key]=frame;self.issued[id(frame)]=frame;return frame

    def tail(self,*,Z,N,count=4,bits=32):return self.subinterval(points.EDGE,s.E,Z=Z,N=N,count=count,bits=bits)
    def integrate(self,*args,**kwargs):raise NotImplementedError('Original full[1,e] patch integral still needs active source')

    def extend(self,frame,history,Y):
        if type(frame) is not OriginalPatchTerminalIntegralFrame or self.issued.get(id(frame)) is not frame or frame.right!=s.E:raise ValueError('Issued tail ending at true Rh required')
        if type(history) is not partial.OriginalPartialHistoryFunctions or history.current.provider is not self.provider.provider:raise ValueError('One genuine live original provider required')
        if history.family!=self.family or history.source_graph_sha256!=self.source_graph_sha256:raise ValueError('Same history source/graph required')
        downstream=history.prefix(Y,Z=str(frame.Z),N=frame.N);c=history.c;forcing={};decay={}
        with mp.workdps(c.dps+40):
            for k in leaves.transport.RATES:
                decay[k]=downstream.decay[k]*history.current.ordinary(frame.decay[k])
                for j in ('C0','Z'):
                    forcing[k,j]=downstream.forcing[k,j]+downstream.decay[k]*history.current.ordinary(frame.forcing[k,j])
        self.bind_hashes(history.hashes)
        return dict(forcing=forcing,decay=decay,record=dict(source_family=self.family,source_graph_sha256=self.source_graph_sha256,
            original_Z_exact=str(frame.Z),explicit_candidate_N=frame.N,exact_patch_left=str(frame.left),reference_O2_right=str(downstream.right),
            ordinary_C0_Z_forcing={k:{j:forcing[k,j] for j in ('C0','Z')} for k in leaves.transport.RATES},
            directed_incoming_multipliers=decay,unknown_incoming_at_patch_left_retained=True,
            actual_original_ten_Rh_incoming_root_roles={str(node):dict(key=k,ordinary_order=j) for node,(k,j) in history.current.boundaries.items()},
            exact_operator_composition='H(Y)=D_ref_to_Y*I_tail+I_ref_to_Y+D_ref_to_Y*D_tail*unknown_H(left)',
            physical_endpoint_identity='x=e <=> original reference offset=-5',
            original_P0_P0_Z_not_reset=True,actual_source_owned_upstream_history_installed=False))


@precise.phase.native.inlet.source_precision
def run(*,return_live=False):
    began=time.monotonic();history=partial.OriginalPartialHistoryFunctions()
    owner=OriginalPatchTerminalIntegrals(points.OriginalPatchTerminalPointLeaves(history.current.provider));frames=[]
    for count,N in ((4,160),(16,160),(4,257)):
        frames.append(owner.tail(Z='.37',N=N,count=count));print('Actual terminal whole-source tail:',count,N,flush=True)
    for a,b in ((s.Rational(2),s.Rational(5,2)),(s.Rational(2),s.Rational(2)),(s.Rational(2),s.Rational(2)+s.Rational(1,10**1000))):
        frames.append(owner.subinterval(a,b,Z='.37',N=160));print('Actual terminal partial width:',str(b-a)[:70],flush=True)
    extended=[owner.extend(frames[0],history,'-.731'),owner.extend(frames[2],history,'.731')]
    result=dict(**{GATE:True},source_family=owner.family,source_graph_sha256=owner.source_graph_sha256,mode=owner.mode,
        actual_terminal_integral_frames=[f.record for f in frames],actual_terminal_to_reference_O2_affine_extensions=[e['record'] for e in extended],
        original_graph_unchanged=leaves.digest_rows(owner.provider.built['graph'].nodes)==owner.provider.graph_digest,
        full_original_patch_integral_installed=False,whole_Z_contract_installed=False,actual_source_owned_upstream_history_installed=False,
        actual_five_controls_installed=False,current_whole_N_selected=False,**dict.fromkeys(precise.phase.packets.OPEN,False),
        input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
        scope='Actual terminal closed whole-source C0/Z partial/tail integrals with dx/x and original signed graph, exact N-dependent coefficients, pressure/phase errors; affine extension to accepted reference/O2 preserves unknown incoming. Active patch/full upstream/control/global-N/recursion/corrected UVW remain open.')
    (HERE/NAME).write_bytes(gzip.compress(json.dumps(precise.encode(result),indent=2).encode()+b'\n',compresslevel=9,mtime=0))
    return (result,owner,frames,history,extended) if return_live else result


if __name__=='__main__':run()
