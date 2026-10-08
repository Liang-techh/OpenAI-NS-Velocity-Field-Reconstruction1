"""Original exact-midplane reference whole-window C0/Z contributions.

The original source-wide pressure odd-derivative identity removes only the
P0_Z error budget at exact Z0. Actual p2_Z and all even pressure errors remain.
All original whole-cell/density/rate programs are reused with reversible AST
projections of two explicit guards. This does not relax the signed service.
"""
import ast
import copy
import gzip
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_reference_whole_cell_integrals as whole

points=whole.points;point=whole.point;base=whole.base;ep=whole.ep
HERE,PREFIX,sha=whole.HERE,whole.PREFIX,whole.sha
NAME=PREFIX+'current_original_reference_midplane_integrals.json.gz'
RECEIPT=PREFIX+'current_original_reference_midplane_integrals_check.json'
GATE='original_exact_midplane_whole_reference_C0_Z_contributions_with_pressure_parity_enclosed'


def project_original_midplane_methods():
    """Certify that all original source math outside named guards is unchanged."""
    tree=ast.parse(Path(whole.__file__).read_bytes())
    cls=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=='OriginalReferenceWholeCells')
    originals={name:next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name==name)
               for name in ('__init__','roots')}
    functions={};proof={}
    for name,original in originals.items():
        fn=copy.deepcopy(original)
        if name=='__init__':
            positions=[i for i,n in enumerate(fn.body) if isinstance(n,ast.If)
                and ast.dump(n.test)==ast.dump(ast.parse('self.Z==0',mode='eval').body)]
            if len(positions)!=1:raise ValueError('Original strict midplane guard changed')
            position=positions[0];guard=fn.body.pop(position)
            restored=copy.deepcopy(fn);restored.body.insert(position,guard)
            changes=['separate exact-Z0 adapter removes original signed-only Z0 rejection']
        else:
            loops=[n for n in ast.walk(fn) if isinstance(n,ast.For)
                and ast.unparse(n.target)=='(order, (sensitivity, bound))']
            if len(loops)!=1:raise ValueError('Original pressure jet error loop changed')
            loop=loops[0]
            expected="enumerate(zip(self.sensitivity[key][index], (5, 10, 44), strict=True))"
            if ast.dump(loop.iter)!=ast.dump(ast.parse(expected,mode='eval').body):
                raise ValueError('Original pressure late factors changed')
            rule=ast.parse("if order==1:\n    if not self.parity['passed']:\n        raise ValueError('Original full-source pressure parity proof required')\n    continue").body[0]
            loop.body.insert(0,rule)
            restored=copy.deepcopy(fn)
            restored_loop=next(n for n in ast.walk(restored) if isinstance(n,ast.For)
                and ast.unparse(n.target)=='(order, (sensitivity, bound))')
            restored_loop.body.pop(0)
            changes=['exact Z0 P0_Z late budget is zero by full source odd derivative']
        if ast.dump(restored)!=ast.dump(original):
            raise ValueError('Original source math changed outside scoped midplane guard')
        env=dict(vars(whole))
        exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),
                     '<original whole-reference method; exact midplane guard only>','exec'),env)
        functions[name]=env[name]
        proof[name]=dict(original_source=Path(whole.__file__).name,source_sha256=sha(Path(whole.__file__).name),
            original_method_ast_sha256=hashlib.sha256(ast.dump(original).encode()).hexdigest(),
            changes=changes,reversed_AST_exactly_equals_original=True)
    return functions,proof


PROJECTED,PROJECTION_PROOF=project_original_midplane_methods()


def original_pressure_midplane_parity(frame):
    """Actual fourteen stages, fixed limits, and finite-end/late remainder."""
    pressure=point.pressure
    receipt=json.loads((HERE/pressure.RECEIPT).read_bytes())
    source=json.loads((HERE/pressure.NAME).read_bytes())
    if not receipt.get('all_passed') or not receipt.get(pressure.GATE) or receipt['source_family']!=frame.family:
        raise ValueError('Accepted original normalized pressure source required')
    theorem=receipt['original_late_source_error_theorem']
    if not theorem['passed'] or theorem!=source['original_late_source_error_theorem']:
        raise ValueError('Same original all-late pressure theorem required')
    operator=frame.owner.pressure;z=operator.partition['z'];stages={}
    for name,density in operator.original_densities.items():
        if s.simplify(density.subs(z,-z)-density)!=0:
            raise ValueError('Original pressure density lost even Z parity: '+name)
        if s.simplify(s.diff(density,z).subs(z,0))!=0:
            raise ValueError('Original midplane first pressure jet not zero: '+name)
        stages[name]=dict(even_Z=True,ordinary_Z_at_zero_exact_zero=True)
    for name,offset in operator.partition['offsets'].items():
        if s.diff(offset,z)!=0:raise ValueError('Original radial stage origin depends on Z: '+name)
    midplanes=[row for row in source['actual_original_pressure_point_queries'] if row['original_Z_exact']=='0']
    if len(midplanes)!=1:raise ValueError('Original accepted exact midplane pressure record required')
    jet=midplanes[0]['normalized_pressure_ordinary_Z_jets'][1]
    if not jet['exact_symmetry_zero'] or not jet['finite_end_and_all_late_source_error']['zero']:
        raise ValueError('Original finite-end/late odd-pressure budget not zero')
    return dict(passed=True,actual_original_even_density_stage_identities=stages,
        original_stage_offsets_Z_independent=True,
        finite_axial_end_and_all_late_stages_in_same_pressure_function=True,
        source_midplane_pressure_first_jet_including_remainder_exact_zero=True,
        only_order1_budget_removed_at_exact_Z0=True,
        order0_and_order2_pressure_late_errors_retained=True,
        no_baseline_only_pressure_symmetry_shortcut=True,
        original_pressure_receipt=pressure.RECEIPT,receipt_sha256=sha(pressure.RECEIPT))


class OriginalReferenceMidplaneWholeCells(whole.OriginalReferenceWholeCells):
    """Separate exact Z0 service; the original signed class remains strict."""
    def __init__(self,*,Z=0):
        if point.pressure.exact_Z(Z)!=0:raise ValueError('This original midplane adapter requires exact Z0')
        PROJECTED['__init__'](self,Z=0)
        self.parity=original_pressure_midplane_parity(self.inputs.frame)
        receipt=json.loads((HERE/whole.RECEIPT).read_bytes())
        if not receipt.get('all_passed') or not receipt.get(whole.GATE) or receipt['source_family']!=self.family:
            raise ValueError('Accepted original whole-reference source/integral service required')
        for name,digest in {**receipt['input_hashes'],whole.RECEIPT:sha(whole.RECEIPT)}.items():
            if sha(name)!=digest:raise ValueError('Original midplane dependency changed: '+name)
            if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Original midplane source families disagree')
            self.hashes[name]=digest
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)

    def roots(self,left,right):
        with mp.workdps(self.ctx.dps+40):
            query=PROJECTED['roots'](self,left,right);roots=query['roots'];kernel=query['kernel']
            if not roots['p2'][(0,0)].zero or not kernel.u.zero or kernel.geometry!='small_r_series':
                raise ArithmeticError('Actual original exact-midplane p2=0 source identity required')
            if roots['p2'][(0,1)].zero:
                raise ArithmeticError('Actual original nonzero midplane p2_Z must survive')
            for name in ('E','a','b','t0'):
                if not roots[name][(0,1)].zero:raise ArithmeticError('Original even/constant first jet changed: '+name)
            if not roots['V'][(0,0)].zero or roots['V'][(0,1)].zero:
                raise ArithmeticError('Original V(0)=0 and V_Z=4/Pstar identities required')
            query['midplane_contract']=dict(actual_p2_C0_exact_zero=True,actual_u_C0_exact_zero=True,
                actual_p2_Z_nonzero_source_row_retained=True,ordinary_E_Z_exact_zero=True,
                exact_L_midplane_one=True,pressure_parity=self.parity,
                formal_p2_Z_amplification_not_zeroed_or_capped=True)
            return query

    def cell(self,left,right,*,N,bits=32):
        with mp.workdps(self.ctx.dps+40):
            result=super().cell(left,right,N=N,bits=bits)
            record=result['record']
            if record['source_geometry']!='small_r_series':raise ArithmeticError('Original regular midplane branch required')
            for piece in record['original_inverse_and_Z_piece_records']:
                if piece['ordinary_Z']['branch']!='exact_midplane_nonzero_p2_Z':
                    raise ArithmeticError('Original exact-midplane slow-Z primitive branch required')
            record['exact_midplane_contract']=dict(actual_p2_C0_exact_zero=True,
                actual_p2_Z_nonzero_source_row_retained=True,actual_u_C0_exact_zero=True,
                original_pressure_parity=self.parity,
                original_native_large_Z_derivative_factors_retained=True)
            return result

    def integrate(self,*,count,N):
        with mp.workdps(self.ctx.dps+40):
            result=super().integrate(count=count,N=N);report=result['report']
            report.update(exact_midplane_whole_window_C0_Z_contributions_installed=True,
                original_exact_midplane_slow_Z_functions_used=True,
                signed_implicit_Z_products_bounded_before_multiplication=False,
                exact_midplane_large_derivative_factors_not_capped=True,
                high_precision_source_and_integral_exports_retained=True,
                whole_Z_functional_or_all_route_closure=False)
            return result


def run():
    begin=time.monotonic();owner=OriginalReferenceMidplaneWholeCells();levels=[]
    for count,N in ((4,160),(16,160),(16,320),(16,16384)):
        levels.append(owner.integrate(count=count,N=N)['report'])
        print('Original exact-midplane reference C0/Z:',count,'cells, N',N,flush=True)
    report=dict(**{GATE:True},source_family=owner.family,
        actual_original_midplane_whole_reference_levels=levels,
        original_pressure_midplane_parity=owner.parity,
        source_method_projections=PROJECTION_PROOF,
        original_source_factor_atlas=owner.atlas.record(),
        exact_midplane_C0_and_true_nonzero_Z_primitives_integrated=True,
        original_nonzero_midplane_p2_Z_not_replaced_by_zero=True,
        original_even_pressure_late_errors_retained=True,
        formal_large_native_Z_factors_retained=True,exact_Z0_only=True,
        all_17_chart_or_24_cell_oracle_installed=False,actual_five_controls_installed=False,current_whole_N_selected=False,
        **dict.fromkeys(point.source.inertial.profiles.loop.OPEN,False),input_hashes=owner.hashes,
        execution_seconds=time.monotonic()-begin,
        scope='Actual original exact-Z0 whole Rh_reference[-5,0] C0/Z own-rate contributions. Original full-pressure odd derivative identity only removes P0_Z tail at exact0; true nonzero p2_Z, source scales, even pressure errors, finite-N coefficients, phase and incoming/P0 separation remain. Large Z derivative factors remain formal. Not a Z neighborhood or whole-Z/all-route terminal closure/global N/recursive corrected field.')
    (HERE/NAME).write_bytes(gzip.compress(json.dumps(base.encoded(report),indent=2).encode()+b'\n',compresslevel=9,mtime=0))
    return report


if __name__=='__main__':run()
