"""Restore exact directed inlet atoms from a committed endpoint receipt."""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_interval_difference import IntervalDifference
from lei_ren_part1_paper_interval_pressure_width_jet import IntervalPressureWidthJet
from lei_ren_part1_paper_interval_axial_second_jet import IntervalAxialSecondJet


DEFAULT_RECEIPT=Path(__file__).with_name('lei_ren_part1_paper_interval_pressure_perturbed_core_check.json')


def read_inlet(path=DEFAULT_RECEIPT):
    path=Path(path);record=json.loads(path.read_text());ctx=MPIntervalContext();ctx.dps=record['precision']
    def interval(row):
        return ctx.mpf([mp.make_mpf(tuple(row['lower_exact_mpf_tuple'])),
                        mp.make_mpf(tuple(row['upper_exact_mpf_tuple']))])
    def ring(row):
        atoms={}
        for key,pair in row.items():
            p,w=map(int,key.strip('()').split(','))
            atoms[(p,w)]=IntervalDifference(ctx,interval(pair['nominal']),interval(pair.get('difference',pair.get('perturbation'))))
        return IntervalPressureWidthJet._from_atoms(ctx,atoms,record['pressure_parameter_order'],record['width_order'])
    def item(row):
        if isinstance(row,dict) and all(name in row for name in ('value','tangent','second')):
            return IntervalAxialSecondJet._from_slots(ctx,*(ring(row[name]) for name in ('value','tangent','second')))
        if isinstance(row,dict) and 'lower_exact_mpf_tuple' in row:return interval(row)
        if isinstance(row,dict):return {name:item(value) for name,value in row.items()}
        return row
    inlet={name:item(value) for name,value in record['inlet'].items()};inlet['ctx']=ctx
    return inlet,dict(input_receipt=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                      conditional_pressure_integral_error_propagated=record['conditional_pressure_integral_error_propagated'])
