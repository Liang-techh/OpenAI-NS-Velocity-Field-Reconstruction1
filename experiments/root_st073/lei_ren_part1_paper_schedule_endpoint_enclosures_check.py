"""Actual stored-schedule endpoints and pre-collar mass upper bounds."""
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_continuous_preheat_pressure_check import source_profile
from lei_ren_part1_paper_schedule_endpoint_enclosures import ScheduleEndpointEnclosures,endpoints


def encode(value):
    if hasattr(value,'_mpi_'):
        lo,hi=endpoints(value)
        return dict(lower=mp.nstr(lo,90),upper=mp.nstr(hi,90),width=mp.nstr(hi-lo,60),
                    lower_exact_mpf_tuple=list(value._mpi_[0]),upper_exact_mpf_tuple=list(value._mpi_[1]))
    if isinstance(value,mp.mpf):return dict(value=mp.nstr(value,90),log_abs=None if value==0 else mp.nstr(mp.log(abs(value)),70),exact_mpf_tuple=list(value._mpf_))
    if isinstance(value,dict):return {k:encode(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [encode(v) for v in value]
    return value


def run():
    s=source_profile().schedule;e=ScheduleEndpointEnclosures(s)
    with mp.workdps(e.precision+30):
        stages=[stage for stage in s._stage_bounds if stage not in ('reference','heat_connection','exact_heat')]
        rows={stage:e.stage_pressure_upper(stage) for stage in stages}
        checkpoints={name:e.log_amplitude_ratio(y) for name,y in s.checkpoint_offsets.items()}
        max_width=max(endpoints(v['interval'])[1]-endpoints(v['interval'])[0] for v in checkpoints.values())
        assert max_width<mp.mpf('1e-200'), 'Endpoint primitive intervals unresolved'
        # Every mass bound is nonnegative; tiny flatten mass stays separate.
        assert all(row['mass_upper']>=0 for row in rows.values())
        report=dict(precision=e.precision,axial_radius='.8',stages=rows,checkpoints=checkpoints,
            maximum_endpoint_log_width=max_width,endpoint_primitive_quadrature_used=False,
            directed_interval_arithmetic=True,stored_schedule_scope=True,
            transcendental_parameter_errors_enclosed=False,heat_collar_included=False,
            full_pressure_error_enclosed=False,source_C2_norm_enclosed=False)
        Path(__file__).with_suffix('.json').write_text(json.dumps(encode(report),indent=2)+'\n',encoding='utf-8')
        print('actual pre-collar stage bounds saved',len(rows),'endpoint max width',mp.nstr(max_width,12),flush=True)


if __name__=='__main__':run()
