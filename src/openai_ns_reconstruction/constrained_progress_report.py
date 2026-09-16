"""Build a development comparison report from saved independent validation data."""
import json
from pathlib import Path


def run():
    root=Path('artifacts/constrained')
    files=sorted(root.glob('*/validation.json'))
    files += [root/'initial_pde_validation.json']
    lines=['# Constrained reconstruction: development results','',
        'These are saved development comparisons, not final acceptance. The validation seed has informed model selection; a fresh frozen-candidate audit remains required.',
        '', '| Experiment | Sampled max momentum residual | Sampled structure | Evidence |',
        '| --- | ---: | --- | --- |']
    count=0
    for path in files:
        if not path.exists():continue
        data=json.loads(path.read_text())
        rows=data.get('rows',[])
        if not rows:continue
        finest=min(r['step'] for r in rows)
        maximum=max(r['residual_sampled_max'] for r in rows if r['step']==finest)
        structural=data.get('structure',{}).get('sampled_constraints_pass')
        label='pass' if structural is True else 'fail' if structural is False else 'not recorded in this artifact'
        name=path.parent.name if path.parent!=root else 'initial'
        lines.append(f'| {name} | {maximum:.6f} | {label} | [{path.name}](../{path.as_posix()}) |')
        count+=1
    lines += ['', '## Acceptance still missing','',
        '- Momentum residual thresholds remain 0.001. None of these development runs meets them.',
        '- Sampled structure checks do not prove uniform-in-space/time constraints.',
        '- Global torque closure is not local momentum closure.',
        '- Five symbolic ansatz checks have documented assumptions; time-dependent extensions must not inherit fixed-profile scaling claims.',
        '- Final fresh-sample audit, current-candidate convergence/sensitivity and clean-environment reproduction remain pending.',
        '', 'Regenerate with `python -m openai_ns_reconstruction.constrained_progress_report`.', '']
    path=Path('reports/CONSTRAINED_PROGRESS.md');path.parent.mkdir(exist_ok=True)
    path.write_text('\n'.join(lines),encoding='utf-8')
    print(f'wrote {count} experiment comparisons')

if __name__=='__main__':run()
