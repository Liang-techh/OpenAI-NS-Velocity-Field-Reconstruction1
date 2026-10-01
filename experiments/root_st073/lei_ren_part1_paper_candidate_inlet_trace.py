"""Directed candidate angular inlet trace without fitted trace resets."""
import argparse
import hashlib
import json
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_candidate_combined_core_budget import DEFAULT_STATE
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode


def run(s_exit='4',state_name=DEFAULT_STATE,output_name=None):
    base=Path(__file__).parent;path=base/state_name
    raw=path.read_bytes();state=json.loads(raw)
    for name,digest in state['source_hashes'].items():
        name=({'driver':state['target'].get('driver_file','lei_ren_part1_paper_candidate_gauge_core.py'),
               'pressure_input':state['target']['pressure_file']}).get(name,name)
        if hashlib.sha256((base/name).read_bytes()).hexdigest()!=digest:
            raise AssertionError('Candidate state dependency changed: '+name)
    ctx=MPIntervalContext();ctx.dps=state['target']['precision']
    with mp.workdps(ctx.dps+40):
        def read(v):
            return ctx.mpf([mp.make_mpf(tuple(v['lower_exact_mpf_tuple'])),
                            mp.make_mpf(tuple(v['upper_exact_mpf_tuple']))])
        s=ctx.mpf(s_exit);z=ctx.mpf(state['target']['Z'])
        if endpoints(s)[0]<=0 or endpoints(s)[1]>mp.mpf('4.1'):
            raise ValueError('Exit must lie in (0,4.1]')
        lam=ctx.mpf(state['target']['Lambda']);eps=1/lam
        dt=ctx.mpf(state['target']['delta']);L=1-dt*z*z;d=1-z*z
        A=[read(row[0])/lam**n for n,row in enumerate(state['A_rows'])]
        Az=[read(row[1])/lam**n for n,row in enumerate(state['A_rows'])]
        U=[read(row[0])/lam**n for n,row in enumerate(state['Uz_rows'])]
        Uz=[read(row[1])/lam**n for n,row in enumerate(state['Uz_rows'])]
        g=-read(state['fixed_jets']['ell_Z_taylor'][0])/lam
        zero=ctx.mpf(0)
        def add(*polys):
            return [sum((p[n] if n<len(p) else zero for p in polys),zero)
                    for n in range(max(map(len,polys)))]
        def scale(p,c):return [v*c for v in p]
        def mul(a,b):
            return [sum((a[i]*b[n-i] for i in range(len(a))
                         if 0<=n-i<len(b)),zero) for n in range(len(a)+len(b)-1)]
        def value(p):
            out=zero
            for v in reversed(p):out=out*s+v
            return out
        W=add([ctx.mpf(1)],[-((1-dt)*z*u+d*uz)/(n+1)
                              for n,(u,uz) in enumerate(zip(U,Uz))])
        H=add([(1-dt)*z/2],scale(U,d))
        # Full original angular RHS divided by Lambda*L*F0.
        rhs=add(scale(mul(W,[(n+1)*v for n,v in enumerate(A)]),eps/L),
                scale(mul(H,Az),eps/L),
                scale(add(A,scale(mul(A,U),-2*z)),eps*dt/(2*L)),
                scale(mul(H,A),-g/L))
        lhs=[2*(n+1)*(n+2)*A[n+1] for n in range(len(A)-1)]
        defect=add(lhs,scale(rhs,-1))
        integrated=sum((v*s**(n+1)/(n+2) for n,v in enumerate(defect)),zero)
        direct=2*sum((n*v*s**n for n,v in enumerate(A)),zero)-sum(
            (v*s**(n+1)/(n+2) for n,v in enumerate(rhs)),zero)
        lo,hi=endpoints(integrated);dl,dh=endpoints(direct)
        if max(lo,dl)>min(hi,dh):
            raise AssertionError('Integral and direct inlet identities disagree')
        angular_value=value(A)
        if endpoints(angular_value)[0]<=0:
            raise AssertionError('Finite inlet angular denominator is not positive')
        ratio=integrated/angular_value
        report=dict(state_sha256=hashlib.sha256(raw).hexdigest(),state_file=path.name,
            accepted_schedule_sha256=state['accepted_schedule_sha256'],
            Lambda=state['target']['Lambda'],Z=state['target']['Z'],
            radial_degree=state['completed_radial_order'],scaled_exit=s,
            finite_Phi_at_exit=angular_value,
            finite_Ttheta_over_F0_integral=integrated,
            finite_Ttheta_over_F0_direct=direct,
            finite_Ttheta_over_F=ratio,
            two_directed_identity_enclosures_overlap=True,
            finite_trace_interval_contains_zero=endpoints(ratio)[0]<=0<=endpoints(ratio)[1],
            interval_containment_is_not_matching_certificate=True,
            analytic_core_angular_trace_zero_by_exact_equation=True,
            finite_trace_reset_or_fit_used=False,whole_axis_trace_certified=False,
            core_to_collar_matching_certified=False,temporal_recursion=False)
        output=base/output_name if output_name else Path(__file__).with_suffix('.json')
        output.write_text(json.dumps(encode(report),indent=2)+'\n',encoding='utf-8')
        print('Candidate inlet degree',state['completed_radial_order'],
              'Ttheta/F interval',mp.nstr(endpoints(ratio)[0],18),
              mp.nstr(endpoints(ratio)[1],18),flush=True)
        return report


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--scaled-exit',default='4')
    parser.add_argument('--state-file',default=DEFAULT_STATE)
    parser.add_argument('--output-name')
    args=parser.parse_args();run(args.scaled_exit,args.state_file,args.output_name)
