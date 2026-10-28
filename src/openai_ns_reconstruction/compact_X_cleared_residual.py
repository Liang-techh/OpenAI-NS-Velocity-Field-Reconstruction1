"""Denominator-cleared ordinary-Z tests for the actual compact-exit X port.

The residuals are exact functions of supplied same-source jets, not bounds
or sampled values. Source-semantic, terminal-zero and positivity gates remain.
"""
import base64
import gzip
import hashlib
import json
from math import comb
import sympy as s

P='lei_ren_part1_paper_compliant_current_original_'
RP='experiments/root_st073/'+P+'Rp_native_C5_pulse_source.json.gz'
RP_SHA='442df70e5ab4051066b241e943c8d0531067d057e774a04dcc341dadbd703ebe'
COMPACT=P+'true_compact_C5_functions.json.gz'
QUIET=P+'true_quiet_Rp_C5_pressure_C3_functions.json.gz'
OUTER=P+'true_outer_C6_C5_functions.json.gz'
COMPACT_SRC='experiments/root_st073/'+P+'true_compact_C5_functions.py'
COMPACT_SHA='f0b0333b19a0e2504d023b97e841033d3cbcc23427c57ad9672264866c01ec7b'
OUTER_SRC='experiments/root_st073/'+P+'true_outer_C6_C5_functions.py'
OUTER_SHA='4464a934d303a95a0da157e9e289f6ed8a04d5ff0adcb0d663857f0d9c57dcd3'
CHECK='experiments/root_st073/'+P+'Rp_native_C5_pulse_source_check.json'
CHECK_SHA='320f4623859c4ad64a4b5c6c4912e43e42651b58d8359b334b96b5b0d6d38d15'
REPORT='experiments/root_st073/lei_ren_part1_paper_compliant_current_frame_quiet_X_quotient_jets.json'
REPORT_SHA='202472d2e002fc06abbdb208a25737ff97e50e6bf278366db33e8273e3f3349f'
RESULT='tmp/frame_quiet_X_quotient_jets_v1/result.json'
RESULT_SHA='7da0f20d6e40046f5779a13c202f0c37520a92007d92fcf3b64c8554ea79ed87'


def cleared_ratio_jets(E,h):
    """Return P_n=E_0^(n+1) d_Z^n(h/E), with no division."""
    if not E or len(E)!=len(h):raise ValueError('matching nonempty ordinary jet arrays required')
    out=[h[0]]
    for n in range(1,len(E)):
        out.append(E[0]**n*h[n]-sum(comb(n,k)*E[k]*E[0]**(k-1)*out[n-k] for k in range(1,n+1)))
    return out


def hashed(root,name,expected):
    raw=(root/name).read_bytes()
    if hashlib.sha256(raw).hexdigest()!=expected:raise ValueError('changed input '+name)
    return raw


def recover(root):
    cs=hashed(root,COMPACT_SRC,COMPACT_SHA).decode('utf-8').replace('\r\n','\n')
    os=hashed(root,OUTER_SRC,OUTER_SHA).decode('utf-8').replace('\r\n','\n')
    anchors=["self.incoming_background=outer.Rc_background['history'];self.incoming_correction=outer.Rc['correction']",
             "self.modified_fields.update(E=f[0]['E'],V=f[0]['V'])",
             "self.modified_fields[k]=a.add(a.trim(fields[k]),self.history[k])",
             "self.terminal=self.at(g.constant(2));self.inlet=self.at(g.one)",
             "absolute_leading_histories_and_pressure_not_zeroed=True"]
    if any(x not in cs for x in anchors):raise ValueError('compact source anchor missing')
    if "canonical_normalized_Rh='E=exp(-1/2)/(1+Z^2); V=m=4Z/Pstar; h=5E/8;" not in os:
        raise ValueError('Rh seed source anchor missing')
    receipt=json.loads(hashed(root,CHECK,CHECK_SHA)); cq=receipt['checks']['compact_quiet']
    terminal=cq['actual_terminal_integral_and_unique_branch_source_proof']
    if not (receipt.get('all_passed') and cq.get('passed') and terminal.get('passed')):
        raise ValueError('inherited compact source receipt not accepted')
    report=json.loads(hashed(root,REPORT,REPORT_SHA)); result_raw=base64.b64decode(report['artifacts'][RESULT]['base64'],validate=True)
    if hashlib.sha256(result_raw).hexdigest()!=RESULT_SHA:raise ValueError('quiet quotient result hash changed')
    prior=json.loads(result_raw)
    if prior.get('status')!='returned':raise ValueError('quiet quotient leaf did not return')
    rp=json.loads(gzip.decompress(hashed(root,RP,RP_SHA))); frame=rp['Rp_native_C5_pulse_source_bindings']
    compact=frame['inputs'][COMPACT];quiet=frame['inputs'][QUIET];outer=frame['inputs'][OUTER]
    exit=compact['terminal'];incoming=quiet['incoming']
    roles={}
    for key in ('E','h'):
        modified=exit['modified_fields'][key]
        if incoming['modified_fields'][key]!=modified:raise ValueError('quiet inlet differs from compact exit: '+key)
        if modified!=prior['result']['source_roles']['incoming_'+('E' if key=='E' else 'h')]:raise ValueError('prior quiet port mismatch: '+key)
        roles[key]={'modified':modified,'background':exit['background_fields'][key],
                    'relative_history':exit['correction'].get(key)}
    seed=outer['seed_identity']
    # Independent polynomial proof: Leibniz construction versus formal D.
    E=s.symbols('E0:6');h=s.symbols('h0:6'); pn=cleared_ratio_jets(E,h)
    formal=[h[0]]
    for n in range(5):
        diff=sum(s.diff(formal[-1],E[j])*E[j+1]+s.diff(formal[-1],h[j])*h[j+1] for j in range(n+1))
        formal.append(s.expand(E[0]*diff-(n+1)*E[1]*formal[-1]))
    checks={'independent_cleared_recurrence_'+str(n):s.expand(pn[n]-formal[n])==0 for n in range(6)}
    Eb,hb,de,dh=s.symbols('Eb hb de dh')
    actual_E=Eb+de;actual_h=hb+dh
    checks['signed_correction_defect']=s.expand(Eb*actual_h-hb*actual_E-(Eb*dh-hb*de))==0
    checks['terminal_defect_if_both_zero']=s.expand((Eb*actual_h-hb*actual_E).subs({de:0,dh:0}))==0
    # Common Z-dependent amplitude a(Z), radial scalars F,G.
    z=s.symbols('z');a=s.Function('a')(z);F,G=s.symbols('F G',nonzero=True)
    ratio=(a*G)/(a*F)
    checks.update({'common_factor_Z'+str(n):s.simplify(s.diff(ratio,z,n))==0 for n in range(1,6)})
    return {'schema':'compact_X_cleared_residual_v1','status':'source_bound_conditional_cleared_test',
            'checks':checks,'checks_passed':all(checks.values()),'compact_exit_roles':roles,
            'current_X_refs':prior['result']['source_roles']['current_X'],
            'cleared_residual_formula':'P0=h0; Pn=E0^n h_n - sum_(k=1)^n binom(n,k) E_k E0^(k-1) P_(n-k)',
            'P4_polynomial':str(s.expand(pn[4])),'P5_polynomial':str(s.expand(pn[5])),
            'signed_correction_defect':'E_bg*h_modified-h_bg*E_modified = E_bg*delta_h-h_bg*delta_E',
            'inherited_source_evidence':{'relative_terminal_zero_flag':terminal.get('actual_relative_terminal_functions_zero_from_integrals_and_unique_branch'),
                'velocity_flatness_flag':cq.get('actual_terminal_velocity_flatness'),
                'absolute_histories_not_zeroed_flag':terminal.get('absolute_histories_and_interior_relative_histories_not_zeroed'),
                'rh_seed_description':seed.get('canonical_normalized_Rh'),
                'rh_whole_Z_identity_receipt':seed.get('original_whole_Z_leading_identity_receipt'),
                'reexecuted':False},
            'conditions':['E,h are same-source ordinary C5 jets and E is nonzero on the full Z domain',
                'background common-factor propagation holds as an all-Z source-function identity',
                'terminal delta_E and delta_h are zero as functions, not merely overlapping bounds',
                'quiet geometry and mu are Z-independent and source semantics are preserved'],
            'actual_P4_zero_proved':False,'actual_P5_zero_proved':False,
            'actual_source_jet_proved':False,'actual_terminal_function_identity_proved':False,
            'actual_X4_X5_compatibility_proved':False,'source_values_computed':False,
            'old_errors_removed':False,'runtime_model_metadata':None}
