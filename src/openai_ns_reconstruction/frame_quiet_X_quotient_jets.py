"""Quiet-frame X=h/E as a same-input quotient jet, through Z order five.

This uses saved port references and exact source definitions; it does not
evaluate source functions, traverse graph edges or call ancestor builders.
The incoming jet and Z-independent quiet geometry remain explicit premises.
"""
import gzip
import hashlib
import json
from math import comb
import sympy as s

P='lei_ren_part1_paper_compliant_current_original_'
RP='experiments/root_st073/'+P+'Rp_native_C5_pulse_source.json.gz'
RP_SHA='442df70e5ab4051066b241e943c8d0531067d057e774a04dcc341dadbd703ebe'
SOURCE='experiments/root_st073/'+P+'true_quiet_Rp_C5_pressure_C3_functions.py'
SOURCE_SHA='0db6ea1a97184662c8e80f011844c12d143eeeb11333510f2fbbb8dcc19de068'
RECEIPT='experiments/root_st073/'+P+'Rp_native_C5_pulse_source_check.json'
RECEIPT_SHA='320f4623859c4ad64a4b5c6c4912e43e42651b58d8359b334b96b5b0d6d38d15'
QUIET=P+'true_quiet_Rp_C5_pressure_C3_functions.json.gz'


def quotient_jets(A,h):
    """Ordinary derivative quotient recurrence, supplied same-source jets."""
    if len(A)!=len(h) or not A:raise ValueError('matching nonempty jet lengths required')
    if A[0]==0:raise ValueError('nonzero inlet E required')
    q=[]
    for n in range(len(A)):
        q.append((h[n]-sum(comb(n,k)*A[k]*q[n-k] for k in range(1,n+1)))/A[0])
    return q


def quiet_X_jets(A,h,T,K):
    """Apply the Z-independent quiet propagation to the quotient jet."""
    q=quotient_jets(A,h)
    return [T*q[0]+K]+[T*x for x in q[1:]]


def hashed(root,path,expected):
    raw=(root/path).read_bytes()
    if hashlib.sha256(raw).hexdigest()!=expected:raise ValueError('hash mismatch: '+path)
    return raw


def recover(root):
    code=hashed(root,SOURCE,SOURCE_SHA).decode('utf-8').replace('\r\n','\n')
    anchors=["E=a.scale(A,g.unary('exp',g.neg(g.mul(alpha,s))))",
             "theta=g.mul(decay('3/2'),s,g.unary('exprel',g.mul(g.sub(g.one,mu),s)))",
             "h=a.add(a.scale(a.trim(own['h']),decay('3/2')),a.scale(A,theta))",
             "X=div(H['h'],E)","self.flow=power_flow(g,self.incoming['modified_fields']['E'],self.incoming['modified_fields'],self.mu,s,5)"]
    if any(x not in code for x in anchors):raise ValueError('quiet source anchor missing')
    receipt=json.loads(hashed(root,RECEIPT,RECEIPT_SHA))
    inherited=receipt['checks']['compact_quiet']
    if receipt.get('all_passed') is not True or inherited.get('passed') is not True:
        raise ValueError('prior quiet receipt not accepted')
    labels=inherited['actual_eight_Rp_frame_C5_rows']
    if not all('X_Z'+str(n) in labels for n in range(6)):
        raise ValueError('prior receipt lacks X row labels')
    archive=json.loads(gzip.decompress(hashed(root,RP,RP_SHA)))
    frame=archive['Rp_native_C5_pulse_source_bindings']; quiet=frame['inputs'][QUIET]
    incoming=quiet['incoming']['modified_fields']; A=incoming['E']; h=incoming['h']
    old_frame=quiet['Rp_frame']['X']; current=frame['current_frame']['X']
    if not all(isinstance(q,list) and len(q)>=6 for q in (A,h,old_frame,current)):
        raise ValueError('missing six saved source rows')
    if old_frame!=current:raise ValueError('quiet X and current frame references differ')
    if current[:3]!=[265772,265776,265782]:raise ValueError('saved current X port differs')
    # Formal Leibniz identities and independent polynomial quotient checks.
    z=s.symbols('z'); av=(s.Symbol('a0',nonzero=True),)+s.symbols('a1:6'); hv=s.symbols('h0:6')
    q=quotient_jets(av,hv)
    checks={}
    for n in range(6):
        checks['leibniz_'+str(n)]=s.cancel(sum(comb(n,k)*av[k]*q[n-k] for k in range(n+1))-hv[n])==0
    ap=sum(av[n]*z**n/s.factorial(n) for n in range(6))
    hp=sum(hv[n]*z**n/s.factorial(n) for n in range(6))
    for n in (0,1,2):
        checks['direct_quotient_'+str(n)]=s.cancel(s.diff(hp/ap,z,n).subs(z,0)-q[n])==0
    mu,quiet_s,G=s.symbols('mu quiet_s G',real=True)
    source_ratio=(hv[0]+av[0]*quiet_s*G)*s.exp(-s.Rational(3,2)*quiet_s)/(av[0]*s.exp(-(s.Rational(1,2)+mu)*quiet_s))
    reduced_ratio=s.exp(-(1-mu)*quiet_s)*(hv[0]/av[0]+quiet_s*G)
    checks['quiet_source_ratio']=s.simplify(source_ratio-reduced_ratio)==0
    T,K=s.symbols('T K',real=True); propagated=quiet_X_jets(av,hv,T,K)
    checks['propagated_Z0']=s.cancel(propagated[0]-T*q[0]-K)==0
    checks['propagated_Z1']=s.cancel(propagated[1]-T*q[1])==0
    checks['propagated_Z2']=s.cancel(propagated[2]-T*q[2])==0
    checks['propagated_Z4']=s.cancel(propagated[4]-T*q[4])==0
    checks['propagated_Z5']=s.cancel(propagated[5]-T*q[5])==0
    return {'schema':'frame_quiet_X_quotient_jets_v1','status':'source_bound_conditional_quiet_quotient_jets',
            'checks':checks,'checks_passed':all(checks.values()),
            'source_roles':{'incoming_E':A,'incoming_h':h,'quiet_X':old_frame,'current_X':current,
                            'mu_ref':quiet.get('mu'),'quiet_coordinate_ref':quiet.get('s'),'length_ref':quiet.get('length')},
            'inherited_receipt_scope':{'compact_quiet_passed':True,'saved_X_labels':['X_Z'+str(n) for n in range(6)],
                                       'receipt_evidence_reexecuted':False},
            'formula':{'lambda':'1-mu','W':'h_in/A','T':'exp(-lambda*s)',
                       'K':'T*s*exprel(lambda*s)','X':'T*W+K',
                       'X_j':'T*W_j for j>=1 when mu and s are Z-independent'},
            'conditions':['A is nonzero and incoming A,h are ordinary C5 jets of the same source functions',
                          'mu and quiet s/endpoint length are Z-independent',
                          'saved source-row meanings and fixed-coordinate substitution are preserved',
                          'exprel is its original entire continuation'],
            'actual_incoming_C5_jet_proved':False,'actual_geometry_Z_independence_proved':False,
            'actual_current_X_C2_jet_proved_by_this_leaf':False,'actual_X4_X5_zero_proved':False,
            'actual_inlet_substitution_transport_proved':False,'physical_endpoint_proved':False,
            'source_values_computed':False,'old_errors_removed':False,'runtime_model_metadata':None}
