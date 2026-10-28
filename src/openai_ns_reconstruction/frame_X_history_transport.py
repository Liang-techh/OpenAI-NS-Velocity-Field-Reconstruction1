"""Source-bound affine Gamma angular-history transport through ordinary Z2.

The defining finite integral and its derivatives remain symbolic functions;
no source functions, integrals, old builders, or numerical bounds are run.
The current X rows must still be proved to form a genuine C2 source jet.
"""
import base64
import hashlib
import json
import sympy as s

P='lei_ren_part1_paper_compliant_current_original_'
SOURCE='experiments/root_st073/'+P+'n1_Gamma_future_functions.py'
SOURCE_SHA='38d3a8e5abd964c42b135410485f6ad58dfe9f7656035b4b5c7b1ed48c2ebb01'
RECEIPT='experiments/root_st073/'+P+'n1_Gamma_future_functions_check.json'
REPORT='experiments/root_st073/lei_ren_part1_paper_compliant_current_frame_axis_X_history_port.json'
REPORT_SHA='ad0136625deff21d15d0ade72a324f43cced95228c363d43c8f100d7585cb04e'
RESULT='tmp/frame_axis_gamma_roles_v1/result.json'
RESULT_SHA='2facb6d999cbbbb386cb8e53945cee571a2a7ab13688bd3c87507aa6cbb67a91'


def hashed(root,name,expected):
    raw=(root/name).read_bytes()
    if hashlib.sha256(raw).hexdigest()!=expected:raise ValueError('hash changed: '+name)
    return raw


def history_jets(X0,X1,X2,J0,J1,J2,P0,P1,P2,C,eq,A,B):
    """Evaluate symbolic jet algebra from supplied function jets, not bounds."""
    D0=C*X0+(1-C)*eq+J0
    D1=C*X1+J1
    D2=C*X2+J2
    F0=P0*D0
    F1=P1*D0+P0*D1
    F2=P2*D0+2*P1*D1+P0*D2
    return {'Xf':(F0,F1,F2),'Xt':(A*F0+B,A*F1,A*F2)}


def finite_J_integrands(z,v,r,sigma):
    """Defining integrands for J and its first two formal Z derivatives."""
    q=(1+z*z)/2
    weight=s.exp(r*v)
    return (weight*q**sigma,
            weight*sigma*z*q**(sigma-1),
            weight*sigma*(q**(sigma-1)+(sigma-1)*z*z*q**(sigma-2)))


def recover(root):
    code=hashed(root,SOURCE,SOURCE_SHA).decode('utf-8').replace('\r\n','\n')
    receipt=json.loads((root/RECEIPT).read_bytes())
    if receipt.get('all_passed') is not True or receipt['input_hashes'].get(SOURCE.split('/')[-1])!=SOURCE_SHA:
        raise ValueError('saved Gamma receipt does not pin source')
    anchors=["q = a.add(f(1), a.mul(z, z))",
             "eq = g.quotient(g.one, r, 'original 0<mu<1')",
             "logratio = a.log(a.scale('1/2', q), 'q/2>=1/2')",
             "F = lambda x: a.exp(a.scale(sigma(g, g.mul(c('1/100'), x)), logratio))",
             "a.add(Xv, flatten_primitive)), F(c(100)), 'positive original flatten factor')",
             "Xrel = a.add(f(eq), a.scale(g.unary('exp', g.neg(g.mul(r, Lrel))), a.sub(Xf, f(eq))))",
             "Xin = a.scale(g.unary('exp', g.mul(c('-1/2'), r)), a.add(Xrel, f(Iin)))",
             "Xs = a.add(Xin, f(Ts))",
             "Xt = a.scale(g.unary('exp', g.mul(c('-1/2'), k)), a.add(Xs, f(Iout)))",
             "variable=z.value.node, value=g.zero.node, Z_independent_substitution=True)"]
    if any(a not in code for a in anchors):raise ValueError('Gamma source definition anchor missing')
    report=json.loads(hashed(root,REPORT,REPORT_SHA)); raw=base64.b64decode(report['artifacts'][RESULT]['base64'],validate=True)
    if hashlib.sha256(raw).hexdigest()!=RESULT_SHA:raise ValueError('saved Gamma-role result changed')
    saved=json.loads(raw)
    if saved.get('status')!='returned':raise ValueError('Gamma-role result did not return')
    evidence=saved['result']; roles={x['history']:x for x in evidence['axis_roles']}
    if not all(roles[n]['same_saved_reference'] for n in ('Xt','Xf')):raise ValueError('saved axis role mismatch')
    if evidence['gamma_native_binding']['imported_X_C1']!=[31655,31678]:raise ValueError('saved Gamma X inlet differs')
    ports=[row for row in evidence['axis_port_roles'] if row['key']=='X' and row['input'].endswith('n1_pulse_C3_functions.json.gz')]
    if len(ports)!=1 or len(ports[0]['rows'])<3:raise ValueError('saved X C2 port rows missing')
    jet_ports=ports[0]['rows'][:3]
    if [row['current_port_ref'] for row in jet_ports]!=[265772,265776,265782] or not all(row['direct_saved_mapping_matches'] for row in jet_ports):
        raise ValueError('saved X port mapping differs')
    # Derive the explicit source formula using sigma(1)=1, q=(1+Z^2)/2.
    z,r,mu,k,L,eq,Iin,Iout,Ts=s.symbols('z r mu k L eq Iin Iout Ts',real=True)
    X=s.Function('X')(z); J=s.Function('J')(z)
    C=s.exp(-13*r/mu); H=s.exp(-r*L); Pz=2*s.exp(-100*r)/(1+z*z)
    Xv=eq+C*(X-eq); Xf=Pz*(Xv+J)
    Xrel=eq+H*(Xf-eq)
    Xt=s.exp(-k/2)*(s.exp(-r/2)*(Xrel+Iin)+Ts+Iout)
    A=s.exp(-k/2)*s.exp(-r/2)*H
    B=s.exp(-k/2)*(s.exp(-r/2)*((1-H)*eq+Iin)+Ts+Iout)
    abstract=history_jets(X,s.diff(X,z),s.diff(X,z,2),J,s.diff(J,z),s.diff(J,z,2),
                          Pz,s.diff(Pz,z),s.diff(Pz,z,2),C,eq,A,B)
    v,sigma=s.symbols('v sigma',real=True)
    integrands=finite_J_integrands(z,v,r,sigma)
    checks={'Xt_affine_Xf':s.simplify(Xt-(A*Xf+B))==0,
            'Xf_value':s.simplify(abstract['Xf'][0]-Xf)==0,
            'Xf_Z1':s.simplify(abstract['Xf'][1]-s.diff(Xf,z))==0,
            'Xf_Z2':s.simplify(abstract['Xf'][2]-s.diff(Xf,z,2))==0,
            'Xt_Z1':s.simplify(abstract['Xt'][1]-s.diff(Xt,z))==0,
            'Xt_Z2':s.simplify(abstract['Xt'][2]-s.diff(Xt,z,2))==0,
            'P_Z1':s.simplify(s.diff(Pz,z)+2*z*Pz/(1+z*z))==0,
            'P_Z2':s.simplify(s.diff(Pz,z,2)-2*(3*z*z-1)*Pz/(1+z*z)**2)==0,
            'axis_P_Z1':s.diff(Pz,z).subs(z,0)==0,
            'J_integrand_Z1':s.simplify(s.diff(integrands[0],z)-integrands[1])==0,
            'J_integrand_Z2':s.simplify(s.diff(integrands[0],z,2)-integrands[2])==0}
    return {'schema':'frame_X_history_transport_v1','status':'source_bound_conditional_history_jets',
            'checks':checks,'checks_passed':all(checks.values()),'axis_role_references':evidence['axis_roles'],
            'native_X_C1':evidence['gamma_native_binding']['imported_X_C1'],'current_X_C2_port_rows':jet_ports,
            'source_definition':{'r':'1-mu','k':'1-delta/2','L':'-30 log(mu)','eq':'1/r',
              'q':'(1+Z^2)/2','C':'exp(-13*r/mu)','P':'2*exp(-100*r)/(1+Z^2)',
              'J':'integral_0^100 exp(r*v)*q(Z)^sigma(v/100) dv; defining integral only',
              'J_Z1_integrand':'exp(r*v)*sigma(v/100)*Z*q^(sigma(v/100)-1)',
              'J_Z2_integrand':'exp(r*v)*sigma(v/100)*(q^(sigma(v/100)-1)+(sigma(v/100)-1)*Z^2*q^(sigma(v/100)-2))',
              'D':'C*X+(1-C)*eq+J','Xf':'P*D','Xt':'A*Xf+B',
              'A':'exp(-k/2)*exp(-r/2)*exp(-r*L)',
              'B':'exp(-k/2)*(exp(-r/2)*((1-exp(-r*L))*eq+Iin)+Ts+Iout)'},
            'conditions':['mu>0 and r!=0','sigma(1)=1 with stated flat source semantics',
              'mu,delta,Iin,Iout,Ts,L independent of Z','current X rows are same-source ordinary C2 jet',
              'J1 and J2 are ordinary derivatives of the defining finite integral',
              'source graphs share the claimed function semantics and axis substitution is valid'],
            'actual_current_X_C2_jet_proved':False,'actual_inlet_substitution_transport_proved':False,
            'actual_Gamma_improper_integral_transport_proved':False,'actual_axis_seam_proved':False,
            'higher_X4_X5_compatibility_proved':False,'source_values_computed':False,
            'integrals_evaluated':False,'old_errors_removed':False,'runtime_model_metadata':None}
