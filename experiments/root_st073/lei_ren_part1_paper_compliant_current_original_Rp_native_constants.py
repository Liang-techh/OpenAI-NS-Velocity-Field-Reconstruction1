"""Current complete Rp functions as source-bound native constructor constants.

The exact function view is separate from the old numerical enclosures.
No shared native owner is mutated and no interval endpoint is selected.
"""
import ast
import json
from pathlib import Path
import time
import sympy as s
import lei_ren_part1_paper_compliant_current_original_C3_Rp_native_frame as frame

HERE,PREFIX,sha=frame.HERE,frame.PREFIX,frame.sha
NAME=PREFIX+'current_original_Rp_native_constants.json'
RECEIPT=PREFIX+'current_original_Rp_native_constants_check.json'
GATES=('current_original_Rp_native_cached_constant_function_identity_installed',)
OPEN=frame.OPEN
KEYS=frozenset(('U','M','K','E_Q','E_Z','C1','C2','C0','C_E','H','Pin','Xp'))
LEGACY_KERNEL_ALIASES={'K_B_at_yd':'K_B_at_cutoff','K_B2_at_yd':'K_B2_at_cutoff'}


def exact_expression(source,bindings):
    """Read checked scalar formulas through a closed arithmetic grammar."""
    tree=ast.parse(source,mode='eval')
    allowed=(ast.Expression,ast.BinOp,ast.UnaryOp,ast.Add,ast.Sub,ast.Mult,
        ast.Div,ast.Pow,ast.USub,ast.UAdd,ast.Call,ast.Name,ast.Load,ast.Constant)
    for node in ast.walk(tree):
        if not isinstance(node,allowed):raise ValueError('Unsupported native constant expression')
        if isinstance(node,ast.Name) and node.id not in bindings and node.id!='exp':
            raise ValueError('Unbound native constant symbol '+node.id)
        if isinstance(node,ast.Call) and (not isinstance(node.func,ast.Name)
                or node.func.id!='exp' or len(node.args)!=1 or node.keywords):
            raise ValueError('Only exact exponential calls admitted')
        if isinstance(node,ast.Constant) and type(node.value) is not int:
            raise ValueError('Only exact integer scalar literals admitted')
    return s.sympify(source,locals={**bindings,'exp':s.exp})


def kernel_functions(owner,q):
    """Same defining continuous primitives at their actual source endpoints."""
    charts=owner.frame.bridge.leading.functions['charts'];out={}
    for chart,names in (('O2_slope',{'J1':'J','I_theta':'theta','I_energy':'energy','I_pressure':'pressure'}),
            ('O3_transition',{'J_transition':'J','K_theta':'theta','K_energy':'energy','K_pressure':'pressure'})):
        source=charts[chart];variable=q.at(source['native_coordinate'])
        for name,key in names.items():out[name]=q.at(source['kernels'][key]).subs(variable,1)
    source=charts['O2_axial'];variable=q.at(source['kernels']['original_y_endpoint'])
    for name,key in (('K_B_at_cutoff','original_B_mass'),('K_B2_at_cutoff','original_B_squared_mass')):
        out[name]=q.at(source['kernels'][key]).subs(variable,s.exp(40))
    if len(out)!=10:raise ValueError('All ten original continuous primitives required')
    return out


class CurrentOriginalRpNativeConstants:
    def __init__(self,source_frame=None,require_checked=True):
        self.frame=source_frame if source_frame is not None else frame.CurrentOriginalC3RpNativeFrame()
        if not self.frame.acceptance_loaded:raise ValueError('Accepted complete current Rp frame required')
        self.identity=self.frame.identity;self.hashes=dict(self.frame.hashes)
        receipt=frame.bridge.outer.rh.read(frame.RECEIPT)
        if not receipt['all_passed'] or not all(receipt[key] for key in frame.GATES):
            raise ValueError('All current Rp function gates required')
        self.hashes.update(receipt['input_hashes'])
        self.hashes.update({frame.NAME:sha(frame.NAME),frame.RECEIPT:sha(frame.RECEIPT),
            Path(__file__).name:sha(Path(__file__).name)})
        for name,digest in self.hashes.items():
            if sha(name)!=digest:raise ValueError('Changed native constant prerequisite '+name)
        self.acceptance_loaded=False
        if require_checked:
            proof=frame.bridge.outer.rh.read(RECEIPT)
            if not proof['all_passed'] or not all(proof[key] for key in GATES):
                raise ValueError('Checked actual cached constant identity required')
            for name,digest in proof['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed actual native constants '+name)
            self.acceptance_loaded=True

    def current_constants(self):
        q,data=self.frame.symbolic_frame();z=q.z;U=data['u'].subs(z,0)
        C1=s.diff(data['m1'],z).subs(z,0);C2=s.diff(data['m2'],z).subs(z,0)
        C0=data['energy'].subs(z,0);CE=s.diff(data['energy'],z,2).subs(z,0)/2
        Xp=data['X'].subs(z,0);Pin=data['Mp'].subs(z,0);Pstar=s.exp(q.logP)
        constants=dict(U=U,M=Pstar*U*C1,K=Pstar*U*U*C2,E_Q=U*U*C0,E_Z=U*U*CE,
            C1=C1,C2=C2,C0=C0,C_E=CE,H=U*Xp,Pin=Pin,Xp=Xp)
        constants={key:s.cancel(value) for key,value in constants.items()}
        if set(constants)!=KEYS or any(z in value.free_symbols for value in constants.values()):
            raise ValueError('All native constructor constants must be Z-independent')
        return q,constants

    def original_native_functions(self,q):
        """Accepted old source formulas, with every primitive bound to its function."""
        old=self.frame.reports[frame.native.NAME]['underlying_function_projection']
        if not old['passed'] or old['identity_count']!=45:
            raise ValueError('Actual old source assignment projection required')
        kernels=kernel_functions(self,q)
        Tw=q.at(self.frame.functions['Tw'])
        bindings={**kernels,'Pstar':s.exp(q.logP),'mu':q.mu,'Tw':Tw,'E':s.E}
        # Historical recipe symbols denote cutoff values, before the exp(-11)
        # suffix. Keep the alias local; public keys name the actual endpoint.
        bindings.update({old:kernels[new] for old,new in LEGACY_KERNEL_ALIASES.items()})
        raw={key:exact_expression(value,bindings) for key,value in old['actual_canonical_coefficients'].items()}
        if set(raw)!=set(('U','M','K','H','Pin','EZ','EQ')):
            raise ValueError('Complete original native coefficient formulas required')
        U=raw['U'];Pstar=s.exp(q.logP)
        values=dict(U=U,M=raw['M'],K=raw['K'],H=raw['H'],Pin=raw['Pin'],
            E_Z=raw['EZ'],E_Q=raw['EQ'],C1=raw['M']/(Pstar*U),
            C2=raw['K']/(Pstar*U*U),C0=raw['EQ']/(U*U),C_E=raw['EZ']/(U*U),Xp=raw['H']/U)
        return values,kernels


def run(source_frame=None):
    began=time.monotonic();owner=CurrentOriginalRpNativeConstants(source_frame,require_checked=False)
    q,constants=owner.current_constants();native,kernels=owner.original_native_functions(q)
    report=dict(candidate_current_native_cached_constant_view_constructed=True,source_family=owner.identity,
        actual_current_native_constants={key:s.sstr(value) for key,value in constants.items()},
        actual_original_native_defining_functions={key:s.sstr(value) for key,value in native.items()},
        actual_original_kernel_endpoint_functions={key:s.sstr(value) for key,value in kernels.items()},
        historical_recipe_kernel_aliases=LEGACY_KERNEL_ALIASES,
        turnoff_cutoff='exp(40)',true_yd='exp(40)+11',
        native_M_K_are_raw_and_current_common_histories_already_divided_by_Pstar=True,
        old_interval_constants_not_replaced_or_selected=True,
        selected_native_pulse_constructor_consumes_current_C3_frame=False,
        native_interval_inlet_callback_installed=False,current_numeric_point_field_oracle_installed=False,
        **dict.fromkeys(GATES+OPEN,False),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_text(json.dumps(report,indent=2)+'\n',encoding='utf8',newline='\n')
    print('Current Rp native constructor constant functions projected',flush=True)
    return owner


if __name__=='__main__':run()
