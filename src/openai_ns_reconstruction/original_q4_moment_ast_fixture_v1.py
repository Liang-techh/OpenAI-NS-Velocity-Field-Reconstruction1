"""Fresh synthetic inputs for source and selected-row witnesses."""
import copy

GRAPH = '''
class FunctionTransportGraph:
 def node(self, operation, **attributes):
  row = dict(operation=operation, **attributes)
  key = json.dumps(row, sort_keys=True, separators=(',', ':'))
  if key not in self.keys:
   self.keys[key] = len(self.nodes)
   self.nodes.append(row)
  return FunctionRef(self, self.keys[key])
 def symbol(self, name):
  return self.node('bound_variable', name=name)
'''
CONTROLS = '''
def integral(g, integrand, variable, lower, upper, **metadata):
 g.ids([integrand, lower, upper])
 return g.node('definite_integral', integrand=integrand.node, variable=variable,
  lower=lower.node, upper=upper.node, exact_function_integral=True,
  numerical_value_not_installed=True, **metadata)
'''
CORE = '''
class JetAlgebra:
 def integral(self,jet,variable,lower,upper,**attrs):
  g=self.g
  return FunctionJet(tuple(controls.integral(g,q,g.nodes[variable.node]['name'],lower,upper,
   endpoint_independent_of_Z=True,ordinary_slow_Z_derivative_order=k,
   true_core_inlet_complete_history=True,numerical_cap_not_defining_integrand=True,**attrs)
   for k,q in enumerate(rows(self.trim(jet)))))
class CurrentOriginalCoreInletAxial6GenericSource:
 def __init__(self):
  t=g.symbol('original_core_inlet_moment_t');sample=g.mul(rho,t)
  phi=a.at(self.Phi,rho,sample);u=a.at(self.U,rho,sample)
  densities=dict(H=a.scale(phi,g.mul(c(2),t)),M=u,K=a.scale(a.mul(phi,u),g.mul(c(2),t)),
   A=a.mul(u,u),B=a.scale(a.mul(phi,phi),t),C=a.mul(phi,phi))
  self.moment_functions={key:a.integral(body,t,g.zero,g.one,actual_core_normalized_moment=key,
   same_unique_nonlinear_core_and_full_density_product=True) for key,body in densities.items()}
'''

def run_fixture(check_source):
    source = {'graph': GRAPH, 'controls': CONTROLS, 'core_inlet': CORE}
    base = check_source(source)
    checks = {'valid_accepted': base.get('status') == 'verified',
              'identity_not_claimed': base.get('scope', {}).get('global_binder_identity_proved') is False}
    for label, key, old, new in (
        ('wrong_symbol', 'graph', "'bound_variable'", "'another'"),
        ('wrong_interning', 'graph', 'self.keys[key] = len(self.nodes)', 'self.keys[key] = 1'),
        ('wrong_transport', 'core_inlet', "g.nodes[variable.node]['name']", 'variable.node'),
        ('wrong_variable', 'controls', 'variable=variable,', 'variable=integrand,'),
        ('wrong_density', 'core_inlet', 'M=u,K=', 'M=phi,K='),
        ('rebound_t', 'core_inlet', "sample=g.mul(rho,t)", "sample=g.mul(rho,t);t=g.zero"),
        ('bad_syntax', 'graph', 'class FunctionTransportGraph:', 'class :')):
        case = dict(source); case[key] = case[key].replace(old, new)
        checks[label + '_rejected'] = check_source(case).get('status') != 'verified'
    return {'passed': all(checks.values()), 'checks': checks}

def run_graph_fixture(bind):
    name = 'original_core_inlet_moment_t'
    nodes = [{'operation': 'bound_variable', 'name': name},
             {'operation': 'density'},
             {'operation': 'definite_integral', 'variable': name, 'integrand': 1}]
    q4 = {'exact_graph_nodes': copy.deepcopy(nodes), 'source_family': {'id': 'synthetic'}}
    core = copy.deepcopy(q4)
    inspection = {'inspection': {'status': 'observed',
                  'evidence': {'integral_node': 2, 'integrand_node': 1, 'binder_name': name,
                               'unshadowed_binder_nodes': [0], 'shadowed_binder_nodes': [],
                               'same_name_node_conflict': False, 'integral_source_row': copy.deepcopy(nodes[2])},
                  'scope': {'integral_binder_identity_proved': False, 'integrand_scan_complete': False}}}
    inputs = (q4, core, [inspection], {'status': 'verified'})
    base = bind(*inputs)
    checks = {'valid_accepted': base.get('status') == 'matched',
              'global_identity_not_claimed': base.get('scope', {}).get('global_string_to_node_identity_proved') is False,
              'capture_safety_not_claimed': base.get('scope', {}).get('capture_safe_substitution_proved') is False}
    for label in ('family', 'row', 'name', 'shadow', 'source'):
        a, b, ins, src = copy.deepcopy(inputs)
        if label == 'family': b['source_family'] = {'id': 'different'}
        if label == 'row': b['exact_graph_nodes'][1]['operation'] = 'other'
        if label == 'name': ins[0]['inspection']['evidence']['binder_name'] = 'different'
        if label == 'shadow': ins[0]['inspection']['evidence']['shadowed_binder_nodes'] = [0]
        if label == 'source': src['status'] = 'mismatch'
        checks[label + '_rejected'] = bind(a, b, ins, src).get('status') != 'matched'
    return {'passed': all(checks.values()), 'checks': checks}
