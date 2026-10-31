"""Bounded source witness and selected archive-row comparison; no evaluation."""
import ast

_GRAPH = '''
class G:
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
_CONTROLS = '''
def integral(g, integrand, variable, lower, upper, **metadata):
 g.ids([integrand, lower, upper])
 return g.node('definite_integral', integrand=integrand.node, variable=variable,
  lower=lower.node, upper=upper.node, exact_function_integral=True,
  numerical_value_not_installed=True, **metadata)
'''
_JET = '''
class J:
 def integral(self,jet,variable,lower,upper,**attrs):
  g=self.g
  return FunctionJet(tuple(controls.integral(g,q,g.nodes[variable.node]['name'],lower,upper,
   endpoint_independent_of_Z=True,ordinary_slow_Z_derivative_order=k,
   true_core_inlet_complete_history=True,numerical_cap_not_defining_integrand=True,**attrs)
   for k,q in enumerate(rows(self.trim(jet)))))
'''
_INIT = '''
t=g.symbol('original_core_inlet_moment_t');sample=g.mul(rho,t)
phi=a.at(self.Phi,rho,sample);u=a.at(self.U,rho,sample)
densities=dict(H=a.scale(phi,g.mul(c(2),t)),M=u,K=a.scale(a.mul(phi,u),g.mul(c(2),t)),
 A=a.mul(u,u),B=a.scale(a.mul(phi,phi),t),C=a.mul(phi,phi))
self.moment_functions={key:a.integral(body,t,g.zero,g.one,actual_core_normalized_moment=key,
 same_unique_nonlinear_core_and_full_density_product=True) for key,body in densities.items()}
'''

def _function(tree, name, class_name=None):
    body = tree.body
    if class_name is not None:
        classes = [n for n in body if isinstance(n, ast.ClassDef) and n.name == class_name]
        if len(classes) != 1:
            raise ValueError('class missing or ambiguous: ' + class_name)
        body = classes[0].body
    found = [n for n in body if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef)) and n.name == name]
    if len(found) != 1:
        raise ValueError('function missing or ambiguous: ' + name)
    return found[0]

def _dump(node):
    return ast.dump(node, include_attributes=False)

def check_source(source_texts):
    scope = {'global_binder_identity_proved': False, 'numeric_field_value_computed': False,
             'old_errors_removed': False, 'archive_constructor_lineage_proved': False}
    checks = {}
    try:
        graph = ast.parse(source_texts['graph'])
        controls = ast.parse(source_texts['controls'])
        core = ast.parse(source_texts['core_inlet'])
        expected_graph = ast.parse(_GRAPH)
        for name in ('node', 'symbol'):
            checks['graph_' + name] = _dump(_function(graph, name, 'FunctionTransportGraph')) == _dump(_function(expected_graph, name, 'G'))
        checks['controls_integral'] = _dump(_function(controls, 'integral')) == _dump(_function(ast.parse(_CONTROLS), 'integral'))
        checks['jet_integral'] = _dump(_function(core, 'integral', 'JetAlgebra')) == _dump(_function(ast.parse(_JET), 'integral', 'J'))
        init = _function(core, '__init__', 'CurrentOriginalCoreInletAxial6GenericSource')
        actual = [_dump(n) for n in init.body]
        wanted = [_dump(n) for n in ast.parse(_INIT).body]
        positions = [[i for i, x in enumerate(actual) if x == item] for item in wanted]
        checks['constructor_statements_unique_and_ordered'] = all(len(p) == 1 for p in positions) and [p[0] for p in positions] == sorted(p[0] for p in positions)
        checks['single_t_store_in_constructor'] = sum(isinstance(n, ast.Name) and n.id == 't' and isinstance(n.ctx, ast.Store) for n in ast.walk(init)) == 1
    except (KeyError, TypeError, SyntaxError, ValueError) as exc:
        return {'status': 'unresolved', 'checks': checks, 'error': str(exc), 'scope': scope}
    return {'status': 'verified' if all(checks.values()) else 'mismatch', 'checks': checks, 'scope': scope}

def bind(q4, core, inspections, source_result):
    scope = {'global_string_to_node_identity_proved': False, 'capture_safe_substitution_proved': False,
             'numeric_field_value_computed': False, 'old_errors_removed': False,
             'archive_constructor_lineage_proved': False}
    checks = {'source_witness_verified': source_result.get('status') == 'verified'}
    evidence = []
    try:
        qnodes, cnodes = q4['exact_graph_nodes'], core['exact_graph_nodes']
        checks['same_source_family'] = q4['source_family'] == core['source_family']
        candidates = set()
        for item in inspections:
            ins = item['inspection']; ev = ins['evidence']; prior_scope = ins['scope']
            idx, body = ev['integral_node'], ev['integrand_node']
            name = ev['binder_name']; unshadowed = ev['unshadowed_binder_nodes']
            candidate = unshadowed[0] if len(unshadowed) == 1 else None
            if candidate is not None: candidates.add(candidate)
            ok = (ins['status'] == 'observed' and candidate is not None and
                  ev['shadowed_binder_nodes'] == [] and ev['same_name_node_conflict'] is False and
                  prior_scope['integral_binder_identity_proved'] is False and
                  name == 'original_core_inlet_moment_t' and
                  qnodes[idx] == cnodes[idx] == ev['integral_source_row'] and
                  qnodes[body] == cnodes[body] and qnodes[candidate] == cnodes[candidate] and
                  qnodes[candidate] == {'operation': 'bound_variable', 'name': name} and
                  qnodes[idx].get('operation') == 'definite_integral' and
                  qnodes[idx].get('variable') == name and qnodes[idx].get('integrand') == body)
            evidence.append({'integral_node': idx, 'integrand_node': body, 'candidate_node': candidate,
                             'selected_rows_match': bool(ok), 'integrand_scan_complete': prior_scope.get('integrand_scan_complete')})
        checks['nonempty_inspections'] = bool(inspections)
        checks['all_selected_rows_match'] = bool(evidence) and all(x['selected_rows_match'] for x in evidence)
        checks['one_candidate_in_inspected_portions'] = len(candidates) == 1
    except (KeyError, IndexError, TypeError, ValueError) as exc:
        return {'status': 'unresolved', 'checks': checks, 'evidence': evidence, 'error': str(exc), 'scope': scope}
    return {'status': 'matched' if all(checks.values()) else 'mismatch', 'checks': checks, 'evidence': evidence, 'scope': scope}
