"""Pinned actual C5 finite-N signed m15 expressions, before numerical boxing.

Only saved source graphs and ASTs are read. The local graph defines the
compact and Rm-patch pressure *increment* differences. It does not replace
their inherited incoming pressures, select cap values, or execute a source.
"""
from __future__ import annotations

import ast
import gzip
import hashlib
import json
from pathlib import Path

from .original_finite_N_density_kernel import ordinary_jet, physical_fields, volume_density_difference, volume_difference_terms
from .original_signed_density_expression import SignedDensityExpressionGraph

PREFIX='lei_ren_part1_paper_compliant_current_original_'
PARENT=PREFIX+'quiet_scalar_m15_numeric.json'
PARENT_SHA256='93f419a4a10809ff15341349e86856a1a191937cc27612891be67842148375f9'
COMPACT=PREFIX+'true_compact_C5_functions.json.gz'
RM=PREFIX+'true_Rm_patch_finite_N_C5_functions.json.gz'
RM_BASE=PREFIX+'true_Rm_patch_C6_functions.json.gz'


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def digest(value):return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def load(path):
    raw=Path(path).read_bytes()
    return json.loads(gzip.decompress(raw) if str(path).endswith('.gz') else raw)
def require(value,message):
    if not value:raise ValueError(message)


def ast_guard(root,relative,name,text):
    source=(root/relative).read_text(encoding='utf-8');tree=ast.parse(source)
    found=[n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name==name]
    require(len(found)==1,'Unique source function required: '+name)
    fn=found[0];body=''.join((ast.get_source_segment(source,fn) or '').split())
    checks={q:''.join(q.split()) in body for q in text}
    require(all(checks.values()),'Changed density/FTC source: '+name)
    return dict(passed=True,path=relative,function=name,required_text=checks,
        AST_sha256=hashlib.sha256(ast.dump(fn,include_attributes=False).encode()).hexdigest(),source_executed=False)


def build_signed_sources(root,parent):
    # Imports are bounded arithmetic only; none of the archived owners run.
    from .original_finite_N_density_kernel import endpoint_density_difference
    root=Path(root);directory=root/'experiments/root_st073'
    require(sha(directory/PARENT)==PARENT_SHA256,'Accepted quiet parent required')
    inputs=dict(parent['input_hashes']);inputs[PARENT]=PARENT_SHA256
    reports={}
    for name in (COMPACT,RM,RM_BASE):
        require(inputs.get(name)==sha(directory/name),'Unpinned actual C5 source: '+name)
        reports[name]=load(directory/name)
        require(all(reports[name]['source_family'].get(k)==v for k,v in parent['source_family'].items()),'Different source family')
    base_nodes=reports[RM_BASE]['exact_graph_nodes'];rm_nodes=reports[RM]['exact_graph_nodes']
    require(rm_nodes[:len(base_nodes)]==base_nodes,'Rm baseline is not the identical source prefix')
    cb=reports[COMPACT]['true_compact_C5_bindings'];rb=reports[RM]['Rm_patch_finite_N_C5_bindings']
    bg=reports[RM_BASE]['true_Rm_patch_C6_bindings']['chart']
    implementation=dict(parent['implementation_hashes'])
    macro='src/openai_ns_reconstruction/original_macro_m15.py'
    phys='experiments/root_st073/'+PREFIX+'modified_Q4_pressure_C3_functions.py'
    outer='experiments/root_st073/'+PREFIX+'true_outer_C6_C5_functions.py'
    rsource='experiments/root_st073/'+PREFIX+'true_Rm_patch_finite_N_C5_functions.py'
    for rel in (macro,phys,outer,rsource):
        expected=implementation.get(rel) if rel.startswith('src/') else inputs.get(Path(rel).name)
        require(expected==sha(root/rel),'Pinned density source changed: '+rel)
    guards={
        'half_S1_coefficients':ast_guard(root,macro,'half_S1_coefficients',('A, D = 1-z*z*delta, 1-z*z',
            'dict(q=one/(2*A), q2=one/4, U2=-z*z/(A*A))','UqZ=D/(2*A)','UUZ=z*D/(A*A)')),
        'physical_recovery':ast_guard(root,phys,'physical_pressure_fields',('U=a.scale(fields[\'V\'],Pstar)',
            'my=a.sub(U,m)','Q=div(numerator(U,m));Qy=div(numerator(Uy,my))')),
        'normalized_baseline_Q':ast_guard(root,phys,'generic_functions',(
            "E,V,m,h,k,e=(a.trim(fields[name]) for name in ('E','V','m','h','k','e'))",
            'Q=div(a.sub(a.scale(a.mul(z,V),2),transport))')),
        'reduced_S':ast_guard(root,phys,'reduced_pressure_S',('linear=div(', 'quadratic=a.mul(Q,', 'axial=div(a.mul(U,')),
        'compact_FTC':ast_guard(root,outer,'pressure_geometry',("('S1',2)","m15=a.sub(a.scale(cumulative['S1'],'1/2'),a.sub(boundary['R_Q_y'],initial['R_Q_y']))")),
        'Rm_FTC':ast_guard(root,rsource,'pressure_continuation',("('S1',2,g.mul(Rm,Rm,x))",'nonzero_preceding_S0_S1_retained=True',
            "m15=a.sub(a.scale(cumulative['S1'],'1/2'),a.sub(boundary['R_Q_y'],initial['R_Q_y']))")),
    }
    specs={
        'compact':dict(report=COMPACT,b=cb,background=cb['background'],x=260181,R=260862,
            delta=25220,Pstar=33572,N=18,Z=20,domain=[1,48]),
        'Rm_patch':dict(report=RM,b=rb,background=bg,x=185676,R=187214,
            delta=25220,Pstar=33572,N=18,Z=20,domain=[1,2948]),
    }
    output={}
    for label,s in specs.items():
        name=s['report'];raw=reports[name];nodes=raw['exact_graph_nodes'];b=s['b'];background=s['background']
        graph=SignedDensityExpressionGraph()
        ref=lambda n:graph.reference(name,inputs[name],n,source_nodes=nodes)
        jet=lambda rows:ordinary_jet(*(ref(n) for n in rows[:3]))
        pair=lambda value:[value.value.node,value.dz.node]
        z=ordinary_jet(ref(s['Z']),graph.one,graph.zero);de=ref(s['delta']);Pstar=ref(s['Pstar'])
        require(nodes[s['Z']].get('operation')=='bound_variable' and nodes[s['Z']].get('name')=='Z','Wrong ordinary slow variable')
        require(b['N']==18 and nodes[18].get('definition')=='N=2^(2^J)' and
            nodes[18].get('one_common_N_for_all_charts_and_control_denominators') is True,'Wrong common exact N')
        require(background['domain']==s['domain'] and background['radius']==s['R'],'Changed chart geometry')
        # generic Q is normalized by Pstar; the same physical recovery applies
        # Pstar once. It retains the entire actual signed background history.
        U0=jet(background['fields']['V'])*Pstar
        Q0=jet(background['generic']['Q'])*Pstar
        Uy0=jet(background['radial']['V_y'])*Pstar
        baseline=physical_fields(U0,Q0)
        full=physical_fields(jet(b['physical']['U']),jet(b['physical']['Q']))
        correction=physical_fields(full.U-U0,full.Q-Q0)
        L=1-z*z*de;D=1-z*z
        Qy0=(z*(1+de)*U0+2*z*Uy0-D*baseline.Uz)/L-Q0
        Qy=jet(b['physical']['Q_y'])
        delta_Qy=Qy-Qy0
        expected={'baseline':baseline,'correction':correction,'Z':z,'delta':de}
        guard=lambda role,value: expected.get(role) is value
        dv=volume_density_difference(baseline,correction,z,de,source_provenance_guard=guard)
        dF=endpoint_density_difference(baseline,correction,z,de,source_provenance_guard=guard)
        terms=volume_difference_terms(baseline,correction,z,de,source_provenance_guard=guard)
        R=ref(s['R']);x=ref(s['x']);one=ref(s['domain'][0])
        delta_boundary=dF*(R*R)-delta_Qy*R
        increment=[];volumes=[]
        for k in (0,1):
            density=R*R*(dv.value if k==0 else dv.dz)/x
            integral=graph.insert(dict(operation='definite_integral',integrand=density.node,variable=x.node,
                lower=one.node,upper=x.node,exact_function_integral=True,numerical_value_not_installed=True,
                physical_measure='R^2 dx/x',ordinary_Z_derivative_order=k))
            B=delta_boundary.value if k==0 else delta_boundary.dz
            inlet=graph.insert(dict(operation='function_substitution',expression=B.node,variable=x.node,value=one.node,
                endpoint_independent_of_Z=True))
            increment.append((B-inlet+integral).node);volumes.append(integral.node)
        output[label]=dict(report=name,report_sha256=inputs[name],source_graph_node_count=len(nodes),
            source_scope='whole chart, not an endpoint-only substitution',native_variable=s['x'],domain=s['domain'],
            physical_radius=s['R'],ordinary_slow_Z_variable=s['Z'],common_exact_N=s['N'],delta=s['delta'],Pstar=s['Pstar'],
            source_ports=dict(full_physical=b['physical'],baseline_V=background['fields']['V'],
                baseline_normalized_Q=background['generic']['Q'],baseline_V_y=background['radial']['V_y'],
                baseline_m=background['fields']['m']),
            signed_local_roots=dict(delta_U_C0_C1=pair(correction.U),delta_Q_C0_C1=pair(correction.Q),
                delta_v_C0_C1=pair(dv),delta_F_C0_C1=pair(dF),delta_Q_y_C0_C1=pair(delta_Qy),
                delta_boundary_C0_C1=pair(delta_boundary),delta_volume_C0_C1=volumes,
                delta_m15_chart_increment_C0_C1=increment),
            factored_six_volume_terms={k:pair(v) for k,v in terms.items()},**graph.packet(),
            baseline_recovery='U0=Pstar*V0, Q0=Pstar*generic.Q0; Qy0=[Z(1+delta)U0+2Z*Uy0-D*Uz0]/L-Q0',
            FTC='DeltaM15(x)-DeltaM15(1)=[R^2*DeltaF-R*DeltaQy](x)-[R^2*DeltaF-R*DeltaQy](1)+integral_1^x R^2*Delta v dx/x',
            incoming_pressure_difference_not_zeroed=True,incoming_pressure_difference_numerically_resolved=False,
            exact_N_materialized=False,extra_N_division=False,numerical_signed_values_installed=False,
            caps_removed_or_narrowed=False,terminal_zero_claim=False)
    return dict(passed=True,source_only=True,source_executed=False,AST=guards,charts=output,
        Rm_baseline_identical_graph_prefix=dict(passed=True,nodes=len(base_nodes),report=RM_BASE,sha256=inputs[RM_BASE]),
        original_function_definitions_bound=True,numerical_signed_values_installed=False,
        finite_N_caps_active=True,global_terminal_closure=False,accepted_ancestors_replayed=False)
