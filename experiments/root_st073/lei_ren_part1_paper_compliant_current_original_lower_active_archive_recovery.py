"""Pinned source-only recovery for the actual lower/active function graph.

No nonlinear oracle, old runtime cache, numerical center, or ancestor run is
created. Restored views expose the exact saved functions required by the
join proof; numerical oracle services remain unavailable in these views.
"""
import copy
import importlib
import json
from pathlib import Path
import sys
from types import SimpleNamespace as NS

sys.setrecursionlimit(max(10000,sys.getrecursionlimit()))
import lei_ren_part1_paper_compliant_current_original_lower_active_common_source_join_functions as current


def restore(require_checked=True):
    global current
    core_source=importlib.import_module(current.PREFIX+'current_original_n1_core_fixed_point_functions')
    importlib.reload(core_source.c3)
    # Refresh in dependency order before creating any typed function view.
    # Reloading a jet class after hydration would invalidate its type contract.
    modules={}
    for stem in ('current_original_n1_core_fixed_point_functions',
        'current_original_n1_lower_collar_functions',
        'current_original_n1_analytic_core_lower_pressure_C3_functions',
        'current_original_core_inlet_axial6_generic_source',
        'current_original_prescribed_bridge_C6_functions',
        'current_original_prescribed_loop_C5_functions',
        'current_original_modified_Q4_pressure_C3_functions'):
        modules[stem]=importlib.reload(importlib.import_module(current.PREFIX+stem))
    current=importlib.reload(current)
    raw=current.parent.linear.read(current.NAME);g=current.source.FunctionTransportGraph()
    assert raw['input_hashes'][Path(current.__file__).name]==current.sha(Path(current.__file__).name)
    assert raw['input_hashes'][current.RECOVERY]==current.sha(current.RECOVERY)
    assert not any(raw[k] for k in current.GATES+current.OPEN)
    g.nodes=raw['exact_graph_nodes']
    g.keys={json.dumps(row,sort_keys=True,separators=(',',':')):i for i,row in enumerate(g.nodes)}
    frozen=copy.deepcopy(g.nodes);family=copy.deepcopy(raw['source_family']);frozen_family=copy.deepcopy(family)
    ref=lambda n:current.source.FunctionRef(g,n)
    def jet(ids):
        assert isinstance(ids,list) and all(type(n) is int and 0<=n<len(g.nodes) for n in ids)
        return current.FunctionJet(tuple(ref(n) for n in ids))
    jets=lambda packet:{k:jet(v) for k,v in packet.items()}
    arrays=lambda packet:{k:[jet(q) for q in v] for k,v in packet.items()}
    registry=[]
    def frame(**values):
        obj=NS(g=g,family=family,numerical_oracle_available=False,**values)
        registry.append((obj,tuple(values),current.encode(values)));return obj
    def load(stem,key):
        module=modules[stem];data=current.parent.linear.read(module.NAME)
        rec=current.parent.linear.read(module.RECEIPT)
        gates=getattr(module,'GATES',(getattr(module,'GATE',None),))
        assert rec['all_passed'] and all(rec[k] for k in gates) and not any(rec[k] for k in module.OPEN)
        assert data['source_family']==family and rec['source_family']==family
        assert data['exact_graph_nodes']==g.nodes[:len(data['exact_graph_nodes'])]
        pins=dict(rec['input_hashes']);archive=pins.pop(module.NAME)
        checker=Path(module.__file__).stem+'_check.py';checker_hash=pins.pop(checker)
        assert pins==data['input_hashes']
        assert archive==current.sha(module.NAME) and checker_hash==current.sha(checker)
        assert pins[Path(module.__file__).name]==current.sha(Path(module.__file__).name)
        for name,digest in rec['input_hashes'].items():assert current.sha(name)==digest,name
        return module,data[key]
    coremod,cb=load('current_original_n1_core_fixed_point_functions','core_fixed_point_bindings')
    core=frame(source_module=coremod.__name__,functions=jets(cb['functions']),
        parameters={k:ref(v) for k,v in cb['parameters'].items()},components=jets(cb['components']),
        definition={**cb['definition'],'source_variables':{k:ref(v) for k,v in cb['definition']['source_variables'].items()}},
        terms={k:[ref(n) for n in v] for k,v in cb['definition']['actual_source_terms'].items()},
        rho=ref(cb['definition']['source_variables']['s']))
    lm,lb=load('current_original_n1_lower_collar_functions','lower_collar_bindings')
    lf={}
    for key,value in lb['fields'].items():
        if key.endswith('_rho_rows'):lf[key]=[jet(q) for q in value]
        elif isinstance(value,list):lf[key]=jet(value)
        elif key=='domain':lf[key]={k:ref(v) for k,v in value.items()}
        elif type(value) is int:lf[key]=ref(value)
        else:lf[key]=value
    lower0=frame(fields=lf,parameters={k:ref(v) for k,v in lb['parameters'].items()},
        core_radial_rows=arrays(lb['core_radial_rows']),sc_bits=lb['geometry_certificate']['sc_exact_mpf_tuple'])
    lower0.core=core
    hm,hb=load('current_original_n1_analytic_core_lower_pressure_C3_functions','core_lower_C5_C3_bindings')
    lower=frame(core_functions=jets(hb['core_functions']),core_radial_rows=arrays(hb['core_radial_rows']),
        fields=jets(hb['fields']),radial_fields=arrays(hb['radial_fields']),pressure=jets(hb['pressure']))
    lower.core=core;lower.lower=lower0
    im,ib=load('current_original_core_inlet_axial6_generic_source','core_inlet_axial6_bindings')
    inlet=frame(core_functions=jets(ib['core_functions']),P0=jet(ib['P0']),
        moment_functions=jets(ib['moments']),inlet_moments=jets(ib['inlet_moments']),
        fields=jets(ib['fields']),Pstar=ref(ib['Pstar']),radius=ref(ib['radius']))
    inlet.previous=lower;inlet.core=core
    bm,bb=load('current_original_prescribed_bridge_C6_functions','prescribed_bridge_C6_bindings')
    def comparison(packet):
        return {**packet,'phase':ref(packet['phase']),'phi':jet(packet['phi']),'U':jet(packet['U']),
            'logphi':jet(packet['logphi']),'moments':jets(packet['moments']),'directions':jets(packet['directions'])}
    def background(packet):
        return {**packet,'phase':ref(packet['phase']),'radius':ref(packet['radius']),
            'physical_measure':ref(packet['physical_measure']),'chi':ref(packet['chi']),
            **{k:jets(packet[k]) for k in ('actual','moments','fields','radial','generic')}}
    bridge=frame(core7=jets(bb['core7']),radial7=arrays(bb['radial7']),initial7=jets(bb['initial7']),
        P0=jet(bb['P0']),hb=ref(bb['hb']),Ra=ref(bb['Ra']),
        comparison_micro=comparison(bb['comparison_micro']),micro=background(bb['micro']))
    bridge.core=core;bridge.inlet=inlet
    native=coremod.c3.C3Function(*current.rows(bridge.P0)[:4])
    core.pulse=NS(own=NS(jets=NS(native={'P0':native})))
    bridge.derivative=hm.HigherCoreZDerivative(lower0)
    for group in (bridge.core7,bridge.initial7):
        for q in group.values():
            # Incoming eighth rows are checked by differentiation below,
            # so do not seed their sixth-to-seventh link here.
            incoming=current.rows(q)[:7]
            for left,right in zip(incoming,incoming[1:]):bridge.derivative.seed(left,right)
    for group in bridge.radial7.values():
        for q in group:
            incoming=current.rows(q)
            for left,right in zip(incoming,incoming[1:]):bridge.derivative.seed(left,right)
    pm,pb=load('current_original_prescribed_loop_C5_functions','prescribed_loop_C5_bindings')
    charts={}
    for name,packet in pb['charts'].items():
        q={**packet}
        for key in ('coordinate','radius','physical_measure','actual_phase'):q[key]=ref(packet[key])
        q['domain']=[ref(n) for n in packet['domain']]
        for key in ('roots','incoming','history','contributions','modified_fields'):q[key]=jets(packet[key])
        q['memory']={k:ref(v) for k,v in packet['memory'].items()}
        loop_packet=packet['loop']
        q['loop']={**loop_packet,**{k:jet(loop_packet[k]) for k in ('q','A','B_over_Pstar','A_phi','B_over_Pstar_phi')},
            'flat_predicate':ref(loop_packet['flat_predicate']),
            'original_collar_upper':ref(loop_packet['original_collar_upper']) if loop_packet['original_collar_upper'] is not None else None}
        finite=packet['finite_N'];q['finite_N']={**finite,
            **{k:jet(finite[k]) for k in ('E_N','V_N','delta_E','delta_V','A_over_N')},
            'five_signed_density_increments':jets(finite['five_signed_density_increments'])}
        q['source_chart']=bridge.micro if name=='first_micro' else None
        charts[name]=q
    loop=frame(charts=charts,N=ref(pb['N']),sc=ref(pb['sc']),zero_inlet=charts['first_micro']['incoming'])
    loop.source_charts={'first_micro':bridge.micro}
    mm,mb=load('current_original_modified_Q4_pressure_C3_functions','modified_Q4_pressure_C3_bindings')
    modified_charts={}
    for name,packet in mb['charts'].items():
        q={**packet,'source_chart':charts[name]}
        q['physical']=jets(packet['physical'])
        q['radial']={k:jet(v) if isinstance(v,list) else v for k,v in packet['radial'].items()}
        modified_charts[name]=q
    modified=frame(charts=modified_charts);modified.loop=loop
    modified.background=NS(bridge=bridge)
    sb=raw['lower_active_common_source_bindings'];owner=current.CurrentOriginalLowerActiveCommonSourceJoinFunctions.__new__(
        current.CurrentOriginalLowerActiveCommonSourceJoinFunctions)
    owner.g=g;owner.family=family;owner.hashes=dict(raw['input_hashes']);owner.lower=lower;owner.core=core
    owner.modified=modified;owner.loop=loop;owner.bridge=bridge;owner.inlet=inlet
    owner.lower_phase=lower0.fields['query'];owner.active_phase=charts['first_micro']['coordinate']
    owner.left_endpoint=lower0.parameters['sc_half'];owner.right_endpoint=charts['first_micro']['domain'][0]
    for key in ('left_fields','right_fields','left_trace','right_trace','lower_histories'):
        setattr(owner,key,jets(sb[key]))
    owner.lower_history_mean=jet(sb['lower_history_mean'])
    owner.lower_history_pressure=jet(sb['lower_history_pressure'])
    owner.inlet_binding=dict(core_radius=core.rho,actual_radius=inlet.radius,
        six_core_moment_functions=inlet.moment_functions,six_rho4_inlet_functions=inlet.inlet_moments,
        comparison_incoming=bridge.initial7,lower_history_incoming=inlet.inlet_moments,
        normalized_P0=bridge.P0,physical_P0=core.functions['P0'],Pstar=inlet.Pstar,
        physical_axis_scale=ref(sb['inlet_binding']['physical_axis_scale']))
    owner.control_binding=dict(comparison=bridge.comparison_micro,actual=bridge.micro,
        lower_ratio=lower.fields['ratio'],lower_phi=lower.radial_fields['Phi_rho_rows'][0],
        lower_psi=lower.radial_fields['Psi_rho_rows'][0],same_hb=bridge.hb,
        actual_chi=lower0.fields['chi'],active_loop={k:charts['first_micro']['loop'][k]
            for k in ('q','A','B_over_Pstar','A_phi','B_over_Pstar_phi','flat_predicate','original_collar_upper')},
        active_finite_N=charts['first_micro']['finite_N'])
    owner.left_radius=ref(sb['left_radius']);owner.right_radius=ref(sb['right_radius'])
    owner.sc_binding=copy.deepcopy(sb['sc_binding']);owner.integrated_theorem=copy.deepcopy(sb['integrated_theorem'])
    owner.primitive_theorem=copy.deepcopy(sb['primitive_theorem']);owner.definition=copy.deepcopy(sb['definition'])
    owner.snapshot=copy.deepcopy(sb);owner.acceptance_loaded=False
    source_coordinates=(owner.lower_phase,owner.active_phase,owner.left_endpoint,owner.right_endpoint)
    owner_hashes=dict(owner.hashes)
    seeded_derivatives=dict(bridge.derivative.derivatives)
    allowed=set(bridge.derivative.allowed)
    def assert_saved():
        assert g.nodes[:len(frozen)]==frozen
        assert owner.g is g and owner.family==frozen_family and owner.hashes==owner_hashes
        for name,digest in owner.hashes.items():assert current.sha(name)==digest,name
        for obj,keys,expected in registry:
            assert obj.g is g and obj.family==frozen_family and obj.numerical_oracle_available is False
            assert current.encode({k:getattr(obj,k) for k in keys})==expected
        assert owner._snapshot()==sb and owner.snapshot==sb
        assert (owner.lower_phase,owner.active_phase,owner.left_endpoint,owner.right_endpoint)==source_coordinates
        assert owner.lower is lower and owner.core is core and owner.modified is modified
        assert owner.loop is loop and owner.bridge is bridge and owner.inlet is inlet
        assert lower.core is core and lower.lower is lower0 and lower0.core is core
        assert inlet.previous is lower and inlet.core is core and bridge.core is core and bridge.inlet is inlet
        assert modified.loop is loop and modified.background.bridge is bridge
        assert loop.source_charts['first_micro'] is bridge.micro
        assert charts['first_micro']['source_chart'] is bridge.micro
        assert modified.charts['first_micro']['source_chart'] is charts['first_micro']
        assert current.rows(core.pulse.own.jets.native['P0'])==current.rows(native)
        assert bridge.derivative.g is g and bridge.derivative.core is core
        assert bridge.derivative.allowed==bridge.derivative._allowed==allowed
        assert all(bridge.derivative.derivatives[k]==v for k,v in seeded_derivatives.items())
        assert owner.numerical_oracle_available is False and owner.archive_only_source_recovery is True
    owner.assert_graph=assert_saved
    if require_checked:
        rec=current.parent.linear.read(current.RECEIPT)
        assert rec['all_passed'] and all(rec[k] for k in current.GATES) and not any(rec[k] for k in current.OPEN)
        assert rec['source_family']==family
        pins=dict(rec['input_hashes']);archive=pins.pop(current.NAME)
        checker=Path(current.__file__).stem+'_check.py';checker_hash=pins.pop(checker)
        assert pins==owner.hashes and archive==current.sha(current.NAME) and checker_hash==current.sha(checker)
        assert len(rec['checks'])==10 and all(c['passed'] for c in rec['checks'].values())
        assert len(rec['checks']['trace_transfer']['actual_registered_two_sided_trace_rows'])==44
        for name,digest in rec['input_hashes'].items():assert current.sha(name)==digest,name
        assert rec['input_hashes'][current.NAME]==current.sha(current.NAME)
        owner.acceptance_loaded=True
    owner.numerical_oracle_available=False;owner.archive_only_source_recovery=True;owner.assert_graph()
    return owner
