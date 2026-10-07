"""Explicit accepted-owner query for the exact formal generic left inlet.

No constructor or saved-cover fallback. The returned objects are genuine
coordinate-query source covers, not selected point function values.
"""
import json
import lei_ren_part1_paper_compliant_current_generic_five_defect_bounds as current

HERE,PREFIX,sha=current.HERE,current.PREFIX,current.sha
packets=current.packets


def original_left_inlet(owner,*,Z=(-1,1)):
    if owner is None:raise ValueError('An existing checked CurrentBridgeBackgroundTensor owner is required; no constructor or cache fallback')
    if type(owner).__name__!='CurrentBridgeBackgroundTensor' or type(owner).__module__!=PREFIX+'current_bridge_background_tensor':
        raise ValueError('Exact existing original bridge tensor owner required')
    if not owner.acceptance_loaded:raise ValueError('The original bridge owner must already have checked acceptance')
    service=packets.CurrentSourcePackets(owners={'bridge':owner})
    receipt=json.loads((HERE/current.RECEIPT).read_bytes())
    if not receipt['all_passed'] or not receipt[current.GATE] or receipt['source_family']!=service.family:
        raise ValueError('Same checked current generic source/domain required')
    service.bind_hashes(receipt['input_hashes']);service.bind_hashes({current.RECEIPT:sha(current.RECEIPT)})
    c=owner.ctx;read=lambda row:packets.interval(c,row);ep=packets.recovery.endpoints
    name=PREFIX+'current_inner_exit_strict_collar.json';left=json.loads((HERE/name).read_bytes())
    proof=left['explicit_current_inner_exit_strict_collar'];sc=read(proof['selected_first_phase_endpoint'])
    phase=sc/2
    if ep(phase)[0]<=0 or ep(phase)[1]>=1:raise ArithmeticError('Same original interior first-bridge phase required')
    domain=json.loads((HERE/(PREFIX+'current_generic_shear_loop_domain.json')).read_bytes())['current_original_generic_loop_domain']
    uniform=json.loads((HERE/(PREFIX+'current_generic_shear_uniform_inputs.json')).read_bytes())['current_actual_logarithmic_loop_scales']
    eta=read(uniform['selected_positive_eta_log'])
    limit=read(domain['allowed_eta_log_upper_for_automatic_constant_edges'])
    if ep(eta)[1]>ep(limit)[0]:raise ValueError('Current actual eta must satisfy the checked left/right flat-edge constraint')
    packet=service.query('bridge_first',Z,phase)
    state=packets.FactoredRecoveryState(packet)
    service.bind_hashes({__file__.replace('\\','/').split('/')[-1]:sha(__file__.replace('\\','/').split('/')[-1])})
    return dict(original_source_packet=packet,original_inlet_history_covers=state.original,
        original_axis_pressure_over_Pstar_squared=packet.P0,generic_inlet_defect_exact_zero=state.defect,
        ordinary_logR_rows_and_physical_half_shift_preserved=True,
        source_family=service.family,source_query_coordinate=phase,
        exact_formal_radius_binding=dict(radius='Ra*exp(hb*s_c/2)',
            source=domain['formal_left_log_radius_offset'],
            positive_microscopic_offset_not_added_to_huge_log_Ra=True),
        left_edge_loop_q_A_B_and_all_derivatives_exact_zero_by_current_collar=True,
        original_P0_and_all_five_incoming_histories_retained=True,
        returned_rows_are_covers_not_selected_point_function_values=True,
        source_graph_ancestor_constructors_called=False,
        full_physical_point_field_or_seam_identity_admitted=False,input_hashes=service.hashes)
