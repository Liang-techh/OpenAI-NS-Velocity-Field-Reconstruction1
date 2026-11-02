# Uncheckpointed seam bridge v1 synthetic execution — supplemental record

The worker was instructed not to execute the new source consumer before freezing. Its initial report asserted that the synthetic fixture passed, then disclosed that it had in fact run multiple uncheckpointed synthetic-only entries, including a failed fixture and failed first four-argument version. It reports no actual receipt or ancestor was accessed. Source bytes at each entry cannot be authenticated retrospectively. Never re-enter the v1 family.

The worker's command prefix (P) was:

```powershell
python -c "import importlib.util; p='src/openai_ns_reconstruction/core_z2_seam_bridge_v1.py'; f='src/openai_ns_reconstruction/core_z2_seam_bridge_fixture_v1.py'; s=importlib.util.spec_from_file_location('bridge',p); m=importlib.util.module_from_spec(s); s.loader.exec_module(m); sf=importlib.util.spec_from_file_location('fixture',f); q=importlib.util.module_from_spec(sf); sf.loader.exec_module(q); 
```

The worker reported these entries, in order. Each Python command after 1 consists of P followed by the shown suffix and a closing double quote:

1. `python -m py_compile src/openai_ns_reconstruction/core_z2_seam_bridge_v1.py src/openai_ns_reconstruction/core_z2_seam_bridge_fixture_v1.py`; exit 0, empty output.
2. `out=q.run_fixture(m.build); print(out['passed']); print(out['checks'])`: `False`, with only `two_fields_per_unknown` false among 13 checks.
3. `old,pro,seam=q.synthetic_inputs(); r=m.build(old,pro,seam); print(r['bindings'][0]['product_ids_by_field']); print(r['bindings'][0]['product_counts_by_field'])`: printed two products each for `F1` and `d_x_F1`.
4. Repeat entry 2 after fixture correction: `True`, 13 checks passed.
5. `out=q.run_fixture(m.build); print(out['passed']); print(out['checks']); old,pro,seam=q.synthetic_inputs(); low=q.synthetic_low_output(); print(m.build(old,pro,seam,low)['product_count'])`: failed with `ValueError: new seam F1[0].input_matrix_coefficient directed endpoints differ`; fixture reported `False` and `{'build_completed': False}`. The later print did not run.
6. Diagnostic suffix `old,pro,seam=q.synthetic_inputs(); low=q.synthetic_low_output(); cells=pro['target_cells'];\nfor fld in ('F1','d_x_F1'):\n for i,t in enumerate(seam['terms_by_order']['Z1'][fld]):\n  try:\n   sid=t['state_profile_term_id']; st=None\n   for c in cells:\n    for f2,ts in c['W1_Z2'].items():\n     for x in ts:\n      if x['term_id']==sid: st=x\n   m._validate_low_payload(low_output=low,product=t,profile_term=st,label=f'{fld}[{i}]')\n  except Exception as e:\n   print('FAIL',fld,i,t['matrix_receipt_input_path'],t['matrix_key'],t['source_cell'],e); raise\nprint('ok')`: reported `FAIL F1 0` and the same endpoint mismatch, and raised `ValueError`. Worker supplied only its own textual traceback; no retained tool result artifact exists.
7. Repeat entry 5 after code correction: `True`, 15 checks passed, product count `96`.
8. `old,pro,seam=q.synthetic_inputs(); low=q.synthetic_low_output(); oldw={'result':{'output':{'result':old}}}; neww={'result':{'output':{'result':pro,'missing_seam_result':seam}}}; loww={'result':{'output':low}}; r=m.build(oldw,neww,neww,loww); print(r['product_count'], r['low_output_consumed'])`: output `96 True`.

The original incident JSON preserves the worker's initial partial account and the observed later source bytes; the four-argument final code was written after earlier calls, so the archive cannot establish the exact source used in each invocation. A separate new v2 will use fresh synthetic inputs and the frozen one-entry protocol.
