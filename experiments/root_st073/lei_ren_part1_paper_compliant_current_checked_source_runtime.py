"""Opt-in runtime reuse of one already accepted symbolic source theorem.

This context does not hydrate numerical owners from report covers. Original
constructors, coordinate queries, acceptance checks and operators still run.
Only the exact collar/Gamma Boolean proof is reused after its entire recorded
source closure has been checked. Digest reuse is invalidated by file metadata.
The normal producer/checker path is untouched outside this serial context.
"""
import copy
import hashlib
import json
import threading
from pathlib import Path

HERE = Path(__file__).parent
ROOT = HERE.parents[1]
PREFIX = 'lei_ren_part1_paper_compliant_'
THEOREM = PREFIX + 'collar_Gamma_C4_check.json'


class CheckedSourceRuntime:
    """Serial, bounded source construction; no field-value or owner cache."""
    _active = None

    def __init__(self):
        self.digests = {}
        self.reused_proofs = 0
        self.digest_reads = 0
        self.digest_hits = 0
        self._entered = False

    def digest(self, name):
        path = (HERE / name).resolve()
        if not path.is_relative_to(ROOT.resolve()):
            raise ValueError('Runtime digest must remain in the original repository')
        st = path.stat()
        stamp = (st.st_dev, st.st_ino, st.st_size, st.st_mtime_ns, st.st_ctime_ns)
        old = self.digests.get(path)
        if old is not None and old[0] == stamp:
            self.digest_hits += 1
            return old[1]
        value = hashlib.sha256(path.read_bytes()).hexdigest()
        after = path.stat()
        if stamp != (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns, after.st_ctime_ns):
            raise ValueError('Source changed while its runtime digest was read: ' + name)
        self.digests[path] = (stamp, value)
        self.digest_reads += 1
        return value

    def validate_proof(self):
        record = json.loads((HERE / THEOREM).read_bytes())
        if not record.get('all_passed') or not record.get('waiting_collar_and_collar_Gamma_joins_certified'):
            raise ValueError('Accepted original collar/Gamma source theorem required')
        closure = record.get('input_hashes', {})
        required = {PREFIX + 'collar_Gamma_C4_check.py', PREFIX + 'collar_Gamma_C4.py'}
        if not required <= set(closure):
            raise ValueError('Original symbolic proof and defining source must be hash-bound')
        for name, value in closure.items():
            if self.digest(name) != value:
                raise ValueError('Accepted original source changed: ' + name)
        proof = record.get('functional_production_source_identities')
        if not isinstance(proof, dict) or not proof or any(value is not True for value in proof.values()):
            raise ValueError('Complete original functional identity result required')
        self.theorem_sha = self.digest(THEOREM)
        self.proof = proof
        self.closure = dict(closure)

    def reuse_proof(self):
        if not self._entered or threading.get_ident() != self.thread:
            raise RuntimeError('Checked theorem reuse is limited to its active serial context')
        # A changed input must fail even after the first successful reuse.
        if self.digest(THEOREM) != self.theorem_sha:
            raise ValueError('Accepted collar/Gamma receipt changed during construction')
        for name, value in self.closure.items():
            if self.digest(name) != value:
                raise ValueError('Accepted original source changed: ' + name)
        self.reused_proofs += 1
        return copy.deepcopy(self.proof)

    def __enter__(self):
        if self._entered or CheckedSourceRuntime._active is not None:
            raise RuntimeError('Checked source runtime contexts cannot overlap or nest')
        import lei_ren_part1_paper_compliant_macro_signed_integrals as macro
        import lei_ren_part1_paper_compliant_collar_Gamma_C4_check as theorem
        import lei_ren_part1_paper_compliant_current_heat_source as heat
        import lei_ren_part1_paper_compliant_current_postpulse_interfaces as postpulse
        if heat.functional_source_identities is not theorem.functional_source_identities:
            raise ValueError('Original symbolic theorem alias differs before runtime entry')
        if postpulse.heat_theorem is not theorem.functional_source_identities:
            raise ValueError('Original postpulse theorem alias differs before runtime entry')
        self.validate_proof()
        self.modules = (macro, theorem, heat, postpulse)
        self.original_digest = macro._sha256
        self.original_proof = theorem.functional_source_identities
        self.thread = threading.get_ident()
        self.digest_callable = self.digest
        self.proof_callable = self.reuse_proof
        macro._sha256 = self.digest_callable
        theorem.functional_source_identities = self.proof_callable
        heat.functional_source_identities = self.proof_callable
        postpulse.heat_theorem = self.proof_callable
        self._entered = True
        CheckedSourceRuntime._active = self
        return self

    def __exit__(self, exc_type, exc, tb):
        macro, theorem, heat, postpulse = self.modules
        changed = (macro._sha256 is not self.digest_callable or
                   theorem.functional_source_identities is not self.proof_callable or
                   heat.functional_source_identities is not self.proof_callable or
                   postpulse.heat_theorem is not self.proof_callable)
        macro._sha256 = self.original_digest
        theorem.functional_source_identities = self.original_proof
        heat.functional_source_identities = self.original_proof
        postpulse.heat_theorem = self.original_proof
        self._entered = False
        CheckedSourceRuntime._active = None
        if changed and exc_type is None:
            raise RuntimeError('Runtime proof/digest bindings changed during source construction')
        return False

    def record(self):
        return dict(accepted_symbolic_theorem=THEOREM, theorem_sha256=self.theorem_sha,
            validated_source_hash_count=len(self.closure), reused_symbolic_proofs=self.reused_proofs,
            source_digest_reads=self.digest_reads, unchanged_source_digest_hits=self.digest_hits,
            numerical_owners_or_field_values_cached=False,
            original_numerical_constructors_and_coordinate_operators_unchanged=True,
            receipt_scope_or_global_completion_gates_promoted=False,
            serial_context_bindings_restored=not self._entered)


def native_bridge_owner(seed=None):
    """One native source seed and its checked bridge, using explicit injection.

    The full downstream tensor chain is unnecessary for original source rows.
    An optional seed must be the real accepted original physical assembly.
    """
    from lei_ren_part1_paper_compliant_current_core_physical_assembly import CurrentCorePhysicalAssembly
    from lei_ren_part1_paper_compliant_current_actual_bridge_mixed_C4 import CurrentActualBridgeMixedC4
    with CheckedSourceRuntime() as runtime:
        if seed is None:
            seed = CurrentCorePhysicalAssembly()
        if (type(seed) is not CurrentCorePhysicalAssembly or
                not seed.core_acceptance_loaded or not seed.physical_acceptance_loaded or
                not all(seed.current_provider_graph().values())):
            raise ValueError('Exact checked original core/heat physical seed required')
        bridge = CurrentActualBridgeMixedC4(owner=seed)
        if not bridge.acceptance_loaded or not all(bridge.current_provider_graph().values()):
            raise ValueError('Original native bridge construction did not load acceptance')
    return bridge, runtime.record()
