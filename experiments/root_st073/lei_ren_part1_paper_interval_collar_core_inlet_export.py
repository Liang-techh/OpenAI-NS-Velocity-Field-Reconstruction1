"""Export trusted accepted finite-core atoms as portable exact JSON input."""
import hashlib
import json
import pickle
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_interval_collar_core_inlet import LOCAL_PICKLE_CACHE,_read_cache


def run():
    data_bytes=Path(LOCAL_PICKLE_CACHE).read_bytes();raw=pickle.loads(data_bytes)
    coefficients=raw['coefficients']
    stored=dict(format='exact-finite-core-atoms-v1',Z=raw['Z'],Z_depth=raw['Z_depth'],
                precision=coefficients['precision'],source_pickle_sha256=hashlib.sha256(data_bytes).hexdigest(),
                Lambda_exact_mpf_tuple=list(raw['Lambda']._mpf_),delta_exact_mpf_tuple=list(raw['delta']._mpf_),
                coefficients={name:[[{str(p):list(v._mpf_) for p,v in entry.atoms.items()} for entry in row]
                                    for row in coefficients[name]] for name in ('F','Uz','P')},
                stored_finite_model_only=True,original_source_errors_enclosed=False)
    path=Path(__file__).with_name('lei_ren_part1_paper_interval_collar_core_inlet_input.json')
    path.write_text(json.dumps(stored,indent=2)+'\n')
    # Confirm exact binary scalar identity before accepting this portable copy.
    decoded,_,_=_read_cache(path)
    for name in ('F','Uz','P'):
        for row_index,(row,encoded_row) in enumerate(zip(coefficients[name],stored['coefficients'][name])):
            for entry_index,(entry,encoded) in enumerate(zip(row,encoded_row)):
                for p,v in entry.atoms.items():
                    if mp.make_mpf(tuple(encoded[str(p)]))!=v:raise AssertionError('Input atom changed')
                    if decoded['coefficients'][name][row_index][entry_index].atoms[p]!=v:
                        raise AssertionError('JSON loader changed input atom at default precision')
    print('portable exact input',path.name,'bytes',path.stat().st_size,flush=True)


if __name__=='__main__':run()
