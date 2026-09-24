"""Check every shipped file hash. Integrity is not PDE admission."""
from pathlib import Path
import hashlib,json,sys
ROOT=Path(__file__).resolve().parent

def verify():
    manifest=json.loads((ROOT/'MANIFEST.json').read_text())
    entries=manifest['files']
    for relative,expected in entries.items():
        p=(ROOT/relative).resolve()
        if not p.is_relative_to(ROOT) or not p.is_file():raise ValueError(f'Missing/escaping path: {relative}')
        actual=hashlib.sha256(p.read_bytes()).hexdigest()
        if actual!=expected:raise ValueError(f'Checksum mismatch: {relative}')
    print(f'Verified {len(entries)} payload files. Original global NS target remains unvalidated.')
    return len(entries)
if __name__=='__main__':
    try:verify()
    except Exception as e:print(str(e),file=sys.stderr);raise SystemExit(1)
