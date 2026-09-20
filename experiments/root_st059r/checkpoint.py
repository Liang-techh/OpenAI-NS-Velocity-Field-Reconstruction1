"""Transactional numerical checkpoints. Each generation stores actual arrays.

POSIX fsync + directory rename commits a complete generation before advancing
LATEST. Resume scans complete generations, not an untrusted latest pointer.
This is local durable storage, not a guarantee about external host retention.
"""
from __future__ import annotations
import hashlib, io, json, os, shutil, tempfile
from pathlib import Path
import numpy as np

class IntegrityError(RuntimeError):
    pass

def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()

def filehash(path: Path) -> str:
    return digest(Path(path).read_bytes())

def syncdir(path: Path) -> None:
    fd = os.open(path, os.O_RDONLY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)

def atomic_bytes(path: Path, data: bytes) -> None:
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix='.'+path.name+'-', dir=path.parent)
    try:
        with os.fdopen(fd, 'wb') as f:
            f.write(data); f.flush(); os.fsync(f.fileno())
        os.replace(tmp, path); syncdir(path.parent)
    finally:
        if os.path.exists(tmp): os.unlink(tmp)

def atomic_json(path: Path, data: dict) -> None:
    atomic_bytes(path, (json.dumps(data, indent=2, allow_nan=False)+'\n').encode())

def npy_bytes(a: np.ndarray) -> bytes:
    b=io.BytesIO(); np.save(b, a, allow_pickle=False); return b.getvalue()

class Store:
    def __init__(self, root: Path):
        self.root=Path(root)
    def initialize(self, binding: dict, arrays: dict[str,np.ndarray], parent: bytes) -> None:
        if self.root.exists() and any(self.root.iterdir()):
            raise IntegrityError('Refuse nonempty checkpoint store')
        self.root.mkdir(parents=True, exist_ok=True)
        atomic_bytes(self.root/'parent.json',parent)
        for name,a in arrays.items():
            if not name.isidentifier(): raise ValueError('Invalid array name')
            atomic_bytes(self.root/f'{name}.npy',npy_bytes(a))
        binding=dict(binding)
        binding['immutable_files']={p.name:filehash(p) for p in self.root.iterdir() if p.is_file()}
        atomic_json(self.root/'binding.json',binding)
    def verify(self) -> dict:
        b=json.loads((self.root/'binding.json').read_text())
        for name,h in b['immutable_files'].items():
            if Path(name).name!=name or filehash(self.root/name)!=h:
                raise IntegrityError('Immutable checkpoint input changed: '+name)
        return b
    def save(self, iteration: int, delta: np.ndarray, candidate: bytes, state: dict,
             _interrupt_before_rename: bool=False) -> Path:
        b=self.verify(); name=f'step-{iteration:06d}'; final=self.root/name
        if final.exists(): raise IntegrityError('Refuse overwrite of committed generation')
        tmp=Path(tempfile.mkdtemp(prefix='.pending-',dir=self.root))
        data={'delta.npy':npy_bytes(delta), 'candidate.json':candidate,
              'state.json':(json.dumps(state,indent=2,allow_nan=False)+'\n').encode()}
        for n,v in data.items(): atomic_bytes(tmp/n,v)
        meta={'iteration':iteration,'binding_sha256':filehash(self.root/'binding.json'),
              'files':{n:digest(v) for n,v in data.items()},'pde_validated':False}
        atomic_json(tmp/'manifest.json',meta); syncdir(tmp)
        if _interrupt_before_rename: raise InterruptedError('Simulated process interruption before commit')
        os.replace(tmp,final); syncdir(self.root)
        atomic_json(self.root/'LATEST.json',{'generation':name,'manifest_sha256':filehash(final/'manifest.json')})
        return final
    def load(self) -> tuple[int,np.ndarray,bytes,dict]:
        self.verify(); generations=sorted(self.root.glob('step-*'))
        if not generations: raise IntegrityError('No committed numerical checkpoint')
        p=generations[-1]; h=json.loads((p/'manifest.json').read_text())
        if h['binding_sha256']!=filehash(self.root/'binding.json'):
            raise IntegrityError('Binding mismatch')
        for n,s in h['files'].items():
            if Path(n).name!=n or filehash(p/n)!=s: raise IntegrityError('Committed generation corrupted: '+n)
        d=np.load(p/'delta.npy',allow_pickle=False)
        if not np.isfinite(d).all(): raise IntegrityError('Nonfinite delta')
        return h['iteration'],d,(p/'candidate.json').read_bytes(),json.loads((p/'state.json').read_text())
