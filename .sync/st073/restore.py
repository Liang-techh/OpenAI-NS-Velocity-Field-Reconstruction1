"""Restore original ST073 bytes; never publish a numerically-only match.
All text and historical reports are transferred verbatim. Deterministic binary
snapshots are recreated from immutable original inputs, then checked against
EVERY original SHA256. Padding bytes are restored without changing values.
No fitting and no new scientific claim.
"""
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import sys,json,hashlib,subprocess,os,tarfile,io,lzma,base64,zipfile,copy,platform
import numpy as np
ROOT=Path(sys.argv[1]).resolve()
TRANSPORT=Path(__file__).resolve().parent
TEXT_SHA='7f59d421deec7207009181aa7911a73cdbc42e1f8160612bc1a3acf57aaba60c'
PAD_SHA='5da7e493c692203e5c44ff8706e44cfb26b946be109884db4c91200f5427546e'
def sha(b):return hashlib.sha256(b).hexdigest()
if not ROOT.exists():
 ROOT.mkdir(parents=True)
 payload=b''.join((TRANSPORT/f'part{i:02d}').read_bytes() for i in range(5))
 if sha(payload)!=TEXT_SHA:raise ValueError('Text transport checksum mismatch')
 with tarfile.open(fileobj=io.BytesIO(payload),mode='r:xz') as archive:
  for member in archive:
   dest=(ROOT/member.name).resolve()
   if not member.isfile() or not dest.is_relative_to(ROOT):raise ValueError('Unsafe archive member')
   dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(archive.extractfile(member).read())
padding_bytes=(TRANSPORT/'padding.xz').read_bytes()
if sha(padding_bytes)!=PAD_SHA:raise ValueError('Padding transport checksum mismatch')
PADDING=json.loads(lzma.decompress(padding_bytes))
sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'upstream'))
from full_radial import FullRadialField
EXPECTED=json.loads((ROOT/'MANIFEST.json').read_text())['files']
def check(path):
 rel=path.relative_to(ROOT).as_posix();got=sha(path.read_bytes());ok=got==EXPECTED[rel]
 print(('MATCH ' if ok else 'MISMATCH ')+rel,flush=True)
 return {'path':rel,'sha256':got,'matches_original':ok}
def save(path,**arrays):
 if not path.exists():np.savez_compressed(path,**arrays)
 return check(path)
def group(task):
 ident,phase,seed,neta,nx=task;f=FullRadialField.load(ROOT/'data'/f'{ident}.json')
 rng=np.random.default_rng(seed);eta=rng.uniform(-.5,.5,neta);X=rng.uniform(0,1/64,(neta,nx));es=np.broadcast_to(eta[:,None],X.shape);kk=np.r_[0,3,6,np.sort(rng.uniform(.02,5.98,5))];records=[]
 for j,k in enumerate(kk):
  path=ROOT/'data'/f'{phase}holdout_{seed}_{j}.npz'
  if path.exists():records.append(check(path));continue
  d=f.evaluate_similarity(X,es,.5*2**(-k))
  if phase:r=save(path,X=X,eta=es,k=k,**d)
  else:r=save(path,X=X,eta=es,k=k,velocity=d['velocity'],pressure=d['pressure'],residual=d['residual'],divergence=d['divergence'])
  records.append(r)
 return records
def restore_padding(path):
 rel=path.relative_to(ROOT).as_posix();record=PADDING[rel];buffer=io.BytesIO()
 with zipfile.ZipFile(path) as old,zipfile.ZipFile(buffer,'w',compression=zipfile.ZIP_DEFLATED) as new:
  for info in old.infolist():
   data=bytearray(old.read(info.filename))
   if info.filename=='coefficient_jets.npy':
    a=np.load(io.BytesIO(data),allow_pickle=False)
    if a.dtype!=np.dtype('longdouble') or a.dtype.itemsize!=16 or len(data)-a.nbytes!=record['head']:raise ValueError('Unexpected longdouble storage')
    stored=np.frombuffer(data,dtype=np.uint8,offset=record['head']).reshape(-1,16)
    stored[:,10:]=np.frombuffer(base64.b64decode(record['padding']),dtype=np.uint8).reshape(-1,6)
    if sha(data)!=record['npy_sha256']:raise ValueError('Numeric interface bytes do not match original')
   with new.open(copy.copy(info),'w',force_zip64=True) as h:h.write(data)
 if sha(buffer.getvalue())!=record['original_sha256']:raise ValueError('Restored NPZ is not byte-identical')
 path.write_bytes(buffer.getvalue())
def interfaces(ident):
 records=[];f=FullRadialField.load(ROOT/'data'/f'{ident}.json')
 for k in [0,3,6]:
  path=ROOT/'data'/f'interface_{ident}_k{k}.npz'
  if not path.exists():
   tau=.5*2**(-k);es=np.linspace(-.5,.5,17);X=np.ones_like(es)/64;d=f.evaluate_similarity(X,es,tau);q=tau/(1-es**2)
   coefs=[f.coefficients(float(e),float(qq)) for e,qq in zip(es,q)]
   np.savez_compressed(path,eta=es,X=X,tau=tau,coefficient_jets=np.asarray(coefs),**d)
  restore_padding(path);records.append(check(path))
 return records
if __name__=='__main__':
 records=[]
 tasks=[('ST073-F','',s,64,8) for s in [9237391,9237392]]+[('ST073-V','V_',s,32,16) for s in [9237395,9237396]]
 with ProcessPoolExecutor(max_workers=4) as pool:
  for result in pool.map(group,tasks):records+=result
  for result in pool.map(interfaces,['ST073-F','ST073-V']):records+=result
 from core_series import build as original
 from general_core import build as extended
 for ident,builder,swirl in [('ST068-I',original,1.),('ST073-V-leading-control',extended,4.)]:
  path=ROOT/'data'/f'{ident}.npz'
  if not path.exists():
   c,_=builder(order=18,eta_degree=64,swirl=swirl)
   if swirl==4.:
    c.meta['id']=ident;c.meta['scope']='Same axis data as V, leading-only comparison; not a new full solution'
   c.save(path)
  records.append(check(path))
 path=ROOT/'figures/full_momentum_radial_convergence.png';path.parent.mkdir(exist_ok=True)
 if not path.exists():subprocess.run([sys.executable,str(ROOT/'render_results.py')],check=True,cwd=ROOT)
 records.append(check(path))
 for name,expected in EXPECTED.items():
  p=ROOT/name
  if not p.is_file() or sha(p.read_bytes())!=expected:raise ValueError('Original payload identity mismatch: '+name)
 if len(EXPECTED)!=173 or not all(x['matches_original'] for x in records):raise ValueError('Incomplete original snapshot')
 import scipy,matplotlib,PIL,zlib
 receipt={'original_archive_sha256':'c9b3089fb727c2834479037f7bbebc2262a35b05479890e392052203f5acc917','payload_files':173,'manifest_files':1,'all_original_hashes_match':True,'binary_files':records,'environment':{'python':platform.python_version(),'platform':platform.platform(),'numpy':np.__version__,'scipy':scipy.__version__,'matplotlib':matplotlib.__version__,'pillow':PIL.__version__,'zlib':zlib.ZLIB_VERSION},'scope':'Transport restoration of exact historical bytes. No fit; no new holdout validation or global PDE acceptance.','pde_validated':False,'global_field_ready':False}
 (ROOT.parent/'st073_transport_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
 print('COMPLETE: every original payload SHA256 matches.',flush=True)
