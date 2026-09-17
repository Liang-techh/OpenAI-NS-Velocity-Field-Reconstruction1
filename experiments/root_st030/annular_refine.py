"""Memory-bounded ST031 quadrature check; same integrands and mode space."""
from __future__ import annotations
import argparse
import numpy as np
import nonaxisymmetric_screen as screen

_original = screen.base_jets

def chunked_base_jets(f, raw, x, t):
    blocks = [_original(f, raw, x[i:i+512], t) for i in range(0, len(x), 512)]
    return tuple(np.concatenate([b[j] for b in blocks], axis=0) for j in range(5))

def main():
    p=argparse.ArgumentParser();p.add_argument('candidate');p.add_argument('--out',required=True);p.add_argument('--order',type=int,default=48);a=p.parse_args()
    screen.base_jets=chunked_base_jets
    screen.run(a.candidate,a.out,order=a.order,energy_order=64)
if __name__=='__main__':main()
