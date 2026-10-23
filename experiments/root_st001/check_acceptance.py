"""Exit nonzero unless every registered numerical gate in the report passes."""
import argparse
import json
from pathlib import Path

def main():
    p=argparse.ArgumentParser();p.add_argument('report',type=Path);a=p.parse_args()
    r=json.loads(a.report.read_text());required=('momentum_max','momentum_L2','divergence_max','divergence_L2','initial_energy','energy_range','energy_quadrature','boundary','core_signs','scaled_core_drift','derived_amplitude')
    missing=[k for k in required if k not in r.get('gates',{})]
    failed=[k for k in required if r.get('gates',{}).get(k) is not True]
    if missing or failed:
        print('FAIL: '+', '.join(failed));raise SystemExit(1)
    print('All reported numerical gates pass; this is not a continuum proof.')
if __name__=='__main__':main()
