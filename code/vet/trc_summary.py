# Summarise TRICERATOPS runs in a directory: FPP/NFPP mean, sd (ddof=1), range, mean prob by star+scenario.
# usage: trc_summary.py DIR
import sys, glob, json, numpy as np, pandas as pd
d=sys.argv[1]; R=[json.load(open(f)) for f in sorted(glob.glob(f'{d}/res_*.json'))]
print(d,'runs',len(R),'n_stars',sorted({r['n_stars'] for r in R}),'n_tdepth>0',sorted({r['n_stars_tdepth_gt0'] for r in R}),'depth_ppm',round(R[0]['depth_ppm']),'err_ppm',round(R[0]['err_ppm']))
for r in R: print(r['run'],'FPP %.4f NFPP %.4f'%(r['FPP'],r['NFPP']))
for k in ['FPP','NFPP']:
    v=np.array([r[k] for r in R]); print(f'{k} mean {v.mean():.4f} sd {v.std(ddof=1):.4f} range {v.min():.4f}-{v.max():.4f}')
P=pd.concat([pd.read_csv(f) for f in sorted(glob.glob(f'{d}/probs_*.csv'))])
m=P.groupby(['ID','scenario'])['prob'].sum()/len(R); print(m.sort_values(ascending=False).head(14).round(4).to_string())
