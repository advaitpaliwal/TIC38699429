# consistency of per-sector PDCSAP depths; SAP/PDCSAP ratio vs CROWDSAP; plot
import pandas as pd, numpy as np, json
from scipy.stats import chi2 as C
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
df=pd.read_csv('evidence/03_depth_by_product_sector.csv'); out={}
for p in ['TESS-SPOC PDCSAP_FLUX','TESS-SPOC SAP_FLUX','TGLC cal_psf_flux']:
    g=df[df['product']==p]; w=1/g.err_ppm**2; m=(g.depth_ppm*w).sum()/w.sum(); c2=float((((g.depth_ppm-m)/g.err_ppm)**2).sum()); n=len(g)-1
    out[p]=dict(mean=m,chi2=c2,dof=n,p_value=float(C.sf(c2,n)))
pd_=df[df['product']=='TESS-SPOC PDCSAP_FLUX'].set_index('sector'); sa=df[df['product']=='TESS-SPOC SAP_FLUX'].set_index('sector')
r=sa.depth_ppm/pd_.depth_ppm; cr=pd_.CROWDSAP
out['SAP_over_PDC_vs_CROWDSAP']=dict(pearson_r=float(np.corrcoef(r,cr)[0,1]),median_ratio_over_crowdsap=float(np.median(r/cr)),rows={int(k):[round(float(r[k]),3),round(float(cr[k]),3)] for k in r.index})
json.dump(out,open('evidence/03b_depth_consistency.json','w'),indent=1); print(json.dumps({k:v for k,v in out.items() if k!='SAP_over_PDC_vs_CROWDSAP'},indent=1)); print({k:v for k,v in out['SAP_over_PDC_vs_CROWDSAP'].items() if k!='rows'})
fig,ax=plt.subplots(1,2,figsize=(13,4.5))
for p,c in [('TESS-SPOC PDCSAP_FLUX','k'),('TESS-SPOC SAP_FLUX','C1'),('TGLC cal_psf_flux','C2'),('eleanor-lite CORR_FLUX','C4')]:
    g=df[df['product']==p]; ax[0].errorbar(g.sector,g.depth_ppm,g.err_ppm,fmt='o',color=c,label=p,ms=4)
ax[0].legend(fontsize=8); ax[0].set_xlabel('sector'); ax[0].set_ylabel('depth (ppm)'); ax[0].set_title('depth per sector and product')
ax[1].plot(cr,r,'ko'); ax[1].plot([0.5,0.95],[0.5,0.95],'r:'); ax[1].set_xlabel('SPOC CROWDSAP (target share of aperture light)'); ax[1].set_ylabel('SAP depth / PDCSAP depth'); ax[1].set_title('SAP vs PDCSAP: ratio tracks crowding')
plt.tight_layout(); plt.savefig('plots/depth_by_product.png',dpi=100)
