# Step 1c + 4: depth per sector and per data product (TESS-SPOC SAP vs PDCSAP, TGLC PSF/aperture, eleanor-lite), local linear baseline per event
import numpy as np, glob, json, re, warnings; warnings.filterwarnings('ignore')
from astropy.io import fits
P=11.718154156303166; T0=1329.529459184984; W=0.05  # full width in days for in-transit
def depth(t,f):
    ok=np.isfinite(t)&np.isfinite(f); t,f=t[ok],f[ok]
    if len(t)<50: return None
    f=f/np.nanmedian(f)
    ev=[];
    for n in range(int(np.floor((t.min()-T0)/P)),int(np.ceil((t.max()-T0)/P))+1):
        dt=t-(T0+n*P); a=np.abs(dt)<W/2; b=(np.abs(dt)>0.08)&(np.abs(dt)<0.5)
        if a.sum()<1 or b.sum()<10 or (dt[b]<0).sum()<3 or (dt[b]>0).sum()<3: continue
        c=np.polyfit(dt[b],f[b],1); r=f[b]-np.polyval(c,dt[b]); s=1.4826*np.median(np.abs(r))
        ev.append(((np.polyval(c,dt[a])-f[a]).mean(), s/np.sqrt(a.sum())))
    if not ev: return None
    e=np.array(ev); w=1/e[:,1]**2
    return dict(depth_ppm=float((e[:,0]*w).sum()/w.sum()*1e6), err_ppm=float(1/np.sqrt(w.sum())*1e6), n_ev=len(ev))
rows=[]
for p in sorted(glob.glob('spoc/**/*.fits',recursive=True)):
    h=fits.open(p); d=h[1].data; nm=p.split('/')[-1]; s=int(re.search(r's00(\d\d)',nm).group(1))
    if 'tess-spoc' in nm:
        q=d['QUALITY']==0; t=d['TIME'][q]
        hd=h[1].header
        for col in ['SAP_FLUX','PDCSAP_FLUX']:
            r=depth(t,d[col][q]); 
            if r: rows.append(dict(product='TESS-SPOC '+col,sector=s,CROWDSAP=hd.get('CROWDSAP'),FLFRCSAP=hd.get('FLFRCSAP'),**r))
    elif 'tglc' in nm:
        q=(d['TESS_flags']==0)&(d['TGLC_flags']==0); t=d['time'][q]
        for col in ['cal_psf_flux','cal_aper_flux']:
            r=depth(t,d[col][q]); 
            if r: rows.append(dict(product='TGLC '+col,sector=s,**r))
    elif 'eleanor' in nm:
        q=d['QUALITY']==0; t=d['TIME'][q]
        for col in ['CORR_FLUX','PCA_FLUX']:
            r=depth(t,d[col][q]); 
            if r: rows.append(dict(product='eleanor-lite '+col,sector=s,**r))
    h.close()
import pandas as pd
df=pd.DataFrame(rows); df.to_csv('evidence/03_depth_by_product_sector.csv',index=False)
def comb(g):
    w=1/g.err_ppm**2; return pd.Series(dict(depth_ppm=(g.depth_ppm*w).sum()/w.sum(),err_ppm=1/np.sqrt(w.sum()),n_sectors=len(g),n_ev=g.n_ev.sum()))
yr=lambda s:'Y1' if s<=13 else 'Y3' if s<=39 else 'Y5'
df['year']=df.sector.map(yr)
summ=df.groupby(['product','year']).apply(comb).reset_index(); tot=df.groupby('product').apply(comb).reset_index(); tot['year']='all'
out=pd.concat([summ,tot]); out.to_csv('evidence/03_depth_by_product_summary.csv',index=False)
print(out.round(0).to_string())
sp=df[df['product'].str.startswith('TESS-SPOC')].pivot(index='sector',columns='product',values=['depth_ppm','err_ppm']).round(0)
sp['CROWDSAP']=df[df['product']=='TESS-SPOC PDCSAP_FLUX'].set_index('sector').CROWDSAP
sp['ratio_SAP_over_PDC']=sp[('depth_ppm','TESS-SPOC SAP_FLUX')]/sp[('depth_ppm','TESS-SPOC PDCSAP_FLUX')]
print(sp.to_string())
