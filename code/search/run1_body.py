# usage: run1.py TIC R M OUTDIR [threads]  -> per-sector best LC (SPOC 2-min > TESS-SPOC > QLP), wotan biweight, TLS 1-60 d, 2 iterations
import sys, os, json, time, warnings, subprocess, shutil, numpy as np, pandas as pd
warnings.filterwarnings('ignore')
import lightkurve as lk
from wotan import flatten
from transitleastsquares import transitleastsquares, transit_mask, cleaned_array
R0=os.path.dirname(os.path.abspath(sys.argv[0]))
tic=int(sys.argv[1]); R=float(sys.argv[2]); M=float(sys.argv[3]); outdir=sys.argv[4]; th=int(sys.argv[5]) if len(sys.argv)>5 else 1
out=f'{outdir}/{tic}.json'
if os.path.exists(out): sys.exit()
t0=time.time()
prod=pd.read_csv(f'{R0}/_sources/mast_lc_products.csv'); prod=prod[prod.target_name.astype(int)==tic]
ts=f'{tic:016d}'; p=f'{ts[0:4]}/{ts[4:8]}/{ts[8:12]}/{ts[12:16]}'
rank={'SPOC':0,'TESS-SPOC':1,'QLP':2}
prod=prod.assign(r=prod.provenance_name.map(rank)).sort_values(['sequence_number','r']).drop_duplicates('sequence_number')
lcd=f'{R0}/lc/{tic}'; os.makedirs(lcd,exist_ok=True)
def get(url,fn):
    r=subprocess.run(['curl','-sf','-m','120','--retry','4','--retry-max-time','240','-o',fn,url])
    if r.returncode!=0:
        if os.path.exists(fn): os.remove(fn)
        return None
    return fn
def urls(row):
    s=int(row.sequence_number)
    if row.provenance_name=='SPOC':
        return [(f'https://archive.stsci.edu/missions/tess/tid/s{s:04d}/{p}/{row.obs_id}_lc.fits','pdcsap_flux','SPOC')]
    if row.provenance_name=='TESS-SPOC':
        return [(f'https://archive.stsci.edu/hlsps/tess-spoc/s{s:04d}/target/{p}/hlsp_tess-spoc_tess_phot_{ts}-s{s:04d}_tess_v1_lc.fits','pdcsap_flux','TESS-SPOC')]
    return [(f'https://archive.stsci.edu/hlsps/qlp/s{s:04d}/{p}/hlsp_qlp_tess_ffi_s{s:04d}-{ts}_tess_v01_llc.fits','kspsap_flux','QLP')]
T=[];F=[];S=[];used=[]
try:
    for _,row in prod.iterrows():
        s=int(row.sequence_number); lc=None
        for u,col,src in urls(row):
            fn=get(u,f'{lcd}/s{s}.fits')
            if not fn: continue
            for c in (col,'det_flux','sap_flux'):
                try: lc=lk.read(fn,flux_column=c,quality_bitmask='hard').remove_nans(); break
                except Exception: lc=None
            os.remove(fn)
        if lc is None or len(lc)<100: continue
        t=lc.time.value; f=np.asarray(lc.flux.value,float); f=f/np.nanmedian(f)
        cad=np.nanmedian(np.diff(t))
        fl=flatten(t,f,method='biweight',window_length=0.75,edge_cutoff=0.,break_tolerance=0.5)
        k=np.isfinite(fl); t,fl=t[k],fl[k]
        # bin to 30 min
        if cad<25/1440:
            b=np.floor((t-t[0])/(30/1440)).astype(int); n=np.bincount(b)
            ok=n>0; tb=np.bincount(b,t)[ok]/n[ok]; fb=np.bincount(b,fl)[ok]/n[ok]; t,fl=tb,fb
        mad=np.nanmedian(np.abs(fl-np.nanmedian(fl)))*1.4826
        k=fl<1+3*mad; t,fl=t[k],fl[k]
        T.append(t);F.append(fl);S.append(np.full(len(t),s)); used.append((s,src,round(cad*1440,1)))
    shutil.rmtree(lcd,ignore_errors=True)
    if len(used)<8: raise RuntimeError(f'only {len(used)} sectors')
    t=np.concatenate(T); f=np.concatenate(F); sid=np.concatenate(S); o=np.argsort(t); t,f,sid=t[o],f[o],sid[o]
    np.savez_compressed(f'{R0}/lc/{tic}.npz',t=t,f=f,s=sid)
    if os.environ.get("LC_ONLY"):sys.exit(0)
    NIT=int(os.environ.get("NIT","3")); res=[]; tt,ff=t.copy(),f.copy()
    for it in range(NIT):
        m=transitleastsquares(tt,ff)
        r=m.power(period_min=1.0,period_max=60.0,n_transits_min=3,R_star=R,M_star=M,R_star_min=0.08,R_star_max=0.8,M_star_min=0.08,M_star_max=0.8,
                  oversampling_factor=int(os.environ.get("OS","1")),duration_grid_step=1.1,use_threads=th,show_progress_bar=False)
        if not np.isfinite(r.SDE): break
        pk=np.argsort(r.power)[::-1]; peaks=[]
        for j in pk:
            if all(abs(r.periods[j]-q)/q>0.01 for q,_ in peaks): peaks.append((float(r.periods[j]),float(r.power[j])))
            if len(peaks)==5: break
        res.append(dict(iter=it,SDE=float(r.SDE),P=float(r.period),T0=float(r.T0),depth=float(1-r.depth),dur=float(r.duration),snr=float(r.snr),
            oe=float(r.odd_even_mismatch),ntr=int(r.distinct_transit_count),ntr_data=int(sum(np.asarray(r.per_transit_count)>0)),rp_rs=float(r.rp_rs),top_peaks=peaks))
        k=~transit_mask(tt,r.period,max(3*r.duration,0.25),r.T0); tt,ff=tt[k],ff[k]
    d=dict(tic=tic,R=R,M=M,used=used,nsec=len(used),npts=len(t),std=float(np.std(f)),results=res,secs=round(time.time()-t0))
except Exception as e:
    shutil.rmtree(lcd,ignore_errors=True)
    d=dict(tic=tic,used=used,err=repr(e)[:300],secs=round(time.time()-t0))
json.dump(d,open(out,'w'))
print(json.dumps(d)[:600])
