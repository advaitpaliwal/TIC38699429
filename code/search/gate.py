# usage: gate.py [--shard I/N] [--only LANE:TIC ...]   (shards only fill the cache; run with no flags to write the CSVs;
#        --from-cache writes the CSVs from cached peaks only and lists the peaks still to vet in pending.csv)
# Runs LEO-vetter 1.2.0 (flux-level tests, published default thresholds) on every v4 peak with ok==True and SDE>=9.
# Inputs per peak: the searched light curve with v4 prep (leo_glue.prep), earlier peaks of that star masked as the
# search did, TLS P and T0 from res_v4, and the v4 expected central-transit duration (TLS's own duration is not stored).
# Caches one JSON per peak in cache/, writes candidates.csv (passing peaks) and all_peaks.csv (every gated peak).
import os,sys,json,glob,time,numpy as np,pandas as pd
H=os.path.dirname(os.path.abspath(__file__));sys.path.insert(0,H)
from leo_glue import prep,star_dict,mask_prev,vet,tdur
B=os.path.abspath(os.path.join(H,'..','..'));LANES=['hz-north','hz-south','hz-near','hz-wide']
def pool(lane):return pd.read_csv(os.path.join(B,lane,'pool.csv')).set_index('ID')
def vet_star(lane,tic,row,res,cache_only=False):
    out=[];lc=None;prev=[]
    R,M,Te=float(res['R']),float(res['M']),float(res['Teff']);lg=float(row.logg) if row is not None and 'logg' in row else None
    for p in res['peaks']:
        gpu=res.get('settings')=='v4gpu'
        if p['ok'] and p['SDE']>=(21.7 if gpu else 9):  # GPU score scale: candidate threshold from inverted-flux calibration (GPU.md)
            cf=os.path.join(H,'cache',f"{lane}_{tic}_{'gpu' if gpu else ''}{p['iter']}.json")
            if os.path.exists(cf):o=json.load(open(cf))
            elif cache_only:o=dict(verdict='PENDING',fails=[],metrics={})
            else:
                if lc is None:lc=prep(os.path.join(B,lane,'lc',f'{tic}.npz'))
                t,f=mask_prev(*lc,prev,R,M);dur=tdur(p['P'],R,M);t0=time.time()
                try:o=vet(tic,t,f,p['P'],p['T0'],dur,star_dict(tic,R,M,Te,lg))
                except Exception as x:o=dict(verdict='ERROR',fails=[repr(x)[:300]],metrics={})
                o.update(lane=lane,tic=tic,iter=p['iter'],P=p['P'],T0=p['T0'],dur_d=dur,secs=round(time.time()-t0,1))
                json.dump(o,open(cf,'w'))
            out.append(dict(lane=lane,tic=tic,iter=p['iter'],P=p['P'],T0=p['T0'],depth_ppm=p['depth_ppm'],rp_re=p['rp_re'],S_earth=p['S_earth'],SDE=p['SDE'],leo=o['verdict'],failing_tests='; '.join(o['fails']),leo_MES=o['metrics'].get('MES'),leo_Rp=o['metrics'].get('Rp')))
        prev.append((p['P'],p['T0']))
    return out
if __name__=='__main__':
    only=set(sys.argv[sys.argv.index('--only')+1:]) if '--only' in sys.argv else None
    sh=tuple(map(int,sys.argv[sys.argv.index('--shard')+1].split('/'))) if '--shard' in sys.argv else None
    rows=[];k=-1
    for lane in LANES:
        pl=pool(lane)
        for fj in sorted(glob.glob(os.path.join(B,lane,'res_v4','*.json'))+glob.glob(os.path.join(B,lane,'res_v4_gpu','*.json'))):
            tic=int(os.path.basename(fj)[:-5])
            if only and f'{lane}:{tic}' not in only:continue
            k+=1
            if sh and k%sh[1]!=sh[0]:continue
            try:res=json.load(open(fj))
            except Exception:continue
            if 'peaks' not in res:continue
            rows+=vet_star(lane,tic,pl.loc[tic] if tic in pl.index else None,res,'--from-cache' in sys.argv)
            if '--from-cache' not in sys.argv:print(lane,tic,len(rows),flush=True)
    df=pd.DataFrame(rows)
    if only is None and sh is None:
        df.to_csv(os.path.join(H,'all_peaks.csv'),index=False)
        df[df.leo=='PC'].to_csv(os.path.join(H,'candidates.csv'),index=False)
        df[df.leo=='PENDING'][['lane','tic','iter','P','SDE']].to_csv(os.path.join(H,'pending.csv'),index=False)
    print(df.groupby(['lane','leo']).size() if len(df) else 'no peaks')
