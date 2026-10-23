# usage: run_star_gpu.py LANE TIC [TIC ...] [--out DIR] [--invert] [--stop X]
# GPU drop-in for run_star.py: same light curve (lane lc/<tic>.npz, else the lane's run1 in LC_ONLY mode), same prep
# (bin30, trim, unlone), same period range (hz-near 0.5 d, others 1 d, to min(60, span/3)), same >=3-good-epochs rule,
# up to 4 passes masking 1.5 x the expected duration. Search = cuvarbase eebls_transit_gpu(use_fast=True) on a
# transit_autofreq grid (bake-off settings), dips only. The fast kernel returns no epoch, so T0 comes from a phase fold
# (bake-off glue fold_t0); snr, depth, rp and odd-even are box-fold estimates, not TLS values. SDE = (power-mean)/std of
# the BLS spectrum, on its own scale: STOP and CAND come from the inverted-flux calibration in GPU.md.
# Writes <lane>/res_v4_gpu/<tic>.json with settings='v4gpu'. Opens the GPU once per process, retrying on out-of-memory.
import sys,os,json,time,csv,subprocess,numpy as np
H=os.path.dirname(os.path.abspath(__file__));B=os.path.abspath(os.path.join(H,'..'));sys.path.insert(0,os.path.join(B,'edgefix'))
from v4_test import tdur,bin30,good_epochs
from fix_test import trim,unlone
STOP=float(os.environ.get('GPU_STOP','20.4'));CAND=float(os.environ.get('GPU_CAND','21.7'))  # inverted-flux calibration, GPU.md: p95 and max
LIVE_PY=os.path.join(B,'.venv','bin','python')  # lane run1 needs lightkurve/wotan, which live in the CPU venv (used read-only)
def star_row(lane,tic):
    for src in ('pool.csv','control.csv'):  # control.csv: lane control stars (hz-near has one)
        p=os.path.join(B,lane,src)
        if not os.path.exists(p):continue
        with open(p) as fh:
            for r in csv.DictReader(fh):
                if int(r['ID'])==tic:return float(r['rad']),float(r['mass']),float(r['Teff'])
    raise KeyError(tic)
def build_lc(lane,tic,R,M,cwd=None):
    # same call as run_star.py; cwd lets a test build into a scratch copy of the lane (hz-north / hz-near only)
    L=cwd or os.path.join(B,lane);npz=os.path.join(L,'lc',f'{tic}.npz')
    env={**os.environ,'LC_ONLY':'1','OMP_NUM_THREADS':'1','NUMBA_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MPLCONFIGDIR':os.path.join(L,'mplcache')}
    for _try in range(2):
        if os.path.exists(npz):break
        old=os.path.join(L,'res',f'{tic}.json')  # lane run1 skips stars that already have res/<tic>.json
        if os.path.exists(old):os.remove(old)
        a={'hz-north':['run1.py','pool.csv',str(tic),'res'],'hz-near':['run1.py',str(tic)],'hz-south':['run1.py',str(tic),str(R),str(M),os.path.join(L,'res'),'1']}[lane]
        subprocess.run(['nice','-n','19',LIVE_PY,'-I',os.path.join(L,a[0])]+a[1:],cwd=L,env=env,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,timeout=3*3600)
    return npz
_cb=None
def gpu():
    # cuvarbase grabs the GPU primary context at import (pycuda.autoprimaryctx). Unified memory is shared with other
    # residents, so that import can fail with out-of-memory: retry it, clearing the half-imported modules each time.
    global _cb
    while _cb is None:
        for i in range(360):
            try:
                import cuvarbase.bls as cb;_cb=cb;break
            except Exception as e:
                for k in [k for k in sys.modules if k.startswith(('cuvarbase','pycuda.autoprimaryctx'))]:del sys.modules[k]
                if i==359:raise
                time.sleep(10)
    return _cb
def bls(t,f,dy,R,M,pmin,pmax,spp):
    cb=gpu();import pycuda.driver as cuda
    fr,q0=cb.transit_autofreq(t,fmin=1/pmax,fmax=1/pmin,samples_per_peak=spp,rho=M/R**3,qmin_fac=0.5,qmax_fac=2.0)
    for i in range(120):
        try:
            fr,pw=cb.eebls_transit_gpu(t,f,dy,freqs=fr,qvals=q0,qmin_fac=0.5,qmax_fac=2.0,use_fast=True,ignore_negative_delta_sols=True)
            return 1/fr,np.asarray(pw)
        except (cuda.MemoryError,cuda.Error) as e:
            if 'out of memory' not in str(e).lower() or i==119:raise
            time.sleep(10)
def fold_t0(t,f,P,D):  # bake-off glue: epoch of the box at period P, width D (deepest mean in phase bins)
    nb=max(8,int(np.ceil(P/(D/4))));ph=((t-t.min())/P)%1;b=(ph*nb).astype(int)
    n=np.bincount(b,minlength=nb);s=np.bincount(b,f,minlength=nb);m=np.where(n>0,s/np.maximum(n,1),np.inf)
    w=max(1,int(round(D/P*nb)));mm=np.array([m.take(range(i,i+w),mode="wrap").mean() for i in range(nb)])
    i=int(np.argmin(mm));return t.min()+(i+w/2)/nb*P
def box_stats(t,f,P,T0,D):
    x=(t-T0)/P;ep=np.round(x);intr=np.abs(x-ep)*P<0.5*D;out=~intr
    base=np.median(f[out]);sd=1.4826*np.median(np.abs(f[out]-base));n=intr.sum()
    if n==0:return 0.0,0.0,0.0
    dep=base-f[intr].mean();snr=dep/sd*np.sqrt(n) if sd>0 else 0.0
    def d(k):
        m=intr&k;return (base-f[m].mean(),sd/np.sqrt(m.sum())) if m.sum() else (np.nan,np.nan)
    (a,sa),(b,sb)=d(ep%2==0),d(ep%2==1);oe=float(abs(a-b)/np.hypot(sa,sb)) if np.isfinite(a+b) and sa+sb>0 else 0.0
    return float(dep),float(snr),oe
def search(t,f,R,M,Teff,lane,stop=None):
    stop=STOP if stop is None else stop;span=t.max()-t.min()
    pmin=0.5 if lane=='hz-near' else 1.0;pmax=min(60.0,span/3);spp=3 if lane=='hz-near' else 2  # spp mirrors TLS oversampling
    f=f/np.median(f);Lum=R**2*(Teff/5772)**4;peaks=[]
    for it in range(4):
        dy=np.full_like(f,1.4826*np.median(np.abs(f-np.median(f))))
        per,pw=bls(t,f,dy,R,M,pmin,pmax,spp);s=(pw-pw.mean())/pw.std();j=int(np.argmax(s))
        P=float(per[j]);SDE=float(s[j]);D=tdur(P,R,M);T0=fold_t0(t,f,P,D);dep,snr,oe=box_stats(t,f,P,T0,D)
        g=good_epochs(t,P,T0,0.4*D);a=(M*(P/365.25)**2)**(1/3)
        peaks.append(dict(iter=it,P=P,T0=float(T0),SDE=SDE,snr=snr,depth_ppm=dep*1e6,rp_re=float(np.sqrt(max(dep,0)))*R*109.1,exp_dur_h=D*24,good_epochs=g,ok=g>=3,oe=oe,S_earth=Lum/a**2,nfreq=int(len(per))))
        if SDE<stop:break
        k=np.abs(((t-T0)/P+0.5)%1-0.5)*P>=0.5*1.5*D  # same mask as TLS transit_mask(t,P,1.5*D,T0)
        t,f=t[k],f[k]
    return peaks,pmin,pmax
def prep(npz,invert=False):
    d=np.load(npz);t,f=d['t'],d['f'];n=len(t)
    if invert:f=2-f
    if np.median(np.diff(t))<25/1440:t,f=bin30(t,f)
    t,f=trim(t,f);t,f=unlone(t,f);return t,f,n
def run_star(lane,tic,outdir=None,invert=False,npz=None,star=None,stop=None,build_cwd=None):
    t0=time.time();outdir=outdir or os.path.join(B,lane,'res_v4_gpu');os.makedirs(outdir,exist_ok=True)
    try:
        R,M,Teff=star or star_row(lane,tic)
        if npz is None:npz=build_lc(lane,tic,R,M,build_cwd)
        tb=time.time()-t0;t,f,n=prep(npz,invert);ts=time.time()
        peaks,pmin,pmax=search(t,f,R,M,Teff,lane,stop)
        D_=dict(tic=tic,lane=lane,settings='v4gpu',R=R,M=M,Teff=Teff,n_pts=n,span=float(t.max()-t.min()),pmin=pmin,pmax=pmax,
                stop=STOP if stop is None else stop,cand=CAND,invert=invert,peaks=peaks,lc_secs=round(tb,1),search_secs=round(time.time()-ts,1),secs=round(time.time()-t0))
    except Exception as e:
        D_=dict(tic=tic,lane=lane,settings='v4gpu',err=repr(e)[:300],secs=round(time.time()-t0))
    json.dump(D_,open(os.path.join(outdir,f'{tic}.json'),'w'));print(json.dumps(D_)[:400],flush=True);return D_
if __name__=='__main__':
    a=sys.argv[1:];out=None;inv='--invert' in a;stop=None
    if '--out' in a:i=a.index('--out');out=a[i+1];del a[i:i+2]
    if '--stop' in a:i=a.index('--stop');stop=float(a[i+1]);del a[i:i+2]
    a=[x for x in a if x!='--invert'];lane=a[0]
    for tic in a[1:]:run_star(lane,int(tic),out,inv,stop=stop)
