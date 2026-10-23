# usage: v4_test.py NPZ R M PMIN PMAX THREADS OS
# v4: bin to 30 min; trim 0.375 d at gap edges; drop lone low points; star R/M +-25%;
# a peak counts only if >=3 epochs each have >=2 points within +-0.4 x central-transit duration (TLS's own duration is not trusted)
import sys,os,json,numpy as np
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from transitleastsquares import transitleastsquares,transit_mask
from fix_test import trim,unlone
def tdur(P,R,M):
    a=(M*(P/365.25)**2)**(1/3)*215.03;return P/np.pi*np.arcsin(min(1,R/a))
def bin30(t,f,w=30/1440):
    b=np.floor((t-t[0])/w).astype(int);n=np.bincount(b);k=n>0
    return np.bincount(b,t)[k]/n[k],np.bincount(b,f)[k]/n[k]
def good_epochs(t,P,T0,half):
    ep=np.round((t-T0)/P);m=np.abs(t-T0-ep*P)<half
    u,c=np.unique(ep[m],return_counts=True);return int((c>=2).sum())
if __name__=='__main__':
    npz,R,M,pmin,pmax=sys.argv[1],*map(float,sys.argv[2:6]);th,os_=int(sys.argv[6]),int(sys.argv[7])
    d=np.load(npz);t,f=d['t'],d['f']
    if np.median(np.diff(t))<25/1440:t,f=bin30(t,f)
    t,f=trim(t,f);t,f=unlone(t,f);pmax=min(pmax,(t.max()-t.min())/3)
    for it in range(5):
        r=transitleastsquares(t,f).power(period_min=pmin,period_max=pmax,R_star=R,M_star=M,R_star_min=0.75*R,R_star_max=1.25*R,M_star_min=0.75*M,M_star_max=1.25*M,oversampling_factor=os_,duration_grid_step=1.1,use_threads=th,show_progress_bar=False)
        P=float(r.period);D=tdur(P,R,M);g=good_epochs(t,P,r.T0,0.4*D)
        print(json.dumps(dict(iter=it,P=P,SDE=float(r.SDE),good_epochs=g,exp_h=D*24,depth_ppm=float((1-r.depth)*1e6),ok=g>=3)),flush=True)
        k=~transit_mask(t,P,1.5*D,r.T0);t,f=t[k],f[k]
