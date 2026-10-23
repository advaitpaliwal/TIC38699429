# usage: v3_test.py NPZ R M PMIN PMAX THREADS OS
# fixes: trim 0.375 d at gap edges; drop lone low points; star R/M limited to +-25%;
# a peak counts only if >=3 transits contain data AND duration >= 0.4 x central-transit duration; otherwise mask and search again
import sys,json,numpy as np
from transitleastsquares import transitleastsquares,transit_mask
sys.path.insert(0,__import__("os").path.dirname(__import__("os").path.abspath(__file__)));from fix_test import trim,unlone
def tdur(P,R,M):  # central transit duration, days (circular, Rp<<R*)
    a=(M*(P/365.25)**2)**(1/3)*215.03;return P/np.pi*np.arcsin(min(1,R/a))
if __name__=='__main__':
    npz,R,M,pmin,pmax,th,os_=sys.argv[1],*map(float,sys.argv[2:6]),int(sys.argv[6]),int(sys.argv[7])
    d=np.load(npz);t,f=d['t'],d['f'];t,f=trim(t,f);t,f=unlone(t,f)
    pmax=min(pmax,(t.max()-t.min())/3)
    for it in range(6):
        r=transitleastsquares(t,f).power(period_min=pmin,period_max=pmax,R_star=R,M_star=M,R_star_min=0.75*R,R_star_max=1.25*R,M_star_min=0.75*M,M_star_max=1.25*M,oversampling_factor=os_,duration_grid_step=1.1,use_threads=th,show_progress_bar=False)
        P=float(r.period);dur=float(r.duration);ep=np.round((t-r.T0)/P);n=len(np.unique(ep[np.abs(t-r.T0-ep*P)<0.5*dur]))
        ok=bool(n>=3 and dur>=0.4*tdur(P,R,M))
        print(json.dumps(dict(iter=it,P=P,SDE=float(r.SDE),ntr_data=n,dur_h=dur*24,exp_h=tdur(P,R,M)*24,depth_ppm=float((1-r.depth)*1e6),ok=ok)),flush=True)
        k=~transit_mask(t,P,3*max(dur,1/24),r.T0);t,f=t[k],f[k]
