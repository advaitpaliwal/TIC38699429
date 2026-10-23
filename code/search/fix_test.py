# usage: trim_test.py NPZ R M PMIN PMAX THREADS ; cuts 0.375 d at both edges of every data stretch (gaps > 0.5 d), then TLS with 3-transit rule
import sys,json,time,numpy as np
from transitleastsquares import transitleastsquares,transit_mask
def trim(t,f,cut=0.375,gap=0.5):
    b=np.r_[0,np.where(np.diff(t)>gap)[0]+1,len(t)];k=np.zeros(len(t),bool)
    for a,z in zip(b[:-1],b[1:]):
        s=t[a:z];k[a:z]=(s>s[0]+cut)&(s<s[-1]-cut)
    return t[k],f[k]
def unlone(t,f,lo=3,nb=1):
    r=f-np.median(f);s=1.4826*np.median(np.abs(r));k=np.ones(len(f),bool)
    k[1:-1]=~((r[1:-1]<-lo*s)&(np.abs(r[:-2])<nb*s)&(np.abs(r[2:])<nb*s));return t[k],f[k]
if __name__=='__main__':
    npz,R,M,pmin,pmax,th=sys.argv[1],float(sys.argv[2]),float(sys.argv[3]),float(sys.argv[4]),float(sys.argv[5]),int(sys.argv[6])
    d=np.load(npz);t,f=d['t'],d['f'];n0=len(t);t,f=trim(t,f);t,f=unlone(t,f);print('kept',len(t),'of',n0,flush=True)
    pmax=min(pmax,(t.max()-t.min())/3)
    for it in range(4):
        r=transitleastsquares(t,f).power(period_min=pmin,period_max=pmax,R_star=R,M_star=M,R_star_min=0.08,R_star_max=1.0,M_star_min=0.08,M_star_max=1.0,oversampling_factor=int(__import__('os').environ.get('OS','1')),duration_grid_step=1.1,use_threads=th,show_progress_bar=False)
        print(json.dumps(dict(iter=it,P=float(r.period),SDE=float(r.SDE),ntr=int(r.distinct_transit_count),dur_h=float(r.duration*24),depth_ppm=float((1-r.depth)*1e6),oe=float(r.odd_even_mismatch))),flush=True)
        k=~transit_mask(t,r.period,3*max(r.duration,1/24),r.T0);t,f=t[k],f[k]
