# usage: run_star.py LANE TIC   (LANE = hz-north | hz-south | hz-near)
# v4 search: reuse lane/lc/TIC.npz (or build it with the lane's run1 in LC_ONLY mode), bin 30 min, trim 0.375 d at gap
# edges, drop lone low points, TLS with star R/M +-25%; a peak counts only if >=3 epochs each have >=2 points within
# +-0.4 x central-transit duration. Up to 4 peaks per star, stopping early once SDE < 9 (measured noise floor).
import sys,os,json,time,subprocess,numpy as np,pandas as pd
H=os.path.dirname(os.path.abspath(__file__));sys.path.insert(0,os.path.join(H,'..','edgefix'))
from transitleastsquares import transitleastsquares,transit_mask
from v4_test import tdur,bin30,good_epochs
from fix_test import trim,unlone
if __name__=='__main__':
    lane,tic=sys.argv[1],int(sys.argv[2]);L=os.path.join(H,'..',lane);t0=time.time()
    row=pd.read_csv(os.path.join(L,'pool.csv')).set_index('ID').loc[tic];R,M,Teff=float(row.rad),float(row.mass),float(row.Teff)
    out=os.path.join(L,'res_v4',f'{tic}.json');npz=os.path.join(L,'lc',f'{tic}.npz');py=sys.executable
    env={**os.environ,'LC_ONLY':'1','OMP_NUM_THREADS':'1','NUMBA_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MPLCONFIGDIR':os.path.join(L,'mplcache')}
    try:
        for _try in range(2):
            if os.path.exists(npz):break
            for old in (os.path.join(L,'res',f'{tic}.json'),):  # lane run1 skips stars that already have res/<tic>.json
                if os.path.exists(old):os.remove(old)
            a={'hz-north':['run1.py','pool.csv',str(tic),'res'],'hz-near':['run1.py',str(tic)],'hz-south':['run1.py',str(tic),str(R),str(M),os.path.join(L,'res'),'1']}[lane]
            subprocess.run([py,'-I',os.path.join(L,a[0])]+a[1:],cwd=L,env=env,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,timeout=3*3600)
        d=np.load(npz);t,f=d['t'],d['f']
        if np.median(np.diff(t))<25/1440:t,f=bin30(t,f)
        t,f=trim(t,f);t,f=unlone(t,f);span=t.max()-t.min()
        pmin=0.5 if lane=='hz-near' else 1.0;pmax=min(60.0,span/3);os_=3 if lane=='hz-near' else 2
        Lum=R**2*(Teff/5772)**4;peaks=[]
        for it in range(4):
            # TLS's results step can build a zero-sample template for sub-day peaks on tiny stars (ValueError,
            # zero-size array). Those peaks are far inside the habitable zone, so retry from 1 d.
            try:r=transitleastsquares(t,f).power(period_min=pmin,period_max=pmax,R_star=R,M_star=M,R_star_min=0.75*R,R_star_max=1.25*R,M_star_min=0.75*M,M_star_max=1.25*M,oversampling_factor=os_,duration_grid_step=1.1,use_threads=1,show_progress_bar=False)
            except ValueError as e:
                if 'zero-size' not in str(e):raise
                if peaks:break  # later passes: keep the peaks already found and stop
                if pmin>=1.0:raise
                pmin=1.0;r=transitleastsquares(t,f).power(period_min=pmin,period_max=pmax,R_star=R,M_star=M,R_star_min=0.75*R,R_star_max=1.25*R,M_star_min=0.75*M,M_star_max=1.25*M,oversampling_factor=os_,duration_grid_step=1.1,use_threads=1,show_progress_bar=False)
            P=float(r.period);D=tdur(P,R,M);g=good_epochs(t,P,r.T0,0.4*D);a=(M*(P/365.25)**2)**(1/3)
            peaks.append(dict(iter=it,P=P,T0=float(r.T0),SDE=float(r.SDE),snr=float(r.snr),depth_ppm=float((1-r.depth)*1e6),rp_re=float(r.rp_rs)*R*109.1,exp_dur_h=D*24,good_epochs=g,ok=g>=3,oe=float(r.odd_even_mismatch),S_earth=Lum/a**2))
            if r.SDE<9:break
            k=~transit_mask(t,P,1.5*D,r.T0);t,f=t[k],f[k]
        D_=dict(tic=tic,lane=lane,settings='v4',R=R,M=M,Teff=Teff,n_pts=int(len(d['t'])),span=float(span),peaks=peaks,secs=round(time.time()-t0))
    except Exception as e:
        D_=dict(tic=tic,lane=lane,settings='v4',err=repr(e)[:300],secs=round(time.time()-t0))
    json.dump(D_,open(out,'w'));print(json.dumps(D_)[:400])
