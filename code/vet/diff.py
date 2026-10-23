# In/out-of-transit difference images from TESScut FFI stamps + PRF fits at each Gaia star.
# star-by-star single-source fits
# free-position grid
# 1-sigma (dchi2<2.3) region extent
# aperture depth on 3x3 around target pixel
import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__))); from common import *
def main():
    import sys, json, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    from PRF import TESS_PRF
    s=int(sys.argv[1]); tpf=tpf_for(s)
    F=tpf.flux.value.astype(float); t=tpf.time.value
    q=(tpf.quality==0)&np.isfinite(t)&np.all(np.isfinite(F.reshape(len(t),-1)),1); t,F=t[q],F[q]
    ny,nx=F.shape[1:]
    cad=np.median(np.diff(t))
    def event_img(tc):
        dt=t-tc; a=np.abs(dt)<0.4*DUR+cad/2; b=(np.abs(dt)>0.75*DUR)&(np.abs(dt)<0.3)
        if a.sum()<max(2,int(0.4*DUR/cad)) or b.sum()<10 or (dt[b]<0).sum()<3 or (dt[b]>0).sum()<3: return None
        X=np.vstack([np.ones(b.sum()),dt[b]]).T; Y=F[b].reshape(b.sum(),-1)
        beta,*_=np.linalg.lstsq(X,Y,rcond=None)
        # robust: sigma-clip out-of-window residuals once
        r=Y-X@beta; s_=1.4826*np.median(np.abs(r),0); ok=np.all(np.abs(r)<5*s_+1e-9,1)
        if ok.sum()>8: beta,*_=np.linalg.lstsq(X[ok],Y[ok],rcond=None)
        Xa=np.vstack([np.ones(a.sum()),dt[a]]).T
        return ((Xa@beta)-F[a].reshape(a.sum(),-1)).mean(0).reshape(ny,nx)
    def stack(off=0.0):
        n0=int(np.floor((t.min()-T0)/P)); n1=int(np.ceil((t.max()-T0)/P)); ims=[]
        for n in range(n0,n1+1):
            im=event_img(T0+n*P+off)
            if im is not None: ims.append(im)
        return (np.mean(ims,0) if ims else None), len(ims)
    D,ne=stack()
    rng=np.random.default_rng(s); nulls=[]
    for o in rng.uniform(0.2,P-0.2,120):
        d,_=stack(o)
        if d is not None: nulls.append(d)
    nulls=np.array(nulls); nulls=nulls-np.median(nulls,axis=(1,2))[:,None,None]; D=D-np.median(D); sig=1.4826*np.median(np.abs(nulls-np.median(nulls,0)),0)
    oot=np.median(F,0)
    g=gaia_field(tpf); g=g[(g.Gmag<20.5)&(g.sep<150)].reset_index(drop=True)
    prf=TESS_PRF(tpf.camera,tpf.ccd,s,tpf.column+nx/2,tpf.row+ny/2,localdatadir=None)
    def model(x,y): return prf.locate(x,y,(ny,nx))
    w=1/sig**2
    def fit_at(x,y):
        m=model(x,y); A=np.vstack([m.ravel(),np.ones(m.size)]).T; W=w.ravel()
        beta=np.linalg.solve(A.T@(A*W[:,None]),A.T@(D.ravel()*W)); r=D.ravel()-A@beta
        cov=np.linalg.inv(A.T@(A*W[:,None]))
        return float((r**2*W).sum()), float(beta[0]), float(np.sqrt(cov[0,0]))
    rows=[]
    for i,r in g.iterrows():
        if not (-1<r.x<nx and -1<r.y<ny): continue
        c2,a,ae=fit_at(r.x,r.y)
        T=r.Gmag-0.43  # rough G->T offset (inferred; replaced by TIC Tmag later where matched)
        rows.append(dict(source=int(r.Source),Gmag=r.Gmag,sep=r.sep,pa=r.pa,x=r.x,y=r.y,chi2=c2,amp=a,amp_err=ae,approx_flux_es=15000*10**(-0.4*((15.27 if int(r.Source)==TARGET_GAIA else r.Gmag-0.43)-10)),depth_needed=a/(15000*10**(-0.4*((15.27 if int(r.Source)==TARGET_GAIA else r.Gmag-0.43)-10)))))
    xs=np.arange(3,12.01,0.1); ys=np.arange(3,12.01,0.1); C=np.full((len(ys),len(xs)),np.nan); Am=np.zeros_like(C)
    for j,yv in enumerate(ys):
        for i,xv in enumerate(xs):
            C[j,i],Am[j,i],_=fit_at(xv,yv)
    j,i=np.unravel_index(np.nanargmin(np.where(Am>0,C,np.inf)),C.shape)
    best=dict(x=float(xs[i]),y=float(ys[j]),chi2=float(C[j,i]),amp=float(Am[j,i]))
    reg=(C-C[j,i])<2.30; yy,xx=np.meshgrid(ys,xs,indexing='ij')
    best['x_lo'],best['x_hi']=float(xx[reg].min()),float(xx[reg].max()); best['y_lo'],best['y_hi']=float(yy[reg].min()),float(yy[reg].max())
    ra,de=tpf.wcs.all_pix2world([[best['x'],best['y']]],0)[0]
    from astropy.coordinates import SkyCoord; import astropy.units as u
    tg=g.iloc[0]; best['sep_from_target_arcsec']=float(SkyCoord(ra*u.deg,de*u.deg).separation(SkyCoord(tg.ra_ep*u.deg,tg.de_ep*u.deg)).arcsec)
    tx,ty=tg.x,tg.y
    cx,cy=int(round(tx)),int(round(ty)); ap=np.zeros((ny,nx),bool); ap[cy-1:cy+2,cx-1:cx+2]=True
    apd=D[ap].sum()/oot[ap].sum(); _n=nulls[:,ap].sum(1); ape=1.4826*np.median(np.abs(_n-np.median(_n)))/oot[ap].sum()
    np.savez(f'out/diff_s{s}.npz',D=D,sig=sig,oot=oot,C=C,xs=xs,ys=ys,nulls=nulls)
    json.dump(dict(sector=s,n_events=ne,cadence_min=cad*1440,ap3x3_depth_ppm=apd*1e6,ap3x3_err_ppm=ape*1e6,target_xy=[tx,ty],best=best,stars=rows),open(f'out/diff_s{s}.json','w'),indent=1,default=float)
    fig,ax=plt.subplots(1,3,figsize=(15,5))
    for a_,im,ti in zip(ax,[oot,D,D/sig],['out-of-transit median (e/s)','in-vs-out difference (e/s, + = dimming)','difference / null std']):
        h=a_.imshow(im,origin='lower',cmap='viridis' if a_ is ax[0] else 'RdBu_r',vmin=None if a_ is ax[0] else -np.nanmax(np.abs(im)),vmax=None if a_ is ax[0] else np.nanmax(np.abs(im)))
        plt.colorbar(h,ax=a_,fraction=.046); a_.set_title(f'S{s} {ti}',fontsize=9)
        b=g[(g.Gmag<17)]; a_.scatter(b.x,b.y,s=np.clip(150/2**((b.Gmag-10)/2),4,150),facecolor='none',edgecolor='k')
        a_.plot(tx,ty,'+',c='lime',ms=14,mew=2); a_.plot(best['x'],best['y'],'x',c='m',ms=12,mew=2)
        a_.set_xlim(-.5,nx-.5); a_.set_ylim(-.5,ny-.5)
    fig.tight_layout(); fig.savefig(f'out/diff_s{s}.png',dpi=80)
    print(s,ne,'ap depth',round(apd*1e6),'+-',round(ape*1e6),'best',{k:round(v,2) for k,v in best.items()})
    
if __name__=='__main__':
    main()