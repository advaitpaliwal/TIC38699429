# Step 1a: box-fit depth per sector/year, odd/even, secondary scan at all phases, 2P and P/2
import numpy as np, sys, json
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
d=np.load(sys.argv[1]); t,f,s=d['t'],d['f'],d['s']
P=11.718124156303166; T0=1329.532459184984
def dtof(P,T0): return (((t-T0)/P+0.5)%1-0.5)*P
def box(P,T0,w,sel=None):
    dt=dtof(P,T0); i=np.abs(dt)<w/2; o=(np.abs(dt)>w*1.5)&(np.abs(dt)<0.5)
    if sel is not None: i&=sel; o&=sel
    if i.sum()<2: return None
    return dict(depth_ppm=(np.median(f[o])-np.mean(f[i]))*1e6, err_ppm=np.std(f[o])/np.sqrt(i.sum())*1e6, n_in=int(i.sum()))
out={'sectors_in_npz':sorted(set(s.tolist()))}
# refine T0/P/w by grid
best=None
for dP in np.linspace(-2e-4,2e-4,41):
  for dT in np.linspace(-0.03,0.03,61):
    for w in [0.04,0.05,0.06,0.07,0.08]:
      r=box(P+dP,T0+dT,w); sn=r['depth_ppm']/r['err_ppm']
      if best is None or sn>best[0]: best=(sn,P+dP,T0+dT,w,r)
_,P,T0,w,r=best; out['refined']=dict(P=P,T0=T0,w_days=w,**r,snr=best[0])
# per-sector & per-year
out['per_sector']={int(k):box(P,T0,w,s==k) for k in sorted(set(s))}
yr=lambda k:1 if k<=13 else 3 if k<=39 else 5 if k<=69 else 7
ys=np.array([yr(k) for k in s])
out['per_year']={int(y):box(P,T0,w,ys==y) for y in sorted(set(ys))}
# epochs list
ep=np.round((t-T0)/P).astype(int); dt=dtof(P,T0); ine=np.abs(dt)<w/2
out['epochs_with_intransit_points']=int(len(set(ep[ine])))
pe={}
for e in sorted(set(ep[ine])):
    m=(ep==e); r=box(P,T0,w,m)
    if r: pe[int(e)]=dict(**r,sector=int(np.median(s[m&ine])))
out['per_epoch']=pe
# odd/even
o=box(P,T0,w,ep%2==1); e_=box(P,T0,w,ep%2==0)
out['odd']=o; out['even']=e_; out['oddeven_sigma']=abs(o['depth_ppm']-e_['depth_ppm'])/np.hypot(o['err_ppm'],e_['err_ppm'])
# 2P: dip at phase 0.5 of 2P?  P/2: dip at phase 0.5 of P
out['at_2P_phase0']=box(2*P,T0,w); out['at_2P_phase0.5']=box(2*P,T0+P,w)
out['P_over_2_phase0.5_of_P']=box(P,T0+P/2,w)
# secondary scan over all phases
ph=np.arange(0.02,0.98,0.0025); sec=[]
for p in ph:
    r=box(P,T0+p*P,w); sec.append((p,r['depth_ppm'],r['err_ppm']))
sec=np.array(sec); k=np.argmax(sec[:,1]/sec[:,2])
out['secondary_best']=dict(phase=sec[k,0],depth_ppm=sec[k,1],err_ppm=sec[k,2],snr=sec[k,1]/sec[k,2])
out['secondary_scan_rms_snr']=float(np.std(sec[:,1]/sec[:,2]))
out['secondary_at_0.5']=box(P,T0+0.5*P,w)
# plots
fig,ax=plt.subplots(2,2,figsize=(14,8))
x=dtof(P,T0)*24; m=np.abs(x)<8
ax[0,0].plot(x[m],(f[m]-1)*1e6,'.',color='0.75',ms=2)
be=np.arange(-8,8.01,0.25); idx=np.digitize(x[m],be); bc=(be[1:]+be[:-1])/2
bm=[np.mean(f[m][idx==j]-1)*1e6 if (idx==j).any() else np.nan for j in range(1,len(be))]
ax[0,0].plot(bc,bm,'ko-',ms=3); ax[0,0].set_title(f'P={P:.5f} fold, depth {out["refined"]["depth_ppm"]:.0f} ppm'); ax[0,0].set_xlabel('h')
out['profile']=list(zip(bc.tolist(),[round(v) if v==v else None for v in bm]))
ph_all=(((t-T0)/P)%1); be=np.linspace(0,1,401); idx=np.digitize(ph_all,be)
ax[0,1].plot((be[1:]+be[:-1])/2,[np.mean(f[idx==j]-1)*1e6 for j in range(1,len(be))],'k.-',ms=2); ax[0,1].set_title('full phase (42-min bins)'); ax[0,1].set_ylim(-3000,1500)
ax[1,0].plot(sec[:,0],sec[:,1]/sec[:,2],'k-'); ax[1,0].set_title('secondary box SNR vs phase'); ax[1,0].axhline(3,color='r',ls=':')
for lab,sel,c in [('odd',ep%2==1,'C0'),('even',ep%2==0,'C3')]:
    mm=m&sel; idx=np.digitize(x[mm],np.arange(-8,8.01,0.25))
    ax[1,1].plot(bc,[np.mean(f[mm][idx==j]-1)*1e6 if (idx==j).any() else np.nan for j in range(1,len(bc)+1)],'o-',color=c,label=lab,ms=3)
ax[1,1].legend(); ax[1,1].set_title('odd vs even')
plt.tight_layout(); plt.savefig('plots/fold_secondary_oddeven.png',dpi=100)
json.dump(out,open('evidence/01a_fold_secondary_oddeven.json','w'),indent=1,default=float)
print(json.dumps({k:v for k,v in out.items() if k not in('per_epoch','profile','per_sector')},indent=1,default=float))
print('per_sector',{k:(round(v['depth_ppm']),round(v['err_ppm']),v['n_in']) if v else None for k,v in out['per_sector'].items()})
print('per_epoch',{k:(round(v['depth_ppm']),round(v['err_ppm']),v['n_in'],v['sector']) for k,v in pe.items()})
print('profile',out['profile'])
