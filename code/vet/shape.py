# Step 1b: batman+emcee transit fit; density free vs TIC density; V vs U shape
import numpy as np, sys, json, batman, emcee
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
d=np.load(sys.argv[1]); t,f=d['t'],d['f']
P=11.718154156303166; T0=1329.529459184984
dt=(((t-T0)/P+0.5)%1-0.5)*P; m=np.abs(dt)<0.3; x,y=dt[m],f[m]
sig=1.4826*np.median(np.abs(y[np.abs(x)>0.06]-np.median(y[np.abs(x)>0.06])))
G=6.674e-8; rhoTIC=0.136941*1.989e33/(4/3*np.pi*(0.166264*6.957e10)**3)
pr=batman.TransitParams(); pr.t0=0; pr.per=P; pr.rp=0.1; pr.a=50; pr.inc=89.9; pr.ecc=0; pr.w=90; pr.u=[0.3,0.3]; pr.limb_dark='quadratic'
TM=batman.TransitModel(pr,x,supersample_factor=9,exp_time=30/1440)
def model(th):
    tc,rp,b,lrho,q1,q2,f0=th
    a=(10**lrho*G*(P*86400)**2/(3*np.pi))**(1/3)
    pr.t0=tc; pr.rp=rp; pr.a=a; pr.inc=np.degrees(np.arccos(b/a)); u1=2*np.sqrt(q1)*q2; u2=np.sqrt(q1)*(1-2*q2); pr.u=[u1,u2]
    return f0*TM.light_curve(pr)
def lp(th,fixrho=None):
    tc,rp,b,lrho,q1,q2,f0=th
    if not(-0.03<tc<0.03 and 0.01<rp<0.8 and 0<=b<1+rp and -1<lrho<2.5 and 0<q1<1 and 0<q2<1 and 0.98<f0<1.02): return -np.inf
    pri=0
    if fixrho is not None: pri=-0.5*((lrho-np.log10(fixrho))/0.05)**2   # ~12% density prior
    return pri-0.5*np.sum(((y-model(th))/sig)**2)
res={'rho_TIC_gcc':rhoTIC,'point_sigma_ppm':sig*1e6,'n_points':int(m.sum())}
fig,ax=plt.subplots(1,2,figsize=(13,4.5))
for k,(name,fx) in enumerate([('density_free',None),('density_TIC_prior',rhoTIC)]):
    nw=28; p0=np.array([0,0.15,0.3,np.log10(rhoTIC),0.4,0.3,1.0])+1e-3*np.random.default_rng(1).standard_normal((nw,7))
    S=emcee.EnsembleSampler(nw,7,lp,kwargs={'fixrho':fx}); S.run_mcmc(p0,5000,progress=False)
    ch=S.get_chain(discard=2000,thin=10,flat=True); lpf=S.get_log_prob(discard=2000,thin=10,flat=True)
    q=lambda v:[float(np.percentile(v,16)),float(np.median(v)),float(np.percentile(v,84))]
    rp=ch[:,1]; b=ch[:,2]; rho=10**ch[:,3]
    a=(rho*G*(P*86400)**2/(3*np.pi))**(1/3)
    T14=P/np.pi*np.arcsin(np.sqrt(np.clip((1+rp)**2-b**2,0,None))/a/np.sin(np.arccos(b/a)))*24
    best=ch[np.argmax(lpf)]
    res[name]=dict(rp_rs=q(rp),Rp_Re=q(rp*0.166264*109.08),b=q(b),rho_gcc=q(rho),a_rs=q(a),T14_h=q(T14),
       frac_grazing_b_gt_1_minus_rp=float(np.mean(b>1-rp)),chi2_best=float(-2*lpf.max()),best=best.tolist(),acc=float(np.mean(S.acceptance_fraction)))
    xs=np.linspace(-0.3,0.3,600); 
    be=np.arange(-0.3,0.3001,0.25/24); idx=np.digitize(x,be); bc=(be[1:]+be[:-1])/2
    ax[k].plot(x*24,(y-1)*1e6,'.',color='0.8',ms=2)
    ax[k].plot(bc*24,[np.mean(y[idx==j]-1)*1e6 if (idx==j).any() else np.nan for j in range(1,len(be))],'ko',ms=3)
    o=np.argsort(x); ax[k].plot(x[o]*24,(model(best)[o]-1)*1e6,'r-')
    ax[k].set_xlim(-5,5); ax[k].set_ylim(-40000,15000); ax[k].set_title(f'{name}: Rp {np.median(rp*0.166264*109.08):.2f} Re, b {np.median(b):.2f}, rho {np.median(rho):.1f}'); ax[k].set_xlabel('h')
res['dchi2_TICprior_minus_free']=res['density_TIC_prior']['chi2_best']-res['density_free']['chi2_best']
plt.tight_layout(); plt.savefig('plots/shape_fit.png',dpi=100)
json.dump(res,open('evidence/01b_shape_fit.json','w'),indent=1); print(json.dumps(res,indent=1))
