# Glue between the v4 search and LEO-vetter 1.2.0 (https://github.com/mkunimoto/LEO-Vetter, used as published).
# Light curve: same prep as run_star.py (30-min bin, 0.375 d gap-edge trim, lone-low-point drop).
# Earlier peaks of the same star are masked (1.5 x expected duration) exactly as the search did before finding peak k.
import os,sys,warnings,numpy as np
H=os.path.dirname(os.path.abspath(__file__));sys.path.insert(0,os.path.join(H,'..','..','edgefix'))
from v4_test import bin30,tdur
from fix_test import trim,unlone
from transitleastsquares import transit_mask
from leo_vetter.main import TCELightCurve
from leo_vetter.stellar import quadratic_ldc
from leo_vetter import thresholds as T
FA_TESTS=[T.weak,T.invalid_transits,T.bad_shape,T.non_unique,T.chases,T.dmm,T.single_event,T.bad_fit,T.sinusoidal,T.unphysical_duration,T.asymmetric,T.chi,T.data_gapped]
FP_TESTS=[T.odd_even,T.vshaped,T.large,T.secondary]  # offset test needs pixel vetting (not run)
def prep(npz):
    d=np.load(npz);t,f=d['t'],d['f']
    if np.median(np.diff(t))<25/1440:t,f=bin30(t,f)
    t,f=trim(t,f);t,f=unlone(t,f);return t,f
def errors(t,f,gap=0.5):
    # no per-point errors are stored in the npz files: use the robust scatter (1.4826 x MAD) of each data stretch
    e=np.empty_like(f);b=np.r_[0,np.where(np.diff(t)>gap)[0]+1,len(t)]
    for a,z in zip(b[:-1],b[1:]):
        x=f[a:z];e[a:z]=1.4826*np.median(np.abs(x-np.median(x)))
    e[~(e>0)]=np.nanmedian(e[e>0]);return e
def star_dict(tic,R,M,Teff,logg=None):
    if logg is None or not np.isfinite(logg):logg=4.438+np.log10(M/R**2)
    u1,u2=quadratic_ldc(float(Teff),float(logg))
    # catalogue errors are not in pool.csv; assume 10% R, M and 100 K (they only set derived-parameter error bars)
    return dict(tic=int(tic),rad=R,mass=M,Teff=Teff,logg=logg,rho=M/R**3,e_rad=0.1*R,e_mass=0.1*M,e_Teff=100.0,e_rho=0.3*M/R**3,u1=float(u1),u2=float(u2))
def mask_prev(t,f,prev,R,M):
    for P,T0 in prev:
        k=~transit_mask(t,P,1.5*tdur(P,R,M),T0);t,f=t[k],f[k]
    return t,f
def vet(tic,t,f,P,T0,dur,star):
    e=errors(t,f)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        tlc=TCELightCurve(tic,t,f,f,e,P,T0,dur)
        tlc.compute_flux_metrics(star,verbose=False)
        m=tlc.metrics
        fails=[]
        for test in FA_TESTS+FP_TESTS:
            try:flag,msg=test(m,T._default_thresholds)
            except Exception as x:flag,msg=True,f'{test.__name__}: error {x!r}'
            if bool(flag):fails.append(f'{test.__name__} ({msg})')
    keep={k:(float(v) if isinstance(v,(int,float,np.floating,np.integer)) else str(v)) for k,v in m.items()}
    return dict(verdict='PC' if not fails else 'FAIL',fails=fails,metrics=keep)
