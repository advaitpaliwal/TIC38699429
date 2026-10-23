# Step 1: TLS refit on the search light curve, period window 5-30 d (covers P/2 and 2P), TIC stellar parameters, 1 thread
import numpy as np, json, sys
from transitleastsquares import transitleastsquares
d=np.load(sys.argv[1]); t,f=d['t'],d['f']
m=transitleastsquares(t,f)
r=m.power(period_min=11.6,period_max=11.85,R_star=0.166264,R_star_min=0.15,R_star_max=0.18,M_star=0.136941,M_star_min=0.12,M_star_max=0.16,use_threads=1,show_progress_bar=False,oversampling_factor=3,duration_grid_step=1.05)
keep=['period','period_uncertainty','T0','depth','duration','SDE','snr','rp_rs','odd_even_mismatch','transit_count','distinct_transit_count','empty_transit_count','depth_mean_even','depth_mean_odd']
out={k:(np.asarray(getattr(r,k)).tolist()) for k in keep}
out['depth_ppm']=(1-r.depth)*1e6; out['duration_h']=r.duration*24; out['Rp_Re_from_rp_rs']=r.rp_rs*0.166264*109.08
pp=np.array(r.periods); pw=np.array(r.power)
for lab,PP in [('P/2',r.period/2),('2P',r.period*2)]:
    j=np.argmin(abs(pp-PP)); out[f'SDE_at_{lab}']=float(pw[max(0,j-3):j+4].max()) if (pp.min()<PP<pp.max()) else None
json.dump(out,open('evidence/01c_tls_refit_narrow.json','w'),indent=1,default=float); print(json.dumps(out,indent=1,default=float))
