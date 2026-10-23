# usage: control.py planets | artifacts
# planets: TLS (v4 settings, narrow period window around the known period) on the v4-prepped light curve gives P, T0, duration;
#          planets found earlier in the v4 search order are masked first, as the search did. Each planet is vetted twice:
#          with TLS's duration and with the v4 expected central duration (the one gate.py must use).
# artifacts: the 23 v4/vet1 stars, through gate.py's own code path on their matching res_v4 peak.
import os,sys,json,numpy as np,pandas as pd
H=os.path.dirname(os.path.abspath(__file__));sys.path.insert(0,H)
from leo_glue import prep,star_dict,mask_prev,vet,tdur
from transitleastsquares import transitleastsquares
B=os.path.abspath(os.path.join(H,'..','..'))
PL=[  # name, tic, npz, R, M, Teff, [(known P) in v4 search order]; target = last entry
 ('TOI-700 c',150428135,'hz-south/control/150428135.npz',0.419457,0.413261,3482.,[16.0512]),
 ('TOI-700 d',150428135,'hz-south/control/150428135.npz',0.419457,0.413261,3482.,[16.0512,9.9773,37.4246]),
 ('TOI-2094 b',356016119,'hz-north/lc/356016119.npz',0.376775,0.365046,3457.,[18.7933]),
 ('TOI-2095 b',235678745,'hz-north/lc/235678745.npz',0.463205,0.460544,3746.,[17.6648]),
]
def tls(t,f,P0,R,M):
    r=transitleastsquares(t,f).power(period_min=P0*0.999,period_max=P0*1.001,R_star=R,M_star=M,R_star_min=0.75*R,R_star_max=1.25*R,M_star_min=0.75*M,M_star_max=1.25*M,oversampling_factor=5,duration_grid_step=1.1,use_threads=1,show_progress_bar=False)
    return float(r.period),float(r.T0),float(r.duration),float(r.SDE),float((1-r.depth)*1e6)
if sys.argv[1]=='planets':
    out=[]
    for name,tic,npz,R,M,Te,seq in PL:
        t,f=prep(os.path.join(B,npz));prev=[]
        for P0 in seq[:-1]:
            P,T0,*_=tls(t,f,P0,R,M);prev.append((P,T0));t,f=mask_prev(t,f,[(P,T0)],R,M)
        P,T0,D,S,dep=tls(t,f,seq[-1],R,M);st=star_dict(tic,R,M,Te)
        for lab,dur in [('tls_dur',D),('v4_exp_dur',tdur(P,R,M))]:
            o=vet(tic,t,f,P,T0,dur,st);m=o['metrics']
            r=dict(name=name,tic=tic,dur_kind=lab,P=P,T0=T0,dur_h=dur*24,tls_SDE=S,tls_depth_ppm=dep,leo=o['verdict'],failing_tests='; '.join(o['fails']),**{k:m.get(k) for k in ['MES','new_MES','N_transit','new_N_transit','max_SES','SHP','DMM','CHI','mean_chases','sine_sig','Rp','sig_dep','transit_sig_dep','N_gap_2.0']})
            print(json.dumps(r),flush=True);out.append(r)
            json.dump(o,open(os.path.join(H,'cache',f'control_{name.replace(" ","_")}_{lab}.json'),'w'))
    pd.DataFrame(out).to_csv(os.path.join(H,'control_planets.csv'),index=False)
else:
    import gate
    s=pd.read_csv(os.path.join(B,'v4','vet1','stars.csv'));out=[]
    for _,x in s.iterrows():
        res=json.load(open(os.path.join(B,x.lane,'res_v4',f'{x.tic}.json')))
        i=int(np.argmin([abs(p['P']-x.P)+abs(p['T0']-x.T0)/1e3 for p in res['peaks']]))
        pk=res['peaks'][i];res2=dict(res,peaks=[dict(p,ok=(j==i) or p['ok'],SDE=p['SDE']) for j,p in enumerate(res['peaks'][:i+1])])
        rows=gate.vet_star(x.lane,int(x.tic),gate.pool(x.lane).loc[int(x.tic)],res2)
        r=[q for q in rows if q['iter']==pk['iter']][0];r['P_vet1']=x.P;print(json.dumps(r),flush=True);out.append(r)
    pd.DataFrame(out).to_csv(os.path.join(H,'control_artifacts.csv'),index=False)
