# One TRICERATOPS FPP/NFPP run with the pixel-ruled-out neighbours removed (evidence/06_trc_removal_match.csv, remove==True), via target.remove_star after calc_depths so dilution is unchanged. Light curve: the search's own npz (TESS-SPOC PDCSAP, flattened, 30-min bins; dilution-corrected).
# usage: trc_run.py RUN_INDEX
import sys, json, pickle, time, numpy as np, warnings; warnings.filterwarnings('ignore')
sys.path.insert(0,'/Users/advaitpaliwal/Scratch/planet-hunt/hunt-ffi/trc'); import gaia_patch
import triceratops.triceratops as tr
k=int(sys.argv[1]); P,T0=11.718154156303166,1329.529459184984
z=np.load('38699429.npz'); t,f=z['t'],z['f']
tg=pickle.load(open('trc/target.pkl','rb'))
ph=t-T0-np.round((t-T0)/P)*P; w=np.abs(ph)<0.25
edges=np.arange(-0.25,0.25+1e-9,10/1440); idx=np.digitize(ph[w],edges)
tb=np.array([ph[w][idx==i].mean() for i in np.unique(idx)]); fb=np.array([f[w][idx==i].mean() for i in np.unique(idx)])
dur=0.0625; oot=np.abs(tb)>dur; err=float(np.std(fb[oot])); depth=float(np.median(fb[oot])-fb[np.abs(tb)<dur/2*0.8].mean())
tg.calc_depths(tdepth=depth,all_ap_pixels=tg._aps,dilution_corrected=True)
import csv; drop=[int(r['tic']) for r in csv.DictReader(open('evidence/06_trc_removal_match.csv')) if r['remove']=='True']; tg.remove_star(drop)
np.random.seed(1000+k); t0=time.time()
tg.calc_probs(time=tb,flux_0=fb,flux_err_0=err,P_orb=P,contrast_curve_file=None,filt='TESS',N=1000000,parallel=True,exptime=30/1440,nsamples=20,verbose=0,dilution_corrected=True)
tg.probs.to_csv(f'trc_rm/probs_{k}.csv',index=False)
json.dump(dict(run=k,n_stars=len(tg.stars),removed=drop,n_stars_tdepth_gt0=int((tg.stars.tdepth>0).sum()),depth_ppm=depth*1e6,err_ppm=err*1e6,nbins=len(tb),FPP=float(tg.FPP),NFPP=float(tg.NFPP),sec=time.time()-t0),open(f'trc_rm/res_{k}.json','w'))
print(k,tg.FPP,tg.NFPP,round(time.time()-t0),flush=True)
