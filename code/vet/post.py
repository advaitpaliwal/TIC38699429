# per-neighbour summary across sectors: free-amplitude PRF fit at each Gaia star vs target-only fit; dchi2 summed over sectors.
# Tmag from Gaia G and BP-RP with the TIC-8 relation (Stassun+2019 eq.1); target uses its TIC Tmag 15.27. Flux zero point 15000 e/s at T=10.
import json,sys,pandas as pd,numpy as np
g=pd.read_csv('evidence/04_gaia_field.csv'); g['sepkey']=None
secs=[int(s) for s in sys.argv[1:]]; rows=[]
for s in secs:
    j=json.load(open(f'out/diff_s{s}.json'))
    st=pd.DataFrame(j['stars']); tgt=st[st.sep<0.5].iloc[0]
    st['dchi2']=st.chi2-tgt.chi2; st['sector']=s; st['ap3x3_depth_ppm']=j['ap3x3_depth_ppm']; st['n_events']=j['n_events']
    st['best_sep_arcsec']=j['best']['sep_from_target_arcsec']; rows.append(st)
df=pd.concat(rows); df['key']=df.Gmag.round(3)
br=dict(zip(g.Gmag.round(2),g['BP-RP']))
df['bprp']=df.Gmag.round(2).map(br)
c=df.bprp.fillna(1.0); df['Tmag']=df.Gmag-0.00522555*c**3+0.0891337*c**2-0.633923*c+0.0324473
df.loc[df.sep<0.5,'Tmag']=15.27
df['flux_es']=15000*10**(-0.4*(df.Tmag-10)); df['depth_needed']=df.amp/df.flux_es
df.to_csv('evidence/04_neighbour_fits_by_sector.csv',index=False)
agg=df.groupby('key').agg(source_float=('source','first'),sep=('sep','mean'),Gmag=('Gmag','first'),Tmag=('Tmag','first'),dchi2_sum=('dchi2','sum'),
    dchi2_by_sector=('dchi2',lambda v:[round(x,1) for x in v]),n_sec=('sector','count'),amp_mean=('amp','mean'),depth_needed_mean=('depth_needed','mean')).sort_values('sep')
agg['verdict']=np.where(agg.sep<0.5,'target',np.where(agg.dchi2_sum>25,'ruled out (dchi2>25)',np.where(agg.depth_needed_mean>1,'ruled out (needs >100%)','not excluded')))
agg.to_csv('evidence/04_neighbour_verdicts.csv'); print(agg.round(3).to_string())
print(df[df.sep<0.5][['sector','n_events','amp','amp_err','depth_needed','ap3x3_depth_ppm','best_sep_arcsec']].round(4).to_string())
