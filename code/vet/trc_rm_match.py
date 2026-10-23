# Cross-match the pixel-fit neighbour verdicts (Gaia) to the TRICERATOPS star table (TIC) by position + G mag.
# Removal rule: dchi2 > 25 in every usable sector.
import csv, ast, numpy as np
g={}
for r in csv.DictReader(open('evidence/04_gaia_field.csv')): g[round(float(r['Gmag']),3)]=r
t=list(csv.DictReader(open('evidence/05_trc_stars.csv')))
out=[]
for r in csv.DictReader(open('evidence/04_neighbour_verdicts.csv')):
    if r['verdict']=='target': continue
    d=ast.literal_eval(r['dchi2_by_sector']); G=float(r['Gmag']); gr=g[round(G,3)]
    ra,de=float(gr['RA_ICRS']),float(gr['DE_ICRS'])
    best=min(t,key=lambda s:np.hypot((float(s['ra'])-ra)*np.cos(np.radians(de)),float(s['dec'])-de))
    sep=3600*np.hypot((float(best['ra'])-ra)*np.cos(np.radians(de)),float(best['dec'])-de)
    gm=float(best['GAIAmag']) if best['GAIAmag'] else np.nan
    ok=sep<3 and abs(gm-G)<0.3
    rule=min(d)>25
    out.append(dict(gaia=gr['Source'],G=G,sep_target=round(float(r['sep']),1),dchi2_by_sector=d,dchi2_min=min(d),dchi2_sum=round(float(r['dchi2_sum']),1),
      tic=best['ID'] if ok else '',match_sep=round(sep,2),tic_G=gm,tic_tdepth=float(best['tdepth']) if ok else '',remove=bool(rule and ok)))
w=csv.DictWriter(open('evidence/06_trc_removal_match.csv','w'),fieldnames=list(out[0])); w.writeheader(); w.writerows(out)
for o in out: print(o)
