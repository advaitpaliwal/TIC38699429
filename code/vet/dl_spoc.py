import lightkurve as lk, warnings; warnings.filterwarnings('ignore')
r=lk.search_lightcurve('TIC 38699429',author=['TESS-SPOC','QLP','SPOC','TASOC','GSFC-ELEANOR-LITE','CDIPS','TGLC'])
print(r)
r=r[[a in('TESS-SPOC','QLP','SPOC','TGLC','GSFC-ELEANOR-LITE') for a in r.author]]
r.download_all(download_dir='spoc')
