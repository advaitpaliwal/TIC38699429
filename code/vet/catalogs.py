import warnings; warnings.filterwarnings('ignore')
import astropy.units as u, json
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier
c=SkyCoord(66.1882793443072*u.deg,-62.0085242002372*u.deg)
V=Vizier(columns=['**'],row_limit=-1)
out={}
for cat,r in [('B/vsx/vsx',120),('I/358/veb',120),('I/358/vclassre',120),('I/355/gaiadr3',150),('J/AJ/156/102',60),('IV/39/tic82',150),('J/ApJS/258/16',60)]:
    try:
        res=V.query_region(c,radius=r*u.arcsec,catalog=cat)
        out[cat]=[t.to_pandas().astype(str).to_dict('records') for t in res]
        print(cat,[len(t) for t in res])
    except Exception as e: out[cat]=repr(e); print(cat,'ERR',e)
json.dump(out,open('evidence/02_catalogs_raw.json','w'),indent=0,default=str)
