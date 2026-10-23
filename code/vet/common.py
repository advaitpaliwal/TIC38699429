import warnings; warnings.filterwarnings('ignore')
import numpy as np, glob, lightkurve as lk, astropy.units as u, os, pandas as pd
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier
P,T0,DUR=11.718154156303166,1329.529459184984,1.5/24
RA0,DE0=66.1882793443072,-62.0085242002372
TARGET_GAIA=4676612253064169344
def tpf_for(s): return lk.read(glob.glob(f'cut/tesscut/tess-s{s:04d}-*_15x15_astrocut.fits')[0])
def gaia_field(tpf, rad=150):
    fn='evidence/04_gaia_field.csv'
    if not os.path.exists(fn):
        v=Vizier(columns=['Source','RA_ICRS','DE_ICRS','pmRA','pmDE','Gmag','BP-RP','Plx','RUWE'],row_limit=-1)
        v.query_region(SkyCoord(RA0*u.deg,DE0*u.deg),radius=rad*u.arcsec,catalog='I/355/gaiadr3')[0].to_pandas().to_csv(fn,index=False)
    t=pd.read_csv(fn)
    dy=(tpf.time[0].jd-2457389.0)/365.25
    ra=(t.RA_ICRS+np.nan_to_num(t.pmRA)*dy/3.6e6/np.cos(np.radians(t.DE_ICRS))).values; de=(t.DE_ICRS+np.nan_to_num(t.pmDE)*dy/3.6e6).values
    x,y=tpf.wcs.all_world2pix(np.vstack([ra,de]).T,0).T
    c=SkyCoord(ra*u.deg,de*u.deg); i0=int(np.where(t.Source.values==TARGET_GAIA)[0][0])
    t=t.assign(ra_ep=ra,de_ep=de,x=x,y=y,sep=c.separation(c[i0]).arcsec,pa=c[i0].position_angle(c).deg)
    return t.sort_values('sep').reset_index(drop=True)
