# TRICERATOPS target: all 31 TESS-SPOC sectors, SPOC optimal apertures from the lc FITS, Gaia via VizieR patch (hunt-ffi/trc/gaia_patch.py). No stars removed.
import sys, glob, pickle, re, numpy as np, warnings; warnings.filterwarnings('ignore')
sys.path.insert(0,'/Users/advaitpaliwal/Scratch/planet-hunt/hunt-ffi/trc')
import gaia_patch
import triceratops.triceratops as tr
from astropy.io import fits
fs=sorted(glob.glob('spoc/**/*tess-spoc*_lc.fits',recursive=True))
secs=np.array(sorted({int(re.search(r'-s00(\d\d)_',f).group(1)) for f in fs}))
tg=tr.target(ID=38699429,sectors=secs); aps=[]
for s in secs:
    f=[x for x in fs if f'-s00{s:02d}_' in x][0]
    with fits.open(f) as h:
        a=h['APERTURE'].data.astype(int); c0=h['APERTURE'].header['CRVAL1P']; r0=h['APERTURE'].header['CRVAL2P']
    rows,cols=np.nonzero(np.bitwise_and(a,2)); aps.append(np.column_stack([cols+c0,rows+r0]))
tg._aps=aps; pickle.dump(tg,open('trc/target.pkl','wb'))
print(len(secs),'sectors; stars in field',len(tg.stars)); print(tg.stars.head(12).to_string())
