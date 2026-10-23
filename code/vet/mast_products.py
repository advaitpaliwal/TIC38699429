import warnings; warnings.filterwarnings('ignore')
from astroquery.mast import Observations
import json
obs=Observations.query_criteria(target_name=['38699429','0000000038699429'],obs_collection=['TESS','HLSP'])
print(len(obs)); 
pr=Observations.get_product_list(obs)
keep=pr['obs_collection','project','productFilename','dataURI','description']
keep.write('_sources/mast_products_38699429.csv',format='csv',overwrite=True)
for r in keep: print(r['project'],r['productFilename'])
