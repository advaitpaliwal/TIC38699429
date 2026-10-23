import sys, lightkurve as lk, warnings; warnings.filterwarnings('ignore')
s=int(sys.argv[1]); r=lk.search_tesscut('TIC 38699429',sector=s)
tpf=r.download(cutout_size=15,download_dir='cut'); print(s,tpf.path,tpf.shape)
