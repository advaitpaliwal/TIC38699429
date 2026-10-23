# usage: drive.py LANE NWORKERS ; resumable: skips stars with res_v4/TIC.json
import sys,os,subprocess,time,pandas as pd
from concurrent.futures import ThreadPoolExecutor
H=os.path.dirname(os.path.abspath(__file__));lane,nw=sys.argv[1],int(sys.argv[2]);nw=min(nw,int(os.environ.get("HUNT_MAX_WORKERS","2")))  # 3 lanes x 5 workers stalled the Mac Mini; 2 per lane keeps it responsive
L=os.path.join(H,'..',lane)
os.makedirs(os.path.join(L,'res_v4'),exist_ok=True);os.makedirs(os.path.join(L,'logs_v4'),exist_ok=True)
pool=pd.read_csv(os.path.join(L,'pool.csv'))
have=lambda i:os.path.exists(os.path.join(L,'lc',f'{i}.npz'))
todo=[int(i) for i in pool.ID if not os.path.exists(os.path.join(L,'res_v4',f'{i}.json'))]
todo.sort(key=lambda i:not have(i))  # stars with a saved light curve first (no download)
print(time.strftime('%F %T'),lane,'todo',len(todo),flush=True)
def run(i):
    with open(os.path.join(L,'logs_v4',f'{i}.log'),'w') as fh:
        try:subprocess.run(["nice","-n","19",sys.executable,'-I',os.path.join(H,'run_star.py'),lane,str(i)],stdout=fh,stderr=subprocess.STDOUT,timeout=5*3600)
        except subprocess.TimeoutExpired:open(os.path.join(L,'res_v4',f'{i}.json'),'w').write('{"tic":%d,"lane":"%s","settings":"v4","err":"timeout 5h"}'%(i,lane))
    print(time.strftime('%F %T'),i,'done',flush=True)
    subprocess.run([sys.executable,'-I',os.path.join(H,'summarize.py')],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
with ThreadPoolExecutor(nw) as ex:list(ex.map(run,todo))
print(time.strftime('%F %T'),lane,'ALL DONE',flush=True)
