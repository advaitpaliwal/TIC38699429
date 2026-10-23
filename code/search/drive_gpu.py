# usage: drive_gpu.py [--prefetch N] LANE [LANE ...]   (run with the GPU venv: gpu_bakeoff/.venv/bin/python)
# One process does every GPU search (one GPU context, opened with out-of-memory retries); N CPU threads build light
# curves ahead of it with the lane run1 (LC_ONLY). Resumable: skips stars with res_v4/<tic>.json or res_v4_gpu/<tic>.json,
# re-checked just before each search. Works from the opposite end of each pool to the CPU drive.py (which takes stars
# with a saved light curve first, in pool order) so the two overlap as little as possible. Lanes are interleaved.
import sys,os,csv,time,itertools
from concurrent.futures import ThreadPoolExecutor
H=os.path.dirname(os.path.abspath(__file__));sys.path.insert(0,H)
import run_star_gpu as g
a=sys.argv[1:];N=4
if '--prefetch' in a:i=a.index('--prefetch');N=int(a[i+1]);del a[i:i+2]
lanes=a or ['hz-north','hz-south','hz-near'];B=g.B
done=lambda L,i:any(os.path.exists(os.path.join(B,L,d,f'{i}.json')) for d in ('res_v4','res_v4_gpu'))
Q=[]
for L in lanes:
    os.makedirs(os.path.join(B,L,'res_v4_gpu'),exist_ok=True)
    ids=[int(r['ID']) for r in csv.DictReader(open(os.path.join(B,L,'pool.csv')))]
    todo=[i for i in ids if not done(L,i)]
    have=lambda i:os.path.exists(os.path.join(B,L,'lc',f'{i}.npz'))
    cpu=sorted(todo,key=lambda i:not have(i))  # drive.py's order
    Q.append([(L,i) for i in reversed(cpu)]);print(time.strftime('%F %T'),L,'todo',len(todo),flush=True)
Q=[x for x in itertools.chain(*itertools.zip_longest(*Q)) if x]
def build(x):
    L,i=x
    if done(L,i):return
    try:R,M,_=g.star_row(L,i);g.build_lc(L,i,R,M)
    except Exception as e:print('lc build failed',L,i,repr(e)[:200],flush=True)
g.gpu()  # open the GPU first (retries while memory is short)
with ThreadPoolExecutor(N) as ex:
    fut={};k=0
    for j,x in enumerate(Q):
        while k<len(Q) and k<j+2*N:fut[k]=ex.submit(build,Q[k]);k+=1  # keep 2N light curves in flight ahead
        fut.pop(j).result();L,i=x
        if done(L,i):continue
        r=g.run_star(L,i);print(time.strftime('%F %T'),L,i,'done',r.get('secs'),flush=True)
print(time.strftime('%F %T'),'ALL DONE',flush=True)
