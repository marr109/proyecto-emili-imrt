import csv,re,glob,os,numpy as np,json
from scipy.stats import wilcoxon,mannwhitneyu,kruskal
B='/Users/marrojasr/Documents/code/emili_imrt/experiments/new'
def ls_runs():
    R=[]
    for t in glob.glob(B+'/local_search/nangshift*/*/data/*/seed*/trajectory.csv'):
        p=t.split('/'); step=p[-6].replace('nangshift',''); var=p[-5]; strat=p[-3]; seed=int(p[-2][4:])
        cat='unrestricted' if var.startswith('unres') else 'restricted'
        order='circular' if 'circular' in var else ('random' if 'rand' in var else 'fixed')
        rows=list(csv.DictReader(open(t)))
        tot=len(rows); nc=sum(r['cached']=='false' for r in rows)
        m=re.search(r'FMO objective f\*\s*:\s*([\d.]+)',open(os.path.dirname(t)+'/run.log').read())
        R.append(dict(fam='LS',step=step,cat=cat,strat=strat,order=order,seed=seed,f=float(m.group(1)),tot=tot,nc=nc,start=rows[0]['angles_deg']))
    return R
LS=ls_runs()
def pr(fam):
    R=[]
    for r in csv.DictReader(open(f'{B}/{fam}/comparisons/per_run_metrics.csv')):
        R.append(dict(fam=r['family'],var=r['variant'],step=r['step'],cat=r['catalog'],strat=r['strategy'],order=r['order'].lower(),seed=int(r['seed']),f=float(r['final_fstar']),tot=int(r['total_evals']),nc=int(r['noncached_evals'])))
    return R
SA=pr('sa'); ILS=[r for r in pr('ils') if r['fam']=='ILS']
print(len(LS),len(SA),len(ILS), set(r['order'] for r in ILS))
def key(r): return (r['step'],r['cat'],r['seed'])
def paired(A,Bb,metric):
    a={key(r):r[metric] for r in A}; b={key(r):r[metric] for r in Bb}
    k=sorted(set(a)&set(b)); x=np.array([a[i] for i in k]); y=np.array([b[i] for i in k])
    d=y-x; rel=np.median(d/x)*100
    p=wilcoxon(x,y).pvalue if np.any(d!=0) else 1.0
    return len(k),np.median(x),np.median(y),rel,p,(d<0).sum(),(d>0).sum()
def unp(A,Bb,metric):
    x=np.array([r[metric] for r in A]);y=np.array([r[metric] for r in Bb])
    rel=(np.median(y)-np.median(x))/np.median(x)*100
    return len(x),len(y),np.median(x),np.median(y),rel,mannwhitneyu(x,y).pvalue
F=lambda R,**kw:[r for r in R if all(r[k]==v for k,v in kw.items())]
out={}
lf=F(LS,strat='first'); lb=F(LS,strat='best')
# start identity
st={}
for r in LS: st.setdefault(key(r),set()).add(r['start'])
print('LS same start',sum(len(v)==1 for v in st.values()),len(st))
for o in ['circular','random']:
  for m in ['f','nc','tot']:
    out[f'LS_first_{o}_vs_fixed_{m}']=paired(F(lf,order='fixed'),F(lf,order=o),m)
# per condition circular vs fixed nc
for s in ['5','10']:
  for c in ['restricted','unrestricted']:
    for o in ['circular','random']:
      for m in ['nc','tot']:
        out[f'LS_cond_{s}_{c}_{o}_{m}']=paired(F(lf,order='fixed',step=s,cat=c),F(lf,order=o,step=s,cat=c),m)
# best invariance
same=0;n=0
for r in F(lb,order='fixed'):
  for o in ['circular','random']:
    q=[x for x in lb if x['order']==o and key(x)==key(r)][0]; n+=1; same+= abs(q['f']-r['f'])<0.01
out['best_inv']=(same,n)
for m in ['f','nc','tot']: out[f'first_vs_best_{m}']=paired(F(lb,order='fixed'),F(lf,order='fixed'),m)
# catalog
for nm,R in [('LSfirst',F(lf,order='fixed')),('LSbest',F(lb,order='fixed')),('SA',SA),('ILS',ILS)]:
  for m in ['f','nc','tot']:
    out[f'cat_{nm}_{m}']=unp(F(R,cat='restricted'),F(R,cat='unrestricted'),m)
    out[f'step_{nm}_{m}']=unp(F(R,step='5'),F(R,step='10'),m)
for m in ['f','nc','tot']:
  out[f'cool_{m}']=paired(F(SA,var='linear'),F(SA,var='geometric'),m)
# perturbation KW pooled within strategy (first fixed order + best)
vars_=sorted(set(r['var'] for r in ILS))
for st_ in ['first','best']:
  for m in ['f','nc','tot']:
    g=[[r[m] for r in ILS if r['var']==v and r['strat']==st_] for v in vars_]
    out[f'pert_{st_}_{m}']=(len(g[0]),[float(np.median(x)) for x in g],kruskal(*g).pvalue)
# ILS order
for m in ['f','nc','tot']:
  g=[[r[m] for r in ILS if r['strat']=='first' and r['order']==o] for o in ['fixed','circular','random']]
  out[f'ilsorder_{m}']=(len(g[0]),[float(np.median(x)) for x in g],kruskal(*g).pvalue)
  # paired circular vs fixed within ILS (same start per family)
# ILS restricted vs unrestricted per (step,strat)
for s in ['5','10']:
  out[f'ilscat_{s}']=unp(F(ILS,cat='restricted',step=s),F(ILS,cat='unrestricted',step=s),'f')
# ILS First vs Best
out['ils_fb_f']=unp(F(ILS,strat='best'),F(ILS,strat='first'),'f')
print(vars_)
for k,v in out.items(): print(k, [round(float(x),4) if not isinstance(x,list) else [round(y,1) for y in x] for x in v])
