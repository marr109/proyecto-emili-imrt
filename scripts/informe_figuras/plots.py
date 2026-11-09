import csv,glob,os,re,numpy as np,matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
exec(open('effects.py').read().split('out={}')[0])  # loads LS, SA, ILS, F
OUT='/Users/marrojasr/Documents/Proyecto de Título/informe/Imagenes/experimentos/'
C=['#2a78d6','#eb6834','#1baf7a','#eda100','#e87ba4','#008300','#6250d6','#e34948']
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,'grid.color':'#e6e6e3','grid.linewidth':0.6,'axes.axisbelow':True,'axes.edgecolor':'#52514e','text.color':'#0b0b0b','axes.labelcolor':'#0b0b0b'})
CONDS=[('5','restricted'),('5','unrestricted'),('10','restricted'),('10','unrestricted')]
CL={('5','restricted'):'Paso 5°\nrestringido',('5','unrestricted'):'Paso 5°\nsin restricción',('10','restricted'):'Paso 10°\nrestringido',('10','unrestricted'):'Paso 10°\nsin restricción'}
def boxes(ax,groups,series,metric,ylabel,conds=CONDS,labels=CL):
    n=len(series); w=0.8/n
    for j,(name,sel) in enumerate(series):
        data=[[r[metric] for r in sel(*c)] for c in conds]
        pos=np.arange(len(conds))+(j-(n-1)/2)*w
        b=ax.boxplot(data,positions=pos,widths=w*0.85,patch_artist=True,showfliers=False,medianprops=dict(color='#0b0b0b',lw=1.5),whiskerprops=dict(color='#52514e'),capprops=dict(color='#52514e'),boxprops=dict(edgecolor='#52514e',lw=0.8))
        for p in b['boxes']: p.set_facecolor(C[j]); p.set_alpha(0.55)
        rng=np.random.default_rng(j)
        for k,d in enumerate(data): ax.scatter(pos[k]+rng.uniform(-w*0.25,w*0.25,len(d)),d,s=9,color=C[j],zorder=3,linewidths=0)
        ax.plot([],[],'s',color=C[j],alpha=0.7,label=name)
    ax.set_xticks(range(len(conds))); ax.set_xticklabels([labels[c] for c in conds]); ax.set_ylabel(ylabel)
def two(series,fn,title=None):
    fig,ax=plt.subplots(1,2,figsize=(11,4.2))
    boxes(ax[0],None,series,'f','Objetivo final del FMO (menor es mejor)')
    boxes(ax[1],None,series,'nc','Resoluciones reales del FMO por corrida')
    ax[0].set_title('(a) Calidad',loc='left'); ax[1].set_title('(b) Costo',loc='left')
    ax[1].legend(frameon=False,loc='upper right',fontsize=9)
    fig.tight_layout(); fig.savefig(OUT+fn,dpi=200); plt.close(fig)
ON={'fixed':'Orden fijo','circular':'Orden circular','random':'Orden aleatorio'}
two([(ON[o],lambda s,c,o=o:F(LS,strat='first',order=o,step=s,cat=c)) for o in ['fixed','circular','random']],'ls_ordenes.png')
two([('Primera mejora (orden fijo)',lambda s,c:F(LS,strat='first',order='fixed',step=s,cat=c)),('Mejor mejora',lambda s,c:F(LS,strat='best',order='fixed',step=s,cat=c))],'ls_first_best.png')
# convergence
def traj(t):
    rows=list(csv.DictReader(open(t))); best=np.inf; xs=[];ys=[];k=0
    for r in rows:
        if r['cached']=='false': k+=1
        best=min(best,float(r['objective'])); xs.append(k); ys.append(best)
    return np.array(xs),np.array(ys)
fig,axs=plt.subplots(1,2,figsize=(11,4))
for a,s in zip(axs,['5','10']):
    for j,(lab,pat) in enumerate([('Primera mejora, orden fijo',f'nangshift{s}/*restricted/data/first'),('Primera mejora, orden circular',f'nangshift{s}/*restricted-shift-circular/data/first'),('Primera mejora, orden aleatorio',f'nangshift{s}/*restricted-shift-rand/data/first'),('Mejor mejora',f'nangshift{s}/*restricted/data/best')]):
        T=[traj(t) for t in glob.glob(B+f'/local_search/{pat}/seed*/trajectory.csv')]
        M=max(x[-1] for x,_ in T); grid=np.arange(1,M+1)
        Y=np.array([[y[np.searchsorted(x,g,side='right')-1] if g>=x[0] else y[0] for g in grid] for x,y in T])
        a.plot(grid,np.median(Y,0),color=C[j],lw=2,label=lab)
        a.fill_between(grid,*np.percentile(Y,[25,75],0),color=C[j],alpha=0.12,lw=0)
    a.set_title(f'Paso {s}° (30 corridas por curva, ambos catálogos)',loc='left'); a.set_xlabel('Resoluciones reales del FMO'); a.set_ylabel('Mejor objetivo encontrado')
axs[1].legend(frameon=False,fontsize=9)
fig.tight_layout(); fig.savefig(OUT+'ls_convergencia.png',dpi=200); plt.close(fig)
# ILS perturbations
PN=[('range3','P1: tres ángulos'),('adaptive','P2: tamaño adaptativo'),('mag2','P3: dos ángulos, ±5°/±10°'),('mag2-stepwise','P4: P3 paso a paso'),('adaptivemag-stepwise','P5: adaptativo paso a paso')]
fig,ax=plt.subplots(figsize=(11,4.4))
boxes(ax,None,[(lab,lambda s,c,v=v:F(ILS,var=v,strat='first',step=s,cat=c)) for v,lab in PN],'f','Objetivo final del FMO')
ax.legend(frameon=False,ncol=3,fontsize=9,loc='upper center',bbox_to_anchor=(0.5,1.16)); fig.tight_layout(); fig.savefig(OUT+'ils_perturbaciones.png',dpi=200); plt.close(fig)
# budget ecdf
fig,ax=plt.subplots(figsize=(7,4.2))
per=list(csv.DictReader(open(B+'/ils/comparisons/per_run_metrics.csv')))+list(csv.DictReader(open(B+'/sa/comparisons/per_run_metrics.csv')))
def frac(rs): return np.sort([int(r['noncached_idx_of_best'])/int(r['noncached_evals']) for r in rs])
seen=set();uniq=[]
for r in per:
    k=tuple(r.values())
    if k not in seen: seen.add(k); uniq.append(r)
for j,(lab,sel) in enumerate([('SA lineal',lambda r:r['family']=='SA' and r['variant']=='linear'),('SA geométrico',lambda r:r['family']=='SA' and r['variant']=='geometric'),('ILS, primera mejora',lambda r:r['family']=='ILS' and r['strategy']=='first'),('ILS, mejor mejora',lambda r:r['family']=='ILS' and r['strategy']=='best')]):
    f=frac([r for r in uniq if sel(r)]); ax.step(f,np.arange(1,len(f)+1)/len(f),where='post',color=C[j],lw=2,label=f'{lab} (n={len(f)})')
ax.axvline(0.5,color='#52514e',lw=0.8,ls='--')
ax.set_xlabel('Fracción del presupuesto (resoluciones reales) al encontrar el mejor plan'); ax.set_ylabel('Proporción acumulada de corridas'); ax.legend(frameon=False,fontsize=9,loc='upper left')
fig.tight_layout(); fig.savefig(OUT+'presupuesto_ecdf.png',dpi=200); plt.close(fig)
# ILS vs SA
fig,ax=plt.subplots(figsize=(11,4.4))
boxes(ax,None,[('Búsqueda local (primera mejora, orden fijo)',lambda s,c:F(LS,strat='first',order='fixed',step=s,cat=c)),('ILS P1, primera mejora, orden fijo',lambda s,c:F(ILS,var='range3',strat='first',order='fixed',step=s,cat=c)),('SA lineal',lambda s,c:F(SA,var='linear',step=s,cat=c)),('SA geométrico',lambda s,c:F(SA,var='geometric',step=s,cat=c))],'f','Objetivo final del FMO')
ax.legend(frameon=False,ncol=2,fontsize=9,loc='upper center',bbox_to_anchor=(0.5,1.17)); fig.tight_layout(); fig.savefig(OUT+'ils_sa.png',dpi=200); plt.close(fig)
# dosimetry scatter
fig,ax=plt.subplots(figsize=(7,4.2))
for j,(lab,fam) in enumerate([('Búsqueda local','local_search'),('SA','sa'),('ILS','ils')]):
    d=[r for r in csv.DictReader(open(f'{B}/{fam}/plans_dosimetry.csv')) if r['status']=='solved' and float(r['objective'])<85000]
    ax.scatter([float(r['objective']) for r in d],[float(r['ptv_d95']) for r in d],s=8,alpha=0.5,color=C[j],linewidths=0,label=f'{lab} (n={len(d)})')
ax.set_xlabel('Objetivo final del FMO'); ax.set_ylabel('D95 del PTV (Gy)'); ax.legend(frameon=False,fontsize=9)
fig.tight_layout(); fig.savefig(OUT+'dosimetria_objetivo_d95.png',dpi=200); plt.close(fig)
