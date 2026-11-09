import csv,sys,glob,numpy as np,matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0,'/Users/marrojasr/Documents/code/emili_imrt/scripts')
from plot_ils_trajectory import detect_perturbations,group_bursts,collapse_bursts
B='/Users/marrojasr/Documents/code/emili_imrt/experiments/new'
OUT='/Users/marrojasr/Documents/Proyecto de Título/informe/Imagenes/experimentos/'
C=['#2a78d6','#eb6834','#1baf7a','#eda100','#e87ba4','#008300','#6250d6','#e34948']
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,'grid.color':'#e6e6e3','grid.linewidth':0.6,'axes.axisbelow':True})
def rd(t):
    R=list(csv.DictReader(open(t)))
    return [int(r['eval']) for r in R],[float(r['objective']) for r in R],[[int(a) for a in r['angles_deg'].split(';')] for r in R],[r['cached']=='true' for r in R]
def cd(a,b): d=abs(a-b)%360; return min(d,360-d)
def nb(a,b,step): 
    d=[j for j in range(len(a)) if a[j]!=b[j]]; return len(d)==1 and cd(a[d[0]],b[d[0]])==step
mode=sys.argv[1]
if mode=='scan':
    # ILS P1 first fixed step10 restricted: count bursts per seed
    for s in ['5','10']:
      for seed in range(1,16):
        t=f'{B}/ils/nangshift{s}/restricted/data/first/seed{seed:02d}/trajectory.csv'
        e,o,a,c=rd(t); rec=collapse_bursts(group_bursts(detect_perturbations(e,o,a,int(s))))
        print(s,seed,len(e),len(rec),round(min(o)))
    sys.exit()
if mode=='plot':
    s,seed=sys.argv[2],int(sys.argv[3])
    # ILS
    t=f'{B}/ils/nangshift{s}/restricted/data/first/seed{seed:02d}/trajectory.csv'
    e,o,a,c=rd(t); rec=collapse_bursts(group_bursts(detect_perturbations(e,o,a,int(s))))
    fig,ax=plt.subplots(figsize=(10,4))
    ax.scatter(e,o,s=10,color='#9a9a96',label='Configuración evaluada',zorder=2)
    ax.step(e,np.minimum.accumulate(o),where='post',color=C[0],lw=2,label='Mejor objetivo encontrado',zorder=3)
    for k,b in enumerate(rec):
        ax.axvspan(b['eval_primer_movimiento']-0.5,b['eval_despues']+0.5,color=C[1],alpha=0.25,lw=0,label='Perturbación (nuevo punto de partida)' if k==0 else None)
        ax.scatter([b['eval_despues']],[b['obj_despues']],color=C[1],s=30,zorder=4,marker='v')
    ax.set_ylim(min(o)-300,np.percentile(o,97)+300)
    cuts=[1]+[b['eval_despues'] for b in rec]+[e[-1]+1]
    for k in range(len(cuts)-1):
        x0,x1=cuts[k],cuts[k+1]
        if x1-x0>3: ax.text((x0+x1)/2,min(o)-150,f'Búsqueda local {k+1}',ha='center',va='bottom',fontsize=8.5,color='#52514e')
    ax.set_xlabel('Número de evaluación'); ax.set_ylabel('Objetivo del FMO'); ax.legend(frameon=False,fontsize=9,loc='upper right')
    fig.tight_layout(); fig.savefig(OUT+'traj_ils.png',dpi=200); plt.close(fig)
    print('ILS bursts',[(b['eval_primer_movimiento'],b['eval_despues']) for b in rec],len(e),min(o))
def nb1(a,b,step=None): return sum(x!=y for x,y in zip(a,b))==1
if mode=='sa':
    s,cat,var,seed=sys.argv[2],sys.argv[3],sys.argv[4],int(sys.argv[5])
    suf='-geometric' if var=='geometric' else ''
    t=f'{B}/sa/nangshift{s}/{cat}-budget{suf}/data/seed{seed:02d}/trajectory.csv'
    e,o,a,c=rd(t); cur=a[0];co=o[0]; curs=[co]; acc=[];bad=0;worse=[]
    for i in range(1,len(a)):
        prev=a[i-1]
        if nb1(a[i],cur) or a[i]==cur: pass
        elif nb1(a[i],prev) or a[i]==prev:
            if o[i-1]>co: worse.append((e[i-1],o[i-1]))
            cur,co=prev,o[i-1]
        else: bad+=1
        curs.append(co)
    print('unexplained',bad,'worse accepted',len(worse),len(e))
    fig,ax=plt.subplots(figsize=(10,4))
    ax.scatter(e,o,s=10,color='#9a9a96',label='Vecino evaluado',zorder=2)
    ax.step(e,curs,where='post',color=C[1],lw=1.6,label='Solución actual (reconstruida)')
    ax.step(e,np.minimum.accumulate(o),where='post',color=C[0],lw=2,label='Mejor objetivo encontrado')
    if worse: ax.scatter(*zip(*worse),marker='^',s=36,color=C[7],zorder=4,label='Empeoramiento aceptado')
    ax.set_ylim(min(o)-300,np.percentile(o,98)+300)
    ax.set_xlabel('Número de evaluación'); ax.set_ylabel('Objetivo del FMO'); ax.legend(frameon=False,fontsize=9,loc='upper right')
    fig.tight_layout(); fig.savefig(OUT+'traj_sa.png',dpi=200); plt.close(fig)
    print([w[0] for w in worse])
if mode=='ls':
    s,cat,seed=sys.argv[2],sys.argv[3],int(sys.argv[4])
    fig,ax=plt.subplots(figsize=(10,4))
    for j,(lab,v,st) in enumerate([('Primera mejora, orden fijo',cat,'first'),('Primera mejora, orden circular',cat+'-shift-circular','first'),('Mejor mejora',cat,'best')]):
        e,o,a,c=rd(f'{B}/local_search/nangshift{s}/{v}/data/{st}/seed{seed:02d}/trajectory.csv')
        k=np.cumsum([not x for x in c])
        ax.step(k,np.minimum.accumulate(o),where='post',color=C[j],lw=2,label=f'{lab} ({k[-1]} resoluciones, final {min(o):.0f})')
        print(lab,k[-1],len(e),min(o))
    ax.set_xlabel('Resoluciones reales del FMO'); ax.set_ylabel('Mejor objetivo encontrado'); ax.legend(frameon=False,fontsize=9)
    fig.tight_layout(); fig.savefig(OUT+'traj_ls.png',dpi=200); plt.close(fig)
