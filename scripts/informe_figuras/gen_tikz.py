# rounds: list of (visit dict slot->order, accepted slot or None)
D={'Orden fijo':[({0:1,1:2,2:3,3:4,4:5},4)]*3+[({0:1,1:2,2:3,3:4,4:5},None)],
 'Orden aleatorio (un sorteo)':[({2:1,0:2,4:3},4),({3:1,1:2,0:3,4:4},4),({4:1},4),({1:1,3:2,0:3,2:4,4:5},None)],
 'Orden circular':[({0:1,1:2,2:3,3:4,4:5},4),({4:1},4),({4:1},4),({4:1,0:2,1:3,2:4,3:5},None)]}
c=0.52;L=[]
L.append(r'\begin{tikzpicture}[font=\footnotesize]')
x0=0
tot={}
for name,R in D.items():
    L.append(rf'\node[font=\small\bfseries,anchor=south] at ({x0+2.5*c},{0.15}) {{{name}}};')
    for s in range(5): L.append(rf'\node[anchor=south,font=\scriptsize] at ({x0+(s+.5)*c},{-0.25}) {{{s}}};')
    tot[name]=[]
    for r,(vis,acc) in enumerate(R):
        y=-(r+1)*c-0.25
        for s in range(5):
            if s in vis: fill='acc' if acc==s else ('fin' if acc is None else 'vis')
            else: fill='none'
            st={'acc':'fill=orange!55','fin':'fill=gray!35','vis':'fill=blue!15','none':'fill=white'}[fill]
            L.append(rf'\draw[{st},draw=black!50] ({x0+s*c},{y}) rectangle ++({c},{c});')
            if s in vis: L.append(rf'\node at ({x0+(s+.5)*c},{y+c/2}) {{{vis[s]}}};')
        ev=2*len(vis); tot[name].append(ev)
        L.append(rf'\node[anchor=west] at ({x0+5*c+0.05},{y+c/2}) {{{ev}}};')
    L.append(rf'\node[anchor=west,font=\scriptsize\itshape] at ({x0+5*c},{0.0}) {{evals.}};')
    x0+=5*c+1.6
# round labels
for r in range(4): L.append(rf'\node[anchor=east,font=\scriptsize] at (-0.05,{-(r+1)*c-0.25+c/2}) {{R{r+1}}};')
# legend
yl=-5*c-0.55
for i,(st,t) in enumerate([('fill=blue!15','casilla visitada sin mejora'),('fill=orange!55','mejora aceptada'),('fill=gray!35','última ronda: certifica el óptimo local')]):
    L.append(rf'\draw[{st},draw=black!50] ({i*4.3},{yl}) rectangle ++(0.3,0.3); \node[anchor=west] at ({i*4.3+0.35},{yl+0.15}) {{{t}}};')
# stacked bars
yb=yl-0.9; sc=0.25
cols=['blue!35','blue!20','blue!35','gray!35']
for k,(name,ev) in enumerate(tot.items()):
    y=yb-k*0.55; x=3.2
    L.append(rf'\node[anchor=east] at ({x-0.1},{y+0.2}) {{{name.replace(" (un sorteo)","")}}};')
    for j,e in enumerate(ev):
        col='gray!35' if j==3 else ('orange!35' if j%2==0 else 'orange!20')
        L.append(rf'\draw[fill={col},draw=white,line width=0.8pt] ({x},{y}) rectangle ++({e*sc},0.4); \node at ({x+e*sc/2},{y+0.2}) {{{e}}};')
        x+=e*sc
    L.append(rf'\node[anchor=west,font=\footnotesize\bfseries] at ({x+0.1},{y+0.2}) {{{sum(ev)}}};')
L.append(r'\end{tikzpicture}')
open('/Users/marrojasr/Documents/Proyecto de Título/informe/Estructura Documento/05.Propuesta de solucion/fig_ordenes.tex','w').write('% Generated figure: scan-order example\n'+'\n'.join(L)+'\n')
print(tot)
