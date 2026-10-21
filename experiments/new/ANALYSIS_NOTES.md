# Análisis de experimentos nuevos (modelo FMO por órgano, K=4, CERR_Prostate)

Generado por `scripts/analyze_new_experiments.py`. Todas las cifras de este documento se calcularon a partir de `trajectory.csv` y `run.log` de cada corrida (15 semillas por condición). No se leyó ningún `plans_dosimetry.csv`.

Convenciones:
- **Resultado final** = `FMO objective f*` de `run.log` (plan final del framework). **Mejor anytime** = mínimo de `trajectory.csv`. En búsqueda local y SA ambos coinciden (se verifica abajo); en ILS pueden diferir.
- `trajectory.csv` **no tiene columna de tiempo** (columnas: eval, angles_deg, objective, cached, términos por estructura). La unidad de todo análisis de 'cuándo' es la **evaluación** (fila), no minutos. Las filas cacheadas no cuestan un solve de FMO, por lo que el índice de evaluación no es proporcional al tiempo de reloj; por eso también se reporta el índice contando solo evaluaciones no cacheadas.
- Tests estadísticos: Mann-Whitney / Kruskal-Wallis / Wilcoxon bilaterales con n=15 por celda (poca potencia). Se realizan muchos tests y **no se aplicó corrección por comparaciones múltiples**; los p-valores son orientativos.

## 0. Verificaciones

| Familia | Corridas leídas | Final (f*) != mínimo trayectoria (tol. 0.01) |
|---|---|---|
| Búsqueda local (first+best, 3 órdenes) | 360 | 0 |
| SA | 120 | 0 |
| ILS (first 3 órdenes + best fijo) | 1200 | 105 |

Corridas sin línea `FMO objective f*` en run.log (se usó el mínimo de la trayectoria): 0.

## A. Búsqueda local

**A1.** Figuras: `local_search/comparisons/A1_first_orders_nangshift{5,10}.png` (solo First, órdenes fijo/circular/random, agrupado por catálogo).

**A2.** Figuras: `A2_first_vs_best_nangshift{5,10}.png`. En cada (paso, catálogo) se compara Best (orden fijo) contra el orden de First con **menor mediana de objetivo final** (elegido con los mismos datos; sesgo de selección leve a favor de First).

**A4.** Estadísticos (valores = objetivo final f*; 'eval. hasta valor final' = índice de la primera fila que alcanza el mínimo de la trayectoria):

| Paso | Catálogo | Mediana First-Fijo | Mediana First-Circular | Mediana First-Random | KW p (3 órdenes) | Orden elegido (menor mediana) | Mediana Best (fijo) | MWU p (First elegido vs Best) | Mediana anytime First | Mediana anytime Best | Mediana eval. hasta valor final (First / Best) | Mediana eval. no cacheadas hasta valor final (First / Best) | Mediana evaluaciones totales (First / Best) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 5 | Restricted | 75326 | 75326 | 75346 | 0.9892 | Fijo | 75395 | 0.4303 | 75326 | 75395 | 39.0 / 57.0 | 33.0 / 51.0 | 47.0 / 65.0 |
| 5 | Unrestricted | 74797 | 74726 | 75356 | 0.3111 | Circular | 75226 | 0.1982 | 74726 | 75226 | 41.0 / 93.0 | 28.0 / 76.0 | 49.0 / 105.0 |
| 10 | Restricted | 75560 | 75560 | 75346 | 0.3441 | Random | 75433 | 0.3947 | 75346 | 75433 | 23.0 / 39.0 | 22.0 / 32.0 | 31.0 / 47.0 |
| 10 | Unrestricted | 74780 | 74780 | 74950 | 0.9155 | Fijo | 75276 | 0.3090 | 74780 | 75276 | 27.0 / 51.0 | 22.0 / 43.0 | 35.0 / 65.0 |

Detalle por orden en `local_search/comparisons/A4_evals_to_final_per_condition.csv`.

**A3.** Figuras de convergencia (mejor objetivo acumulado vs. índice de evaluación; línea = mediana de las 15 semillas, banda = IQR; las corridas más cortas se rellenan con su valor final): `A3_conv_orders_nangshift{5,10}.png`, `A3_conv_first_vs_best_nangshift{5,10}.png` (x = número de evaluación, **incluye filas cacheadas**) y las mismas con sufijo `_noncached` (x = evaluaciones no cacheadas, es decir solves reales). First y Best requieren un número distinto de evaluaciones por iteración, por eso la comparación en x es de costo, no de iteraciones.

## B. Adecuación del presupuesto (18 min = 1080 s de reloj)

Para cada corrida: evaluaciones totales (filas), no cacheadas, e índice de evaluación en que se alcanzó por primera vez el mínimo de la trayectoria, como fracción del total. Como `trajectory.csv` no tiene tiempos, la fracción es **en evaluaciones**; la fracción en evaluaciones no cacheadas es una mejor aproximación del tiempo de reloj (los solves de FMO dominan el costo). Tabla agregada: SA por calendario (120 corridas = 2 pasos x 2 catálogos x 15 semillas x 2 calendarios, 60 por calendario); ILS por perturbación y estrategia, agrupando pasos, catálogos y órdenes (First: 180 corridas por perturbación, Best fijo: 60).

| Variante | n corridas | Mediana evals totales | Mediana evals no cacheadas | Fracción (evals) mediana / Q1 / Q3 | Fracción p90 | Fracción (no cacheadas) mediana / Q1 / Q3 | Fracción (no cach.) p90 | % corridas con mejor en <=50% (evals) | % corridas con mejor en <=50% (no cach.) |
|---|---|---|---|---|---|---|---|---|---|
| SA lineal | 60 | 76 | 73 | 0.76 / 0.52 / 0.93 | 0.98 | 0.76 / 0.51 / 0.94 | 0.98 | 23.3% | 25.0% |
| SA geométrico | 60 | 78 | 73 | 0.82 / 0.72 / 0.94 | 0.99 | 0.83 / 0.73 / 0.94 | 0.99 | 10.0% | 11.7% |
| ILS range3 (3 ang.), first | 180 | 84 | 72 | 0.75 / 0.55 / 0.93 | 0.98 | 0.75 / 0.55 / 0.93 | 0.99 | 20.0% | 19.4% |
| ILS range3 (3 ang.), best | 60 | 87 | 73 | 0.91 / 0.70 / 0.97 | 0.98 | 0.91 / 0.70 / 0.97 | 0.99 | 10.0% | 10.0% |
| ILS adaptive, first | 180 | 83 | 72 | 0.79 / 0.56 / 0.94 | 0.98 | 0.78 / 0.56 / 0.94 | 0.98 | 17.2% | 19.4% |
| ILS adaptive, best | 60 | 86 | 73 | 0.89 / 0.78 / 0.96 | 0.98 | 0.89 / 0.78 / 0.96 | 0.99 | 3.3% | 3.3% |
| ILS mag2 (2 ang.), first | 180 | 95 | 73 | 0.74 / 0.45 / 0.90 | 0.97 | 0.75 / 0.54 / 0.90 | 0.97 | 30.0% | 20.6% |
| ILS mag2 (2 ang.), best | 60 | 86 | 73 | 0.92 / 0.79 / 0.97 | 0.98 | 0.93 / 0.81 / 0.97 | 0.98 | 10.0% | 6.7% |
| ILS mag2-stepwise, first | 180 | 98 | 73 | 0.74 / 0.44 / 0.90 | 0.98 | 0.76 / 0.55 / 0.90 | 0.97 | 30.0% | 20.0% |
| ILS mag2-stepwise, best | 60 | 87 | 73 | 0.94 / 0.81 / 0.98 | 0.99 | 0.93 / 0.83 / 0.97 | 1.00 | 6.7% | 5.0% |
| ILS adaptivemag-stepwise, first | 180 | 113 | 73 | 0.64 / 0.37 / 0.92 | 0.98 | 0.71 / 0.50 / 0.91 | 0.97 | 38.9% | 26.1% |
| ILS adaptivemag-stepwise, best | 60 | 89 | 73 | 0.91 / 0.69 / 0.96 | 0.98 | 0.92 / 0.77 / 0.96 | 0.98 | 6.7% | 0.0% |

Desglose por paso y catálogo en `ils/comparisons/B_budget_adequacy.csv` y `sa/comparisons/B_budget_adequacy.csv`; métricas por corrida en `per_run_metrics.csv` de cada carpeta.

Figuras: `sa/comparisons/B_budget_fraction_{box,ecdf}.png`, `ils/comparisons/B_budget_fraction_{box,ecdf}.png`.

## C. Consistencia de evaluaciones con el mismo presupuesto de 18 min

Mediana [mín – máx] por corrida, agrupando catálogos (y órdenes en ILS First). 'Cacheadas' = filas con `cached=true` (reevaluaciones de un plan ya visto). 'No cacheadas' = solves FMO reales.

| Familia / variante | Paso | n | Total evals | Cacheadas | No cacheadas | Mediana % cacheadas |
|---|---|---|---|---|---|---|
| SA lineal | 5 | 30 | 77 [69 – 101] | 3 [0 – 14] | 73 [67 – 96] | 4.3% |
| SA geométrico | 5 | 30 | 78 [69 – 107] | 4 [1 – 18] | 73 [67 – 96] | 5.6% |
| ILS range3 (3 ang.), first | 5 | 90 | 85 [75 – 102] | 12 [4 – 26] | 72 [65 – 86] | 14.4% |
| ILS range3 (3 ang.), best | 5 | 30 | 87 [76 – 100] | 13 [10 – 21] | 72 [66 – 83] | 15.3% |
| ILS adaptive, first | 5 | 90 | 84 [74 – 101] | 12 [4 – 30] | 72 [66 – 86] | 14.2% |
| ILS adaptive, best | 5 | 30 | 87 [75 – 100] | 13 [8 – 21] | 73 [65 – 83] | 15.3% |
| ILS mag2 (2 ang.), first | 5 | 90 | 108 [76 – 294] | 34 [6 – 224] | 72 [66 – 84] | 31.9% |
| ILS mag2 (2 ang.), best | 5 | 30 | 87 [76 – 102] | 14 [10 – 27] | 72 [66 – 84] | 16.4% |
| ILS mag2-stepwise, first | 5 | 90 | 109 [77 – 309] | 35 [6 – 239] | 72 [66 – 85] | 32.0% |
| ILS mag2-stepwise, best | 5 | 30 | 88 [75 – 103] | 15 [10 – 29] | 73 [65 – 83] | 16.8% |
| ILS adaptivemag-stepwise, first | 5 | 90 | 126 [79 – 309] | 52 [6 – 234] | 73 [67 – 84] | 41.9% |
| ILS adaptivemag-stepwise, best | 5 | 30 | 89 [76 – 145] | 16 [10 – 72] | 73 [66 – 84] | 17.8% |
| SA lineal | 10 | 30 | 74 [67 – 99] | 2 [0 – 6] | 73 [67 – 93] | 2.0% |
| SA geométrico | 10 | 30 | 76 [68 – 104] | 3 [0 – 10] | 73 [67 – 95] | 4.2% |
| ILS range3 (3 ang.), first | 10 | 90 | 84 [72 – 98] | 12 [4 – 21] | 72 [66 – 79] | 14.6% |
| ILS range3 (3 ang.), best | 10 | 30 | 86 [76 – 101] | 13 [10 – 20] | 73 [66 – 84] | 15.1% |
| ILS adaptive, first | 10 | 90 | 83 [72 – 96] | 10 [4 – 26] | 73 [67 – 80] | 12.0% |
| ILS adaptive, best | 10 | 30 | 86 [75 – 100] | 12 [8 – 20] | 73 [65 – 83] | 13.7% |
| ILS mag2 (2 ang.), first | 10 | 90 | 88 [77 – 114] | 16 [4 – 41] | 73 [67 – 79] | 17.0% |
| ILS mag2 (2 ang.), best | 10 | 30 | 86 [76 – 125] | 12 [10 – 52] | 73 [66 – 83] | 14.0% |
| ILS mag2-stepwise, first | 10 | 90 | 90 [78 – 119] | 18 [6 – 47] | 73 [67 – 80] | 19.3% |
| ILS mag2-stepwise, best | 10 | 30 | 86 [75 – 126] | 13 [10 – 55] | 73 [65 – 82] | 15.1% |
| ILS adaptivemag-stepwise, first | 10 | 90 | 106 [80 – 183] | 34 [9 – 113] | 73 [67 – 79] | 31.1% |
| ILS adaptivemag-stepwise, best | 10 | 30 | 88 [75 – 116] | 14 [10 – 40] | 74 [65 – 85] | 16.3% |

Efecto del orden de escaneo en ILS First (mediana de evals totales / no cacheadas, agrupando perturbaciones, catálogos y ambos pasos):

| Orden | n | Mediana evals totales | Mediana no cacheadas |
|---|---|---|---|
| Fijo | 300 | 88 | 73 |
| Circular | 300 | 96 | 72 |
| Random | 300 | 85 | 72 |

Referencia: búsqueda local (termina por convergencia, no por presupuesto):

| Familia | Paso | n | Mediana total | Mediana cacheadas | Mediana no cacheadas |
|---|---|---|---|---|---|
| Búsqueda local first (fijo) | 5 | 30 | 55 | 8 | 45 |
| Búsqueda local best (fijo) | 5 | 30 | 77 | 12 | 64 |
| Búsqueda local first (fijo) | 10 | 30 | 34 | 6 | 28 |
| Búsqueda local best (fijo) | 10 | 30 | 52 | 8 | 44 |

Resumen global: mediana SA = 77 evals totales (73 no cacheadas); mediana ILS = 88 totales (73 no cacheadas); razón de medianas ILS/SA: 1.14 (totales), 1.00 (no cacheadas).

## D. ILS vs SA con el mismo presupuesto

**Elección del mejor ILS:** para cada (paso, catálogo) se toma, entre las 20 configuraciones ILS (5 perturbaciones x [First con 3 órdenes + Best fijo]), la de **menor mediana de objetivo final f***. La elección usa los mismos datos de la comparación, por lo que el 'mejor ILS' está sesgado a favor (selección entre 20 configuraciones); los p-valores contra SA son optimistas para ILS. Referencia adicional: búsqueda local First con el orden elegido en A.

| Paso | Catálogo | Mejor ILS (pert/orden/estrategia) | Mediana ILS mejor | Mediana SA lineal | Mediana SA geom. | Mediana LS First | Mediana ILS range3-First-fijo | MWU p ILS mejor vs SA lin | MWU p ILS mejor vs SA geom | MWU p SA lin vs geom | MWU p ILS mejor vs LS First | MWU p ILS range3-fijo vs SA lin | MWU p ILS range3-fijo vs SA geom |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 5 | Restricted | ILS adaptivemag-stepwise/Circular/first | 74380 | 75750 | 75282 | 75326 | 75395 | 0.0310 | 0.1466 | 0.2628 | 0.2900 | 0.3401 | 0.9339 |
| 5 | Unrestricted | ILS range3/Circular/first | 74332 | 74992 | 74651 | 74726 | 74880 | 0.1249 | 0.3837 | 0.4553 | 0.5069 | 0.7089 | 0.9669 |
| 10 | Restricted | ILS adaptive/Random/first | 74246 | 74278 | 74819 | 75346 | 74597 | 1.0000 | 0.9669 | 0.8195 | 0.4306 | 0.5614 | 0.7089 |
| 10 | Unrestricted | ILS range3/Circular/first | 73958 | 75198 | 75141 | 74780 | 74130 | 0.0344 | 0.0620 | 0.6187 | 0.1150 | 0.0620 | 0.1150 |

Figuras: `ils/comparisons/ils_vs_sa_nangshift{5,10}.png`. Se añade la columna con ILS range3 (valor de literatura), First orden fijo, sin selección, como comparación libre de sesgo de selección.

## E. Tamaño y tipo de perturbación en ILS

**E1.** Medianas del objetivo final por perturbación y Kruskal-Wallis entre las 5 perturbaciones:

| Paso | Catálogo | Estrategia | Orden | Mediana range3 (3 ang.) | Mediana adaptive | Mediana mag2 (2 ang.) | Mediana mag2-stepwise | Mediana adaptivemag-stepwise | KW p (5 pert.) |
|---|---|---|---|---|---|---|---|---|---|
| 5 | Restricted | first | Fijo | 75395 | 75460 | 75460 | 75460 | 75337 | 0.9991 |
| 5 | Restricted | first | Circular | 74422 | 74699 | 74649 | 74649 | 74380 | 0.9869 |
| 5 | Restricted | first | Random | 74760 | 74649 | 74649 | 74649 | 74649 | 0.9992 |
| 5 | Restricted | best | Fijo | 76208 | 76208 | 76208 | 76208 | 76208 | 1.0000 |
| 5 | Unrestricted | first | Fijo | 74880 | 74756 | 75248 | 75248 | 75188 | 0.9532 |
| 5 | Unrestricted | first | Circular | 74332 | 75047 | 74995 | 74995 | 74469 | 0.5984 |
| 5 | Unrestricted | first | Random | 75319 | 75373 | 75373 | 75373 | 75265 | 0.9491 |
| 5 | Unrestricted | best | Fijo | 76147 | 76147 | 76147 | 76147 | 75984 | 0.9999 |
| 10 | Restricted | first | Fijo | 74597 | 74623 | 74640 | 74640 | 74424 | 0.9998 |
| 10 | Restricted | first | Circular | 74559 | 74470 | 74550 | 74550 | 75097 | 0.9802 |
| 10 | Restricted | first | Random | 74605 | 74246 | 74844 | 74844 | 74882 | 0.9458 |
| 10 | Restricted | best | Fijo | 75117 | 75041 | 75117 | 75117 | 75117 | 0.9995 |
| 10 | Unrestricted | first | Fijo | 74130 | 74378 | 74346 | 74346 | 74329 | 0.9981 |
| 10 | Unrestricted | first | Circular | 73958 | 74130 | 74315 | 74315 | 74129 | 0.9368 |
| 10 | Unrestricted | first | Random | 74701 | 74615 | 74677 | 74677 | 74610 | 0.9750 |
| 10 | Unrestricted | best | Fijo | 74772 | 74772 | 74794 | 74794 | 74755 | 0.9999 |

**E2.** Comparaciones entre perturbaciones (Mann-Whitney no pareado). 'stepwise vs no stepwise' se mide con mag2 vs mag2-stepwise (misma magnitud, 2 ángulos) y adaptive vs adaptivemag-stepwise (esta última cambia también el tipo de magnitud, por lo que no es un contraste limpio de stepwise):

| Paso | Catálogo | Estrategia | Orden | MWU p range3 vs mag2 | MWU p mag2 vs mag2-stepwise | MWU p adaptive vs adaptivemag-stepwise | MWU p range3 vs adaptive |
|---|---|---|---|---|---|---|---|
| 5 | Restricted | first | Fijo | 1.0000 | 0.9669 | 0.9504 | 1.0000 |
| 5 | Restricted | first | Circular | 0.8357 | 1.0000 | 0.6783 | 0.6187 |
| 5 | Restricted | first | Random | 1.0000 | 1.0000 | 0.8845 | 0.9834 |
| 5 | Restricted | best | Fijo | 0.9834 | 0.9669 | 0.9669 | 0.9503 |
| 5 | Unrestricted | first | Fijo | 0.5066 | 1.0000 | 0.8519 | 0.7714 |
| 5 | Unrestricted | first | Circular | 0.2716 | 0.8032 | 0.4552 | 0.2211 |
| 5 | Unrestricted | first | Random | 0.7557 | 1.0000 | 0.5200 | 0.6630 |
| 5 | Unrestricted | best | Fijo | 1.0000 | 0.9008 | 1.0000 | 1.0000 |
| 10 | Restricted | first | Fijo | 0.9339 | 1.0000 | 0.9009 | 0.8357 |
| 10 | Restricted | first | Circular | 0.9669 | 1.0000 | 0.4807 | 0.7874 |
| 10 | Restricted | first | Random | 0.8195 | 1.0000 | 0.3951 | 0.8357 |
| 10 | Restricted | best | Fijo | 1.0000 | 1.0000 | 1.0000 | 0.9009 |
| 10 | Unrestricted | first | Fijo | 0.8035 | 1.0000 | 0.9669 | 0.8519 |
| 10 | Unrestricted | first | Circular | 0.5614 | 1.0000 | 0.9010 | 0.5614 |
| 10 | Unrestricted | first | Random | 0.6782 | 1.0000 | 0.9835 | 0.6782 |
| 10 | Unrestricted | best | Fijo | 1.0000 | 0.9834 | 0.9009 | 0.8681 |

Semillas en que mag2 y mag2-stepwise terminan con el mismo f* (|dif|<0.01): 227 de 240 (94.6%); las medianas y p-valores casi idénticos entre ambas columnas reflejan esto.

**E3.** Efecto del orden de escaneo dentro de cada perturbación (ILS First):

| Paso | Catálogo | Perturbación | Mediana Fijo | Mediana Circular | Mediana Random | KW p (3 órdenes) |
|---|---|---|---|---|---|---|
| 5 | Restricted | range3 (3 ang.) | 75395 | 74422 | 74760 | 0.2994 |
| 5 | Restricted | adaptive | 75460 | 74699 | 74649 | 0.5192 |
| 5 | Restricted | mag2 (2 ang.) | 75460 | 74649 | 74649 | 0.4526 |
| 5 | Restricted | mag2-stepwise | 75460 | 74649 | 74649 | 0.4258 |
| 5 | Restricted | adaptivemag-stepwise | 75337 | 74380 | 74649 | 0.4807 |
| 5 | Unrestricted | range3 (3 ang.) | 74880 | 74332 | 75319 | 0.4541 |
| 5 | Unrestricted | adaptive | 74756 | 75047 | 75373 | 0.8956 |
| 5 | Unrestricted | mag2 (2 ang.) | 75248 | 74995 | 75373 | 0.6279 |
| 5 | Unrestricted | mag2-stepwise | 75248 | 74995 | 75373 | 0.6934 |
| 5 | Unrestricted | adaptivemag-stepwise | 75188 | 74469 | 75265 | 0.3813 |
| 10 | Restricted | range3 (3 ang.) | 74597 | 74559 | 74605 | 0.8977 |
| 10 | Restricted | adaptive | 74623 | 74470 | 74246 | 0.8901 |
| 10 | Restricted | mag2 (2 ang.) | 74640 | 74550 | 74844 | 0.9192 |
| 10 | Restricted | mag2-stepwise | 74640 | 74550 | 74844 | 0.9192 |
| 10 | Restricted | adaptivemag-stepwise | 74424 | 75097 | 74882 | 0.9397 |
| 10 | Unrestricted | range3 (3 ang.) | 74130 | 73958 | 74701 | 0.3868 |
| 10 | Unrestricted | adaptive | 74378 | 74130 | 74615 | 0.6496 |
| 10 | Unrestricted | mag2 (2 ang.) | 74346 | 74315 | 74677 | 0.7499 |
| 10 | Unrestricted | mag2-stepwise | 74346 | 74315 | 74677 | 0.7499 |
| 10 | Unrestricted | adaptivemag-stepwise | 74329 | 74129 | 74610 | 0.9276 |

Figuras: `ils/comparisons/E_perturbation_{first,best}_nangshift{5,10}.png`, `E_order_by_perturbation_first_nangshift{5,10}.png`.

De 16 tests KW entre perturbaciones, 0 tienen p<0.05 (sin corrección por comparaciones múltiples).

**E4. Pareo ILS vs búsqueda local por semilla.**

Control: la solución inicial (fila 1) es la misma entre todas las configuraciones de una misma familia para cada (paso, catálogo, semilla): búsqueda local 60/60, ILS 60/60. Es decir, la semilla fija el inicio dentro de cada familia, pero no entre familias.

Comparación del objetivo de la fila 1 de `trajectory.csv` (solución inicial) entre ILS y búsqueda local para igual (paso, catálogo, semilla, estrategia, orden): coincidencias 0 de 1200 (0.0%).

El pareo NO es válido en todas las comparaciones; se usa Mann-Whitney no pareado donde no coincide el inicio (columna 'Test'). Filas de orden fijo; todas en `E4_ils_vs_ls_paired.csv`.

| Paso | Catálogo | Estrategia | Orden | Perturbación | Semillas con igual inicio | Mediana LS final | Mediana ILS final | Mediana (ILS - LS) | Semillas con ILS<LS | Test | p |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 5 | Restricted | first | Fijo | range3 (3 ang.) | 0/15 | 75326 | 75395 | 423 | 5/15 | MWU (no pareado) | 0.7400 |
| 5 | Restricted | first | Fijo | adaptive | 0/15 | 75326 | 75460 | 84 | 6/15 | MWU (no pareado) | 0.6333 |
| 5 | Restricted | first | Fijo | mag2 (2 ang.) | 0/15 | 75326 | 75460 | 423 | 5/15 | MWU (no pareado) | 0.5897 |
| 5 | Restricted | first | Fijo | mag2-stepwise | 0/15 | 75326 | 75460 | 423 | 5/15 | MWU (no pareado) | 0.5897 |
| 5 | Restricted | first | Fijo | adaptivemag-stepwise | 0/15 | 75326 | 75337 | 186 | 7/15 | MWU (no pareado) | 0.8519 |
| 5 | Restricted | best | Fijo | range3 (3 ang.) | 0/15 | 75395 | 76208 | 1282 | 4/15 | MWU (no pareado) | 0.1249 |
| 5 | Restricted | best | Fijo | adaptive | 0/15 | 75395 | 76208 | 1282 | 4/15 | MWU (no pareado) | 0.1057 |
| 5 | Restricted | best | Fijo | mag2 (2 ang.) | 0/15 | 75395 | 76208 | 1180 | 4/15 | MWU (no pareado) | 0.1057 |
| 5 | Restricted | best | Fijo | mag2-stepwise | 0/15 | 75395 | 76208 | 1180 | 4/15 | MWU (no pareado) | 0.1057 |
| 5 | Restricted | best | Fijo | adaptivemag-stepwise | 0/15 | 75395 | 76208 | 1282 | 4/15 | MWU (no pareado) | 0.1057 |
| 5 | Unrestricted | first | Fijo | range3 (3 ang.) | 0/15 | 74797 | 74880 | 432 | 6/15 | MWU (no pareado) | 0.8357 |
| 5 | Unrestricted | first | Fijo | adaptive | 0/15 | 74797 | 74756 | 862 | 7/15 | MWU (no pareado) | 0.6783 |
| 5 | Unrestricted | first | Fijo | mag2 (2 ang.) | 0/15 | 74797 | 75248 | 862 | 7/15 | MWU (no pareado) | 0.3195 |
| 5 | Unrestricted | first | Fijo | mag2-stepwise | 0/15 | 74797 | 75248 | 862 | 7/15 | MWU (no pareado) | 0.3195 |
| 5 | Unrestricted | first | Fijo | adaptivemag-stepwise | 0/15 | 74797 | 75188 | 862 | 7/15 | MWU (no pareado) | 0.3837 |
| 5 | Unrestricted | best | Fijo | range3 (3 ang.) | 0/15 | 75226 | 76147 | 854 | 6/15 | MWU (no pareado) | 0.1466 |
| 5 | Unrestricted | best | Fijo | adaptive | 0/15 | 75226 | 76147 | 854 | 6/15 | MWU (no pareado) | 0.1466 |
| 5 | Unrestricted | best | Fijo | mag2 (2 ang.) | 0/15 | 75226 | 76147 | 854 | 6/15 | MWU (no pareado) | 0.1466 |
| 5 | Unrestricted | best | Fijo | mag2-stepwise | 0/15 | 75226 | 76147 | 1565 | 5/15 | MWU (no pareado) | 0.1354 |
| 5 | Unrestricted | best | Fijo | adaptivemag-stepwise | 0/15 | 75226 | 75984 | 1565 | 6/15 | MWU (no pareado) | 0.1466 |
| 10 | Restricted | first | Fijo | range3 (3 ang.) | 0/15 | 75560 | 74597 | -601 | 10/15 | MWU (no pareado) | 0.1354 |
| 10 | Restricted | first | Fijo | adaptive | 0/15 | 75560 | 74623 | -601 | 11/15 | MWU (no pareado) | 0.0890 |
| 10 | Restricted | first | Fijo | mag2 (2 ang.) | 0/15 | 75560 | 74640 | -506 | 12/15 | MWU (no pareado) | 0.0890 |
| 10 | Restricted | first | Fijo | mag2-stepwise | 0/15 | 75560 | 74640 | -425 | 12/15 | MWU (no pareado) | 0.0890 |
| 10 | Restricted | first | Fijo | adaptivemag-stepwise | 0/15 | 75560 | 74424 | -522 | 12/15 | MWU (no pareado) | 0.0890 |
| 10 | Restricted | best | Fijo | range3 (3 ang.) | 0/15 | 75433 | 75117 | -37 | 8/15 | MWU (no pareado) | 0.8035 |
| 10 | Restricted | best | Fijo | adaptive | 0/15 | 75433 | 75041 | -37 | 8/15 | MWU (no pareado) | 0.5614 |
| 10 | Restricted | best | Fijo | mag2 (2 ang.) | 0/15 | 75433 | 75117 | -37 | 8/15 | MWU (no pareado) | 0.5897 |
| 10 | Restricted | best | Fijo | mag2-stepwise | 0/15 | 75433 | 75117 | -37 | 8/15 | MWU (no pareado) | 0.5897 |
| 10 | Restricted | best | Fijo | adaptivemag-stepwise | 0/15 | 75433 | 75117 | -37 | 8/15 | MWU (no pareado) | 0.5338 |
| 10 | Unrestricted | first | Fijo | range3 (3 ang.) | 0/15 | 74780 | 74130 | -449 | 9/15 | MWU (no pareado) | 0.2134 |
| 10 | Unrestricted | first | Fijo | adaptive | 0/15 | 74780 | 74378 | -416 | 10/15 | MWU (no pareado) | 0.1844 |
| 10 | Unrestricted | first | Fijo | mag2 (2 ang.) | 0/15 | 74780 | 74346 | -651 | 10/15 | MWU (no pareado) | 0.2628 |
| 10 | Unrestricted | first | Fijo | mag2-stepwise | 0/15 | 74780 | 74346 | -651 | 10/15 | MWU (no pareado) | 0.2628 |
| 10 | Unrestricted | first | Fijo | adaptivemag-stepwise | 0/15 | 74780 | 74329 | -345 | 10/15 | MWU (no pareado) | 0.1585 |
| 10 | Unrestricted | best | Fijo | range3 (3 ang.) | 0/15 | 75276 | 74772 | -39 | 8/15 | MWU (no pareado) | 0.9669 |
| 10 | Unrestricted | best | Fijo | adaptive | 0/15 | 75276 | 74772 | -106 | 8/15 | MWU (no pareado) | 0.9010 |
| 10 | Unrestricted | best | Fijo | mag2 (2 ang.) | 0/15 | 75276 | 74794 | 85 | 7/15 | MWU (no pareado) | 0.9669 |
| 10 | Unrestricted | best | Fijo | mag2-stepwise | 0/15 | 75276 | 74794 | 85 | 7/15 | MWU (no pareado) | 0.9339 |
| 10 | Unrestricted | best | Fijo | adaptivemag-stepwise | 0/15 | 75276 | 74755 | 64 | 7/15 | MWU (no pareado) | 0.8357 |

## Anexo: resultado final f* vs mejor anytime en ILS

Corridas ILS donde f* es peor que el mínimo de la trayectoria (>0.01): 105 de 1200. Mediana de la diferencia entre ellas: 225.4.

Por estrategia y perturbación:

| Estrategia | Perturbación | n | f* peor que anytime | Mediana de la diferencia |
|---|---|---|---|---|
| first | range3 (3 ang.) | 180 | 0 | n/a |
| first | adaptive | 180 | 0 | n/a |
| first | mag2 (2 ang.) | 180 | 0 | n/a |
| first | mag2-stepwise | 180 | 10 | 109.9 |
| first | adaptivemag-stepwise | 180 | 2 | 40.5 |
| best | range3 (3 ang.) | 60 | 19 | 225.4 |
| best | adaptive | 60 | 16 | 421.8 |
| best | mag2 (2 ang.) | 60 | 17 | 236.6 |
| best | mag2-stepwise | 60 | 26 | 210.8 |
| best | adaptivemag-stepwise | 60 | 15 | 471.9 |

## Perturbaciones por corrida de ILS (tope de intensidad 100)

Medido sobre ILS base (3 ángulos), `first`, con los tres órdenes juntos (45 corridas por celda), usando el detector de perturbaciones de `scripts/plot_ils_trajectory.py` (válido solo para perturbaciones por rango y `first`):

| Paso | Catálogo | Perturbaciones por corrida (mediana, mín-máx) | Corridas sin ninguna | Evaluación en que termina la 1ª búsqueda local (mediana) |
|---|---|---|---|---|
| 5 | restringido | 2 (0-3) | 1 de 45 | 40,5 |
| 5 | sin restricción | 2 (0-4) | 5 de 45 | 34,5 |
| 10 | restringido | 3 (1-4) | 0 de 45 | 24 |
| 10 | sin restricción | 3 (1-4) | 0 de 45 | 25 |

- En total, mediana de 2 perturbaciones por corrida y 6 de 180 corridas sin ninguna.
- Con ~73 solves reales por corrida, la primera búsqueda local consume entre 24 y 40 evaluaciones y quedan 2 o 3 perturbaciones; el diseño anterior hacía 10.
- No incluye `best` ni las perturbaciones de magnitud discreta (el detector no es válido para ellas).
- Consecuencia: la ausencia de diferencias entre perturbaciones (sección E) se lee como "con este presupuesto no se distinguen", no como "la perturbación no influye".
