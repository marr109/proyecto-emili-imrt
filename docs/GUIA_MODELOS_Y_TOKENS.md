# EMILI IMRT: guía de modelos FMO y uso de tokens

Esta guía está pensada para quien quiere (A) agregar un modelo distinto al
modelo FMO cuadrático por órgano y (B) usar la línea de comandos de `emili`.
Todo lo que se afirma aquí fue leído en el código del repositorio. Lo que no se
pudo comprobar está en la sección final "No verificado".

Rutas relativas a la raíz del repositorio `emili_imrt`.

---

# PARTE A. Cómo agregar un modelo distinto

## A.1 Dónde se conecta el modelo FMO

El problema BAO elige K ángulos. Para cada conjunto de ángulos, el FMO calcula
las intensidades de los beamlets. El FMO se resuelve con AMPL + Gurobi.

| Pieza | Archivo | Rol |
|---|---|---|
| Modelo baseline | `ampl_gurobi/fmo.mod` | Un solo `w_over` para todos los OAR. Congelado: no editar. |
| Modelo por órgano | `ampl_gurobi/fmo_organs.mod` | `w_over_organ[o]` por OAR y `normalize` opcional. |
| Copia para compartir | `ampl_gurobi/fmo_organs_envio.mod` | Mismo código, sin comentarios largos. |
| Lado C++ | `imrt/imrt_fmo.cpp` (`ImrtFmoSolver`) | Elige el `.mod`, carga parámetros, resuelve, lee resultados. |
| Contrato de datos | `imrt/imrt_fmo_source.h` (`IFmoDataSource`) | Qué datos entrega la instancia al solver. |
| Instancia CERR | `imrt/cerr_instance.cpp` (`CerrFmoSource`) | Lee `instances/CERR_Prostate` y fija Dmin, Dmax, pesos por defecto. |
| Caché y evaluación | `imrt/imrt_bao.cpp` (`BaoProblem`) | Llama a `fmo_.solve(...)` y guarda el resultado. |
| Driver Python | `ampl_gurobi/cerr_fmo_ampl.py`, `cerr_fmo_gurobi.py`, `cerr_fmo_osqp.py` | Mismo problema fuera de C++ (validación, dosimetría). |

### Selección en tiempo de ejecución

`ImrtFmoSolver::initAmpl()` (`imrt/imrt_fmo.cpp`, líneas 140 a 147) lee
`EMILI_FMO_MODEL`:

| Valor | Archivo que carga | Parámetros específicos |
|---|---|---|
| sin definir o `baseline` | `ampl_gurobi/fmo.mod` | `w_over` (viene de `source_.w_over()`, 0.5) |
| `organs` | `ampl_gurobi/fmo_organs.mod` | `EMILI_W_OVER_ORGANS`, `EMILI_FMO_NORMALIZE` |
| otro valor | nada | lanza `runtime_error`; el constructor lo captura y `emili` termina con "FMO solver initialization failed" |

El modelo se lee una sola vez (el objeto `ampl::AMPL` es persistente). Solo se
refrescan los datos por solve.

### Qué recibe el modelo

Parámetros fijos, cargados en `initAmpl()`:

| Parámetro AMPL | Origen |
|---|---|
| `n_ptv`, `n_oar` | suma de `n_boxets` de cada estructura |
| `max_intensity`, `w_under`, `w_ptv_over` | `IFmoDataSource` (`EMILI_MAX_INTENSITY`, `EMILI_W_PTV_OVER`) |
| `dmin`, `dmax_ptv` (por fila PTV), `dmax` (por fila OAR) | `FmoOrganRef` de cada órgano |
| `n_organs`, `w_over_organ`, `oar_organ`, `normalize` | solo con `organs` |

Datos por solve, cargados en `solve()` tras `reset data DIMLETS, PTV_DOSE, OAR_DOSE, d_ptv, d_oar`:

| Conjunto o parámetro | Contenido |
|---|---|
| `DIMLETS` | ids globales de los beamlets de los ángulos activos |
| `PTV_DOSE`, `OAR_DOSE` | pares (fila de voxel, beamlet) con dosis no nula |
| `d_ptv`, `d_oar` | tasa de dosis de cada par |

### Qué espera de vuelta el C++

Después de `ampl_->solve()` el código exige, por nombre:

| Requisito | Dónde | Si falta |
|---|---|---|
| `solve_result == "solved"` | `imrt_fmo.cpp` ~365 | excepción; el solve devuelve objetivo `1e30` |
| objetivo llamado `fmo_objective` | `getObjective("fmo_objective")` | excepción, `1e30` |
| variable `x` indexada por `DIMLETS` | intensidades `x_full` | excepción, `1e30` |
| variables `u` (filas PTV) y `v` (filas OAR) | desglose por órgano | excepción, `1e30` |

`FmoResult` (en `imrt_fmo.h`) contiene: `intensities`, `objective`,
`ptv_underdose_sq` (suma de u² por PTV) y `oar_overdose_sq` (suma de v² por OAR).
El desglose siempre se calcula como suma de cuadrados de `u` y `v`, sea cual sea la
fórmula del modelo.

`BaoProblem::calcObjectiveFunctionValue` guarda esos cuatro campos en la
solución. El informe clínico (`imrt/imrt_report.cpp`) recalcula las dosis desde
las intensidades, así que no depende de la fórmula del modelo.

## A.2 Receta: agregar una formulación FMO nueva

Ejemplo ilustrativo: una variante `organs_l1` (penalización lineal en lugar de
cuadrática). El nombre y la fórmula son un ejemplo, no existen en el repositorio.

### Paso 1. Crear el `.mod`

Copiar `fmo_organs.mod` a `ampl_gurobi/fmo_organs_l1.mod`. Conservar los nombres que
exige A.1 (`x`, `u`, `v`, `fmo_objective`, y los parámetros y conjuntos de datos).
Cambiar solo el objetivo.

```ampl
# ampl_gurobi/fmo_organs_l1.mod (fragmento; el resto igual a fmo_organs.mod)
param p_scale default 1;     # nuevo parámetro; el default reproduce el comportamiento previo

minimize fmo_objective:
    w_under * scale_ptv * sum {b in PTV_B} u[b]
  + sum {b in OAR_B} w_over_organ[oar_organ[b]] * scale_organ[oar_organ[b]] * v[b]
  + w_ptv_over * scale_ptv * sum {b in PTV_B} w[b];
```

Si el objetivo ya no es cuadrático, comprobar que el solver elegido lo admite
(ver "No verificado").

### Paso 2. Registrar el valor del selector (`imrt/imrt_fmo.cpp`)

Archivos y lugares que hay que tocar:

| # | Archivo | Contexto | Cambio |
|---|---|---|---|
| 1 | `imrt/imrt_fmo.cpp` | líneas 140 a 147 (`model_name`, `per_organ`, `ampl_->read`) | aceptar el nuevo valor y elegir el archivo |
| 2 | `imrt/imrt_fmo.cpp` | líneas 153 a 215 (rama `if (!per_organ)` / `else`) | si el modelo comparte parámetros con `organs`, reutilizar la rama; si no, agregar una rama propia |
| 3 | `imrt/imrt_fmo.h` | comentario de clase (líneas 32 a 36) | documentar el valor nuevo |
| 4 | `README.md` | tabla "Per-organ FMO model" | agregar la variable y el valor |
| 5 | `ampl_gurobi/README.md` | tabla de archivos | agregar el `.mod` |

Forma mínima del cambio 1 (la lógica actual usa un booleano `per_organ`; hay que
generalizarlo):

```cpp
const std::string model_name = envOr("EMILI_FMO_MODEL", "baseline");
if (model_name != "baseline" && model_name != "organs" && model_name != "organs_l1")
    throw std::runtime_error("EMILI_FMO_MODEL must be 'baseline', 'organs' or 'organs_l1', got '"
                             + model_name + "'");
const bool per_organ = (model_name == "organs" || model_name == "organs_l1");
const char* mod_file = (model_name == "organs_l1") ? "ampl_gurobi/fmo_organs_l1.mod"
                     : per_organ                   ? "ampl_gurobi/fmo_organs.mod"
                                                   : "ampl_gurobi/fmo.mod";
ampl_->read(repoPath(mod_file));
```

### Paso 3. Exponer parámetros como variables `EMILI_*`

Patrón existente: `envOr("EMILI_X", "default")` dentro de `initAmpl()`
(ejemplos: `EMILI_FMO_NORMALIZE`, `EMILI_W_OVER_ORGANS`), seguido de
`ampl_->getParameter("nombre").set(valor)`.

```cpp
const double p_scale = std::stod(envOr("EMILI_FMO_PSCALE", "1"));
ampl_->getParameter("p_scale").set(p_scale);
```

Reglas:

- El valor por defecto debe reproducir el comportamiento previo.
- Validar el texto y lanzar `runtime_error` con mensaje claro, como hace `EMILI_W_OVER_ORGANS`.
- Los parámetros de la instancia (Dmin, `EMILI_W_PTV_OVER`, ...) se leen en
  `imrt/cerr_instance.cpp`, no en `initAmpl()`. Para uno nuevo de ese tipo,
  agregar el accessor en `IFmoDataSource`.

### Paso 4. Mantener los scripts coherentes

| Archivo | Por qué |
|---|---|
| `scripts/run_catalog_job.sh` línea 24 (y los `run_*_organs.sh`) | exportan `EMILI_FMO_MODEL`, con `organs` por defecto |
| `irace/target-runner` línea 79 | fija el entorno del experimento con `EMILI_FMO_MODEL="organs"` |
| `ampl_gurobi/cerr_fmo_ampl.py` línea 33 | `MOD` fija `fmo_organs.mod` para todos los scripts Python |
| `ampl_gurobi/cerr_fmo_gurobi.py` | usa ese `MOD` |
| `ampl_gurobi/cerr_fmo_osqp.py` | reconstruye el QP a mano; no lee el `.mod`, hay que replicar la fórmula |
| `ampl_gurobi/dosimetry_batch.py` (`metrics`) | las columnas `term_*` asumen la forma cuadrática |

### Paso 5. Validar

1. Prueba mínima del `.mod` con datos sintéticos, siguiendo `ampl_gurobi/test_mod_normalize.py`
   (usa `fmo_organs_example.dat` y compara con un cálculo a mano):
   ```bash
   python ampl_gurobi/test_mod_normalize.py   # con el Python del entorno virtual de ampl_gurobi
   ```
2. Equivalencia con el comportamiento previo: con los defaults del paso 3 el objetivo
   debe coincidir con el de `organs`. El README de `ampl_gurobi` hace esta misma
   comparación entre `fmo.mod` y `fmo_organs.mod`.
3. Cruce C++ contra Python: `ampl_gurobi/cpp_fmo_crosscheck.cpp` imprime el objetivo
   para los ángulos {0, 90, 180, 270}. No está en `CMakeLists.txt`; hay que compilarlo
   aparte. Comparar con `cerr_fmo_ampl.py` usando las mismas variables `EMILI_*`.
4. Datos exactos del C++ para un `ampl` externo: `scripts/export_fmo_dat.cpp` escribe
   el `.dat` de un conjunto de ángulos (instrucciones de compilación en su cabecera).
5. Dosimetría de los planes finales: `ampl_gurobi/dosimetry_batch.py <carpeta> --out x.csv`
   (agregar `--geud` para columnas gEUD). Debe correrse con las mismas variables `EMILI_*`.
6. Corrida corta de humo (ver B.3). Un solve tarda unos 8 s en el equipo donde se
   probó esta guía, así que usar `-it 30` o más. Revisar en el log:
   `fmo_model: ...`, `FMO objective f*` y que no aparezca `[FMO] ERROR`.

## A.3 Agregar un PROBLEMA distinto (no solo un modelo)

Plantilla: las clases BAO de `imrt/imrt_bao.h` y `imrt/imrt_bao.cpp`.

| Clase base (`emilibase.h`) | Métodos a implementar | Ejemplo BAO |
|---|---|---|
| `Problem` | `calcObjectiveFunctionValue`, `evaluateSolution`, `problemSize` | `BaoProblem` |
| `Solution` | `clone`, `getSolutionRepresentation`, `isFeasible`, `getRawData`, `setRawData` | `BaoSolution` |
| `InitialSolution` | `generateSolution`, `generateEmptySolution` | `RandomKAnglesInit` |
| `Neighborhood` | `computeStep`, `reverseLastMove`, `begin`, `reset`, `random`, `size` | `AngleShiftNeighborhood` |
| `Perturbation` | `perturb` | `AngleShiftMultiPerturbation` |
| `Acceptance` | `accept` | `BaoImproveAccept` |
| `TabuMemory` | `tabu_check`, `forbid`, `reset` | `BaoTabuMemory` |
| `Shake` | `shake` | `MultiScaleAngleShake` |

Registro de tokens, en la clase `ImrtBuilder` (`imrt/imrt_builder.h/.cpp`):

1. Definir el nombre con `#define` (`imrt_builder.cpp`, líneas 12 a 57), por ejemplo `#define NEIGH_MIO "nmio"`.
2. Agregar una rama `else if (tm.checkToken(NEIGH_MIO))` en la función `buildX` que corresponda:
   `buildInitialSolution`, `buildNeighborhood`, `buildPerturbation`, `buildAcceptance`,
   `buildTermination`, `buildTabuTenure`, `buildShake`, `buildAlgo`.
3. Leer argumentos con `tm.getInteger()`, `tm.getDecimal()` y opciones con `tm.checkToken("...")`.
4. Un problema nuevo se identifica en `isCompatibleWith` / `canOpenInstance` / `openInstance`
   (hoy aceptan `imrt` y `baoimrt`). Lo más simple es una subclase nueva de `Builder`,
   registrada en `main.cpp` con `ps.addBuilder(&nuevo)`, junto a `emb` e `imrtb`.
5. `CMakeLists.txt` usa `aux_source_directory(./imrt SRC_LIST)`: un `.cpp` nuevo dentro de
   `imrt/` entra solo, pero hay que volver a ejecutar `cmake ..`.

Gramática e irace para un token nuevo:

- `grammars/imrt_bao.xml`: agregar la alternativa dentro de `bao:neighbourhood`,
  `bao:perturbation` o la derivación que corresponda; si tiene parámetros numéricos,
  declararlos arriba (`type`, `min`, `max`, `stepIfEnumerated`).
- `irace/parameters.irace` y `irace/target-runner`: agregar el parámetro y su traducción
  a texto (diccionarios `NB` y funciones `perturbation()`/`neighbourhoods()`).
- `scripts/check_irace_params.py`: actualizar `NUMERIC` (rangos), `TRACK` (parámetros
  que definen dependencias) y `ALGORITHMS`; documentar desvíos intencionales en `KNOWN`
  o `NOTES`. Ejecutarlo al final:
  ```bash
  python3 scripts/check_irace_params.py     # imprime OK / DIFF / KNOWN
  python3 scripts/sample_grammar.py 5 1     # imprime 5 comandos al azar
  ```

## A.4 Advertencias observadas en el código

| # | Advertencia | Evidencia |
|---|---|---|
| 1 | Los objetivos no son comparables entre modelos ni entre valores de `EMILI_DMIN`, `EMILI_W_PTV_OVER`, `EMILI_EXCLUDE_OVERLAP`. Cambian el problema. | comentarios en `cerr_instance.cpp` |
| 2 | Los pesos son por estructura. `EMILI_W_OVER_ORGANS` debe tener un valor por OAR, en el orden de `oarOrgans()` (BLADDER, RECTUM). Si la cuenta no coincide, error. | `imrt_fmo.cpp` ~180 |
| 3 | `baseline` ignora `EMILI_W_OVER_ORGANS` y `EMILI_FMO_NORMALIZE`: solo se leen en la rama `organs`. | `imrt_fmo.cpp` ~153 a 215 |
| 4 | La caché está indexada por el conjunto de ángulos ordenado, sin considerar slots. Nunca se vacía. | `imrt_bao.cpp` ~84 |
| 5 | Un solve fallido devuelve `1e30` y esa entrada también se guarda en la caché. | `imrt_bao.cpp` ~96 a 106 |
| 6 | Si el modelo no define `u` y `v`, todo solve devuelve `1e30` (la lectura está dentro del mismo `try`). | `imrt_fmo.cpp` ~381 a 410 |
| 7 | El desglose por órgano (`ptv_*_u2`, `oar_*_v2` en el csv) es suma de cuadrados aunque la fórmula del modelo sea otra. | `imrt_fmo.cpp`, `BaoProblem::openCsvLog` |
| 8 | El generador es `std::mt19937` con `std::uniform_int_distribution`. La distribución no está fijada por el estándar, así que la misma semilla puede dar trayectorias distintas entre plataformas o librerías. | `emilibase.cpp` líneas 43 a 49 |
| 9 | El presupuesto `-it` llega por `SIGALRM`, bloqueada durante cada solve. El tiempo real puede pasarse hasta un solve. | `AlarmGuard` en `imrt_fmo.cpp` |
| 10 | Los scripts Python fijan 2 PTV y 2 OAR (`PTV_ORGANS`, `OAR_ORGANS`) y la ruta de instancia por defecto está fijada en el código, por lo que conviene pasarla siempre de forma explícita. | `cerr_fmo_ampl.py` líneas 31 a 36 |
| 11 | `ampl_gurobi/README.md` dice que C++ no pasa `oar_organ` / `w_over_organ`. Está desactualizado: `initAmpl()` ya los pasa. | `imrt_fmo.cpp` ~190 a 197 |

---

# PARTE B. Guía de uso de los tokens

## B.1 Sintaxis general

```text
./build/emili <instancia> <problema> [K] [verbose] [csv <archivo>] <descripcion del algoritmo> [-it <seg>] [rnds <semilla>] [ps]
```

| Parte | Significado |
|---|---|
| `<instancia>` | carpeta CERR, p. ej. `instances/CERR_Prostate` (debe contener `beamletIndex.txt`) |
| `<problema>` | `baoimrt` (búsqueda de ángulos + FMO) o `imrt` (FMO clásico) |
| `K` | solo `baoimrt`: cantidad de ángulos activos (entero, 1 a 360) |
| `verbose` | opcional; imprime cada evaluación. Va antes de `csv` |
| `csv <archivo>` | opcional; una fila por evaluación (`eval,angles_deg,objective,cached` + columnas por órgano) |
| descripción | algoritmo y componentes, en el orden de B.2 |
| `-it <seg>` | presupuesto de reloj (segundos). `-ro <rho>` da `rho * tamaño del problema` |
| `rnds <semilla>` | semilla; sin este token es 0 |
| `ps` | imprime la solución al terminar |

Restricciones de orden, verificadas en `generalParser.cpp`:

- `-it` debe ir **antes** de `rnds`. El parser lee `-it` justo después de la descripción;
  si `rnds` aparece primero, `-it` no se consume y no hay límite de tiempo.
- `rnds` se busca en toda la línea (`seek`), pero conviene dejarlo al final.
- `verbose` va antes de `csv` (el builder los revisa en ese orden).
- Los componentes se leen en secuencia: cada algoritmo consume sus componentes en el orden
  de la tabla B.2.4.

## B.2 Tablas de tokens

### B.2.1 Solución inicial

| Token | Sintaxis | Significado |
|---|---|---|
| `irandomk` | `irandomk [unrestricted]` | K ángulos al azar. Por defecto, catálogo restringido a múltiplos de 5. |
| `ifirstk` | `ifirstk [unrestricted]` | los primeros K ángulos |
| `unrestricted` | sufijo | usa los 360 ángulos. También vale tras `nangswapsample`, `nangswapwindow`, `nangflip`, `pgreedysample`. |
| `izero`, `iuniform x`, `irandom max` | solo `imrt` | intensidades en 0, uniformes o aleatorias |

### B.2.2 Vecindarios

| Token | Sintaxis | Significado |
|---|---|---|
| `nangswap` | `nangswap` | cambia un ángulo activo por uno inactivo (vecindario completo K x (n-K)) |
| `nangshift` | `nangshift <step> [fixed\|circular\|randorder]` | mueve un ángulo activo +-step posiciones (grados) |
| `nangshiftmulti` | `nangshiftmulti <n> <s1> ... <sn>` | acumula los desplazamientos +-s_i de la lista |
| `nangswapsample` | `nangswapsample <n> [unrestricted]` | n swaps distintos sorteados en cada exploración (n >= 1) |
| `nangswapwindow` | `nangswapwindow <delta> [unrestricted]` | reemplazos a +-delta grados (1 a 180) |
| `nangflip` | `nangflip <delta> [unrestricted]` | reemplazos a +-delta del ángulo opuesto (0 a 180) |

Orden de recorrido de `nangshift`:

| Valor | Efecto |
|---|---|
| `fixed` (por defecto) | slots 0..K-1 |
| `circular` | parte del slot modificado en la ronda anterior y rota |
| `randorder` | permutación aleatoria de slots en cada exploración |

Dar `fixed` junto con otro orden termina el programa con error. Según el comentario
de `grammars/imrt_bao.xml`, el orden distinto de `fixed` solo tiene sentido con
`first`; con `best` el grammar no lo genera.

### B.2.3 Perturbaciones y shake

| Token | Sintaxis | Significado |
|---|---|---|
| `prangswap` | `prangswap <p>` | p swaps aleatorios |
| `prangshift` | `prangshift <step> <numSteps> [stepwise] [distinct\|replace]` | mueve numSteps ángulos una magnitud al azar en (step, 2 step) |
| `prangshiftadaptive` | `prangshiftadaptive <step> <maxNumSteps> <stagnation> [stepwise] [distinct\|replace]` | empieza con 1 ángulo; sube un nivel tras `stagnation` rondas sin mejora |
| `prangshiftmag` | `prangshiftmag <numSteps> <count> <m1> ... <mc> [stepwise] [distinct\|replace]` | como `prangshift`, magnitud sorteada de la lista, signo aleatorio |
| `prangshiftadaptivemag` | `prangshiftadaptivemag <maxNumSteps> <stagnation> <count> <m1> ... <mc> [stepwise] [distinct\|replace]` | versión adaptativa de la anterior |
| `pgreedy` | `pgreedy <D>` | destruye y reconstruye D ángulos de forma voraz (Iterated Greedy) |
| `pgreedysample` | `pgreedysample <D> <sample> [unrestricted]` | igual, evaluando solo `sample` candidatos por paso |
| `bangshake` | `bangshake <p_max>` | shake multiescala para VNS (nivel k reemplaza k ángulos) |
| `pershake` | `pershake <pert> <pert> ...` | shake formado por una lista de perturbaciones |

Opciones de las perturbaciones de desplazamiento:

- `stepwise`: evalúa tras cada movimiento individual.
- `distinct` (por defecto): mueve ángulos distintos. `replace`: sorteo con reemplazo.
- Si se dan ambas, `stepwise` va antes que `distinct|replace`.
- Con K = 4, un nivel de shake 4 reemplaza todos los ángulos (el grammar lo limita a 3).

### B.2.4 Aceptación, terminación y memorias tabú

| Token | Sintaxis | Significado |
|---|---|---|
| `improve` | `improve` | acepta solo si mejora (genérico de EMILI) |
| `aimprove` | `aimprove` | igual, versión del módulo IMRT |
| `baoimprove` | `baoimprove [rejectrepeated]` | mejora; con `rejectrepeated` rechaza conjuntos de ángulos ya visitados |
| `saacc` | `saacc <T_ini> <T_fin> <beta> <intervalo> <alpha>` | Metropolis con enfriamiento geométrico (en los scripts: `beta = 0`, `intervalo = 1`) |
| `sa_metropolis` | `sa_metropolis <T_ini> <T_fin> <beta>` | Metropolis, otra forma de enfriamiento |
| `accng` | `accng <aceptación>` | cambio de vecindario de VNS; toma la aceptación siguiente |
| `locmin` | `locmin` | termina en el óptimo local |
| `tmaxiter` | `tmaxiter <n>` | máximo de iteraciones |
| `ttime` | `ttime <seg>` | máximo de segundos |
| `maxstep` | `maxstep <n>` | máximo de pasos (genérico) |
| `TBao_fixed` | `TBao_fixed <tenure>` | prohíbe los últimos `tenure` conjuntos de ángulos visitados |
| `TBao_adaptive` | `TBao_adaptive <min> <max>` | igual, con tenure adaptativo (parte en (min+max)/2) |
| `TBao_attr` | `TBao_attr <tenure>` | prohíbe reinsertar cualquiera de los últimos `tenure` ángulos retirados |

### B.2.5 Algoritmos

| Token | Sintaxis (orden de componentes) |
|---|---|
| `first`, `best` | `first <inicial> <terminación> <vecindario>` |
| `ils` | `ils <búsqueda local completa> <terminación> <perturbación> <aceptación>` |
| `vnd` | `vnd first\|best <inicial> <terminación> <vecindario>+` (lista de vecindarios, en orden) |
| `vns` | `vns <búsqueda local o vnd> <terminación> <shake> accng <aceptación>` |
| `sa` | `sa <inicial> <terminación> <vecindario> <aceptación>` |
| `tabubao` | `tabubao [first\|best] [stall <n>] <inicial> <terminación> <vecindario> <memoria TBao_*>` |
| `nols` | sin búsqueda local (solo solución inicial) |

Notas de `tabubao`: sin `first|best` usa `best`. `stall <n>` detiene tras n iteraciones
sin un solve FMO nuevo (por defecto 100; 0 lo desactiva). Exige una memoria `TBao_*`.

## B.3 Ejemplos completos

Variables comunes (valores de los experimentos; ver `scripts/run_catalog_job.sh`):

| Variable | Valor del experimento | Significado |
|---|---|---|
| `EMILI_FMO_MODEL` | `organs` | usa `fmo_organs.mod` |
| `EMILI_W_OVER_ORGANS` | `0.5,0.8` | pesos de sobredosis: vejiga, recto |
| `EMILI_DMIN` | `68` | prescripción del PTV (Gy); el techo es 1.07 x Dmin |
| `EMILI_W_PTV_OVER` | `0.5` | penalización de punto caliente en el PTV (por defecto 0) |
| `EMILI_EXCLUDE_OVERLAP` | `1` | un voxel PTV+OAR cuenta solo como PTV |
| `EMILI_MAX_INTENSITY` | `100` | cota superior de cada beamlet (por defecto 15000) |
| `EMILI_FMO_NORMALIZE` | `0` | `1` divide cada término por el número de voxeles |
| `EMILI_GUROBI_OPTIONS` | vacío | opciones de Gurobi, p. ej. `threads=1` |
| `EMILI_AMPL_BIN_DIR`, `EMILI_GUROBI_BIN` | ruta del venv | binarios de AMPL y Gurobi |
| `EMILI_TABU_LOG`, `EMILI_TABU_TRACE`, `EMILI_PERT_TRACE` | opcional | depuración de tabú y perturbaciones |

```bash
export EMILI_FMO_MODEL=organs EMILI_W_OVER_ORGANS=0.5,0.8 EMILI_DMIN=68 \
       EMILI_W_PTV_OVER=0.5 EMILI_EXCLUDE_OVERLAP=1 EMILI_MAX_INTENSITY=100
# Si no existe ampl_gurobi/.venv, apuntar al entorno virtual donde estén instalados AMPL y Gurobi:
SP=<ruta-al-venv>/lib/python3.X/site-packages
export EMILI_AMPL_BIN_DIR=$SP/ampl_module_base/bin EMILI_GUROBI_BIN=$SP/ampl_module_gurobi/bin/gurobi
B="./build/emili instances/CERR_Prostate baoimrt 4 csv out.csv"
```

1. Búsqueda local, primera mejora, orden circular:

   ```bash
   $B first irandomk locmin nangshift 10 circular -it 1080 rnds 1
   ```

2. ILS con perturbación de desplazamiento (2 ángulos, con reemplazo, evaluación por paso):

   ```bash
   $B ils first irandomk locmin nangshift 10 tmaxiter 1000000 \
      prangshift 10 2 stepwise replace baoimprove rejectrepeated -it 1080 rnds 1
   ```

3. VND con cinco vecindarios:

   ```bash
   $B vnd first irandomk locmin nangshift 10 nangshift 5 nangswapsample 10 \
      nangswapwindow 20 nangflip 0 -it 1080 rnds 1
   ```

4. VNS general con shake multiescala:

   ```bash
   $B vns first irandomk locmin nangshift 10 tmaxiter 1000000 bangshake 3 \
      accng baoimprove rejectrepeated -it 1080 rnds 1
   ```

5. Recocido simulado (T de 1000 a 10, alpha 0.99):

   ```bash
   $B sa irandomk tmaxiter 1000000 nangshift 10 saacc 1000 10 0 1 0.99 -it 1080 rnds 1
   ```

6. Tabú con memoria por atributo:

   ```bash
   $B tabubao first irandomk tmaxiter 1000000 nangshift 10 TBao_attr 10 -it 1080 rnds 1
   ```

7. Extra: ILS con Iterated Greedy muestreado:

   ```bash
   $B ils first irandomk locmin nangshift 10 tmaxiter 1000000 pgreedysample 2 10 baoimprove -it 1080 rnds 1
   ```

Prueba de humo realizada para esta guía (los 7 comandos, con `-it 3`, csv fuera del repositorio):
los 7 terminaron con código 0 y el parser reconoció todos los componentes (el eco del
log lista cada uno con sus argumentos). Con `-it 3` no se alcanza a terminar el primer
solve (unos 8 s): el informe sale vacío (`K=0 beams`, objetivo `1e30`). Sirve para
validar la sintaxis, no para evaluar resultados. Para una corrida real usar `-it 1080`
(valor de los experimentos) o al menos `-it 30`.

Detalle observado: como el catálogo es 0..359 con paso de 1 grado, `step` en posiciones
equivale a grados.

## B.4 Gramática XML, irace y ejecutores

| Elemento | Ruta | Relación con los tokens |
|---|---|---|
| Gramática | `grammars/imrt_bao.xml` | genera solo la parte "descripción del algoritmo". El llamador antepone `instances/CERR_Prostate baoimrt 4 csv <traj>` y agrega `-it <seg> rnds <semilla>` al final. |
| Espacio de irace | `irace/parameters.irace` | subconjunto de la gramática escrito a mano: lo que se ejecutó |
| Traductor | `irace/target-runner` | arma la línea de comandos desde los parámetros, fija las variables `EMILI_*`, ejecuta `build/emili` y lee `FMO objective f*` del log. No lee el XML. |
| Escenario | `irace/scenario.irace` | presupuesto (300 experimentos) y paralelismo |
| Verificador | `scripts/check_irace_params.py` | compara límites, valores y dependencias de XML contra `parameters.irace` |
| Muestreo | `scripts/sample_grammar.py [n] [semilla]` | imprime n comandos de algoritmo al azar |

Se ejecutaron ambos scripts en esta revisión: `check_irace_params.py` informó `OK` en
las dependencias de cada algoritmo, y `sample_grammar.py` produjo líneas válidas
(por ejemplo `best irandomk locmin nangshift 5`).

Notas de la gramática:

- Es un superconjunto del espacio de irace: cubre todo lo implementado y usado en los experimentos.
- Usa `step` en {1, 5, 10, 20}, mientras que irace usa {5, 10}.
- Expone `distinct|replace`; la corrida de irace ejecutada usó siempre `distinct`.
- Ofrece las memorias tabú `TBao_adaptive`, `TBao_fixed` y `TBao_attr` (irace usó solo `TBao_adaptive`), las perturbaciones
  `pgreedy` y `pgreedysample` (sin evaluación experimental) y el catálogo `irandomk unrestricted`.
- Deja fuera a propósito: `nangswap` completo (272 vecinos no caben en el presupuesto), `nangshiftmulti`, las perturbaciones
  heredadas `prangswap` y `prandom`, el problema `imrt` solo de FMO, la opción `stall` del tabú y `unrestricted` dentro de `pgreedysample`.
- Validación: 150 líneas muestreadas de la gramática pasaron el analizador de EMILI sin errores (`-it 2`, solo valida sintaxis).
- `sa_metropolis` y `saacc` requieren T_ini > T_fin; no se expresa en el XML.

Ejecutores de experimentos:

| Script | Uso |
|---|---|
| `scripts/run_catalog_job.sh <familia> <config> <restricted\|unrestricted> <semilla>` | una corrida por (familia, config, catálogo, semilla). Familias: `vnd`, `vns`, `tabu`, `tabuattr`. Presupuesto en `BUDGET_SECONDS` (1080 por defecto); salida bajo `${OUT_ROOT:-experiments/new}`. |
| `scripts/run_*_organs.sh`, `run_*_batch.sh` | lotes por familia de algoritmo (ILS, SA, VND, VNS, tabú, búsqueda local) |
| `scripts/run_irace_validation.sh` | validación posterior a irace |

Ejemplo:

```bash
scripts/run_catalog_job.sh tabuattr t10 restricted 1
```

Este script define las mismas variables `EMILI_*` de B.3 con esos valores por defecto,
y en Linux resuelve los binarios de AMPL y Gurobi desde `ampl_gurobi/.venv`.

---

## No verificado

- Que un objetivo no cuadrático (`organs_l1` u otro) lo resuelva el driver `gurobi` de
  AMPL: es una consideración general sobre solvers, no se probó el `.mod` de ejemplo.
  El `.mod` y los cambios en C++ de A.2 son ilustrativos y no se compilaron ni ejecutaron.
- El comando exacto para compilar `cpp_fmo_crosscheck.cpp`. No hay objetivo en
  `CMakeLists.txt`; solo se comprobó que el archivo existe y qué imprime.
- El comportamiento de `imrt` (FMO clásico: `izero`, `iuniform`, `irandom`, `nshift`,
  `nswap`, `prandom`, `nactive`). Se leyó el builder; no se ejecutó.
- Que `-it` colocado después de `rnds` se ignore: sale de leer `getTime` y `getSeed`
  en `generalParser.cpp` (y del README), pero no se ejecutó esa combinación.
- Que `circular`/`randorder` sean rechazados con `best`: solo se cita el comentario del
  XML; el builder de `nangshift` no distingue `first` de `best`.
- Los tokens genéricos de EMILI que no son del módulo IMRT (`ig`, `tabu`, `gtabu`, `als`,
  `lspipe`, `timero`, `pmetro`, `nshake`, `mshake`...): solo se listaron los que usa el
  proyecto.
- La ejecución completa de `dosimetry_batch.py`, `validate_seed.py` y `test_mod_normalize.py`
  (no se corrieron; requieren licencia AMPL vigente y datos de experimentos).
- Resultados de corridas de más de 3 s: no se hicieron.
