# EMILI IMRT — Optimización de Ángulos de Haz y Mapas de Fluencia

<p align="center">
  <pre>
   ______ __  __ _____ _      _____ 
  |  ____|  \/  |_   _| |    |_   _|
  | |__  | \  / | | | | |      | |  
  |  __| | |\/| | | | | |      | |  
  | |____| |  | |_| |_| |____ _| |_ 
  |______|_|  |_|_____|______|_____|
  </pre>
</p>

**EMILI** (*Easily Modifiable Iterated Local search Implementation*) es un framework metaheurístico extensible en C++ para optimización combinatoria, extendido aquí con un **módulo de planificación de tratamientos IMRT** que resuelve el problema acoplado de **Optimización de Ángulos de Haz (BAO)** y **Optimización de Mapas de Fluencia (FMO)** mediante AMPL + Gurobi y búsqueda local iterada.

---

## Tabla de Contenidos

- [Arquitectura](#arquitectura)
- [Formulación del Problema IMRT](#formulación-del-problema-imrt)
- [Inicio Rápido](#inicio-rápido)
- [Uso de la Línea de Comandos](#uso-de-la-línea-de-comandos)
- [Ajuste de Hiperparámetros (irace)](#ajuste-de-hiperparámetros-irace)
- [Informe Clínico (ICRU-83)](#informe-clínico-icru-83)
- [Competencia de Algoritmos](#competencia-de-algoritmos)
- [Estructura del Proyecto](#estructura-del-proyecto)
- [Opciones de Compilación](#opciones-de-compilación)
- [Contribución](#contribución)
- [Licencia y Créditos](#licencia-y-créditos)

---

## Arquitectura

El código sigue un patrón de **framework + extensión**. La biblioteca núcleo de EMILI provee bloques metaheurísticos abstractos; el directorio `imrt/` los extiende para radioterapia.

```
┌──────────────────────────────────────────────────┐
│  NÚCLEO DEL FRAMEWORK  (emilibase.h/.cpp)        │
│  Problem · Solution · InitialSolution            │
│  Neighborhood · Perturbation · Acceptance        │
│  Termination · LocalSearch · TabuMemory          │
├──────────────────────────────────────────────────┤
│  PARSER / BUILDER  (generalParser.h/.cpp)        │
│  Gramática CLI → composición en ejecución        │
├──────────────────────────────────────────────────┤
│  MÓDULO IMRT  (imrt/)                            │
│  CerrFmoSource · ImrtProblem · ImrtFmoSolver     │
│  BaoProblem · AngleSwapNeighborhood              │
│  RandomAnglesPerturbation · GreedyAngles...      │
└──────────────────────────────────────────────────┘
```

### Metaheurísticas disponibles de fábrica

| Algoritmo | Token | BL interna | Componentes especiales |
|---|---|---|---|
| Búsqueda Local Iterada | `ils` | Mejor/Primera mejora | `prangswap`, `improve` |
| Búsqueda de Vecindario Variable | `vns` | Mejor/Primera mejora | `bangshake`, `accng` |
| Búsqueda Tabú (BAO) | `tabubao` | ninguna (ciclo propio, vecinos a partir de la solución actual) | `TBao_fixed`, `TBao_adaptive`, `TBao_attr` |
| Descenso de Vecindario Variable | `vnd` | lista de vecindarios | `nangshift`, `nangswapsample`, `nangswapwindow`, `nangflip` |
| Greedy Iterado | `ils` + `pgreedy` / `pgreedysample` | Mejor/Primera mejora | destrucción-reconstrucción (variante con muestreo) |
| Recocido Simulado | `ils nols` + `saacc` | Vacía (aceptación de Metropolis) | `saacc s e β i α` |

*Todos los algoritmos se componen en tiempo de ejecución mediante tokens de la CLI, sin necesidad de recompilar.*

---

## Formulación del Problema IMRT

### Contexto clínico

La Radioterapia de Intensidad Modulada entrega una dosis conformada a un Volumen Blanco de Planificación (PTV) mientras protege los Órganos en Riesgo (OAR). El problema consta de dos subproblemas acoplados:

1. **BAO (Optimización de Ángulos de Haz)** — seleccionar *K* ángulos del gantry entre *N* candidatos
2. **FMO (Optimización de Mapas de Fluencia)** — calcular las intensidades óptimas de los beamlets para los ángulos elegidos

### Objetivo de penalización cuadrática (FMO)

El subproblema FMO se resuelve de forma **exacta** mediante AMPL + Gurobi:

```
min  w_under · Σ || u_b ||²  +  w_over · Σ || v_b ||²  +  w_ptv_over · Σ || w_b ||²
s.t. D_ptv·x + u ≥ Dmin            (PTV underdose slack)
     D_oar·x  − v ≤ Dmax           (OAR overdose slack)
     D_ptv·x  − w ≤ 1.07·Dmin      (PTV hot-spot constraint)
     u, v, w ≥ 0
     0 ≤ x_j ≤ max_intensity
```

| Peso | Controla | Efecto clínico al aumentarlo |
|---|---|---|
| `w_under` | Penalización de subdosis en el PTV | ↑ fuerza la cobertura V95% |
| `w_over` | Penalización de sobredosis en OAR | ↑ protege los órganos sanos |
| `w_ptv_over` | Penalización de punto caliente en el PTV | ↓ reduce D2 / mejora la homogeneidad |

### Matriz de dosis dispersa

Los datos de deposición de dosis se leen desde el **formato de exportación CERR** sin procesar (`instances/CERR_Prostate`) mediante `CerrFmoSource`: una lista de vóxeles por órgano más archivos de dosis por ángulo, leyendo la cantidad de beamlets por ángulo en lugar de suponerla uniforme. El solver FMO ensambla un **QP compacto** usando únicamente los beamlets de los *K* ángulos activos, lo que mantiene el tamaño del problema en O(K) independientemente del tamaño del conjunto de ángulos candidatos.

---

## Inicio Rápido

### Prerrequisitos

- Compilador C++11 (gcc ≥ 4.7 o Clang)
- CMake ≥ 2.8
- AMPL + Gurobi, instalados localmente mediante `ampl_gurobi/.venv` (véase `ampl_gurobi/` para la configuración; la API C++ de AMPL se incluye integrada en `third_party/ampl_cppapi/`, por lo que no requiere una instalación separada para esa parte)

### Compilación

```bash
mkdir build && cd build
cmake ..
make -j$(nproc)
```

El solver FMO (`imrt/imrt_fmo.cpp`) resuelve el binario de AMPL y el controlador de Gurobi
en tiempo de ejecución mediante `EMILI_AMPL_BIN_DIR`/`EMILI_GUROBI_BIN` (por defecto, las
rutas dentro de `ampl_gurobi/.venv`). Si ese entorno virtual no está configurado, `baoimrt`
e `imrt` se compilarán y ejecutarán igualmente, pero cada resolución de FMO fallará de forma explícita
y devolverá un valor objetivo de peor caso en lugar de un resultado real.

El binario es `build/emili`.

---

## Uso de la Línea de Comandos

```bash
./build/emili <instance_dir> <problem> <algorithm_description> [options]
```

### BAO+FMO (K=4 ángulos, ILS con primera mejora)

```bash
./build/emili instances/PROSTATE_36ang baoimrt 4 \
    ils first irandomk tmaxiter 1 nangswap \
    tmaxiter 3 prangswap 2 improve \
    rnds 42
```

### Solo FMO (todos los ángulos, vecindario de desplazamiento de un único beamlet)

```bash
./build/emili instances/TINY_test imrt 36 \
    bi izero nshift 0.5 tmaxiter 100 \
    rnds 42
```

### Referencia rápida de tokens de algoritmo

| Categoría | Token | Significado |
|---|---|---|
| Problema | `imrt <dir> [nactive K] [verbose]` | FMO sobre todos los ángulos |
| Problema | `baoimrt <K> <dir> [verbose]` | BAO+FMO |
| Inicialización | `izero` / `iuniform X` / `irandom M` | Intensidades iniciales de FMO |
| Inicialización | `ifirstk` / `irandomk` | Conjunto inicial de ángulos de BAO |
| Vecindario | `nshift Δ` / `nswap` | Movimientos de beamlets de FMO |
| Vecindario | `nangswap` | Intercambio de ángulos de BAO |
| Perturbación | `prangswap P` / `pgreedy D` | Perturbaciones de BAO |
| Aceptación | `improve` / `saacc s e β i α` | Criterio de aceptación |
| Terminación | `tmaxiter N` / `tfeasible` | Condición de parada |

---

### Referencia de tokens de BAO

Para una guía completa de uso de los tokens y de cómo agregar un modelo FMO distinto, ver [docs/GUIA_MODELOS_Y_TOKENS.md](docs/GUIA_MODELOS_Y_TOKENS.md).

Vecindarios (todos modifican uno de los K ángulos activos):

| Token | Vecinos | Notas |
|---|---|---|
| `nangshift <step> [fixed\|circular\|randorder]` | 2K (ángulo -step, +step) | step en grados (1, 5, 10, 20); el orden de exploración solo importa para la primera mejora, por defecto `fixed` |
| `nangswap` | K x (catálogo - K) | intercambio completo, costoso (272 vecinos con el catálogo de 5 grados) |
| `nangswapsample <n> [unrestricted]` | n intercambios muestreados | un plan sin muestras que mejoren solo es probablemente un óptimo local |
| `nangswapwindow <deg> [unrestricted]` | intercambios dentro de +-deg del ángulo removido | |
| `nangflip 0` | como máximo K (ángulo +180) | |

Perturbaciones (ILS y la sacudida `pershake` del VNS general):

| Token | Significado |
|---|---|
| `prangshift <step> <n> [stepwise] [distinct\|replace]` | desplaza n ángulos en una magnitud aleatoria en (step, 2 step); `distinct` (por defecto) mueve n ángulos diferentes, `replace` sortea la posición de cada movimiento de forma independiente, de modo que un ángulo puede moverse más de una vez; `stepwise` evalúa después de cada movimiento individual |
| `prangshiftadaptive <step> <max> <stagnation> [distinct\|replace]` | comienza con 1 ángulo y agrega uno por cada `stagnation` rondas sin mejora |
| `prangshiftmag` / `prangshiftadaptivemag` | igual, con magnitudes discretas |
| `pgreedy <D>` / `pgreedysample <D> <sample>` | destruye D ángulos y reconstruye de forma voraz (la variante con muestreo evalúa `sample` candidatos por cada posición vacía) |

Búsqueda tabú:

```
tabubao [first|best] [stall <n>] irandomk tmaxiter <N> nangshift <step> <memory>
```

- `first|best` es opcional (por defecto `best`, la estrategia de los textos clásicos); `stall` detiene la búsqueda tras n iteraciones sin una nueva resolución de FMO (por defecto 100).
- Memorias: `TBao_fixed <tenure>` y `TBao_adaptive <min> <max>` prohíben conjuntos de ángulos completos ya visitados; `TBao_attr <tenure>` prohíbe reinsertar cualquiera de los últimos `tenure` ángulos que salieron del plan.
- Salida de depuración: `EMILI_TABU_LOG=<file>` escribe una fila csv por evento tabú (start, move, tabu_skip) con la lista de prohibidos antes y después; `EMILI_PERT_TRACE=1` imprime la posición y el cambio de ángulo de cada movimiento de perturbación.

Ejemplo:

```bash
./build/emili instances/CERR_Prostate baoimrt 4 csv trajectory.csv \
    tabubao first irandomk tmaxiter 1000000 nangshift 10 TBao_adaptive 3 8 \
    -it 1080 rnds 1
```

`-it <seconds>` debe ir antes de `rnds`. La gramática XML `grammars/imrt_bao.xml` genera estas líneas de comandos para la configuración automática; `scripts/check_irace_params.py` la verifica contra el espacio de parámetros de irace y `scripts/sample_grammar.py` imprime muestras aleatorias.

---

## Modelo de FMO por órgano y variables de experimento

El modelo de FMO por defecto es `ampl_gurobi/fmo.mod` (un único peso de sobredosis para todos los órganos en riesgo). El modelo por órgano
`ampl_gurobi/fmo_organs.mod` asigna a cada órgano su propio peso. Se selecciona en tiempo de ejecución, sin recompilar:

| Variable | Por defecto | Significado |
|---|---|---|
| `EMILI_FMO_MODEL` | `baseline` | `baseline` lee `fmo.mod`; `organs` lee `fmo_organs.mod` |
| `EMILI_W_OVER_ORGANS` | `w_over` para cada órgano | Pesos de sobredosis, separados por comas, en el orden de los órganos en riesgo (p. ej., `0.5,0.8`) |
| `EMILI_FMO_NORMALIZE` | `0` | `1` divide cada término por el número de vóxeles de su estructura |
| `EMILI_DMIN` | `65` | Prescripción del PTV en Gy (el techo del PTV es 1.07 x este valor) |
| `EMILI_W_PTV_OVER` | `0` | Peso de la penalización de punto caliente en el PTV |
| `EMILI_EXCLUDE_OVERLAP` | `0` | `1` cuenta los vóxeles compartidos por un PTV y un órgano solo como PTV |
| `EMILI_MAX_INTENSITY` | `15000` | Cota superior de la intensidad de cada beamlet |
| `EMILI_GUROBI_OPTIONS` | vacío | Se pasan al solver, p. ej., `threads=1` |
| `EMILI_AMPL_BIN_DIR`, `EMILI_GUROBI_BIN` | ruta del `.venv` del repositorio | Binarios de AMPL y del solver (necesarios en Linux) |

Los valores por defecto reproducen el comportamiento original. Los ejecutores en `scripts/` establecen las variables por órgano utilizadas en los experimentos
(por ejemplo, `scripts/run_catalog_job.sh`). La licencia del solver no forma parte del repositorio (`*.lic` se ignora);
véase `ampl_gurobi/README.md`.

## Ajuste de Hiperparámetros (irace)

`irace/` contiene la configuración de [irace](https://iridia.ulb.ac.be/irace/) utilizada para configurar las metaheurísticas de BAO:

| Archivo | Contenido |
|---|---|
| `irace/parameters.irace` | espacio de parámetros: algoritmo, estrategia, paso, orden de exploración, vecindarios, perturbaciones, tenencias tabú, enfriamiento |
| `irace/scenario.irace` | presupuesto e instancias |
| `irace/target-runner` | construye la línea de comandos de EMILI a partir de los parámetros muestreados y devuelve el objetivo |
| `irace/results/irace.Rdata` | resultado final de la ejecución realizada |

```bash
irace --scenario irace/scenario.irace
```

La gramática `grammars/imrt_bao.xml` es un superconjunto de este espacio: ofrece además los pasos 1 y 20, las posiciones `replace`,
las memorias tabú `TBao_fixed` y `TBao_attr`, las perturbaciones `pgreedy` y `pgreedysample` y el catálogo sin restricción
(`irandomk unrestricted`). `scripts/check_irace_params.py` lista esas diferencias y falla ante cualquier otra desviación.

---

## Informe Clínico (ICRU-83)

Cada ejecución de BAO+FMO imprime al final un resumen clínico del plan, con estas secciones:

- **Resultado del FMO:** ángulos seleccionados, valor del objetivo y dosis de prescripción.
- **Estadísticas de dosis por estructura (Gy):** número de vóxeles, Dmin, Dmean, Dmax, D95, D5 y D2.
- **Restricciones de DVH:** meta, valor alcanzado y estado de cada una (por ejemplo D95 del PTV frente a 0,95·Rx, D2 del PTV frente a 1,07·Rx y Dmax de los órganos en riesgo).
- **Índices de calidad del plan:** índice de conformidad (CI), índice de homogeneidad (HI) y cobertura V95% del PTV.

También se genera un **CSV de DVH** (`dvh.csv`) en el directorio de ejecución. La trayectoria de búsqueda de una ejecución (opción `csv <archivo>`) puede graficarse con `scripts/plot_trajectory.py`.

### Métricas de calidad del plan

| Métrica | Definición | Objetivo |
|---|---|---|
| **V95%** | % del volumen del PTV que recibe ≥ 95% de Rx | ≥ 95% |
| **CI** | Índice de Conformidad: V_Rx / V_PTV | → 1.0 |
| **HI** | Índice de Homogeneidad: (D2 − D98) / Rx | → 0 |
| **D95** | Dosis al 95% del PTV | ≥ 0.95·Rx |
| **D2** | Dosis casi máxima (2% del volumen) | ≤ 1.07·Rx |

---

## Competencia de Algoritmos

`run_comparison.sh` evalúa comparativamente múltiples algoritmos con semillas idénticas:

```bash
./run_comparison.sh
```

Salida por algoritmo en `comparison_results/`:

| Algoritmo | Archivo de registro | Componentes clave |
|---|---|---|
| ILS | `ILS.log` | ILS + inicialización aleatoria + perturbación de 2 intercambios |
| Tabú | `Tabu.log` | Tabú con memoria de conjuntos de ángulos de BAO |
| VNS | `VNS.log` | rVNS con sacudida multiescala |
| SA | `SA.log` | ILS + enfriamiento de Metropolis |
| Greedy | `Greedy.log` | ILS + destrucción-reconstrucción |

Los gráficos de trayectoria y de comparación se encuentran en `scripts/` (`plot_trajectory.py`, `plot_aggregate_comparison.py`, `plot_boxplot_comparison.py`).

Análisis por semilla de las ejecuciones de tabú: `scripts/tabu_replay_events.py` reconstruye los eventos de la lista tabú a partir de cada `trajectory.csv`
(validando la reproducción contra las evaluaciones registradas) y `scripts/tabu_seed_convergence_csv.py` escribe `convergence.csv`
(mejor objetivo hasta el momento frente a resoluciones reales de FMO).

---

## Estructura del Proyecto

```
emili_imrt/
├── main.cpp                   # Punto de entrada
├── emilibase.h/.cpp           # Framework núcleo (Problem, Solution, LS, ILS, VNS, Tabu...)
├── generalParser.h/.cpp       # Parser de tokens de la CLI + patrón Builder
├── setup.h                    # Banderas de tiempo de compilación
├── CMakeLists.txt             # Sistema de compilación
│
├── imrt/                      # Módulo de extensión IMRT
│   ├── cerr_instance.h/.cpp   # CerrFmoSource: cargador de exportaciones CERR sin procesar
│   ├── imrt_fmo_source.h/.cpp # Interfaz IFmoDataSource compartida por FMO/BAO
│   ├── imrt.h/.cpp            # FMO: ImrtProblem, vecindarios de beamlets
│   ├── imrt_fmo.h/.cpp        # Solver FMO con AMPL+Gurobi
│   ├── imrt_bao.h/.cpp        # BAO: vecindarios de ángulos, perturbaciones
│   └── imrt_builder.h/.cpp    # Builder de la CLI para tokens IMRT
│
├── grammars/imrt_bao.xml      # Gramática de los algoritmos de BAO (formato grammar2code)
├── irace/                     # Configuración de irace (parámetros, escenario, target-runner)
├── ampl_gurobi/               # Modelos de FMO, controladores del solver y herramientas de dosimetría
├── scripts/                   # Ejecutores de experimentos, herramientas de análisis y graficación
├── experiments/               # Salida de los ejecutores de experimentos (no se distribuye)
└── run_comparison.sh          # Ejecutor de benchmark multialgoritmo
```

---

## Opciones de Compilación

| Bandera de CMake | Por defecto | Descripción |
|---|---|---|
| `DEBUG_FLAGS` | ON | Símbolos de depuración (`-g`) |
| `O3_FLAGS` | OFF | Optimización `-O3` |
| `IRACE_OPTIMISED_FLAGS` | OFF | Banderas de GCC ajustadas con irace |
| `OPT_OPTIMISED_FLAGS` | OFF | Banderas de GCC ajustadas con OpenTuner |
| `WITH_STATS` | ON | Salida de estadísticas de la solución |
| `USE_NEW_MAIN` | ON | Carga de algoritmos basada en componentes |
| `PORTABLE` | OFF | Enlazado estático |

---

## Contribución

Este repositorio sigue un **flujo de trabajo basado en PR** gestionado por Gentle AI. Todos los cambios pasan por:

1. **Rama de funcionalidad** a partir de `develop`
2. **Cambios de código** con commits convencionales
3. **Pull Request** con revisión automatizada
4. **Fusión en `develop`** tras la aprobación

### Prefijos de commits convencionales

- `feat:` — nueva funcionalidad
- `fix:` — corrección de errores
- `perf:` — mejora de rendimiento
- `refactor:` — reestructuración de código
- `docs:` — documentación
- `ci:` — configuración de CI/irace

### Datos de instancias

Los archivos de instancias clínicas (matrices Dij, VOILists) **no se versionan** en este repositorio: `instances/` está en el `.gitignore`, de modo que el código se pueda reproducir con cualquier instancia propia. Contacte al mantenedor para obtener acceso a las instancias de entrenamiento.

#### Usar otra instancia

1. Copie la instancia a un directorio bajo `instances/` (por ejemplo `instances/MI_CASO`) con el formato de exportación CERR que lee `CerrFmoSource` (`imrt/cerr_instance.h`):
   - `<ORGANO>.txt`: ids globales de vóxel; el orden de las líneas es el índice local.
   - `<ORGANO>_<indice_angulo>.txt`: filas `id_global_voxel id_local_beamlet(desde 1) tasa_de_dosis`.
   - `beamletIndex.txt`: filas `indice_angulo inicio_global fin_global` (desde 1, inclusivo).
2. Ejecute indicando la instancia, por ejemplo `./build/emili instances/MI_CASO baoimrt 4 csv trajectory.csv <algoritmo> -it 1080 rnds 1`. Los lanzadores `scripts/run_catalog_job.sh` y `scripts/run_tabu_organs.sh` la toman de la variable `INSTANCE` (por defecto `instances/CERR_Prostate`); el resto de los scripts de `scripts/` todavía fija esa ruta y hay que editarla.
3. Ajuste lo que depende de las estructuras de la instancia: `EMILI_W_OVER_ORGANS` debe tener un peso por cada órgano en riesgo, en el mismo orden; y los scripts de Python de `ampl_gurobi/` asumen 2 PTV y 2 OAR (`PTV_ORGANS` y `OAR_ORGANS` en `cerr_fmo_ampl.py`), por lo que hay que cambiarlos si la instancia tiene otras estructuras.
4. Los valores del objetivo de instancias distintas no son comparables entre sí; compare algoritmos dentro de una misma instancia.

---

## Licencia y Créditos

- **Framework EMILI**: Federico Pagnozzi (BSD 2-Clause) — `federico.pagnozzi@ulb.ac.be`
- **Extensión IMRT** (BAO+FMO, integración AMPL+Gurobi, ajuste con irace): Marco Rojas
- **AMPL / Gurobi**: AMPL Optimization Inc. / Gurobi Optimization, LLC

```
EMILI IMRT — BAO & FMO optimization for radiation therapy planning
Copyright (c) 2014 Federico Pagnozzi. All rights reserved.
This file is distributed under the BSD 2-Clause License. See LICENSE.txt
```
