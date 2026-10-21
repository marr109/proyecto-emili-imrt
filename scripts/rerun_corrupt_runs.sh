#!/usr/bin/env bash
set -uo pipefail

# Re-ejecuta las corridas cuyo FMO devolvio valores invalidos.
#
# Que paso
# --------
# ImrtFmoSolver::solve() llamaba a ampl_->solve() sin verificar el resultado y
# leia el objetivo inmediatamente despues. Cuando el solve no se completaba
# --Gurobi sin licencia, modelo infactible, solver interrumpido-- AMPL conservaba
# el objetivo de la instancia ANTERIOR y lo devolvia sin senalar nada. Con el
# solver caido por completo el valor emitido era 0, que en un problema de
# minimizacion es el mejor posible: la busqueda lo habria tomado como optimo.
#
# Sintoma: decenas de configuraciones de angulos distintas con el mismo objetivo
# hasta el decimal, y busquedas que cortaban por falso minimo local.
# Corregido en imrt_fmo.cpp verificando solve_result.
#
# Alcance: 10 corridas, TODAS anteriores a la matriz de presupuesto de tiempo.
# experiments/ils/*/*budget1080s* (0 de 300) y experiments/sa (0 de 30) estan
# limpios. Este script cubre las 6 que sustentan tablas o figuras del informe;
# quedan fuera tabu/ y _legacy/, que nada cita.
#
# Por que invoca los comandos directamente
# ----------------------------------------
# Una version previa delegaba en los scripts de lote confiando en que saltearian
# las corridas ya completas. Esa suposicion es falsa en este repositorio: los
# scripts detectan lo completo buscando "Found solution" en el run.log, y NINGUNO
# de los run.log del servidor tiene esa marca --se generaron con una version
# anterior del binario que no la imprimia--. El resultado fue que el lote empezo
# a regenerar desde seed01 y sobrescribio una corrida sana antes de detenerlo.
# Por eso aqui cada corrida corrupta se invoca por separado, con su comando
# exacto reconstruido de los scripts de lote correspondientes, y no se toca
# ninguna otra semilla.

BACKUP="experiments/_corrupt_backup_$(date +%Y%m%d)"
DRY="${DRY:-0}"

say() { echo "[$(date '+%F %H:%M:%S')] $*"; }

# grep -a en vez de strings: strings vive en /usr/bin y no siempre esta en el
# PATH que hereda un screen desatendido, y entonces la tuberia sale vacia y el
# guard rechaza un binario que si tiene el fix.
if ! grep -aq "solve_result=" ./build/emili 2>/dev/null; then
  echo "ERROR: ./build/emili no contiene la verificacion de solve_result." >&2
  echo "       Compilar primero, o instalar build-fix/emili en build/emili." >&2
  exit 1
fi

if [ -z "${EMILI_AMPL_BIN_DIR:-}" ] || [ -z "${EMILI_GUROBI_BIN:-}" ]; then
  SP="$(echo ampl_gurobi/.venv/lib/python*/site-packages)"
  [ -d "$SP/ampl_module_base" ] || { echo "No se encontro ampl_module_base bajo $SP" >&2; exit 1; }
  export EMILI_AMPL_BIN_DIR="$PWD/$SP/ampl_module_base/bin"
  export EMILI_GUROBI_BIN="$PWD/$SP/ampl_module_gurobi/bin/gurobi"
fi
export EMILI_GUROBI_OPTIONS="${EMILI_GUROBI_OPTIONS:-threads=4}"

# dir | semilla | tokens del algoritmo, tal como los arma su script de lote
#   ils/*            -> run_ils_batch.sh          (prangshift STEP 3)
#   ils/*-adaptive/* -> run_ils_adaptive_batch.sh (prangshiftadaptive STEP 4 2)
#   local_search/*   -> corrida manual            (locmin, sin perturbacion)
RUNS=(
  "experiments/ils/nangshift5/restricted/data/first/seed12|12|ils first irandomk locmin nangshift 5 tmaxiter 10 prangshift 5 3 baoimprove rejectrepeated"
  "experiments/ils/nangshift5/restricted-adaptive/data/first/seed03|3|ils first irandomk locmin nangshift 5 tmaxiter 10 prangshiftadaptive 5 4 2 baoimprove rejectrepeated"
  "experiments/ils/nangshift5/unrestricted-adaptive/data/first/seed03|3|ils first irandomk unrestricted locmin nangshift 5 tmaxiter 10 prangshiftadaptive 5 4 2 baoimprove rejectrepeated"
  "experiments/ils/nangshift10/restricted-adaptive/data/first/seed03|3|ils first irandomk locmin nangshift 10 tmaxiter 10 prangshiftadaptive 10 4 2 baoimprove rejectrepeated"
  "experiments/ils/nangshift10/unrestricted-adaptive/data/first/seed03|3|ils first irandomk unrestricted locmin nangshift 10 tmaxiter 10 prangshiftadaptive 10 4 2 baoimprove rejectrepeated"
  "experiments/local_search/nangshift10/unrestricted/data/best/seed02|2|best irandomk unrestricted locmin nangshift 10"
)

say "regenerando ${#RUNS[@]} corridas invalidas (ninguna otra semilla se toca)"
for entry in "${RUNS[@]}"; do
  IFS='|' read -r dir seed tokens <<< "$entry"

  # Respalda si todavia existe: una ejecucion previa pudo haberla removido ya.
  if [ -d "$dir" ]; then
    dest="${BACKUP}/${dir#experiments/}"
    if [ "$DRY" = "1" ]; then
      say "  [dry] respaldaria $dir"
    else
      mkdir -p "$(dirname "$dest")" && cp -r "$dir" "$dest" && rm -rf "$dir"
      say "  respaldada: ${dir#experiments/}"
    fi
  fi

  if [ "$DRY" = "1" ]; then
    say "  [dry] ./build/emili ... csv ${dir}/trajectory.csv ${tokens} rnds ${seed}"
    continue
  fi

  mkdir -p "$dir"
  say "  corriendo ${dir#experiments/}"
  # shellcheck disable=SC2086
  ./build/emili instances/CERR_Prostate baoimrt 4 csv "${dir}/trajectory.csv" \
    $tokens rnds "$seed" > "${dir}/run.log" 2>&1

  if grep -q "Conformity Index\|Found solution" "${dir}/run.log" 2>/dev/null; then
    say "    OK"
  else
    say "    SIN MARCA DE CORRIDA COMPLETA -- revisar ${dir}/run.log"
  fi
done

say "COMPLETO. Verificar con la deteccion de objetivos repetidos antes de reanalizar."
