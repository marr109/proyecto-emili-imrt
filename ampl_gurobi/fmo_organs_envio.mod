# =============================================================================
#  FMO con penalización de sobredosis por órgano
#  Dados K ángulos fijos, calcula las intensidades de los beamlets que minimizan
#  la penalización cuadrática (QP convexo). Su óptimo es el f(B) del BAO.
#  Con w_over_organ igual para todos y normalize = 0 equivale a fmo.mod.
# =============================================================================

# --- Tamaños y conjuntos -----------------------------------------------------
param n_ptv integer >= 0;        # vóxeles del PTV
param n_oar integer >= 0;        # vóxeles de todos los OAR
param n_organs integer >= 1;     # número de OAR (ej. vejiga y recto)

set DIMLETS;                     # beamlets de los ángulos elegidos
set PTV_B := 0 .. n_ptv - 1;
set OAR_B := 0 .. n_oar - 1;
set ORGANS := 0 .. n_organs - 1;

# --- Parámetros de la penalización -------------------------------------------
param max_intensity > 0;         # intensidad máxima por beamlet
param w_under >= 0;              # peso de subdosis del PTV
param w_ptv_over >= 0 default 0; # peso de puntos calientes del PTV (0 = sin penalizar)
param oar_organ {OAR_B} integer >= 0, < n_organs;   # órgano de cada vóxel de OAR
param w_over_organ {ORGANS} >= 0;                   # peso de sobredosis por órgano

# normalize = 1: cada término se divide por el número de vóxeles de su estructura.
# Cambia la escala del objetivo (no comparable con normalize = 0).
param normalize binary default 0;
param organ_size {o in ORGANS} := card {b in OAR_B: oar_organ[b] = o};
param scale_ptv := if normalize = 1 and n_ptv > 0 then 1 / n_ptv else 1;
param scale_organ {o in ORGANS} :=
    if normalize = 1 and organ_size[o] > 0 then 1 / organ_size[o] else 1;

# --- Límites de dosis por vóxel ----------------------------------------------
param dmin {PTV_B} >= 0;                        # mínimo exigido al PTV
param dmax {OAR_B} >= 0;                        # máximo tolerado en OAR
param dmax_ptv {PTV_B} >= 0 default Infinity;   # techo del PTV

# --- Matriz de dosis dispersa: d[b,j] = dosis del beamlet j al vóxel b -------
set PTV_DOSE within {PTV_B, DIMLETS};
set OAR_DOSE within {OAR_B, DIMLETS};
param d_ptv {PTV_DOSE} >= 0;
param d_oar {OAR_DOSE} >= 0;

# --- Variables ---------------------------------------------------------------
var x {DIMLETS} >= 0, <= max_intensity;   # intensidad de cada beamlet
var u {PTV_B} >= 0;                       # subdosis del PTV (falta para dmin)
var v {OAR_B} >= 0;                       # sobredosis del OAR (exceso sobre dmax)
var w {PTV_B} >= 0;                       # exceso del PTV sobre dmax_ptv

# --- Objetivo: tres penalizaciones cuadráticas -------------------------------
minimize fmo_objective:
    w_under * scale_ptv * sum {b in PTV_B} u[b]^2
  + sum {b in OAR_B} w_over_organ[oar_organ[b]] * scale_organ[oar_organ[b]] * v[b]^2
  + w_ptv_over * scale_ptv * sum {b in PTV_B} w[b]^2;

# --- Restricciones (dosis del vóxel b = suma_j d[b,j] * x[j]) ----------------
subject to ptv_floor {b in PTV_B}:
    sum {(b,j) in PTV_DOSE} d_ptv[b,j] * x[j] + u[b] >= dmin[b];

subject to oar_ceiling {b in OAR_B}:
    sum {(b,j) in OAR_DOSE} d_oar[b,j] * x[j] - v[b] <= dmax[b];

subject to ptv_ceiling {b in PTV_B}:
    sum {(b,j) in PTV_DOSE} d_ptv[b,j] * x[j] - w[b] <= dmax_ptv[b];
