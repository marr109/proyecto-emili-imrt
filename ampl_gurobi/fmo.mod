# FMO QP -- mirrors imrt/imrt_fmo.cpp::ImrtFmoSolver::solve() exactly.
# Variables: x = active beamlet intensities, u = PTV underdose slack,
#            v = OAR overdose slack, w = PTV overdose slack.
# Setting w_ptv_over = 0 and dmax_ptv = Infinity reproduces the
# "w disabled" behavior of the C++ solver without needing conditional
# constraint activation (the ptv_ceiling row becomes non-binding and w
# has zero objective weight, so it is driven to 0 at no cost).

param n_ptv integer >= 0;
param n_oar integer >= 0;

set DIMLETS;
set PTV_B := 0 .. n_ptv - 1;
set OAR_B := 0 .. n_oar - 1;

param max_intensity > 0;
param w_under >= 0;
param w_over  >= 0;
param w_ptv_over >= 0 default 0;

param dmin {PTV_B} >= 0;
param dmax {OAR_B} >= 0;
param dmax_ptv {PTV_B} >= 0 default Infinity;

set PTV_DOSE within {PTV_B, DIMLETS};
set OAR_DOSE within {OAR_B, DIMLETS};
param d_ptv {PTV_DOSE} >= 0;
param d_oar {OAR_DOSE} >= 0;

var x {DIMLETS} >= 0, <= max_intensity;
var u {PTV_B} >= 0;
var v {OAR_B} >= 0;
var w {PTV_B} >= 0;

minimize fmo_objective:
    w_under    * sum {b in PTV_B} u[b]^2
  + w_over     * sum {b in OAR_B} v[b]^2
  + w_ptv_over * sum {b in PTV_B} w[b]^2;

subject to ptv_floor {b in PTV_B}:
    sum {(b,j) in PTV_DOSE} d_ptv[b,j] * x[j] + u[b] >= dmin[b];

subject to oar_ceiling {b in OAR_B}:
    sum {(b,j) in OAR_DOSE} d_oar[b,j] * x[j] - v[b] <= dmax[b];

subject to ptv_ceiling {b in PTV_B}:
    sum {(b,j) in PTV_DOSE} d_ptv[b,j] * x[j] - w[b] <= dmax_ptv[b];
