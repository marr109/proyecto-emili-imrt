param n_ptv integer >= 0;   
param n_oar integer >= 0;    
param n_organs integer >= 1;    

set DIMLETS;                        

set PTV_B := 0 .. n_ptv - 1;         
set OAR_B := 0 .. n_oar - 1;         
set ORGANS := 0 .. n_organs - 1;     

param max_intensity > 0;           
param w_under >= 0;                  
param w_ptv_over >= 0 default 0;    
                                    
param oar_organ {OAR_B} integer >= 0, < n_organs;
param w_over_organ {ORGANS} >= 0;

param normalize binary default 0;
param organ_size {o in ORGANS} := card {b in OAR_B: oar_organ[b] = o};  
param scale_ptv := if normalize = 1 and n_ptv > 0 then 1 / n_ptv else 1;
param scale_organ {o in ORGANS} :=
    if normalize = 1 and organ_size[o] > 0 then 1 / organ_size[o] else 1;

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
    w_under * scale_ptv * sum {b in PTV_B} u[b]^2
  + sum {b in OAR_B} w_over_organ[oar_organ[b]] * scale_organ[oar_organ[b]] * v[b]^2
  + w_ptv_over * scale_ptv * sum {b in PTV_B} w[b]^2;

subject to ptv_floor {b in PTV_B}:
    sum {(b,j) in PTV_DOSE} d_ptv[b,j] * x[j] + u[b] >= dmin[b];

subject to oar_ceiling {b in OAR_B}:
    sum {(b,j) in OAR_DOSE} d_oar[b,j] * x[j] - v[b] <= dmax[b];

subject to ptv_ceiling {b in PTV_B}:
    sum {(b,j) in PTV_DOSE} d_ptv[b,j] * x[j] - w[b] <= dmax_ptv[b];
