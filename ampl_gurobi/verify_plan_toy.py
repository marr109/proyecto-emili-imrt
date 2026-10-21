"""Independent dosimetric check of the optimal plan from fi_search_toy.py.
Recomputes D = A x and the holgura from their definitions and compares with AMPL."""
import os, tempfile
from amplpy import AMPL
import fi_search_toy as fs

ACTIVE = [1, 2]          # optimal angle set found by first improvement
K = len(ACTIVE)
TOL = 1e-3               # solver tolerance for barrier solutions (residuals ~1e-4)

def build_dat(active):
    K = len(active)
    d_ptv = " ".join(f"{p} {j} {fs.A_PTV[p][a]}" for p in range(2) for j, a in enumerate(active))
    d_oar = " ".join(f"{o} {j} {fs.A_OAR[o][a]}" for o in range(2) for j, a in enumerate(active))
    return f"""param n_ptv := 2; param n_oar := 2; param n_organs := 2;
set DIMLETS := {' '.join(str(j) for j in range(K))};
param max_intensity := {fs.M}; param w_under := {fs.W_UNDER}; param w_ptv_over := {fs.W_PTV_OVER};
param oar_organ := 0 0 1 1;
param w_over_organ := 0 {fs.W_OVER_ORGAN[0]} 1 {fs.W_OVER_ORGAN[1]};
param dmin := 0 {fs.DMIN} 1 {fs.DMIN}; param dmax := 0 {fs.DMAX} 1 {fs.DMAX};
param dmax_ptv := 0 {fs.DMAX_PTV} 1 {fs.DMAX_PTV};
set PTV_DOSE := {' '.join(f'({p},{j})' for p in range(2) for j in range(K))};
set OAR_DOSE := {' '.join(f'({o},{j})' for o in range(2) for j in range(K))};
param d_ptv := {d_ptv}; param d_oar := {d_oar};
"""

ampl = AMPL(); ampl.option["solver"] = "gurobi"; ampl.option["show_stats"] = 0
ampl.read(fs.MOD)
with tempfile.NamedTemporaryFile("w", suffix=".dat", delete=False) as f:
    f.write(build_dat(ACTIVE)); name = f.name
ampl.read_data(name); os.unlink(name)
ampl.solve()
assert ampl.solve_result == "solved", ampl.solve_result
x = [ampl.get_variable("x")[j].value() for j in range(K)]
u = [ampl.get_variable("u")[b].value() for b in range(2)]
v = [ampl.get_variable("v")[b].value() for b in range(2)]
obj_solver = ampl.get_objective("fmo_objective").value()
ampl.close()

# Recalculation from the definitions, no AMPL involved.
D_ptv = [sum(fs.A_PTV[p][a] * x[j] for j, a in enumerate(ACTIVE)) for p in range(2)]
D_oar = [sum(fs.A_OAR[o][a] * x[j] for j, a in enumerate(ACTIVE)) for o in range(2)]
u_chk = [max(0.0, fs.DMIN - d) for d in D_ptv]
v_chk = [max(0.0, d - fs.DMAX) for d in D_oar]
w_chk = [max(0.0, d - fs.DMAX_PTV) for d in D_ptv]
obj_chk = (fs.W_UNDER * sum(t * t for t in u_chk)
           + sum(fs.W_OVER_ORGAN[o] * v_chk[o] ** 2 for o in range(2))
           + fs.W_PTV_OVER * sum(t * t for t in w_chk))

ok = True
print(f"x = {[round(t, 4) for t in x]} for angles {ACTIVE}")
for name_, s, c in [("PTV 0", u[0], u_chk[0]), ("PTV 1", u[1], u_chk[1]),
                    ("bladder", v[0], v_chk[0]), ("rectum", v[1], v_chk[1])]:
    good = abs(s - c) < TOL; ok &= good
    print(f"  {name_:<8} slack solver={s:9.4f}  recalc={c:9.4f}  {'OK' if good else 'MISMATCH'}")
good = abs(obj_solver - obj_chk) < TOL; ok &= good
print(f"  objective solver={obj_solver:.6f}  recalc={obj_chk:.6f}  {'OK' if good else 'MISMATCH'}")
print(f"\nDosimetric review (provisional dmin/dmax = {fs.DMIN}/{fs.DMAX} Gy-equiv):")
print(f"  PTV min dose = {min(D_ptv):.2f}  -> {'covered' if min(D_ptv) >= fs.DMIN - TOL else 'UNDERDOSE by %.2f' % (fs.DMIN - min(D_ptv))}")
print(f"  OAR max dose = {max(D_oar):.2f}  -> {'within Dmax' if max(D_oar) <= fs.DMAX + TOL else 'ABOVE Dmax by %.2f' % (max(D_oar) - fs.DMAX)}")
clinical_ok = min(D_ptv) >= fs.DMIN - TOL and max(D_oar) <= fs.DMAX + TOL
print("\nSOLVER CONSISTENCY:", "PASS" if ok else "FAIL")
print("CLINICAL ACCEPTABILITY:", "PASS" if clinical_ok else "FAIL (plan violates provisional limits)")
