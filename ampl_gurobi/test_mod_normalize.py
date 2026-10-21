"""Check fmo_organs.mod on the toy instance: normalize=0 vs normalize=1 against a hand calculation."""
import os
from amplpy import AMPL
HERE = os.path.dirname(os.path.abspath(__file__))
def run(norm):
    a = AMPL(); a.option["solver"] = "gurobi"; a.option["show_stats"] = 0
    a.read(os.path.join(HERE, "fmo_organs.mod"))
    a.read_data(os.path.join(HERE, "fmo_organs_example.dat"))
    a.get_parameter("normalize").set(norm)
    a.solve()
    r = (a.solve_result, a.get_objective("fmo_objective").value(),
         [a.get_variable(n).to_dict() for n in ("u", "v", "w")])
    a.close(); return r
r0, r1 = run(0), run(1)
print("normalize=0:", r0[0], round(r0[1], 9))
print("normalize=1:", r1[0], round(r1[1], 9))
# toy sizes: n_ptv=2, organ_size = {0:2, 1:1}; weights 1.0 / [0.5, 0.8] / 0.3
u, v, w = r1[2]
chk = (1.0 / 2) * sum(t * t for t in u.values()) \
    + 0.5 / 2 * (v[0] ** 2 + v[1] ** 2) + 0.8 / 1 * v[2] ** 2 \
    + 0.3 / 2 * sum(t * t for t in w.values())
print("normalize=1 recomputed from solution:", round(chk, 9), "match:", abs(chk - r1[1]) < 1e-6)
