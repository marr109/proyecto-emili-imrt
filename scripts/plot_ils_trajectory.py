#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Grafica una trayectoria de ILS resaltando en qué evaluaciones ocurrió una
perturbación (AngleShiftMultiPerturbation) en vez de un paso normal de
búsqueda local (AngleShiftNeighborhood).

Distingue ambas por la magnitud del cambio: la búsqueda local siempre mueve
un ángulo por exactamente ±step; la perturbación lo mueve por una magnitud
aleatoria en (step, 2*step) — nunca step exacto (ver
AngleShiftMultiPerturbation, imrt/imrt_bao.cpp). Comparando cada fila contra
la anterior (mismo slot, cuánto cambió), esa diferencia de magnitud es
suficiente para clasificar cada evaluación sin necesitar ninguna columna
extra en el CSV.

Cada perturbación se reporta AGRUPADA: `prangshift <mag> <numSteps>` mueve
numSteps ángulos y el evaluador registra un eval por movimiento, pero eso es
UNA sola perturbación. Tanto el CSV (--perturbations-csv) como el gráfico
tratan la ráfaga completa como un único salto: s_p aceptada -> solución
perturbada. El detalle movimiento a movimiento sigue disponible con
--moves-csv para auditar la clasificación.

CLI example:
  python scripts/plot_ils_trajectory.py \
      --csv experiments/ils/nangshift10/grid5/data/first/seed01/trajectory.csv \
      --step 10 \
      --output experiments/ils/nangshift10/grid5/data/first/seed01/ils_perturbations.png
"""

import argparse
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_trajectory import read_trajectory, running_best


def circular_delta(a, b, n_angles=360):
    d = abs(a - b) % n_angles
    return min(d, n_angles - d)


def detect_perturbations(evals, objectives, angle_sets, step, reject_repeated=True):
    """Simula la semántica real de ILS con First Improvement — no solo el
    rebase dentro de una ronda de búsqueda local, sino también la decisión
    de aceptar/rechazar el resultado de cada iteración externa
    (BaoImproveAccept, imrt_bao.cpp), porque con `rejectrepeated` activo el
    punto desde el que se perturba en la iteración siguiente puede REVERTIR
    a un punto anterior si el resultado de la ronda no mejoró o repitió una
    solución ya visitada — y si no se modela ese revert, la primera
    perturbación de esa ráfaga parece (incorrectamente) cambiar varios
    slots a la vez en vez de uno solo.

    Máquina de estados con dos fases:
      "search": recorre vecinos ±step de `base` (ronda de búsqueda local
        vigente), rebasando cada vez que encuentra una mejora — igual que
        FirstImprovementSearch.
      "perturbing": encadena movimientos de magnitud != step a partir de
        `pert_chain` (arranca en s_p, el punto ACEPTADO por el ILS, no en
        el último base de búsqueda local si esa ronda fue rechazada).

    La transición search -> perturbing resuelve el accept/reject de la
    ronda que acaba de terminar contra s_p (con el mismo criterio de
    BaoImproveAccept, incluyendo el registro de visitados si
    reject_repeated=True) ANTES de clasificar la fila que disparó la
    transición como el primer paso de la perturbación.

    Solo válido para `first` (Best no rebasa a mitad de ronda).
    """
    s_p, s_p_obj, s_p_eval = angle_sets[0], objectives[0], evals[0]
    visited = set()
    base, base_obj, base_eval = s_p, s_p_obj, s_p_eval
    phase = "search"
    pert_chain, pert_chain_obj, pert_chain_eval = None, None, None
    records = []

    for i in range(1, len(angle_sets)):
        cur, cur_obj = angle_sets[i], objectives[i]

        if phase == "search":
            diffs = [j for j in range(len(cur)) if cur[j] != base[j]]
            if len(diffs) == 1 and circular_delta(cur[diffs[0]], base[diffs[0]]) == step:
                if cur_obj < base_obj:
                    base, base_obj, base_eval = cur, cur_obj, evals[i]
                continue

            # Ronda terminada: resuelve accept/reject de `base` (candidato)
            # contra s_p, igual que BaoImproveAccept::accept.
            candidate, candidate_obj, candidate_eval = base, base_obj, base_eval
            visited.add(tuple(s_p))
            if candidate_obj < s_p_obj and not (reject_repeated and tuple(candidate) in visited):
                visited.add(tuple(candidate))
                s_p, s_p_obj, s_p_eval = candidate, candidate_obj, candidate_eval
            pert_chain, pert_chain_obj, pert_chain_eval = s_p, s_p_obj, s_p_eval
            phase = "perturbing"
            # sigue abajo: clasifica esta misma fila ya en fase "perturbing"

        # phase == "perturbing"
        diffs = [j for j in range(len(cur)) if cur[j] != pert_chain[j]]
        if len(diffs) == 1 and circular_delta(cur[diffs[0]], pert_chain[diffs[0]]) == step:
            # La búsqueda local ya retomó: esta fila es un vecino normal de
            # la base perturbada — vuelve a fase "search" y la reprocesa.
            phase = "search"
            base, base_obj, base_eval = pert_chain, pert_chain_obj, pert_chain_eval
            if cur_obj < base_obj:
                base, base_obj, base_eval = cur, cur_obj, evals[i]
            continue

        slot = diffs[0] if len(diffs) == 1 else None
        records.append({
            "eval_antes": pert_chain_eval,
            "eval_despues": evals[i],
            "angles_antes": ";".join(str(a) for a in pert_chain),
            "angles_despues": ";".join(str(a) for a in cur),
            "slots_distintos": len(diffs),
            "slot_movido": slot if slot is not None else "",
            "delta_deg": circular_delta(cur[slot], pert_chain[slot]) if slot is not None else "",
            "obj_antes": pert_chain_obj,
            "obj_despues": cur_obj,
        })
        pert_chain, pert_chain_obj, pert_chain_eval = cur, cur_obj, evals[i]

    return records


def group_bursts(records):
    """Numera cada perturbación detectada dentro de su ráfaga consecutiva
    (evals correlativos = misma perturbación de numSteps movimientos)."""
    burst_id = 0
    step_in_burst = 0
    prev_eval = None
    for r in records:
        if prev_eval is not None and r["eval_despues"] == prev_eval + 1:
            step_in_burst += 1
        else:
            burst_id += 1
            step_in_burst = 1
        r["rafaga"] = burst_id
        r["paso_en_rafaga"] = step_in_burst
        prev_eval = r["eval_despues"]
    return records


def collapse_bursts(records, n_angles=360):
    """Colapsa cada rafaga a UNA sola perturbacion: el cambio completo desde
    la solucion aceptada por el ILS (s_p, antes del primer movimiento de la
    rafaga) hasta la solucion ya perturbada (despues del ultimo movimiento).

    Una perturbacion de numSteps movimientos aparece en la trayectoria como
    numSteps filas consecutivas porque el evaluador registra cada movimiento
    por separado, pero conceptualmente es UN solo salto. Los slots se
    comparan extremo contra extremo, asi que un slot movido y devuelto a su
    valor original no cuenta como movido.
    """
    bursts = []
    for r in records:
        if bursts and bursts[-1]["rafaga"] == r["rafaga"]:
            b = bursts[-1]
            b["eval_despues"] = r["eval_despues"]
            b["angles_despues"] = r["angles_despues"]
            b["obj_despues"] = r["obj_despues"]
            b["n_movimientos"] += 1
        else:
            bursts.append({
                "rafaga": r["rafaga"],
                "n_movimientos": 1,
                "eval_antes": r["eval_antes"],
                # eval_antes es donde se ACEPTÓ s_p, que puede quedar muy
                # atrás en la ronda; el primer movimiento de la perturbación
                # marca dónde empieza realmente la ráfaga.
                "eval_primer_movimiento": r["eval_despues"],
                "eval_despues": r["eval_despues"],
                "angles_antes": r["angles_antes"],
                "angles_despues": r["angles_despues"],
                "obj_antes": r["obj_antes"],
                "obj_despues": r["obj_despues"],
            })

    for b in bursts:
        antes = [int(a) for a in b["angles_antes"].split(";")]
        despues = [int(a) for a in b["angles_despues"].split(";")]
        moved = [j for j in range(len(antes)) if antes[j] != despues[j]]
        b["slots_movidos"] = ";".join(str(j) for j in moved)
        b["n_slots_movidos"] = len(moved)
        b["deltas_deg"] = ";".join(
            str(circular_delta(despues[j], antes[j], n_angles)) for j in moved)
        b["delta_deg_total"] = sum(
            circular_delta(despues[j], antes[j], n_angles) for j in moved)
        b["delta_obj"] = b["obj_despues"] - b["obj_antes"]
    return bursts


def write_csv(path, rows, fieldnames):
    out_dir = os.path.dirname(path) or "."
    os.makedirs(out_dir, exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--csv", required=True)
    p.add_argument("--step", type=int, required=True,
                    help="step de nangshift usado en esta corrida (10 o 5)")
    p.add_argument("--output", required=True)
    p.add_argument("--title", default=None)
    p.add_argument("--perturbations-csv", default=None,
                    help="Si se pasa, escribe UNA fila por perturbación completa: "
                         "la ráfaga de movimientos colapsada al cambio total "
                         "(s_p aceptada -> solución perturbada).")
    p.add_argument("--moves-csv", default=None,
                    help="Si se pasa, escribe el detalle movimiento a movimiento "
                         "dentro de cada ráfaga (para auditar la clasificación).")
    p.add_argument("--no-reject-repeated", action="store_true",
                    help="Pasar si la corrida NO usó `baoimprove rejectrepeated` "
                         "(por defecto se asume que sí, que es lo que usan las "
                         "corridas actuales de ILS).")
    args = p.parse_args(argv)

    evals, objectives, angle_sets = read_trajectory(args.csv)
    if not evals:
        raise SystemExit("error: empty trajectory CSV: {0}".format(args.csv))

    best = running_best(objectives)
    records = group_bursts(detect_perturbations(
        evals, objectives, angle_sets, args.step,
        reject_repeated=not args.no_reject_repeated))
    bursts = collapse_bursts(records)

    # El punto marcado es el FINAL de cada perturbación completa: el estado
    # ya perturbado desde el que arranca la próxima búsqueda local.
    obj_by_eval = dict(zip(evals, objectives))
    pert_xs = [b["eval_despues"] for b in bursts]
    pert_objs = [obj_by_eval[e] for e in pert_xs]

    if args.perturbations_csv:
        write_csv(args.perturbations_csv, bursts,
                  ["rafaga", "n_movimientos", "eval_antes",
                   "eval_primer_movimiento", "eval_despues",
                   "angles_antes", "angles_despues", "n_slots_movidos",
                   "slots_movidos", "deltas_deg", "delta_deg_total",
                   "obj_antes", "obj_despues", "delta_obj"])
        print("perturbations CSV written: {0} ({1} perturbaciones completas)".format(
            args.perturbations_csv, len(bursts)))

    if args.moves_csv:
        write_csv(args.moves_csv, records,
                  ["rafaga", "paso_en_rafaga", "eval_antes", "eval_despues",
                   "angles_antes", "angles_despues", "slots_distintos",
                   "slot_movido", "delta_deg", "obj_antes", "obj_despues"])
        print("moves CSV written: {0} ({1} movimientos)".format(
            args.moves_csv, len(records)))

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(evals, objectives, color="#9aa5b1", linewidth=0.7, alpha=0.5, label="evaluado")
    ax.step(evals, best, where="post", color="#2f6f9f", linewidth=2, label="mejor hasta el momento")

    # Cada perturbación se dibuja como UNA banda (desde la solución aceptada
    # hasta la solución ya perturbada), no como una marca por movimiento.
    for b in bursts:
        ax.axvspan(b["eval_primer_movimiento"] - 0.5, b["eval_despues"] + 0.5,
                   color="#d1495b", alpha=0.35, linewidth=0, zorder=1)
    if pert_xs:
        ax.scatter(pert_xs, pert_objs, color="#d1495b", s=28, zorder=3,
                   label="perturbación completa ({0})".format(len(bursts)))

    ax.set_xlabel("evaluación")
    ax.set_ylabel("FMO objective")
    ax.set_title(args.title or "Trayectoria ILS — perturbaciones resaltadas")
    ax.legend(loc="upper right")
    ax.grid(alpha=0.3)

    fig.tight_layout()
    out_dir = os.path.dirname(args.output) or "."
    os.makedirs(out_dir, exist_ok=True)
    fig.savefig(args.output, dpi=150)
    plt.close(fig)
    print("ILS trajectory plot written: {0} ({1} perturbaciones completas, "
          "{2} movimientos, de {3} evals)".format(
              args.output, len(bursts), len(records), len(evals)))


if __name__ == "__main__":
    sys.exit(main())
