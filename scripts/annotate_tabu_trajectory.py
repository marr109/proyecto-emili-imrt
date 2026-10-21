#!/usr/bin/env python3
"""Agrega columnas 'ronda' y 'tenure' a un trajectory.csv de tabú (Best o First),
reconstruyendo fielmente la lógica real de BestTabuSearch/FirstTabuSearch +
AdaptiveBaoTabuMemory (emilibase.cpp:937-1010, imrt/imrt_bao.cpp:921-984).

Estas columnas NO las escribe el binario -- son derivadas, reconstruidas por
este script a partir de los valores ya evaluados en el CSV, replicando la
misma secuencia de decisiones (bestSoFar/incumbent/tabu_check/forbid/resize)
que el algoritmo real. Verificado contra la salida real ("Found solution")
de varias semillas antes de usarse (ver docs/analisis-tabu-adaptativo.md).

Uso: python3 scripts/annotate_tabu_trajectory.py <trajectory.csv> <first|best> [--out OUT.csv]
"""
import csv
import sys
import argparse


def simulate(rows, strategy, tenure_min=3, tenure_max=8):
    """Replica BestTabuSearch/FirstTabuSearch + AdaptiveBaoTabuMemory.

    Devuelve una lista de dicts por fila (en el mismo orden que `rows`,
    empezando en la fila 2 del CSV -- la fila 1 es la solución inicial, no
    un candidato) con: ronda, tenure_en_esa_ronda, prohibido, aceptado.
    """
    tenure = (tenure_min + tenure_max) // 2
    memory = []  # buffer circular simulado como lista, más reciente al final
    since_last_revisit = 0

    def tabu_check(angles):
        return angles in memory

    def forbid(angles):
        nonlocal tenure, since_last_revisit, memory
        if tabu_check(angles):
            if tenure < tenure_max:
                tenure += 1
            since_last_revisit = 0
        else:
            since_last_revisit += 1
            if since_last_revisit >= 2 * tenure and tenure > tenure_min:
                tenure -= 1
                since_last_revisit = 0
        memory.append(angles)
        if len(memory) > tenure:
            memory = memory[-tenure:]

    best_so_far_val = float(rows[0]['objective'])
    best_so_far_ang = rows[0]['angles_deg']
    incumbent_val = best_so_far_val
    incumbent_ang = best_so_far_ang

    annotations = []  # una entrada por fila de candidato (rows[1:])
    idx = 1
    round_num = 0

    while idx < len(rows):
        round_num += 1
        if incumbent_val < best_so_far_val:
            best_so_far_val, best_so_far_ang = incumbent_val, incumbent_ang

        tenure_this_round = tenure

        # El incumbente se fuerza al primer candidato de la ronda (línea real
        # del código: *incumbent = *ithSolution antes de entrar al for).
        first_row = rows[idx]
        incumbent_val = float(first_row['objective'])
        incumbent_ang = first_row['angles_deg']
        annotations.append({'ronda': round_num, 'tenure': tenure_this_round,
                             'prohibido': tabu_check(incumbent_ang), 'aceptado': True})
        idx += 1

        # Candidatos siguientes de la misma ronda (ya van 1 de 8 contado arriba).
        candidates_seen_this_round = 1
        NEIGHBORHOOD_SIZE = 8  # 2*K, K=4 ángulos activos
        while idx < len(rows) and candidates_seen_this_round < NEIGHBORHOOD_SIZE:
            row = rows[idx]
            val = float(row['objective'])
            ang = row['angles_deg']
            is_tabu = tabu_check(ang)
            beats_best = val < best_so_far_val
            eligible = (not is_tabu) or beats_best
            improves = incumbent_val > val
            accepted = improves and eligible
            annotations.append({'ronda': round_num, 'tenure': tenure_this_round,
                                 'prohibido': is_tabu, 'aceptado': accepted})
            idx += 1
            candidates_seen_this_round += 1
            if accepted:
                incumbent_val, incumbent_ang = val, ang
                if strategy == 'first':
                    break  # First corta apenas encuentra un candidato elegible que mejora
            # Best sigue hasta agotar los 8 candidatos de la ronda (sin cortar).

        tenure_before_forbid = tenure
        forbid(incumbent_ang)

    return annotations, best_so_far_ang, best_so_far_val


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('csv_path')
    ap.add_argument('strategy', choices=['first', 'best'])
    ap.add_argument('--out', default=None)
    ap.add_argument('--tenure-min', type=int, default=3)
    ap.add_argument('--tenure-max', type=int, default=8)
    args = ap.parse_args()

    with open(args.csv_path) as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        rows = list(reader)

    annotations, best_ang, best_val = simulate(rows, args.strategy, args.tenure_min, args.tenure_max)

    # Se descartan las columnas de desglose por órgano (ptv_*/oar_*) -- no
    # aportan al análisis de mecánica de tabú, solo agregan ruido visual.
    keep_fields = [f for f in fieldnames if not (f.startswith('ptv_') or f.startswith('oar_'))]

    out_path = args.out or args.csv_path.replace('.csv', '_annotated.csv')
    new_fields = keep_fields + ['ronda', 'tenure', 'prohibido', 'aceptado']
    with open(out_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=new_fields, extrasaction='ignore')
        writer.writeheader()
        # Primera fila (solución inicial) no tiene ronda/tenure -- se marca vacía.
        first = dict(rows[0])
        first.update({'ronda': 0, 'tenure': '', 'prohibido': '', 'aceptado': ''})
        writer.writerow(first)
        for row, ann in zip(rows[1:], annotations):
            merged = dict(row)
            merged.update(ann)
            writer.writerow(merged)

    print(f"escrito: {out_path}")
    print(f"bestSoFar final (reconstruido): {best_ang}  f={best_val:.2f}")
    print(f"rondas totales: {annotations[-1]['ronda']}")
    print(f"tenure final: {annotations[-1]['tenure']}")


if __name__ == '__main__':
    main()
