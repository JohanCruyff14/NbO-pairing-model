"""DFT-energy-driven finite reduced pairing prototype; no bulk quadrature claim.

Requires numpy and scipy. Run from any directory. All output is local.
Pair level p comprises (band,k,up) and (band,-k,down), with disjoint
fermionic orbitals between levels. Both k orientations are retained.
"""
from pathlib import Path
from itertools import combinations
import csv
import hashlib
import json
import math

import numpy as np
from scipy.sparse import coo_matrix, diags, save_npz
from scipy.sparse.linalg import eigsh

BASE = Path(__file__).resolve().parent
OUT = BASE / "results"


def write_csv(name, rows, fields=None):
    rows = list(rows)
    with (OUT / name).open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields or list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def wrap(k):
    return tuple(round(float(x) % 1.0, 10) % 1.0 for x in k)


def choose_levels(config):
    with (BASE / "data/candidate_bands_25_27.csv").open(
        encoding="utf-8-sig", newline=""
    ) as stream:
        rows = list(csv.DictReader(stream))
    levels = []
    for band in config["candidate_bands"]:
        for side in (-1, 1):
            options = [r for r in rows if int(r["band_index"]) == band
                       and side * float(r["xi_eV"]) > 0]
            if not options:
                raise ValueError(f"No state on side {side} for band {band}")
            row = min(options, key=lambda r: (abs(float(r["xi_eV"])), int(r["k_index"])))
            k = wrap([row[f"k{i}"] for i in (1, 2, 3)])
            minus_k = wrap([-x for x in k])
            # At a TR-invariant momentum these are the same pair, so include once.
            orientations = [k] if k == minus_k else [k, minus_k]
            for orientation, momentum in enumerate(orientations):
                partner = wrap([-x for x in momentum])
                item = dict(
                    p=len(levels), band_index=band, selected_side="below" if side < 0 else "above",
                    source_k_index=int(row["k_index"]), orientation=orientation,
                    k1=momentum[0], k2=momentum[1], k3=momentum[2],
                    partner_k1=partner[0], partner_k2=partner[1], partner_k3=partner[2],
                    xi_eV=float(row["xi_eV"]), energy_eV=float(row["energy_eV"]),
                    source_irreducible_weight=float(row["weight"]),
                    source_occupation_per_spin=float(row["occupation_per_spin"]),
                )
                levels.append(item)
    # Spin orbitals must be disjoint, including at Brillouin-zone boundaries.
    orbitals = set()
    for level in levels:
        for prefix, spin in (("", "up"), ("partner_", "down")):
            key = (level["band_index"], tuple(level[f"{prefix}k{i}"] for i in (1, 2, 3)), spin)
            assert key not in orbitals, f"Duplicate spin orbital: {key}"
            orbitals.add(key)
    return levels


def basis_and_pair_operators(m, pair_number):
    basis = sorted(sum(1 << p for p in occupied)
                   for occupied in combinations(range(m), pair_number))
    indices = {state: i for i, state in enumerate(basis)}
    occupation = np.array([[(state >> p) & 1 for p in range(m)] for state in basis], dtype=float)
    # A_pq = b_p^dagger b_q; even-parity pair operators commute on distinct levels.
    operators = {}
    for p in range(m):
        for q in range(m):
            if p == q:
                operators[p, q] = diags(occupation[:, p], format="csr")
                continue
            rr, cc = [], []
            for col, state in enumerate(basis):
                if ((state >> q) & 1) and not ((state >> p) & 1):
                    rr.append(indices[state ^ (1 << p) ^ (1 << q)])
                    cc.append(col)
            operators[p, q] = coo_matrix((np.ones(len(rr)), (rr, cc)),
                                         shape=(len(basis), len(basis))).tocsr()
    return basis, occupation, operators


def hamiltonian(xi, coupling, occupation, operators):
    h = diags(2 * occupation @ xi, format="csr")
    for (p, q), op in operators.items():
        h = h - coupling[p, q] * op
    h.eliminate_zeros()
    return h


def qubit_terms(xi, coupling):
    # |0>=empty, |1>=pair; Z|0>=|0>. Qubit p is bit p (least-significant first).
    m = len(xi)
    terms = [dict(coefficient_eV=float(sum(xi - np.diag(coupling) / 2)), pauli="I", qubits=[])]
    for p in range(m):
        terms.append(dict(coefficient_eV=float(-xi[p] + coupling[p, p] / 2), pauli="Z", qubits=[p]))
    for p in range(m):
        for q in range(p + 1, m):
            for pauli in ("XX", "YY"):
                terms.append(dict(coefficient_eV=float(-coupling[p, q] / 2), pauli=pauli, qubits=[p, q]))
    return terms


def check_qubit_action(terms, h, basis):
    # Independent Pauli action, including cancellation of number-changing XX/YY.
    index = {state: i for i, state in enumerate(basis)}
    rr, cc, vv = [], [], []
    leakage = 0.0
    for col, state in enumerate(basis):
        outputs = {}
        for term in terms:
            target, coefficient = state, complex(term["coefficient_eV"])
            if term["pauli"] != "I":
                for pauli, p in zip(term["pauli"], term["qubits"]):
                    bit = (state >> p) & 1
                    if pauli == "Z":
                        coefficient *= 1 - 2 * bit
                    else:
                        target ^= 1 << p
                        if pauli == "Y":
                            coefficient *= 1j * (1 - 2 * bit)
            outputs[target] = outputs.get(target, 0) + coefficient
        for target, coefficient in outputs.items():
            if target not in index:
                leakage = max(leakage, abs(coefficient))
            elif coefficient:
                rr.append(index[target]); cc.append(col); vv.append(coefficient)
    restricted = coo_matrix((vv, (rr, cc)), shape=h.shape).tocsr()
    difference = restricted - h
    error = float(np.max(np.abs(difference.data))) if difference.nnz else 0.0
    assert leakage < 1e-12 and error < 1e-12
    return error, leakage


def main():
    config = json.loads((BASE / "model_config.json").read_text(encoding="utf-8"))
    assert config["selection"] == "nearest_strictly_below_and_above_EF_per_band"
    assert config["include_both_pair_orientations"] is True
    assert config["pair_number"] == "zero_coupling_filling"
    levels = choose_levels(config)
    xi = np.array([r["xi_eV"] for r in levels])
    m = len(levels)
    n_pairs = int(np.count_nonzero(xi < 0))
    basis, occupation, operators = basis_and_pair_operators(m, n_pairs)
    assert len(basis) == math.comb(m, n_pairs)
    fs_state = sum(1 << p for p in range(m) if xi[p] < 0)
    fs_index = basis.index(fs_state)
    reference_strength = float(config["reference_strength_eV"])
    strengths = sorted(set([0.0, reference_strength] + list(config["scan_strengths_eV"])))
    assert all(float(s) >= 0 for s in strengths)
    OUT.mkdir(exist_ok=True)
    write_csv("active_levels.csv", levels)
    write_csv("basis.csv", [dict(index=i, bitmask=state,
              ket_pMminus1_to_p0=format(state, f"0{m}b")) for i, state in enumerate(basis)])
    scan, reference = [], None
    for strength in strengths:
        # gamma/M is a convention of THIS finite model, not a DFT-derived BZ integral.
        coupling = np.full((m, m), float(strength) / m)
        h = hamiltonian(xi, coupling, occupation, operators)
        assert (h - h.T).nnz == 0
        if strength == 0:
            psi = np.zeros(len(basis)); psi[fs_index] = 1.0
            energy = float(h[fs_index, fs_index])
            assert abs(energy - 2 * sum(xi[xi < 0])) < 1e-12
        else:
            values, vectors = eigsh(h, k=1, which="SA", tol=1e-12,
                                   v0=np.ones(len(basis)) / np.sqrt(len(basis)))
            energy, psi = float(values[0]), vectors[:, 0]
            if psi[np.argmax(abs(psi))] < 0:
                psi = -psi
        residual = float(np.linalg.norm(h @ psi - energy * psi))
        assert residual < 1e-9
        correlation = np.array([[float(psi @ (operators[p, q] @ psi))
                                 for q in range(m)] for p in range(m)])
        assert abs(np.trace(correlation) - n_pairs) < 1e-9
        assert np.max(abs(correlation - correlation.T)) < 1e-10
        assert np.linalg.eigvalsh(correlation)[0] > -1e-9
        e_fs = float(h[fs_index, fs_index])
        item = dict(strength_eV=float(strength), G_each_element_eV=float(strength) / m,
                    ground_energy_eV=energy, reference_FS_energy_same_H_eV=e_fs,
                    correlation_energy_eV=energy - e_fs,
                    pair_matrix_largest_eigenvalue=float(np.linalg.eigvalsh(correlation)[-1]),
                    offdiagonal_pair_sum=float(correlation.sum() - np.trace(correlation)),
                    reference_FS_probability=float(psi[fs_index] ** 2), residual_norm_eV=residual)
        assert energy <= e_fs + 1e-10
        scan.append(item)
        if strength == reference_strength:
            reference = (h, coupling, psi, correlation, item)
    write_csv("strength_scan.csv", scan)
    h, coupling, psi, correlation, result = reference
    terms = qubit_terms(xi, coupling)
    pauli_error, leakage = check_qubit_action(terms, h, basis)
    save_npz(OUT / "hamiltonian_fixed_pairs.npz", h)
    write_csv("G_matrix_eV.csv", [dict(p=p, **{f"q{q}": coupling[p, q] for q in range(m)})
                                  for p in range(m)])
    write_csv("pair_correlation.csv", [dict(p=p, **{f"q{q}": correlation[p, q] for q in range(m)})
                                       for p in range(m)])
    write_csv("ground_state.csv", [dict(basis_index=i, bitmask=state,
              amplitude=psi[i], probability=psi[i] ** 2) for i, state in enumerate(basis)])
    write_csv("occupations.csv", [dict(p=p, band_index=levels[p]["band_index"], xi_eV=xi[p],
              pair_occupation=correlation[p, p], electron_occupation=2 * correlation[p, p])
              for p in range(m)])
    (OUT / "qubit_hamiltonian.json").write_text(json.dumps(dict(
        n_qubits=m, basis_convention="qubit p = bit p; 0 empty, 1 doubly occupied",
        target_pair_number=n_pairs, coefficient_unit="eV", terms=terms), indent=2) + "\n")
    record = dict(
        model="DFT-energy-driven phenomenological reduced pairing prototype v1",
        schema_version=1,
        numpy_version=np.__version__,
        config=config, n_pair_levels=m, n_spin_orbitals=2*m,
        n_pairs=n_pairs, n_active_electrons=2*n_pairs, dimension=len(basis),
        frozen_electron_contribution_included=False,
        irreducible_weights_used_in_hamiltonian=False,
        full_crystal_point_symmetry_preserved=False,
        interaction_is_first_principles=False,
        hamiltonian="H=2 sum_p xi_p N_p - sum_pq G_pq bdag_p b_q; G_pq=gamma/M incl. p=q",
        electron_number_rule="fill all selected negative-xi pair levels at gamma=0; keep this number",
        checks=dict(pauli_restriction_max_error_eV=pauli_error, number_sector_leakage=leakage,
                    ground_state_residual_norm_eV=result["residual_norm_eV"]),
        reference_result=result,
        inputs_sha256={name: hashlib.sha256((BASE / name).read_bytes()).hexdigest()
                       for name in ("model_config.json", "data/candidate_bands_25_27.csv")},
    )
    (OUT / "model_record.json").write_text(json.dumps(record, indent=2) + "\n")
    text = ["# 第一版有限配对模型：已计算结果", "",
            "DFT 只提供单粒子能量；相互作用为参数化吸引。不是 NbO 体相超导预测。", "",
            f"- 配对能级：{m}；自旋轨道：{2*m}；配对编码量子比特：{m}。",
            f"- 固定 {n_pairs} 对电子；矩阵维度 {len(basis)}。",
            "- 选态：每带取最近负能态和正能态，补齐两个 k 配对方向。",
            "- 粒子数：按这些选定能级的零作用填充确定，不等于体相电子数下折叠。",
            "- G_pq = gamma/M，包括对角项。gamma 是有限模型参数，单位 eV。", "",
            "| band | source k | xi / meV |", "| --- | --- | --- |"]
    for level in levels:
        if level["orientation"] == 0:
            text.append(f"| {level['band_index']} | ({level['k1']}, {level['k2']}, {level['k3']}) | {1000*level['xi_eV']:.6f} |")
    text += ["", "每行另含反向 k 的配对能级，因此能量出现两次。", "",
             "| gamma / meV | ground energy / eV | correlation energy / meV | largest eigenvalue of C |",
             "| --- | --- | --- | --- |"]
    for row in scan:
        text.append(f"| {1000*row['strength_eV']:.1f} | {row['ground_energy_eV']:.10f} | {1000*row['correlation_energy_eV']:.6f} | {row['pair_matrix_largest_eigenvalue']:.6f} |")
    text += ["", "关联能比较的是同一个相互作用哈密顿量下的精确基态与未关联填充态；它不是材料凝聚能。",
             "最大配对关联本征值的增加只描述本有限模型，不构成体相超导证据。", "",
             f"参考输出取 gamma={1000*reference_strength:g} meV，每个 G 元素={1000*reference_strength/m:.6f} meV。",
             f"基态残差 {result['residual_norm_eV']:.3e} eV；Pauli 映射误差 {pauli_error:.3e} eV。", "",
             "能量是减去参考费米能后的有限空间能量，不是 VASP 总能量。",
             "详见上一级 MODEL.md：模型约定、局限与读取方法。", ""]
    (OUT / "RESULTS.md").write_text("\n".join(text), encoding="utf-8")
    print("Pair levels:", m, "pairs:", n_pairs, "dimension:", len(basis))
    print("Reference gamma/eV:", reference_strength)
    print("Ground energy/eV:", result["ground_energy_eV"])
    print("Correlation energy/eV:", result["correlation_energy_eV"])
    print("Residual/eV:", result["residual_norm_eV"])
    print("Qubit mapping error/eV:", pauli_error)
    print("All outputs saved in results/.")


if __name__ == "__main__":
    main()
