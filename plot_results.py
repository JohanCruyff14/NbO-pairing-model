"""Plot physical observables from the existing finite-pairing calculation.

Run build_pairing_model.py first, then run this script from any directory.
No eigensolver, DFT calculation, or parameter changes are performed here.
"""
from pathlib import Path
import csv
import hashlib
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BASE = Path(__file__).resolve().parent
RESULTS = BASE / "results"
FIGURES = BASE / "figures"


def read_csv(name):
    with (RESULTS / name).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(name, rows):
    with (FIGURES / name).open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def save_figure(fig, name):
    # Vector files preserve editable text; PNG is the convenient screen preview.
    for extension in ("png", "svg", "pdf"):
        fig.savefig(FIGURES / f"{name}.{extension}", dpi=300,
                    facecolor="white", bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)


def plot_occupations(xi, baseline, occupation, strength):
    order = np.argsort(xi, kind="stable")
    x = np.arange(len(xi))
    fig, (ax, delta_ax) = plt.subplots(
        2, 1, sharex=True, figsize=(183 / 25.4, 120 / 25.4),
        gridspec_kw={"height_ratios": [2.3, 1]}, layout="constrained")
    fig.suptitle("Pair occupation redistribution", fontsize=11)
    ax.set_title(rf"$\gamma={strength * 1000:g}$ meV; fixed pair number", fontsize=8)
    ax.scatter(x, baseline[order], facecolors="none", edgecolors="0.4",
               marker="D", s=38, linewidths=0.9, label=r"$\gamma=0$")
    ax.scatter(x, occupation[order], color="#286B9B", marker="o", s=16,
               label=rf"$\gamma={strength * 1000:g}$ meV", zorder=3)
    ax.set_ylabel(r"Pair occupation $\langle N_p\rangle$")
    ax.set_ylim(-0.07, 1.14)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
    ax.legend(loc="center right", fontsize=7)
    delta = occupation - baseline
    delta_ax.bar(x, delta[order], width=0.62, color="#286B9B")
    delta_ax.axhline(0, color="0.4", linewidth=0.7)
    delta_ax.set_ylabel(r"$\Delta\langle N_p\rangle$")
    delta_ax.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
    delta_ax.set_xticks(x, [f"p{p}\n{xi[p] * 1000:.2f}" for p in order], fontsize=6.8)
    delta_ax.set_xlabel(r"Pair level $p$ / $\xi_p$ (meV), ordered by energy; equal level spacing")
    boundary = int(baseline.sum()) - 0.5
    for panel, label in ((ax, "a"), (delta_ax, "b")):
        panel.axvline(boundary, color="0.65", linestyle="--", linewidth=0.8)
        panel.set_xlim(-0.65, len(xi) - 0.35)
        panel.text(-0.075, 1.02, label, transform=panel.transAxes,
                   fontweight="bold", fontsize=9)
        panel.grid(axis="y", color="0.92", linewidth=0.5)
        panel.set_axisbelow(True)
    ax.text(boundary, 1.055, r"$\xi=0$ boundary", ha="center", fontsize=7)
    save_figure(fig, "occupations")


def plot_correlations(correlation, bands, strength):
    m = len(correlation)
    offdiag = np.ma.array(correlation, mask=np.eye(m, dtype=bool))
    maximum = float(np.max(np.abs(offdiag)))
    fig, axes = plt.subplots(1, 2, figsize=(183 / 25.4, 94 / 25.4),
                             layout="constrained")
    fig.suptitle(rf"Pair correlations: $\gamma={strength * 1000:g}$ meV", fontsize=11)
    panels = [(correlation, "Full matrix", 1.0,
               r"$C_{pq}$ (dimensionless)"),
              (offdiag, "Transfer coherence (diagonal masked)", maximum or 1.0,
               r"$C_{pq}$, $p\ne q$ (dimensionless)")]
    for ax, (values, title, vmax, color_label), letter in zip(axes, panels, ("a", "b")):
        # Do not discard a sign if a future real-valued model has negative correlations.
        negative = float(np.min(values)) < -1e-12
        cmap = plt.get_cmap("RdBu_r" if negative else "Blues").copy()
        cmap.set_bad("#E5E5E5")
        image = ax.imshow(values, origin="upper", interpolation="nearest", cmap=cmap,
                          vmin=-vmax if negative else 0, vmax=vmax)
        ax.set_title(title, fontsize=8, pad=8)
        ax.set_xticks(range(m), labels=range(m), fontsize=6)
        ax.set_yticks(range(m), labels=range(m), fontsize=6)
        ax.set_xlabel(r"Source pair level $q$")
        ax.set_ylabel(r"Destination pair level $p$")
        ax.text(-0.14, 1.08, letter, transform=ax.transAxes,
                fontweight="bold", fontsize=9)
        for boundary in np.flatnonzero(np.diff(bands)) + 0.5:
            ax.axhline(boundary, color="0.6", linewidth=0.5)
            ax.axvline(boundary, color="0.6", linewidth=0.5)
        colorbar = fig.colorbar(image, ax=ax, orientation="horizontal",
                               fraction=0.065, pad=0.13, shrink=0.95)
        colorbar.set_label(color_label, fontsize=7)
        colorbar.ax.tick_params(labelsize=6)
    save_figure(fig, "pair_correlations")


def main():
    inputs = ("active_levels.csv", "occupations.csv", "pair_correlation.csv", "model_record.json")
    missing = [name for name in inputs if not (RESULTS / name).is_file()]
    if missing:
        raise SystemExit("Missing results: " + ", ".join(missing)
                         + ". Run python build_pairing_model.py first.")
    levels = sorted(read_csv("active_levels.csv"), key=lambda row: int(row["p"]))
    occupations = sorted(read_csv("occupations.csv"), key=lambda row: int(row["p"]))
    correlations = sorted(read_csv("pair_correlation.csv"), key=lambda row: int(row["p"]))
    record = json.loads((RESULTS / "model_record.json").read_text(encoding="utf-8"))
    m, n_pairs = record["n_pair_levels"], record["n_pairs"]
    for rows in (levels, occupations, correlations):
        if [int(row["p"]) for row in rows] != list(range(m)):
            raise ValueError("Result level indices are missing, repeated, or inconsistent.")
    xi = np.array([float(row["xi_eV"]) for row in levels])
    occupation = np.array([float(row["pair_occupation"]) for row in occupations])
    correlation = np.array([[float(row[f"q{q}"]) for q in range(m)] for row in correlations])
    bands = np.array([int(row["band_index"]) for row in levels])
    baseline = (xi < 0).astype(float)
    strength = float(record["reference_result"]["strength_eV"])
    if np.any(xi == 0) or int(baseline.sum()) != n_pairs:
        raise ValueError("Zero-coupling filling is ambiguous or differs from the recorded pair number.")
    checks = {
        "all_values_finite": bool(np.isfinite(correlation).all() and np.isfinite(xi).all()
                                  and np.isfinite(occupation).all()),
        "occupation_matches_matrix_diagonal": bool(np.allclose(occupation, np.diag(correlation), atol=1e-10, rtol=0)),
        "occupation_matches_level_energies": bool(np.allclose(
            xi, [float(row["xi_eV"]) for row in occupations], atol=1e-12, rtol=0)),
        "pair_number_conserved": bool(abs(occupation.sum() - n_pairs) < 1e-9),
        "occupation_in_unit_interval": bool(np.all(occupation >= -1e-10) and np.all(occupation <= 1 + 1e-10)),
        "correlation_symmetric": bool(np.allclose(correlation, correlation.T, atol=1e-10, rtol=0)),
        "correlation_positive_semidefinite": bool(np.linalg.eigvalsh(correlation)[0] >= -1e-9),
    }
    if not all(checks.values()):
        raise ValueError(f"Inconsistent physical observables: {checks}")
    FIGURES.mkdir(exist_ok=True)
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": 0.7, "legend.frameon": False,
        "svg.fonttype": "none", "pdf.fonttype": 42,
    })
    plot_occupations(xi, baseline, occupation, strength)
    plot_correlations(correlation, bands, strength)
    write_csv("occupation_comparison.csv", [dict(
        p=p, band_index=int(bands[p]), orientation=int(levels[p]["orientation"]),
        xi_meV=1000 * xi[p], pair_occupation_gamma_zero=baseline[p],
        pair_occupation_reference=occupation[p], delta_pair_occupation=occupation[p] - baseline[p])
        for p in range(m)])
    write_csv("correlation_comparison.csv", [dict(
        p=p, q=q, C_gamma_zero=baseline[p] if p == q else 0.0,
        C_reference=correlation[p, q], is_diagonal=p == q)
        for p in range(m) for q in range(m)])
    offdiag = correlation.copy()
    np.fill_diagonal(offdiag, 0)
    summary = dict(
        reference_strength_meV=1000 * strength, G_each_element_meV=1000 * strength / m,
        n_pairs=n_pairs, pair_occupation_sum=float(occupation.sum()),
        pairs_transferred_above_reference_EF=float(occupation[xi > 0].sum()),
        max_absolute_occupation_change=float(np.max(abs(occupation - baseline))),
        max_absolute_offdiagonal_correlation=float(np.max(abs(offdiag))),
        offdiagonal_correlation_sum=float(offdiag.sum()),
        checks=checks, numpy_version=np.__version__, matplotlib_version=matplotlib.__version__,
        inputs_sha256={name: hashlib.sha256((RESULTS / name).read_bytes()).hexdigest() for name in inputs},
        baseline="Exact gamma=0 product state of the same finite model, not VASP smearing occupations.",
        statistics="Deterministic ground-state observables; no sampling or statistical error bars.",
        offdiagonal_display="Diagonal masked, not set to zero; separate linear color scale.",
    )
    (FIGURES / "plot_record.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    notes = f"""# 有限配对模型：两张图的读取方法

两张图由 `plot_results.py` 读取现有计算输出生成。PNG 供查看，SVG/PDF 供放大或编辑；CSV 是实际绘图数据。图中使用英文标签以便跨平台复现。

默认参考参数：gamma = {1000 * strength:g} meV；每个 G_pq = {1000 * strength / m:.6f} meV；固定 {n_pairs} 对电子。

## 1. 占据数对比

![占据数对比](occupations.png)

- a：空心菱形为 gamma=0 的填充态，蓝点为参考作用强度下的基态。纵轴是电子对占据数 `<N_p>`，取值在 0 到 1；平均电子数是其两倍。
- 横轴逐个展示按 xi 排序的配对能级；每个标签上行为原始编号 p，下行为 xi（meV）。横向间距只是等距排列，不是实际能量间隔。虚线分开负能级和正能级。
- 两个配对方向即使能量相同也保留为两个点，未平均或合并。编号、能带及方向对应 `occupation_comparison.csv` 和计算输出的 `active_levels.csv`。
- b：`Delta<N_p> = <N_p>_reference - <N_p>_gamma=0`。负能级失去部分占据，正能级获得部分占据；总变化为零。
- 正能级的总电子对占据为 **{summary['pairs_transferred_above_reference_EF']:.6f}**；总占据仍为 **{occupation.sum():.6f}**。分数占据来自多体叠加，不代表恢复了单电子占据态。
- 零作用基线是同一有限模型的精确填充态，不是 VASP 展宽后的单粒子占据数。

## 2. 配对关联矩阵

![配对关联矩阵](pair_correlations.png)

- 定义：`C_pq = <Psi_0 | b_p^dagger b_q | Psi_0>`，无量纲。行 p 是电子对到达的能级，列 q 是电子对移出的能级。
- a：完整矩阵，`C_pp = <N_p>`；色标范围固定为 0 到 1。
- b：同一矩阵的非对角部分，灰色对角线被遮罩，并非数值为零。采用独立的线性色标，显式展示较小的转移相干性；不能直接跨面板比较颜色深浅。
- p 按原始编号排列，细线划分来自不同候选能带的能级块。
- 零作用填充态的所有非对角元为零；参考作用下最大非对角元绝对值为 **{summary['max_absolute_offdiagonal_correlation']:.6f}**。逐项对比见 `correlation_comparison.csv`。
- 非对角元是相干期望值，不是转移概率或转移速率，也不是具有空间距离坐标的长程关联。固定粒子数下 `<b_p>` 为零，不妨碍 `<b_p^dagger b_q>` 非零。

## 复现与检查

在项目目录执行：

```bash
python build_pairing_model.py
python plot_results.py
```

这是确定性的有限模型基态计算，没有随机抽样，不画统计误差条。数值检查包括占据范围、总电子对数、关联矩阵对称性与半正定性，以及占据数与矩阵对角元的一致性；检查和输入哈希记录在 `plot_record.json`。

这些图说明参数化吸引作用下，当前受限配对空间中的占据重分配和能级间关联；不能单独证明 NbO 体相超导。DFT 提供 xi，配对作用仍为人为设置。
"""
    (FIGURES / "README.md").write_text(notes, encoding="utf-8")
    print(f"Saved two figures (PNG/SVG/PDF), source CSVs and interpretation to {FIGURES}")
    print(f"Pairs above reference EF: {summary['pairs_transferred_above_reference_EF']:.6f}")
    print(f"Largest |C_pq| for p != q: {summary['max_absolute_offdiagonal_correlation']:.6f}")
    print("Physical observable checks:", checks)


if __name__ == "__main__":
    main()
