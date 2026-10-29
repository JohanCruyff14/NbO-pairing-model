# 第一版：使用 NbO 能量的有限配对哈密顿量

## 定义

使用 data 中第 25–27 带的真实 xi = epsilon - EF。每带选距离 EF 最近的一个负能态、一个正能态。没有对不可约网格做等权积分，也没有将权重乘入能量；这是选定的离散轨道原型，尚未对舍弃态做受控下折叠。

对六个源态分别保留两个配对方向（此处没有时间反演不变点）：

- p=(n,k)：b†_p = d†_(n,k,up) d†_(n,-k,down)。
- 反向能级：b†_(n,-k) = d†_(n,-k,up) d†_(n,k,down)。

两种方向用到不同自旋轨道，不能合并成同一个硬核配对能级；它们的对称组合可构成自旋单态。反向能量依据非磁性、无 SOC 的时间反演对称性补齐，而非声称 CSV 已包含所有这些点。保留时间反演相关结构，但未补齐整个晶体点群的 k 点星。

共有 M=12 个配对能级（24 个自旋轨道）。每级只保留空态和双占据态，即 seniority-zero 限制，排除单占据和破对激发。这是特定配对空间的基态，不保证是完整电子问题的基态。

```text
H = 2 sum_p xi_p N_p - sum_pq G_pq b†_p b_q
G_pq = gamma/M  (所有 p,q，包括 p=q)
```

gamma 是人为选择的有限模型作用强度，未从电子声子作用导出。除以 M 是本原型的归一化约定，不等于从原 16^3 网格直接得到的相互作用。不能将其数值解释成 NbO 材料参数。

默认 gamma=0.02 eV。另算 0、0.01、0.05、0.1 eV 作对照，不研究温度或相变。固定粒子数下，对角作用为常数 -gamma*N_pairs/M，报告关联能时用同一个 H 比较，避免将这一平移误认为关联效果。

六个负能配对能级在零作用时填满，所以固定 N_pairs=6（12 个活性电子），D=choose(12,6)=924。这个粒子数来自选定轨道的零作用填充规则，不是将体相每晶胞的三电子直接映射到小模型，更不是已经完成冻结态计数。

## 量子比特形式

令 |0>=空、|1>=双占据，N_p=(1-Z_p)/2，b†=(X-iY)/2：

```text
Hq = C + sum_p (-xi_p + G_pp/2) Z_p
       - (1/2) sum_(p<q) G_pq (X_p X_q + Y_p Y_q)
C = sum_p (xi_p - G_pp/2)
```

配对编码需要 12 个量子比特，完整配对 Hilbert 空间维度 4096；当前只解其中 pair_number=6 的 924 维子空间。JSON 中 qubit p 对应整数位 p（最低位是 p=0）。初态可以是六个 xi<0 的位为 1 的填充态。量子变分线路需要保持总激发数，否则可能跑到不同粒子数扇区。

## 运行

使用 Python 3.12，安装 requirements.txt 中的 NumPy、SciPy。进入本目录：

```bash
python -m pip install -r requirements.txt
python build_pairing_model.py
```

修改 model_config.json 中 reference_strength_eV 和 scan_strengths_eV 后重跑；会覆盖 results 中自动生成的本地结果，不修改源数据。现阶段选择规则和粒子数规则固定，其它选择需要明确修改脚本。输出不提交到 Git，默认参数的轻量快照保存在 reference_results 中。

## 输出与复用

| 文件 | 内容 |
| --- | --- |
| results/active_levels.csv | 12 个配对能级、自旋上下方向的 k 与 -k、DFT 能量及源索引 |
| results/G_matrix_eV.csv | 12×12 作用矩阵，单位 eV |
| results/hamiltonian_fixed_pairs.npz | 924×924 稀疏矩阵，按 basis.csv 的顺序排列 |
| results/basis.csv | 固定六对电子的基底，bitmask 明确标出占据 |
| results/qubit_hamiltonian.json | 包括常数项的 Pauli 项，可用于量子算法 |
| results/ground_state.csv | 基态振幅与概率，与 basis.csv 对应 |
| results/occupations.csv | 每级电子对占据与电子占据 |
| results/pair_correlation.csv | C_pq = <b†_p b_q> |
| results/strength_scan.csv | 作用强度、基态能、关联能和配对关联 |
| results/model_record.json | 假设、参数、输入哈希与数值验证记录 |
| results/RESULTS.md | 可直接阅读的结果摘要 |

加载矩阵：

```python
from scipy.sparse import load_npz
from scipy.sparse.linalg import eigsh
H = load_npz('results/hamiltonian_fixed_pairs.npz')
energies, states = eigsh(H, k=1, which='SA')
print(energies[0])
```

脚本检查：自旋轨道无重复、零作用极限、厄米性、基态残差、配对关联矩阵的迹与半正定性，并通过独立 Pauli 作用核对量子哈密顿量的矩阵及粒子数守恒。

下一阶段应改变活性空间和配对作用形式，检查结论依赖性。不能用这一有限原型单独证明 NbO 体相超导，也不直接预测材料凝聚能。
