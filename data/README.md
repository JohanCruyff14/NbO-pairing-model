# DFT 输入数据

材料结构参考：[Materials Project mp-2311](https://legacy.materialsproject.org/materials/mp-2311/)。本目录的能量与投影来自用户独立执行的 VASP PBE 计算，并非直接下载的网站能带。

六原子立方晶胞（Nb3 O3），a=4.261457 Å；非磁性 ISPIN=1、无 SOC；ENCUT=520 eV；Gamma-centered 16×16×16 网格，165 个不可约点，32 带，51 个电子；高斯展宽 0.05 eV。EF=5.97117024 eV。这些参数未经过系统收敛测试。

| 文件 | 内容 |
| --- | --- |
| candidate_bands_25_27.csv | 三带 × 165 个不可约 k 点，共 495 行，模型的必要输入 |
| band_summary.csv | 全 32 带的范围、电子数、费米能穿越情况，供选带依据检查 |
| metadata.json | 提取约定、核对结果和原始文件 SHA256 |

CSV 为 UTF-8 BOM，索引从 1 开始，能量单位 eV。

- `k_index,band_index`：原始不可约 k 点和 VASP 能带序号。
- `k1,k2,k3`：倒格矢的分数坐标。
- `weight`：布里渊区积分权重，打印精度下总和约为 1；不是有限配对能级的振幅。
- `energy_eV`：原始 KS 能量；`xi_eV=energy_eV-EF`。
- `occupation_per_spin`：当前 ISPIN=1 XML 的 0–1 占据数；积分电子数需要乘 2。高斯展宽占据不等于零温阶跃。
- `paw_Nb_s` 等：PROCAR 的物种/角动量原始投影；除以 `paw_sum` 可得到投影内部比例，但不是完整波函数概率。

原始 VASP 输出未随仓库提供。用自有输出重新生成上述表格：

```bash
python extract_scf.py --source /path/to/02_scf
```

所需文件：`vasprun.xml`、`EIGENVAL`、`PROCAR`、`POSCAR`；来源哈希还会记录该目录的其它文件。提取器针对本 NbO 非磁性结构编写，不是通用 VASP 解析器。默认会更新 `data/`，若只想核对，请添加 `--output /path/to/check_export`。

哈希用于标记原始快照，原始文件不在仓库内，不能仅凭这些表格重新完成 DFT。表格不含配对作用，也不足以计算电子声子矩阵元。
