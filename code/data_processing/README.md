# 📁 data_processing 文件夹说明 / Folder Overview

本文件夹包含对 SNP 数据进行预处理、GWAS 分析、以及后续特征筛选的核心代码。  
This folder contains essential scripts for preprocessing SNP data, conducting GWAS, and filtering features for downstream modeling.

---

## 🔧 工具路径 / Tool Paths

| 工具名称 / Tool Name       | 使用路径 / Path                                                           |
|----------------------------|---------------------------------------------------------------------------|
| `plink`                    | `/nfs/my/Huang/lzm/default/plink/plink`                                  |
| `PopLDdecay`               | `/nfs/my/Huang/lzm/default/PopLDdecay-3.42/bin/PopLDdecay`               |
| `Plot_OnePop.plink`        | `/nfs/my/Huang/lzm/default/PopLDdecay-3.42/bin/Plot_OnePop.pl`           |
| `admixture`                | `/nfs/my/Huang/lzm/admixture_linux-1.3.0/admixture`                      |
| `S2_gemma.sh`              | `/nfs/my/Huang/lzm/soybean_C219/code/data_processing/S2_gemma.sh`        |
| `S3_assoc2CMplot.R`        | `/nfs/my/Huang/lzm/soybean_C219/code/data_processing/S3_assoc2CMplot.R`  |
| `S4_CMplot.R`              | `/nfs/my/Huang/lzm/soybean_C219/code/data_processing/S4_CMplot.R`        |
| `S5_p_wald_filter.py`      | `/nfs/my/Huang/lzm/soybean_C219/code/data_processing/S5_p_wald_filter.py`|
| `S6_p_wald_shap.py`        | `/nfs/my/Huang/lzm/soybean_C219/code/data_processing/S6_p_wald_shap.py`  |

---

## 📜 分析流程 / Analysis Workflow

### 中文说明：

1. `S1_hmp2ped.py`：将 `219_snp.hmp` 转换为 `.ped` 格式。
2. 根据表型样本划分数据集为五折。
3. 对每一折训练集进行质量控制（QC），参数为 `--geno 0.05 --maf 0.05`。
4. 进行连锁不平衡（LD）分析：`--indep-pairwise 50 5 0.2`，筛选独立 SNP。
5. `S2_gemma.sh`：使用 GEMMA 对每个性状做 GWAS 分析（LD前后都做），可用 `run_gemma.sh` 批量运行。
6. `S3_assoc2CMplot.py`：整合连续三年 GWAS p 值，生成 CMplot 格式文件，建议用 `run_assoc2CMplot.sh`。
7. `S4_CMplot.R`：使用 R 包 `CMplot` 画出曼哈顿图、Q-Q 图、密度图等。
8. `S5_p_wald_filter.py`：筛选 `p_wald < 0.01` 的显著 SNP 位点，取 2015 年与三年的交集 rs。
9. 用 plink 根据 rs 提取 SNP 数据，保存为 `.ped` 和 `.raw` 格式用于后续 ML / DL 模型。
10. `S6_p_wald_shap.py`：从最优折中选出用于 SHAP 分析的 SNP 位点，并附加 SHAP 值。

---

### English Description:

1. `S1_hmp2ped.py`: Converts `219_snp.hmp` to `.ped` format.
2. Split dataset into 5 folds based on samples with phenotypes.
3. Apply quality control (QC) to each training fold with `--geno 0.05 --maf 0.05`.
4. Perform linkage disequilibrium (LD) filtering: `--indep-pairwise 50 5 0.2`.
5. `S2_gemma.sh`: Run GWAS for each trait using GEMMA before and after LD filtering. Use `run_gemma.sh` for automation.
6. `S3_assoc2CMplot.py`: Merge GWAS p-values from 3 consecutive years into CMplot format. Run with `run_assoc2CMplot.sh`.
7. `S4_CMplot.R`: Generate Manhattan, Q-Q, and density plots using the `CMplot` R package.
8. `S5_p_wald_filter.py`: Select SNPs with `p_wald < 0.01`, and find the intersection of 2015 and 3-year sets.
9. Use `plink` to extract SNPs by rs IDs, output as `.ped` and `.raw` for ML/DL modeling.
10. `S6_p_wald_shap.py`: Select SNPs from the best fold for SHAP analysis and compute SHAP values.

---

📌 **Note**: Please make sure to update file and tool paths to your own environment before running the scripts.
