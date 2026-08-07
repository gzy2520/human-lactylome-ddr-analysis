# 新对话项目启动 Prompt

你正在继续处理项目：

`/Users/gzy2520/Desktop/Research/kla`

这是一个人源蛋白乳酸化质谱重分析项目。开始任何工作前，先阅读：

1. `/Users/gzy2520/Desktop/Research/kla/PROJECT_INDEX.md`
2. `/Users/gzy2520/Desktop/Research/kla/reanalysis/README.md`
3. 与当前任务相关的 `reanalysis/reports/*.md`

## 必须遵守的分析原则

- 生物信息学任务优先使用 R；只有在 Python/Node 明显更适合时才使用其他语言。
- 不得使用 Gene Symbol 作为分析主键、蛋白去重键、GO 交集键或回退匹配键。
- 优先使用稳定 ID：UniProt BaseAccession、Ensembl protein ID 或 Entrez ID。
- 当前项目的蛋白主键是去除 isoform 后缀的 UniProt BaseAccession。
- Ensembl protein ID 必须通过明确映射转换后再与 UniProt GO 表匹配。
- `GeneSymbol` 只能作为来源元数据或显示标签，不能参与分析判断。
- 如果来源文件只有 Gene Symbol，必须先转换为稳定 ID，并输出转换审计、歧义记录和未映射记录。
- 不能用 Kla 信号代替普通全蛋白组强度。
- 不要删除任何数据；旧结果移动到 archive，不覆盖旧结果。
- 新结果写入现有 `reanalysis` 目录或新的子目录，并记录运行命令、日志和文件清单。
- 热图使用白色到暖色；低值接近白色，高值为浅黄、橙色到深红色；不要使用冷色作为低值主色。

## 当前数据范围

主 Kla 数据：

`PXD014870, PXD028488, PXD050470, PXD053474, PXD060185, PXD078013, PXD078736`

暂不分析：

`PXD038880, PXD050906`

PXD038880/PXD050906 只能保留原始下载文件、元数据和 excluded 说明，不能进入
Kla 蛋白集合、GO-DDR 交集或 Venn。

## 当前结果状态

- 扩展比较表有 40 个细胞系/组织样本组。
- 37 个已完成 Kla 与普通全蛋白参照的 DDR 占比计算。
- PXD037371 的 3 组因 TMT 通道映射问题排除。
- PXD050470 的 Kla 侧为 853 个蛋白，DDR 为 29 个。
- PXD050470 的普通参照是 PXD043880 正常 CA1 海马全蛋白组：
  - 74 个供体
  - 原始 2,092 个蛋白特征
  - symbol 转换后 2,105 个 reviewed UniProt BaseAccession
  - DDR 为 83 个
- PXD043880 的 symbol 转换审计：
  `/Users/gzy2520/Desktop/Research/kla/reanalysis/results/tables/PXD043880_hippocampus_symbol_to_reviewed_uniprot_mapping_audit.csv`
- 人源 reviewed UniProt 映射快照：
  `/Users/gzy2520/Desktop/Research/kla/reanalysis/config/uniprot_human_reviewed_2026-08-05.tsv`

## 两张热图的含义

1. `kla_regulator_cross_study_relative_intensity_heatmap.png`
   - 调控蛋白在各自 Kla 数据中的相对信号百分位。
   - 37/40 个样本组有可审计 Kla 定量。
   - 不是普通表达量，也不是 Log FC。

2. `kla_regulator_whole_proteome_relative_intensity_heatmap.png`
   - 调控蛋白在对应普通全蛋白组中的相对信号百分位。
   - 16/40 个样本组有普通全蛋白强度。
   - 没有强度的组显示 `?`，不能用 Kla 信号替代。

## 关键输出

- DDR 占比图：
  `/Users/gzy2520/Desktop/Research/kla/reanalysis/results/figures/cell_type_kla_vs_reference_ddr_fraction.png`
- 中文 DDR 统计表：
  `/Users/gzy2520/Desktop/Research/kla/reanalysis/results/tables/cell_type_kla_vs_reference_ddr_statistics_accession_only_zh.csv`
- Kla 调控因子热图：
  `/Users/gzy2520/Desktop/Research/kla/reanalysis/results/figures/kla_regulator_cross_study_relative_intensity_heatmap.png`
- 普通全蛋白调控因子热图：
  `/Users/gzy2520/Desktop/Research/kla/reanalysis/results/figures/kla_regulator_whole_proteome_relative_intensity_heatmap.png`
- 普通全蛋白可用性审计：
  `/Users/gzy2520/Desktop/Research/kla/reanalysis/results/tables/kla_regulator_whole_proteome_intensity_availability_audit.csv`
- 最终文件清单：
  `/Users/gzy2520/Desktop/Research/kla/reanalysis/reports/final_file_manifest_sha256.csv`

## 推荐先做的事情

先运行：

```bash
cd /Users/gzy2520/Desktop/Research/kla
sed -n '1,260p' PROJECT_INDEX.md
```

然后根据用户的新要求选择相关脚本。修改代码前先检查当前文件和结果，不要假设
旧报告中的数字仍然是最新数字。完成后运行相关 R 测试，并更新：

```bash
python3 reanalysis/scripts/build_final_manifest.py
```

回复用户时，要明确区分：

- Kla 蛋白检测/定量；
- 普通全蛋白参照；
- GO-DDR 交集；
- Kla 信号热图；
- 普通全蛋白信号热图；
- 哪些样本缺少普通强度而只能显示 `?`。
