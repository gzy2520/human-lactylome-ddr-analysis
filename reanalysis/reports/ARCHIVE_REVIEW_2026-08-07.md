# 项目归档审计

审计日期：2026-08-07

## 已归档

- `reanalysis/scripts/build_final_manifest.py`：R 版与 Python 版均生成 552 行清单，
  两份 CSV 逐字节一致，SHA256 均为
  `842a6c82255bc7022b631d999e7207a172dafead5844b6ab78f97406f19c020a`。
  因项目生物信息学流程优先使用 R，正式入口保留为
  `reanalysis/scripts/build_final_manifest.R`。
- `reanalysis/reports/.Rhistory`：交互式历史不是可重复分析输入，移入归档。
- `reanalysis/scripts/analyze_regulators.R`：旧整合入口没有同步严格普通参照策略，
  已于 2026-08-08 移入
  `archive/2026-08-08_strict_reference_policy/scripts/`。
- `reanalysis/reports/REANALYSIS_37GROUP_REFERENCE_UPDATE_2026-08-07.md`：
  旧报告错误声称 37 组都具有普通全蛋白参照，已由
  `STRICT_REFERENCE_UPDATE_2026-08-08.md` 替代并移入本轮归档。
- `reanalysis/logs/workbook_inspection/*.inspect.ndjson`：旧工作簿检查缓存共约
  191 MB，当前工作簿已重建并重新检查，缓存移入本轮归档以保持主目录清晰。

迁移记录：

`archive/2026-08-07_git_baseline_cleanup/MIGRATION_MANIFEST.csv`

## 暂时保留

下列脚本与新的 R 主流程存在功能重叠，但仍被旧报告、回归测试或历史工作簿引用。
本轮不移动，待做对应输出的逐文件回归后再归档：

| 文件或类别 | 当前判断 |
|---|---|
| `run_pipeline.py`、`extractors.py`、`common.py` | 旧 Python 主流程；仍被 `test_pipeline.py` 使用 |
| `analyze_reference_proteome_ddr.py` | 旧参照蛋白组 DDR 流程；历史报告仍引用 |
| 旧 `build_*workbook.mjs` | 多数功能已并入 `build_workbooks.R`，但旧交付仍可追溯 |
| `build_project_metadata.py`、`audit_target_sources.py` | 旧项目盘点与来源审计；保留用于历史复现 |

## 保持活动

- 当前 Kla 景观、Kla 强度、普通全蛋白强度、DDR 占比和四分类 Venn 的拆分 R 脚本。
- 所有自动测试、配置、ID 映射表和方法报告。
- `data` 中原始数据、检索结果、补充表与 metadata。
- 已有 `archive` 内容。归档目录本身不再重复迁移。

## Git 策略

Git 跟踪代码、配置、测试、报告、GO 注释表和调控因子 identifier 表。
大型原始质谱数据、历史 archive、中间文件、日志、生成图片和生成结果表继续保留在本机，
但由 `.gitignore` 排除；最终结果完整性由
`reanalysis/reports/final_file_manifest_sha256.csv` 记录。
