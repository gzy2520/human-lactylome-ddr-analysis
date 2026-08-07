### Task 9: 并入目标蛋白审计 → `analyze_ddr.R`（第四部分）

**Files:**
- Modify: `reanalysis/scripts/analyze_ddr.R`
- 参照（只读）：`reanalysis/scripts/audit_target_sources.py`

**Interfaces:**
- Produces: 无新输出文件（Task 6 的 pipeline 段已翻译 `build_target_source_audit`）；本任务把 Python 原版与 R 版输出做交叉验证。

- [ ] **Step 1: 交叉验证 R 版审计输出**

Run:
```bash
cd /Users/gzy2520/Desktop/Research/kla
Rscript reanalysis/scripts/analyze_ddr.R . --stage pipeline
diff <(Rscript -e 'cat(readLines("reanalysis/results/tables/target_protein_source_level_audit_MRE11_XLF_NBS1.csv"), sep="\n")') \
     <(PYTHONPATH=reanalysis/scripts python3 -c '
import pandas as pd
print(pd.read_csv("reanalysis/results/tables/target_protein_source_level_audit_MRE11_XLF_NBS1.csv", dtype=str).to_csv(index=False))')
```
Expected: 无输出差异（同一文件在 R 版 pipeline 重跑后与 Python 版基线字节一致，由 verify_outputs 在 Task 6 已保证；本任务重点确认 audit 相关测试点 6/17 通过）。

- [ ] **Step 2: 检查点**

---

