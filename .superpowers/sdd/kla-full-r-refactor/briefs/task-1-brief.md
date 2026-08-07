### Task 1: 建立回归基线（输出 sha256 快照 + 现有测试确认）

**Files:**
- Create: `archive/refactor_baseline_2026-08-06/output_sha256_baseline.csv`
- Create: `archive/refactor_baseline_2026-08-06/PRE_REFACTOR_TEST_RESULTS.txt`

**Interfaces:**
- Produces: 基线表 `output_sha256_baseline.csv`（列：`RelativePath,SizeBytes,SHA256`），后续所有任务的验证都 diff 它。

- [ ] **Step 1: 确认现有测试全过**

Run:
```bash
cd /Users/gzy2520/Desktop/Research/kla
Rscript reanalysis/tests/test_expanded_ddr_fraction_by_accession.R .
Rscript reanalysis/tests/test_kla_regulator_intensity.R .
Rscript reanalysis/tests/test_kla_regulator_landscape.R .
Rscript reanalysis/tests/test_kla_regulator_whole_proteome_intensity.R .
Rscript reanalysis/tests/test_lactylome_acquisition.R .
PYTHONPATH=reanalysis/scripts python3 -m unittest discover -s reanalysis/tests -p 'test_*.py' -v 2>&1 | tail -5
```
Expected: 6 个测试文件全部通过（含 Python 侧 19 项），输出保存到 `PRE_REFACTOR_TEST_RESULTS.txt`。若任一失败，先停下报告，不要带病重构。

- [ ] **Step 2: 生成输出文件 sha256 基线**

```bash
cd /Users/gzy2520/Desktop/Research/kla
mkdir -p archive/refactor_baseline_2026-08-06
# 用现有 manifest 覆盖的输出全集做基线（reanalysis/ + previous_umap/ + 根文档，排除自身）
python3 - <<'EOF'
import hashlib, csv
from pathlib import Path
root = Path(".")
paths = set()
for f in (Path("PROJECT_INDEX.md"), Path("NEW_CHAT_PROJECT_PROMPT.md")):
    if f.exists(): paths.add(f)
for d in (Path("reanalysis"), Path("previous_umap")):
    paths.update(p for p in d.rglob("*") if p.is_file())
rows = []
for p in sorted(paths, key=lambda x: x.as_posix()):
    if p.name == ".DS_Store" or p.suffix == ".pyc": continue
    rows.append({"RelativePath": p.as_posix(), "SizeBytes": p.stat().st_size,
                 "SHA256": hashlib.sha256(p.read_bytes()).hexdigest()})
with open("archive/refactor_baseline_2026-08-06/output_sha256_baseline.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["RelativePath","SizeBytes","SHA256"]); w.writeheader(); w.writerows(rows)
print(f"Baseline files hashed: {len(rows)}")
EOF
```
Expected: 输出 "Baseline files hashed: N"（N 与当前 manifest 行数接近）。

- [ ] **Step 3: 写基线 diff 工具脚本（后续任务复用）**

Create: `reanalysis/scripts/lib/verify_outputs.R`（保留在 lib/ 顶层，供各任务手工调用；本脚本不属于合并目标，是重构期间的工具）：
```r
#!/usr/bin/env Rscript
# 用法: Rscript lib/verify_outputs.R <project_root> <baseline.csv> [排除正则...]
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(args[[1]])
baseline <- read.csv(args[[2]], check.names = FALSE, stringsAsFactors = FALSE)
exclude <- if (length(args) >= 3) args[[3]] else NULL
sha1 <- function(p) {
  con <- file(p, "rb"); on.exit(close(con))
  openssl::sha256(con)
}
changed <- list()
for (i in seq_len(nrow(baseline))) {
  rel <- baseline$RelativePath[i]
  p <- file.path(root, rel)
  if (!is.null(exclude) && grepl(exclude, rel, perl = TRUE)) next
  if (!file.exists(p)) { changed[[length(changed) + 1]] <- data.frame(RelativePath = rel, Status = "MISSING"); next }
  if (baseline$SizeBytes[i] != file.info(p)$size ||
      baseline$SHA256[i] != sha1(p)) {
    changed[[length(changed) + 1]] <- data.frame(RelativePath = rel, Status = "CHANGED")
  }
}
if (length(changed)) {
  out <- do.call(rbind, changed)
  write.csv(out, file.path(dirname(args[[2]]), "diff_report.csv"), row.names = FALSE)
  cat("DIFF FOUND:", nrow(out), "files\n")
  print(out)
  quit(status = 1)
}
cat("OK: all", nrow(baseline), "baseline files unchanged\n")
```
Expected: 运行 `Rscript reanalysis/scripts/lib/verify_outputs.R . archive/refactor_baseline_2026-08-06/output_sha256_baseline.csv` 输出 "OK: all N baseline files unchanged"。

- [ ] **Step 4: 检查点（代替 commit）**

将本任务产出写入检查点文件 `docs/superpowers/plans/CHECKPOINTS.md`（每任务完成后追加一行"Task N 完成: <验证命令> <结果摘要>"）。

---

