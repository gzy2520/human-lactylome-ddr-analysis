### Task 11: 迁移 manifest → `build_final_manifest.R`

**Files:**
- Create: `reanalysis/scripts/build_final_manifest.R`
- 参照（只读）：`reanalysis/scripts/build_final_manifest.py`

**Interfaces:**
- Produces: CLI `Rscript build_final_manifest.R <project_root>`；输出 `reanalysis/reports/final_file_manifest_sha256.csv`，表头 `RelativePath,SizeBytes,SHA256`，包含范围与 Python 版完全一致（根两文档 + reanalysis/ + previous_umap/ + data/PXD*/metadata/，排除 .DS_Store/.pyc/自身）。

- [ ] **Step 1: 翻译为 R**

用 `openssl::sha256(file)` 或 `digest::digest(file, algo = "sha256", serialize = FALSE)`；目录遍历用 `list.files(recursive = TRUE)`；排序键 `path.as_posix()` 与 Python 一致。

- [ ] **Step 2: 对比新旧 manifest（同文件集）**

Run:
```bash
cd /Users/gzy2520/Desktop/Research/kla
Rscript reanalysis/scripts/build_final_manifest.R . && mv reanalysis/reports/final_file_manifest_sha256.csv /tmp/manifest_r.csv
python3 reanalysis/scripts/build_final_manifest.py && mv reanalysis/reports/final_file_manifest_sha256.csv /tmp/manifest_py.csv
diff /tmp/manifest_r.csv /tmp/manifest_py.csv && echo MANIFEST_IDENTICAL
```
Expected: `MANIFEST_IDENTICAL`（同一时刻同文件集，两个实现输出逐字节一致）。

- [ ] **Step 3: 检查点**

---

