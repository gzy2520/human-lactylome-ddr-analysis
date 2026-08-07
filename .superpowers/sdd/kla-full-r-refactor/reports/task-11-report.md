# Task 11 报告：manifest → `build_final_manifest.R`

**状态**: DONE

**产出文件**:
- 新建 `/Users/gzy2520/Desktop/Research/kla/reanalysis/scripts/build_final_manifest.R`（77 行）
- 更新 `/Users/gzy2520/Desktop/Research/kla/reanalysis/reports/final_file_manifest_sha256.csv`（499 行 = 表头 + 498 个文件）
- 未修改任何 `.py` 文件

## R 实现说明

接口：`Rscript build_final_manifest.R <project_root>`，与 Python 版（`Path(__file__).resolve().parents[2]` 推导根目录）等效——按 brief 约定取命令行参数，验证命令从项目根传入 `.`。

实现要点（逐一对应 Python 版 68 行逻辑）：

1. **包含范围**：根目录 `PROJECT_INDEX.md` + `NEW_CHAT_PROJECT_PROMPT.md`（存在才加入）+ `reanalysis/` 全树 + `previous_umap/` 全树 + `data/PXD*/metadata/` 全树（PXD 目录按排序访问）。与 Python 版完全一致。
2. **目录遍历**：`list.files(dir, recursive = TRUE, all.files = TRUE, full.names = TRUE)`。**`all.files = TRUE` 是关键**——Python 的 `pathlib.rglob("*")` 包含点文件（如 `reanalysis/reports/.Rhistory` 确实在旧 manifest 中），而 R 默认 `all.files = FALSE` 会漏掉它们；不设此项会导致逐字节不一致。
3. **sha256**：`digest::digest(file = path, algo = "sha256", serialize = FALSE)`，返回 64 位小写 hex，与 `hashlib.sha256().hexdigest()` 一致。digest 对 file 走 C 层 `digest_impl`（分块流式读取），实测 5 GB 文件无内存问题（全树约 0.7 GB，R 版全程约 4 秒）。
4. **排除规则**：
   - `basename != ".DS_Store"`；
   - `.pyc` 排除用**模拟 pathlib 的 suffix 语义**（`has_pyc_suffix`）：只有"最后一个点不在首字符位"且后缀为 `pyc` 才排除——pathlib 中名字 `.pyc`（纯前导点）的 suffix 是空串、会被包含，R 版复刻了此边界行为；
   - 排除 manifest 自身（`paths != OUTPUT`）；
   - 仅保留文件（`file.exists & !dir.exists`，对应 `is_file()`）。
5. **排序键**：Python 按绝对 posix 路径字符串排序；因所有路径共享同一 `ROOT` 前缀，按相对路径字符串排序与之等价。R 用 `order(rel, method = "radix")`——radix 排序是字节序（对 UTF-8 即码点序），与 Python 字符串比较一致，且不受 locale 影响。PXD 目录排序同样用 radix。
6. **CSV 输出**：`write.csv(df, row.names = FALSE, quote = FALSE, eol = "\r\n", fileEncoding = "UTF-8")`，表头 `RelativePath,SizeBytes,SHA256`。
   - `eol = "\r\n"` 匹配 Python csv 模块默认行结束符（`newline=""` 下恒为 CRLF）；
   - `quote = FALSE` 匹配 `QUOTE_MINIMAL`（实测全树无逗号/引号/换行文件名，无需引用；已用 find 核查）；
   - 末尾行后同样有 `\r\n`（Python csv writer 每记录后都写行结束符），已用 `od -c` 逐字节确认表头与结尾。
7. **SizeBytes 预处理为字符串**（`as.character(file.info()$size)`）：树中最大文件 4.98 GB 超 2^31，R 整型会溢出；且 write.table 默认 `digits` 会对 double 做有效数字截断（如 4984509491 → `1.235e+09` 之类）。全字符串矩阵绕开这两个坑，数值逐位原样输出。
8. **去重**：`unique()` 对应 Python 的 `set()`（本树三个来源互不重叠，双保险）。

## diff 输出（验证结果）

两次背靠背完整对比（R 生成 → Python 生成 → diff）：

```bash
$ Rscript reanalysis/scripts/build_final_manifest.R . && mv .../final_file_manifest_sha256.csv /tmp/manifest_r.csv
Final files hashed: 498
Manifest: reanalysis/reports/final_file_manifest_sha256.csv

$ python3 reanalysis/scripts/build_final_manifest.py && mv .../final_file_manifest_sha256.csv /tmp/manifest_py.csv
Final files hashed: 498
Manifest: reanalysis/reports/final_file_manifest_sha256.csv

$ diff /tmp/manifest_r.csv /tmp/manifest_py.csv && echo MANIFEST_IDENTICAL
MANIFEST_IDENTICAL
```

复验：第二次完整运行同样 `SECOND_RUN_IDENTICAL`；两次运行间 R 版输出与 Python 版输出各自逐字节稳定（`cmp` 通过），排除文件集并发漂移干扰。额外用 `od -c` 确认：表头 `RelativePath,SizeBytes,SHA256\r\n`、行尾 CRLF、文件末尾 `\r\n`。

## 文件数

- manifest 数据行：**498 个文件**（+ 1 行表头 = 499 行）
- 旧 manifest 为 492：多出的 6 个是其他 agent 并行新建的文件（`analyze_kla_regulator_intensity.R`、`analyze_kla_regulator_whole_proteome_intensity.R` 等）及本任务自身的 `build_final_manifest.R`——属预期行为，manifest 本就应反映运行时刻的文件集。
- 样本抽检：`reanalysis/scripts/build_final_manifest.R,3604,...`（自身含入，3604 字节与实文件一致）；`reanalysis/reports/.Rhistory,0,e3b0c442...`（点文件正确含入）；`NEW_CHAT_PROJECT_PROMPT.md` 排在 `PROJECT_INDEX.md` 之前（码点序 `N` < `P`，排序正确）。

## 问题与处理

| 问题 | 处理 |
|---|---|
| R `list.files` 默认漏掉点文件（Python rglob 包含，如 `.Rhistory`） | `all.files = TRUE`；已抽检 `.Rhistory` 在 manifest 中 |
| `write.csv` 默认 `eol = "\n"` 与 Python csv 的 `\r\n` 不一致 | 显式 `eol = "\r\n"`；`od -c` 确认 CRLF 与结尾换行 |
| `write.csv` 默认会引用所有字符列，Python `QUOTE_MINIMAL` 只在必要时引用 | `quote = FALSE`（已核查全树无需引用字段） |
| 最大文件 4.98 GB > 2^31，R 整型溢出 / write.table `digits` 截断 double | SizeBytes 列预转 `as.character()`，全字符串矩阵输出 |
| `.pyc` 判断的 pathlib 语义边界（`suffix` 对前导点名字为空） | 自定义 `has_pyc_suffix()` 精确复刻 |
| 排序 locale 依赖 | `method = "radix"`（字节序 = 码点序，与 Python str 排序一致） |
| 并发 agent 可能改变文件集 | 两次背靠背运行均 `MANIFEST_IDENTICAL`，且 R/Python 各次输出相互 cmp 一致，确认无漂移 |

## 最终状态

`reanalysis/reports/final_file_manifest_sha256.csv` 已恢复为最新正确内容（R 版与 Python 版输出逐字节相同，任一份等价）。`reanalysis/scripts/build_final_manifest.py` 未被修改。
