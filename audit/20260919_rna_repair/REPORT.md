# RNA 修复与验收（2026-09-19 启动，2026-09-20 继续）

状态：数值重建、92 项验收、18 套图视觉检查和本地发布均已完成。

## 范围

按用户要求，保留原单因素与双因素 ANOVA，不调整显著性检验或配对模型。修复输入换算、缓存/输出安全、来源与样本元数据、qsmooth 诊断说明、图形排版、版本混用及发布可追溯性。旧冻结蛋白组数据与其分析不改。

## 已实施的修复

- FPKM/RPKM 只按列重标至 TPM，不二次除基因长度；有限/非负/正列和检查。
- Stage 2 必须使用空目录，仅同一次运行中共享 ReferenceKey 可复用；来源失败非零退出。部分提取不冒充完整汇总。
- Stage 3 在全部 28 对象、对应矩阵和合同样本数齐备后统一生成完整 manifest/QC/index。
- 新合同 `config/rnaseq_expression_extraction_contract_20260919.csv` 修正原海马 ExpectedN=200 为一直实际使用的 197，将 MCF10A 的文本 ~3 规范为整数 3；选择器与样本范围未变。
- GSE181540 实际以 counts 计算的路径正确记录为 counts；n=1 的相关性为 NA。
- 目录来源不再对目录本身取 MD5；单独记录 28 来源共 1,039 个实际原始文件的 MD5。
- qsmooth 算法保持不变；权重边车改用 QuantileRank，分开报告 A/B 权重、MeanAbsDelta 与 MaxAbsDelta。
- 核心 renderer 显式接收 expression/qsmooth/panel 路径，使用冻结 196 基因 Hallmark 映射；保留 28 材料/1,898 样本口径。取消写旧交付和桌面的隐式副作用。
- PDF 标题按画布换行，NA 副标题修复实际进入新输出；现有配色、面板顺序和计算方式保持。
- 以发布清单、SHA-256 和独立当前输入清单替换过时 preflight；保留旧冻结包哈希检查和旧脚本归档。
- 旧方法文档归档，当前文档据代码重写，去除重复长度公式、错误的误差上界、无依据的批次/QC 确证说法。

## 可复核证据

- `server_rebuild.log`：28 个参考矩阵从源重新提取，failures=0。最初 finalizer 被历史 ExpectedN 文本类型阻断；不是源数据处理失败。
- `server_finalize.log`：使用修订样本合同重新汇总/验收。
- `corrected_fpkm_verification.log`：独立源复算，BPH/TALL-104 新矩阵与正确公式的 log2 最大差及中位差均为 0。
- `source_fingerprints.log`：28 来源 / 1,039 个唯一原始文件。
- `regression_tests.log`：FPKM 长度无关性、非法输入、非空输出保护、不完整汇总阻断。
- `render_regression.log`：相同旧输入下，7 张核心表与原 28 材料版本一致；统计模型调用未变。
- `render_preview.log`：隔离预览，无隐式交付复制；预览不作为修正后正式结果。

初次 Stage 2 按历史 20260916 合同提取，原合同与部署代码指纹保存在 `extraction_input_fingerprints.csv`；汇总采用 20260919 修订合同并另存 `finalization_contract.csv` 及指纹。两个合同的选择器相同，仅修正 ExpectedN/说明字段。后续完整重跑入口统一使用新合同。

## 保留的边界

本次不改变统计模型，也没有重做每个来源的供体、建库和模型身份审查。来源批次、外部 RNA 与蛋白组非同一样本等局限仍需披露。数值修复后，按原算法得到的统计数值可能改变；这不是更换检验方法。

## 最终完成记录

- 当前发布：`outputs/20260919_repair_release/`；核心图在 `core/`，联合状态描述图在 `reference/`。单一入口为 `config/rna_current_release.csv`。
- 92 项验收全部通过，18 套 PDF/PNG（14 套核心图 + 4 套描述图）齐备，渲染的 PDF 完成视觉检查。
- 28 个源矩阵逐元素比对：仅 BPH/TALL-104 改变（最大绝对 log2 变化 6.501166 / 5.219513），其余 26 份与旧版保持一致；共同基因 18,332，样本 1,898。
- 修正后 A/B 最大绝对差 1.38658615；14,375 个基因的跨材料最大差超过 0.5，而没有基因的跨材料平均绝对差超过 0.5。两种指标分别披露。
- 除核心标题换行和 NA 修复外，视觉 QA 修复联合柱图低表达小区块标签错放/裁切，以及 RNA_1 小通路条和窄类别标题裁切；3 张描述表逐行保持不变。
- 当前 preflight 通过：163 个原冻结文件、17 个当前输入、174 个发布/依赖文件哈希核验；28 材料、1,898 样本、31 行对应关系一致。
- 桌面当前交付：`/Users/gzy2520/Desktop/renew/kla/RNA/`；旧桌面 RNA 目录及旧顶层 RNA 图移入同目录下的 `archive/rna_before_20260919/`。蛋白组图未动。
- 原 `results/rna_teacher_questions` 已移入 `results/archive/rna_teacher_questions_before_20260919/`，该熟悉入口现在链接到新 core；旧 dated delivery 保留并加历史提示。
- 服务器保留新 expression 对象及全部源文件指纹；新发布和验证后的入口脚本同步到隔离项目。原入口脚本备份在 `workflow/releases/20260920_previous_entrypoints/`，验证版在 `workflow/releases/20260920_validated/`。
- 探索脚本的 Spearman ties 提示属于本次明确不调整的统计部分；图形的字体度量及失效 nudge 参数警告已修复。
- 同步后逐文件复核：服务器 60 个交付文件及 30 个 qsmooth/panel/assisted 支撑文件与本地 MD5 全部一致；证据分别见 `server_delivery_verification.log` 和 `server_supporting_verification.log`。
