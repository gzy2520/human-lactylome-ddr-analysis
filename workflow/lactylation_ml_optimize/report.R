library(data.table)
out <- 'outputs/20260922_br_v5'; audit <- 'audit/20260922_br_v5'
s<-fread(file.path(out,'performance_summary.csv'));u<-fread(file.path(out,'conditional_uncertainty.csv'));z<-fread(file.path(out,'s2_summary.csv'))
t<-fread(file.path(audit,'acceptance_tests.csv'))
lines<-c('# B-R v5：有限范围的概率优化', '',
'结论：有小幅平均预测改善，但没有解决材料之间的检出率差异。将 v5 保留为探索性候选，v4 继续保留为已验收的历史基线；不宣称独立验证优越性。', '',
'## 数据与设计',
'沿用 v4 的 6,744 行、2,708 个阳性、28 份 RNA 参考、19 个来源连通组。标签为 Ref 检出背景中的 Kla 流程检出；任务只适用于已知 DDR 蛋白跨材料预测，不能解释为生理乳酸化程度或实测乳酸。',
'全部使用 R，种子 25。没有增加 RNA 特征或修改标签。候选网格在本轮性能计算前写入 config/lactylation_ml_optimize/v5_contract.md：9 种频率平滑配置；分别使用原始平滑概率（SHRINK）和训练内交叉拟合的逻辑概率校准（SHRINK_CAL）。',
'外层留一连通组，内层 5 折按连通组等权 LogLoss 选参。每个内层训练集内部再做 4 折交叉拟合，以拟合校准器；最终校准使用外层训练集的 5 折交叉拟合概率。测试标签不能进入频率、校准或候选选择。校准固定惩罚 0.01*(slope-1)^2，不按外层结果调节。',
'19 个外层训练集、两种模型均选择 row_8：保留逐行蛋白历史频率，使用训练连通组等权阳性率作平滑中心，先验强度 8。候选组等权计票未被选中。选择一致只说明这组数据中的内层偏好，不表示强度 8 为普遍最优。', '',
'## 留出结果（19 组等权）',
'| 模型 | LogLoss ↓ | AP ↑ | ROC ↑ | Brier ↓ | 28材料校准MAE ↓ |',
'|---|---:|---:|---:|---:|---:|')
for(i in seq_len(nrow(s)))lines<-c(lines,sprintf('| %s | %.6f | %.6f | %.6f | %.6f | %.6f |',s$Model[i],s$LogLoss[i],s$AP[i],s$ROC[i],s$Brier[i],s$MaterialMAE[i]))
gain<-s[Model=='MP_cal',LogLoss]-s[Model=='SHRINK_CAL',LogLoss]
lines<-c(lines,'',sprintf('SHRINK_CAL 相对旧最佳概率基线 MP_cal 的平均 LogLoss 下降 %.6f（相对下降 %.2f%%）。AP 小幅上升，但仍低于旧 ME_P/ME2P 的 AP；材料平均概率校准 MAE 反而升高。因此优化只体现在部分预测指标上。',gain,100*gain/s[Model=='MP_cal',LogLoss]),'',
'## 稳健性与局限')
for(i in seq_len(nrow(u)))lines<-c(lines,sprintf('- %s：相对 MP_cal 有 %d 组改善、%d 组变差；平均 ΔLogLoss=%.6f。对已保存的 19 组差值重采样 10,000 次，95%% 分位范围 [%.6f, %.6f]；删一评价组的平均差值范围 [%.6f, %.6f]。',u$Model[i],u$Better[i],u$Worse[i],u$MeanDelta[i],u$BootstrapLow[i],u$BootstrapHigh[i],u$DeleteOneLow[i],u$DeleteOneHigh[i]))
lines<-c(lines,'以上重采样仅描述已保存预测的组间差异，外层训练集重叠且同数据已多次开发，不能视作独立样本置信区间或正式显著性检验。',
'在看到主结果后增加了删组重训的稳定性检查；未更改模型网格、惩罚或选参规则。对每个删除组，剩余 18 个测试组均重新完成整个内层选择与校准；共 342 次外层拟合、684 条模型评价。S2 的 MP 对照同样使用删除后的训练集重新计算，不复用原始 MP_cal 来冒充重训对照。')
for(m in unique(z$Model)){q<-z[Model==m];lines<-c(lines,sprintf('- S2 %s：19 个删组方案平均 LogLoss 范围 [%.6f, %.6f]；对同期重算 MP 的平均 ΔLogLoss 范围 [%.6f, %.6f]；相对原模型在相同剩余测试组上的平均变化范围 [%.6f, %.6f]。',m,min(q$MeanLogLoss),max(q$MeanLogLoss),min(q$MeanDeltaVsMP),max(q$MeanDeltaVsMP),min(q$MeanChangeVsOriginal),max(q$MeanChangeVsOriginal)))}
lines<-c(lines,'',
'## 可以怎么用',
'可将 SHRINK_CAL 作为本任务平均 LogLoss 较好的探索候选，同时保留更简单的 SHRINK 和 v4 MP/MP_cal。当前收益小且跨组不一致，不将新模型升级为生理乳酸化预测器，不继续按已看见的外层分数扩展搜索。',
'材料检出率的极端差异仍无法预测，说明当前候选没有解决来源间差异；仅凭此结果不能判断差异源自生物学还是检测流程。要回答老师提出的乳酸来源/去向与乳酸化关系，后续需要同平台可比定量与同一样本的代谢/乳酸、RNA、总蛋白及 Kla 数据。此轮没有重新检验或否定代谢机制。', '',
'## 验收与复现',sprintf('%d 项验收全部通过：包括外层标签扰动、内层隔离、交叉拟合标签隔离、独立 base R 逐组指标复算、v4 对照保持一致、S2 完整性。103 个 v4 清单文件哈希全部保持不变。两张图的 PDF 已渲染检查，无截字或遮挡。',nrow(t)),
'主结果：outputs/20260922_br_v5/；代码：workflow/lactylation_ml_optimize/。输入和交付哈希见本审计目录。',
'在隔离副本且移除副本内本轮输出目录后，运行 bash workflow/lactylation_ml_optimize/run_all.sh。入口会拒绝覆盖已存在的主预测；不要删除正式目录来原地重跑。输入与交付哈希、人工图件检查另见本审计目录。')
writeLines(lines,file.path(audit,'REPORT.md'))
