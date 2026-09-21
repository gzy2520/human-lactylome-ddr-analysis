library(data.table)
O<-'outputs/20260921_br_v4'; A<-'audit/20260921_br_v4'
versions<-list(v3='outputs/20260921_br_v3',weight_only='outputs/20260921_br_v4_weight_only',corrected_L=O)
comparisons<-list(); deltas<-list()
for(v in names(versions)){
 dir<-versions[[v]];tag<-if(v=='v3') 'BRv3' else 'BRv4'
 p<-fread(file.path(dir,paste0(tag,'_oos_preds.csv.gz')))
 d<-p[,.(ME_P=-mean(y*log(p_ME_P)+(1-y)*log1p(-p_ME_P)),ME2P=-mean(y*log(p_ME2P)+(1-y)*log1p(-p_ME2P))),by=ConnGroup]
 d[,`:=`(delta=ME2P-ME_P,Version=v)];deltas[[v]]<-d
 perf<-fread(file.path(dir,paste0(tag,'_performance_summary.csv')))[Summary=='equal_weight_19folds']
 perf[,Version:=v];comparisons[[v]]<-perf
}
fwrite(rbindlist(comparisons,fill=TRUE),file.path(O,'version_performance_comparison.csv'))
d<-rbindlist(deltas)
fwrite(d,file.path(O,'module_paired_fold_deltas.csv'))
cnt<-d[,.(MeanDelta=mean(delta),Worse=sum(delta>1e-12),Better=sum(delta< -1e-12),Tied=sum(abs(delta)<=1e-12)),by=Version]
fwrite(cnt,file.path(O,'module_delta_counts.csv'))
stopifnot(cnt[Version=='v3',Worse]==9,all(cnt$Worse+cnt$Better+cnt$Tied==19))
md<-function(x){paste(c(paste('|',paste(names(x),collapse=' | '),'|'),paste('|',paste(rep('---',ncol(x)),collapse=' | '),'|'),apply(x,1,function(z)paste('|',paste(z,collapse=' | '),'|'))),collapse='\n')}
tab<-comparisons$corrected_L[,.(Model,LogLoss,AP,Brier,ROC)];tab[,2:5]<-lapply(tab[,2:5],round,6)
s2<-fread(file.path(O,'S2_train_sensitivity_summary.csv'))
s2brief<-s2[,.(MeanPairedDelta=mean(diff_eqw_LL),MinPairedDelta=min(diff_eqw_LL),MaxPairedDelta=max(diff_eqw_LL)),by=Model]
qa<-fread(file.path(O,'acceptance_tests_record.csv'))
extras<-fread(file.path(O,'additional_acceptance.csv'))
stopifnot(all(qa$Status=='PASS'),all(extras$Status=='PASS'))
report<-paste0('# B-R v4 修复与结果（2026-09-21）\n\n',
'本轮修复内层评价权重、模块生物学定义和报告计数。保留全部v2/v3原始结果。当前结论仅适用于Ref检出背景下的已知DDR蛋白流程检出预测，同数据集探索性分析，不是独立外部确认。\n\n',
'## 修复与勘误\n\n',
'- v3报告的“15/19折变差”为手工写入错误，直接预测复算为9/19。v4所有计数来自逐折差值表；正数表示变差。\n',
'- 内层改为每个验证连通组各计算LogLoss，再对全部训练连通组等权汇总。不是先按蛋白行平均再平均5个内层折。保留每组损失供独立重算。\n',
'- 源矩阵独立复算确认6744行目标表达和28个旧模块均值一致，无错位。备用LDHA ID ENSG00000288299在全部28个源矩阵中缺失；旧“6个可用成员”误将NA占位行当作可用。实际旧均值由5个成员构成。\n',
'- 新模块包含LDHA/LDHB/LDHC/LDHAL6B四个不同基因，名称为L-乳酸/丙酮酸转换表达，不能解释为单向生成量。LDHD仅另列D-乳酸氧化表达，不参与本轮新增模型。\n',
'- 独立ID证据：两个LDHA ID对应同一Entrez 3939 / UniProt P00338。采用NCBI参考ID ENSG00000134333，未依据预测性能选ID。\n\n',
'## 三版本受控比较\n\n',md(cnt), '\n\n',
'v3为旧混合模块与旧内层权重；weight_only仅修正权重、保留旧特征；corrected_L同时修正权重及L/D模块定义。完整指标见version_performance_comparison.csv。\n\n',
'## 当前外层指标（19连通组等权）\n\n',md(tab),'\n\n',
'模型定义沿用v3，ME2/ME2P/XEN中的模块现为四基因L-乳酸/丙酮酸转换。模块增量必须比较ME_P与ME2P，不能用系数符号替代逐折增量。\n\n',
'## 结果解释\n\n',
sprintf('当前LogLoss最低模型为%s（%.6f）；ME_P与ME2P分别为%.6f和%.6f，平均变化%+.6f。修正模块的逐折LogLoss为%d组变差、%d组改善，不能概括为跨组一致。AP从%.6f至%.6f，排序与概率误差指标应分开解释。\n\n',tab[which.min(LogLoss),Model],min(tab$LogLoss),tab[Model=='ME_P',LogLoss],tab[Model=='ME2P',LogLoss],cnt[Version=='corrected_L',MeanDelta],cnt[Version=='corrected_L',Worse],cnt[Version=='corrected_L',Better],tab[Model=='ME_P',AP],tab[Model=='ME2P',AP]),
'在本任务和预定特征表示下，该转换表达模块未改善平均LogLoss；不能据此断言乳酸代谢无作用、模型达到上限或负系数揭示了因果抑制。\n\n',
'## S2同剩余组配对变化（各删除组的等权LogLoss变化）\n\n',md(s2brief),'\n\n',
'## 验证和范围\n\n',
'- 17项继承验收及新增验收通过；完成实际标签扰动测试，仅说明所测试的隔离条件通过，不声称已证明所有可能的信息泄漏均不存在。\n',
'- S1为删一测试组后重新汇总；S2对MP、ME_P、ME2P完成19个删除组×18个外层折，使用同一内层调参引擎，与原结果在相同剩余组上比较。其它模型没有本轮S2，不扩大宣称。\n',
'- 调参与评价连通组等权；模型训练似然仍按蛋白行等权，未改变v3拟合目标。仍可能受研究流程及蛋白可检出倾向影响。\n',
'- 只评估这个预定义转换表达模块，不能概括全部乳酸代谢、真实乳酸浓度或乳酸化调控机制；不根据负系数推断消耗、抑制或因果机制。\n',
'- 旧稿“高度稳健”“15折变差”“因为无法泛化导致校准恶化”等超出计算证据的表述不再沿用。\n\n',
'## 证据来源\n\n',
'- [LDHD原始酶学研究](https://pubmed.ncbi.nlm.nih.gov/37863926/)：D-乳酸氧化。\n',
'- [人类LDHD缺陷研究](https://pubmed.ncbi.nlm.nih.gov/30931947/)。\n',
'- [NCBI LDHA](https://www.ncbi.nlm.nih.gov/gene/3939)与项目冻结UniProt映射。\n',
'- [LDHAL6B UniProt](https://www.uniprot.org/uniprotkb/Q9BYZ2/entry)：功能注释，非本项目直接酶活测量。\n')
writeLines(trimws(report, which = 'right'),file.path(A,'REPORT.md'))
writeLines(c('# Current B-R v4 delivery','Primary results: corrected four-gene L-lactate/pyruvate conversion module.','See audit/20260921_br_v4/REPORT.md and config/lactylation_ml_review/model_config_v4.md.','Historical v3 figures/report are retained as prior versions; do not use the old 15/19 claim.'),file.path(O,'README.md'))
cat(md(cnt),'\n',md(tab),'\n')
