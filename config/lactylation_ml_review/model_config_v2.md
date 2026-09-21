# 冻结：B-R v2 模型比较配置（2026-09-22，在查看任何新模型外层结果之前写定）
# 任务 B-R：行=(GroupID, DDR蛋白)|RefDetected；y=该流程下Kla检出；19连通分组留一。种子25。
# 标签为流程检出结果；不称浓度/位点占有率/生物学真值。
Model,FeatureSpec,Formula,Rationale,PrimaryMetric,Note
M0,空,常数=训练折总体正例率,下界基线,logloss,"常数预测AP=正例率；ROC=0.5"
MP,蛋白检出倾向,训练折Ref背景下该蛋白Kla检出频率(Beta先验a=1,b=1,即在Ref背景二项计数上加1/2),竞争基线:蛋白固有可检出倾向,logloss,"外层测试行用全部外层训练数据；训练行频率用内层5折连通分组交叉拟合防自泄漏；蛋白未见于训练→回退M0并记录"
ME,目标基因表达,train-z(target_expr),d11220e M1对应,logloss,
ME2,目标基因+乳酸生成模块,"train-z(target_expr, mact_generation_lactate_mean)",机制增量:乳酸生成,L-logloss,模块均值材料级
ME2P,目标基因+乳酸生成+MP,train-z(target_expr, mact_generation_lactate_mean)+MP频率作特征,MP与机制组合,logloss,MP频率仍隔离构造
MP_only,MP频率,同MP,MP单独检验,logloss,
XEN,交互探索(仅1项),"train-z(target_expr) + train-z(mact_generation_lactate_mean) + 乘积交互",预先声明唯一交互;系数正则化由内层5选1完成,logloss,"内层不足时固定lambda=0.1并声明;不得因外层结果追加交互"
# 预固定正则化候选(内层5折连通分组选择)：ridge lambda ∈ {0.05,0.1,0.3,1.0}（按n缩放: lambda*n/1000 用于glm目标)
# 特征构造规则（冻结）：材料级基因值=UNSMOOTHED log2TPM逐参考材料内样本中位数（该步不依赖训练折，可全量计算——仅是材料内聚合，非跨材料统计）
# 标准化：z中心/尺度仅由训练折的“不重复参考材料”集合估计（每材料计1次，不按蛋白行重复加权）；应用至验证折。
# 缺失：材料级特征NA时用训练折参考材料中位数填补（填补值亦为训练折估计）；行级target_expr NA行已在02记录（7目标基因无矩阵值→对应行剔除，同一集合全体模型可比）。
# 模块均值：模块成员log2TPM材料中位数的均值；方向不加正负号（双向反应）；成员及方向已于本轮开跑前冻结于 modules_frozen_v2.csv。
# 评价：主指标=外层logloss（材料级模块只平移材料概率，主要影响概率校准而非折内排序，故主指标不设ROC/AP）；同报Brier、AP（阈值聚合修正版）、ROC-AUC（材料内不敏感，仅列参考）。
# 汇报口径：每折、19折等权、行数加权、pooled（pooled仅参考，禁用跨折常数比较的pooled ROC）。
# 两类敏感性：(S1)评价汇总对单个测试组的敏感性=删一测试折重算均值；(S2)训练影响=删一连通组后重新执行全部外层验证（重算预处理/频率/调参），仅在剩余组上与原模型对比。
