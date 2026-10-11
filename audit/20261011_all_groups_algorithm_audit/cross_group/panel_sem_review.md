# 来源观察单位及processed matrix面板处理复核

## 审核结论

拟议处理符合本轮原表证据：修正4个细胞Kla位置7个样本中的反向序列分母；将4个已处理完整量化矩阵各保留为1个作者panel；保留来源列ID/数目但不计算样本SEM。当前325条展示记录据此变成307条，62个柱位不变，有SEM位置40→36。独立计数核验见 `panel_plan_independent_check.csv`，复核脚本为 `panel_plan_check.R`。

## 四个panel的处理

| GroupID及组学 | 原来源列数 | panel分子/分母 | 处理依据 |
|---|---:|---:|---|
| 05 Hippocampus Kla | 3 | 29/853 | S3为作者筛选quantifiable Class I位点panel；各列全正值，同一最小值重复73/70/1279次；论文所述原始检出并非全部3样本共有 |
| 05 Hippocampus whole | 3 | 219/6082 | S4表头注明相对量化值；三列6082行全正值，共同最小值分别重复84/69/594次 |
| 09 Adjacent liver whole | 8 | 194/4251 | ANTs为完整iTRAQ相对量化矩阵，Non-rec1列4254行全部为1，其他列全正值；它保留比值而非独立检出mask |
| 12 HCC whole | 8 | 271/5188 | CISs为完整iTRAQ相对量化矩阵，Non-rec1列5198行全部为1，其他列全正值 |

原始行数与去重后稳定ID数不同：ANTs4254行→4251个accession，CISs5198行→5188个accession；保留原稳定ID去重规则。四个位置原来源列的蛋白集合和占比各自完全相同，collapse不会改变柱高。325−[(3−1)+(3−1)+(8−1)+(8−1)]=307。

建议每个panel明确设 `UnitType=processed_matrix_panel`、`SemEligible=FALSE`、`N=1`、`SD/SEM/PlotSEM=NA`。`SourceObservationIDs`完整保留旧列名，`SourceObservationN`为3/3/8/8；它记录矩阵来源列数，不用来继续计算SEM。原始列值和旧样本membership留在历史审计文件，panel membership仅按稳定ID保留一套，不能把同一集合复制3或8次再对外记为样本重复。图注为panel增加一种简短标记并注明作者已处理量化矩阵即可。

不能仅凭重复最小值推定具体插补软件或算法；panel处理的依据是来源表性质以及其未保留原始样本检出mask。此次四panel处理不需要重新搜索原始质谱数据。

## 反向记录与条件标签

细胞独立审计支持先按原始代表Protein/Leading proteins排除REV__记录，再执行原样本检出、稳定ID和去重规则。HepG2、HMC3、pretreated HK2、HK2 ctr/man共7个样本各少1个分母蛋白，DDR分子不变；其他正向原始行未提供被删ID的合格支持。故这是原始记录修正，不是全局删除某个合法人蛋白。

HK2 man的正式来源只支持ctrl/man命名。显示材料名和当前Condition应使用原始ctr/man，原SampleGroup中的mannitol仅作为历史关联键保留；不应在新增显示说明中继续赋予mannitol或hypoxia-reoxygenation释义。这样不影响既有稳定分组ID及旧版对照。

## 其他目前无需强制重拆的来源单位

- Tendon的reporter正值表有明确零值和非零值；RCT2/3蛋白集合相同由原表检测结果产生，不能单凭相同集合把所有TMT来源折成panel。当前3个reporter通道保持；___1/2/3是位点量化multiplicity，不作为3个独立重复。原全鉴定union与通道正值集合的统计对象变化已单独记录。
- HMC3和Roseburia每条件仅1个来源观测，继续不画SEM；合并均值仅为对应既定材料柱的来源均值，Condition仍保存两个条件。
- NSC/GBM多个模型保持模型单位，不画样本SEM；不同模型不得折作技术重复。
- HEK293T碰撞能量、HCT116亚细胞组分与fraction已合并；无需再次拆成生物样本。
- 全蛋白组technical_preparation保留其技术制备标记t，SEM描述这些制备观察值；没有新的生物样本证据时不更名为患者重复。
- ESCC部分样本membership相同不足以判为复制panel；目前已审计原始T样本列有逐列不同检出集合与数值，保留94个来源列。
- Adjacent liver/HCC Kla仍是同一患者的两个组织条件，各自一个单例；这次collapse仅改独立参考队列的whole proteome位置，不把两种组学误当同一患者配对。

除已确认的反向记录、四panel和HK2条件释义外，本轮现有证据未要求进一步改变展示位置或样本拆分。最终release生成后仍需重新核验307条记录、62位置、36处SEM、panel来源ID保留、修正membership及RNA数值不变；本文件为处理方案审核，尚不宣称最终文件通过。
