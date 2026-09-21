# 乳酸来源与去向探索验收

开始于2026-09-20，完成于2026-09-21。接受用户“先按计划执行”的授权。本轮不改变旧显著性检验、不训练模型。

## 交付

- `outputs/20260920_lactate_metabolism/`：128个实测稳定ID基因、28材料、1898样本标识，23套PDF/PNG、中文汇报、方法、证据与校验清单。
- 6个完整的人类GO查询集合，10个人工机制模块，GO与人工成员独立记录。1275条请求范围内的注释记录；全部分页和is_a/part_of关系检查通过。
- 51个人工机制基因的实测Ensembl映射均经当前reviewed UniProt复核；全部蛋白GO映射补查了当前UniProt，未把ComplexPortal当作基因复制。
- 140个候选/替代Ensembl ID中128个进入共同矩阵。12个未覆盖ID逐来源核查；缺失未填0。
- 保留RNA实际对照条件标签及原31条蛋白组材料映射。没有用历史跳号RowOrder生成GroupID。

## 验证

1. 分析入口10条结构/数值检查全部通过，见交付tables/validation.csv。
2. GO下载、分页和关系检查6个集合均通过，见tables/go_download_validation.csv。
3. 独立空目录重跑后，17张数据表完全一致；压缩表按解压内容比较。见reproduction_comparison.log。
4. 全部23个PDF已渲染并查看，包括更正实际RNA标签后的最终版本。visual_qa.csv绑定PDF哈希；visual/保存预览与联系表。
5. 原冻结163文件、17当前输入、174既有发布依赖的preflight在复现前后均通过。未变动既有标准化、旧统计或旧交付。
6. input_sha256.csv记录脚本/配置/矩阵等输入；release_sha256.csv记录新交付文件。

## 开发中发现并修正

- 初稿运输GO号码经官方名称检查确认为胆汁酸运输，拒绝纳入，改用GO:0015727；原查询证据留在rejected_term/。
- 计划中的+1与当前代码不符，实际沿用log2(TPM+0.5)。
- GroupID不能从RowOrder派生：改为PXD+SampleGroup验证连接，完整31/28断言通过。
- G6PC旧别名改按安装注释中的当前G6PC1解析，并用当前UniProt复核。
- UniProt映射查询排除非UniProt复合物ID；复合物注释仍保留在原注释及缺失映射表。
- 图中文字使用RNA实际对照条件，避免暗示RNA来自KO/感染样本；示意图注明双向LDH/MCT和非通量性质。

这些问题均在新分析构建/验收过程中处理，未改变上一轮验收数据。

## 解释边界

重点对比为HCC/癌旁肝和TCGA-PRAD/BPH。BPH为有背景用药的良性病变队列，跨研究比较不是匹配正常对照。RNA不是与质谱同样本的实测；KlaProteinCount是材料蛋白并集检出数。未将其叫作乳酸化强度或位点占有率，未计算虚假的样本级相关。

没有乳酸浓度或可比独立Kla定量标签，不能证实“竞争分流”或“写擦稳态”；模型训练条件尚不满足。这些是本轮可行性结论。

## 发布副本

桌面新目录为 `/Users/gzy2520/Desktop/renew/kla/RNA_lactate_metabolism_20260921/`，原RNA目录未覆盖。桌面交付及实验室服务器的新交付逐文件SHA256复核通过；服务器同时保存本阶段脚本、冻结配置和方法，见desktop_verification.log及server_verification.log。新输入129个文件、交付88个文件（清单本身另存）在本地再次核验通过。
