# 重复数缺口的替代数据落实（2026-09-15）

本轮是老师意见审核的后续解决方案。优先读取 `config/rnaseq_replicate_resolution_gate_20260915.csv`，原老师颜色和原始审核保留。服务器已实测通过 Tailscale 100.121.229.123 连接；工作目录 `/home/user/gzy/kla31-rnaseq-20260914`，复用原服务器已知主机密钥，不关闭主机身份验证。

## 已下载并验收的完整计数矩阵

|材料|来源|矩阵样本数|仍须保留的条件|
|---|---|---:|---|
|HK2（27/28）|GSE240748|3|WT未处理24h；独立培养细节待核，不代表mannitol处理|
|MCF7（13）|GSE157383|3|DMSO 7天，浓度未注明|
|MDA-MB-468（20）|GSE157383|3|DMSO 7天，浓度未注明|
|HepG2（18）|GSE158552|3|DMSO 24h，亲本不是SIRT1/3 KO|
|HUVEC（31）|GSE203551|4|独立供体32/36/46/47；未刺激未加药，P3|
|精子（07）|GSE65683|4|Group I，临床纳入及采样批次需进一步核对|
|非糖尿病撕裂肌腱（01）|GSE236746|3|撕裂肌腱，不是健康肌腱|
|PC-3M（22）|GSE235595|3|DMSO 4天|
|HCT116（14/15/16）|GSE253699|3|WT control；排除FSTL3；三个清单行共享这套参照，不是三个独立队列|
|RKO（24）|GSE318640|3|DMSO时长/浓度缺失；亲本不是GSK3B KO|

10套矩阵共32个样本列，覆盖13个清单行。前8套为NCBI生成的39,376个Entrez基因counts；后2套是作者提交的58,735个Ensembl基因counts。所有提取矩阵通过稳定ID格式、唯一性、有限非负整数、非空文库检查。记录原始文件MD5、精确GSM及原始列名，见 `resolved_matrices/`。未按Symbol做分析，也未填补缺失值。

NCBI入口此前漏检：肌腱与PC-3M并非必须重算FASTQ；现在已有所选样本的基因计数矩阵。原FASTQ保留作可追溯备份。官方说明：[NCBI RNA-seq counts](https://www.ncbi.nlm.nih.gov/geo/info/rnaseqcounts.html)。HK2作者原表存在长坐标字段断裂和非稳定ID，本轮改用NCBI计数，不再修补缺失为零。

精子原拟选7例，仅GSM1602978、GSM1602981、GSM1602982、GSM1602983在NCBI矩阵中，不能报n=7。其全基因Spearman为0.792–0.861，其他来源为0.908–0.967；这只是初步诊断，零表达并列会影响相关系数，不能据此直接判定生物学重复成立或删掉供体。尚需表达过滤后PCA、供体/培养重复与原蛋白材料核对。所有AnalysisReady仍为FALSE，技术验收通过与生物学最终通过分开记录。

## 仍未闭合的来源及具体处理

- A549：GSE171750所选GSM5233375–77是untreated三例，但作者RSEM表只有Symbol。NCBI计数列表两次未取得；不能断言不存在。下一步限定精确样本的原始定量或原注释映射；禁止直接以Symbol拼表。
- T47D：进一步找到GSE283812（属于GSE243418），GSM8672552–55四个DMSO生物学重复，RPMI1640+10%FBS。样本写24h，但通用处理描述写48h，时间有冲突。下载到完整计数表，行标仍是Symbol，不能直接验收。GSE310291四个4h对照保留次选，不混时间点凑重复。
- TALL-104：当前bulk n=1；三重复线索为芯片。解决方案是继续限定bulk来源；如无合格来源，该组只作描述性结果或从主要带误差线比较中排除，不伪造重复。
- MES28：原GSE266884同一shMCT1对照只有2例，不能拼入不同shKU70实验。新查GSE119834只有一个MES28库GSM3384841（标题GSC44），也不解决重复数。保留n=2描述性用途；独立来源不能直接当同批第三例。

以上四行仍有矩阵/重复数缺口；其他已准备矩阵的条件限制不等于已消失。此前13行不足n>=3的清单中，9行已取得至少三个样本列的稳定ID矩阵候选；剩余4行是A549、T47D、TALL-104、MES28。全31行仍需依次完成生物学资格审查，尤其癌旁组织、未标色GDC和NSC，不能以本轮13行替代全部31行结论。

## 复现与连接

脚本在 `workflow/`：`resolve_replicate_sources_20260915.R`获取SOFT元数据；`acquire_ncbi_counts_20260915.R`与`acquire_replicate_matrices_20260915.R`下载；两份`prepare_resolved_*_20260915.R`校验和提取；`update_replicate_resolution_gate_20260915.R`生成新门控表。矩阵单位保留counts，不混入TPM/FPKM，也不把不同研究合并成同批强度表。

原下载状态命令仍是 `bash workflow/check_server_rnaseq_download_20260914.sh`，已默认走Tailscale。实测FASTQ complete=12、md5_match=12、failed=0，下载进程STOPPED；GDC 1017/1017、DepMap 3/3完成。这是下载完整性，不是FASTQ比对或正式分析完成。
