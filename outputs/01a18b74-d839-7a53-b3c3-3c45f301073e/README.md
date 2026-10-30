# 图a、b与DDR七通路逐点数据

本目录保存2026年10月30日提供的图a、b及通路图d–j对应的历史作图数值。

- 图a：272个数据点；Whole proteome四类为60、107、13、20，Kla为31、11、9、21。
- 图b：31个GroupID的增殖得分，四类为9、3、7、12；保留原图共享RNA ReferenceKey的条目。
- 七通路：BER、NER、MMR、FA、HR、AEJ、NHEJ，每通路72个样本，每样本同时保留Pro、Inh。
- 总计1311个数值。四类依次为非肿瘤组织、肿瘤组织、正常细胞系、癌细胞系。

主工作簿为`GraphPad_图a_b及七通路_逐点数据_20261030.xlsx`。图a及通路采用Grouped布局：四类为列组，组内每个子列对应一个观测，行分别为Whole/Kla或Pro/Inh。图b为四列Column表。空白位置保持为空，保留实际零值。`逐点明细`提供单元格位置、样本标识、来源及计数分子分母。

图a数值为DDR检出蛋白数除以该来源观测的检出蛋白数，再乘100；通路数值为通路Pro或Inh蛋白数除以该样本Kla-DDR蛋白数，再乘100。图b直接使用历史源表的`ProliferationMedian`：每RNA样本对196个G2M_CHECKPOINT基因的log2(qsmooth TPM+0.5)取均值，再按历史GroupID取样本中位数。

`source_snapshot/`保留三份作图源CSV。比例与通路来源为`corrected_final_result_20260905_sample_only/Corrected_Data/sample_only_inputs/`；增殖得分来源为`outputs/20260917_delivery_31group/rna/`。源文件SHA256列在工作簿使用说明中。

复现步骤（在项目根目录执行，Python与Node路径按运行环境指定）：

```bash
python3 workflow/export_teacher_prism_screenshot_20261030.py --input-dir outputs/01a18b74-d839-7a53-b3c3-3c45f301073e/source_snapshot
node workflow/build_teacher_prism_screenshot_20261030.mjs
```

工作簿生成依赖`@oai/artifact-tool`。数值整理使用Python标准库；未进行新的显著性检验。脚本检查原图各组样本数、观测唯一性、计数分子分母与百分比的一致性，并逐一核对工作簿内1311个数值。交付前读取导出的xlsx进行独立核对，最大比例计算误差为4.98×10⁻¹⁴个百分点。
