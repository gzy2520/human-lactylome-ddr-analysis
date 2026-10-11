from pathlib import Path
import csv, json, hashlib, zipfile
root=Path('/Users/gzy2520/Desktop/Research/kla/data');out=Path(__file__).resolve().parent
rows=[
('PXD075377','adjacent liver;HCC','tissue_sample','1 per tissue','No','Intensity Control;Intensity HCC','metadata Description: paired tissues from one clinical patient; sample.txt only HCC/Control','search_results/extracted/2-Basic_analysis/MS_identified_information.txt','旧pool分类改为一名患者的两类组织样本；每类n=1'),
('PXD028488','HCT116','single_sample','1','No','Source File=20211127-HCT116-LAC-Enrichment-1.raw','One source raw file in DB search psm.csv','search_results/Enrichment-Search files/HCT116-Enrichment-Search files/DB search psm.csv','逐PSM重新提取，有一个原始样本'),
('PXD028488','HEK293T','technical_sourceunit','1 merged sourceunit','No','HCD24;HCD28;HCD32;HCD36;HCD40','Five collision energies; raw names contain mix-HCDxx; no biological replicate column','search_results/Enrichment-Search files/HEK293T-Enrichment-all HCD-Search files/DB search psm.csv','各独立HCD目录解析后合并；all HCD目录为技术合集仅审计；生物学重复不能拆分'),
('PXD053474','HCT116','technical_sourceunit','1 merged sourceunit','No','DIA F1/F2/F3/F4; compartments; S3 Subcellular and Whole-cell lysates','All fractions/acquisition modes collapsed; author S3 plus DDA and DIA site corroboration','supplementary/js4c00366_si_003.xlsx','原始数据重建590蛋白；亚细胞结构与色谱分级不能作为重复；S2是既往文献汇总'),
('PXD058534','pretreated HK-2','single_sample','1','No','Experiment HK2; raw FE458LPLa_HK2_Slot1-35_1_14667','summary.txt one non-Total raw row','search_results/extracted_pairing/txt/La (K)Sites.txt','一个样本可计算比例，无SEM'),
('PXD070007','neural stem cells;glioblastoma stem cells','model','2 NSC;6 GSC models','Model variation only','ENSA;HMP1;G2907;G3028;G3264;GSC23;MES28;RKI','Annotation_Combine has model intensities; no replicate/raw columns; local ZIP contains same single workbook','search_results/SA206LPLaB1_Annotation.xlsx','可计算各模型比例；未找到模型内生物学样本列；如展示误差应明确是模型间变异'),
('PXD073311','HUVEC control and Pg infection','sample','3 A0h;3 A6h','Yes within condition','A0h_1.raw;A0h_2.raw;A0h_3.raw;A6h_1.raw;A6h_2.raw;A6h_3.raw','Distinct sample raw-file intensity columns in pr_matrix','search_results/extracted_pairing/IPX0015307003_Database_search_result/Database_search_result/report.pr_matrix.tsv','保留全部6样本及条件字段；旧condition标签错误遮蔽重复信息'),
('PXD066351','HCT116 control and Roseburia co-culture','single_sample_per_condition','1 NC116;1 R116','No within condition','R.Condition;R.FileName;R.Fraction;R.Replicate','Two raw files; R.Fraction NA; R.Replicate 0; DIA_result confirms only 2 acquisition columns','search_results/XB01472B1DPLa-MSstats_Input.csv','两条件各单样本；保留各自比例，不把条件差异当同条件重复'),
('PXD078013','RKO WT and GSK3B KO','sample','2 WT;2 KO','Yes within condition','Experiment sg_1;sg_2;sg737_1;sg737_2','metadata states two replicates each; SDRF lists 4 sources with biological replicate blank','search_results/evidence.txt','保留4样本；evidence Kla-positive site IDs与proteinGroups连接，蛋白组代表ID遵循历史脚本')]
keys=['PXD','SampleGroup','UnitType','AvailableUnits','ErrorbarSupport','ColumnsOrRawFiles','Evidence','RelativeSource','RecoveryConclusionZh']
with (out/'source_ledger.csv').open('w') as f:
 w=csv.writer(f);w.writerow(keys);w.writerows(rows)
# Inventory all local alternative source files, including archives.
with (out/'source_file_inventory.csv').open('w') as f:
 w=csv.writer(f);w.writerow(['PXD','RelativePath','Bytes'])
 for px in sorted(set(x[0] for x in rows)):
  for p in sorted((root/px).rglob('*')):
   if p.is_file() and p.name!='.DS_Store':w.writerow([px,str(p.relative_to(root)),p.stat().st_size])
with (out/'source_headers.json').open('w') as f:
 d={}
 for row in rows:
  p=root/row[0]/row[7]
  if p.suffix in ['.csv','.tsv','.txt']:
   with p.open(errors='replace') as src:d[str(p.relative_to(root))]=next(csv.reader(src,delimiter=',' if p.suffix=='.csv' else '\t'))
 json.dump(d,f,ensure_ascii=False,indent=2)
