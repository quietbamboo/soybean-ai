plink               ：  /nfs/my/Huang/lzm/default/plink/plink  
PopLDdecay          :   /nfs/my/Huang/lzm/default/PopLDdecay-3.42/bin/PopLDdecay  
Plot_OnePop.plink   :   (perl) /nfs/my/Huang/lzm/default/PopLDdecay-3.42/bin/Plot_OnePop.pl  
admixture           :   /nfs/my/Huang/lzm/default/admixture_linux-1.3.0/admixture  
S2_gemma.sh         :   /nfs/my/Huang/lzm/soybean_C219/code/data_processing/S2_gemma.sh  
S3_assoc2CMplot     :   /nfs/my/Huang/lzm/soybean_C219/code/data_processing/S3_assoc2CMplot.py  
S4_CMplot.R         :   /nfs/my/Huang/lzm/soybean_C219/code/data_processing/S4_CMplot.R  
S5_p_wald_filter.py :   /nfs/my/Huang/lzm/soybean_C219/code/data_processing/S5_p_wald_filter.py  
S6_p_wald_shap.py   :   /nfs/my/Huang/lzm/soybean_C219/code/data_processing/S6_p_wald_shap.py  

01、S1_hmp2ped.py 将219_snp.hmp转换为plink格式  
02、根据有表型的样本数据，将数据集使用五折划分  
03、分别对每一折训练集中的样本进行质量控制，QC：--geno 0.05 --maf 0.05  
04、分别对每一折训练集中的样本进行LD分析，并筛选出一批相对独立的SNP样本，LD：--indep-pairwise 50 5 0.2  
05、S2_gemma.sh 使用gemma对每种性状进行GWAS分析，LD前后都进行分析，可使用run_gemma.sh脚本  
06、S3_assoc2CMplot.py 整合每个性状连续3年的 GWAS p值到CMplot格式，使用run_assoc2CMplot.sh脚本  
07、S4_CMplot.R 使用R包 CMplot 绘制曼哈顿图、Q-Q图、密度图等  
08、S5_p_wald_filter 筛选出 p_wald < 0.01 的SNP位点的 rs，2015年与3年的交集  
09、使用 plink 通过 rs 提取数据集，每一折的训练集、测试集，15年和3年的交集，分别保存为 .ped 和 .raw 格式  
10、S6_p_wald_shap.py 筛选用来 SHAP 分析的最优折的符合条件的位点，保存到一个文件后继续添加 SHAPLEY 值  
