#!/bin/bash
# ./gemma.sh train_avg.txt num plink pheno
# train_avg.txt 保存为所有的表型数据，第一行为列名，第一列为样本名字，后面的列为每种性状的表型值
# num 为计算的性状在train_avg.txt的第几列（从第2列开始）
# plink 为输入SNP的plink二进制格式的前缀（plink.bed plink.fam plink.bim）
# pheno 为最终输出文件名字的前缀，建议使用性状名字命名

# 检查输入参数的数量
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <trait_avg> <num> <plink> <pheno>"
    exit 1
fi




# 定义输入和输出文件
trait_avg=$1
num=$2
plink_prefix=$3
output_prefix=$4

# 检查输入文件是否存在
if [ ! -f "$trait_avg" ]; then
    echo "trait_avg file not found: $trait_avg"
    exit 1
fi

if [ ! -f "$plink_prefix.bed" ]; then
    echo "bed file not found: $plink_prefix.bed"
    exit 1
fi

if [ ! -f "$plink_prefix.fam" ]; then
    echo "fam file not found: $plink_prefix.fam"
    exit 1
fi

if [ ! -f "$plink_prefix.bim" ]; then
    echo "bim file not found: $plink_prefix.bim"
    exit 1
fi

# 工具路径
plink="/nfs/my/Huang/lzm/default/plink/plink"
gemma="/nfs/my/Huang/lzm/default/gemma-0.98.5-linux-static-AMD64"

# 使用 awk 筛选第 num 列表型不为0的样本ID
# awk 'NR==1 {next} $num != 0 {print $1}' $train_avg > temp_sample.txt
# awk -v col="$num" 'NR==1 {next} $col != 0 {print $1}' "$trait_avg" > temp_sample.txt
awk -v col="$num" 'NR==1 {next} {
    val = $col + 0;  # 强制将值转换为数字类型
    if (val != 0) print $1
}' "$trait_avg" > temp_sample.txt

# 复制样本ID，生成第二列
paste temp_sample.txt temp_sample.txt > sample.txt

# 删除临时文件
rm temp_sample.txt

# plink 提取有表型数据的样本的 SNP
$plink --bfile "$plink_prefix" --keep sample.txt --make-bed --out trait

# 将表型值放入 .fam 文件中
# awk 'NR==FNR {sample[$1]=$num; next} $1 in sample {$6=sample[$2]} {print $0}' $train_avg trait.fam > updated_trait.fam
awk -v col="$num" 'NR==FNR {sample[$1]=$col; next} $1 in sample {$6=sample[$1]} {print $0}' "$trait_avg" trait.fam > updated_trait.fam
mv updated_trait.fam trait.fam

# gemma 跑亲缘关系
$gemma -bfile trait -gk 1 -o $output_prefix

# 获取生成的亲缘关系矩阵文件路径
kinship_matrix="./output/${output_prefix}.cXX.txt"

# 进行关联分析
$gemma -bfile trait -k "$kinship_matrix" -lmm 1 -o $output_prefix

# 删除中间文件
rm trait.* sample.txt