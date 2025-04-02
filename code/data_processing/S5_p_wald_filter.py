import sys
import pandas as pd

# input_dataset         输入数据集
# p_value_threshold     p值
# output_15_years_file  保存15年筛选数据的文件名
# output_3_years_file   保存3年交叉结果的文件名

# 检查命令行参数数量是否正确
if len(sys.argv) != 5:
    print("Usage: python script.py <input_dataset> <p_value_threshold> <output_15_years_file> <output_3_years_file>")
    sys.exit(1)

arg1 = sys.argv[1]
arg2 = sys.argv[2]
arg3 = sys.argv[3]
arg4 = sys.argv[4]

try:
    # 将阈值参数转换为浮点型
    thrd = float(arg2)
except ValueError:
    print("Error: <p_value_threshold> must be a numeric value.")
    sys.exit(1)

assoc = pd.read_csv(arg1, sep='\t')

assoc_15 = assoc[(assoc.iloc[:,3] < thrd)][['rs']]
assoc_3 = assoc[(assoc.iloc[:,3] < thrd) & (assoc.iloc[:,4] < thrd) & (assoc.iloc[:,5] < thrd)][['rs']]

assoc_15.to_csv(arg3, header=None, index=None, sep='\t')
assoc_3.to_csv(arg4, header=None, index=None, sep='\t')
