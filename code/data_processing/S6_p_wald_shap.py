import sys
import pandas as pd

# input_dataset         输入数据集
# p_value_threshold     p值
# output_3_years_file   保存3年交叉结果的文件名

# 检查命令行参数数量是否正确
if len(sys.argv) != 4:
    print("Usage: python script.py <input_dataset> <p_value_threshold> <output_file>")
    sys.exit(1)

input_dataset = sys.argv[1]
p_value_thrd  = sys.argv[2]
output_file   = sys.argv[3]

try:
    # 将阈值参数转换为浮点型
    thrd = float(p_value_thrd)
except ValueError:
    print("Error: <p_value_threshold> must be a numeric value.")
    sys.exit(1)

assoc = pd.read_csv(input_dataset, sep='\t')
assoc_3 = assoc[(assoc.iloc[:,3] < thrd) & (assoc.iloc[:,4] < thrd) & (assoc.iloc[:,5] < thrd)]

assoc_3.to_csv(output_file, index=None, sep='\t')
