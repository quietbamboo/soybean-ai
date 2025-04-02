import sys
import pandas as pd

# 检查是否提供了足够的参数
# arg1-3为3个性状的assoc文件, arg4为输出文件的名字, arg5为性状名
# 例 python assoc2CMplot.py oil15.assoc.txt oil16.assoc.txt oil17.assoc.txt train.QC.oil.txt oil
if len(sys.argv) != 6:
    print("Usage: python script.py <arg1> <arg2> <arg3> <arg4> <arg5>")
    sys.exit(1)

# 获取传入的参数
arg1 = sys.argv[1]
arg2 = sys.argv[2]
arg3 = sys.argv[3]
arg4 = sys.argv[4]
arg5 = sys.argv[5]

trait15 = pd.read_csv(arg1, sep='\t')
trait15 = trait15[['rs', 'chr', 'ps', 'p_wald']]
trait15.rename(columns={'p_wald': arg5+'15'}, inplace=True)

trait16 = pd.read_csv(arg2, sep='\t')
trait16 = trait16[['rs', 'chr', 'ps', 'p_wald']]
trait16.rename(columns={'p_wald': arg5+'16'}, inplace=True)

trait17 = pd.read_csv(arg3, sep='\t')
trait17 = trait17[['rs', 'chr', 'ps', 'p_wald']]
trait17.rename(columns={'p_wald': arg5+'17'}, inplace=True)

trait = pd.merge(trait15, trait16, on=['rs', 'chr', 'ps'], how='outer')
trait = pd.merge(trait, trait17, on=['rs', 'chr', 'ps'], how='outer')
# trait = trait.apply(lambda x: x.fillna(x.max()), axis=0)      # 用每一列的最大值填充
trait = trait.fillna(1.0)

trait.to_csv(arg4, index=None, sep='\t')