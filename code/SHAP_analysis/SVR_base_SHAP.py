import argparse
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from sklearn.svm import SVR
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_absolute_error, r2_score
from scipy.stats import pearsonr

import shap

# 增加解析命令行参数
# python SVR_base_SHAP.py --fold 4 --trait Oil --output_file SVR.base.LD.best_fold.txt
parser = argparse.ArgumentParser(description="Train and evaluate different p_wald with SVR model and base type")
parser.add_argument('--fold', type=int, required=True, help="Best fold number")
parser.add_argument('--trait', type=str, required=True, help="Trait to predict (Oil/Protein/WSPC)")
parser.add_argument('--output_file', type=str, required=True, help="Name of output results file")
parser.add_argument('--base_path', type=str, default=r"../../Data_Train_Test", 
                    help="Base path for SNP data (default: ../../Data_Train_Test)")
parser.add_argument('--output_file_path', type=str, default=r"../../results", 
                    help="Base path for SNP data (default: ../../results)")

args = parser.parse_args()
fold = args.fold
trait = args.trait
output_file = args.output_file
base_path = args.base_path
output_file_path = args.output_file_path


def ped2onehot(ped_file_path, Line):
    base_to_index = {'A':0, 'T':1, 'C':2, 'G':3}
    SNP = []
    with open(ped_file_path, 'r') as ped:
        for line in ped:
            encodes = []
            snp = line.split()
            if snp[0] not in Line.values:
                continue
            for i in range(6, len(snp), 2):
                encode = [0,0,0,0]
                if snp[i] != '0':
                    encode[base_to_index[snp[i]]] += 1
                    encode[base_to_index[snp[i+1]]] += 1
                encodes.append(encode)
            SNP.append(encodes)
    return np.array(SNP)

def load_data(ped_file_path, phe_file_path, trait_column_num):
    phe = pd.read_csv(phe_file_path, header=0, sep='\s+')
    phe = phe.iloc[:,[0, trait_column_num]]
    phe = phe[phe.iloc[:, 1] != 0]
    Line = phe.iloc[:,0]
    SNP = ped2onehot(ped_file_path, Line)
    phe = phe.iloc[:, 1].values
    return SNP.reshape(SNP.shape[0], -1), phe.astype(np.float32)

# 定义评估指标的函数
def calculate_metrics(y_true, y_pred, dec = 4):
    mae = round(mean_absolute_error(y_true, y_pred), dec)
    pcc = round(pearsonr(y_true, y_pred)[0], dec)
    r2 = round(r2_score(y_true, y_pred), dec)
    return mae, pcc, r2

# 设置路径模板
train_test_sets = [ fr"{base_path}/LD/{trait}/Train/train{fold}.LD.{trait}3.p_2.ped",
                    fr"{base_path}/LD/{trait}/Test/test{fold}.LD.{trait}3.p_2.ped",
                    fr"{base_path}/traits/train/train_{fold}.txt",
                    fr"{base_path}/traits/test/test_{fold}.txt"]
trait_to_index = {"Protein":1, "Oil":4, "WSPC":7}
column_num = trait_to_index[trait]

# 检查保存结果的文件是否存在，如果不存在则创建
os.makedirs(output_file_path, exist_ok=True)
results_saved = os.path.join(output_file_path, output_file)

with open(results_saved, 'a') as file:
    file.write(f"Type: base\tML_model : SVR\tTrait : {trait}\tFold : {fold}\t\n")


# 训练最好一折的模型
train_ped, test_ped, train_phe, test_phe = train_test_sets

SNP_train, phe_train = load_data(train_ped, train_phe, column_num)
SNP_test, phe_test = load_data(test_ped, test_phe, column_num)

scaler = StandardScaler()
phe_train = scaler.fit_transform(phe_train.reshape(-1, 1)).reshape(-1)
phe_test = scaler.transform(phe_test.reshape(-1, 1)).reshape(-1)

model = SVR()
model.fit(SNP_train, phe_train)

phe_train_pred = model.predict(SNP_train)
phe_test_pred = model.predict(SNP_test)

results = []
train_mae, train_pcc, train_r2 = calculate_metrics(phe_train, phe_train_pred)
results.append(["Train", train_mae, train_pcc, train_r2])

test_mae,  test_pcc,  test_r2  = calculate_metrics(phe_test,  phe_test_pred)
results.append(["Test", test_mae, test_pcc, test_r2])

# SNP_all = np.vstack((SNP_train, SNP_test))
# phe_all = np.hstack((phe_train, phe_test))
SNP_all = SNP_test
phe_all = phe_test

phe_all_pred = model.predict(SNP_all)
all_mae, all_pcc, all_r2  = calculate_metrics(phe_all, phe_all_pred)
results.append(["All", all_mae, all_pcc, all_r2])

print(results)
# 将结果转换为DataFrame
results_df = pd.DataFrame(results, columns=['Data', 'MAE', 'PCC', 'R^2'])

# 将结果输出到文件，并确保列对齐
with open(results_saved, 'a') as file:
    file.write(results_df.to_string(index=False, col_space=12))
    file.write('\n\n')

# 使用SHAP计算特征的重要性
print("SHAP挖掘显著位点……")
# 初始化一个解释器（explainer）
batch_size = 30
explainer = shap.KernelExplainer(model.predict, shap.sample(SNP_all, batch_size))
# 计算 SHAP 值
shap_values = explainer.shap_values(SNP_all)
print(f"shap_values shape : {shap_values.shape}")


# 创建一个布尔数组，表示哪些元素是非零的
non_zero_mask = (SNP_all != 0)
# 沿着列的方向对布尔数组进行求和，得到每列非零元素的数量
non_zero_count_per_column = np.sum(non_zero_mask, axis=0)
# 找到具有非零值的列的索引
non_zero_columns = np.where(non_zero_count_per_column > 0)[0]
# 打印具有非零值的列的索引
print(f"具有非零值的列的索引：{len(non_zero_columns)}")

# filter 将样本中位点为0的shap值也置零
# 样本中位点为0表示该特征并未出现在样本中，0是一个占位符，所以需要删除
shap_values_filter = non_zero_mask * shap_values


# 保留非0列的SHAP值，再计算特征的平均SHAP值
# 先求绝对值再平均，SHAP值重要性
# 理解为什么先求绝对值：一个特征受到其他特征的影响
# 比如A特征和B特征在一起对模型起到正向影响，但和C特征在一起就会对模型起到负向影响
# 所有要对shap值先求平均，从而评估这个特征的重要性，负责会正负抵消
mean_abs_shap_values_filter = np.abs(shap_values_filter[:, non_zero_columns]).mean(axis=0)

# 求shap_values的平均值，查看每个特征整体起到正向影响还是负向影响
mean_shap_values_filter = np.mean(shap_values_filter[:, non_zero_columns], axis=0)

p_value_file_name = f"train{fold}.LD.p_2.{trait}.txt"
p_value_file = pd.read_csv(p_value_file_name, header=0, sep='\t')

# 与 p 值保存在一起
p_value_file['mean_abs_shap_filter_0'] = mean_abs_shap_values_filter[0::2]
p_value_file['mean_abs_shap_filter_1'] = mean_abs_shap_values_filter[1::2]

p_value_file['mean_shap_filter_0'] = mean_abs_shap_values_filter[0::2]
p_value_file['mean_shap_filter_1'] = mean_abs_shap_values_filter[1::2]

output_file_name = f"fold{fold}.{trait}.p_wald.shap_values.txt"
p_value_file.to_csv(output_file_name, index=None, sep='\t')


shap_values_png = f"shap_values.fold{fold}.{trait}.png"
shap.summary_plot(shap_values, SNP_all, max_display=40, show=False)
plt.savefig(shap_values_png)
plt.clf()

shap_values_filter_png = f"shap_values_filter.fold{fold}.{trait}.png"
shap.summary_plot(shap_values_filter, SNP_all, max_display=40, show=False)
plt.savefig(shap_values_filter_png)
plt.clf()
