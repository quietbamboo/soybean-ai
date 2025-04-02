import argparse
import pandas as pd
import numpy as np

from sklearn.linear_model import Ridge
from sklearn.svm import SVR
from sklearn.ensemble import RandomForestRegressor
from xgboost import XGBRegressor
from sklearn.model_selection import GridSearchCV
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_absolute_error, r2_score
from scipy.stats import pearsonr
import os
import sys

# 增加解析命令行参数
parser = argparse.ArgumentParser(description="Train and evaluate model with SNP data")
parser.add_argument('--ML_model', type=str, required=True, help="ML model (SVR/RF/XGBoost)")
parser.add_argument('--SNP_type', type=str, required=True, help="SNP data type (QC/LD)")
parser.add_argument('--trait', type=str, required=True, help="Trait to predict (Oil/Protein/WSPC)")
parser.add_argument('--year', type=int, required=True, help="Year of trait (e.g. 3, 15)")
parser.add_argument('--output_file', type=str, required=True, help="Name of output results file")
parser.add_argument('--base_path', type=str, default=r"../../Data_Train_Test", 
                    help="Base path for SNP data (default: ../../Data_Train_Test)")
parser.add_argument('--output_file_path', type=str, default=r"../../results/ML_results", 
                    help="Base path for SNP data (default: ../../results/ML_results)")

args = parser.parse_args()
ML_model = args.ML_model
SNP_type = args.SNP_type
trait = args.trait
year = args.year
output_file = args.output_file
base_path = args.base_path
output_file_path = args.output_file_path

print(f"Type: gene\t\tML_model: {ML_model}\t\tSNP_type: {SNP_type}\t\ttrait: {trait}\t\tyear:{year}")

def raw2onehot(raw_file_path, Line):
    raw_to_index = {'0':0, '1':1, '2':2, 'NA':3}
    SNP = []
    with open(raw_file_path, 'r') as ped:
        for line in ped:
            encodes = []
            snp = line.split()
            if snp[0] not in Line.values:
                continue
            for i in range(6, len(snp)):
                encode = [0,0,0,0]
                encode[raw_to_index[snp[i]]] += 1
                encodes.append(encode)
            SNP.append(encodes)
    return np.array(SNP)

def load_data(raw_file_path, phe_file_path, trait_column_num):
    phe = pd.read_csv(phe_file_path, header=0, sep='\s+')
    phe = phe.iloc[:,[0, trait_column_num]]
    phe = phe[phe.iloc[:, 1] != 0]
    Line = phe.iloc[:,0]
    SNP = raw2onehot(raw_file_path, Line)
    phe = phe.iloc[:, 1].values
    return SNP.reshape(SNP.shape[0], -1), phe.astype(np.float32)

# 设置路径模板
train_test_sets = [(fr"{base_path}/{SNP_type}/{trait}/Train/train{i+1}.{SNP_type}.{trait}{year}.p_2.raw",
                    fr"{base_path}/{SNP_type}/{trait}/Test/test{i+1}.{SNP_type}.{trait}{year}.p_2.raw",
                    fr"{base_path}/traits/train/train_{i+1}.txt",
                    fr"{base_path}/traits/test/test_{i+1}.txt") for i in range(5)]
trait_to_index = {"Protein":1, "Oil":4, "WSPC":7}
column_num = trait_to_index[trait]
results = []

# 检查保存结果的文件是否存在，如果不存在则创建
os.makedirs(output_file_path, exist_ok=True)
results_saved = os.path.join(output_file_path, output_file)

with open(results_saved, 'a') as file:
    file.write(f"ML_model : {ML_model}\tSNP_type : {SNP_type}\tTrait : {trait}\tYear : {year} \n")

fold = 0
for train_ped, test_ped, train_phe, test_phe in train_test_sets:
    fold = fold + 1
    SNP_train, phe_train = load_data(train_ped, train_phe, column_num)
    SNP_test, phe_test = load_data(test_ped, test_phe, column_num)

    scaler = StandardScaler()
    phe_train = scaler.fit_transform(phe_train.reshape(-1, 1)).reshape(-1)
    phe_test = scaler.transform(phe_test.reshape(-1, 1)).reshape(-1)

    if ML_model == "SVR":
        model = SVR()
    elif ML_model == "SVR_linear":
        model = SVR(kernel='linear')
    elif ML_model == "RF":
        model = RandomForestRegressor()
    elif ML_model == "XGBoost":
        model = XGBRegressor()
    elif ML_model == "XGBoost_linear":
        model = XGBRegressor(booster='gblinear')
    elif ML_model == "Ridge":
        model = Ridge(alpha=1.0)
    else:
        print("ML_model error!")
        sys.exit(1)
    model.fit(SNP_train, phe_train)

    phe_train_pred = model.predict(SNP_train)
    phe_test_pred = model.predict(SNP_test)
    
    dec = 4
    train_mae = round(mean_absolute_error(phe_train, phe_train_pred), dec)
    train_pcc = round(pearsonr(phe_train, phe_train_pred)[0], dec)
    train_r2  = round(r2_score(phe_train, phe_train_pred), dec)

    test_mae = round(mean_absolute_error(phe_test, phe_test_pred), dec)
    test_pcc = round(pearsonr(phe_test, phe_test_pred)[0], dec)
    test_r2  = round(r2_score(phe_test, phe_test_pred), dec)
    
    results.append([fold, train_mae, test_mae, train_pcc, test_pcc, train_r2, test_r2])
  
# 将结果转换为DataFrame并计算前5行的平均值
results_df = pd.DataFrame(results, columns=['Fold', 'Train_MAE', 'Test_MAE', 'Train_PCC', 'Test_PCC', 'Train_R^2', 'Test_R^2'])
# 计算平均值（忽略'Fold'列）
mean_values = results_df.iloc[:, 1:].mean().round(dec)
# 添加一行表示平均值
mean_row = ['Avg'] + mean_values.tolist()
results_df.loc[len(results_df)] = mean_row

# 将结果输出到文件，并确保列对齐
with open(results_saved, 'a') as file:
    file.write(results_df.to_string(index=False, col_space=12))
    file.write('\n\n')
