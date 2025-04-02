from PredModel.dualCNN_torch_version import DualCNN
from PredModel.DNNGP import DNNGP

import argparse
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_absolute_error, r2_score
from scipy.stats import pearsonr

# 增加解析命令行参数
parser = argparse.ArgumentParser(description="Train and evaluate model with SNP data")
parser.add_argument('--DL_model', type=str, required=True, help="DL model (DualCNN/DNNGP)")
parser.add_argument('--SNP_type', type=str, required=True, help="SNP data type (QC/LD)")
parser.add_argument('--trait', type=str, required=True, help="Trait to predict (Oil/Protein/WSPC)")
parser.add_argument('--year', type=int, required=True, help="Year of trait (e.g. 3, 15)")
parser.add_argument('--output_file', type=str, required=True, help="Name of output results file")
parser.add_argument('--base_path', type=str, default=r"../../Data_Train_Test", 
                    help="Base path for SNP data (default: ../../Data_Train_Test)")
parser.add_argument('--output_file_path', type=str, default=r"../../results/DL_results", 
                    help="Base path for SNP data (default: ../../results/DL_results)")

args = parser.parse_args()
DL_model = args.DL_model
SNP_type = args.SNP_type
trait = args.trait
year = args.year
output_file = args.output_file
base_path = args.base_path
output_file_path = args.output_file_path

'''
pytorch 设置gpu设备
torch.cuda.is_available() 查看gpu是否可用
'''
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(device)

print(f"Type: base\t\tML_model: {DL_model}\t\tSNP_type: {SNP_type}\t\ttrait: {trait}\t\tyear:{year}")

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

def load_data(ped_file_path, phe_file_path, trait_column_num):
    phe = pd.read_csv(phe_file_path, header=0, sep='\s+')
    phe = phe.iloc[:, [0, trait_column_num]]
    phe = phe[phe.iloc[:, 1] != 0]
    Line = phe.iloc[:, 0]
    SNP = raw2onehot(ped_file_path, Line)
    phe = phe.iloc[:, 1].values
    return SNP, phe.astype(np.float32)

# 设置路径模板
train_test_sets = [(fr"{base_path}/{SNP_type}/{trait}/Train/train{i+1}.{SNP_type}.{trait}{year}.p_2.raw",
                    fr"{base_path}/{SNP_type}/{trait}/Test/test{i+1}.{SNP_type}.{trait}{year}.p_2.raw",
                    fr"{base_path}/traits/train/train_{i+1}.txt",
                    fr"{base_path}/traits/test/test_{i+1}.txt") for i in range(5)]
trait_to_index = {"Protein":1, "Oil":4, "WSPC":7}
column_num = trait_to_index[trait]

# 检查保存结果的文件是否存在，如果不存在则创建
os.makedirs(output_file_path, exist_ok=True)
results_saved = os.path.join(output_file_path, output_file)

with open(results_saved, 'a') as file:
    file.write(f"DL_model : {DL_model}\tSNP_type : {SNP_type}\tTrait : {trait}\tYear : {year} \n")

def train_model(model, train_loader, val_loader, criterion, optimizer, epochs=1000):
    best_val_mae = float('inf')
    best_model_state_dict = None
    no_improvement_counter = 0
    early_stop_patience = 20

    train_mae_recodes = []
    val_mae_recodes = []

    for epoch in range(epochs):
        model.train()
        train_loss = 0.0
        train_mae = 0.0
        for inputs, targets in train_loader:
            inputs, targets = inputs.to(device), targets.to(device)
            optimizer.zero_grad()
            outputs = model(inputs)
            loss = criterion(outputs.squeeze(), targets)
            loss.backward()
            optimizer.step()
            
            train_loss += loss.item() * inputs.size(0)
            train_mae += torch.mean(torch.abs(outputs.squeeze() - targets)).item()
        train_loss /= len(train_loader.dataset)
        train_mae /= len(train_loader)
        train_mae_recodes.append(train_mae)

        model.eval()
        val_mae = 0.0
        with torch.no_grad():
            for inputs, targets in val_loader:
                inputs, targets = inputs.to(device), targets.to(device)
                outputs = model(inputs)
                val_mae += torch.mean(torch.abs(outputs.squeeze() - targets)).item()
            val_mae /= len(val_loader)
            val_mae_recodes.append(val_mae)
    
        print(f"Epoch {epoch + 1}:\tTrain Loss:{train_loss:.4f}\tTrain MAE:{train_mae:.4f}\tVal MAE:{val_mae:.4f}")

        if val_mae < best_val_mae:
            best_val_mae = val_mae
            best_model_state_dict = model.state_dict()
            no_improvement_counter = 0
        else:
            no_improvement_counter += 1
            if no_improvement_counter >= early_stop_patience:
                print("Validation MAE no longer decreasing. Early stopping.")
                break

    return best_model_state_dict, train_mae_recodes, val_mae_recodes

def calculate_metrics(y_true, y_pred):
    mae = mean_absolute_error(y_true, y_pred)
    pcc, _ = pearsonr(y_true, y_pred)
    r2 = r2_score(y_true, y_pred)
    return mae, pcc, r2

def calculate_test(model, test_loader):
    model.eval()
    predictions = []
    true_labels = []
    with torch.no_grad():
        for inputs, targets in test_loader:
            inputs = inputs.to(device)
            outputs = model(inputs).squeeze().cpu()

            if outputs.dim() == 0:
                predictions.append(outputs.item())
            else:
                predictions.extend(outputs.tolist())
            true_labels.extend(targets.tolist())

        test_mae, test_pcc, test_r2 = calculate_metrics(predictions, true_labels)
    return test_mae, test_pcc, test_r2

batch_size = 10
learning_rate = 0.0005

results = {"fold":[], "train_mae":[], "test_mae":[], "train_pcc":[], "test_pcc":[], "train_r2":[], "test_r2":[]}
fold = 0
for train_ped, test_ped, train_phe, test_phe in train_test_sets:
    fold += 1
    print(f"Fold {fold}")
    SNP_train, phe_train = load_data(train_ped, train_phe, column_num)
    SNP_test, phe_test = load_data(test_ped, test_phe, column_num)

    label_scaler = StandardScaler()
    label_scaler.fit(phe_train.reshape(-1, 1))
    phe_train = label_scaler.transform(phe_train.reshape(-1, 1)).flatten()
    phe_test = label_scaler.transform(phe_test.reshape(-1, 1)).flatten()

    SNP_train = torch.tensor(SNP_train, dtype=torch.float32).permute(0, 2, 1)
    phe_train = torch.tensor(phe_train, dtype=torch.float32)
    SNP_test = torch.tensor(SNP_test, dtype=torch.float32).permute(0, 2, 1)
    phe_test = torch.tensor(phe_test, dtype=torch.float32)

    train_dataset = TensorDataset(SNP_train, phe_train)
    test_dataset = TensorDataset(SNP_test, phe_test)

    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
    test_loader = DataLoader(test_dataset, batch_size=batch_size, shuffle=False)

    if DL_model == "DualCNN":
        model = DualCNN(in_dim=SNP_train.shape[2])
    elif DL_model == "DNNGP":
        model = DNNGP(in_dim=SNP_train.shape[2])
    else:
        print("DL_model error!")
        sys.exit(1)
    model.to(device)
    criterion = nn.MSELoss()
    l2_regularization = 0.01
    optimizer = optim.Adam(model.parameters(), lr=learning_rate, weight_decay=l2_regularization)

    best_model_state_dict, train_mae_recodes, test_mae_recodes = train_model(model, train_loader, test_loader, criterion, optimizer, epochs=100)
    model_saved_file_name = os.path.join(output_file_path, f'Fold_{fold}_best_model.pth')
    torch.save(best_model_state_dict, model_saved_file_name)

    model.load_state_dict(best_model_state_dict)

    results["fold"].append(fold)
    train_mae, train_pcc, train_r2 = calculate_test(model, train_loader)
    results["train_mae"].append(train_mae)
    results["train_pcc"].append(train_pcc)
    results["train_r2"].append(train_r2)
    test_mae, test_pcc, test_r2 = calculate_test(model, test_loader)
    results["test_mae"].append(test_mae)
    results["test_pcc"].append(test_pcc)
    results["test_r2"].append(test_r2)

    plt.plot(test_mae_recodes, label=f'Fold {fold} Validation Loss')

plt.legend()
# plt.ylim(None, 1)
plt.xlabel('Epochs')
plt.ylabel('Loss MAE')
plt.title('Training and Validation Convergence Curves for 5 Folds')
img_saved_path = os.path.join(output_file_path, f'fold5_curve_{SNP_type}_{trait}_{year}.png')
plt.savefig(img_saved_path)

results = pd.DataFrame(results)
average_row = results.mean()
results = results.append(average_row, ignore_index=True)
results.iloc[-1, 0] = "avg"
print(results)

with open(results_saved, 'a') as file:
    file.write(results.to_string(index=False, col_space=12))
    file.write('\n\n')
