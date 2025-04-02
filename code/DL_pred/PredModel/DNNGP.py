import torch
import torch.nn as nn
import math

class DNNGP(nn.Module):
    def __init__(self, in_dim, dp1_rate=0.5, dp2_rate=0.5, nb_classes=4):
        super(DNNGP, self).__init__()

        self.in_dim = in_dim
        self.in_dim = int(self.in_dim-3)
        self.in_dim = int(self.in_dim-3)
        self.in_dim = int(self.in_dim-3)

        # 第一层卷积
        self.conv1 = nn.Conv1d(nb_classes, 64, 4)
        self.act1 = nn.ReLU(inplace=True)
        self.dp1 = nn.Dropout(dp1_rate)
        self.norm = nn.BatchNorm1d(64)

        # 第二层卷积
        self.conv2 = nn.Conv1d(64, 64, 4)
        self.act2 = nn.ReLU(inplace=True)
        self.dp2 = nn.Dropout(dp2_rate)

        # 第三层卷积
        self.conv3 = nn.Conv1d(64, 64, 4)
        self.act3 = nn.ReLU(inplace=True)

        # 全连接层
        self.flatten = nn.Flatten()
        self.fc = nn.Linear(self.in_dim*64, 1)

    def forward(self, x):
        # 第一层卷积
        x = self.conv1(x)
        x = self.act1(x)
        x = self.dp1(x)
        x = self.norm(x)

        # 第二层卷积
        x = self.conv2(x)
        x = self.act2(x)
        x = self.dp2(x)

        # 第三层卷积
        x = self.conv3(x)
        x = self.act3(x)

        # 全连接层
        # x = x.reshape(x.size(0), -1)
        x = self.flatten(x)
        x = self.fc(x)
        
        return x
    
if __name__ == '__main__':
    model = DNNGP(100)
    model = model.cuda()

    from torchsummary import summary
    summary(model, (4, 100))