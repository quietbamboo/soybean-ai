import torch
import torch.nn as nn
import torch.nn.functional as F

class ISRU(nn.Module):
    def __init__(self, a=0.03):
        super(ISRU, self).__init__()
        self.a = a

    def forward(self, x):
        return x / torch.sqrt(1 + self.a * x**2)

class DualCNN(nn.Module):
    def __init__(self, in_dim, nb_classes=4):
        super(DualCNN, self).__init__()

        # First part
        self.conv1 = nn.Conv1d(nb_classes, 10, kernel_size=4, padding='same')
        self.conv2 = nn.Conv1d(10, 10, kernel_size=20, padding='same')
        self.dropout1 = nn.Dropout(0.75)

        # Second part
        self.shortcut = nn.Conv1d(nb_classes, 10, kernel_size=4, padding='same')

        # Third part
        self.conv3 = nn.Conv1d(10, 10, kernel_size=4, padding='same')
        self.dropout2 = nn.Dropout(0.75)
        self.flatten = nn.Flatten()
        self.dropout3 = nn.Dropout(0.75)
        self.out = nn.Linear(10*in_dim, 1)

    def forward(self, x):

        # First part
        out1 = F.relu(self.conv1(x))
        out1 = F.relu(self.conv2(out1))
        out1 = self.dropout1(out1)

        # Second part
        shortcut = self.shortcut(x)

        # Third part
        merged = out1 + shortcut

        out2 = F.relu(self.conv3(merged))
        out2 = self.dropout2(out2)
        out2 = self.flatten(out2)
        out2 = self.dropout3(out2)

        # Output layer
        out = self.out(out2)
        out = ISRU()(out)

        return out

# # Create an instance of the model
# input_length =  # provide the input length
# model = DualCNN(input_length)

# # Print the model architecture
# print(model)
