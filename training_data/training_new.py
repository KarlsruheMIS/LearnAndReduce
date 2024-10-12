import csv
import torch
import random
import sys
import math

import torch.nn.functional as F

from torch_geometric.data import Data

dataset = []
path = 'csv/'

reduction = int(sys.argv[1])
features = 8

train = val = test = 0
train_p = val_p = test_p = 0

with open('instances.txt', newline='') as instancefile:
    for instance in instancefile.read().splitlines():
        edge_list = []
        field_names = []
        fields = []

        original_graph = instance + '_original_graph.csv'
        original_reduction_data = instance + '_original_reduction_data.csv'

        with open(path + original_graph, newline="") as csvfile:
            csvfile.readline()
            reader = csv.reader(csvfile, delimiter=';', quoting = csv.QUOTE_NONNUMERIC)
            edge_list = list(reader)

        with open(path + original_reduction_data, newline="") as csvfile:
            field_names = csvfile.readline().split(';')
            reader = csv.reader(csvfile, delimiter=';', quoting = csv.QUOTE_NONNUMERIC)
            fields = list(reader)

        edge_index = torch.tensor(edge_list, dtype=torch.long)

        weight = [row[1] for row in fields]
        degree = [0.0 for _ in range(len(fields))]
        for u, v in edge_list:
            degree[int(u)] += 1.0

        neighborhood_degree = [0.0 for _ in range(len(fields))]
        neighborhood_weight = [0.0 for _ in range(len(fields))]
        
        neighborhood_min_weight = [math.inf if degree[i] > 0 else 0.0 for i in range(len(fields))]
        neighborhood_max_weight = [0.0 for _ in range(len(fields))]
        neighborhood_min_degree = [math.inf if degree[i] > 0 else 0.0 for i in range(len(fields))]
        neighborhood_max_degree = [0.0 for _ in range(len(fields))]

        max_degree = 0.0
        max_weight = 0.0
        
        for u, v in edge_list:
            neighborhood_degree[int(u)] += degree[int(v)]
            neighborhood_weight[int(u)] += weight[int(v)]

            neighborhood_min_weight[int(u)] = min(neighborhood_min_weight[int(u)], weight[int(v)])
            neighborhood_max_weight[int(u)] = max(neighborhood_max_weight[int(u)], weight[int(v)])
            neighborhood_min_degree[int(u)] = min(neighborhood_min_degree[int(u)], degree[int(v)])
            neighborhood_max_degree[int(u)] = max(neighborhood_max_degree[int(u)], degree[int(v)])

            max_degree = max(max_degree, degree[int(u)])
            max_weight = max(max_weight, weight[int(u)])

        for i in range(len(fields)):
            if degree[i] == 0.0:
                continue

            neighborhood_weight[i] = weight[i] - neighborhood_weight[i]
            neighborhood_degree[i] = (neighborhood_degree[i] / degree[i])

        x0 = torch.tensor(weight, dtype=torch.float).unsqueeze(1)
        x1 = torch.tensor(neighborhood_weight, dtype=torch.float).unsqueeze(1)
        x2 = torch.tensor(neighborhood_min_weight, dtype=torch.float).unsqueeze(1)
        x3 = torch.tensor(neighborhood_max_weight, dtype=torch.float).unsqueeze(1)
        x4 = torch.tensor(degree, dtype=torch.float).unsqueeze(1)
        x5 = torch.tensor(neighborhood_degree, dtype=torch.float).unsqueeze(1)
        x6 = torch.tensor(neighborhood_min_degree, dtype=torch.float).unsqueeze(1)
        x7 = torch.tensor(neighborhood_max_degree, dtype=torch.float).unsqueeze(1)

        x = torch.cat([x0, x1, x2, x3, x4, x5, x6, x7], 1)

        
        y = torch.tensor([[row[reduction]] for row in fields], dtype=torch.float)

        data_split = [10 if int(row[reduction]) == 2 else random.randint(0, 9) for row in fields]

        train_mask = torch.tensor([p < 6 for p in data_split], dtype=torch.bool)
        val_mask = torch.tensor([p >= 6 and p < 8 for p in data_split], dtype=torch.bool)
        test_mask = torch.tensor([p >= 8 and p < 10 for p in data_split], dtype=torch.bool)

        train += train_mask.sum().item()
        train_p += y[train_mask].sum().item()
        val += val_mask.sum().item()
        val_p += y[val_mask].sum().item()
        test += test_mask.sum().item()
        test_p += y[test_mask].sum().item()

        data = Data(x=x, y=y, edge_index=edge_index.t().contiguous(),
                    train_mask=train_mask, val_mask=val_mask, test_mask=test_mask)

        data.validate(raise_on_error=True)

        dataset.append(data)

print("Using:", field_names[reduction])
print("Size of dataset:", len(dataset))
print(train, int(train_p), val, int(val_p), test, int(test_p))

from torch_geometric.nn import MessagePassing, GCNConv, SAGEConv, GENConv, GINEConv, TransformerConv, PNAConv

class LRConv(MessagePassing):
    def __init__(self, in_channels, out_channels):
        super().__init__(flow='target_to_source', aggr='mean')
        self.lin = torch.nn.Linear(2 * in_channels, out_channels)
        self.reset_parameters()

    def reset_parameters(self):
        self.lin.reset_parameters()

    def forward(self, x, edge_index, edge_attr):
        out = self.propagate(edge_index, x=x)
        return out

    def message(self, x_i, x_j):
        x = torch.cat([x_i, x_j], dim=-1)
        x = self.lin(x)
        return x
        
class LR_GCN(torch.nn.Module):
    def __init__(self, in_channels, hidden_channels, out_channels):
        super().__init__()
        self.conv1 = LRConv(in_channels, hidden_channels)
        self.conv2 = LRConv(hidden_channels, hidden_channels)
        self.lin1 = torch.nn.Linear(2 * hidden_channels + in_channels, hidden_channels)
        self.lin2 = torch.nn.Linear(hidden_channels, out_channels)

    def forward(self, x, edge_index, edge_weight=None):
        x1 = F.relu(self.conv1(x, edge_index, edge_weight))
        x1 = F.dropout(x1, p=0.2, training=self.training)
        x2 = F.relu(self.conv2(x1, edge_index, edge_weight))
        x2 = F.dropout(x2, p=0.2, training=self.training)
        x = torch.cat([x, x1, x2], dim=-1)
        x = F.relu(self.lin1(x))
        x = F.sigmoid(self.lin2(x))
        return x

class GCN(torch.nn.Module):
    def __init__(self, in_channels, hidden_channels, out_channels):
        super().__init__()
        self.conv1 = GCNConv(in_channels, hidden_channels)
        self.conv2 = GCNConv(hidden_channels, hidden_channels)
        self.lin1 = torch.nn.Linear(2 * hidden_channels + in_channels, hidden_channels)
        self.lin2 = torch.nn.Linear(hidden_channels, out_channels)

    def forward(self, x, edge_index, edge_weight=None):
        x1 = F.relu(self.conv1(x, edge_index))
        x1 = F.dropout(x1, p=0.2, training=self.training)
        x2 = F.relu(self.conv2(x1, edge_index))
        x2 = F.dropout(x2, p=0.2, training=self.training)
        x = torch.cat([x, x1, x2], dim=-1)
        x = F.relu(self.lin1(x))
        x = F.sigmoid(self.lin2(x))
        return x

class SAGE(torch.nn.Module):
    def __init__(self, in_channels, hidden_channels, out_channels):
        super().__init__()
        self.conv1 = SAGEConv(in_channels, hidden_channels)
        self.conv2 = SAGEConv(hidden_channels, hidden_channels)
        self.lin1 = torch.nn.Linear(2 * hidden_channels + in_channels, hidden_channels)
        self.lin2 = torch.nn.Linear(hidden_channels, out_channels)

    def forward(self, x, edge_index, edge_weight=None):
        x1 = F.relu(self.conv1(x, edge_index))
        x1 = F.dropout(x1, p=0.2, training=self.training)
        x2 = F.relu(self.conv2(x1, edge_index))
        x2 = F.dropout(x2, p=0.2, training=self.training)
        x = torch.cat([x, x1, x2], dim=-1)
        x = F.relu(self.lin1(x))
        x = F.sigmoid(self.lin2(x))
        return x

def train():
    model.train()

    total_loss = total_examples = 0
    for data in dataset[0:1]:
        data = data.to(device)
        optimizer.zero_grad()

        out = model(data.x, data.edge_index)
        loss = F.binary_cross_entropy(out[data.train_mask], data.y[data.train_mask]) # pos_weight=torch.tensor([10.0])

        loss.backward()
        optimizer.step()
        total_loss += loss.item() * data.num_nodes
        total_examples += data.num_nodes
        
    return total_loss / total_examples

from torcheval.metrics.functional import binary_confusion_matrix

@torch.no_grad()
def test():
    model.eval()

    accs = torch.tensor([0 for _ in range(12)])
    for data in dataset:
        out = model(data.x.to(device), data.edge_index.to(device))
        pred = out.round()
        target = data.y.to(device)

        train = binary_confusion_matrix(pred[data.train_mask].flatten(), target[data.train_mask].flatten().long())
        val = binary_confusion_matrix(pred[data.val_mask].flatten(), target[data.val_mask].flatten().long())
        test = binary_confusion_matrix(pred[data.test_mask].flatten(), target[data.test_mask].flatten().long())

        accs += torch.cat([train.flatten(), val.flatten(), test.flatten()])
    
    return accs

def store_model(model, path):
    with open(path, "w") as f:
        f.write(str(torch.nn.utils.parameters_to_vector(model.parameters()).size()[0])+'\n')
        for p in torch.nn.utils.parameters_to_vector(model.parameters()):
            f.write(f'{p.item():.9f}'+'\n')

def train_model(model, it):
    print('epoch;tn_t;fp_t;fn_t;tp_t;tn_v;fp_v;fn_v;tp_v;tn_e;fp_e;fn_e;tp_e')
    for epoch in range(1, it + 1):
        loss = train()
        accs = test()
        print(f'{epoch:d};{loss:.4f}', end='')
        for p in range(12):
            print(f';{accs[p].item():d}', end='')
        print()

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

print('Training GCN')
model = GCN(features, 16, 1).to(device)
optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
train_model(model, 400)

print('Training SAGE')
model = SAGE(features, 16, 1).to(device)
optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
train_model(model, 400)

print('Training LR_CONV')
model = LR_GCN(features, 16, 1).to(device)
optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
train_model(model, 400)

store_model(model,  field_names[reduction] + '.lr_gcn')