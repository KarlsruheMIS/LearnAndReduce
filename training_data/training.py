import torch
from torch_geometric.data import Data

import csv
import os
import sys

reduction = sys.argv[1]
#reduction = 'domination'
print('Using', reduction)

dataset = []

def parse_data(filename):
    (path, ext) = os.path.splitext(filename)
    with open(path + '_meta' + ext) as f:
        reader = csv.DictReader(f, delimiter=';')
        if reduction not in reader.fieldnames:
            return
        y = torch.tensor([[int(r[reduction])] for r in reader], dtype=torch.float)

    with open(path + '_meta' + ext) as f:
        reader = csv.DictReader(f, delimiter=';')
        x = torch.tensor([[float(r['d']), float(r['w']), float(r['nw']), float(r['l'])] for r in reader], dtype=torch.float)
    
    with open(filename) as f:
        reader = csv.DictReader(f, delimiter=';')
        edge_index = torch.tensor([[int(r['source']), int(r['target'])] for r in reader], dtype=torch.long)

    with open(filename) as f:
        reader = csv.DictReader(f, delimiter=';')
        edge_attr = torch.tensor([[float(r['d']), float(r['w']), float(r['dr']), 
                                   float(r['wr']), float(r['e0']), float(r['e1']),
                                   float(r['e2']), float(r['e4'])] for r in reader], dtype=torch.float)

    if x.size()[0] > 100:
        dataset.append(Data(x=x, y=y, edge_index=edge_index.t().contiguous(), edge_attr=edge_attr)) #, edge_attr=edge_attr

path = sys.argv[2]
#path = 'csv/'

directory = os.fsencode(path)
for file in os.listdir(directory):
    filename = os.fsdecode(file)
    if not(filename.endswith('_meta.csv')):
        parse_data(path + filename)

weight_scale = float(sys.argv[3])

print('Found', len(dataset))

def model_fit(model, epoch):
    optimizer = torch.optim.Adam(model.parameters(), lr=0.01)
    loss = torch.nn.BCEWithLogitsLoss(pos_weight=torch.tensor([weight_scale]))
    scheduler = torch.optim.lr_scheduler.ExponentialLR(optimizer, gamma=0.99)
    
    N = len(dataset)
    M = int(N * 0.25)
    N -= M
    
    print(N, M)

    for e in range(epoch):
        running_loss = 0.0
        model.train()
        optimizer.zero_grad()
        for data in dataset[0:N]:
            out = model(data)
            l = loss(out, data.y)
            running_loss += l.item()
            l.backward()
        optimizer.step()
        if (e % 10) == 0:
            if scheduler.get_last_lr()[0] > 0.0005:
                scheduler.step()
            model.eval()
            tp, fp, tn, fn = 0, 0, 0, 0
            for data in dataset[N:N+M]:
                out = model(data)
                for i in range(data.y.size()[0]):
                    if data.y[i,0].item() > 0.5:
                        if out[i,0].item() > 0.0:
                            tp += 1
                        else:
                            fn += 1
                    else:
                        if out[i,0].item() > 0.0:
                            fp += 1
                        else:
                            tn += 1
            print(e, "%.5f" % scheduler.get_last_lr()[0], "%.5f" % (running_loss / N), tp, fp, tn, fn)

def store_model(model, path):
    with open(path, "w") as f:
        f.write(str(torch.nn.utils.parameters_to_vector(model.parameters()).size()[0])+'\n')
        for p in torch.nn.utils.parameters_to_vector(model.parameters()):
            f.write(str(p.item())+'\n')

from torch.nn import Linear, Parameter, Sequential, ReLU, Sigmoid
from torch_geometric.nn import MessagePassing

class LRConv(MessagePassing):
    def __init__(self, in_channels, out_channels, edge_channels):
        super().__init__(flow='target_to_source', aggr='max')  # "Add" aggregation (Step 5).
        self.seq = Sequential(
          Linear(in_channels * 2 + edge_channels, 16),
          ReLU(),
          Linear(16, out_channels),
          ReLU()
        )
        self.reset_parameters()

    def reset_parameters(self):
        for layer in self.seq:
            if hasattr(layer, 'reset_parameters'):
                layer.reset_parameters()

    def forward(self, x, edge_index, edge_attr):
        out = self.propagate(edge_index, x=x, edge_attr=edge_attr)
        #out = torch.cat([out, x], dim=-1)
        #out = out + self.bias

        return out

    def message(self, x_i, x_j, edge_attr):
        return self.seq(torch.cat([x_i, x_j, edge_attr], dim=-1))

from torch_geometric.nn import GCNConv, SAGEConv, GENConv, GINEConv, TransformerConv, PNAConv

class LR_GCN(torch.nn.Module):
    def __init__(self):
        super().__init__()
        self.conv1 = LRConv(4, 16, 8)
        self.conv2 = LRConv(16, 16, 8)
        self.lin1 = Linear(16, 16)
        self.a1 = ReLU()
        self.lin2 = Linear(16, 1)

    def forward(self, data):
        x, edge_index, edge_attr = data.x, data.edge_index, data.edge_attr
        x = self.conv1(x, edge_index, edge_attr)
        x = self.conv2(x, edge_index, edge_attr)
        x = self.lin1(x)
        x = self.a1(x)
        x = self.lin2(x)
        return x

class GCN(torch.nn.Module):
    def __init__(self):
        super().__init__()
        self.conv1 = GCNConv(4, 16)
        self.a1 = ReLU()
        self.conv2 = GCNConv(16, 1)

    def forward(self, data):
        x, edge_index, edge_attr = data.x, data.edge_index, data.edge_attr
        x = self.conv1(x, edge_index)
        x = self.a1(x)
        x = self.conv2(x, edge_index)
        return x

class SAGE(torch.nn.Module):
    def __init__(self):
        super().__init__()
        self.conv1 = SAGEConv(4, 16)
        self.a1 = ReLU()
        self.conv2 = SAGEConv(16, 1)

    def forward(self, data):
        x, edge_index, edge_attr = data.x, data.edge_index, data.edge_attr
        x = self.conv1(x, edge_index)
        x = self.a1(x)
        x = self.conv2(x, edge_index)
        return x
    

it = 2000

lr = LR_GCN()
print(torch.nn.utils.parameters_to_vector(lr.parameters()).size())
print('Learn and Reduce')
model_fit(lr, it)
store_model(lr, reduction+'_lr'+'.gnn')

# gcn = GCN()
# print(torch.nn.utils.parameters_to_vector(gcn.parameters()).size())
# print('GCN')
# model_fit(gcn, it)
# store_model(gcn, reduction+'_gcn'+'.gnn')

# sage = GCN()
# print(torch.nn.utils.parameters_to_vector(sage.parameters()).size())
# print('SAGE')
# model_fit(sage, it)
# store_model(gcn, reduction+'_sage'+'.gnn')