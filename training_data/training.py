import torch
from torch_geometric.data import Data

import csv
import os
import sys

reduction = sys.argv[1]
#reduction = 'domination'
print('Using', reduction)
scale = 100.0

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
        x = torch.tensor([[float(r['w']) / scale, float(r['d']) / scale, float(r['c'])] for r in reader], dtype=torch.float)
    
    with open(filename) as f:
        reader = csv.DictReader(f, delimiter=';')
        edge_index = torch.tensor([[int(r['source']), int(r['target'])] for r in reader], dtype=torch.long)

    with open(filename) as f:
        reader = csv.DictReader(f, delimiter=';')
        edge_attr = torch.tensor([[float(r['uc']) / scale, float(r['vc']) / scale, float(r['ic']) / scale, 
                                   float(r['uw']) / scale, float(r['vw']) / scale, float(r['iw']) / scale] for r in reader], dtype=torch.float)

    dataset.append(Data(x=x, y=y, edge_index=edge_index.t().contiguous(), edge_attr=edge_attr))

path = sys.argv[2]

directory = os.fsencode(path)
for file in os.listdir(directory):
    filename = os.fsdecode(file)
    if not(filename.endswith('_meta.csv')):
        parse_data(path + filename)

print('Found', len(dataset))

from torch.nn import Linear, Parameter, Sequential, ReLU, Sigmoid
from torch_geometric.nn import MessagePassing

class LRConv(MessagePassing):
    def __init__(self, in_channels, out_channels, edge_channels):
        super().__init__(aggr='max')  # "Add" aggregation (Step 5).
        self.seq = Sequential(
          Linear(in_channels * 2 + edge_channels, out_channels * 2),
          ReLU(),
          Linear(out_channels * 2, out_channels),
          ReLU()
        )
        self.reset_parameters()

    def reset_parameters(self):
        for layer in self.seq:
            if hasattr(layer, 'reset_parameters'):
                layer.reset_parameters()

    def forward(self, x, edge_index, edge_attr):
        out = self.propagate(edge_index, x=x, edge_attr=edge_attr)
        #out = out + self.bias

        return out

    def message(self, x_i, x_j, edge_attr):
        return self.seq(torch.cat([x_i, x_j, edge_attr], dim=-1))
    

def model_fit(model, epoch):
    optimizer = torch.optim.Adam(model.parameters())
    loss = torch.nn.MSELoss()
    
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
            optimizer.zero_grad()
        if e < 50 or (e % 10) == 0:
            model.eval()
            tp, fp, tn, fn = 0, 0, 0, 0
            for data in dataset[N:N+M]:
                out = model(data)
                for i in range(data.y.size()[0]):
                    if data.y[i,0].item() > 0.5:
                        if out[i,0].item() > 0.5:
                            tp += 1
                        else:
                            fn += 1
                    else:
                        if out[i,0].item() > 0.5:
                            fp += 1
                        else:
                            tn += 1
            print(e, running_loss / N, tp, fp, tn, fn)

from torch_geometric.nn import GCNConv, SAGEConv, GENConv, GINEConv, TransformerConv, PNAConv

class LR_GCN(torch.nn.Module):
    def __init__(self):
        super().__init__()
        self.conv1 = LRConv(3, 32, 6)
        self.lin1 = Linear(32, 16)
        self.a1 = ReLU()
        self.lin2 = Linear(16, 1)
        self.a2 = Sigmoid()

    def forward(self, data):
        x, edge_index, edge_attr = data.x, data.edge_index, data.edge_attr
        x = self.conv1(x, edge_index, edge_attr)
        x = self.lin1(x)
        x = self.a1(x)
        x = self.lin2(x)
        return self.a2(x)

class GCN(torch.nn.Module):
    def __init__(self):
        super().__init__()
        self.conv1 = GCNConv(3, 16)
        self.a1 = ReLU()
        self.conv2 = GCNConv(16, 1)
        self.a2 = Sigmoid()

    def forward(self, data):
        x, edge_index, edge_attr = data.x, data.edge_index, data.edge_attr
        x = self.conv1(x, edge_index)
        x = self.a1(x)
        x = self.conv2(x, edge_index)
        return self.a2(x)

class SAGE(torch.nn.Module):
    def __init__(self):
        super().__init__()
        self.conv1 = SAGEConv(3, 16)
        self.a1 = ReLU()
        self.conv2 = SAGEConv(16, 1)
        self.a2 = Sigmoid()

    def forward(self, data):
        x, edge_index, edge_attr = data.x, data.edge_index, data.edge_attr
        x = self.conv1(x, edge_index)
        x = self.a1(x)
        x = self.conv2(x, edge_index)
        return self.a2(x)

class GEN(torch.nn.Module):
    def __init__(self):
        super().__init__()
        self.conv1 = GENConv(3, 16, edge_dim=6)
        self.a1 = ReLU()
        self.conv2 = GENConv(16, 1, edge_dim=6)
        self.a2 = Sigmoid()

    def forward(self, data):
        x, edge_index, edge_attr = data.x, data.edge_index, data.edge_attr
        x = self.conv1(x, edge_index, edge_attr)
        x = self.a1(x)
        x = self.conv2(x, edge_index, edge_attr)
        return self.a2(x)

class GINE(torch.nn.Module):
    def __init__(self):
        super().__init__()
        self.seq1 = Sequential(Linear(3, 16), ReLU())
        self.conv1 = GINEConv(self.seq1, edge_dim=6)
        self.a1 = ReLU()
        self.seq2 = Sequential(Linear(16, 1), ReLU())
        self.conv2 = GINEConv(self.seq2, edge_dim=6)
        self.a2 = Sigmoid()

    def forward(self, data):
        x, edge_index, edge_attr = data.x, data.edge_index, data.edge_attr
        x = self.conv1(x, edge_index, edge_attr)
        x = self.a1(x)
        x = self.conv2(x, edge_index, edge_attr)
        return self.a2(x)

class Transformer(torch.nn.Module):
    def __init__(self):
        super().__init__()
        self.conv1 = TransformerConv(3, 16, edge_dim=6)
        self.a1 = ReLU()
        self.conv2 = TransformerConv(16, 1, edge_dim=6)
        self.a2 = Sigmoid()

    def forward(self, data):
        x, edge_index, edge_attr = data.x, data.edge_index, data.edge_attr
        x = self.conv1(x, edge_index, edge_attr)
        x = self.a1(x)
        x = self.conv2(x, edge_index, edge_attr)
        return self.a2(x)
    

it = 1000

lr = LR_GCN()
print('Learn and Reduce')
model_fit(lr, it)

gcn = GCN()
print('GCN')
model_fit(gcn, it)

sage = SAGE()
print('Sage')
model_fit(sage, it)

gen = GEN()
print('GENConv')
model_fit(gen, it)

gine = GINE()
print('GINEConv')
model_fit(gine, it)

transformer = Transformer()
print('TransformerConv')
model_fit(transformer, it)