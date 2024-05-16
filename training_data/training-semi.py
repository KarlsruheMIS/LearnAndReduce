import torch
from torch_geometric.data import Data

import csv
import os
import sys

reduction = sys.argv[1]
print('Using', reduction)

dataset = []

def parse_data(filename):
    (path, ext) = os.path.splitext(filename)
    with open(path + '_meta' + ext) as f:
        reader = csv.DictReader(f, delimiter=';')
        y_in = torch.tensor([int(r['i']) for r in reader], dtype=torch.float)
    with open(path + '_meta' + ext) as f:
        reader = csv.DictReader(f, delimiter=';')
        y_out = torch.tensor([int(r['e']) for r in reader], dtype=torch.float)
        
    mask = torch.tensor([(y_in[i].item() + y_out[i].item() == 1.0) for i in range(y_in.size()[0])], dtype=torch.bool)
    y = torch.tensor([0.5 for i in range(y_in.size()[0])], dtype=torch.float)
    for i in range(y_in.size()[0]):
        if y_in[i].item() > 0.5 and y_out[i].item() < 0.5:
            y[i] = 1.0
        elif y_in[i].item() < 0.5 and y_out[i].item() > 0.5:
            y[i] = 0.0

    with open(path + '_meta' + ext) as f:
        reader = csv.DictReader(f, delimiter=';')
        x = torch.tensor([[float(r['d']), float(r['w']), float(r['nw']), float(r['l'])] for r in reader], dtype=torch.float)
    
    with open(filename) as f:
        reader = csv.DictReader(f, delimiter=';')
        edge_index = torch.tensor([[int(r['source']), int(r['target'])] for r in reader], dtype=torch.long)

    if x.size()[0] > 100:
        dataset.append(Data(x=x, y=y, edge_index=edge_index.t().contiguous(), mask=mask)) #, edge_attr=edge_attr

    for data in dataset:
        print(data)

path = sys.argv[2]
#path = 'csv/'

directory = os.fsencode(path)
for file in os.listdir(directory):
    filename = os.fsdecode(file)
    if not(filename.endswith('_meta.csv')):
        parse_data(path + filename)

print('Found', len(dataset))

for data in dataset:
    print(data.mask.sum(), data.y[data.mask].sum())

def model_fit(model, epoch):
    optimizer = torch.optim.Adam(model.parameters(), lr=0.01)
    loss = torch.nn.MSELoss()
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
            l = loss(out[data.mask], data.y[data.mask])
            running_loss += l.item()
            l.backward()
        optimizer.step()
        if (e % 10) == 0:
            if scheduler.get_last_lr()[0] > 0.0005:
                scheduler.step()
            model.eval()
            tp, fp, tn, fn = 0, 0, 0, 0
            new_p, new_n, total = 0, 0, 0
            for data in dataset[N:N+M]:
                out = model(data)
                for i in range(data.y.size()[0]):
                    if not(data.mask[i].item()):
                        if out[i].item() > 0.95:
                            new_p += 1
                        elif out[i].item() < 0.05:
                            new_n += 1
                        total += 1
                        continue
                    if data.y[i].item() > 0.5:
                        if out[i].item() > 0.5:
                            tp += 1
                        else:
                            fn += 1
                    elif data.y[i].item() < 0.5:
                        if out[i].item() > 0.5:
                            fp += 1
                        else:
                            tn += 1
            print(e, "%.5f" % scheduler.get_last_lr()[0], "%.5f" % (running_loss / N), tp, fp, tn, fn, new_p, new_n, total)

def store_model(model, path):
    with open(path, "w") as f:
        f.write(str(torch.nn.utils.parameters_to_vector(model.parameters()).size()[0])+'\n')
        for p in torch.nn.utils.parameters_to_vector(model.parameters()):
            f.write(str(p.item())+'\n')

from torch.nn import Linear, Parameter, Sequential, ReLU, Sigmoid
from torch_geometric.nn import MessagePassing

class LRConv(MessagePassing):
    def __init__(self, in_channels, out_channels):
        super().__init__(flow='target_to_source', aggr='max')  # "Add" aggregation (Step 5).
        self.seq = Sequential(
          Linear(in_channels * 2, 16),
          ReLU(),
          Linear(16, out_channels),
          ReLU()
        )
        self.reset_parameters()

    def reset_parameters(self):
        for layer in self.seq:
            if hasattr(layer, 'reset_parameters'):
                layer.reset_parameters()

    def forward(self, x, edge_index):
        out = self.propagate(edge_index, x=x)
        #out = torch.cat([out, x], dim=-1)
        #out = out + self.bias

        return out

    def message(self, x_i, x_j):
        return self.seq(torch.cat([x_i, x_j], dim=-1))

from torch_geometric.nn import GCNConv, SAGEConv, GENConv, GINEConv, TransformerConv, PNAConv

class LR_GCN(torch.nn.Module):
    def __init__(self):
        super().__init__()
        self.conv1 = LRConv(4, 16)
        self.conv2 = LRConv(16, 16)
        self.lin1 = Linear(16, 16)
        self.a1 = ReLU()
        self.lin2 = Linear(16, 1)
        self.a2 = Sigmoid()

    def forward(self, data):
        x, edge_index = data.x, data.edge_index
        x = self.conv1(x, edge_index)
        x = self.conv2(x, edge_index)
        x = self.lin1(x)
        x = self.a1(x)
        x = self.lin2(x)
        x = self.a2(x).squeeze()
        return x
    

it = 2000

lr = LR_GCN()
print(torch.nn.utils.parameters_to_vector(lr.parameters()).size())
print('Learn and Reduce')
model_fit(lr, it)
store_model(lr, reduction+'_lr'+'.gnn')
