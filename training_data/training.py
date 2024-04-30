import torch
from torch_geometric.data import Data

import csv
import os

reduction = 'domination'

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
        x = torch.tensor([[int(r['w']), int(r['d']), float(r['c'])] for r in reader], dtype=torch.float)
    
    with open(filename) as f:
        reader = csv.DictReader(f, delimiter=';')
        edge_index = torch.tensor([[int(r['source']), int(r['target'])] for r in reader], dtype=torch.long)

    with open(filename) as f:
        reader = csv.DictReader(f, delimiter=';')
        edge_attr = torch.tensor([[int(r['uc']), int(r['vc']), int(r['ic']), int(r['uw']), int(r['vw']), int(r['iw'])] for r in reader], dtype=torch.float)

    dataset.append(Data(x=x, y=y, edge_index=edge_index.t().contiguous(), edge_attr=edge_attr))

directory = os.fsencode('/scratch/problems/MWIS_learn_and_reduce/training_data/csv')
for file in os.listdir(directory):
    filename = os.fsdecode(file)
    if not(filename.endswith('_meta.csv')):
        parse_data('/scratch/problems/MWIS_learn_and_reduce/training_data/csv/' + filename)

from torch_geometric.nn import TransformerConv
import torch.nn.functional as F

class GCN(torch.nn.Module):
    def __init__(self):
        super().__init__()
        self.conv1 = TransformerConv(3, 16, edge_dim=6)
        self.conv2 = TransformerConv(16, 1, edge_dim=6)

    def forward(self, data):
        x, edge_index, edge_attr = data.x, data.edge_index, data.edge_attr

        x = self.conv1(x, edge_index, edge_attr)
        x = F.relu(x)
        x = self.conv2(x, edge_index, edge_attr)

        return F.sigmoid(x)
    
model = GCN()
optimizer = torch.optim.Adam(model.parameters(), lr=0.01)
loss = torch.nn.MSELoss()

model.train()
for epoch in range(200):
    loss_v = 0.0
    for data in dataset[0:5]:
        optimizer.zero_grad()
        out = model(data)
        output = loss(out, data.y)
        output.backward()
        loss_v += output.item()
        optimizer.step()
    print(loss_v / 5.0)