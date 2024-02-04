import numpy as np
import pandas as pd
import numpy as np
import torch.nn as nn
import torch.nn.functional as F
import scanpy as sc
import random
import torch
from sklearn.semi_supervised import LabelPropagation, LabelSpreading


class VAE(nn.Module):
    def __init__(self, input_size, output_size, h_dim=400, z_dim=20):
        super(VAE, self).__init__()
        self.output_size=output_size
        self.fc1 = nn.Linear(input_size, h_dim)
        self.batch_norm1 = torch.nn.BatchNorm1d(h_dim)
        self.fc2 = nn.Linear(h_dim, z_dim)
        self.fc3 = nn.Linear(h_dim, z_dim)
        self.fc4 = nn.Linear(z_dim, h_dim)
        self.batch_norm2 = torch.nn.BatchNorm1d(h_dim)
        self.fc5 = nn.Linear(h_dim, output_size)

    # 编码过程
    def encode(self, x):
        h = F.relu(self.fc1(x))
        h = self.batch_norm1(h)
        return self.fc2(h), self.fc3(h)

    # 随机生成隐含向量
    def reparameterize(self, mu, log_var):
        std = torch.exp(log_var / 2)
        eps = torch.randn_like(std)
        return mu + eps * std

    # 解码过程
    def decode(self, z):
        h = F.relu(self.fc4(z))
        h = self.batch_norm2(h)
        return torch.sigmoid(self.fc5(h))

    # 整个前向传播过程：编码-》解码
    def forward(self, x):
        mu, log_var = self.encode(x)
        z = self.reparameterize(mu, log_var)
        x_reconst = self.decode(z)
        return x_reconst, mu, log_var


def Gen_TrainSet(h5adFile='resources/V1-Mouse-Brain.h5ad'):
    adata = sc.read_h5ad(h5adFile)

    # generate leiden cluster for resolution=10
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, flavor='seurat', inplace=True)
    sc.pp.pca(adata, n_comps=50, use_highly_variable=True, svd_solver='arpack')
    sc.pp.neighbors(adata)
    # use umap and leiden for clustering
    sc.tl.umap(adata)

    candidates = pd.DataFrame(index=adata.obs_names)
    print("For small resolutions")
    can_num = 0
    for res in range(1, 10, 2):
        lores = res*0.1
        sc.tl.leiden(adata, key_added='clusters', resolution=lores)
        tmp_can = np.array(list(map(lambda x: np.array(adata.obs['clusters'] == x), adata.obs['clusters'].unique()))).T
        tmp_can = pd.DataFrame(tmp_can, index=adata.obs_names)
        candidates = pd.concat([candidates, tmp_can], axis=1, ignore_index=True)
        uni_num = len(adata.obs['clusters'].unique())
        print("unique candidates:{}".format(uni_num))
        can_num += uni_num
    print("For large resolutions")
    for res in [1, 2, 5, 10]:
        sc.tl.leiden(adata, key_added='clusters', resolution=res)
        tmp_can = np.array(list(map(lambda x: np.array(adata.obs['clusters'] == x), adata.obs['clusters'].unique()))).T
        tmp_can = pd.DataFrame(tmp_can, index=adata.obs_names)
        candidates = pd.concat([candidates, tmp_can], axis=1, ignore_index=True)
        uni_num = len(adata.obs['clusters'].unique())
        print("unique candidates:{}".format(uni_num))
        can_num += uni_num
    c = list(candidates.columns)
    random.shuffle(c)
    candidates = candidates.reindex(columns=c)
    print("All candidates Number: {}".format(can_num))
    return adata, candidates


def Gen_maskSet(candidate: pd.DataFrame, errRate=0.20):
    sele_can = candidate[candidate==True]
    cell_len = len(sele_can)
    mask_can = candidate.copy()
    errArray = random.sample(range(cell_len), int(cell_len*errRate))
    for cell in errArray:
        print(sele_can.index[cell])
        mask_can.loc[sele_can.index[cell]] = not mask_can.loc[sele_can.index[cell]]
    return mask_can


def train_VAE(candidates, h_dim=400, z_dim=40, num_epochs=50, learning_rate=1e-3, batch_size=128):
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    model = VAE(input_size=len(candidates), output_size=len(candidates), h_dim=h_dim, z_dim=z_dim).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)

    for epoch in range(num_epochs):
        kl_divs = []
        mask_x = torch.from_numpy(np.array([Gen_maskSet(candidates.iloc[:, i], errRate=0.05) for i in range(len(candidates.columns))], dtype=np.float32))
        print(mask_x.shape)
        origin_x = torch.from_numpy(np.array(candidates, dtype=np.float32).T)
        mask_x = mask_x.to(device)

        reconst_x, mu, log_var = model(mask_x)
        reconst_loss = F.binary_cross_entropy(reconst_x, origin_x, reduction='mean')
        kl_div = - 0.5 * torch.sum(1 + log_var - mu.pow(2) - log_var.exp())

        loss = reconst_loss + kl_div
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        kl_divs.append(kl_div)

        print("Epoch[{}/{}] Reconst Loss: {}, KL Div: {}"
                  .format(epoch, num_epochs, reconst_loss.item(), kl_div.item()))

    with torch.no_grad():

        test_x = torch.from_numpy(np.array([Gen_maskSet(candidates.iloc[:, i], errRate=0.05) for i in range(len(candidates.columns))], dtype=np.float32))
        print(test_x)
        out, _, _ = model(test_x)
        npout = out.cpu().numpy().T
        npout = pd.DataFrame(npout)
        prob = candidates.sum(axis=0) / len(candidates)
        for i in range(len(npout.columns)):
            sorted_arr = np.sort(np.array(npout[i]))
            index = int((len(sorted_arr) - 1) * (1-prob.iloc[i])) + 1
            flag = sorted_arr[index]
            print(flag)
            npout[i] = npout[i].apply(lambda x: 1 if x >= flag else 0)
        return npout


def train_LPA(candidate, adata, use_model=LabelPropagation, errRate=0.05):
    X = adata.obsm['X_pca']
    y = Gen_maskSet(candidate, errRate)
    model = use_model().fit(X, y)
    y_pred = model.predict(X)

    def acc(y_true, y_pred):
        return np.sum(np.equal(y_true, y_pred))/len(y_true)
    print("acc:{}".format(acc(candidate, y_pred)))


def do():
    candidates = Gen_TrainSet()
    out = train_VAE(candidates)



















