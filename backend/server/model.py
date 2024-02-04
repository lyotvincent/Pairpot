from typing import Dict, Tuple
from tqdm import tqdm
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import DataLoader
from torchvision import models, transforms
from torchvision.datasets import MNIST
from torchvision.utils import save_image, make_grid
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from sklearn.semi_supervised import LabelPropagation, LabelSpreading
# Single Cell Referenecs Down Sampling Through

device = "mps"
print(f"Using {device} device")
# VAE model

class VAE(nn.Module):
    def __init__(self, input_size, output_size, h_dim=400, z_dim=20):
        super(VAE, self).__init__()
        self.fc1 = nn.Linear(input_size, h_dim)
        self.fc2 = nn.Linear(h_dim, z_dim)  # 均值 向量
        self.fc3 = nn.Linear(h_dim, z_dim)  # 保准方差 向量
        self.fc4 = nn.Linear(z_dim, h_dim)
        self.fc5 = nn.Linear(h_dim, output_size)

    # 编码过程
    def encode(self, x):
        h = F.relu(self.fc1(x))
        return self.fc2(h), self.fc3(h)

    # 随机生成隐含向量
    def reparameterize(self, mu, log_var):
        std = torch.exp(log_var / 2)
        eps = torch.randn_like(std)
        return mu + eps * std

    # 解码过程
    def decode(self, z):
        h = F.relu(self.fc4(z))
        return torch.sigmoid(self.fc5(h))

    # 整个前向传播过程：编码-》解码
    def forward(self, x):
        mu, log_var = self.encode(x)
        z = self.reparameterize(mu, log_var)
        x_reconst = self.decode(z)
        return x_reconst, mu, log_var


def conv1d_5(inplanes, outplanes, stride=1):

    return nn.Conv1d(inplanes, outplanes, kernel_size=5, stride=stride,
                     padding=2, bias=False)


class Block(nn.Module):
    def __init__(self, inplanes, planes, stride=1, downsample=None, bn=False):
        super(Block, self).__init__()
        self.bn = bn
        self.conv1 = conv1d_5(inplanes, planes, stride)
        self.bn1 = nn.BatchNorm1d(planes)
        self.relu = nn.ReLU(inplace=False)
        self.conv2 = conv1d_5(planes, planes)
        self.bn2 = nn.BatchNorm1d(planes)
        self.downsample = downsample
        self.stride = stride

    def forward(self, x):
        residual = x

        out = self.conv1(x)
        if self.bn:
            out = self.bn1(out)
        out = self.relu(out)

        out = self.conv2(out)
        if self.bn:
            out = self.bn2(out)
        out = self.relu(out)

        if self.downsample is not None:
            residual = self.downsample(x)

        out += residual
        out = self.relu(out)

        return out


class Decoder_block(nn.Module):

    def __init__(self, inplanes, outplanes, kernel_size=5, stride=5):
        super(Decoder_block, self).__init__()
        self.upsample = nn.ConvTranspose1d(in_channels=inplanes, out_channels=outplanes,
                           padding=0, kernel_size=kernel_size, stride=stride, bias=False)
        self.conv1 = conv1d_5(inplanes, outplanes)
        self.relu = nn.ReLU()
        self.conv2 = conv1d_5(outplanes, outplanes)

    def forward(self, x1, x2):
        x1 = self.upsample(x1)
        out = torch.cat((x1, x2), dim=1)

        out = self.conv1(out)
        out = self.relu(out)

        out = self.conv2(out)
        out = self.relu(out)

        return out


class Uresnet(nn.Module):

    def __init__(self, inplanes=5, layers=[2,2,2,2]):
        super(Uresnet, self).__init__()
        self.inplanes = inplanes
        self.encoder1 = self._make_encoder(Block, 32, layers[0], 5)
        self.encoder2 = self._make_encoder(Block, 64, layers[1], 5)
        self.encoder3 = self._make_encoder(Block, 128, layers[2], 5)
        self.encoder4 = self._make_encoder(Block, 256, layers[3], 4)
        self.decoder3 = Decoder_block(256, 128, stride=4, kernel_size=4)
        self.decoder2 = Decoder_block(128, 64)
        self.decoder1 = Decoder_block(64, 32)
        # 如果要修改输出的通道，在这里修改
        self.conv1x1 = nn.ConvTranspose1d(32, 1, kernel_size=5, stride=5, bias=False)

    def _make_encoder(self, block, planes, blocks, stride=1):
        downsample = None
        if self.inplanes != planes or stride != 1:
            downsample = nn.Conv1d(self.inplanes, planes, stride=stride, kernel_size=1, bias=False)
        layers = []
        layers.append(block(self.inplanes, planes, stride, downsample))
        self.inplanes = planes
        for i in range(1, blocks):
            layers.append(block(self.inplanes, planes))
        return nn.Sequential(*layers)


    def forward(self, x):
        down1 = self.encoder1(x)
        down2 = self.encoder2(down1)
        down3 = self.encoder3(down2)
        down4 = self.encoder4(down3)

        up3 = self.decoder3(down4, down3)
        up2 = self.decoder2(up3, down2)
        up1 = self.decoder1(up2, down1)
        out = self.conv1x1(up1)

        return out


class ResidualConvBlock(nn.Module):
    def __init__(
            self, in_channels: int, out_channels: int, is_res: bool = False
    ) -> None:
        super().__init__()
        '''
        standard ResNet style convolutional block
        '''
        self.same_channels = in_channels == out_channels
        self.is_res = is_res
        self.conv1 = nn.Sequential(
            nn.Conv1d(in_channels, out_channels, 3, 1, 1),
            nn.BatchNorm1d(out_channels),
            nn.GELU(),
        )
        self.conv2 = nn.Sequential(
            nn.Conv1d(out_channels, out_channels, 3, 1, 1),
            nn.BatchNorm1d(out_channels),
            nn.GELU(),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        if self.is_res:
            x1 = self.conv1(x)
            x2 = self.conv2(x1)

            if self.same_channels:
                out = x + x2
            else:
                out = x1 + x2
            return out / 1.414
        else:
            x1 = self.conv1(x)
            x2 = self.conv2(x1)
            return x2


class UnetDown(nn.Module):
    def __init__(self, in_channels, out_channels):
        super(UnetDown, self).__init__()
        '''
        process and downscale the image feature maps
        '''
        layers = [ResidualConvBlock(in_channels, out_channels), nn.MaxPool1d(2)]
        self.model = nn.Sequential(*[ResidualConvBlock(in_channels, out_channels), nn.MaxPool1d(2)])

    def forward(self, x):
        return self.model(x)


class UnetUp(nn.Module):
    def __init__(self, in_channels, out_channels):
        super(UnetUp, self).__init__()
        '''
        process and upscale the image feature maps
        '''
        layers = [
            nn.ConvTranspose1d(in_channels, out_channels, 4, 2, padding=1),
            ResidualConvBlock(out_channels, out_channels),
            ResidualConvBlock(out_channels, out_channels),
        ]
        self.model = nn.Sequential(*layers)

    def forward(self, x, skip):
        x = torch.cat((x, skip), 1)
        x = self.model(x)
        return x


class EmbedFC(nn.Module):
    def __init__(self, input_dim, emb_dim):
        super(EmbedFC, self).__init__()
        '''
        generic one layer FC NN for embedding things  
        '''
        self.input_dim = input_dim
        layers = [
            nn.Linear(input_dim, emb_dim),
            nn.GELU(),
            nn.Linear(emb_dim, emb_dim),
        ]
        self.model = nn.Sequential(*layers)

    def forward(self, x):
        x = x.view(-1, self.input_dim)
        return self.model(x)


class ContextUnet(nn.Module):
    def __init__(self, in_channels, n_feat=256, n_classes=10):
        super(ContextUnet, self).__init__()

        self.in_channels = in_channels
        self.n_feat = n_feat
        self.n_classes = n_classes

        self.init_conv = ResidualConvBlock(in_channels, n_feat, is_res=True)

        self.down1 = UnetDown(n_feat, n_feat)
        self.down2 = UnetDown(n_feat, 2 * n_feat)

        self.to_vec = nn.Sequential(nn.AvgPool1d(7), nn.GELU())

        self.timeembed1 = EmbedFC(1, 2 * n_feat)
        self.timeembed2 = EmbedFC(1, 1 * n_feat)
        self.contextembed1 = EmbedFC(n_classes, 2 * n_feat)
        self.contextembed2 = EmbedFC(n_classes, 1 * n_feat)

        self.up0 = nn.Sequential(
            # nn.ConvTranspose2d(6 * n_feat, 2 * n_feat, 7, 7), # when concat temb and cemb end up w 6*n_feat
            nn.ConvTranspose1d(2 * n_feat, 2 * n_feat, 7, 7),  # otherwise just have 2*n_feat
            nn.GroupNorm(8, 2 * n_feat),
            nn.ReLU(),
        )

        self.up1 = UnetUp(4 * n_feat, n_feat)
        self.up2 = UnetUp(2 * n_feat, n_feat)
        self.out = nn.Sequential(
            nn.Conv1d(2 * n_feat, n_feat, 3, 1, 1),
            nn.GroupNorm(8, n_feat),
            nn.ReLU(),
            nn.Conv1d(n_feat, self.in_channels, 3, 1, 1),
        )

    def forward(self, x, c, t, context_mask):
        x = self.init_conv(x)
        # print("x:",x.shape)
        down1 = self.down1(x)
        # print("down1:", down1.shape)
        down2 = self.down2(down1)
        # print("down2:", down2.shape)
        hiddenvec = self.to_vec(down2)
        # print("hiddenvec:", hiddenvec.shape)

        c = nn.functional.one_hot(c, num_classes=self.n_classes).type(torch.float)

        context_mask = context_mask[:, None]
        context_mask = context_mask.repeat(1, self.n_classes)
        context_mask = (-1 * (1 - context_mask))  # need to flip 0 <-> 1
        c = c * context_mask

        cemb1 = self.contextembed1(c).view(-1, self.n_feat * 2, 1)
        # print("cemb1", cemb1.shape)
        temb1 = self.timeembed1(t).view(-1, self.n_feat * 2, 1)
        # print("temb1", temb1.shape)
        cemb2 = self.contextembed2(c).view(-1, self.n_feat, 1)
        # print("cemb2", cemb2.shape)
        temb2 = self.timeembed2(t).view(-1, self.n_feat, 1)
        # print("temb2", temb2.shape)

        up1 = self.up0(hiddenvec)
        # print("up1", up1.shape)
        # up2 = self.up1(up1, down2) # if want to avoid add and multiply embeddings
        up2 = self.up1(cemb1 * up1 + temb1, down2)  # add and multiply embeddings
        # print("up2", up2.shape)
        up3 = self.up2(cemb2 * up2 + temb2, down1)
        out = self.out(torch.cat((up3, x), 1))
        return out


def ddpm_schedules(beta1, beta2, T):
    assert beta1 < beta2 < 1.0, "beta1 and beta2 must be in (0, 1)"

    beta_t = (beta2 - beta1) * torch.arange(0, T + 1, dtype=torch.float32) / T + beta1
    sqrt_beta_t = torch.sqrt(beta_t)
    alpha_t = 1 - beta_t
    log_alpha_t = torch.log(alpha_t)
    alphabar_t = torch.cumsum(log_alpha_t, dim=0).exp()

    sqrtab = torch.sqrt(alphabar_t)
    oneover_sqrta = 1 / torch.sqrt(alpha_t)

    sqrtmab = torch.sqrt(1 - alphabar_t)
    mab_over_sqrtmab_inv = (1 - alpha_t) / sqrtmab

    return {
        "alpha_t": alpha_t,  # \alpha_t
        "oneover_sqrta": oneover_sqrta,  # 1/\sqrt{\alpha_t}
        "sqrt_beta_t": sqrt_beta_t,  # \sqrt{\beta_t}
        "alphabar_t": alphabar_t,  # \bar{\alpha_t}
        "sqrtab": sqrtab,  # \sqrt{\bar{\alpha_t}}
        "sqrtmab": sqrtmab,  # \sqrt{1-\bar{\alpha_t}}
        "mab_over_sqrtmab": mab_over_sqrtmab_inv,  # (1-\alpha_t)/\sqrt{1-\bar{\alpha_t}}
    }


class DDPM(nn.Module):
    def __init__(self, nn_model, betas, n_T, device, drop_prob=0.1):
        super(DDPM, self).__init__()
        self.nn_model = nn_model.to(device)

        # register_buffer allows accessing dictionary produced by ddpm_schedules
        # e.g. can access self.sqrtab later
        for k, v in ddpm_schedules(betas[0], betas[1], n_T).items():
            self.register_buffer(k, v)

        self.n_T = n_T
        self.device = device
        self.drop_prob = drop_prob
        self.loss_mse = nn.MSELoss()

    def forward(self, x, c):

        _ts = torch.randint(1, self.n_T, (x.shape[0],)).to(self.device)  # t ~ Uniform(0, n_T)
        noise = torch.randn_like(x)  # eps ~ N(0, 1)
        x_t = (
                self.sqrtab[_ts, None, None] * x
                + self.sqrtmab[_ts, None, None] * noise
        )
        context_mask = torch.bernoulli(torch.zeros_like(c) + self.drop_prob).to(self.device)
        # print(context_mask.shape)
        return self.loss_mse(noise, self.nn_model(x_t, c, _ts / self.n_T, context_mask))

    def sample(self, n_sample, size, device, guide_w=0.0):

        x_i = torch.randn(n_sample, *size).to(device)  # x_T ~ N(0, 1), sample initial noise
        c_i = torch.arange(0, 10).to(device)  # context for us just cycles throught the mnist labels
        c_i = c_i.repeat(int(n_sample / c_i.shape[0]))

        # don't drop context at test time
        context_mask = torch.zeros_like(c_i).to(device)

        # double the batch
        c_i = c_i.repeat(2)
        context_mask = context_mask.repeat(2)
        context_mask[n_sample:] = 1.  # makes second half of batch context free

        x_i_store = []  # keep track of generated steps in case want to plot something
        print()
        for i in range(self.n_T, 0, -1):
            print(f'sampling timestep {i}', end='\r')
            t_is = torch.tensor([i / self.n_T]).to(device)
            t_is = t_is.repeat(n_sample, 1, 1)

            # double batch
            x_i = x_i.repeat(2, 1, 1)
            t_is = t_is.repeat(2, 1, 1)

            z = torch.randn(n_sample, *size).to(device) if i > 1 else 0

            eps = self.nn_model(x_i, c_i, t_is, context_mask)
            eps1 = eps[:n_sample]
            eps2 = eps[n_sample:]
            eps = (1 + guide_w) * eps1 - guide_w * eps2
            x_i = x_i[:n_sample]
            x_i = (
                    self.oneover_sqrta[i] * (x_i - eps * self.mab_over_sqrtmab[i])
                    + self.sqrt_beta_t[i] * z
            )
            if i % 20 == 0 or i == self.n_T or i < 8:
                x_i_store.append(x_i.detach().cpu().numpy())

        x_i_store = np.array(x_i_store)
        return x_i, x_i_store


n_epoch = 5
batch_size = 64
n_T = 400
device = torch.device('cuda' if torch.cuda.is_available() else 'mps')  # 设备
n_classes = 10
n_feat = 128
lrate = 1e-4

ws_test = [0.0, 0.5, 2.0]
ddpm = DDPM(nn_model=ContextUnet(in_channels=1, n_feat=n_feat, n_classes=n_classes), betas=(1e-4, 0.02), n_T=n_T,
            device=device, drop_prob=0.1)
ddpm = ddpm.to(device)

tf = transforms.Compose([transforms.ToTensor()])  # mnist is already normalised 0 to 1

dataset = MNIST("MNIST", train=True, download=True, transform=tf)
# dataloader = DataLoader(dataset, batch_size=batch_size, shuffle=True, num_workers=5)
dataloader = DataLoader(dataset, batch_size=batch_size, shuffle=True)
optim = torch.optim.Adam(ddpm.parameters(), lr=lrate)

pp = None
for ep in range(n_epoch):
    print(f'epoch {ep}')
    pp = ddpm.train()
    # linear lrate decay
    optim.param_groups[0]['lr'] = lrate*(1-ep/n_epoch)

    pbar = tqdm(dataloader)
    loss_ema = None
    for x, c in pbar:
        x = np.array([tensor.view(tensor.size(0), -1) for tensor in x])
        x = torch.tensor(x)
        optim.zero_grad()
        x = x.to(device)
        c = c.to(device)
        loss = ddpm(x, c)
        loss.backward()
        if loss_ema is None:
            loss_ema = loss.item()
        else:
            loss_ema = 0.95 * loss_ema + 0.05 * loss.item()
        pbar.set_description(f"loss: {loss_ema:.4f}")
        optim.step()


ddpm = ddpm.eval()
with torch.no_grad():
    n_sample = 4*n_classes
    for w_i, w in enumerate(ws_test):
        x_gen, x_gen_store = ddpm.sample(n_sample, (1, 784), device, guide_w=w)
        print("output shape:",x_gen.shape)
        x_gen = x_gen.view(-1, 1, 28, 28)
        print("2d shape:", x_gen.shape)
        x_real = torch.Tensor(x_gen.shape).to(device)
        # for k in range(n_classes):
        #     for j in range(int(n_sample/n_classes)):
        #         try:
        #             idx = torch.squeeze((c == k).nonzero())[j]
        #         except:
        #             idx = 0
        #         x_real[k+(j*n_classes)] = x[idx]

        x_all = x_gen
        print(x_all)
        grid = make_grid(x_all*-1 + 1, nrow=10).cpu().detach().numpy().transpose(1, 2, 0)
        plt.imshow(grid)
        plt.show()

