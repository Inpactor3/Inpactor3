"""
CNN residual 1D con cabezas múltiples estilo YORO/YOLO para Inpactor3.

Entrada:  (B, 4, W)      W = 50_000
Salida:   (B, S, 4+C)    S = W // cell_size = 500
    canal 0 : objectness (logit)
    canal 1 : offset dentro de celda (sigmoide en decode)
    canal 2 : log-longitud (real)
    canal 3+: logits de clase (C = n_lineajes + 1 background)

El stride total del backbone iguala `cell_size`. Con cell_size=100 usamos
stride 100 = 2 * 5 * 5 * 2 mediante bloques con stride.
"""
from __future__ import annotations

import torch
import torch.nn as nn
import torch.nn.functional as F


class ResBlock1D(nn.Module):
    def __init__(self, c_in: int, c_out: int, stride: int = 1):
        super().__init__()
        self.conv1 = nn.Conv1d(c_in, c_out, 7, stride=stride, padding=3, bias=False)
        self.bn1 = nn.BatchNorm1d(c_out)
        self.conv2 = nn.Conv1d(c_out, c_out, 7, padding=3, bias=False)
        self.bn2 = nn.BatchNorm1d(c_out)
        self.short = (
            nn.Identity()
            if (c_in == c_out and stride == 1)
            else nn.Conv1d(c_in, c_out, 1, stride=stride, bias=False)
        )

    def forward(self, x):
        h = F.relu(self.bn1(self.conv1(x)))
        h = self.bn2(self.conv2(h))
        return F.relu(h + self.short(x))


class Yoro1D(nn.Module):
    def __init__(self, n_classes: int, cell_size: int = 100, base_ch: int = 32):
        super().__init__()
        assert cell_size == 100, "Demo pensada para cell_size=100 (stride total)."
        self.n_classes = n_classes
        self.cell_size = cell_size
        c = base_ch
        # strides que multiplican a 100: 2, 5, 5, 2
        self.stem = ResBlock1D(4, c, stride=2)
        self.b1 = ResBlock1D(c, c * 2, stride=5)
        self.b2 = ResBlock1D(c * 2, c * 4, stride=5)
        self.b3 = ResBlock1D(c * 4, c * 4, stride=2)
        self.b4 = ResBlock1D(c * 4, c * 4, stride=1)
        self.head = nn.Conv1d(c * 4, 3 + n_classes, 1)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        h = self.stem(x)
        h = self.b1(h)
        h = self.b2(h)
        h = self.b3(h)
        h = self.b4(h)
        out = self.head(h)  # (B, 3+C, S)
        return out.transpose(1, 2).contiguous()  # (B, S, 3+C)


def detection_loss(
    pred: torch.Tensor,
    target: dict,
    lambda_reg: float = 1.0,
    lambda_cls: float = 1.0,
    pos_weight: float = 200.0,
) -> tuple[torch.Tensor, dict]:
    """BCE(objectness) con pos_weight (desbalance ~250:1) + SmoothL1 + CE(clase)."""
    obj_logit = pred[..., 0]
    off = pred[..., 1]
    lenp = pred[..., 2]
    cls_logits = pred[..., 3:]

    y_obj = target["y_obj"]
    y_off = target["y_off"]
    y_len = target["y_len"]
    y_cls = target["y_cls"]

    pw = torch.tensor(pos_weight, device=pred.device)
    loss_obj = F.binary_cross_entropy_with_logits(obj_logit, y_obj, pos_weight=pw)

    pos = y_obj > 0.5
    if pos.any():
        # y_off puede caer fuera de [0,1) por elementos que empiezan antes de la celda;
        # se lo comparamos al offset "crudo" (sin sigmoide) via SmoothL1.
        loss_off = F.smooth_l1_loss(off[pos], y_off[pos])
        loss_len = F.smooth_l1_loss(lenp[pos], y_len[pos])
        loss_cls = F.cross_entropy(cls_logits[pos], y_cls[pos])
    else:
        loss_off = loss_len = loss_cls = torch.zeros((), device=pred.device)

    loss = loss_obj + lambda_reg * (loss_off + loss_len) + lambda_cls * loss_cls
    return loss, {
        "obj": loss_obj.item(),
        "off": loss_off.detach().item(),
        "len": loss_len.detach().item(),
        "cls": loss_cls.detach().item(),
    }


@torch.no_grad()
def decode(pred: torch.Tensor, cell_size: int, thr: float = 0.5):
    """Devuelve, por muestra del batch, lista de (start, end, cls_id, score)."""
    B, S, _ = pred.shape
    obj = torch.sigmoid(pred[..., 0])
    off = pred[..., 1]  # crudo, coherente con la loss SmoothL1
    lenp = pred[..., 2]
    cls = pred[..., 3:].argmax(-1)
    outs = []
    for b in range(B):
        keep = (obj[b] > thr).nonzero(as_tuple=True)[0]
        boxes = []
        for s in keep.tolist():
            cell_start = s * cell_size
            start = cell_start + off[b, s].item() * cell_size
            length = float(torch.exp(lenp[b, s]).item()) * cell_size
            boxes.append((start, start + length, int(cls[b, s].item()), float(obj[b, s])))
        outs.append(boxes)
    return outs
