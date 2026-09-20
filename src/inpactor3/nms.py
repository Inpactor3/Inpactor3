"""1D Non-Maximum Suppression."""
from __future__ import annotations


def iou_1d(a: tuple[float, float], b: tuple[float, float]) -> float:
    s = max(a[0], b[0])
    e = min(a[1], b[1])
    inter = max(0.0, e - s)
    union = (a[1] - a[0]) + (b[1] - b[0]) - inter
    return inter / union if union > 0 else 0.0


def nms(boxes: list[tuple[float, float, int, float]], iou_thr: float = 0.3):
    """boxes = [(start, end, cls, score)]; conserva mayor score por cluster."""
    kept: list[tuple[float, float, int, float]] = []
    for b in sorted(boxes, key=lambda x: -x[3]):
        if all(iou_1d((b[0], b[1]), (k[0], k[1])) < iou_thr for k in kept):
            kept.append(b)
    return kept
