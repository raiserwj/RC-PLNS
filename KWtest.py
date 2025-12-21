total = 0.0

# with open("record.txt", "r", encoding="utf-8") as f:
#     for line in f:
#         line = line.strip()
#         if not line:
#             continue
#         total += float(line)
# print("Sum =", total)
import numpy as np
from scipy.stats import kruskal
import scikit_posthocs as sp
from statsmodels.stats.multitest import multipletests

def kw_dunn_holm(ours, brkga, brkga_bssf, alpha=0.05):
    """
    ours/brkga/brkga_bssf: 1D array-like, length=30 (or any n), independent runs.
    return: dict with kw p-value and Holm-adjusted pairwise p-values.
    """
    data=[]
    with open("record_.txt", "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            data.append(line)
    out_=[]
    for num in range(12):
        ours = np.asarray(data[num*30:(num+1)*30], dtype=float)
        brkga = np.asarray(data[num*30+360:(num+1)*30+360], dtype=float)
        brkga_bssf = np.asarray(data[num*30+720:(num+1)*30+720], dtype=float)

        # 1) Kruskal–Wallis
        H, p_kw = kruskal(ours, brkga, brkga_bssf)

        out = {"H": float(H), "p_kw": float(p_kw), "pairwise_raw": None, "pairwise_holm": None}

        # 2) Dunn post hoc only if KW significant (常见做法；你也可以不加这个 if)
        if p_kw < alpha:
            groups = [ours, brkga, brkga_bssf]
            labels = ["Ours", "BRKGA", "BRKGA+BSSF"]

            # Dunn raw p-values matrix (unadjusted)
            pmat = sp.posthoc_dunn(groups, p_adjust=None)
            pmat.index = labels
            pmat.columns = labels

            # Holm correction (对所有两两比较：3 pairs)
            pairs = [("Ours", "BRKGA"), ("Ours", "BRKGA+BSSF"), ("BRKGA", "BRKGA+BSSF")]
            raw_ps = np.array([pmat.loc[a, b] for a, b in pairs], dtype=float)
            _, holm_ps, _, _ = multipletests(raw_ps, alpha=alpha, method="holm")

            holm = {f"{a} vs {b}": float(p) for (a, b), p in zip(pairs, holm_ps)}
            raw = {f"{a} vs {b}": float(p) for (a, b), p in zip(pairs, raw_ps)}

            out["pairwise_raw"] = raw
            out["pairwise_holm"] = holm

        out_.append(out)
    return out_

# 示例（把下面三行替换成你的 30-run 结果）
ours = np.random.rand(30)
brkga = np.random.rand(30)
brkga_bssf = np.random.rand(30)

res = kw_dunn_holm(ours, brkga, brkga_bssf, alpha=0.05)
print(res)

print("Sum =", total)