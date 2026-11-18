"""Auditoria do experimento TCGA: procura vazamento e confirma se a tarefa é genuinamente fácil.

Uso: python scripts/auditoria_tcga.py   (precisa dos .gz em data/raw/tcga)
"""
import gzip
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.metrics import accuracy_score, f1_score, roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.tree import DecisionTreeClassifier

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
TCGA = ROOT / "data" / "raw" / "tcga"
RS = 42
COORTES = {"COAD": "cólon", "STAD": "estômago", "LUAD": "pulmão", "LIHC": "fígado", "BRCA": "mama", "KIRC": "rim"}

# --- recarrega a matriz SEM o filtro de variância (mesma amostragem do notebook) ---
pm = pd.read_csv(TCGA / "probeMap_450k.tsv", sep="\t", index_col=0)
auto = set(pm.index[pm["chrom"].isin({f"chr{i}" for i in range(1, 23)})])
rng = np.random.default_rng(RS)
mats, origem = {}, {}
for c in COORTES:
    path = TCGA / f"TCGA.{c}.HumanMethylation450.gz"
    with gzip.open(path, "rt") as f:
        cols = f.readline().rstrip("\n").split("\t")[1:]
    for s in cols:
        origem[s] = c
    normais = [s for s in cols if s[13:15] == "11"]
    tum = [s for s in cols if s[13:15] == "01"]
    pn = {s[:12] for s in normais}
    par = [s for s in tum if s[:12] in pn]
    resto = [s for s in tum if s[:12] not in pn]
    extra = list(rng.choice(resto, size=min(max(0, 150 - len(par)), len(resto)), replace=False))
    sel = normais + par + extra
    blocos = [ch[ch.index.isin(auto)] for ch in pd.read_csv(path, sep="\t", index_col=0, usecols=["sample"] + sel,
                                                            dtype={s: np.float32 for s in sel}, chunksize=50_000)]
    mats[c] = pd.concat(blocos)
    print("lido", c, flush=True)
M = pd.concat(mats.values(), axis=1, join="inner").T
del mats
M = M.loc[:, M.isna().mean() == 0]
print("matriz sem filtro de variância:", M.shape, flush=True)

meta = pd.DataFrame(index=M.index)
meta["paciente"] = meta.index.str[:12]
meta["tipo"] = np.where(meta.index.str[13:15] == "11", "normal", "tumor")
meta["tecido"] = meta.index.map(origem).map(COORTES)
tum = meta.index[meta.tipo == "tumor"]

# --- 1. integridade: IDs duplicados e pacientes cruzando folds ---
print("\n[1] IDs de amostra duplicados:", M.index.duplicated().sum())
cv = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=RS)
y1, g1 = meta["tipo"], meta["paciente"]
cruz = sum(len(set(g1.iloc[tr]) & set(g1.iloc[te])) for tr, te in cv.split(M, y1, g1))
print("[1] pacientes presentes em treino E teste do mesmo fold (etapa 1):", cruz)


class TopVariance(BaseEstimator, TransformerMixin):
    """Filtro de variância ajustado só no treino (versão 'paranoica' do filtro não supervisionado)."""
    def __init__(self, n=20_000):
        self.n = n
    def fit(self, X, y=None):
        self.idx_ = np.argsort(np.asarray(X).var(axis=0))[::-1][: self.n]
        return self
    def transform(self, X):
        return np.asarray(X)[:, self.idx_]


def rf_pipe(pre=None, k=500, model=None):
    steps = ([("var", pre)] if pre is not None else []) + [
        ("sel", SelectKBest(f_classif, k=k)),
        ("clf", model or RandomForestClassifier(500, class_weight="balanced", random_state=RS, n_jobs=-1))]
    return Pipeline(steps)


def avaliar(nome, X, y, g, pipe):
    p = cross_val_predict(pipe, X, y, groups=g, cv=cv, method="predict_proba")
    classes = sorted(y.unique())
    pred = np.array(classes)[p.argmax(1)]
    auc = roc_auc_score(y == classes[1], p[:, 1]) if len(classes) == 2 else roc_auc_score(y, p, multi_class="ovr")
    print(f"  {nome:55s} acc={accuracy_score(y, pred):.3f}  f1_macro={f1_score(y, pred, average='macro'):.3f}  auc={auc:.3f}", flush=True)


y2, g2 = meta.loc[tum, "tecido"], meta.loc[tum, "paciente"]
top_global = M.var().sort_values(ascending=False).index[:20_000]
aleat = M.columns[np.random.default_rng(RS).choice(M.shape[1], 20_000, replace=False)]

for etapa, X, y, g in [("ETAPA 1 (tumor vs normal)", M, y1, g1), ("ETAPA 2 (tecido de origem)", M.loc[tum], y2, g2)]:
    print(f"\n=== {etapa} ===")
    avaliar("[ref] top20k variância global (como no notebook)", X[top_global], y, g, rf_pipe())
    avaliar("[2] top20k variância ajustado SÓ no treino de cada fold", X, y, g, rf_pipe(TopVariance(20_000)))
    avaliar("[3] 20k CpGs ALEATÓRIOS", X[aleat], y, g, rf_pipe())
    y_perm = pd.Series(np.random.default_rng(RS).permutation(y.values), index=y.index)
    avaliar("[4] rótulos EMBARALHADOS (deve dar ~acaso)", X[top_global], y_perm, g, rf_pipe())
    avaliar("[5] árvore de decisão com 3 níveis (modelo trivial)", X[top_global], y, g,
            rf_pipe(model=DecisionTreeClassifier(max_depth=3, class_weight="balanced", random_state=RS)))
    avaliar("[5] 1 único CpG + árvore de 1 nível", X[top_global], y, g,
            rf_pipe(k=1, model=DecisionTreeClassifier(max_depth=1, class_weight="balanced", random_state=RS)))
